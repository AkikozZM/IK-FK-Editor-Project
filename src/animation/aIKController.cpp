#include "aIKController.h"
#include "aActor.h"
#include <Eigen\Dense>
#include "aMatrix.h"

#pragma warning (disable : 4018)

int IKController::gIKmaxIterations = 5;
double IKController::gIKEpsilon = 0.1;
double lastAngle = 0;

// AIKchain class functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////
AIKchain::AIKchain()
{
	mWeight0 = 0.1;
}

AIKchain::~AIKchain()
{

}

AJoint* AIKchain::getJoint(int index) 
{ 
	return mChain[index]; 
}

void AIKchain::setJoint(int index, AJoint* pJoint) 
{ 
	mChain[index] = pJoint; 
}

double AIKchain::getWeight(int index) 
{ 
	return mWeights[index]; 
}

void AIKchain::setWeight(int index, double weight) 
{ 
	mWeights[index] = weight; 
}

int AIKchain::getSize() 
{ 
	return mChain.size(); 
}

std::vector<AJoint*>& AIKchain::getChain() 
{ 
	return mChain; 
}

std::vector<double>& AIKchain::getWeights() 
{ 
	return mWeights; 
}

void AIKchain::setChain(std::vector<AJoint*> chain) 
{
	mChain = chain; 
}

void AIKchain::setWeights(std::vector<double> weights) 
{ 
	mWeights = weights; 
}

// AIKController class functions
/////////////////////////////////////////////////////////////////////////////////////////////////////////

IKController::IKController()
{
	m_pActor = NULL;
	m_pSkeleton = NULL;
	mvalidLimbIKchains = false;
	mvalidCCDIKchains = false;

	// Limb IK
	m_pEndJoint = NULL;
	m_pMiddleJoint = NULL;
	m_pBaseJoint = NULL;
	m_rotationAxis = vec3(0.0, 1.0, 0.0);

	ATransform desiredTarget = ATransform();
	mTarget0.setLocal2Parent(desiredTarget);  // target associated with end joint
	mTarget1.setLocal2Parent(desiredTarget);  // optional target associated with middle joint - used to specify rotation of middle joint about end/base axis
	mTarget0.setLocal2Global(desiredTarget);
	mTarget1.setLocal2Global(desiredTarget);

	//CCD IK
	mWeight0 = 0.1;  // default joint rotation weight value

}

IKController::~IKController()
{
}

ASkeleton* IKController::getSkeleton()
{
	return m_pSkeleton;
}

const ASkeleton* IKController::getSkeleton() const
{
	return m_pSkeleton;
}

ASkeleton* IKController::getIKSkeleton()
{
	return &mIKSkeleton;
}

const ASkeleton* IKController::getIKSkeleton() const
{
	return &mIKSkeleton;
}

AActor* IKController::getActor()
{
	return m_pActor;
}

void IKController::setActor(AActor* actor)

{
	m_pActor = actor;
	m_pSkeleton = m_pActor->getSkeleton();
}


AIKchain IKController::createIKchain(int endJointID, int desiredChainSize, ASkeleton* pSkeleton)
{
	// TODO: given the end joint ID and the desired length of the IK chain, 
	// add the corresponding skeleton joint pointers to the AIKChain "chain" data member starting with the end joint
	// desiredChainSize = -1 should create an IK chain of maximum length (where the last chain joint is the joint before the root joint)
	// also add weight values to the associated AIKChain "weights" data member which can be used in a CCD IK implemention
	AIKchain ikChain;
	AJoint* currentJoint = pSkeleton->getJointByID(endJointID);
	int chainLength = 0;
	// Traverse the parent Joints
	while (currentJoint != nullptr && (chainLength < desiredChainSize || desiredChainSize == -1)) {
		ikChain.getChain().push_back(currentJoint);
		ikChain.getWeights().push_back(0.1);
		if (desiredChainSize == -1 && currentJoint->getParent() == nullptr) {
			break;
		}
		currentJoint = currentJoint->getParent();
		chainLength++;
	}
	
	return ikChain;
}



bool IKController::IKSolver_Limb(int endJointID, const ATarget& target)
{
	// Implements the analytic/geometric IK method assuming a three joint limb  

	// copy transforms from base skeleton
	mIKSkeleton.copyTransforms(m_pSkeleton);

	if (!mvalidLimbIKchains)
	{
		mvalidLimbIKchains = createLimbIKchains();
		if (!mvalidLimbIKchains) { return false; }
	}

	vec3 desiredRootPosition;

	// The joint IDs are not const, so we cannot use switch here
	if (endJointID == mLhandID)
	{
		mLhandTarget = target;
		computeLimbIK(mLhandTarget, mLhandIKchain, -axisY, &mIKSkeleton);
	}
	else if (endJointID == mRhandID)
	{
		mRhandTarget = target;
		computeLimbIK(mRhandTarget, mRhandIKchain, axisY, &mIKSkeleton);
	}
	else if (endJointID == mLfootID)
	{
		mLfootTarget = target;
		computeLimbIK(mLfootTarget, mLfootIKchain, axisX, &mIKSkeleton);
	}
	else if (endJointID == mRfootID)
	{
		mRfootTarget = target;
		computeLimbIK(mRfootTarget, mRfootIKchain, axisX, &mIKSkeleton);
	}
	else if (endJointID == mRootID)
	{
		desiredRootPosition = target.getGlobalTranslation();
		mIKSkeleton.getJointByID(mRootID)->setLocalTranslation(desiredRootPosition);
		mIKSkeleton.update();
		computeLimbIK(mLhandTarget, mLhandIKchain, -axisY, &mIKSkeleton);
		computeLimbIK(mRhandTarget, mRhandIKchain, axisY, &mIKSkeleton);
		computeLimbIK(mLfootTarget, mLfootIKchain, axisX, &mIKSkeleton);
		computeLimbIK(mRfootTarget, mRfootIKchain, axisX, &mIKSkeleton);
	}
	else
	{
		mIKchain = createIKchain(endJointID, 3, &mIKSkeleton);
		computeLimbIK(target, mIKchain, axisY, &mIKSkeleton);
	}

	// update IK Skeleton transforms
	mIKSkeleton.update();

	// copy IK skeleton transforms to main skeleton
	m_pSkeleton->copyTransforms(&mIKSkeleton);

	return true;
}



int IKController::createLimbIKchains()
{
	bool validChains = false;
	int desiredChainSize = 3;

	// create IK chains for Lhand, Rhand, Lfoot and Rfoot 
	mLhandIKchain = createIKchain(mLhandID, desiredChainSize, &mIKSkeleton);
	mRhandIKchain = createIKchain(mRhandID, desiredChainSize, &mIKSkeleton);
	mLfootIKchain = createIKchain(mLfootID, desiredChainSize, &mIKSkeleton);
	mRfootIKchain = createIKchain(mRfootID, desiredChainSize, &mIKSkeleton);
	
	if (mLhandIKchain.getSize() == 3 && mRhandIKchain.getSize() == 3 && mLfootIKchain.getSize() == 3 && mRfootIKchain.getSize() == 3)
	{
		validChains = true;
		
		// initalize end joint target transforms for Lhand, Rhand, Lfoot and Rfoot based on current position and orientation of joints
		mIKSkeleton.copyTransforms(m_pSkeleton);
		mLhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mLhandID)->getLocal2Global());
		mRhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mRhandID)->getLocal2Global());
		mLfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mLfootID)->getLocal2Global());
		mRfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mRfootID)->getLocal2Global());
	}

	return validChains;
}


int IKController::computeLimbIK(ATarget target, AIKchain& IKchain, const vec3 midJointAxis, ASkeleton* pIKSkeleton)
{
	// TODO: Implement the analytic/geometric IK method assuming a three joint limb  
	// The actual position of the end joint should match the target position within some episilon error 
	// the variable "midJointAxis" contains the rotation axis for the middle joint
	if (IKchain.getSize() != 3) {
		return false;
	}
	// basejoint = shoulder, endjoint = hand
	AJoint* baseJoint = IKchain.getJoint(2);
	AJoint* midJoint = IKchain.getJoint(1);
	AJoint* endJoint = IKchain.getJoint(0);

	vec3 basePos = baseJoint->getGlobalTranslation();
	vec3 midPos = midJoint->getGlobalTranslation();
	vec3 endPos = endJoint->getGlobalTranslation();
	vec3 targetPos = target.getGlobalTranslation();
	// base to mid
	double length1 = (midPos - basePos).Length();
	// mid to end
	double length2 = (endPos - midPos).Length();
	// base to target
	double length3 = (targetPos - basePos).Length();

	double cosAngle1 = (pow(length1, 2) + pow(length2, 2) - pow(length3, 2)) / (2 * length1 * length2);
	cosAngle1 = std::max(-1.0, std::min(1.0, cosAngle1));
	double angle1 = acos(cosAngle1);
	angle1 = M_PI - angle1;
	mat3 rotation;
	rotation.FromAxisAngle(midJointAxis, angle1);

	// Apply the rotation to the base joint
	midJoint->setLocalRotation(rotation);
	pIKSkeleton->update();

	// Recalculate positions after the rotation
	vec3 rd = (midPos - basePos).Normalize();
	vec3 dir = (targetPos - basePos).Normalize();
	vec3 rotationAxis = rd.Cross(dir).Normalize();
	double cosAngle2 = Dot(rd, dir);
	double rotationAngle = acos(cosAngle2);

	quat rot2;
	rot2.FromAxisAngle(rotationAxis, rotationAngle);

	// Apply the rotation to the mid joint
	baseJoint->setLocalRotation(baseJoint->getLocalRotation() * rot2.ToRotation());
	pIKSkeleton->update();

	return true;
}

bool IKController::IKSolver_CCD(int endJointID, const ATarget& target)
{
	// Implements the CCD IK method assuming a three joint limb 

	bool validChains = false;

	if (!mvalidCCDIKchains)
	{
		mvalidCCDIKchains = createCCDIKchains();
		assert(mvalidCCDIKchains);
	}

	// copy transforms from base skeleton
	mIKSkeleton.copyTransforms(m_pSkeleton);

	vec3 desiredRootPosition;


	if (endJointID == mLhandID)
	{
		mLhandTarget = target;
		computeCCDIK(mLhandTarget, mLhandIKchain, &mIKSkeleton);
	}
	else if (endJointID == mRhandID)
	{
		mRhandTarget = target;
		computeCCDIK(mRhandTarget, mRhandIKchain, &mIKSkeleton);
	}
	else if (endJointID == mLfootID)
	{
		mLfootTarget = target;
		computeCCDIK(mLfootTarget, mLfootIKchain, &mIKSkeleton);
	}
	else if (endJointID == mRfootID)
	{
		mRfootTarget = target;
		computeCCDIK(mRfootTarget, mRfootIKchain, &mIKSkeleton);
	}
	else if (endJointID == mRootID)
	{
		desiredRootPosition = target.getGlobalTranslation();
		mIKSkeleton.getJointByID(mRootID)->setLocalTranslation(desiredRootPosition);
		mIKSkeleton.update();
		computeCCDIK(mLhandTarget, mLhandIKchain, &mIKSkeleton);
		computeCCDIK(mRhandTarget, mRhandIKchain, &mIKSkeleton);
		computeCCDIK(mLfootTarget, mLfootIKchain, &mIKSkeleton);
		computeCCDIK(mRfootTarget, mRfootIKchain, &mIKSkeleton);
	}
	else
	{
		mIKchain = createIKchain(endJointID, -1, &mIKSkeleton);
		computeCCDIK(target, mIKchain, &mIKSkeleton);
	}

	// update IK Skeleton transforms
	mIKSkeleton.update();

	// copy IK skeleton transforms to main skeleton
	m_pSkeleton->copyTransforms(&mIKSkeleton);

	return true;
}

int IKController::createCCDIKchains()
{
	bool validChains = false;

	int desiredChainSize = -1;  // default of -1 creates IK chain of maximum length from end joint to child joint of root


	// create IK chains for Lhand, Rhand, Lfoot and Rfoot 
	mLhandIKchain = createIKchain(mLhandID, desiredChainSize, &mIKSkeleton);
	mRhandIKchain = createIKchain(mRhandID, desiredChainSize, &mIKSkeleton);
	mLfootIKchain = createIKchain(mLfootID, desiredChainSize, &mIKSkeleton);
	mRfootIKchain = createIKchain(mRfootID, desiredChainSize, &mIKSkeleton);

	if (mLhandIKchain.getSize() > 1 && mRhandIKchain.getSize() > 1 && mLfootIKchain.getSize() > 1 && mRfootIKchain.getSize() > 1)
	{
		validChains = true;

		// initalize end joint target transforms for Lhand, Rhand, Lfoot and Rfoot based on current position and orientation of joints
		mIKSkeleton.copyTransforms(m_pSkeleton);
		mLhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mLhandID)->getLocal2Global());
		mRhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mRhandID)->getLocal2Global());
		mLfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mLfootID)->getLocal2Global());
		mRfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mRfootID)->getLocal2Global());
	}

	return validChains;
}


int IKController::computeCCDIK(ATarget target, AIKchain& IKchain, ASkeleton* pIKSkeleton)
{

	// TODO: Implement CCD IK  
	// The actual position of the end joint should match the desiredEndPos within some episilon error 
	vec3 targetPos = target.getGlobalTranslation();
	AJoint* endJoint = IKchain.getJoint(0);
	vec3 endPos = endJoint->getGlobalTranslation();
	int size = IKchain.getSize();
	// Check if we're already close enough
	if ((endPos - targetPos).Length() < gIKEpsilon)
	{
		return true;
	}

	// 1. compute axis and angle for a joint in the IK chain (distal to proximal) in global coordinates
	for (int i = 0; i < size; i++) {
		AJoint* joint = IKchain.getJoint(i);
		vec3 jointPos = joint->getGlobalTranslation();
		endPos = endJoint->getGlobalTranslation();
		vec3 toEnd = (endPos - jointPos).Normalize();
		vec3 toTarget = (targetPos - jointPos).Normalize();
		vec3 rotationAxis = toEnd.Cross(toTarget).Normalize();
		double cosAngle = Dot(toEnd, toTarget);
		double angle = acos(std::min(std::max(cosAngle, -1.0), 1.0));
		// 2. once you have the desired axis and angle, convert axis to local joint coords 
		vec3 localAxis = joint->getGlobalRotation().Inverse() * rotationAxis;
		// 3. compute desired change to local rotation matrix
		quat deltaRot;
		deltaRot.FromAxisAngle(localAxis, angle);
		mat3 deltaRotMat;
		deltaRotMat.FromQuaternion(deltaRot);
		// 4. set local rotation matrix to new value
		mat3 newLocalRot = joint->getLocalRotation() * deltaRotMat;
		joint->setLocalRotation(newLocalRot);
		// 5. update transforms for joint and all children
		pIKSkeleton->update();
		// Check if we've reached the target
		endPos = endJoint->getGlobalTranslation();
		if ((endPos - targetPos).Length() < gIKEpsilon)
		{
			return true;
		}
	}
	return false;
}

int IKController::createPseudoInvIKchains()
{
	bool validChains = false;

	int desiredChainSize = -1;  // default of -1 creates IK chain of maximum length from end joint to child joint of root


	// create IK chains for Lhand, Rhand, Lfoot and Rfoot 
	mLhandIKchain = createIKchain(mLhandID, desiredChainSize, &mIKSkeleton);
	mRhandIKchain = createIKchain(mRhandID, desiredChainSize, &mIKSkeleton);
	mLfootIKchain = createIKchain(mLfootID, desiredChainSize, &mIKSkeleton);
	mRfootIKchain = createIKchain(mRfootID, desiredChainSize, &mIKSkeleton);

	if (mLhandIKchain.getSize() > 1 && mRhandIKchain.getSize() > 1 && mLfootIKchain.getSize() > 1 && mRfootIKchain.getSize() > 1)
	{
		validChains = true;

		// initalize end joint target transforms for Lhand, Rhand, Lfoot and Rfoot based on current position and orientation of joints
		mIKSkeleton.copyTransforms(m_pSkeleton);
		mLhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mLhandID)->getLocal2Global());
		mRhandTarget.setLocal2Global(mIKSkeleton.getJointByID(mRhandID)->getLocal2Global());
		mLfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mLfootID)->getLocal2Global());
		mRfootTarget.setLocal2Global(mIKSkeleton.getJointByID(mRfootID)->getLocal2Global());
	}

	return validChains;
}

int IKController::computePseudoInvIK(ATarget target, AIKchain& IKchain, ASkeleton* pIKSkeleton)
{
	bool result = false;
	// TODO: Implement Pseudo Inverse-based IK  
	// The actual position of the end joint should match the target position after the skeleton is updated with the new joint angles
	//TODO: YOUR PSEUDO INVERSE CODE GOES HERE
	AJoint* endJoint = IKchain.getJoint(0);
	vec3 endPos = endJoint->getGlobalTranslation();
	vec3 targetPos = target.getGlobalTranslation();
	vec3 errorVector = targetPos - endPos;
	//aMatrix Jacobian(3, IKchain.getChain().size())

	return result;
}


bool IKController::IKSolver_PseudoInv(int endJointID, const ATarget& target)
{
	bool validChains = false;
	
		if (!mvalidPseudoInvIKchains)
		{
			mvalidPseudoInvIKchains = createPseudoInvIKchains();
			//assert(mvalidCCDIKchains);
		}

		// copy transforms from base skeleton
		mIKSkeleton.copyTransforms(m_pSkeleton);

		vec3 desiredRootPosition;

		if (endJointID == mLhandID)
		{
			mLhandTarget = target;
			computePseudoInvIK(mLhandTarget, mLhandIKchain, &mIKSkeleton);
		}
		else if (endJointID == mRhandID)
		{
			mRhandTarget = target;
			computePseudoInvIK(mRhandTarget, mRhandIKchain, &mIKSkeleton);
		}
		else if (endJointID == mLfootID)
		{
			mLfootTarget = target;
			computePseudoInvIK(mLfootTarget, mLfootIKchain, &mIKSkeleton);
		}
		else if (endJointID == mRfootID)
		{
			mRfootTarget = target;
			computePseudoInvIK(mRfootTarget, mRfootIKchain, &mIKSkeleton);
		}
		else if (endJointID == mRootID)
		{
			desiredRootPosition = target.getGlobalTranslation();
			mIKSkeleton.getJointByID(mRootID)->setLocalTranslation(desiredRootPosition);
			mIKSkeleton.update();
			computePseudoInvIK(mLhandTarget, mLhandIKchain, &mIKSkeleton);
			computePseudoInvIK(mRhandTarget, mRhandIKchain, &mIKSkeleton);
			computePseudoInvIK(mLfootTarget, mLfootIKchain, &mIKSkeleton);
			computePseudoInvIK(mRfootTarget, mRfootIKchain, &mIKSkeleton);
		}
		else
		{
			mIKchain = createIKchain(endJointID, -1, &mIKSkeleton);
			computePseudoInvIK(target, mIKchain, &mIKSkeleton);
		}

			// update IK Skeleton transforms
		mIKSkeleton.update();

		// copy IK skeleton transforms to main skeleton
		m_pSkeleton->copyTransforms(&mIKSkeleton);
	
	return true;
}


bool IKController::IKSolver_Other(int endJointID, const ATarget& target)
{
	// TODO: Put Optional IK implementation or enhancements here
	 
	return true;
}