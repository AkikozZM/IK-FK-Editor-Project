#include "aActor.h"

#pragma warning(disable : 4018)



/****************************************************************
*
*    	    Actor functions
*
****************************************************************/

AActor::AActor() 
{
	m_pInternalSkeleton = new ASkeleton();
	m_pSkeleton = m_pInternalSkeleton;

	m_BVHController = new BVHController();
	m_BVHController->setActor(this);

	m_IKController = new IKController();
	m_IKController->setActor(this);

	m_BehaviorController = new BehaviorController();
	m_BehaviorController->setActor(this);

	// code to update additional Actor data goes here
	resetGuide();

}

AActor::AActor(const AActor* actor)
{
	*this = *actor;
}

AActor& AActor::operator = (const AActor& actor)
{
	// Performs a deep copy
	if (&actor == this)
	{
		return *this;
	}
	m_pSkeleton = actor.m_pSkeleton;

	// code to update additional Actor data goes here


	return *this;
}

AActor::~AActor()
{
	 delete m_BehaviorController;
	 delete m_IKController;
	 delete m_BVHController;
	 delete m_pInternalSkeleton;

}

void AActor::clear()
{
	// looks like it is clearing more times than the number of actors.  as a result, m_pSkeleton is not defined for last case.
	m_pSkeleton->clear();  

	// code to update additional Actor data goes here
}

void AActor::update()
{
	if (!m_pSkeleton->getRootNode() )
		 return; // Nothing loaded
	else m_pSkeleton->update();

	// code to update additional Actor data goes here

}

ASkeleton* AActor::getSkeleton()
{
	return m_pSkeleton;
}

void AActor::setSkeleton(ASkeleton* pExternalSkeleton)
{
	m_pSkeleton = pExternalSkeleton;
}

void AActor::resetSkeleton()
{
	m_pSkeleton = m_pInternalSkeleton;
}

BVHController* AActor::getBVHController()
{
	return m_BVHController;
}

IKController* AActor::getIKController()
{
	return m_IKController;
}

BehaviorController* AActor::getBehaviorController()
{
	return m_BehaviorController;
}

quat FromTwoVectors(const vec3& v1, const vec3& v2) {
	vec3 cross = v1.Cross(v2);
	double dot = Dot(v1, v2);
	double s = sqrt((1 + dot) * 2);
	double invs = 1 / s;

	quat q(s * 0.5, cross[0] * invs, cross[1] * invs, cross[2] * invs);
	return q;
}

void AActor::updateGuideJoint(vec3 guideTargetPos)
{
	AJoint* rootJoint = m_pSkeleton->getRootNode();
	if (!rootJoint) { 
		return; 
	}
	// TODO: 
	// 1.	Set the global position of the guide joint to the global position of the root joint
	vec3 rootGlobalPos = m_Guide.getLocal2Global() * rootJoint->getGlobalTranslation();
	m_Guide.setGlobalTranslation(rootGlobalPos);
	// 2.	Set the y component of the guide position to 0
	vec3 pos = m_Guide.getGlobalTranslation();
	pos[1] = 0;
	m_Guide.setGlobalTranslation(pos);
	// 3.	Set the global rotation of the guide joint towards the guideTarget
	// Hint: Return of getGlobalRotation() here is in target space but not in world/global space
	vec3 guidePos = m_Guide.getGlobalTranslation();
	vec3 target = guideTargetPos - guidePos;
	target[1] = 0;
	target.Normalize();

	vec3 forward = vec3(0,0,1);
	vec3 cross = forward.Cross(target);
	double dot = Dot(forward, target);
	double s = sqrt((1 + dot) * 2);
	double invs = 1 / s;
	quat q(s * 0.5, cross[0] * invs, cross[1] * invs, cross[2] * invs);

	mat3 targetRotation = q.ToRotation();

	// convert to global
	m_Guide.setGlobalRotation(targetRotation);
}

void AActor::solveFootIK(float leftHeight, float rightHeight, bool rotateLeft, bool rotateRight, vec3 leftNormal, vec3 rightNormal)
{
	if (!m_pSkeleton->getRootNode()) { 
		return; 
	}
	AJoint* rootJoint = m_pSkeleton->getRootNode();
	AJoint* leftFoot = m_pSkeleton->getJointByID(m_IKController->mLfootID);
	AJoint* rightFoot = m_pSkeleton->getJointByID(m_IKController->mRfootID);

	// TODO: 
	// The normal and the height given are in the world space
	// Hint: Return of getGlobal***() here is in target space but not in world/global space

	// 1.	Update the local translation of the root based on the left height and the right height
	vec3 rootPos = rootJoint->getLocalTranslation();
	float avgHeight = (leftHeight + rightHeight) / 2.0f;
	rootPos[1] += avgHeight;
	rootJoint->setLocalTranslation(rootPos);

	m_pSkeleton->update();

	// 2.	Update the character with Limb-based IK 
	vec3 leftFootTargetPos = leftFoot->getGlobalTranslation();
	leftFootTargetPos[1] += leftHeight;
	ATarget leftFootTarget;
	leftFootTarget.setGlobalTranslation(leftFootTargetPos);
	m_IKController->IKSolver_Limb(m_IKController->mLfootID, leftFootTarget);

	vec3 rightFootTargetPos = rightFoot->getGlobalTranslation();
	rightFootTargetPos[1] += rightHeight;
	ATarget rightFootTarget;
	rightFootTarget.setGlobalTranslation(rightFootTargetPos);
	m_IKController->IKSolver_Limb(m_IKController->mRfootID, rightFootTarget);
	
	m_pSkeleton->update();



	// Rotate Foot
	if (rotateLeft)
	{
		// Update the local orientation of the left foot based on the left normal
		ATransform& Linverse = leftFoot->getLocal2Global().Inverse();
		leftNormal = Linverse.Rotate(leftNormal);
		double angle = acos(Dot(vec3(0, 1, 0), leftNormal.Normalize()));
		vec3 axis = vec3(0, 1, 0).Cross(leftNormal).Normalize();
		mat3 rot;
		rot.FromAxisAngle(axis, angle);
		leftFoot->setLocalRotation(leftFoot->getLocalRotation() * rot);
	}
	if (rotateRight)
	{
		// Update the local orientation of the right foot based on the right normal
		ATransform& Rinverse = rightFoot->getLocal2Global().Inverse();
		rightNormal = Rinverse.Rotate(rightNormal);
		double angle = acos(Dot(vec3(0, 1, 0), rightNormal.Normalize()));
		vec3 axis = vec3(0, 1, 0).Cross(rightNormal).Normalize();
		mat3 rot;
		rot.FromAxisAngle(axis, angle);
		rightFoot->setLocalRotation(rightFoot->getLocalRotation() * rot);
	}
	m_pSkeleton->update();
}
