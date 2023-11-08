using System.Collections;
using System.Collections.Generic;
using System.IO;
using UnityEditor;
using UnityEngine;

public class FKIKCharacterController : MonoBehaviour
{
    public enum ControlMode { KEYBOARD, TARGET, BEHAVIOR };
    public ControlMode m_controlMode = ControlMode.KEYBOARD;

    [Header("Target Joints")]
    public GameObject m_root;   // Hip joint
    public GameObject m_leftHand;
    public GameObject m_rightHand;
    public GameObject m_leftFoot;
    public GameObject m_rightFoot;

    // Guide is the gameObject.transform
    public GameObject m_guideTarget;    // The target of the guide 

    [Header("FK")]
    [Range(0.0f, 2.0f)] public float m_speedScale = 1.0f;
    [SerializeField] private float m_t = 0;
    [SerializeField] private float m_maxt = 0;
    [SerializeField] private int m_keyID = 0;
    [SerializeField] private int m_keySize = 0;
    public Object m_defaultBVH;

    public bool m_fixedRoot = false;    // If true, the root local position x and z will be zero (used in keyboard mode)

    // IK constraint
    [Header("IK")]
    public bool m_enableFootIK = false;
    [Range(0, 1)] public float m_rootSpeedY = 0.25f;
    [Range(0, 1)] public float m_footSpeedY = 0.5f;
    private float m_lastRootOffsetY = 0, m_lastLeftFootY = 0, m_lastRightFootY = 0;
    // Record the foot-ground contact state
    public bool m_leftFootGroundContact = false;
    public bool m_rightFootGroundContact = false;
    public float m_footGroundOffsetMin = 13;
    public float m_footGroundOffsetMax = 15;

    // Raycast parameters
    public float m_footAboveOffset = 1000;
    public LayerMask m_environmentLayer;

    // Containers
    private List<GameObject> m_joints = new List<GameObject>();
    private FKIKPlugin.JointData[] m_jointDataArray;
    private SortedDictionary<int, GameObject> m_idJointMap = new SortedDictionary<int, GameObject>();   // <Joint ID, Joint object>

    // States
    private bool m_canMove = false;
    private bool m_bvhLoaded = false;
    private bool m_bvhNotMatch = true;
    private int m_id;

    // Control the virtual joints
    private JointPainter m_jointPainter;

    // Help functions
    // Cast the y component of a vector3 to 0
    Vector3 Vector3ToGround(Vector3 pos)
    {
        return new Vector3(pos.x, 0, pos.z);
    }

    // Start is called before the first frame update
    void Start()
    {
        m_id = FKIKPlugin.CreateActor();
        m_jointPainter = this.GetComponent<JointPainter>();

        // Initialize the plugin and virtual joints        
        TraverseJoints(m_root.transform);   // Add joints to the m_joints and m_idJointMap
        LinkJoints();   // Set the parent child relationship
        m_jointPainter.InitVirtualJoints(this.name + "joints", m_joints);       
        ConstructJointDataArray();  // Initialize the joint data array  
        ConstructIKChains();    // Construct IK chains for hands and feet
        GameObject[] endJoints = new GameObject[] { m_root, m_leftHand, m_rightHand, m_leftFoot, m_rightFoot };
        m_jointPainter.InitializeTargetJoints(endJoints);

        if (m_defaultBVH)
        {
            string projectPath = Application.dataPath;
            string filePath = AssetDatabase.GetAssetPath(m_defaultBVH);
            string filename = Path.Combine(projectPath, "../", filePath);
            Debug.Log(filename);
            LoadBVHFile(filename);
            StartMove();         
        }
        
    }

    private void Update()
    {
        InputProcess();
        m_jointPainter.DrawJointsAndSkeletons();
    }

    private void FixedUpdate()
    {
        if (!m_bvhNotMatch) { return; }
          
        if (m_canMove)
        {
            m_t += Time.fixedDeltaTime * m_speedScale;
            if (m_t >= m_maxt)
            {
                m_t -= m_maxt;
                switch(m_controlMode)
                {
                    case ControlMode.TARGET:
                        UpdateGuideByTarget(true);  // Update Guide to the new root position
                        break;
                    case ControlMode.KEYBOARD:
                        UpdateGuideToRoot();
                        break;
                    case ControlMode.BEHAVIOR:
                        break;
                }
            }
            
        }
        UpdateJointTransform(m_t);   
    }

    private void InputProcess()
    {
        if (m_controlMode != ControlMode.KEYBOARD) { return; }
        if (Input.GetKey(KeyCode.W))
        {
            m_speedScale += 0.1f;
            m_speedScale = Mathf.Min(m_speedScale, 2);
        }
        if (Input.GetKey(KeyCode.S))
        {
            m_speedScale -= 0.1f;
            m_speedScale = Mathf.Max(m_speedScale, 0);
        }
        if (!m_bvhLoaded) { return; }
        if (Input.GetKey(KeyCode.A))
        {
            Vector3 guidePos = m_root.transform.InverseTransformPoint(gameObject.transform.position);
            guidePos = Quaternion.AngleAxis(-5, Vector3.up) * guidePos;
            guidePos = m_root.transform.TransformPoint(guidePos);
            gameObject.transform.position = guidePos;
            gameObject.transform.rotation = gameObject.transform.rotation * Quaternion.AngleAxis(-5, Vector3.up);
        }
        if (Input.GetKey(KeyCode.D))
        {
            Vector3 guidePos = m_root.transform.InverseTransformPoint(gameObject.transform.position);
            guidePos = Quaternion.AngleAxis(5, Vector3.up) * guidePos;
            guidePos = m_root.transform.TransformPoint(guidePos);
            gameObject.transform.position = guidePos;
            gameObject.transform.rotation = gameObject.transform.rotation * Quaternion.AngleAxis(5, Vector3.up);
        }
    }

    public void SwitchControlMode(ControlMode newMode)
    {
        m_controlMode = newMode;
        switch (newMode)
        {
            case ControlMode.TARGET:
            case ControlMode.KEYBOARD:
                m_fixedRoot = false;
                break;
            case ControlMode.BEHAVIOR:
                m_controlMode = newMode;
                m_fixedRoot = true;
                UpdateGuideToRoot();
                break;
        }
    }

    // Set joints' localRotation and localPosition with JointData
    private void ApplyJointDataArray()
    {
        for (int i = 0; i < m_joints.Count; ++i)
        {
            FKIKPlugin.JointData data = m_jointDataArray[i];
            GameObject joint = m_idJointMap[data.id];
            joint.transform.localPosition = FKIKPlugin.BVHToUnityTranslation(data.localTranslation);
            joint.transform.localRotation = FKIKPlugin.BVHToUnityQuaternion(data.localRotation);
        }

        if (m_fixedRoot)
        {
            // Set root local position x and z to zero
            m_root.transform.localPosition = new Vector3(0, m_root.transform.localPosition.y, 0);
        }
    }

    private void FootIK()
    {
        // Get foot position in world space
        Vector3 leftFootPos = m_leftFoot.transform.position;
        Vector3 rightFootPos = m_rightFoot.transform.position;
        /// C++ Implementation
        RaycastFootHeight(leftFootPos, out float leftHeight, out Vector3 leftNormal);
        RaycastFootHeight(rightFootPos, out float rightHeight, out Vector3 rightNormal);
        Debug.Log("L: " + leftHeight);
        Debug.Log("R: " + rightHeight);
        FKIKPlugin.SolveFootIK(m_id, leftHeight, rightHeight, m_leftFootGroundContact, m_rightFootGroundContact,
            FKIKPlugin.UnityToBVHTranslation(leftNormal), FKIKPlugin.UnityToBVHTranslation(rightNormal));

        FKIKPlugin.GetJointData(m_id, m_jointDataArray, m_joints.Count);
        ApplyJointDataArray();

        /// Unity Implementation

        //if (RaycastFootHeight(leftFootPos, out float leftHeight, out Vector3 leftNormal))
        //{
        //    leftFootPos.y = Mathf.Lerp(leftHeight + leftFootPos.y, m_lastLeftFootY, m_footSpeedY);
        //}
        //if (RaycastFootHeight(rightFootPos, out float rightHeight, out Vector3 rightNormal))
        //{
        //    rightFootPos.y = Mathf.Lerp(rightHeight + rightFootPos.y, m_lastRightFootY, m_footSpeedY);
        //}

        //m_lastLeftFootY = leftFootPos.y;
        //m_lastRightFootY = rightFootPos.y;
        //// Transform to guide space
        //leftFootPos = gameObject.transform.InverseTransformPoint(leftFootPos); // World space to guide space
        //rightFootPos = gameObject.transform.InverseTransformPoint(rightFootPos); // World space to guide space

        //// Move root height
        //float rootOffsetY = leftFootPos.y < rightFootPos.y ? leftHeight : rightHeight;
        //rootOffsetY = Mathf.Lerp(m_lastRootOffsetY, rootOffsetY, m_rootSpeedY);
        //m_lastRootOffsetY = rootOffsetY;
        //Vector3 rootPos = m_root.transform.position + new Vector3(0, rootOffsetY, 0);
        //rootPos = gameObject.transform.InverseTransformPoint(rootPos);

        //FKIKPlugin.SetRootJointTranslation(m_id, FKIKPlugin.UnityToBVHTranslation(rootPos));
        //FKIKPlugin.UpdateSkeleton(m_id);
        //FKIKPlugin.SolveLimbIK(m_id, FindIdByJoint(m_leftFoot), FKIKPlugin.UnityToBVHTranslation(leftFootPos));
        //FKIKPlugin.SolveLimbIK(m_id, FindIdByJoint(m_rightFoot), FKIKPlugin.UnityToBVHTranslation(rightFootPos));
        //FKIKPlugin.GetJointData(m_id, m_jointDataArray, m_joints.Count);

        //ApplyJointDataArray();

        //// Foot rotation 
        //if (m_leftFootGroundContact)
        //{
        //    Vector3 leftForward = Vector3.Cross(m_leftFoot.transform.right, leftNormal); // Calculate new forward vector
        //    m_leftFoot.transform.up = leftNormal;
        //    m_leftFoot.transform.forward = leftForward;
        //}
        //if (m_rightFootGroundContact)
        //{
        //    Vector3 rightForward = Vector3.Cross(m_rightFoot.transform.right, rightNormal);
        //    m_rightFoot.transform.up = rightNormal;
        //    m_rightFoot.transform.forward = rightForward;
        //}
    }

    public void CheckFootGroundContact()
    {
        if (!m_leftFootGroundContact && m_leftFoot.transform.position.y < m_footGroundOffsetMin)
        {
            m_leftFootGroundContact = true;
        }
        if (m_leftFootGroundContact && m_leftFoot.transform.position.y > m_footGroundOffsetMax)
        {
            m_leftFootGroundContact = false;
        }
        if (!m_rightFootGroundContact && m_rightFoot.transform.position.y < m_footGroundOffsetMin)
        {
            m_rightFootGroundContact = true;
        }
        if (m_rightFootGroundContact && m_rightFoot.transform.position.y > m_footGroundOffsetMax)
        {
            m_rightFootGroundContact = false;
        }
    }

   private void UpdateJointTransform(float t)
    {
        // FK
        if (m_bvhLoaded)
        {
            FKIKPlugin.UpdateBVHSkeleton(m_id, t);
            FKIKPlugin.GetJointData(m_id, m_jointDataArray, m_joints.Count);
        }
        else
        {
            FKIKPlugin.UpdateSkeleton(m_id);
            FKIKPlugin.GetJointData(m_id, m_jointDataArray, m_joints.Count);
        }
        ApplyJointDataArray();

        // Check foot-ground contact
        CheckFootGroundContact();

        // IK
        if (m_enableFootIK)
        {
            FootIK();
        }
    }

    #region Initialization
    // Construct the jointData array to transfer data between C# and C++
    void ConstructJointDataArray()
    {
        m_jointDataArray = new FKIKPlugin.JointData[m_joints.Count];
        for (int i = 0; i < m_joints.Count; ++i)
        {
            m_jointDataArray[i] = new FKIKPlugin.JointData(i);
        }
    }

    // Construct a map between joints and indices
    bool ConstructIdJointMap()
    {
        m_idJointMap.Clear();
        foreach (GameObject joint in m_joints)
        {
            int id = -1;
            // End Joint: can only be found by its parent's name
            if (joint.transform.childCount == 0)
            {
                id = FKIKPlugin.GetJointIdByParentName(m_id, joint.transform.parent.name);
            }
            else
            {
                id = FKIKPlugin.GetJointIdByName(m_id, joint.transform.name);
            }

            // Not found
            if (id == -1)
            {
                m_idJointMap.Clear();
                Debug.Log("Joint Not Found:" + joint.name);
                return false;
            }
            m_idJointMap.Add(id, joint);
        }
        return true;
    }

    // Set id of end joints and construct the IK chains
    void ConstructIKChains()
    {
        foreach(KeyValuePair<int, GameObject> pair in m_idJointMap)
        {
            if (pair.Value == m_root)
            {
                //Debug.Log("Root:" + pair.Key.ToString());
                FKIKPlugin.SetRootID(m_id, pair.Key);
            }
            else if (pair.Value == m_leftFoot)
            {
                //Debug.Log("Left Foot:" + pair.Key.ToString());
                FKIKPlugin.SetLeftFootID(m_id, pair.Key);
            }
            else if (pair.Value == m_rightFoot)
            {
                //Debug.Log("Right Foot:" + pair.Key.ToString());
                FKIKPlugin.SetRightFootID(m_id, pair.Key);
            }
            else if (pair.Value == m_leftHand)
            {
                //Debug.Log("Left Hand:" + pair.Key.ToString());
                FKIKPlugin.SetLeftHandID(m_id, pair.Key);
            }
            else if (pair.Value == m_rightHand)
            {
                //Debug.Log("Right Hand:" + pair.Key.ToString());
                FKIKPlugin.SetRightHandID(m_id, pair.Key);
            }
        }
        FKIKPlugin.CreateLimbIKChains(m_id);
    }

    void TraverseJoints(Transform transform)
    {
        m_joints.Add(transform.gameObject);

        // Initialize joint in plugin
        int id = FKIKPlugin.CreateJoint(m_id, transform.gameObject.name, transform.gameObject == m_root);
        m_idJointMap.Add(id, transform.gameObject);

        // Travser childeren joints
        foreach (Transform child in transform)
        {
            TraverseJoints(child);
        }
    }

    void LinkJoints()
    {
        foreach (GameObject joint in m_joints)
        {
            FKIKPlugin.JointData data;
            data.id = FindIdByJoint(joint);
            data.localTranslation = FKIKPlugin.UnityToBVHTranslation(joint.transform.localPosition);
            data.localRotation = FKIKPlugin.UnityToBVHQuaternion(joint.transform.localRotation);
            int parentID = FindIdByJoint(joint.transform.parent.gameObject);
            FKIKPlugin.SetJointData(m_id, data, parentID);
        }
        FKIKPlugin.UpdateSkeleton(m_id);
    }


    // Return false if the BVH file and the FBX model do not match
    public bool LoadBVHFile(string filepath)
    {
        m_idJointMap.Clear();
        m_bvhLoaded = false;
        // m_canMove = false;
        m_bvhNotMatch = true;
        if (FKIKPlugin.LoadBVH(m_id, filepath))
        {
            // Reconstrcut idJoint Map. If files do not match, return false
            if (!ConstructIdJointMap())
            {
                m_canMove = false;
                m_bvhNotMatch = false;
                return false;
            }

            m_bvhLoaded = true;
            // m_canMove = true;
            m_t = 0;
            m_maxt = FKIKPlugin.GetDuration(m_id);
            m_keySize = FKIKPlugin.GetKeySize(m_id);
            m_lastRootOffsetY = m_root.transform.position.y;
            m_lastLeftFootY = m_leftFoot.transform.position.y;
            m_lastRightFootY = m_rightFoot.transform.position.y;

            // Reset end joint ids
            ConstructIKChains();
            return true;
        }
        return false;
    }

    #endregion

    #region Getter And Setter
    public int GetKeySize()
    {
        return m_keySize;
    }

    public void SetTimeByProgressBar(float id)
    {
        if (id >= m_keySize)
        {
            Debug.LogError("Key id out of range");
            return;
        }
        m_t = FKIKPlugin.GetKeyTime(m_id, (int)id);
        m_keyID = (int)id;
    }

    public bool StartMove()
    {
        if (m_bvhLoaded)
        {
            m_canMove = !m_canMove;
            return true;
        }
        return false;
    }

    public void EditKeyRotation(GameObject virtualJoint)
    {
        int jointID = FindIdByJoint(virtualJoint.GetComponent<VirtualJointController>().GetBindedJoint());
        FKIKPlugin.SetJointRotation(m_id, jointID, FKIKPlugin.UnityToBVHQuaternion(virtualJoint.transform.localRotation));
    }

    public void SetSpeedScale(float scale)
    {
        m_speedScale = scale;
    }

    // [TODO]: Cache the id
    public int FindIdByJoint(GameObject joint)
    {
        foreach (KeyValuePair<int, GameObject> pair in m_idJointMap)
        {
            if (pair.Value == joint)
            {
                return pair.Key;
            }
        }
        return -1;
    }

    public void SetEnableFootIK(bool enable)
    {
        m_enableFootIK = enable;
    }

    public bool BVHLoaded()
    {
        return m_bvhLoaded;
    }
    #endregion

    #region IKSolver
    public void SolveIK(GameObject targetJoint, Vector3 pos)
    {
        int jointID = FindIdByJoint(targetJoint.GetComponent<VirtualJointController>().GetBindedJoint());
        FKIKPlugin.SolveLimbIK(m_id, jointID, FKIKPlugin.UnityToBVHTranslation(pos));
    }

    // Cast a ray down to get the height of the environment
    private bool RaycastFootHeight(Vector3 footPos, out float height, out Vector3 normal)
    {
        Vector3 origin = footPos + Vector3.up * m_footAboveOffset;
        RaycastHit hit = new RaycastHit();
        if (Physics.Raycast(origin, Vector3.down, out hit, Mathf.Infinity, m_environmentLayer))
        {
            height = hit.point.y;
            normal = hit.normal;
            return true;
        }
        height = 0;
        normal = new Vector3(0, 1, 0);
        return false;
    }
    #endregion

    #region Target
    void UpdateGuideByTarget(bool rotateToTarget)
    {
        float[] newPos = new float[3];
        float[] newQuat = new float[4];
        FKIKPlugin.UpdateGuideJointByTarget(m_id, FKIKPlugin.UnityToBVHTranslation(m_guideTarget.transform.position),
            newPos, newQuat);
        // Update guide
        gameObject.transform.position = FKIKPlugin.BVHToUnityTranslation(newPos);
        gameObject.transform.rotation = FKIKPlugin.BVHToUnityQuaternion(newQuat);
    }

    void UpdateGuideToRoot()
    {
        gameObject.transform.localRotation = Quaternion.LookRotation(
            Vector3ToGround(m_root.transform.position) - Vector3ToGround(gameObject.transform.position), Vector3.up);
        gameObject.transform.position = Vector3ToGround(m_root.transform.position);
    }

    #endregion


}
