using RTG;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class JointPainter : MonoBehaviour
{

    public GameObject m_jointPrefab;
    public GameObject m_skeletonPrefab;
    public GameObject m_targetPrefab;
    public GameObject m_root;

    public GameObject m_jointSkeletonParent = null;
    public GameObject m_targetJointParent = null;
    // Map between virtual joint and the joint game object
    public Dictionary<GameObject, GameObject> m_jointMap = new Dictionary<GameObject, GameObject>();
    public List<GameObject> m_targetJoints = new List<GameObject>();

    public bool m_drawSkeleton = false;

    // Start is called before the first frame update
    void Start()
    {
        
    }

    // Update is called once per frame
    void Update()
    {
    }

    public void InitVirtualJoints(string name, List<GameObject> joints)
    {
        if (m_jointSkeletonParent == null)
        {
            m_jointSkeletonParent = new GameObject(name);   // An empty game object
        }
        foreach (GameObject joint in joints)
        { 
            GameObject virtualJoint = Instantiate(m_jointPrefab, joint.transform.position, joint.transform.rotation, 
                m_jointSkeletonParent.transform);
            VirtualJointController controller = virtualJoint.AddComponent<VirtualJointController>();
            
            m_jointMap.Add(virtualJoint, joint);

            GameObject virtualSkeleton = Instantiate(m_skeletonPrefab, m_jointSkeletonParent.transform);

            controller.Initialize(joint, joint == m_root ? null : joint.transform.parent.gameObject, virtualSkeleton);
        }

        // Draw skeletons or not
        if (m_drawSkeleton)
        {
            DrawJointsAndSkeletons();
        }
        else
        {
            SetDrawSkeletons(false);
        }
    }

    public void DrawJointsAndSkeletons()
    {
        if (!m_jointSkeletonParent.activeSelf)
        {
            return;
        }
        // Update the transform of the virtual joints
        foreach (KeyValuePair<GameObject, GameObject> pair in m_jointMap)
        {
            pair.Key.transform.position = pair.Value.transform.position;
            pair.Key.transform.localRotation = pair.Value.transform.localRotation;
            pair.Key.GetComponent<VirtualJointController>().DrawSkeleton();
        }

    }

    public void SetDrawSkeletons(bool draw)
    {
        m_jointSkeletonParent.SetActive(draw);
        DrawJointsAndSkeletons();
    }

    public void InitializeTargetJoints(GameObject[] endJoints)
    {
        if (m_targetJointParent == null)
        {
            m_targetJointParent = new GameObject("Target");   // An empty game object
        }
        foreach (GameObject endJoint in endJoints)
        {
            GameObject targetJoint = Instantiate(m_targetPrefab, endJoint.transform.position, endJoint.transform.rotation,
                m_targetJointParent.transform);
            VirtualJointController controller = targetJoint.AddComponent<VirtualJointController>();
            controller.Initialize(endJoint, null, null);
            m_targetJoints.Add(targetJoint);
        }
        // Do not draw targets at first
        SetDrawTargets(false);
    }

    public void SetDrawTargets(bool draw)
    {
        m_targetJointParent.SetActive(draw);
    }

    // Reset the position of the target joints
    public void ResetTargetPosition()
    {
        foreach (GameObject targetJoint in m_targetJoints)
        {
            GameObject bindedJoint = targetJoint.GetComponent<VirtualJointController>().GetBindedJoint();
            targetJoint.transform.position = bindedJoint.transform.position;
            targetJoint.transform.localRotation = bindedJoint.transform.localRotation;
        }
    }

}
