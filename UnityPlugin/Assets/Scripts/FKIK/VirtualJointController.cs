using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class VirtualJointController : MonoBehaviour
{
    private GameObject m_joint;
    private GameObject m_parentJoint;

    private LineRenderer m_skeletonRenderer;

    private Color m_selectedColor = Color.gray;
    private Color m_initialColor;
    private Color m_hoveredColor = Color.cyan;
    private Renderer m_renderer;

    // Start is called before the first frame update
    void Start()
    {
        m_renderer = this.GetComponent<Renderer>();
        m_initialColor = m_renderer.material.GetColor("_Color");
    }

    // Update is called once per frame
    void Update()
    {

    }

    public void Initialize(GameObject joint, GameObject parentJoint, GameObject skeleton)
    {
        m_joint = joint;
        if (parentJoint == null)
        {
            return;
        }
        m_parentJoint = parentJoint;
        m_skeletonRenderer = skeleton.GetComponent<LineRenderer>();
    }

    public void DrawSkeleton()
    {
        if (m_parentJoint == null)
        {
            return;
        }
        m_skeletonRenderer.positionCount = 2;
        m_skeletonRenderer.SetPosition(0, m_parentJoint.transform.position);
        m_skeletonRenderer.SetPosition(1, m_joint.transform.position);
    }

    public void SelectJoint()
    {
        m_renderer.material.SetColor("_Color", m_selectedColor);
    }

    public void HoverJoint()
    {
        m_renderer.material.SetColor("_Color", m_hoveredColor);
    }

    public void ExitJoint()
    {
        m_renderer.material.SetColor("_Color", m_initialColor);
    }

    public void UpdateBindedJointRotation()
    {
        m_joint.transform.localRotation = this.transform.localRotation;
    }

    public GameObject GetBindedJoint()
    {
        return m_joint;
    }

}
