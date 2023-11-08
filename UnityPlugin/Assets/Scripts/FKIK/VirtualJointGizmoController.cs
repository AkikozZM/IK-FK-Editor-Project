using RTG;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class VirtualJointGizmoController : MonoBehaviour
{
    public LayerMask m_layer;
    public FKIKCharacterController m_jointController;
    private GameObject m_selectedJoint;
    private GameObject m_hoveredJoint;
    private ObjectTransformGizmo m_objectRotationGizmo;
    private ObjectTransformGizmo m_objectTranslationGizmo;

    public enum Mode { FK, IK };
    public Mode m_mode = Mode.FK;

    // Start is called before the first frame update
    void Start()
    {
        m_objectRotationGizmo = RTGizmosEngine.Get.CreateObjectRotationGizmo();
        m_objectRotationGizmo.SetTransformSpace(GizmoSpace.Local);
        m_objectRotationGizmo.Gizmo.SetEnabled(false);
        m_objectRotationGizmo.Gizmo.PostUpdateEnd += OnRotationGizmoPostUpdateEnd;

        m_objectTranslationGizmo = RTGizmosEngine.Get.CreateObjectMoveGizmo();
        m_objectTranslationGizmo.SetTransformSpace(GizmoSpace.Local);
        m_objectTranslationGizmo.Gizmo.SetEnabled(false);
        m_objectTranslationGizmo.Gizmo.PostUpdateEnd += OnTranslationGizmoPostUpdateEnd;
    }

    // Update is called once per frame
    void Update()
    {
        // Check selected joint
        if (Input.GetMouseButtonDown(0) && RTGizmosEngine.Get.HoveredGizmo == null)
        {
            GameObject pickedObject = PickVirtualJoint();
            if (pickedObject != m_selectedJoint)
            {
                ChangeSelectedJoint(pickedObject);
            }
        }
        else if (RTGizmosEngine.Get.HoveredGizmo == null)
        {
            GameObject hoveredObject = PickVirtualJoint();
            if (hoveredObject != m_hoveredJoint)
            {
                ChangeHoveredJoint(hoveredObject);
            }
        }
    }

    private GameObject PickVirtualJoint()
    {
        // Build a ray using the current mouse cursor position
        Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);

        // Check if the ray intersects a game object. If it does, return it
        RaycastHit hit;

        // First check if the ray hit the virtual joint
        // Bit shift the index of the layer (9) to get a bit mask
        if (Physics.Raycast(ray, out hit, Mathf.Infinity, m_layer))
        {
            return hit.collider.gameObject;
        }

        // No object is intersected by the ray. Return null.
        return null;
    }

    public void SetFKIKMode(Mode m)
    {
        m_mode = m;
    }

    public void ChangeSelectedJoint(GameObject newSeletcedObject)
    {
        // Recover the color of the previous seleced joint
        if (m_selectedJoint != null)
        {
            m_selectedJoint.GetComponent<VirtualJointController>().ExitJoint();
        }

        // Store the new target object
        m_selectedJoint = newSeletcedObject;

        m_objectRotationGizmo.Gizmo.SetEnabled(false);
        m_objectTranslationGizmo.Gizmo.SetEnabled(false);

        // Check if the new target object is valid
        if (m_selectedJoint != null)
        {
            if (m_mode == Mode.FK)
            {
                m_objectRotationGizmo.SetTargetObject(m_selectedJoint);
                m_objectRotationGizmo.Gizmo.SetEnabled(true);
                m_objectTranslationGizmo.SetTargetObject(m_selectedJoint);
                m_objectTranslationGizmo.Gizmo.SetEnabled(true);
            }
            else if (m_mode == Mode.IK && m_selectedJoint.CompareTag("TargetJoint"))
            {
                m_objectTranslationGizmo.SetTargetObject(m_selectedJoint);
                m_objectTranslationGizmo.Gizmo.SetEnabled(true);
            }

            // Change the color of the new selected joint
            m_selectedJoint.GetComponent<VirtualJointController>().SelectJoint();
        }

    }

    private void ChangeHoveredJoint(GameObject newHoveredObject)
    {
        if (m_hoveredJoint != null && m_hoveredJoint != m_selectedJoint)
        {
            m_hoveredJoint.GetComponent<VirtualJointController>().ExitJoint();
        }
        m_hoveredJoint = newHoveredObject;
        
        if (m_hoveredJoint != null)
        {
            m_hoveredJoint.GetComponent<VirtualJointController>().HoverJoint();
        }
    }

    private void OnRotationGizmoPostUpdateEnd(RTG.Gizmo gizmo)
    {
        if (m_selectedJoint == null)
        {
            return;
        }
        
        m_jointController.EditKeyRotation(m_selectedJoint);
        m_objectTranslationGizmo.SetTargetObject(m_selectedJoint);

    }

    private void OnTranslationGizmoPostUpdateEnd(RTG.Gizmo gizmo)
    {
        if (m_selectedJoint == null || !m_selectedJoint.CompareTag("TargetJoint"))
        {
            return;
        }
        m_jointController.SolveIK(m_selectedJoint, m_selectedJoint.transform.position);
    }
}
