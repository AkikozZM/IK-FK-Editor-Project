using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class TargetController : MonoBehaviour
{
    public LayerMask m_environmentLayer;
    public GameObject m_target;
    // Start is called before the first frame update
    void Start()
    {
        
    }

    // Update is called once per frame
    void Update()
    {
        // Raycast to set new target position
        if (Input.GetKey(KeyCode.LeftControl) && Input.GetMouseButtonDown(0))
        {
            Ray ray = Camera.main.ScreenPointToRay(Input.mousePosition);
            if (Physics.Raycast(ray, out RaycastHit hit, Mathf.Infinity, m_environmentLayer))
            {
                m_target.transform.position = hit.point;
                //m_jointController.SetGuideTarget(hit.point);
            }
        }
    }

}
