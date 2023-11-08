using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class FKIKSpeedScaleController : MonoBehaviour
{
    public FKIKCharacterController m_characterController;
    private Vector3 m_lastPosition;

    // Start is called before the first frame update
    void Start()
    {
        m_lastPosition = this.transform.position;
    }

    // Update is called once per frame
    void FixedUpdate()
    {
        Vector3 vel = (m_lastPosition - this.transform.position) / Time.fixedDeltaTime;
        m_lastPosition = this.transform.position;
        m_characterController.SetSpeedScale(vel.magnitude / 164.0f);
    }

}
