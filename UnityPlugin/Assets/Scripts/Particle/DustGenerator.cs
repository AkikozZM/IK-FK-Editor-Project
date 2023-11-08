using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class DustGenerator : MonoBehaviour
{
    public bool m_enableDust = false;
    public FKIKCharacterController m_jointController;
    public GameObject m_dustParticleSystem;

    // Copy the states from jointController
    [SerializeField] private bool m_leftFootGroundContact = true;
    [SerializeField] private bool m_rightFootGroundContact = true;

    [SerializeField] private List<AParticleSystem> m_particleSystems = new List<AParticleSystem>();
    // Start is called before the first frame update
    void Start()
    {
        
    }

    // Update is called once per frame
    void Update()
    {
        // Clear particle systems that are not alive
        foreach (AParticleSystem p in m_particleSystems)
        {
            if (!p.IsAlive())
            {
                p.Destroy();
            }
        }
        m_particleSystems.RemoveAll(item => item == null);

        if (!m_enableDust) { return; }
        if (m_leftFootGroundContact != m_jointController.m_leftFootGroundContact)
        {
            m_leftFootGroundContact = m_jointController.m_leftFootGroundContact;
            // Left foot starts to contact the ground
            if (m_leftFootGroundContact)
            {
                GameObject dustObject = Instantiate(m_dustParticleSystem);
                AParticleSystem dust = dustObject.GetComponent<AParticleSystem>();
                dust.m_position = m_jointController.m_leftFoot.transform.position;
                m_particleSystems.Add(dust);
            }
        }
        if (m_rightFootGroundContact != m_jointController.m_rightFootGroundContact)
        {
            m_rightFootGroundContact = m_jointController.m_rightFootGroundContact;
            // Right foot starts to contact the ground
            if (m_rightFootGroundContact)
            {
                GameObject dustObject = Instantiate(m_dustParticleSystem);
                AParticleSystem dust = dustObject.GetComponent<AParticleSystem>();
                dust.m_position = m_jointController.m_rightFoot.transform.position;
                m_particleSystems.Add(dust);
            }
        }
        
    }
}
