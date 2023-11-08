using System;
using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class AParticleSystem : MonoBehaviour
{
    private int m_id;     // The id of the particle system
    private ParticlePlugin.ParticleData[] m_particleDataArray;  // The array of the particle data
    public GameObject m_particlePrefab;
    private List<GameObject> m_particles = new List<GameObject>();

    public bool m_infinite = true;
    private int m_maxParticleObjects = 50;  // The max number of maxParticleNum
    [Range (1, 50)] public int m_maxParticles = 20;
    [Range (0, 2.0f)] public float m_lifeTime = 0.4f;
    public Vector3 m_gravity = new Vector3(0, -20.0f, 0);
    public Vector3 m_velocity = new Vector3(0, 40, 0);
    public Vector3 m_position = new Vector3(0, 0, 0);
    public Color m_startColor = new Color(1, 0.8f, 0);
    public Color m_endColor = new Color(0, 0, 0);
    [Range(0, 1)] public float m_startAlpha = 0.4f;
    [Range(0, 10)]public float m_startSize = 4.0f;
    [Range(0, 10)]public float m_endSize = 2.0f;
    public Vector2 m_positionJitter = new Vector2(-10.0f, 10.0f); // min and max
    public Vector2 m_velocityJitter = new Vector2(-20f, 20f); // min and max

    // Start is called before the first frame update
    void Start()
    {
        m_id = ParticlePlugin.CreateParticleSystem();
        m_particleDataArray = new ParticlePlugin.ParticleData[m_maxParticleObjects];    // Initialize the array
        // Initialize particle objects
        for (int i = 0; i < m_maxParticleObjects; i++)
        {
            m_particleDataArray[i] = new ParticlePlugin.ParticleData();
            GameObject particle = GameObject.Instantiate(m_particlePrefab, this.transform);
            m_particles.Add(particle);
            particle.SetActive(false);
        }
        UpdateParameters();
    }

    // Update is called once per frame
    void Update()
    {
        ParticlePlugin.UpdateParticleSystem(m_id, Time.deltaTime);
        int num = ParticlePlugin.GetParticleNum(m_id);
        ParticlePlugin.GetParticleData(m_particleDataArray, m_id, num);
        for (int i = 0; i < m_maxParticleObjects; ++i)
        {
            GameObject particle = m_particles[i];
            if ((i < num && !m_particleDataArray[i].isAlive) || i >= num) 
            {
                particle.SetActive(false);
            }
            else
            {
                particle.SetActive(true);
                // Set position
                particle.transform.position = ParticlePlugin.FloatArrayToVector3(m_particleDataArray[i].position);
                // Set scale
                float scale = m_particleDataArray[i].scale;
                particle.transform.localScale = new Vector3(scale, scale, scale);
                // Set color
                Vector3 color = ParticlePlugin.FloatArrayToVector3(m_particleDataArray[i].color);
                float alpha = m_particleDataArray[i].alpha;
                Material material = particle.GetComponent<Renderer>().material;
                material.color = new Color(color.x, color.y, color.z, alpha);
            }
        }
    }

    void OnValidate()
    {
        try { if (Application.isPlaying) { UpdateParameters(); } }
        catch (DllNotFoundException) { }
    }

    void UpdateParameters()
    {
        ParticlePlugin.ParticleSystemParameter parm = new ParticlePlugin.ParticleSystemParameter(m_infinite,
            m_maxParticles, m_lifeTime, m_gravity, m_position, m_velocity, m_startColor, m_endColor,
            m_startAlpha, m_startSize, m_endSize, m_positionJitter, m_velocityJitter);
        ParticlePlugin.SetParticleSystemParameters(m_id, parm);
    }

    public bool IsAlive()
    {
        return ParticlePlugin.IsAlive(m_id);
    }

    public void Destroy()
    {
        ParticlePlugin.RemoveParticleSystem(m_id);
        Destroy(this.gameObject);
    }
}
