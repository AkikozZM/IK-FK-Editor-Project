using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class BehaviorUIController : MonoBehaviour
{
    public BehaviorPluginManager m_behaviorPluginManager;
    public Toggle m_obstacleToggle;
    // Start is called before the first frame update
    void Start()
    {
        m_obstacleToggle.isOn = false;
    }

    // Update is called once per frame
    void Update()
    {
        
    }

    public void SetShowObstacle()
    {
        if (m_obstacleToggle.isOn)
        {
            m_behaviorPluginManager.GenerateObstacles();
        }
        else
        {
            m_behaviorPluginManager.ClearObstacles();
        }
    }
}
