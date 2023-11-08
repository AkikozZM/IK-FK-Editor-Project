using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;

[System.Serializable]
public class Gain
{
    public float gMaxSpeed = 1000.0f;
    public float gMaxAngularSpeed = 200.0f;
    public float gMaxForce = 2000.0f;
    public float gMaxTorque = 2000.0f;
    public float gKNeighborhood = 500.0f;
    public float gOriKv = 32.0f;
    public float gOriKp = 256.0f;
    public float gVelKv = 10.0f;
    public float gAgentRadius = 80.0f;

    public float gMass = 1.0f;
    public float gIntertia = 1.0f;
    public float KArrival = 1.0f;
    public float KDepartucre = 12000.0f;
    public float KNoise = 15.0f;
    public float KWander = 80.0f;
    public float KAvoid =  600.0f;
    public float TAvoid = 1000.0f;
    public float KSeparation = 12000.0f;
    public float KAlignment = 1.0f;
    public float KCohesion = 1.0f;
    public float leaderFollowDistance = 200.0f;

    public float[] GetGainFloatArray(int gainNum)
    {
        float[] gains = new float[gainNum];
        gains[0] = gMaxSpeed;
        gains[1] = gMaxAngularSpeed;
        gains[2] = gMaxForce;
        gains[3] = gMaxTorque;
        gains[4] = gKNeighborhood;
        gains[5] = gOriKv;
        gains[6] = gOriKp;
        gains[7] = gVelKv;
        gains[8] = gAgentRadius;

        gains[9] = gMass;
        gains[10] = gIntertia;
        gains[11] = KArrival;
        gains[12] = KDepartucre;
        gains[13] = KNoise;
        gains[14] = KWander;
        gains[15] = KAvoid;
        gains[16] = TAvoid;
        gains[17] = KSeparation;
        gains[18] = KAlignment;
        gains[19] = KCohesion;
        gains[20] = leaderFollowDistance;
        return gains;
    }
}

public class BehaviorPluginManager : MonoBehaviour
{
    public enum BehaviorType { SEEK, FLEE, ARRIVAL, DEPARTURE, AVOID, WANDER,
        ALIGNMENT, SEPARATION, COHESION, FLOCKING, LEADER };

    private int m_actorNum;
    
    public BehaviorType m_behaviorType = BehaviorType.SEEK;

    public List<GameObject> m_actorList = new List<GameObject>();
    private List<GameObject> m_obstacleList = new List<GameObject>();

    public Gain m_gain;
    public GameObject m_target;
    public GameObject m_leader;
    public GameObject m_leaderIndicator;

    [Header("Obstacle")]
    public int m_obstacleNum = 25;
    public GameObject m_obstaclePrefab;
    public Vector2 m_obstaclePosMinMax = new Vector2(-100, 100);
    public Vector2 m_obstacleRadiusMinMax = new Vector2(10.0f, 20.0f);
    public LayerMask m_environmentLayer;

    // Data Array
    private BehaviorPlugin.ActorData[] m_actorDataArray;
    private BehaviorPlugin.ObstacleData[] m_obstacleDataArray;

    // Start is called before the first frame update
    void Start()
    {
        //ResetPlugin();
        SwitchBehavior((int)m_behaviorType);
    }

    void OnValidate()
    {
        try
        {
            if (Application.isPlaying)
            {
                ResetPlugin() ;
            }
        }
        catch (System.DllNotFoundException) { }
    }

    public void ResetPlugin()
    {
        m_actorNum = m_actorList.Count;
        BehaviorPlugin.InitializeBehaviorPlugin(m_actorNum, m_obstacleNum, (int)m_behaviorType);
        int gainNum = BehaviorPlugin.GetGainNum();
        BehaviorPlugin.SetControllerGains(m_gain.GetGainFloatArray(gainNum), gainNum);

        m_actorDataArray = new BehaviorPlugin.ActorData[m_actorNum];
        for (int i = 0; i < m_actorNum; ++i)
        {
            m_actorDataArray[i] = new BehaviorPlugin.ActorData(i);
            m_actorDataArray[i].globalPosition = PluginHelpFunction.Vector3ToFloatArray(m_actorList[i].transform.position);
            m_actorDataArray[i].globalRotation = PluginHelpFunction.QuaternionToFloatArray(m_actorList[i].transform.rotation);
        }
        SetLeader(m_leader);
        ClearObstacles();
        if (m_behaviorType == BehaviorType.AVOID)
        {
            GenerateObstacles();
        }
    }

    public void GenerateObstacles()
    {
        Vector3[] positions = new Vector3[] { new Vector3(500, 0, 500) ,
            new Vector3(500, 0, -500), new Vector3(-500, 0, -500), new Vector3(-500, 0, 500)};
        ClearObstacles();
        m_obstacleDataArray = new BehaviorPlugin.ObstacleData[m_obstacleNum];
        for (int i = 0; i < m_obstacleNum; ++i)
        {
            // Generates random obstacles
            GameObject obstacle = GameObject.Instantiate(m_obstaclePrefab);
            Vector3 position;
            if (i <= 3)
            {
                position = positions[i];
            }
            else
            {
                position = new Vector3(Random.Range(m_obstaclePosMinMax[0], m_obstaclePosMinMax[1]),
                0, Random.Range(m_obstaclePosMinMax[0], m_obstaclePosMinMax[1]));
            }
            Physics.Raycast(position + new Vector3(0, 1000, 0), Vector3.down, out RaycastHit hit, Mathf.Infinity, m_environmentLayer);
            position.y = hit.point.y;

            float randomRadius = Random.Range(m_obstacleRadiusMinMax[0], m_obstacleRadiusMinMax[1]);
            Vector3 scale = new Vector3(randomRadius, randomRadius, randomRadius);
            obstacle.transform.position = position;
            obstacle.transform.localScale = scale;
            obstacle.transform.rotation = Random.rotation;
            m_obstacleList.Add(obstacle);

            // Set obstacle data
            m_obstacleDataArray[i].globalPosition = PluginHelpFunction.Vector3ToFloatArray(position);
            m_obstacleDataArray[i].radius = randomRadius;
        }
        BehaviorPlugin.SetObstacleData(m_obstacleDataArray, m_obstacleNum);

    }

    public void ClearObstacles()
    {
        BehaviorPlugin.ClearObstacles();
        foreach(GameObject obstacle in m_obstacleList)
        {
            Destroy(obstacle);
        }
        m_obstacleList.Clear();
    }

    private void Update()
    {
        if (m_leader)
        {
            m_leaderIndicator.SetActive(true);
            FKIKCharacterController character = m_leader.GetComponent<FKIKCharacterController>();
            if (!character) { Debug.LogError("No character controller!"); }
            Vector3 leaderPos = character.m_root.transform.position;
            m_leaderIndicator.transform.position = new Vector3(leaderPos.x, leaderPos.y + 100.0f, leaderPos.z);
        }
        else
        {
            m_leaderIndicator.SetActive(false);
        }
    }

    private void FixedUpdate()
    {
        //1.  send position of target/goal to Plugin
        Vector3 targetPos = m_target.transform.position;
        BehaviorPlugin.SetTarget(new float[] { targetPos[0], targetPos[1], targetPos[2] });

        //2. transfer all agent data necessary for force and torque calculations from Unity to Plugin
        //3. update agent position and orientation
        //4. transfer all data necessary for Unity agent update from Plugin 
        WriteActorDataArray();
        BehaviorPlugin.UpdateActorData(Time.fixedDeltaTime, m_actorDataArray, m_actorNum);
        WriteActorTransform();
    }


    //Copy transform and rigidbody data to m_actorDataArray
    public void WriteActorDataArray()
    {
        for (int i = 0; i < m_actorNum; ++i)
        {
            ref BehaviorPlugin.ActorData actorData = ref m_actorDataArray[i];

            // Get the leader's root joint transform because its controlled by keyboard
            if (m_actorList[i] == m_leader)
            {
                GameObject root = m_leader.GetComponent<FKIKCharacterController>().m_root;
                actorData.globalPosition = PluginHelpFunction.Vector3ToFloatArray(root.transform.position);
                actorData.globalRotation = PluginHelpFunction.QuaternionToFloatArray(root.transform.rotation);
            }
        }
    }

    // Copy transform to rigidbody data
    public void WriteActorTransform()
    {
        for (int i = 0; i < m_actorNum; ++i)
        {
            ref BehaviorPlugin.ActorData actorData = ref m_actorDataArray[i];
            if (m_actorList[i] != m_leader)
            {
                GameObject character = m_actorList[i];
                character.transform.position = PluginHelpFunction.FloatArrayToVector3(actorData.globalPosition);
                character.transform.rotation = PluginHelpFunction.FloatArrayToQuaternion(actorData.globalRotation);
                //GameObject root = m_actorList[i].GetComponent<FKIKCharacterController>().m_root;
                //root.transform.position = PluginHelpFunction.FloatArrayToVector3(actorData.globalPosition);
                //root.transform.rotation = PluginHelpFunction.FloatArrayToQuaternion(actorData.globalRotation);
            }
        }
    }

    public void SetLeader(GameObject actor)
    {
        int index = m_actorList.IndexOf(actor);
        if (index == -1) { return; }
        BehaviorPlugin.SetLeaderIndex(index);
        m_leader = actor;
    }

    public GameObject GetLeader()
    {
        return m_leader;
    }

    public void SwitchBehavior(int num)
    {
        BehaviorType newType = (BehaviorType)num;
        m_behaviorType = newType;
        ResetPlugin();
        if (newType != BehaviorType.LEADER && m_leader)
        {
            m_leader.GetComponent<FKIKCharacterController>().SwitchControlMode(FKIKCharacterController.ControlMode.BEHAVIOR);
            m_leader.GetComponent<FKIKSpeedScaleController>().enabled = true;
            m_leader = null;
        }
        else if (newType == BehaviorType.LEADER)
        {
            m_leader = m_actorList[0];
            m_leader.GetComponent<FKIKCharacterController>().SwitchControlMode(FKIKCharacterController.ControlMode.KEYBOARD);
            m_leader.GetComponent<FKIKSpeedScaleController>().enabled = false;
        }
    }
}
