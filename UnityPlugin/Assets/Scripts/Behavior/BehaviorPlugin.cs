using System.Collections;
using System.Collections.Generic;
using System.Runtime.InteropServices;
using UnityEngine;

public class BehaviorPlugin
{
    // Initialize the plugin
    [DllImport("BehaviorPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void InitializeBehaviorPlugin(int agentNum, int obstacleNum, int activeBehavior);

    // Set the target position
    [DllImport("BehaviorPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void SetTarget(float[] target);

    // Set all obstacle data
    [DllImport("BehaviorPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void SetObstacleData([In, Out] ObstacleData[] obstacleDataArray, int obstacleNum);

    // Set all obstacle data
    [DllImport("BehaviorPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void ClearObstacles();

    // Set all gains 
    [DllImport("BehaviorPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern int GetGainNum();

    // Set all gains 
    [DllImport("BehaviorPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void SetControllerGains(float[] gains, int num);

    // Update actor states
    [DllImport("BehaviorPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void UpdateActorData(float timestep, [In, Out] ActorData[] ActorDataArray, int actorNum);

    // Set the index of the leader actor
    [DllImport("BehaviorPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void SetLeaderIndex(int index);

    // Get the index of the leader actor
    [DllImport("BehaviorPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern int GetLeaderIndex();

    [StructLayout(LayoutKind.Sequential, Pack = 1)]
    public unsafe struct ActorData
    {
        public ActorData(int index)
        {
            id = index;
            globalPosition = new float[3] { 0, 0, 0 };
            globalRotation = new float[4] { 1, 0, 0, 0 };
            linearVelocity = new float[3] { 0, 0, 0 };
            angularVelocity = new float[3] { 0, 0, 0 };
            globalInertialTensor = new float[4] { 1, 0, 0, 0 };
            mass = 1;
        }
        public int id;
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3)]
        public float[] globalPosition;
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 4)]
        public float[] globalRotation;      // quat (w, x, y, z)
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3)]
        public float[] linearVelocity;
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3)]
        public float[] angularVelocity;
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 4)]
        public float[] globalInertialTensor;
        public float mass;
    }


    [StructLayout(LayoutKind.Sequential, Pack = 1)]
    public unsafe struct ObstacleData
    {
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3)]
        public float[] globalPosition;
        public float radius;
    }

}
