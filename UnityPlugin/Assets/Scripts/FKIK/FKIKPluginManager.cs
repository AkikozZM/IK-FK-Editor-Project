using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Runtime.InteropServices;
using System;
using System.Text;

public class FKIKPlugin
{
    // Return the actor ID
    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern int CreateActor();

    // Remove the actor
    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void RemoveActor(int id);

    // Return the joint ID
    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern int CreateJoint(int id, string name, bool isRoot);

    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void SetJointData(int id, JointData data, int parentID);

    // Update all joints' global transform
    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void UpdateSkeleton(int id);

    // Update IK joints
    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void UpdateIKSkeleton(int id);

    // Use absolute file path
    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern bool LoadBVH(int id, string filepath);

    // Get the number of joints
    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern int GetJointSize(int id);

    // Return the id of the joint by its name
    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern int GetJointIdByName(int id, string name);

    // Return the id of the end joint by its parent name.
    // This function is used to get the end joint id 
    // beacause in .BVH files the end joint does not have name
    // If the joint is not the end joint, it will return -1
    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern int GetJointIdByParentName(int id, string name);

    // Get Joint Data in each frame
    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern int GetJointData(int id, [In, Out] JointData[] jointDataArray, int size);

    // Update BVH at time t
    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern int UpdateBVHSkeleton(int id, float t);

    // Get the duration of the bvh file
    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern float GetDuration(int id);

    // Get the number of keys
    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern int GetKeySize(int id);

    // Get time of the key
    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern float GetKeyTime(int id, int keyID);

    // Set joint rotation
    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void SetJointRotation(int id, int jointID, float[] q);

    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void SetRootJointTranslation(int id, float[] pos);


    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void SetRootJointRotation(int id, float[] q);

    /* IK */
    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void SolveLimbIK(int id, int jointID, float[] pos);

    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void SetLeftHandID(int id, int jointID);

    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void SetRightHandID(int id, int jointID);

    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void SetLeftFootID(int id, int jointID);

    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void SetRightFootID(int id, int jointID);

    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void SetRootID(int id, int jointID);

    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void CreateLimbIKChains(int id);

    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void SolveFootIK(int id, float leftHeight, float rightHeight, 
        bool leftRotate, bool rightRotate, float[] leftNormal, float[] rightNormal);

    [DllImport("FKIKPlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void UpdateGuideJointByTarget(int id, float[] targetPos, float[] newPos, float[] newQuat);

    [StructLayout(LayoutKind.Sequential, Pack = 1)]
    public unsafe struct JointData
    {
        public JointData(int index)
        {
            id = index;
            localRotation = new float[4];   // Right-hand rotation (w, x, y, z)
            localTranslation = new float[3];    // Right-hand translation
        }
        public int id;
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 4)]
        public float[] localRotation;      // Quaternion
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3)]
        public float[] localTranslation;    // Vector3        
    }


    // Convert between the right-hand coordinate system (BVH) and the left-hand coordinate system (Unity)
    public static Vector3 BVHToUnityTranslation(float[] vec3)
    {
        if (vec3.Length != 3)
        {
            Debug.LogError("Wrong array size");
        }
        return new Vector3(-vec3[0], vec3[1], vec3[2]);
    }

    public static Quaternion BVHToUnityQuaternion(float[] vec4)
    {
        if (vec4.Length != 4)
        {
            Debug.LogError("Wrong array size");
        }
        return new Quaternion(-vec4[1], vec4[2], vec4[3], -vec4[0]);
    }

    public static float[] UnityToBVHQuaternion(Quaternion quat)
    {
        float[] q = new float[] {-quat.w, -quat.x, quat.y, quat.z };
        return q;
    }

    public static float[] UnityToBVHTranslation(Vector3 vec)
    {
        float[] v = new float[] { -vec.x, vec.y, vec.z };
        return v;
    }

}

//public class FKIKPluginManager : MonoBehaviour
//{
//    // Start is called before the first frame update
//    void Start()
//    {
//    }

//    // Update is called once per frame
//    void Update()
//    {
        
//    }
//}
