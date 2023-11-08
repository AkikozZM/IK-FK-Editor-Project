using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System.Runtime.InteropServices;

public class ParticlePlugin
{
    [StructLayout(LayoutKind.Sequential, Pack = 1)]
    public unsafe struct ParticleData
    {
        public ParticleData(int index)
        {
            isAlive = false;
            position = new float[3];
            color = new float[3];
            scale = 0;
            alpha = 1;
        }
        public bool isAlive;
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3)]
        public float[] position;
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3)]
        public float[] color;
        public float scale;
        public float alpha;
    }

    [StructLayout(LayoutKind.Sequential, Pack = 1)]
    public unsafe struct ParticleSystemParameter
    {
        public ParticleSystemParameter(bool _infinite, int _maxParticles, float _lifeTime, Vector3 _gravity,
            Vector3 _position, Vector3 _velocity, Color _startColor, Color _endColor, float _startAlpha, 
            float _startSize, float _endSize, Vector2 _positionJitter, Vector2 _velocityJitter)
        {
            infinite = _infinite;
            maxParticles = _maxParticles;
            lifeTime = _lifeTime;
            gravity = new float[3] { _gravity[0], _gravity[1], _gravity[2] };
            position = new float[3] { _position[0], _position[1], _position[2] };
            velocity = new float[3] { _velocity[0], _velocity[1], _velocity[2] };
            startColor = new float[3] { _startColor.r, _startColor.g, _startColor.b };
            endColor = new float[3] { _endColor.r, _endColor.g, _endColor.b };
            startAlpha = _startAlpha;
            startSize = _startSize;
            endSize = _endSize;
            positionJitter = new float[2] { _positionJitter[0], _positionJitter[1] };
            velocityJitter = new float[2] { _velocityJitter[0], _velocityJitter[1] };
        }
        public bool infinite;
        public int maxParticles;
        public float lifeTime;
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3)]
        public float[] gravity;
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3)]
        public float[] position;
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3)]
        public float[] velocity;
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3)]
        public float[] startColor;
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 3)]
        public float[] endColor;
        public float startAlpha;
        public float startSize;
        public float endSize;
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 2)]
        public float[] positionJitter;
        [MarshalAs(UnmanagedType.ByValArray, SizeConst = 2)]
        public float[] velocityJitter;
    }

    [DllImport("ParticlePlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern int CreateParticleSystem();

    [DllImport("ParticlePlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void RemoveParticleSystem(int id);

    [DllImport("ParticlePlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern int GetParticleNum(int id);

    [DllImport("ParticlePlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern int GetMaxParticleNum(int id);

    [DllImport("ParticlePlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void UpdateParticleSystem(int id, float deltaT);

    [DllImport("ParticlePlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void SetParticleSystemParameters(int id, ParticleSystemParameter parm);

    [DllImport("ParticlePlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern void GetParticleData([In, Out] ParticleData[] particleDataArray, int id, int size);

    [DllImport("ParticlePlugin", CallingConvention = CallingConvention.Cdecl)]
    public static extern bool IsAlive(int id);

    public static Vector3 FloatArrayToVector3(float[] vec3)
    {
        if (vec3.Length != 3) { Debug.LogError("Wrong array size"); }
        return new Vector3(vec3[0], vec3[1], vec3[2]);
    }


}

    public class ParticlePluginManager : MonoBehaviour
{
    // Start is called before the first frame update
    void Start()
    {
        
    }

    // Update is called once per frame
    void Update()
    {
        
    }
}
