using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class TargetCurveController : MonoBehaviour
{
    public CurveController m_rootCurve;
    public FKIKCharacterController m_character;
    public GameObject m_target;
    public float m_threshold = 150;
    private float m_curveT = 0;

    // Start is called before the first frame update
    void Start()
    {
    }

    // Update is called once per frame
    void Update()
    {
        if (m_rootCurve.GetKeyNumber() < 2) { return; }
        Vector3 rootPos = m_character.m_root.transform.position;
        rootPos.y = 0;

        m_rootCurve.GetCurveValue(m_curveT, out Vector3 curvePos, out Quaternion curveQuat);
        Vector3 targetPos = curvePos;
        m_target.transform.position = targetPos;
        // Rotate the root towards the root target on the curve
        if (Vector3.Distance(rootPos, targetPos) < m_threshold && m_curveT <= m_rootCurve.GetDuration())
        {
            m_curveT += Time.deltaTime;
        }
    }

    public void ResetCurve()
    {
        m_curveT = 0;
        m_rootCurve.ClearCurve();
    }
}
