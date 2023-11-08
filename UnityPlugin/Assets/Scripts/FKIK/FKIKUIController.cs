using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

//public class FKIKSceneGlobal
//{
//    public enum Mode { FK, IK };
//    public static Mode mode = Mode.FK;
//}

public class FKIKUIController : MonoBehaviour
{

    private JointPainter m_jointPainter;
    private FKIKCharacterController m_characterController;
    public GameObject m_model;
    public VirtualJointGizmoController m_gizmoController;

    public Slider m_progressBar;
    public TMPro.TMP_InputField m_progressBarValue;
    public PlayButtonController m_playButtonController;

    public Toggle m_rootCurveToggle;
    public Toggle m_FKToggle;

    public TargetController m_targetController;
    public TargetCurveController m_targetCurveController;

    // Start is called before the first frame update
    void Start()
    {
        m_jointPainter = m_model.GetComponent<JointPainter>();
        m_characterController = m_model.GetComponent<FKIKCharacterController>();
    }

    // Update is called once per frame
    void Update()
    {
        //if (m_jointController.BVHLoaded())
        //{
        //    FKIKSceneGlobal.mode = FKIKSceneGlobal.Mode.FK;
        //    if (m_FKToggle) { m_FKToggle.isOn = true; }
        //}
    }

    public void SetupProgressBar()
    {
        if (!m_progressBar) { return; }
        m_progressBar.maxValue = m_characterController.GetKeySize() - 1;
        m_progressBar.value = 0;
    }

    public void SetProgressBarValue(float value)
    {
        m_progressBarValue.text = value.ToString();
    }

    public void SetFK()
    {
        m_jointPainter.SetDrawTargets(false);
        m_gizmoController.ChangeSelectedJoint(null);
        m_gizmoController.SetFKIKMode(VirtualJointGizmoController.Mode.FK);
        m_jointPainter.ResetTargetPosition();
    }

    public void SetIK()
    {
        m_jointPainter.SetDrawTargets(true);
        m_gizmoController.ChangeSelectedJoint(null);
        m_gizmoController.SetFKIKMode(VirtualJointGizmoController.Mode.IK);
        m_jointPainter.ResetTargetPosition();
    }


    public void OnClickPlayButton()
    {
        bool move = m_characterController.StartMove();
        if (move) { m_playButtonController.ChangeText(); }
    }
    
    public void SetTargetControl(int mode)
    {
        
        if (mode == 0)  // Control the target by mouse button
        {
            m_targetCurveController.ResetCurve();
            m_targetController.enabled = true;
            m_targetCurveController.enabled = false;
        }
        else
        {
            m_targetController.enabled = false;
            m_targetCurveController.enabled = true;
        }
    }


}
