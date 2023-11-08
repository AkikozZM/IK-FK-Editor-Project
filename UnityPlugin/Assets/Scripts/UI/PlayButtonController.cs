using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class PlayButtonController : MonoBehaviour
{
    private TMPro.TMP_Text m_text;
    // Start is called before the first frame update
    void Start()
    {
        m_text = this.GetComponentInChildren<TMPro.TMP_Text>();
    }

    // Update is called once per frame
    void Update()
    {
        
    }

    public void ChangeText()
    {
        if (m_text.text.Equals("Play")) { m_text.text = "Pause"; }
        else { m_text.text = "Play"; }
    }

    public void Pause()
    {
        m_text.text = "Play";
    }
}
