using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using SimpleFileBrowser;

public class BVHFileBrowserController : MonoBehaviour
{

    public FKIKCharacterController m_jointController;

    public FKIKUIController m_uiController;
    // Start is called before the first frame update
    void Start()
    {
        FileBrowser.SetFilters(false, new FileBrowser.Filter("BVH", ".bvh"));


    }

    // Update is called once per frame
    void Update()
    {
        
    }

    public void ShowFileBrowser()
    {
        if (FileBrowser.IsOpen)
            return;

        string defaultPath = Application.dataPath + "/../../motions";
        // Show a select folder dialog 
        // onSuccess event: print the selected folder's path
        // onCancel event: print "Canceled"
        // Load file/folder: folder, Allow multiple selection: false
        // Initial path: default (Documents), Title: "Select Folder", submit button text: "Select"
        FileBrowser.ShowLoadDialog(SelectFile, () => { Debug.Log("Canceled"); },
                                   false, false, defaultPath, "Load", "Select");
    }

    void SelectFile(string[] paths)
    {
        Debug.Log(paths[0]);
        if (!m_jointController.LoadBVHFile(paths[0]))
        {
            Debug.LogError("BVH file does not match!");
            return;
        }
        m_uiController.SetupProgressBar();
    }
}
