/*
 * MGUIEFileSelector.cxx
 *
 *
 * Copyright (C) by Andreas Zoglauer.
 * All rights reserved.
 *
 *
 * This code implementation is the intellectual property of
 * Andreas Zoglauer.
 *
 * By copying, distributing or modifying the Program (or any work
 * based on the Program) you indicate your acceptance of this statement,
 * and all its terms.
 *
 */


////////////////////////////////////////////////////////////////////////////////
//
// MGUIEFileSelector
//
//
// This class is an elementary GUI-widget:
// It contains a text-label and a input-field.
// It can be checked, if the input is within a preselected range of values.
//
//
// Example:
//
// 
//
////////////////////////////////////////////////////////////////////////////////


// Include the header:
#include "MGUIEFileSelector.h"

// Standard libs:
#include "MStreams.h"

// ROOT libs:
#include "TGMsgBox.h"
#include "TSystem.h"
#include "TGFileDialog.h"

// MEGAlib libs:
#include "MFile.h"


////////////////////////////////////////////////////////////////////////////////


#ifdef ___CLING___
ClassImp(MGUIEFileSelector)
#endif


////////////////////////////////////////////////////////////////////////////////


MGUIEFileSelector::MGUIEFileSelector(const TGWindow* Parent, MString Label, 
                                     MString FileName) :
  MGUIElement(Parent, kVerticalFrame)
{
  // Creates a frame containing two entry-boxes 
  //
  // Label:    The text of the label
  // FileName: Default name of the file

  m_Label = Label;
  m_FileName = FileName;
  m_FileTypes = new const char*[4];
  m_FileTypes[0] = "All files";
  m_FileTypes[1] = "*";
  m_FileTypes[2] = 0;
  m_FileTypes[3] = 0;
  m_NFileTypes = 0;

  Create();
}


////////////////////////////////////////////////////////////////////////////////


MGUIEFileSelector::~MGUIEFileSelector()
{
  // default destructor

  if (MustCleanup() == kNoCleanup) {
    delete m_TextLabel;
    delete m_TextLabelLayout;

    delete m_InputLayout;
    delete m_Input;  

    delete m_ButtonFolderLayout;
    delete m_ButtonFolder;

    delete m_InputFrame;
    delete m_InputFrameLayout;
  }
  
  // We do only delete the first 2*m_NFileTypes. The ones added in the constructor have not been created by new!
  for (unsigned int f = 0; f < 2*m_NFileTypes; ++f) {
    delete [] m_FileTypes[f];
  }
  delete [] m_FileTypes;
}


////////////////////////////////////////////////////////////////////////////////


void MGUIEFileSelector::Create()
{
  // Create the label and the input-field.

  // label:
  if (m_Label != "") {
    m_TextLabel = new TGLabel(this, new TGString(m_Label));
    m_TextLabelLayout =
      new TGLayoutHints(kLHintsLeft | kLHintsTop | kLHintsCenterY, 0, 30, 0, 0);
    AddFrame(m_TextLabel, m_TextLabelLayout);
  } else {
    m_TextLabel = 0;
    m_TextLabelLayout = 0;
  }

  m_InputFrame = new TGHorizontalFrame(this, 100, 40);
  if (m_Label != "") {
    m_InputFrameLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 5, 0);
  } else {
    m_InputFrameLayout = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX, 0, 0, 0, 0);
  }
  AddFrame(m_InputFrame, m_InputFrameLayout);

  // input-field
  m_InputLayout = new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 0, 0, 0, 0);
  m_Input = new TGTextEntry(m_InputFrame, m_FileName);
  m_Input->SetAlignment(kTextRight);
  m_InputFrame->AddFrame(m_Input, m_InputLayout);
  m_Input->Layout();

  // input label:
  MString Icon = g_MEGAlibPath + "/resource/icons/folder.xpm";
  MFile::ExpandFileName(Icon);

  m_ButtonFolderLayout = new TGLayoutHints(kLHintsRight | kLHintsCenterY, 5, 0, 0, 0);
  m_ButtonFolder = new TGPictureButton(m_InputFrame, fClient->GetPicture(Icon), 99);
  m_ButtonFolder->Associate(this);
  m_InputFrame->AddFrame(m_ButtonFolder, m_ButtonFolderLayout);


  return;
}


////////////////////////////////////////////////////////////////////////////////


MString MGUIEFileSelector::GetFileName()
{
  // Return the name of the selected file

  return m_Input->GetText();
}


////////////////////////////////////////////////////////////////////////////////


void MGUIEFileSelector::SetFileName(MString FileName)
{
  // Set the name of the selected file

  m_Input->SetText(FileName);
}


////////////////////////////////////////////////////////////////////////////////


void MGUIEFileSelector::SetFileType(MString Name, MString Suffix)
{
  // Add a new file type to the back of the list
  
  ++m_NFileTypes;
  const char** Types = new const char*[4 + 2*m_NFileTypes]; 
 
  // Copy the old added content
  for (unsigned int f = 0; f < 2*(m_NFileTypes-1); ++f) {
    Types[f] = m_FileTypes[f];
  }

  // Add the new content
  Types[2*m_NFileTypes-2] = StrDup(Name.Data());
  Types[2*m_NFileTypes-1] = StrDup(Suffix.Data());
  
  // Add the bottom (all files)
  for (unsigned int f = 2*m_NFileTypes; f < 4 + 2*m_NFileTypes; ++f) {
    Types[f] = m_FileTypes[f-2];
  }
  
  // Cleanup and copy
  delete [] m_FileTypes; // Although we have a char** we do not want to delete the strings, since we copies them! 
  m_FileTypes = Types;
}


////////////////////////////////////////////////////////////////////////////////


void MGUIEFileSelector::SetEnabled(bool flag)
{
  m_IsEnabled = flag;
  
  m_Input->SetEnabled(m_IsEnabled);

  if (flag == true) {
    m_ButtonFolder->SetState(kButtonUp);
  } else {
    m_ButtonFolder->SetState(kButtonDisabled);
  }
}


////////////////////////////////////////////////////////////////////////////////


bool MGUIEFileSelector::ProcessMessage(long Message, long Parameter1, 
                                       long Parameter2)
{
  // Process the messages for this application
  
  switch (GET_MSG(Message)) {
  case kC_COMMAND:
    switch (GET_SUBMSG(Message)) {
    case kCM_BUTTON:
      switch (Parameter1) {
      case 99:
        SelectFileName();
        break;
      default:
        break;
      }
    default:
      break;
    }
  default:
    break;
  }

  return true;
}


////////////////////////////////////////////////////////////////////////////////


bool MGUIEFileSelector::SelectFileName()
{
  TGFileInfo Info;
  Info.fFilename = StrDup(gSystem->BaseName(m_FileName));
  Info.fIniDir = StrDup(gSystem->DirName(m_FileName));
  Info.fFileTypes = m_FileTypes;
  
  new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &Info);
  
  // Get the filename ...
  if ((char *) Info.fFilename != 0) {
    m_FileName = MString((char *) Info.fFilename);
    if (m_FileName.IsEmpty()) {
      return false;
    }
  } 
  // ... or return when cancel has been pressed
  else {
    return false;
  }  
  
  m_Input->SetText(m_FileName);
  Layout();
  
  return true;
}


// MGUIEFileSelector.cxx: the end...
////////////////////////////////////////////////////////////////////////////////
