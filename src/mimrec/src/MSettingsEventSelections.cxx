/*
 * MSettingsEventSelections.cxx
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
// MSettingsEventSelections.cxx
//
////////////////////////////////////////////////////////////////////////////////


// Include the header:
#include "MSettingsEventSelections.h"

// Standard libs:
#include <limits>
using namespace std;
#include <iomanip>

// ROOT libs:

// MEGAlib libs:
#include "MStreams.h"
#include "MEarthHorizon.h"
#include "MProjection.h"
#include "MLMLAlgorithms.h"


////////////////////////////////////////////////////////////////////////////////


#ifdef ___CINT___
ClassImp(MSettingsEventSelections)
#endif


///////////////////////////////////////////////////////////////////////////////


MSettingsEventSelections::MSettingsEventSelections() : MSettingsInterface()
{
  // default constructor

  m_EventSelectionModified = true;


  // Event selector
  m_EventSelectorTab = 0;

  m_EventIdRangeMin = 0;
  m_EventIdRangeMax = numeric_limits<int>::max();

  m_TrackLengthRangeMin = 1;
  m_TrackLengthRangeMax = 100;

  m_SequenceLengthRangeMin = 2;
  m_SequenceLengthRangeMax = 10;

  m_ClusteringQualityFactorRangeMin = 0;
  m_ClusteringQualityFactorRangeMax = 100;

  m_ComptonQualityFactorRangeMin = 0;
  m_ComptonQualityFactorRangeMax = 100;

  m_TrackQualityFactorRangeMin = 0;
  m_TrackQualityFactorRangeMax = 100;

  m_FirstEnergyRangeMin = 0;
  m_FirstEnergyRangeMax = 2000;

  m_SecondEnergyRangeMin = 0;
  m_SecondEnergyRangeMax = 0;

  m_ThirdEnergyRangeMin = 0;
  m_ThirdEnergyRangeMax = 0;

  m_FourthEnergyRangeMin = 0;
  m_FourthEnergyRangeMax = 0;

  m_EnergyRangeGammaMin = 0;
  m_EnergyRangeGammaMax = 100000;

  m_EnergyRangeElectronMin = 0;
  m_EnergyRangeElectronMax = 100000;

  m_ComptonAngleRangeMin = 0;
  m_ComptonAngleRangeMax = 180;

  m_EHCType = 0;
  m_EHCProbability = 1.0;
  m_EHCComptonProbabilityFileName = g_StringNotDefined;
  m_EHCPairProbabilityFileName = g_StringNotDefined;
  m_EHCEarthPosition = MVector(0, 0, -1E+20);
  m_EHCAngle = 0;

  m_SourceUsePointSource = false;
  m_SourceCoordinates = MProjection::c_Galactic;
  m_SourcePhi = 0.0;
  m_SourceTheta = 0.0;
  m_SourceLongitude = 184.56;
  m_SourceLatitude = -5.78;
  m_SourceX = 0.0;
  m_SourceY = 0.0; 
  m_SourceZ = 100.0;
  m_SourceARMMin = 0.0;
  m_SourceARMMax = 180.0;
  m_SourceSPDMin = 0.0;
  m_SourceSPDMax = 180.0;

  m_BeamUse = false;
  m_BeamStartX = 0.0;
  m_BeamStartY = 0.0; 
  m_BeamStartZ = 10000.0;
  m_BeamFocalSpotX = 0.0;
  m_BeamFocalSpotY = 0.0; 
  m_BeamFocalSpotZ = 0.0;
  m_BeamRadius = 1000.0;
  m_BeamDepth = 1000.0;

  m_ThetaDeviationMax = 180;

  m_FirstDistanceRangeMin = 0;
  m_FirstDistanceRangeMax = numeric_limits<int>::max();

  m_DistanceRangeMin = 0;
  m_DistanceRangeMax = numeric_limits<int>::max();

  m_TimeRangeMin = 0.0;
  m_TimeRangeMax = double(numeric_limits<int>::max());

  m_TimeWalkRangeMin = -1000;
  m_TimeWalkRangeMax = numeric_limits<int>::max();

  m_CoincidenceWindowRangeMin = 0;
  m_CoincidenceWindowRangeMax = 1;

  m_OpeningAnglePairMin = 0;
  m_OpeningAnglePairMax = 180;

  m_InitialEnergyDepositPairMin = 0;
  m_InitialEnergyDepositPairMax = 10000;

  m_EventTypeCompton = 1;
  m_EventTypeDoubleCompton = 1;
  m_EventTypeComptonNotTracked = 1;
  m_EventTypeComptonTracked = 1;
  m_EventTypePair = 1;
  m_EventTypePhoto = 0;
  m_EventTypeUnidentifiable = 0;
  m_EventTypeDecay = 0;

  m_FlaggedAsBad = false;

  m_ExcludedDetectors.clear();

  // The special GUI mode
  m_SpecialMode = false;
}


////////////////////////////////////////////////////////////////////////////////


MSettingsEventSelections::~MSettingsEventSelections()
{
  // default destructor
}


////////////////////////////////////////////////////////////////////////////////


bool MSettingsEventSelections::WriteXml(MXmlNode* Node)
{
   // Write content to an XML tree

  MXmlNode* aNode = 0;
  MXmlNode* bNode = 0;

  aNode = new MXmlNode(Node, "EventSelections");

  // Menu Eventselection:
  new MXmlNode(aNode, "EventSelectorTab", m_EventSelectorTab);
  new MXmlNode(aNode, "UseEventTypeCompton", m_EventTypeCompton);
  new MXmlNode(aNode, "UseEventTypeComptonNotTracked", m_EventTypeComptonNotTracked);
  new MXmlNode(aNode, "UseEventTypeComptonTracked", m_EventTypeComptonTracked);
  new MXmlNode(aNode, "UseEventTypePair", m_EventTypePair);
  new MXmlNode(aNode, "UseEventTypePhoto", m_EventTypePhoto); 
  new MXmlNode(aNode, "UseEventTypeUnidentifiable", m_EventTypeUnidentifiable); 
  new MXmlNode(aNode, "UseEventTypeDecay", m_EventTypeDecay);
  new MXmlNode(aNode, "FlaggedAsBad", m_FlaggedAsBad);

  new MXmlNode(aNode, "EventID", m_EventIdRangeMin, m_EventIdRangeMax);
  new MXmlNode(aNode, "TrackLength", m_TrackLengthRangeMin, m_TrackLengthRangeMax);
  new MXmlNode(aNode, "SequenceLength", m_SequenceLengthRangeMin, m_SequenceLengthRangeMax);
  new MXmlNode(aNode, "ClusteringQualityFactor", m_ClusteringQualityFactorRangeMin, m_ClusteringQualityFactorRangeMax);
  new MXmlNode(aNode, "ComptonQualityFactor", m_ComptonQualityFactorRangeMin, m_ComptonQualityFactorRangeMax);
  new MXmlNode(aNode, "TrackQualityFactor", m_TrackQualityFactorRangeMin, m_TrackQualityFactorRangeMax);
  new MXmlNode(aNode, "FirstEnergyWindow", m_FirstEnergyRangeMin, m_FirstEnergyRangeMax);
  new MXmlNode(aNode, "SecondEnergyWindow", m_SecondEnergyRangeMin, m_SecondEnergyRangeMax);
  new MXmlNode(aNode, "ThirdEnergyWindow", m_ThirdEnergyRangeMin, m_ThirdEnergyRangeMax);
  new MXmlNode(aNode, "FourthEnergyWindow", m_FourthEnergyRangeMin, m_FourthEnergyRangeMax);
  new MXmlNode(aNode, "EnergyGamma", m_EnergyRangeGammaMin, m_EnergyRangeGammaMax);
  new MXmlNode(aNode, "EnergyElectron", m_EnergyRangeElectronMin, m_EnergyRangeElectronMax);
  new MXmlNode(aNode, "ComptonAngle", m_ComptonAngleRangeMin, m_ComptonAngleRangeMax);

  new MXmlNode(aNode, "ThetaDeviationMax", m_ThetaDeviationMax);
  new MXmlNode(aNode, "FirstDistance", m_FirstDistanceRangeMin, m_FirstDistanceRangeMax);
  new MXmlNode(aNode, "AnyDistanceRange", m_DistanceRangeMin, m_DistanceRangeMax);
  new MXmlNode(aNode, "Time", m_TimeRangeMin, m_TimeRangeMax);
  new MXmlNode(aNode, "TimeWalk", m_TimeWalkRangeMin, m_TimeWalkRangeMax);
  new MXmlNode(aNode, "CoincidenceWindow", m_CoincidenceWindowRangeMin, m_CoincidenceWindowRangeMax);
  new MXmlNode(aNode, "OpeningAnglePair", m_OpeningAnglePairMin, m_OpeningAnglePairMax);
  new MXmlNode(aNode, "InitialEnergyDepositPair", m_InitialEnergyDepositPairMin, m_InitialEnergyDepositPairMax);
  
  bNode = new MXmlNode(aNode, "EarthHorizonCut");
  new MXmlNode(bNode, "Type", m_EHCType); 
  new MXmlNode(bNode, "Probability", m_EHCProbability);
  new MXmlNode(bNode, "ComptonProbabilityFile", m_EHCComptonProbabilityFileName); 
  new MXmlNode(bNode, "PairProbabilityFile", m_EHCPairProbabilityFileName);
  new MXmlNode(bNode, "Angle", m_EHCAngle);
  new MXmlNode(bNode, "EarthPosition", m_EHCEarthPosition); 

  bNode = new MXmlNode(aNode, "Source");  
  new MXmlNode(bNode, "UsePointSource", m_SourceUsePointSource);
  new MXmlNode(bNode, "Coordinates", m_SourceCoordinates);
  new MXmlNode(bNode, "Phi", m_SourcePhi); 
  new MXmlNode(bNode, "Theta", m_SourceTheta);
  new MXmlNode(bNode, "Longitude", m_SourceLongitude); 
  new MXmlNode(bNode, "Latitude", m_SourceLatitude);
  new MXmlNode(bNode, "X", m_SourceX); 
  new MXmlNode(bNode, "Y", m_SourceY); 
  new MXmlNode(bNode, "Z", m_SourceZ); 
  new MXmlNode(bNode, "ARM", m_SourceARMMin, m_SourceARMMax); 
  new MXmlNode(bNode, "SPD", m_SourceSPDMin, m_SourceSPDMax);

  bNode = new MXmlNode(aNode, "Beam");
  new MXmlNode(bNode, "UseBeam", m_BeamUse);
  new MXmlNode(bNode, "StartX", m_BeamStartX); 
  new MXmlNode(bNode, "StartY", m_BeamStartY); 
  new MXmlNode(bNode, "StartZ", m_BeamStartZ); 
  new MXmlNode(bNode, "FocalSpotX", m_BeamFocalSpotX); 
  new MXmlNode(bNode, "FocalSpotY", m_BeamFocalSpotY); 
  new MXmlNode(bNode, "FocalSpotZ", m_BeamFocalSpotZ); 
  new MXmlNode(bNode, "Radius", m_BeamRadius); 
  new MXmlNode(bNode, "BeamDepth", m_BeamDepth);

  bNode = new MXmlNode(aNode, "ExcludedDetectors");
  for (unsigned int i = 0; i < m_ExcludedDetectors.size(); ++i) {
    new MXmlNode(bNode, "ExcludedDetector",  m_ExcludedDetectors[i]);
  }

  return true;
}


////////////////////////////////////////////////////////////////////////////////


bool MSettingsEventSelections::ReadXml(MXmlNode* Node)
{
  // Retrieve the content from an XML tree
  
  MXmlNode* aNode = 0;
  MXmlNode* bNode = 0;
  MXmlNode* cNode = 0;
 
 
  if ((aNode = Node->GetNode("EventSelections")) != 0) {
    if ((bNode = aNode->GetNode("EventSelectorTab")) != 0) {
      m_EventSelectorTab = bNode->GetValueAsInt();
    }
    if ((bNode = aNode->GetNode("UseEventTypeCompton")) != 0) {
      m_EventTypeCompton = bNode->GetValueAsInt();
    }
    if ((bNode = aNode->GetNode("UseEventTypeComptonNotTracked")) != 0) {
      m_EventTypeComptonNotTracked = bNode->GetValueAsInt();
    }
    if ((bNode = aNode->GetNode("UseEventTypeComptonTracked")) != 0) {
      m_EventTypeComptonTracked = bNode->GetValueAsInt();
    }
    if ((bNode = aNode->GetNode("UseEventTypePair")) != 0) {
      m_EventTypePair = bNode->GetValueAsInt();
    }
    if ((bNode = aNode->GetNode("UseEventTypePhoto")) != 0) {
      m_EventTypePhoto = bNode->GetValueAsInt();
    }
    if ((bNode = aNode->GetNode("UseEventTypeUnidentifiable")) != 0) {
      m_EventTypeUnidentifiable = bNode->GetValueAsInt();
    }
    if ((bNode = aNode->GetNode("UseEventTypeDecay")) != 0) {
      m_EventTypeDecay = bNode->GetValueAsInt();
    }
    if ((bNode = aNode->GetNode("FlaggedAsBad")) != 0) {
      m_FlaggedAsBad = bNode->GetValueAsInt();
    }
    if ((bNode = aNode->GetNode("EventID")) != 0) {
      m_EventIdRangeMin = bNode->GetMinValueAsLong();
      m_EventIdRangeMax = bNode->GetMaxValueAsLong();
    }
    if ((bNode = aNode->GetNode("TrackLength")) != 0) {
      m_TrackLengthRangeMin = bNode->GetMinValueAsInt();
      m_TrackLengthRangeMax = bNode->GetMaxValueAsInt();
    }
    if ((bNode = aNode->GetNode("SequenceLength")) != 0) {
      m_SequenceLengthRangeMin = bNode->GetMinValueAsInt();
      m_SequenceLengthRangeMax = bNode->GetMaxValueAsInt();
    }
    if ((bNode = aNode->GetNode("ClusteringQualityFactor")) != 0) {
      m_ClusteringQualityFactorRangeMin = bNode->GetMinValueAsDouble();
      m_ClusteringQualityFactorRangeMax = bNode->GetMaxValueAsDouble();
    }
    if ((bNode = aNode->GetNode("ComptonQualityFactor")) != 0) {
      m_ComptonQualityFactorRangeMin = bNode->GetMinValueAsDouble();
      m_ComptonQualityFactorRangeMax = bNode->GetMaxValueAsDouble();
    }
    if ((bNode = aNode->GetNode("TrackQualityFactor")) != 0) {
      m_TrackQualityFactorRangeMin = bNode->GetMinValueAsDouble();
      m_TrackQualityFactorRangeMax = bNode->GetMaxValueAsDouble();
    }
    if ((bNode = aNode->GetNode("FirstEnergyWindow")) != 0) {
      m_FirstEnergyRangeMin = bNode->GetMinValueAsDouble();
      m_FirstEnergyRangeMax = bNode->GetMaxValueAsDouble();
    }
    if ((bNode = aNode->GetNode("SecondEnergyWindow")) != 0) {
      m_SecondEnergyRangeMin = bNode->GetMinValueAsDouble();
      m_SecondEnergyRangeMax = bNode->GetMaxValueAsDouble();
    }
    if ((bNode = aNode->GetNode("ThirdEnergyWindow")) != 0) {
      m_ThirdEnergyRangeMin = bNode->GetMinValueAsDouble();
      m_ThirdEnergyRangeMax = bNode->GetMaxValueAsDouble();
    }
    if ((bNode = aNode->GetNode("FourthEnergyWindow")) != 0) {
      m_FourthEnergyRangeMin = bNode->GetMinValueAsDouble();
      m_FourthEnergyRangeMax = bNode->GetMaxValueAsDouble();
    }
    if ((bNode = aNode->GetNode("EnergyGamma")) != 0) {
      m_EnergyRangeGammaMin = bNode->GetMinValueAsDouble();
      m_EnergyRangeGammaMax = bNode->GetMaxValueAsDouble();
    }
    if ((bNode = aNode->GetNode("EnergyElectron")) != 0) {
      m_EnergyRangeElectronMin = bNode->GetMinValueAsDouble();
      m_EnergyRangeElectronMax = bNode->GetMaxValueAsDouble();
    }
    if ((bNode = aNode->GetNode("ComptonAngle")) != 0) {
      m_ComptonAngleRangeMin = bNode->GetMinValueAsDouble();
      m_ComptonAngleRangeMax = bNode->GetMaxValueAsDouble();
    }
    if ((bNode = aNode->GetNode("ThetaDeviationMax")) != 0) {
      m_ThetaDeviationMax = bNode->GetValueAsDouble();
    }
    if ((bNode = aNode->GetNode("FirstDistance")) != 0) {
      m_FirstDistanceRangeMin = bNode->GetMinValueAsDouble();
      m_FirstDistanceRangeMax = bNode->GetMaxValueAsDouble();
    }
    if ((bNode = aNode->GetNode("AnyDistanceRange")) != 0) {
      m_DistanceRangeMin = bNode->GetMinValueAsDouble();
      m_DistanceRangeMax = bNode->GetMaxValueAsDouble();
    }
    if ((bNode = aNode->GetNode("Time")) != 0) {
      m_TimeRangeMin = bNode->GetMinValueAsDouble();
      m_TimeRangeMax = bNode->GetMaxValueAsDouble();
    }
    if ((bNode = aNode->GetNode("TimeWalk")) != 0) {
      m_TimeWalkRangeMin = bNode->GetMinValueAsDouble();
      m_TimeWalkRangeMax = bNode->GetMaxValueAsDouble();
    }
    if ((bNode = aNode->GetNode("CoincidenceWindow")) != 0) {
      m_CoincidenceWindowRangeMin = bNode->GetMinValueAsDouble();
      m_CoincidenceWindowRangeMax = bNode->GetMaxValueAsDouble();
    }
    if ((bNode = aNode->GetNode("OpeningAnglePair")) != 0) {
      m_OpeningAnglePairMin = bNode->GetMinValueAsDouble();
      m_OpeningAnglePairMax = bNode->GetMaxValueAsDouble();
    }
    if ((bNode = aNode->GetNode("InitialEnergyDepositPair")) != 0) {
      m_InitialEnergyDepositPairMin = bNode->GetMinValueAsDouble();
      m_InitialEnergyDepositPairMax = bNode->GetMaxValueAsDouble();
    }

    if ((bNode = aNode->GetNode("EarthHorizonCut")) != 0) {
      if ((cNode = bNode->GetNode("Type")) != 0) {
        m_EHCType = cNode->GetValueAsInt();
      }
      if ((cNode = bNode->GetNode("Probability")) != 0) {
        m_EHCProbability = cNode->GetValueAsDouble();
      }
      if ((cNode = bNode->GetNode("ComptonProbabilityFile")) != 0) {
        m_EHCComptonProbabilityFileName = cNode->GetValueAsString();
      }
      if ((cNode = bNode->GetNode("PairProbabilityFile")) != 0) {
        m_EHCPairProbabilityFileName = cNode->GetValueAsString();
      }
      if ((cNode = bNode->GetNode("Angle")) != 0) {
        m_EHCAngle = cNode->GetValueAsDouble();
      }
      if ((cNode = bNode->GetNode("EarthPosition")) != 0) {
        m_EHCEarthPosition = cNode->GetValueAsVector();
      }
    }

    if ((bNode = aNode->GetNode("Source")) != 0) {
      if ((cNode = bNode->GetNode("UsePointSource")) != 0) {
        m_SourceUsePointSource = cNode->GetValueAsBoolean();
      }
      if ((cNode = bNode->GetNode("Coordinates")) != 0) {
        m_SourceCoordinates = cNode->GetValueAsInt();
      }
      if ((cNode = bNode->GetNode("Phi")) != 0) {
        m_SourcePhi = cNode->GetValueAsDouble();
      }
      if ((cNode = bNode->GetNode("Theta")) != 0) {
        m_SourceTheta = cNode->GetValueAsDouble();
      }
      if ((cNode = bNode->GetNode("Longitude")) != 0) {
        m_SourceLongitude = cNode->GetValueAsDouble();
      }
      if ((cNode = bNode->GetNode("Latitude")) != 0) {
        m_SourceLatitude = cNode->GetValueAsDouble();
      }
      if ((cNode = bNode->GetNode("X")) != 0) {
        m_SourceX = cNode->GetValueAsDouble();
      }
      if ((cNode = bNode->GetNode("Y")) != 0) {
        m_SourceY = cNode->GetValueAsDouble();
      }
      if ((cNode = bNode->GetNode("Z")) != 0) {
        m_SourceZ = cNode->GetValueAsDouble();
      }
      if ((cNode = bNode->GetNode("ARM")) != 0) {
        m_SourceARMMin = cNode->GetMinValueAsDouble(); 
        m_SourceARMMax = cNode->GetMaxValueAsDouble();
      }
      if ((cNode = bNode->GetNode("SPD")) != 0) {
        m_SourceSPDMin = cNode->GetMinValueAsDouble(); 
        m_SourceSPDMax = cNode->GetMaxValueAsDouble();
      }
    }

    if ((bNode = aNode->GetNode("Beam")) != 0) {
      if ((cNode = bNode->GetNode("UseBeam")) != 0) {
        m_BeamUse = cNode->GetValueAsBoolean();
      }
      if ((cNode = bNode->GetNode("StartX")) != 0) {
        m_BeamStartX = cNode->GetValueAsDouble();
      }
      if ((cNode = bNode->GetNode("StartY")) != 0) {
        m_BeamStartY = cNode->GetValueAsDouble();
      }
      if ((cNode = bNode->GetNode("StartZ")) != 0) {
        m_BeamStartZ = cNode->GetValueAsDouble();
      }
      if ((cNode = bNode->GetNode("FocalSpotX")) != 0) {
        m_BeamFocalSpotX = cNode->GetValueAsDouble();
      }
      if ((cNode = bNode->GetNode("FocalSpotY")) != 0) {
        m_BeamFocalSpotY = cNode->GetValueAsDouble();
      }
      if ((cNode = bNode->GetNode("FocalSpotZ")) != 0) {
        m_BeamFocalSpotZ = cNode->GetValueAsDouble();
      }
      if ((cNode = bNode->GetNode("Radius")) != 0) {
        m_BeamRadius = cNode->GetValueAsDouble();
      }
      if ((cNode = bNode->GetNode("BeamDepth")) != 0) {
        m_BeamDepth = cNode->GetValueAsDouble();
      }
    }

    if ((bNode = aNode->GetNode("ExcludedDetectors")) != 0) {
      m_ExcludedDetectors.clear();
      for (unsigned int n = 0; n < bNode->GetNNodes(); ++n) {
        if ((cNode = bNode->GetNode(n)) != 0) {
          m_ExcludedDetectors.push_back(cNode->GetValue());
        }
      }
    }
  }


  return true;
}


// MSettingsEventSelections.cxx: the end...
////////////////////////////////////////////////////////////////////////////////
