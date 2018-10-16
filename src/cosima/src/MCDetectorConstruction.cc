/*
 * MCDetectorConstruction.cxx
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


// Standard:
#include <sstream>
using namespace std;

// Cosima:
#include "MCDetectorConstruction.hh"
#include "MC2DStripSD.hh"
#include "MCScintillatorSD.hh"
#include "MCCalorBarSD.hh"
#include "MCDriftChamberSD.hh"
#include "MCAngerCameraSD.hh"
#include "MCVoxel3DSD.hh"
#include "MCRunManager.hh"
#include "MCSource.hh"
#include "MCRun.hh"
#include "MCEventAction.hh"
#include "MCRegion.hh"
#include "MCGeometryConverter.hh"

// Geant4:
#include "G4SystemOfUnits.hh"
#include "globals.hh"
#include "G4Element.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
#include "G4Cons.hh"
#include "G4Polycone.hh"
#include "G4Polyhedra.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SDManager.hh"
#include "G4UserLimits.hh"
#include "G4Region.hh"
#include "G4VisAttributes.hh"
#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4NistManager.hh"

// MEGAlib:
#include "MAssert.h"
#include "MStreams.h"
#include "MDVolume.h"
#include "MDShapeBRIK.h"
#include "MDShapeTUBS.h"
#include "MDShapeSPHE.h"
#include "MDShapeTRD1.h"
#include "MDShapeTRD2.h"
#include "MDShapeTRAP.h"
#include "MDShapeCONE.h"
#include "MDShapeCONS.h"
#include "MDShapePGON.h"
#include "MDShapePCON.h"
#include "MDShapeIntersection.h"
#include "MDShapeUnion.h"
#include "MDShapeSubtraction.h"
#include "MDMaterial.h"
#include "MDMaterialComponent.h"
#include "MDDetector.h"
#include "MDStrip2D.h"
#include "MDCalorimeter.h"
#include "MDACS.h"
#include "MDDriftChamber.h"
#include "MDVoxel3D.h"
#include "MDGuardRing.h"
#include "MDAngerCamera.h"
#include "MVector.h"
#include "MRotation.h"

// Root:
#include "TRotMatrix.h"
#include "TMatrix.h" 
#include "TMath.h" 
//#include "TColor.h"


/******************************************************************************
 * Default constructor
 */
MCDetectorConstruction::MCDetectorConstruction(MCParameterFile& RunParameters) 
  : m_RunParameters(RunParameters), m_WorldVolume(0)
{
  // Intentionally left blank
}


/******************************************************************************
 * Default destructor
 */
MCDetectorConstruction::~MCDetectorConstruction()
{
  // Intentionally left blank

  delete m_Geometry;
}


/******************************************************************************
 * Default Geant4 Construct Method, returns the world volume 
 */
G4VPhysicalVolume* MCDetectorConstruction::Construct()
{
  if (m_WorldVolume == 0) {
    mout<<"Geometry is not correctly initialized"
        <<"this will end in a fatal system crash..."<<endl;
  }

  return m_WorldVolume;
}


/******************************************************************************
 * Should be directly called after constructor:
 * Initializes the geometry
 */
bool MCDetectorConstruction::Initialize()
{
  // Cleanup old geometry

  G4GeometryManager::GetInstance()->OpenGeometry();
  G4PhysicalVolumeStore::GetInstance()->Clean();
  G4LogicalVolumeStore::GetInstance()->Clean();
  G4SolidStore::GetInstance()->Clean();

  m_Geometry = new MDGeometryQuest(); 

  // Read and initilize the geometry from MEGAlib
  if (m_Geometry->ScanSetupFile(m_RunParameters.GetGeometryFileName().c_str(), false, false, false) == false) {
    return false;
  }


  // Initialize start area for all sources, but only if it has not been set (parameter file preceeds geometry file):
  vector<MCRun>& Runs = MCRunManager::GetMCRunManager()->GetRuns();
  for (unsigned int r = 0; r < Runs.size(); ++r) {
    vector<MCSource*>& Sources = Runs[r].GetSourceList();
    for (unsigned t = 0; t < Sources.size(); ++t) {
      if (Sources[t]->GetStartAreaType() == MCSource::c_StartAreaUnknown) {
        Sources[t]->SetStartAreaType(MCSource::c_StartAreaSphere);
        if (Sources[t]->SetStartAreaParameters(m_Geometry->GetStartSpherePosition()[0]*cm,
                                               m_Geometry->GetStartSpherePosition()[1]*cm,
                                               m_Geometry->GetStartSpherePosition()[2]*cm,
                                               0.0,
                                               0.0,
                                               0.0,
                                               m_Geometry->GetStartSphereRadius()*cm) == false) {
          return false;
        }
      }
    }
  }

  if (ConstructMaterials() == false) return false;
  if (ConstructVolumes() == false) return false;

  // Recursively position all volumes
  if (PositionVolumes(m_Geometry->GetWorldVolume()) == false) {
    cout<<"Error during positioning volumes. Aborting"<<endl;
    return false;
  }
    
  if (m_WorldVolume == 0) return false;

  if (ConstructDetectors() == false) return false;

  if (ConstructRegions() == false) return false;

  // Modify the trigger unit:
  // In order to be able to modify energies, thresholds etc. we only pre-trigger here
  if (m_Geometry->GetTriggerUnit() == 0) {
    mout<<"Error: Geometry has no trigger unit. Please define a trigger criteria in the geometry file."<<endl;
    return false;
  }
  // That means:
  // (a) No vetoes
  m_Geometry->GetTriggerUnit()->IgnoreVetoes(true);
  // (b) Thresholds are at zero
  m_Geometry->GetTriggerUnit()->IgnoreThresholds(true);
  
  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  return true;
}


/******************************************************************************
 * Create a rotation matrix 
 */
G4RotationMatrix* MCDetectorConstruction::CreateRotation(MDVolume* Volume)
{
  MVector xvcolumn, yvcolumn, zvcolumn;

  xvcolumn = MVector(1.,0.,0.);
  yvcolumn = MVector(0.,1.,0.);
  zvcolumn = MVector(0.,0.,1.);

  MRotation RotMatrix = Volume->GetRotationMatrix();

  xvcolumn = RotMatrix*xvcolumn;
  yvcolumn = RotMatrix*yvcolumn;
  zvcolumn = RotMatrix*zvcolumn;   

  G4ThreeVector xColumn(xvcolumn[0], xvcolumn[1], xvcolumn[2]);
  G4ThreeVector yColumn(yvcolumn[0], yvcolumn[1], yvcolumn[2]);
  G4ThreeVector zColumn(zvcolumn[0], zvcolumn[1], zvcolumn[2]);

  return new G4RotationMatrix(xColumn, yColumn, zColumn);
}


/******************************************************************************
 * Create a rotation matrix 
 */
G4RotationMatrix* MCDetectorConstruction::CreateRotation(MDOrientation* Orientation)
{
  MVector xvcolumn, yvcolumn, zvcolumn;

  xvcolumn = MVector(1.,0.,0.);
  yvcolumn = MVector(0.,1.,0.);
  zvcolumn = MVector(0.,0.,1.);

  MRotation RotMatrix = Orientation->GetRotationMatrix();

  xvcolumn = RotMatrix*xvcolumn;
  yvcolumn = RotMatrix*yvcolumn;
  zvcolumn = RotMatrix*zvcolumn;   

  G4ThreeVector xColumn(xvcolumn[0], xvcolumn[1], xvcolumn[2]);
  G4ThreeVector yColumn(yvcolumn[0], yvcolumn[1], yvcolumn[2]);
  G4ThreeVector zColumn(zvcolumn[0], zvcolumn[1], zvcolumn[2]);

  return new G4RotationMatrix(xColumn, yColumn, zColumn);
}


/******************************************************************************
 * Create detectors
 */
bool MCDetectorConstruction::ConstructDetectors()
{
  // Define detectors and sensitive material stuff...
  G4SDManager* SDManager = G4SDManager::GetSDMpointer();
  MCEventAction* EventAction = 
    (MCEventAction *) (G4EventManager::GetEventManager()->GetUserEventAction());

  int Type;
  MString Name;
  G4ThreeVector Vector;

  MDDetector* Detector = 0;
  MCSD* SD = 0;
  for (unsigned int d = 0; d < m_Geometry->GetNDetectors(); ++d) {

    Type = m_Geometry->GetDetectorAt(d)->GetType();
    Name = m_Geometry->GetDetectorAt(d)->GetName() + "SD";

    Detector = m_Geometry->GetDetectorAt(d);

    if (Type == MDDetector::c_Strip2D || Type == MDDetector::c_Strip3D) {
      MDStrip2D* Strip = dynamic_cast<MDStrip2D*>(Detector);
      
      MC2DStripSD* TwoDStripSD = new MC2DStripSD(Name.Data());
      SD = TwoDStripSD;
      
      if (Type == MDDetector::c_Strip3D) {
        TwoDStripSD->SetIs3D(true);
      } else {
        TwoDStripSD->SetIs3D(false);
      }

      // Add the dimensions:
      TwoDStripSD->SetSize(0.5*Strip->GetWidthX()*cm, 
                           0.5*Strip->GetWidthY()*cm);
      TwoDStripSD->SetOffset(Strip->GetOffsetX()*cm, 
                             Strip->GetOffsetY()*cm);
      TwoDStripSD->SetPitch(Strip->GetPitchX()*cm, 
                            Strip->GetPitchY()*cm);
      TwoDStripSD->SetNStrips(Strip->GetNStripsX(), 
                              Strip->GetNStripsY());
      if (Strip->HasGuardRing() == true) {
        TwoDStripSD->SetUniqueGuardringPosition(G4ThreeVector(Strip->GetGuardRing()->GetUniquePosition().X()*cm,
                                                              Strip->GetGuardRing()->GetUniquePosition().Y()*cm,
                                                              Strip->GetGuardRing()->GetUniquePosition().Z()*cm));
      }

      SDManager->AddNewDetector(TwoDStripSD);
      mdebug<<"Adding Strip detector for "<<Name<<endl;
      
      // Get the name of the sensitive volume and set its SD:
      for (unsigned int sv = 0; 
           sv < Detector->GetNSensitiveVolumes(); ++sv) {
        MString SenName = 
          Detector->GetSensitiveVolume(sv)->GetName() + 
          "Log";

        G4LogicalVolumeStore* SS = G4LogicalVolumeStore::GetInstance();
        for (unsigned int vs = 0; vs < SS->size(); ++vs) {
          if (SenName == SS->at(vs)->GetName().c_str()) {
            SS->at(vs)->SetSensitiveDetector(TwoDStripSD);
          }
        }
      }
    } 

    else if (Type == MDDetector::c_Voxel3D) {
      MDVoxel3D* Voxler = dynamic_cast<MDVoxel3D*>(Detector);
      
      MCVoxel3DSD* Voxel3DSD = new MCVoxel3DSD(Name.Data());
      SD = Voxel3DSD;

      // Add the dimensions:
      Voxel3DSD->SetSize(0.5*Voxler->GetWidthX()*cm, 
                         0.5*Voxler->GetWidthY()*cm,
                         0.5*Voxler->GetWidthZ()*cm);
      Voxel3DSD->SetOffset(Voxler->GetOffsetX()*cm, 
                           Voxler->GetOffsetY()*cm,
                           Voxler->GetOffsetZ()*cm);
      Voxel3DSD->SetVoxelSize(Voxler->GetVoxelSizeX()*cm, 
                              Voxler->GetVoxelSizeY()*cm,
                              Voxler->GetVoxelSizeZ()*cm);
      Voxel3DSD->SetNVoxels(Voxler->GetNVoxelsX(), 
                            Voxler->GetNVoxelsY(),
                            Voxler->GetNVoxelsZ());
      if (Voxler->HasGuardRing() == true) {
        Voxel3DSD->SetUniqueGuardringPosition(G4ThreeVector(Voxler->GetGuardRing()->GetUniquePosition().X()*cm,
                                                            Voxler->GetGuardRing()->GetUniquePosition().Y()*cm,
                                                            Voxler->GetGuardRing()->GetUniquePosition().Z()*cm));
      }

      SDManager->AddNewDetector(Voxel3DSD);
      mdebug<<"Adding Voxler detector for "<<Name<<endl;
      
      // Get the name of the sensitive volume and set its SD:
      for (unsigned int sv = 0; 
           sv < Detector->GetNSensitiveVolumes(); ++sv) {
        MString SenName = 
          Detector->GetSensitiveVolume(sv)->GetName() + 
          "Log";

        G4LogicalVolumeStore* SS = G4LogicalVolumeStore::GetInstance();
        for (unsigned int vs = 0; vs < SS->size(); ++vs) {
          if (SenName == SS->at(vs)->GetName().c_str()) {
            SS->at(vs)->SetSensitiveDetector(Voxel3DSD);
          }
        }
      }
    } 

    else if (Type == MDDetector::c_DriftChamber) {
      MDDriftChamber* Chamber = dynamic_cast<MDDriftChamber*>(Detector);
      
      MCDriftChamberSD* ChamberSD = new MCDriftChamberSD(Name.Data());
      SD = ChamberSD;

      // Add the dimensions:
      ChamberSD->SetSize(0.5*Chamber->GetWidthX()*cm, 
                         0.5*Chamber->GetWidthY()*cm);
      ChamberSD->SetOffset(Chamber->GetOffsetX()*cm, 
                           Chamber->GetOffsetY()*cm);
      ChamberSD->SetPitch(Chamber->GetPitchX()*cm, 
                          Chamber->GetPitchY()*cm);
      ChamberSD->SetNStrips(Chamber->GetNStripsX(), 
                            Chamber->GetNStripsY());

      ChamberSD->SetLightSpeed(Chamber->GetLightSpeed()*cm/s);
      ChamberSD->SetLightDetectorPosition(Chamber->GetLightDetectorPosition());
      ChamberSD->SetDriftConstant(Chamber->GetDriftConstant());
      ChamberSD->SetEnergyPerElectron(Chamber->GetEnergyPerElectron()*keV);


      SDManager->AddNewDetector(ChamberSD);
      mdebug<<"Adding drift chamber detector for "<<Name<<endl;
      
      // Get the name of the sensitive volume and set its SD:
      for (unsigned int sv = 0; 
           sv < Detector->GetNSensitiveVolumes(); ++sv) {
        MString SenName = 
          Detector->GetSensitiveVolume(sv)->GetName() + 
          "Log";

        G4LogicalVolumeStore* SS = G4LogicalVolumeStore::GetInstance();
        for (unsigned int vs = 0; vs < SS->size(); ++vs) {
          if (SenName == SS->at(vs)->GetName().c_str()) {
            SS->at(vs)->SetSensitiveDetector(ChamberSD);
          }
        }
      }
    } 

    else if (Type == MDDetector::c_Calorimeter) {
      MCCalorBarSD* CalorimeterSD = new MCCalorBarSD(Name.Data());
      SD = CalorimeterSD;
      SDManager->AddNewDetector(CalorimeterSD);
      mdebug<<"Adding calorimeter for "<<Name<<endl;

      // Check if we have a 3D calorimeter:
      MDCalorimeter* Calorimeter = dynamic_cast<MDCalorimeter*>(Detector);
      if (Calorimeter->HasDepthResolution() == true) {
        CalorimeterSD->SetIs3D(true);
      }
    
      // Get the name of the sensitive volume and set its SD:
      for (unsigned int sv = 0; 
           sv < Detector->GetNSensitiveVolumes(); ++sv) {
        MString SenName = 
          Detector->GetSensitiveVolume(sv)->GetName() 
          + "Log";

        G4LogicalVolumeStore* SS = G4LogicalVolumeStore::GetInstance();
        for (unsigned vs = 0; vs < SS->size(); ++vs) {
          if (SenName == SS->at(vs)->GetName().c_str()) {
            SS->at(vs)->SetSensitiveDetector(CalorimeterSD);
          }
        }
      }
    } 

    else if (Type == MDDetector::c_ACS) {
      MCScintillatorSD* S = new MCScintillatorSD(Name.Data());
      SD = S;
      SDManager->AddNewDetector(S);
      mout<<"Adding szintillator for "<<Name<<endl;
     
      S->SetCommonVolumeName(Detector->GetCommonVolume()->GetName().Data());

      MVector UniquePositionInCommon(0.0, 0.0, 0.0);

      // Get the name of the sensitive volumes and set its SD:
      for (unsigned int sv = 0; 
           sv < Detector->GetNSensitiveVolumes(); ++sv) {
        MString SenName = Detector->GetSensitiveVolume(sv)->GetName() + "Log";
        MVector Pos = Detector->GetSensitiveVolume(sv)->GetShape()->GetUniquePosition();
        // The position of the hits is the unique position in the first volume
        if (sv == 0) {
          S->SetUniquePosition(SenName.Data(), G4ThreeVector(Pos[0]*cm, Pos[1]*cm, Pos[2]*cm));
          // Get the position in the common volume:
          MDVolume* V = Detector->GetSensitiveVolume(sv);
          while (V != 0 && V != Detector->GetCommonVolume()) {
            UniquePositionInCommon = V->GetPositionInMotherVolume(UniquePositionInCommon);
            V = V->GetMother();
          }
        } else {
          // Build a volume tree:
          vector<MDVolume* > Vs;
          Vs.push_back(Detector->GetSensitiveVolume(sv));
          while (Vs.back() != 0 && Vs.back() != Detector->GetCommonVolume()) {
            Vs.push_back(Vs.back()->GetMother());
          }
          // Rotate and translate the unique position in common volume into this sensitive volume:
          MVector NewPos = UniquePositionInCommon;
          for (int i = int(Vs.size()) - 2; i >= 0; --i) {
            NewPos -= Vs[i]->GetPosition();           // translate 
            if (Vs[i]->IsRotated() == true) {
              NewPos = Vs[i]->GetRotationMatrix() * NewPos;   // rotate
            }
          }
          S->SetUniquePosition(SenName.Data(), G4ThreeVector(NewPos[0]*cm, NewPos[1]*cm, NewPos[2]*cm));
        }
        G4LogicalVolumeStore* SS = G4LogicalVolumeStore::GetInstance();
        for (unsigned vs = 0; vs < SS->size(); ++vs) {
          if (SenName == SS->at(vs)->GetName().c_str()) {
            //cout<<"SenName: "<<SenName<<endl;
            SS->at(vs)->SetSensitiveDetector(S);
          }
        }
      }
    }

    else if (Type == MDDetector::c_AngerCamera) {
      MCAngerCameraSD* AngerCameraSD = new MCAngerCameraSD(Name.Data());
      SD = AngerCameraSD;
      SDManager->AddNewDetector(AngerCameraSD);
      mout<<"Adding Anger camera for "<<Name<<endl;
     
      // Check if we have a 3D calorimeter:
      MDAngerCamera* AC = dynamic_cast<MDAngerCamera*>(Detector);
      if (AC->GetPositioning() == MDAngerCamera::c_PositionResolutionXYZ) {
        AngerCameraSD->SetIs3D(true);
      }
      
      // Get the name of the sensitive volumes and set its SD:
      // Get the name of the sensitive volume and set its SD:
      for (unsigned int sv = 0; sv < Detector->GetNSensitiveVolumes(); ++sv) {
        MString SenName = Detector->GetSensitiveVolume(sv)->GetName() + "Log";

        G4LogicalVolumeStore* SS = G4LogicalVolumeStore::GetInstance();
        for (unsigned vs = 0; vs < SS->size(); ++vs) {
          if (SenName == SS->at(vs)->GetName().c_str()) {
            SS->at(vs)->SetSensitiveDetector(AngerCameraSD);
          }
        }
      }
    }
    
    else if (Type == MDDetector::c_GuardRing) {
      // This is part of the strip & voxel detectors!
    } else {
      merr<<"DetectorType not yet implemented ("
          <<Detector->GetName()<<")!"<<endl;
      return false;
    }

    EventAction->SetCollectionName(Name.Data(), Type);

    // Add general information about the detector structure:
    Name = Detector->GetDetectorVolume()->GetName() + "Log";
    SD->SetDetectorVolumeName(Name.Data()); 
    Vector = G4ThreeVector(2*Detector->GetStructuralSize().X()*cm,
                           2*Detector->GetStructuralSize().Y()*cm,
                           2*Detector->GetStructuralSize().Z()*cm);
    SD->SetDetectorStructuralSize(Vector);
    Vector = G4ThreeVector(Detector->GetStructuralPitch().X()*cm,
                           Detector->GetStructuralPitch().Y()*cm,
                           Detector->GetStructuralPitch().Z()*cm);
    SD->SetDetectorStructuralPitch(Vector);
    Vector = G4ThreeVector(Detector->GetStructuralOffset().X()*cm,
                           Detector->GetStructuralOffset().Y()*cm,
                           Detector->GetStructuralOffset().Z()*cm);
    SD->SetDetectorStructuralOffset(Vector);

    MDShapeBRIK* Brik = (MDShapeBRIK*) (Detector->GetDetectorVolume()->GetShape());
    Vector = G4ThreeVector(Brik->GetSizeX()*cm,
                           Brik->GetSizeY()*cm,
                           Brik->GetSizeZ()*cm);
    SD->SetDetectorStructuralDimension(Vector);
    SD->SetDiscretizeHits(m_RunParameters.DiscretizeHits());

    if (Detector->GetEnergyLossType() == MDDetector::c_EnergyLossTypeMap) {
      SD->SetEnergyLossMap(Detector->GetEnergyLossMap());
    }

    SD->SetHasTimeResolution(Detector->HasTimeResolution());
  }

  return true;
}


/******************************************************************************
 * Create volumes
 */
bool MCDetectorConstruction::ConstructVolumes()
{
  MString Type; 
  MString VolumeName;
  MString MaterialName;
  MString LogName;
  G4VSolid* Solid = 0;
  //bool ROOTColorNotInitialize = false;
  
  // Step 1: Convert all shapes
  // We have to ensure to convert boolean volumes correctly since they need to 
  // there internal shapes to be alreday converted!
  
  list<MDShape*> AllShapes;
  for (unsigned int i = 0; i < m_Geometry->GetNShapes(); ++i) {
    AllShapes.push_back(m_Geometry->GetShapeAt(i));
    //cout<<"Input shapes: "<<AllShapes.back()->GetName()<<endl;
  }
  list<G4VSolid*> AllSolids;
  for (list<MDShape*>::iterator I = AllShapes.begin(); I != AllShapes.end(); ++I) {
    MDShape* Shape = (*I);
    // A. If not all sub shapes of this shape are defined, then push it back to the end to take care of it later
    bool AllSubShapesDefined = true;
    for (unsigned int sh = 0; sh < Shape->GetNSubShapes(); ++sh) {
      G4String Name = Shape->GetSubShape(sh)->GetName().GetString();
      bool Found = false;
      for (list<G4VSolid*>::iterator J = AllSolids.begin(); J != AllSolids.end(); ++J) {
        if ((*J)->GetName() == Name) {
          Found = true;
          break;
        }
      }
      if (Found == false) {
        AllSubShapesDefined = false;
        break;
      }
    }
    if (AllSubShapesDefined == false) {
      AllShapes.push_back(Shape);
      continue;
    }
    
    // B. Convert it
    Type = Shape->GetType();
    
    if (Type == "BRIK") {
      MDShapeBRIK* BRIK = dynamic_cast<MDShapeBRIK*>(Shape);
      Solid = new G4Box(Shape->GetName().GetString(), BRIK->GetSizeX()*cm, BRIK->GetSizeY()*cm, BRIK->GetSizeZ()*cm);
    } else if (Type == "TUBS") {
      MDShapeTUBS* TUBS = dynamic_cast<MDShapeTUBS*>(Shape);
      Solid = new G4Tubs(Shape->GetName().GetString(), 
                     TUBS->GetRmin()*cm, 
                     TUBS->GetRmax()*cm, 
                     TUBS->GetHeight()*cm,
                     TUBS->GetPhi1()*deg,
                     (TUBS->GetPhi2()-TUBS->GetPhi1())*deg);
    } else if (Type == "SPHE") {
      MDShapeSPHE* SPHE = dynamic_cast<MDShapeSPHE*>(Shape);
      Solid = new G4Sphere(Shape->GetName().GetString(), 
                       SPHE->GetRmin()*cm, 
                       SPHE->GetRmax()*cm, 
                       SPHE->GetPhimin()*deg,
                       (SPHE->GetPhimax()-SPHE->GetPhimin())*deg,
                       SPHE->GetThetamin()*deg,
                       (SPHE->GetThetamax()-SPHE->GetThetamin())*deg);
    } else if (Type == "TRD1") {
      MDShapeTRD1* TRD1 = dynamic_cast<MDShapeTRD1*>(Shape);
      Solid = new G4Trd(Shape->GetName().GetString(), 
                    TRD1->GetDx1()*cm, 
                    TRD1->GetDx2()*cm, 
                    TRD1->GetY()*cm,
                    TRD1->GetY()*cm,
                    TRD1->GetZ()*cm);
    } else if (Type == "TRD2") {
      MDShapeTRD2* TRD2 = dynamic_cast<MDShapeTRD2*>(Shape);
      Solid = new G4Trd(Shape->GetName().GetString(), 
                    TRD2->GetDx1()*cm, 
                    TRD2->GetDx2()*cm, 
                    TRD2->GetDy1()*cm,
                    TRD2->GetDy2()*cm,
                    TRD2->GetZ()*cm);
    }else if (Type == "CONE") {
      MDShapeCONE* CONE = dynamic_cast<MDShapeCONE*>(Shape);
      Solid = new G4Cons(Shape->GetName().GetString(), 
                     CONE->GetRminBottom()*cm, 
                     CONE->GetRmaxBottom()*cm, 
                     CONE->GetRminTop()*cm,
                     CONE->GetRmaxTop()*cm, 
                     CONE->GetHalfHeight()*cm, 
                     0.0*deg, 360.0*deg);
    } else if (Type == "CONS") {
      MDShapeCONS* CONS = dynamic_cast<MDShapeCONS*>(Shape);
      Solid = new G4Cons(Shape->GetName().GetString(), 
                     CONS->GetRminBottom()*cm, 
                     CONS->GetRmaxBottom()*cm, 
                     CONS->GetRminTop()*cm,
                     CONS->GetRmaxTop()*cm, 
                     CONS->GetHalfHeight()*cm, 
                     CONS->GetPhiMin()*deg, CONS->GetPhiMax()*deg);
    } else if (Type == "PCON") {
      MDShapePCON* PCON = dynamic_cast<MDShapePCON*>(Shape);
      double* z = new double[PCON->GetNSections()];
      double* rmin = new double[PCON->GetNSections()];
      double* rmax = new double[PCON->GetNSections()];
      for (unsigned int i = 0; i < PCON->GetNSections(); ++i) {
        z[i] = PCON->GetZ(i)*cm;
        rmin[i] = PCON->GetRmin(i)*cm;
        rmax[i] = PCON->GetRmax(i)*cm;
      }
      Solid = new G4Polycone(Shape->GetName().GetString(), 
                         PCON->GetPhi()*deg, 
                         PCON->GetDPhi()*deg, 
                         PCON->GetNSections(),
                         z,
                         rmin,
                         rmax);
    } else if (Type == "PGON") {
      MDShapePGON* PGON = dynamic_cast<MDShapePGON*>(Shape);
      double* z = new double[PGON->GetNSections()];
      double* rmin = new double[PGON->GetNSections()];
      double* rmax = new double[PGON->GetNSections()];
      for (unsigned int i = 0; i < PGON->GetNSections(); ++i) {
        z[i] = PGON->GetZ(i)*cm;
        rmin[i] = PGON->GetRmin(i)*cm;
        rmax[i] = PGON->GetRmax(i)*cm;
      }
      Solid = new G4Polyhedra(Shape->GetName().GetString(),
              PGON->GetPhi()*deg,
              PGON->GetDPhi()*deg,
              PGON->GetNSides(),
              PGON->GetNSections(),
              z,
              rmin,
              rmax);
    } else if (Type == "TRAP") {
      MDShapeTRAP* TRAP = dynamic_cast<MDShapeTRAP*>(Shape);
      // For Geant4 no value is allowed to be zero:
      double dz = TRAP->GetDz()*cm;
      if (dz <= 0) dz = 1.0E-10;  
      double theta = TRAP->GetTheta()*deg;
      double phi = TRAP->GetPhi()*deg;
      double h1 = TRAP->GetH1()*cm;
      double bl1 = TRAP->GetBl1()*cm;
      double tl1 = TRAP->GetTl1()*cm;
      double a1 = TRAP->GetAlpha1()*deg;
      double h2 = TRAP->GetH2()*cm;
      double bl2 = TRAP->GetBl2()*cm;
      double tl2 = TRAP->GetTl2()*cm;
      double a2 = TRAP->GetAlpha2()*deg;
      // Geant4 does not accept zeros...
      // So we need some kind of logic, that if we add very small values,
      // that everything still is planar...
      // This does not take into account each possible  case...

      int NZeros = 0;
      if (h1 <= 0) NZeros++;  
      if (bl1 <= 0) NZeros++;  
      if (tl1 <= 0) NZeros++;  
      if (h2 <= 0) NZeros++;  
      if (bl2 <= 0) NZeros++;  
      if (tl2 <= 0) NZeros++;  
        
      if (NZeros != 0) {
        mout<<endl;
        mout<<"SEVERE WARNING!"<<endl;
        mout<<"... for Shape TRAP:"<<endl;
        mout<<"One of the parameters is zero & Geant4 does not allow this!"<<endl;
        mout<<"Trying to estimate other parameters:"<<endl;
        mout<<"Start: "<<endl;
        mout<<"h1: "<<h1<<"  bl1: "<<bl1<<"  tl1: "<<tl1<<endl;
        mout<<"h2: "<<h2<<"  bl2: "<<bl2<<"  tl2: "<<tl2<<endl;

        const double Exp = 1E-10;
        if (NZeros == 1) {
          if (h1 <= 0) h1 = Exp*cm;  
          if (bl1 <= 0) bl1 = Exp*cm;  
          if (tl1 <= 0) tl1 = Exp*cm;  
          if (h2 <= 0) h2 = Exp*cm;  
          if (bl2 <= 0) bl2 = Exp*cm;  
          if (tl2 <= 0) tl2 = Exp*cm;  
        } else {
          double Ratio = 0.0;
          if (h1 != 0 && h2 != 0) Ratio = h1/h2;
          else if (tl1 != 0 && tl2 != 0) Ratio = tl1/tl2;
          else if (bl1 != 0 && bl2 != 0) Ratio = tl1/tl2;
          
          if (Ratio == 0.0) {
            if (h1 != 0) Ratio = 1/Exp;
            else Ratio = Exp;
          }

          if (h1 != 0 && h2 == 0) h2=h1/Ratio;
          if (h1 == 0 && h2 != 0) h1=h2*Ratio;
            
          if (tl1 == 0 && tl2 == 0) {
            if (h1 > h2) tl1 = Exp*cm;
            else tl2 = Exp*cm;
          }
          if (tl1 != 0 && tl2 == 0) tl2=tl1/Ratio;
          if (tl1 == 0 && tl2 != 0) tl1=tl2*Ratio;
            
          if (bl1 == 0 && bl2 == 0) {
            if (h1 > h2) bl1 = Exp*cm;
            else bl2 = Exp*cm;
          }
          if (bl1 != 0 && bl2 == 0) bl2=bl1/Ratio;
          if (bl1 == 0 && bl2 != 0) bl1=bl2*Ratio;
        }
          
        // If still something is missing
        if (h1 <= 0) h1 = 1.0E-6*cm;  
        if (bl1 <= 0) bl1 = 1.0E-6*cm;  
        if (tl1 <= 0) tl1 = 1.0E-6*cm;  
        if (h2 <= 0) h2 = 1.0E-6*cm;  
        if (bl2 <= 0) bl2 = 1.0E-6*cm;  
        if (tl2 <= 0) tl2 = 1.0E-6*cm;  
         
        mout<<"Final:"<<endl;
        mout<<"h1: "<<h1<<"  bl1: "<<bl1<<"  tl1: "<<tl1<<endl;
        mout<<"h2: "<<h2<<"  bl2: "<<bl2<<"  tl2: "<<tl2<<endl;
        mout<<"This estimation is a very risky thing..."<<endl;
        mout<<"Better correct your geometry by yourself!!!!"<<endl;
        mout<<"END SEVERE WARNING!"<<endl;
        mout<<endl;
      } 

      Solid = new G4Trap(Shape->GetName().GetString(), 
                     dz, 
                     theta, 
                     phi,
                     h1,
                     bl1,
                     tl1,
                     a1,
                     h2,
                     bl2,
                     tl2,
                     a2);
    } else if (Type == "Subtraction") {
      MDShapeSubtraction* Subtraction = dynamic_cast<MDShapeSubtraction*>(Shape);
      
      // a) Find the solids
      G4String MinuendName = Subtraction->GetMinuend()->GetName().GetString();
      G4String SubtrahendName = Subtraction->GetSubtrahend()->GetName().GetString();

      G4VSolid* MinuendSolid = 0;
      G4VSolid* SubtrahendSolid = 0;
      for (list<G4VSolid*>::iterator J = AllSolids.begin(); J != AllSolids.end(); ++J) {
        if ((*J)->GetName() == MinuendName) {
          MinuendSolid = *J;
        }
        if ((*J)->GetName() == SubtrahendName) {
          SubtrahendSolid = *J;
        }
      }
      if (MinuendSolid == 0) {
        merr<<"Fatal error: No solid found with name "<<MinuendName<<endl;
        return false;
      }
      if (SubtrahendSolid == 0) {
        merr<<"Fatal error: No solid found with name "<<SubtrahendName<<endl;
        return false;
      }
      
      // b) Create the solid
      MDOrientation* O = Subtraction->GetOrientationSubtrahend() ;
      G4ThreeVector Pos(O->GetPosition().X()*cm, O->GetPosition().Y()*cm, O->GetPosition().Z()*cm);
      Solid = new G4SubtractionSolid(Shape->GetName().GetString(), MinuendSolid, SubtrahendSolid, CreateRotation(O), Pos);
    } else if (Type == "Intersection") {
      MDShapeIntersection* Intersection = dynamic_cast<MDShapeIntersection*>(Shape);
      
      // a) Find the solids
      G4String MinuendName = Intersection->GetShapeA()->GetName().GetString();
      G4String SubtrahendName = Intersection->GetShapeB()->GetName().GetString();

      G4VSolid* MinuendSolid = 0;
      G4VSolid* SubtrahendSolid = 0;
      for (list<G4VSolid*>::iterator J = AllSolids.begin(); J != AllSolids.end(); ++J) {
        if ((*J)->GetName() == MinuendName) {
          MinuendSolid = *J;
        }
        if ((*J)->GetName() == SubtrahendName) {
          SubtrahendSolid = *J;
        }
      }
      if (MinuendSolid == 0) {
        merr<<"Fatal error: No solid found with name "<<MinuendName<<endl;
        return false;
      }
      if (SubtrahendSolid == 0) {
        merr<<"Fatal error: No solid found with name "<<SubtrahendName<<endl;
        return false;
      }
      
      // b) Create the solid
      MDOrientation* O = Intersection->GetOrientationShapeB() ;
      G4ThreeVector Pos(O->GetPosition().X()*cm, O->GetPosition().Y()*cm, O->GetPosition().Z()*cm);
      Solid = new G4IntersectionSolid(Shape->GetName().GetString(), MinuendSolid, SubtrahendSolid, CreateRotation(O), Pos);
    } else if (Type == "Union") {
      MDShapeUnion* Union = dynamic_cast<MDShapeUnion*>(Shape);
      
      // a) Find the solids
      G4String MinuendName = Union->GetAugend()->GetName().GetString();
      G4String SubtrahendName = Union->GetAddend()->GetName().GetString();

      G4VSolid* MinuendSolid = 0;
      G4VSolid* SubtrahendSolid = 0;
      for (list<G4VSolid*>::iterator J = AllSolids.begin(); J != AllSolids.end(); ++J) {
        if ((*J)->GetName() == MinuendName) {
          MinuendSolid = *J;
        }
        if ((*J)->GetName() == SubtrahendName) {
          SubtrahendSolid = *J;
        }
      }
      if (MinuendSolid == 0) {
        merr<<"Fatal error: No solid found with name "<<MinuendName<<endl;
        return false;
      }
      if (SubtrahendSolid == 0) {
        merr<<"Fatal error: No solid found with name "<<SubtrahendName<<endl;
        return false;
      }
      
      // b) Create the solid
      MDOrientation* O = Union->GetOrientationAddend() ;
      G4ThreeVector Pos(O->GetPosition().X()*cm, O->GetPosition().Y()*cm, O->GetPosition().Z()*cm);
      Solid = new G4UnionSolid(Shape->GetName().GetString(), MinuendSolid, SubtrahendSolid, CreateRotation(O), Pos);
    } else {
      merr<<"Unknown volume type: "<<Type<<endl;
      return false;
    }
    AllSolids.push_back(Solid);
  }
  
  //for (list<G4VSolid*>::iterator J = AllSolids.begin(); J != AllSolids.end(); ++J) {
  //  cout<<"Solid: "<<(*J)->GetName()<<endl;
  //}
  
  
  // Step 2: All shapes are now defined, so we can build the volumes
  
  
  for (unsigned int v = 0; v < m_Geometry->GetNVolumes(); ++v) {
    if (m_Geometry->GetVolumeAt(v)->GetCloneTemplate() == 0) {
      
      G4String ShapeName = m_Geometry->GetVolumeAt(v)->GetShape()->GetName().GetString();
      VolumeName = m_Geometry->GetVolumeAt(v)->GetName();
      MaterialName = m_Geometry->GetVolumeAt(v)->GetMaterial()->GetName();
      LogName = VolumeName + "Log";

      Solid = 0;
      for (list<G4VSolid*>::iterator J = AllSolids.begin(); J != AllSolids.end(); ++J) {
        if ((*J)->GetName() == ShapeName) {
          Solid = *J;
          break;
        }
      }
      if (Solid == 0) {
        merr<<"Fatal error: No solid found with name "<<ShapeName<<" for volume "<<VolumeName<<endl;
        return false;
      }
      
      G4LogicalVolume* LV = 
        new G4LogicalVolume(Solid, 
                            G4Material::GetMaterial(MaterialName.Data()), 
                            LogName.Data(), 0, 0, 0);
      G4VisAttributes* Vis = new G4VisAttributes();
      if (m_Geometry->GetVolumeAt(v)->GetVisibility() == 0) {
        Vis->SetVisibility(false);
      } else {
        Vis->SetVisibility(true);       
      }
      /*
      // Convert ROOT color to Geant4 color
      TColor* C = gROOT->GetColor(m_Geometry->GetVolumeAt(v)->GetColor());
      if (C != 0) {
        Vis->SetColor(G4Color(C->GetRed(), C->GetGreen(), C->GetBlue(), C->GetAlpha()));
      } else {
        if (ROOTColorNotInitialize == false) {
          mout<<"ROOT colors not initialized. Using black and white..."<<endl;
          ROOTColorNotInitialize = true;
        }
      }
      */
      LV->SetVisAttributes(Vis);
   
      mdebug<<"Generating Log: "<<LogName.Data()<<" mat="<<MaterialName.Data()<<endl;
    } // no clone 
  }

  return true;
}


/******************************************************************************
 * Create materials
 */
bool MCDetectorConstruction::ConstructMaterials()
{
  mdebug<<"Defining materials..."<<endl;

  if (m_Geometry->IsScanned() == false) return false;
 
  G4NistElementBuilder ElementBuilder(0);
  
  // List of used materials:
  vector<MDMaterial*> UnusedMaterials = m_Geometry->GetListOfUnusedMaterials();
  
  for (unsigned int mat = 0; mat < m_Geometry->GetNMaterials(); ++mat) {
    // Ignore unused materials:
    if (m_RunParameters.CreateCrossSectionFiles() == false) {
      if (find(UnusedMaterials.begin(), UnusedMaterials.end(), m_Geometry->GetMaterialAt(mat)) != UnusedMaterials.end()) {
        mdebug<<"Removing unused material: "<<m_Geometry->GetMaterialAt(mat)->GetName()<<endl;
        continue;
      }
    }
    
    G4Material* Material = 
      new G4Material(m_Geometry->GetMaterialAt(mat)->GetName().Data(), 
                     m_Geometry->GetMaterialAt(mat)->GetDensity()*g/cm3,
                     m_Geometry->GetMaterialAt(mat)->GetNComponents());
    
    
    for (unsigned int c = 0; 
         c < m_Geometry->GetMaterialAt(mat)->GetNComponents(); ++c) {
      MDMaterialComponent* Component = 
        m_Geometry->GetMaterialAt(mat)->GetComponentAt(c);

      ostringstream LongName;
      LongName<<"LongName_"<<m_Geometry->GetMaterialAt(mat)->GetName().Data()
              <<"_El"<<c+1<<endl;

      ostringstream ShortName;
      ShortName<<m_Geometry->GetMaterialAt(mat)->GetName().Data()
               <<"_El"<<c+1<<endl;

      double A = Component->GetAtomicWeight();
      double Z = Component->GetAtomicNumber();

      if (A < 1 && Z < 1) {
        mout<<m_Geometry->GetMaterialAt(mat)->GetName().Data()<<": Probably found Geant3 vaccum: upgrading to Geant4 vacuum"<<endl; 
        A = 1;
        Z = 1;
      }

      G4Element* Element = 0;
      if (Component->HasNaturalIsotopeComposition() == true) {
        Element = ElementBuilder.FindOrBuildElement(Z, true);
      } else {
        Element =  new G4Element(LongName.str(), ShortName.str(), Z, A*g/mole);
      }
      
      if (Element == 0) {
        merr<<"Couldn't find all elements of: "<<m_Geometry->GetMaterialAt(mat)->GetName()<<" Missing: Z="<<Z<<endl;
        return false;
      }
      
      if (Component->GetWeightingType() == MDMaterialComponentWeightingType::c_ByAtoms) {
        Material->AddElement(Element, Component->GetWeightingByAtoms());
      } else {
        Material->AddElement(Element, Component->GetWeightingByMass());
      }
    }
    mdebug<<Material<<endl;
  }

  return true;
}


/******************************************************************************
 * Create regions
 */
bool MCDetectorConstruction::ConstructRegions()
{
  mdebug<<"Defining regions..."<<endl;

  if (m_Geometry->IsScanned() == false) return false;
 
  const vector<MCRegion>& Regions = m_RunParameters.GetRegionList();

  G4LogicalVolumeStore* LVS = G4LogicalVolumeStore::GetInstance();
  
  // Loop over all regions from the input file and define them 
  for (unsigned int r = 0; r < Regions.size(); ++r) {
    // Check if volume exists in the geomega volume tree
    if (m_Geometry->GetWorldVolume()->ContainsVolume(Regions[r].GetVolumeName(), true) == true) {
      MString Name = Regions[r].GetVolumeName() + "Log";
      bool Found = false;
      // Find the volume in the logical volume store
      for (unsigned lv = 0; lv < LVS->size(); ++lv) {
        if (Name == LVS->at(lv)->GetName().c_str()) {
          // Define region. The cuts are set in the physics list
          G4Region* Region = new G4Region(Regions[r].GetName().Data());
          Region->AddRootLogicalVolume(LVS->at(lv));
          Found = true;
          break;
        }
      }
      if (Found == false) {
        mout<<"Volume "<<Regions[r].GetVolumeName()<<" defining region ";
        mout<<Regions[r].GetName()<<" is not found in the LogicalVolumeStore."<<endl;
        mout<<"Available volumes are: ";
        for (unsigned lv = 0; lv < LVS->size(); ++lv) mout<<LVS->at(lv)->GetName().c_str()<<" ";
        mout<<endl;
        return false;
      }
    } else {
      mout<<"Volume "<<Regions[r].GetVolumeName()<<" defining region ";
      mout<<Regions[r].GetName()<<" is not part of the world volume."<<endl;
      mout<<"Make sure that the spelling is correct and that the volume is not a virtual volume."<<endl;
      mout<<"If the geometry contains virtual volumes, the names of the volume might have changed."<<endl;
      mout<<"In this case use geomega to determine its new name."<<endl;
      return false;
    }
  }

  return true;
}


/******************************************************************************
 * Recursively position all volumes - in the it smother volumes
 */
bool MCDetectorConstruction::PositionVolumes(MDVolume* Volume)
{
  // This is far too complex to REALLY understand without investing time ...
  // But since it works: "Don't f&*k with it!" 
  
  MDVolume* VC;
  MString Name, MotherName, CopyName;

  if ((VC = Volume->GetCloneTemplate()) != 0) {
      
    // If the volume is a copy, then reposition the CopyOf volume:
    
    Name = Volume->GetName();
    CopyName = Volume->GetCloneTemplate()->GetName() + "Log";
    MotherName = Volume->GetMother()->GetName() + "Log";
 

    G4LogicalVolume* ThisLog = 0;
    G4LogicalVolume* MotherLog = 0;

    G4LogicalVolumeStore* SS = G4LogicalVolumeStore::GetInstance();
    for (unsigned lv = 0; lv < SS->size(); ++lv) {
      if (ThisLog == 0 && CopyName == SS->at(lv)->GetName().c_str()) {
        ThisLog = SS->at(lv);
      }
      if (MotherLog == 0 && MotherName == SS->at(lv)->GetName().c_str()) {
        MotherLog = SS->at(lv);
      }
      if (ThisLog != 0 && MotherLog != 0) break;
    }

    if (ThisLog == 0) {
      mout<<"Unable to find copy-log: "<<CopyName<<endl;
      return false;
    }

    if (MotherLog == 0) {
      mout<<"Unable to find mother-log: "<<MotherName<<endl;
      return false;
    }
    G4PVPlacement* P = 
      new G4PVPlacement(CreateRotation(Volume), 
                        G4ThreeVector(Volume->GetPosition().X()*cm, 
                                      Volume->GetPosition().Y()*cm, 
                                      Volume->GetPosition().Z()*cm), 
                        ThisLog, 
                        Volume->GetName().Data(), 
                        MotherLog, 
                        false, 0);
    if (m_RunParameters.PerformOverlapCheck() == true) {
      P->CheckOverlaps(m_RunParameters.ResolutionOverlapCheck(), m_RunParameters.ToleranceOverlapCheck(), false);
    }

    mdebug<<"Positioning clone: "<<CopyName<<" in "<<MotherName
          <<" at ("<<Volume->GetPosition().X()<<","
          <<Volume->GetPosition().Y()<<","<<Volume->GetPosition().Z()
          <<") with Name "<<Volume->GetName().Data()<<endl;
    
    if (Volume->GetCloneTemplate()->AreCloneTemplateDaughtersWritten() == 
        kFALSE) {
      for (unsigned int i = 0; i < Volume->GetNDaughters(); i++) {
        if (PositionVolumes(Volume->GetDaughterAt(i)) == false) return false;
      }
      Volume->GetCloneTemplate()->SetCloneTemplateDaughtersWritten(kTRUE);
    }
  } else {
    // If it is no copy then position it *only*
    // when it has a mother or is the root volume
    if (Volume->GetMother() !=  0 || Volume->IsWorldVolume() == true) {
      
      Name = Volume->GetName() + "Log";
      if (Volume->IsWorldVolume() == false) {
        MotherName = Volume->GetMother()->GetName() + "Log";
      }

      G4LogicalVolume* ThisLog = 0;
      G4LogicalVolume* MotherLog = 0;

      G4LogicalVolumeStore* SS = G4LogicalVolumeStore::GetInstance();
      for (unsigned lv = 0; lv < SS->size(); ++lv) {
        if (ThisLog == 0 && Name == SS->at(lv)->GetName().c_str()) {
          ThisLog = SS->at(lv);
        }
        if (MotherLog == 0 && MotherName == SS->at(lv)->GetName().c_str()) {
          MotherLog = SS->at(lv);
        }
        if (ThisLog != 0 && MotherLog != 0) break;
      }

      if (ThisLog == 0) {
        merr<<"Unable to find this-log: "<<Name<<endl;
        return false;
      }
      
      if (MotherLog == 0 && Volume->IsWorldVolume() == false) {
        merr<<"Unable to find mother-log: "<<MotherName<<endl;
        return false;
      }


      if (Volume->IsWorldVolume() == true) {
        m_WorldVolume = 
          new G4PVPlacement(0, 
                            G4ThreeVector(Volume->GetPosition().X()*cm, 
                                          Volume->GetPosition().Y()*cm, 
                                          Volume->GetPosition().Z()*cm), 
                            Volume->GetName().Data(), ThisLog, 0, false, 0);
        if (m_RunParameters.PerformOverlapCheck() == true) {
          m_WorldVolume->CheckOverlaps(m_RunParameters.ResolutionOverlapCheck(), m_RunParameters.ToleranceOverlapCheck(), false);
        }
      } else {
        G4PVPlacement* P = 
          new G4PVPlacement(CreateRotation(Volume),
                            G4ThreeVector(Volume->GetPosition().X()*cm, 
                                          Volume->GetPosition().Y()*cm, 
                                          Volume->GetPosition().Z()*cm), 
                            ThisLog, 
                            Volume->GetName().Data(), 
                            MotherLog, 
                            false, 0);
        if (m_RunParameters.PerformOverlapCheck() == true) {
          P->CheckOverlaps(m_RunParameters.ResolutionOverlapCheck(), m_RunParameters.ToleranceOverlapCheck(), false);
        }
      }

      mdebug<<"Positioning original: "<<Name<<" in "<<MotherName
            <<" at ("<<Volume->GetPosition().X()<<","<<Volume->GetPosition().Y()
            <<","<<Volume->GetPosition().Z()<<") with Name "
            <<Volume->GetName().Data()<<endl;

    } else {
      //cout<<"   is ignored..."<<endl;
    }

    // Do the same for all daughters:
    for (unsigned int i = 0; i < Volume->GetNDaughters(); i++) {
      if (PositionVolumes(Volume->GetDaughterAt(i)) == false) return false;
    }
  }
  
  return true;
}


/******************************************************************************
 * Return a random position in the given volume
 */
G4ThreeVector MCDetectorConstruction::GetRandomPosition(MString VolumeName)
{
  MVector V = m_Geometry->GetRandomPositionInVolume(VolumeName);

  return G4ThreeVector(V[0]*cm, V[1]*cm, V[2]*cm);
}


/******************************************************************************
 * Return true if the volume exists:
 */
bool MCDetectorConstruction::HasVolume(const MString& VolumeName) const
{
  return (m_Geometry->GetVolume(VolumeName) == nullptr) ? false: true;
}


/******************************************************************************
 * Return true if the volume is valid volume in the geometry
 */
bool MCDetectorConstruction::IsValidVolume(MString VolumeName)
{
  // No clue why I need a cast here....
  if (dynamic_cast<MDGeometry*>(m_Geometry)->GetVolume(VolumeName) == 0) {
    return false;
  }

  return true;
}


/******************************************************************************
 * Return the (Geomega) hash of the given material 
 */
unsigned long MCDetectorConstruction::GetMaterialHash(const G4Material* Material)
{
  // No clue why I need a cast here....
  MDMaterial* M = dynamic_cast<MDGeometry*>(m_Geometry)->GetMaterial(Material->GetName());
  if (M == 0) {
    merr<<"Material "<<Material->GetName()<<" not found in Geomega geometry!"<<show;
    return 0;
  }

  return M->GetHash();
}


/*
 * MCDetectorConstruction.cc: the end...
 ******************************************************************************/
