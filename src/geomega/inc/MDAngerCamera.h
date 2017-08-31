/*
 * MDAngerCamera.h
 *
 * Copyright (C) by Andreas Zoglauer.
 * All rights reserved.
 *
 * Please see the source-file for the copyright-notice.
 *
 */


#ifndef __MDAngerCamera__
#define __MDAngerCamera__


////////////////////////////////////////////////////////////////////////////////


// Standard libs:
#include <map>
using namespace std;

// ROOT libs:

// MEGAlib libs:
#include "MGlobal.h"
#include "MDDetector.h"
#include "MDVolume.h"
#include "MVector.h"

// Forward declarations:


////////////////////////////////////////////////////////////////////////////////


class MDAngerCamera : public MDDetector
{
  // public interface:
 public:
  MDAngerCamera(MString String = "");
  MDAngerCamera(const MDAngerCamera& A);
  virtual ~MDAngerCamera();

  virtual MDDetector* Clone();
  //! Copy data to named detectors
  virtual bool CopyDataToNamedDetectors();

  virtual bool Validate();

  virtual void Noise(MVector& Pos, double& Energy, double& Time, MDVolume* Volume) const;
  virtual vector<MDGridPoint> Discretize(const MVector& Pos, 
                                         const double& Energy, 
                                         const double& Time,
                                         MDVolume* Volume) const;
  //! Return the Grid point of this position
  virtual MDGridPoint GetGridPoint(const MVector& Pos) const;
  //! Return a position in detector volume coordinates
  virtual MVector GetPositionInDetectorVolume(const unsigned int xGrid, 
                                              const unsigned int yGrid,
                                              const unsigned int ZGrid,
                                              const MVector PositionInGrid,
                                              const unsigned int Type,
                                              MDVolume* Volume);
  virtual MVector GetPositionResolution(const MVector& Pos, const double Energy) const;

  virtual MString GetGeomega() const;
  virtual MString ToString() const;

  //! Set the type XY-only or XYZ-only
  void SetPositioning(int Positioning) { m_PositionResolutionType = Positioning; }
  void SetPositionResolution(const double Energy, const double Resolution);

  static const int c_PositionResolutionUnknown;
  static const int c_PositionResolutionXY;
  static const int c_PositionResolutionXYZ;


  // protected methods:
 protected:


  // private methods:
 private:



  // protected members:
 protected:


  // private members:
 private:
  //! Type of the positioning (XY only or XYZ)
  int m_PositionResolutionType;
  //! Position solution
  MFunction m_PositionResolution; 


#ifdef ___CINT___
 public:
  ClassDef(MDAngerCamera, 0) // a basic AngerCamera detector
#endif

};

#endif


////////////////////////////////////////////////////////////////////////////////
