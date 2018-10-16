/*
 * MImageSpheric.h
 *
 * Copyright (C) by Andreas Zoglauer.
 * All rights reserved.
 *
 * Please see the source-file for the copyright-notice.
 *
 */


#ifndef __MImageSpheric__
#define __MImageSpheric__


////////////////////////////////////////////////////////////////////////////////


// ROOT libs:
#include <TGaxis.h>

// MEGAlib libs:
#include "MGlobal.h"
#include "MImage2D.h"

// Forward declarations:


////////////////////////////////////////////////////////////////////////////////


class MImageSpheric : public MImage2D
{
  // public interface:
 public:
  //! Default constructor
  MImageSpheric();
  //! Standard constructor
  MImageSpheric(MString Title, double *IA,
                MString xTitle, double xMin, double xMax, int xNBins, 
                MString yTitle, double yMin, double yMax, int yNBins, 
                MString vTitle, int Spectrum = c_Rainbow, int DrawOption = c_COLCONT4Z);
  //! Standard destructor
  virtual ~MImageSpheric();

  //! Clone this image
  virtual MImage* Clone();

  //! Set the image array and redisplay it
  virtual void SetImageArray(double* Array);

  //! Display the histogram in the given canvas
  virtual void Display(TCanvas* Canvas = nullptr);
 
  // protected methods:
 protected:


  // private methods:
 private:



  // protected members:
 protected:

  // private members:
 private:
  //! The new y axis
  TGaxis* m_YAxis;


#ifdef ___CLING___
 public:
  ClassDef(MImageSpheric, 0) // class containing an astropysical image
#endif

};

#endif


////////////////////////////////////////////////////////////////////////////////
