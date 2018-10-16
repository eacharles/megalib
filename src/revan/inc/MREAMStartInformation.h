/*
 * MREAMStartInformation.h
 *
 * Copyright (C) by Andreas Zoglauer.
 * All rights reserved.
 *
 * Please see the source-file for the copyright-notice.
 *
 */


#ifndef __MREAMStartInformation__
#define __MREAMStartInformation__


////////////////////////////////////////////////////////////////////////////////


// Standard libs:

// ROOT libs:

// MEGAlib libs:
#include "MGlobal.h"
#include "MREAM.h"
#include "MRESE.h"

// Forward declarations:


////////////////////////////////////////////////////////////////////////////////


class MREAMStartInformation : public MREAM
{
  // public interface:
 public:
  //! Default constructor for a charged particle direction as measured by a Strip3DStartInformation detector
  MREAMStartInformation();
  //! Copy constructor for a charged particle direction
  MREAMStartInformation(const MREAMStartInformation& REAMStartInformation);
  //! Default destructor for a charged particle direction
  virtual ~MREAMStartInformation();

  //! Default assignment constructor
  const MREAMStartInformation& operator=(const MREAMStartInformation&);

  //! In almost all cases when we need a replica of this object we do not know what it is -
  //! Thus we need a general duplication mechanism
  virtual MREAM* Clone();

  //! Dump the essential content of this REAM
  virtual MString ToString(int Level = 0);


  //! Set original position
  void SetPosition(const MVector Position) { m_Position = Position; }
  //! Set original direction
  void SetDirection(const MVector Direction) { m_Direction = Direction; }
  //! Set original polarization
  void SetPolarization(const MVector Polarization) { m_Polarization = Polarization; }
  //! Set original energy
  void SetEnergy(const double Energy) { m_Energy = Energy; }
  //! Return original position
  MVector GetPosition() const { return m_Position; }
  //! Return original position
  MVector GetPolarization() const { return m_Polarization; }
  //! Return original direction 
  MVector GetDirection() const { return m_Direction; }
  //! Return original energy 
  double GetEnergy() const { return m_Energy; }

  // protected methods:
 protected:

  // private methods:
 private:

  // protected members:
 protected:

  // private members:
 private:
  //! Original position
  MVector m_Position;
  //! Original direction
  MVector m_Direction;
  //! Original polarization
  MVector m_Polarization;
  //! Original energy
  double m_Energy;

#ifdef ___CLING___
 public:
  ClassDef(MREAMStartInformation, 0) // no description
#endif

};

#endif


////////////////////////////////////////////////////////////////////////////////
