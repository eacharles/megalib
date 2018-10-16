/*
 * MExceptions.h
 *
 * Copyright (C) by Andreas Zoglauer.
 * All rights reserved.
 *
 * Please see the source-file for the copyright-notice.
 *
 */


#ifndef __MExceptions__
#define __MExceptions__


////////////////////////////////////////////////////////////////////////////////


// Standard libs:

// ROOT libs:
#include <exception>
#include <iostream>
#include <sstream>
#include <cstdlib>
using namespace std;

// MEGAlib libs:
#include "MGlobal.h"
#include "MString.h"

// Forward declarations:


////////////////////////////////////////////////////////////////////////////////


//! Exception base class
class MException : public exception
{
public:
  // Default constructor
  MException() {
    if (m_Abort == true) {
      cout<<what()<<endl<<flush;
      abort();
    }
  }
  //! Default destructor
  virtual ~MException() throw() {}
  
  //! Abort instead of throwing anexception
  static void UseAbort(bool Abort = true) { m_Abort = Abort; }
  
protected:
  //! Flag storing the global variable if exception should abort the program instead of being thrown
  static bool m_Abort;
  
#ifdef ___CLING___
public:
  ClassDef(MException, 1)
#endif  
};


////////////////////////////////////////////////////////////////////////////////


//! This exception is thrown when a parmeter is out of range and 
//! the error cannot be reconvered gracefully
class MExceptionTestFailed : public MException
{
public:
  //! Constructor giving the minium array index, its size, and the accessing index to the array 
  MExceptionTestFailed(MString Description, double Value1, MString Test, double Value2) : MException() {
    ostringstream out;
    out<<"Test failed: "<<Description<<": "<<Value1<<" "<<Test<<" "<<Value2<<endl;
    m_Description = out.str();
  }
  /*
  //! Constructor giving the minium array index, its size, and the accessing index to the array 
  MExceptionTestFailed(MString Description, int Value1, int Value2, MString Test) : MException() {
    ostringstream out;
    out<<"Test failed: "<<Description<<": "<<Value1<<" "<<Test<<" "<<Value2<<endl;
    m_Description = out.str();
  }
  //! Constructor giving the minium array index, its size, and the accessing index to the array 
  MExceptionTestFailed(MString Description, unsigned int Value1, unsigned int Value2, MString Test) : MException() {
    ostringstream out;
    out<<"Test failed: "<<Description<<": "<<Value1<<" "<<Test<<" "<<Value2<<endl;
    m_Description = out.str();
  }
  */
  //! Default destructor
  virtual ~MExceptionTestFailed() throw() {}
  //! The error message
  virtual const char* what() const throw() {
    return m_Description.Data();
  }
  
private:
  //! Description of the test
  MString m_Description;
  
#ifdef ___CLING___
public:
  ClassDef(MExceptionTestFailed, 1)
#endif  
};


////////////////////////////////////////////////////////////////////////////////


//! This exception is thrown when a parmeter is out of range and 
//! the error cannot be reconvered gracefully
class MExceptionParameterOutOfRange : public MException
{
public:
  //! Default constructor
  MExceptionParameterOutOfRange() : MException(), m_IsEmpty(true), m_Value(0), m_Minimum(0), m_Maximum(0), m_Name("") {
  }
  //! Constructor giving the minium array index, its size, and the accessing index to the array 
  MExceptionParameterOutOfRange(double Value, double Min, double Max, MString Name) : MException() {
    SetMinSizeIndex(Value, Min, Max, Name);
  }
  //! Default destructor
  virtual ~MExceptionParameterOutOfRange() throw() {}
  //! Set the data minimum array index, its size, and the accessing index to the array
  void SetMinSizeIndex(double Value, double Min, double Max, MString Name) {
    m_Value = Value; m_Minimum = Min; m_Maximum = Max; m_Name = Name; m_IsEmpty = false;
  }
  //! The error message
  virtual const char* what() const throw() {
    if (m_IsEmpty == false) {
      ostringstream stream;
      stream<<"Parameter \""<<m_Name<<"\" out of range - allowed: ["<<m_Minimum<<".."<<m_Maximum<<"] - you have: "<<m_Value<<endl;
      return stream.str().c_str();
    } else {
      return "Index out of bounds!"; 
    }
  }
  
private:
  //! True if SetValueMinMaxName() has never been called
  bool m_IsEmpty;
  //! The current value
  double m_Value;
  //! The minimum value
  double m_Minimum;
  //! The maximum value
  double m_Maximum;
  //! The name of the parameter
  MString m_Name;
  
  
#ifdef ___CLING___
public:
  ClassDef(MExceptionParameterOutOfRange, 1)
#endif  
};


////////////////////////////////////////////////////////////////////////////////


//! This exception is thrown when an index is out-of-bounds and 
//! the error cannot be recovered gracefully
class MExceptionIndexOutOfBounds : public MException
{
public:
  //! Default constructor
  MExceptionIndexOutOfBounds() : MException(), m_IsEmpty(true), m_Min(0), m_Size(0), m_Index(0) {
  }
  //! Constructor giving the minium array index, its size, and the accessing index to the array 
  MExceptionIndexOutOfBounds(unsigned int Min, unsigned int Size, unsigned int Index) : MException() {
    SetMinSizeIndex(Min, Size, Index);
  }
  //! Default destructor
  virtual ~MExceptionIndexOutOfBounds() throw() {}
  //! Set the data minimum array index, its size, and the accessing index to the array
  void SetMinSizeIndex(unsigned int Min, unsigned int Size, unsigned int Index) {
    m_Min = Min; m_Size = Size; m_Index = Index; m_IsEmpty = false;
  }
  //! The error message
  virtual const char* what() const throw() {
    if (m_IsEmpty == false) {
      ostringstream stream;
      stream<<"Index out of bounds - allowed: ["<<m_Min<<".."<<m_Size<<"[ - you are trying to access with index: "<<m_Index<<endl; 
      return stream.str().c_str();
    } else {
      return "Index out of bounds!"; 
    }
  }
  
private:
  //! True if SetMinSizeIndex() has never been called
  bool m_IsEmpty;
  //! The minimum index
  unsigned int m_Min;
  //! The size of the array, vector, etc.
  unsigned int m_Size;
  //! The index with which we wanted to access the array, vector, etc.
  unsigned int m_Index;
  
  
#ifdef ___CLING___
public:
  ClassDef(MExceptionIndexOutOfBounds, 1)
#endif  
};


////////////////////////////////////////////////////////////////////////////////


//! This exception is thrown when an index is out-of-bounds and 
//! the error cannot be recovered gracefully
class MExceptionEmptyArray : public MException
{
public:
  //! Default constructor
  MExceptionEmptyArray() : MException(), m_Name("") {
  }
  //! Constructor giving the minium array index, its size, and the accessing index to the array 
  MExceptionEmptyArray(MString Name) : MException(), m_Name(Name) {
  }
  //! Default destructor
  virtual ~MExceptionEmptyArray() throw() {}
  //! The error message
  virtual const char* what() const throw() {
    if (m_Name.IsEmpty() == true) {
      return "The array has no elements!"; 
    } else {
      ostringstream stream;
      stream<<"The array \""<<m_Name<<"\" has no elements."<<endl; 
      return stream.str().c_str();
    }
  }
  
private:
  //! Name of the array
  MString m_Name;
  
#ifdef ___CLING___
public:
  ClassDef(MExceptionEmptyArray, 1)
#endif  
};


////////////////////////////////////////////////////////////////////////////////


//! This exception is thrown when an value is out-of-bounds and 
//! the error cannot be recovered gracefully
class MExceptionValueOutOfBounds : public MException
{
public:
  //! Default constructor
  MExceptionValueOutOfBounds() : MException() {
    m_Description = "Value out of bounds!"; 
  }
  //! Constructor giving the minimum and maximum and the given value 
  MExceptionValueOutOfBounds(double Value) : MException() {
    m_Description = "Value out of bounds: ";
    m_Description += Value;
  }
  //! Constructor giving the minimum and maximum and the given value 
  MExceptionValueOutOfBounds(double Min, double Max, double Value) : MException() {
    ostringstream stream;
    stream<<"Value out of bounds - allowed: ["<<Min<<".."<<Max<<"] - your's: "<<Value<<endl; 
    m_Description = stream.str();
  }
  //! Default destructor
  virtual ~MExceptionValueOutOfBounds() throw() {}
  //! The error message
  virtual const char* what() const throw() {
    return m_Description.Data();
  }
  
private:
  //! The description
  MString m_Description;
  
#ifdef ___CLING___
public:
  ClassDef(MExceptionValueOutOfBounds, 1)
#endif  
};


////////////////////////////////////////////////////////////////////////////////


//! This exception is thrown when a division by zero would have occured
class MExceptionDivisionByZero : public MException
{
public:
  //! Default constructor
  MExceptionDivisionByZero() : MException() {}
  //! Default destructor
  virtual ~MExceptionDivisionByZero() throw() {}
  //! The error message
  virtual const char* what() const throw() {
    return "Division by zero!";
  }
  
#ifdef ___CLING___
public:
  ClassDef(MExceptionDivisionByZero, 1)
#endif  
};


////////////////////////////////////////////////////////////////////////////////


//! This exception is thrown when a number is not finite, i.e. nan or inf
class MExceptionNumberNotFinite : public MException
{
public:
  //! Default constructor
  MExceptionNumberNotFinite() : MException() {}
  //! Default destructor
  virtual ~MExceptionNumberNotFinite() throw() {}
  //! The error message
  virtual const char* what() const throw() {
    return "Number not finite!";
  }
  
#ifdef ___CLING___
public:
  ClassDef(MExceptionNumberNotFinite, 1)
#endif  
};


////////////////////////////////////////////////////////////////////////////////


//! This exception is thrown when an index is out-of-bounds and 
//! the error cannot be recovered gracefully
class MExceptionObjectDoesNotExist : public MException
{
public:
  //! Default constructor
  MExceptionObjectDoesNotExist() : MException(), m_IsEmpty(true), m_Name("") {
  }
  //! Standard constructor
  MExceptionObjectDoesNotExist(const MString& Name) : MException(), m_IsEmpty(true), m_Name(Name) {
  }
  //! Default destructor
  virtual ~MExceptionObjectDoesNotExist() throw() {}
  //! Set the data minimum array index, its size, and the accessing index in the array
  void SetName(const MString& Name) { m_Name = Name; }
  //! The error message
  virtual const char* what() const throw() {
    if (m_IsEmpty == false) {
      ostringstream stream;
      stream<<"Object/Key "<<m_Name<<" not found in the list/vector/array!"<<endl; 
      return stream.str().c_str();
    } else {
      return "Object/Key not found in the list/vector/array!"; 
    }
  }
  
private:
  //! True if SetName() has never been called
  bool m_IsEmpty;
  //! The Name of the not found object
  MString m_Name;
  
#ifdef ___CLING___
public:
  ClassDef(MExceptionObjectDoesNotExist, 1)
#endif  
};


////////////////////////////////////////////////////////////////////////////////


//! This exception is thrown when an index is out-of-bounds and 
//! the error cannot be recovered gracefully
class MExceptionObjectsNotIdentical : public MException
{
public:
  //! Default constructor
  MExceptionObjectsNotIdentical() : MException(), m_Name1(""), m_Name2("") {
  }
  //! Standard constructor
  MExceptionObjectsNotIdentical(const MString& Name1, const MString& Name2) : MException(), m_Name1(Name1), m_Name2(Name2) {
  }
  //! Default destructor
  virtual ~MExceptionObjectsNotIdentical() throw() {}
  //! The error message
  virtual const char* what() const throw() {
    if (m_Name1 != "" && m_Name2 != "") {
      ostringstream stream;
      stream<<"The object "<<m_Name1<<" and "<<m_Name2<<" are not identical!"<<endl; 
      return stream.str().c_str();
    } else {
      return "Objects not identical!"; 
    }
  }
  
private:
  //! The name of the first object
  MString m_Name1;
  //! The name of the second object
  MString m_Name2;
  
#ifdef ___CLING___
public:
  ClassDef(MExceptionObjectsNotIdentical, 1)
#endif  
};


////////////////////////////////////////////////////////////////////////////////


//! This exception is thrown when an index is out-of-bounds and 
//! the error cannot be recovered gracefully
class MExceptionUnknownMode : public MException
{
public:
  //! Default constructor
  MExceptionUnknownMode() : MException(), m_IsEmpty(true), m_Type(""), m_Mode("") {
  }
  //! Standard constructor
  MExceptionUnknownMode(const MString& Mode) : MException(), m_IsEmpty(true), m_Type(""), m_Mode(Mode) {
  }
  //! Standard constructor
  MExceptionUnknownMode(const MString& Type, int i) : MException(), m_IsEmpty(true), m_Type(Type) {
    m_Mode += i;
    abort();
  }
  //! Default destructor
  virtual ~MExceptionUnknownMode() throw() {}
  //! Set a name for the unknown mode
  void SetName(const MString& Mode) { m_Mode = Mode; }
  //! The error message
  virtual const char* what() const throw() {
    if (m_IsEmpty == false) {
      ostringstream stream;
      if (m_Type == "") {
        stream<<"Unknown mode "<<m_Mode<<"!"<<endl;
      } else {
        stream<<"Unknown "<<m_Type<<" mode "<<m_Mode<<"!"<<endl;        
      }
      return stream.str().c_str();
    } else {
      return "Unknown mode!"; 
    }
  }

private:
  //! True if SetName() has never been called
  bool m_IsEmpty;
  //! The type of the mode
  MString m_Type;
  //! The mode which is unknown
  MString m_Mode;
  
#ifdef ___CLING___
public:
  ClassDef(MExceptionUnknownMode, 1)
#endif  
};


////////////////////////////////////////////////////////////////////////////////


//! This exception is thrown when an index is out-of-bounds and 
//! the error cannot be recovered gracefully
class MExceptionNeverReachThatLineOfCode : public MException
{
public:
  //! Default constructor
  MExceptionNeverReachThatLineOfCode() : MException(), m_IsEmpty(true) {
  }
  //! Standard constructor
  MExceptionNeverReachThatLineOfCode(const MString& Description) : MException(), m_IsEmpty(false), m_Description(Description) {
  }
  //! Default destructor
  virtual ~MExceptionNeverReachThatLineOfCode() throw() {}
  //! Set a name for the unknown mode
  void SetDescription(const MString& Description) { m_Description = Description; m_IsEmpty = false; }
  //! The error message
  virtual const char* what() const throw() {
    if (m_IsEmpty == false) {
      ostringstream stream;
      if (m_Description == "") {
        stream<<"We should have never reached that line of code: "<<endl;
        stream<<m_Description<<endl;
      } else {
        stream<<"We should have never reached that line of code!"<<endl;        
      }
      return stream.str().c_str();
    } else {
      return "We should have never reached that line of code!"; 
    }
  }
  
private:
  //! True if SetName() has never been called
  bool m_IsEmpty;
  //! The description of what went wrong
  MString m_Description;
  
#ifdef ___CLING___
public:
  ClassDef(MExceptionNeverReachThatLineOfCode, 1)
#endif  
};



////////////////////////////////////////////////////////////////////////////////


//! This exception is thrown when an index is out-of-bounds and 
//! the error cannot be recovered gracefully
class MExceptionArbitrary : public MException
{
public:
  //! Default constructor
  MExceptionArbitrary() : MException(), m_IsEmpty(true) {
  }
  //! Standard constructor
  MExceptionArbitrary(const MString& Description) : MException(), m_IsEmpty(false), m_Description(Description) {
  }
  //! Default destructor
  virtual ~MExceptionArbitrary() throw() {}
  //! Set a name for the unknown mode
  void SetDescription(const MString& Description) { m_Description = Description; m_IsEmpty = false; }
  //! The error message
  virtual const char* what() const throw() {
    if (m_IsEmpty == false) {
      ostringstream stream;
      if (m_Description != "") {
        stream<<"An exception was triggered: "<<endl;
        stream<<m_Description<<endl;
      } else {
        stream<<"An exception was triggered!"<<endl;        
      }
      return stream.str().c_str();
    } else {
      return "An exception was triggered!"; 
    }
  }
  
private:
  //! True if SetName() has never been called
  bool m_IsEmpty;
  //! The description of what went wrong
  MString m_Description;
  
#ifdef ___CLING___
public:
  ClassDef(MExceptionArbitrary, 1)
#endif  
};


#endif


////////////////////////////////////////////////////////////////////////////////
