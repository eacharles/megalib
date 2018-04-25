/*
 * BackgroundGenerator.cxx
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

// Standard
#include <iostream>
#include <string>
#include <sstream>
#include <csignal>
#include <iomanip>
#include <cmath>
#include <cstdlib>
using namespace std;

// ROOT
#include <TROOT.h>
#include <TEnv.h>
#include <TSystem.h>
#include <TApplication.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH2.h>

// MEGAlib
#include "MGlobal.h"
#include "MStreams.h"


/******************************************************************************/

class BackgroundGenerator
{
public:
  //! Default constructor
  BackgroundGenerator();
  //! Default destructor
  ~BackgroundGenerator();

  //! Parse the command line
  bool ParseCommandLine(int argc, char** argv);
  //! Analyze what eveer needs to be analyzed...
  bool Analyze();
  //! Interrupt the analysis
  void Interrupt() { m_Interrupt = true; }

protected:
  //! Generate a cosmic photo spectrum according to Gruber et al.
  bool GenerateCosmicPhotonGruber();
  //! Generate a cosmic proton spectrum from SPENVIS files
  bool GenerateCosmicProtonsSpenvis();
  //! Generate a cosmic alpha spectrum from SPENVIS files
  bool GenerateCosmicAlphasSpenvis();
  //! Generate a cosmic electron spectrum according to Mizuno
  bool GenerateCosmicElectronsMizuno();
  //! Generate a cosmic positron spectrum according to Mizuno
  bool GenerateCosmicPositronsMizuno();


    //! Generate a cosmic proton spectrum from AMS data
  bool GenerateCosmicProtonsAguilar();
  //! Generate a cosmic alpha spectrum from AMS data
  bool GenerateCosmicAlphasAguilar();
  //! Generate a cosmic electron spectrum AMS data
  bool GenerateCosmicElectronsAguilar();
  //! Generate a cosmic positron spectrum AMS data
  bool GenerateCosmicPositronsAguilar();
  //! Generate a albedo photo spectrum according to SazonovChurazovMizunoAbdo
  bool GenerateCosmicPhotonTurlerMizunoAckermann();
  //! Generate a albedo photo spectrum according to SazonovChurazovMizunoAbdo
  bool GenerateAlbedoPhotonsSazonovChurazovMizunoAbdo();

  //! Generate a trapped proton and electron spectrum from SPENVIS
  bool GenerateTrappedProtonsElectronsSpenvis();

  //! Generate the albedo photons according to Ajello and Mizuno et al
  bool GenerateAlbedoPhotonsAjelloMizuno();
  //! Generate the albedo photons according to Tuerler, Mizuno & Abdo et al
  bool GenerateAlbedoPhotonsTuerlerMizunoAbdo();
  //! Generate the albedo photon lines according to Harris et al
  bool GenerateAlbedoPhotonLinesHarris();
  //! Generate the albedo electron spectrum according to Alcaracz 2000 (AMS) and Mizuno
  bool GenerateAlbedoElectronsAlcarazMizuno();
  //! Generate the albedo positron spectrum according to Alcaracz 2000 (AMS) and Mizuno
  bool GenerateAlbedoPositronsAlcarazMizuno();
  //! Generate the albedo proton spectrum according to Alcaracz 2000 (AMS)
  bool GenerateAlbedoProtonsAlcaraz();
  //! Generate the albedo proton spectrum according to Mizuno 2004 (AMS)
  bool GenerateAlbedoProtonsAveragedMizuno();
  //! Generate the albedo proton spectrum according to Alcaracz 2000 (AMS)
  bool GenerateAlbedoProtonsAveragedAlcaraz();
  //! Generate the albedo neutron spectrum according to Kole+ 2014 and Morris
  bool GenerateAlbedoNeutronsMorrisKole();
  //! Generate the albedo neutron spectrum according to Kole+ 2014
  bool GenerateAlbedoNeutronsKole();


  //!
  bool IsValid(double ModelMin, double ModelMax);

  bool WriteSourceFile(MString SourceName, double Flux, int ParticleID, MString Comments = "", MString SpectralString = "");
  bool WriteEnergyFile(MString SourceName, vector<double> Spectrum);
  bool WriteAngleFile(MString SourceName, vector<double> Angle);

  //! Write the final summary files
  bool WriteSummaryFiles();

private:
  //! True, if the analysis needs to be interrupted
  bool m_Interrupt;

  //! True, if we are in MGGPOD mode
  bool m_IsMGGPOD;

  //! True, if we want to use different models
  bool m_PlanB;

  //! True if we are in orbit around Earth
  bool m_IsEarthOrbit;

  //! The orbit altitude in kilometers
  double m_Altitude;
  //! The orbit inclination in degrees (equatorial = 0 degree)
  double m_Inclination;

  //! Angle from zenith to earth horizon in degree
  double m_HorizonAngle;
  //! Average Earth radius in km
  double m_EarthRadius;
  //! Height of atmosphere in km
  double m_AtmosphereHeight;
  //! Average geomagnetic cutoff
  double m_AverageGeomagneticCutOff;

  //! Minimum energy range in keV
  double m_EnergyMin;
  //! Maximum energy range in keV
  double m_EnergyMax;

  //! The number of energy Bins
  unsigned int m_NEnergyBins;
  //! The energy binning
  vector<double> m_EnergyBins;

  //! The number of angle Bins
  unsigned int m_NAngleBins;
  //! The angle binning
  vector<double> m_AngleBins;

  //! A trapped proton and electrons file from Spenvis
  MString m_AverageTrappedProtonsElectronsSpenvis;

  //! A cosmic proton file from Spenvis
  MString m_AverageCosmicProtonsSpenvis;

  //! A cosmic alphas file from Spenvis
  MString m_AverageCosmicAlphasSpenvis;

  //! The photonic summary files
  vector<MString> m_PhotonicSummaryFiles;

  //! The leptonic summary files
  vector<MString> m_LeptonicSummaryFiles;

  //! The hadronic summary files
  vector<MString> m_HadronicSummaryFiles;

  //! The trapped hadronic summary files
  vector<MString> m_TrappedHadronicSummaryFiles;

};

/******************************************************************************/


/******************************************************************************
 * Default constructor
 */
BackgroundGenerator::BackgroundGenerator() : m_Interrupt(false)
{
  gStyle->SetPalette(1, 0);

  m_IsEarthOrbit = false;
  m_Altitude = 575.0; // km
  m_Inclination = 6.0; // deg

  m_HorizonAngle = 0.0; // Calculated later
  m_EarthRadius = 6370; // km
  m_AtmosphereHeight = 40; // km
  m_AverageGeomagneticCutOff = 12.6; // GV

  m_EnergyMin = 5; // keV
  m_EnergyMax = 100000000; // keV

  m_NEnergyBins = 100;
  m_NAngleBins = 100;

  m_IsMGGPOD = false;

  m_PlanB = false;
}


/******************************************************************************
 * Default destructor
 */
BackgroundGenerator::~BackgroundGenerator()
{
  // Intentionally left blanck
}


/******************************************************************************
 * Parse the command line
 */
bool BackgroundGenerator::ParseCommandLine(int argc, char** argv)
{
  ostringstream Usage;
  Usage<<endl;
  Usage<<"  Usage: BackgroundGenerator <options>"<<endl;
  Usage<<"    General options:"<<endl;
  Usage<<"         -e <d> <d>:   energy range (min/max) as two doubles"<<endl;
  Usage<<"         -b <i>:       number of energy and angle bins"<<endl;
  Usage<<"         -o <d> <d>:   orbit altitude (in km), inclination (in degree), and average geomagnetic cutoff (in GV) as three doubles - not giving these value assumes no orbit arround Earth"<<endl;
  Usage<<"         -tps <file>:  A trapped proton file from SPENVIS"<<endl;
  Usage<<"         -cps <file>:  A cosmic proton file from SPENVIS"<<endl;
  Usage<<"         -cas <file>:  A cosmic alpha particle file from SPENVIS"<<endl;
  Usage<<"         -pb: Generating background using a different model"<<endl;
  Usage<<"         -h:           print this help"<<endl;
  Usage<<endl;

  string Option;

  // Check for help
  for (int i = 1; i < argc; i++) {
    Option = argv[i];
    if (Option == "-h" || Option == "--help" || Option == "?" || Option == "-?") {
      cout<<Usage.str()<<endl;
      return false;
    }
  }

  // Now parse the command line options:
  for (int i = 1; i < argc; i++) {
    Option = argv[i];

    // First check if each option has sufficient arguments:
    // Single argument
    if (Option == "-f" || Option == "-b" || Option == "-tps" || Option == "-cps" || Option == "-cas") {
      if (!((argc > i+1) &&
            (argv[i+1][0] != '-' || isalpha(argv[i+1][1]) == 0))){
        cout<<"Error: Option "<<argv[i][1]<<" needs a second argument!"<<endl;
        cout<<Usage.str()<<endl;
        return false;
      }
    }
    // Multiple arguments template
    else if (Option == "-e") {
      if (!((argc > i+2) &&
            (argv[i+1][0] != '-' || isalpha(argv[i+1][1]) == 0) &&
            (argv[i+2][0] != '-' || isalpha(argv[i+2][1]) == 0))){
        cout<<"Error: Option "<<argv[i][1]<<" needs two arguments!"<<endl;
        cout<<Usage.str()<<endl;
        return false;
      }
    } else if (Option == "-o") {
      if (!((argc > i+3) &&
            (argv[i+1][0] != '-' || isalpha(argv[i+1][1]) == 0) &&
            (argv[i+2][0] != '-' || isalpha(argv[i+2][1]) == 0) &&
            (argv[i+3][0] != '-' || isalpha(argv[i+3][1]) == 0))) {
        cout<<"Error: Option "<<argv[i][1]<<" needs three arguments!"<<endl;
        cout<<Usage.str()<<endl;
        return false;
      }
    }

    // Then fulfill the options:
    if (Option == "-f") {
      //m_FileName = argv[++i];GenerateAlbedoNeutronsArmstrongMorris
      //cout<<"Accepting file name: "<<m_FileName<<endl;
    } else if (Option == "-b") {
      m_NEnergyBins = atoi(argv[++i]);
      m_NAngleBins = m_NEnergyBins;
      cout<<"Accepting "<<m_NEnergyBins<<" bins"<<endl;
    } else if (Option == "-tps") {
      m_AverageTrappedProtonsElectronsSpenvis = argv[++i];
      cout<<"Accepting trappeds proton and electrons file from SPENVIS: "<<m_AverageTrappedProtonsElectronsSpenvis<<endl;
    } else if (Option == "-cps") {
      m_AverageCosmicProtonsSpenvis = argv[++i];
      cout<<"Accepting cosmic proton file from SPENVIS: "<<m_AverageCosmicProtonsSpenvis<<endl;
    } else if (Option == "-cas") {
      m_AverageCosmicAlphasSpenvis = argv[++i];
      cout<<"Accepting cosmic alphas file from SPENVIS: "<<m_AverageCosmicAlphasSpenvis<<endl;
    } else if (Option == "-e") {
      m_EnergyMin = atof(argv[++i]);
      m_EnergyMax = atof(argv[++i]);
      cout<<"Accepting minimum energy ("<<m_EnergyMin<<") and maximum energy ("<<m_EnergyMax<<")"<<endl;
    } else if (Option == "-o") {
      m_Altitude = atof(argv[++i]);
      m_Inclination = atof(argv[++i]);
      m_AverageGeomagneticCutOff = atof(argv[++i]);
      m_IsEarthOrbit = true;
      cout<<"Accepting altitude ("<<m_Altitude<<" km), inclination ("<<m_Inclination<<"), and average geomagnetic cutoff ("<<m_AverageGeomagneticCutOff<<")"<<endl;
    } else if (Option == "-m") {
      m_IsMGGPOD = true;
      cout<<"Generating MGGPOD output format"<<endl;
    } else if (Option == "-pb") {
      m_PlanB = true;
      cout<<"Generating background using a different model for:"<<endl;
      cout<<"* Primary Protons"<<endl;
      cout<<"* Primary Alphas"<<endl;
      cout<<"* Primary Electrons"<<endl;
      cout<<"* Primary Positrons"<<endl;
      cout<<"* Cosmic Photons"<<endl;
      cout<<"* Albedo Photons"<<endl;
    } else {
      cout<<"Error: Unknown option \""<<Option<<"\"!"<<endl;
      cout<<Usage.str()<<endl;
      return false;
    }
  }

  return true;
}


/******************************************************************************
 * Do whatever analysis is necessary
 */
bool BackgroundGenerator::Analyze()
{
  if (m_Interrupt == true) return false;

  if (m_IsEarthOrbit == true) {
    m_HorizonAngle = 90.0 + acos((m_EarthRadius + m_AtmosphereHeight)/(m_EarthRadius+m_Altitude))*c_Deg;
  } else {
    m_HorizonAngle = 180.0;
  }

  double Min = m_EnergyMin;
  if (Min <= 0) Min = 1;
  Min = log(Min);
  double Max = log(m_EnergyMax);
  double Dist = (Max-Min)/(m_NEnergyBins);
  for (unsigned int i = 0; i <= m_NEnergyBins; ++i) {
    m_EnergyBins.push_back(exp(Min + i*Dist));
  }
  m_NEnergyBins = m_EnergyBins.size();

  Min = 0;
  Max = 180;
  Dist = (Max-Min)/(m_NAngleBins);
  for (unsigned int i = 0; i <= m_NAngleBins; ++i) {
    m_AngleBins.push_back(Min + i*Dist);
  }

  if (m_PlanB == false){
    if (GenerateCosmicPhotonGruber() == false) return false;
    if (GenerateCosmicElectronsMizuno() == false) return false;
    if (GenerateCosmicPositronsMizuno() == false) return false;


    if (m_AverageCosmicProtonsSpenvis != "") {
      if (GenerateCosmicProtonsSpenvis() == false) return false;
    }
    if (m_AverageCosmicAlphasSpenvis != "") {
      if (GenerateCosmicAlphasSpenvis() == false) return false;
    }
    if (m_AverageTrappedProtonsElectronsSpenvis != "") {
      if (GenerateTrappedProtonsElectronsSpenvis() == false) return false;
    }

    if (m_IsEarthOrbit == true) {
      if (GenerateAlbedoPhotonsTuerlerMizunoAbdo() == false) return false;
      if (GenerateAlbedoPhotonLinesHarris() == false) return false;
      if (GenerateAlbedoElectronsAlcarazMizuno() == false) return false;
      if (GenerateAlbedoPositronsAlcarazMizuno() == false) return false;
      if (GenerateAlbedoProtonsAveragedMizuno() == false) return false;
      if (GenerateAlbedoNeutronsKole() == false) return false;
    }
  } else {
    if (GenerateCosmicPhotonTurlerMizunoAckermann() == false) return false;
    if (GenerateCosmicElectronsAguilar() == false) return false;
    if (GenerateCosmicPositronsAguilar() == false) return false;

    if (GenerateCosmicProtonsAguilar() == false) return false;

    if (GenerateCosmicAlphasAguilar() == false) return false;

    if (m_AverageTrappedProtonsElectronsSpenvis != "") {
      if (GenerateTrappedProtonsElectronsSpenvis() == false) return false;
    }

    if (m_IsEarthOrbit == true) {
      if (GenerateAlbedoPhotonsSazonovChurazovMizunoAbdo() == false) return false;
      if (GenerateAlbedoPhotonLinesHarris() == false) return false;
      if (GenerateAlbedoElectronsAlcarazMizuno() == false) return false;
      if (GenerateAlbedoPositronsAlcarazMizuno() == false) return false;
      if (GenerateAlbedoProtonsAveragedMizuno() == false) return false;
      if (GenerateAlbedoNeutronsKole() == false) return false;
    }
  }

  WriteSummaryFiles();

  return true;
}


/******************************************************************************
 * Generate a cosmic photon spectrum after
 * Türler et al. 2010
 *   doi:10.1051/0004-6361/200913072
 * Mizuno et al. 2004
 *   http://stacks.iop.org/0004-637X/614/i=2/a=1113
 * Ackermann et al. 2015
 *  doi:10.1088/0004-637X/799/1/86
 */
bool BackgroundGenerator::GenerateCosmicPhotonTurlerMizunoAckermann()
{
  mout<<endl;
  mout<<"Generating a cosmic photon spectrum according to Turler, Mizuno, Ackermann..."<<endl;

  if (IsValid(1, 1000000000) == false) {
    return false;
  }

  auto Turler = [this] (double EnergykeV) {
    double Flux = 0.109 / (pow(EnergykeV/28,1.4)+pow(EnergykeV/28,2.88));
    return Flux;
  };

  auto Mizuno = [this] (double EnergykeV) {
    double Flux = 40.*pow(EnergykeV/1000, -2.15)/(10000000);
    return Flux;
  };

  auto Ackermann = [this] (double EnergykeV) {
    double I100 = 0.95*pow(10,-7)/1000;
    double gamma = 2.32;
    double Ecut = 279*1000000;

    double Flux = I100 * pow(EnergykeV/(100*1000),-gamma)*exp(EnergykeV/Ecut);
    return Flux;
  };

  auto TurlerMizunoAckermann = [Turler, Ackermann, Mizuno] (double EnergykeV) {
    if (EnergykeV < 890) return Turler(EnergykeV);
    else if (EnergykeV > 1200) return Ackermann(EnergykeV);
    else return Mizuno(EnergykeV);
  };


  vector<double> Spectrum;
  for (unsigned int b = 0; b < m_NEnergyBins; ++b) {
    Spectrum.push_back(TurlerMizunoAckermann(m_EnergyBins[b]));
  }

  vector<double> Angle;
  for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
    if (m_AngleBins[b] < m_HorizonAngle) {
      Angle.push_back(1.0);
    } else {
      Angle.push_back(0.0);
    }
  }

  // Calculate flux:

  // Integration factor over sphere
  double AngleFactor = 2*c_Pi * (1 - cos(m_HorizonAngle*c_Rad));

  // Determine the flux in ph/s/cm2 via numerical integration...
  double Flux = 0;
  int Bins = 10000;
  double Min = log(m_EnergyMin);
  double Max = log(m_EnergyMax);
  double Dist = (Max-Min)/Bins;
  for (double e = 0; e <= Bins-1; ++e) {
    double EMin = exp(Min + e*Dist);
    double EMax = exp(Min + (e+1)*Dist);
    double Average = 0.5*(TurlerMizunoAckermann(EMin) + TurlerMizunoAckermann(EMax));
    Flux += Average*(EMax-EMin);
  }
  cout<<"Anglefactor: "<<AngleFactor<<endl;
  cout<<"Flux: "<<Flux<<" ph/cm2/s/sr"<<endl;
  Flux *= AngleFactor;
  cout<<"Average flux: "<<Flux<<" ph/cm2/s"<<endl;




  WriteEnergyFile("CosmicPhotonTurlerMizunoAckermann", Spectrum);
  WriteAngleFile("CosmicPhotonTurlerMizunoAckermann", Angle);
  WriteSourceFile("CosmicPhotonTurlerMizunoAckermann", Flux, 1);

  m_PhotonicSummaryFiles.push_back("CosmicPhotonTurlerMizunoAckermann");

  return true;
}

/******************************************************************************
 * Generate an albedo photon spectrum after
 * Sazonov et al. 2007
 *   doi:10.1111/j.1365-2966.2007.11746.x
 * Churazov et al. 2006
 *   doi:10.1111/j.1365-2966.2008.12918.x
 * Türler et al. 2010
 *   doi:10.1051/0004-6361/200913072
 * Mizuno et al. 2004
 *   http://stacks.iop.org/0004-637X/614/i=2/a=1113
 * Abdo et al. 2009
 *   doi:10.1103/PhysRevD.80.122004
 * Mizuno is used as the absolute normalization
 */
bool BackgroundGenerator::GenerateAlbedoPhotonsSazonovChurazovMizunoAbdo()
{
  // Sazonov low energy component
  auto Sazonov = [this](double EnergykeV) {

    if (m_IsEarthOrbit == true) {
      m_HorizonAngle = 90.0 + acos((m_EarthRadius + m_AtmosphereHeight)/(m_EarthRadius+m_Altitude))*c_Deg;
    } else {
      m_HorizonAngle = 180.0;
    }

    double thetamax = 180 - m_HorizonAngle;
    double cosomega = cos(thetamax*c_Rad);
    double phi = 650. / 1000.;

    double num = 1.47*0.0178/(pow(phi/2.8,0.4)+pow(phi/2.8,1.5));
    double den = sqrt(1+pow(m_AverageGeomagneticCutOff/(1.3*pow(phi,0.25)*(1+2.5*pow(phi,0.4))),2));
    double fac = 3*cosomega*(1+cosomega)/5*c_Pi;

    double c = fac*num/den;

    return c / (pow(EnergykeV/44,-5)+pow(EnergykeV/44,1.4));
  };

  // Tuerler primary low energy component
  auto TuerlerCosmicPhotons = [](double EnergykeV) {
    return 0.109 / (pow(EnergykeV/28,1.4)+pow(EnergykeV/28, 2.88));
  };
  // Churazov low energy component
  auto Churazov = [this, TuerlerCosmicPhotons](double EnergykeV) {

    if (m_IsEarthOrbit == true) {
      m_HorizonAngle = 90.0 + acos((m_EarthRadius + m_AtmosphereHeight)/(m_EarthRadius+m_Altitude))*c_Deg;
    } else {
      m_HorizonAngle = 180.0;
    }

    double thetamax = 180 - m_HorizonAngle; // deg max polar angle wrt nadir
    double omega = 2*c_Pi*(1-cos(thetamax*c_Rad));

    double first = 1.22/(pow(EnergykeV/28.5, -2.54)+pow(EnergykeV/51.3, 1.57)-0.37);
    double second = (2.93+pow(EnergykeV/3.08, 4))/(1+pow(EnergykeV/3.08,4));
    double third = (0.123+pow(EnergykeV/91.83,3.44))/(1+pow(EnergykeV/91.83, 3.44));

    return omega*TuerlerCosmicPhotons(EnergykeV)*first*second*third;
  };

  auto ChurazovSazonov = [Churazov, Sazonov](double EnergykeV) {
    return Churazov(EnergykeV)+Sazonov(EnergykeV);
  };

  // Mizuno downward component
  auto Mizuno = [](double EnergykeV) {
    // Valid above ~ 1 MeV
    if (EnergykeV < 20000) {
      return 1010.0*pow(EnergykeV/1000, -1.34) / 1000.0 / 10000.0;
    } else if (EnergykeV >= 20000 && EnergykeV < 1000000) {
      return 7290.0*pow(EnergykeV/1000, -2.0) / 1000.0 / 10000.0;
    } else if (EnergykeV >= 1000000) {
      return 29000*pow(EnergykeV/1000, -2.2) / 1000.0 / 10000.0;
    }
    return 0.0;
  };

  // Abdo FERMI
  auto Abdo = [](double EnergykeV) {
    // Valid from 200 MeV & up
    return 1.823e-8*pow(EnergykeV/200000, -2.8);
  };

  auto ChurazovSazonovMizunoAbdo = [ChurazovSazonov, Mizuno, Abdo](double EnergykeV, double ScalerChurazovSazonov, double ScalerMizuno, double ScalerAbdo) {
    if (EnergykeV < 1850) {
      return ChurazovSazonov(EnergykeV)*ScalerChurazovSazonov;
    } else if (EnergykeV > 200000) {
      return Abdo(EnergykeV)*ScalerAbdo;
    } else {
      return Mizuno(EnergykeV)*ScalerMizuno;
    }
  };

  if (m_IsEarthOrbit == false) return true;

  mout<<endl;
  mout<<"Generating an albedo photon spectrum according to Churazov Sazonov, Mizuno, Abdo..."<<endl;

  if (IsValid(1, 1000000000) == false) {
    return false;
  }


  vector<double> Angle;
  for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
    if (m_AngleBins[b] >= m_HorizonAngle) {
      Angle.push_back(1.0);
    } else {
      Angle.push_back(0.0);
    }
  }

  // Scaling from Mizuno et al.
  double Rcut_desired = m_AverageGeomagneticCutOff;
  double Rcut_Mizuno = 4.5;
  double ScalerMizuno = pow(Rcut_desired/Rcut_Mizuno, -1.13);

  // Make sure to scale the Churazov & Sazonov result to the Mizuno result at 1.85 MeV:
  double MizunoValue = ScalerMizuno*Mizuno(1850);
  double ChurazovSazonovValue = ChurazovSazonov(1850);
  double ScalerChurazovSazonov= MizunoValue/ChurazovSazonovValue;

  MizunoValue = ScalerMizuno*Mizuno(200000);
  double AbdoValue = Abdo(200000);
  double ScalerAbdo = MizunoValue/AbdoValue;

  vector<double> Spectrum;
  for (unsigned int b = 0; b < m_NEnergyBins; ++b) {
    Spectrum.push_back(ChurazovSazonovMizunoAbdo(m_EnergyBins[b], ScalerChurazovSazonov, ScalerMizuno, ScalerAbdo));
    //cout<<m_EnergyBins[b]<<":"<<ChurazovSazonovMizunoAbdo(m_EnergyBins[b], ScalerChurazovSazonov, ScalerMizuno, ScalerAbdo)<<endl;
  }


  // Determine the flux in ph/s/cm2 via numerical integration...
  double Flux = 0;
  int Bins = 10000;
  double Min = log(m_EnergyMin);
  double Max = log(m_EnergyMax);
  double Dist = (Max-Min)/Bins;
  for (double e = 0; e <= Bins-1; ++e) {
    double EMin = exp(Min + e*Dist);
    double EMax = exp(Min + (e+1)*Dist);
    double Average = 0.5*(ChurazovSazonovMizunoAbdo(EMin, ScalerChurazovSazonov, ScalerMizuno, ScalerAbdo) + ChurazovSazonovMizunoAbdo(EMax, ScalerChurazovSazonov, ScalerMizuno, ScalerAbdo));
    Flux += Average*(EMax-EMin);
  }
  // Integration factor over sphere
  double AngleFactor = 2*c_Pi * (cos(m_HorizonAngle*c_Rad) - (-1));

  cout<<"Angle-factor: "<<AngleFactor<<endl;
  cout<<"Flux: "<<Flux<<" ph/cm2/s/sr"<<endl;
  Flux *= AngleFactor;
  cout<<"Flux: "<<Flux<<" ph/cm2/s"<<endl;

  MString Comments;
  Comments += "# Albedo photons determinined after Sazonov et al. 2007 & Churazov et al. 2006 (Turler+ 2010), Mizuno+ 2004, Abdo+ 2009\n";
  Comments += "# Assumptions: \n";
  Comments += "# (1) Use Churazov+Sazonov below 1850 keV, Abdo+ above 200 MeV, Mizuno in between \n";
  Comments += "# (2) Everything is normalized to Mizuno \n";
  Comments += "# (3) Geomagnetic cutoff difference scaled as in Mizuno \n";
  Comments += "# (4) Angular distribution is flat out to the Earth-horizon (no limb brightening!) \n";
  Comments += "# (5) Limb brightening and darkening is NOT included\n";
  Comments += "# (6) There should be a bump form the 511-keV-Compton scatters belog ~450 keV which is NOT included \n";

  WriteEnergyFile("AlbedoPhotonsChurazovSazonovMizunoAbdo", Spectrum);
  WriteAngleFile("AlbedoPhotonsChurazovSazonovMizunoAbdo", Angle);
  WriteSourceFile("AlbedoPhotonsChurazovSazonovMizunoAbdo", Flux, 1, Comments);

  m_PhotonicSummaryFiles.push_back("AlbedoPhotonsChurazovSazonovMizunoAbdo");

  return true;
}

/******************************************************************************
 * Generate a cosmic alpha particle spectrum from AMS data
 * Aguilar et al. 2015b
 *   doi:10.1103/PhysRevLett.115.211101
 */
bool BackgroundGenerator::GenerateCosmicAlphasAguilar()
{
  mout<<endl;
  mout<<"Generating a cosmic alpha particle spectrum using AMS data..."<<endl;

  if (IsValid(1, 1000000000) == false) {
    return false;
  }

  ifstream in;
  MString m_AverageCosmicAlphasAguilar = "./AguilarAlphas.dat";
  in.open(m_AverageCosmicAlphasAguilar);
  if (in.is_open() == false) {
    mout<<"Error: Unable to open file "<<m_AverageCosmicAlphasAguilar<<endl;
    return false;
  }

  vector<double> Energies;
  vector<double> Fluxes;

  double E0 = 2*(0.938 + 0.940);

  MString Line;
  bool Start = false;
  while (in.good() == true) {
    Line.ReadLine(in);
    if (Line.Length() < 2) continue;
    if (Line.Contains("Flux") == true) {
      Start = true;
      continue;
    }
    if (Line.Contains("End of File") == true) break;
    if (Start == false) continue;
    vector<MString> Tokens = Line.Tokenize(" ");
    if (Tokens.size() != 2) {
      cout<<"Error: The line should have two tokens, rigidity, and differential flux"<<endl;
      cout<<Line<<endl;
      continue;
    }

    double Flux = atof(Tokens[1])*atof(Tokens[0]);
    double Energy = 4*(sqrt(pow(E0,2)+pow((atof(Tokens[0])/2),2))-E0)*1000000;
    Flux = Flux/Energy/10000;

    if (Flux > 0) {
      //cout<<Line<<endl;
      Energies.push_back(Energy);
      Fluxes.push_back(Flux);
    }
  }

  auto Interp = [this, Energies, Fluxes] (double EnergykeV) {

    //solar modulation potential (MV): ~550 for solar minimum ~1100 for solar maximum
    double solmod = 650.;

    double E0 = 2*(0.938 + 0.940);

    // Geomagnetic modulation factor from Mizuno et al. 2004
    double Rigidity = sqrt(EnergykeV*EnergykeV/1000000/1000000 + 2*EnergykeV*E0/1000000);
    double redfac = 1/(1+pow(Rigidity/m_AverageGeomagneticCutOff,-12.0));

    //Solar modulation factor from Gleeson & Axford 1968
    double solmodfac = ((EnergykeV/1000000+E0)*(EnergykeV/1000000+E0)-E0*E0)/(
                (EnergykeV/1000000+E0+solmod/1000)*(EnergykeV/1000000+E0+solmod/1000)-E0*E0);

    double modfac = solmodfac * redfac;

    double Flux = 0.0;
    if (EnergykeV < Energies[0]) {
      double m = (log(modfac*Fluxes[1]) - log(modfac*Fluxes[0]))/(log(Energies[1]) - log(Energies[0]));
      double t = log(modfac*Fluxes[1]) - m * log(Energies[1]);
      double logy = m * log(EnergykeV) + t;
      Flux = exp(logy);
    } else if (EnergykeV > Energies.back()) {
      int Last = Energies.size() - 1;
      double m = (log(modfac*Fluxes[Last]) - log(modfac*Fluxes[Last-1]))/(log(Energies[Last]) - log(Energies[Last-1]));
      double t = log(modfac*Fluxes[Last]) - m * log(Energies[Last]);
      double logy = m * log(EnergykeV) + t;
      Flux = exp(logy);
    } else {
      for (unsigned int e = 1; e < Energies.size(); ++e) {
        if (EnergykeV < Energies[e]) {
          double m = (log(modfac*Fluxes[e]) - log(modfac*Fluxes[e-1]))/(log(Energies[e]) - log(Energies[e-1]));
          double t = log(modfac*Fluxes[e]) - m * log(Energies[e]);
          double logy = m * log(EnergykeV) + t;
          Flux = exp(logy);
          break;
        }
      }
    }

    return Flux;
  };



  vector<double> Angle;
  for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
    if (m_AngleBins[b] <= m_HorizonAngle) {
      Angle.push_back(1.0);
    } else {
      Angle.push_back(0.0);
    }
  }

  vector<double> Spectrum;
  for (unsigned int b = 0; b < m_NEnergyBins; ++b) {
    Spectrum.push_back(Interp(m_EnergyBins[b]));
  }


  // Determine the flux in ph/s/cm2 via numerical integration...
  double Flux = 0;
  int Bins = 10000;
  double Min = log(m_EnergyMin);
  double Max = log(m_EnergyMax);
  double Dist = (Max-Min)/Bins;
  for (double e = 0; e <= Bins-1; ++e) {
    double EMin = exp(Min + e*Dist);
    double EMax = exp(Min + (e+1)*Dist);
    double Average = 0.5*(Interp(EMin) + Interp(EMax));
    Flux += Average*(EMax-EMin);
  }
  // Integration factor over sphere
  double AngleFactor = 2*c_Pi * (1 - cos(m_HorizonAngle*c_Rad));

  cout<<"Angle-factor: "<<AngleFactor<<endl;
  cout<<"Flux: "<<Flux<<" particle/cm2/s/sr"<<endl;
  Flux *= AngleFactor;
  cout<<"Flux: "<<Flux<<" particle/cm2/s"<<endl;

  MString Comments;
  Comments += "# Cosmic alpha particles using AMS data\n";
  Comments += "# Assumptions: \n";
  Comments += "# * log-log extrapolation from AMS data\n";
  Comments += "# * Angular distribution is flat out to the Earth-horizon \n";

  WriteEnergyFile("CosmicAlphasAguilar", Spectrum);
  WriteAngleFile("CosmicAlphasAguilar", Angle);
  WriteSourceFile("CosmicAlphasAguilar", Flux, 21, Comments);

  m_HadronicSummaryFiles.push_back("CosmicAlphasAguilar");

  return true;
}

/******************************************************************************
 * Generate a cosmic proton particle spectrum from AMS data
 * Aguilar et al. 2015
 *   doi:10.1103/PhysRevLett.114.171103
 */
bool BackgroundGenerator::GenerateCosmicProtonsAguilar()
{
  mout<<endl;
  mout<<"Generating a cosmic protons spectrum using AMS data..."<<endl;

  if (IsValid(1, 1000000000) == false) {
    return false;
  }

  ifstream in;
  MString m_AverageCosmicProtonsAguilar = "./AguilarProton.dat";
  in.open(m_AverageCosmicProtonsAguilar);
  if (in.is_open() == false) {
    mout<<"Error: Unable to open file "<<m_AverageCosmicProtonsAguilar<<endl;
    return false;
  }

  vector<double> Energies;
  vector<double> Fluxes;

  double E0 = 0.938;

  MString Line;
  bool Start = false;
  while (in.good() == true) {
    Line.ReadLine(in);
    if (Line.Length() < 2) continue;
    if (Line.Contains("Flux") == true) {
      Start = true;
      continue;
    }
    if (Line.Contains("End of File") == true) break;
    if (Start == false) continue;
    vector<MString> Tokens = Line.Tokenize(" ");
    if (Tokens.size() != 2) {
      cout<<"Error: The line should have two tokens, rigidity, and differential flux"<<endl;
      cout<<Line<<endl;
      continue;
    }

    double Flux = atof(Tokens[1])*atof(Tokens[0]);
    double Energy = (sqrt(pow(E0,2)+pow((atof(Tokens[0])/2),2))-E0)*1000000;
    Flux = Flux/Energy/10000;

    if (Flux > 0) {
      //cout<<Line<<endl;
      Energies.push_back(Energy);
      Fluxes.push_back(Flux);
    }
  }

  auto Interp = [this, Energies, Fluxes] (double EnergykeV) {

    //solar modulation potential (MV): ~550 for solar minimum ~1100 for solar maximum
    double solmod = 650.;

    double E0 = 0.938;

    // Geomagnetic modulation factor from Mizuno et al. 2004
    double Rigidity = sqrt(EnergykeV*EnergykeV/1000000/1000000 + 2*EnergykeV*E0/1000000);
    double redfac = 1/(1+pow(Rigidity/m_AverageGeomagneticCutOff,-12.0));

    //Solar modulation factor from Gleeson & Axford 1968
    double solmodfac = ((EnergykeV/1000000+E0)*(EnergykeV/1000000+E0)-E0*E0)/(
                (EnergykeV/1000000+E0+solmod/1000)*(EnergykeV/1000000+E0+solmod/1000)-E0*E0);

  double   modfac = solmodfac * redfac;

    double Flux = 0.0;
    if (EnergykeV < Energies[0]) {
      double m = (log(modfac*Fluxes[1]) - log(modfac*Fluxes[0]))/(log(Energies[1]) - log(Energies[0]));
      double t = log(modfac*Fluxes[1]) - m * log(Energies[1]);
      double logy = m * log(EnergykeV) + t;
      Flux = exp(logy);
    } else if (EnergykeV > Energies.back()) {
      int Last = Energies.size() - 1;
      double m = (log(modfac*Fluxes[Last]) - log(modfac*Fluxes[Last-1]))/(log(Energies[Last]) - log(Energies[Last-1]));
      double t = log(modfac*Fluxes[Last]) - m * log(Energies[Last]);
      double logy = m * log(EnergykeV) + t;
      Flux = exp(logy);
    } else {
      for (unsigned int e = 1; e < Energies.size(); ++e) {
        if (EnergykeV < Energies[e]) {
          double m = (log(modfac*Fluxes[e]) - log(modfac*Fluxes[e-1]))/(log(Energies[e]) - log(Energies[e-1]));
          double t = log(modfac*Fluxes[e]) - m * log(Energies[e]);
          double logy = m * log(EnergykeV) + t;
          Flux = exp(logy);
          break;
        }
      }
    }

    return Flux;
  };



  vector<double> Angle;
  for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
    if (m_AngleBins[b] <= m_HorizonAngle) {
      Angle.push_back(1.0);
    } else {
      Angle.push_back(0.0);
    }
  }

  vector<double> Spectrum;
  for (unsigned int b = 0; b < m_NEnergyBins; ++b) {
    Spectrum.push_back(Interp(m_EnergyBins[b]));
  }


  // Determine the flux in ph/s/cm2 via numerical integration...
  double Flux = 0;
  int Bins = 10000;
  double Min = log(m_EnergyMin);
  double Max = log(m_EnergyMax);
  double Dist = (Max-Min)/Bins;
  for (double e = 0; e <= Bins-1; ++e) {
    double EMin = exp(Min + e*Dist);
    double EMax = exp(Min + (e+1)*Dist);
    double Average = 0.5*(Interp(EMin) + Interp(EMax));
    Flux += Average*(EMax-EMin);
  }
  // Integration factor over sphere
  double AngleFactor = 2*c_Pi * (1 - cos(m_HorizonAngle*c_Rad));

  cout<<"Angle-factor: "<<AngleFactor<<endl;
  cout<<"Flux: "<<Flux<<" particle/cm2/s/sr"<<endl;
  Flux *= AngleFactor;
  cout<<"Flux: "<<Flux<<" particle/cm2/s"<<endl;

  MString Comments;
  Comments += "# Cosmic protons using AMS data\n";
  Comments += "# Assumptions: \n";
  Comments += "# * log-log extrapolation from AMS data\n";
  Comments += "# * Angular distribution is flat out to the Earth-horizon \n";

  WriteEnergyFile("CosmicProtonsAguilar", Spectrum);
  WriteAngleFile("CosmicProtonsAguilar", Angle);
  WriteSourceFile("CosmicProtonsAguilar", Flux, 21, Comments);

  m_HadronicSummaryFiles.push_back("CosmicProtonsAguilar");

  return true;
}


/******************************************************************************
 * Generate a cosmic electron spectrum from AMS data
 * Aguilar et al. 2014
 *   doi:10.1103/PhysRevLett.113.121102
 */
bool BackgroundGenerator::GenerateCosmicElectronsAguilar()
{
  mout<<endl;
  mout<<"Generating a cosmic electrons spectrum using AMS data..."<<endl;

  if (IsValid(1, 1000000000) == false) {
    return false;
  }

  ifstream in;
  MString m_AverageCosmicEPAguilar = "./AguilarElectronPositron.dat";
  in.open(m_AverageCosmicEPAguilar);
  if (in.is_open() == false) {
    mout<<"Error: Unable to open file "<<m_AverageCosmicEPAguilar<<endl;
    return false;
  }

  vector<double> Energies;
  vector<double> Fluxes;

  MString Line;
  bool Start = false;
  while (in.good() == true) {
    Line.ReadLine(in);
    if (Line.Length() < 2) continue;
    if (Line.Contains("EkeV") == true) {
      Start = true;
      continue;
    }
    if (Start == false) continue;
    vector<MString> Tokens = Line.Tokenize(" ");
    if (Tokens.size() != 3) {
      cout<<"Error: The line should have two tokens, rigidity, and differential flux"<<endl;
      cout<<Line<<endl;
      continue;
    }

    double Flux = atof(Tokens[1])/pow(10,10);
    double Energy = atof(Tokens[0]);

    if (Flux > 0) {
      //cout<<Energy<<" "<<Flux<<endl;
      Energies.push_back(Energy);
      Fluxes.push_back(Flux);
    }
  }

  auto Interp = [this, Energies, Fluxes] (double EnergykeV) {

    //solar modulation potential (MV): ~550 for solar minimum ~1100 for solar maximum
    double solmod = 650.;

    double E0 = 0.511/1000;

    // Geomagnetic modulation factor from Mizuno et al. 2004
    double Rigidity = sqrt(EnergykeV*EnergykeV/1000000/1000000 + 2*EnergykeV*E0/1000000);
    double redfac = 1/(1+pow(Rigidity/m_AverageGeomagneticCutOff,-6.0));

    //Solar modulation factor from Gleeson & Axford 1968
    double solmodfac = ((EnergykeV/1000000+E0)*(EnergykeV/1000000+E0)-E0*E0)/(
                (EnergykeV/1000000+E0+solmod/1000)*(EnergykeV/1000000+E0+solmod/1000)-E0*E0);

    double modfac = solmodfac * redfac;

    double Flux = 0.0;
    if (EnergykeV < Energies[0]) {
      double m = (log(modfac*Fluxes[1]) - log(modfac*Fluxes[0]))/(log(Energies[1]) - log(Energies[0]));
      double t = log(modfac*Fluxes[1]) - m * log(Energies[1]);
      double logy = m * log(EnergykeV) + t;
      Flux = exp(logy);
    } else if (EnergykeV > Energies.back()) {
      int Last = Energies.size() - 1;
      double m = (log(modfac*Fluxes[Last]) - log(modfac*Fluxes[Last-1]))/(log(Energies[Last]) - log(Energies[Last-1]));
      double t = log(modfac*Fluxes[Last]) - m * log(Energies[Last]);
      double logy = m * log(EnergykeV) + t;
      Flux = exp(logy);
    } else {
      for (unsigned int e = 1; e < Energies.size(); ++e) {
        if (EnergykeV < Energies[e]) {
          double m = (log(modfac*Fluxes[e]) - log(modfac*Fluxes[e-1]))/(log(Energies[e]) - log(Energies[e-1]));
          double t = log(modfac*Fluxes[e]) - m * log(Energies[e]);
          double logy = m * log(EnergykeV) + t;
          Flux = exp(logy);
          break;
        }
      }
    }

    return Flux;
  };



  vector<double> Angle;
  for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
    if (m_AngleBins[b] <= m_HorizonAngle) {
      Angle.push_back(1.0);
    } else {
      Angle.push_back(0.0);
    }
  }

  vector<double> Spectrum;
  for (unsigned int b = 0; b < m_NEnergyBins; ++b) {
    Spectrum.push_back(Interp(m_EnergyBins[b]));
  }


  // Determine the flux in ph/s/cm2 via numerical integration...
  double Flux = 0;
  int Bins = 10000;
  double Min = log(m_EnergyMin);
  double Max = log(m_EnergyMax);
  double Dist = (Max-Min)/Bins;
  for (double e = 0; e <= Bins-1; ++e) {
    double EMin = exp(Min + e*Dist);
    double EMax = exp(Min + (e+1)*Dist);
    double Average = 0.5*(Interp(EMin) + Interp(EMax));
    Flux += Average*(EMax-EMin);
  }
  // Integration factor over sphere
  double AngleFactor = 2*c_Pi * (1 - cos(m_HorizonAngle*c_Rad));

  cout<<"Angle-factor: "<<AngleFactor<<endl;
  cout<<"Flux: "<<Flux<<" particle/cm2/s/sr"<<endl;
  Flux *= AngleFactor;
  cout<<"Flux: "<<Flux<<" particle/cm2/s"<<endl;

  MString Comments;
  Comments += "# Cosmic electrons using AMS data\n";
  Comments += "# Assumptions: \n";
  Comments += "# * log-log extrapolation from AMS data\n";
  Comments += "# * Angular distribution is flat out to the Earth-horizon \n";

  WriteEnergyFile("CosmicElectronsAguilar", Spectrum);
  WriteAngleFile("CosmicElectronsAguilar", Angle);
  WriteSourceFile("CosmicElectronsAguilar", Flux, 21, Comments);

  m_HadronicSummaryFiles.push_back("CosmicElectronsAguilar");

  return true;
}

/******************************************************************************
 * Generate a cosmic positron spectrum from AMS data
 * Aguilar et al. 2014
 *   doi:10.1103/PhysRevLett.113.121102
 */
bool BackgroundGenerator::GenerateCosmicPositronsAguilar()
{
  mout<<endl;
  mout<<"Generating a cosmic positrons spectrum using AMS data..."<<endl;

  if (IsValid(1, 1000000000) == false) {
    return false;
  }

  ifstream in;
  MString m_AverageCosmicEPAguilar = "./AguilarElectronPositron.dat";
  in.open(m_AverageCosmicEPAguilar);
  if (in.is_open() == false) {
    mout<<"Error: Unable to open file "<<m_AverageCosmicEPAguilar<<endl;
    return false;
  }

  vector<double> Energies;
  vector<double> Fluxes;

  MString Line;
  bool Start = false;
  while (in.good() == true) {
    Line.ReadLine(in);
    if (Line.Length() < 2) continue;
    if (Line.Contains("EkeV") == true) {
      Start = true;
      continue;
    }
    if (Line.Contains("End of File") == true) break;
    if (Start == false) continue;
    vector<MString> Tokens = Line.Tokenize(" ");
    if (Tokens.size() != 3) {
      cout<<"Error: The line should have two tokens, rigidity, and differential flux"<<endl;
      cout<<Line<<endl;
      continue;
    }

    double Flux = atof(Tokens[2])/pow(10,10);
    double Energy = atof(Tokens[0]);

    if (Flux > 0) {
      //cout<<Line<<endl;
      Energies.push_back(Energy);
      Fluxes.push_back(Flux);
    }
  }

  auto Interp = [this, Energies, Fluxes] (double EnergykeV) {

    //solar modulation potential (MV): ~550 for solar minimum ~1100 for solar maximum
    double solmod = 650.;

    double E0 = 0.511/1000;

    // Geomagnetic modulation factor from Mizuno et al. 2004
    double Rigidity = sqrt(EnergykeV*EnergykeV/1000000/1000000 + 2*EnergykeV*E0/1000000);
    double redfac = 1/(1+pow(Rigidity/m_AverageGeomagneticCutOff,-6.0));

    //Solar modulation factor from Gleeson & Axford 1968
    double solmodfac = ((EnergykeV/1000000+E0)*(EnergykeV/1000000+E0)-E0*E0)/(
                (EnergykeV/1000000+E0+solmod/1000)*(EnergykeV/1000000+E0+solmod/1000)-E0*E0);

    double modfac = solmodfac * redfac;

    double Flux = 0.0;
    if (EnergykeV < Energies[0]) {
      double m = (log(modfac*Fluxes[1]) - log(modfac*Fluxes[0]))/(log(Energies[1]) - log(Energies[0]));
      double t = log(modfac*Fluxes[1]) - m * log(Energies[1]);
      double logy = m * log(EnergykeV) + t;
      Flux = exp(logy);
    } else if (EnergykeV > Energies.back()) {
      int Last = Energies.size() - 1;
      double m = (log(modfac*Fluxes[Last]) - log(modfac*Fluxes[Last-1]))/(log(Energies[Last]) - log(Energies[Last-1]));
      double t = log(modfac*Fluxes[Last]) - m * log(Energies[Last]);
      double logy = m * log(EnergykeV) + t;
      Flux = exp(logy);
    } else {
      for (unsigned int e = 1; e < Energies.size(); ++e) {
        if (EnergykeV < Energies[e]) {
          double m = (log(modfac*Fluxes[e]) - log(modfac*Fluxes[e-1]))/(log(Energies[e]) - log(Energies[e-1]));
          double t = log(modfac*Fluxes[e]) - m * log(Energies[e]);
          double logy = m * log(EnergykeV) + t;
          Flux = exp(logy);
          break;
        }
      }
    }

    return Flux;
  };



  vector<double> Angle;
  for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
    if (m_AngleBins[b] <= m_HorizonAngle) {
      Angle.push_back(1.0);
    } else {
      Angle.push_back(0.0);
    }
  }

  vector<double> Spectrum;
  for (unsigned int b = 0; b < m_NEnergyBins; ++b) {
    Spectrum.push_back(Interp(m_EnergyBins[b]));
  }


  // Determine the flux in ph/s/cm2 via numerical integration...
  double Flux = 0;
  int Bins = 10000;
  double Min = log(m_EnergyMin);
  double Max = log(m_EnergyMax);
  double Dist = (Max-Min)/Bins;
  for (double e = 0; e <= Bins-1; ++e) {
    double EMin = exp(Min + e*Dist);
    double EMax = exp(Min + (e+1)*Dist);
    double Average = 0.5*(Interp(EMin) + Interp(EMax));
    Flux += Average*(EMax-EMin);
  }
  // Integration factor over sphere
  double AngleFactor = 2*c_Pi * (1 - cos(m_HorizonAngle*c_Rad));

  cout<<"Angle-factor: "<<AngleFactor<<endl;
  cout<<"Flux: "<<Flux<<" particle/cm2/s/sr"<<endl;
  Flux *= AngleFactor;
  cout<<"Flux: "<<Flux<<" particle/cm2/s"<<endl;

  MString Comments;
  Comments += "# Cosmic positrons using AMS data\n";
  Comments += "# Assumptions: \n";
  Comments += "# * log-log extrapolation from AMS data\n";
  Comments += "# * Angular distribution is flat out to the Earth-horizon \n";

  WriteEnergyFile("GenerateCosmicPositronsAguilar", Spectrum);
  WriteAngleFile("GenerateCosmicPositronsAguilar", Angle);
  WriteSourceFile("GenerateCosmicPositronsAguilar", Flux, 21, Comments);

  m_HadronicSummaryFiles.push_back("GenerateCosmicPositronsAguilar");

  return true;
}

/******************************************************************************
 * Generate an albedo photon spectrum after Ajello et al. 2008 & Mizuno 2004
 * Mizuno is used as the absolute normalization and the Ajello shape is scaled
 * to the Mizuno normalization
 */
bool BackgroundGenerator::GenerateAlbedoPhotonsAjelloMizuno()
{
  auto Ajello = [](double EnergykeV) {
    // Valid up to ~ 1MeV
    return 0.0148/(pow(EnergykeV/33.7, -5) + pow(EnergykeV/33.7, 1.72));
  };

  /* Mizuno downward:
  auto Mizuno = [](double EnergykeV) {
    // Valid above ~ 1 MeV
    if (EnergykeV < 1000000) {
      return (250.0*pow(EnergykeV/1000, -1.7) + 114000*pow(EnergykeV/1000, -2.5)*exp(-pow(EnergykeV/1000/120, -1.5))) / 1000.0 / 10000.0;
    } else if (EnergykeV >= 1000000) {
      return 21500*pow(EnergykeV/1000, -2.2) / 1000.0 / 10000.0;
    }
    return 0.0;
  };
  */

  // Mizuno upward
  auto Mizuno = [](double EnergykeV) {
    // Valid above ~ 1 MeV
    if (EnergykeV < 20000) {
      return 1010.0*pow(EnergykeV/1000, -1.34) / 1000.0 / 10000.0;
    } else if (EnergykeV >= 20000 && EnergykeV < 1000000) {
      return 7290.0*pow(EnergykeV/1000, -2.0) / 1000.0 / 10000.0;
    } else if (EnergykeV >= 1000000) {
      return 29000*pow(EnergykeV/1000, -2.2) / 1000.0 / 10000.0;
    }
    return 0.0;
  };

  auto AjelloMizuno = [Ajello, Mizuno](double EnergykeV, double ScalerAjello, double ScalerMizuno) {
    if (EnergykeV < 1000) {
      return Ajello(EnergykeV)*ScalerAjello;
    } else {
      return Mizuno(EnergykeV)*ScalerMizuno;
    }
  };

  if (m_IsEarthOrbit == false) return true;

  mout<<endl;
  mout<<"Generating an albedo photon spectrum according to Ajello and Mizuno..."<<endl;

  if (IsValid(1, 1000000000) == false) {
    return false;
  }


  vector<double> Angle;
  for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
    if (m_AngleBins[b] >= m_HorizonAngle) {
      Angle.push_back(1.0);
    } else {
      Angle.push_back(0.0);
    }
  }

  // Scaling from Mizuno et al.
  double Rcut_desired = m_AverageGeomagneticCutOff;
  double Rcut_Mizuno = 4.5;
  double ScalerMizuno = pow(Rcut_desired/Rcut_Mizuno, -1.13);

  // Make sure to scale the Ajello result to the Mizuno result at 1 MeV:
  double MizunoValue = Mizuno(1000)*pow(Rcut_desired/Rcut_Mizuno, -1.13);
  double AjelloValue = Ajello(1000);
  double ScalerAjello = MizunoValue/AjelloValue;

  vector<double> Spectrum;
  for (unsigned int b = 0; b < m_NEnergyBins; ++b) {
    Spectrum.push_back(AjelloMizuno(m_EnergyBins[b], ScalerAjello, ScalerMizuno));
    cout<<m_EnergyBins[b]<<":"<<AjelloMizuno(m_EnergyBins[b], ScalerAjello, ScalerMizuno)<<endl;
  }


  // Determine the flux in ph/s/cm2 via numerical integration...
  double Flux = 0;
  int Bins = 10000;
  double Min = log(m_EnergyMin);
  double Max = log(m_EnergyMax);
  double Dist = (Max-Min)/Bins;
  for (double e = 0; e <= Bins-1; ++e) {
    double EMin = exp(Min + e*Dist);
    double EMax = exp(Min + (e+1)*Dist);
    double Average = 0.5*(AjelloMizuno(EMin, ScalerAjello, ScalerMizuno) + AjelloMizuno(EMax, ScalerAjello, ScalerMizuno));
    Flux += Average*(EMax-EMin);
  }
  // Integration factor over sphere
  double AngleFactor = 2*c_Pi * (cos(m_HorizonAngle*c_Rad) - (-1));

  cout<<"Angle-factor: "<<AngleFactor<<endl;
  cout<<"Flux: "<<Flux<<" ph/cm2/s/sr"<<endl;
  Flux *= AngleFactor;
  cout<<"Flux: "<<Flux<<" ph/cm2/s"<<endl;

  MString Comments;
  Comments += "# Albedo photons determinined after Ajello 2008 & Mizuno 2004\n";
  Comments += "# Assumptions: \n";
  Comments += "# (1) Use Ajello below 1 MeV, Mizuno above \n";
  Comments += "# (2) Geomegnetic cutoff difference scaled as in Mizuno \n";
  Comments += "# (3) Ajello normalized to Mizuno \n";
  Comments += "# (4) Angular distribution is flat out to the Earth-horizon (no limb brightening!) \n";

  WriteEnergyFile("AlbedoPhotonsAjelloMizuno", Spectrum);
  WriteAngleFile("AlbedoPhotonsAjelloMizuno", Angle);
  WriteSourceFile("AlbedoPhotonsAjelloMizuno", Flux, 1, Comments);

  m_PhotonicSummaryFiles.push_back("AlbedoPhotonsAjelloMizuno");

  return true;
}



/******************************************************************************
 * Generate an albedo photon spectrum after Tuerler+ (A&A 512, A49, 2010) & Mizuno 2004 & Abdo (PHYSICAL REVIEW D80, 122004, 2009)
 * Mizuno is used as the absolute normalization and the Ajello shape is scaled
 * to the Mizuno normalization
 */
bool BackgroundGenerator::GenerateAlbedoPhotonsTuerlerMizunoAbdo()
{
  // Tuerler low energy component
  auto Tuerler = [](double EnergykeV) {
    if (EnergykeV < 49.4) {
      return 0.05*pow(EnergykeV, -0.37);
    } else {
      return 0.05*pow(49.4, 1.78-0.37)*pow(EnergykeV, -1.78);
    }
  };

  // Mizuno downward component
  auto Mizuno = [](double EnergykeV) {
    // Valid above ~ 1 MeV
    if (EnergykeV < 20000) {
      return 1010.0*pow(EnergykeV/1000, -1.34) / 1000.0 / 10000.0;
    } else if (EnergykeV >= 20000 && EnergykeV < 1000000) {
      return 7290.0*pow(EnergykeV/1000, -2.0) / 1000.0 / 10000.0;
    } else if (EnergykeV >= 1000000) {
      return 29000*pow(EnergykeV/1000, -2.2) / 1000.0 / 10000.0;
    }
    return 0.0;
  };

  // Abdo FERMI
  auto Abdo = [](double EnergykeV) {
    // Valid from 200 MeV & up
    return 1.823e-8*pow(EnergykeV/200000, -2.8);
  };

  auto TuerlerMizunoAbdo = [Tuerler, Mizuno, Abdo](double EnergykeV, double ScalerTuerler, double ScalerMizuno, double ScalerAbdo) {
    if (EnergykeV < 750) {
      return Tuerler(EnergykeV)*ScalerTuerler;
    } else if (EnergykeV > 200000) {
      return Abdo(EnergykeV)*ScalerAbdo;
    } else {
      return Mizuno(EnergykeV)*ScalerMizuno;
    }
  };

  if (m_IsEarthOrbit == false) return true;

  mout<<endl;
  mout<<"Generating an albedo photon spectrum according to Tuerler and Mizuno..."<<endl;

  if (IsValid(1, 1000000000) == false) {
    return false;
  }


  vector<double> Angle;
  for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
    if (m_AngleBins[b] >= m_HorizonAngle) {
      Angle.push_back(1.0);
    } else {
      Angle.push_back(0.0);
    }
  }

  // Scaling from Mizuno et al.
  double Rcut_desired = m_AverageGeomagneticCutOff;
  double Rcut_Mizuno = 4.5;
  double ScalerMizuno = pow(Rcut_desired/Rcut_Mizuno, -1.13);

  // Make sure to scale the Tuerler result to the Mizuno result at 1 MeV:
  double MizunoValue = ScalerMizuno*Mizuno(750);
  double TuerlerValue = Tuerler(750);
  double ScalerTuerler = MizunoValue/TuerlerValue;

  MizunoValue = ScalerMizuno*Mizuno(200000);
  double AbdoValue = Abdo(200000);
  double ScalerAbdo = MizunoValue/AbdoValue;

  vector<double> Spectrum;
  for (unsigned int b = 0; b < m_NEnergyBins; ++b) {
    Spectrum.push_back(TuerlerMizunoAbdo(m_EnergyBins[b], ScalerTuerler, ScalerMizuno, ScalerAbdo));
    cout<<m_EnergyBins[b]<<":"<<TuerlerMizunoAbdo(m_EnergyBins[b], ScalerTuerler, ScalerMizuno, ScalerAbdo)<<endl;
  }


  // Determine the flux in ph/s/cm2 via numerical integration...
  double Flux = 0;
  int Bins = 10000;
  double Min = log(m_EnergyMin);
  double Max = log(m_EnergyMax);
  double Dist = (Max-Min)/Bins;
  for (double e = 0; e <= Bins-1; ++e) {
    double EMin = exp(Min + e*Dist);
    double EMax = exp(Min + (e+1)*Dist);
    double Average = 0.5*(TuerlerMizunoAbdo(EMin, ScalerTuerler, ScalerMizuno, ScalerAbdo) + TuerlerMizunoAbdo(EMax, ScalerTuerler, ScalerMizuno, ScalerAbdo));
    Flux += Average*(EMax-EMin);
  }
  // Integration factor over sphere
  double AngleFactor = 2*c_Pi * (cos(m_HorizonAngle*c_Rad) - (-1));

  cout<<"Angle-factor: "<<AngleFactor<<endl;
  cout<<"Flux: "<<Flux<<" ph/cm2/s/sr"<<endl;
  Flux *= AngleFactor;
  cout<<"Flux: "<<Flux<<" ph/cm2/s"<<endl;

  MString Comments;
  Comments += "# Albedo photons determinined after Turler+ 2010, Mizuno+ 2004, Abdo+ 2009\n";
  Comments += "# Assumptions: \n";
  Comments += "# (1) Use Tuerler below 750 keV, Abdo+ above 200 MeV, Mizuno in between \n";
  Comments += "# (2) Everything is normalized to Mizuno \n";
  Comments += "# (3) Geomagnetic cutoff difference scaled as in Mizuno \n";
  Comments += "# (4) Angular distribution is flat out to the Earth-horizon (no limb brightening!) \n";
  Comments += "# (5) Tuerler includes the reflected CXB component, not just the Albedo component thus is larger at lower energies than Ajello \n";
  Comments += "# (6) The lower energy limit is ~20 keV \n";
  Comments += "# (7) Limb brightening and darkening is NOT included\n";
  Comments += "# (8) There should be a bump form the 511-keV-Compton scatters belog ~450 keV which is NOT included \n";

  WriteEnergyFile("AlbedoPhotonsTuerlerMizunoAbdo", Spectrum);
  WriteAngleFile("AlbedoPhotonsTuerlerMizunoAbdo", Angle);
  WriteSourceFile("AlbedoPhotonsTuerlerMizunoAbdo", Flux, 1, Comments);

  m_PhotonicSummaryFiles.push_back("AlbedoPhotonsTuerlerMizunoAbdo");

  return true;
}



/******************************************************************************
 * Generate a cosmic photon spectrum after Gruber et al. 1999
 */
bool BackgroundGenerator::GenerateAlbedoPhotonLinesHarris()
{
  if (m_IsEarthOrbit == false) return true;

  mout<<endl;
  mout<<"Generating a Cosima spectrum for the 511 keV emission according to Harris..."<<endl;

  if (IsValid(1, 1000000000) == false) {
    return false;
  }

  cout<<"Energy: Mono"<<endl;
  cout<<"Beam: Assume symmetric"<<endl;

  // SMM average orbit: 400-570 km
  double SMMAltitude = 500;
  // SMM fluxes: in 10-3 ph/cm2/s
  // 41.3 ± 0.1 - 45.7 ± 0.1 (< 7 GV)
  // 29.9 ± 0.1 - 31.6 ± 0.1 (7 - 11 GV)
  // 23.3 ± 0.1 - 22.2 ± 0.1 (> 11 GV)
  double SMMFlux = 22.7E-3;

  vector<double> Angle;
  for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
    if (m_AngleBins[b] < m_HorizonAngle) {
      Angle.push_back(1.0);
    } else {
      Angle.push_back(0.0);
    }
  }

  // Scale Flux with altitude
  double Flux = SMMFlux * SMMAltitude*SMMAltitude / (m_Altitude*m_Altitude);

  WriteAngleFile("AnnihilationLineHarris", Angle);


  MString Comments;
  Comments += "# Annihilation line determined after Harris, JGR v.108, 2003\n";
  Comments += "# Assumptions: Average geomagnetic cut-off > 11 GV \n";
  Comments += "#              Average SMM altitude 500 \n";

  WriteSourceFile("AnnihilationLineHarris", Flux, 1, Comments, "Mono 511");

  m_PhotonicSummaryFiles.push_back("AnnihilationLineHarris");

  return true;
}


/******************************************************************************
 * Generate an albedo electron spectrum after Alcaraz et al. 2000 & Mizuno 2004
 */
bool BackgroundGenerator::GenerateAlbedoElectronsAlcarazMizuno()
{
  auto AlcarazMizuno = [] (double EnergyGeV) {
    double F0 = 0.3;
    double a = 2.2;
    double b = 4.0;
    double Ebreak = 3.0;

    double Flux = 0.0;
    if (EnergyGeV < 0.001) return 0.0;
    if (EnergyGeV < Ebreak) {
      if (EnergyGeV < 0.1) {
        Flux = F0 * pow(EnergyGeV/0.1, -1.0); // e/m2/s/sr/MeV
      } else {
        Flux = F0 * pow(EnergyGeV/0.1, -a); // e/m2/s/sr/MeV
      }
    } else {
      Flux = F0 * pow(Ebreak/0.1, -a) * pow(EnergyGeV/Ebreak, -b); // e/m2/s/sr/MeV
    }

    Flux /= (1000*10000); // p/cm2/s/sr/keV
    return Flux;
  };

  if (m_IsEarthOrbit == false) return true;

  mout<<endl;
  mout<<"Generating an albedo electron spectrum according to Alcaraz and Mizuno..."<<endl;

  if (IsValid(1, 1000000000) == false) {
    return false;
  }


  vector<double> Angle;
  for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
    if (m_AngleBins[b] >= m_HorizonAngle) {
      Angle.push_back(1.0);
    } else {
      Angle.push_back(0.0);
    }
  }

  vector<double> Spectrum;
  for (unsigned int b = 0; b < m_NEnergyBins; ++b) {
    Spectrum.push_back(AlcarazMizuno(m_EnergyBins[b]/1000000));
  }


  // Determine the flux in ph/s/cm2 via numerical integration...
  double Flux = 0;
  int Bins = 10000;
  double Min = log(m_EnergyMin);
  double Max = log(m_EnergyMax);
  double Dist = (Max-Min)/Bins;
  for (double e = 0; e <= Bins-1; ++e) {
    double EMin = exp(Min + e*Dist);
    double EMax = exp(Min + (e+1)*Dist);
    double Average = 0.5*(AlcarazMizuno(EMin/1000000) + AlcarazMizuno(EMax/1000000));
    Flux += Average*(EMax-EMin);
  }
  // Integration factor over sphere
  double AngleFactor = 2*c_Pi * (cos(m_HorizonAngle*c_Rad) - (-1));

  cout<<"Angle-factor: "<<AngleFactor<<endl;
  cout<<"Flux: "<<Flux<<" ph/cm2/s/sr"<<endl;
  Flux *= AngleFactor;
  cout<<"Flux: "<<Flux<<" ph/cm2/s"<<endl;

  MString Comments;
  Comments += "# Albedo electrons determinined after Alcaraz 2000 & Mizuno 2004\n";
  Comments += "# Assumptions: \n";
  Comments += "# * Only valid near the equator! \n";
  Comments += "# * Angular distribution is flat out to the Earth-horizon (no limb brightening!) \n";

  WriteEnergyFile("AlbedoElectronsAlcarazMizuno", Spectrum);
  WriteAngleFile("AlbedoElectronsAlcarazMizuno", Angle);
  WriteSourceFile("AlbedoElectronsAlcarazMizuno", Flux, 3, Comments);

  m_LeptonicSummaryFiles.push_back("AlbedoElectronsAlcarazMizuno");

  return true;
}


/******************************************************************************
 * Generate an albedo positron spectrum after Alcaraz et al. 2000 & Mizuno 2004
 */
bool BackgroundGenerator::GenerateAlbedoPositronsAlcarazMizuno()
{
  auto AlcarazMizuno = [] (double EnergyGeV) {
    double F0 = 3.33*0.3; // 3.33 from electron positron difference in equatorial orbit
    double a = 2.2;
    double b = 4.0;
    double Ebreak = 3.0;

    double Flux = 0.0;
    if (EnergyGeV < 0.001) return 0.0;
    if (EnergyGeV < Ebreak) {
      if (EnergyGeV < 0.1) {
        Flux = F0 * pow(EnergyGeV/0.1, -1.0); // p/m2/s/sr/MeV
      } else {
        Flux = F0 * pow(EnergyGeV/0.1, -a); // p/m2/s/sr/MeV
      }
    } else {
      Flux = F0 * pow(Ebreak/0.1, -a) * pow(EnergyGeV/Ebreak, -b); // p/m2/s/sr/MeV
    }

    Flux /= (1000*10000); // p/cm2/s/sr/keV
    return Flux;
  };

  if (m_IsEarthOrbit == false) return true;

  mout<<endl;
  mout<<"Generating an albedo positron spectrum according to Alcaraz and Mizuno..."<<endl;

  if (IsValid(1, 1000000000) == false) {
    return false;
  }


  vector<double> Angle;
  for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
    if (m_AngleBins[b] >= m_HorizonAngle) {
      Angle.push_back(1.0);
    } else {
      Angle.push_back(0.0);
    }
  }

  vector<double> Spectrum;
  for (unsigned int b = 0; b < m_NEnergyBins; ++b) {
    Spectrum.push_back(AlcarazMizuno(m_EnergyBins[b]/1000000));
  }


  // Determine the flux in ph/s/cm2 via numerical integration...
  double Flux = 0;
  int Bins = 10000;
  double Min = log(m_EnergyMin);
  double Max = log(m_EnergyMax);
  double Dist = (Max-Min)/Bins;
  for (double e = 0; e <= Bins-1; ++e) {
    double EMin = exp(Min + e*Dist);
    double EMax = exp(Min + (e+1)*Dist);
    double Average = 0.5*(AlcarazMizuno(EMin/1000000) + AlcarazMizuno(EMax/1000000));
    Flux += Average*(EMax-EMin);
  }
  // Integration factor over sphere
  double AngleFactor = 2*c_Pi * (cos(m_HorizonAngle*c_Rad) - (-1));

  cout<<"Angle-factor: "<<AngleFactor<<endl;
  cout<<"Flux: "<<Flux<<" ph/cm2/s/sr"<<endl;
  Flux *= AngleFactor;
  cout<<"Flux: "<<Flux<<" ph/cm2/s"<<endl;

  MString Comments;
  Comments += "# Albedo positrons determinined after Alcaraz 2000 & Mizuno 2004\n";
  Comments += "# Assumptions: \n";
  Comments += "# * Don't remember... \n";
  Comments += "# * Angular distribution is flat out to the Earth-horizon (no limb brightening!) \n";

  WriteEnergyFile("AlbedoPositronsAlcarazMizuno", Spectrum);
  WriteAngleFile("AlbedoPositronsAlcarazMizuno", Angle);
  WriteSourceFile("AlbedoPositronsAlcarazMizuno", Flux, 2, Comments);

  m_LeptonicSummaryFiles.push_back("AlbedoPositronsAlcarazMizuno");

  return true;
}


/******************************************************************************
 * Generate an albedo proton spectrum after Alcaraz et al. 2000
 */
bool BackgroundGenerator::GenerateAlbedoProtonsAlcaraz()
{
  auto AlcarazDownward = [] (double EnergykeV) {
    double EnergyGeV = 0.000001*EnergykeV;

    vector<double> Energies;
    vector<double> Fluxes;

    // Extrapolated to lower energies
    Energies.push_back(0.001); Fluxes.push_back(1.07);
    Energies.push_back(0.002); Fluxes.push_back(0.7819);
    Energies.push_back(0.005); Fluxes.push_back(0.51533);
    Energies.push_back(0.01); Fluxes.push_back(0.3759);
    Energies.push_back(0.02); Fluxes.push_back(0.2742);
    Energies.push_back(0.05); Fluxes.push_back(0.180722);
    // Original measurements near equator
    Energies.push_back(0.07); Fluxes.push_back(1.67E-01);
    Energies.push_back(0.10); Fluxes.push_back(1.21E-01);
    Energies.push_back(0.15); Fluxes.push_back(9.79E-02);
    Energies.push_back(0.22); Fluxes.push_back(8.62E-02);
    Energies.push_back(0.31); Fluxes.push_back(7.01E-02);
    Energies.push_back(0.44); Fluxes.push_back(5.04E-02);
    Energies.push_back(0.62); Fluxes.push_back(3.28E-02);
    Energies.push_back(0.85); Fluxes.push_back(2.06E-02);
    Energies.push_back(1.15); Fluxes.push_back(1.16E-02);
    Energies.push_back(1.54); Fluxes.push_back(6.69E-03);
    Energies.push_back(2.02); Fluxes.push_back(2.86E-03);
    Energies.push_back(2.62); Fluxes.push_back(1.10E-03);
    Energies.push_back(3.38); Fluxes.push_back(4.43E-04);
    Energies.push_back(4.31); Fluxes.push_back(1.57E-04);
    Energies.push_back(5.45); Fluxes.push_back(6.10E-05);
    // Extrapolated to higher energies
    Energies.push_back(6.86); Fluxes.push_back(1.0E-05);
    Energies.push_back(8.60); Fluxes.push_back(2.0E-06);
    Energies.push_back(10.73); Fluxes.push_back(1.0E-07);
    Energies.push_back(13.34);

    for (unsigned int e = 0; e < Energies.size() - 1; ++e) {
      Energies[e] = 0.5*(Energies[e] + Energies[e+1]);
    }
    Energies.resize(Fluxes.size());

    if (EnergyGeV < Energies[0] || EnergyGeV >= Energies.back()) return 0.0;

    double Flux = 0.0;
    for (unsigned int e = 1; e < Energies.size()-1; ++e) {
      if (EnergyGeV < Energies[e]) {
        // Interpolate - double log:
        double m = (log(Fluxes[e]) - log(Fluxes[e-1]))/(log(Energies[e]) - log(Energies[e-1]));
        double t = log(Fluxes[e]) - m * log(Energies[e]);
        double logy = m * log(EnergyGeV) + t;
        Flux = exp(logy);
        break;
      }
    }

    return Flux/1000.0/10000.0;
  };

  auto AlcarazUpward = [] (double EnergykeV) {
    double EnergyGeV = 0.000001*EnergykeV;

    vector<double> Energies;
    vector<double> Fluxes;

    // Extrapolated up
    Energies.push_back(0.001); Fluxes.push_back(2.52);
    Energies.push_back(0.002); Fluxes.push_back(1.62);
    Energies.push_back(0.005); Fluxes.push_back(0.91);
    Energies.push_back(0.01); Fluxes.push_back(0.59);
    Energies.push_back(0.02); Fluxes.push_back(0.38);
    Energies.push_back(0.05); Fluxes.push_back(0.21);
    // Original
    Energies.push_back(0.07); Fluxes.push_back(16.4E-2);
    Energies.push_back(0.10); Fluxes.push_back(10.9E-2);
    Energies.push_back(0.15); Fluxes.push_back(85.3E-3);
    Energies.push_back(0.22); Fluxes.push_back(84.8E-3);
    Energies.push_back(0.31); Fluxes.push_back(66.8E-3);
    Energies.push_back(0.44); Fluxes.push_back(48.4E-3);
    Energies.push_back(0.62); Fluxes.push_back(32.7E-3);
    Energies.push_back(0.85); Fluxes.push_back(20.2E-3);
    Energies.push_back(1.15); Fluxes.push_back(12.4E-3);
    Energies.push_back(1.54); Fluxes.push_back(62.0E-4);
    Energies.push_back(2.02); Fluxes.push_back(25.9E-4);
    Energies.push_back(2.62); Fluxes.push_back(10.7E-4);
    Energies.push_back(3.38); Fluxes.push_back(29.7E-5);
    Energies.push_back(4.31); Fluxes.push_back(11.2E-5);
    Energies.push_back(5.45); Fluxes.push_back(37.0E-6);
    // Extrapolated down
    Energies.push_back(6.86); Fluxes.push_back(1.0E-05);
    Energies.push_back(8.60); Fluxes.push_back(2.0E-06);
    Energies.push_back(10.73); Fluxes.push_back(1.0E-07);
    Energies.push_back(13.34);

    for (unsigned int e = 0; e < Energies.size() - 1; ++e) {
      Energies[e] = 0.5*(Energies[e] + Energies[e+1]);
    }
    Energies.resize(Fluxes.size());

    if (EnergyGeV < Energies[0] || EnergyGeV >= Energies.back()) return 0.0;

    double Flux = 0.0;
    for (unsigned int e = 1; e < Energies.size()-1; ++e) {
      if (EnergyGeV < Energies[e]) {
        // Interpolate - double log:
        double m = (log(Fluxes[e]) - log(Fluxes[e-1]))/(log(Energies[e]) - log(Energies[e-1]));
        double t = log(Fluxes[e]) - m * log(Energies[e]);
        double logy = m * log(EnergyGeV) + t;
        Flux = exp(logy);
        break;
      }
    }

    return Flux/1000.0/10000.0;
  };

  if (m_IsEarthOrbit == false) return true;

  mout<<endl;
  mout<<"Generating an albedo proton spectrum according to Alcaraz..."<<endl;

  if (IsValid(1, 1000000000) == false) {
    return false;
  }


  vector<double> Angle;
  vector<double> Spectrum;
  double Flux = 0;
  int Bins = 10000;
  double AngleFactor = 2*c_Pi * (cos(m_HorizonAngle*c_Rad) - (-1));
  MString Comments;


  // Downward component

  Angle.clear();
  for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
    if (m_AngleBins[b] <= 90) {
      Angle.push_back(1.0);
    } else {
      Angle.push_back(0.0);
    }
  }

  Spectrum.clear();
  for (unsigned int b = 0; b < m_NEnergyBins; ++b) {
    Spectrum.push_back(AlcarazDownward(m_EnergyBins[b]));
  }


  // Determine the flux in ph/s/cm2 via numerical integration...
  Flux = 0;
  Bins = 10000;
  double Min = log(m_EnergyMin);
  double Max = log(m_EnergyMax);
  double Dist = (Max-Min)/Bins;
  for (double e = 0; e <= Bins-1; ++e) {
    double EMin = exp(Min + e*Dist);
    double EMax = exp(Min + (e+1)*Dist);
    double Average = 0.5*(AlcarazDownward(EMin) + AlcarazDownward(EMax));
    Flux += Average*(EMax-EMin);
  }
  // Integration factor over half sphere
  AngleFactor = 2*c_Pi;

  cout<<"Angle-factor: "<<AngleFactor<<endl;
  cout<<"Flux: "<<Flux<<" ph/cm2/s/sr"<<endl;
  Flux *= AngleFactor;
  cout<<"Flux: "<<Flux<<" ph/cm2/s"<<endl;

  Comments = "";
  Comments += "# Albedo proton spectrum downward component determinined after Alcaraz 2000\n";
  Comments += "# Assumptions: \n";
  Comments += "# * Equatorial low-earth orbit \n";
  Comments += "# * Angular distribution is flat \n";
  Comments += "# * Everything below 70 MeV is a wild guess extrapolation \n";

  WriteEnergyFile("AlbedoProtonAlcarazDownward", Spectrum);
  WriteAngleFile("AlbedoProtonAlcarazDownward", Angle);
  WriteSourceFile("AlbedoProtonAlcarazDownward", Flux, 4, Comments);

  m_HadronicSummaryFiles.push_back("AlbedoProtonAlcarazDownward");


  // Upward component

  Angle.clear();
  for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
    if (m_AngleBins[b] >= 90) {

      Angle.push_back(1.0);
    } else {
      Angle.push_back(0.0);
    }
  }

  Spectrum.clear();
  for (unsigned int b = 0; b < m_NEnergyBins; ++b) {
    Spectrum.push_back(AlcarazUpward(m_EnergyBins[b]));
  }


  // Determine the flux in ph/s/cm2 via numerical integration...
  Flux = 0;
  Bins = 10000;
  Min = log(m_EnergyMin);
  Max = log(m_EnergyMax);
  Dist = (Max-Min)/Bins;
  for (double e = 0; e <= Bins-1; ++e) {
    double EMin = exp(Min + e*Dist);
    double EMax = exp(Min + (e+1)*Dist);
    double Average = 0.5*(AlcarazUpward(EMin) + AlcarazUpward(EMax));
    Flux += Average*(EMax-EMin);
  }
  // Integration factor over half sphere
  AngleFactor = 2*c_Pi;

  cout<<"Angle-factor: "<<AngleFactor<<endl;
  cout<<"Flux: "<<Flux<<" ph/cm2/s/sr"<<endl;
  Flux *= AngleFactor;
  cout<<"Flux: "<<Flux<<" ph/cm2/s"<<endl;

  Comments = "";
  Comments += "# Albedo proton spectrum upward component determinined after Alcaraz 2000\n";
  Comments += "# Assumptions: \n";
  Comments += "# * Equatorial low-earth orbit \n";
  Comments += "# * Angular distribution is flat out to the Earth-horizon \n";

  WriteEnergyFile("AlbedoProtonAlcarazUpward", Spectrum);
  WriteAngleFile("AlbedoProtonAlcarazUpward", Angle);
  WriteSourceFile("AlbedoProtonAlcarazUpward", Flux, 4, Comments);

  m_HadronicSummaryFiles.push_back("AlbedoProtonAlcarazUpward");

  return true;
}


/******************************************************************************
 * Generate an albedo proton spectrum after Alcaraz et al. 2000
 * Combines and averages up and downward component
 */
bool BackgroundGenerator::GenerateAlbedoProtonsAveragedAlcaraz()
{
  auto AlcarazDownward = [] (double EnergykeV) {
    double EnergyGeV = 0.000001*EnergykeV;

    vector<double> Energies;
    vector<double> Fluxes;

    // Extrapolated to lower energies from last two values
    Energies.push_back(0.001); Fluxes.push_back(7.753);
    Energies.push_back(0.002); Fluxes.push_back(4.15);
    Energies.push_back(0.005); Fluxes.push_back(1.812);
    Energies.push_back(0.010); Fluxes.push_back(0.969);
    Energies.push_back(0.020); Fluxes.push_back(0.518);
    Energies.push_back(0.050); Fluxes.push_back(0.226);
    // Original measurements near equator, units: p/m^2/s/sr/MeV
    Energies.push_back(0.07); Fluxes.push_back(1.67E-01);
    Energies.push_back(0.10); Fluxes.push_back(1.21E-01);
    Energies.push_back(0.15); Fluxes.push_back(9.79E-02);
    Energies.push_back(0.22); Fluxes.push_back(8.62E-02);
    Energies.push_back(0.31); Fluxes.push_back(7.01E-02);
    Energies.push_back(0.44); Fluxes.push_back(5.04E-02);
    Energies.push_back(0.62); Fluxes.push_back(3.28E-02);
    Energies.push_back(0.85); Fluxes.push_back(2.06E-02);
    Energies.push_back(1.15); Fluxes.push_back(1.16E-02);
    Energies.push_back(1.54); Fluxes.push_back(6.69E-03);
    Energies.push_back(2.02); Fluxes.push_back(2.86E-03);
    Energies.push_back(2.62); Fluxes.push_back(1.10E-03);
    Energies.push_back(3.38); Fluxes.push_back(4.43E-04);
    Energies.push_back(4.31); Fluxes.push_back(1.57E-04);
    Energies.push_back(5.45); Fluxes.push_back(6.10E-05);
    // Extrapolated to higher energies
    Energies.push_back(6.86); Fluxes.push_back(1.0E-05);
    Energies.push_back(8.60); Fluxes.push_back(2.0E-06);
    Energies.push_back(10.73); Fluxes.push_back(1.0E-07);
    Energies.push_back(13.34);

    // Re-extrapolate down:
    double Gradient = (log10(Fluxes[7]) - log10(Fluxes[6])) / (log10(Energies[7]) - log10(Energies[6]));
    for (unsigned int i = 5; i <= 5; --i) {
      Fluxes[i] = pow(10, log10(Fluxes[i+1]) - Gradient*(log10(Energies[i+1]) - log10(Energies[i])));
    }

    for (unsigned int e = 0; e < Energies.size() - 1; ++e) {
      Energies[e] = 0.5*(Energies[e] + Energies[e+1]);
    }
    Energies.resize(Fluxes.size());

    if (EnergyGeV < Energies[0] || EnergyGeV >= Energies.back()) return 0.0;

    double Flux = 0.0;
    for (unsigned int e = 1; e < Energies.size()-1; ++e) {
      if (EnergyGeV < Energies[e]) {
        // Interpolate - double log:
        double m = (log(Fluxes[e]) - log(Fluxes[e-1]))/(log(Energies[e]) - log(Energies[e-1]));
        double t = log(Fluxes[e]) - m * log(Energies[e]);
        double logy = m * log(EnergyGeV) + t;
        Flux = exp(logy);
        break;
      }
    }

    return Flux/1000.0/10000.0; // return Flux in p/cm2/s/sr/keV
  };

  auto AlcarazUpward = [] (double EnergykeV) {
    double EnergyGeV = 0.000001*EnergykeV;

    vector<double> Energies;
    vector<double> Fluxes;

    // Extrapolated up
    Energies.push_back(0.001); Fluxes.push_back(21.29);
    Energies.push_back(0.002); Fluxes.push_back(9.62);
    Energies.push_back(0.005); Fluxes.push_back(3.369);
    Energies.push_back(0.01); Fluxes.push_back(1.523);
    Energies.push_back(0.02); Fluxes.push_back(0.689);
    Energies.push_back(0.05); Fluxes.push_back(0.241);
    // Original
    Energies.push_back(0.07); Fluxes.push_back(16.4E-2);
    Energies.push_back(0.10); Fluxes.push_back(10.9E-2);
    Energies.push_back(0.15); Fluxes.push_back(85.3E-3);
    Energies.push_back(0.22); Fluxes.push_back(84.8E-3);
    Energies.push_back(0.31); Fluxes.push_back(66.8E-3);
    Energies.push_back(0.44); Fluxes.push_back(48.4E-3);
    Energies.push_back(0.62); Fluxes.push_back(32.7E-3);
    Energies.push_back(0.85); Fluxes.push_back(20.2E-3);
    Energies.push_back(1.15); Fluxes.push_back(12.4E-3);
    Energies.push_back(1.54); Fluxes.push_back(62.0E-4);
    Energies.push_back(2.02); Fluxes.push_back(25.9E-4);
    Energies.push_back(2.62); Fluxes.push_back(10.7E-4);
    Energies.push_back(3.38); Fluxes.push_back(29.7E-5);
    Energies.push_back(4.31); Fluxes.push_back(11.2E-5);
    Energies.push_back(5.45); Fluxes.push_back(37.0E-6);
    // Extrapolated down
    Energies.push_back(6.86); Fluxes.push_back(1.0E-05);
    Energies.push_back(8.60); Fluxes.push_back(2.0E-06);
    Energies.push_back(10.73); Fluxes.push_back(1.0E-07);
    Energies.push_back(13.34);

    // Re-extrapolate down:
    double Gradient = (log10(Fluxes[7]) - log10(Fluxes[6])) / (log10(Energies[7]) - log10(Energies[6]));
    for (unsigned int i = 5; i <= 5; --i) {
      Fluxes[i] = pow(10, log10(Fluxes[i+1]) - Gradient*(log10(Energies[i+1]) - log10(Energies[i])));
    }

    for (unsigned int e = 0; e < Energies.size() - 1; ++e) {
      Energies[e] = 0.5*(Energies[e] + Energies[e+1]);
    }
    Energies.resize(Fluxes.size());

    if (EnergyGeV < Energies[0] || EnergyGeV >= Energies.back()) return 0.0;

    double Flux = 0.0;
    for (unsigned int e = 1; e < Energies.size()-1; ++e) {
      if (EnergyGeV < Energies[e]) {
        // Interpolate - double log:
        double m = (log(Fluxes[e]) - log(Fluxes[e-1]))/(log(Energies[e]) - log(Energies[e-1]));
        double t = log(Fluxes[e]) - m * log(Energies[e]);
        double logy = m * log(EnergyGeV) + t;
        Flux = exp(logy);
        break;
      }
    }

    return Flux/1000.0/10000.0; // return Flux in p/cm2/s/sr/keV
  };

  if (m_IsEarthOrbit == false) return true;

  mout<<endl;
  mout<<"Generating an albedo proton spectrum according to Alcaraz..."<<endl;

  if (IsValid(1, 1000000000) == false) {
    return false;
  }


  vector<double> Angle;
  vector<double> Spectrum;
  double Flux = 0;
  int Bins = 10000;


  // Upward and Downward component combined

  Angle.clear();
  for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
    Angle.push_back(1.0);
  }

  Spectrum.clear();
  for (unsigned int b = 0; b < m_NEnergyBins; ++b) {
    Spectrum.push_back(0.5*(AlcarazUpward(m_EnergyBins[b]) + AlcarazDownward(m_EnergyBins[b])) );
  }


  // Determine the flux in ph/s/cm2/sr via numerical integration...
  Flux = 0;
  Bins = 10000;
  double Min = log(m_EnergyMin);
  double Max = log(m_EnergyMax);
  double Dist = (Max-Min)/Bins;
  for (double e = 0; e <= Bins-1; ++e) {
    double EMin = exp(Min + e*Dist);
    double EMax = exp(Min + (e+1)*Dist);
    double Average = 0.5*(0.5*(AlcarazUpward(EMin) + AlcarazDownward(EMin)) + 0.5*(AlcarazUpward(EMax) + AlcarazDownward(EMax)));
    Flux += Average*(EMax-EMin);
  }

  // Integration factor over full sphere
  double AngleFactor = 4*c_Pi;

  cout<<"Angle-factor: "<<AngleFactor<<endl;
  cout<<"Flux: "<<Flux<<" ph/cm2/s/sr"<<endl;
  Flux *= AngleFactor;
  cout<<"Flux: "<<Flux<<" ph/cm2/s"<<endl;

  MString Comments;
  Comments = "";
  Comments += "# Albedo proton spectrum upward and downward component combined and averaged determinined after Alcaraz 2000\n";
  Comments += "# Assumptions: \n";
  Comments += "# * Equatorial low-earth orbit \n";
  Comments += "# * Angular distribution is isotropic \n";
  Comments += "# * Alcaraz's data goes from 70 MeV to 5 GeV. Below that we extrapolate from gradient from last data point.\n";
  Comments += "# * We arbitrarily cut a 1 MeV .\n";

  WriteEnergyFile("AlbedoProtonAlcaraz", Spectrum);
  WriteAngleFile("AlbedoProtonAlcaraz", Angle);
  WriteSourceFile("AlbedoProtonAlcaraz", Flux, 4, Comments);

  m_HadronicSummaryFiles.push_back("AlbedoProtonAlcaraz");

  return true;
}


/******************************************************************************
 * Generate an albedo proton spectrum after Mizuno+ 2004
 * Combines and averages up and downward component
 */
bool BackgroundGenerator::GenerateAlbedoProtonsAveragedMizuno()
{
  // 0.0 ≤ θM ≤ 0.2 down/upward cutoff PL 0.136/0.123/0.155/0.51  F0 /F1 /a/Ec (GeV) f

  auto MizunoAveraged = [] (double EnergykeV) {
    double EnergyMeV = 0.001*EnergykeV;

    if (EnergyMeV < 1) return 0.0;

    double Flux = 0.0;
    if (EnergyMeV > 100) {
      Flux = 0.123*pow(EnergyMeV/1000, -0.155) * exp(-pow(EnergyMeV/511, -0.155+1)); // p/m2/s/sr/MeV
    } else {
      Flux = 0.136*pow(EnergyMeV/100, -1);
    }

    return 2*Flux/1000.0/10000.0; // return Flux in p/cm2/s/sr/keV Factor 2 to consider upward and downward component
  };

  if (m_IsEarthOrbit == false) return true;

  mout<<endl;
  mout<<"Generating an albedo proton spectrum according to Mizuno..."<<endl;

  if (IsValid(1, 1000000000) == false) {
    return false;
  }


  vector<double> Angle;
  vector<double> Spectrum;
  double Flux = 0;
  int Bins = 10000;


  // Upward and Downward component combined

  Angle.clear();
  for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
    Angle.push_back(1.0);
  }

  Spectrum.clear();
  for (unsigned int b = 0; b < m_NEnergyBins; ++b) {
    Spectrum.push_back(MizunoAveraged(m_EnergyBins[b]));
  }


  // Determine the flux in ph/s/cm2/sr via numerical integration...
  Flux = 0;
  Bins = 10000;
  double Min = log(m_EnergyMin);
  double Max = log(m_EnergyMax);
  double Dist = (Max-Min)/Bins;
  for (double e = 0; e <= Bins-1; ++e) {
    double EMin = exp(Min + e*Dist);
    double EMax = exp(Min + (e+1)*Dist);
    double Average = 0.5*(MizunoAveraged(EMin) + MizunoAveraged(EMax));
    Flux += Average*(EMax-EMin);
  }
  cout<<"Flux: "<<Flux<<" ph/cm2/s/sr"<<endl;

  // Integration factor over full sphere
  double AngleFactor = 4*c_Pi;
  cout<<"Angle-factor: "<<AngleFactor<<endl;
  Flux *= AngleFactor;
  cout<<"Flux: "<<Flux<<" ph/cm2/s"<<endl;

  MString Comments;
  Comments = "";
  Comments += "# Albedo proton spectrum upward and downward component combined and averaged determinined after Mizuno 2004\n";
  Comments += "# Assumptions: \n";
  Comments += "# * Equatorial low-earth orbit \n";
  Comments += "# * Angular distribution is isotropic \n";

  WriteEnergyFile("AlbedoProtonMizuno", Spectrum);
  WriteAngleFile("AlbedoProtonMizuno", Angle);
  WriteSourceFile("AlbedoProtonMizuno", Flux, 4, Comments);

  m_HadronicSummaryFiles.push_back("AlbedoProtonMizuno");

  return true;
}



/******************************************************************************
 * Generate an albedo neutron spectrum after Morris and Kole
 */
bool BackgroundGenerator::GenerateAlbedoNeutronsMorrisKole()
{
  auto Morris = [this] (double EnergyMeV) {
    double Angle = 0;
    double M = 1.0;
    double b = 0.00255;
    double alpha = 0.152;

    double Flux = 0.0;
    if (EnergyMeV < 70) {
      Flux = 0.036*M*(1-b*Angle)*exp(-alpha*m_AverageGeomagneticCutOff)*pow(EnergyMeV, -0.6);
    } else {
      Flux = 8.7*M*(1-b*Angle)*exp(-alpha*m_AverageGeomagneticCutOff)*pow(EnergyMeV, -1.89);
    }

    return Flux / 1000.0; // Switch from n/MeV/cm2/s to n/keV/cm2/s
  };

  auto Kole = [this] (double EnergyMeV) {
    double Pressure = 0;
    double MagLat = m_Inclination; // Orig: 10
    double SolarActivity = 0.5; // Solar Minimum (0), Solar Maximum (1)

    double a = 0.0003 + (7.0-5.0*SolarActivity)*0.001*(1-tanh(c_Rad*(180-4.0*MagLat)));
    double b = 0.0140 + (1.4-0.9*SolarActivity)*  0.1*(1-tanh(c_Rad*(180-3.5*MagLat)));
    double c = 180 -                             42*(1-tanh(c_Rad*(180-5.5*MagLat)));
    double d = -0.008 + (6.0-1.0*SolarActivity)*0.001*(1-tanh(c_Rad*(180-4.4*MagLat)));

    double Slope1 = -0.29 * exp(-Pressure/7.5) + 0.735;
    double Norm1 = (a*Pressure + b)*exp(-Pressure/c) + d;
    double Slope2 = -0.247 * exp(-Pressure/36.5) + 1.4;
    double Norm2 = Norm1*pow(0.9, -Slope1+Slope2);
    double Slope3 = -0.40 * exp(-Pressure/40.0) + 0.9;
    double Norm3 = Norm2*pow(15, -Slope2+Slope3);
    double Slope4 = -0.46 * exp(-Pressure/100.0) + 2.53;
    double Norm4 = Norm3*pow(70, -Slope3+Slope4);

    double Flux = 0.0;
    if (EnergyMeV < 0.9) {
      Flux = Norm1 * pow(EnergyMeV, -Slope1);
    } else if (EnergyMeV >= 0.9 && EnergyMeV < 15) {
      Flux = Norm2 * pow(EnergyMeV, -Slope2);
    } else if (EnergyMeV >= 15 && EnergyMeV < 70) {
      Flux = Norm3 * pow(EnergyMeV, -Slope3);
    } else if (EnergyMeV >= 70) {
      Flux = Norm4 * pow(EnergyMeV, -Slope4);
    }

    return Flux / 1000.0; // Switch from n/MeV/cm2/s to n/keV/cm2/s
  };

  if (m_IsEarthOrbit == false) return true;

  mout<<endl;
  mout<<"Generating an albedo neutron spectrum according to Morris and Kole..."<<endl;

  if (IsValid(1, 1000000000) == false) {
    return false;
  }


  vector<double> Angle;
  for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
    if (m_AngleBins[b] >= m_HorizonAngle) {
      Angle.push_back(1.0);
    } else {
      Angle.push_back(0.0);
    }
  }
  // Integration factor over sphere
  double AngleFactor = 2*c_Pi * (cos(m_HorizonAngle*c_Rad) - (-1));
  cout<<"Angle-factor: "<<AngleFactor<<endl;

  cout<<"Scaler: "<<Morris(150)/Kole(150)<<endl;
  double Scaler = 1.0; //Morris(150)/Kole(150);
  cout<<"Normalization at 30 MeV: Morris: "<<Morris(30)<<" vs. Kole "<<Kole(30)<<endl;
  cout<<"Normalization at 150 MeV: Morris: "<<Morris(150)<<" vs. Kole "<<Kole(150)<<endl;

  vector<double> Spectrum;
  for (unsigned int b = 0; b < m_NEnergyBins; ++b) {
    Spectrum.push_back(Kole(m_EnergyBins[b]/1000)*Scaler/AngleFactor); // We need the angle factor here to switch from ph/cm2/s/keV to ph/cm2/s/keV/sr
  }


  // Determine the flux in ph/s/cm2 via numerical integration...
  double Flux = 0;
  int Bins = 10000;
  double Min = log(m_EnergyMin);
  double Max = log(m_EnergyMax);
  double Dist = (Max-Min)/Bins;
  for (double e = 0; e <= Bins-1; ++e) {
    double EMin = exp(Min + e*Dist);
    double EMax = exp(Min + (e+1)*Dist);
    double Average = 0.5*(Kole(EMin/1000) + Kole(EMax/1000))*Scaler;
    Flux += Average*(EMax-EMin);
  }

  cout<<"Flux: "<<Flux<<" n/cm2/s/sr"<<endl;
  //Flux *= AngleFactor;
  //cout<<"Flux: "<<Flux<<" n/cm2/s"<<endl;

  MString Comments;
  Comments += "# Albedo neutrons determinined after Morris 1995 & Kole 2014\n";
  Comments += "# Assumptions: \n";
  Comments += "# * Use the Morris absolute flux value @ 150 MeV and the shape of Kole \n";
  Comments += "# * Angular distribution is flat out to the Earth-horizon (no limb brightening!) \n";
  Comments += "# * Solar activity is assumed to be half way between solar minimum and maximum \n";
  Comments += "# * The inclination was used to approximate the average Magnetic Latitude \n";

  WriteEnergyFile("AlbedoNeutronsMorrisKole", Spectrum);
  WriteAngleFile("AlbedoNeutronsMorrisKole", Angle);
  WriteSourceFile("AlbedoNeutronsMorrisKole", Flux, 6, Comments);

  m_HadronicSummaryFiles.push_back("AlbedoNeutronsMorrisKole");

  return true;
}


/******************************************************************************
 * Generate an albedo neutron spectrum after Kole
 */
bool BackgroundGenerator::GenerateAlbedoNeutronsKole()
{

  auto Kole = [this] (double EnergyMeV) {
    double Pressure = 0;
    double MagLat = m_Inclination; // Orig: 10
    double SolarActivity = 0.5; // Solar Minimum (0), Solar Maximum (1)

    double a = 0.0003 + (7.0-5.0*SolarActivity)*0.001*(1-tanh(c_Rad*(180-4.0*MagLat)));
    double b = 0.0140 + (1.4-0.9*SolarActivity)*  0.1*(1-tanh(c_Rad*(180-3.5*MagLat)));
    double c = 180 -                             42*(1-tanh(c_Rad*(180-5.5*MagLat)));
    double d = -0.008 + (6.0-1.0*SolarActivity)*0.001*(1-tanh(c_Rad*(180-4.4*MagLat)));

    double Slope1 = -0.29 * exp(-Pressure/7.5) + 0.735;
    double Norm1 = (a*Pressure + b)*exp(-Pressure/c) + d;
    double Slope2 = -0.247 * exp(-Pressure/36.5) + 1.4;
    double Norm2 = Norm1*pow(0.9, -Slope1+Slope2);
    double Slope3 = -0.40 * exp(-Pressure/40.0) + 0.9;
    double Norm3 = Norm2*pow(15, -Slope2+Slope3);
    double Slope4 = -0.46 * exp(-Pressure/100.0) + 2.53;
    double Norm4 = Norm3*pow(70, -Slope3+Slope4);

    double Flux = 0.0;
    if (EnergyMeV < 0.9) {
      Flux = Norm1 * pow(EnergyMeV, -Slope1);
    } else if (EnergyMeV >= 0.9 && EnergyMeV < 15) {
      Flux = Norm2 * pow(EnergyMeV, -Slope2);
    } else if (EnergyMeV >= 15 && EnergyMeV < 70) {
      Flux = Norm3 * pow(EnergyMeV, -Slope3);
    } else if (EnergyMeV >= 70) {
      Flux = Norm4 * pow(EnergyMeV, -Slope4);
    }

    return Flux / 1000.0; // Switch from n/MeV/cm2/s to n/keV/cm2/s
  };

  if (m_IsEarthOrbit == false) return true;

  mout<<endl;
  mout<<"Generating an albedo neutron spectrum according to Morris and Kole..."<<endl;

  if (IsValid(1, 1000000000) == false) {
    return false;
  }


  vector<double> Angle;
  for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
    if (m_AngleBins[b] >= m_HorizonAngle) {
      Angle.push_back(1.0);
    } else {
      Angle.push_back(0.0);
    }
  }
  // Integration factor over sphere
  double AngleFactor = 2*c_Pi * (cos(m_HorizonAngle*c_Rad) - (-1));

  vector<double> Spectrum;
  for (unsigned int b = 0; b < m_NEnergyBins; ++b) {
    Spectrum.push_back(Kole(m_EnergyBins[b]/1000)/AngleFactor); // We need the angle factor here to switch from ph/cm2/s/keV to ph/cm2/s/keV/sr
  }


  // Determine the flux in ph/s/cm2 via numerical integration...
  double Flux = 0;
  int Bins = 10000;
  double Min = log(m_EnergyMin);
  double Max = log(m_EnergyMax);
  double Dist = (Max-Min)/Bins;
  for (double e = 0; e <= Bins-1; ++e) {
    double EMin = exp(Min + e*Dist);
    double EMax = exp(Min + (e+1)*Dist);
    double Average = 0.5*(Kole(EMin/1000) + Kole(EMax/1000));
    Flux += Average*(EMax-EMin);
  }

  MString Comments;
  Comments += "# Albedo neutrons determinined after Kole+ 2014\n";
  Comments += "# Assumptions: \n";
  Comments += "# * Angular distribution is flat out to the Earth-horizon (no limb brightening!) \n";
  Comments += "# * Solar activity is assumed to be half way between solar minimum and maximum \n";
  Comments += "# * The inclination was used to approximate the average Magnetic Latitude \n";

  WriteEnergyFile("AlbedoNeutronsKole", Spectrum);
  WriteAngleFile("AlbedoNeutronsKole", Angle);
  WriteSourceFile("AlbedoNeutronsKole", Flux, 6, Comments);

  m_HadronicSummaryFiles.push_back("AlbedoNeutronsKole");

  return true;
}


/******************************************************************************
 * Generate a cosmic photon spectrum after Gruber et al. 1999
 */
bool BackgroundGenerator::GenerateCosmicPhotonGruber()
{
  mout<<endl;
  mout<<"Generating a cosmic photon spectrum according to Gruber..."<<endl;

  if (IsValid(1, 1000000000) == false) {
    return false;
  }

  auto Gruber = [this] (double EnergykeV) {
    double Flux = 0.0;
    if (EnergykeV < 60) {
      Flux = 7.877*pow(EnergykeV, -1.29)*exp(-EnergykeV/41.13);
    } else {
      Flux = 0.0259*pow(60, 5.5)*pow(EnergykeV, -6.5) +
             0.504*pow(60, 1.58)*pow(EnergykeV, -2.58) +
             0.0288*pow(60, 1.05)*pow(EnergykeV, -2.05);
    }
    return Flux;
  };


  vector<double> Spectrum;
  for (unsigned int b = 0; b < m_NEnergyBins; ++b) {
    // Gruber 1999, Formula (1) in ph/cm2/s/sr/keV
    if (m_EnergyBins[b] < 60) {
      Spectrum.push_back(7.877*pow(m_EnergyBins[b], -1.29)*exp(-m_EnergyBins[b]/41.13));
    } else {
      Spectrum.push_back(0.0259*pow(60, 5.5)*pow(m_EnergyBins[b], -6.5) +
                         0.504*pow(60, 1.58)*pow(m_EnergyBins[b], -2.58) +
                        0.0288*pow(60, 1.05)*pow(m_EnergyBins[b], -2.05));
    }
  }

  vector<double> Angle;
  for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
    if (m_AngleBins[b] < m_HorizonAngle) {
      Angle.push_back(1.0);
    } else {
      Angle.push_back(0.0);
    }
  }

  // Calculate flux:

  // Integration factor over sphere
  double AngleFactor = 2*c_Pi * (1 - cos(m_HorizonAngle*c_Rad));

  // Integration over spectrum:
  double Flux = 0.0;
  if (m_EnergyMax <= 60) {
    // We need to perform a numerical integration...
    int Bins = 10000;
    double Dist = (m_EnergyMax-m_EnergyMin)/Bins;

    // We have fine enough binning that this approach is ok:
    for (double e = 0; e <= Bins-1; ++e) {
      double EMin = m_EnergyMin + e*Dist;
      double EMax = m_EnergyMin + (e+1)*Dist;
      double Average = 0.5*(7.877*pow(EMin, -1.29)*exp(-EMin/41.13) + 7.877*pow(EMax, -1.29)*exp(-EMax/41.13));
      Flux += Average*(EMax-EMin);
    }

  } else if (m_EnergyMin >= 60) {
    Flux += 0.0259/(pow(60, -5.5)*(-6.5 + 1)) * (pow(m_EnergyMax, -5.5) - pow(m_EnergyMin, -5.5));
    Flux += 0.504/(pow(60, -1.58)*(-2.58 + 1)) * (pow(m_EnergyMax, -1.58) - pow(m_EnergyMin, -1.58));
    Flux += 0.0288/(pow(60, -1.05)*(-2.05 + 1)) * (pow(m_EnergyMax, -1.05) - pow(m_EnergyMin, -1.05));
  } else {
    // We need to perform a numerical integration...
    int Bins = 10000;
    double Dist = (60-m_EnergyMin)/Bins;

    // We have fine enough binning that this approach is ok:
    for (double e = 0; e <= Bins-1; ++e) {
      double EMin = m_EnergyMin + e*Dist;
      double EMax = m_EnergyMin + (e+1)*Dist;
      double Average = 0.5*(7.877*pow(EMin, -1.29)*exp(-EMin/41.13) + 7.877*pow(EMax, -1.29)*exp(-EMax/41.13));
      Flux += Average*(EMax-EMin);
    }

    Flux += 0.0259/(pow(60, -5.5)*(-6.5 + 1)) * (pow(m_EnergyMax, -5.5) - pow(60, -5.5));
    Flux += 0.504/(pow(60, -1.58)*(-2.58 + 1)) * (pow(m_EnergyMax, -1.58) - pow(60, -1.58));
    Flux += 0.0288/(pow(60, -1.05)*(-2.05 + 1)) * (pow(m_EnergyMax, -1.05) - pow(60, -1.05));
  }
  cout<<"Anglefactor: "<<AngleFactor<<endl;
  cout<<"Flux: "<<Flux<<" ph/cm2/s/sr"<<endl;
  Flux *= AngleFactor;
  cout<<"Average flux: "<<Flux<<" ph/cm2/s"<<endl;




  WriteEnergyFile("CosmicPhotonsGruber", Spectrum);
  WriteAngleFile("CosmicPhotonsGruber", Angle);
  WriteSourceFile("CosmicPhotonsGruber", Flux, 1);

  m_PhotonicSummaryFiles.push_back("CosmicPhotonsGruber");

  // Determine the flux in ph/s/cm2 via numerical integration...
  double Flux2 = 0;
  int Bins = 10000;
  double Min = log(m_EnergyMin);
  double Max = log(m_EnergyMax);
  double Dist = (Max-Min)/Bins;
  for (double e = 0; e <= Bins-1; ++e) {
    double EMin = exp(Min + e*Dist);
    double EMax = exp(Min + (e+1)*Dist);
    double Average = 0.5*(Gruber(EMin) + Gruber(EMax));
    Flux2 += Average*(EMax-EMin);
  }

  cout<<"Angle-factor: "<<AngleFactor<<endl;
  cout<<"Flux sanity check: "<<Flux2<<" ph/cm2/s/sr"<<endl;
  Flux2 *= AngleFactor;
  cout<<"Flux sanity check:: "<<Flux2<<" ph/cm2/s"<<endl;

  return true;
}


/******************************************************************************
 * Generate a cosmic electron file from Mizuno 2004
 */
bool BackgroundGenerator::GenerateCosmicElectronsMizuno()
{
  auto Mizuno = [this] (double EnergykeV) {
    double A = 0.65; // counts / s / m2 / sr / MeV
    double Index = 3.3;

    double EnergyGeV = 0.000001*EnergykeV;

    double Rigidity = sqrt(EnergyGeV*EnergyGeV + 2*EnergyGeV*0.000511); // GV!

    // The unmodulated flux:
    double Flux = A*pow(Rigidity, -Index); // Formulas (1) & (13) in Mizuno are for a rigidity in GV

    // The flux due to geomagnetic cut-off:
    Flux *= 1.0 / (1.0 + pow(Rigidity/m_AverageGeomagneticCutOff, -6));

    return Flux / 1000.0 / 10000.0; // Switch from n/MeV/m2/s to n/keV/cm2/s
  };

  if (m_IsEarthOrbit == false) return true;

  mout<<endl;
  mout<<"Generating a cosmic electron spectrum using Mizuno 2004 approximations..."<<endl;

  if (IsValid(1, 1000000000) == false) {
    return false;
  }


  vector<double> Angle;
  for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
    if (m_AngleBins[b] >= m_HorizonAngle) {
      Angle.push_back(1.0);
    } else {
      Angle.push_back(0.0);
    }
  }

  vector<double> Spectrum;
  for (unsigned int b = 0; b < m_NEnergyBins; ++b) {
    Spectrum.push_back(Mizuno(m_EnergyBins[b]));
  }


  // Determine the flux in ph/s/cm2 via numerical integration...
  double Flux = 0;
  int Bins = 10000;
  double Min = log(m_EnergyMin);
  double Max = log(m_EnergyMax);
  double Dist = (Max-Min)/Bins;
  for (double e = 0; e <= Bins-1; ++e) {
    double EMin = exp(Min + e*Dist);
    double EMax = exp(Min + (e+1)*Dist);
    double Average = 0.5*(Mizuno(EMin) + Mizuno(EMax));
    Flux += Average*(EMax-EMin);
  }
  // Integration factor over sphere
  double AngleFactor = 2*c_Pi * (1 - cos(m_HorizonAngle*c_Rad));

  cout<<"Angle-factor: "<<AngleFactor<<endl;
  cout<<"Flux: "<<Flux<<" n/cm2/s/sr"<<endl;
  Flux *= AngleFactor;
  cout<<"Flux: "<<Flux<<" n/cm2/s"<<endl;

  MString Comments;
  Comments += "# Cosmic electrons after Mizuno 2004\n";
  Comments += "# Assumptions: \n";
  Comments += "# * Solar modulation has been ignored! \n";
  Comments += "# * Angular distribution is flat out to the Earth-horizon \n";

  WriteEnergyFile("CosmicElectronsMizuno", Spectrum);
  WriteAngleFile("CosmicElectronsMizuno", Angle);
  WriteSourceFile("CosmicElectronsMizuno", Flux, 3, Comments);

  m_LeptonicSummaryFiles.push_back("CosmicElectronsMizuno");

  return true;
}


/******************************************************************************
 * Generate a cosmic positron file from Mizuno 2004
 */
bool BackgroundGenerator::GenerateCosmicPositronsMizuno()
{
  auto Mizuno = [this] (double EnergykeV) {
    double A = 0.055; // counts / s / m2 / sr / MeV
    double Index = 3.3;

    double EnergyGeV = 0.000001*EnergykeV;

    double Rigidity = sqrt(EnergyGeV*EnergyGeV + 2*EnergyGeV*0.000511); // GV!

    // The unmodulated flux:
    double Flux = A*pow(Rigidity, -Index); // Formulas (1) & (13) in Mizuno are for a rigidity in GV

    // The flux due to geomagnetic cut-off:
    Flux *= 1.0 / (1.0 + pow(Rigidity/m_AverageGeomagneticCutOff, -6));

    return Flux / 1000.0 / 10000.0; // Switch from n/MeV/m2/s to n/keV/cm2/s
  };

  if (m_IsEarthOrbit == false) return true;

  mout<<endl;
  mout<<"Generating a cosmic positron spectrums using Mizuno 2004 approximations..."<<endl;

  if (IsValid(1, 1000000000) == false) {
    return false;
  }


  vector<double> Angle;
  for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
    if (m_AngleBins[b] >= m_HorizonAngle) {
      Angle.push_back(1.0);
    } else {
      Angle.push_back(0.0);
    }
  }

  vector<double> Spectrum;
  for (unsigned int b = 0; b < m_NEnergyBins; ++b) {
    Spectrum.push_back(Mizuno(m_EnergyBins[b]));
  }


  // Determine the flux in ph/s/cm2 via numerical integration...
  double Flux = 0;
  int Bins = 10000;
  double Min = log(m_EnergyMin);
  double Max = log(m_EnergyMax);
  double Dist = (Max-Min)/Bins;
  for (double e = 0; e <= Bins-1; ++e) {
    double EMin = exp(Min + e*Dist);
    double EMax = exp(Min + (e+1)*Dist);
    double Average = 0.5*(Mizuno(EMin) + Mizuno(EMax));
    Flux += Average*(EMax-EMin);
  }
  // Integration factor over sphere
  double AngleFactor = 2*c_Pi * (1 - cos(m_HorizonAngle*c_Rad));

  cout<<"Angle-factor: "<<AngleFactor<<endl;
  cout<<"Flux: "<<Flux<<" n/cm2/s/sr"<<endl;
  Flux *= AngleFactor;
  cout<<"Flux: "<<Flux<<" n/cm2/s"<<endl;

  MString Comments;
  Comments += "# Cosmic positrons after Mizuno 2004\n";
  Comments += "# Assumptions: \n";
  Comments += "# * Solar modulation has been ignored! \n";
  Comments += "# * Angular distribution is flat out to the Earth-horizon \n";

  WriteEnergyFile("CosmicPositronsMizuno", Spectrum);
  WriteAngleFile("CosmicPositronsMizuno", Angle);
  WriteSourceFile("CosmicPositronsMizuno", Flux, 2, Comments);

  m_LeptonicSummaryFiles.push_back("CosmicPositronsMizuno");

  return true;
}


/******************************************************************************
 * Generate a cosmic proton spectrum from SPENVIS files
 */
bool BackgroundGenerator::GenerateCosmicProtonsSpenvis()
{
  mout<<endl;
  mout<<"Generating a cosmic proton file using SPENVIS data..."<<endl;

  if (IsValid(1, 1000000000) == false) {
    return false;
  }

  ifstream in;
  in.open(m_AverageCosmicProtonsSpenvis);
  if (in.is_open() == false) {
    mout<<"Error: Unable to open file "<<m_AverageCosmicProtonsSpenvis<<endl;
    return false;
  }

  vector<double> Energies;
  vector<double> Fluxes;

  MString Line;
  bool Start = false;
  while (in.good() == true) {
    Line.ReadLine(in);
    if (Line.Length() < 2) continue;
    if (Line.Contains("DFlux") == true) {
      Start = true;
      continue;
    }
    if (Line.Contains("End of File") == true) break;
    if (Start == false) continue;
    vector<MString> Tokens = Line.Tokenize(",");
    double Energy = 1000*atof(Tokens[0]);
    double Flux = atof(Tokens[2])/1000/10000;
    if (Flux > 0) {
      //cout<<Line<<endl;
      Energies.push_back(Energy);
      Fluxes.push_back(Flux);
    }
  }

  /*
  cout<<"Data: "<<endl;
  for (unsigned int e = 0; e < Energies.size(); ++e) {
    cout<<Energies[e]<<":"<<Fluxes[e]<<endl;
  }
  */

  auto Spenvis = [Energies, Fluxes] (double EnergykeV) {

    double Flux = 0.0;
    if (EnergykeV < Energies[0]) {
      double m = (log(Fluxes[1]) - log(Fluxes[0]))/(log(Energies[1]) - log(Energies[0]));
      double t = log(Fluxes[1]) - m * log(Energies[1]);
      double logy = m * log(EnergykeV) + t;
      Flux = exp(logy);
    } else if (EnergykeV > Energies.back()) {
      int Last = Energies.size() - 1;
      double m = (log(Fluxes[Last]) - log(Fluxes[Last-1]))/(log(Energies[Last]) - log(Energies[Last-1]));
      double t = log(Fluxes[Last]) - m * log(Energies[Last]);
      double logy = m * log(EnergykeV) + t;
      Flux = exp(logy);
    } else {
      for (unsigned int e = 1; e < Energies.size(); ++e) {
        if (EnergykeV < Energies[e]) {
          double m = (log(Fluxes[e]) - log(Fluxes[e-1]))/(log(Energies[e]) - log(Energies[e-1]));
          double t = log(Fluxes[e]) - m * log(Energies[e]);
          double logy = m * log(EnergykeV) + t;
          Flux = exp(logy);
          break;
        }
      }
    }

    return Flux;
  };



  vector<double> Angle;
  for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
    if (m_AngleBins[b] <= m_HorizonAngle) {
      Angle.push_back(1.0);
    } else {
      Angle.push_back(0.0);
    }
  }

  vector<double> Spectrum;
  for (unsigned int b = 0; b < m_NEnergyBins; ++b) {
    Spectrum.push_back(Spenvis(m_EnergyBins[b]));
  }


  // Determine the flux in ph/s/cm2 via numerical integration...
  double Flux = 0;
  int Bins = 10000;
  double Min = log(m_EnergyMin);
  double Max = log(m_EnergyMax);
  double Dist = (Max-Min)/Bins;
  for (double e = 0; e <= Bins-1; ++e) {
    double EMin = exp(Min + e*Dist);
    double EMax = exp(Min + (e+1)*Dist);
    double Average = 0.5*(Spenvis(EMin) + Spenvis(EMax));
    Flux += Average*(EMax-EMin);
  }
  // Integration factor over sphere
  double AngleFactor = 2*c_Pi * (1 - cos(m_HorizonAngle*c_Rad));

  cout<<"Angle-factor: "<<AngleFactor<<endl;
  cout<<"Flux: "<<Flux<<" n/cm2/s/sr"<<endl;
  Flux *= AngleFactor;
  cout<<"Flux: "<<Flux<<" n/cm2/s"<<endl;

  MString Comments;
  Comments += "# Cosmic protons using SPENVIS\n";
  Comments += "# Assumptions: \n";
  Comments += "# * log-log extrapolation beyond upper energy limit from SPENVIS\n";
  Comments += "# * Angular distribution is flat out to the Earth-horizon \n";

  WriteEnergyFile("CosmicProtonsSpenvis", Spectrum);
  WriteAngleFile("CosmicProtonsSpenvis", Angle);
  WriteSourceFile("CosmicProtonsSpenvis", Flux, 4, Comments);

  m_HadronicSummaryFiles.push_back("CosmicProtonsSpenvis");

  return true;
}


/******************************************************************************
 * Generate a cosmic alpha particle spectrum from SPENVIS files
 */
bool BackgroundGenerator::GenerateCosmicAlphasSpenvis()
{
  mout<<endl;
  mout<<"Generating a cosmic alpha particle spectrum using SPENVIS data..."<<endl;

  if (IsValid(1, 1000000000) == false) {
    return false;
  }

  ifstream in;
  in.open(m_AverageCosmicAlphasSpenvis);
  if (in.is_open() == false) {
    mout<<"Error: Unable to open file "<<m_AverageCosmicAlphasSpenvis<<endl;
    return false;
  }

  vector<double> Energies;
  vector<double> Fluxes;

  MString Line;
  bool Start = false;
  while (in.good() == true) {
    Line.ReadLine(in);
    if (Line.Length() < 2) continue;
    if (Line.Contains("DFlux") == true) {
      Start = true;
      continue;
    }
    if (Line.Contains("End of File") == true) break;
    if (Start == false) continue;
    vector<MString> Tokens = Line.Tokenize(",");
    if (Tokens.size() != 3) {
      cout<<"Error: The line should have three tokens, energy, integrated flux, and differential flux"<<endl;
      cout<<Line<<endl;
      continue;
    }
    double Energy = 1000*atof(Tokens[0]) * 4; // from MeV to keV and x4 for 4 nucleons
    double Flux = atof(Tokens[2])/1000/10000 / 4; // Flux is per nucleus thus 1/4th
    if (Flux > 0) {
      //cout<<Line<<endl;
      Energies.push_back(Energy);
      Fluxes.push_back(Flux);
    }
  }

  /*
  cout<<"Data: "<<endl;
  for (unsigned int e = 0; e < Energies.size(); ++e) {
    cout<<Energies[e]<<":"<<Fluxes[e]<<endl;
  }
  */

  auto Spenvis = [Energies, Fluxes] (double EnergykeV) {

    double Flux = 0.0;
    if (EnergykeV < Energies[0]) {
      double m = (log(Fluxes[1]) - log(Fluxes[0]))/(log(Energies[1]) - log(Energies[0]));
      double t = log(Fluxes[1]) - m * log(Energies[1]);
      double logy = m * log(EnergykeV) + t;
      Flux = exp(logy);
    } else if (EnergykeV > Energies.back()) {
      int Last = Energies.size() - 1;
      double m = (log(Fluxes[Last]) - log(Fluxes[Last-1]))/(log(Energies[Last]) - log(Energies[Last-1]));
      double t = log(Fluxes[Last]) - m * log(Energies[Last]);
      double logy = m * log(EnergykeV) + t;
      Flux = exp(logy);
    } else {
      for (unsigned int e = 1; e < Energies.size(); ++e) {
        if (EnergykeV < Energies[e]) {
          double m = (log(Fluxes[e]) - log(Fluxes[e-1]))/(log(Energies[e]) - log(Energies[e-1]));
          double t = log(Fluxes[e]) - m * log(Energies[e]);
          double logy = m * log(EnergykeV) + t;
          Flux = exp(logy);
          break;
        }
      }
    }

    return Flux; // p/cm2/s/keV/(MeV/n)  <-- switching from MeV/n to MeV in energy dimension does not change flux since we just change scale on x-axis!
  };



  vector<double> Angle;
  for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
    if (m_AngleBins[b] <= m_HorizonAngle) {
      Angle.push_back(1.0);
    } else {
      Angle.push_back(0.0);
    }
  }

  vector<double> Spectrum;
  for (unsigned int b = 0; b < m_NEnergyBins; ++b) {
    Spectrum.push_back(Spenvis(m_EnergyBins[b]));
  }


  // Determine the flux in ph/s/cm2 via numerical integration...
  double Flux = 0;
  int Bins = 10000;
  double Min = log(m_EnergyMin);
  double Max = log(m_EnergyMax);
  double Dist = (Max-Min)/Bins;
  for (double e = 0; e <= Bins-1; ++e) {
    double EMin = exp(Min + e*Dist);
    double EMax = exp(Min + (e+1)*Dist);
    double Average = 0.5*(Spenvis(EMin) + Spenvis(EMax));
    Flux += Average*(EMax-EMin);
  }
  // Integration factor over sphere
  double AngleFactor = 2*c_Pi * (1 - cos(m_HorizonAngle*c_Rad));

  cout<<"Angle-factor: "<<AngleFactor<<endl;
  cout<<"Flux: "<<Flux<<" n/cm2/s/sr"<<endl;
  Flux *= AngleFactor;
  cout<<"Flux: "<<Flux<<" n/cm2/s"<<endl;

  MString Comments;
  Comments += "# Cosmic alpha particles using SPENVIS\n";
  Comments += "# Assumptions: \n";
  Comments += "# * log-log extrapolation beyond upper energy limit from SPENVIS\n";
  Comments += "# * Angular distribution is flat out to the Earth-horizon \n";

  WriteEnergyFile("CosmicAlphasSpenvis", Spectrum);
  WriteAngleFile("CosmicAlphasSpenvis", Angle);
  WriteSourceFile("CosmicAlphasSpenvis", Flux, 21, Comments);

  m_HadronicSummaryFiles.push_back("CosmicAlphasSpenvis");

  return true;
}


/******************************************************************************
 * Generate an average trapped protons and electrons spectrum from SPENVIS files
 */
bool BackgroundGenerator::GenerateTrappedProtonsElectronsSpenvis()
{
  mout<<endl;
  mout<<"Generating an average trapped proton and electron spectrum based on SPENVIS data..."<<endl;

  if (IsValid(1, 1000000000) == false) {
    return false;
  }

  ifstream in;
  in.open(m_AverageTrappedProtonsElectronsSpenvis);
  if (in.is_open() == false) {
    mout<<"Error: Unable to open file "<<m_AverageTrappedProtonsElectronsSpenvis<<endl;
    return false;
  }

  vector<double> ProtonEnergies;
  vector<double> ProtonFluxes;

  vector<double> ElectronEnergies;
  vector<double> ElectronFluxes;

  MString Line;
  bool IsProtons = true;
  bool Start = false;
  while (in.good() == true) {
    Line.ReadLine(in);
    if (Line.Length() < 2) continue;
    if (Line.Contains("DFlux") == true) {
      Start = true;
      continue;
    }
    if (Line.Contains("End of Block") == true) {
      IsProtons = false;
      Start = false;
    }
    if (Line.Contains("End of File") == true) break;
    if (Start == false) continue;
    vector<MString> Tokens = Line.Tokenize(",");
    if (Tokens.size() != 3) {
      cout<<"Error: The line should have three tokens, energy, integrated flux, differential flux"<<endl;
      cout<<Line<<endl;
      continue;
    }
    double Energy = 1000*atof(Tokens[0]); //
    double Flux = atof(Tokens[2])/1000/4/c_Pi; // <-- Assuming 4 pi
    if (Flux > 0) {
      //cout<<Line<<endl;
      if (IsProtons == true) {
        ProtonEnergies.push_back(Energy);
        ProtonFluxes.push_back(Flux);
      } else {
        ElectronEnergies.push_back(Energy);
        ElectronFluxes.push_back(Flux);
      }
    }
  }

  /*
  cout<<"Data protons: "<<endl;
  for (unsigned int e = 0; e < ProtonEnergies.size(); ++e) {
    cout<<ProtonEnergies[e]<<":"<<ProtonFluxes[e]<<endl;
  }

  cout<<"Data electrons: "<<endl;
  for (unsigned int e = 0; e < ElectronEnergies.size(); ++e) {
    cout<<ElectronEnergies[e]<<":"<<ElectronFluxes[e]<<endl;
  }
  */

  auto ProtonsSpenvis = [ProtonEnergies, ProtonFluxes] (double EnergykeV) {

    if (ProtonEnergies.size() < 2) return 0.0;

    double Flux = 0.0;
    if (EnergykeV < ProtonEnergies[0]) {
      double m = (log(ProtonFluxes[1]) - log(ProtonFluxes[0]))/(log(ProtonEnergies[1]) - log(ProtonEnergies[0]));
      double t = log(ProtonFluxes[1]) - m * log(ProtonEnergies[1]);
      double logy = m * log(EnergykeV) + t;
      Flux = exp(logy);
    } else if (EnergykeV > ProtonEnergies.back()) {
      int Last = ProtonEnergies.size() - 1;
      double m = (log(ProtonFluxes[Last]) - log(ProtonFluxes[Last-1]))/(log(ProtonEnergies[Last]) - log(ProtonEnergies[Last-1]));
      double t = log(ProtonFluxes[Last]) - m * log(ProtonEnergies[Last]);
      double logy = m * log(EnergykeV) + t;
      Flux = exp(logy);
    } else {
      for (unsigned int e = 1; e < ProtonEnergies.size(); ++e) {
        if (EnergykeV < ProtonEnergies[e]) {
          double m = (log(ProtonFluxes[e]) - log(ProtonFluxes[e-1]))/(log(ProtonEnergies[e]) - log(ProtonEnergies[e-1]));
          double t = log(ProtonFluxes[e]) - m * log(ProtonEnergies[e]);
          double logy = m * log(EnergykeV) + t;
          Flux = exp(logy);
          break;
        }
      }
    }

    return Flux; // p/cm2/s/keV/sr
  };

  auto ElectronsSpenvis = [ElectronEnergies, ElectronFluxes] (double EnergykeV) {

    if (ElectronEnergies.size() < 2) return 0.0;

    double Flux = 0.0;
    if (EnergykeV < ElectronEnergies[0]) {
      double m = (log(ElectronFluxes[1]) - log(ElectronFluxes[0]))/(log(ElectronEnergies[1]) - log(ElectronEnergies[0]));
      double t = log(ElectronFluxes[1]) - m * log(ElectronEnergies[1]);
      double logy = m * log(EnergykeV) + t;
      Flux = exp(logy);
    } else if (EnergykeV > ElectronEnergies.back()) {
      int Last = ElectronEnergies.size() - 1;
      double m = (log(ElectronFluxes[Last]) - log(ElectronFluxes[Last-1]))/(log(ElectronEnergies[Last]) - log(ElectronEnergies[Last-1]));
      double t = log(ElectronFluxes[Last]) - m * log(ElectronEnergies[Last]);
      double logy = m * log(EnergykeV) + t;
      Flux = exp(logy);
    } else {
      for (unsigned int e = 1; e < ElectronEnergies.size(); ++e) {
        if (EnergykeV < ElectronEnergies[e]) {
          double m = (log(ElectronFluxes[e]) - log(ElectronFluxes[e-1]))/(log(ElectronEnergies[e]) - log(ElectronEnergies[e-1]));
          double t = log(ElectronFluxes[e]) - m * log(ElectronEnergies[e]);
          double logy = m * log(EnergykeV) + t;
          Flux = exp(logy);
          break;
        }
      }
    }

    return Flux; // p/cm2/s/keV/sr
  };



  vector<double> Angle;
  for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
    Angle.push_back(1.0);
  }
  // Integration factor over sphere
  double AngleFactor = 4*c_Pi;

  vector<double> ProtonSpectrum;
  for (unsigned int b = 0; b < m_NEnergyBins; ++b) {
    ProtonSpectrum.push_back(ProtonsSpenvis(m_EnergyBins[b]));
  }
  vector<double> ElectronSpectrum;
  for (unsigned int b = 0; b < m_NEnergyBins; ++b) {
    ElectronSpectrum.push_back(ElectronsSpenvis(m_EnergyBins[b]));
  }


  // Determine the flux in ph/s/cm2 via numerical integration...
  double ProtonFlux = 0;
  double ElectronFlux = 0;
  int Bins = 10000;
  double Min = log(m_EnergyMin);
  double Max = log(m_EnergyMax);
  double Dist = (Max-Min)/Bins;
  for (double e = 0; e <= Bins-1; ++e) {
    double EMin = exp(Min + e*Dist);
    double EMax = exp(Min + (e+1)*Dist);
    double Average = 0.5*(ProtonsSpenvis(EMin) + ProtonsSpenvis(EMax));
    ProtonFlux += Average*(EMax-EMin);
    Average = 0.5*(ElectronsSpenvis(EMin) + ElectronsSpenvis(EMax));
    ElectronFlux += Average*(EMax-EMin);
  }

  cout<<"Angle-factor: "<<AngleFactor<<endl;
  cout<<"Proton Flux: "<<ProtonFlux<<" p/cm2/s/sr"<<endl;
  ProtonFlux *= AngleFactor;
  cout<<"Proton Flux: "<<ProtonFlux<<" p/cm2/s"<<endl;
  cout<<"Electron Flux: "<<ElectronFlux<<" e/cm2/s/sr"<<endl;
  ElectronFlux *= AngleFactor;
  cout<<"Electron Flux: "<<ElectronFlux<<" e/cm2/s"<<endl;

  MString Comments;
  Comments += "# Trapped protons particles using SPENVIS\n";
  Comments += "# Assumptions: \n";
  Comments += "# * The values are orbit averaged - thus do NOT simulate alongside other hadrons: \n";
  Comments += "#   Since we are OFF during SAA passages, do a separate simulation of the trapped protons \n";
  Comments += "#   and then just include the activations not the prompt decays for normal background calculations. \n";
  Comments += "# * Angular distribution is assumed isotropic \n";

  WriteEnergyFile("TrappedProtonsSpenvis", ProtonSpectrum);
  WriteAngleFile("TrappedProtonsSpenvis", Angle);
  WriteSourceFile("TrappedProtonsSpenvis", ProtonFlux, 4, Comments);

  m_TrappedHadronicSummaryFiles.push_back("TrappedProtonsSpenvis");

  Comments = "";
  Comments += "# Trapped electron particles using SPENVIS\n";
  Comments += "# Assumptions: \n";
  Comments += "# * The values are orbit averaged. \n";
  Comments += "# * DO NOT include in normal background calculations since we are OFF during SAA passages.\n";
  Comments += "# * Angular distribution is assumed isotropic \n";

  WriteEnergyFile("TrappedElectronsSpenvis", ElectronSpectrum);
  WriteAngleFile("TrappedElectronsSpenvis", Angle);
  WriteSourceFile("TrappedElectronsSpenvis", ElectronFlux, 3, Comments);

  return true;
}



/******************************************************************************
 * Check if the chosen spectrum is within the validity range model
 */
bool BackgroundGenerator::WriteEnergyFile(MString SourceName, vector<double> Spectrum)
{
  MString EnergyFileName = SourceName;
  EnergyFileName += ".spectrum.dat";

  ofstream fout;
  fout.open(EnergyFileName);
  if (fout.is_open() == false) {
    mout<<"Unable to open: "<<EnergyFileName<<endl;
    return false;
  }

  if (m_IsMGGPOD == false) {
    fout<<"# Spectrum file for cosima"<<endl;
    fout<<endl;
    fout<<"IP LIN"<<endl;
    fout<<endl;

    for (unsigned int b = 0; b < m_EnergyBins.size(); ++b) {
      fout<<"DP "<<m_EnergyBins[b]<<" "<<Spectrum[b]<<endl;
    }
    fout<<endl;
  } else {
    fout.setf(ios_base::scientific, ios_base::floatfield);
    fout.precision(6);
    for (unsigned int b = 0; b < m_EnergyBins.size(); ++b) {
      fout<<" "<<setw(11)<<m_EnergyBins[b]/1000000<<"  "<<setw(11)<<Spectrum[b]*1000000*10000<<endl;
    }
    fout<<endl;
  }

  fout.close();

  return true;
}


/******************************************************************************
 * Check if the chosen spectrum is within the validity range model
 */
bool BackgroundGenerator::WriteAngleFile(MString SourceName, vector<double> Angle)
{
  MString AngleFileName = SourceName;
  AngleFileName += ".beam.dat";

  ofstream fout;
  fout.open(AngleFileName);
  if (fout.is_open() == false) {
    mout<<"Unable to open: "<<AngleFileName<<endl;
    return false;
  }

  if (m_IsMGGPOD == false) {
    fout<<"# Beam file for cosima"<<endl;
    fout<<endl;
    fout<<"IP LIN"<<endl;
    fout<<endl;

    for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
      fout<<"DP "<<m_AngleBins[b]<<" "<<Angle[b]<<endl;
    }
    fout<<endl;
  } else {
    fout.setf(ios_base::scientific, ios_base::floatfield);
    fout.precision(6);
    for (unsigned int b = 0; b < m_AngleBins.size(); ++b) {
      fout<<" "<<setw(11)<<m_AngleBins[b]<<" "<<setw(11)<<Angle[b]<<endl;
    }
    fout<<endl;
  }

  fout.close();

  return true;
}


/******************************************************************************
 * Check if the chosen spectrum is within the validity range model
 */
bool BackgroundGenerator::WriteSourceFile(MString SourceName, double Flux, int ParticleID, MString Comments, MString SpectralString )
{
  MString FileName = SourceName;
  FileName += ".partial.source";

  MString EnergyFileName = SourceName;
  EnergyFileName += ".spectrum.dat";

  MString AngleFileName = SourceName;
  AngleFileName += ".beam.dat";

  ofstream fout;
  fout.open(FileName);
  if (fout.is_open() == false) {
    mout<<"Unable to open: "<<FileName<<endl;
    return false;
  }

  fout<<"# Partial beam file for cosima"<<endl;
  fout<<endl;
  //fout<<"Include Common.partial.source"<<endl;
  //fout<<endl;
  //fout<<"SpaceSim.FileName      "<<SourceName<<endl;
  //fout<<"SpaceSim.Time          10.0"<<endl;
  //fout<<endl;
  fout<<Comments;
  //fout<<endl;
  fout<<"SpaceSim.Source "<<SourceName<<endl;
  fout<<SourceName<<".ParticleType           "<<ParticleID<<endl;
  if (m_IsEarthOrbit == true) {
    fout<<SourceName<<".Beam                   FarFieldFileZenithDependent "<<AngleFileName<<endl;
  } else {
    fout<<SourceName<<".Beam                   FarFieldAreaSource 0 180 0 360"<<endl;
  }
  if (SpectralString != "") {
    fout<<SourceName<<".Spectrum               "<<SpectralString<<endl;
  } else {
    fout<<SourceName<<".Spectrum               File "<<EnergyFileName<<endl;
  }
  fout<<SourceName<<".Flux                   "<<Flux<<endl;
  fout<<endl;
  fout.close();

  return true;
}


/******************************************************************************
 * Write the final summary files
 */
bool BackgroundGenerator::WriteSummaryFiles()
{
  ofstream fout;
  fout.open("Common.partial.source");
  if (fout.is_open() == false) {
    mout<<"Unable to open: Common.partial.source"<<endl;
    return false;
  }

  fout<<"# Common setup option for all simulations"<<endl;
  fout<<endl;
  fout<<"# Version of the source file"<<endl;
  fout<<"Version 1"<<endl;
  fout<<endl;
  fout<<"# Geometry"<<endl;
  fout<<"Geometry             > You geometry here <"<<endl;
  fout<<"DetectorTimeConstant 0.000005"<<endl;
  fout<<endl;
  fout<<"# Output formats"<<endl;
  fout<<"StoreSimulationInfo  init-only"<<endl;
  fout<<endl;

  fout.close();


  ofstream pout;
  pout.open("PhotonicComponents.source");
  if (pout.is_open() == false) {
    mout<<"Unable to open: PhotonicComponents.source"<<endl;
    return false;
  }

  pout<<"# Components file for the photon simulations"<<endl;
  pout<<endl;
  pout<<"Include Common.partial.source"<<endl;
  pout<<endl;
  pout<<"# Physics list"<<endl;
  pout<<"PhysicsListEM LivermorePol"<<endl;
  pout<<endl;
  pout<<"# The simulation run"<<endl;
  pout<<"Run SpaceSim"<<endl;
  pout<<"SpaceSim.FileName  PhotonicBackground"<<endl;
  pout<<"SpaceSim.Events    50000000"<<endl;
  pout<<endl;
  for (MString S: m_PhotonicSummaryFiles) {
    pout<<"Include "<<S<<".partial.source"<<endl;
  }
  pout<<endl;

  pout.close();


  ofstream lout;
  lout.open("LeptonicComponents.source");
  if (lout.is_open() == false) {
    mout<<"Unable to open: LeptonicComponents.source"<<endl;
    return false;
  }

  lout<<"# Components file for the lepton simulations"<<endl;
  lout<<endl;
  lout<<"Include Common.partial.source"<<endl;
  lout<<endl;
  lout<<endl;
  lout<<"# Physics list"<<endl;
  lout<<"PhysicsListHD   qgsp-bic-hp"<<endl;
  lout<<"PhysicsListEM   LivermorePol"<<endl;
  lout<<endl;
  lout<<endl;
  lout<<"Run SpaceSim"<<endl;
  lout<<"SpaceSim.FileName  LeptonicBackground"<<endl;
  lout<<"SpaceSim.Events    50000000"<<endl;
  lout<<endl;
  for (MString S: m_LeptonicSummaryFiles) {
    lout<<"Include "<<S<<".partial.source"<<endl;
  }
  lout<<endl;

  lout.close();



  ofstream hiout;
  hiout.open("HadronicComponentsPrompt.source");
  if (hiout.is_open() == false) {
    mout<<"Unable to open: HadronicComponentsPrompt.source"<<endl;
    return false;
  }

  hiout<<"# Components file for the hadron simulations - prompt"<<endl;
  hiout<<endl;
  hiout<<"Include Common.partial.source"<<endl;
  hiout<<endl;
  hiout<<endl;
  hiout<<"# Physics list"<<endl;
  hiout<<"PhysicsListHD   qgsp-bic-hp"<<endl;
  hiout<<"PhysicsListEM   LivermorePol"<<endl;
  hiout<<"DecayMode       ActivationBuildup"<<endl;
  hiout<<endl;
  hiout<<endl;
  hiout<<"Run SpaceSim"<<endl;
  hiout<<"SpaceSim.FileName                  HadronicBackground"<<endl;
  hiout<<"SpaceSim.Events                    5000000"<<endl;
  hiout<<"SpaceSim.IsotopeProductionFile     HadronicBackgroundIsotopes"<<endl;
  hiout<<endl;
  for (MString S: m_HadronicSummaryFiles) {
    hiout<<"Include "<<S<<".partial.source"<<endl;
  }
  hiout<<endl;

  hiout.close();



  ofstream hbout;
  hbout.open("HadronicComponentsBuildUp.source");
  if (hbout.is_open() == false) {
    mout<<"Unable to open: HadronicComponentsBuildUp.source"<<endl;
    return false;
  }

  hbout<<"# Components file for the hadron simulations - build up"<<endl;
  hbout<<endl;
  hbout<<"Include Common.partial.source"<<endl;
  hbout<<endl;
  hbout<<"# Physics list"<<endl;
  hbout<<"PhysicsListHD                      qgsp-bic-hp"<<endl;
  hbout<<"PhysicsListEM                      LivermorePol"<<endl;
  hbout<<endl;
  hbout<<endl;
  hbout<<"# The activation calculation for one year const irradiation"<<endl;
  hbout<<"Activator A"<<endl;
  hbout<<"A.ActivationMode          ConstantIrradiation  31556736"<<endl;
  hbout<<"A.ActivationFile          HadronicBackgroundActivation.dat"<<endl;
  hbout<<endl;
  hbout<<"# Now generate a source file with all the HadronicBackgroundIsotopes.* files, via"<<endl;
  hbout<<"# User > for i in `ls HadronicBackgroundIsotopes.p1.inc*.dat`; do echo \"A.IsotopeProductionFile $i\" >> HadronicBackgroundIsotopes.source; done"<<endl;
  hbout<<"Include HadronicBackgroundIsotopes.source"<<endl;
  hbout<<endl;

  hbout.close();



  ofstream hdout;
  hdout.open("HadronicComponentsDecay.source");
  if (hdout.is_open() == false) {
    mout<<"Unable to open: HadronicComponentsDecay.source"<<endl;
    return false;
  }

  hdout<<"# Components file for the hadron simulations - decay"<<endl;
  hdout<<endl;
  hdout<<"Include Common.partial.source"<<endl;
  hdout<<endl;
  hdout<<"# Physics list"<<endl;
  hdout<<"PhysicsListHD                      qgsp-bic-hp"<<endl;
  hdout<<"PhysicsListEM                      LivermorePol"<<endl;
  hdout<<"DecayMode                          ActivationDelayedDecay"<<endl;
  hdout<<endl;
  hdout<<"Run ActivationStep3"<<endl;
  hdout<<"ActivationStep3.FileName                         HadronicBackgroundDecay"<<endl;
  hdout<<"ActivationStep3.Events                           5000000"<<endl;
  hdout<<endl;
  hdout<<"ActivationStep3.ActivationSources                HadronicBackgroundActivation.dat"<<endl;
  hdout<<endl;

  hdout.close();





  ofstream tiout;
  tiout.open("TrappedHadronicComponentsPrompt.source");
  if (tiout.is_open() == false) {
    mout<<"Unable to open: TrappedHadronicComponentsPrompt.source"<<endl;
    return false;
  }

  tiout<<"# Components file for the trapped hadron simulations - prompt"<<endl;
  tiout<<endl;
  tiout<<"Include Common.partial.source"<<endl;
  tiout<<endl;
  tiout<<endl;
  tiout<<"# Physics list"<<endl;
  tiout<<"PhysicsListHD   qgsp-bic-hp"<<endl;
  tiout<<"PhysicsListEM   LivermorePol"<<endl;
  tiout<<"DecayMode       ActivationBuildup"<<endl;
  tiout<<endl;
  tiout<<endl;
  tiout<<"Run SpaceSim"<<endl;
  tiout<<"SpaceSim.FileName                  TrappedHadronicBackground"<<endl;
  tiout<<"SpaceSim.Events                    5000000"<<endl;
  tiout<<"SpaceSim.IsotopeProductionFile     TrappedHadronicBackgroundIsotopes"<<endl;
  tiout<<endl;
  for (MString S: m_TrappedHadronicSummaryFiles) {
    tiout<<"Include "<<S<<".partial.source"<<endl;
  }
  tiout<<endl;

  tiout.close();



  ofstream tbout;
  tbout.open("TrappedHadronicComponentsBuildUp.source");
  if (tbout.is_open() == false) {
    mout<<"Unable to open: TrappedHadronicComponentsBuildUp.source"<<endl;
    return false;
  }

  tbout<<"# Components file for the trapped hadron simulations - build up"<<endl;
  tbout<<endl;
  tbout<<"Include Common.partial.source"<<endl;
  tbout<<endl;
  tbout<<"# Physics list"<<endl;
  tbout<<"PhysicsListHD                      qgsp-bic-hp"<<endl;
  tbout<<"PhysicsListEM                      LivermorePol"<<endl;
  tbout<<endl;
  tbout<<endl;
  tbout<<"# The activation calculation for one year const irradiation"<<endl;
  tbout<<"Activator A"<<endl;
  tbout<<"A.ActivationMode          ConstantIrradiation  31556736"<<endl;
  tbout<<"A.ActivationFile          TrapopedHadronicBackgroundActivation.dat"<<endl;
  tbout<<endl;
  tbout<<"# Now generate a source file with all the TrappedHadronicBackgroundIsotopes.* files, via"<<endl;
  tbout<<"# User > for i in `ls TrappedHadronicBackgroundIsotopes.p1.inc*.dat`; do echo \"A.IsotopeProductionFile $i\" >> TrappedHadronicBackgroundIsotopes.source; done"<<endl;
  tbout<<"Include TrappedHadronicBackgroundIsotopes.source"<<endl;
  tbout<<endl;

  tbout.close();



  ofstream tdout;
  tdout.open("TrappedHadronicComponentsDecay.source");
  if (tdout.is_open() == false) {
    mout<<"Unable to open: TrappedHadronicComponentsDecay.source"<<endl;
    return false;
  }

  tdout<<"# Components file for the trapped hadron simulations - decay"<<endl;
  tdout<<endl;
  tdout<<"Include Common.partial.source"<<endl;
  tdout<<endl;
  tdout<<"# Physics list"<<endl;
  tdout<<"PhysicsListHD                      qgsp-bic-hp"<<endl;
  tdout<<"PhysicsListEM                      LivermorePol"<<endl;
  tdout<<"DecayMode                          ActivationDelayedDecay"<<endl;
  tdout<<endl;
  tdout<<"Run ActivationStep3"<<endl;
  tdout<<"ActivationStep3.FileName             TrappedHadronicBackgroundDecay"<<endl;
  tdout<<"ActivationStep3.Events               5000000"<<endl;
  tdout<<endl;
  tdout<<"ActivationStep3.ActivationSources    TrappedHadronicBackgroundActivation.dat"<<endl;
  tdout<<endl;

  tdout.close();



  return false;
}


/******************************************************************************
 * Check if the chosen spectrum is within the validity range model
 */
bool BackgroundGenerator::IsValid(double ModelMin, double ModelMax)
{
  if (m_EnergyMin < ModelMin) {
    mout<<"Choosen minimum energy ("<<m_EnergyMin<<" keV) is below the minimum allowed by the model ("<<ModelMin<<" keV)!"<<endl;
    return false;
  }
  if (m_EnergyMax > ModelMax) {
    mout<<"Choosen maximum energy ("<<m_EnergyMax<<" keV) is above the maximum allowed by the model ("<<ModelMax<<" keV)!"<<endl;
    return false;
  }

  return true;
}



/******************************************************************************/

BackgroundGenerator* g_Prg = 0;
int g_NInterruptCatches = 1;

/******************************************************************************/


/******************************************************************************
 * Called when an interrupt signal is flagged
 * All catched signals lead to a well defined exit of the program
 */
void CatchSignal(int a)
{
  if (g_Prg != 0 && g_NInterruptCatches-- > 0) {
    cout<<"Catched signal Ctrl-C (ID="<<a<<"):"<<endl;
    g_Prg->Interrupt();
  } else {
    abort();
  }
}


/******************************************************************************
 * Main program
 */
int main(int argc, char** argv)
{
  //void (*handler)(int);
  //handler = CatchSignal;
  //(void) signal(SIGINT, CatchSignal);

  // Initialize global MEGALIB variables, especially mgui, etc.
  MGlobal::Initialize();

  TApplication BackgroundGeneratorApp("BackgroundGeneratorApp", 0, 0);

  g_Prg = new BackgroundGenerator();

  if (g_Prg->ParseCommandLine(argc, argv) == false) {
    cerr<<"Error during parsing of command line!"<<endl;
    return -1;
  }
  if (g_Prg->Analyze() == false) {
    cerr<<"Error during analysis!"<<endl;
    return -2;
  }

  //BackgroundGeneratorApp.Run();

  cout<<"Program exited normally!"<<endl;

  return 0;
}

/*
 * Cosima: the end...
 ******************************************************************************/
