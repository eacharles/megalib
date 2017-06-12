/*
* MERTrackKalman3D.cxx
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
// MERTrackKalman3D
//
// Base class for event reconstruction tasks, e.g. find clusters, tracks,
// Compton sequence
//
////////////////////////////////////////////////////////////////////////////////


// Include the header:
#include "MERTrackKalman3D.h"

// Standard libs:
#include <list>
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <functional>
using namespace std;

// ROOT libs:
#include "TMath.h"
#include "TH1F.h"


// MEGAlib libs:
#include "MAssert.h"
#include "MStreams.h"
#include "MRESE.h"
#include "MRESEIterator.h"
#include "MRETrack.h"
#include "MREHit.h"
#include "MMuonEvent.h"
#include "MTimer.h"
#include "MFile.h"

////////////////////////////////////////////////////////////////////////////////


#ifdef ___CINT___
ClassImp(MERTrackKalman3D)
#endif


////////////////////////////////////////////////////////////////////////////////


MERTrackKalman3D::MERTrackKalman3D() : MERTrack()
{
  // Construct an instance of MERTrackKalman3D
  w=pow(tan(60*(TMath::Pi())/180.),2);
  sfg_factor=24;//1.5/sigma;
  max_offset=30.;
  start_nrg_fit=20;
  n_fit = 7;

}


////////////////////////////////////////////////////////////////////////////////


MERTrackKalman3D::~MERTrackKalman3D()
{
  // Delete this instance of MERTrackKalman3D
}


////////////////////////////////////////////////////////////////////////////////


bool MERTrackKalman3D::EvaluateTracks(MRERawEvent* RE)
{
  // Evaluate the given tracks...

  if (RE->GetVertex() != 0) {
    return EvaluatePairs(RE);
  }

  // Mark: Add you evaluations here:
  mout<<"MRawEventList* MERTrackKalman3D::EvaluateTracks(MRERawEvent* ER): nyi!"<<endl;
  //return MERTrack::EvaluateTracks(RE);

  // Store your results:
  RE->SetTrackQualityFactor(MRERawEvent::c_NoQualityFactor);
  RE->SetPairQualityFactor(MRERawEvent::c_NoQualityFactor);

  return true;
}



////////////////////////////////////////////////////////////////////////////////


bool MERTrackKalman3D::SetSpecialParameters(double HeightX0,
  double SigmaHitPos,
  unsigned int NLayersForVertexSearch)
  {
    m_NLayersForVertexSearch = NLayersForVertexSearch;
    m_heightX0 = HeightX0;
    m_sigma = SigmaHitPos;

    if (m_NLayersForVertexSearch < 4) {
      merr<<"Error: Revan (tracking): NLayersForVertexSearch (="<<m_NLayersForVertexSearch<<") must be >= 4. Setting it to 4."<<show;
      m_NLayersForVertexSearch = 4;
    }
    return true;
  }




  ////////////////////////////////////////////////////////////////////////////////


  MString MERTrackKalman3D::ToString(bool CoreOnly) const
  {
    // Dump an options string gor the tra file:

    ostringstream out;

    if (CoreOnly == false) {
      out<<"# Tracking - Kalman options:"<<endl;
      out<<"# "<<endl;
    }
    out<<MERTrack::ToString(true);
    if (CoreOnly == false) {
      out<<"# "<<endl;
    }

    return out.str().c_str();
  }


  ////////////////////////////////////////////////////////////////////////////////


  MRawEventList* MERTrackKalman3D::CheckForPair(MRERawEvent* RE)
  {
    // Check if this event could be a pair:
    //
    // The typical pattern of a pair is:
    // (a) There is a vertex
    // (b) In two layers above/below are at least two hits

    // Search the vertex:
    // The vertex is at this point a (non-ambiguous) track or a single (clustered)
    // hit followed by at least two layers with two hits per layer:

    MRESE* Vertex = 0;

    MRawEventList *List = 0;
    bool OnlyHitInLayer = false;
    unsigned int MaximumLayerJump = 2;
    //MRESE* RESE = 0;


    // Create a list of RESEs dorted by depth in tracker
    vector<MRESE*> ReseList;
    for (int h = 0; h < RE->GetNRESEs(); h++) {
      ReseList.push_back(RE->GetRESEAt(h));
    }
    sort(ReseList.begin(), ReseList.end(), CompareRESEByZ());

    mdebug<<"RESE's sorted by depth: "<<endl;
    vector<MRESE*>::iterator Iterator1;
    vector<MRESE*>::iterator Iterator2;
    for (Iterator1 = ReseList.begin(); Iterator1 != ReseList.end(); Iterator1++) {
      if (IsInTracker(*Iterator1) == false) continue;
      mdebug<<(*Iterator1)->GetID()<<": "<<(*Iterator1)->GetPosition().Z()<<endl;
    }

    // For each of the RESE's in the list check if it could be the first of the vertex
    for (Iterator1 = ReseList.begin(); Iterator1 != ReseList.end(); Iterator1++) {
      if (IsInTracker(*Iterator1) == false) continue;

      // If it is a single hit, and if it is the only one in its layer:
      OnlyHitInLayer = true;
      for (Iterator2 = ReseList.begin(); Iterator2 != ReseList.end(); Iterator2++) {
        if ((*Iterator1) == (*Iterator2)) continue;
        if (IsInTracker(*Iterator2) == false) continue;
        if (m_Geometry->AreInSameLayer((*Iterator1), (*Iterator2)) == true) {
          OnlyHitInLayer = false;
          break;
        }
      }
      if (OnlyHitInLayer == false) continue;
      mdebug<<"Search vertex: Only hit in layer:"<<endl;
      mdebug<<(*Iterator1)->ToString()<<endl;

      // We only have one hit:
      vector<int> NBelow(m_NLayersForVertexSearch, 0);
      vector<int> NAbove(m_NLayersForVertexSearch, 0);

      int Distance;
      for (Iterator2 = ReseList.begin(); Iterator2 != ReseList.end(); Iterator2++) {
        if ((*Iterator1) == (*Iterator2)) continue;
        if (IsInTracker(*Iterator2) == false) continue;

        Distance = m_Geometry->GetLayerDistance((*Iterator1), (*Iterator2));
        if (Distance > 0 && Distance < int(m_NLayersForVertexSearch)) NAbove[Distance]++;
        if (Distance < 0 && abs(Distance) < int(m_NLayersForVertexSearch)) NBelow[abs(Distance)]++;
        massert(Distance != 0); // In this case the algorithm is broken...
        mdebug<<"Distance "<<(*Iterator1)->GetID()<<" - "<<(*Iterator2)->GetID()<<": "<<Distance<<endl;
      }


      // Under the following conditions we do have a pair:

      // Pair starting from top
      Vertex = 0;
      //int VertexDirection = 0;

      unsigned int NLGE = 0;    // number of layers hit
      unsigned int NLGETwo = 0; // number of layers with hits greater or equal 2
      unsigned int NLGEOne = 0; // number of layers with exact one hit

      mdebug<<"Search vertex ("<<(*Iterator1)->GetPosition().Z()<<"): Above: "
      <<NAbove[1]<<" "<<NAbove[2]<<" "<<NAbove[3]<<" "<<NAbove[4]<<" "<<NAbove[5]<<endl;
      mdebug<<"Search vertex ("<<(*Iterator1)->GetPosition().Z()<<"): Below: "
      <<NBelow[1]<<" "<<NBelow[2]<<" "<<NBelow[3]<<" "<<NBelow[4]<<" "<<NBelow[5]<<endl;

      // Check for vertex below
      if (NAbove[1] == 0) {
        NLGETwo = 0;
        NLGEOne = 1; // We count the first one too...
        NLGE = 1; // We count the first one two...

        for (unsigned int d = 1; d < m_NLayersForVertexSearch; ++d) {
          if (NBelow[d] >= 1) NLGE++;
          if (NBelow[d] == 1) NLGEOne++;
          if (NBelow[d] >= 2) NLGETwo++;
        }
        mdebug<<"Vertec statistics (max: "<<m_NLayersForVertexSearch<<"): layers=:"<<NLGE<<" with one hit: "<<NLGEOne<<" 2 or more: "<<NLGETwo<<endl;
        if (NLGE >= m_NLayersForVertexSearch-MaximumLayerJump && NLGETwo >= 2 && NLGEOne >= 1) {
          Vertex = (*Iterator1);
          //VertexDirection = -1;
        }
      }

      if (Vertex != 0) {
        if (List == 0) List = new MRawEventList();
        break; // Only take first right now
      }
    }

    if (Vertex != 0){
      MRESE* PrelVertex = Vertex;
      vector < Float_t > Chi2;
      vector < Int_t > ID;

      Int_t nplanes = 0;
      Double_t PlaneZ = (*ReseList.begin())->GetPosition().Z();

      // Pair starting from top

      int VertexDirection = 0;

      for (Iterator1 = ReseList.begin(); Iterator1 != ReseList.end() && nplanes<2; Iterator1++) {
        if (IsInTracker(*Iterator1) == false) continue;

        if (abs(PlaneZ - (*Iterator1)->GetPosition().Z())>  + 0.0001){
          nplanes++;
          PlaneZ = (*Iterator1)->GetPosition().Z();
        }
        Float_t ChiT1=0.;
        Float_t ChiT2=0.;

        vector < TMatrix > Track1;
        vector < TMatrix > Track2;

        Float_t E1=0.;
        Float_t E2=0.;

        MRESEList Chosen1;
        MRESEList Chosen2;

        RE->SetVertex(*Iterator1);

        tie(Track1, E1, ChiT1, Track2, E2, ChiT2, Chosen1, Chosen2) = SearchTracks(RE);

        if ((ChiT1==0. && ChiT2==0.) || !(Chosen1.GetNRESEs()>2 && Chosen2.GetNRESEs()>2)) continue;
        if((*Iterator1)==PrelVertex) Chi2.push_back(ChiT1+ChiT2-sfg_factor);
        else Chi2.push_back(ChiT1+ChiT2);
        ID.push_back((*Iterator1)->GetID());
      }

      for (UInt_t i = 0; i < Chi2.size(); i++) {
        if(Chi2[i]==*min_element(Chi2.begin(), Chi2.end())) {
          Vertex=RE->GetRESE(ID[i]);
        }
      }

      if (Vertex != 0) {
        if (List == 0) List = new MRawEventList();
        RE->SetVertex(Vertex);
        RE->SetVertexDirection(VertexDirection);
        MRERawEvent *New = RE->Duplicate();
        RE->SetVertex(0);
        List->AddRawEvent(New);
        mdebug<<"Search vertex: Found vertex: "<<Vertex->GetID()<<endl;
      }
      //}
    }
    return List;
  }


  ////////////////////////////////////////////////////////////////////////////////

  Float_t MERTrackKalman3D::Kalman(MRERawEvent* RE, Int_t n, Float_t en, MRESEList Previous, vector < TMatrix > &trk, MRESEList &Chosen)
  {

    vector < TMatrix > Ck;
    vector < TMatrix > Csmoothk;
    vector < TMatrix > Fk;
    vector < TMatrix > FTk;
    vector < TMatrix > Cprojk;
    vector < TMatrix > CprojInvk;
    //vector < TMatrix > Qk;
    vector < TMatrix > xk;
    vector < TMatrix > xprojk;
    vector < TMatrix > xsmoothk;

    Float_t _planedistance = 0.;
    Float_t height = m_heightX0*9.37;

    TMatrix C(4,4);
    TMatrix C_proj(4,4);
    Cprojk.push_back(C_proj);
    CprojInvk.push_back(C_proj);
    TMatrix F(4,4);
    TMatrix Q(4,4);

    Chosen.RemoveAllRESEs();
    Chosen.Compress();

    Float_t chi2 = 0.;
    MRESE vertex = RE->GetVertex();
    Chosen.AddRESE(RE->GetVertex());
    //Filtering
    Float_t t;// tan of the RMS of the scattering angle

    Float_t V = (m_sigma*m_sigma);  //Weight Matrix

    TMatrix G(3,3);
    G(0,0)=1/(m_sigma*m_sigma);
    G(0,1)=0.;
    G(0,2)=0.;
    G(1,0)=0.;
    G(1,1)=0.;
    G(1,2)=0.;
    G(2,0)=0.;
    G(2,1)=0.;
    G(2,2)=1/(m_sigma*m_sigma);


    C(0,0)=m_sigma*m_sigma;
    C(0,1)=0.;
    C(1,0)=0.;
    C(0,2)=0.;
    C(0,3)=0.;
    C(1,1)=w;
    C(1,2)=0.;
    C(1,3)=0.;
    C(2,0)=0.;
    C(2,1)=0.;
    C(3,0)=0.;
    C(3,1)=0.;
    C(2,2)=m_sigma*m_sigma;
    C(3,3)=w;
    C(3,2)=0.;
    C(2,3)=0.;

    Ck.push_back(C);

    TMatrix x(4,1);
    x(0,0)=vertex.GetPositionX();
    x(1,0)=0.;
    x(2,0)=vertex.GetPositionY();
    x(3,0)=0.;

    xk.push_back(x);

    vector < TMatrix > mk;

    TMatrix m_temp(3,1);
    m_temp(0,0) = x(0,0);
    m_temp(1,0) = 0.;
    m_temp(2,0) = x(2,0);

    mk.push_back(m_temp);

    TMatrix x_proj(4,1);
    xprojk.push_back(x_proj);

    Int_t i=0;

    MRESEList List;
    List.RemoveAllRESEs();
    List.Compress();

    // Create a list of RESEs sorted by depth in tracker
    vector<MRESE*> ReseList;
    for (int h = 0; h < RE->GetNRESEs(); h++) {
      ReseList.push_back(RE->GetRESEAt(h));
    }
    sort(ReseList.begin(), ReseList.end(), CompareRESEByZ());

    //mdebug<<"RESE's sorted by depth: "<<endl;
    vector<MRESE*>::iterator Iterator1;
    for (Iterator1 = ReseList.begin(); Iterator1 != ReseList.end(); Iterator1++) {
      if (IsInTracker(*Iterator1) == false) continue;
      List.AddRESE(*Iterator1);
    }
    ReseList.clear();

    Float_t PlaneZ = vertex.GetPosition().Z();

    for (int r = 0; r < List.GetNRESEs() && i<n-1; r++) {

      if (PlaneZ <= (List.GetRESEAt(r)->GetPosition().Z() + 0.0001)) continue;

      PlaneZ = List.GetRESEAt(r)->GetPosition().Z();

      for (int k = 0; k < r; k++) {
        if(List.GetRESEAt(k)->GetPosition().Z()-PlaneZ < 0.0001) break;
        _planedistance = List.GetRESEAt(k)->GetPosition().Z()-PlaneZ;
        if(k==r-1 && _planedistance>1.5) {
          return chi2;
        }
      }

      if(_planedistance>1.5) return chi2; // we accept only hits on subsequent planes
      t=(13.6/en)*sqrt(m_heightX0)*(1+0.038*log(m_heightX0));

      F(0,0)=1.;
      F(0,1)=_planedistance;
      F(0,2)=0.;
      F(0,3)=0.;
      F(1,0)=0.;
      F(1,1)=1.;
      F(1,2)=0.;
      F(1,3)=0.;
      F(2,0)=0.;
      F(2,1)=0.;
      F(2,2)=1.;
      F(2,3)=_planedistance;
      F(3,0)=0.;
      F(3,1)=0.;
      F(3,2)=0.;
      F(3,3)=1.;

      Fk.push_back(F);

      Q(0,0)=(t*t*height*height)/3.;
      Q(2,2)=(t*t*height*height)/3.;
      Q(0,1)=t*t*height/2.;
      Q(1,0)=t*t*height/2.;
      Q(2,3)=t*t*height/2.;
      Q(3,2)=t*t*height/2.;
      Q(1,1)=t*t;
      Q(3,3)=t*t;
      Q(0,2)=0.;
      Q(0,3)=0.;
      Q(1,2)=0.;
      Q(1,3)=0.;
      Q(2,0)=0.;
      Q(2,1)=0.;
      Q(3,0)=0.;
      Q(3,1)=0.;

      //Qk.push_back(Q);

      TMatrix H(3,4);
      H(0,0)=1.;
      H(0,1)=0.;
      H(0,2)=0.;
      H(0,3)=0.;
      H(1,0)=0.;
      H(1,1)=0.;
      H(1,2)=0.;
      H(1,3)=0.;
      H(2,0)=0.;
      H(2,1)=0.;
      H(2,2)=1.;
      H(2,3)=0.;

      TMatrix HT=H;
      HT.T();

      x_proj = (F*x);

      xprojk.push_back(x_proj);

      TMatrix FT = F;
      FT.T();
      FTk.push_back(FT);

      C_proj = F*C*FT + Q;

      Cprojk.push_back(C_proj);

      TMatrix C_projInv=C_proj;

      C_projInv.Invert();

      CprojInvk.push_back(C_projInv);

      C = C_projInv;

      C+=HT*G*H;

      C.Invert();

      Ck.push_back(C);

      vector < Float_t > weight;

      for (int k = 0; k < List.GetNRESEs(); k++) {

        if(abs(List.GetRESEAt(k)->GetPosition().Z()-PlaneZ) > 0.0001) continue;

        Bool_t previous = false;

        for (int j = 0; j < Previous.GetNRESEs(); j++) {
          if(abs(Previous.GetRESEAt(j)->GetPosition().Z()-PlaneZ)  > 0.0001) continue;

          if(Previous.GetRESEAt(j)->ComputeMinDistance(List.GetRESEAt(k)) <0.001) previous = true;
        }
        Float_t distx = (List.GetRESEAt(k)->GetPosition().X()-x_proj(0,0));
        Float_t disty = (List.GetRESEAt(k)->GetPosition().Y()-x_proj(2,0));
        Float_t dist = distx*distx + disty * disty;
        if (previous) weight.push_back(sfg_factor+sqrt(dist/C_proj(0,0)));
        else weight.push_back(sqrt(dist/C_proj(0,0)));
      }

      Int_t kmin=0;

      for(int kk=0; kk<(int) weight.size(); kk++){
        if(weight[kk]==*min_element(weight.begin(), weight.end())) kmin=kk;
        if(*min_element(weight.begin(), weight.end())>max_offset && i>2) kmin=-1;
      }

      if (kmin>=0)
      {
        Int_t kind=0;
        for (int kk = 0; kk < List.GetNRESEs(); kk++) {
          if(abs(List.GetRESEAt(kk)->GetPosition().Z()-PlaneZ) > 0.0001) continue;
          else {
            if (kind==kmin){
              m_temp(0,0)=List.GetRESEAt(kk)->GetPosition().X();
              m_temp(1,0)=0.;
              m_temp(2,0)=List.GetRESEAt(kk)->GetPosition().Y();

              Chosen.AddRESE(List.GetRESEAt(kk));
              break;
            }
            kind ++;
          }
        }
      } else {
        m_temp(0,0)=x_proj(0,0);
        m_temp(1,0)=0.;
        m_temp(2,0)=x_proj(2,0);
      }

      mk.push_back(m_temp);
      x = C * (C_projInv*x_proj+HT*G*m_temp);
      xk.push_back(x);

      i++;
    }


    if ((Int_t)xk.size()< n) n = xk.size();

    //Smoothing

    TMatrix xsmooth(4,1);

    xsmoothk.resize(xk.size(), xsmooth);
    xsmoothk[xsmoothk.size()-1] = xk[xk.size()-1];
    trk.push_back(xsmoothk[xsmoothk.size()-1]);

    TMatrix Csmooth(4,4);
    Csmoothk.resize(Ck.size(), Csmooth);
    Csmoothk[Ck.size()-1] = C;


    for (Int_t i = n-2; i >= 0; i--) {

      TMatrix A(4,4);

      A = Ck[i]*FTk[i]*CprojInvk[i+1];

      TMatrix AT = A;
      AT.T();

      xsmoothk[i] = xk[i] + A * (xsmoothk[i+1] - xprojk[i+1]);

      trk.push_back(xsmoothk[i]);

      Csmooth = Ck[i] + A*(Csmoothk[i+1] - Cprojk[i+1])*AT;
      Csmoothk[i] = Csmooth;

      Float_t rkx = mk[i](0,0) - xsmoothk[i](0,0);
      Float_t rky = mk[i](2,0) - xsmoothk[i](2,0);
      Float_t rk = sqrt(rkx*rkx+rky*rky);

      Float_t Rkx = V -  Csmoothk[i](0,0);
      Float_t Rky = V -  Csmoothk[i](2,2);
      Float_t Rk = sqrt(Rkx*Rkx+Rky*Rky);

      chi2 += rk * rk / Rk;
    }

    chi2 /= n-3;
    reverse(trk.begin(),trk.end());

    return chi2;
  }

  tuple<vector < TMatrix >, Float_t, Float_t, vector < TMatrix >, Float_t, Float_t, MRESEList, MRESEList> MERTrackKalman3D::SearchTracks(MRERawEvent* RE)
  {
    MRESEList PrevChoice;
    Float_t en = start_nrg_fit;
    Float_t nrg_est1;
    Float_t nrg_est2;
    Float_t chi1=0.;
    Float_t chi2=0.;
    Int_t loop1=0;
    Bool_t stop=false;
    Int_t n = n_fit;

    MRESEList Chosen1;
    MRESEList Chosen2;

    vector < TMatrix > trk1;
    vector < TMatrix > trk2;
    mdebug<<"First track: "<<endl;
    PrevChoice.RemoveAllRESEs();
    PrevChoice.Compress();

    //track 1
    while(!stop){

      chi1 = Kalman(RE, n, en, PrevChoice, trk1, Chosen1);

      if (trk1.size()<3){
        return std::make_tuple(trk1, 0, 0, trk2, 0, 0, Chosen1, Chosen2);
      }

      nrg_est1=0.;

      //Energy Estimation
      //nrg_est1 = MultipleScattering(trk1);
      nrg_est1 = EnergyEst3D(trk1);

      if (nrg_est1> 50000.) nrg_est1=50000.;
      if ((nrg_est1/en)<1/3){
        //loop1++;
        en=en/3.;
        trk1.clear();
      } else if ((nrg_est1/en)>3. && loop1<9){
        loop1++;
        en=en*3.;
        trk1.clear();
      } else stop=true;
    }
    for (size_t i = 0; i < trk1.size(); i++) {
      mdebug<<"Track1 after Kalman "<<trk1.size()<<" "<<trk1[i](0,0)<<" "<<trk1[i](1,0)<<endl;
    }
    stop=false;
    en = start_nrg_fit;
    Int_t loop2=0;
    mdebug<<"Second track: "<<endl;

    while(!stop){

      chi2 = Kalman(RE, n, en, Chosen1, trk2, Chosen2);

      if (trk2.size()<3){
        return std::make_tuple(trk1, 0, 0, trk2, 0, 0, Chosen1, Chosen2);
      }


      nrg_est2=0;

      //Energy Estimation
      //nrg_est2 = MultipleScattering(trk2);
      nrg_est2 = EnergyEst3D(trk2);

      if (nrg_est2> 50000.) nrg_est2=50000.;
      if ((nrg_est2/en)<1/3){
        //loop1++;
        en=en/3.;
        trk2.clear();
      } else if ((nrg_est2/en)>3. && loop2<9){
        loop2++;
        en=en*3.;
        trk2.clear();
      } else stop=true;
    }

    for (size_t i = 0; i < trk2.size(); i++) {
      mdebug<<"Track2 after Kalman "<<trk2.size()<<" "<<trk2[i](0,0)<<" "<<trk2[i](1,0)<<endl;
    }
    return std::make_tuple(trk1, nrg_est1, chi1, trk2, nrg_est2, chi2, Chosen1, Chosen2);
  }

  Float_t MERTrackKalman3D::MultipleScattering(vector < TMatrix > trk)
  {
    Int_t n = trk.size();
    TVector3 vi0(1,0,0);
    TVector3 vi1(1,0,0);

    Float_t E = 0;

    for (Int_t i = 0; i < n-1; i++) {
      if(trk[i](3,0) == 0) vi0.SetTheta(0.);
      else vi0.SetTheta(atan(trk[i](3,0)/sin(atan2(trk[i](3,0),trk[i](1,0)))));

      vi0.SetPhi(atan2(trk[i](3,0),trk[i](1,0)));

      if (i != 0) {
        Float_t angle = vi0.Angle(vi1);
        E = (13.6/angle)*sqrt(m_heightX0/cos(vi0.Theta()))*(1+0.038*log(m_heightX0/cos(vi0.Theta())));
        //break;
      }
      vi1=vi0;
    }

    return E;
  }

  Float_t MERTrackKalman3D::MultipleScattering3D(vector < TMatrix > trk, Double_t e0)
  {
    Double_t epsilon = 1.5;
    Int_t n = trk.size();
    TVector3 vi0(1,0,0);
    TVector3 vi1(1,0,0);

    vector< Float_t > theta(n);
    vector< Float_t > thickness;
    vector< Float_t > en;
    vector< Float_t > angle;

    for (Int_t i = 0; i < n-1; i++) {

      if(trk[i](3,0) == 0) vi0.SetTheta(0.);
      else vi0.SetTheta(atan(trk[i](3,0)/sin(atan2(trk[i](3,0),trk[i](1,0)))));
      vi0.SetPhi(atan2(trk[i](3,0),trk[i](1,0)));

      theta[i] = vi0.Theta();

      if (i != 0) {
        angle.push_back(vi0.Angle(vi1));

      }
      vi1=vi0;
    }
    Double_t sigma = 0;
    Double_t chi = 0;

    for (Int_t i = 0; i < n-1; i++) {
      thickness.push_back(m_heightX0/cos(theta[i]));
      if ((e0 - (i + 1) * epsilon * cos(theta[i])) > 1)
      en.push_back(e0 - (i + 1) * epsilon * cos( theta[i]));
      else
      en.push_back(1);
    }

    for (Int_t i = 0; i < n - 1; i++) {
      sigma = 13.6 / en[i] * sqrt(2*thickness[i]);
      chi += (angle[i] / sigma) * (angle[i] / sigma);
    }
    chi /= n - 1;

    return chi;
  }

  Float_t MERTrackKalman3D::EnergyEst3D(vector < TMatrix > trk)
  {
    Int_t jjmax = 13;
    Int_t jjmin = 1;
    Int_t jmin = 0;
    Double_t ejmin = pow(2., jjmax);
    Double_t ej;

    Double_t absm[5];
    for (Int_t i = jjmax; i >= jjmin; i--) {
      Double_t m[5] = {10000, 10000, 10000, 10000, 10000};

      for (Int_t j = 0; j < 5; j++) {
        ej = ejmin - pow(2., i) + j * pow(2., (i - 1));
        m[j] = MultipleScattering3D(trk, ej);
        absm[j] = abs(m[j] - 1);
      }

      for (Int_t j = 0; j < 5; j++) {
        if (absm[j] == *min_element(absm, absm+5)) jmin = j;
      }
      if (ejmin > (pow(2., i) - jmin * pow(2., (i - 1)))) {
        ejmin = ejmin - pow(2., i) + jmin * pow(2., (i - 1));
      } else {
        ejmin = 0;
      }

    }
    Double_t energy = ejmin;
    return energy;
  }

  void MERTrackKalman3D::TrackPairs(MRERawEvent* RE)
  {
    if (RE->GetVertex() == 0) {
      Fatal("void MRERawEvent::TrackPairs()",
      "No vertex found!");
    }

    Float_t ChiT1=0.;
    Float_t ChiT2=0.;

    vector < TMatrix > Track1;
    vector < TMatrix > Track2;

    Float_t E1=0.;
    Float_t E2=0.;
    Float_t Egamma=0.;


    MRESEList Chosen1;
    MRESEList Chosen2;

    MRETrack* Electron = new MRETrack();
    MRETrack* Positron = new MRETrack();

    // Now create the vertex:
    Electron->AddRESE(RE->GetVertex());
    Electron->SetStartPoint(RE->GetVertex());
    Positron->AddRESE(RE->GetVertex());
    Positron->SetStartPoint(RE->GetVertex());

    mdebug<<"Vertex : "<<" "<<RE->GetVertex()->GetPosition().X()<<" "<<RE->GetVertex()->GetPosition().Y()<<" "<<RE->GetVertex()->GetPosition().Z()<<endl;

    tie(Track1, E1, ChiT1, Track2, E2, ChiT2, Chosen1, Chosen2) = SearchTracks(RE);

    if ((ChiT1>500. && ChiT2>500.) || isinf(ChiT1+ChiT2) || (ChiT1==0. && ChiT2==0.) || !(Chosen1.GetNRESEs()>2 && Chosen2.GetNRESEs()>2)) return;

    for (int r = 0; r < Chosen1.GetNRESEs(); r++) {
      Electron->AddRESE(Chosen1.GetRESEAt(r));

    }

    for (int r = 0; r < Chosen2.GetNRESEs(); r++) {
      Positron->AddRESE(Chosen2.GetRESEAt(r));
    }

    Egamma= E1+E2;

    Float_t theta1x=Track1[0](1,0);
    Float_t theta1y=Track1[0](3,0);

    Float_t theta2x=Track2[0](1,0);
    Float_t theta2y=Track2[0](3,0);

    Float_t phi1 = atan2(theta1y,theta1x);
    Float_t theta1 = atan(theta1y/sin(phi1));
    if(theta1y == 0) theta1 =0.;
    else theta1 = atan(theta1y/sin(phi1));
    Float_t phi2 = atan2(theta2y,theta2x);
    Float_t theta2 = atan(theta2y/sin(phi2));
    if(theta2y == 0) theta2 =0.;
    else theta2 = atan(theta2y/sin(phi2));

    Float_t x1= sin(theta1)*cos(phi1);
    Float_t y1= sin(theta1)*sin(phi1);
    Float_t z1= cos(theta1);

    Float_t x2= sin(theta2)*cos(phi2);
    Float_t y2= sin(theta2)*sin(phi2);
    Float_t z2= cos(theta2);

    Float_t x= (x1*E1+x2*E2)/Egamma;
    Float_t y= (y1*E1+y2*E2)/Egamma;
    Float_t z= (z1*E1+z2*E2)/Egamma;

    Float_t mag=sqrt(x*x+y*y+z*z);
    Float_t mag1=sqrt(x1*x1+y1*y1+z1*z1);
    Float_t mag2=sqrt(x2*x2+y2*y2+z2*z2);

    MVector direction(x/mag,y/mag,-z/mag);
    MVector direction1(-x1/mag1,-y1/mag1,-z1/mag1);
    MVector direction2(-x2/mag2,-y2/mag2,-z2/mag2);

    Electron->SetFixedDirection(direction1);
    Positron->SetFixedDirection(direction2);

    Electron->SetEnergy(E1*1000); //So far the energy was calculated in MeV, MEGAlib works in keV
    Positron->SetEnergy(E2*1000);

    mdebug<<"Energy of the two particles: "<<E1<<" "<<E2<<endl;
    mdebug<<"Reconstructed theta, phi and energy of the gamma: "<<c_Deg*direction.Theta()<<" "<<c_Deg*direction.Phi()<<" "<<Egamma<<endl;
    RE->SetElectronTrack(Electron);
    RE->SetPositronTrack(Positron);

    RE->SetTrackQualityFactor((ChiT1+ChiT2)/100);
    RE->SetPairQualityFactor((ChiT1+ChiT2)/100);
    RE->SetGoodEvent(true);
    RE->SetEventType(MRERawEvent::c_PairEvent);

  }
  // MERTrackKalman3D.cxx: the end...
  ////////////////////////////////////////////////////////////////////////////////
