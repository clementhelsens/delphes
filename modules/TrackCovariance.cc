/*
 *  Delphes: a framework for fast simulation of a generic collider experiment
 *  Copyright (C) 2020  Universite catholique de Louvain (UCLouvain), Belgium
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/** \class TrackCovariance
 *
 *  Smears track parameters according to appropriate covariance matrix.
 *
 *  \authors P. Demin - UCLouvain, Louvain-la-Neuve
 *           M. Selvaggi - CERN
 *
 */

//FIXME add reference to Bedeschi-code
//FIXME make sure about units of P, X
//FIXME fix pt > 200 GeV issue and angle > 6.41

#include "modules/TrackCovariance.h"

#include "classes/DelphesClasses.h"

#include "TrackCovariance/SolGeom.h"
#include "TrackCovariance/SolGridCov.h"
#include "TrackCovariance/ObsTrk.h"

#include "TLorentzVector.h"
#include "TMath.h"
#include "TObjArray.h"

#include <iostream>
#include <sstream>

using namespace std;

//------------------------------------------------------------------------------

TrackCovariance::TrackCovariance() :
  fGeometry(0), fCovariance(0), fItInputArray(0)
{
  fGeometry = new SolGeom();
  fCovariance = new SolGridCov();
}

//------------------------------------------------------------------------------

TrackCovariance::~TrackCovariance()
{
  if(fGeometry) delete fGeometry;
  if(fCovariance) delete fCovariance;
}

//------------------------------------------------------------------------------

void TrackCovariance::Init()
{
  fBz = GetDouble("Bz", 0.0);
  fGeometry->Read(GetString("DetectorGeometry", ""));

  fCovariance->Calc(fGeometry);

  // import input array

  fInputArray = ImportArray(GetString("InputArray", "TrackMerger/tracks"));
  fItInputArray = fInputArray->MakeIterator();

  // create output array

  fOutputArray = ExportArray(GetString("OutputArray", "tracks"));
}

//------------------------------------------------------------------------------

void TrackCovariance::Finish()
{
  if(fItInputArray) delete fItInputArray;
}

//------------------------------------------------------------------------------

void TrackCovariance::Process()
{
  Candidate *candidate, *mother;
  Double_t mass, p, pt, q, ct;
  Double_t dd0, ddz, dphi, dct, dp, dpt, dC;
  

  fItInputArray->Reset();
  while((candidate = static_cast<Candidate *>(fItInputArray->Next())))
  {
    const TLorentzVector &candidatePosition = candidate->InitialPosition;
    const TLorentzVector &candidateMomentum = candidate->Momentum;

    //std::cout << "  x  " << candidatePosition.Vect().X() << "  y  " << candidatePosition.Vect().Y() << "  z  " << candidatePosition.Vect().Z() << std::endl;

    mass = candidateMomentum.M();
    TVector3 poscorrect;
    poscorrect(0)=candidatePosition.Vect().X()/1000.;
    poscorrect(1)=candidatePosition.Vect().Y()/1000.;
    poscorrect(2)=candidatePosition.Vect().Z()/1000.;

    ObsTrk track(poscorrect, candidateMomentum.Vect(), candidate->Charge, fBz, fCovariance);

    mother    = candidate;
    candidate = static_cast<Candidate *>(candidate->Clone());

    candidate->Momentum.SetVectM(track.GetObsP(), mass);
    candidate->InitialPosition.SetXYZT(track.GetObsX().X(),track.GetObsX().Y(),track.GetObsX().Z(),candidatePosition.T());

    /*float p_pre = candidate->Momentum.P();
    float p_pre_px = candidate->Momentum.Px();
    float p_pre_py = candidate->Momentum.Py();
    float p_pre_pz = candidate->Momentum.Pz();
    float p_after = candidate->Momentum.P();
    float p_after_px = candidate->Momentum.Px();
    float p_after_py = candidate->Momentum.Py();
    float p_after_pz = candidate->Momentum.Pz();

    if (p_pre>2.0){
      std::cout << "p_pre " << p_pre << "  p_after  " << p_after << std::endl;
      std::cout << "p_pre_px " << p_pre_px << "  p_after_px  " << p_after_px << "  track obs x " << track.GetObsP()[0]<< std::endl;
      std::cout << "p_pre_py " << p_pre_py << "  p_after_py  " << p_after_py << "  track obs y " << track.GetObsP()[1]<<std::endl;
      std::cout << "p_pre_pz " << p_pre_pz << "  p_after_pz  " << p_after_pz << "  track obs z " << track.GetObsP()[2]<<std::endl;
      std::cout << "phi_pre  " << TMath::ATan(p_pre_py/p_pre_px) << "  phi_after  " << track.GetObsPar()[1] << std::endl;
      float lambda =  TMath::ATan(p_pre/candidate->Momentum.Pt());
      float theta = TMath::Pi()/2.-lambda;
      float theta_after = TMath::ATan(1./track.GetObsPar()[4]);
      std::cout << "theta_pre  " << theta << "  theta_after  " << theta_after << std::endl;
      }*/

    // save full covariance 5x5 matrix internally (D0, phi, Curvature, dz, ctg(theta))
    candidate->TrackCovariance = track.GetCov();

    pt = candidate->Momentum.Pt();
    p  = candidate->Momentum.P();
    q  = track.GetObsQ();
    ct = track.GetObsPar()[4];

    candidate->Xd = track.GetObsX().X();
    candidate->Yd = track.GetObsX().Y();
    candidate->Zd = track.GetObsX().Z();
    
    candidate->D0 = track.GetObsPar()[0];
    candidate->Phi = track.GetObsPar()[1];
    candidate->C = track.GetObsPar()[2];
    candidate->DZ = track.GetObsPar()[3];
    candidate->CtgTheta = track.GetObsPar()[4];
    candidate->P  = track.GetObsP().Mag();
    candidate->PT = pt;
    candidate->Charge = q;

    dd0       = TMath::Sqrt(track.GetCov()(0, 0)); 
    ddz       = TMath::Sqrt(track.GetCov()(3, 3)); 
    dphi      = TMath::Sqrt(track.GetCov()(1, 1)); 
    dct       = TMath::Sqrt(track.GetCov()(4, 4)); 
    dpt       = 2 * TMath::Sqrt( track.GetCov()(2, 2))*pt*pt / (0.2998*fBz);
    dp        = TMath::Sqrt((1.+ct*ct)*dpt*dpt + 4*pt*pt*ct*ct*dct*dct/(1.+ct*ct)/(1.+ct*ct));
    dC        = TMath::Sqrt(track.GetCov()(2, 2));

    candidate->ErrorD0 = dd0;
    candidate->ErrorDZ = ddz;
    candidate->ErrorP = dp;
    candidate->ErrorC = dC;
    candidate->ErrorCtgTheta = dct;
    candidate->ErrorPhi = dphi;
    candidate->ErrorPT = dpt;
    //candidate->TrackResolution = dpt / pt;
    candidate->TrackResolution = dp / p;


    candidate->AddCandidate(mother);

    fOutputArray->Add(candidate);
  }
}

//------------------------------------------------------------------------------
