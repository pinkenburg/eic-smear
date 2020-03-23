// Adapted from https://gitlab.com/eic/ejana/-/blob/master/src/plugins/reco/eic_smear/JleicSmear.cc
// Work in progress - no guarantees for correctness!

#include "eicsmear/erhic/VirtualParticle.h"
#include "eicsmear/smear/Acceptance.h"
#include "eicsmear/smear/Device.h"
#include "eicsmear/smear/Detector.h"
#include "eicsmear/smear/Smearer.h"
#include "eicsmear/smear/ParticleMCS.h"
#include "eicsmear/smear/PerfectID.h"
#include <eicsmear/smear/Smear.h>
#include <eicsmear/erhic/ParticleMC.h>
#include "Math/Vector4D.h"

using std::cout;
using std::endl;

double ThetaFromEta( const double eta ) {
  // The default value of -19 is used in the main eRHIC code,
  // so use that for consistency.
  if ( !isnan(eta) && !isinf(eta)   ) {
    return 2.0 * atan( exp( -eta ));
  }

  throw std::runtime_error("ThetaFromEta called with NaN or Inf");
  return -1;
}

double EtaFromTheta( const double theta ) {
  // The default value of -19 is used in the main eRHIC code,
  // so use that for consistency.
  double eta = -19.;
  if (theta > 0. && theta < TMath::Pi() && !TMath::IsNaN(theta)) {
    eta = -log(tan(theta / 2.));
  }
  return eta;
}


Smear::Detector BuildJLEIC() {

  gSystem->Load("libeicsmear");

  // EM Calorimeters
  // ---------------
  // eta = -3 --  -1.1
  // sigma^2 = (0.02*sqrt(E))^2 + (0.001*E)^2 + 0.005^2;
  Smear::Device emcalBck(Smear::kE, "sqrt( pow(0.02*sqrt(E),2) + pow(0.001*E,2) + pow(0.005,2) )" );
  auto thMaxBck = ThetaFromEta ( -3 );  auto thMinBck = ThetaFromEta ( -1.1 );
  //   if (thMinBck >  thMaxBck ) std::swap(thMinBck,thMaxBck);
  Smear::Acceptance::Zone emcalBckZone(thMinBck, thMaxBck);
  emcalBck.Accept.AddZone(emcalBckZone);
  emcalBck.Accept.SetGenre(Smear::kElectromagnetic);
    
  // eta = -1.1 --  3.5
  // sigma^2 = (0.1*sqrt(E))^2 + (0.01*E)^2;
  Smear::Device emcalMidFwd(Smear::kE, "sqrt( pow(0.1*sqrt(E),2) + pow(0.01*E,2) )" );
  auto thMaxMidFwd = ThetaFromEta ( -3 );  auto thMinMidFwd = ThetaFromEta ( -1.1 );
  Smear::Acceptance::Zone emcalMidFwdZone(thMinMidFwd, thMaxMidFwd);
  emcalMidFwd.Accept.AddZone(emcalMidFwdZone);
  emcalMidFwd.Accept.SetGenre(Smear::kElectromagnetic);

  // Hadronic Calorimeters
  // ---------------------
  Smear::Device hcal(Smear::kE, "sqrt(E)");
  Smear::Acceptance::Zone hcalZone; // Without arguments, accept 4pi
  hcal.Accept.AddZone(hcalZone);
  hcal.Accept.SetGenre(Smear::kHadronic);
  
  // Note: Original implementation has separate zdc for neutrons:
  //   // Neutron coming to Zero Degree calorimeter
  // if(abs(particle->pdg) == 2112 && p.Theta() < 0.01){
  //        // zdc
  //     new_e = gRandom->Gaus(p.E(), 1. * sqrt(p.E()));
  //     is_smeared_e = true;
  // }
  // But those should already be captured by the infinite hcal above.
  // Should it be doubly smeared?
    
  // TODO: separated smearing for e-endcap and h-endcap
    

  // Tracking
  // --------
  // Note: Original implementation has vertex smearing,
  // which (at least for now) is not supported here
  
  // eta = -3.5 --  3.5
  auto thMaxTrk = ThetaFromEta ( -3 );  auto thMinTrk = ThetaFromEta ( -1.1 );
  Smear::Acceptance::Zone trk(thMinTrk,thMaxTrk);
  Smear::Device TrackPt(Smear::kPt, "sqrt ( pow( 0.01 * pow(pT,2), 2) + pow(0.005 * pT,2))");
  Smear::Device TrackTheta(Smear::kTheta, "0.001");
  Smear::Device TrackPhi(Smear::kPhi, "0.001");

  TrackPt.Accept.AddZone(trk);
  TrackTheta.Accept.AddZone(trk);
  TrackPhi.Accept.AddZone(trk);
  TrackPt.Accept.SetGenre(Smear::kHadronic);
  TrackTheta.Accept.SetGenre(Smear::kHadronic);
  TrackPhi.Accept.SetGenre(Smear::kHadronic);
  TrackPt.Accept.SetCharge(Smear::kCharged);
  TrackTheta.Accept.SetCharge(Smear::kCharged);
  TrackPhi.Accept.SetCharge(Smear::kCharged);

  // TODO: separated smearing for e-endcap and h-endcap

  // Roman Pots, Theta < 0.01
  Smear::Acceptance::Zone trkRP(0,0.01);
  Smear::Device TrackPtRP(Smear::kPt, "0.02");
  TrackPtRP.Accept.AddZone(trkRP);
  TrackPtRP.Accept.SetCharge(Smear::kCharged);

  // D1, 0.01 < Theta < 0.05
  Smear::Acceptance::Zone trkD1(0.01,0.05);
  Smear::Device TrackPtD1(Smear::kPt, "0.01");
  TrackPtD1.Accept.AddZone(trkD1);
  TrackPtD1.Accept.SetCharge(Smear::kCharged);

  // Create a detector and add the devices
  Smear::Detector det;
  det.AddDevice(TrackPt);
  // det.AddDevice(TrackTheta);
  // det.AddDevice(TrackPhi);

  // det.AddDevice(TrackPtRP);
  // det.AddDevice(TrackPtD1);

  // det.AddDevice(emcalBck);
  // det.AddDevice(emcalMidFwd);
  // det.AddDevice(hcal);


  det.SetEventKinematicsCalculator("NM JB DA"); // The detector will calculate event kinematics from smeared values

  return det;
}
