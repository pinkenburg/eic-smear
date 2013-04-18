//
// ParticleIdentifier.cxx
//
// Created by TB on 6/10/11.
// Copyright 2011 BNL. All rights reserved.
//

#include <algorithm>
#include <iostream>

#include "eicsmear/erhic/EventMC.h"
#include "eicsmear/erhic/ParticleIdentifier.h"

//	=============================================================================
// Constructor
//	=============================================================================
ParticleIdentifier::ParticleIdentifier(const int leptonBeamPdgCode)
: mChargedCurrent(false)
, mLeptonBeamPdgCode(leptonBeamPdgCode)
, mScatteredPdgCode(leptonBeamPdgCode)
{ }

//	=============================================================================
// Identify the lepton beam particle
//	=============================================================================
bool ParticleIdentifier::isBeamLepton(
   const erhic::VirtualParticle& particle) const {
   // Test against status 201 for SOPHIA events
   return (21 == particle.GetStatus() or 201 == particle.GetStatus()) and
      GetLeptonBeamPdgCode() == particle.Id() and
      particle.GetParentIndex() == 0;
}

//	=============================================================================
// Identify the scattered lepton
//	=============================================================================
bool ParticleIdentifier::isScatteredLepton(
      const erhic::VirtualParticle& particle) const {
   return 1 == particle.GetStatus() and mScatteredPdgCode == particle.Id();
}

//	=============================================================================
// Return true if the particle matches any of the following:
//  - Outgoing electron with KS 21 - we instead pick up the
//                                   repeat entry, when KS < 21.
//  - Repeat of incoming nucleon.
//  - Intermediate gluon.
//  - Intermediate quark.
//	=============================================================================
bool ParticleIdentifier::SkipParticle(
   const erhic::VirtualParticle& particle) const {
   int kI1 = particle.GetStatus();
   int pdgCode = particle.Id();
   int parent = particle.GetParentIndex();
   // Remove duplicate outgoing electron, info is picked up later when KS < 21
   if(21 == kI1 and pdgCode == GetLeptonBeamPdgCode() and parent == 1) {
      return true;
   }	//	if
   // Remove repeat of incoming nucelon
   if(21 == kI1 and (pdgCode == 2112 or pdgCode == 2212) and parent == 2) {
      return true;
   }	//	if
   // Remove intermediate gluons
   if(21 == kI1 and pdgCode == 21) {
      return true;
   }	//	if
   // Remove intermediate (anti)quarks
   if(21 == kI1 and ::abs(pdgCode) < 10) {
      return true;
   }	//	if
   return false;
}

//	=============================================================================
// Identify the virtual photon based on its PDG code and PYTHIA K(I,1)
//	=============================================================================
bool ParticleIdentifier::IsVirtualPhoton(
      const erhic::VirtualParticle& particle) const {
   const int pdg = abs(particle.Id());
   return pdg > 21 and pdg < 25 and 21 == particle.GetStatus();
}

//	=============================================================================
// Identify the nucleon beam particle
//	=============================================================================
bool ParticleIdentifier::isBeamNucleon(
   const erhic::VirtualParticle& particle) const {
   // Test against status 201 for SOPHIA events
   return (21  == particle.GetStatus() or 201 == particle.GetStatus()) and
      (2112 == particle.Id() or 2212 == particle.Id()) and
      particle.GetParentIndex() == 0;
}

//	=============================================================================
// beam electron (11) --> scattered nu_e (12)
// beam positron (-11) --> scattered nu_e_bar (-12)
//	=============================================================================
Int_t ParticleIdentifier::DetermineScatteredType(Int_t beamType) {
   if(not mChargedCurrent) {
      return beamType;
   } // if
   int sign = beamType / abs(beamType);
   return sign * (abs(beamType) + 1);
}

//	=============================================================================
//	=============================================================================
void ParticleIdentifier::SetLeptonBeamPdgCode(const int pdgCode) {
   mLeptonBeamPdgCode = pdgCode;
   if(mChargedCurrent) {
      mScatteredPdgCode = DetermineScatteredType(pdgCode);
   } // if
   else {
      mScatteredPdgCode = pdgCode;
   } // else
}

bool ParticleIdentifier::SetChargedCurrent(bool cc) {
   // If charged current status is changed we need to update
   // the expected scattered lepton type. But we need to call
   // DetermineScatteredType() *after* updating mChargedCurrent.
   bool changed = not mChargedCurrent and cc;
   mChargedCurrent = cc;
   if(changed) {
      mScatteredPdgCode = DetermineScatteredType(mLeptonBeamPdgCode);
   } // if
   return cc;
}

//	=============================================================================
// Identify the beams from an event and store their properties in a
// BeamParticles object.
// See BeamParticles.h for the quantities stored.
// Returns true if all beams are found, false if not.
//	=============================================================================
bool ParticleIdentifier::IdentifyBeams(const erhic::VirtualEvent& event,
                                       BeamParticles& beams) {
   beams.Reset();
   std::vector<const erhic::VirtualParticle*> particles;
   bool foundAll = IdentifyBeams(event, particles);
   if(particles.at(0)) {
      beams.SetBeamLepton(particles.at(0)->Get4Vector());
   } // if
   if(particles.at(1)) {
      beams.SetBeamHadron(particles.at(1)->Get4Vector());
   } // if
   if(particles.at(2)) {
      beams.SetBoson(particles.at(2)->Get4Vector());
   } // if
   if(particles.at(3)) {
      beams.SetScatteredLepton(particles.at(3)->Get4Vector());
   } // if
   return foundAll;
}

//	=============================================================================
//	=============================================================================
bool ParticleIdentifier::IdentifyBeams(const erhic::VirtualEvent& event,
                           std::vector<const erhic::VirtualParticle*>& beams) {
   // Initialise a vector with four NULL pointers,
   // one beam particle of interest.
   const erhic::VirtualParticle* const null(NULL);
   beams.assign(4, null);
   // Configure the particle finder.
   ParticleIdentifier finder;
   if(event.GetTrack(0)) {
      finder.SetLeptonBeamPdgCode(event.GetTrack(0)->Id());
   } // if
   // Count leptons so we don't overwrite the beam and scattered lepton with
   // subsequent leptons of the same type.
   int leptonCount(0);
   // Set to true once we find the first virtual photon, so we can skip
   // subsequent virtual photons.
   bool foundExchangeBoson(false);
   bool foundAll(false);
   for(unsigned n(0); n < event.GetNTracks(); ++n ) {
      const erhic::VirtualParticle* particle = event.GetTrack(n);
      if(not particle) {
         continue;
      } // if
      // Test for beam lepton/hadron, exchange boson and scattered lepton.
      if(finder.isBeamNucleon(*particle)) {
         beams.at(1) = particle;
      }	//	if...
      else if(finder.isBeamLepton(*particle) and 0 == leptonCount) {
         beams.at(0) = particle;
         ++leptonCount;
      }	//	else if...
      else if(finder.isScatteredLepton(*particle) and 1 == leptonCount) {
         beams.at(3) = particle;
         // Protect against additional KS == 1 leptons following this
         ++leptonCount;
      }	//	if
      else if(finder.IsVirtualPhoton(*particle) and not foundExchangeBoson) {
         beams.at(2) = particle;
         foundExchangeBoson = true;
         // Check for charged current events, in which the scattered lepton
         // ID will not be the same as the incident lepton ID.
         finder.SetChargedCurrent(abs(particle->Id()) == 24);
      }	//	if
      // Break out if we've found all four beams (should happen at/near the
      // start of the event record, so checking the other particles is a
      // waste of time).
      foundAll = std::find(beams.begin(), beams.end(), null) == beams.end();
      if(foundAll) {
         break;
      } // if
   } // for
   return foundAll;
}
