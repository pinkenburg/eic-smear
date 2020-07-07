/**
   \file
   Implementation of class erhic::EventFactoryHepMC.
 
   \author    Chris Pinkenburg, Kolja Kauder
   \date      2020-07-07
   \copyright 2020 Brookhaven National Lab
*/

#include "eicsmear/erhic/EventFactoryHepMC.h"

#include <memory>
#include <stdexcept>
#include <string>

#include <TClass.h>
#include <TProcessID.h>

#include "eicsmear/erhic/BeamParticles.h"
#include "eicsmear/erhic/EventPythia.h"
#include "eicsmear/erhic/EventHepMC.h"
#include "eicsmear/erhic/EventMilou.h"
#include "eicsmear/erhic/EventDjangoh.h"
#include "eicsmear/erhic/EventDpmjet.h"
#include "eicsmear/erhic/EventRapgap.h"
#include "eicsmear/erhic/EventPepsi.h"
#include "eicsmear/erhic/EventGmcTrans.h"
#include "eicsmear/erhic/EventSimple.h"
#include "eicsmear/erhic/EventSartre.h"
#include "eicsmear/functions.h"  // For getFirstNonBlank()
#include "eicsmear/erhic/Kinematics.h"
#include "eicsmear/erhic/ParticleIdentifier.h"
#include "eicsmear/erhic/ParticleMC.h"

#include <TVector3.h>
#include <TParticlePDG.h>
#include <TLorentzVector.h>
#include <TDatabasePDG.h>

#include<map>

using std::cout;
using std::cerr;
using std::endl;
using std::map;

namespace erhic {
  
  // Use this struct to automatically reset TProcessID object count.
  struct TProcessIdObjectCount {
    // Initialse object with current TProcessID object count.
    TProcessIdObjectCount() {
      count = TProcessID::GetObjectCount();
    }
    // Restore object count to the value at initialisation.
    // See example in $ROOTSYS/test/Event.cxx
    // To save space in the table keeping track of all referenced objects
    // we assume that our events do not address each other.
    ~TProcessIdObjectCount() {
      TProcessID::SetObjectCount(count);
    }
    int count;
  };

  EventFromAsciiFactory<erhic::EventHepMC>::EventFromAsciiFactory(std::istream& is):
    mInput(&is),
    mEvent(nullptr) 
  {
    adapter2 = std::make_shared<HepMC3::ReaderAsciiHepMC2>(is);
  }

  std::string EventFromAsciiFactory<erhic::EventHepMC>::EventName() const {
    return erhic::EventHepMC::Class()->GetName();
  }

  erhic::EventHepMC* EventFromAsciiFactory<erhic::EventHepMC>::Create()
  {
    TProcessIdObjectCount objectCount;
    mEvent.reset(new erhic::EventHepMC());
    if (!AddParticle()) {
      mEvent.reset(nullptr);
    }  // if
    
    return mEvent.release();
  }

  bool EventFromAsciiFactory<erhic::EventHepMC>::AddParticle() {
    try {
      if (mEvent.get()) {
        HepMC3::GenEvent evt(HepMC3::Units::GEV,HepMC3::Units::MM);
	adapter2->read_event(evt);
	if ( adapter2->failed() )    return false;
	
	// std::cout << evt.cross_section()->get_attempted_events() << std::endl;
	// std::cout << evt.cross_section()->get_accepted_events() << std::endl;
	// std::cout << evt.cross_section()->xsec() << std::endl;

	particleindex = 1;
	// Can't use GenParticle::children() because they don't have indices assigned yet
	// map each HepMC particle onto its corresponding particleindex
	hepmcp_index.clear();

	// start with the beam plus gamma* and scattered lepton
	// ParticleIdentifier expects
	//   beams.SetBeamLepton(particles.at(0)->Get4Vector());
	//   beams.SetBeamHadron(particles.at(1)->Get4Vector());
	//   beams.SetBoson(particles.at(2)->Get4Vector());
	//   beams.SetScatteredLepton(particles.at(3)->Get4Vector());
	// While we don't _have_ to use that class, it makes sense to follow the convention
	// and not reinvent the wheel
	// for (auto& p : evt.beams() ) {
	//   cout << p->pid() << endl;
	// }




	// cout << " ===== "  << endl;
	// for (auto& p : evt.vertices().at(0)->particles_out() ) {
	//   cout << p->pid() << endl;
	// }
	// cout << " - "  << endl;
	// for (auto& p : evt.vertices().at(1)->particles_out() ) {
	//   cout << p->pid() << endl;
	// }
	// cout << " - "  << endl;
	// for (auto& p : evt.vertices().at(2)->particles_out() ) {
	//   cout << p->pid() << endl;
	// }
	// cout << " - "  << endl;
	// for (auto& p : evt.vertices().at(3)->particles_out() ) {
	//   cout << p->pid() << endl;
	// }
	// cout << " - "  << endl;
	// for (auto& p : evt.vertices().at(4)->particles_out() ) {
	//   cout << p->pid() << endl;
	// }
	// cout << " - "  << endl;
	// for (auto& p : evt.vertices().at(5)->particles_out() ) {
	//   cout << p->pid() << endl;
	// }
	// cout << " =============== "  << endl;

	
	for (auto& v : evt.vertices() ){	  
	  for (auto& p : v->particles_out() ) {	    
	    HandleHepmcParticle( p, v );
	    // std::cout << particle.GetIndex()<<std::endl;
	  }
	}

	for (auto& v : evt.vertices() ){	  
	  for (auto& p : v->particles_out() ) {
	    for (auto& parent : p->parents() ) {
	      
	      if ( hepmcp_index[parent] == 0 ){
		// cerr << hepmcp_index[p] << "  " << parent->pid() << "  " << hepmcp_index[parent] << endl;
	      }
	      // explicitly assumes that particleid = Event entry # +1
	      // cout << hepmcp_index[parent] -1 << endl;
	      // cout << "  " << parent->pid() << endl;
	      // cout << "  " << hepmcp_index[parent] << endl;
	    }
	  }
	}
	// for ( auto& hepmcp : hepmcps ){
	//   cout << hepmcp->parents().size() << endl;
	// }
	// for ( auto& pmc_gp : PMC_HepMCGP_lookup ){
	//   auto& pmc = pmc_gp.first;
	//   auto& gp = pmc_gp.second;
	//   std::cout << "PMC index: " << pmc->GetIndex() << std::endl;
	//   std::cout << "gp " << gp << std::endl;
	//   // std::cout << "Parent size: " << gp->parents().size() << std::endl;
	//   // std::cout << "Children size: " << gp->children().size() << std::endl;
	// }



	//ParticleMCeA *particle = new ParticleMCeA(mLine);  // Throws if the string is bad
	//particle->SetEvent(mEvent.get());
	//mEvent->AddLast(particle);
	//delete particle;
      }  // if
      return true;
    }  // try
    catch(std::exception& error) {
      std::cerr << "Exception building particle: " << error.what() << std::endl;
      return false;
    }
  }

  void EventFromAsciiFactory<erhic::EventHepMC>::HandleHepmcParticle( const HepMC3::GenParticlePtr& p, const HepMC3::GenVertexPtr& v){
    ParticleMC particle;
    
    TVector3 vertex(v->position().x(),v->position().y(),v->position().z());
    particle.SetVertex(vertex);
    // takes care of  xv, yv, zv;
    
    TLorentzVector lovec(p->momentum().x(),p->momentum().y(),p->momentum().z(),p->momentum().e());
    particle.Set4Vector(lovec);
    // takes care of  E, px, py, pz, m, derived quantities
    // derived quantities: pt, p, rapidity, eta, theta, phi
    
    // fill up the missing parts
    particle.SetId( p->pid() );
    particle.SetStatus(p->status());
    
    // Index: Runs from 1 to N.
    particle.SetIndex ( particleindex );
    
    // - parent and daughter questions should return 0 if none are found.
    // - logic:
    // --- fill Index up in this round
    // --- build up information
    // --- make the parent/child connection in a second loop
    // (I think) we cannot rely on parents always being read before their children.
    // So the afterburner needs to do the following
    // loop again over GenParticle's, say gp
    //    find gp's corresponding Index
    //    find gp's parents (vector<GenParticlePtr>)
    //    foreach parent
    //        find corresponding ParticleMC pmc
    //        add gp's index to the list of pmc's children
    //    done
    // loop over pmc's and condense the list of children
    
    // remember this HepMC3::GenParticlePtr <-> index connection
    // hepmcps.push_back ( p );
    hepmcp_index[p] = particleindex;
    
    // PMC_HepMCGP_lookup.push_back( PMC_HepMCGP_pair ( &particle, p ) );
    particleindex++;
    
    // std::cout << p->parents().size()<<std::endl;
    // std::cout << p->children().size()<<std::endl;
    // std::cout << std::endl;
    
    //
    // I(-1)
    // , orig(-1)
    // , daughter(-1)
    // , ldaughter(-1)
    // , parentId(std::numeric_limits<Int_t>::min())
    
    // These should come from FinishEvent()
    // , z(0.)
    // , xFeynman(0.)
    // , thetaGamma(0.)
    // , ptVsGamma(0.)
    // , phiPrf(0.)
    
    particle.SetEvent(mEvent.get());
    mEvent->AddLast(&particle);
  }
    
}  // namespace erhic

namespace {

  // Need this to generate the CINT code for each version
  erhic::EventFromAsciiFactory<erhic::EventHepMC> eh;
}  // namespace
