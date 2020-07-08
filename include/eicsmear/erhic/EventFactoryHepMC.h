#ifndef INCLUDE_EICSMEAR_ERHIC_EVENTFACTORYHEPMC_H_
#define INCLUDE_EICSMEAR_ERHIC_EVENTFACTORYHEPMC_H_

#include "eicsmear/erhic/EventFactory.h"
#include "eicsmear/erhic/EventPythia.h"
#include "eicsmear/smear/EventSmear.h"

#include <HepMC3/ReaderAsciiHepMC2.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/GenParticle.h>

#include<map>

namespace HepMC3
{
  class ReaderAsciiHepMC2;
}

namespace erhic {

/* template<> */
/* EventFromAsciiFactory<erhic::EventHepMC> EventHepMC* Create(); */
/*   /\* /\\** *\/ */
/*   /\*  Perform end-of-event operations. *\/ */
/*   /\*  *\\/ *\/ */
/*   /\* Int_t FinishEvent() {return 0;} *\/ */

/*   /\** */
/*    Create a new particle from the last data read from the input stream. */
/*    *\/ */
/* template<> */
/* EventFromAsciiFactory<erhic::EventHepMC> bool AddParticle(); */

 /* int particleindex; */
 /* std::map < HepMC3::GenParticlePtr, int > hepmcp_index; */
  void HandleHepmcParticle( const HepMC3::GenParticlePtr& p, std::map < HepMC3::GenParticlePtr, int >& hepmcp_index, int& particleindex, std::unique_ptr<erhic::EventHepMC>& mEvent );
      

  // Warning: explicitly putting the erhic:: namespace before the class
  // name doesn't seen to work for template classes.
  //   ClassDef(EventFromAsciiFactory<erhic::EventHepMC>, 1)
}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_EVENTFACTORY_H_
