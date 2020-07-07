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


#include <HepMC3/ReaderAsciiHepMC2.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/GenParticle.h>

template<>
class EventFromAsciiFactory<erhic::EventHepMC> : public VirtualEventFactory {
 public:
  /**
   Constructor.
   */
  EventFromAsciiFactory() { }

  /**
   Destructor.
   */
  virtual ~EventFromAsciiFactory() { }

  /**
   Initialise the factory from an input stream. 
   */
  explicit EventFromAsciiFactory(std::istream& is);

  /**
   Returns a new event instance.
   */
  virtual erhic::EventHepMC* Create();

  /**
   Returns the name of the event class created by this factory.
   */
  virtual std::string EventName() const;

  std::istream* mInput;  //!
  std::string mLine;  //!
  std::unique_ptr<erhic::EventHepMC> mEvent;  //!

 protected:

  std::shared_ptr<HepMC3::ReaderAsciiHepMC2> adapter2;
  /**
   Returns true when an end-of-event marker is encountered in the input stream.
   */
  bool AtEndOfEvent() const {return false;}

  /**
   Perform end-of-event operations.
   */
  Int_t FinishEvent() {return 0;}

  /**
   Create a new particle from the last data read from the input stream.
   */
  bool AddParticle();

 private:
  int particleindex;
  std::map < HepMC3::GenParticlePtr, int > hepmcp_index;
  void HandleHepmcParticle( const HepMC3::GenParticlePtr& p, const HepMC3::GenVertexPtr& v);
      

  // Warning: explicitly putting the erhic:: namespace before the class
  // name doesn't seen to work for template classes.
  ClassDef(EventFromAsciiFactory<erhic::EventHepMC>, 1)
};
}  // namespace erhic

#endif  // INCLUDE_EICSMEAR_ERHIC_EVENTFACTORY_H_
