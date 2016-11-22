#include "MilliQActionInitialization.hh"

#include "MilliQPrimaryGeneratorAction.hh"

#include "MilliQRunAction.hh"
#include "MilliQEventAction.hh"
#include "MilliQTrackingAction.hh"
#include "MilliQSteppingAction.hh"
#include "MilliQStackingAction.hh"
#include "MilliQSteppingVerbose.hh"

#include "MilliQRecorderBase.hh"

#include "G4Types.hh"

void MilliQActionInitialization::BuildForMaster() const {
  SetUserAction(new MilliQRunAction(fRecorder, fPTree));
}

void MilliQActionInitialization::Build() const {

  SetUserAction(new MilliQPrimaryGeneratorAction(fPTree));

  SetUserAction(new MilliQStackingAction());

  SetUserAction(new MilliQRunAction(fRecorder, fPTree));
  SetUserAction(new MilliQEventAction(fRecorder, fPTree));
  SetUserAction(new MilliQTrackingAction(fRecorder));
  SetUserAction(new MilliQSteppingAction(fRecorder, fAlternate));
}

G4VSteppingVerbose* MilliQActionInitialization::InitializeSteppingVerbose() const {
  return new MilliQSteppingVerbose();
}
