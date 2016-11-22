
#include "MilliQPhysicsList.hh"
#include "FTFP_BERT.hh"
#include "G4StepLimiterPhysics.hh"


#include "MilliQGeneralPhysics.hh"
#include "MilliQEMPhysics.hh"
#include "MilliQMonopolePhysics.hh"
#include "MilliQMonopole.hh"
#include "MilliQMuonPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4EmStandardPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"
//#include "HadronPhysicsFTFP_BERT.hh"


#include "G4OpticalPhysics.hh"
#include "G4OpticalProcessIndex.hh"

#include "G4SystemOfUnits.hh"

MilliQPhysicsList::MilliQPhysicsList(const boost::property_tree::ptree pt) : G4VModularPhysicsList()
{
  // default cut value  (1.0mm)
  defaultCutValue = 1.0*mm;

  // General Physics
  //  RegisterPhysics( new MilliQGeneralPhysics("general") );

  // Radioactive decay
  RegisterPhysics(new G4RadioactiveDecayPhysics());

  MilliQMonopolePhysics* MonopolePhys = new MilliQMonopolePhysics(pt);
  //Monopole Physics
  RegisterPhysics( MonopolePhys );

  //These are the processes contained in the FTFP_BERT PhysList contained in the monopole.cc example
  G4int ver = 0;
  // EM physics
  RegisterPhysics( new G4EmStandardPhysics(ver));
  // Synchroton Radiation & GN Physics
  RegisterPhysics( new G4EmExtraPhysics(ver) );
  //Decays
  RegisterPhysics( new G4DecayPhysics(ver) );
  // Hadron Elastic scattering
  RegisterPhysics( new G4HadronElasticPhysics(ver) );
  // Hadron Physics
  //    RegisterPhysics(  new HadronPhysicsFTFP_BERT(ver));
  // Stopping Physics
  RegisterPhysics( new G4StoppingPhysics(ver) );
  // Ion Physics
  RegisterPhysics( new G4IonPhysics(ver));
  // Neutron tracking cut
  RegisterPhysics( new G4NeutronTrackingCut(ver));

  // EM Physics
  RegisterPhysics( new MilliQEMPhysics("standard EM"));

  // Muon Physics
  RegisterPhysics( new MilliQMuonPhysics("muon"));

  // Optical Physics
  G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();
  RegisterPhysics( opticalPhysics );


  opticalPhysics->SetScintillationYieldFactor(1.0);
  opticalPhysics->SetScintillationExcitationRatio(0.0);

  opticalPhysics->SetMaxNumPhotonsPerStep(100);
  opticalPhysics->SetMaxBetaChangePerStep(10.0);

  opticalPhysics->SetTrackSecondariesFirst(kCerenkov,true);
  opticalPhysics->SetTrackSecondariesFirst(kScintillation,true);

}


MilliQPhysicsList::~MilliQPhysicsList() {}


void MilliQPhysicsList::SetCuts(){
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets
  //   the default cut value for all particle types
  SetCutsWithDefault();
}
