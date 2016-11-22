#include "MilliQRunAction.hh"
#include "MilliQRecorderBase.hh"
#include "MilliQDetectorConstruction.hh"
#include "MilliQPrimaryGeneratorAction.hh"
#include "MilliQRunActionMessenger.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleDefinition.hh"
#include "g4root.hh"
#include "G4ios.hh"

#include "G4EmCalculator.hh"

MilliQRunAction::MilliQRunAction(MilliQRecorderBase* r, const boost::property_tree::ptree pt) : fRecorder(r), fMilliQRunActionMessenger(0) {

  fDetector = new MilliQDetectorConstruction(pt);
  fKinematic = new MilliQPrimaryGeneratorAction(pt);

  fVerboseLevel = 0;
  fMilliQRunActionMessenger = new MilliQRunActionMessenger(this);

}

MilliQRunAction::~MilliQRunAction() {
  delete G4AnalysisManager::Instance();
  delete fMilliQRunActionMessenger;
  delete dataFormat;
}

void MilliQRunAction::BeginOfRunAction(const G4Run* aRun) {

  // Book histograms, ntuple
  dataFormat = new MilliQDataFormat();
  dataFormat->CreateDataMembers();


  if(fRecorder) fRecorder->RecordBeginOfRun(aRun);

  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  // Get analysis manager
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  // Open an output file
  // The default file name is set in B5RunAction::B5RunAction(),
  // it can be overwritten in a macro
  analysisManager->OpenFile();

}

void MilliQRunAction::EndOfRunAction(const G4Run* aRun) {

  if(fRecorder) fRecorder->RecordEndOfRun(aRun);

  // save histograms & ntuple
  //
  G4AnalysisManager* analysisManager = G4AnalysisManager::Instance();

  ////////////////////////////

  G4int nEvents = aRun->GetNumberOfEvent();
  if(nEvents == 0) {
    return;
  }

  if (fDetector->GetAlternateGeometry() == 1) {
    //run conditions
    //
    fDetector->DefineMaterials();
    const G4Material* material = fDetector->GetScintMaterial();
    G4String matName = material->GetName();

    const G4ParticleDefinition* particle = fKinematic->GetParticleGun()->GetParticleDefinition();
    G4String particleName = particle->GetParticleName();

    G4double ekin[100], dedx[100];
    G4EmCalculator calc;
    //  calc.SetVerbose(0);

    //
    G4int i;
    for(i = 0; i < 100; ++i) {
      ekin[i] = 0.2 * std::pow(10., 0.07 * G4double(i)) * std::pow(10., G4double(3)) * MeV;
      dedx[i] = calc.ComputeElectronicDEDX(ekin[i], particleName, matName);
    }
    G4cout << ekin[0] << G4endl;
    for(i = 0; i < 100; i++) G4cout << " E(MeV)= " << ekin[i] << " dedx= " << dedx[i] << G4endl;

    G4cout << "### End of dEdX table" << G4endl;

    for(i = 0; i < 100; ++i) {
      analysisManager->FillNtupleDColumn(MilliQDataFormat::kDEDX, 0, ekin[i]);
      analysisManager->FillNtupleDColumn(MilliQDataFormat::kDEDX, 1, dedx[i]);
      analysisManager->AddNtupleRow(MilliQDataFormat::kDEDX);
    }
  }
  //	CLHEP::HepRandom::showEngineStatus();


  /////////////////////////
  analysisManager->Write();
  analysisManager->CloseFile();

}
