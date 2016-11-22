//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: MilliQEventAction.cc 68752 2013-04-05 10:23:47Z gcosmo $
//
/// \file optical/MilliQ/src/MilliQEventAction.cc
/// \brief Implementation of the MilliQEventAction class
//
//
#include "MilliQEventAction.hh"
#include "MilliQPMTHit.hh"
#include "MilliQScintHit.hh"
#include "MilliQUserEventInformation.hh"
#include "MilliQDetectorConstruction.hh"
#include "MilliQTrajectory.hh"
#include "MilliQRecorderBase.hh"
#include "MilliQAnalysis.hh"

#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "g4root.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4UImanager.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "vector"
#include <algorithm>
#include <numeric>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
MilliQEventAction::MilliQEventAction(MilliQRecorderBase* r, const boost::property_tree::ptree pt)
  : fRecorder(r), fSaveThreshold(0), fPMTCollID(-1), fPMTAllCollID(-1), fScintCollID(-1), fVerbose(0),
    fPMTThreshold(1), fForcedrawphotons(false), fForcenophotons(false), fPTree(pt)
{
  fEventMessenger = new MilliQEventMessenger(this);

  try {
    boost::property_tree::ini_parser::read_ini(fPTree.get<std::string>("Configuration.GeometryConfigFile"), fGeometryPTree);
  }
  catch(boost::property_tree::ptree_error &e) {
    G4ExceptionDescription msg;
    msg << G4endl << "Configuration file " << e.what() << G4endl;
    G4Exception("MilliQAnalysis::ReadConfiguration()", "MilliQAnalysis::ConfigFileReadError", FatalException, msg);
  }


  // Calculates the CAEN IDs for which we need to check if there has been activation
  NBlocks = fGeometryPTree.get<G4int>("DetectorGeometry.NBlocks_X")*fGeometryPTree.get<G4int>("DetectorGeometry.NBlocks_Y")*fGeometryPTree.get<G4int>("DetectorGeometry.NBlocks_Z");
  NStacks = fGeometryPTree.get<G4int>("DetectorGeometry.NStacks");


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MilliQEventAction::~MilliQEventAction(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MilliQEventAction::BeginOfEventAction(const G4Event* anEvent){

  //New event, add the user information object
  G4EventManager::GetEventManager()->SetUserInformation(new MilliQUserEventInformation);

  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  if(fPMTCollID < 0) fPMTCollID = SDman->GetCollectionID("pmtHitCollection");
  if(fPMTAllCollID < 0) fPMTAllCollID = SDman->GetCollectionID("pmtAllHitCollection");
  if(fScintCollID < 0) fScintCollID = SDman->GetCollectionID("scintCollection");

  if(fRecorder) fRecorder->RecordBeginOfEvent(anEvent);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MilliQEventAction::EndOfEventAction(const G4Event * anEvent){

  G4AnalysisManager * analysisManager = G4AnalysisManager::Instance();
  MilliQUserEventInformation * eventInformation = (MilliQUserEventInformation*)anEvent->GetUserInformation();
  G4TrajectoryContainer * trajectoryContainer=anEvent->GetTrajectoryContainer();

  G4int n_trajectories = 0;
  if(trajectoryContainer) n_trajectories = trajectoryContainer->entries();

  // extract the trajectories and draw them
  if(G4VVisManager::GetConcreteInstance()) {
    for(G4int i = 0; i < n_trajectories; i++) {
      MilliQTrajectory * trj = (MilliQTrajectory*)((*(anEvent->GetTrajectoryContainer()))[i]);
      if(trj->GetParticleName()=="opticalphoton"){
        trj->SetForceDrawTrajectory(fForcedrawphotons);
        trj->SetForceNoDrawTrajectory(fForcenophotons);
      }
      trj->DrawTrajectory();
    }
  }

  // Get all the hit collections stored in the event from G4SDManager

  MilliQScintHitsCollection * scintHC = 0;
  MilliQPMTHitsCollection * pmtHC = 0;
  MilliQPMTHitsCollection * pmtAllHC = 0;
  G4HCofThisEvent * hitsCE = anEvent->GetHCofThisEvent();

  if(!hitsCE) {
    G4ExceptionDescription msg;
    msg << "No hits collection of this event found." << G4endl;
    G4Exception("MilliQEventAction::EndOfEventAction()", "B5Code001", JustWarning, msg);
    return;
  }

  // Now get the different hit collections

  if(hitsCE) {
    if(fScintCollID >= 0) scintHC = (MilliQScintHitsCollection*)(hitsCE->GetHC(fScintCollID));
    if(fPMTCollID >= 0) pmtHC = (MilliQPMTHitsCollection*)(hitsCE->GetHC(fPMTCollID));
    if(fPMTAllCollID >= 0) pmtAllHC = (MilliQPMTHitsCollection*)(hitsCE->GetHC(fPMTAllCollID));
  }

  // Hits in scintillator
  if(scintHC) {

    int n_hit = scintHC->entries();
    G4ThreeVector eWeightPos(0.);
    G4double edep;
    G4double edepMax = 0;

    //gather info on hits in scintillator
    for(int i = 0; i < n_hit; i++) {
      edep = (*scintHC)[i]->GetEdep();
      eventInformation->IncEDep(edep); // track the total energy deposition for entire event
      eWeightPos += (*scintHC)[i]->GetPos() * edep; // calculate energy-weighted position of hits
      if(edep > edepMax) {
	      edepMax = edep;//store max energy deposit
	      G4ThreeVector posMax = (*scintHC)[i]->GetPos();
	      eventInformation->SetPosMax(posMax, edep); // track the position and energy of largest energy deposit
      }
    }

    if(eventInformation->GetEDep() == 0.) {
      if(fVerbose > 0) G4cout << "No hits in the scintillator this event." << G4endl;
    }
    else {
      // Finish calculation of energy weighted position
      eWeightPos /= eventInformation->GetEDep();
      eventInformation->SetEWeightPos(eWeightPos);
      if(fVerbose > 0) G4cout << "\tEnergy weighted position of hits in MilliQ : " << eWeightPos / mm << G4endl;
    }

  } // End scintillator hits collection

  // Hits in PMTs
  if(pmtHC) {
    G4ThreeVector reconPos(0., 0., 0.);
    G4int pmts = pmtHC->entries();

    //Gather info from all PMTs
    for(G4int i = 0; i < pmts; i++) {
      eventInformation->IncHitCount((*pmtHC)[i]->GetPhotonCount()); // track number of PMT hits
      reconPos += (*pmtHC)[i]->GetPMTPos() * (*pmtHC)[i]->GetPhotonCount(); // calculate NHit-weighted hit position

      if((*pmtHC)[i]->GetPhotonCount() >= fPMTThreshold) eventInformation->IncPMTSAboveThreshold();
      else (*pmtHC)[i]->SetDrawit(false); //wasnt above the threshold, turn it back off
    }

    // dont bother unless there were hits
    if(eventInformation->GetHitCount() > 0) {
      reconPos /= eventInformation->GetHitCount(); // NHit-weighted hit position
      if(fVerbose > 0) G4cout << "\tReconstructed position of hits in MilliQ : " << reconPos / mm << G4endl;
      eventInformation->SetReconPos(reconPos);
    }
    pmtHC->DrawAllHits();
  }

  std::vector< std::vector<G4double> >  pmtTime(NBlocks*NStacks); // we want the first scintillator hit, and the first pmt hit
  std::vector< std::vector<G4double> >  scintTime(NBlocks*NStacks); // we want the first scintillator hit, and the first pmt hit
  std::vector< std::vector<G4double> >  scintEnergy(NBlocks*NStacks); // we want the first scintillator hit, and the first pmt hit

  // durp: what is the difference between pmtHC and pmtAllHC? pmtHC is used for eventInformation max energy and positions of hits

  bool PrintStats = false;
  if(scintHC && pmtAllHC){

    // For all hits in PMTs, record the time of hit in pmtTimes(PMT position, hit number)
    for(G4int j = 0; j < pmtAllHC->entries(); j++) {
      if((*pmtAllHC)[j]->GetPMTNumber() > -1) { // It was hit!
	      PrintStats = true;
	      pmtTime[(*pmtAllHC)[j]->GetPMTNumber()].push_back((*pmtAllHC)[j]->GetTime()) ;
      }
    }

    // For all hits in scintillators, record the time and energy in scintTime/scintEnergy
    for(int i = 0;i < scintHC->entries(); i++) {
      if((*scintHC)[i]->GetCpNum() > -1) { // It was hit!
	     scintTime[(*scintHC)[i]->GetCpNum()].push_back((*scintHC)[i]->GetTime());
	     scintEnergy[(*scintHC)[i]->GetCpNum()].push_back((*scintHC)[i]->GetEdep());
      }
    }

  }

  // Order hit times within each block/module
  for(unsigned int i = 0; i < pmtTime.size(); i++) sort(pmtTime[i].begin(), pmtTime[i].end());
  for(unsigned int i = 0; i < scintTime.size(); i++) sort(scintTime[i].begin(), scintTime[i].end());



  MilliQAnalysis* mcpan = new MilliQAnalysis(pmtTime, scintTime, scintEnergy,fPTree);

  if(mcpan->IsActive()){
      // Output information to ROOT MilliQAll
	  ActivePMT = mcpan->GetActivePMT();
	  NumberPMTHits = mcpan->GetNumberPMTHits();
	  PmtAllHitTimes = mcpan->GetPmtAllHitTimes();
	  PmtMedianHitTimes = mcpan->GetPmtMedianHitTimes();
	  TotalEnergyDeposit = mcpan->GetTotalEnergyDeposit();
	  TimeOfFlight = mcpan->GetTimeOfFlight();
	  analysisManager->FillNtupleIColumn(MilliQDataFormat::kAll, 6, mcpan->GetFirstHitScintillator() );
	  analysisManager->FillNtupleIColumn(MilliQDataFormat::kAll, 7, eventInformation->GetPhotonCount_Scint() );
	  analysisManager->AddNtupleRow(MilliQDataFormat::kAll);

      // Output information to ROOT MilliQCAEN
	  TriggerCount = mcpan->GetTriggerCount();
	  TimeCount = mcpan->GetTimeCount();
	  EventID = mcpan->GetEventID();
	  TDC = mcpan->GetTDC();
	  WaveMax = mcpan->GetWaveMax();
	  WaveMin = mcpan->GetWaveMin();
	  ChargeIntegral = mcpan->GetChargeIntegral();
	  analysisManager->AddNtupleRow(MilliQDataFormat::kCAEN);



  }

  // Output All information

  delete(mcpan);




  if(fVerbose>0){//PrintStats){
    //End of event output. later to be controlled by a verbose level
    G4cout << "\tNumber of photons that hit PMTs in this event : "
           << eventInformation->GetHitCount() << G4endl;
    G4cout << "\tNumber of PMTs above threshold("<<fPMTThreshold<<") : "
           << eventInformation->GetPMTSAboveThreshold() << G4endl;
    G4cout << "\tNumber of photons produced by scintillation in this event : "
           << eventInformation->GetPhotonCount_Scint() << G4endl;
    G4cout << "\tNumber of photons produced by cerenkov in this event : "
           << eventInformation->GetPhotonCount_Ceren() << G4endl;
    G4cout << "\tNumber of photons absorbed (OpAbsorption) in this event : "
           << eventInformation->GetAbsorptionCount() << G4endl;
    G4cout << "\tNumber of photons absorbed at boundaries (OpBoundary) in "
           << "this event : " << eventInformation->GetBoundaryAbsorptionCount()
           << G4endl;
    G4cout << "Unacounted for photons in this event : "
           << (eventInformation->GetPhotonCount_Scint() +
               eventInformation->GetPhotonCount_Ceren() -
               eventInformation->GetAbsorptionCount() -
               eventInformation->GetHitCount() -
               eventInformation->GetBoundaryAbsorptionCount())
           << G4endl;
  }
  //If we have set the flag to save 'special' events, save here
  if(fSaveThreshold&&eventInformation->GetPhotonCount() <= fSaveThreshold)
    G4RunManager::GetRunManager()->rndmSaveThisEvent();

  if(fRecorder)fRecorder->RecordEndOfEvent(anEvent);

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MilliQEventAction::SetSaveThreshold(G4int save){
  /*Sets the save threshold for the random number seed. If the number of photons
    generated in an event is lower than this, then save the seed for this event
    in a file called run###evt###.rndm
  */
  fSaveThreshold=save;
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  G4RunManager::GetRunManager()->SetRandomNumberStoreDir("random/");
  //  G4UImanager::GetUIpointer()->ApplyCommand("/random/setSavingFlag 1");
}
