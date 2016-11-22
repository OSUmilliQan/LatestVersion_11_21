#include "MilliQDataFormat.hh"
#include "G4RunManager.hh"

#include "G4ios.hh"

MilliQDataFormat::MilliQDataFormat() {

  fVerboseLevel = 0;

}

MilliQDataFormat::~MilliQDataFormat() {

}

void MilliQDataFormat::CreateAnalysisManager() {

  analysisManager = G4AnalysisManager::Instance();

  // Choice of analysis technology is done by including g4root.hh
  G4cout << "Using G4AnalysisManager type: " << analysisManager->GetType() << G4endl;

  // Default settings
  analysisManager->SetVerboseLevel(fVerboseLevel);
  analysisManager->SetFileName("MilliQ");

}

void MilliQDataFormat::CreateHistograms() {

  // Creating 1D histograms
  analysisManager->CreateH1("PMT", "PMT # Hits", 50, 0., 50); // h1 Id = 0

  // Creating 2D histograms
  // analysisManager->CreateH2("Chamber1 XY","Drift Chamber 1 X vs Y",           // h2 Id = 0
  //

}

void MilliQDataFormat::CreateNtuples() {

  const MilliQEventAction* constEventAction =
		  static_cast<const MilliQEventAction*>(G4RunManager::GetRunManager()->GetUserEventAction());
  MilliQEventAction* eventAction = const_cast<MilliQEventAction*>(constEventAction);

  // kAll
  analysisManager->CreateNtuple("MilliQAll", "All");
  analysisManager->CreateNtupleIColumn("ActivePMT", eventAction->GetActivePMT());						// column Id = 0
  analysisManager->CreateNtupleDColumn("PmtMedianHitTimes", eventAction->GetPmtMedianHitTimes());		// column Id = 1
  analysisManager->CreateNtupleDColumn("PmtAllHitTimes", eventAction->GetPmtAllHitTimes());				// column Id = 2
  analysisManager->CreateNtupleDColumn("TimeOfFlight", eventAction->GetTimeOfFlight());  				// column Id = 3
  analysisManager->CreateNtupleDColumn("TotalEnergyDeposit", eventAction->GetTotalEnergyDeposit()); 	// column Id = 4
  analysisManager->CreateNtupleIColumn("NumberPMTHits", eventAction->GetNumberPMTHits());				// column Id = 5
  analysisManager->CreateNtupleIColumn("FirstHitScintillator");											// column Id = 6
  analysisManager->CreateNtupleIColumn("PhotonCountAllScintillators");									// column Id = 7
//analysisManager->CreateNtupleDColumn("ScintillatorEnergyDeposits", eventAction->Get()); 				// column Id = 8
//analysisManager->CreateNtupleDColumn("ScintillatorTimes", eventAction->Get()); 						// column Id = 9
  analysisManager->FinishNtuple();


  // kCAEN
  analysisManager->CreateNtuple("MilliQCAEN", "All");
  analysisManager->CreateNtupleFColumn("TriggerCount", eventAction->GetTriggerCount()); // column Id = 0
  analysisManager->CreateNtupleFColumn("TimeCount", eventAction->GetTimeCount()); // column Id = 0
  analysisManager->CreateNtupleIColumn("EventID", eventAction->GetEventID()); // column Id = 0
  analysisManager->CreateNtupleIColumn("TDC", eventAction->GetTDC()); // column Id = 0
  analysisManager->CreateNtupleFColumn("WaveMax", eventAction->GetWaveMax()); // column Id = 0
  analysisManager->CreateNtupleFColumn("WaveMin", eventAction->GetWaveMin()); // column Id = 0
  analysisManager->CreateNtupleFColumn("ChargeIntegral", eventAction->GetChargeIntegral()); // column Id = 0
  analysisManager->FinishNtuple();







  // kDEDX
  analysisManager->CreateNtuple("MilliQDedx", "DEDX");
  analysisManager->CreateNtupleDColumn("EKinMeV"); // column Id = 0
  analysisManager->CreateNtupleDColumn("MeVpermm"); // column Id = 1
  analysisManager->FinishNtuple();







}
