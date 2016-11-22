#ifndef MilliQAnalysis_H
#define MilliQAnalysis_H 1

#include "MilliQDataFormat.hh"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "MilliQWaveform.hh"



#include <vector>

class MilliQAnalysis {
 public:
  MilliQAnalysis(std::vector< std::vector<G4double> > pmtTime,
		 std::vector< std::vector<G4double> > scintTime,
		 std::vector< std::vector<G4double> > scintEn,
		 const boost::property_tree::ptree pt);

  void NearestN();
  void CalculateUsefulInfo();
  void CAENID();
  void OnlineTrigger();
  void Waveform();
  bool IsActive() { return fIsActive; }

  // Output All info
  std::vector<G4int> GetActivePMT() { return ActivePMT; }
  std::vector<G4int> GetNumberPMTHits() { return NumberPMTHits; }
  std::vector<G4double> GetPmtMedianHitTimes() { return PmtMedianHitTimes; }
  std::vector<G4double> GetPmtAllHitTimes() { return PmtAllHitTimes; }
  std::vector<G4double> GetTotalEnergyDeposit() { return TotalEnergyDeposit; }
  std::vector<G4double> GetTimeOfFlight() { return TimeOfFlight; }
  G4int GetFirstHitScintillator() { return FirstHitScintillator; }

  // Output MilliQCAEN info
  std::vector<G4float> GetTriggerCount() { return TriggerCount;}
  std::vector<G4float> GetTimeCount() { return TimeCount;}
  std::vector<G4int> GetEventID() {return EventID;}
  std::vector<G4int> GetTDC() { return TDC;}
  std::vector<G4float> GetWaveMax() {return WaveMax;}
  std::vector<G4float> GetWaveMin() { return WaveMin;}
  std::vector<G4float> GetChargeIntegral() { return ChargeIntegral;}


  inline void SetVerbose(bool v) { fVerbose = v; }
  
 private:

  std::vector< std::vector<G4double> > fpmtTime;
  std::vector< std::vector<G4double> > fscintTime;
  std::vector< std::vector<G4double> > fscintEn;

  // Output MilliQAll info
  std::vector<G4int> ActivePMT;
  std::vector<G4int> NumberPMTHits;
  std::vector<G4double> PmtAllHitTimes;
  std::vector<G4double> PmtMedianHitTimes;
  std::vector<G4double> TotalEnergyDeposit;
  std::vector<G4double> TimeOfFlight;
  G4int FirstHitScintillator;

  // Output MilliQCAEN info
  std::vector<G4float> TriggerCount;
  std::vector<G4float> TimeCount;
  std::vector<G4int> EventID;
  std::vector<G4int> TDC;
  std::vector<G4float> WaveMax;
  std::vector<G4float> WaveMin;
  std::vector<G4float> ChargeIntegral;


  std::vector< std::vector<G4int> > fModules;
  G4ThreeVector fNblock;
  G4int fNstack;
  unsigned int fDetectorSize;
  bool fIsActive;
  G4double fcoincidenceThreshold;

  const boost::property_tree::ptree fPTree;
  boost::property_tree::ptree fGeometryPTree;
  boost::property_tree::ptree fParticlePTree;
  boost::property_tree::ptree fPMTPTree;


  bool fVerbose;
};

#endif
