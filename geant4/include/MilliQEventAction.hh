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
// $Id: MilliQEventAction.hh 68752 2013-04-05 10:23:47Z gcosmo $
//
/// \file optical/MilliQ/include/MilliQEventAction.hh
/// \brief Definition of the MilliQEventAction class
//

#ifndef MilliQEventAction_h
#define MilliQEventAction_h 1

#include "MilliQDataFormat.hh"
#include "MilliQEventMessenger.hh"
#include "G4UserEventAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

class G4Event;
class MilliQRecorderBase;

class MilliQEventAction : public G4UserEventAction {
 public:

  MilliQEventAction(MilliQRecorderBase*, const boost::property_tree::ptree pt);
  virtual ~MilliQEventAction();

  virtual void BeginOfEventAction(const G4Event*);
  virtual void EndOfEventAction(const G4Event*);

  void SetSaveThreshold(G4int);

  void SetEventVerbose(G4int v){ fVerbose = v;}

  void SetPMTThreshold(G4int t) { fPMTThreshold = t;}

  void SetForceDrawPhotons(G4bool b) { fForcedrawphotons = b;}
  void SetForceDrawNoPhotons(G4bool b) { fForcenophotons = b;}

  // Output All info
  std::vector<G4int>& GetActivePMT() { return ActivePMT; }
  std::vector<G4int>& GetNumberPMTHits() { return NumberPMTHits; }
  std::vector<G4double>& GetPmtMedianHitTimes() { return PmtMedianHitTimes; }
  std::vector<G4double>& GetPmtAllHitTimes() { return PmtAllHitTimes; }
  std::vector<G4double>& GetTotalEnergyDeposit() { return TotalEnergyDeposit; }
  std::vector<G4double>& GetTimeOfFlight() { return TimeOfFlight; }

  // Output MilliQCAEN info
  std::vector<G4float>& GetTriggerCount() { return TriggerCount;}
  std::vector<G4float>& GetTimeCount() { return TimeCount;}
  std::vector<G4int>& GetEventID() {return EventID;}
  std::vector<G4int>& GetTDC() { return TDC;}
  std::vector<G4float>& GetWaveMax() {return WaveMax;}
  std::vector<G4float>& GetWaveMin() { return WaveMin;}
  std::vector<G4float>& GetChargeIntegral() { return ChargeIntegral;}


 private:

  MilliQRecorderBase*   fRecorder;
  MilliQEventMessenger* fEventMessenger;

  G4int fSaveThreshold;

  G4int fPMTCollID;
  G4int fPMTAllCollID;
  G4int fScintCollID;

  G4int fVerbose;

  G4int fPMTThreshold;

  G4bool fForcedrawphotons;
  G4bool fForcenophotons;

  G4int NBlocks;
  G4int NStacks;

  // Output all info
  std::vector<G4int> ActivePMT;
  std::vector<G4int> NumberPMTHits;
  std::vector<G4double> PmtAllHitTimes;
  std::vector<G4double> PmtMedianHitTimes;
  std::vector<G4double> TotalEnergyDeposit;
  std::vector<G4double> TimeOfFlight;

  // Output MilliQCAEN info
  std::vector<G4float> TriggerCount;
  std::vector<G4float> TimeCount;
  std::vector<G4int> EventID;
  std::vector<G4int> TDC;
  std::vector<G4float> WaveMax;
  std::vector<G4float> WaveMin;
  std::vector<G4float> ChargeIntegral;


  const boost::property_tree::ptree fPTree;
  boost::property_tree::ptree fGeometryPTree;


};

#endif
