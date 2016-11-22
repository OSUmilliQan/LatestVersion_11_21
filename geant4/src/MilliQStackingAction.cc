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
// $Id: MilliQStackingAction.cc 68752 2013-04-05 10:23:47Z gcosmo $
//
/// \file optical/MilliQ/src/MilliQStackingAction.cc
/// \brief Implementation of the MilliQStackingAction class
//
//
#include "MilliQStackingAction.hh"
#include "MilliQUserEventInformation.hh"
#include "MilliQSteppingAction.hh"

#include "G4ios.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MilliQStackingAction::MilliQStackingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MilliQStackingAction::~MilliQStackingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack MilliQStackingAction::ClassifyNewTrack(const G4Track * aTrack) {
 
  MilliQUserEventInformation* eventInformation = (MilliQUserEventInformation*)G4EventManager::GetEventManager()->GetConstCurrentEvent()->GetUserInformation();
 
  //Count what process generated the optical photons
  if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition()) {
    // particle is optical photon
    if(aTrack->GetParentID() > 0) {
      // particle is secondary
      if(aTrack->GetCreatorProcess()->GetProcessName() == "Scintillation") eventInformation->IncPhotonCount_Scint();
      else if(aTrack->GetCreatorProcess()->GetProcessName() == "Cerenkov") eventInformation->IncPhotonCount_Ceren();
    }
  }
  else {
    //    
  }
  return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MilliQStackingAction::NewStage() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MilliQStackingAction::PrepareNewEvent() {}
