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
// $Id: MilliQScintSD.cc 68752 2013-04-05 10:23:47Z gcosmo $
//
//
//
#include "../include/MilliQScintHit.hh"
#include "../include/MilliQScintSD.hh"
#include "MilliQDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MilliQScintSD::MilliQScintSD(G4String name, const G4int numBlocks)
  : G4VSensitiveDetector(name), NBlocks(numBlocks)
{
  fScintCollection = NULL;
  collectionName.insert("scintCollection");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MilliQScintSD::~MilliQScintSD() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MilliQScintSD::Initialize(G4HCofThisEvent* hitsCE){
  fScintCollection = new MilliQScintHitsCollection
                      (SensitiveDetectorName,collectionName[0]);
  //A way to keep all the hits of this event in one place if needed
  static G4int hitsCID = -1;
  if(hitsCID<0){
    hitsCID = GetCollectionID(0);
  }
  hitsCE->AddHitsCollection( hitsCID, fScintCollection );

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool MilliQScintSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ){
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep==0.) return false; //No edep so dont count as hit

  G4StepPoint* thePrePoint = aStep->GetPreStepPoint();
  G4TouchableHistory* theTouchable =
    (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  G4VPhysicalVolume* thePrePV = theTouchable->GetVolume();



  G4int scintNumberBlock = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(2);//GetReplicaNumber(1);
  G4int scintNumberStack = aStep->GetPreStepPoint()->GetTouchable()->GetReplicaNumber(3);//GetReplicaNumber(1);
  G4int scintNumber = scintNumberStack*NBlocks+scintNumberBlock;
 // G4cout<<"scintNumberBlock "<<scintNumberBlock<<" scintNumberStack "<<scintNumberStack<<" scintNumber "<<scintNumber<<G4endl;

  G4double hitTime = thePrePoint->GetGlobalTime();

  G4StepPoint* thePostPoint = aStep->GetPostStepPoint();

  //Get the average position of the hit
  G4ThreeVector pos = thePrePoint->GetPosition() + thePostPoint->GetPosition();
  pos/=2.;

  MilliQScintHit* scintHit = new MilliQScintHit(thePrePV);
  scintHit->SetEdep(edep);
  scintHit->SetPos(pos);
  scintHit->SetCpNum(scintNumber);
  scintHit->SetTime(hitTime);


  fScintCollection->insert(scintHit);


  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MilliQScintSD::EndOfEvent(G4HCofThisEvent* ) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MilliQScintSD::clear() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MilliQScintSD::DrawAll() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MilliQScintSD::PrintAll() {}
