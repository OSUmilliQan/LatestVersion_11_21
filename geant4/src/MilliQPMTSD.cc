
#include "MilliQPMTSD.hh"
#include "MilliQPMTHit.hh"
#include "MilliQDetectorConstruction.hh"
#include "MilliQUserTrackInformation.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4ParticleTypes.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"

MilliQPMTSD::MilliQPMTSD(G4String name, const G4int numBlocks, const G4int numStacks)
  : G4VSensitiveDetector(name),fPMTHitCollection(0), fPMTAllHitCollection(0), NBlocks(numBlocks), NStacks(numStacks)
{
  collectionName.insert("pmtHitCollection");
  collectionName.insert("pmtAllHitCollection");

}

MilliQPMTSD::~MilliQPMTSD() {}

void MilliQPMTSD::Initialize(G4HCofThisEvent* hitsCE) {

  fPMTHitCollection = new MilliQPMTHitsCollection(SensitiveDetectorName,
						  collectionName[0]);
  fPMTAllHitCollection = new MilliQPMTHitsCollection(SensitiveDetectorName,
						     collectionName[1]);

  //Store collection with event and keep ID
  static G4int hitCID = -1;
  if(hitCID<0){
    // hitCID = GetCollectionID(0);
    hitCID = G4SDManager::GetSDMpointer()->GetCollectionID(fPMTHitCollection);
  }
  hitsCE->AddHitsCollection( hitCID, fPMTHitCollection );

  static G4int AllhitCID = -1;
  if(AllhitCID<0){
    // hitCID = GetCollectionID(0);
    AllhitCID = G4SDManager::GetSDMpointer()->GetCollectionID(fPMTAllHitCollection);
  }
  hitsCE->AddHitsCollection( AllhitCID, fPMTAllHitCollection );

  // fill calorimeter hits with zero energy deposition
  for(G4int iStack = 0; iStack < NStacks; iStack++) {
    for(G4int iBlock = 0; iBlock < NBlocks; iBlock++) { // Total blocks in 1 stack
      MilliQPMTHit* hit = new MilliQPMTHit();
      fPMTHitCollection->insert(hit);
    }
  }
}

//Calling this does the regular Sensitive detector business
//It calculates the energy deposed in the PMT by particles
G4bool MilliQPMTSD::ProcessHits(G4Step* step,G4TouchableHistory* ){

  G4double edep = step->GetTotalEnergyDeposit();
  if (edep==0.) return true;

  G4TouchableHistory* touchable
    = (G4TouchableHistory*)(step->GetPreStepPoint()->GetTouchable());

  G4VPhysicalVolume* cellPhysical = touchable->GetVolume(3); // Block Parameterization
  G4int blockNo = cellPhysical->GetCopyNo();

  G4VPhysicalVolume* stackPhysical = touchable->GetVolume(4); //Stack Parameterization
  G4int stackNo = stackPhysical->GetCopyNo();

  /*
    G4cout<<"Volume() "<<touchable->GetVolume()->GetName() <<G4endl;
    G4cout<<"Volume 0 "<<touchable->GetVolume(0)->GetName() <<G4endl;
    G4cout<<"Volume 1 "<<touchable->GetVolume(1)->GetName() <<G4endl;
    G4cout<<"Volume 2 "<<touchable->GetVolume(2)->GetName() <<G4endl;
    G4cout<<"Volume 3 "<<touchable->GetVolume(3)->GetName() <<G4endl;
    G4cout<<"Volume 4 "<<touchable->GetVolume(4)->GetName() <<G4endl;
    G4cout<<"Volume 5 "<<touchable->GetVolume(5)->GetName() <<G4endl;
  */

  G4int hitID = stackNo*NBlocks+blockNo;
  MilliQPMTHit* hit = (*fPMTHitCollection)[hitID];

  //	    G4cout<<"Block Number "<<blockNo<<" Stack Number "<<stackNo<<" COPY NUMBER "<<hitID<<G4endl;

  // check if it is first touch
  if (hit->GetStackID()<0)
    {
      hit->SetStackID(stackNo);
      hit->SetBlockID(blockNo);

      G4int depth = touchable->GetHistory()->GetDepth();
      G4AffineTransform transform
	= touchable->GetHistory()->GetTransform(depth-2);
      transform.Invert();
      hit->SetRot(transform.NetRotation());
      hit->SetPos(transform.NetTranslation());
    }

  // add energy deposition
  hit->AddEdep(edep);

  return true;
}


//Generates a hit and uses the postStepPoint's mother volume replica number
//PostStepPoint because the hit is generated manually when the photon is
//absorbed by the photocathode (uses skin settings defined in the detector
G4bool MilliQPMTSD::ProcessHits_constStep(const G4Step* aStep,
                                          G4TouchableHistory*)
{

  //need to know if this is an optical photon
  if(aStep->GetTrack()->GetDefinition() != G4OpticalPhoton::OpticalPhotonDefinition()) {
    return false;
  }

  G4int pmtNumberBlock = aStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(3);//GetReplicaNumber(1);
  G4int pmtNumberStack = aStep->GetPostStepPoint()->GetTouchable()->GetReplicaNumber(4);//GetReplicaNumber(1);
  G4int pmtNumber = pmtNumberStack*NBlocks + pmtNumberBlock;

  G4double hitTime = aStep->GetPostStepPoint()->GetGlobalTime();

  G4VPhysicalVolume* physVol = aStep->GetPostStepPoint()->GetTouchable()->GetVolume(3);

  //    G4cout<<"pmtNumberBlock "<<pmtNumberBlock<<" pmtNumberStack "<<pmtNumberStack<<" pmtNumber "<<pmtNumber<<G4endl;

  //Find the correct hit collection
  G4int n = fPMTHitCollection->entries();
  MilliQPMTHit* hit = NULL;
  for(G4int i = 0; i < n; i++) {
    if((*fPMTHitCollection)[i]->GetPMTNumber() == pmtNumber) {
      hit = (*fPMTHitCollection)[i];
      break;
    }
  }

  if(hit == NULL) { // this pmt wasnt previously hit in this event
    hit = new MilliQPMTHit(); //so create new hit
    hit->SetPMTNumber(pmtNumber);
    hit->SetPMTPhysVol(physVol);
    fPMTHitCollection->insert(hit);
  }

  hit->IncPhotonCount(); //increment hit for the selected pmt

  hit->SetDrawit(true);

  MilliQPMTHit* Allhit = new MilliQPMTHit();
  Allhit->SetPMTPhysVol(physVol);
  Allhit->SetPMTNumber(pmtNumber);
  Allhit->SetTime(hitTime);

  fPMTAllHitCollection->insert(Allhit);

  return true;
}


void MilliQPMTSD::EndOfEvent(G4HCofThisEvent* ) {}


void MilliQPMTSD::clear() {}


void MilliQPMTSD::DrawAll() {}


void MilliQPMTSD::PrintAll() {}
