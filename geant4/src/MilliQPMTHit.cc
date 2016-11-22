
#include "MilliQPMTHit.hh"
#include "MilliQDetectorConstruction.hh"

#include "G4ios.hh"
#include "G4VVisManager.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Box.hh"
#include "G4AttDefStore.hh"
#include "G4AttDef.hh"
#include "G4AttValue.hh"
#include "G4UIcommand.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


G4ThreadLocal G4Allocator<MilliQPMTHit>* MilliQPMTHitAllocator=0;


MilliQPMTHit::MilliQPMTHit()
  : G4VHit(), fPmtNumber(-1), fPhotons(0), fPos(0), fPhysVol(0), fDrawit(false), fStackID(-1), fBlockID(-1), fEdep(0.), fTime(0.) {}

//MilliQPMTHit::MilliQPMTHit(G4int iCol, G4int iBlock)
//  : G4VHit(), fPmtNumber(-1), fPhotons(0), fPos(0), fPhysVol(0), fDrawit(false), fStackID(-1), fBlockID(-1), fEdep(0.), fTime(0.) {}

MilliQPMTHit::~MilliQPMTHit() {}

MilliQPMTHit::MilliQPMTHit(const MilliQPMTHit &right) : G4VHit()
{
  fPmtNumber = right.fPmtNumber;
  fPhotons = right.fPhotons;
  fPhysVol = right.fPhysVol;
  fDrawit = right.fDrawit;
  fStackID = right.fStackID;
  fBlockID = right.fBlockID;
  fEdep = right.fEdep;
  fTime = right.fTime;
  fPos = right.fPos;
  fRot = right.fRot;
}

const MilliQPMTHit& MilliQPMTHit::operator=(const MilliQPMTHit &right){
  fPmtNumber = right.fPmtNumber;
  fPhotons = right.fPhotons;
  fPhysVol = right.fPhysVol;
  fDrawit = right.fDrawit;
  fStackID = right.fStackID;
  fBlockID = right.fBlockID;
  fEdep = right.fEdep;
  fTime = right.fTime;
  fPos = right.fPos;
  fRot = right.fRot;
  return *this;
}


G4int MilliQPMTHit::operator==(const MilliQPMTHit &right) const {
  return (fPmtNumber == right.fPmtNumber &&
          fStackID == right.fStackID &&
          fBlockID == right.fBlockID);
}


void MilliQPMTHit::Draw(){
  if(fDrawit&&fPhysVol){ //ReDraw only the PMTs that have hit counts > 0
    //Also need a physical volume to be able to draw anything
    G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
    if(pVVisManager){//Make sure that the VisManager exists
      G4VisAttributes attribs(G4Colour(1.,0.,0.));
      attribs.SetForceSolid(true);
      G4RotationMatrix rot;
      if(fPhysVol->GetRotation())//If a rotation is defined use it
        rot=*(fPhysVol->GetRotation());
      G4Transform3D trans(rot,fPhysVol->GetTranslation());//Create transform
      pVVisManager->Draw(*fPhysVol,attribs,trans);//Draw it
    }
  }
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

const std::map<G4String,G4AttDef>* MilliQPMTHit::GetAttDefs() const
{
    G4bool isNew;
    std::map<G4String,G4AttDef>* store
    = G4AttDefStore::GetInstance("MilliQPMTHit",isNew);

    if (isNew) {
        (*store)["HitType"]
          = G4AttDef("HitType","Hit Type","Physics","","G4String");

        (*store)["Stack"]
          = G4AttDef("Stack","Stack ID","Physics","","G4int");

        (*store)["Block"]
          = G4AttDef("Block","Block ID","Physics","","G4int");

        (*store)["Energy"]
          = G4AttDef("Energy","Energy Deposited","Physics","G4BestUnit",
                     "G4double");

        (*store)["Time"]
          = G4AttDef("Time","Time of Hit","Physics","G4BestUnit",
                     "G4double");

        (*store)["Pos"]
          = G4AttDef("Pos", "Position", "Physics","G4BestUnit",
                     "G4ThreeVector");
    }
    return store;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

std::vector<G4AttValue>* MilliQPMTHit::CreateAttValues() const
{
    std::vector<G4AttValue>* values = new std::vector<G4AttValue>;

    values
      ->push_back(G4AttValue("HitType","MilliQPMTHit",""));
    values
      ->push_back(G4AttValue("Stack",G4UIcommand::ConvertToString(fStackID),
                             ""));
    values
      ->push_back(G4AttValue("Block",G4UIcommand::ConvertToString(fBlockID),""));
    values
      ->push_back(G4AttValue("Energy",G4BestUnit(fEdep,"Energy"),""));
    values
      ->push_back(G4AttValue("Time",G4BestUnit(fTime,"Time"),""));
    values
      ->push_back(G4AttValue("Pos",G4BestUnit(fPos,"Length"),""));

    return values;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MilliQPMTHit::Print()
{
    G4cout << "  Cell[" << fBlockID << ", " << fStackID << "] "
    << fEdep/MeV << " (MeV) " << fPos << G4endl;
    G4cout << "*************HIT****************"<<G4endl;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
