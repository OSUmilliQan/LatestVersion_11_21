
#ifndef MilliQPMTHit_h
#define MilliQPMTHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"


#include "tls.hh"

class G4VTouchable;
class G4AttDef;
class G4AttValue;

class MilliQPMTHit : public G4VHit
{
  public:

    MilliQPMTHit();
//    MilliQPMTHit(G4int iCol,G4int iBlock);
    virtual ~MilliQPMTHit();
    MilliQPMTHit(const MilliQPMTHit &right);

    const MilliQPMTHit& operator=(const MilliQPMTHit &right);
    G4int operator==(const MilliQPMTHit &right) const;

    inline void *operator new(size_t);
    inline void operator delete(void *aHit);

    virtual void Draw();
    virtual const std::map<G4String,G4AttDef>* GetAttDefs() const;
    virtual std::vector<G4AttValue>* CreateAttValues() const;
    virtual void Print();

    inline void SetDrawit(G4bool b) { fDrawit = b; }
    inline G4bool GetDrawit() { return fDrawit; }

    inline void IncPhotonCount() { fPhotons++; }
    inline G4int GetPhotonCount() { return fPhotons; }

    inline void SetPMTNumber(G4int n) { fPmtNumber = n; }
    inline G4int GetPMTNumber() { return fPmtNumber; }

    inline void SetPMTPhysVol(G4VPhysicalVolume* physVol) { this->fPhysVol = physVol; }
    inline G4VPhysicalVolume* GetPMTPhysVol() { return fPhysVol; }

//    inline void SetPMTPos(G4double x,G4double y,G4double z){
 //     fPos=G4ThreeVector(x,y,z);

//    }
     inline G4ThreeVector GetPMTPos() { return fPos; }


     //I added these:

      void SetStackID(G4int z) { fStackID = z; }
      G4int GetStackID() const { return fStackID; }

      void SetBlockID(G4int z) { fBlockID = z; }
      G4int GetBlockID() const { return fBlockID; }

      void SetEdep(G4double de) { fEdep = de; }
      void AddEdep(G4double de) { fEdep += de; }
      G4double GetEdep() const { return fEdep; }

      void SetTime(G4double t) { fTime = t; }
      G4double GetTime() const { return fTime; }

      void SetPos(G4ThreeVector xyz) { fPos = xyz; }
      G4ThreeVector GetPos() const { return fPos; }

      void SetRot(G4RotationMatrix rmat) { fRot = rmat; }
      G4RotationMatrix GetRot() const { return fRot; }


  private:

    G4int fPmtNumber;
    G4int fPhotons;
    G4ThreeVector fPos;
    G4VPhysicalVolume* fPhysVol;
    G4bool fDrawit;

    G4int fStackID;
    G4int fBlockID;
    G4double fEdep;
    G4double fTime;

    G4RotationMatrix fRot;

};

typedef G4THitsCollection<MilliQPMTHit> MilliQPMTHitsCollection;

extern G4ThreadLocal G4Allocator<MilliQPMTHit>* MilliQPMTHitAllocator;

inline void* MilliQPMTHit::operator new(size_t){
  if(!MilliQPMTHitAllocator)
      MilliQPMTHitAllocator = new G4Allocator<MilliQPMTHit>;
  return (void *) MilliQPMTHitAllocator->MallocSingle();
}

inline void MilliQPMTHit::operator delete(void *aHit){
  MilliQPMTHitAllocator->FreeSingle((MilliQPMTHit*) aHit);
}

#endif
