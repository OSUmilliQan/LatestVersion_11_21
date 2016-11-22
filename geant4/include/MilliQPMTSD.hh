
#ifndef MilliQPMTSD_h
#define MilliQPMTSD_h 1

#include "G4DataVector.hh"
#include "G4VSensitiveDetector.hh"
#include "MilliQPMTHit.hh"

#include <vector>

class G4Step;
class G4HCofThisEvent;

class MilliQPMTSD : public G4VSensitiveDetector {

public:

  MilliQPMTSD(G4String name, const G4int numBlocks, const G4int numStacks);
  virtual ~MilliQPMTSD();

  virtual void Initialize(G4HCofThisEvent* );
  virtual G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* );

  //A version of processHits that keeps aStep constant
  G4bool ProcessHits_constStep(const G4Step* ,
			       G4TouchableHistory* );
  virtual void EndOfEvent(G4HCofThisEvent* );
  virtual void clear();
  void DrawAll();
  void PrintAll();

private:

  MilliQPMTHitsCollection* fPMTHitCollection;
  MilliQPMTHitsCollection* fPMTAllHitCollection; //Record details of all PMT hits

  const G4int NBlocks;
  const G4int NStacks;

};

#endif
