
#include "G4VPVParameterisation.hh"
#include "G4Box.hh"
#include "G4VPhysicalVolume.hh"

// For G4ThreeVector storage of dimensional data:
//      x dimension is treated as length
//      y dimension is treated as height
//      z dimension is treated as width

class MilliQDetectorStackParameterisation : public G4VPVParameterisation
{
public:
    
    MilliQDetectorStackParameterisation(G4int  pN,
                                        G4ThreeVector  pStackDimensions,
                                        G4ThreeVector  pDistributionUnitVector,
										G4ThreeVector  pOffset,
                                        G4double       pStartDepth,
                                        G4double       pEndDepth);
    
    virtual ~MilliQDetectorStackParameterisation();
    
    
    virtual void ComputeTransformation(const G4int pId,
                                       G4VPhysicalVolume* pPV) const;

    
    G4int GetNumberOfBlocks(){return fN;}
    G4ThreeVector GetStackDimensions();
    
private:
    
    G4double    fStartCenterDepth;
    G4double    fDepthGap;
    
    G4int fN;
    G4ThreeVector fStackDimensions;
    G4ThreeVector fDistributionUnitVector;
    G4ThreeVector fOffset;
};
