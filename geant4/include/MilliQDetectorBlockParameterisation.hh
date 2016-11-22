
#include "G4VPVParameterisation.hh"
#include "G4Box.hh"
#include "G4VPhysicalVolume.hh"
#include "MilliQDetectorBlockLV.hh"

// For G4ThreeVector storage of dimensional data:
//      x dimension is treated as width
//      y dimension is treated as height
//      z dimension is treated as length

class MilliQDetectorBlockParameterisation : public G4VPVParameterisation
{
public:
    
    MilliQDetectorBlockParameterisation(G4ThreeVector  pN,
                                        G4ThreeVector  pBlockDimensions,
                                        G4ThreeVector  pGapDimensions);
    
    virtual ~MilliQDetectorBlockParameterisation();
    
    
    virtual void ComputeTransformation(const G4int pId,
                                       G4VPhysicalVolume* pPV) const;
    
    void GetIdDimensionsFromId(G4int           pId,
                               G4ThreeVector&  pIdDimensions) const;
    
    G4int GetNumberOfBlocks(){return (G4int)(fN.x()*fN.y()*fN.z());}
    G4ThreeVector GetStackDimensions();
    
private:
    
    G4ThreeVector fN;
    G4ThreeVector fBlockDimensions, fGapDimensions;
};
