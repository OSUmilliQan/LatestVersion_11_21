
#ifndef MilliQDetectorStack_H
#define MilliQDetectorStack_H 1

#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"
#include "G4PVParameterised.hh"

#include "MilliQDetectorBlockParameterisation.hh"
#include "MilliQDetectorBlockLV.hh"
#include "MilliQPMTSD.hh"

class MilliQDetectorStackLV : public G4LogicalVolume
{
public:

    MilliQDetectorStackLV(G4VSolid*               pSolid,
                          G4Material*             pMaterial,
                          const G4String&         pName,
                          G4FieldManager*         pFieldManager,
                          G4VSensitiveDetector*   pSD,
                          G4UserLimits*           pUserLimits,
                          G4bool                  pOptimise,

                          G4ThreeVector           pNumberOfBlocks,
                          G4ThreeVector           pBetweenBlockSpacing,

                          G4ThreeVector           pScintillatorDimensions,
                          G4double                pScintillatorHousingThickness,
						              G4double				        pLightGuideLength,
                          G4MaterialPropertiesTable * pScintillatorHousingPT,

                          G4double                pPmtPhotocathodeRadius,
                          G4double                pPmtPhotocathodeHeight,
                          G4MaterialPropertiesTable * pLGHousingPT,
                          G4MaterialPropertiesTable * pPmtPT,
                          G4VSensitiveDetector*   pPmtSD,
						              G4VSensitiveDetector*   pScintSD);

    G4ThreeVector GetDimensions() { return fStackDimensions; }

private:
    MilliQDetectorBlockParameterisation* fDetectorBlockParameterisation;
    G4ThreeVector fStackDimensions;
    void VisAttributes();
    void AdjustVolume();
};

#endif
