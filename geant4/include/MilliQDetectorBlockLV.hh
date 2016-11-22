
#ifndef MilliQDetectorBlock_H
#define MilliQDetectorBlock_H 1

#include "G4PVPlacement.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Trd.hh"
#include "G4Material.hh"
#include "G4LogicalVolume.hh"
#include "G4OpticalSurface.hh"

#include "MilliQDetectorConstruction.hh"
#include "MilliQPMTSD.hh"

class MilliQDetectorBlockLV : public G4LogicalVolume
{
  public:

    MilliQDetectorBlockLV(G4VSolid*              pSolid,
                          G4Material*            pMaterial,
                          const G4String&        pName,
                          G4FieldManager*        pFieldManager,
                          G4VSensitiveDetector*  pSD,
                          G4UserLimits*          pUserLimits,
                          G4bool                 pOptimise,

                          G4ThreeVector  pScintillatorDimensions,
                          G4double       pScintillatorHousingThickness,
						              G4double		   pLightGuideLength,
                          G4MaterialPropertiesTable * pScintillatorHousingPT,

                          G4double  pPmtPhotocathodeRadius,		// PMT radius
                          G4double  pPmtPhotocathodeHeight,     // Height of PMT
                          G4MaterialPropertiesTable * pLGHousingPT,
                          G4MaterialPropertiesTable * pPmtPT,
                          G4VSensitiveDetector*   pPmtSD,
						              G4VSensitiveDetector*   pScintSD);

    G4ThreeVector GetDimensions() { return fDimensions; };

  private:

    void VisAttributes();
    void SurfaceProperties(G4MaterialPropertiesTable * pScintillatorHousingPT,
                           G4MaterialPropertiesTable * pLGHousingPT,
                           G4MaterialPropertiesTable * pPmtPT);

    G4ThreeVector    fDimensions;

    G4LogicalVolume* fScintillatorHousingLV;
    G4LogicalVolume* fScintillatorLGHousingLV;
    G4LogicalVolume* fScintillatorLV;
    G4LogicalVolume* fLightGuideLV;

    G4LogicalVolume* fPmtGlassSectionLV;
    G4LogicalVolume* fPmtPhotocathodeSectionLV;
};

#endif
