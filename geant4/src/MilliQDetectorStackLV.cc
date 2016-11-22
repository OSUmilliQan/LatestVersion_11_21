
#include "globals.hh"

#include "MilliQDetectorStackLV.hh"

#include "G4LogicalSkinSurface.hh"
#include "G4LogicalBorderSurface.hh"

#include "G4SystemOfUnits.hh"


MilliQDetectorStackLV::MilliQDetectorStackLV(G4VSolid*               pSolid,
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
										                         G4double				         pLightGuideLength,
                                             G4MaterialPropertiesTable * pScintillatorHousingPT,

                                             G4double                pPmtPhotocathodeRadius,
                                             G4double                pPmtPhotocathodeHeight,
                                             G4MaterialPropertiesTable * pLGHousingPT,
                                             G4MaterialPropertiesTable * pPmtPT,
                                             G4VSensitiveDetector*   pPmtSD,
										                         G4VSensitiveDetector*   pScintSD)
                    : G4LogicalVolume(pSolid,
                                      pMaterial,
                                      pName,
                                      pFieldManager,
                                      pSD,
                                      pUserLimits,
                                      pOptimise)
{
    //
    // Detector Blocks
    //

    // Detector Blocks - Volume
    G4Box* houseingV = new G4Box("Detector Block Housing Volume" ,
                                 1*m , 1*m , 1*m );//temp dimension

    // Detector Blocks - Logical Volume
    MilliQDetectorBlockLV* aDetectorBlock
    = new MilliQDetectorBlockLV(houseingV, //pSolid
                                G4Material::GetMaterial("Air"), //pMaterial
                                "Detector Block Housing Logical Volume", //pName
                                0, //pFieldManager
                                0, //pSD
                                0, //pUserLimits
                                true, //pOptimise

                                pScintillatorDimensions, //pScintillatorDimensions
                                pScintillatorHousingThickness, //pScintillatorHousingThickness
								                pLightGuideLength,
                                pScintillatorHousingPT,

								                pPmtPhotocathodeRadius, //pPmtHeight
								                pPmtPhotocathodeHeight, //pPmtPhotocathodeDepth
                                pLGHousingPT,
                                pPmtPT,
                                pPmtSD,
								                pScintSD);

    // Detector Blocks - Parameterisation
    fDetectorBlockParameterisation
    = new MilliQDetectorBlockParameterisation(pNumberOfBlocks,                  //n
                                              aDetectorBlock->GetDimensions(),  //block dimensions
                                              pBetweenBlockSpacing);            //gap dimensions

    // Detector Blocks - Physical Volume
    new G4PVParameterised( "Detector Block Housing Physical Volume",
                            aDetectorBlock,
                            this,
                            kZAxis,
                            fDetectorBlockParameterisation->GetNumberOfBlocks(),
                            fDetectorBlockParameterisation);


//    G4cout<<"The number of block parametrizations "<<fDetectorBlockParameterisation->GetNumberOfBlocks()<<G4endl;


    AdjustVolume();
    VisAttributes();
}


void MilliQDetectorStackLV::VisAttributes()
{
    this->SetVisAttributes(G4VisAttributes::Invisible);
}


void MilliQDetectorStackLV::AdjustVolume()
{
    //set final dimensions of stack housing and set this logical volume
    //
    fStackDimensions = fDetectorBlockParameterisation->GetStackDimensions();

    G4Box* stackHouseingV = (G4Box*)this->GetSolid();
    stackHouseingV->SetXHalfLength(fStackDimensions.x()/2.);
    stackHouseingV->SetYHalfLength(fStackDimensions.y()/2.);
    stackHouseingV->SetZHalfLength(fStackDimensions.z()/2.);
    this->SetSolid(stackHouseingV);
}
