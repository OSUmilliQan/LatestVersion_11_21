#include "globals.hh"

#include "MilliQDetectorBlockLV.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"

#include "G4LogicalSkinSurface.hh"
#include "G4RotationMatrix.hh"
#include "G4Transform3D.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4PhysicalConstants.hh"

#include "G4SystemOfUnits.hh"

MilliQDetectorBlockLV::MilliQDetectorBlockLV(G4VSolid* pSolid,
					     															 G4Material* pMaterial,
																						 const G4String& pName,
					     												 			 G4FieldManager* pFieldManager,
																						 G4VSensitiveDetector* pSD,
					     												 			 G4UserLimits* pUserLimits,
																						 G4bool pOptimise,

					     												 			 G4ThreeVector pScintillatorDimensions,
					     												 			 G4double pScintillatorHousingThickness,
																						 G4double pLightGuideLength,
																						 G4MaterialPropertiesTable * pScintillatorHousingPT,

					     												 			 G4double pPmtPhotocathodeRadius,
					     												 			 G4double pPmtPhotocathodeHeight,
					     												 			 G4MaterialPropertiesTable * pLGHousingPT,
																						 G4MaterialPropertiesTable * pPmtPT,
																						 G4VSensitiveDetector* pPmtSD,
					     												 			 G4VSensitiveDetector* pScintSD) :
  G4LogicalVolume(pSolid, pMaterial, pName, pFieldManager, pSD, pUserLimits, pOptimise) {

  if (pLightGuideLength > 0.9 * pScintillatorDimensions.x())
    pLightGuideLength = 0.9 * pScintillatorDimensions.x();

  if (pPmtPhotocathodeHeight > pScintillatorHousingThickness)
    pPmtPhotocathodeHeight = pScintillatorHousingThickness;

  if (pPmtPhotocathodeRadius
      > std::min(pScintillatorDimensions.y() / 2.,
		 pScintillatorDimensions.z() / 2.)) {
    pPmtPhotocathodeRadius = std::min(pScintillatorDimensions.y() / 2.,
				      pScintillatorDimensions.z() / 2.);
  }

  //
  // Detector Block(this)
  //

  fDimensions = G4ThreeVector(
			      pScintillatorDimensions.x() + 2*pScintillatorHousingThickness,
			      pScintillatorDimensions.y() + 2*pScintillatorHousingThickness,
			      pScintillatorDimensions.z() + 2*pScintillatorHousingThickness);

  G4Box* detectorBlockV = new G4Box("Detector Block Volume",
				    												fDimensions.x() / 2., fDimensions.y() / 2., fDimensions.z() / 2.);

  this->SetSolid(detectorBlockV);

  //****************************************************************************************************************************

  // Scintillator Housing - Volume
  G4Box* scintillatorHousingV = new G4Box("Scintillator Housing Volume", //name
					  															pScintillatorDimensions.x() / 2. - pLightGuideLength/2+pScintillatorHousingThickness/2, //half x dimension
					  															pScintillatorDimensions.y() / 2. + pScintillatorHousingThickness, //half y dimension
					  															pScintillatorDimensions.z() / 2. + pScintillatorHousingThickness); //half z dimension

  // Scintillator Housing - Logical Volume
  fScintillatorHousingLV = new G4LogicalVolume(scintillatorHousingV, //volume
					       															 G4Material::GetMaterial("Aluminium"), //material
					       											 				 "Scintillator Housing Logical Volume"); //name
  // Scintillator Housing - Physical Volume
  new G4PVPlacement(0,                                              //rotation
		    						G4ThreeVector(-pLightGuideLength/2 - pScintillatorHousingThickness/2, 0, 0),     //translation
		    						fScintillatorHousingLV,                             //logical volume
		    						"Scintillator Housing Physical Volume",                       //name
		    						this,                                        //mother logical volume
		    						false,                                                        //many
		    						0);                                                         //copy n

  // Scintillator - Volume
  G4Box* scintillatorV = new G4Box("Scintillator Volume",             //name
				   												 pScintillatorDimensions.x() / 2. - pLightGuideLength/2, //half x dimension
				   										 		 pScintillatorDimensions.y() / 2.,    //half y dimension
				   										 		 pScintillatorDimensions.z() / 2.);   //half z dimension

  // Scintillator - Logical Volume
  fScintillatorLV = new G4LogicalVolume(scintillatorV,                //volume
																				G4Material::GetMaterial("Scintillator"),  //material
																				"Scintillator Logical Volume",			//name
																				0, 										//field manager
																				pScintSD, 								//sensitive detector
																				0, 										//userlimits
																				true); 									//optimise
  // Scintillator - Physical Volume
  new G4PVPlacement(0,                                    //rotation
		    						G4ThreeVector(pScintillatorHousingThickness/2, 0, 0),               //translation
		    						fScintillatorLV,                      //logical volume
		    						"Scintillator Physical Volume",       //name
		    						fScintillatorHousingLV,               //mother logical volume
		    						false,                                //many
		    						0);                                   //copy n

  //////

  // Scintillator LG Housing - Volume
  G4Box* scintillatorLGHousingV = new G4Box("Scintillator LG Housing Volume", //name
					    															pLightGuideLength/2 + pScintillatorHousingThickness/2, //half x dimension
					    															pScintillatorDimensions.y() / 2. + pScintillatorHousingThickness, //half y dimension
					    															pScintillatorDimensions.z() / 2. + pScintillatorHousingThickness); //half z dimension

  // Scintillator LG Housing - Logical Volume
  fScintillatorLGHousingLV = new G4LogicalVolume(scintillatorLGHousingV, //volume
						 																		 G4Material::GetMaterial("Aluminium"), //material
						 																	 	 "Scintillator LG Housing Logical Volume"); //name
  // Scintillator LG Housing - Physical Volume
  new G4PVPlacement(0,                                              //rotation
		    						G4ThreeVector(pScintillatorDimensions.x() / 2. - pLightGuideLength/2 + pScintillatorHousingThickness/2, 0, 0),     //translation
		    						fScintillatorLGHousingLV,                             //logical volume
		    						"Scintillator LG Housing Physical Volume",                       //name
		    						this,                                        //mother logical volume
		    						false,                                                        //many
		    						0);                                                         //copy n


  G4RotationMatrix* rotm = new G4RotationMatrix(); // This orientation makes the sensitive part farthest away from beam
  rotm->rotateX(90 * deg);
  rotm->rotateY(270 * deg);
  rotm->rotateZ(0 * deg);

  //Light Guide
  // Light Guide - Volume

  G4Trd* lightGuideV = new G4Trd("Light Guide Volume",
				 												 pScintillatorDimensions.y() / 2, pScintillatorDimensions.y() / 2,
				 											 	 pScintillatorDimensions.z() / 2, 3*pScintillatorDimensions.z()/6,
				 											 	 pLightGuideLength/2);

  // Light Guide - Logical Volume
  fLightGuideLV = new G4LogicalVolume(lightGuideV,                    //volume
				      												G4Material::GetMaterial("Scintillator"),  //material
				      												"Light Guide Logical Volume",			//name
				      												0, 										//field manager
				      												pScintSD, 								//sensitive detector
				      												0, 										//userlimits
				      												true); 									//optimise
  // Light Guide - Physical Volume
  new G4PVPlacement(rotm,                                    //rotation
		    						G4ThreeVector(-pScintillatorHousingThickness/2, 0, 0),                 //translation
		    						fLightGuideLV,                       //logical volume
		    						"Light Guide Physical Volume",       //name
		    						fScintillatorLGHousingLV,               //mother logical volume
		    						false,                                //many
		    						0);                                   //copy n
  //////////////////////

  //
  // PMT
  // Composed of "PMT Glass", which contains (on the farthest side from the scintillator)
  // the "PMT Photocathode Section".


  // PMT Glass Section - Volume
  G4Tubs* pmtGlassSectionV = new G4Tubs("PMT Glass Section Volume", //name
																				0. * cm,                                   //inner radius
																				pPmtPhotocathodeRadius,           //outer radius
																				pPmtPhotocathodeHeight / 2.,          //half height
																				0. * deg,                                  //start angle
																				360. * deg);                               //end angle

  // PMT Glass Section - Logical Volume
  fPmtGlassSectionLV = new G4LogicalVolume(pmtGlassSectionV, //volume
					   															 G4Material::GetMaterial("Glass"),  //material
					   													 		 "PMT Glass Section Logical Volume"); //name

  // PMT Photocathode Section - Physical Volume
  new G4PVPlacement(rotm,                                            //rotation
		    						G4ThreeVector(pLightGuideLength/2 - pScintillatorHousingThickness/2 + pPmtPhotocathodeHeight/2, 0, 0), //translation
		    						fPmtGlassSectionLV,                    //logical volume
		    						"PMT Glass Section Physical Volume",        //name
		    						fScintillatorLGHousingLV,							//mother logical volume
		    						false,                                        //many
		    						0);

  // PMT Photocathode Section - Volume
  G4Tubs* pmtPhotocathodeSectionV = new G4Tubs("PMT Photocathode Section Volume", //name
					       															 0. * cm,                                   //inner radius
					       											 				 pPmtPhotocathodeRadius,           //outer radius
					       											 				 pPmtPhotocathodeHeight / 2.,          //half height
					       											 				 0. * deg,                                  //start angle
					       											 				 360. * deg);                               //end angle

  // PMT Photocathode Section - Logical Volume
  fPmtPhotocathodeSectionLV = new G4LogicalVolume(pmtPhotocathodeSectionV, //volume
						  																		G4Material::GetMaterial("Aluminium"),  //material
						  																		"PMT Photocathode Section Logical Volume", //name
						  																		0, //field manager
						  																		pPmtSD, //sensitive detector
						  																		0, //userlimits
						  																		true); //optimise

  // PMT Photocathode Section - Physical Volume
  new G4PVPlacement(0,                                            //rotation
		    						G4ThreeVector(0, 0,	0), //translation
		    						fPmtPhotocathodeSectionLV,                    //logical volume
		    						"PMT Photocathode Section Physical Volume",        //name
		    						fPmtGlassSectionLV,							//mother logical volume
		    						false,                                        //many
		    						0);

  VisAttributes();
  SurfaceProperties(pScintillatorHousingPT, pLGHousingPT, pPmtPT);
}

void MilliQDetectorBlockLV::VisAttributes() {
  fScintillatorHousingLV->SetVisAttributes(G4Colour::Red());
  fScintillatorLGHousingLV->SetVisAttributes(G4Colour::Red());
  fScintillatorLV->SetVisAttributes(G4Colour::Blue());
  fLightGuideLV->SetVisAttributes(G4Colour::Blue()); //(0.94, 0.5, 0.5)); // Pink
  fPmtGlassSectionLV->SetVisAttributes(G4Colour::Cyan());
  //   fPmtVacuumSectionLV->SetVisAttributes(G4Colour::Brown());
  fPmtPhotocathodeSectionLV->SetVisAttributes(G4Colour::Magenta());

  this->SetVisAttributes(G4Colour::Gray());
}

//Settings for BC408
void MilliQDetectorBlockLV::SurfaceProperties(G4MaterialPropertiesTable * pScintillatorHousingPT,
																							G4MaterialPropertiesTable * pLGHousingPT,
																							G4MaterialPropertiesTable * pPmtPT) {

  //**Scintillator housing properties
  G4OpticalSurface * OpScintHousingSurface = new G4OpticalSurface("HousingSurface");
  OpScintHousingSurface->SetType(dielectric_metal);
  OpScintHousingSurface->SetFinish(polished);
  OpScintHousingSurface->SetModel(glisur);
  OpScintHousingSurface->SetMaterialPropertiesTable(pScintillatorHousingPT);

  //**Light guide housing properties
  G4OpticalSurface* OpLGHousingSurface = new G4OpticalSurface("LGHousingSurface");
  OpLGHousingSurface->SetType(dielectric_metal);
  OpLGHousingSurface->SetFinish(polished);
  OpLGHousingSurface->SetModel(glisur);
  OpLGHousingSurface->SetMaterialPropertiesTable(pLGHousingPT);

  //**Photocathode surface properties
  //reference: http://www.hamamatsu.com/resources/pdf/etd/R329-02_TPMH1254E.pdf
  G4OpticalSurface * photocath_opsurf = new G4OpticalSurface("photocath_opsurf", glisur, polished, dielectric_metal);
  photocath_opsurf->SetMaterialPropertiesTable(pPmtPT);

  //**Create logical skin surfaces
  new G4LogicalSkinSurface("photocath_surf", fScintillatorHousingLV, OpScintHousingSurface);
  new G4LogicalSkinSurface("lg_surf", fScintillatorLGHousingLV, OpLGHousingSurface);
  new G4LogicalSkinSurface("photocath_surf", fPmtPhotocathodeSectionLV, photocath_opsurf);
}
