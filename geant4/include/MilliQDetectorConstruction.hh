#ifndef MilliQDetectorConstruction_H
#define MilliQDetectorConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Box;
class G4Tubs;
class MilliQDetectorStack;
class G4UniformMagField;
class G4GlobalMagFieldMessenger;
//class MilliQMonopoleFieldSetup;

#include "G4Cache.hh"
#include "G4Material.hh"
#include "MilliQDetectorMessenger.hh"
#include "G4VisAttributes.hh"
#include "G4RotationMatrix.hh"
#include "MilliQMonopoleFieldSetup.hh"

#include "MilliQPMTSD.hh"
#include "MilliQScintSD.hh"
#include "MilliQDetectorStackParameterisation.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4Cache.hh"
#include "globals.hh"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

class MilliQDetectorConstruction : public G4VUserDetectorConstruction {

public:

  MilliQDetectorConstruction(const boost::property_tree::ptree pt);
  virtual ~MilliQDetectorConstruction() {};

  virtual G4VPhysicalVolume* Construct();

  virtual void ConstructSDandField();
  void ConstructShield(G4LogicalVolume*, G4double, G4double);

  void ConstructCheckGeometry();

  void ReadConfiguration();
  void ReadGeometryConfiguration();
  void ReadScintillatorConfiguration();
  void ReadScintillatorHousingConfiguration();
  void ReadPMTConfiguration();

  //Functions to modify the geometry
  void SetMagField(G4double);
  void SetHousingThickness(G4double );
  void SetPMTRadius(G4double );
  void SetDefaults();

  //Get values
  G4double GetHousingThickness(){return fD_mtl;}
  G4double GetPMTRadius(){return fOuterRadius_pmt;}
  G4int GetNblocksPerStack() { return NBlocks.x() * NBlocks.y() * NBlocks.z(); }
  G4int GetNstacks(){return NStacks;}
  G4int GetAlternateGeometry(){return fAlternate;}
  inline const G4Material* GetScintMaterial() {return fScintillatorMaterial;};

  void SetHousingReflectivity(G4double );
  G4double GetHousingReflectivity(){return fRefl;}

  void SetMainScintYield(G4double );
  void DefineMaterials();

private:

  G4VPhysicalVolume* ConstructDetector();

  MilliQDetectorMessenger* fDetectorMessenger;

  G4VPhysicalVolume* fWorldPV;
  G4VPhysicalVolume* fWorldPVCheckPhysics;
  G4LogicalVolume* worldLV;
  G4LogicalVolume* fMagneticVolume;

  //Materials & Elements
  G4Material* fScintillatorMaterial;
  G4Material* fAluminiumMaterial;
  G4Material* fAirMaterial;
  G4Material* fGlassMaterial;
  G4Material* fVacuumMaterial;
  G4Material* fConcreteMaterial;
  G4Material* led;
  G4Material* polyethylene;

  //Geometry
  G4double fScint_x;
  G4double fScint_y;
  G4double fScint_z;
  G4double fD_mtl;
  G4double fScintHouseThick;
  G4double fLightGuideLength;

  G4double fScintDensity;
  G4int fScintCarbonContent;
  G4int fScintHydrogenContent;

  G4double fScintillationYield;
  G4double fScinitResolutionScale;
  G4double fScintFastTimeConstant;
  G4double fScintFastRiseTime;
  G4double fScintSlowTimeConstant;
  G4double fScintYieldRatio;
  G4double fScintBirksConstant;

  G4double fScintHousingRefl;
  G4double fScintHousingEff;

  G4double fLGHousingRefl;
  G4double fLGHousingEff;

  G4double fPmtRad;
  G4double fPmtPhotoRad;
  G4double fPmtPhotoHeight;
  //G4double fLGHouseRefl;
  G4double fOuterRadius_pmt;
  G4double fRefl;

  G4double fPmtReR;
  G4double fPmtImR;

  G4double NStacks;
  G4ThreeVector fOffset;
  G4ThreeVector shield1Thick;
  G4ThreeVector shield2Thick;
  G4ThreeVector detShieldGap;
  G4ThreeVector NBlocks;
  G4ThreeVector fBetweenBlockSpacing;
  G4int fAlternate;
  static G4ThreadLocal G4GlobalMagFieldMessenger*  fMagFieldMessenger;

//  G4UniformMagField* fMagField;
  MilliQDetectorStack* fDetectorStack;
  MilliQMonopoleFieldSetup* fMonFieldSetup;

  G4Cache<MilliQMonopoleFieldSetup*> fEmFieldSetup;

  G4MaterialPropertyVector * fScintRelativeOutput;
  G4MaterialPropertyVector * fScintRIndex;
  G4MaterialPropertyVector * fScintAbsLength;

  G4MaterialPropertyVector * fScintHousingReflectivity;
  G4MaterialPropertyVector * fScintHousingEfficiency;

  G4MaterialPropertyVector * fLGHousingReflectivity;
  G4MaterialPropertyVector * fLGHousingEfficiency;

  G4MaterialPropertyVector * fPmtEfficiency;
  G4MaterialPropertyVector * fPmtRealRIndex;
  G4MaterialPropertyVector * fPmtImaginaryRIndex;

  G4MaterialPropertyVector * fGlassRIndex;
  G4MaterialPropertyVector * fGlassAbsLength;

  G4MaterialPropertyVector * fVacuumRIndex;

  G4MaterialPropertiesTable* fScintillator_mt;

  G4MaterialPropertiesTable * fScintillatorHousingPT;
  G4MaterialPropertiesTable * fLightGuideHousingPT;
  G4MaterialPropertiesTable * fPhotoCathodePT;

  //Sensitive Detectors
  G4Cache<MilliQPMTSD*> fPmt_SD;
  G4Cache<MilliQScintSD*> fScint_SD;

  // Configuration
  const boost::property_tree::ptree fPTree;
  boost::property_tree::ptree fGeometryPTree;
  boost::property_tree::ptree fScintillatorPTree;
  boost::property_tree::ptree fPMTPTree;

};

#endif
