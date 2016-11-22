
#include "MilliQDetectorMessenger.hh"
#include "MilliQDetectorConstruction.hh"

#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"
#include "G4UIcmdWithAnInteger.hh"
#include "G4UIcommand.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4Scintillation.hh"


MilliQDetectorMessenger::MilliQDetectorMessenger(MilliQDetectorConstruction* detector)
 : fMilliQDetector(detector)
{
  //Setup a command directory for detector controls with guidance
  fDetectorDir = new G4UIdirectory("/MilliQ/detector/");
  fDetectorDir->SetGuidance("Detector geometry control");

  fVolumesDir = new G4UIdirectory("/MilliQ/detector/volumes/");
  fVolumesDir->SetGuidance("Enable/disable volumes");


  fMilliQCmd = new G4UIcmdWithABool("/MilliQ/detector/volumes/MilliQ",this);
  fMilliQCmd->SetGuidance("Enable/Disable the main detector volume.");
  fMilliQCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fMilliQCmd->SetToBeBroadcasted(false);

  fDefaultsCmd = new G4UIcommand("/MilliQ/detector/defaults",this);
  fDefaultsCmd->SetGuidance("Set all detector geometry values to defaults.");
  fDefaultsCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
  fDefaultsCmd->SetToBeBroadcasted(false);

  fMainScintYield=new G4UIcmdWithADouble("/MilliQ/detector/MainScintYield",this);
  fMainScintYield->SetGuidance("Set scinitillation yield of main volume.");
  fMainScintYield->SetGuidance("Specified in photons/MeV");
  fMainScintYield->AvailableForStates(G4State_PreInit,G4State_Idle);
  fMainScintYield->SetToBeBroadcasted(false);
/*
  fStepSizeCmd = new G4UIcmdWithADoubleAndUnit("/testex/det/setStepSize",this);
  fStepSizeCmd->SetGuidance("Set maxStepSize in the absorber");
  fStepSizeCmd->SetParameterName("StepSize",false);
  fStepSizeCmd->SetRange("StepSize>0.");
  fStepSizeCmd->SetUnitCategory("Length");
  fStepSizeCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fMagFieldCmd = new G4UIcmdWithADoubleAndUnit("/testex/det/setField",this);
  fMagFieldCmd->SetGuidance("Define magnetic field.");
  fMagFieldCmd->SetGuidance("Magnetic field will be in Z direction.");
  fMagFieldCmd->SetParameterName("Bz",false);
  fMagFieldCmd->SetUnitCategory("Magnetic flux density");
  fMagFieldCmd->AvailableForStates(G4State_PreInit,G4State_Idle);
*/
}


MilliQDetectorMessenger::~MilliQDetectorMessenger()
{
  delete fDetectorDir;
  delete fVolumesDir;
  delete fDefaultsCmd;
  delete fMilliQCmd;
  delete fMainScintYield;
}


void MilliQDetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == fDefaultsCmd){
    fMilliQDetector->SetDefaults();
  }
  else if (command == fMainScintYield){
   fMilliQDetector->SetMainScintYield(fMainScintYield->GetNewDoubleValue(newValue));
  }
}
