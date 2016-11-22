#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

#include "G4UImanager.hh"
#include "G4String.hh"
#include "G4Types.hh"

#include "MilliQPhysicsList.hh"
#include "MilliQDetectorConstruction.hh"
#include "G4PhysListFactory.hh"
#include "G4VModularPhysicsList.hh"

#include "MilliQActionInitialization.hh"
#include "MilliQRecorderBase.hh"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

int main(int argc, char** argv) {

#ifdef G4MULTITHREADED
  G4MTRunManager * runManager = new G4MTRunManager;
#else
  G4RunManager * runManager = new G4RunManager;
#endif

  boost::property_tree::ptree pt;
  boost::property_tree::ptree pt_geometry;

  // Read primary configuration file
  std::string configFile = (argc >= 3) ? argv[2] : "config/default.ini";
  try {
    boost::property_tree::ini_parser::read_ini(configFile, pt); // std::string, ptree
    boost::property_tree::ini_parser::read_ini(pt.get<std::string>("Configuration.GeometryConfigFile"), pt_geometry);
  }
  catch(boost::property_tree::ptree_error &e) {
    G4ExceptionDescription msg;
    msg << G4endl << "Configuration file " << e.what() << G4endl;
    G4Exception("MilliQ::main()", "MilliQ::ConfigFileReadError", FatalException, msg);
  }

  const G4int geometryVersion = pt.get<G4int>("Configuration.Version");

  runManager->SetUserInitialization(new MilliQDetectorConstruction(pt));
  runManager->SetUserInitialization(new MilliQPhysicsList(pt));

  MilliQRecorderBase* recorder = NULL; //No recording is done in this example
  runManager->SetUserInitialization(new MilliQActionInitialization(recorder, geometryVersion, pt));

#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive("Quiet");
#endif

  runManager->Initialize();
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if(argc == 1){
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    //UImanager->ApplyCommand("/control/macroPath /Users/JamesLondon/geant/geant-builds/MilliQ");
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute vis.mac");
#endif
    if(ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;
#endif
  }

  else {
    G4String command = "/control/execute ";
    G4String filename = argv[1];
    UImanager->ApplyCommand(command+filename);
  }

#ifdef G4VIS_USE
  delete visManager;
#endif

  delete runManager;
  return 0;
}
