#ifndef MilliQRunAction_h
#define MilliQRunAction_h 1

#include "MilliQDataFormat.hh"

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4Types.hh"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

class MilliQRecorderBase;
class MilliQDetectorConstruction;
class MilliQPrimaryGeneratorAction;
class G4Run;
class MilliQRunActionMessenger;

class MilliQRunAction : public G4UserRunAction {
public:

  MilliQRunAction(MilliQRecorderBase*, const boost::property_tree::ptree pt);
  virtual ~MilliQRunAction();

  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);

  inline void SetVerbose(G4int verbose) { fVerboseLevel = verbose;}
  inline G4int GetVerbose()             { return fVerboseLevel;}

  void InitializeDetectorParameters(const std::string filename);
  
private:

  MilliQRecorderBase*           fRecorder;
  MilliQDetectorConstruction*   fDetector;
  MilliQPrimaryGeneratorAction* fKinematic;
  MilliQRunActionMessenger*     fMilliQRunActionMessenger;

  MilliQDataFormat* dataFormat;
  
  G4int fVerboseLevel;

};


#endif
