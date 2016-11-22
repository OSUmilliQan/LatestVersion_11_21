//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: MilliQSteppingAction.hh 68752 2013-04-05 10:23:47Z gcosmo $
//
/// \file optical/MilliQ/include/MilliQSteppingAction.hh
/// \brief Definition of the MilliQSteppingAction class
//
#ifndef MilliQSteppingAction_H
#define MilliQSteppingACtion_H 1

#include "globals.hh"
#include "G4UserSteppingAction.hh"

#include "G4OpBoundaryProcess.hh"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>

class MilliQRecorderBase;
class MilliQEventAction;
class MilliQDetectorConstruction;
class MilliQTrackingAction;
class MilliQSteppingMessenger;

class MilliQSteppingAction : public G4UserSteppingAction {

public:

  MilliQSteppingAction(MilliQRecorderBase*, const G4int geometryVersion);
  virtual ~MilliQSteppingAction();
  virtual void UserSteppingAction(const G4Step*);

  void SetOneStepPrimaries(G4bool b) { fOneStepPrimaries = b; }
  G4bool GetOneStepPrimaries() { return fOneStepPrimaries; }

private:

  MilliQRecorderBase* fRecorder;
  G4bool fOneStepPrimaries;
  MilliQSteppingMessenger* fSteppingMessenger;

  const G4int fAlternate;

  G4OpBoundaryProcessStatus fExpectedNextStatus;
};

#endif
