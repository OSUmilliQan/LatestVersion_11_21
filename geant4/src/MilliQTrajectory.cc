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
// $Id: MilliQTrajectory.cc 72349 2013-07-16 12:13:16Z gcosmo $
//
/// \file optical/MilliQ/src/MilliQTrajectory.cc
/// \brief Implementation of the MilliQTrajectory class
//
//
#include "MilliQTrajectory.hh"
#include "G4TrajectoryPoint.hh"
#include "G4Trajectory.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleTypes.hh"
#include "G4ThreeVector.hh"
#include "G4Polyline.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4VVisManager.hh"
#include "G4Polymarker.hh"

G4ThreadLocal G4Allocator<MilliQTrajectory>* MilliQTrajectoryAllocator = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MilliQTrajectory::MilliQTrajectory() : G4Trajectory(), fDrawit(false), fForceNoDraw(false), fForceDraw(false) {
  fParticleDefinition = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MilliQTrajectory::MilliQTrajectory(const G4Track* aTrack) : G4Trajectory(aTrack),fDrawit(false) {
  fParticleDefinition = aTrack->GetDefinition();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MilliQTrajectory::MilliQTrajectory(MilliQTrajectory &right) : G4Trajectory(right), fDrawit(right.fDrawit) {
  fParticleDefinition = right.fParticleDefinition;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MilliQTrajectory::~MilliQTrajectory() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MilliQTrajectory::DrawTrajectory() const
{
  // i_mode is no longer available as an argument of G4VTrajectory.
  // In this exampple it was always called with an argument of 50.
  const G4int i_mode = 50;
  // Consider using commands /vis/modeling/trajectories.
  
  //Taken from G4VTrajectory and modified to select colours based on particle
  //type and to selectively eliminate drawing of certain trajectories.

  if(!fForceDraw && (!fDrawit || fForceNoDraw))
    return;

  // If i_mode>=0, draws a trajectory as a polyline and, if i_mode!=0,
  // adds markers - yellow circles for step points and magenta squares
  // for auxiliary points, if any - whose screen size in pixels is
  // given by std::abs(i_mode)/1000.  E.g: i_mode = 5000 gives easily
  // visible markers.
 
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if (!pVVisManager) return;
 
  const G4double markerSize = std::abs(i_mode)/1000;
  G4bool lineRequired (i_mode >= 0);
  G4bool markersRequired (markerSize > 0.);
 
  G4Polyline trajectoryLine;
  G4Polymarker stepPoints;
  G4Polymarker auxiliaryPoints;
 
  for (G4int i = 0; i < GetPointEntries() ; i++) {
    G4VTrajectoryPoint* aTrajectoryPoint = GetPoint(i);
    const std::vector<G4ThreeVector>* auxiliaries
      = aTrajectoryPoint->GetAuxiliaryPoints();
    if (auxiliaries) {
      for (size_t iAux = 0; iAux < auxiliaries->size(); ++iAux) {
        const G4ThreeVector pos((*auxiliaries)[iAux]);
        if (lineRequired) {
          trajectoryLine.push_back(pos);
        }
        if (markersRequired) {
          auxiliaryPoints.push_back(pos);
        }
      }
    }
    const G4ThreeVector pos(aTrajectoryPoint->GetPosition());
    if (lineRequired) {
      trajectoryLine.push_back(pos);
    }
    if (markersRequired) {
      stepPoints.push_back(pos);
    }
  }
 
  if (lineRequired) {
    G4Colour colour;
 
    if(fParticleDefinition==G4OpticalPhoton::OpticalPhotonDefinition()){
        //Scintillation and Cerenkov photons are green
        colour = G4Colour(0.,1.,0.);
    }
    else if(fParticleDefinition==G4Electron::ElectronDefinition()){
        //red
        colour = G4Colour(1.,0.,0.);
    }
    else if(fParticleDefinition==G4Gamma::GammaDefinition()){
        //cyan
        colour = G4Colour(0.,1.,1.);
    }
    else //All other particles are yellow
        colour = G4Colour(1.,1.,0.);
 
    G4VisAttributes trajectoryLineAttribs(colour);
    trajectoryLine.SetVisAttributes(&trajectoryLineAttribs);
    pVVisManager->Draw(trajectoryLine);
  }
  if (markersRequired) {
    auxiliaryPoints.SetMarkerType(G4Polymarker::squares);
    auxiliaryPoints.SetScreenSize(markerSize);
    auxiliaryPoints.SetFillStyle(G4VMarker::filled);
    G4VisAttributes auxiliaryPointsAttribs(G4Colour(0.,1.,1.));  // Magenta
    auxiliaryPoints.SetVisAttributes(&auxiliaryPointsAttribs);
    pVVisManager->Draw(auxiliaryPoints);

    stepPoints.SetMarkerType(G4Polymarker::circles);
    stepPoints.SetScreenSize(markerSize);
    stepPoints.SetFillStyle(G4VMarker::filled);
    G4VisAttributes stepPointsAttribs(G4Colour(1.,1.,0.));  // Yellow
    stepPoints.SetVisAttributes(&stepPointsAttribs);
    pVVisManager->Draw(stepPoints);
  }
}
