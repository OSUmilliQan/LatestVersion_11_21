//
// ********************************************************************
// * License and Disclaimer                                          *
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
// $Id: MilliQPrimaryGeneratorAction.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file MilliQPrimaryGeneratorAction.cc
/// \brief Implementation of the MilliQPrimaryGeneratorAction class

#include "MilliQPrimaryGeneratorAction.hh"
#include "MilliQPrimaryGeneratorMessenger.hh"
#include "MilliQMonopole.hh"
#include "MilliQMonopolePhysics.hh"

#include "G4PhysicalConstants.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4Box.hh"
#include "G4LorentzVector.hh"
#include "G4Event.hh"
#include <vector>
#include <math.h>
#include "G4LorentzVector.hh"
#include "fstream"
#include "iostream"
#include "sstream"
#include <string>
#include "vector"
#include "globals.hh"
#include "G4ParticleGun.hh"

#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4Electron.hh"

#include "Randomize.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int MilliQPrimaryGeneratorAction::neventLHE = 0;
G4bool MilliQPrimaryGeneratorAction::firstPass = true;
std::vector<std::vector<G4double> >	MilliQPrimaryGeneratorAction::vertexList;
std::vector<std::vector<G4double> >	MilliQPrimaryGeneratorAction::momentumList;
std::vector<std::vector<G4double> >	MilliQPrimaryGeneratorAction::qmeList;

MilliQPrimaryGeneratorAction::MilliQPrimaryGeneratorAction(const boost::property_tree::ptree pt) :
  G4VUserPrimaryGeneratorAction(), fParticleGun(0), fGunMessenger(0), fRndmFlag("off"),
  fXVertex(0.), fYVertex(0.), fZVertex(0.), fEnergy(0.), fVertexDefined(false), fCalibDefined(false), fPTree(pt) {

  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);

  // create a messenger for this class
  fGunMessenger = new MilliQPrimaryGeneratorMessenger(this);

  // default particle kinematic

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  std::string particleName;

  G4ParticleDefinition* particle = particleTable->FindParticle(particleName = "monopole");

  fParticleGun->SetParticleDefinition(particle);
  fgPrimaryParticle = particle;

  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1., 0., 0.));

  fXVertex = -990. * cm;
  fYVertex = 0.;
  fZVertex = 0.;
  fEnergy = 10. * GeV;
  fParticleGun->SetParticlePosition(G4ThreeVector(fXVertex, fYVertex, fZVertex));
  fParticleGun->SetParticleEnergy(fEnergy);

  try {
    boost::property_tree::ini_parser::read_ini(fPTree.get<std::string>("Configuration.ParticleConfigFile"), fParticlePTree);
  }
  catch(boost::property_tree::ptree_error &e) {
    G4ExceptionDescription msg;
    msg << G4endl << "Configuration file " << e.what() << G4endl;
    G4Exception("MilliQPrimaryGeneratorAction::PrimaryGeneratorAction()", "ConfigFileReadError", FatalException, msg);
  }

  fPathname = fParticlePTree.get<G4String>("ParticleProperties.PathName");
  fFilename = fParticlePTree.get<G4String>("ParticleProperties.FileName");


}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

MilliQPrimaryGeneratorAction::~MilliQPrimaryGeneratorAction() {
  delete fParticleGun;
  delete fGunMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MilliQPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent) {

  //Load in the data
  if(firstPass == true){
    GetLHEFourVectors();
    firstPass = false;
  }

  G4double x0, y0, z0, xMo, yMo, zMo, En, MoNorm;
  if (fVertexDefined) {
    x0 = fXVertex;
    y0 = fYVertex;
    z0 = fZVertex;
    En = fEnergy;
    xMo = 1.;
    yMo = -0.1 + 0.2 * G4UniformRand();
    zMo = 0.;
  }
  else {
    xMo = momentumList[neventLHE][0];
    yMo = momentumList[neventLHE][1];
    zMo = momentumList[neventLHE][2];
    x0 = 0.00*m;//vertexList[neventLHE][0]*m;
    y0 = vertexList[neventLHE][1]*m;
    z0 = vertexList[neventLHE][2]*m;
    En = qmeList[neventLHE][2]*GeV;
  }

  MoNorm = sqrt(pow(xMo, 2) + pow(yMo, 2) + pow(zMo, 2));

  if(neventLHE % 100 == 0) G4cout <<"neventLHE "<<neventLHE<<G4endl;

  //	G4cout << "xMoGun " << xMo << " yMoGun " << yMo << " zMoGun " << zMo << " xGun " << x0 << " yGun " << y0 << " zGun " << z0 << G4endl;
  //	G4cout << "nevent " << neventLHE << " Q " << qmeList[neventLHE][0] << " M " << qmeList[neventLHE][1] << " E "<< qmeList[neventLHE][2] << G4endl;

  fParticleGun->SetParticleMomentumDirection(	G4ThreeVector(xMo / MoNorm, yMo / MoNorm, zMo / MoNorm));
  fParticleGun->SetParticleEnergy(En);
  fParticleGun->SetParticlePosition(G4ThreeVector(x0, y0, z0));
  fParticleGun->GeneratePrimaryVertex(anEvent);
  fParticleGun->SetParticleDefinition(fgPrimaryParticle);
  neventLHE++;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4String MilliQPrimaryGeneratorAction::GetPrimaryName() {
  return fgPrimaryParticle->GetParticleName();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MilliQPrimaryGeneratorAction::SetXVertex(G4double x) {
  fVertexDefined = true;
  fXVertex = x;
  G4cout << " X coordinate of the primary vertex = " << fXVertex / mm << " mm." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MilliQPrimaryGeneratorAction::SetYVertex(G4double y) {
  fVertexDefined = true;
  fYVertex = y;
  G4cout << " Y coordinate of the primary vertex = " << fYVertex / mm << " mm." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MilliQPrimaryGeneratorAction::SetZVertex(G4double z) {
  fVertexDefined = true;
  fZVertex = z;
  G4cout << " Z coordinate of the primary vertex = " << fZVertex / mm << " mm." << G4endl;
}

void MilliQPrimaryGeneratorAction::SetCalibEnergy(G4double e) {
  fCalibDefined = true;
  fEnergy = e;
  G4cout << " Energy of the primary vertex = " << fEnergy / GeV << " GeV." << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void MilliQPrimaryGeneratorAction::GetLHEFourVectors() {
  std::ifstream infile;

  std::cout << std::endl << "opening " << fPathname << fFilename << std::endl;

  infile.open(fPathname.append(fFilename).c_str());
  G4String line;
  G4double fe, fq, fm, fx, fy, fz, fpx, fpy, fpz;
  std::vector<G4double> Tver(3), Tmo(3), Tqme(3);
  while (std::getline(infile, line)) {

    std::istringstream iss(line);
    iss.clear();
    iss.str(line);

    if (!(iss >> fq >> fm >> fx >> fy >> fz >> fpx >> fpy >> fpz))
      break;

    fe = std::sqrt(std::pow(fpx,2)+std::pow(fpy,2)+std::pow(fpz,2)+std::pow(fm,2))-fm;

    Tver[0]=fx;	Tver[1]=fy;	Tver[2]=fz;
    Tmo[0]=fpx;	Tmo[1]=fpy;	Tmo[2]=fpz;
    Tqme[0]=fq;	Tqme[1]=fm;	Tqme[2]=fe;
    MilliQPrimaryGeneratorAction::vertexList.push_back( Tver );
    momentumList.push_back( Tmo );
    qmeList.push_back( Tqme );

  }

  std::cout << std::endl << "found " << momentumList.size() << " events" << std::endl;

}
