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
// $Id: RunAction.cc 75216 2013-10-29 16:08:11Z gcosmo $
//
/// \file RunAction.cc
/// \brief Implementation of the RunAction class

#include "RunAction.hh"
#include "Run.hh"
#include "PrimaryGeneratorGun2.hh"
#include "PrimaryGeneratorAction.hh"
#include "DetectorConstruction.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4ProcessManager.hh"

//#include <fstream>
//using namespace std;
//extern ofstream myfile;

RunAction::RunAction(DetectorConstruction* det, PrimaryGeneratorGun2* prim): G4UserRunAction(),fDetector(det), fPrimary(prim),fRun(0)
//(DetectorConstruction* det, PrimaryGeneratorAction* kin):G4UserRunAction(),fDetector(det),fPrimary(kin),fRun(0),fHistoManager(0)
{ 
G4RunManager::GetRunManager()->SetPrintProgress(1);         
}

RunAction::~RunAction()
{
}
G4Run* RunAction::GenerateRun()
{ 
  fRun = new Run(fDetector); 
  return fRun;
}
void RunAction::BeginOfRunAction(const G4Run* run)
{ 
 // myfile << G4endl << G4endl<< "RUNACTION.Begin()" << G4endl << G4endl;
  //inform the runManager to save random number seed
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  G4int runID = run->GetRunID();
 // myfile << "RunID=" <<runID << endl << endl;
}

void RunAction::EndOfRunAction(const G4Run* run)
{
 //myfile<< " runID " << runID  << " particle " << particle  << " energy " << energy << G4endl;
  G4int NbOfEvents = run->GetNumberOfEvent();
  if (NbOfEvents == 0) return;

 //G4Material* material = fDetector->GetMaterial();  
 //G4ParticleDefinition* particle = fPrimary->GetParticleGun()->GetParticleDefinition(); 
 //G4double energy = fPrimary->GetParticleGun()->GetParticleEnergy();
 //G4double charge   = particle->GetPDGCharge();

 // myfile << "END OF RUN " << run->GetRunID() << ". FINISH" << G4endl;
}
