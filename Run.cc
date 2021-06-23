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
/// \brief Implementation of the Run class
//
// $Id: Run.cc 71376 2013-06-14 07:44:50Z maire $
#include "Run.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorGun2.hh"
#include "PrimaryGeneratorAction.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "G4EmCalculator.hh"

#include <iomanip>
#include <fstream>
  using namespace std;
  extern ofstream myfile;

Run::Run(DetectorConstruction* det): G4Run(),fDetector(det),fParticle(0),fEkin(0.),fNbOfTraks0(0),fNbOfTraks1(0),fNbOfSteps0(0),fNbOfSteps1(0),
  fEdep(0),fTrueRange(0.), fTrueRange2(0.),fProjRange(0.),fProjRange2(0.),fTransvDev(0.), fTransvDev2(0.)
{
}

Run::~Run()
{
}

void Run::SetPrimary(G4ParticleDefinition* particle, G4double energy)
{ 
  fParticle = particle;
  fEkin = energy;
} 

void Run::PrintSummary() const
{    
  G4int nbOfEvents     = GetNumberOfEvent();
  G4String partName    = fParticle->GetParticleName(); 
  G4Material* material = fDetector->GetMaterial();
  G4double density     = material->GetDensity();
     
   // 24/10/2016  myfile << " The run was: " << nbOfEvents << " " << partName << " of "<< G4BestUnit(fEkin,"Energy") << " through " << " of "
          // 24/10/2016  << material->GetName() << " (density: " << G4BestUnit(density,"Volumic Mass") << ")" << G4endl;           
}   
void Run::Merge(const G4Run* run)
{
  const Run* localRun = static_cast<const Run*>(run);

  // pass information about primary particle
  fParticle = localRun->fParticle;
  fEkin     = localRun->fEkin;

  // accumulate sums
  //
  fNbOfTraks0 += localRun->fNbOfTraks0;  
  fNbOfTraks1 += localRun->fNbOfTraks1;  
  fNbOfSteps0 += localRun->fNbOfSteps0;
  fNbOfSteps1 += localRun->fNbOfSteps1;   
  fEdep       += localRun->fEdep;  
  fTrueRange  += localRun->fTrueRange;
  fTrueRange2 += localRun->fTrueRange2;
  fProjRange  += localRun->fProjRange;
  fProjRange2 += localRun->fProjRange2;
  fTransvDev  += localRun->fTransvDev;
  fTransvDev2 += localRun->fTransvDev2; 
  
  G4Run::Merge(run); 
} 

