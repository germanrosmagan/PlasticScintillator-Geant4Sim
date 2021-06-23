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
/// \file runAndEvent//src/TrackingAction.cc
/// \brief Implementation of the TrackingAction class
//
// $Id: $
//
#include "PrimaryGeneratorGun2.hh"
#include "PrimaryGeneratorAction.hh"
#include "Run.hh"
#include "G4RunManager.hh"
#include "G4Track.hh"
#include "G4UnitsTable.hh"
#include "G4SteppingManager.hh"
#include "G4ProcessManager.hh"
#include "TrackingAction.hh"
//#include "Trajectory.hh"

#include <fstream>
using namespace std;
extern ofstream myfile;
extern ofstream myfile2;


TrackingAction::TrackingAction(PrimaryGeneratorGun2* prim):G4UserTrackingAction(),fPrimary(prim)
{
}

void TrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
        // myfile << G4endl << G4endl << "PRETRACKINGACTION " << G4endl;
	fpTrackingManager->SetStoreTrajectory(true);
}

void TrackingAction::PostUserTrackingAction(const G4Track* aTrack)
{
G4ThreeVector vertex = fPrimary->GetParticleGun()->GetParticlePosition();
//G4ThreeVector position = aTrack->GetPosition() - vertex;  
G4ThreeVector position = aTrack->GetPosition();   

//myfile5 << "TrackID " << aTrack->GetTrackID() << " process " << aTrack->GetCreatorProcess()->GetProcessName() << " secundarios " << fpTrackingManager->GimmeSecondaries()->size()<< G4endl;
	//geometrical information in G4Track is identical to "PostStepPoint"

	//true and projected ranges for primary particle
	
/*	G4ThreeVector vertex = fPrimary->GetParticleGun()->GetParticlePosition();
	G4ThreeVector position = aTrack->GetPosition() - vertex;     
	G4double energytot=  aTrack->GetTotalEnergy(); 
	G4double Kinenergy=aTrack->GetKineticEnergy();
	G4String particle2= aTrack->GetParticleDefinition()->GetParticleName();  		
	G4TrackVector * secondaries = fpTrackingManager->GimmeSecondaries();// Número de sencundarios en el Track
	myfile << "TrackID " << aTrack->GetTrackID()<< " proveniente de " << aTrack->GetParentID()<< " vertex " << vertex << "   position " << position  << " part. " << particle2  << "   KINenergy " << G4BestUnit(Kinenergy, "Energy") << " energytot " << G4BestUnit(energytot, "Energy") << " Num de sec en el Track " << secondaries->size() << G4endl;*/

	
/*	if (aTrack->GetParentID() != 0)
	{
		G4String process= aTrack->GetCreatorProcess()->GetProcessName();		
		myfile << " process " << process  <<  G4endl; 
		// if (process=="Scintillation") myfile  <<  G4endl <<  G4endl << " CENTELLEO  " <<  G4endl <<  G4endl;      
	}*/
}

