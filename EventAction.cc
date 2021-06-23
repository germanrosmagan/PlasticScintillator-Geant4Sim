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
// $Id: EventAction.cc 75117 2013-10-28 09:38:37Z gcosmo $
//
/// \file EventAction.cc
/// \brief Implementation of the EventAction class

#include "EventAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

#include <fstream>
using namespace std;
extern ofstream myfile;


EventAction::EventAction(): G4UserEventAction()
{
} 

EventAction::~EventAction()
{}

void EventAction::BeginOfEventAction(const G4Event*)
{ 
  // myfile << G4endl<< G4endl << "EVENTACTION.Begin() " << G4endl;
}

void EventAction::EndOfEventAction(const G4Event* event)
{   
  // myfile << G4endl<< "EVENTACTION.End() " << G4endl;
  G4TrajectoryContainer* trajectoryContainer = event->GetTrajectoryContainer();
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
  // periodic printing
  //#G4Step* theStep;
  //#G4int* isPrimary = theStep->GetTrack()->GetParrentID();
  //theStep->GetTrack()->GetParentID() == 0  means the partcile is secondary
  G4int eventID = event->GetEventID();
  //if ( eventID < 100 || eventID % 100 == 0) {
     // 24/10/2016 myfile << "END OF eventID= "  << eventID  << "   n_trajectories= "  << n_trajectories << endl;
     // 24/10/2016  myfile << endl << endl;
  //if ( trajectoryContainer ) {
  //    G4cout << "    " << n_trajectories << " trajectories stored in this event." << G4endl;
  //  }
  //}
}
