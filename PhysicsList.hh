#ifndef PhysicsList_h
#define PhysicsList_h 1

// Defines mandatory user class PhysicsList

#include "G4VUserPhysicsList.hh"
#include "globals.hh"
#include "TrackingAction.hh"
#include "SteppingAction.hh"
#include "PrimaryGeneratorGun2.hh"
#include "G4Track.hh"
//Class prototype
class G4VPhysicsConstructor;

/*
   This mandatory user class provides the physics

   It is responsible for
    - Definition of particles
    - Construction of physics processes
    - setting of user cuts (limits where you throw out physics below a certain energy)
*/

class PhysicsList: public G4VUserPhysicsList
{
protected:
  // Construct particles
  void ConstructParticle();
  void ConstructBosons();
  void ConstructLeptons();
  void ConstructMesons();
  void ConstructBaryons();
  void ConstructIons();
  void ConstructNuclei();         
  // Construct physics processes
  void ConstructProcess();
  void ConstructEM();
  void ConstructGeneral();  
  //void ConstructOp(TrackingAction(fPrimary));
  void ConstructScintillation();
  // Define user cuts
  void SetCuts();

private:
  G4VPhysicsConstructor*  emPhysicsList;
  SteppingAction*  fstep;
 

public:
  //Constructor
  //PhysicsList();
  PhysicsList():G4VUserPhysicsList(){}
  PhysicsList(SteppingAction*);
  //Destructor
  ~PhysicsList();
 virtual void ConstructOp();

  G4Step* aStep;
 //PrimaryGeneratorGun2*  fPrimary;
};
#endif

 
