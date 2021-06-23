#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh" //This is so that you can use functions in the UI
#include "DetectorConstruction.hh" 
#include "G4StepLimiterPhysics.hh"
//#include "G4GeneralParticleSource.hh"
#include "RunAction.hh"
#include "EventAction.hh"
//#include "PrimaryGeneratorAction.hh"
#include "PrimaryGeneratorGun2.hh"
#include "SteppingAction.hh"
#include "TrackingAction.hh"
#include "PhysicsList.hh"
//#include "ActionInitialization.hh"
#include "G4TrajectoryDrawByParticleID.hh"
#include "QGSP_BERT.hh"//A short hand list of the physics things to simulate

#include <fstream>
using namespace std;
ofstream myfile("Output_surface.dat");
ofstream myfile2("Output_detector.dat");
ofstream myfile3("Output_track.dat");
ofstream myfile4("Output_ProcessNULL.dat");
ofstream myfile5("Output_stepping.dat");// Very heavy. Only for small tests. 
ofstream myfile6("Output_OpAbsorption.dat");
ofstream myfile7("Output_PrimaryParticle.dat");
ofstream myfile8("Output_TrueEvent.dat");

//The start of any C++ program is the main() function
int main(int argc,char** argv)
{

CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
//set random seed with system time
G4long seed = SSEEDD;
CLHEP::HepRandom::setTheSeed(seed);

/*The only manager which should be explicitly constructed in the user's main(). The run manager
 must be given all the information necessary:detector construction, particle and physics, primary 
partcile production and others (see **)*/
  G4RunManager* runManager = new G4RunManager();

//**// detector construction		
  DetectorConstruction* detector = new DetectorConstruction();
  runManager->SetUserInitialization(detector);

//**// particles and physics involved in the simulation
//G4VModularPhysicsList* physicsList = new QGSP_BERT;
//physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  //G4VUserPhysicsList* physics = new PhysicsList();
  //runManager->SetUserInitialization(physics);
 runManager->SetUserInitialization(new PhysicsList());

//**// how the primary particle is started?
//Start GEANT4 with the configuration the user has specficed
  runManager->Initialize();

  PrimaryGeneratorGun2* gen_action = new PrimaryGeneratorGun2();
  runManager->SetUserAction(gen_action);
  RunAction* runAct = new RunAction(detector,gen_action);
  runManager->SetUserAction(runAct);
  EventAction* evtAct = new EventAction();
  runManager->SetUserAction(evtAct);
  SteppingAction* stpAct = new SteppingAction(detector,runAct,evtAct);
  runManager->SetUserAction(stpAct);
  TrackingAction* trackact=new TrackingAction(gen_action);
  runManager->SetUserAction(trackact);

 // ActionInitialization* actionin  = new ActionInitialization(detector,gen_action);
 //runManager->SetUserInitialization(actionin);



//This is the visualization manager that allows graphics to be displayed
  G4VisManager* visManager = new G4VisExecutive();
  visManager->Initialize();
  G4TrajectoryDrawByParticleID* model = new G4TrajectoryDrawByParticleID;
  G4TrajectoryDrawByParticleID* model2 = new G4TrajectoryDrawByParticleID("test");

  //model->SetDefault("blue");
  model->Set("gamma", "yellow");
  model->Set("e+", "magenta");
  //model->Set("e-", G4Colour(0.3, 0.3, 0.3));
  model->Set("e-", "red");
  model->Set("mu+", "blue");
  model->Set("mu-", "brown");
  model->Set("pi+", "white");
  model->Set("pi-", "cyan");
  model->Set("proton", "grey");
  model->Set("alpha", "grey");
  model->Set("He3", "grey");
  model->Set("GenericIon", "grey");
 visManager->RegisterModel(model);
  visManager->RegisterModel(model2);

  visManager->SelectTrajectoryModel(model->Name());

//The User Interface (UI) manager so you may start a user interface
  G4UImanager* UImanager = G4UImanager::GetUIpointer();
 
  if (argc!=1) //if argc != 1 then we are in batch mode
  {
    //command line contains name of the macro to execute
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    
    //execute the command
    UImanager->ApplyCommand(command+fileName);
  }
  else //here we are ing interactive mode, so we need to define the UI session
  {
    G4UIExecutive * ui = new G4UIExecutive(argc,argv);

    //execute the visual macro to initialize the visuals you desire
    UImanager->ApplyCommand("/control/execute vis.mac");//VIEWER POPS-UP WHEN THE SCRIPT RUNS/**********************

//UImanager->ApplyCommand("/vis/open OGLSQt");
UImanager->ApplyCommand("/run/verbose 1");
UImanager->ApplyCommand("/event/verbose 1");
UImanager->ApplyCommand("/tracking/verbose 1");
UImanager->ApplyCommand("/vis/drawVolume");
UImanager->ApplyCommand("/vis/scene/add/trajectories");
//UImanager->ApplyCommand("/vis/modeling/trajectories/create/drawByCharge drawByParticleID");
/*UImanager->ApplyCommand("/vis/modeling/trajectories/create/drawByAttribute-0");
UImanager->ApplyCommand("/vis/modeling/trajectories/drawByAttribute-0/verbose true");
UImanager->ApplyCommand("/vis/modeling/trajectories/drawByAttribute-0/setAttribute CPN");
UImanager->ApplyCommand("/vis/modeling/trajectories/drawByAttribute-0/addValue brem_key  eBrem");
UImanager->ApplyCommand("/vis/modeling/trajectories/drawByAttribute-0/addValue annihil_key annihil");
UImanager->ApplyCommand("/vis/modeling/trajectories/drawByAttribute-0/addValue decay_key Decay");
UImanager->ApplyCommand("/vis/modeling/trajectories/drawByAttribute-0/addValue muIon_key muIoni");
UImanager->ApplyCommand("/vis/modeling/trajectories/drawByAttribute-0/addValue eIon_key  eIoni");

UImanager->ApplyCommand("/vis/modeling/trajectories/drawByAttribute-0/brem_key/setLineColour     red");
UImanager->ApplyCommand("/vis/modeling/trajectories/drawByAttribute-0/annihil_key/setLineColour  green");
UImanager->ApplyCommand("/vis/modeling/trajectories/drawByAttribute-0/decay_key/setLineColour    cyan");
UImanager->ApplyCommand("/vis/modeling/trajectories/drawByAttribute-0/eIon_key/setLineColour     yellow");
UImanager->ApplyCommand("/vis/modeling/trajectories/drawByAttribute-0/muIon_key/setLineColour magenta");*/
UImanager->ApplyCommand("/vis/scene/endOfEventAction accumulate");

/*
UImanager->ApplyCommand("/gps/particle  proton");
UImanager->ApplyCommand("/gps/pos/type Surface");
UImanager->ApplyCommand("/gps/pos/shape Sphere");
UImanager->ApplyCommand("/gps/pos/radius 10 cm");
UImanager->ApplyCommand("/gps/pos/centre 0. 0. 0. cm");
UImanager->ApplyCommand("/gps/ang/type cos");
UImanager->ApplyCommand("/gps/ang/mintheta 0 deg");
UImanager->ApplyCommand("/gps/ang/maxtheta 90 deg");
UImanager->ApplyCommand("/gps/ene/mono 10 MeV");
*/

UImanager->ApplyCommand("/run/beamOn 1");

    //start the user interface and keep it running until you exit
    // COMMENTED SO THAT THE VIEWER DOES NOT EXIT WHEN THE SCRIPT IS EXECUTED / *********************
    //ui->SessionStart(); 

    delete ui;
  }
///////////////////////////////////////////////////////////////////////////////

//Delete runManager to delete user actions, physics_list and detector_description
  delete runManager; 

  return 0; 
}
