//Include the prototypes of the members
#include "PhysicsList.hh"
#include "globals.hh"
//units
#include "G4SystemOfUnits.hh"

//Manages all the processes that you active/allow to happen
#include "G4ProcessManager.hh"

//the list of all particle types in Geant4
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh" 
#include "G4ParticleTypes.hh"
#include "G4LeptonConstructor.hh"
#include "G4BosonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"   
#include "G4GenericIon.hh"
#include "G4Deuteron.hh"
#include "G4Triton.hh"
#include "G4He3.hh"
#include "G4Alpha.hh"
#include "G4OpticalPhoton.hh"

//the list of processes

#include "G4ProcessManager.hh"
#include "G4ProcessVector.hh"
#include "G4VProcess.hh"

#include "G4EmStandardPhysics.hh"//Electromagnetic physics methods

#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4RayleighScattering.hh"
#include "G4CoulombScattering.hh"
#include "G4GammaConversionToMuons.hh"
#include "G4GammaConversion.hh"
#include "G4eMultipleScattering.hh"
#include "G4MuMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"
#include "G4MuIonisation.hh"
#include "G4MuBremsstrahlung.hh"
#include "G4MuPairProduction.hh"
#include "G4hIonisation.hh"
#include "G4hBremsstrahlung.hh"
#include "G4hPairProduction.hh"
#include "G4ionIonisation.hh"
#include "G4NuclearStopping.hh"
#include "G4IonParametrisedLossModel.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"

#include "G4Decay.hh"
#include "G4SynchrotronRadiation.hh"
#include "G4Cerenkov.hh"
#include "G4Scintillation.hh"
#include "G4EmSaturation.hh"

//Helper and options
#include "G4PhysicsListHelper.hh"
#include "G4PhysicsListOrderingParameter.hh" 
#include "G4EmProcessOptions.hh"
#include "G4LossTableManager.hh"//IDK what loss table man does

#include "SteppingAction.hh"
#include "G4Step.hh" 


//Constructor
//two mandatory virtual methods must be implemented: ConstructParticle() and ConstructProcess()
PhysicsList::PhysicsList(SteppingAction *stepp):G4VUserPhysicsList(), fstep(stepp)
{

	ConstructParticle();
	ConstructProcess();

	//Set the Cut value to 0.1 micrometers this item is in the G4VUserPhysicsList so it is inherited
	defaultCutValue = 0.1*um;

	SetVerboseLevel(1);//0:silent. 1:only the valid commands are shown. 2:comment lines are also shown (default).

	// Only allow EM physics.
	//emPhysicsList = new G4EmStandardPhysics();//Lo comento, luego no se usa

}

//Generic Destructor
PhysicsList::~PhysicsList() {}

void PhysicsList::ConstructParticle()
{
	//IMPORTANT: we must define ALL particles ussed in our application, not only primary particles (also secondaries
	//created by physical processes)
	ConstructBosons();
	ConstructLeptons();
	ConstructMesons();
	ConstructBaryons();
	ConstructIons();         
	ConstructNuclei(); 

}

void PhysicsList::ConstructBosons()//photons, phonons, bosons W and Z, gluons, Higgs. We only need photons.
{
	// gamma and optical photon
	G4Gamma::GammaDefinition(); 
	G4OpticalPhoton::OpticalPhotonDefinition(); 
}

void PhysicsList::ConstructLeptons()
{
	// leptons
	G4Electron::ElectronDefinition();
	G4Positron::PositronDefinition();
	G4NeutrinoE::NeutrinoEDefinition();
	G4AntiNeutrinoE::AntiNeutrinoEDefinition();
	G4MuonPlus::MuonPlusDefinition();
	G4MuonMinus::MuonMinusDefinition();
	G4NeutrinoMu::NeutrinoMuDefinition();
	G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();
}

void PhysicsList::ConstructMesons()
{
	//  mesons
	G4PionPlus::PionPlusDefinition();
	G4PionMinus::PionMinusDefinition();
	G4PionZero::PionZeroDefinition();
	G4KaonZero::KaonZeroDefinition();
	G4KaonPlus::KaonPlusDefinition();
	G4KaonMinus::KaonMinusDefinition();
}

void PhysicsList::ConstructBaryons()
{
	//  barions
	G4Proton::ProtonDefinition();
	G4AntiProton::AntiProtonDefinition();
	G4Neutron::NeutronDefinition();
	G4AntiNeutron::AntiNeutronDefinition();
}

void PhysicsList::ConstructIons()
{
	// ions
	G4IonConstructor iConstructor;
	iConstructor.ConstructParticle(); 
	G4GenericIon::GenericIonDefinition();
}

void PhysicsList::ConstructNuclei()
{
	G4Deuteron::Deuteron();
	G4Triton::Triton();
	G4He3::He3();
	G4Alpha::Alpha();
}

// This is required to specify which physical processes you want to simulate
void PhysicsList::ConstructProcess()
{
	AddTransportation();
	ConstructGeneral();
	ConstructEM();
	ConstructOp();
	ConstructScintillation();
	SetCuts();
}

void PhysicsList::ConstructEM()
{
	//Activate atomic deexcitation 
	G4EmProcessOptions* theParameters = new G4EmProcessOptions();
	theParameters -> SetFluo(true);
	theParameters -> SetAuger(true);
	theParameters -> SetPIXE(true);

	theParticleIterator->reset();

	while( (*theParticleIterator)() ){

		G4ParticleDefinition* particle = theParticleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		G4String particleName = particle->GetParticleName();
		//   G4PhysicsListHelper* plhelper = G4PhysicsListHelper::GetPhysicsListHelper();    //Get pointer to G4PhysicsListHelper

		if (particleName == "gamma"){// gamma
			// Construct processes for gamma

			pmanager->AddDiscreteProcess(new G4GammaConversion());
			pmanager->AddDiscreteProcess(new G4ComptonScattering());
			pmanager->AddDiscreteProcess(new G4RayleighScattering());
			pmanager->AddDiscreteProcess(new G4PhotoElectricEffect());
			pmanager->AddDiscreteProcess(new G4GammaConversionToMuons());

			/*	
			 plhelper->RegisterProcess(new G4GammaConversion(), particle);
			 plhelper->RegisterProcess(new G4ComptonScattering(), particle);
			 plhelper->RegisterProcess(new G4RayleighScattering(), particle);
			 plhelper->RegisterProcess(new G4PhotoElectricEffect(), particle);   
			 plhelper->RegisterProcess(new G4GammaConversionToMuons(), particle);
			 */
		} 
		else if (particleName == "e-"){//electron
			// Construct processes for electron			
			pmanager->AddProcess(new G4eMultipleScattering(),-1,1,1); // only elastic scattering. 
			pmanager->AddProcess(new G4eIonisation(),-1,2,2);
			pmanager->AddProcess(new G4eBremsstrahlung(),-1,3,3);///OJO, en varios -1, -3, 3; en otros cmo está aquí, otros -1, .1 ,3
			pmanager->AddProcess(new G4CoulombScattering()); 
			/*
			 plhelper->RegisterProcess(new G4eMultipleScattering(), particle);
			 plhelper->RegisterProcess(new G4eIonisation(), particle);
			 plhelper->RegisterProcess(new G4eBremsstrahlung(), particle);
			 plhelper->RegisterProcess(new G4CoulombScattering(), particle);
			 */
			G4Cerenkov*   theCerenkovProcess = new G4Cerenkov("Cerenkov");
			if (theCerenkovProcess->IsApplicable(*particle)) {
				pmanager->AddProcess(theCerenkovProcess);
				pmanager->SetProcessOrdering(theCerenkovProcess,idxPostStep);   
			}
		} 
		else if (particleName == "e+") {//positron
			// Construct processes for positron

			pmanager->AddProcess(new G4eMultipleScattering(),-1,1,1);
			pmanager->AddProcess(new G4eIonisation(),-1,2,2);
			pmanager->AddProcess(new G4eBremsstrahlung(),-1,3,3);///OJO, en varios -1, -3, 3; en otros cmo está aquí, otros -1, -1 ,3
			pmanager->AddProcess(new G4eplusAnnihilation(),0,-1,4);//Positron annihilation into two gammas
			pmanager->AddProcess(new G4CoulombScattering()); 
			/*
			 plhelper->RegisterProcess(new G4eMultipleScattering(), particle);
			 plhelper->RegisterProcess(new G4eIonisation(), particle);
			 plhelper->RegisterProcess(new G4eBremsstrahlung(), particle);
			 plhelper->RegisterProcess(new G4eplusAnnihilation(), particle);
			 plhelper->RegisterProcess(new G4CoulombScattering(), particle);
			 */
			G4Cerenkov*   theCerenkovProcess = new G4Cerenkov("Cerenkov");
			if (theCerenkovProcess->IsApplicable(*particle)) {
				pmanager->AddProcess(theCerenkovProcess);
				pmanager->SetProcessOrdering(theCerenkovProcess,idxPostStep);
			}
		} 
		else if(particleName == "mu+" || particleName == "mu-") {//muon
			// Construct processes for muon

			pmanager->AddProcess(new G4MuMultipleScattering(),-1,1,1);
			pmanager->AddProcess(new G4MuIonisation(),-1,2,2);
			pmanager->AddProcess(new G4MuBremsstrahlung(),-1,3,3);///OJO, en otros -1, -1, 3
			pmanager->AddProcess(new G4MuPairProduction(),-1,4,4);///OJO, en otros -1, -1, 4
				/*
				 plhelper->RegisterProcess(new G4MuMultipleScattering(), particle);
				 plhelper->RegisterProcess(new G4MuIonisation(), particle);
				 plhelper->RegisterProcess(new G4MuBremsstrahlung(), particle);
				 plhelper->RegisterProcess(new G4MuPairProduction(), particle);
				 */
			G4Cerenkov*   theCerenkovProcess = new G4Cerenkov("Cerenkov"); 
			//DUDA: diferencia entre esto y hacer directamentepmanager->AddProcess(new G4Cerenkov())????
			if (theCerenkovProcess->IsApplicable(*particle)) {
				pmanager->AddProcess(theCerenkovProcess);
				pmanager->SetProcessOrdering(theCerenkovProcess,idxPostStep);
			}
		}
		else if ((!particle->IsShortLived()) &&
		         (particle->GetPDGCharge() != 0.0) && 
		         (particle->GetParticleName() != "chargedgeantino")) {
					 //all others charged particles except geantino

					 pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
					 pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
					 /*
					  plhelper->RegisterProcess(new G4hMultipleScattering, particle);
					  plhelper->RegisterProcess(new G4hIonisation, particle);
					  */
				 }

	}//end while


}//end PhysicsList::ConstructEM()


void PhysicsList::ConstructOp()
{   

  theParticleIterator->reset();
  while( (*(theParticleIterator))() )
    {
      G4ParticleDefinition* particle = theParticleIterator->value();
      G4ProcessManager* pmanager = particle->GetProcessManager();
      G4String particleName = particle->GetParticleName();
		/*
		if (particleName == "opticalphoton") { 	         	
			 pmanager->AddProcess(new G4OpBoundaryProcess());
			 pmanager->AddProcess(new G4OpAbsorption());
			 pmanager->AddProcess(new G4OpRayleigh());
		}	
		*/
		if (particleName == "opticalphoton") {		
	    //http://geant4.cern.ch/G4UsersDocuments/UsersGuides/ForApplicationDeveloper/html/TrackingAndPhysics/physicsProcess.html  (ver 5.2.5)				
		  pmanager->AddDiscreteProcess(new G4OpAbsorption());//In output it is indicated as OpAbsorption. In DetectorConstruction.cc material abosoption length needs to be indicated.
		  pmanager->AddDiscreteProcess(new G4OpRayleigh());//no Rayleigh scattering attenutation length specified by the user, the program automatically calls the RayleighAttenuationLengthGenerator
		  // 24/10/2016 pmanager->AddProcess(new G4OpBoundaryProcess());
		  pmanager->AddDiscreteProcess(new G4OpBoundaryProcess());

		}
		
		/*
		 G4PhysicsListHelper* helper = G4PhysicsListHelper::GetPhysicsListHelper();
		 helper->RegisterProcess(new G4OpBoundaryProcess(), OpPhoton);
		 helper->RegisterProcess(new G4OpAbsorption(), OpPhoton);
		 helper->RegisterProcess(new G4OpRayleigh(), OpPhoton);
		 */

	}//end while
}


void PhysicsList::ConstructGeneral()//Aquí inlcuyo el decay y protón, alpha, ion.
{

	G4Decay* theDecayProcess = new G4Decay();
	theParticleIterator->reset();
	while( (*theParticleIterator)() ){
		G4ParticleDefinition* particle = theParticleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		if (theDecayProcess->IsApplicable(*particle)) {
			pmanager->AddDiscreteProcess(theDecayProcess); 
			// set ordering for PostStepDoIt and AtRestDoIt
			pmanager ->SetProcessOrdering(theDecayProcess, idxPostStep);
			pmanager ->SetProcessOrdering(theDecayProcess, idxAtRest);
		}
	}

	theParticleIterator->reset();

	while( (*theParticleIterator)() ){

		G4ParticleDefinition* particle = theParticleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		G4String particleName = particle->GetParticleName();
		//G4PhysicsListHelper* plhelper = G4PhysicsListHelper::GetPhysicsListHelper();    //Get pointer to G4PhysicsListHelper

		if( particleName == "proton" ||
		   particleName == "pi-" ||
		   particleName == "pi+"    ) {
			   //proton  

			   pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
			   pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
			   pmanager->AddProcess(new G4hBremsstrahlung,     -1, 3, 3);
			   pmanager->AddProcess(new G4hPairProduction,     -1, 4, 4);              
			   /*
				plhelper->RegisterProcess(new G4hMultipleScattering, particle);
				plhelper->RegisterProcess(new G4hIonisation, particle);
				plhelper->RegisterProcess(new G4hBremsstrahlung, particle);
				plhelper->RegisterProcess(new G4hPairProduction, particle);
				*/
		   } 

		else if( particleName == "alpha" || 
		        particleName == "He3" ) {
					//alpha and He3

					pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
					pmanager->AddProcess(new G4ionIonisation,       -1, 2, 2);
					pmanager->AddProcess(new G4NuclearStopping,     -1, 3,-1);             

					/*	plhelper->RegisterProcess(new G4hMultipleScattering, particle);
					plhelper->RegisterProcess(new G4ionIonisation, particle);
					plhelper->RegisterProcess(new G4NuclearStopping, particle);
					*/
				} 

		else if( particleName == "GenericIon" ) {
			//Ions 

			//de esta forma es necesario descomentar la definición de pmanager arriba.
			G4ionIonisation* ionIoni = new G4ionIonisation();
			ionIoni->SetEmModel(new G4IonParametrisedLossModel());      
			pmanager->AddProcess(ionIoni,       -1, 2, 2);      
			pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
			pmanager->AddProcess(new G4NuclearStopping,     -1, 3,-1);

			/*
			 plhelper->RegisterProcess(new G4hMultipleScattering, particle);
			 plhelper->RegisterProcess(new G4ionIonisation, particle);
			 plhelper->RegisterProcess(new G4NuclearStopping, particle);	
			 */                
		} 

		//Añado aquí también para que al resto le tenga en cuenta estos procesos
		else if ((!particle->IsShortLived()) &&
		         (particle->GetPDGCharge() != 0.0) && 
		         (particle->GetParticleName() != "chargedgeantino")) {
					 //all others charged particles except geantino
					 pmanager->AddProcess(new G4hMultipleScattering, -1, 1, 1);
					 pmanager->AddProcess(new G4hIonisation,         -1, 2, 2);
					 /*
					  plhelper->RegisterProcess(new G4hMultipleScattering, particle);
					  plhelper->RegisterProcess(new G4hIonisation, particle);
					  */
				 }
		else if (particle->GetParticleName() == "chargedgeantino") { // OK, no sale ninguno

			G4cout << G4endl << "hay chargedgeantino ? " << particle->GetParticleName() << G4endl << G4endl;

			}

	}//end while

}//end ConstructGeneral



void PhysicsList::ConstructScintillation()
{
 //G4bool theScintProcessDefNeverUsed = true;

	G4Scintillation* Scint = new G4Scintillation("Scintillation");
	Scint->SetTrackSecondariesFirst(true);

  //allows for different scintillation yields depending on the particle type – in such case, separate scintillation processes must be attached 		to the various particles. http://hypernews.slac.stanford.edu/HyperNews/geant4/get/opticalphotons/389/1.html
	//ExcitationRatio = 1.0; is the default, where the ratio is the intensity of fast component to the total scintillation intensity . http://hypernews.slac.stanford.edu/HyperNews/geant4/get/opticalphotons/219.html
	Scint->SetScintillationYieldFactor(1.);
 	Scint->SetScintillationExcitationRatio(1.);
        
	//Use Birks Correction in the Scintillation process. 
	G4EmSaturation* emSaturation = G4LossTableManager::Instance()->EmSaturation();
	Scint->AddSaturation(emSaturation);


 	Scint->SetVerboseLevel(1);
	theParticleIterator->reset();
	while ((*theParticleIterator)()) {  
		G4ParticleDefinition* particle = theParticleIterator->value();
		G4ProcessManager* pmanager = particle->GetProcessManager();
		if (Scint->IsApplicable(*particle)){
		//CONTINUO SIN PROCESSORDERING NO CENTELLEA. CON LAS DOS LÍNEAS, CONTINUO Y DISCRETO ES IGUAL. EN DISCRETO SIN LAS DOS LINEAS SALEN MENOS PARTICULAS PERO LA DIFERENCIA ES ÍNFIMA
		pmanager->AddDiscreteProcess(Scint);
		//pmanager->AddProcess(Scint);
	        pmanager->SetProcessOrderingToLast(Scint, idxAtRest); 
	        pmanager->SetProcessOrderingToLast(Scint, idxPostStep);
		}
	}
	
	/*
	 G4Scintillation* Scint = new G4Scintillation("Scintillation");
	 Scint->SetTrackSecondariesFirst(true);   
	 G4PhysicsListHelper* helper = G4PhysicsListHelper::GetPhysicsListHelper();
	 theParticleIterator->reset();
	 while ((*theParticleIterator)()) {
		 G4ParticleDefinition* particle = theParticleIterator->value();
		 if (Scint->IsApplicable(*particle))helper->RegisterProcess(Scint, particle);
	}
	*/

}


void PhysicsList::SetCuts() //to avoid infrared divergence, threshold below no secondary are generated
{

	//G4VUserPhysicsList::SetCutsWithDefault method sets the default cut value for all particle types
	SetCutsWithDefault();

	// set cut values for gamma at first and for e- second and next for e+, because some processes for e+/e- need cut values for gamma 
	//SetCutValue(defaultCutValue, "gamma");
	//SetCutValue(defaultCutValue, "e-");
	//SetCutValue(defaultCutValue, "e+");

	if(verboseLevel>0) DumpCutValuesTable(); //Request to print out information of cut values
}

