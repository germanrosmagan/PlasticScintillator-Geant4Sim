
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
/// \file electromagnetic/TestEm3/src/SteppingAction.cc
/// \brief Implementation of the SteppingAction class
//

#include "SteppingAction.hh"
#include "G4SystemOfUnits.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "G4UnitsTable.hh"
#include "G4Step.hh"
#include "G4Positron.hh"
#include "G4RunManager.hh"
#include "G4PhysicalConstants.hh"
#include "G4DynamicParticle.hh"
#include "G4SystemOfUnits.hh"

#include <fstream>
using namespace std;
extern ofstream myfile2;//detector
extern ofstream myfile3;//track
extern ofstream myfile4;//processNull
extern ofstream myfile5;//detalles de cada step.
extern ofstream myfile6;//OpAbsoprtion
extern ofstream myfile7;//Primary particle information
extern ofstream myfile8;//TrueEvent

SteppingAction::SteppingAction(DetectorConstruction* det, RunAction* run, EventAction* evt)
:G4UserSteppingAction(),fDetector(det),fRunAct(run),fEventAct(evt) 
{
}

SteppingAction::~SteppingAction()
{
}
//********************************************************************************************************//
//All secondaries generated along a track are accumulated in G4SteppingManager and there is no direct way of 
//accessing to the secondaries produced by a particular step. From your UserSteppingAction, you can get all 
//such information by accessing appropriately to G4SteppingManager.
//********************************************************************************************************//

void SteppingAction::UserSteppingAction(const G4Step* aStep)
{

 G4bool SacarFile5 = false;//Genera archivos enormes y tarda. Solo para verificar cosas. No lo saco aunque debajo saco este file simplificado

 if(SacarFile5) 
 {

   myfile5 << endl << endl << endl << "Nuevo Step" << endl << endl;


   //http://geant4.slac.stanford.edu/Tips/event/6.html

   myfile5 << "Step is limited by " << aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName()  << G4endl;
   myfile5 << "StepStatus " << aStep->GetPostStepPoint()->GetStepStatus() << G4endl;
   myfile5 << "Processes involved to the step" << G4endl;
  G4StepStatus stepStatus = fpSteppingManager->GetfStepStatus();

  if(stepStatus==fAtRestDoItProc)
  {
    G4ProcessVector* procAtRest = fpSteppingManager->GetfAtRestDoItVector();
    G4SelectedAtRestDoItVector* selProcAtRest
     = fpSteppingManager->GetfSelectedAtRestDoItVector();
    size_t MAXofAtRestLoops = fpSteppingManager->GetMAXofAtRestLoops();
    for(size_t i1=0;i1<MAXofAtRestLoops;i1++)
    {
      if((*selProcAtRest)[MAXofAtRestLoops-i1-1]==2)
      {  myfile5 << "  At rest : " << (*procAtRest)[i1]->GetProcessName() << " (forced)" << G4endl; }
      else if((*selProcAtRest)[MAXofAtRestLoops-i1-1]==1)
      {  myfile5 << "  At rest : " << (*procAtRest)[i1]->GetProcessName() << G4endl; }
    }
  }

  if(stepStatus!=fExclusivelyForcedProc && stepStatus!=fAtRestDoItProc)
  {
    G4ProcessVector* procAlong = fpSteppingManager->GetfAlongStepDoItVector();
    size_t MAXofAlongStepLoops = fpSteppingManager->GetMAXofAlongStepLoops();
    for(size_t i2=0;i2<MAXofAlongStepLoops;i2++)
    {
      if((*procAlong)[i2]!=0)
      myfile5 << "  Along step : " << (*procAlong)[i2]->GetProcessName() << G4endl;
    }
  }

  if(stepStatus!=fAtRestDoItProc)
  {
    G4ProcessVector* procPost = fpSteppingManager->GetfPostStepDoItVector();
    G4SelectedPostStepDoItVector* selProcPost
     = fpSteppingManager->GetfSelectedPostStepDoItVector();
    size_t MAXofPostStepLoops = fpSteppingManager->GetMAXofPostStepLoops();
    for(size_t i3=0;i3<MAXofPostStepLoops;i3++)
    {
      if((*selProcPost)[MAXofPostStepLoops-i3-1]==2)
      {  myfile5 << "  Post step : " << (*procPost)[i3]->GetProcessName() << " (forced)" << G4endl; }
      else if((*selProcPost)[MAXofPostStepLoops-i3-1]==1)
      {  myfile5 << "  Post step : " << (*procPost)[i3]->GetProcessName() << G4endl; }
    }
  }

  G4int nSecAtRest2 = fpSteppingManager->GetfN2ndariesAtRestDoIt();
  G4int nSecAlong2  = fpSteppingManager->GetfN2ndariesAlongStepDoIt();
  G4int nSecPost2   = fpSteppingManager->GetfN2ndariesPostStepDoIt();
  G4int nSecTotal2  = nSecAtRest2+nSecAlong2+nSecPost2;
  //G4TrackVector* secVec = fpSteppingManager->GetfSecondary();

  if(nSecTotal2>0)
  {
     myfile5 << "  :----- List of 2ndaries - " << std::setw(3) << nSecTotal2
           << " (Rest=" << std::setw(2) << nSecAtRest2
           << ",Along=" << std::setw(2) << nSecAlong2
           << ",Post="  << std::setw(2) << nSecPost2 << ")" << G4endl;

	/* 
	    for(size_t lp1=(*secVec).size()-nSecTotal2; lp1<(*secVec).size(); lp1++)
	    {
	      myfile5 << "    : "
		     << G4BestUnit((*secVec)[lp1]->GetPosition(), "Length") << " "
		     << std::setw( 9) << G4BestUnit((*secVec)[lp1]->GetKineticEnergy() , "Energy") << " "
		     << std::setw(18) << (*secVec)[lp1]->GetDefinition()->GetParticleName()
		     << " generated by " << (*secVec)[lp1]->GetCreatorProcess()->GetProcessName() << G4endl;
	    }
	*/ 
  }

 } //end if file5

	//******************************************************************//
	
	const G4StepPoint* prePoint = aStep->GetPreStepPoint();
	const G4StepPoint* endPoint = aStep->GetPostStepPoint();
	//To get their positions in the global coordinate system:
	G4ThreeVector preStepPosition = prePoint->GetPosition()/cm;
	G4ThreeVector postStepPosition = endPoint->GetPosition()/cm;
	//To get energies
    	G4double energypre = aStep->GetPreStepPoint()->GetTotalEnergy() ;//MeV
	G4double energypos = aStep->GetPostStepPoint()->GetTotalEnergy();//MeV
	//To get momentums
	G4ThreeVector momentumpre = aStep->GetPreStepPoint()->GetMomentum();
    	G4ThreeVector momentumpos = aStep->GetPostStepPoint()->GetMomentum();
    	//To get materials
	G4String Material1 = aStep->GetPreStepPoint()->GetPhysicalVolume()->GetName();
    


	if(Material1 != "Galactic") myfile8 << "YES" << endl;
	//const G4ParticleDefinition* particle1 = aStep->GetTrack()->GetDefinition();  //= dynamic particle
    	const G4DynamicParticle* dynparticle = aStep->GetTrack()->GetDynamicParticle(); // G4DynamicParticle tiene momento, masa, carga....
	const G4ParticleDefinition* particle = dynparticle->GetDefinition();  
        pname=particle->GetParticleName();
	const G4String aProcess = endPoint->GetProcessDefinedStep()->GetProcessName();		
	if (aStep->GetPostStepPoint()->GetPhysicalVolume() != NULL) 
		if (aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName() != "Galactic") 
			myfile8 << "YES" << endl;

	
	//To get the process which has limited the current step: 
	//A step is a delta information to a track. track is not a collection of steps. Instead, track is updated by steps. 
     

	//******************************************************************//
	//G4TouchableHandle theTouchable = preStepPoint->GetTouchableHandle();
	//G4ThreeVector worldPosition = preStepPoint->GetPosition();
	//G4ThreeVector localPosition = theTouchable->GetHistory()-> GetTopTransform().TransformPoint(worldPosition);
	//where worldPosition here stays for the position related to the world volume, while localPosition refers to the coordinates
	// local to the volume where the particle is currently placed. 
	//******************************************************************//
	

	G4int   nbStep = aStep->GetTrack()->GetCurrentStepNumber();
    	G4double Trleng = aStep->GetTrack()->GetTrackLength()/cm;//se va actualizando, va sumando lo que se recorre en cada Step

	//Numero total de secundarios: "  No. of secodaries = "	
	G4int nSecAtRest = fpSteppingManager->GetfN2ndariesAtRestDoIt();
	G4int nSecAlong  = fpSteppingManager->GetfN2ndariesAlongStepDoIt();
	G4int nSecPost   = fpSteppingManager->GetfN2ndariesPostStepDoIt();
	G4int nSecTotal  = nSecAtRest+nSecAlong+nSecPost; //numero total de secundarios en el Step
	//const std::vector<const G4Track*>* secondaries = aStep->GetSecondaryInCurrentStep();

    	G4double edep = aStep->GetTotalEnergyDeposit()/MeV;
    	//Total energy deposited during the step - this is the sum of: the energy deposited by the energy loss process, and the energy lost by secondaries which have NOT been generated because each of their energies was below the cut threshold
	// collect energy deposit taking into account track weight
    	G4double edep2 = aStep->GetTotalEnergyDeposit()*aStep->GetTrack()->GetWeight()/MeV; //sale igual que edep (weight es 1). Quizá los track sólo tienen pesos en los casos como las cascadas de partículas
    	if(edep!=edep2) cout<< endl << endl << "OJO: edep!=edep2 " << endl << endl;
	
	G4double Kineticenergy=dynparticle->GetKineticEnergy();
	G4double Totalenergy=dynparticle->GetTotalEnergy();    // kinetic energy + rest energy -> = energypos
	G4double Time=aStep->GetDeltaTime(); //delta entre pre y post
	G4double TimeTr=aStep->GetTrack()->GetGlobalTime(); //it is updated in every step
    	G4TrackVector *fSecondary=fpSteppingManager->GetfSecondary();// Numner of secondarys generated along the Track. It is updated in Steps
	
    /*angle between PRESTEP and POSSTEP
    G4double momentumproduct =  momentumpre(0)*momentumpos(0)+momentumpre(1)*momentumpos(1)+momentumpre(2)*momentumpos(2);
    G4double momentumpremodule = std::sqrt( momentumpre(0)*momentumpre(0) +  momentumpre(1)*momentumpre(1) + momentumpre(2)*momentumpre(2) );
    G4double momentumposmodule = std::sqrt( momentumpos(0)*momentumpos(0) +  momentumpos(1)*momentumpos(1) + momentumpos(2)*momentumpos(2) );
    G4double anglebetmom = std::acos(momentumproduct/(momentumpremodule*momentumposmodule))*deg;
    */

    //http://geant4.slac.stanford.edu/Tips/event/6.html
	G4String Material2="NULL";
 	if((aStep->GetTrack()->GetTrackID()>0)&&(aStep->GetPostStepPoint()->GetPhysicalVolume()!= NULL))	
	Material2 = aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();	

	if(SacarFile5) 
        myfile5 << "  " << aStep->GetTrack()->GetTrackID() << "  " << Material1 << "  " << Material2 << "  " << 
 	postStepPosition[0] << "  " <<  postStepPosition[1]<< "  " <<  postStepPosition[2] << "  " << pname << "  " << aProcess << " " << (*fSecondary).size() << endl;


   //if(pname=="opticalphoton") myfile << aStep->GetTrack()->GetTrackID();

   //if(aProcess!="Transportation" && aProcess!="OpAbsorption")
  
   if(aProcess == "OpAbsorption")
	 myfile6 << "TRACK ID = " << aStep->GetTrack()->GetTrackID() << " que viene del TRACK = " << aStep->GetTrack()->GetParentID() 
    		 << "	  "  << aProcess  << "   " << pname << "  AcomulatedTrackLength(cm) = " <<  Trleng << " time(ns) = " <<TimeTr/ns
		 << " Before step Tot.Energy(eV) = " << energypre/eV << "    After step Tot.Energy(eV) " << energypos/eV 
		 << "	PreStepPosition(cm) = " << preStepPosition << " in " << Material1
		 << "	PostStepPosition(cm) = " << postStepPosition 
		 << endl;

		

    if(aStep->GetTrack()->GetTrackID()==1)  //primary info   
     myfile7 << aProcess  << "   " << pname << " " << Material1 << " " << Material2 
	 	<< " 	AcomulatedTrackLength(cm)= " <<  Trleng << " time(ns)= " <<TimeTr/ns
       		<< "	BeforeStepTotEnergy(eV)= " << energypre/eV << "    AfterStepTotEnergy(eV)= " << energypos/eV 
  	     	<< "    EnergyDep(eV)= "  << edep/eV 
     		<< "	NumSecundariosEnStep= "  <<  nSecTotal 
	        << "	PreStepPosition(cm)= " << preStepPosition[0] << "  " <<  preStepPosition[1]<< "  " <<  preStepPosition[2]
                << "	PostStepPosition(cm)= " <<  postStepPosition[0] << "  " <<  postStepPosition[1]<< "  " <<  postStepPosition[2]
		<< endl << endl;
  
	
	if(aStep->GetTrack()->GetTrackID()==1 && Material2=="Plastic") aStep->GetTrack()->SetTrackStatus(fStopAndKill);	


	//secondaries
	
    G4double releasedenergy = 0;
    G4double EnergyScint=0, EnergyCer=0;
    int Count=0;

    for(size_t lp=(*fSecondary).size()-nSecTotal; lp<(*fSecondary).size(); lp++)
	//(*fSecondary).size() -> es el número de secundarios acumulados en el track
	//nsecTotal -> es el número de secundarios del Step
	{
			
	    if( (*fSecondary)[lp]->GetCreatorProcess()->GetProcessName() == "Cerenkov")  EnergyCer += (*fSecondary)[lp]->GetTotalEnergy()/eV;
	    
	    if( (*fSecondary)[lp]->GetCreatorProcess()->GetProcessName() == "Scintillation") Count++;
	    if( (*fSecondary)[lp]->GetCreatorProcess()->GetProcessName() == "Scintillation") EnergyScint += (*fSecondary)[lp]->GetTotalEnergy()/eV;

	    releasedenergy+=(*fSecondary)[lp]->GetTotalEnergy()/MeV;//cinética+reposo

	    if(SacarFile5) 
            {
		 myfile5 << "	Parent ID " <<  aStep->GetTrack()->GetTrackID() << " " 
		 << (*fSecondary)[lp]->GetCreatorProcess()->GetProcessName() << "	Part. name "  
		 << (*fSecondary)[lp]->GetDefinition()->GetParticleName() << " Pos(cm) = " 
		 << (*fSecondary)[lp]->GetPosition()/cm << "   Tot.Energy(eV) = "  << (*fSecondary)[lp]->GetTotalEnergy()/eV << "    Kin.Energy(eV) "  
		 << (*fSecondary)[lp]->GetKineticEnergy()/eV << "     Time(ns) "  << (*fSecondary)[lp]->GetGlobalTime()/ns
		 << "  energypre  " << energypre/eV << "   energypost  "  << energypos/eV << "   Count " << Count << G4endl;
	     }
	


		myfile3 <<  aStep->GetTrack()->GetTrackID() << " " << Count << "  " << (*fSecondary)[lp] ->GetDefinition()->GetParticleName() << "  " 				<< (*fSecondary)[lp]->GetCreatorProcess()->GetProcessName() << "  " << (*fSecondary)[lp]->GetKineticEnergy() << "  " 
			<< (*fSecondary)[lp]->GetTotalEnergy() << "  " 
			<< (*fSecondary)[lp]->GetPosition().x()/cm  << " " << (*fSecondary)[lp]->GetPosition().y()/cm  << " " 
			<< (*fSecondary)[lp]->GetPosition().z()/cm << endl;


   	}//end for secs del step


    if( energypos+edep+releasedenergy != energypre &&  nSecTotal!=0) 
       //cout << G4endl<< G4endl<< G4endl<< "FALLA LA SUMA DE ENERGÍAS!!!!!!!" << " " << aProcess << " " << energypre << " " << energypos << " " << edep << " " << releasedenergy << endl<< endl<< endl<< endl;


   
   if(aStep->GetTrack()->GetCreatorProcess()== NULL)			
	myfile4 << "ID " << aStep->GetTrack()->GetTrackID()   << "  ParentID " << aStep->GetTrack()->GetParentID() 
	<< "    " << pname << "  AcomulatedTrackLength(cm) " <<  Trleng << "  time(ns) " <<TimeTr/ns
	<< "    Position " << postStepPosition[0]<< " " << postStepPosition[1] << " " << postStepPosition[2]
	<< "    Energy(eV) "  << Totalenergy/eV << "    KinEnergy(eV) "  << Kineticenergy/eV << G4endl;

	
	
    //aStep->GetPostStepPoint()->GetPhysicalVolume()!= NULL

    //aStep->GetTrack()->GetCreatorProcess()!= NULL -> para evita break con el primario

    if( (aStep->GetPostStepPoint()->GetPhysicalVolume()!= NULL) && (aStep->GetTrack()->GetCreatorProcess()!= NULL) ) 
    {
	       
		G4String Material2 = aStep->GetPostStepPoint()->GetPhysicalVolume()->GetName();	
		
	
		
		if( (pname=="opticalphoton" || pname=="gamma") &&
		    ( (Material1=="Plastic" && Material2=="Sensor") || (Material1=="Sensor" && Material2=="Plastic") ) )
		{	
			//aProcess = endPoint->GetProcessDefinedStep()->GetProcessName() es siempre Transportation hasta el detector
			//aStep->GetTrack()->GetCreatorProcess()->GetProcessName() me da el proceso que creo esta partícula,
			       
			myfile2 << "ID " << aStep->GetTrack()->GetTrackID()   << "  ParentID " << aStep->GetTrack()->GetParentID() 
			<< "    " << pname << " "<< aStep->GetTrack()->GetCreatorProcess()->GetProcessName()  << "  AcomulatedTrackLength(cm) " 
			<<  Trleng << "  time(ns) " <<TimeTr/ns
			<< "    Position(cm) " << postStepPosition[0] << " " << postStepPosition[1] << " " << postStepPosition[2]
			<< "    Energy(eV) "  << Totalenergy/eV << "    KinEnergy(eV) "  << Kineticenergy/eV << G4endl;

				aStep->GetTrack()->SetTrackStatus(fStopAndKill);	
		
		}//end if 
	       
     }//end if evitar NULL

}//END 
