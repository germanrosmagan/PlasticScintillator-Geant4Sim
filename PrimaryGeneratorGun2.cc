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
/// \file eventgenerator/Gun/src/PrimaryGeneratorGun2.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
//
// $Id: PrimaryGeneratorGun2.cc 68734 2013-04-05 09:47:02Z gcosmo $
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "PrimaryGeneratorGun2.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
//#include "G4Isotope.hh"
//#include "G4Element.hh"
//#incluede "G4NucleiProperties.hh"
//#include "G4Deuteron.hh" 
//#include "G4IonConstructor.hh"
#include "G4IonTable.hh"

#include <fstream>
using namespace std;

#define G4UniformRandom() CLHEP::HepRandom::getTheEngine()->flat()

PrimaryGeneratorGun2::PrimaryGeneratorGun2()
 : G4VUserPrimaryGeneratorAction(), fParticleGun(0)
{
  // default particle kinematic
  //
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
 
  //Set variables for Particle energy. Variables defined in PrimaryGeneratorGun2.hh
  //Esto es para lanzar energías según E-gamma.
  //fgamma=2.85;//AMS02, PRL114, 2015
  //fEmin=1*MeV;//Less than 1GeV no particles due to geomagnetic cutoff, shield by Earths magnetic field (Gustavo: no necesariamente)
  //fEmax=100*MeV;

  // vertex -> it is used for isotropic or uniform. Ratio of the sphere 
  fRadius = 20*cm;
  
  //opening angle (uniform vs isotropic; it is used below)
  //Lanzar uniform
  //fAlphaMin =  180.*deg;
  //fAlphaMax = 360.*deg;
  //Lanzar isotropic	
  fAlphaMin =  0.*deg;
  fAlphaMax = 90.*deg;

  
}

PrimaryGeneratorGun2::~PrimaryGeneratorGun2()
{
  delete fParticleGun;
}

void PrimaryGeneratorGun2::GeneratePrimaries(G4Event* anEvent)
{

  //Set primary particle--------------------------------------------------------------------
  //
  G4ParticleDefinition* particle;

  G4bool ShootPimary=0;
  //0->partícula. Especificar "e-" que es su nombre
  //1->nucleo. Especificar NUMEROATOMICO y NUMEROMASICO
  

   if(ShootPimary==0){
    particle = G4ParticleTable::GetParticleTable()->FindParticle("e-");
   } 

   if(ShootPimary==1){
      G4int Z = 26, A = 55;
      G4double excitEnergy= 0.*eV;      
      particle = G4IonTable::GetIonTable()->GetIon(Z, A, excitEnergy);//Z, A, excit. energy      
   }
  
  
  fParticleGun->SetParticleDefinition(particle);
  G4cout << G4endl << "Particle " << particle->GetParticleName() << G4endl << G4endl;


  //Set Particle energy-----------------------------------------------------------------

  //E^(-gamma) entre Emin y Emax
  //G4double Eminim = std::pow(fEmin,1.-fgamma);
  //G4double Emaxim = std::pow(fEmax,1.-fgamma);
  //G4double energy = std::pow( Eminim + G4UniformRandom()*(Emaxim-Eminim) , 1./(1.-fgamma) );//in MeV

  G4double energy = 1*MeV;

  fParticleGun->SetParticleEnergy(energy);

  //Launch----------------------------------------------- ------------------
  // 0-> is fired uniformly or isotropically (now isotropically) from a sphere surrounding the detector
  // 1 -> it is thrown diagonally to the plastic
  // 2-4 -> is thrown perpendicular to the faces of the detector (3 is the one in front, 2 and 4 are equivalent)
  // 5-> is like -3, towards the detector going through wheels
  G4int ShootOpt=0;	

  if(ShootOpt==0){//ISOTROPIC OR UNIFORM, SELECT BELOW
	  	
	  //http://geant4.web.cern.ch/geant4/UserDocumentation/Doxygen/examples_doc/html/ExampleparticleGun.html
	  //http://www.apc.univ-paris7.fr/~franco/g4doxy4.10/html/_primary_generator_action4_8cc_source.html
	  //http://geant4.web.cern.ch/geant4/geant4_public/source/geant4/examples/extended/eventgenerator/particleGun/README

	  //vertex position uniform in sphere	
	  G4double cosTheta = 2*G4UniformRand() - 1;  //cosTheta uniform in [0, pi]
	  G4double sinTheta = std::sqrt(1. - cosTheta*cosTheta);
	  G4double phi      = twopi*G4UniformRand();  //phi uniform in [0, 2*pi]
	  G4ThreeVector ur(sinTheta*std::cos(phi),sinTheta*std::sin(phi),cosTheta);

	  fParticleGun->SetParticlePosition(fRadius*ur);

	  //Option 1: particle direction uniform around ur.  
          //dN=A sin(alpha) dalpha dphi donde A es una constante de normalizacion  
	  //cosAlpha uniform in [cos(alphaMin), cos(alphaMax)]
	  //G4double CosAlphaMin = std::cos(fAlphaMin);
  	  //G4double CosAlphaMax = std::cos(fAlphaMax);  
	  //G4double cosAlpha = CosAlphaMin-G4UniformRand()*(CosAlphaMin-CosAlphaMax);
	  //CosAlphaMin+G4UniformRand()*(CosAlphaMax-CosAlphaMin) pero da lo mismo
	  //G4double sinAlpha = std::sqrt(1. - cosAlpha*cosAlpha);

	  //Option 2: particle direction isotropic around ur. 
          //dN=A cos(alpha) sin(alpha) dalpha dphi
          //sin2Alpha uniform in [sin^2(alphaMin), sin^2(alphaMax)]
	  G4double Sin2AlphaMin = std::sin(fAlphaMin)*std::sin(fAlphaMin);
	  G4double Sin2AlphaMax = std::sin(fAlphaMax)*std::sin(fAlphaMax);
	  G4double sin2Alpha = Sin2AlphaMin + G4UniformRand()*(Sin2AlphaMax-Sin2AlphaMin);
 	  G4double sinAlpha = std::sqrt(sin2Alpha);
	  G4double cosAlpha = std::sqrt(1. - sin2Alpha);


          //For both options psi uniform in (0,2*pi)
	  G4double psi      = twopi*G4UniformRand();
	  G4ThreeVector dir(sinAlpha*std::cos(psi),sinAlpha*std::sin(psi),cosAlpha); //rotate dir (rotateUz transforms uz to ur) It is composition of two simple rotations: theta around oy, then phi around oz (non commutative).
	  dir.rotateUz(ur);           

	  fParticleGun->SetParticleMomentumDirection(-1.0*dir);;//-1.0 para lanzar hacia dentro en lugar de hacia afuera
 
  }

 if(ShootOpt==1){

    G4double nn1=2.5*3/2.;//2.5cm is detector length
    G4double nn2=2.5*3/2.*fRadius;
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(nn1,-1.0*nn1,-1.0*nn1));//(CentX/2,-CentY/2,-CentZ/2)//Cent mide 2.5*3
    fParticleGun->SetParticlePosition(G4ThreeVector(-1.0*nn2,nn2,nn2));//Radius*(-CentX/2,CentY/2,CentZ/2)

  }


  if(ShootOpt==2){

    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,-1.,0.));
	fParticleGun->SetParticlePosition(G4ThreeVector(0.,1.0*fRadius,0.));

  }

  if(ShootOpt==3){

    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
    fParticleGun->SetParticlePosition(G4ThreeVector(-1.0*fRadius,0.,0.));

  }
  if(ShootOpt==4){

    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,-1.));
    fParticleGun->SetParticlePosition(G4ThreeVector(0.,0.,1.0*fRadius));

  }
  if(ShootOpt==5){

    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(-1.,0.,0.));
    fParticleGun->SetParticlePosition(G4ThreeVector(1.0*fRadius,0.,0.));

  }


  //create vertex
  fParticleGun->GeneratePrimaryVertex(anEvent);
 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
