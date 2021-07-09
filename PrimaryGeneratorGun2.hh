#ifndef PrimaryGeneratorGun2_h
#define PrimaryGeneratorGun2_h 1
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ParticleGun.hh"
#include "DetectorConstruction.hh"
#include "globals.hh"

#include "G4VUserPrimaryGeneratorAction.hh"
class G4ParticleGun;
class G4Event;
class G4VPrimaryGenerator;
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PrimaryGeneratorGun2 : public G4VUserPrimaryGeneratorAction
{
private:  
  G4VPrimaryGenerator* InitializeGPS();
  G4VPrimaryGenerator* gun;
 std::ofstream* outfile;
  public:
    PrimaryGeneratorGun2();
   ~PrimaryGeneratorGun2();

  public:
    virtual void GeneratePrimaries(G4Event*);

    G4double GetRadius() {return fRadius;};
    G4double GetAlphaMin() {return fAlphaMin;};
    G4double GetAlphaMax() {return fAlphaMax;};
    G4double GetGamma() {return fgamma;};
    G4double GetEmin() {return fEmin;};
    G4double GetEmax() {return fEmax;};
    G4ParticleGun* GetParticleGun() {return fParticleGun;};

  private:
    G4double fRadius;             //vertex surface
    G4double fAlphaMin;
    G4double fAlphaMax;        //opening angle
    G4double fgamma;		  //variables for energy
    G4double fEmin;               
    G4double fEmax;
    //G4ParticleDefinition* fParticle;

   G4ParticleGun*             fParticleGun;
    DetectorConstruction*      fDetector;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
