#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

/*
  In GEANT4 the object which you are going to pass particles through
  is called "the Detector" in GEANT4.  You build a detector
  by specifying it's shape (eg Cube, Pyrimid, Sphere, Cyclinder, etc.)
  it's material it is made of, in GEANT4 this is refered to
  as a "Logical Volume."  Finally you place the objects that
  make up the detector in a base volume known in GEANT4 as 
  "the World."

  The DetectorConstruction.hh and .cc
  defines the mandatory user class "DetectorConstruction."
  Setup of the detector geometry, definition of dectector properties,
  and location are all handled here.
*/

//Standard to all simulations
#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

//to be able to specify 3-vectors in GEANT4
#include "G4ThreeVector.hh"
#include "TargetSD.hh"
//Class prototypes, so you can use them without having to define them yet
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4Material;
class G4VSensitiveDetector;
class G4LogicalBorderSurface;

// The inherited G4VUserDetectorConstruction class is the base class for constructing the detector
class DetectorConstruction : public G4VUserDetectorConstruction
{
private:
  // This function defines needed materials
  void DefineMaterials();
  // This function initializes geometry parameters (lengths, dimensions, etc.)
  void ComputeGeometry();
 
  // Materials to be defined
  G4Material* silicon;
  G4Material* galactic;
  G4Material* TiO2;
  G4Material* vehicule;
  G4Material* paint;
  G4Material* coestruded;
  G4Material* plastic;
  G4Material* fMaterial;
  G4Material* epoxy;
  G4Material* aluminium;

  //Logical Volumes

  G4double halfWorldLength;
  G4LogicalVolume* logicWorld; //global mother volume which all things will be placed
  G4LogicalVolume* DetectorLogic;
  G4LogicalVolume* PlasticLogic;
  G4LogicalVolume* CoatingLogic;	
  G4LogicalVolume* GlueLogic;

  G4LogicalVolume* logicbox_x1;
  G4LogicalVolume* logicbox_x2;
  G4LogicalVolume* logicbox_x3;
  G4LogicalVolume* logicbox_x4;
  G4LogicalVolume* logicbox_x5;
  G4LogicalVolume* logicbox_x6;
  G4LogicalVolume* logicbox_x7;
  G4LogicalVolume* logicbox_x8;	

  G4LogicalVolume* logicbox_y1;
  G4LogicalVolume* logicbox_y2;
  G4LogicalVolume* logicbox_y3;
  G4LogicalVolume* logicbox_y4;
  G4LogicalVolume* logicbox_y5;
  G4LogicalVolume* logicbox_y6;	

  G4LogicalVolume* logicbox_z1;
  G4LogicalVolume* logicbox_z2;
  G4LogicalVolume* logicbox_z3;
  G4LogicalVolume* logicbox_z4;
  G4LogicalVolume* logicbox_z5;
  G4LogicalVolume* logicbox_z6;

  G4LogicalVolume* logiccyl_x;
  G4LogicalVolume* logiccyl_y;
  G4LogicalVolume* logiccyl_z;

  //Phyisical volumes

  G4VPhysicalVolume* physiWorld;
  G4VPhysicalVolume* PlasticPhys;
  G4VPhysicalVolume* DetectorPhys;	
  G4VPhysicalVolume* CoatingPhys;
  G4VPhysicalVolume* GluePhys;

  G4VPhysicalVolume* physicbox_x1;
  G4VPhysicalVolume* physicbox_x2;
  G4VPhysicalVolume* physicbox_x3;
  G4VPhysicalVolume* physicbox_x4;
  G4VPhysicalVolume* physicbox_x5;
  G4VPhysicalVolume* physicbox_x6;
  G4VPhysicalVolume* physicbox_x7;
  G4VPhysicalVolume* physicbox_x8; 

  G4VPhysicalVolume* physicbox_y1;
  G4VPhysicalVolume* physicbox_y2;
  G4VPhysicalVolume* physicbox_y3;
  G4VPhysicalVolume* physicbox_y4;
  G4VPhysicalVolume* physicbox_y5;
  G4VPhysicalVolume* physicbox_y6;

  G4VPhysicalVolume* physicbox_z1;
  G4VPhysicalVolume* physicbox_z2;
  G4VPhysicalVolume* physicbox_z3;
  G4VPhysicalVolume* physicbox_z4;
  G4VPhysicalVolume* physicbox_z5;
  G4VPhysicalVolume* physicbox_z6;

  G4VPhysicalVolume* physiccyl_x;
  G4VPhysicalVolume* physiccyl_y;
  G4VPhysicalVolume* physiccyl_z;
	

  //Surfaces
  G4LogicalBorderSurface* PlasticSurface;
  G4LogicalBorderSurface* ExternalSurface; 
  G4LogicalBorderSurface* GlueSurface;
 
  //detector
  G4ThreeVector posDet;
  G4double DetLength;
  G4double DetX;
  G4double DetY;
  G4double DetZ;
  
  //plastic scintillator
  G4ThreeVector posCent;
  G4double CentX;
  G4double CentY;
  G4double CentZ;

//

G4double fractionmass;

TargetSD* fTargetSD;
  
public:
  DetectorConstruction();// Constructor
  ~DetectorConstruction();// Destructor

  
  // REQUIRED Construct the entire geometry of the setup REQUIRED
  G4VPhysicalVolume* Construct();
  G4Material* GetMaterial()   {return fMaterial;};
  //virtual void ConstructSDandField();
};
#endif


