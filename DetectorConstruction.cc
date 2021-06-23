//Include all the prototypes and data members we have defined in DetectorConstruction.hh
#include "DetectorConstruction.hh"
//Included to make a some geometrical objects
#include "G4GeometryTolerance.hh"
#include "G4GeometryManager.hh"
//A list of a box, a tube, logical volume, world placement, and being able make copies of objects, or replicating objects
#include "G4Box.hh"
#include "G4SubtractionSolid.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4OpticalSurface.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
//Allows use of units in the calculations
#include "G4SystemOfUnits.hh"
//Materials manager to specify materials
#include "G4Material.hh"
#include "G4NistManager.hh" //materials from the NIST database 
#include "G4Element.hh"
#include "G4Isotope.hh"
#include "G4MaterialTable.hh"
//Allows us to specficy what the objects will look like in the visualization
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "globals.hh"
#include "G4SDManager.hh"
//sensitive detectors
//#include "B2TrackerSD.hh"
#include "G4VSensitiveDetector.hh"
#include "TargetSD.hh"
//#include "B2TrackerHit.hh"
#include "G4MultiUnion.hh"
//To create ascii file
#include <fstream>

//Constructor with a default value of silicon
DetectorConstruction::DetectorConstruction():fMaterial(0), fTargetSD(0)
{
  //Execute the following on the creation of the detector
  ComputeGeometry();	
  DefineMaterials();
 
}

//Generic destructor
DetectorConstruction::~DetectorConstruction(){}


//*********************************************GEOMETRIA**********************************************************
//Executed each time a detector is constructed
void DetectorConstruction::ComputeGeometry() 
{//This function defines the geometry dimensions. Variables are defined in .hh

  //The coordinates origin is located at the center of the scintillator (world)
  
  halfWorldLength = 5*m; 
  
  // The detector measures dxd (and arbitrary thickness)
  // The plastic will be a cube with side L, with L = 3 * d.
    
  DetLength=2*cm;

  DetX = 1*cm;
  DetY = DetLength;
  DetZ = DetLength;

    CentX = 3*DetLength;
  CentY = 3*DetLength; 
  CentZ = 3*DetLength;
 
  posCent = G4ThreeVector(0,0,0); //plastic in center of coordinates
  posDet  = G4ThreeVector(CentX/2.,0.,0.);// detector in the center of one detector  face (plano x=CentX/2)


}

//*********************************************MATERIALS**********************************************************
//Executed each time a detector is constructed
void DetectorConstruction::DefineMaterials()
{
 
  G4double a; // mass of a mole
  G4double z; // mean number of protons
  G4String name, symbol;
 

  //Get Materials from a NIST database (DB)
  //En ./source/materials/src/G4NistMaterialBuilder.cc
  G4NistManager* man = G4NistManager::Instance();
  man->SetVerbose(1); //display nothing about the materials in the UI  
	

  //PLASTIC SCINTILLATOR ------------------------------------------

  plastic =  man->FindOrBuildMaterial("G4_POLYSTYRENE");//C8H8
  G4double densityplastic=1.08*g/cm3;
	
  //BC418 "G4_POLYVINYL_TOLUENE" also very used
	
  //cOULD ALSO BE DEFINED
  /*  	
  //BC404 http://hypernews.slac.stanford.edu/HyperNews/geant4/get/AUX/2011/06/01/03.27-76314-2DetectorConstruction.txt
  G4Element* C  = new G4Element("Carbon"  ,symbol="C" , z= 6., a= 12.01*g/mole);
  G4Element* H  = new G4Element("Hydrogen",symbol="H" , z= 1., a= 1.01*g/mole);	
  density = 1.032*g/cm3;
  G4Material* Sci404 = new G4Material("Scintillator", density = 1.032*g/cm3, ncomponents=2);
  Sci404->AddElement(C, natoms=9);
  Sci404->AddElement(H, natoms=10);
  */	
	

	
  //galactic empty----------------------------------------------------
 
   galactic  = man->FindOrBuildMaterial("G4_Galactic");
 
   //or defined
   /*
   G4double temperature, pressure;
   density     = 1.e-25*g/cm3; // data from CLHEP
   pressure    = 3.e-18*pascal;
   temperature = 2.73*kelvin;
   galactic =   new G4Material(name="Galactic", z=1., a=1.01*g/mole, density, kStateGas,temperature,pressure);  
   */
	
  //TiO2--------------------------------------------------------------
 
  // Definition of element Oxygen
  G4Element* elO = new G4Element(name="Oxygen", symbol="O", z=8., a=16.00*g/mole);
  // Definition of element Titanium
  G4Element* elTi = new G4Element(name="Titanium", symbol="Ti", z=22., a=47.87*g/mole);
    
  G4double densityTiO2 = 4.23*g/cm3;

  TiO2 = new G4Material("TiO2",  //its name 
                       densityTiO2,    //its density
                       2);         //number of components

  //Add Element for Material "TiO2" specifiyng the number of each element
  TiO2->AddElement(elO,2);
  TiO2->AddElement(elTi,1); 


  //PAINT (TiO2+vehiculE(C4H6O2 -> Polyvinyl acetate)) ----------------------------------

  // Definition of element Hidrogen
  G4Element* elH = new G4Element(name="Hidrogen", symbol="H", z=1., a=1*g/mole);
  // Definition of element Carbon
  G4Element* elC = new G4Element(name="Carbon", symbol="C", z=6., a=12.00*g/mole);

  // http://scientificpolymer.com/technical-library/refractive-index-of-polymers-by-index/
  G4double densityvehicule = 1.19*g/cm3;
  vehicule = new G4Material("Vehicule",densityvehicule, 3); 
  vehicule->AddElement(elC,4);  
  vehicule->AddElement(elH,6);  
  vehicule->AddElement(elO,2);   

  //paint and vehicule + TiO2 
    //https://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides/ForApplicationDeveloper/html/ch04s02.html
  G4double densitypaint = densityTiO2*15/100.+densityvehicule*(100.-15)/100.;
  //http://publications.rwth-aachen.de/record/667646/files/667646.pdf pag 203
  paint = new G4Material("Paint",densitypaint, 2);   
  paint->AddMaterial(TiO2, fractionmass=15*perCent);
  paint->AddMaterial(vehicule, fractionmass=(100.-15)*perCent);


  //coestruded  -----------------------------------------------------------

  //https://geant4.web.cern.ch/geant4/UserDocumentation/UsersGuides/ForApplicationDeveloper/html/ch04s02.html
  G4double densitycoestruded = densityTiO2*15/100.+densityplastic*(100.-15)/100.;
  coestruded = new G4Material("Coestruded",densitycoestruded, 2);   
  coestruded->AddMaterial(TiO2, fractionmass=15*perCent);
  coestruded->AddMaterial(plastic, fractionmass=(100.-15)*perCent);


  //Detector -------------------------------------------------------------------------
  silicon = man->FindOrBuildMaterial("G4_Si");
  

  //Glue----------------------------------------------------------------------------
  //Epoxy definition C21H25ClO5
  //http://www.ashland.com/file_source/Ashland/Documents/Sustainability/rc_bisphenol_a_epoxy_diacrylate_pss.pdf
  //http://mstg.unile.it/compositi/line1.html

  G4Element* elCl = new G4Element(name="Cloro", symbol="Cl", z=17., a=35.45*g/mole);
  G4double densityepoxy = 1.195*g/cm3;
  epoxy = new G4Material("epoxy",  densityepoxy, 4);
  epoxy->AddElement(elC,21);
  epoxy->AddElement(elH,25); 
  epoxy->AddElement(elCl,1);
  epoxy->AddElement(elO,5);

  //Inertial wheels Al----------------------------------------------------------------------------
  // Definition of material Aluminimum
  G4double densityAl = 2.7*g/cm3;
  aluminium = new G4Material(name="Aluminum", z=13.,a=26.98*g/mole, densityAl);

  //-----------------------------------------------------------------------------------------------------------------
  //Listar materiales 
  //G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;

}
 
 
//Construct() This function is called by G4RunManager during the initialization. The user must implement the detector and world geometry in the 
//  Construct() function. The function must return the pointer to the physical volume of the world. 

G4VPhysicalVolume* DetectorConstruction::Construct()
{

  //Set limits for the simulation (World volume)
  G4bool checkOverlaps = true;
  G4GeometryManager::GetInstance()->SetWorldMaximumExtent(2.*halfWorldLength);
  

    //0-> no coating
  //1-> coestrudado (plastic and dopant + TiO2);
  //2-> pintura (vehículo + TiO2); -> should be revised, it is not correct now
  G4int SelectCoating=1;

  //With or Wo Glue
  G4bool WithGlue=0;

  //With or Wo Telescope (bars and wheels)
  G4bool WithTelescope=0;


  //World********************************************************************************************	

  //Make a solid box, or the shape of the world
  G4Box* solidWorld = new G4Box("world",halfWorldLength,halfWorldLength,halfWorldLength);  
  
  //Make a logical volume filled with galactic vaccuum called world located at (0,0,0)INTEGERS
  logicWorld = new G4LogicalVolume(solidWorld, galactic, "World", 0, 0, 0);
  
  //Must place the World Physical volume unrotated at (0,0,0)
  physiWorld = new G4PVPlacement(0,// no rotation INTEGERS
			    G4ThreeVector(),// at (0,0,0)
			    logicWorld,// its logical volume
       			    "Galactic",// its name
			    0,// its mother volume doesnt exsist
			    false,// no boolean operations
			    0, // copy number 
			    checkOverlaps);


  //Scintillator and detector dimensions*************************************************************

  //define Scintillator dimensions
  G4double halfCentSizeX = CentX/2.;
  G4double halfCentSizeY = CentY/2.;
  G4double halfCentSizeZ = CentZ/2.;

  // define detector dimensions
  G4double halfDetSizeX = DetX/2.;
  G4double halfDetSizeY = DetY/2.;
  G4double halfDetSizeZ = DetZ/2.;


 //Material properties**********************************************************************************
	
 //G4double h_Planck = 4.13566733*1e-15;//eV*s
 //G4double c_light = 299792458;//m/s

 //Propiedades opticas del Plástico-----------------------------------------------------------------
  
 //optical properties of BC418 (en el enlace vienen otros valores para el BC444).
 //http://hypernews.slac.stanford.edu/HyperNews/geant4/get/AUX/2011/06/01/03.27-76314-2DetectorConstruction.txt   
 //Explanation (ver 5.2.5.2)
 //http://geant4.cern.ch/G4UsersDocuments/UsersGuides/ForApplicationDeveloper/html/TrackingAndPhysics/physicsProcess.html
 //http://wiki.opengatecollaboration.org/index.php/Users_Guide_V7.2:Generating_and_tracking_optical_photons
 //http://www.crystals.saint-gobain.com/uploadedFiles/SG-Crystals/Documents/SGC%20BC418-420-422%20Data%20Sheet.pdf
 //BC418 ENERGY EMISSION SPECTRUM
 //http://www.eljentechnology.com/products/plastic-scintillators/ej-228-ej-230
 //http://www.crystals.saint-gobain.com/sites/imdf.crystals.com/files/documents/sgc-bc418-420-422-data-sheet_69699.pdf 


  //Scintillation spectrum of POLYSTYERNE DOPADO with PPO and POPOP (1 % and 0.03 % by weight, respectively).
  //Measurement of neutrino interaction and three flavour neutrino interactions in the T2K experiment fIG 3.8
   //http://hypernews.slac.stanford.edu/HyperNews/geant4/get/opticalphotons/442/1.html?inline=-1	
  const G4int numEntries =8;   
  G4double photonEnergies[numEntries] = {3.01*eV,2.92*eV,2.76*eV,2.61*eV,2.48*eV,2.25*eV,2.01*eV,1.91*eV};   
  //These values specify the probability that a photon with the given energy is emitted.
  //http://wiki.opengatecollaboration.org/index.php/Users_Guide_V7.2:Generating_and_tracking_optical_photon
  G4double scintweight[numEntries] = {0.90,0.96,0.85,0.60,0.31,0.08,0.05,0.05};
 
  //OPTICAL PROPERTIES
  //http://www.crystals.saint-gobain.com/uploadedFiles/ SG-Crystals/Documents/SGC%20BC418-420-422%20Data%20Sheet.pdf
  G4double refindexplastic=1.58;
  //valores tipicos 150-400 cm http://www.southernscientific.co.uk/data/file/9/5/Table%20of%20Physical%20Constants%20of%20Scintillators.1438855290.pdf
  G4double abslengthplastic=250.0*cm;
  //Rayleigh. Sólo aqúi dice 95m para 400nm at room temperature. G4OpRayleigh.cc lo calcula si no se le da.
  //http://hypernews.slac.stanford.edu/HyperNews/geant4/get/emprocess/333.html?inline=-1   
  G4double raylengthplastic=95.0*m;
 
 G4double PlasticRefractiveIndex[numEntries] = {refindexplastic, refindexplastic, refindexplastic, refindexplastic, refindexplastic, refindexplastic, refindexplastic, refindexplastic};
 G4double PlasticAbsorptionLength[numEntries] = {abslengthplastic, abslengthplastic, abslengthplastic, abslengthplastic, abslengthplastic, abslengthplastic, abslengthplastic, abslengthplastic}; 
  G4double PlasticRayleighLength[numEntries] = {raylengthplastic, raylengthplastic, raylengthplastic, raylengthplastic, raylengthplastic, raylengthplastic, raylengthplastic, raylengthplastic};

  G4MaterialPropertiesTable* mptPlastic = new G4MaterialPropertiesTable();// defining an object of the G4MaterialPropertiesTable
  mptPlastic->AddProperty("RINDEX",photonEnergies,PlasticRefractiveIndex,numEntries);
  mptPlastic->AddProperty("ABSLENGTH",photonEnergies,PlasticAbsorptionLength,numEntries);//necesario para que ocurra G4OpAbsorption     
  mptPlastic->AddProperty("RAYLEIGH",photonEnergies,PlasticRayleighLength,numEntries); 
  
  //SCINTILLATION PROPERTIES
  //the number of photons that is emitted per amount of energy absorbed, or, more precisely, it gives the expectation value of this number, since the real number of emitted photons follows a normal distribution whose mean is SCINTILLATIONYIELD and sigma is RESOLUTIONSCALE
  //the actual number of emitted photons during a step fluctuates around the mean number of photons with a width given by ResolutionScale*sqrt(MeanNumberOfPhotons). The average light yield, MeanNumberOfPhotons, has a linear dependence on the local energy  
  //RESOLUTIONSCALE = 1 energy resolution follows a Poisson distribution. To broaden the statistical distribution of generated photons. the variance of this normal distribution is RESOLUTION-SCALE times this expectation value. 1 gives a Poisson law, 0 is "perfect" resolution
  mptPlastic->AddConstProperty("SCINTILLATIONYIELD", 8000./MeV); //POLYSTYERENE VALUE
  mptPlastic->AddConstProperty("RESOLUTIONSCALE", 1.0);// with

  //Organic plastic fast component dominates
  //ONLY ONE COMPONENT
  mptPlastic->AddProperty("FASTCOMPONENT", photonEnergies, scintweight, numEntries);  
  mptPlastic->AddConstProperty("FASTTIMECONSTANT", 3.6*ns); //POLYSTYERENE VALUE
  //Two components
  //Then comment in PhysicsList.cc: Scint->SetScintillationYieldFactor(1.) Scint->SetScintillationExcitationRatio(1.); ya que sino fija solo componente rápida. 
  //mptPlastic->AddProperty("SLOWCOMPONENT", photonEnergies, scintweight, numEntries);
  //mptPlastic->AddConstProperty("SLOWTIMECONSTANT", 14.2*ns); //Se construyen POLYSTYERENE hasta con 300-400 ns (http://www.amcrys.com/details.html?cat_id=146&id=4286); 
  //YIELDRATIO->The relative strength of the fast component as a fraction of total scintillation yield (fast vs slow)
  //ONLY ONE COMPONENT
  mptPlastic->AddConstProperty("YIELDRATIO", 1.0);
  //TWO COMPONENTS
  //mptPlastic->AddConstProperty("YIELDRATIO", 0.73);

  //Birks law
  //Es un corrección que se hace al light yield per path length. The constant depends on the scintillator type
  //La luz producida en un centelleador no es directamente proporcional a la Edep por la partıcula que lo atraviesa, sino que sigue la ley de Birk.
  //Birks' Law has mostly been tested for organic scintillators. Its applicability to inorganic scintillators is debated 
  //Aquí no lo usa: //http://hypernews.slac.stanford.edu/HyperNews/geant4/get/AUX/2011/06/01/03.27-76314-2DetectorConstruction.txt
  //Aquí sí lo usa: https://www-zeuthen.desy.de/geant4/g4course2011/day2/6_opticalphysics/DetectorConstruction_8cc-source.html
  //Aqui el código de cómo calcularlo y lo meten en SteppingAction.cc: http://hypernews.slac.stanford.edu/HyperNews/geant4/get/opticalphotons/311.html?inline=-1
  //https://arxiv.org/pdf/0911.3041v1.pdf C8H8 vale 9.0 × 10−3 g MeV−1 cm-2, C16H18 vale 6.8 × 10−3 g MeV−1 cm-2, C9H12 vale 35 × 10−3 g MeV−1 cm-2
  plastic->GetIonisation()->SetBirksConstant(0.126*mm/MeV);//for polystyerene (0.0126*mm/MeV para agua)

// set these to material properties
  plastic->SetMaterialPropertiesTable(mptPlastic);


  //Propiedades opticas del pintura ----------------------------------------------------------------
  
  //(Polyvinyl acetate)
  //https://www.chemours.com/Titanium_Technologies/en_US/assets/downloads/Ti-Pure-for-coatings-overview.pdf
  G4double refindexvehicule=1.47;
  //indice TiO2
  //https://www.chemours.com/Titanium_Technologies/en_US/assets/downloads/Ti-Pure-for-coatings-overview.pdf
  G4double refindexTiO2=2.55;//Titanium dioxide (anatase) 2.55 ;Titanium dioxide (rutile) 2.73
  //http://publications.rwth-aachen.de/record/667646/files/667646.pdf pag 203
  G4double refindexpaint = refindexvehicule*(100.-15)/100.+refindexTiO2*15/100.;
  
  G4double RefractiveIndexPaint[numEntries] = {refindexpaint, refindexpaint, refindexpaint, refindexpaint, refindexpaint, refindexpaint, refindexpaint, refindexpaint};
  
  G4MaterialPropertiesTable* mptPaint = new G4MaterialPropertiesTable();	
  mptPaint->AddProperty("RINDEX",photonEnergies,RefractiveIndexPaint,numEntries);    	  
  paint->SetMaterialPropertiesTable(mptPaint);

 //oPTICAL PROPERTIES OF STRUDED ----------------------------------------------------------------

 // http://publications.rwth-aachen.de/record/667646/files/667646.pdf pag 203
  G4double refindexcoestruded = refindexplastic*(100.-15)/100.+refindexTiO2*15/100.;

  G4double RefractiveIndexCoestruded[numEntries] ={refindexcoestruded, refindexcoestruded, refindexcoestruded, refindexcoestruded, refindexcoestruded, refindexcoestruded, refindexcoestruded, refindexcoestruded};

  G4MaterialPropertiesTable* mptCoestruded = new G4MaterialPropertiesTable();	
  mptCoestruded->AddProperty("RINDEX",photonEnergies,RefractiveIndexCoestruded,numEntries);    	
  coestruded->SetMaterialPropertiesTable(mptCoestruded);


  //oPTICAL PROPOERTIES OF Si	----------------------------------------------------------------
  
  //http://refractiveindex.info/?shelf=main&book=Si&page=Green
  //En el enlace viene también la reflactancia R(angulo de incidencia y polarización)
  //Opcion 1: valor medio de n.
  G4double refractiveIndexSi[numEntries] = {3.4860,3.4860,3.4860,3.4860,3.4860,3.4860,3.4860,3.4860};
  //Opcion 2: n(ldo)->cambia mucho a ldo<1micra. 
  //http://www.filmetrics.com/refractive-index-databas
  //G4double refractiveIndexSi[numEntries] = {6.089, 6.810, 5.940, 5.570, 5.085, 4.790, 4.580, 4.416};
 
  G4MaterialPropertiesTable* mptSi = new G4MaterialPropertiesTable();
  mptSi->AddProperty("RINDEX",photonEnergies,refractiveIndexSi,numEntries);	  
  silicon->SetMaterialPropertiesTable(mptSi);

	
  // optical properties of galactic vacuum----------------------------------------------------------------

  //http://hypernews.slac.stanford.edu/HyperNews/geant4/get/AUX/2011/06/01/03.27-76314-2DetectorConstruction.txt
  G4double refractiveIndexVacuum[numEntries] = {1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00};
  // defining, an object of the G4MaterialPropertiesTable specific to Vacuum
  G4MaterialPropertiesTable* mptVacuum = new G4MaterialPropertiesTable();
  mptVacuum->AddProperty("RINDEX",photonEnergies,refractiveIndexVacuum,numEntries);
  galactic->SetMaterialPropertiesTable(mptVacuum);
  


  //optical properties of GLUE----------------------------------------------------------------

  //http://www.amstechnologies.com/fileadmin/amsmedia/downloads/5187_nttatopticaladhesives.pdf
  //Tomo el de un epoxy. De otro viene n(lambda), pero como este valor fijo par el epoxy que dan.
  G4double refractiveIndexEpoxy[numEntries] = {1.627,1.627,1.627,1.627,1.627,1.627,1.627,1.627};
  G4MaterialPropertiesTable* mptGlue = new G4MaterialPropertiesTable();
  mptGlue->AddProperty("RINDEX",photonEnergies,refractiveIndexEpoxy,numEntries);	  
  epoxy->SetMaterialPropertiesTable(mptGlue);




  //Make solids***********************************************************************************************************


  //da igual que a todos estos que son G4Box o G4SubtractionSolid, les defina al inicio como G4VSolid	
  G4Box* DetectorSolid = new G4Box("SensorSolid", halfDetSizeX, halfDetSizeY, halfDetSizeZ);
  G4Box* PlasticSolid = new G4Box("PlasticSolid",halfCentSizeX,halfCentSizeY,halfCentSizeZ);

  //https://books.google.es/books?id=JJV71U3ZIw8C&pg=PA241&lpg=PA241&dq=espesor+pintura+tio2&source=bl&ots=3V9Gsa7Khx&sig=CFB_-45iXXQMv4r-wvyS40AkP7I&hl=es&sa=X&ved=0ahUKEwjtjZDp8L3NAhXDVhoKHceADB8Q6AEIKTAB#v=onepage&q=espesor%20pintura%20tio2&f=false
  G4double layerwidth=0.025*cm;
  G4double numlayers=2; //2 se fija desde el script
  G4double coatingwidth=layerwidth*numlayers; //grosor total pintura/coestrudado

  G4Box* CoatingSolid0=new G4Box("CoatingSolid0",halfCentSizeX,halfCentSizeY,halfCentSizeZ);
  G4Box* CoatingSolid1 = new G4Box("CoatingSolid1",halfCentSizeX+coatingwidth,halfCentSizeY+coatingwidth,halfCentSizeZ+coatingwidth);
  G4SubtractionSolid* CoatingSolid2 = new G4SubtractionSolid("CoatingSolid2", CoatingSolid1, CoatingSolid0);
  G4Box* CoatingSolid3=new G4Box("CoatingSolid3",coatingwidth/2,halfDetSizeY,halfDetSizeZ);	
  G4SubtractionSolid* CoatingSolid = new G4SubtractionSolid("CoatingSolid", CoatingSolid2, CoatingSolid3, 0, G4ThreeVector(halfCentSizeX+coatingwidth/2,0.,0.));//el 0 es porque no hay rotación, el G4ThreeVector es la traslación del CoatingSolid3



  //GLUE SHOULD BE REVISED
  G4double halfgluewidth=0.0025*cm; 
  if(2*halfgluewidth > coatingwidth) G4cout << G4endl <<G4endl << G4endl << G4endl << "ERROR: PEGAMENTO DEMASIADO GRUESO" <<G4endl << G4endl << G4endl;
    G4Box * GlueSolid = new G4Box("GlueSolid",halfgluewidth, halfDetSizeY,halfDetSizeZ);

  
  //tELESCOPE STRUCTURE
  G4double halfwidthbar = 0.25*cm;
  G4double halflongbar = 5*cm;
  G4Box* box_x1= new G4Box("Box_x1",halflongbar,halfwidthbar,halfwidthbar);
  G4Box* box_x2= new G4Box("Box_x2",halflongbar,halfwidthbar,halfwidthbar);
  G4Box* box_x3= new G4Box("Box_x3",halflongbar,halfwidthbar,halfwidthbar);
  G4Box* box_x4= new G4Box("Box_x4",halflongbar,halfwidthbar,halfwidthbar);
  G4Box* box_x5= new G4Box("Box_x5",halflongbar,halfwidthbar,halfwidthbar);
  G4Box* box_x6= new G4Box("Box_x6",halflongbar,halfwidthbar,halfwidthbar);
  G4Box* box_x7= new G4Box("Box_x7",halflongbar,halfwidthbar,halfwidthbar);
  G4Box* box_x8= new G4Box("Box_x8",halflongbar,halfwidthbar,halfwidthbar);

  G4Box* box_y1= new G4Box("Box_y1",halfwidthbar,halflongbar,halfwidthbar);
  G4Box* box_y2= new G4Box("Box_y2",halfwidthbar,halflongbar,halfwidthbar);
  G4Box* box_y3= new G4Box("Box_y3",halfwidthbar,halflongbar,halfwidthbar);
  G4Box* box_y4= new G4Box("Box_y4",halfwidthbar,halflongbar,halfwidthbar);
  G4Box* box_y5= new G4Box("Box_y5",halfwidthbar,halflongbar,halfwidthbar);
  G4Box* box_y6= new G4Box("Box_y6",halfwidthbar,halflongbar,halfwidthbar);

  G4Box* box_z1= new G4Box("Box_z1",halfwidthbar,halfwidthbar,halflongbar);
  G4Box* box_z2= new G4Box("Box_z2",halfwidthbar,halfwidthbar,halflongbar);
  G4Box* box_z3= new G4Box("Box_z3",halfwidthbar,halfwidthbar,halflongbar);
  G4Box* box_z4= new G4Box("Box_z4",halfwidthbar,halfwidthbar,halflongbar);
  G4Box* box_z5= new G4Box("Box_z5",halfwidthbar,halfwidthbar,halflongbar);
  G4Box* box_z6= new G4Box("Box_z6",halfwidthbar,halfwidthbar,halflongbar);

  G4Tubs* cyl_x = new G4Tubs("Cylinder1",0,1.5*cm,0.5*cm,0.,360.*deg); //name, inner radius, outer radius, halflengthz, startangle, endangle
  G4Tubs* cyl_y = new G4Tubs("Cylinder2",0,1.5*cm,0.5*cm,0.,360.*deg);
  G4Tubs* cyl_z = new G4Tubs("Cylinder3",0,1.5*cm,0.5*cm,0.,360.*deg);
 
  //Make logical volumes***************************************************************************************

  DetectorLogic = new G4LogicalVolume(DetectorSolid,// its solid
				   silicon,//its material
				  "DetectorLogic",//its name
				   0,0,0);//estos ceros son los valores por defecto.

  PlasticLogic = new G4LogicalVolume(PlasticSolid,//its solid
				 plastic,//its material
				 "PlasticLogic");//its name

  if(SelectCoating==1) CoatingLogic = new G4LogicalVolume(CoatingSolid, coestruded, "CoatingLogic");

  if(SelectCoating==2) CoatingLogic = new G4LogicalVolume(CoatingSolid, paint, "CoatingLogic");

 
  if(WithGlue) GlueLogic = new G4LogicalVolume(GlueSolid, epoxy, "GlueLogic");

  if(WithTelescope) {

	logicbox_x1 = new G4LogicalVolume(box_x1, aluminium, "logicbox_x1");
	logicbox_x2 = new G4LogicalVolume(box_x2, aluminium, "logicbox_x2");
	logicbox_x3 = new G4LogicalVolume(box_x3, aluminium, "logicbox_x3");
	logicbox_x4 = new G4LogicalVolume(box_x4, aluminium, "logicbox_x4");
	logicbox_x5 = new G4LogicalVolume(box_x5, aluminium, "logicbox_x5");
	logicbox_x6 = new G4LogicalVolume(box_x6, aluminium, "logicbox_x6");
	logicbox_x7 = new G4LogicalVolume(box_x7, aluminium, "logicbox_x7");
	logicbox_x8 = new G4LogicalVolume(box_x8, aluminium, "logicbox_x8");
	
	logicbox_y1 = new G4LogicalVolume(box_y1, aluminium, "logicbox_y1");
	logicbox_y2 = new G4LogicalVolume(box_y2, aluminium, "logicbox_y2");
	logicbox_y3 = new G4LogicalVolume(box_y3, aluminium, "logicbox_y3");
	logicbox_y4 = new G4LogicalVolume(box_y4, aluminium, "logicbox_y4");
	logicbox_y5 = new G4LogicalVolume(box_y5, aluminium, "logicbox_y5");
	logicbox_y6 = new G4LogicalVolume(box_y6, aluminium, "logicbox_y6");

	logicbox_z1 = new G4LogicalVolume(box_z1, aluminium, "logicbox_z1");
	logicbox_z2 = new G4LogicalVolume(box_z2, aluminium, "logicbox_z2");
	logicbox_z3 = new G4LogicalVolume(box_z3, aluminium, "logicbox_z3");
	logicbox_z4 = new G4LogicalVolume(box_z4, aluminium, "logicbox_z4");
	logicbox_z5 = new G4LogicalVolume(box_z5, aluminium, "logicbox_z5");
	logicbox_z6 = new G4LogicalVolume(box_z6, aluminium, "logicbox_z6");

	logiccyl_x = new G4LogicalVolume(cyl_x, aluminium, "logiccyl_x");
	logiccyl_y = new G4LogicalVolume(cyl_y, aluminium, "logiccyl_y");
	logiccyl_z = new G4LogicalVolume(cyl_z, aluminium, "logiccyl_z");

   }

  //Make the physical volumes ***************************************************************************************

  PlasticPhys = new G4PVPlacement(0, //no rotation
			     posCent,//location
			     PlasticLogic, //its logical volume
			     "Plastic",//its name
			     logicWorld,//its mother volume 
			     false, //no boolean operation
			     0); //copy number

//Nota: Otra opción que da exactamente lo mismo es arriba poner como CoatingSolid lo que ahora es CoatingSolid1 (cubo todo lleno de coating, bueno con el hueco del detector además) y aquí poner al PlasticPhys como madre CoatingLogic. Al intersectar plastic entonces al coating lo sustituye. El hijo intersecta a la madre.
 
  if(SelectCoating!=0) CoatingPhys = new G4PVPlacement(0, //no rotation
			     posCent,//location
			     CoatingLogic, //its logical volume
			     "Coating",//its name
			     logicWorld,//its mother volume
			     false, //no boolean operation
			     0); //copy number
 

  if(WithGlue) GluePhys = new G4PVPlacement(0, 
                                  G4ThreeVector(halfCentSizeX+halfgluewidth,0.,0.), //next to the plastic
                                  GlueLogic, //its logical volume
                                  "Glue", 
                                  logicWorld, 
                                  false, 
                                  checkOverlaps,
  				  0);


 //Detector -> no coating between detector and plastic

 if(!WithGlue) DetectorPhys = new G4PVPlacement(0, //no rotation
                             G4ThreeVector(halfCentSizeX+halfDetSizeX,0.,0.), //location (next to the plastic)
                             DetectorLogic, //its logical volume
                             "Sensor", //its name
                             logicWorld, //its mother volume 
                             false, //no boolean operation
                             checkOverlaps,//G4bool checkOverlaps = true;
  			     0);

 if(WithGlue) DetectorPhys = new G4PVPlacement(0, //no rotation
                             G4ThreeVector(halfCentSizeX+2*halfgluewidth+halfDetSizeX,0.,0.), //location (next to the plastic)
                             DetectorLogic, //its logical volume
                             "Sensor", //its name
                             logicWorld, //its mother volume 
                             false, //no boolean operation
                             checkOverlaps,//G4bool checkOverlaps = true;
  			     0);


  if(WithTelescope) {

         //detloc -> design parameter
	 G4double detloc = 2.*cm;

	 physicbox_x1 = new G4PVPlacement(0, G4ThreeVector(halflongbar+detloc,halflongbar,halflongbar),
  			     logicbox_x1,"physicbox_x1",logicWorld,false,0);
	 physicbox_x2 = new G4PVPlacement(0,G4ThreeVector(halflongbar+detloc,-halflongbar,halflongbar),
  			     logicbox_x2,"physicbox_x2",logicWorld,false,0); 
	 physicbox_x3 = new G4PVPlacement(0,G4ThreeVector(halflongbar+detloc,halflongbar,-halflongbar),
  			     logicbox_x3,"physicbox_x3",logicWorld,false,0); 
	 physicbox_x4 = new G4PVPlacement(0,G4ThreeVector(halflongbar+detloc,-halflongbar,-halflongbar),
  			     logicbox_x4,"physicbox_x4",logicWorld,false,0); 
	 physicbox_x5 = new G4PVPlacement(0,G4ThreeVector(-halflongbar+detloc,halflongbar,halflongbar),
  			     logicbox_x5,"physicbox_x5",logicWorld,false,0); 
	 physicbox_x6 = new G4PVPlacement(0,G4ThreeVector(-halflongbar+detloc,-halflongbar,halflongbar),
  			     logicbox_x6,"physicbox_x6",logicWorld,false,0); 
 	 physicbox_x7 = new G4PVPlacement(0,G4ThreeVector(-halflongbar+detloc,halflongbar,-halflongbar),
  			     logicbox_x7,"physicbox_x7",logicWorld,false,0); 
	 physicbox_x8 = new G4PVPlacement(0,G4ThreeVector(-halflongbar+detloc,-halflongbar,-halflongbar),
  			     logicbox_x8,"physicbox_x8",logicWorld,false,0); 


	 physicbox_y1 = new G4PVPlacement(0, G4ThreeVector(2.*halflongbar+detloc,0.,halflongbar),
  			     logicbox_y1,"physicbox_y1",logicWorld,false,0);
	 physicbox_y2 = new G4PVPlacement(0, G4ThreeVector(2.*halflongbar+detloc,0.,-halflongbar),
  			     logicbox_y2,"physicbox_y2",logicWorld,false,0);
	 physicbox_y3 = new G4PVPlacement(0, G4ThreeVector(-2.*halflongbar+detloc,0.,halflongbar),
  			     logicbox_y1,"physicbox_y3",logicWorld,false,0);
	 physicbox_y4 = new G4PVPlacement(0, G4ThreeVector(-2.*halflongbar+detloc,0.,-halflongbar),
  			     logicbox_y2,"physicbox_y4",logicWorld,false,0);
	 physicbox_y5 = new G4PVPlacement(0, G4ThreeVector(detloc,0.,halflongbar),
  			     logicbox_y1,"physicbox_y5",logicWorld,false,0);
	 physicbox_y6 = new G4PVPlacement(0, G4ThreeVector(detloc,0.,-halflongbar),
  			     logicbox_y2,"physicbox_y6",logicWorld,false,0);

	 physicbox_z1 = new G4PVPlacement(0, G4ThreeVector(2.*halflongbar+detloc,halflongbar,0.),
  			     logicbox_z1,"physicbox_z1",logicWorld,false,0);
	 physicbox_z2 = new G4PVPlacement(0, G4ThreeVector(2.*halflongbar+detloc,-halflongbar,0.),
  			     logicbox_z2,"physicbox_z2",logicWorld,false,0);
	 physicbox_z3 = new G4PVPlacement(0, G4ThreeVector(detloc,halflongbar,0.),
  			     logicbox_z3,"physicbox_z3",logicWorld,false,0);
	 physicbox_z4 = new G4PVPlacement(0, G4ThreeVector(detloc,-halflongbar,0.),
  			     logicbox_z4,"physicbox_z4",logicWorld,false,0);
	 physicbox_z5 = new G4PVPlacement(0, G4ThreeVector(-2.*halflongbar+detloc,halflongbar,0.),
  			     logicbox_z5,"physicbox_z5",logicWorld,false,0);
	 physicbox_z6 = new G4PVPlacement(0, G4ThreeVector(-2.*halflongbar+detloc,-halflongbar,0.),
  			     logicbox_z6,"physicbox_z6",logicWorld,false,0);

	physiccyl_x = new G4PVPlacement(0, G4ThreeVector(halflongbar+detloc,0.,-detloc),
  			     logiccyl_x,"physiccyl_x",logicWorld,false,0);

	G4double phi = 90*deg;
	G4ThreeVector u = G4ThreeVector(0, 0, -1);
	G4ThreeVector v = G4ThreeVector(-std::sin(phi), std::cos(phi),0.);
	G4ThreeVector w = G4ThreeVector( std::cos(phi), std::sin(phi),0.);
	G4RotationMatrix rotm1 = G4RotationMatrix(u, v, w);
	G4RotationMatrix* rotm1Inv = new G4RotationMatrix(rotm1.inverse()); 
	
	physiccyl_y = new G4PVPlacement(rotm1Inv, G4ThreeVector(halflongbar+detloc,-detloc,0.),
			     logiccyl_y,"physiccyl_y",logicWorld,false,0);

	G4double phi2 = 180*deg;
	G4ThreeVector u2 = G4ThreeVector(0, -1, 0);
	G4ThreeVector v2 = G4ThreeVector(-std::sin(phi2), std::cos(phi2),0.);
	G4ThreeVector w2 = G4ThreeVector( std::cos(phi2), std::sin(phi2),0.);
	G4RotationMatrix rotm2 = G4RotationMatrix(u2, v2, w2);
	G4RotationMatrix* rotm2Inv = new G4RotationMatrix(rotm2.inverse()); 

	physiccyl_z = new G4PVPlacement(rotm2Inv, G4ThreeVector(halflongbar,0.,0.),
  			     logiccyl_z,"physiccyl_z",logicWorld,false,0);

  }


 //Make surfaces ******************************************************************
        
 G4double polish_value = 0.5;// same value for all the surfaces 

 if(SelectCoating!=0){

  //Superficie plástico - coestrudado/pintura -------------------------------------------------------------------

 //See Section Boundary Process                         
 //http://geant4.slac.stanford.edu/UsersWorkshop/PDF/Peter/OpticalPhoton.pdf
 //http://geant4.cern.ch/G4UsersDocuments/UsersGuides/ForApplicationDeveloper/html/TrackingAndPhysics/physicsProcess.html 	 

  //http://publications.rwth-aachen.de/record/667646/files/667646.pdf

 
  G4OpticalSurface* OpPlasticSurface = new G4OpticalSurface("PlasticSurface"); 
  OpPlasticSurface -> SetType(dielectric_metal);
  OpPlasticSurface -> SetModel(glisur);//dos modelos, glisur o unified. Default: glisur (más moderno)
  OpPlasticSurface -> SetFinish(ground);//polished by defalut
   OpPlasticSurface -> SetPolish(polish_value);//0 is maximum roughness, 1 perfectamente pulida

  //COESTRUDADO
  //http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.540.5022&rep=rep1&type=pdf
  //http://publications.rwth-aachen.de/record/667646/files/667646.pdf Fig D.3 la R del Tyvek es 97% aprox, lo toma ref[178] Janecek
  G4double reflectivity_coestruded[numEntries] = {0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90};

  //http://publications.rwth-aachen.de/record/667646/files/667646.pdf
    G4double specularlobe_coestruded[numEntries] = {0.85,0.85,0.85,0.85,0.85,0.85,0.85,0.85};
  
  //Paint
  //Reflectivity TiO2 paint Chemours 
  //https://www.chemours.com/Titanium_Technologies/en_US/assets/downloads/Ti-Pure-for-coatings-overview.pdf
  G4double reflectivity_paint[numEntries] = {0.08, 0.18, 0.61, 0.73, 0.90, 0.93, 0.95, 0.96}; 
  
  G4MaterialPropertiesTable* SurfProp = new G4MaterialPropertiesTable();  

  if(SelectCoating==1) SurfProp-> AddProperty("REFLECTIVITY",photonEnergies,reflectivity_coestruded,numEntries);     
  if(SelectCoating==1) SurfProp -> AddProperty("SPECULARLOBECONSTANT",photonEnergies,specularlobe_coestruded,numEntries);
  if(SelectCoating==2) SurfProp -> AddProperty("REFLECTIVITY",photonEnergies,reflectivity_paint,numEntries);     


  //G4double efficiency[numEntries] = {0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5}; 
  //SurfProp->AddProperty("EFFICIENCY",photonEnergies,efficiency,numEntries);  


  OpPlasticSurface -> SetMaterialPropertiesTable(SurfProp);


  PlasticSurface = new G4LogicalBorderSurface("PlasticSurface",PlasticPhys,CoatingPhys,OpPlasticSurface);//

  //Superficie pintura/coestrudado->vacio -------------------------------------------------------------------

  G4OpticalSurface* OpExternalSurface = new G4OpticalSurface("ExternalSurface");
  OpExternalSurface->SetType(dielectric_dielectric);//el paso es al vacio que es dielectric. Con la anteior metal esta debedar igual
  OpExternalSurface->SetModel(glisur);
  OpExternalSurface->SetFinish(ground);
  OpExternalSurface->SetPolish(polish_value); 
  
  G4MaterialPropertiesTable* ExternalProp = new G4MaterialPropertiesTable();   
  OpExternalSurface-> SetMaterialPropertiesTable(ExternalProp);

  ExternalSurface = new G4LogicalBorderSurface("ExternalSurface",CoatingPhys,physiWorld,OpExternalSurface);

 }

 if(SelectCoating==0){

  G4OpticalSurface* OpSurf = new G4OpticalSurface("OpSurf");
  OpSurf->SetType(dielectric_dielectric);
  OpSurf->SetFinish(ground);
  OpSurf->SetModel(glisur);
  OpSurf->SetPolish(polish_value); 
  ExternalSurface = new G4LogicalBorderSurface("ExternalSurface",PlasticPhys,physiWorld,OpSurf);

 }



 if(WithGlue){

  //plastic-glue surface -------------------------------------------------------------------

  G4OpticalSurface* OpGlueSurf = new G4OpticalSurface("OpGlueSurf");
  OpGlueSurf->SetType(dielectric_dielectric);
  OpGlueSurf->SetFinish(ground);
  OpGlueSurf->SetModel(glisur);
  OpGlueSurf->SetPolish(polish_value); 

  G4MaterialPropertiesTable* GlueProp = new G4MaterialPropertiesTable(); 
  G4double transmittanceGlue[numEntries] = {0.90,0.90,0.90,0.90,0.90,0.90,0.90,0.90};
  //http://www.4spepro.org/view.php?source=003410-2010-12-13
  //http://www.4spepro.org/previewimage.php?imgx=images%2F003410%2F003410_08_fig2.jpg Fig2
  GlueProp->AddProperty("TRANSMITTANCE",photonEnergies,transmittanceGlue,numEntries);
  OpGlueSurf -> SetMaterialPropertiesTable(GlueProp);

  GlueSurface = new G4LogicalBorderSurface("GlueSurface",PlasticPhys,GluePhys,OpGlueSurf);

  

  }


  // visualization attributes ------------------------------------------------
  //make some colors
  G4Color red(1.0,0.0,0.0),yellow(1.0,1.0,0.0),green(0.0,1.0,0.0), blue(0.0,0.0,1.0),brown(0.4,0.4,0.1),white(1.0,1.0,1.0);  
  logicWorld->SetVisAttributes(G4VisAttributes::Invisible);//set the world volume to invisible borders
  PlasticLogic->SetVisAttributes(new G4VisAttributes(red)); //scintillator is painted in red 
  DetectorLogic->SetVisAttributes(new G4VisAttributes(yellow));//detector is painted in yellow
  if(SelectCoating!=0) CoatingLogic->SetVisAttributes(new G4VisAttributes(green));//recubrimiento is painted in green
  if(WithGlue) GlueLogic->SetVisAttributes(new G4VisAttributes(brown));//glue is painted in brown
  
 
  //always return the physical World
  return physiWorld;
	
} //END G4VPhysicalVolume* DetectorConstruction::Construct()



