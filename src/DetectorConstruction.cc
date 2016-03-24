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
/// \file electromagnetic/TestEm8/src/DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// $Id: DetectorConstruction.cc 89039 2015-03-18 09:28:24Z gcosmo $
//
/////////////////////////////////////////////////////////////////////////
//
// TestEm8: Gaseous detector
//
// Created: 31.08.2010 V.Ivanchenko ob base of V.Grichine code
//
// Modified:
//
////////////////////////////////////////////////////////////////////////
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"
#include "DetectorMessenger.hh"
#include "TargetSD.hh"

#include "G4Material.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"

#include "G4SDManager.hh"
#include "G4GeometryManager.hh"
#include "G4RunManager.hh"
#include "G4NistManager.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"
#include "G4ProductionCuts.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4UnitsTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"
#include "TestParameters.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//http://www.slac.stanford.edu/pubs/icfa/summer98/paper3/paper3.pdf
DetectorConstruction::DetectorConstruction()
  : G4VUserDetectorConstruction(),
    fGasMat(0), fWindowMat(0), fWorldMaterial(0), fNoLayers(1),
    fPhysWorld(0), fLogicWorld(0), fLogicWind(0), fLogicDet(0),
    fDetectorMessenger(0), fGasDetectorCuts(0), fRegGasDet(0)
{
  fGasThickness = 42.0*mm;
  fGasRadius    = 10.*cm;

  fWindowThick  = 8.0*micrometer;
  fAbsThick = 0;

  DefineMaterials();

  fDetectorMessenger = new DetectorMessenger(this);

  G4double cut = 10*mm; //original 23 mm
  G4double solidCut = 1*um;
  fGasDetectorCuts   = new G4ProductionCuts();
  fGasDetectorCuts->SetProductionCut(cut,"gamma");
  fGasDetectorCuts->SetProductionCut(cut*10,"e-");
  fGasDetectorCuts->SetProductionCut(cut,"e+");
  fGasDetectorCuts->SetProductionCut(cut,"proton");

  fSolidDetectorCuts   = new G4ProductionCuts();
  fSolidDetectorCuts->SetProductionCut(solidCut,"gamma");
  fSolidDetectorCuts->SetProductionCut(cut,"e-");
  fSolidDetectorCuts->SetProductionCut(cut,"e+");
  fSolidDetectorCuts->SetProductionCut(cut,"proton");

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{
  delete fDetectorMessenger;
  //delete fGasDetectorCuts;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* water();

void DetectorConstruction::DefineMaterials()
{
  //This function illustrates the possible ways to define materials
  G4String name, symbol ;
  G4double density;
  G4int nel;
  G4int ncomponents;
  G4double fractionmass;

  G4NistManager* manager = G4NistManager::Instance();
  //
  // define Elements
  //

  //Solids
  manager->FindOrBuildElement(26); //iron
  manager->FindOrBuildElement(13); //aluminium

  //Gases
  G4Element* elH  = manager->FindOrBuildElement(1);
  G4Element* elC  = manager->FindOrBuildElement(6);
  G4Element* elO  = manager->FindOrBuildElement(8);
  G4Element* elF  = manager->FindOrBuildElement(9);
  G4Element* elNe = manager->FindOrBuildElement(10);
  G4Element* elXe = manager->FindOrBuildElement(54);
  //
  // simple gases at STP conditions
  //
  G4Material* Argon = manager->FindOrBuildMaterial("G4_Ar");
  Argon->GetIonisation()->SetMeanEnergyPerIonPair(26*eV);
  G4Material* Kr = manager->FindOrBuildMaterial("G4_Kr");
  G4Material* Xe     = manager->FindOrBuildMaterial("G4_Xe");
  //
  // gases at STP conditions
  //
  G4Material* CarbonDioxide =
    manager->FindOrBuildMaterial("G4_CARBON_DIOXIDE");
  G4Material* Mylar  = manager->FindOrBuildMaterial("G4_MYLAR");
  //G4Material* Silicon = manager->FindOrBuildMaterial("G4_SILICON_DIOXIDE");
  G4Material* Methane= manager->FindOrBuildMaterial("G4_METHANE");
  G4Material* Propane= manager->FindOrBuildMaterial("G4_PROPANE");
  G4Material* empty  = manager->FindOrBuildMaterial("G4_Galactic");

  fWater = water();

  //Water->GetMaterialPropertiesTable()->DumpTable();

  // 93% Kr + 7% CH4, STP
  density = 3.491*mg/cm3 ;
  G4Material* Kr7CH4 = new G4Material(name="Kr7CH4"  , density,
                                      ncomponents=2);
  Kr7CH4->AddMaterial( Kr,       fractionmass = 0.986 ) ;
  Kr7CH4->AddMaterial( Methane,  fractionmass = 0.014 ) ;

  G4double TRT_Xe_density = 5.485*mg/cm3;
  G4Material* TRT_Xe = new G4Material(name="TRT_Xe", TRT_Xe_density, nel=1,
                                      kStateGas,293.15*kelvin,1.*atmosphere);
  TRT_Xe->AddElement(elXe,1);
  TRT_Xe->GetIonisation()->SetMeanEnergyPerIonPair(22*eV);

  G4double TRT_CO2_density = 1.842*mg/cm3;
  G4Material* TRT_CO2 = new G4Material(name="TRT_CO2", TRT_CO2_density, nel=2,
                                       kStateGas,293.15*kelvin,1.*atmosphere);
  TRT_CO2->AddElement(elC,1);
  TRT_CO2->AddElement(elO,2);
  TRT_CO2->GetIonisation()->SetMeanEnergyPerIonPair(33*eV);

  G4double TRT_CF4_density = 3.9*mg/cm3;
  G4Material* TRT_CF4 = new G4Material(name="TRT_CF4", TRT_CF4_density, nel=2,
                                       kStateGas,293.15*kelvin,1.*atmosphere);
  TRT_CF4->AddElement(elC,1);
  TRT_CF4->AddElement(elF,4);
  TRT_CF4->GetIonisation()->SetMeanEnergyPerIonPair(54*eV);

  // ATLAS TRT straw tube gas mixture (20 C, 1 atm)
  G4double XeCO2CF4_density = 4.76*mg/cm3;
  G4Material* XeCO2CF4 =
    new G4Material(name="XeCO2CF4", XeCO2CF4_density,
                   ncomponents=3,
                   kStateGas,293.15*kelvin,1.*atmosphere);
  XeCO2CF4->AddMaterial(TRT_Xe,0.807);
  XeCO2CF4->AddMaterial(TRT_CO2,0.039);
  XeCO2CF4->AddMaterial(TRT_CF4,0.154);
  XeCO2CF4->GetIonisation()->SetMeanEnergyPerIonPair(27.357*eV);

  // C3H8,20 C, 2 atm
  density = 3.758*mg/cm3 ;
  G4Material* C3H8 = new G4Material(name="C3H8",density,nel=2) ;
  C3H8->AddElement(elC,3) ;
  C3H8->AddElement(elH,8) ;

  // 87.5% Xe + 7.5% CH4 + 5% C3H8, 20 C, 1 atm
  density = 4.9196*mg/cm3 ;
  G4Material* XeCH4C3H8 = new G4Material(name="XeCH4C3H8"  ,
                                  density,  ncomponents=3);
  XeCH4C3H8->AddMaterial( Xe,      fractionmass = 0.971 ) ;
  XeCH4C3H8->AddMaterial( Methane, fractionmass = 0.010 ) ;
  XeCH4C3H8->AddMaterial( Propane, fractionmass = 0.019 ) ;

  // 93% Ar + 7% CH4, STP
  density = 1.709*mg/cm3 ;
  G4Material* Ar7CH4 = new G4Material(name="Ar7CH4", density, ncomponents=2);
  Ar7CH4->AddMaterial( Argon,    fractionmass = 0.971 ) ;
  Ar7CH4->AddMaterial( Methane,  fractionmass = 0.029 ) ;

  // 80% Ar + 20% CO2, STP
  density = 1.8223*mg/cm3 ;
  G4Material* Ar_80CO2_20 = new G4Material(name="ArCO2"  , density,
                                           ncomponents=2);
  Ar_80CO2_20->AddMaterial( Argon,           fractionmass = 0.783 ) ;
  Ar_80CO2_20->AddMaterial( CarbonDioxide,   fractionmass = 0.217 ) ;

  // 80% Xe + 20% CO2, STP
  density = 5.0818*mg/cm3 ;
  G4Material* Xe20CO2 = new G4Material(name="Xe20CO2", density,
                                       ncomponents=2);
  Xe20CO2->AddMaterial( Xe,            fractionmass = 0.922 ) ;
  Xe20CO2->AddMaterial( CarbonDioxide, fractionmass = 0.078 ) ;

  // 80% Kr + 20% CO2, STP
  density = 3.601*mg/cm3 ;
  G4Material* Kr20CO2 = new G4Material(name="Kr20CO2"  , density,
                                       ncomponents=2);
  Kr20CO2->AddMaterial( Kr,            fractionmass = 0.89 ) ;
  Kr20CO2->AddMaterial( CarbonDioxide, fractionmass = 0.11 ) ;

  // ALICE mixture TPC_Ne-CO2-2
  density = 0.939*mg/cm3 ;
  G4Material* NeCO2 = new G4Material(name="TPC_Ne-CO2-2", density,
                                            ncomponents=3);
  NeCO2->AddElement( elNe, fractionmass = 0.8039 ) ;
  NeCO2->AddElement( elO,  fractionmass = 0.1426 ) ;
  NeCO2->AddElement( elC,  fractionmass = 0.0535 ) ;

  fGasMat = TRT_CF4;
  fWindowMat = Mylar;

/*  G4double gasDensity = fGasMat->GetDensity();
  G4double windowDensity = fWindowMat->GetDensity();
  G4double ratio = (2*fWindowThick)/(2*fWindowThick + fGasThickness);
  G4double totDensity = gasDensity*(1.0-ratio) + windowDensity*(ratio);
  G4double fractionWindow = ratio*windowDensity/totDensity;
  G4double fractionGas = (1.0-ratio)*gasDensity/totDensity;
  G4Material* radiatorMat = new G4Material("radiatorMat", totDensity, ncomponents=2);
  radiatorMat->AddMaterial(fWindowMat, fractionWindow);
  radiatorMat->AddMaterial(fGasMat, fractionGas);

  fRadiatorMat = radiatorMat;*/
  fWorldMaterial = empty;

  //G4cout << *(G4Material::GetMaterialTable()) << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4double contThick = (fWindowThick*2 + fGasThickness);
  G4double contR     = fWindowThick*2 + fGasRadius;

  G4double strawsThick = contThick*fNoLayers;

  G4double worldSizeZ = (strawsThick+fAbsThick)*2;
  G4double worldSizeR = contR*1.2;

  TestParameters::GetPointer()->SetPositionZ(-0.55*contThick);

  // Printout parameters
  G4cout << "\n The  WORLD   is made of "
         << worldSizeZ/mm << "mm of " << fWorldMaterial->GetName() ;
  G4cout << ", the transverse size (R) of the world is " << worldSizeR/mm
         << " mm. " << G4endl;
  G4cout << " The CONTAINER is made of "
         << fWindowThick/mm << "mm of " << fWindowMat->GetName() << G4endl;
  G4cout << " The TARGET is made of "
         << fGasThickness/mm << "mm of " << fGasMat->GetName() ;
  G4cout << ", the transverse size (R) is " << fGasRadius/mm << " mm. "
         << G4endl;
  G4cout << G4endl;

  // World
  G4Tubs* SolidWorld =
    new G4Tubs("World",0.,worldSizeR,worldSizeZ/2.,0.,CLHEP::twopi);

  fLogicWorld = new G4LogicalVolume(SolidWorld, fWorldMaterial, "World");

  fPhysWorld = new G4PVPlacement(0,
                                   G4ThreeVector(0.,0.,0.),
                                 "World",
                                 fLogicWorld,
                                 0,
                                 false,
                                 0);

  //Front Piece
 if(fAbsThick > 0){
  G4Tubs* buffer = new G4Tubs("Buffer",
          	  0.,contR,fAbsThick/2.,0.,CLHEP::twopi);

  G4LogicalVolume* fLogicBuff = new G4LogicalVolume(buffer, fWater, "FrontPiece");
  /*G4PVPlacement* PhysBuff =*/ new G4PVPlacement(0, G4ThreeVector(0.,0.,fAbsThick/2),
                                              "FrontPiece",  fLogicBuff,
                                              fPhysWorld, false, 0);
 }

 //
 // Straws
 //
 G4VSolid* straws
   = new G4Tubs("Straws",     // its name
               0., contR, strawsThick/2.,0.,CLHEP::twopi); // its size

 //G4Material* usedMat = fNoLayers == 1? fWorldMaterial : fRadiatorMat;
  G4LogicalVolume* strawLV = new G4LogicalVolume(
                straws,     // its solid
                fWorldMaterial,  // its material
                "Straws");   // its name

  new G4PVPlacement(
                 0,                // no rotation
                 G4ThreeVector(0.,0.,fAbsThick+strawsThick/2),
                 strawLV,          // its logical volume
                 "Straws",    // its name
                 fLogicWorld,          // its mother  volume
                 false,            // no boolean operation
                 0);  // checking overlaps

  // Window
  G4Tubs* wind = new G4Tubs("Absorber",
		  0.,contR,contThick/2.,
		  0.,CLHEP::twopi);

  fLogicWind = new G4LogicalVolume(wind, fWindowMat, "Window");
/*  G4PVPlacement* PhysWind2 = new G4PVPlacement(0, G4ThreeVector(0.,0.,fAbsThick+contThick/2),
		  "Window",  fLogicWind,
		  PhysWind, false, 0);*/

  G4PVReplica* PhysWind = new G4PVReplica(
                 "Absorber",          // its name
                 fLogicWind,          // its logical volume
                 strawLV,          // its mother
                 kZAxis,           // axis of replication
                 fNoLayers,        // number of replica
                 contThick);  // witdth of replica

  // Detector volume
  G4Tubs* det = new G4Tubs("Gas", 0., fGasRadius, fGasThickness/2.,
		  0., CLHEP::twopi);

  fLogicDet = new G4LogicalVolume(det, fGasMat, "Gas");
  new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), "Gas", fLogicDet,
		  PhysWind, false, 0);

  // defined gas detector region
  fRegGasDet = new G4Region("GasDetector");
  fRegGasDet->SetProductionCuts(fGasDetectorCuts);
  fRegGasDet->AddRootLogicalVolume(fLogicDet);

  fRegWindow = new G4Region("WindowArea");
  fRegWindow->SetProductionCuts(fSolidDetectorCuts);
  fRegWindow->AddRootLogicalVolume(fLogicWind);

  // visualisation
  fLogicWorld->SetVisAttributes(G4VisAttributes::Invisible);
  G4VisAttributes* color1 = new G4VisAttributes(G4Colour(0.3, 0.3, 0.3));
  fLogicWind->SetVisAttributes(color1);
  G4VisAttributes* color2 = new G4VisAttributes(G4Colour(0.0, 0.3, 0.7));
  fLogicDet->SetVisAttributes(color2);

  if(0.0 == fGasMat->GetIonisation()->GetMeanEnergyPerIonPair()) {
    SetPairEnergy(20*eV);
  }
  return fPhysWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::ConstructSDandField()
{
  SetSensitiveDetector(fLogicDet, new TargetSD("GasSD"));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGasMaterial(const G4String& name)
{
  // get the pointer to the existing material
  G4Material* mat = G4Material::GetMaterial(name, false);

  // create the material by its name
  if(!mat) { mat = G4NistManager::Instance()->FindOrBuildMaterial(name); }

  if (mat && mat != fGasMat) {
    G4cout << "### New target material: " << mat->GetName() << G4endl;
    fGasMat = mat;
    if(fLogicDet) {
      fLogicDet->SetMaterial(mat);
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetContainerMaterial(const G4String& name)
{
  // get the pointer to the existing material
  G4Material* mat = G4Material::GetMaterial(name, false);

  // create the material by its name
  if(!mat) { mat = G4NistManager::Instance()->FindOrBuildMaterial(name); }

  if (mat && mat != fWindowMat) {
    G4cout << "### New material for container: " << mat->GetName() << G4endl;
    fWindowMat = mat;
    if(fLogicWind) {
      fLogicWind->SetMaterial(mat);
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetWorldMaterial(const G4String& name)
{
  // get the pointer to the existing material
  G4Material* mat = G4Material::GetMaterial(name, false);

  // create the material by its name
  if(!mat) { mat = G4NistManager::Instance()->FindOrBuildMaterial(name); }

  if (mat && mat != fWorldMaterial) {
    G4cout << "### New World material: " << mat->GetName() << G4endl;
    fWorldMaterial = mat;
    if(fLogicWorld) {
      fLogicWorld->SetMaterial(mat);
      G4RunManager::GetRunManager()->PhysicsHasBeenModified();
    }
  }
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetContainerMaterial(){
	return fWindowMat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Material* DetectorConstruction::GetGasMaterial(){
	return fGasMat;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetContainerThickness(){
	return fWindowThick;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetGasThickness(){
	return fGasThickness;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int DetectorConstruction::GetNumberOfLayers(){
	return fNoLayers;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4LogicalVolume* DetectorConstruction::GetLogicalRadiator(){
	return fLogicWind;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetAbsorberThickness(){
	return fAbsThick;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4double DetectorConstruction::GetGasRadius(){
	return fGasRadius;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGasThickness(G4double val)
{
  if(fPhysWorld) {
    G4Exception ("DetectorConstruction::SetGasThickness()", "test005",
                 JustWarning,
                 "Attempt to change already constructed geometry is ignored");
  } else {
    fGasThickness = val;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetGasRadius(G4double val)
{
  if(fPhysWorld) {
    G4Exception ("DetectorConstruction::SetGasRadius()", "test005",
                 JustWarning,
                 "Attempt to change already constructed geometry is ignored");
  } else {
    fGasRadius = val;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetNumberOfLayers(G4int val)
{
  if(fPhysWorld) {
    G4Exception ("DetectorConstruction::SetNumberOfLayers()", "test005",
                 JustWarning,
                 "Attempt to change already constructed geometry is ignored");
  } else {
    fNoLayers = val;
  }
}

void DetectorConstruction::SetAbsorberThickness(G4double val)
{
  if(fPhysWorld) {
    G4Exception ("DetectorConstruction::SetAbsorberThickness()", "test005",
                 JustWarning,
                 "Attempt to change already constructed geometry is ignored");
  } else {
    fAbsThick = val;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetContainerThickness(G4double val)
{
  if(fPhysWorld) {
    G4Exception ("DetectorConstruction::SetContainerThickness()", "test005",
                 JustWarning,
                 "Attempt to change already constructed geometry is ignored");
  } else {
    fWindowThick = val;
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::SetPairEnergy(G4double val)
{
  if(val > 0.0) {
    fGasMat->GetIonisation()->SetMeanEnergyPerIonPair(val);
  }
}

G4Material* water(){
	  G4double a, z, density;
	  G4int nelements;

	  G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.01*g/mole);
	  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);

	  G4Material* water = new G4Material("Water", density= 1.0*g/cm3, nelements=2);
	  water->AddElement(H, 2);
	  water->AddElement(O, 1);

	//
	// ------------ Generate & Add Material Properties Table ------------
	//
	  G4double photonEnergy[] =
	            { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
	              2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
	              2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
	              2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
	              2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
	              3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
	              3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
	              3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };

	  const G4int nEntries = sizeof(photonEnergy)/sizeof(G4double);

	//
	// Water
	//
	  G4double refractiveIndex1[] =
	            { 1.3435, 1.344,  1.3445, 1.345,  1.3455,
	              1.346,  1.3465, 1.347,  1.3475, 1.348,
	              1.3485, 1.3492, 1.35,   1.3505, 1.351,
	              1.3518, 1.3522, 1.3530, 1.3535, 1.354,
	              1.3545, 1.355,  1.3555, 1.356,  1.3568,
	              1.3572, 1.358,  1.3585, 1.359,  1.3595,
	              1.36,   1.3608};

	  assert(sizeof(refractiveIndex1) == sizeof(photonEnergy));

	  G4double absorption[] =
	           {3.448*m,  4.082*m,  6.329*m,  9.174*m, 12.346*m, 13.889*m,
	           15.152*m, 17.241*m, 18.868*m, 20.000*m, 26.316*m, 35.714*m,
	           45.455*m, 47.619*m, 52.632*m, 52.632*m, 55.556*m, 52.632*m,
	           52.632*m, 47.619*m, 45.455*m, 41.667*m, 37.037*m, 33.333*m,
	           30.000*m, 28.500*m, 27.000*m, 24.500*m, 22.000*m, 19.500*m,
	           17.500*m, 14.500*m };

	  assert(sizeof(absorption) == sizeof(photonEnergy));

	  G4double scintilFast[] =
	            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
	              1.00, 1.00, 1.00, 1.00 };

	  assert(sizeof(scintilFast) == sizeof(photonEnergy));

	  G4double scintilSlow[] =
	            { 0.01, 1.00, 2.00, 3.00, 4.00, 5.00, 6.00,
	              7.00, 8.00, 9.00, 8.00, 7.00, 6.00, 4.00,
	              3.00, 2.00, 1.00, 0.01, 1.00, 2.00, 3.00,
	              4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 8.00,
	              7.00, 6.00, 5.00, 4.00 };

	  assert(sizeof(scintilSlow) == sizeof(photonEnergy));

	  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();

	  myMPT1->AddProperty("RINDEX",       photonEnergy, refractiveIndex1,nEntries)
	        ->SetSpline(true);
	  myMPT1->AddProperty("ABSLENGTH",    photonEnergy, absorption,     nEntries)
	        ->SetSpline(true);
	  myMPT1->AddProperty("FASTCOMPONENT",photonEnergy, scintilFast,     nEntries)
	        ->SetSpline(true);
	  myMPT1->AddProperty("SLOWCOMPONENT",photonEnergy, scintilSlow,     nEntries)
	        ->SetSpline(true);

	  myMPT1->AddConstProperty("SCINTILLATIONYIELD",50./MeV);
	  myMPT1->AddConstProperty("RESOLUTIONSCALE",1.0);
	  myMPT1->AddConstProperty("FASTTIMECONSTANT", 1.*ns);
	  myMPT1->AddConstProperty("SLOWTIMECONSTANT",10.*ns);
	  myMPT1->AddConstProperty("YIELDRATIO",0.8);

	  G4double energy_water[] = {
	     1.56962*eV, 1.58974*eV, 1.61039*eV, 1.63157*eV,
	     1.65333*eV, 1.67567*eV, 1.69863*eV, 1.72222*eV,
	     1.74647*eV, 1.77142*eV, 1.7971 *eV, 1.82352*eV,
	     1.85074*eV, 1.87878*eV, 1.90769*eV, 1.93749*eV,
	     1.96825*eV, 1.99999*eV, 2.03278*eV, 2.06666*eV,
	     2.10169*eV, 2.13793*eV, 2.17543*eV, 2.21428*eV,
	     2.25454*eV, 2.29629*eV, 2.33962*eV, 2.38461*eV,
	     2.43137*eV, 2.47999*eV, 2.53061*eV, 2.58333*eV,
	     2.63829*eV, 2.69565*eV, 2.75555*eV, 2.81817*eV,
	     2.88371*eV, 2.95237*eV, 3.02438*eV, 3.09999*eV,
	     3.17948*eV, 3.26315*eV, 3.35134*eV, 3.44444*eV,
	     3.54285*eV, 3.64705*eV, 3.75757*eV, 3.87499*eV,
	     3.99999*eV, 4.13332*eV, 4.27585*eV, 4.42856*eV,
	     4.59258*eV, 4.76922*eV, 4.95999*eV, 5.16665*eV,
	     5.39129*eV, 5.63635*eV, 5.90475*eV, 6.19998*eV
	  };

	  const G4int numentries_water = sizeof(energy_water)/sizeof(G4double);

	  //assume 100 times larger than the rayleigh scattering for now.
	  G4double mie_water[] = {
	     167024.4*m, 158726.7*m, 150742  *m,
	     143062.5*m, 135680.2*m, 128587.4*m,
	     121776.3*m, 115239.5*m, 108969.5*m,
	     102958.8*m, 97200.35*m, 91686.86*m,
	     86411.33*m, 81366.79*m, 76546.42*m,
	     71943.46*m, 67551.29*m, 63363.36*m,
	     59373.25*m, 55574.61*m, 51961.24*m,
	     48527.00*m, 45265.87*m, 42171.94*m,
	     39239.39*m, 36462.50*m, 33835.68*m,
	     31353.41*m, 29010.30*m, 26801.03*m,
	     24720.42*m, 22763.36*m, 20924.88*m,
	     19200.07*m, 17584.16*m, 16072.45*m,
	     14660.38*m, 13343.46*m, 12117.33*m,
	     10977.70*m, 9920.416*m, 8941.407*m,
	     8036.711*m, 7202.470*m, 6434.927*m,
	     5730.429*m, 5085.425*m, 4496.467*m,
	     3960.210*m, 3473.413*m, 3032.937*m,
	     2635.746*m, 2278.907*m, 1959.588*m,
	     1675.064*m, 1422.710*m, 1200.004*m,
	     1004.528*m, 833.9666*m, 686.1063*m
	  };

	  assert(sizeof(mie_water) == sizeof(energy_water));

	  // gforward, gbackward, forward backward ratio
	  G4double mie_water_const[3]={0.99,0.99,0.8};

	  myMPT1->AddProperty("MIEHG",energy_water,mie_water,numentries_water)
	        ->SetSpline(true);
	  myMPT1->AddConstProperty("MIEHG_FORWARD",mie_water_const[0]);
	  myMPT1->AddConstProperty("MIEHG_BACKWARD",mie_water_const[1]);
	  myMPT1->AddConstProperty("MIEHG_FORWARD_RATIO",mie_water_const[2]);

	  G4cout << "Water G4MaterialPropertiesTable" << G4endl;
	  //myMPT1->DumpTable();

	  water->SetMaterialPropertiesTable(myMPT1);

	  // Set the Birks Constant for the Water scintillator

	  water->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
	  return water;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
