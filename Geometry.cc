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
//
/// \file B4/B4c/src/DetectorConstruction.cc
/// \brief Implementation of the B4c::DetectorConstruction class

#include "DetectorConstruction.hh"
#include "CalorimeterSD.hh"
#include "G4Material.hh"
#include "G4NistManager.hh"
#include "G4RotationMatrix.hh"

#include "G4Box.hh" 
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4PVReplica.hh"
#include "G4GlobalMagFieldMessenger.hh"
#include "G4AutoDelete.hh"

#include "G4SDManager.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

namespace B4c
{

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ThreadLocal
G4GlobalMagFieldMessenger* DetectorConstruction::fMagFieldMessenger = nullptr;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  // Define materials
  DefineMaterials();
  // Define volumes
  return DefineVolumes();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void DetectorConstruction::DefineMaterials()
{
  // Lead material defined using NIST Manager
  auto nistManager = G4NistManager::Instance();
  nistManager->FindOrBuildMaterial("G4_Pb");

  // Air
  G4double z,a,fractionmass,density;  // mass of a mole z=mean number of protons;Scintillator
  G4String name, symbol;
  G4int ncomponents;
  a = 14.01*g/mole;
  G4Element*e1N = new G4Element(name="Nitrogen",symbol="N",z=7.,a);
  a =16.00*g/mole;
  G4Element*e1O = new G4Element(name="Oxygen",symbol="O",z=8.,a);
  density = 1.290*mg/cm3;
  G4Material*Air = new G4Material(name="Air", density, ncomponents=2);
  Air -> AddElement(e1N, fractionmass=70*perCent);
  Air -> AddElement(e1O, fractionmass=30*perCent);
        

  // Vacuum
  new G4Material("Galactic", z=1., a=1.01*g/mole,density= universe_mean_density,
                  kStateGas, 2.73*kelvin, 3.e-18*pascal);

  // Print materials
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::DefineVolumes()
{
  // Geometry parameters
//  fNofLayers = 10;
  G4int fNofSlab = 16;
//  G4double absoThickness = 10.*mm;
  G4double SlabLength = 160.*cm;
  G4double SlabWidth = 10.*cm;
  G4double SlabThickness = 3.*cm;
  G4double gapThickness = 100.*cm;
//  G4double gapThickness =  5.*mm;
//  G4double calorSizeXY  = 10.*cm;

//  auto layerThickness = absoThickness + gapThickness;
  auto TotalThickness = 2 * SlabThickness + gapThickness;
  auto TotalWidth = fNofSlab * SlabWidth;
//  auto calorThickness = fNofLayers * layerThickness;
//  auto worldSizeXY = 1.2 * calorSizeXY;
  auto worldSizeX = 1.2 * SlabLength;
  auto worldSizeY = 1.2 * SlabWidth;
  auto worldSizeZ = 1.2 * TotalThickness;
  
//  auto worldSizeZ  = 1.2 * calorThickness;

  // Get materials
  auto defaultMaterial = G4Material::GetMaterial("Galactic");
  //auto absorberMaterial = G4Material::GetMaterial("G4_ANTHRACENE");
  auto absorberMaterial = G4Material::GetMaterial("Galactic");
  auto gapMaterial = G4Material::GetMaterial("Air");

  if ( ! defaultMaterial || ! absorberMaterial || ! gapMaterial ) {
    G4ExceptionDescription msg;
    msg << "Cannot retrieve materials already defined.";
    G4Exception("DetectorConstruction::DefineVolumes()",
      "MyCode0001", FatalException, msg);
  }

  //
  // World
  //
  auto worldS
    = new G4Box("World",           // its name
//                 worldSizeXY/2, worldSizeXY/2, worldSizeZ/2); // its size
                 worldSizeX/2, worldSizeY/2, worldSizeZ/2);

  auto worldLV
    = new G4LogicalVolume(
                 worldS,           // its solid
                 defaultMaterial,  // its material
                 "World");         // its name

  auto worldPV = new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),                         // at (0,0,0)
    worldLV,                                 // its logical volume
    "World",                                 // its name
    nullptr,                                 // its mother  volume
    false,                                   // no boolean operation
    0,                                       // copy number
    fCheckOverlaps);                         // checking overlaps

  //
  // Calorimeter
  //
/*
  auto calorimeterS
    = new G4Box("Calorimeter",     // its name
//                 calorSizeXY/2, calorSizeXY/2, calorThickness/2); // its size
                 SlabLength/2, SlabWidth/2, TotalThickness/2);

  auto calorLV
    = new G4LogicalVolume(
                 calorimeterS,     // its solid
                 defaultMaterial,  // its material
                 "Calorimeter");   //                  calorSizeXY/2, calorSizeXY/2, gapThickness/2); // its sizeits name

  new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),          // at (0,0,0)
    calorLV,                  // its logical volume
    "Calorimeter",            // its name
    worldLV,                  // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    fCheckOverlaps);          // checking overlaps
*/
  //
  // Layer
  //
//  auto layerS
//    = new G4Box("Layer",           // its name
//                 calorSizeXY/2, calorSizeXY/2, layerThickness/2); //its size

  //
  //Slab
  //

  auto Slab
    = new G4Box("Slab",
                SlabLength/2, SlabWidth/2, SlabThickness/2);

//  auto SlabLV
    auto SlabLV
    = new G4LogicalVolume(
                 Slab,
//                 SLabS,           // its solid
                 absorberMaterial,  // its material
                 "Slab");         // its name

  //
  //Panel1
  //
  auto Panel1
    = new G4Box("Panel1",
                SlabLength/2, TotalWidth/2, SlabThickness/2);


    auto Panel1LV
    = new G4LogicalVolume(
                 Panel1,
//                 layerS,           // its solid
                 defaultMaterial,  // its material
                 "Panel1");         // its name

    
  new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(),          // at (0,0,0)
    Panel1LV,                  // its logical volume
    "Panel1",            // its name
    worldLV,                  // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    fCheckOverlaps); 
  

  new G4PVReplica(
                 "Panel1",          // its name
                 SlabLV,          // its logical volume
                 Panel1LV,          // its mother
                 kYAxis,           // axis of replication
                 fNofSlab,        // number of replica
                 SlabWidth);  // witdth of replica

/// Panel2
  auto Panel2
    = new G4Box("Panel2",
                SlabLength/2, TotalWidth/2, SlabThickness/2);
    
    auto Panel2LV
    = new G4LogicalVolume(
                 Panel2,
//                 layerS,           // its solid
                 defaultMaterial,  // its material
                 "Panel2");         // its name

  
  G4RotationMatrix* r = new G4RotationMatrix;
  r ->rotateZ(M_PI/2.*rad);  
    
  new G4PVPlacement(r,  // no rotation
    G4ThreeVector(0,0,SlabThickness/2),          // at (0,0,0)
    Panel2LV,                  // its logical volume
    "Panel2",            // its name
    worldLV,                  // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    fCheckOverlaps); 
  

  new G4PVReplica(
                 "Panel2",          // its name
                 SlabLV,          // its logical volume
                 Panel2LV,          // its mother
                 kYAxis,           // axis of replication
                 fNofSlab,        // number of replica
                 SlabWidth);  // witdth of replica

  
  ////Panel3
  auto Panel3
    = new G4Box("Panel3",
                SlabLength/2, TotalWidth/2, SlabThickness/2);


    auto Panel3LV
    = new G4LogicalVolume(
                 Panel3,
//                 layerS,           // its solid
                 defaultMaterial,  // its material
                 "Panel3");         // its name

    
  new G4PVPlacement(nullptr,  // no rotation
    G4ThreeVector(0,0,(SlabThickness+gapThickness)/2),          // at (0,0,0)
    Panel3LV,                  // its logical volume
    "Panel3",            // its name
    worldLV,                  // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    fCheckOverlaps); 
  


  new G4PVReplica(
                 "Panel3",          // its name
                 SlabLV,          // its logical volume
                 Panel3LV,          // its mother
                 kYAxis,           // axis of replication
                 fNofSlab,        // number of replica
                 SlabWidth);  // witdth of replica

/// Panel4
  auto Panel4
    = new G4Box("Panel4",
                SlabLength/2, TotalWidth/2, SlabThickness/2);
    
  auto Panel4LV
    = new G4LogicalVolume(
                 Panel4,
//                 layerS,           // its solid
                 defaultMaterial,  // its material
                 "Panel4");         // its name


    
  new G4PVPlacement(r,  // no rotation
    G4ThreeVector(0,0,(3*SlabThickness+gapThickness)/2),          // at (0,0,0)
    Panel4LV,                  // its logical volume
    "Panel4",            // its name
    worldLV,                  // its mother  volume
    false,                    // no boolean operation
    0,                        // copy number
    fCheckOverlaps); 
  

  new G4PVReplica(
                 "Panel4",          // its name
                 SlabLV,          // its logical volume
                 Panel4LV,          // its mother
                 kYAxis,           // axis of replication
                 fNofSlab,        // number of replica
                 SlabWidth);  // witdth of replica

  
  

  // Gap
  //
/*
  auto gapScintillator
    = new G4Box("Gap",             // its name
                 SlabLength/2, SlabWidth/2, gapThickness/2); // its size

  auto gapScintillatorLV
    = new G4LogicalVolume(
                 gapScintillator,             // its solid
                 gapMaterial,      // its material
                 "GapLV");         // its name

  new G4PVPlacement(nullptr,                   // no rotation
    G4ThreeVector(0., 0.,  SlabThickness / 2),  // its position
    gapLV,                                     // its logical volume
    "Gap",                                     // its name
    layerLV,                                   // its mother  volume
    false,                                     // no boolean operation
    0,                                         // copy number
    fCheckOverlaps);                           // checking overlaps

  */
  //
  //Scintillator 2
  //

  //
  // Absorber
  //
/*
  auto absorberS
    = new G4Box("Abso",            // its name
                 calorSizeXY/2, calorSizeXY/2, Thickness/2); //its size


  auto absorberLV
    = new G4LogicalVolume(
                 absorberS,        // its solid
                 absorberMaterial, // its material
                 "AbsoLV");        // its name

  new G4PVPlacement(nullptr,                   // no rotationfile:///home/form10/Documents/geant4-master/examples/basic/B4/B4c
    G4ThreeVector(0., 0., -gapThickness / 2),  // its position
    absorberLV,                                // its logical volume
    "Abso",                                    // its name
    layerLV,                                   // its mother  volume
    false,                                     // no boolean operation
    0,                                         // copy number
    fCheckOverlaps);                           // checking overlaps
*/
  //

  //
  // print parameters
  //
/*
  G4cout
    << G4endl
    << "------------------------------------------------------------" << G4endl
    << "---> The calorimeter is " << fNofLayers << " layers of: [ "
    << absoThickness/mm << "mm of " << absorberMaterial->GetName()
    << " + "
    << gapThickness/mm << "mm of " << gapMaterial->GetName() << " ] " << G4endl
    << "------------------------------------------------------------" << G4endl;
*/
  //
  // Visualization attributes
  //
  worldLV->SetVisAttributes (G4VisAttributes::GetInvisible());

  auto simpleBoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  simpleBoxVisAtt->SetVisibility(true);
  SlabLV->SetVisAttributes(simpleBoxVisAtt);

  //
  // Always return the physical World
  //
  return worldPV;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
/*
void DetectorConstruction::ConstructSDandField()
{
  // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  //
  // Sensitive detectors
  //
  auto absoSD
    = new CalorimeterSD("AbsorberSD", "AbsorberHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(absoSD);
  SetSensitiveDetector("AbsoLV",absoSD);

  auto gapSD
    = new CalorimeterSD("GapSD", "GapHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(gapSD);
  SetSensitiveDetector("GapLV",gapSD);

  //
  // Magnetic field
  //
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}
*/

void DetectorConstruction::ConstructSDandField()
{
  // G4SDManager::GetSDMpointer()->SetVerboseLevel(1);

  //
  // Sensitive detectors
  //
/*
  auto layerScintillatorSD
    = new CalorimeterSD("AbsorberSD", "AbsorberHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(absoSD);
  SetSensitiveDetector("AbsoLV",absoSD);

  auto gapSD
    = new CalorimeterSD("GapSD", "GapHitsCollection", fNofLayers);
  G4SDManager::GetSDMpointer()->AddNewDetector(gapSD);
  SetSensitiveDetector("GapLV",gapSD);
*/

  //
  // Magnetic field
  //
  // Create global magnetic field messenger.
  // Uniform magnetic field is then created automatically if
  // the field value is not zero.
  G4ThreeVector fieldValue;
  fMagFieldMessenger = new G4GlobalMagFieldMessenger(fieldValue);
  fMagFieldMessenger->SetVerboseLevel(1);

  // Register the field messenger for deleting
  G4AutoDelete::Register(fMagFieldMessenger);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

}
