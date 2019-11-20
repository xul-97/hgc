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
/// \file B1DetectorConstruction.cc
/// \brief Implementation of the B1DetectorConstruction class

#include "B1DetectorConstruction.hh"

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Torus.hh"
#include "G4Tubs.hh"
#include "G4Hype.hh"
#include "G4Para.hh"
#include "G4EllipticalTube.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4IntersectionSolid.hh"
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::B1DetectorConstruction()
: G4VUserDetectorConstruction(),
  fScoringVolume(0)
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

B1DetectorConstruction::~B1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* B1DetectorConstruction::Construct()
{  
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Envelope parameters
  //
  G4double env_sizeXY = 40*cm, env_sizeZ = 200*cm;
  G4Material* env_mat = nist->FindOrBuildMaterial("G4_WATER");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeXY = 1.2*env_sizeXY;
  G4double world_sizeZ  = 1.2*env_sizeZ;
  G4Material* world_mat = nist->FindOrBuildMaterial("G4_AIR");
  
  G4Box* solidWorld =    
    new G4Box("World",                       //its name
       0.5*world_sizeXY, 0.5*world_sizeXY, 0.5*world_sizeZ);     //its size
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Envelope
  //  
  G4Box* solidEnv =    
    new G4Box("Envelope",                    //its name
        0.5*env_sizeXY, 0.5*env_sizeXY, 0.5*env_sizeZ); //its size
      
  G4LogicalVolume* logicEnv =                         
    new G4LogicalVolume(solidEnv,            //its solid
                        env_mat,             //its material
                        "Envelope");         //its name
               
  new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicEnv,                //its logical volume
                    "Envelope",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  G4RotationMatrix* rmx1 = new G4RotationMatrix();
  rmx1->rotateX(90*deg); 
  G4RotationMatrix* rmx2 = new G4RotationMatrix();
  rmx2->rotateX(-90*deg); 
  G4RotationMatrix* rmy1 = new G4RotationMatrix();
  rmy1->rotateY(90*deg); 
  G4RotationMatrix* rmy2 = new G4RotationMatrix();
  rmy2->rotateY(-90*deg); 

  //
  // letter A
  //

G4Material* mat = nist->FindOrBuildMaterial("G4_BONE_COMPACT_ICRU");
  G4ThreeVector posA = G4ThreeVector(0, 0*cm, -85*cm);

  G4Trd* solidA1 =    
    new G4Trd("A1",             //its name
              2*cm, 2*cm, 
              1*cm, 6*cm, 8.5*cm); //its size

  G4Trd* solidA2 =    
    new G4Trd("A2",             //its name
              2.01*cm, 2.01*cm, 
              0.1*cm, 1.5*cm, 3*cm); //its size  
  G4SubtractionSolid* solidA3 =
   new G4SubtractionSolid("A3", solidA1, solidA2,
                           0, G4ThreeVector(0, 0,-2*cm));

    G4Trd* solidA4 =    
    new G4Trd("A4",                      //its name
              2.01*cm, 2.01*cm, 
              2*cm, 4*cm, 3*cm); //its size  
    G4SubtractionSolid* solidA5 =
   new G4SubtractionSolid("A5", solidA3, solidA4,
                           0, G4ThreeVector(0,0,6*cm));
   G4LogicalVolume* logicA5 =                         
    new G4LogicalVolume(solidA5,         //its solid
                        mat,          //its material
                        "A5");           //its name
    

  new G4PVPlacement(rmx2,                       // rotation
                    posA,                    //at position
                    logicA5,             //its logical volume
                    "A5",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
    fScoringVolume = logicA5;

//
// letter D
  //

  G4ThreeVector posD = G4ThreeVector(0, 0*cm, -70*cm);

 
  G4Box* solidD1 =    
    new G4Box("D1",             //its name
              1.5*cm,  
              8.5*cm,
              2*cm); //its size

  G4Torus* solidD2 =    
    new G4Torus("D2",             //its name
              0.1*cm, 2*cm, 
              6.5*cm, -90.*deg, 180.*deg); //its size  
    
  G4UnionSolid* solidD3 =
   new G4UnionSolid("D3", solidD1, solidD2,
                           0, G4ThreeVector(1.5*cm,0,0));

   G4LogicalVolume* logicD3 =                         
    new G4LogicalVolume(solidD3,         //its solid
                        mat,          //its material
                        "D3");           //its name
   
  new G4PVPlacement(rmy1,                       // rotation
                    posD,                    //at position
                    logicD3,             //its logical volume
                    "D3",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
    fScoringVolume = logicD3;

//
  // letter E
  //

  G4ThreeVector posE1 = G4ThreeVector(0, 0*cm, -55*cm);
  G4ThreeVector posE2 = G4ThreeVector(0, 7.5*cm, -50.25*cm);
  G4ThreeVector posE3 = G4ThreeVector(0, 0*cm, -50.75*cm);
  G4ThreeVector posE4 = G4ThreeVector(0, -7.5*cm, -50.25*cm);

  G4Box* solidE1 =    
    new G4Box("E1",             //its name
              2*cm,  
              8.5*cm,
              1.5*cm); //its size

    G4LogicalVolume* logicE1 =                         
    new G4LogicalVolume(solidE1,         //its solid
                        mat,          //its material
                        "E1");           //its name
    
  new G4PVPlacement(0,                       //no rotation
                    posE1,                    //at position
                    logicE1,             //its logical volume
                    "E1",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  G4Box* solidE2 =    
    new G4Box("E2",             //its name
              2*cm,  
              1*cm,
              4.25*cm); //its size
    G4LogicalVolume* logicE2 =                         
    new G4LogicalVolume(solidE2,         //its solid
                        mat,          //its material
                        "E2");           //its name
    
  new G4PVPlacement(0,                       //no rotation
                    posE2,                    //at position
                    logicE2,             //its logical volume
                    "E2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

    G4Box* solidE3 =    
    new G4Box("E3",             //its name
              2*cm,  
              1*cm,
              3.75*cm); //its size
    G4LogicalVolume* logicE3 =                         
    new G4LogicalVolume(solidE3,         //its solid
                        mat,          //its material
                        "E1");           //its name
    
  new G4PVPlacement(0,                       //no rotation
                    posE3,                    //at position
                    logicE3,             //its logical volume
                    "E3",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

    G4Box* solidE4 =    
    new G4Box("E4",             //its name
              2*cm,  
              1*cm,
              4.25*cm); //its size
   G4LogicalVolume* logicE4 =                         
    new G4LogicalVolume(solidE4,         //its solid
                        mat,          //its material
                        "E4");           //its name
    

  new G4PVPlacement(0,                       //no rotation
                    posE4,                    //at position
                    logicE4,             //its logical volume
                    "E4",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
  fScoringVolume = logicE1;

  //
  // letter F
  //

  G4ThreeVector posF1 = G4ThreeVector(0, 0*cm, -40*cm);
  G4ThreeVector posF2 = G4ThreeVector(0, 7.5*cm, -35.25*cm);
  G4ThreeVector posF3 = G4ThreeVector(0, 1*cm, -35.75*cm);
  

  G4Box* solidF1 =    
    new G4Box("F1",             //its name
              2*cm,  
              8.5*cm,
              1.5*cm); //its size

    G4LogicalVolume* logicF1 =                         
    new G4LogicalVolume(solidF1,         //its solid
                        mat,          //its material
                        "F1");           //its name
    
  new G4PVPlacement(0,                       //no rotation
                    posF1,                    //at position
                    logicF1,             //its logical volume
                    "F1",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  G4Box* solidF2 =    
    new G4Box("F2",             //its name
              2*cm,  
              1*cm,
              4.25*cm); //its size
    G4LogicalVolume* logicF2 =                         
    new G4LogicalVolume(solidF2,         //its solid
                        mat,          //its material
                        "F2");           //its name
    
  new G4PVPlacement(0,                       //no rotation
                    posF2,                    //at position
                    logicF2,             //its logical volume
                    "F2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

    G4Box* solidF3 =    
    new G4Box("F3",             //its name
              2*cm,  
              1*cm,
              3.75*cm); //its size
    G4LogicalVolume* logicF3 =                         
    new G4LogicalVolume(solidF3,         //its solid
                        mat,          //its material
                        "F3");           //its name
    
  new G4PVPlacement(0,                       //no rotation
                    posF3,                    //at position
                    logicF3,             //its logical volume
                    "F3",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

 
  fScoringVolume = logicF1;

//
  // letter H
  //

  G4ThreeVector posH1 = G4ThreeVector(0, 0*cm, -25*cm);
  G4ThreeVector posH2 = G4ThreeVector(0, 0*cm, -16.25*cm);
  G4ThreeVector posH3 = G4ThreeVector(0, 0*cm, -20.25*cm);
  

  G4Box* solidH1 =    
    new G4Box("H1",             //its name
              2*cm,  
              8.5*cm,
              1.5*cm); //its size

    G4LogicalVolume* logicH1 =                         
    new G4LogicalVolume(solidH1,         //its solid
                        mat,          //its material
                        "H1");           //its name
    
  new G4PVPlacement(0,                       //no rotation
                    posH1,                    //at position
                    logicH1,             //its logical volume
                    "H1",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  G4Box* solidH2 =    
    new G4Box("H2",             //its name
              2*cm,  
              8.5*cm,
              1.5*cm); //its size
    G4LogicalVolume* logicH2 =                         
    new G4LogicalVolume(solidH2,         //its solid
                        mat,          //its material
                        "H2");           //its name
    
  new G4PVPlacement(0,                       //no rotation
                    posH2,                    //at position
                    logicH2,             //its logical volume
                    "H2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

    G4Box* solidH3 =    
    new G4Box("H3",             //its name
              2*cm,  
              1*cm,
              4.25*cm); //its size
    G4LogicalVolume* logicH3 =                         
    new G4LogicalVolume(solidH3,         //its solid
                        mat,          //its material
                        "H3");           //its name
    
  new G4PVPlacement(0,                       //no rotation
                    posH3,                    //at position
                    logicH3,             //its logical volume
                    "H3",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

 
  fScoringVolume = logicH1;

  //
   // letter J
   //

  G4ThreeVector posJ1 = G4ThreeVector(0, 3*cm, 0*cm);
  G4ThreeVector posJ2 = G4ThreeVector(0, -3.5*cm, -4*cm);

 
  G4Tubs* solidJ1 =    
    new G4Tubs("J1",             //its name
              0.1*cm,  
              2*cm,
              6.5*cm,
              0,
              360.*deg); //its size
    G4LogicalVolume* logicJ1 =                         
    new G4LogicalVolume(solidJ1,         //its solid
                        mat,          //its material
                        "J1");           //its name
 
new G4PVPlacement(rmx2,                       // no rotation
                    posJ1,                    //at position
                    logicJ1,             //its logical volume
                    "J1",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  G4Torus* solidJ2 =    
    new G4Torus("J2",             //its name
              0.1*cm, 2*cm, 
              4*cm, -180.*deg, 180.*deg); //its size  
    
 
   G4LogicalVolume* logicJ2 =                         
    new G4LogicalVolume(solidJ2,         //its solid
                        mat,          //its material
                        "J2");           //its name
                        
  new G4PVPlacement(rmy2,                       // rotation
                    posJ2,                    //at position
                    logicJ2,             //its logical volume
                    "J2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
    fScoringVolume = logicJ1;

//
  // letter L
  //

  G4ThreeVector posL1 = G4ThreeVector(0, 0.5*cm, 10*cm);
  G4ThreeVector posL2 = G4ThreeVector(0, -7*cm, 14.75*cm);
  

  G4Box* solidL1 =    
    new G4Box("L1",             //its name
              2*cm,  
              8.5*cm,
              1.5*cm); //its size

    G4LogicalVolume* logicL1 =                         
    new G4LogicalVolume(solidL1,         //its solid
                        mat,          //its material
                        "L1");           //its name
    
  new G4PVPlacement(0,                       //no rotation
                    posL1,                    //at position
                    logicL1,             //its logical volume
                    "L1",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  G4Box* solidL2 =    
    new G4Box("L2",             //its name
              2*cm,  
              1*cm,
              4.25*cm); //its size
    G4LogicalVolume* logicL2 =                         
    new G4LogicalVolume(solidL2,         //its solid
                        mat,          //its material
                        "L2");           //its name
    
  new G4PVPlacement(0,                       //no rotation
                    posL2,                    //at position
                    logicL2,             //its logical volume
                    "L2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  fScoringVolume = logicL1;

//
// letter N
  //

  G4ThreeVector posN1 = G4ThreeVector(0, 0*cm, 25.5*cm);
  G4ThreeVector posN2 = G4ThreeVector(0, 0*cm, 34.5*cm);
   G4ThreeVector posN3 = G4ThreeVector(0, 0*cm, 30*cm); 

  G4Box* solidN1 =    
    new G4Box("N1",             //its name
              2*cm,  
              8.5*cm,
              1.5*cm); //its size

    G4LogicalVolume* logicN1 =                         
    new G4LogicalVolume(solidN1,         //its solid
                        mat,          //its material
                        "N1");           //its name
    
  new G4PVPlacement(0,                       //no rotation
                    posN1,                    //at position
                    logicN1,             //its logical volume
                    "N1",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  G4Box* solidN2 =    
    new G4Box("N2",             //its name
              2*cm,  
              8.5*cm,
              1.5*cm); //its size
    G4LogicalVolume* logicN2 =                         
    new G4LogicalVolume(solidN2,         //its solid
                        mat,          //its material
                        "N2");           //its name
    
  new G4PVPlacement(0,                       //no rotation
                    posN2,                    //at position
                    logicN2,             //its logical volume
                    "N2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

  G4Para* solidN3 =
   new G4Para("N3",          //its name 
              1.5*cm,
              8.5*cm,
              2*cm,
              25.*deg, 0,0); // its size
   
 G4LogicalVolume* logicN3 =                         
    new G4LogicalVolume(solidN3,         //its solid
                        mat,          //its material
                        "N3");           //its name


  new G4PVPlacement(rmy2,                       //no rotation
                    posN3,                    //at position
                    logicN3,             //its logical volume
                    "N3",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  fScoringVolume = logicN3;

//
// letter O
  //

  G4ThreeVector posO = G4ThreeVector(0, 0*cm, 45*cm);

 
  G4EllipticalTube* solidO1 =    
    new G4EllipticalTube("O1",             //its name
                         5.5*cm,  
                         9*cm,
                         2*cm); //its size

  G4EllipticalTube* solidO2 =    
    new G4EllipticalTube("O2",             //its name
                         3*cm,
                         5.5*cm, 
                         2.01*cm); //its size  
    
  G4SubtractionSolid* solidO3 =
   new G4SubtractionSolid("O3", solidO1, solidO2,
                          0, G4ThreeVector(0,0,0));

   G4LogicalVolume* logicO3 =                         
    new G4LogicalVolume(solidO3,         //its solid
                        mat,          //its material
                        "O3");           //its name
   
  new G4PVPlacement(rmy1,                       // rotation
                    posO,                    //at position
                    logicO3,             //its logical volume
                    "O3",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
    fScoringVolume = logicO3;

//
// letter S
  //

  G4ThreeVector posS1 = G4ThreeVector(0, 4*cm, 60*cm);
  G4ThreeVector posS2 = G4ThreeVector(0, -4*cm, 60*cm);

G4Torus* solidS1 =    
    new G4Torus("S1",             //its name
              0.1*cm, 1.5*cm, 
              4*cm, -90.*deg, 270.*deg); //its size  
    
 
   G4LogicalVolume* logicS1 =                         
    new G4LogicalVolume(solidS1,         //its solid
                        mat,          //its material
                        "S1");           //its name

   
  new G4PVPlacement(rmy2,                       // rotation
                    posS1,                    //at position
                    logicS1,             //its logical volume
                    "S1",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

G4Torus* solidS2 =    
    new G4Torus("S2",             //its name
              0.1*cm, 1.5*cm, 
              4*cm, 90.*deg, 270.*deg); //its size  
    
 
   G4LogicalVolume* logicS2 =                         
    new G4LogicalVolume(solidS2,         //its solid
                        mat,          //its material
                        "S2");           //its name
   
  new G4PVPlacement(rmy2,                       // rotation
                    posS2,                    //at position
                    logicS2,             //its logical volume
                    "S2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
    fScoringVolume = logicS1;

//
// letter Y
  //

  G4ThreeVector posY1 = G4ThreeVector(0, -4*cm, 75*cm);
  G4ThreeVector posY2 = G4ThreeVector(0, 4*cm, 77.5*cm);
   G4ThreeVector posY3 = G4ThreeVector(0, 4*cm, 72.5*cm); 

  G4Box* solidY1 =    
    new G4Box("Y1",             //its name
              2*cm,  
              4*cm,
              1.5*cm); //its size

    G4LogicalVolume* logicY1 =                         
    new G4LogicalVolume(solidY1,         //its solid
                        mat,          //its material
                        "Y1");           //its name
    
  new G4PVPlacement(0,                       //no rotation
                    posY1,                    //at position
                    logicY1,             //its logical volume
                    "Y1",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking

G4Para* solidY2 =
   new G4Para("Y2",          //its name 
              1.5*cm,
              4*cm,
              2*cm,
              30.*deg, 0,0); // its size
   
 G4LogicalVolume* logicY2 =                         
    new G4LogicalVolume(solidY2,         //its solid
                        mat,          //its material
                        "Y2");           //its name


  new G4PVPlacement(rmy1,                       //no rotation
                    posY2,                    //at position
                    logicY2,             //its logical volume
                    "Y2",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  G4Para* solidY3 =
   new G4Para("Y3",          //its name 
              1.5*cm,
              4*cm,
              2*cm,
              -30.*deg, 0,0); // its size
   
 G4LogicalVolume* logicY3 =                         
    new G4LogicalVolume(solidY3,         //its solid
                        mat,          //its material
                        "Y3");           //its name
 
  new G4PVPlacement(rmy1,                       //no rotation
                    posY3,                    //at position
                    logicY3,             //its logical volume
                    "Y3",                //its name
                    logicEnv,                //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking


  fScoringVolume = logicY3;

  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
