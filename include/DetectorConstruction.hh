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
// $Id: DetectorConstruction.hh,v 1.1 2010/10/18 15:56:17 maire Exp $
// GEANT4 tag $Name: geant4-09-04 $
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4Tubs.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Material.hh"
#include "G4UniformMagField.hh"
//#include "DetectorMessenger.hh"
//#include "globals.hh"
// class G4Box;
// class G4LogicalVolume;
// class G4VPhysicalVolume;
// class G4Material;
// class G4UniformMagField;
// class DetectorMessenger;
//--------------------------------------------------------------------------------

class DetectorConstruction : public G4VUserDetectorConstruction {
public:
  
  DetectorConstruction();
  ~DetectorConstruction(){};
  
public:
  void SetTargetMaterial (G4String); 
  void SetMagField(G4double);
  G4VPhysicalVolume* Construct();
public:
  
  G4Material* GetTargetMaterial()  {return TargetMaterial;};  
  const G4VPhysicalVolume* GetphysiWorld() {return physiWorld;};           
  const G4VPhysicalVolume* GetTargetPV()   {return physiTarget;};
  
private:
  bool               background;
  G4Material*        TargetMaterial;
  G4Material*        AlMaterial;
  G4Material*        CMaterial;
  G4Material*        CuMaterial;
  G4Material*        PbMaterial;
  G4Material*        TiMaterial;
  G4Material*        StainlessSteel; 
  
  G4double           TargetSizeX;
  G4double           TargetSizeY;
  G4double           TargetThickness;

  G4Material*        defaultMaterial;
  G4double           WorldSizeYZ;
  G4double           WorldSizeX;
  
  G4Box*             solidWorld;    //pointer to the solid World 
  G4LogicalVolume*   logicWorld;    //pointer to the logical World
  G4VPhysicalVolume* physiWorld;    //pointer to the physical World
  
  G4Box*             solidTarget; //pointer to the solid target
  G4LogicalVolume*   logicTarget; //pointer to the logical target
  G4VPhysicalVolume* physiTarget; //pointer to the physical target
  
  G4UniformMagField* magField;      //pointer to the magnetic field
  
  
private:
  void DefineMaterials();
  
};

//--------------------------------------------------------------------------------

#endif

