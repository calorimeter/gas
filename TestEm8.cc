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
/// \file electromagnetic/TestEm8/TestEm8.cc
/// \brief Main program of the electromagnetic/TestEm8 example
//
// $Id: TestEm8.cc 85243 2014-10-27 08:22:42Z gcosmo $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


#include "G4RunManager.hh"

#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "FTFP_BERT.hh"
#include "LBE.hh"

#include "Randomize.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

//simulation-specific files
#include "DetectorConstruction.hh"
#include "ActionInitialization.hh"
#include "PhysicsList.hh"
#include "OpNovicePhysicsList.hh"

namespace {
  void PrintUsage() {
    G4cerr << " Wrong Usage! Program Terminated [placeholder] " << G4endl;
  }
}

int main(int argc,char** argv)
{
  // Evaluate arguments
  //
  int amtArgs = 2;
  if ( argc-1 > amtArgs*2 ) {
    PrintUsage();
    return 1;
  }

  G4String macro;
  G4String session;

  for ( G4int i=1; i<argc; i=i+2 ) {
    if      ( G4String(argv[i]) == "-m" ) macro = argv[i+1];
    else if ( G4String(argv[i]) == "-u" ) session = argv[i+1];
    else {
      PrintUsage();
      return 1;
    }
  }

  // Detect interactive mode (if no macro provided) and define UI session
  //
  G4UIExecutive* ui = 0;
  if ( ! macro.size() ) {
    ui = new G4UIExecutive(argc, argv);
  }

  // Choose the Random engine
  //
  G4Random::setTheEngine(new CLHEP::RanecuEngine);

  // Construct the default run manager
  //
  G4RunManager * runManager = new G4RunManager;

  // Set mandatory initialization classes
  //
  DetectorConstruction* detConstruction = new DetectorConstruction(); //Add onto this
  runManager->SetUserInitialization(detConstruction);

  //G4VModularPhysicsList* physicsList = new FTFP_BERT;
  G4VUserPhysicsList* physicsList = new OpNovicePhysicsList();
  //G4VUserPhysicsList* physicsList = new PhysicsList();
  runManager->SetUserInitialization(physicsList);

  ActionInitialization* actionInitialization
     = new ActionInitialization();
  runManager->SetUserInitialization(actionInitialization);

  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  // Process macro or start UI session
  //
  if ( macro.size() ) {
    // batch mode
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command+macro);
  } else  {
    // interactive mode : define UI session
    UImanager->ApplyCommand("/control/execute init_vis.mac");
    if (ui->IsGUI()) {
      UImanager->ApplyCommand("/control/execute gui.mac");
    }
    ui->SessionStart();
    delete ui;
  }

  delete visManager;
  delete runManager;
}

/*
#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "G4UIcommand.hh"
#include "Randomize.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "ActionInitialization.hh"

//#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
//#endif

//#ifdef G4UI_USE
#include "G4UIExecutive.hh"
//#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv) 
{
  //choose the Random engine
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  
//#ifdef G4MULTITHREADED
//  G4MTRunManager* runManager = new G4MTRunManager;
//  // Number of threads can be defined via 3rd argument
//  G4int nThreads = 2;
//  if (argc==3) {
//    if(G4String(argv[2]) == "NMAX") {
//      nThreads = G4Threading::G4GetNumberOfCores();
//    } else {
//      nThreads = G4UIcommand::ConvertToInt(argv[2]);
//    }
//  }
//  if (nThreads > 0) { runManager->SetNumberOfThreads(nThreads); }
//  G4cout << "===== TestEm8 is started with "
//         <<  runManager->GetNumberOfThreads() << " threads =====" << G4endl;
//#else
  G4RunManager* runManager = new G4RunManager;
//#endif

  // set mandatory initialization classes
  runManager->SetUserInitialization(new PhysicsList);
  runManager->SetUserInitialization(new DetectorConstruction());
 
  // set user action classes
  runManager->SetUserInitialization(new ActionInitialization());
  
  G4UImanager* UI = G4UImanager::GetUIpointer();  

  if (argc!=1)   // batch mode  
    {
      G4String command = "/control/execute ";
      G4String fileName = argv[1];
      UI->ApplyCommand(command+fileName);
    }
    
  else           //define visualization and UI terminal for interactive mode
    { 
#ifdef G4VIS_USE
      G4VisManager* visManager = new G4VisExecutive;
      visManager->Initialize();
#endif    
     
#ifdef G4UI_USE
      G4UIExecutive * ui = new G4UIExecutive(argc,argv);      
      ui->SessionStart();
      delete ui;
#endif
     
#ifdef G4VIS_USE
      delete visManager;
#endif     
    } 
   
  // job termination
  //
  delete runManager;

  return 0;
}
*/
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

