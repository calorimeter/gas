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
/// \file electromagnetic/TestEm5/src/StackingAction.cc
/// \brief Implementation of the StackingAction class
//
// $Id: StackingAction.cc,v 1.8 2009-03-06 18:04:23 maire Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "StackingAction.hh"

#include "Run.hh"
#include "StackingMessenger.hh"

#include "G4Track.hh"
#include "G4Step.hh"
#include "G4RunManager.hh"
#include "G4Gamma.hh"
#include "G4OpticalPhoton.hh"
#include "DetectorConstruction.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::StackingAction()
  : G4UserStackingAction(),
    fKillSecondary(false),
    fStackMessenger(0)
{
  fStackMessenger = new StackingMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

StackingAction::~StackingAction()
{
  delete fStackMessenger;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void StackingAction::PrepareNewEvent()
{


}

void StackingAction::NewStage()
{

}

G4ClassificationOfNewTrack
StackingAction::ClassifyNewTrack(const G4Track* aTrack)
{

	if(aTrack->GetCreatorProcess() && aTrack->GetCreatorProcess()->GetProcessName() != "eIoni")
		G4cout << "created a new " << aTrack->GetParticleDefinition()->GetParticleName() << " of " << aTrack->GetKineticEnergy() << " MeV through " << aTrack->GetCreatorProcess()->GetProcessName() << G4endl;
	if( aTrack->GetDefinition()==G4OpticalPhoton::OpticalPhotonDefinition() )
	{
		DetectorConstruction* c = (DetectorConstruction*) G4RunManager::GetRunManager()->GetUserDetectorConstruction();
		//G4cout << " photon with " << aTrack->GetKineticEnergy()<< ", "x
		//		<< aTrack->GetMomentumDirection().cosTheta()
		//		<< ", " << aTrack->GetCreatorProcess()->GetProcessName() <<G4endl;
		G4double z_max = ( c->GetAbsorberThickness());
		//G4double z_min = aTrack->GetVolume()->GetLogicalVolume()->GetSolid()->
		G4double z = aTrack->GetPosition().z();

		Run* run = (Run*) G4RunManager::GetRunManager()->GetCurrentRun();
		if(z > z_max){

			G4double x = aTrack->GetPosition().x();
			G4double y = aTrack->GetPosition().y();
			G4double r = sqrt(x*x+y*y);
			//G4cout << z_max << " a " << z << " and " << r << G4endl;
			if(r >=0 && r <= c->GetGasRadius()){
				run->fGammaCounter++;
				run->fGammaEnergy += aTrack->GetKineticEnergy();
				run->fAngle += r;
			}
		} else {
			G4double r = tan(aTrack->GetMomentumDirection().theta()) * (z_max - z);
			//G4cout << z_max << " b " << z << " and " << r << G4endl;
			if(r >=0 && r <= c->GetGasRadius()){
				run->fGammaCounter++;
				run->fGammaEnergy += aTrack->GetKineticEnergy();
				run->fAngle += r;
			}
		}
	}

  //stack or delete secondaries
  G4ClassificationOfNewTrack status = fUrgent;

  //keep primary particle
  if (aTrack->GetParentID() == 0 || !fKillSecondary) { return status; }
  
  Run* run
    = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());  

  // charged tracks are killed only inside sensitive volumes
  if(aTrack->GetVolume()->GetLogicalVolume()->GetSensitiveDetector() &&
     aTrack->GetDefinition()->GetPDGCharge() != 0.0) 
    {
      run->AddEnergy(aTrack->GetKineticEnergy(), 0); 
      status = fKill;    
    }
  return status;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
