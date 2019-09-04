#include "PrtPixelSD_dummy.h"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4RunManager.hh"


#include "PrtEvent.h"
#include "PrtPrizmHit.h"

#include "PrtRunAction.h"
#include "PrtManager.h"
#include "PrtPrimaryGeneratorAction.h"

#include "TH2F.h"
#include "TVector3.h"

PrtPixelSD_dummy::PrtPixelSD_dummy( const G4String& name,
                                   const G4String& hitsCollectionName,
                                   G4int nofCells
                                   )
: G4VSensitiveDetector(name){
    collectionName.insert(hitsCollectionName);
}

PrtPixelSD_dummy::~PrtPixelSD_dummy(){
}

void PrtPixelSD_dummy::Initialize(G4HCofThisEvent* hce){
    
}

G4bool PrtPixelSD_dummy::ProcessHits(G4Step* step, G4TouchableHistory* hist){
    
    
    if(step == 0) return false;
    
    //G4ThreeVector translation = hist->GetTranslation();
    //G4ThreeVector localpos = step->GetPreStepPoint()->GetPhysicalVolume()->GetObjectTranslation();
    G4TouchableHistory* touchable = (G4TouchableHistory*)(step->GetPostStepPoint()->GetTouchable());
    
    // Get cell id
    G4int layerNumber = touchable->GetReplicaNumber(0);
    //G4cout<< "###### PixelId = "<<layerNumber << G4endl;
    G4Track* track = step->GetTrack();

    const G4DynamicParticle* dynParticle = track->GetDynamicParticle();
    G4ParticleDefinition* particle = dynParticle->GetDefinition();
    G4String ParticleName = particle->GetParticleName();
    
    G4ThreeVector globalpos = step->GetPostStepPoint()->GetPosition();
    G4ThreeVector localpos = touchable->GetHistory()->GetTopTransform().TransformPoint(globalpos);
    
    G4ThreeVector g4pos = track->GetVertexPosition();
    
    
    TVector3 localPos(localpos.x(),localpos.y(),localpos.z());
    
    Double_t x(localPos.x()), y(localPos.y()), z(localPos.z());
    
    G4String VolumName = step->GetPostStepPoint()->GetPhysicalVolume()->GetName();
    
    if (VolumName == "wPixel") track->SetTrackStatus(fStopAndKill);
    
    //G4cout<< VolumName << G4endl;
    G4double stepLength = step->GetStepLength();
    //G4cout<< "###### x = "<<x << " y = "<<y <<" z = "<<z << G4endl;
    //G4cout<< "stepLength = "<<stepLength << G4endl;
    
    
    
    
    return true;
}

void PrtPixelSD_dummy::EndOfEvent(G4HCofThisEvent*){
    
}

