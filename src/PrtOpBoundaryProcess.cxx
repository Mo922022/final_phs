#include "PrtOpBoundaryProcess.h"

#include "PrtManager.h"

PrtOpBoundaryProcess::PrtOpBoundaryProcess()
  : G4OpBoundaryProcess(){
  fLensId = PrtManager::Instance()->GetLens();
}

G4VParticleChange* PrtOpBoundaryProcess::PostStepDoIt(const G4Track& aTrack, const G4Step& aStep){
  G4StepPoint* pPreStepPoint  = aStep.GetPreStepPoint();
  G4StepPoint* pPostStepPoint = aStep.GetPostStepPoint();
  G4VParticleChange* particleChange = G4OpBoundaryProcess::PostStepDoIt(aTrack, aStep); 
  
  // int parentId = aTrack.GetParentID();
  // std::cout<<"parentId   "<<parentId <<std::endl;
  // if(parentId==1) particleChange->ProposeTrackStatus(fStopAndKill);
  

  // ideal focusing
  if(PrtManager::Instance()->GetLens() == 10){
    G4String ParticleName = aTrack.GetDynamicParticle()->GetParticleDefinition()->GetParticleName();  
    if (ParticleName == "opticalphoton"){
      double endofbar = 1250/2.;
      G4ThreeVector theGlobalPoint1 = pPostStepPoint->GetPosition();
      G4TouchableHistory* touchable = (G4TouchableHistory*)(pPostStepPoint->GetTouchable());
      G4ThreeVector lpoint =  touchable->GetHistory()->GetTransform( 1 ).TransformPoint(theGlobalPoint1);
     
      if(lpoint.getZ() < endofbar+0.0001 && lpoint.getZ() > endofbar-0.0001){
	G4ThreeVector ww  = pPreStepPoint->GetTouchableHandle()->GetHistory()->
	  GetTopTransform().Inverse().TransformPoint(G4ThreeVector(0,0,endofbar));
	if(aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()!="wBar") 
	  particleChange->ProposeTrackStatus(fStopAndKill);
	else{
	  G4Navigator* theNavigator = G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
	  theNavigator->LocateGlobalPointWithinVolume(ww);
	  aParticleChange.ProposePosition(ww.getX(), ww.getY(),ww.getZ());
	}
	return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
      }
    }
  }

  if(PrtManager::Instance()->GetRunType() == 1 && pPostStepPoint->GetPosition().z()<pPreStepPoint->GetPosition().z()){
    particleChange->ProposeTrackStatus(fStopAndKill);
  }
    
    

    
//    // wo prizm Reflections phs
//    if(PrtManager::Instance()->GetRunType() == 13 && pPreStepPoint->GetPhysicalVolume()->GetName()=="wPrizm" && pPostStepPoint->GetPhysicalVolume()->GetName()!="wOpticalGrease" && pPostStepPoint->GetPhysicalVolume()->GetName()!="wBar"&& pPostStepPoint->GetPhysicalVolume()->GetName()!="wLens1" && pPostStepPoint->GetPhysicalVolume()->GetName()!="wLens2" && pPostStepPoint->GetPhysicalVolume()->GetName()!="wLens3"
//      /*     && pPostStepPoint->GetMomentum()== pPreStepPoint->GetMomentum() */
//    ){
//        particleChange->ProposeTrackStatus(fStopAndKill);
//    }
//    
//    // wo prizm Reflections cherenkov
//    if(PrtManager::Instance()->GetRunType() != 13 && pPreStepPoint->GetPhysicalVolume()->GetName()=="wPrizm" && pPostStepPoint->GetPhysicalVolume()->GetName()!="wMcp"  ){
//        //std::cout<< "@@@@@@@@@ geo step  "<< pPreStepPoint->GetPhysicalVolume()->GetName()<<"   "<<  pPostStepPoint->GetPhysicalVolume()->GetName()<<std::endl;
//        particleChange->ProposeTrackStatus(fStopAndKill);
//    }
    
    
    // filter
    if (PrtManager::Instance()->GetRunType() == 13 && pPostStepPoint->GetPhysicalVolume()->GetName()=="wcylinder_filter"){
    
    //std::cout<< "@@@@@@@@@ geo step  "<< pPreStepPoint->GetPhysicalVolume()->GetName()<<"   "<<  pPostStepPoint->GetPhysicalVolume()->GetName()<<std::endl;
        particleChange->ProposeTrackStatus(fStopAndKill);
    }
    

  // if(PrtManager::Instance()->GetRunType() == 0 && pPostStepPoint->GetPhysicalVolume()->GetName()=="wPrizm" && pPostStepPoint->GetPosition().z()<pPreStepPoint->GetPosition().z()){
  //   particleChange->ProposeTrackStatus(fStopAndKill);
  // }

  if(PrtManager::Instance()->GetRunType() == 5 &&  pPreStepPoint->GetPhysicalVolume()->GetName()=="wDirc" && pPostStepPoint->GetPhysicalVolume()->GetName()=="wPrizm" && GetStatus() == FresnelRefraction){
    particleChange->ProposeTrackStatus(fStopAndKill);
  }
  
  // kill photons outside bar and prizm
  if(GetStatus() == FresnelRefraction 
     && aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wDirc"){
    //rd for air gap
    if(PrtManager::Instance()->GetLens()!=4) particleChange->ProposeTrackStatus(fStopAndKill);
  }

  if((aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens1" 
      || aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens2"
      || aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wLens3") 
     &&  aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wDirc"){
    if(PrtManager::Instance()->GetLens()!=4) particleChange->ProposeTrackStatus(fStopAndKill);
  }

  if((aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wOpticalGreased"
      || aStep.GetPreStepPoint()->GetPhysicalVolume()->GetName()=="wOpticalGrease")
     &&  aStep.GetPostStepPoint()->GetPhysicalVolume()->GetName()=="wDirc"){
    if(PrtManager::Instance()->GetLens()==2) particleChange->ProposeTrackStatus(fStopAndKill);
  }

  return particleChange;

}
