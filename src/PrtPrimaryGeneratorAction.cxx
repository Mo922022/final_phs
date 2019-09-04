#include "PrtPrimaryGeneratorAction.h"
#include "PrtPrimaryGeneratorMessenger.h"

#include "Randomize.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include "PrtManager.h"
#include "G4PhysicalVolumeStore.hh" // phs


PrtPrimaryGeneratorAction::PrtPrimaryGeneratorAction():G4VUserPrimaryGeneratorAction(),fParticleGun(0){
    G4int n_particle = 1;
    fParticleGun = new G4ParticleGun(n_particle);
    
    //create a messenger for this class
    fGunMessenger = new PrtPrimaryGeneratorMessenger(this);
    
    //default kinematic
    G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
    fParticleP = particleTable->FindParticle("proton");
    fParticlePi = particleTable->FindParticle("pi+");
    
    fParticleGun->SetParticleDefinition(fParticleP);
    fParticleGun->SetParticleTime(0.0*ns);
    fParticleGun->SetParticlePosition(G4ThreeVector(0.0*cm,0.0*cm,0.0*cm));
    fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));
    fParticleGun->SetParticleEnergy(7*MeV);
    // phs
    
    auto store = G4PhysicalVolumeStore::GetInstance();
    std::cout<<"############ store->size() "<< store->size() <<std::endl; // here
    for (size_t i=0;i<store->size();i++){
        std::cout<<"########################### Geaometry =  "<< (*store)[i]->GetName()<<std::endl; // here
        if((*store)[i]->GetName()=="wDirc")  vdirc =(*store)[i]->GetTranslation();
        if((*store)[i]->GetName()=="wScan") {  vscan =(*store)[i]->GetTranslation();
            //std::cout<<"################## vscan GetCopyNo=  "<< (*store)[i]->GetCopyNo()<<std::endl; // here
            std::cout<<"################## vscan coordinate global =  "<< vdirc+(vscan).rotateY(PrtManager::Instance()->GetAngle()*deg-180*deg)<<std::endl; // here
        }
        if((*store)[i]->GetName()=="wMcp") {  vmcp[(*store)[i]->GetCopyNo()  ] =(*store)[i]->GetTranslation();
            std::cout<<"################## MCP GetCopyNo=  "<< (*store)[i]->GetCopyNo()<<std::endl; // here
        }
        if((*store)[i]->GetName()=="wPixel") {vpix[(*store)[i]->GetCopyNo() -1] =(*store)[i]->GetTranslation();
            std::cout<<"############ Pix GetCopyNo=  "<< (*store)[i]->GetCopyNo() -1<<std::endl; // here
        }
    }
    //    //std::cout<<"############ m["<<mid<<"]" <<" pid = "<< pid<<std::endl; // here
    //    for(auto m=0; m<12; m++){
    //        for(auto p=0; p<64; p++){
    //              vpixminus[p]= G4ThreeVector(vpix[p].x(),vpix[p].y(), vpix[p].z()- 0.6);
    //            //vpixminus[p]= G4ThreeVector(vpix[p].x(),vpix[p].y(), vpix[p].z()- 0.6);
    //            //std::cout<<"###### m["<<m<<"]" <<" vpixminus[p] = "<< vpixminus[p]<<std::endl; // here
    //            gpix[m][p] =vdirc+(vmcp[m]+vpixminus[p]).rotateY(PrtManager::Instance()->GetAngle()*deg-180*deg);
    //        }
    //    }
    ftest1 = PrtManager::Instance()->GetTest1();
    ftest2 = PrtManager::Instance()->GetTest2();
}

PrtPrimaryGeneratorAction::~PrtPrimaryGeneratorAction(){
    delete fParticleGun;
    delete fGunMessenger;
}

void PrtPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){
    G4double x,y,z;
    G4double radiatorL = PrtManager::Instance()->GetRadiatorL();
    G4double radiatorW = PrtManager::Instance()->GetRadiatorW();
    G4double radiatorH = PrtManager::Instance()->GetRadiatorH();
    
    if(PrtManager::Instance()->GetMixPiP()){
        if(PrtManager::Instance()->GetParticle()==211 || PrtManager::Instance()->GetParticle()==0){
            fParticleGun->SetParticleDefinition(fParticleP);
            PrtManager::Instance()->SetParticle(2212);
        }else{
            fParticleGun->SetParticleDefinition(fParticlePi);
            PrtManager::Instance()->SetParticle(211);
        }
    }
    
    PrtManager::Instance()->AddEvent(PrtEvent());
    //std::cout<<"################### PrtManager::Instance()->GetBeamDinsion() "<<PrtManager::Instance()->GetBeamDinsion() <<std::endl; // here
    //std::cout<<"################### PrtManager::Instance()->GetRunType() "<<PrtManager::Instance()->GetRunType() <<std::endl; // here
    
    
    if(PrtManager::Instance()->GetBeamDinsion() == -1){ // random momentum
        fParticleGun->SetParticleMomentum(G4ThreeVector(0, 0, 4.0*GeV*G4UniformRand()));
    }
    
    if(PrtManager::Instance()->GetBeamDinsion() > 0){ // smearing and divergence
        G4double sigma = PrtManager::Instance()->GetBeamDinsion()*mm;
        z = fParticleGun->GetParticlePosition().z();
        
        // // gaussian smearing
        // x = G4RandGauss::shoot(0,sigma);
        // y = G4RandGauss::shoot(0,sigma);
        
        // box smearing
        x = (0.5-G4UniformRand())*sigma;
        y = (0.5-G4UniformRand())*sigma;
        
        fParticleGun->SetParticlePosition(G4ThreeVector(x,y,z));
        PrtManager::Instance()->Event()->SetPosition(TVector3(x,y,z));
        G4double angle = -G4UniformRand()*M_PI;
        G4ThreeVector vec(0,0,1);
        vec.setTheta(G4RandGauss::shoot(0,0.0025)); //beam divergence
        vec.setPhi(2*M_PI*G4UniformRand());
        
        fParticleGun->SetParticleMomentumDirection(vec);
    }
    if(PrtManager::Instance()->GetRunType() == 1){ // LUT generation
        //fParticleGun->SetParticlePosition(G4ThreeVector(radiatorH*(0.5-G4UniformRand()),radiatorW*(0.5-G4UniformRand()),radiatorL/2.-0.1));
        fParticleGun->SetParticlePosition(G4ThreeVector(PrtManager::Instance()->GetRStepY(),//+5-10*G4UniformRand(),
                                                        PrtManager::Instance()->GetRStepX(),//+10-20*G4UniformRand(),
                                                        radiatorL/2.-0.1));
        G4double angle = -G4UniformRand()*M_PI;
        G4ThreeVector vec(0,0,1);
        vec.setTheta(acos(G4UniformRand()));
        vec.setPhi(2*M_PI*G4UniformRand());
        
        //    vec.rotateY(-M_PI/2.);
        fParticleGun->SetParticleMomentumDirection(vec);
    }
    
    // phs space lookup table generation
    if(PrtManager::Instance()->GetRunType() == 13){
        //std::cout<<"############ m["<<mid<<"]" <<" pid = "<< pid<<std::endl; // here
        G4ThreeVector rand(0,0,0);
        for(auto m=0; m<12; m++){
            for(auto p=0; p<64; p++){
                vpixminus[p]= G4ThreeVector(vpix[p].x(),vpix[p].y(), vpix[p].z() ); // -0.6
                //vpixminus[p]= G4ThreeVector(vpix[p].x(),vpix[p].y(), vpix[p].z()- 0.6);
                //std::cout<<"###### m["<<m<<"]" <<" vpixminus[p] = "<< vpixminus[p]<<std::endl; // here
                gpix[m][p] =vdirc+(vmcp[m]+vpixminus[p]).rotateY(PrtManager::Instance()->GetAngle()*deg-180*deg);
                // square randmization
                G4double x_pix_length(3.3125*2), y_pix_length(3.3125*2), z_pix_length(0.05*2);
                G4double sq_diff_length_x(x_pix_length), sq_diff_length_y(y_pix_length), sq_rand_varx((G4UniformRand()-0.5)*x_pix_length), sq_rand_vary ((G4UniformRand()-0.5)*y_pix_length);
                G4ThreeVector rand= G4ThreeVector(sq_rand_varx,sq_rand_vary, 0).rotateY(PrtManager::Instance()->GetAngle()*deg-180*deg);
                gpix[m][p]=gpix[m][p]+ rand;
                vpixminus[p]=vpixminus[p]+ rand;
            }
        }
        PrtManager::Instance()->AddEvent(PrtEvent());
        ///////////////
        // pos dir ////
        ///////////////
        fParticleGun->SetParticlePosition(gpix[ftest1][ftest2]);
        
        //PrtManager::Instance()->setStartVertix(vpixminus[ftest2]); // add to stepping Acttion
        ///////////////
        // mom dir ////
        ///////////////
        G4double angle = -G4UniformRand()*M_PI;
        G4ThreeVector vec(0,0,1);
        vec.setTheta(acos(G4UniformRand()));
        vec.setPhi(2*M_PI*G4UniformRand());
        vec.rotateY(PrtManager::Instance()->GetAngle()*deg-180*deg);
        fParticleGun->SetParticleMomentumDirection(-vec);
    }
    
    if(PrtManager::Instance()->GetRunType() == 5){ // calibration light
        G4double shift = PrtManager::Instance()->GetShift();
        
        fParticleGun->SetParticlePosition(G4ThreeVector(-radiatorL/2.+0.1-shift,0,5+tan(45*M_PI/180.)*shift+25));
        G4double angle = -G4UniformRand()*M_PI;
        G4ThreeVector vec(0,0,1);
        vec.setTheta(acos(G4UniformRand()));
        vec.setPhi(2*M_PI*G4UniformRand());
        
        vec.rotateY(-M_PI/2.);
        fParticleGun->SetParticleMomentumDirection(vec);
    }
    if(PrtManager::Instance()->GetRunType() == 6){ // for determining focal plane of the lens
        G4double shiftx = -radiatorL/2.+0.1;
        G4double shifty = radiatorW/2. - G4UniformRand()*radiatorW;
        G4double shiftz = radiatorH/2. - G4UniformRand()*radiatorH;
        
        G4double angle = 0.7*(M_PI/2.-G4UniformRand()*M_PI);
        G4ThreeVector vec(0,0,1);
        vec.setTheta(angle);
        ///vec.setTheta(acos(G4UniformRand()));
        //vec.setPhi(2*M_PI*G4UniformRand());
        //std::cout<<"angle "<<angle*180/M_PI <<std::endl;
        
        G4double lensThickness=15;
        G4double separation=PrtManager::Instance()->GetBeamDinsion();
        if(separation<0.001) separation=10;
        G4double rotShiftX=0.5*separation*std::cos(angle)+(0.5*lensThickness+0.1)*std::tan(angle);
        G4double rotShiftY=-0.5*radiatorL +0.1;
        
        //fParticleGun->SetParticlePosition(G4ThreeVector(shiftx,shifty,shiftz));
        //fParticleGun->SetParticlePosition(G4ThreeVector(rotShiftY, 0,-rotShiftX));
        fParticleGun->SetParticlePosition(G4ThreeVector(rotShiftY, 0,0.5*separation));
        
        
        vec.rotateY(-M_PI/2.);
        fParticleGun->SetParticleMomentumDirection(vec);
        
        fParticleGun->GeneratePrimaryVertex(anEvent);
        rotShiftX=-0.5*separation*std::cos(angle)+(0.5*lensThickness+0.1)*std::tan(angle);
        shiftx = -radiatorL/2.+0.1;
        shifty = radiatorW/2. - G4UniformRand()*radiatorW;
        shiftz = radiatorH/2. - G4UniformRand()*radiatorH;
        //fParticleGun->SetParticlePosition(G4ThreeVector(shiftx,shifty,shiftz));
        //fParticleGun->SetParticlePosition(G4ThreeVector(rotShiftY,0,-rotShiftX));
        fParticleGun->SetParticlePosition(G4ThreeVector(rotShiftY, 0,-0.5*separation));
    }
    
    fParticleGun->GeneratePrimaryVertex(anEvent);
    
    G4ThreeVector dir = fParticleGun->GetParticleMomentumDirection();
    dir *= fParticleGun->GetParticleMomentum();
    PrtManager::Instance()->SetMomentum(TVector3(dir.x(),dir.y(),dir.z()));
}

void PrtPrimaryGeneratorAction::SetOptPhotonPolar(){
    G4double angle = G4UniformRand() * 360.0*deg;
    SetOptPhotonPolar(angle);
}

void PrtPrimaryGeneratorAction::SetOptPhotonPolar(G4double angle){
    if (fParticleGun->GetParticleDefinition()->GetParticleName()!="opticalphoton"){
        G4cout << "--> warning from PrimaryGeneratorAction::SetOptPhotonPolar() :"
        "the particleGun is not an opticalphoton " <<
        fParticleGun->GetParticleDefinition()->GetParticleName()<< G4endl;
        return;
    }
    
    G4ThreeVector normal (1., 0., 0.);
    G4ThreeVector kphoton = fParticleGun->GetParticleMomentumDirection();
    G4ThreeVector product = normal.cross(kphoton);
    G4double modul2       = product*product;
    
    G4ThreeVector e_perpend (0., 0., 1.);
    if (modul2 > 0.) e_perpend = (1./std::sqrt(modul2))*product;
    G4ThreeVector e_paralle    = e_perpend.cross(kphoton);
    
    G4ThreeVector polar = std::cos(angle)*e_paralle + std::sin(angle)*e_perpend;
    fParticleGun->SetParticlePolarization(polar);
}


