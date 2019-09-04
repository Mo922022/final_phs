// -----------------------------------------
// PrtLutReco.cpp
//
// Created on: 13.07.2013
// Author: R.Dzhygadlo at gsi.de
// -----------------------------------------

#include "PrtLutReco.h"

#include "PrtManager.h"

#include "PrtLutNode.h"
#include "PrtTrackInfo.h"
#include "PrtPhotonInfo.h"
#include "PrtAmbiguityInfo.h"

#include "TROOT.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TRotation.h"
#include "TGraph.h"
#include <TVirtualFitter.h>
#include <TArc.h>
#include <TLegend.h>


#include "TH3F.h"
#define prt__sim
#include "../../prttools/datainfo.C"
#include "../../prttools/prttools.C"

using std::cout;
using std::endl;


////

TH1F* hist_total_time_0bar= new TH1F("hist_total_time_0bar",";t_{total} [ns];entries [#]", 250,0,50);
TH1F* hist_bar_time_0bar= new TH1F("hist_bar_time_0bar",";t_{in bar} [ns];entries [#]", 250,0,50);
TH1F* hist_hit_time_0bar= new TH1F("hist_hit_time_0bar",";t_{measured} [ns];entries [#]", 250,0,50);
TH1F* hist_hit_phs_0bar= new TH1F("hist_hit_phs_0bar",";t_{phs 0bar} [ns];entries [#]", 250,0,50);
TH1F* hist_timeDiff_angle_cut_0bar= new TH1F("hist_timeDiff_angle_cut_0bar",";t_{diff angle cut} [ns];entries [#]", 500,-15,15);

TH2F*  hist_time_angle_all_0bar = new TH2F("hist_time_angle_all_0bar",";time [ns];solution angle [rad]", 250,0,50,800,0,4 );
TH2F*  hist_time_angle_all_correlation_0bar = new TH2F("hist_time_angle_all_correlation_0bar","; #Delta time [ns]; Angle [rad]", 1000,-25,25,2000,0,4 ); //1000,-25,25,2000,0,4


TH2F*  hist_time_angle_all_correlation_0bar_zoom = new TH2F("hist_time_angle_all_correlation_0bar","; #Delta time [ns]; Angle [rad]", 80,-2,2,200,0.6,1 );



TH1F* hist_timeDiff_0bar= new TH1F("hist_timeDiff_0bar",";t_{diff} [ns];entries [#]", 500,-15,15);
TH1F* hist_cherenkov_0bar= new TH1F("hist_cherenkov_0bar",  "cherenkov angle;#theta_{C} [rad];entries [#]",150,0.6,1); //150  //80 ,300,0,5 );//

////////

TH1F*  hist_cherenkov_phs = new TH1F("RecoCherenkov_phs",  "cherenkov angle;#theta_{C} [rad];entries [#]",150,0.6,1); //150  //80 ,300,0,5 );//
TH1F*  hist_cherenkov_phs_dir_cut = new TH1F("hist_cherenkov_phs_dir_cut",  "cherenkov angle;#theta_{C} [rad];entries [#]",150,0.6,1); //150  //80 ,300,0,5 );//
TH1F*  hist_cherenkov_phs_bg = new TH1F("hist_cherenkov_phs_bg",  "cherenkov angle;#theta_{C} [rad];entries [#]",150,0.6,1); //150  //80 ,300,0,5 );//

TH1F*  hist_cherenkov_lut = new TH1F("RecoCherenkov_lut",  "cherenkov angle;#theta_{C} [rad];entries [#]",150,0.6,1); //150  //80 ,300,0,5 );//

TH1F*  hist_timeDiff_phs = new TH1F("timeDiff_phs",";t_{phs}-t_{measured} [ns];entries [#]", 500,-15,15);
TH1F*  hist_timeDiff_lut = new TH1F("timeDiff_lut",";t_{lut}-t_{measured} [ns];entries [#]", 500,-15,15);

TH1F*  hist_phsTime = new TH1F("phsTime",";t_{phs} [ns];entries [#]", 250,0,50);
TH1F*  hist_hitTimed = new TH1F("hitTimed",";t_{measured} [ns];entries [#]", 250,0,50);
TH2F*  hist_time_angle_all_phs = new TH2F("hist_time_angle_all_phs",";time [ns];solution angle [rad]", 250,0,50,800,0,4 );

TH2F*  hist_time_angle_all_correlation_phs = new TH2F("hist_time_angle_all_correlation_phs","; #Delta time [ns]; Angle [rad]", 1000,-25,25,2000,0,4 ); //1000,-25,25,2000,0,4
TH2F*  hist_time_angle_all_correlation_phs_zoom = new TH2F("hist_time_angle_all_correlation_phs_zoom","; #Delta time [ns]; Angle [rad]", 80,-2,2,200,0.6,1 );

TH2F*  hist_time_angle_all_correlation_phs_center = new TH2F("hist_time_angle_all_correlation_phs_center","; #Delta time [ns]; Angle [rad]", 1000,-25,25,2000,0,4 );

TH2F*  hist_time_angle_all_lut = new TH2F("hist_time_angle_all_lut",";time [ns];solution angle [rad]", 250,0,50,800,0,4 );
TH2F*  hist_time_angle_all_correlation_lut = new TH2F("hist_time_angle_all_correlation_lut","; #Delta time [ns]; Angle [rad]", 1000,-25,25,2000,0,4 );
TH2F*  hist_time_angle_all_correlation_lut_zoom = new TH2F("hist_time_angle_all_correlation_lut_zoom","; #Delta time [ns]; Angle [rad]", 80,-2,2,200,0.6,1 );

TH1F*  hist_phs_solution_number = new TH1F("hist_phs_solution_number",";number of solutions [#];entries [#]", 500,0,500);
TH1F*  hist_lut_solution_number = new TH1F("hist_lut_solution_number",";number of solutions [#];entries [#]", 500,0,500);

TH1F*  hit_b4 = new TH1F("hit_b4",";t_{measured} [ns];entries [#]", 250,0,50);
TH1F*  hit_phs = new TH1F("hit_phs",";t_{measured} [ns];entries [#]", 250,0,50);
TH1F*  hit_lut = new TH1F("hit_lut",";t_{measured} [ns];entries [#]", 250,0,50);

TH1F*  fHistPhotonEnergy = new TH1F("fHistPhotonEnergy",";|#alpha|[degree];entries [#]", 200, 0, 8);

/*
 TH1F*  hist_dir_x = new TH1F("hist_dir_x",";dir x component ;entries [#]", 100,-1.0,1.0);
 TH1F*  hist_dir_y = new TH1F("hist_dir_y",";dir y component ;entries [#]", 100,-1.0,1.0);
 TH1F*  hist_dir_z = new TH1F("hist_dir_z",";dir z component;entries [#]", 100,-1.0,1.0);
 */
TH2F*  hist_dir_x = new TH2F("hist_dir_x",";dir x component ;Reco. angel [#]", 400,-1.0,1.0, 800, 0, 4);
TH2F*  hist_dir_y = new TH2F("hist_dir_y",";dir y component ;Reco. angel [#]", 400,-1.0,1.0, 800, 0, 4);
TH2F*  hist_dir_z = new TH2F("hist_dir_z",";dir z component;Reco. angel [#]",  400,-1.0,1.0, 800, 0, 4);
TH2F*  hist_dir_xy = new TH2F("hist_dir_xy",";dir y component;dir x component",  400,-1.0,1.0, 400,-1.0,1.0);
TH2F*  hist_dir_xy_test = new TH2F("hist_dir_xy_test",";dir y component;dir x component",  400,-1.0,1.0, 400,-1.0,1.0);
TH2F*  hist_pos_phs_x_angle = new TH2F("hist_pos_phs_x_angle",";pos x component [mm];Reco angel [rad]",  800,-30,30,400, 0, 4);
TH2F*  hist_pos_phs_y_angle = new TH2F("hist_pos_phs_y_angle",";pos y component [mm];Reco angel [rad]",  800,-30,30,400, 0, 4);
TH3F *hist_dir_xyz = new TH3F("hist_dir_xyz","hist_dir_xyz fitz",100,-1.0,1.0, 100,-1.0,1.0,400, 0, 4);


TH2F* hist_pos_xy= new TH2F("hist_pos_xy",";pos x [mm];pos y [mm]",  100,-30,30,100,-30,30);
TH2F* hist_pos_xy_test= new TH2F("hist_pos_xy_test",";pos x [mm];pos y [mm]",  100,-30,30,100,-30,30);

TH2F*  hist_dir_x_lut= new TH2F("hist_dir_x_lut",";dir x component ;Reco. angel [#]", 400,-1.0,1.0, 800, 0, 4);
TH2F*  hist_dir_y_lut = new TH2F("hist_dir_y_lut",";dir y component ;Reco. angel [#]", 400,-1.0,1.0, 800, 0, 4);
TH2F*  hist_dir_z_lut = new TH2F("hist_dir_z_lut",";dir z component;Reco. angel [#]",  400,-1.0,1.0, 800, 0, 4);
TH2F*  hist_dir_xy_lut = new TH2F("hist_dir_xy_lut",";dir y component;dir x component",  400,-1.0,1.0, 400,-1.0,1.0);


TH2F* hist_dir_xy_angle = new TH2F("hist_dir_xy_angle",";dir y component;dir x component",  400,-1.0,1.0, 400,-1.0,1.0);
TH2F* hist_pos_phs_xy_angle= new TH2F("hist_pos_phs_xy_angle",";pos x [mm];pos y [mm]",  100,-30,30,100,-30,30);

TH2F* hist_dir_xy_time = new TH2F("hist_dir_xy_time",";dir y component;dir x component",  400,-1.0,1.0, 400,-1.0,1.0);
TH2F* hist_dir_xy_time_time = new TH2F("hist_dir_xy_time_time",";dir y component;dir x component",  400,-1.0,1.0, 400,-1.0,1.0);


TF1 *fcut = new TF1("fcut", ".76*x*x-.77", -1, 1);
TF1 *fcut2 = new TF1("fcut2", ".9*x*x-.7", -1, 1);

TF1 *fcut3 = new TF1("fcut3", "exp(-[0]/(x-[1])+[2])", 0.0, 50);
//TF1 *fcut = new TF1("fcut", ".76*x*x-.82", -1, 1);
//TF1 *fcut2 = new TF1("fcut2", ".9*x*x-.6", -1, 1);


TH1F*  fHist0 = new TH1F("timediff",";t_{calc}-t_{measured} [ns];entries [#]", 500,-10,10);
TH1F*  fHist0i = new TH1F("timediffi",";t_{calc}-t_{measured} [ns];entries [#]", 500,-10,10);
TH1F*  fhNph = new TH1F("fhNph",";detected photons [#];entries [#]", 150,0,150);
TH1F*  fhNph_pi = new TH1F("fhNph_pi",";detected photons [#];entries [#]", 150,0,150);
TH1F*  fhNph_p = new TH1F("fhNph_p",";detected photons [#];entries [#]", 150,0,150);
TH1F*  fHist1 = new TH1F("time1",";measured time [ns];entries [#]",   1000,0,100);
TH1F*  fHist2 = new TH1F("time2",";calculated time [ns];entries [#]", 1000,0,100);
TH1F*  fHist6 = new TH1F("time6",";measured time [ns];entries [#]", 1000,0,100);

TH2F*  fHist3 = new TH2F("time3",";calculated time [ns];measured time [ns]", 500,0,80, 500,0,40);
TH2F*  fHist4 = new TH2F("time4",";#theta_{c}sin(#varphi_{c});#theta_{c}cos(#varphi_{c}", 100,-1,1, 100,-1,1);
TH2F*  fHist5 = new TH2F("time5",";#theta_{c}sin(#varphi_{c});#theta_{c}cos(#varphi_{c}", 100,-1,1, 100,-1,1);

//TH1F *hLnDiffP = new TH1F("hLnDiffP",  ";ln L(p) - ln L(#pi);entries [#]",200,-30,30);
//TH1F *hLnDiffPi = new TH1F("hLnDiffPi",";ln L(p) - ln L(#pi);entries [#]",200,-30,30);

TH1F *hLnDiffP = new TH1F("hLnDiffP",  ";ln L(p) - ln L(#pi);entries [#]",200,-300,300);
TH1F *hLnDiffPi = new TH1F("hLnDiffPi",";ln L(p) - ln L(#pi);entries [#]",200,-300,300);

TH1F *hLnDiffP_phs = new TH1F("hLnDiffP_phs",  ";ln L(p) - ln L(#pi);entries [#]",200,-300,300);
TH1F *hLnDiffPi_phs = new TH1F("hLnDiffPi_phs",";ln L(p) - ln L(#pi);entries [#]",200,-300,300);

//TH1F *hLnDiffP_phs = new TH1F("hLnDiffP_phs",  ";ln L(p) - ln L(#pi);entries [#]",200,-2000,2000);
//TH1F *hLnDiffPi_phs = new TH1F("hLnDiffPi_phs",";ln L(p) - ln L(#pi);entries [#]",200,-2000,2000);

TH1F *hLnDiffP_pdf = new TH1F("hLnDiffP_pdf",  ";ln L(p) - ln L(#pi);entries [#]",200,-30,30);// 30000
TH1F *hLnDiffPi_pdf = new TH1F("hLnDiffPi_pdf",";ln L(p) - ln L(#pi);entries [#]",200,-30,30);// 30000

TF1 *gF1 = new TF1("gaus0","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0.7,0.9);
TF1 *gF2= new TF1("gaus0","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0.7,0.9);

Int_t gg_i(0), gg_ind(0);
TGraph gg_gr;
PrtLutNode *fLutNode[5000];
PrtLutNode *fLutNode_phs[5000];

TH1F*  fHistMcp[15];
TH1F*  fHistCh[960];

TH2F* hist_dir_x_tcut_test[21],
*hist_dir_y_tcut_test[21],
*hist_dir_xy_tcut_test[21],
*hist_dir_xy_time_tcut_test[21],
*hist_dir_xy_occy_tcut_test[21],
*hist_dir_xy_angle_tcut_test[21],
*hist_dir_xy_tdiff_tcut_test[21];


TH1F* hist_cherenkov_phs_tcut_test[21];








// -----   Default constructor   -------------------------------------------
PrtLutReco::PrtLutReco(TString infile, TString lutfile, Int_t verbose) {
    
    //TString phs_file_path = "/Users/ahmed/Desktop/std/prtdirc/build/1m_l3_20_sub_nobug.root";//final_1m_lut_all_20_447_l6_randPix.root";
    //TString phs_file_path = "/Users/ahmed/Desktop/std/prtdirc/build/10m_l0_1mf.root";//final_1m_lut_all_20_447_l6_randPix.root";
    //TString phs_file_path = "/Users/ahmed/Desktop/std/prtdirc/build/10m_l0_1mf_sub.root";//final_1m_lut_all_20_447_l6_randPix.root";
    //TString phs_file_path = "/Users/ahmed/Desktop/std/prtdirc/build/10m_l0_1mf_sub_0.25mf_4entries.root";
    
    //1m_l3_180_inbar_final.root
    //1m_l3_180_0bar_final.root
    //10m_l3_180_0bar_f_final.root
    //1m_l3_20_final.root
    //10m_l3_20_f_final.root
    //1m_l3_20_rot_final.root
    //10m_l3_20_f_rot_final.root
    
    
    //TString phs_file_path = "/Users/ahmed/final_phs/build/10m_phs_l3_2017_90_f0.2.root";
    TString phs_file_path = "/Users/ahmed/final_phs/build/10m_phs_l3_2017_90_fixpos_avr_thr10.root";
    
    
    //TString phs_file_path = "/Users/ahmed/final_phs/build/10m_phs_l3_2017_90_f0.2_avr.root";
    //TString phs_file_path = "/Users/ahmed/final_phs/build/1m_phs_l3_2017_90_avr.root";
    
    //TString phs_file_path = "/Users/ahmed/final_phs/build/10m_phs_l3_2017_90_f0.2_avr_thr10.root";
    //TString phs_file_path = "/Users/ahmed/final_phs/build/1m_phs_l3_2017_90_avr_thr10.root";
    
    fFile_phs = new TFile(phs_file_path);
    fTree_phs=(TTree *) fFile_phs->Get("prtlut") ;
    fLut_phs = new TClonesArray("PrtLutNode");
    fTree_phs->SetBranchAddress("LUT",&fLut_phs);
    fTree_phs->GetEntry(0);
    
    
    
    fVerbose = verbose;
    fChain = new TChain("data");
    fChain->Add(infile);
    fChain->SetBranchAddress("PrtEvent", &fEvent);
    
    fFile = new TFile(lutfile);
    fTree=(TTree *) fFile->Get("prtlut") ;
    fLut = new TClonesArray("PrtLutNode");
    fTree->SetBranchAddress("LUT",&fLut);
    fTree->GetEntry(0);
    
    fHist = new TH1F("chrenkov_angle_hist",  "chrenkov angle;#theta_{C} [rad];entries [#]", 150,0.6,1); //150
    fHistPi = new TH1F("chrenkov_angle_hist_Pi",  "chrenkov angle pi;#theta_{C} [rad];entries [#]", 150,0.6,1); //150
    fHisti = new TH1F("chrenkov_angle_histi","chrenkov angle;#theta_{C} [rad];entries [#]", 80,0.6,1); //150
    fFit = new TF1("fgaus","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
    fSpect = new TSpectrum(10);
    fRadiator=1;
    
    if(infile.Contains("beam_")) {
        TString fileid(infile);
        fileid.Remove(0,fileid.Last('/')+1);
        fileid.Remove(fileid.Last('.')-1);
        prt_data_info = getDataInfo(fileid);
        fRadiator =  prt_data_info.getRadiatorId();
        
        TString opath(infile);
        opath.Remove(opath.Last('/'));
        if(infile.Contains("C.root")) {
            prt_savepath = opath+Form("/%dr/%d",prt_data_info.getStudyId(),prt_data_info.getFileId());
        } else {
            prt_savepath = opath+Form("/%ds/%d",prt_data_info.getStudyId(),prt_data_info.getFileId());
        }
        
    } else prt_savepath="data/sim";
    std::cout<<"prt_savePath  "<< prt_savepath <<std::endl;
    
    for(Int_t i=0; i<5000; i++) {
        fLutNode[i] = (PrtLutNode*) fLut->At(i);
        fLutNode_phs[i] = (PrtLutNode*) fLut_phs->At(i);
    }
    cout << "-I- PrtLutReco: Intialization successfull" << endl;
    
    for(Int_t i=0; i<prt_nmcp; i++) {
        fHistMcp[i] = new TH1F(Form("fHistMcp_%d",i),Form("fHistMcp_%d;#theta_{C} [rad];entries [#]",i), 160,0.6,1); //150
    }
    
    for(Int_t i=0; i<960; i++) {
        fHistCh[i] = new TH1F(Form("fHistCh_%d",i),Form("fHistCh_%d;#theta_{C} [rad];entries [#]",i), 100,0.6,1); //150
    }
    
    
    
    for(Int_t i=0; i<21; i++) {
        //hist_dir_x_tcut_test[i] = new TH2F(Form("hist_dir_x_tcut_test_%d",i),Form("hist_dir_x_tcut_test_%d;dir x component ;Reco. angel [#]",i), 400,-1.0,1.0, 800, 0, 4);
        //hist_dir_y_tcut_test[i] = new TH2F(Form("hist_dir_y_tcut_test_%d",i),Form("hist_dir_y_tcut_test_%d;dir y component ;Reco. angel [#]",i), 400,-1.0,1.0, 800, 0, 4);
        
        hist_dir_xy_tcut_test[i] = new TH2F(Form("hist_dir_xy_tcut_test_%d",i),Form("hist_dir_xy_tcut_test_%d;dir y component;dir x component",i),  400,-1.0,1.0, 400,-1.0,1.0);
        hist_dir_xy_time_tcut_test[i] = new TH2F(Form("hist_dir_xy_time_tcut_test_%d",i),Form("hist_dir_xy_time_tcut_test_%d;dir y component;dir x component",i),  400,-1.0,1.0, 400,-1.0,1.0);
        
        hist_dir_xy_occy_tcut_test[i] = new TH2F(Form("hist_dir_xy_test_tcut_test_%d",i),Form("hist_dir_xy_time_tcut_test_%d;dir y component;dir x component",i),  400,-1.0,1.0, 400,-1.0,1.0);
        //hist_dir_xy_angle_tcut_test[i] = new TH2F(Form("hist_dir_xy_angle_tcut_test_%d",i),Form("hist_dir_xy_angle_tcut_test_%d;dir y component;dir x component",i),  400,-1.0,1.0, 400,-1.0,1.0);
        hist_dir_xy_tdiff_tcut_test[i] = new TH2F(Form("hist_dir_xy_time_time_tcut_test_%d",i),Form("hist_dir_xy_time_time_tcut_test_%d;dir y component;dir x component",i),  400,-1.0,1.0, 400,-1.0,1.0);
        
        
        hist_cherenkov_phs_tcut_test[i] = new TH1F(Form("hist_cherenkov_phs_tcut_test_%d",i),Form("hist_cherenkov_phs_tcut_test_%d;#theta_{C} [rad];entries [#]",i),150,0.6,1);
        
        
    }
    
}

// -----   Destructor   ----------------------------------------------------
PrtLutReco::~PrtLutReco() {
    
}

Int_t mcpdata[15][65];
Int_t cluster[15][65];
Int_t lneighbours[65];
Int_t lsize(0);

Int_t getneighbours(Int_t m, Int_t p) {
    for(Int_t i=0; i<65; i++) if(p==lneighbours[i]) return -1;
    lneighbours[lsize]=p;
    lsize++;
    for(Int_t t=0; t<65; t++) {
        if(mcpdata[m][t]) {
            for(Int_t i=0; i<65; i++) if(t==lneighbours[i]) continue;
            if((t==p-1 && p%8!=0) || (t==p+1 && p%8!=7) ||
               (t==p+8 && p<57) || (t==p-8 && p>8)) getneighbours(m,t);
        }
    }
    return lsize;
}

void getclusters() {
    for(Int_t m=0; m<prt_nmcp; m++) {
        for(Int_t p=0; p<65; p++) {
            if(mcpdata[m][p])  cluster[m][p] = getneighbours(m,p);
            lsize=0;
            for(Int_t i=0; i<65; i++) lneighbours[i]=0;
        }
    }
}

//-------------- Loop over tracks ------------------------------------------
void PrtLutReco::Run(Int_t start, Int_t end) {
    
    ///////////////
    // pdf histo///
    ///////////////
    TH1F*  histo_t_pdf_p_read[12][64], *histo_t_pdf_pi_read[12][64];
    TH1F*  histo_t_pdf_p[12][64], *histo_t_pdf_pi[12][64];
    
    TFile *ffile_cherenkov_pdf_t;
    TString cherenkov_pdf_t_path;
    
    for(Int_t m=0; m<12; m++) {
        for(Int_t p=0; p<64; p++) {
            histo_t_pdf_p[m][p] = new TH1F(Form("histo_t_pdf_p_%dm_%dp",m,p),Form("histo_t_pdf_p_%dm_%dp;#theta_{C} [rad];entries [#]",m,p),   15,0,50); // 250
            histo_t_pdf_pi[m][p] = new TH1F(Form("histo_t_pdf_pi_%dm_%dp",m,p),Form("histo_t_pdf_pi_%dm_%dp;#theta_{C} [rad];entries [#]",m,p),   15,0,50); // 250
        }
    }
    
    
    cherenkov_pdf_t_path ="/Users/ahmed/final_phs/build/outFile_pdf_t.root";
    cout<<"cherenkov_pdf_t_path= " <<cherenkov_pdf_t_path<<endl;
    ffile_cherenkov_pdf_t  = new TFile(cherenkov_pdf_t_path, "READ");
    for(Int_t m=0; m<12; m++) {
        for(Int_t p=0; p<64; p++) {
            histo_t_pdf_p_read[m][p] = (TH1F*)ffile_cherenkov_pdf_t->Get(Form("histo_t_pdf_p_%dm_%dp",m,p));
            histo_t_pdf_pi_read[m][p] = (TH1F*)ffile_cherenkov_pdf_t->Get(Form("histo_t_pdf_pi_%dm_%dp",m,p));
        }
    }
    
    
    
    
    TVector3 direction, direction2  ;
    Double_t time_phs, pos_x, pos_y, pos_z ;
    TVector3 momInBar_phs(0,0,-1), momInBar_phs_0bar(0,0,-1);
    Double_t ref_point(9.8);
    Int_t sensorId_phs(-1);
    Int_t size_phs(-1);
    Double_t dir_phs_x(-10),dir_phs_y(-10),dir_phs_z(-10);
    
    Int_t kt(-1), ka(-1), count_hit_histo(0);
    Double_t kt_center(-1), ka_center(-1), average_bin(0), content_hist_dir_xy(0),content_hist_dir_xy_time(0),average_bin_time(0), content_hist_dir_xy_test(0);
    
    
    Int_t kt_pos(-1), ka_pos(-1);
    Double_t average_bin_pos(0), content_hist_pos_phs_xy(0), content_hist_pos_phs_xy_test(0);
    Double_t tcut_test;
    
    TVector3 dird,dird_phs, dird_phs_0bar, dir,dir_0bar, dir_norm, momInBar(0,0,1),posInBar,cz;
    Double_t mom, cangle,spr,tangle,tangle_0bar,tangle_phs,likelihood(0),boxPhi,weight,evtime,bartime,bartime_0bar, lenz,dirz,luttheta,luttheta_0bar, barHitTime, hitTime;
    Int_t  tofPid(0),distPid(0),likePid(0),pdgcode, evpointcount=0;
    Bool_t reflected = kFALSE;
    gStyle->SetOptFit(111);
    
    TVector3 fnX1 = TVector3 (1,0,0);
    TVector3 fnY1 = TVector3( 0,1,0);
    bool testTrRes = false;
    Double_t angdiv,dtheta,dtphi,prtangle;
    
    TString outFile = PrtManager::Instance()->GetOutName()+"_spr.root";
    Double_t theta(0),prtphi(0), trr(0),  nph(0),nph_err(0),
    par1(0), par2(0), par3(0), par4(0), par5(0), par6(0), test1(0), test2(0), test3(0),
    separation(0),separation_phs(0),separation_pdf(0),beamx(0),beamz(0),nnratio(0),nnratio_p(0),nnratio_pi(0);
    Double_t minChangle(0);
    Double_t maxChangle(1);
    Double_t deg = TMath::Pi()/180.;
    Double_t criticalAngle = asin(1.00028/1.47125); // n_quarzt = 1.47125; //(1.47125 <==> 390nm)
    
    prt_setRootPalette(1);
    prt_createMap();
    prt_initDigi();
    
    TFile file(outFile,"recreate");
    TTree tree("dirc","SPR");
    tree.Branch("mom", &mom,"mom/D");
    tree.Branch("tofPid", &tofPid,"tofPid/I");
    tree.Branch("distPid", &distPid,"distPid/I");
    tree.Branch("likePid", &likePid,"likePid/I");
    tree.Branch("spr", &spr,"spr/D");
    tree.Branch("trr", &trr,"trr/D");
    tree.Branch("nph",&nph,"nph/D");
    tree.Branch("nph_err",&nph_err,"nph_err/D");
    tree.Branch("cangle",&cangle,"cangle/D");
    tree.Branch("likelihood",&likelihood,"par3/D");
    tree.Branch("separation",&separation,"separation/D");
    tree.Branch("par5",&par5,"par5/D");
    tree.Branch("par6",&par6,"par6/D");
    tree.Branch("test1",&test1,"test1/D");
    tree.Branch("test2",&test2,"test2/D");
    tree.Branch("test3",&test3,"test3/D");
    tree.Branch("nnratio",&nnratio,"nnratio/D");
    tree.Branch("nnratio_p",&nnratio_p,"nnratio_p/D");
    tree.Branch("nnratio_pi",&nnratio_pi,"nnratio_pi/D");
    tree.Branch("theta",&theta,"theta/D");
    tree.Branch("beamx",&beamx,"beamx/D");
    tree.Branch("beamz",&beamz,"beamz/D");
    tree.Branch("prtphi",&prtphi,"prtphi/D");
    
    test1 = PrtManager::Instance()->GetTest1();
    test2 = PrtManager::Instance()->GetTest2();
    test3 = PrtManager::Instance()->GetTest3();
    beamx = PrtManager::Instance()->GetBeamX();
    beamz = PrtManager::Instance()->GetBeamZ();
    par5 = PrtManager::Instance()->GetPrismStepX();
    par6 = PrtManager::Instance()->GetPrismStepY();
    
    Double_t timeRes = PrtManager::Instance()->GetTimeRes();
    fMethod = PrtManager::Instance()->GetRunType();
    
    Int_t nEvents = fChain->GetEntries();
    if(end==0) end = nEvents;
    
    std::cout<<"Run started for ["<<start<<","<<end <<"]"<<std::endl;
    Int_t nsHits(0),nsEvents(0),studyId(0), nHits(0), ninfit(1);
    
    if(start<0) {
        ninfit=abs(start);
        start=0;
    }
    
    for (Int_t ievent=start; ievent<start+end; ievent++) { //&& ievent<end
        Int_t nhhits(0);
        fChain->GetEntry(ievent);
        nHits = fEvent->GetHitSize();
        if(ievent%1000==0) std::cout<<"Event # "<< ievent << " has "<< nHits <<" hits"<<std::endl;
        
        if(ievent-start==0) {
            tree.SetTitle(fEvent->PrintInfo());
            prtangle =fEvent->GetAngle();// prt_data_info.getAngle();// fEvent->GetAngle();
            prtphi = fEvent->GetPhi(); //prt_data_info.getPhi(); //fEvent->GetPhi();
            studyId = fEvent->GetGeometry();
            mom=fEvent->GetMomentum().Mag();
            std::cout<<"prtangle++  "<<prtangle<< " phi "<<prtphi<<std::endl;
            
            if(fEvent->GetType()==0) {
                momInBar.RotateY(TMath::Pi()-prtangle*deg);
                momInBar.RotateZ(prtphi*deg);
            } else {
                momInBar.RotateY(TMath::Pi()-prtangle*deg);
                momInBar.RotateZ(prtphi*deg);
                
                momInBar_phs_0bar.RotateY(TMath::Pi()-prtangle*deg);
                momInBar_phs_0bar.RotateZ(prtphi*deg);
            }
            
            if(fVerbose==3) {
                cz = momInBar.Unit();
                cz = TVector3(-cz.X(),cz.Y(),cz.Z());
            }
        }
        Double_t momentum=fEvent->GetMomentum().Mag();
        if( fEvent->GetType()==1) momentum /= 1000;
        tofPid=fEvent->GetParticle();
        
        Int_t pdg[]= {11,13,211,321,2212};
        Double_t mass[] = {0.000511,0.1056584,0.139570,0.49368,0.9382723};
        Double_t angle1(0), angle2(0),sum1(0),sum2(0),sum1_phs(0),sum2_phs(0),sum1_pdf(0),sum2_pdf(0), sigma(0.009),range(5*sigma),noise(0.3);
        
        fAngleP = acos(sqrt(momentum*momentum+ mass[4]*mass[4])/momentum/1.4738)-0.00; //1.4738 = 370 = 3.35
        fAnglePi= acos(sqrt(momentum*momentum + mass[2]*mass[2])/momentum/1.4738)-0.00; //-0.0014 for 160 25deg
        
        gF1->SetParameter(0,1);
        gF2->SetParameter(0,1);
        
        gF1->SetParameter(1,fAngleP);
        gF2->SetParameter(1,fAnglePi);
        
        gF1->SetParameter(2,sigma);
        gF2->SetParameter(2,sigma);
        
        
        //if(fMethod==2 && tofPid!=2212) continue;
        
        if(fEvent->GetType()==0) {
            
            Int_t gch, ndirc(0), t2(0), t3h(0), t3v(0);
            Int_t hodo1(0), hodo2(0);
            for(auto h=0; h<nHits; h++) {
                gch = fEvent->GetHit(h).GetChannel();
                if(gch<prt_maxdircch) ndirc++;
                
                if(gch==817) t2++;
                if(gch==818) t3h++;
                if(gch==819) t3v++;
                if(gch>=1350 && gch<=1351) hodo1++;
                //if(gch>=1350 && gch<=1351) hodo1++;
                //if(gch>=1369 && gch<=1370)
                hodo2++;
            }
            if(ndirc<5) continue;
            if(!(t2 && t3h && t3v && hodo1 && hodo2)) continue;
            // if(!(t3h && t3v)) continue;
        }
        
        
        //   //clusters search
        //   for(Int_t h=0; h<nHits; h++) {
        //     Int_t mid=fEvent->GetHit(h).GetMcpId();
        //     Int_t pid=fEvent->GetHit(h).GetPixelId()-1;
        //     mcpdata[mid][pid]=1;
        //   }
        //   getclusters();
        
        //   Int_t bad(0);
        //   for(Int_t j=0; j<prt_nmcp; j++){
        //     for(Int_t i=0; i<65; i++){
        //   	if(cluster[j][i]>7){
        // 	  bad=cluster[j][i];
        // 	  goto lllll;
        // 	}
        //     }
        //   }
        
        // lllll:
        //   std::cout<<"bad  "<<bad <<std::endl;
        
        //   if(bad){
        //     for(Int_t j=0; j<prt_nmcp; j++){
        //   	for(Int_t i=0; i<65; i++){
        //   	  mcpdata[j][i]=0;
        //   	  cluster[j][i]=0;
        //   	}
        //     }
        
        //     continue;
        //   }
        
        for(Int_t h=0; h<nHits; h++) {
            fHit = fEvent->GetHit(h);
            hitTime = fHit.GetLeadTime();
            fHistPhotonEnergy->Fill(fHit.GetMomentum().Mag());
            
            
            hit_b4->Fill(hitTime);
            if(fEvent->GetType()!=0) hitTime+=fRand.Gaus(0,test1); // 0.2time resol. in case it was not simulated
            

            
            //======================================== dynamic cuts
            {
                Double_t cut1(7);
                {   //time cuts
                    if(prtangle<=80) {
                        if(hitTime<cut1 || hitTime>45) continue;
                        reflected = kTRUE;
                    } else if(prtangle>94) {
                        if(hitTime<3 || hitTime>20) continue;
                        reflected = kFALSE;
                    } else {
                        if(hitTime<14)  reflected = kFALSE; //13.5
                        else reflected = kTRUE;
                    }
                }
            }
            //==================================================
            Double_t radiatorL = (fRadiator==2)? 1224.9 : 1200; //plate : bar
            
            if(fEvent->GetType()==1) lenz = radiatorL/2.-fHit.GetPosition().Z();
            else lenz = fHit.GetPosition().Z();
            
            if(fVerbose==3) {
                TVector3 cd = fHit.GetMomentum();
                fHist5->Fill(cd.Theta()*TMath::Sin(cd.Phi()),cd.Theta()*TMath::Cos(cd.Phi()));
            }
            
            // TVector3 vv = fHit.GetMomentum();
            // vv.RotateY(prtangle*deg);
            // dirz = vv.Z();
            // if(dirz<0) reflected = kTRUE;
            // else reflected = kFALSE;
            
            Int_t pixid=fHit.GetPixelId()-1;
            Int_t mcpid=fHit.GetMcpId();
            if(reflected) lenz = 2*radiatorL - lenz;
            Int_t ch = map_mpc[mcpid][pixid];
            if(prt_isBadChannel(ch)) continue;
            
            Int_t kp = histo_t_pdf_p_read[mcpid][pixid]->GetXaxis()->FindBin(hitTime);
            Int_t kpi = histo_t_pdf_pi_read[mcpid][pixid]->GetXaxis()->FindBin(hitTime);
            sum1_pdf += TMath::Log(histo_t_pdf_p_read[mcpid][pixid]->GetBinContent(kp));
            sum2_pdf += TMath::Log(histo_t_pdf_pi_read[mcpid][pixid]->GetBinContent(kpi));
            //if(cluster[mcpid][pixid]>8) continue;
            
            // Int_t x(0),y(0), piid(pixid) , nedge(0); //new
            // for(Int_t h=0; h<nHits; h++) {
            // 	Int_t pid=fEvent->GetHit(h).GetPixelId();
            // 	Int_t mid=fEvent->GetHit(h).GetMcpId();
            // 	Double_t tdif=fabs(hitTime-fEvent->GetHit(h).GetLeadTime());
            // 	if(mid!=mcpid || pid==piid || tdif>0.3) continue;
            // 	if(pid==piid-1 && piid%8!=0) y-=1;
            // 	if(pid==piid+1 && piid%8!=7) y+=1;
            
            // 	if(pid==piid+8 && piid<57) x-=1;
            // 	if(pid==piid-8 && piid>8)  x+=1;
            // }
            
            Int_t x(0),y(0), piid(pixid+1) , nedge(0); //old
            for(Int_t h=0; h<nHits; h++) {
                Int_t pid=fEvent->GetHit(h).GetPixelId();
                Int_t mid=fEvent->GetHit(h).GetMcpId();
                if(mid!=mcpid || pid==piid) continue;
                if(pid==piid-1 && piid%8!=1) x-=1;
                if(pid==piid+1 && piid%8!=0) x+=1;
                
                if(pid==piid+8 && piid<57) y+=1;
                if(pid==piid-8 && piid>8)  y-=1;
            }
            
            if(x== 0 && y== 0) nedge=0;
            if(x==-1 && y== 0) nedge=1;
            if(x==-1 && y== 1) nedge=2;
            if(x== 0 && y== 1) nedge=3;
            if(x== 1 && y== 1) nedge=4;
            if(x== 1 && y== 0) nedge=5;
            if(x== 1 && y==-1) nedge=6;
            if(x== 0 && y==-1) nedge=7;
            if(x==-1 && y==-1) nedge=8;
            
            //std::cout<< pixid << " nedge "<<nedge <<" x " <<x << "  y  "<<y<<std::endl;
            
            Int_t sensorId = 100*mcpid+fHit.GetPixelId();
            if(sensorId==1) continue;
            
            Bool_t isGoodHit(0);
            
            if(false) {
                
                
                Bool_t bool_inbar(0);
                
                //ref_point =9.8;
                ref_point =0;
                sensorId_phs = 100*mcpid+pixid;
                size_phs =fLutNode_phs[sensorId_phs]->Entries();
                for(Int_t i=0; i<1000; i++) { // size_phs  5000
                    weight = 1; //fLutNode[sensorId]->GetWeight(i);
                    dird_phs   = fLutNode_phs[sensorId_phs]->GetEntry(i).Unit();
                    //dird_phs   = fLutNode_phs[sensorId_phs]->GetEntryCs(i,nedge); // nedge=0
                    
                    
                    
                    time_phs = fLutNode_phs[sensorId_phs]->GetTime(i) ;//
                    dir_phs_x=dird_phs.X();
                    dir_phs_y=dird_phs.Y();
                    dir_phs_z=dird_phs.Z();
                    
                    
                    if (dir_phs_z > 0 && !bool_inbar ) continue;
                    if (dir_phs_z < 0 &&  bool_inbar ) continue;
                    
                    //hist_time_angle_all_phs->Fill(time_phs,momInBar_phs.Angle(dird_phs));
                    
                    for(int u=0; u<4; u++) {
                        // if((pathid==190000 || pathid==210000) && u == 0) continue; //one from left-right
                        // if((pathid==290000 || pathid==310000) && u == 0) continue; //two from left-right
                        // if((pathid==130000 || pathid==199000) && u == 0) continue; //from up-bottom
                        
                        if(u == 0) dir_0bar = dird_phs;
                        if(u == 1) dir_0bar.SetXYZ( -dird_phs.X(), dird_phs.Y(), dird_phs.Z());
                        if(u == 2) dir_0bar.SetXYZ(  dird_phs.X(),-dird_phs.Y(), dird_phs.Z()); //no need when no divergence in vertical plane
                        if(u == 3) dir_0bar.SetXYZ( -dird_phs.X(),-dird_phs.Y(), dird_phs.Z()); //no need when no divergence in vertical plane
                        if(reflected) dir_0bar.SetXYZ( dir_0bar.X(), dir_0bar.Y(), -dir_0bar.Z());
                        if(dir_0bar.Angle(fnX1) < criticalAngle || dir_0bar.Angle(fnY1) < criticalAngle) continue;
                        //hist_time_angle_all_phs->Fill(time_phs,momInBar_phs_0bar.Angle(dir_0bar));
                        
                        
                        if (bool_inbar)dir_0bar.SetXYZ( dir_0bar.X(), dir_0bar.Y(), -dir_0bar.Z());
                        
                        luttheta_0bar = dir_0bar.Theta();
                        if(luttheta_0bar > TMath::PiOver2()) luttheta_0bar = TMath::Pi()-luttheta_0bar;
                        
                        bartime_0bar = fabs(lenz/cos(luttheta_0bar)/198.);
                        double totaltime_0bar = bartime_0bar+time_phs; // phs at the end of the bar
                        if (bool_inbar)totaltime_0bar = time_phs-0.8168;
                        
                        
                        double timediff_0bar = totaltime_0bar-hitTime;
                        
                        tangle_0bar = momInBar_phs_0bar.Angle(dir_0bar);
                        
                        
                        if(fabs(tangle_0bar-fAngleP)<0.04) {
                            hist_total_time_0bar->Fill(totaltime_0bar);
                            hist_bar_time_0bar->Fill(bartime_0bar);
                            hist_hit_time_0bar->Fill(hitTime);
                            hist_timeDiff_angle_cut_0bar->Fill(timediff_0bar);
                            
                            hist_hit_phs_0bar->Fill(time_phs);
                        }
                        
                        hist_time_angle_all_0bar->Fill(totaltime_0bar,tangle_0bar);
                        hist_time_angle_all_correlation_0bar->Fill(timediff_0bar, tangle_0bar);
                        hist_time_angle_all_correlation_0bar_zoom->Fill(timediff_0bar, tangle_0bar);
                        hist_timeDiff_0bar->Fill(timediff_0bar);
                        
                        if(fabs(timediff_0bar)>test2) continue;
                        hist_cherenkov_0bar->Fill(tangle_0bar ,weight);
                    }
                    
                }// end of phs lut loop
            }
            

            
            if(true) {
                //ref_point =9.8;
                ref_point =0;
                sensorId_phs = 100*mcpid+pixid;
                size_phs =fLutNode_phs[sensorId_phs]->Entries();
                Int_t phs_solution_counter=0;
                for(Int_t i=0; i<size_phs; i++) { // size_phs  5000
                    weight = 1; //fLutNode[sensorId]->GetWeight(i);
                    dird_phs   = fLutNode_phs[sensorId_phs]->GetEntry(i).Unit();
                    //dird_phs.Rotate(TMath::Pi()/2,TVector3(0,0,1));
                    
                    
                    time_phs = fLutNode_phs[sensorId_phs]->GetTime(i);//
                    dir_phs_x=dird_phs.X();
                    dir_phs_y=dird_phs.Y();
                    dir_phs_z=dird_phs.Z();
                    
                    
                    
                    
                    //if (dir_phs_x > 0 ) continue; // for 20 deg
                    //if (time_phs<ref_point) continue;
                    //if (dir_phs_z > 0 ) continue;
                    
                    //test
                    //if (dir_phs_y > 0 ) continue;
                    
                    
                    // test
                    //if (dir_phs_x > -0.71 ||dir_phs_x < -0.75  ) continue;
                    //if(fabs(dir_phs_y)>0.4) continue;
                    //if(fabs(dir_phs_x)>-0.7) continue;
                    
                    
                    /*
                     hist_dir_x->Fill(dir_phs_x);
                     hist_dir_y->Fill(dir_phs_y);
                     hist_dir_z->Fill(dir_phs_z);
                     */
                    
                    //if (dird_phs.X() > 0 && hitTime> ref_point) continue;
                    //if (dird_phs.X() < 0 && hitTime< ref_point) continue;
                    //prt_hdigi[mcpid]->Fill(pixid%8, pixid/8);
                    pos_x = fLutNode_phs[sensorId_phs]->GetHitPos(i).X();
                    pos_y = fLutNode_phs[sensorId_phs]->GetHitPos(i).Y();
                    
                    //if(fabs(pos_x)>2) continue;
                    //if(fabs(pos_y)>4) continue;
                    
                    
                    //hist_phs_xy->Fill(pos_x,pos_y);
                    //cout<< "@@@@@@@@@@@@@@@@@@  time_phs= "<<time_phs<<" dirx= " <<dird_phs.X()<<endl;
                    
                    
                    if (hitTime >14)momInBar_phs=TVector3(0,0,1);
                    tangle_phs = momInBar_phs.Angle(dird_phs);
                    
                    //////////////////////////////////////////////
                    
                    
                    // read pdf per pmt
                    
//                    Int_t kp = histo_t_pdf_p_read[mcpid][pixid]->GetXaxis()->FindBin(hitTime);
//                    Int_t kpi = histo_t_pdf_pi_read[mcpid][pixid]->GetXaxis()->FindBin(hitTime);
//                    sum1_pdf += TMath::Log(histo_t_pdf_p_read[mcpid][pixid]->GetBinContent(kp));
//                    sum2_pdf += TMath::Log(histo_t_pdf_pi_read[mcpid][pixid]->GetBinContent(kpi));
                    
                    // use histograms

                    
                    if(false){
                        
                        //////////////////////////////////////////////
                        
                        //fcut3->SetParameters(2,8.5,0.48);
                        //if (tangle_phs < fcut3->Eval(time_phs)) continue;
                        
                        hist_time_angle_all_phs->Fill(time_phs,tangle_phs);
                        hist_time_angle_all_correlation_phs->Fill(hitTime - time_phs, tangle_phs);
                        hist_time_angle_all_correlation_phs_zoom->Fill(hitTime - time_phs, tangle_phs);
                        /*
                         kt = hist_time_angle_all_correlation_phs->GetXaxis()->FindBin(time_phs);
                         ka = hist_time_angle_all_correlation_phs->GetYaxis()->FindBin(tangle_phs);
                         kt_center=hist_time_angle_all_correlation_phs->GetXaxis()->GetBinCenter(kt);
                         ka_center=hist_time_angle_all_correlation_phs->GetYaxis()->GetBinCenter(ka);
                         //std::cout<< "########################## bin "<<kt<<"time_phs = "<< time_phs<<" "<<kt_center<<std::endl;
                         //std::cout<< "##########################  tangle_phs = "<< tangle_phs<<" "<<ka_center<<std::endl;
                         hist_time_angle_all_correlation_phs_center->Fill(hitTime - kt_center, ka_center);
                         */
                        //hist_dir_xyz->Fill(dir_phs_x, dir_phs_y, tangle_phs );
                        //if(tangle_phs > TMath::PiOver2()) tangle_phs = TMath::Pi()-tangle_phs;
                        if(fabs(tangle_phs-fAngleP)<0.04) {
                            hist_timeDiff_phs->Fill(hitTime-time_phs);
                            hist_phsTime->Fill(time_phs);
                            hist_hitTimed->Fill(hitTime);
                            hit_phs->Fill(hitTime);
                            //std::cout<< "##########################  diff =  "<< hitTime - time_phs<<std::endl;
                        }
                        
                        
                        for(auto i=0; i<21; i++) {
                            
                            tcut_test=(Double_t)i/10;
                            //std::cout<< "##########################  tcut_test "<< tcut_test<<std::endl;
                            
                            if(fabs(hitTime - time_phs)<tcut_test ) {
                                
                                
                                //hist_dir_x_tcut_test[i]->Fill(dir_phs_x,tangle_phs);
                                //hist_dir_y_tcut_test[i]->Fill(dir_phs_y,tangle_phs);
                                
                                hist_dir_xy_tcut_test[i]->Fill(dir_phs_y, dir_phs_x,tangle_phs);
                                hist_dir_xy_time_tcut_test[i]->Fill(dir_phs_y, dir_phs_x,fabs(hitTime - time_phs));// time_phs
                                hist_dir_xy_occy_tcut_test[i]->Fill(dir_phs_y, dir_phs_x);
                                
                                kt = hist_dir_xy_tcut_test[i]->GetXaxis()->FindBin(dir_phs_x);
                                ka = hist_dir_xy_tcut_test[i]->GetYaxis()->FindBin(dir_phs_y);
                                content_hist_dir_xy=hist_dir_xy_tcut_test[i]->GetBinContent(ka,kt);
                                content_hist_dir_xy_time=hist_dir_xy_time_tcut_test[i]->GetBinContent(ka,kt);
                                content_hist_dir_xy_test=hist_dir_xy_occy_tcut_test[i]->GetBinContent(ka,kt);
                                average_bin=content_hist_dir_xy/content_hist_dir_xy_test;
                                average_bin_time=content_hist_dir_xy_time/content_hist_dir_xy_test;
                                //hist_dir_xy_angle_tcut_test[i]->SetBinContent(ka,kt,average_bin);
                                hist_dir_xy_tdiff_tcut_test[i]->SetBinContent(ka,kt,average_bin_time);
                                
                                hist_cherenkov_phs_tcut_test[i]->Fill(tangle_phs,weight);
                                
                            }
                            
                            
                        }
                        
                        
                        if(fabs(hitTime - time_phs)>test2 ) continue;
                        
                        //std::cout<< "##########################  tangle_phs "<< tangle_phs<<std::endl;
                        //fHist->Fill(tangle_phs ,weight);
                        
                        
                        
                        hist_cherenkov_phs->Fill(tangle_phs,weight);
                        
                        if (dir_phs_x < fcut->Eval(dir_phs_y) || dir_phs_x > fcut2->Eval(dir_phs_y)) {
                            
                            hist_cherenkov_phs_bg->Fill(tangle_phs,weight);
                            
                            
                        } else {
                            hist_cherenkov_phs_dir_cut->Fill(tangle_phs,weight);
                        }
                        
                        
                        hist_dir_x->Fill(dir_phs_x,tangle_phs);
                        hist_dir_y->Fill(dir_phs_y,tangle_phs);
                        hist_dir_z->Fill(dir_phs_z,tangle_phs);
                        
                        hist_dir_xy->Fill(dir_phs_y, dir_phs_x,tangle_phs);
                        hist_dir_xy_time->Fill(dir_phs_y, dir_phs_x,fabs(hitTime - time_phs));// time_phs
                        hist_dir_xy_test->Fill(dir_phs_y, dir_phs_x);
                        
                        kt = hist_dir_xy->GetXaxis()->FindBin(dir_phs_x);
                        ka = hist_dir_xy->GetYaxis()->FindBin(dir_phs_y);
                        content_hist_dir_xy=hist_dir_xy->GetBinContent(ka,kt);
                        content_hist_dir_xy_time=hist_dir_xy_time->GetBinContent(ka,kt);
                        content_hist_dir_xy_test=hist_dir_xy_test->GetBinContent(ka,kt);
                        average_bin=content_hist_dir_xy/content_hist_dir_xy_test;
                        average_bin_time=content_hist_dir_xy_time/content_hist_dir_xy_test;
                        hist_dir_xy_angle->SetBinContent(ka,kt,average_bin);
                        hist_dir_xy_time_time->SetBinContent(ka,kt,average_bin_time);
                        
                        
                        //////////////
                        
                        hist_pos_phs_x_angle->Fill(pos_x,tangle_phs);
                        hist_pos_phs_y_angle->Fill(pos_y,tangle_phs);
                        
                        hist_pos_xy->Fill(pos_x,pos_y,tangle_phs);
                        hist_pos_xy_test->Fill(pos_x,pos_y);
                        
                        kt_pos = hist_pos_xy_test->GetXaxis()->FindBin(pos_x);
                        ka_pos = hist_pos_xy_test->GetYaxis()->FindBin(pos_y);
                        
                        content_hist_pos_phs_xy=hist_pos_xy->GetBinContent(kt_pos,ka_pos);
                        content_hist_pos_phs_xy_test=hist_pos_xy_test->GetBinContent(kt_pos,ka_pos);
                        average_bin_pos=content_hist_pos_phs_xy/content_hist_pos_phs_xy_test;
                        // std::cout<< kt_pos <<" "<<ka_pos<<" "<<content_hist_pos_phs_xy<<" "<<content_hist_pos_phs_xy_test<<" "<<average_bin_pos<<std::endl;
                        hist_pos_phs_xy_angle->SetBinContent(kt_pos,ka_pos,average_bin_pos);
                        
                        // prt_hdigi[mcpid]->Fill(pixid%8, pixid/8);
                        
                        if( true && tangle_phs>0.6 && tangle_phs<1 ) {
                            sum1_phs += TMath::Log(gF1->Eval(tangle_phs)+noise);
                            sum2_phs += TMath::Log(gF2->Eval(tangle_phs)+noise);
                            
                            //sum1_phs += TMath::Log(gF1->Eval(ka_center)+noise);
                            //sum2_phs += TMath::Log(gF2->Eval(ka_center)+noise);
                        }
                        ++phs_solution_counter;
                        //std::cout<< "##########################  phs_solution_counter "<< phs_solution_counter<<std::endl;
                    }// end of phs lut loop
                    hist_phs_solution_number->Fill(phs_solution_counter);
                    //if(phs_solution_counter% 2 == 0)hist_cherenkov_phs->Fill(tangle_phs ,weight);
                    
                }
            }
            
            if(false) {
                
                Int_t size =fLutNode[sensorId]->Entries();
                Int_t lut_solution_counter=0;
                for(Int_t i=0; i<size; i++) {
                    weight = 1; //fLutNode[sensorId]->GetWeight(i);
                    //dird   = fLutNode[sensorId]->GetEntryCs(i,nedge); // nedge=0
                    dird   = fLutNode[sensorId]->GetEntry(i);
                    evtime = fLutNode[sensorId]->GetTime(i);
                    Int_t pathid = fLutNode[sensorId]->GetPathId(i);
                    Bool_t samepath(false);
                    if(pathid==fHit.GetPathInPrizm()) samepath=true;
                    //if(fLutNode[sensorId]->GetNRefl(i)!=1 ) continue;
                    //if(pathid != 130000 && pathid != 199000) continue;
                    //std::cout<<"pathid "<< pathid <<std::endl;
                    //if(!samepath) continue;
                    
                    //hist_time_angle_all_lut->Fill(evtime,momInBar.Angle(dird));
                    for(int u=0; u<4; u++) {
                        // if((pathid==190000 || pathid==210000) && u == 0) continue; //one from left-right
                        // if((pathid==290000 || pathid==310000) && u == 0) continue; //two from left-right
                        // if((pathid==130000 || pathid==199000) && u == 0) continue; //from up-bottom
                        
                        if(u == 0) dir = dird;
                        if(u == 1) dir.SetXYZ( -dird.X(), dird.Y(), dird.Z());
                        if(u == 2) dir.SetXYZ( dird.X(),-dird.Y(),  dird.Z()); //no need when no divergence in vertical plane
                        if(u == 3) dir.SetXYZ( -dird.X(),-dird.Y(), dird.Z()); //no need when no divergence in vertical plane
                        if(reflected) dir.SetXYZ( dir.X(), dir.Y(), -dir.Z());
                        if(dir.Angle(fnX1) < criticalAngle || dir.Angle(fnY1) < criticalAngle) continue;
                        
                        //hist_time_angle_all_lut->Fill(evtime,momInBar.Angle(dir));
                        
                        luttheta = dir.Theta();
                        if(luttheta > TMath::PiOver2()) luttheta = TMath::Pi()-luttheta;
                        
                        bartime = fabs(lenz/cos(luttheta)/198.);
                        double totaltime = bartime+evtime;
                        double timediff = totaltime-hitTime;
                        fHist0->Fill(totaltime-hitTime);
                        if(samepath)  fHist0i->Fill(timediff);
                        //	  fHist1->Fill(hitTime);
                        fHist2->Fill(totaltime);
                        tangle = momInBar.Angle(dir);
                        
                        //dir_norm=dir.Unit();
                        //tangle = momInBar.Angle(dir_norm.X());
                        
                        hist_time_angle_all_lut->Fill(totaltime,tangle);
                        hist_time_angle_all_correlation_lut->Fill(timediff, tangle);
                        hist_time_angle_all_correlation_lut_zoom->Fill(timediff, tangle);
                        
                        
                        
                        if(fabs(tangle-fAngleP)<0.04) {
                            hist_timeDiff_lut->Fill(timediff);
                            
                            //std::cout<< "##########################  diff =  "<< hitTime - time_phs<<std::endl;
                        }
                        
                        if(fabs(timediff)>test3) continue;
                        
                        fHist3->Fill(fabs(totaltime),hitTime);
                        //tangle = momInBar.Angle(dir)-0.002;
                        
                        
                        
                        if(fabs(tangle-fAngleP)<0.04) {
                            hit_lut->Fill(hitTime);
                        }
                        
                        if(tangle > minChangle && tangle < maxChangle && tangle < 1.85) {
                            if(tofPid==211 && fMethod==2) fHistPi->Fill(tangle ,weight);
                            else {
                                fHist->Fill(tangle ,weight);
                                
                                
                            }
                            hist_cherenkov_lut->Fill(tangle ,weight);
                            hist_dir_x_lut->Fill(dir_norm.X(),tangle);
                            hist_dir_y_lut->Fill(dir_norm.Y(),tangle);
                            hist_dir_z_lut->Fill(dir_norm.Z(),tangle);
                            hist_dir_xy_lut->Fill(dir_norm.Y(), dir_norm.X());
                            
                            
                            
                            if(tofPid==2212) fHistMcp[mcpid]->Fill(tangle ,weight);
                            fHistCh[ch]->Fill(tangle ,weight);
                            
                            if(true && tangle>0.6 && tangle<1 ) {
                                sum1 += TMath::Log(gF1->Eval(tangle)+noise);
                                sum2 += TMath::Log(gF2->Eval(tangle)+noise);
                            }
                            
                            // //if(samepath) fHist->Fill(tangle ,weight);
                            if(fRadiator==1 && fabs(tangle-0.815)<0.05) isGoodHit=true;
                            if(fRadiator==2 && fabs(tangle-0.815)<0.2)  isGoodHit=true;
                            
                            if(fVerbose==3) {
                                TVector3 rdir = TVector3(-dir.X(),dir.Y(),dir.Z());
                                rdir.RotateUz(cz);
                                Double_t phi = rdir.Phi();
                                Double_t tt =  rdir.Theta();
                                fHist4->Fill(tt*TMath::Sin(phi),tt*TMath::Cos(phi));
                                
                                //for cherenckov circle fit
                                gg_gr.SetPoint(gg_i,tt*TMath::Sin(phi),tt*TMath::Cos(phi));
                                gg_i++;
                            }
                        }
                        ++lut_solution_counter;
                    }
                    
                    
                }// end of lut loop
                hist_lut_solution_number->Fill(lut_solution_counter);
            }
            fHist1->Fill(hitTime);
            if(isGoodHit) {
                fHist6->Fill(hitTime);
                nhhits++;
                nsHits++;
                prt_hdigi[mcpid]->Fill(pixid%8, pixid/8);
            }
        }
        
        if(nhhits>35)
            fhNph->Fill(nhhits);
        
        
        if(tofPid==2212) fhNph_p->Fill(nhhits);
        if(tofPid==211) fhNph_pi->Fill(nhhits);
        
        // for(Int_t j=0; j<prt_nmcp; j++){
        //   for(Int_t i=0; i<65; i++){
        // 	mcpdata[j][i]=0;
        // 	cluster[j][i]=0;
        //   }
        // }
        
        Double_t sum = sum1-sum2;
        if(sum!=0) {
            if(tofPid==2212) hLnDiffP->Fill(sum);
            if(tofPid==211) hLnDiffPi->Fill(sum);//211
            
            /// std::cout<< " tofPid "<<tofPid <<" sum "<<   sum  <<std::endl;
            //likelihood=sum;
        }
        
        Double_t sum_phs = sum1_phs-sum2_phs;
        if(sum_phs!=0) {
            if(tofPid==2212) hLnDiffP_phs->Fill(sum_phs);
            if(tofPid==211) hLnDiffPi_phs->Fill(sum_phs);//211
            
            //std::cout<< " sum_phs "<<sum_phs <<std::endl;
            
            //likelihood=sum_phs;
        }
        
        Double_t sum_pdf = sum1_pdf-sum2_pdf;
        if(sum_pdf!=0) {
            if(tofPid==2212) {
                hLnDiffP_pdf->Fill(sum_pdf);
                //std::cout<< " ####### hLnDiffP_pdf "<<hLnDiffP_pdf->GetEntries() <<std::endl;
                
            }
            if(tofPid==211){
                hLnDiffPi_pdf->Fill(sum_pdf);
                //std::cout<< " ####### hLnDiffPi_pdf "<<hLnDiffPi_pdf->GetEntries() <<std::endl;
                
            }
            
            //std::cout<< " ####### sum_pdf "<<sum_pdf <<std::endl;
            
            likelihood=sum_pdf;
        }
        
        // if(fVerbose==1){
        //   prt_canvasAdd("ff",800,400);
        //   gF1->Draw();
        //   gF2->SetLineColor(4);
        //   gF2->Draw("same");
        
        //   prt_waitPrimitive("ff");
        //   prt_canvasDel("ff");
        //   //prt_canvasSave(1,0);
        //   //prt_canvasDel(Form("lh_%d",gg_ind));
        // }
        
        if(fVerbose>0 &&  fMethod==3 && nsEvents%ninfit==0) {
            if(nsHits>10) {
                // if(tofPid==2212 && sum > 0){
                //   std::cout<<"p  "<<sum1 << "   pi "<<sum2 << "  s "<< sum<<std::endl;
                //   if(fVerbose>0)  if(!FindPeak(cangle,spr, prtangle, tofPid)) continue;
                // }
                
                FindPeak(cangle,spr, prtangle, tofPid);
                distPid = FindPdg(momentum,cangle);
                nph = nsHits/(Double_t)ninfit;
                spr = spr*1000;
                trr = spr/sqrt(nph);
                theta = fEvent->GetAngle();
                par3 = fEvent->GetTest1();
                tree.Fill();
            }
            ResetHists();
            nsHits=0;
        }
        
        if(++nsEvents>=end) break;
    }
    
    nnratio = fhNph->GetEntries()/(double)end;
    nnratio_pi = fhNph_pi->GetEntries()/(double)end;
    nnratio_p = fhNph_p->GetEntries()/(double)end;
    
    std::cout<<"nnratio "<<nnratio<<" "<<end <<"  "<< fhNph->GetEntries()<<std::endl;
    
    TF1 *ff;
    if(fMethod==2) {
        gROOT->SetBatch(1);
        if(fhNph->GetEntries()>20) {
            fhNph->Fit("gaus","","MQN",20,120);
            ff = fhNph->GetFunction("gaus");
            nph=ff->GetParameter(1);
            nph_err=ff->GetParError(1);
        }
        //nph = prt_fit(fhNph,40,10,50,1).X();
        gROOT->SetBatch(0);
        FindPeak(cangle,spr, prtangle);
        //nph = nsHits/(Double_t)nsEvents;
        spr = spr*1000;
        trr = spr/sqrt(nph);
        theta = fEvent->GetAngle();
        par3 = fEvent->GetTest1();
        if(fVerbose) std::cout<<Form("SPR=%2.2F N=%2.2f +/- %2.2f",spr,nph,nph_err)<<std::endl;
        tree.Fill();
    } else {
        if(fVerbose<2) gROOT->SetBatch(1);
        prt_canvasAdd("r_lhood",800,400);
        prt_normalize(hLnDiffP,hLnDiffPi);
        hLnDiffP->SetLineColor(2);
        
        Double_t m1,m2,s1,s2;
        if(hLnDiffP->GetEntries()>10) {
            hLnDiffP->Fit("gaus","S");
            ff = hLnDiffP->GetFunction("gaus");
            m1=ff->GetParameter(1);
            s1=ff->GetParameter(2);
        }
        if(hLnDiffPi->GetEntries()>10) {
            hLnDiffPi->Fit("gaus","S");
            ff = hLnDiffPi->GetFunction("gaus");
            m2=ff->GetParameter(1);
            s2=ff->GetParameter(2);
        }
        separation = (fabs(m2-m1))/(0.5*(s1+s2));
        std::cout<<"######### separation "<< separation <<std::endl;
        
        //gStyle->SetOptFit(0);
        //gStyle->SetOptStat(0);
        
        hLnDiffP->SetName(Form("s_%2.2f",separation));
        hLnDiffP->Draw();
        hLnDiffPi->SetLineColor(4);
        hLnDiffPi->Draw("same");
        
        
        ////////
        prt_canvasAdd("r_lhood_phs",800,400);
        prt_normalize(hLnDiffP_phs,hLnDiffPi_phs);
        hLnDiffP_phs->SetLineColor(2);
        
        Double_t m1_phs,m2_phs,s1_phs,s2_phs;
        if(hLnDiffP_phs->GetEntries()>10) {
            hLnDiffP_phs->Fit("gaus","S");
            ff = hLnDiffP_phs->GetFunction("gaus");
            m1_phs=ff->GetParameter(1);
            s1_phs=ff->GetParameter(2);
        }
        if(hLnDiffPi_phs->GetEntries()>10) {
            hLnDiffPi_phs->Fit("gaus","S");
            ff = hLnDiffPi_phs->GetFunction("gaus");
            m2_phs=ff->GetParameter(1);
            s2_phs=ff->GetParameter(2);
        }
        separation_phs = (fabs(m2_phs-m1_phs))/(0.5*(s1_phs+s2_phs));
        std::cout<<"######### separation_phs "<< separation_phs <<std::endl;
        
        //gStyle->SetOptFit(0);
        //gStyle->SetOptStat(0);
        
        hLnDiffP_phs->SetName(Form("s_%2.2f",separation_phs));
        hLnDiffP_phs->Draw();
        hLnDiffPi_phs->SetLineColor(4);
        hLnDiffPi_phs->Draw("same");
        
        
        ////////
        prt_canvasAdd("r_lhood_pdf",800,400);
        prt_normalize(hLnDiffP_pdf,hLnDiffPi_pdf);
        hLnDiffP_pdf->SetLineColor(2);
        
        Double_t m1_pdf,m2_pdf,s1_pdf,s2_pdf;
        if(hLnDiffP_pdf->GetEntries()>10) {
            hLnDiffP_pdf->Fit("gaus","S");
            ff = hLnDiffP_pdf->GetFunction("gaus");
            m1_pdf=ff->GetParameter(1);
            s1_pdf=ff->GetParameter(2);
        }
        if(hLnDiffPi_pdf->GetEntries()>10) {
            hLnDiffPi_pdf->Fit("gaus","S");
            ff = hLnDiffPi_pdf->GetFunction("gaus");
            m2_pdf=ff->GetParameter(1);
            s2_pdf=ff->GetParameter(2);
        }
        separation_pdf = (fabs(m2_pdf-m1_pdf))/(0.5*(s1_pdf+s2_pdf));
        std::cout<<"######### separation_pdf "<< separation_pdf <<std::endl;
        
        //gStyle->SetOptFit(0);
        //gStyle->SetOptStat(0);
        
        hLnDiffP_pdf->SetName(Form("s_%2.2f",separation_pdf));
        hLnDiffP_pdf->Draw();
        hLnDiffPi_pdf->SetLineColor(4);
        hLnDiffPi_pdf->Draw("same");
        
        
        prt_canvasSave(1,0);
        //prt_waitPrimitive("r_lhood","w");
        if(fVerbose) gROOT->SetBatch(0);
        tree.Fill();
    }
    
    if(fVerbose) ResetHists();
    
    
    
    tree.Write();
    hist_cherenkov_phs->Write();
    hist_cherenkov_phs_dir_cut->Write();
    hist_cherenkov_phs_bg->Write();
    hist_cherenkov_lut->Write();
    
    //hist_hitTimed->Write();
    hist_time_angle_all_phs->Write();
    hist_time_angle_all_lut->Write();
    
    hist_time_angle_all_correlation_phs->Write();
    hist_time_angle_all_correlation_phs_zoom->Write();
    hist_time_angle_all_correlation_lut->Write();
    //hist_time_angle_all_correlation_phs_center->Write();
    hist_time_angle_all_correlation_lut_zoom->Write();
    
    
    hist_timeDiff_phs->Write();
    hist_timeDiff_lut->Write();
    hist_phsTime->Write();
    hit_b4->Write();
    hit_phs->Write();
    hit_lut->Write();
    
    //fHistPhotonEnergy->Write();
    
    hist_dir_x->Write();
    hist_dir_y->Write();
    hist_dir_z->Write();
    hist_dir_xy->Write();
    hist_dir_xy_test->Write();
    hist_dir_xy_angle->Write();
    
    
    hist_pos_phs_x_angle->Write();
    hist_pos_phs_y_angle->Write();
    hist_pos_phs_xy_angle->Write();
    hist_pos_xy->Write();
    hist_pos_xy_test->Write();
    
    hist_dir_xy_time->Write();
    hist_dir_xy_time_time->Write();
    
    /*
     hist_dir_x_lut->Write();
     hist_dir_y_lut->Write();
     hist_dir_z_lut->Write();
     hist_dir_xy_lut->Write();
     */
    
    
    hist_phs_solution_number->Write();
    hist_lut_solution_number->Write();
    
    
    
    for(auto i=0; i<21; i++) {
        hist_cherenkov_phs_tcut_test[i]->Write();
    }
    
    for(auto i=0; i<21; i++) {
        hist_dir_xy_occy_tcut_test[i]->Write();
    }
    for(auto i=0; i<21; i++) {
        hist_dir_xy_tdiff_tcut_test[i]->Write();
    }
    
    
    if (false) {
        hLnDiffP->Write();
        hLnDiffPi->Write();
        hLnDiffP_phs->Write();
        hLnDiffPi_phs->Write();
        
        hist_total_time_0bar->Write();
        hist_hit_time_0bar->Write();
        
        hist_bar_time_0bar->Write();
        hist_hit_phs_0bar->Write();
        
        hist_time_angle_all_0bar->Write();
        hist_time_angle_all_correlation_0bar->Write();
        hist_time_angle_all_correlation_0bar_zoom->Write();
        
        hist_timeDiff_0bar->Write();
        hist_timeDiff_angle_cut_0bar->Write();
        
        hist_cherenkov_0bar->Write();
    }
    
    Double_t binContent_phs = hist_phs_solution_number->GetBinContent(1);
    Double_t Entries_phs = hist_phs_solution_number->GetEntries();
    
    Double_t binContent_lut = hist_lut_solution_number->GetBinContent(1);
    Double_t Entries_lut = hist_lut_solution_number->GetEntries();
    
    
    Double_t ratio_phs = binContent_phs*100/Entries_phs;
    Double_t ratio_lut = binContent_lut*100/Entries_lut;
    
    std::cout<<"ratio_phs "<<ratio_phs<<" "<<binContent_phs <<"  "<< binContent_phs<<std::endl;
    std::cout<<"ratio_lut "<<ratio_lut<<" "<<binContent_lut <<"  "<< binContent_lut<<std::endl;
    
    
    file.Write();
}

Int_t g_num =0;
Bool_t PrtLutReco::FindPeak(Double_t& cangle, Double_t& spr, Double_t a, Int_t tofpdg) {
    cangle=0;
    spr=0;
    //  gStyle->SetCanvasPreferGL(kTRUE);
    //if(fHist->GetEntries()>20 || fHistPi->GetEntries()>20) {
    if(false) {
        gROOT->SetBatch(1);
        Int_t nfound = fSpect->Search(fHist,1,"",0.9); //0.6
        if(nfound>0) cangle = fSpect->GetPositionX()[0];
        else cangle =  fHist->GetXaxis()->GetBinCenter(fHist->GetMaximumBin());
        cangle =  fHist->GetXaxis()->GetBinCenter(fHist->GetMaximumBin());
        
        if(cangle>0.85) cangle=0.82;
        fFit->SetParameters(100,cangle,0.010);
        fFit->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
        fFit->SetParLimits(0,0.1,1E6);
        fFit->SetParLimits(1,cangle-0.04,cangle+0.04);
        fFit->SetParLimits(2,0.005,0.016); // width 7-10
        // fFit->FixParameter(2,0.01);
        // fFit->FixParameter(3,0);
        // fFit->FixParameter(4,0);
        
        Int_t status(0);
        if(fMethod==3) status = fHist->Fit("fgaus","lq","",0.6,1);
        else status =fHist->Fit("fgaus","M","",cangle-0.06,cangle+0.06);
        Double_t chi = fFit->GetChisquare()/fFit->GetNDF();
        
        // if(fFit->GetParError(1)>0.0035){
        // //   // if(fFit->GetParameter(2)>0.011){
        // //   // if(fabs(chi-1<0.3 ){
        //   spr=0;
        //   cangle=0;
        //   fTest=0;
        //   return false;
        // }else{
        //   fTest=chi;
        // }
        
        cangle = fFit->GetParameter(1);
        spr = fFit->GetParameter(2);
        
        if(fVerbose>1) gROOT->SetBatch(0);
        
        if(fMethod==2 && fVerbose>0) {
            
            TString nid = "";//Form("_%2.0f",a);
            prt_canvasAdd("r_tangle"+nid,800,400);
            
            // fFit->SetParLimits(2,0.004,0.008); // width 7-10
            // for(Int_t i=0; i<prt_nmcp; i++){
            // 	prt_canvasAdd(Form("r_tangle_%d",i),800,400);
            // 	fHistMcp[i]->Fit("fgaus","lq","",fAngleP-0.03,fAngleP+0.03);
            // 	std::cout<<"if(mcpid=="<< i<<") tangle += "<<fAngleP-fFit->GetParameter(1)<<";" <<std::endl;
            // 	fHistMcp[i]->Draw();
            // 	drawTheoryLines();
            // }
            
            // for(Int_t i=0; i<960; i++){
            // 	prt_canvasAdd(Form("r_tangle_ch_%d",i),800,400);
            // 	fHistCh[i]->Fit("fgaus","lq","",fAngleP-0.03,fAngleP+0.03);
            // 	std::cout<<"if(ch=="<< i<<") tangle += "<<fAngleP-fFit->GetParameter(1)<<";" <<std::endl;
            // 	fHistCh[i]->Draw();
            // }
            
            //      TString name = Form("r_tangle_%3.1f",test3);
            fHist->SetTitle(Form("theta %3.1f , TOF PID = %d", a, tofpdg));
            fHist->SetMinimum(0);
            //fHist->Scale(1/fHist->GetMaximum());
            
            prt_normalize(fHist,fHistPi);
            fHistPi->SetLineColor(4);
            fHist->SetLineColor(1);
            
            fHist->Draw();
            fHistPi->Draw("same");
            // gF1->Draw("same");
            // gF2->Draw("same");
            fHisti->SetLineColor(kRed+2);
            if(fHisti->GetEntries()>5) fHisti->Draw("same");
            
            drawTheoryLines();
            
            prt_canvasAdd("r_time",800,400);
            prt_normalize(fHist1,fHist6);
            fHist1->SetTitle(Form("theta %3.1f", a));
            fHist1->SetLineColor(2);
            fHist1->Draw();
            fHist6->Draw("same");
            //fHist2->Draw("same");
            
            prt_canvasAdd("r_nph"+nid,800,400);
            fhNph->Draw();
            
            prt_canvasAdd("r_diff"+nid,800,400);
            fHist0->SetTitle(Form("theta %3.1f", a));
            fHist0->Draw();
            fHist0i->SetLineColor(kRed+2);
            if(fHist0i->GetEntries()>5)  fHist0i->Draw("same");
            
            // prt_canvasAdd("r_cm"+nid,800,400);
            // fHist3->SetTitle(Form("theta %3.1f", a));
            // fHist3->Draw("colz");
            
            if(false) {
                Int_t tmax, max=0;
                for(Int_t m=0; m<prt_nmcp; m++) {
                    prt_hdigi[m]->Rebin2D(8,8);
                    prt_hdigi[m]->GetXaxis()->SetNdivisions(0);
                    prt_hdigi[m]->GetYaxis()->SetNdivisions(0);
                    prt_hdigi[m]->GetXaxis()->SetTickLength(0);
                    prt_hdigi[m]->GetYaxis()->SetTickLength(0);
                    prt_hdigi[m]->GetXaxis()->SetAxisColor(1);
                    prt_hdigi[m]->GetYaxis()->SetAxisColor(1);
                    prt_hdigi[m]->SetMarkerSize(10);
                    tmax = prt_hdigi[m]->GetMaximum();
                    if(max<tmax) max = tmax;
                }
                for(Int_t m=0; m<prt_nmcp; m++) {
                    prt_hdigi[m]->Scale(1/(Double_t)max);
                }
            }
            
            prt_drawDigi("m,p,v\n",2017);
            prt_cdigi->SetName("r_hp"+nid);
            prt_canvasAdd(prt_cdigi);
            
            if(fVerbose>1)  prt_waitPrimitive("r_time"+nid);
            prt_canvasSave(1,0);
            prt_canvasDel("*");
            
            if(fVerbose==3) {
                TCanvas* c2 = new TCanvas("c2","c2",0,0,800,400);
                c2->Divide(2,1);
                c2->cd(1);
                
                fHist4->SetStats(0);
                fHist4->SetTitle(Form("Calculated from LUT, #theta = %3.1f#circ", a));
                fHist4->Draw("colz");
                Double_t x0(0), y0(0), theta(cangle);
                FitRing(x0,y0,theta);
                TVector3 corr(x0,y0,1-TMath::Sqrt(x0*x0+y0*y0));
                std::cout<<"Tcorr "<< corr.Theta()*1000<< "  Pcorr "<< corr.Phi() <<std::endl;
                
                TLegend *leg = new TLegend(0.5,0.7,0.85,0.87);
                //      leg->SetFillColor(0);
                //leg->SetFillColorAlpha(0,0.8);
                leg->SetFillStyle(0);
                //leg->SetFillStyle(4000);
                leg->SetBorderSize(0);
                leg->AddEntry((TObject*)0,Form("Entries %0.0f",fHist4->GetEntries()),"");
                leg->AddEntry((TObject*)0,Form("#Delta#theta_{c} %f [mrad]",corr.Theta()*1000),"");
                leg->AddEntry((TObject*)0,Form("#Delta#varphi_{c} %f [mrad]",corr.Phi()),"");
                leg->Draw();
                
                TArc *arc = new TArc(x0,y0,theta);
                arc->SetLineColor(kRed);
                arc->SetLineWidth(1);
                arc->SetFillStyle(0);
                arc->Draw();
                gg_i=0;
                gg_gr.Set(0);
                
                c2->cd(2);
                gStyle->SetOptStat(1110);
                fHist5->SetTitle(Form("True from MC, #theta = %d#circ", a));
                fHist5->Draw("colz");
                
                c2->Print(Form("spr/tcorr_%d.png", a));
                c2->Modified();
                c2->Update();
                c2->WaitPrimitive("");
            }
        }
    }
    
    if(fVerbose<2) gROOT->SetBatch(0);
    
    return (cangle>0 && cangle<1);
}

void PrtLutReco::ResetHists() {
    fHist->Reset();
    fHisti->Reset();
    fHist0->Reset();
    fHist0i->Reset();
    fHist1->Reset();
    fHist2->Reset();
    fHist3->Reset();
    fHist4->Reset();
    for(Int_t m=0; m<prt_nmcp; m++) prt_hdigi[m]->Reset();
}

TF1 *lFit = new TF1("lgaus","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.6,0.9);
TF1 *lFitPi = new TF1("lgausPi","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.6,0.9);

Double_t PrtLutReco::fillLnDiffPPi(Double_t cangle, Int_t tofPid, Double_t mom) {
    if(fHist->GetEntries()>20 ) {
        Int_t pdg[]= {11,13,211,321,2212};
        Double_t mass[] = {0.000511,0.1056584,0.139570,0.49368,0.9382723};
        
        Double_t angle1(0), angle2(0), sigma(0.006),range(0.015);
        
        // //fHist->Scale(1/fHist->GetMaximum());
        
        // Double_t d1,d2, sum1(0),sum2(0);
        // Int_t sbin = fHist->FindBin(fAngleP-range);
        // Int_t ebin = fHist->FindBin(fAngleP+range);
        // // fHist->GetXaxis()->GetNbins()
        // for(Int_t i=sbin; i< ebin; i++){
        //   if(fHist->GetBinContent(i) < 0.01 ) continue;
        //   d1 = gF1->Eval(fHist->GetBinCenter(i))- fHist->GetBinContent(i);
        //   d2 = gF1->Eval(fHist->GetBinCenter(i))- fHist->GetBinContent(i);
        
        //   std::cout<<"f1 "<< gF1->Eval(fHist->GetBinCenter(i)) << "   f2 "<<gF2->Eval(fHist->GetBinCenter(i)) << "    v "<< fHist->GetBinContent(i) <<std::endl;
        
        //   // if(d1>0) sum1+=TMath::Log(d1);
        //   // if(d2>0) sum2+=TMath::Log(d2);
        //   sum1+=TMath::Log(fabs(d1));
        //   sum2+=TMath::Log(fabs(d2));
        
        // }
        // Double_t amin(sum1),amin2(sum2);
        
        // lFit->SetRange(fAngleP-range,fAngleP+range);
        // lFit->FixParameter(0,fFit->GetParameter(0));
        // lFit->FixParameter(1,fAngleP);
        // if(fFit->GetParameter(2)>sigma) sigma=fFit->GetParameter(2);
        // lFit->FixParameter(2,sigma);
        // lFit->FixParameter(3,fFit->GetParameter(3));
        // lFit->FixParameter(4,fFit->GetParameter(4));
        
        lFit->SetRange(fAngleP-range,fAnglePi+range);
        lFit->FixParameter(0,fHist->GetMaximum()-0.5);
        lFit->FixParameter(1,fAngleP);
        lFit->FixParameter(2,0.01);
        lFit->FixParameter(3,0);
        lFit->FixParameter(4,0.5);
        
        
        fHist->Fit("lgaus","lq","",fAngleP-range,fAnglePi+range);
        Double_t amin,amin2,edm,errdef;
        Int_t nvpar,nparx;
        TVirtualFitter *fitter = TVirtualFitter::Fitter(fHist);
        fitter->GetStats(amin,edm,errdef,nvpar,nparx);
        
        // lFitPi->SetRange(fAnglePi-range,fAnglePi+range);
        // lFitPi->SetLineColor(4);
        // lFitPi->FixParameter(0,fFit->GetParameter(0));
        // lFitPi->FixParameter(1,fAnglePi);
        // lFitPi->FixParameter(2,sigma);
        // lFitPi->FixParameter(3,fFit->GetParameter(3));
        // lFitPi->FixParameter(4,fFit->GetParameter(4));
        
        lFitPi->SetRange(fAngleP-range,fAnglePi+range);
        lFitPi->SetLineColor(4);
        lFitPi->FixParameter(0,fHist->GetMaximum()-0.5);
        lFitPi->FixParameter(1,fAnglePi);
        lFitPi->FixParameter(2,0.01);
        lFitPi->FixParameter(3,0);
        lFitPi->FixParameter(4,0.5);
        
        fHist->Fit("lgausPi","lq","",fAngleP-range,fAnglePi+range);
        fitter = TVirtualFitter::Fitter(fHist);
        fitter->GetStats(amin2,edm,errdef,nvpar,nparx);
        
        if(fVerbose) printf("tofPid %04d | %1.4f (%1.4f/%1.4f) likelihood is %1.2f/%1.2f \n",tofPid,cangle,fAngleP,fAnglePi, amin, amin2);
        gg_ind++;
        
        if(fVerbose==1) {
            prt_canvasAdd("ff",800,400);
            //prt_canvasAdd(Form("lh_%d",gg_ind),800,400);
            fHist->SetTitle(Form("%d",tofPid));
            fHist->Draw();
            lFit->SetLineColor(2);
            lFit->Draw("same");
            // gF1->Draw("same");
            // gF2->SetLineColor(4);
            // gF2->Draw("same");
            
            //if(fabs(amin-amin2)<5)
            prt_waitPrimitive("ff");
            prt_canvasDel("ff");
            //prt_canvasSave(1,0);
            //prt_canvasDel(Form("lh_%d",gg_ind));
        }
        
        return amin-amin2;
    }
    return 1000;
}

Double_t PrtLutReco::fillLnDiffPPi2(Double_t cangle, Int_t tofPid, Double_t mom) {
    if(fHist->GetEntries()>20 ) {
        Int_t pdg[]= {11,13,211,321,2212};
        Double_t mass[] = {0.000511,0.1056584,0.139570,0.49368,0.9382723};
        
        Double_t angle1(0), angle2(0), sigma(0.006),range(0.03);
        
        Double_t d1,d2, sum1(0),sum2(0);
        Int_t sbin = fHist->FindBin(fAnglePi-range);
        Int_t ebin = fHist->FindBin(fAngleP+range);
        for(Int_t i=sbin; i< ebin; i++) {
            if(fHist->GetBinContent(i)<1 ) continue;
            d1 = 10*fabs(fHist->GetBinContent(i) *(fAngleP  - fHist->GetBinCenter(i)));
            d2 = 10*fabs(fHist->GetBinContent(i) *(fAnglePi - fHist->GetBinCenter(i)));
            if(d1>0 && d2>0) {
                std::cout<<"d1  "<<d1 << "   d2    "<< d2 <<std::endl;
                sum1+=TMath::Log(d1);
                sum2+=TMath::Log(d2);
            }
        }
        
        if(fVerbose) printf("tofPid %04d | %1.4f (%1.4f/%1.4f) likelihood is %1.2f/%1.2f \n",tofPid,cangle,fAngleP,fAnglePi, sum1, sum2);
        gg_ind++;
        
        if(fVerbose==1) {
            prt_canvasAdd("ff",800,400);
            //prt_canvasAdd(Form("lh_%d",gg_ind),800,400);
            fHist->SetTitle(Form("%d",tofPid));
            fHist->Draw();
            lFit->SetLineColor(2);
            lFit->Draw("same");
            // gFp->Draw("same");
            // gFpi->SetLineColor(4);
            // gFpi->Draw("same");
            
            //if(fabs(amin-amin2)<5)
            prt_waitPrimitive("ff");
            prt_canvasDel("ff");
            //prt_canvasSave(1,0);
            //prt_canvasDel(Form("lh_%d",gg_ind));
        }
        
        return sum1-sum2;
    }
    return 1000;
}

void circleFcn(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
    Int_t np = gg_gr.GetN();
    f = 0;
    Double_t *x = gg_gr.GetX();
    Double_t *y = gg_gr.GetY();
    for (Int_t i=0; i<np; i++) {
        Double_t u = x[i] + par[0];
        Double_t v = y[i] + par[1];
        Double_t dr = par[2] - TMath::Sqrt(u*u+v*v);
        f += dr*dr;
    }
    std::cout<<"fcn  "<< f<<std::endl;
    
}

void circleFcn2(Int_t &, Double_t *, Double_t &f, Double_t *par, Int_t) {
    Int_t np = gg_gr.GetN();
    f = 0;
    Double_t *x = gg_gr.GetX();
    Double_t *y = gg_gr.GetY();
    for (Int_t i=0; i<np; i++) {
        Double_t u = x[i] + par[0];
        Double_t v = y[i] + par[1];
        Double_t dr = par[2] - TMath::Sqrt(u*u+v*v);
        if(dr>0.07) f += dr*dr;
        else f += fabs(dr);
    }
}

void PrtLutReco::FitRing(Double_t& x0, Double_t& y0, Double_t& theta) {
    TGraph ff_gr;
    Int_t ff_i(0);
    Int_t np = gg_gr.GetN();
    Double_t *x = gg_gr.GetX();
    Double_t *y = gg_gr.GetY();
    for (Int_t i=0; i<np; i++) {
        if( fabs(theta - TMath::Sqrt(x[i]*x[i]+y[i]*y[i]))<0.05) {
            ff_gr.SetPoint(ff_i,x[i],y[i]);
            ff_i++;
        }
    }
    gg_gr = ff_gr;
    
    //Fit a circle to the graph points
    TVirtualFitter::SetDefaultFitter("Minuit");  //default is Minuit
    TVirtualFitter *fitter = TVirtualFitter::Fitter(0, 3);
    fitter->SetPrecision(0.00000001);
    fitter->SetMaxIterations(1000);
    
    fitter->SetFCN(circleFcn);
    fitter->SetParameter(0, "x0",   0.03, 0.01, -0.05,0.05);
    fitter->SetParameter(1, "y0",   0, 0.01, -0.05,0.05);
    fitter->SetParameter(2, "R",    theta, 0.01, theta-0.05,theta+0.05);
    
    //fitter->FixParameter(0);
    //fitter->FixParameter(1);
    fitter->FixParameter(2);
    Double_t arglist[1] = {0};
    fitter->ExecuteCommand("MINIMIZE", arglist, 0);
    
    // fitter->SetFCN(circleFcn2);
    // fitter->ExecuteCommand("MINIMIZE", arglist, 0);
    
    x0 = fitter->GetParameter(0);
    y0 = fitter->GetParameter(1);
    theta = fitter->GetParameter(2);
}

Int_t PrtLutReco::FindPdg(Double_t mom, Double_t cangle) {
    // Int_t pdg[]={11,13,211,321,2212};
    // Double_t mass[] = {0.000511,0.1056584,0.139570,0.49368,0.9382723};
    // Int_t pdg[]={211,321,2212};
    // Double_t mass[] = {0.139570,0.49368,0.9382723};
    cangle =  fHist->GetXaxis()->GetBinCenter(fHist->GetMaximumBin());
    
    Int_t pdg[]= {211,2212};
    Double_t mass[] = {0.139570,0.9382723};
    Double_t tdiff, diff=100;
    Int_t minid=0;
    for(Int_t i=0; i<2; i++) {
        tdiff = fabs(cangle - acos(sqrt(mom*mom + mass[i]*mass[i])/mom/1.46907)); //1.46907 - fused silica
        if(tdiff<diff) {
            diff = tdiff;
            minid = i;
        }
    }
    return pdg[minid];
}

void PrtLutReco::drawTheoryLines() {
    gPad->Update();
    TLine *line = new TLine(0,0,0,1000);
    line->SetX1(fAngleP);
    line->SetX2(fAngleP);
    line->SetY1(gPad->GetUymin());
    line->SetY2(gPad->GetUymax());
    line->SetLineColor(kRed);
    line->Draw();
    
    TLine *line1 = new TLine(0,0,0,1000);
    line1->SetX1(fAnglePi);
    line1->SetX2(fAnglePi);
    line1->SetY1(gPad->GetUymin());
    line1->SetY2(gPad->GetUymax());
    line1->SetLineColor(kBlue);
    line1->Draw();
}



























