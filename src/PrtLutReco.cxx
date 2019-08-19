//2018  ok// ----------w-------------------------------
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
#define prt__sim
#include "../../prttools/datainfo.C"
#include "../../prttools/prttools.C"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "Randomize.hh"
using std::cout;
using std::endl;

//TH2F*  hist_phs_xy = new TH2F("hist_phs_xy",";pos x [mm];pos y [mm]", 800,-30,30, 800,-30,30);
TH1F*  hist_phs_solution_number = new TH1F("hist_phs_solution_number",";number of solutions [#];entries [#]", 500,0,500);
TH2F*  hist_time_angle_all = new TH2F("hist_time_angle_all",";time [ns];solution angle [m rad]", 250,0,50,800,0,4 );
TH2F*  hist_time_angle_all_cut = new TH2F("hist_time_angle_all_cut",";time [ns];solution angle [m rad]", 250,0,50,800,0,4 );
TH2F*  hist_time_angle_all_tcut = new TH2F("hist_time_angle_all_tcut",";time [ns];solution angle [m rad]", 250,0,50,800,0,4 );
//TH2F*  hist_time_angle_all_correlation = new TH2F("hist_time_angle_correlation","; #Delta time [ns]; Angle [m rad]", 1000,-25,25,2000,0,4 );
//TH2F*  hist_time_angle_all_correlation_tcut = new TH2F("hist_time_angle_correlation_tcut","; #Delta time [ns]; Angle [m rad]", 250,-25,25,2000,0,4 );
//TH2F*  hist_time_angle_all_correlation_acut = new TH2F("hist_time_angle_correlation_acut","; #Delta time [ns]; Angle [m rad]", 250,-25,25,2000,0,4 );

TH1F*  fHist0 = new TH1F("timeDiff",";t_{phs}-t_{measured} [ns];entries [#]", 500,-50,50);
TH1F*  fHist0i = new TH1F("phsTime",";t_{phs} [ns];entries [#]", 250,0,50);
TH1F*  fHist0i_bg = new TH1F("hitTimed",";t_{measured} [ns];entries [#]", 250,0,50);
TH1F*  hit_b4_phs = new TH1F("hit_b4_phs",";t_{measured} [ns];entries [#]", 250,0,50);

//TH1F*  fHist1 = new TH1F("time1",";measured time [ns];entries [#]",   1000,-500,500);  // test
//TH1F*  fHist2 = new TH1F("time2",";calculated time [ns];entries [#]", 1000,-500,500); //test
TH1F*  hist_ambiguity = new TH1F("hist_ambiguity",";Pixel number [#]; Ambiguity [#]", 770,0,770); //test

TH1F*  fHist1 = new TH1F("time1",";measured time [ns];entries [#]",   500,0,50);  // test
TH1F*  fHist2 = new TH1F("time2",";calculated time [ns];entries [#]", 500,0,50); //test



TH2F*  lut_pix_pos_xy = new TH2F("lut_pix_pos_xy",";X pos []; Y pos []", 400,-100,300, 400,-200,200);


TH2F*  fHist3 = new TH2F("time3",";calculated time [ns];measured time [ns]", 500,0,80, 500,0,40);
TH2F*  fHist4 = new TH2F("time4",";#theta_{c}sin(#varphi_{c});#theta_{c}cos(#varphi_{c}", 100,-1,1, 100,-1,1);
TH2F*  fHist5 = new TH2F("time5",";#theta_{c}sin(#varphi_{c});#theta_{c}cos(#varphi_{c}", 100,-1,1, 100,-1,1);
TH1F*  falpha = new TH1F("alpha",";|#alpha|[degree];entries [#]",   1000,-360,360);
TH1F*  falphai = new TH1F("alphai",";|#alphaf| [degree];entries [#]",   1000,-360,360);
TH1F*  fHistPhotonEnergy = new TH1F("fHistPhotonEnergy",";|#alpha|[degree];entries [#]", 200, 0, 8);

TH1F*  fnHits = new TH1F("fnHits",";number of photons per track [#];entries [#]",   100,0,200);
TH1F*  fnHits_p = new TH1F("fnHits_p",";number of photons per proton track [#];entries [#]",   100,0,200);
TH1F*  fnHits_p_good = new TH1F("fnHits_p_good",";number of Good photons per proton track [#];entries [#]",   100,0,200);

TH1F*  nHits_dac = new TH1F("nHits_dac",";number of photons solutions per proton track [#];entries [#]",   100,0,200);
TH1F*  nHits_dac_syscut_p = new TH1F("nHits_dac_syscut_p",";number of photons per proton track after sys cut [#];entries [#]",   100,0,200);
TH1F*  fnHits_true_sim = new TH1F("fnHits_true_sim",";number of photons per track [#];entries [#]",   100,0,200);

//TH1F *hLnDiffP = new TH1F("hLnDiffP",  ";ln L(p) - ln L(#pi);entries [#]",200,-20000,20000);
//TH1F *hLnDiffPi = new TH1F("hLnDiffPi",";ln L(p) - ln L(#pi);entries [#]",200,-20000,20000);

// in case of standared method choose smaller range for the liklhood histo

//top3

//TH1F *hLnDiffP = new TH1F("hLnDiffP",  ";ln L(p) - ln L(#pi);entries [#]",200,-100, 100);
//TH1F *hLnDiffPi = new TH1F("hLnDiffPi",";ln L(p) - ln L(#pi);entries [#]",200,-100, 100);

//TH1F *hLnDiffP = new TH1F("hLnDiffP",  ";ln L(p) - ln L(#pi);entries [#]",200,-500, 500);
//TH1F *hLnDiffPi = new TH1F("hLnDiffPi",";ln L(p) - ln L(#pi);entries [#]",200,-500, 500);


//TH1F *hLnDiffP = new TH1F("hLnDiffP",  ";ln L(p) - ln L(#pi);entries [#]",200,-2000, 2000);
//TH1F *hLnDiffPi = new TH1F("hLnDiffPi",";ln L(p) - ln L(#pi);entries [#]",200,-2000, 2000);

TH1F *hLnDiffP = new TH1F("hLnDiffP",  ";ln L(p) - ln L(#pi);entries [#]",200,-5000, 5000);
TH1F *hLnDiffPi = new TH1F("hLnDiffPi",";ln L(p) - ln L(#pi);entries [#]",200,-5000, 5000);



//TH1F *hLnDiffP = new TH1F("hLnDiffP",  ";ln L(p) - ln L(#pi);entries [#]",1000,-1000,1000);
//TH1F *hLnDiffPi = new TH1F("hLnDiffPi",";ln L(p) - ln L(#pi);entries [#]",1000,-1000,1000);

TF1 *fcut = new TF1("fcut", "exp(-[0]/(x-[1])+[2])", 0.0, 50);
TF1 *fcut2 = new TF1("fcut2", "exp(-[0]/(x-[1])+[2])", 0.0, 50);




TH1F*  test_hist = new TH1F("test_hist",";t_{calc}-t_{measured} [ns];entries [#]", 500,-1000,1000);

TF1 *gF1 = new TF1("gaus0","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0.7,0.9);
TF1 *gF2= new TF1("gaus0","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2])",0.7,0.9);
Int_t gg_i(0), gg_ind(0);
TGraph gg_gr;
PrtLutNode *fLutNode[5000];
PrtLutNode *fLutNode_phs[5000];


TH1F*  fHistMcp[12]; // changed
TH1F*  fHistMcp_same_path[12]; // changed
TH1F*  fHistCh[960], *fHistCh_read_p[960], *fHistCh_read_pi[960], *hist_nph_wtc_p, *hist_nph_wtc_pi;
TH2F* hist_time_angle_ch[960], *hist_time_angle_ch2[960], *hist_time_angle_ch3[960];

TH2F*  fHistCh_2D[800], *fHistCh_2D_read_p[800], *fHistCh_2D_read_pi[800];



TGraph *fHistCh_graph_p[960], *fHistCh_graph_pi[960];

TH1F*  hit_time_ch[960],*  hit_time_phs_ch[960], * hit_time_diff_ch[960];
TH1F*  hit_time_ch_ccut[960],*  hit_time_phs_ch_ccut[960], * hit_time_diff_ch_ccut[960];

// tof
Double_t TOT_tof1, TOT_tof2, TOT_trg1, TOT_trg2, TOT_trgmzH, TOT_trgmzV, delta_tof2tof1;
//tof2tof1
Int_t bin_d_tof2tof1=100; //400;
Double_t x_d_tof2tof1=30.7; // 30;//25;//65
Double_t y_d_tof2tof1=34; // 35;//85;//75



TAxis *x_axis_data_p[960], *x_axis_data_pi[960];
Int_t x_bmin_data_p[960],  x_bmax_data_p[960], x_bmin_data_pi[960], x_bmax_data_pi[960];

TAxis *y_axis_data_p[960], *y_axis_data_pi[960];
Int_t y_bmin_data_p[960],  y_bmax_data_p[960], y_bmin_data_pi[960], y_bmax_data_pi[960];



Double_t integral_data[960], integral_data_pi[960], scale_p[960], scale_pi[960] ;


Double_t xmin_data = 0;
Double_t xmax_data = 4.0;

Double_t ymin_data = 0;
Double_t ymax_data = 50;

Double_t corrected_ch(0), non_corrected_ch(0);



/*
 // LE
 TH1F * htof1_le =       new TH1F("htof1_le","",400, -500, 500);
 TH1F * htof2_le =       new TH1F("htof2_le","",300, -275, -125); // 150
 TH1F * htrg1_le =       new TH1F("htrg1_le","",100, -160, -110); //50
 TH1F * htrg2_le =       new TH1F("htrg2_le","",400, -500, 500);
 TH1F * htrgmzH_le =      new TH1F("htrgmzH_le","",400, -500, 500);
 TH1F * htrgmzV_le =      new TH1F("htrgmzV_le","",400, -500, 500);
 // TOT histo
 TH1F * htof1_tot =      new TH1F("htof1_tot","",1000, 110, 150);
 TH1F * htof2_tot =      new TH1F("htof2_tot","",1000, 110, 150);
 TH1F * htrg1_tot =      new TH1F("htrg1_tot","",100, 110, 125);
 TH1F * htrg2_tot =      new TH1F("htrg2_tot","",1000, 110, 150);
 TH1F * htrgmzH_tot =     new TH1F("htrgmzH_tot","",1000, 110, 150);
 TH1F * htrgmzV_tot =     new TH1F("htrgmzV_tot","",1000, 110, 150);
 */

// delta
TH1F *hdelta_tof2tof1= new TH1F("hdelta_tof2tof1","",bin_d_tof2tof1, x_d_tof2tof1, y_d_tof2tof1);
TH1F *hdelta_tof2tof1_isproton= new TH1F("hdelta_tof2tof1_isproton","",bin_d_tof2tof1, x_d_tof2tof1, y_d_tof2tof1);
TH1F *histo_photon_ambiguity_wo= new TH1F("histo_photon_ambiguity_wo","",100, 0, 100);
TH1F *histo_photon_ambiguity_wt= new TH1F("histo_photon_ambiguity_wt","",100, 0, 100);
TH1F *histo_photon_ambiguity_wtc= new TH1F("histo_photon_ambiguity_wtc","",100, 0, 100);

/*
 // hodo
 const int nHodoPixel=16;
 TH2F * hHodo= new TH2F("hHodo","",16,0,16,16,0,16); // 9 sep 2017
 TH1F * hHodoPixelH[nHodoPixel],hHodoPixelV[nHodoPixel] ;
 */
TH1D * hodoV = new TH1D("1","vertical fibers",16,0,16);
TH1D * hodoH = new TH1D("2","horizontal fibers",16,0,16);
TH2D * hodoF = new TH2D("3","4",16,0,16,16,0,16);

TH1D * hodoV_select = new TH1D("1_select","vertical fibers selection",16,0,16);
TH1D * hodoH_select = new TH1D("2_select","horizontal fibers selection ",16,0,16);
TH2D *hodo_afterCut = new TH2D("5","6",16,0,16,16,0,16);
TH2D *hodo_multi_withmedHfiber = new TH2D("7","8",16,0,16,16,0,16);
TH2D *hodo_multi_withmedVfiber = new TH2D("9","10",16,0,16,16,0,16);

// 2D histo
// LE TOT histo
TH2F *htof1_le_tot = new TH2F("htof1_le_tot","",500, -300, -50, 100, 10, 60);
TH2F *htof2_le_tot = new TH2F("htof2_le_tot","",500, -300, -50, 100, 10, 60);
TH2F *htrg1_le_tot = new TH2F("htrg1_le_tot","",500, -300, -50, 200, 105, 126);
TH2F *htrg2_le_tot = new TH2F("htrg2_le_tot","",500, -300, -50, 200, 125, 140);
TH2F *htrgmzH_le_tot = new TH2F("htrgmzH_le_tot","",500,-300, -50, 200, 128, 145);
TH2F *htrgmzV_le_tot = new TH2F("htrgmzV_le_tot","",500,-300, -50, 200, 133, 150);

// 1D histo
TH1F * countmulti_tof1 =        new TH1F("countmulti_tof1","",7, 0, 7);
TH1F * countmulti_tof2 =        new TH1F("countmulti_tof2","",7, 0, 7);
TH1F * countmulti_trg1 =        new TH1F("countmulti_trg1","",7, 0, 7);
TH1F * countmulti_trg2 =        new TH1F("countmulti_trg2","",7, 0, 7);
TH1F * countmulti_trgmzH =       new TH1F("countmulti_trgmzH","",7, 0, 7);
TH1F * countmulti_trgmzV =       new TH1F("countmulti_trgmzV","",7, 0, 7);

TH1F * countmulti_hodoV =        new TH1F("countmulti_hodoV","",7, 0, 7);
TH1F * countmulti_hodoH =        new TH1F("countmulti_hodoH","",7, 0, 7);

// -----   Default constructor   -------------------------------------------
PrtLutReco::PrtLutReco(TString infile, TString lutfile, Int_t verbose) {


    //    TFile *ffile_error_polar_sim_p,
    //    ffile_error_tcut_sim_p = TFile::Open(calc_error_tcut_sim_p_path,"read");
    //phs_event =0 ;
    //mcpid_phs=0;
    //pixid_phs=0;
    //for ( auto m=0 ; m<12 ; ++m ){
    // for ( auto p=0 ; p<64 ; ++p ){
    //chain_phs[m][p] = new TChain("data");
    //TString phs_file_path = Form("/lustre/nyx/panda/aali/phs/build/phs/1m_phs_mcp_%d_pix_%d_all.root", m, p);
    //TString phs_file_path = Form("/Users/ahmed/phs/build/phs_space/phs_mcp_%d_pix_%d_all.root", m, p);
    //chain_phs[m][p]->Add(phs_file_path);
    //chain_phs[m][p]->SetBranchAddress("PrtEvent", &phs_event);
    //entries_phs[m][p] = chain_phs[m][p]->GetEntries();
    //std::cout<<"phs_file_path= "<<phs_file_path<<"Entries in chain:  "<<entries_phs[m][p] <<std::endl;
    /*
     fFile_phs[m][p] = new TFile(phs_file_path);
     fTree_phs[m][p]=(TTree *) fFile_phs[m][p]->Get("prtlut") ;
     fLut_phs[m][p] = new TClonesArray("PrtLutNode");
     fTree_phs[m][p]->SetBranchAddress("LUT",&fLut_phs[m][p]);
     fTree_phs[m][p]->GetEntry(0);
     fLutNode_phs[m][p] = (PrtLutNode*) fLut_phs[m][p]->At(0);
     */
    //  }
    //}

    //TString phs_file_path = "/lustre/nyx/panda/aali/phs_v2/prtdirc/build/1m_lut_all_20_447_l3.root";


    //top3
    //TString phs_file_path = "/lustre/nyx/panda/aali/phs_v2/prtdirc/build/final_1m_lut_all_20_447_l6_w11.root";
    //TString phs_file_path = "/lustre/nyx/panda/aali/phs_v2/prtdirc/build/pdf/final_1m_lut_all_20_447_l6_org.root";
    // TString phs_file_path = "/lustre/nyx/panda/aali/phs_v2/prtdirc/build/pdf/final_10m_lut_all_20_447_l6_sub_250_2000.root";

    TString phs_file_path = "/lustre/nyx/panda/aali/phs_v2/prtdirc/build/100k_lut_all_20_447_l6_randPix.root";//final_1m_lut_all_20_447_l6_randPix.root";
    //TString phs_file_path = "/lustre/nyx/panda/aali/phs_v2/prtdirc/build/100k_lut_all_20_447_l6_randPix.root";

    //TString phs_file_path = "/lustre/nyx/panda/aali/phs_v2/prtdirc/build/pdf/final_10m_lut_all_20_447_l6_sub_250_2000_f1mm.root";
    //TString phs_file_path = "/lustre/nyx/panda/aali/phs_v2/prtdirc/build/pdf/final_10m_lut_all_20_447_l6_f1mm.root";

    // TString phs_file_path = "/lustre/nyx/panda/aali/phs_v2/prtdirc/build/1m_lut_all_20_447.root";
    //TString phs_file_path = "/lustre/nyx/panda/aali/phs/build/1m_lut_all_20_447.root";
    //TString phs_file_path = "/lustre/nyx/panda/aali/sim/prtdirc/build/1m_lut_all_150.root";
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
    fHist = new TH1F("fHist",  "chrenkov angle;#theta_{C} [rad];entries [#]", 80,0.6,1); //150  //80     // warning
    fHist_copy = new TH1F("fHist_copy",  "chrenkov angle;#theta_{C} [rad];entries [#]", 80,0.6,1); //150  //80
    fHist_correction = new TH1F("RecoCherenkov",  "cherenkov angle;#theta_{C} [rad];entries [#]",500,0.6,1); //150  //80 ,300,0,5 );//
    //fHist_correction = new TH1F("fHist_correction",  "chrenkov angle;#theta_{C} [rad];entries [#]",300,0,5); //,80,0.6,1); //150  //80 ,300,0,5 );//


    fHist_same_path = new TH1F("fHist_same_path",  "chrenkov angle;#theta_{C} [rad];entries [#]", 80,0.6,1); //150  //80
    fHist_same_path_wotc = new TH1F("fHist_same_path_wotc",  "chrenkov angle;#theta_{C} [rad];entries [#]", 80,0.6,1); //150  //80


    fHist_bg = new TH1F("fHist_bg",  "chrenkov angle;#theta_{C} [rad];entries [#]", 80,0.6,1); //150  //80
    fHistPi = new TH1F("chrenkov_angle_hist_Pi",  "chrenkov angle pi;#theta_{C} [rad];entries [#]", 80,0.6,1); //150
    fHisti = new TH1F("chrenkov_angle_histi","chrenkov angle;#theta_{C} [rad];entries [#]", 80,0.6,1); //150


    fFit = new TF1("fgaus","[0]*exp(-0.5*((x-[1])/[2])*(x-[1])/[2]) +x*[3]+[4]",0.35,0.9);
    fSpect = new TSpectrum(10);
    if(infile.Contains("beam_")) {
        TString fileid(infile);
        fileid.Remove(0,fileid.Last('/')+1);
        fileid.Remove(fileid.Last('.')-1);
        prt_data_info = getDataInfo(fileid);
        TString opath(infile);
        opath.Remove(opath.Last('/'));
        if(infile.Contains("C.root")) {
            prt_savepath = opath+Form("/%dr/%d",prt_data_info.getStudyId(),prt_data_info.getFileId());
        } else {
            prt_savepath = opath+Form("/%ds/%d",prt_data_info.getStudyId(),prt_data_info.getFileId());
        }
    } else prt_savepath="data/sim";
    std::cout<<"fSavePath  "<< prt_savepath <<std::endl;
    for(Int_t i=0; i<5000; i++) {
        fLutNode[i] = (PrtLutNode*) fLut->At(i);
        fLutNode_phs[i] = (PrtLutNode*) fLut_phs->At(i);
    }


    for(Int_t i=0; i<770; i++) {
        Int_t direction_lut =fLutNode[i]->Entries();
        hist_ambiguity->Fill(i, direction_lut);
    }


    Double_t pos_x,pos_y, pos_z;
    for (Int_t mcpid_int=0; mcpid_int<12; mcpid_int++) {
        for (Int_t pixid_int=1; pixid_int<65; pixid_int++) {
            Int_t sensorId_int = 100*mcpid_int + pixid_int ;
            pos_x   = fLutNode[sensorId_int]->GetDigiPos().X();
            pos_y   = fLutNode[sensorId_int]->GetDigiPos().Y();
            pos_z   = fLutNode[sensorId_int]->GetDigiPos().Z();
            lut_pix_pos_xy->Fill(pos_x ,pos_y);
        }
    }


    cout << "-I- PrtLutReco: Intialization successfull" << endl;
    for(Int_t i=0; i<prt_nmcp; i++) {
        fHistMcp[i] = new TH1F(Form("fHistMcp_%d",i),Form("fHistMcp_%d;#theta_{C} [rad];entries [#]",i), 80,0.6,1); //150
    }

    for(Int_t i=0; i<prt_nmcp; i++) {
        fHistMcp_same_path[i] = new TH1F(Form("fHistMcp_same_path_%d",i),Form("fHistMcp_same_path_%d;#theta_{C} [rad];entries [#]",i), 80,0.6,1); //150
    }

    for(Int_t i=0; i<800; i++) {

        // fHistCh[i] = new TH1F(Form("fHistCh_%d",i),Form("fHistCh_%d;#theta_{C} [rad];entries [#]",i), 2000,0.6,1); //150
        fHistCh_2D[i] = new TH2F(Form("fHistCh_2D_%d",i),Form("fHistCh_2D_%d; time [ns]; solution angle [m rad]",i),1000,0,25,2000,0.6,1 );

        hit_time_ch[i] = new TH1F(Form("hit_time_ch_%d",i),Form("hit_time_ch_%d;t_{measured} [ns];entries [#]",i), 500,-50,50);
        hit_time_phs_ch[i] = new TH1F(Form("hit_time_phs_ch_%d",i),Form("hit_time_phs_ch_%d;t_{reco} [ns];entries [#]",i), 500,-50,50);
        hit_time_diff_ch[i] = new TH1F(Form("hit_time_diff_ch_%d",i),Form("hit_time_diff_ch_%d;|t_{reco}-t_{measured}| [ns];entries [#]",i), 500,-50,50);

        hit_time_ch_ccut[i] = new TH1F(Form("hit_time_ch_ccut_%d",i),Form("hit_time_ch_ccut_%d;t_{measured} [ns];entries [#]",i), 500,-50,50);
        hit_time_phs_ch_ccut[i] = new TH1F(Form("hit_time_phs_ch_ccut_%d",i),Form("hit_time_phs_ch_ccut_%d;t_{reco} [ns];entries [#]",i), 500,-50,50);
        hit_time_diff_ch_ccut[i] = new TH1F(Form("hit_time_diff_ch_ccut_%d",i),Form("hit_time_diff_ch_ccut_%d;|t_{reco}-t_{measured}| [ns];entries [#]",i), 500,-50,50);

        hist_time_angle_ch[i] = new TH2F(Form("hist_time_angle_ch_%d",i),Form("hist_time_angle_ch_%d; time [ns]; solution angle [m rad]",i), 250,0,50,40,0.6,1 );
        hist_time_angle_ch2[i] = new TH2F(Form("hist_time_angle_ch2_%d",i),Form("hist_time_angle_ch2_%d; time [ns]; solution angle [m rad]",i), 250,0,50,40,0.6,1 );
        hist_time_angle_ch3[i] = new TH2F(Form("hist_time_angle_ch3_%d",i),Form("hist_time_angle_ch3_%d; time [ns]; solution angle [m rad]",i), 250,0,50,40,0.6,1 );
    }


    /*
     // hodo
     for(Int_t p=0; p<nHodoPixel; p++){
     hHodoPixelH[p]  = new TH1F(Form("hHodoPixelH_%d",p),Form("pixelH %d", p),400,0,25);
     hHodoPixelV[p]  = new TH1F(Form("hHodoPixelV_%d",p),Form("pixelV %d", p),400,0,25);
     }
     */

}

// -----   Destructor   ----------------------------------------------------
PrtLutReco::~PrtLutReco() {
}

Int_t mcpdata[12][65]; // changed
Int_t cluster[12][65]; // changed
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
    cout << "No Problem 2" << endl;

    TVector3 dird, dir, momInBar(0,0,1),posInBar,cz, momInBar_phs(0,0,0);
    Double_t mom, cangle,spr,tangle,likelihood(0),xlikelihood(0),ylikelihood(0),boxPhi,weight,evtime,bartime, lenz,dirz,luttheta, barHitTime, hitTime, photonEnergy, alphAngle;
    Int_t  tofPid(0),distPid(0),likePid(0),pdgcode, evpointcount=0;
    Bool_t reflected = kFALSE;
    gStyle->SetOptFit(111);
    TVector3 fnX1 = TVector3 (1,0,0);
    TVector3 fnY1 = TVector3( 0,1,0);
    bool testTrRes = false;
    Double_t angdiv,dtheta,dtphi,prtangle, counter(0);

    /////////////////////
    /////////////////////
    /////////////////////

    Double_t theta(0),phi(0), trr(0),  nph(0),par1(0), par2(0), par3(0), par4(0), par5(0), par6(0), test1(0), test2(0), test3(0),separation(0),recoAngle(0), chAngleCut(0), timeRes(0);
    Double_t recoP(0), recoPi(0), gPDF(0), method_type(-1),xseparation(0),yseparation(0), all_separation(0);
    Int_t solution_number_approach_selection(0),solution_number(0);
    Int_t openChCorr(0);

    chAngleCut = 0.04; //PrtManager::Instance()->GetchAngleCut();
    //recoAngle = PrtManager::Instance()->GetrecoAngle();
    timeRes = PrtManager::Instance()->GetTimeRes();

    //top3

    // 1 1 0 spr sim , 1 1 1  pdf creation, 0 0 2 separation using pdf creation , 003 separation std
    recoP = 1;//PrtManager::Instance()->GetrecoP();
    recoPi= 0; // PrtManager::Instance()->GetrecoPi();
    gPDF = 0; // PrtManager::Instance()->GetgPDF();
    openChCorr =0 ;// PrtManager::Instance()->GetopenChCorr();
    method_type = PrtManager::Instance()->GetRunType();

    cout<<"@@@@@@@@@@@@ gPDF="<< gPDF << endl;
    cout<<"@@@@@@@@@@@@ recoP="<< recoP << endl;
    cout<<"@@@@@@@@@@@@ recoPi="<< recoPi << endl;
    cout<<"@@@@@@@@@@@@ method_type="<< method_type << endl;


    TFile *ffile_data_p, *ffile_data_pi;
    Double_t prtangle_pdf, pdf_nph_p, pdf_nph_pi;
    TString cherenkov_data_p_path, cherenkov_data_pi_path;



    TVector3 direction, direction2  ;
    Double_t time_phs, pos_x, pos_y, pos_z ;







    if (method_type == 3) {
        prtangle_pdf = 20;

        cout<<"@@@@@@@@@@@@ prtangle_pdf="<< prtangle_pdf << endl;

        //top3

        cherenkov_data_p_path = "/lustre/nyx/panda/aali/phs_v2/prtdirc/build/pdf/histo_p_final_1m_lut_all_20_447_l6_randPix_w10_cherenkovPDF.root";
        cherenkov_data_pi_path ="/lustre/nyx/panda/aali/phs_v2/prtdirc/build/pdf/histo_pi_final_1m_lut_all_20_447_l6_randPix_w10_cherenkovPDF.root";

        //cherenkov_data_p_path = "/lustre/nyx/panda/aali/phs_v2/prtdirc/build/pdf/histo_p_final_1m_lut_all_20_447_l6_org_cherenkovPDF.root";
        //cherenkov_data_pi_path ="/lustre/nyx/panda/aali/phs_v2/prtdirc/build/pdf/histo_pi_final_1m_lut_all_20_447_l6_org_cherenkovPDF.root";

        //cherenkov_data_p_path = "/lustre/nyx/panda/aali/phs_v2/prtdirc/build/pdf/histo_p_final_10m_lut_all_20_447_l6_sub_250_2000_cherenkovPDF.root";
        //cherenkov_data_pi_path ="/lustre/nyx/panda/aali/phs_v2/prtdirc/build/pdf/histo_pi_final_10m_lut_all_20_447_l6_sub_250_2000_cherenkovPDF.root";

        //cherenkov_data_p_path = "/lustre/nyx/panda/aali/phs_v2/prtdirc/build/pdf/histo_p_final_10m_lut_all_20_447_l6_sub_250_2000_f1mm_cherenkovPDF.root";
        //cherenkov_data_pi_path ="/lustre/nyx/panda/aali/phs_v2/prtdirc/build/pdf/histo_pi_final_10m_lut_all_20_447_l6_sub_250_2000_f1mm_cherenkovPDF.root";

        //cherenkov_data_p_path = "/lustre/nyx/panda/aali/phs_v2/prtdirc/build/pdf/histo_p_final_10m_lut_all_20_447_l6_f1mm_cherenkovPDF.root";
        //cherenkov_data_pi_path ="/lustre/nyx/panda/aali/phs_v2/prtdirc/build/pdf/histo_pi_final_10m_lut_all_20_447_l6_f1mm_cherenkovPDF.root";



        //cherenkov_data_p_path = Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/data/332/ambiguit_pdf/histo_%g_sph_p_data_cherenkovPDF.root", prtangle_pdf);
        //cherenkov_data_pi_path = Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/data/332/ambiguit_pdf/histo_%g_sph_pi_data_cherenkovPDF.root", prtangle_pdf);
        //cherenkov_data_p_path = Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/data/332/ambiguit_pdf/histo_2BarRefl_%g_sph_p_data_cherenkovPDF.root", prtangle_pdf);
        //cherenkov_data_pi_path = Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/data/332/ambiguit_pdf/histo_2BarRefl_%g_sph_pi_data_cherenkovPDF.root", prtangle_pdf);
        //data
        //cherenkov_data_p_path = Form("/lustre/nyx/panda/aali/phs_v2/prtdirc/build/pdf/histo_%g_cyl_p_sim_hd_cherenkovPDF.root", prtangle_pdf);
        //cherenkov_data_pi_path = Form("/lustre/nyx/panda/aali/phs_v2/prtdirc/build/pdf/histo_%g_cyl_pi_sim_hd_cherenkovPDF.root", prtangle_pdf);




        //sim
        //cherenkov_data_p_path = Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/sim/332/pdf/histo_%g_sph_p_sim_cherenkovPDF.root", prtangle_pdf);
        //cherenkov_data_pi_path = Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/sim/332/pdf/histo_%g_sph_pi_sim_cherenkovPDF.root", prtangle_pdf);
        //test  number of event sgenerate PDF
        //cherenkov_data_p_path =  Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/sim/332/histo_sim_%d_%g_sph_p_data_cherenkovPDF.root", openChCorr,prtangle_pdf);
        //cherenkov_data_pi_path = Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/sim/332/histo_sim_%d_%g_sph_pi_data_cherenkovPDF.root",openChCorr,prtangle_pdf);

        //cherenkov_data_p_path = Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/data/332/histo_%d_sph_p_data_cherenkovPDF.root", 40);
        //cherenkov_data_pi_path = Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/workspace/testbeam/recon/data/332/histo_%d_sph_pi_data_cherenkovPDF.root", 40);
        cout<<"cherenkov_data_p_path= " <<cherenkov_data_p_path<<endl;
        cout<<"cherenkov_data_pi_path= " <<cherenkov_data_pi_path<<endl;
        ffile_data_p  = new TFile(cherenkov_data_p_path, "READ");
        ffile_data_pi  = new TFile(cherenkov_data_pi_path, "READ");
        for(Int_t pix=0; pix<800; pix++) {
            //fHistCh_graph_p[pix] =new TGraph(   (TH1F*)ffile_data_p->Get(Form("fHistCh_%d",pix))   );
            //i[pix] =new TGraph(   (TH1F*)ffile_data_pi->Get(Form("fHistCh_%d",pix))   );
            fHistCh_2D_read_p[pix] = (TH2F*)ffile_data_p->Get(Form("fHistCh_2D_%d",pix));
            fHistCh_2D_read_pi[pix] = (TH2F*)ffile_data_pi->Get(Form("fHistCh_2D_%d",pix));
        }
        cout<<"@@@@@@@@@@@@ no problem after read histo"<< endl;
        //hist_nph_wtc_p=(TH1F*)ffile_data_p->Get("fnHits_p_good");
        //hist_nph_wtc_pi=(TH1F*)ffile_data_pi->Get("fnHits_p_good");
        //Double_t pdf_nph_p=42556;
        //Double_t pdf_nph_pi=45334;
        //pdf_nph_p=hist_nph_wtc_p->GetEntries();
        //pdf_nph_pi=hist_nph_wtc_pi->GetEntries();
        //cout<<"@@@@@@@@@@@@ pdf_nph_p="<< pdf_nph_p << endl;
        //cout<<"@@@@@@@@@@@@ pdf_nph_pi="<< pdf_nph_pi << endl;
        /*
        for(Int_t pix=0; pix<8000; pix++) {
            x_axis_data_p[pix] = fHistCh_2D_read_p[pix]->GetXaxis();
            x_bmin_data_p[pix] = x_axis_data_p[pix]->FindBin(xmin_data);
            x_bmax_data_p[pix] = x_axis_data_p[pix]->FindBin(xmax_data);

            y_axis_data_p[pix] = fHistCh_2D_read_p[pix]->GetXaxis();
            y_bmin_data_p[pix] = y_axis_data_p[pix]->FindBin(xmin_data);
            y_bmax_data_p[pix] = y_axis_data_p[pix]->FindBin(xmax_data);



            //integral_data[pix] = fHistCh_read_p[pix]->Integral(x_bmin_data_p[pix],x_bmax_data_p[pix]);
            //fHistCh_read_p[pix]->Scale(1/integral_data[pix]);
            //fHistCh_read_p[pix]->Scale(1/pdf_nph_p);

            x_axis_data_pi[pix] = fHistCh_read_pi[pix]->GetXaxis();
            x_bmin_data_pi[pix] = x_axis_data_pi[pix]->FindBin(xmin_data);
            x_bmax_data_pi[pix] = x_axis_data_pi[pix]->FindBin(xmax_data);

            y_axis_data_pi[pix] = fHistCh_read_pi[pix]->GetXaxis();
            y_bmin_data_pi[pix] = y_axis_data_pi[pix]->FindBin(xmin_data);
            y_bmax_data_pi[pix] = y_axis_data_pi[pix]->FindBin(xmax_data);


            //integral_data_pi[pix] = fHistCh_read_pi[pix]->Integral(x_bmin_data_pi[pix],x_bmax_data_pi[pix]);
            //fHistCh_read_pi[pix]->Scale(1/integral_data_pi[pix]);
            //fHistCh_read_pi[pix]->Scale(1/pdf_nph_pi);

            //normalize the histogram to 1
            //scale_p[pix] = 1/fHistCh_read_p[pix]->Integral();
            //scale_pi[pix] = 1/fHistCh_read_pi[pix]->Integral();
            //or
            //normalize histogram per entry and per unit of X axis
            //scale_p[pix] = x_axis_data_p[pix]->GetBinWidth(1)/fHistCh_read_p[pix]->GetIntegral(); // problem
            //scale_pi[pix]= x_axis_data_pi[pix]->GetBinWidth(1)/fHistCh_read_pi[pix]->GetIntegral(); // problem

            //fHistCh_read_p[pix]->Scale(scale_p[pix]);
            //fHistCh_read_pi[pix]->Scale(scale_pi[pix]);

            fHistCh_read_pi[pix]->SetLineColor(kRed);
            // pdf graphs not used
            fHistCh_graph_p[pix] =new TGraph(fHistCh_read_p[pix]);
            fHistCh_graph_pi[pix] =new TGraph(fHistCh_read_pi[pix]);
        }
        */
    }
    TString outFile;
    if(method_type==3) {
        outFile = PrtManager::Instance()->GetOutName()+"_separation.root";
    } else if(gPDF==1) {
        outFile = PrtManager::Instance()->GetOutName()+"_cherenkovPDF.root";
    } else {
        outFile = PrtManager::Instance()->GetOutName()+"_spr.root";
    }

    //TString outFile =PrtManager::Instance()->GetOutName()+"_spr.root" ;
    Double_t minChangle(0);
    Double_t maxChangle(1);
    Double_t rad = TMath::Pi()/180.;
    Double_t criticalAngle = asin(1.00028/1.47125); // n_quarzt = 1.47125; //(1.47125 <==> 390nm)
    Int_t nsEvents(0);

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
    tree.Branch("cangle",&cangle,"cangle/D");
    tree.Branch("likelihood",&likelihood,"par3/D");
    tree.Branch("separation",&separation,"separation/D");
    tree.Branch("par5",&par5,"par5/D");
    tree.Branch("par6",&par6,"par6/D");
    tree.Branch("test1",&test1,"test1/D");
    tree.Branch("test2",&test2,"test2/D");
    tree.Branch("test3",&test3,"test3/D");
    tree.Branch("theta",&theta,"theta/D");
    tree.Branch("phi",&phi,"phi/D");
    tree.Branch("chAngleCut",&chAngleCut,"chAngleCut/D");
    tree.Branch("recoAngle",&recoAngle,"recoAngle/D");
    tree.Branch("timeRes",&timeRes,"timeRes/D");
    tree.Branch("openChCorr",&openChCorr,"openChCorr/I");
    tree.Branch("end",&end,"end/I");
    tree.Branch("solution_number_approach_selection",&end,"solution_number_approach_selection/I");
    tree.Branch("solution_number",&end,"solution_number/I");
    tree.Branch("nsEvents",&nsEvents,"nsEvents/I");
    tree.Branch("corrected_ch",&corrected_ch,"corrected_ch/D");
    tree.Branch("non_corrected_ch",&non_corrected_ch,"non_corrected_ch/D");

    test1 = PrtManager::Instance()->GetTest1();
    test2 = PrtManager::Instance()->GetTest2();
    test3 = PrtManager::Instance()->GetTest3();
    Int_t radiator = PrtManager::Instance()->GetRadiator();
    fMethod = PrtManager::Instance()->GetRunType();
    Int_t nEvents = fChain->GetEntries();
    if(end==0) end = nEvents;
    //if (gPDF==1 )start = 500001; // for data
    if (gPDF==1 )start = 0;   // for simulation
    cout<<"@@@@@@@@@@@@ test1="<< test1 << endl;
    cout<<"@@@@@@@@@@@@ test2="<< test2 << endl;
    cout<<"@@@@@@@@@@@@  sim ?="<< fEvent->GetType()<< endl;
    cout<<"@@@@@@@@@@@@ gPDF="<< gPDF << "    start= "<< start << endl;
    std::cout<<"Run started for ["<<start<<","<<end <<"]"<<std::endl;
    Int_t nsHits(0),studyId(0), nHits(0), ninfit(1);
    if(start<0) {
        ninfit=abs(start);
        start=0;
    }
    Int_t counter_event_loop =0;
    for (Int_t ievent=start; ievent<end; ievent++) { //&& ievent<end
        std::cout<< "##########################  event = "<< counter_event_loop << std::endl;
        fChain->GetEntry(ievent);
        nHits = fEvent->GetHitSize();
        nHits_dac->Fill(nHits);
        if(ievent%1==0) std::cout<<"Event # "<< ievent << " has "<< nHits <<" hits"<<std::endl;
        //if (ievent > 1) if (nHits < 40 || nHits > 50) continue;
        if (ievent >1 && nHits < 17) continue;
        if(ievent-start==0) {
            tree.SetTitle(fEvent->PrintInfo());
            prtangle = fEvent->GetAngle(); //changed 2017
            // here
            std::cout<<"@@@@@@@@@@@@@@@@@@@"<<" prtangle= "<<prtangle<<std::endl;
            // 5 sigma of true path Cherenkov and time diff. devided by 5
            if( prtangle==20) {
                chAngleCut=0.038775/5.0*chAngleCut;
                timeRes=1.862127/5.0*timeRes;
            }
            if( prtangle==30) {
                chAngleCut=0.041868/5.0*chAngleCut;
                timeRes=1.910943/5.0*timeRes;
            }
            if( prtangle==40) {
                chAngleCut=0.042588/5.0*chAngleCut;
                timeRes=1.878375/5.0*timeRes;
            }
            if( prtangle==50) {
                chAngleCut=0.041037/5.0*chAngleCut;
                timeRes=1.895868/5.0*timeRes;
            }
            if( prtangle==60) {
                chAngleCut=0.045299/5.0*chAngleCut;
                timeRes=1.894080/5.0*timeRes;
            }
            if( prtangle==70) {
                chAngleCut=0.046581/5.0*chAngleCut;
                timeRes=1.966618/5.0*timeRes;
            }
            if( prtangle==80) {
                chAngleCut=0.045421/5.0*chAngleCut;
                timeRes=2.132993/5.0*timeRes;
            }
            if( prtangle==90) {
                chAngleCut=0.054653/5.0*chAngleCut;
                timeRes=1.823477/5.0*timeRes;
            }
            if( prtangle==100) {
                chAngleCut=0.047778/5.0*chAngleCut;
                timeRes=1.402801/5.0*timeRes;
            }
            if( prtangle==110) {
                chAngleCut=0.046473/5.0*chAngleCut;
                timeRes=1.383924/5.0*timeRes;
            }
            if( prtangle==120) {
                chAngleCut=0.044891/5.0*chAngleCut;
                timeRes=1.385515/5.0*timeRes;
            }
            if( prtangle==130) {
                chAngleCut=0.043477/5.0*chAngleCut;
                timeRes=1.412302/5.0*timeRes;
            }
            if( prtangle==140) {
                chAngleCut=0.041921/5.0*chAngleCut;
                timeRes=1.437529/5.0*timeRes;
            }
            if( prtangle==150) {
                chAngleCut=0.042944/5.0*chAngleCut;
                timeRes=1.529106/5.0*timeRes;
            }
            if( prtangle==160) {
                chAngleCut=0.042944/5.0*chAngleCut;
                timeRes=1.529106/5.0*timeRes;
            }

            studyId = fEvent->GetGeometry();
            if(studyId==152 || studyId==153 || studyId==161 || studyId==162 || studyId==171 || studyId==172 || studyId==173 || studyId==175 || studyId==176 || studyId==177 || studyId==178) {
                radiator=2;
            }
            mom=fEvent->GetMomentum().Mag(); // changed 2017
            //mom=7.0; // here
            std::cout<<"No Problem1  "<< " mom ="<< mom <<std::endl; // changed
            Double_t beam_corr(0);
            //if(studyId==151) beam_corr = 0.0045; //125 deg //!
            // if(studyId==151) beam_corr = 0.001; // 20 deg
            // if(studyId==151) beam_corr = -0.003; // 25 deg
            //beam_corr = 0.002; // 125 deg s160
            if(fEvent->GetType()==0 && test1==99) {

                // 332 study ID beam correction
                if (prtangle==20) {
                    test1=0.0/1000  ;
                    test2= -6.5/1000 ;
                }
                if (prtangle==30) {
                    test1=1.0/1000 ;
                    test2= -6.0/1000 ;
                }
                if (prtangle==40) {
                    test1=-1.0/1000  ;
                    test2= 5.5/1000 ;
                }
                if (prtangle==50) {
                    test1=-2.0/1000  ;
                    test2= 5.5/1000 ;
                }
                if (prtangle==60) {
                    test1=-1.0/1000  ;
                    test2= 0.0/1000 ;
                }
                if (prtangle==70) {
                    test1=0/1000    ;
                    test2= 8.0/1000 ;
                }
                if (prtangle==80) {
                    test1=-2.0/1000   ;
                    test2= 5/1000 ;
                }
                if (prtangle==90) {
                    test1=1/1000 ;
                    test2= 0.0/1000 ;
                }
                if (prtangle==100) {
                    test1=0.0/1000 ;
                    test2= 7.0/1000 ;
                }
                if (prtangle==110) {
                    test1=1.0/1000 ;
                    test2= 0.0/1000 ;
                }
                if (prtangle==120) {
                    test1=1/1000 ;
                    test2=4.0/1000 ;
                }
                if (prtangle==130) {
                    test1=0.0/1000  ;
                    test2=1.0/1000 ;
                }
                if (prtangle==140) {
                    test1=0.0/1000 ;
                    test2=-1.0/1000 ;
                }
                if (prtangle==150) {
                    test1=0.0/1000 ;
                    test2=7.0/1000 ;
                }
                if (prtangle==160) {
                    test1=0.0/1000 ;
                    test2=7.0/1000 ;
                }



                momInBar.RotateY(TMath::Pi()-prtangle*rad-0);// 0  test1


                momInBar.RotateX(0);//0 test2

                //momInBar.RotateY(TMath::Pi()-prtangle*rad+0.0);
                //momInBar.RotateX(0.0);
                // -+, --, +-, ++

            } else {
                momInBar.RotateY(TMath::Pi()-prtangle*rad);
            }
            // momInBar = fEvent->GetMomentum().Unit();
            if(fVerbose==3) {
                cz = momInBar.Unit();
                cz = TVector3(-cz.X(),cz.Y(),cz.Z());
            }
        }
        //if(nHits<5) continue;  // changed
        //std::cout<<"@@@@@@@@"<<" prtangle= "<<prtangle<<std::endl; // changed
        Double_t momentum=fEvent->GetMomentum().Mag(); // changed 2017
        //Double_t momentum=7.0; // here
        if( fEvent->GetType()==1) momentum /= 1000;
        tofPid=fEvent->GetParticle();
        if(tofPid==212) tofPid=211;
        //std::cout<<"$$$$$$$$$$$$$$$$$$$ up $$$$$$ tofPid "<<tofPid <<std::endl;
        //std::cout<<"$$$$$$$$$$$$$$$$$$$$$$$$$  "<< fEvent->GetType() <<std::endl;


        Int_t pdg[]= {11,13,211,321,2212};
        Double_t mass[] = {0.000511,0.1056584,0.139570,0.49368,0.9382723};
        Double_t angle1(0), angle2(0),sum1(0),sum2(0), sigma(0.009),range(5*sigma),noise(0.3);

        fAngleP = acos(sqrt(momentum*momentum+ mass[4]*mass[4])/momentum/1.4738)-0.00; //1.4738 = 370 = 3.35
        fAnglePi= acos(sqrt(momentum*momentum + mass[2]*mass[2])/momentum/1.4738)-0.00; //-0.0014 for 160 25deg
        //std::cout<<" ####### ###### chrencove Angle  "<< " proton ="<< fAngleP <<std::endl; // changed
        //std::cout<<" ####### ###### chrencove Angle  "<< " pi  ="<< fAnglePi <<std::endl;   // changed
        gF1->SetParameter(0,1);
        gF2->SetParameter(0,1);

        // move the model

        gF1->SetParameter(1,fAngleP);
        gF2->SetParameter(1,fAnglePi);
        //std::cout<<"No Problem2  " <<std::endl;
        /*
         if(prtangle == 20)gF1->SetParameter(1,0.818334);
         if(prtangle == 30)gF1->SetParameter(1,0.816051);
         if(prtangle == 40)gF1->SetParameter(1,0.816986);
         if(prtangle == 50)gF1->SetParameter(1,0.817087);
         if(prtangle == 60)gF1->SetParameter(1,0.815137);
         if(prtangle == 70)gF1->SetParameter(1,0.813078);
         if(prtangle == 80)gF1->SetParameter(1,0.812915);
         if(prtangle == 90)gF1->SetParameter(1,0.810881);
         if(prtangle == 100)gF1->SetParameter(1,0.811012);
         if(prtangle == 110)gF1->SetParameter(1,0.813848);
         if(prtangle == 120)gF1->SetParameter(1,0.814565);
         if(prtangle == 130)gF1->SetParameter(1,0.817052);
         if(prtangle == 140)gF1->SetParameter(1,0.816796);
         if(prtangle == 150)gF1->SetParameter(1,0.815613);


         if(prtangle == 20)gF2->SetParameter(1,0.826920);
         if(prtangle == 30)gF2->SetParameter(1,0.824533);
         if(prtangle == 40)gF2->SetParameter(1,0.825388);
         if(prtangle == 50)gF2->SetParameter(1,0.824601);
         if(prtangle == 60)gF2->SetParameter(1,0.823733);
         if(prtangle == 70)gF2->SetParameter(1,0.820671);
         if(prtangle == 80)gF2->SetParameter(1,0.820424);
         if(prtangle == 90)gF2->SetParameter(1,0.819132);
         if(prtangle == 100)gF2->SetParameter(1,0.820097);
         if(prtangle == 110)gF2->SetParameter(1,0.821585);
         if(prtangle == 120)gF2->SetParameter(1,0.822756);
         if(prtangle == 130)gF2->SetParameter(1,0.824813);
         if(prtangle == 140)gF2->SetParameter(1,0.825208);
         if(prtangle == 150)gF2->SetParameter(1,0.823934);
         */

        gF1->SetParameter(2,sigma);
        gF2->SetParameter(2,sigma);
        //////////////////////////////
        // select proton sim
        if(fMethod==2 && tofPid!=2212 && fEvent->GetType()==1&& recoP==1 && recoPi==0) continue; //2212, 211 will remove pion , proton respectively in simulation
        if(fMethod==2 && tofPid!=211 && fEvent->GetType()==1&& recoP==1 && recoPi==0) goto goon;

        //select pi sim
        if(fMethod==2 && tofPid!=211 && fEvent->GetType()==1&& recoP==0 && recoPi==1) continue; //2212, 211 will remove pion , proton respectively in simulation
        if(fMethod==2 && tofPid!=2212 && fEvent->GetType()==1&& recoP==0 && recoPi==1) goto goon;

        if(fMethod==3 && fEvent->GetType()==1) goto goon;

        if(fEvent->GetType()==0) { // Pions and protons will enter the loop the filltering will be applyied from the tof system
            Bool_t tof1_bool (false);
            Bool_t tof2_bool (false);
            Bool_t trg1_bool (false);
            Bool_t trg2_bool = false;
            Bool_t trgmzH_bool = false;
            Bool_t trgmzV_bool = false;
            Bool_t multi_bool = false;
            Bool_t multi_bool_tof1 = false;
            Bool_t multi_bool_tof2 = false;
            Bool_t multi_bool_trg1 = false;
            Bool_t multi_bool_trg2 = false;
            Bool_t multi_bool_trgmzH = false;
            Bool_t multi_bool_trgmzV = false;
            Bool_t bad_bool_hodoH = false;
            Bool_t bad_bool_hodoV = false;


            //Vertical
            Bool_t hodo_bool_1349= false;
            Bool_t hodo_bool_1350= false;
            Bool_t hodo_bool_1351= false;
            Bool_t hodo_bool_1352= false;

            //Horizental
            Bool_t hodo_bool_1367= false;
            Bool_t hodo_bool_1368= false;
            Bool_t hodo_bool_1369= false;
            Bool_t hodo_bool_1370= false;
            Bool_t hodo_bool_1371= false;

            Double_t LE_tof1=0;
            Double_t LE_tof2=0;
            Double_t LE_trg1=0;
            Double_t LE_trg2=0;
            Double_t LE_trgmzH=0;
            Double_t LE_trgmzV=0;

            TOT_tof1=0;
            TOT_tof2=0;
            TOT_trg1=0;
            TOT_trg2=0;
            TOT_trgmzH=0;
            TOT_trgmzV=0;
            delta_tof2tof1=0;

            Int_t count_tof1 =0;
            Int_t count_tof2 =0;
            Int_t count_trg1 =0;
            Int_t count_trg2 =0;
            Int_t count_trgmzH =0;
            Int_t count_trgmzV =0;

            Int_t count_hodoH_1368 =0;
            Int_t count_hodoH_1369 =0;
            Int_t count_hodoH_1370 =0;

            Int_t count_hodoV_1350 =0;
            Int_t count_hodoV_1351 =0;
            Int_t count_hodoV_1352 =0;

            for(Int_t h=0; h<fEvent->GetHitSize(); h++) {
                fHit = fEvent->GetHit(h);
                Int_t gch=fHit.GetChannel();
                //std::cout<<"starting hits loop "<<std::endl;
                //std::cout<<"ievent"<<ievent<< " 	ihits= "<< h <<std::endl;
                // TOF 1
                if (gch == 1392) { //1392
                    LE_tof1 = fHit.GetLeadTime();
                    TOT_tof1= fHit.GetTotTime();
                    tof1_bool = true;
                    ++count_tof1;
                }
                if (count_tof1 > 1) {
                    multi_bool = true;
                    multi_bool_tof1= true;
                }
                // TOF 2
                if (gch == 1398) { //1398
                    LE_tof2 = fHit.GetLeadTime();
                    TOT_tof2= fHit.GetTotTime();
                    tof2_bool = true;
                    ++count_tof2;
                }
                if (count_tof2 > 1) {
                    multi_bool = true;
                    multi_bool_tof2= true;
                }
                // Trg1
                if (gch == 816) {
                    LE_trg1 = fHit.GetLeadTime();
                    TOT_trg1= fHit.GetTotTime();
                    trg1_bool = true;
                    ++count_trg1;
                }
                if (count_trg1 > 1) {
                    multi_bool = true;
                    multi_bool_trg1= true;
                }
                // trg2
                if (gch == 817) {
                    LE_trg2 = fHit.GetLeadTime();
                    TOT_trg2= fHit.GetTotTime();
                    trg2_bool = true;
                    ++count_trg2;
                }
                if (count_trg2 > 1) {
                    // multi_bool = true;
                    multi_bool_trg2= true;
                }
                // TrgmzH
                if (gch == 818) {
                    LE_trgmzH = fHit.GetLeadTime();
                    TOT_trgmzH= fHit.GetTotTime();
                    trgmzH_bool = true;
                    ++count_trgmzH;
                }
                if (count_trgmzH > 1) {
                    //  multi_bool = true;
                    multi_bool_trgmzH= true;
                }
                // TrgmzV
                if (gch == 819) {
                    LE_trgmzV = fHit.GetLeadTime();
                    TOT_trgmzV= fHit.GetTotTime();
                    trgmzV_bool = true;
                    ++count_trgmzV;
                }
                if (count_trgmzV > 1) {
                    //  multi_bool = true;
                    multi_bool_trgmzV= true;
                }

                if(0<=(gch-1344) && (gch-1344)<16) { // hodoscope V all
                    hodoV->Fill(gch-1344);
                }

                if(16<=(gch-1344) && (gch-1344)<32) { // hodoscope H all
                    hodoH->Fill(gch-1344-16);
                }


                // select bad fibers
                if (gch==1344||gch==1345||gch==1346||gch==1347||gch==1348||gch==1349||gch==1353||gch==1354||gch==1355||gch==1356||gch==1357||gch==1358||gch==1359) bad_bool_hodoV=true;
                if (gch==1360||gch==1361||gch==1362||gch==1363||gch==1364||gch==1365||gch==1366||gch==1367||gch==1371||gch==1372||gch==1373||gch==1374||gch==1375) bad_bool_hodoH=true;
                // fill good fibers

                // fill good fibers
                if (gch==1350 && bad_bool_hodoV==false && bad_bool_hodoH== false) { // Vertical
                    hodo_bool_1350=true;
                    for(Int_t h4=0; h4<fEvent->GetHitSize(); h4++) {
                        fHit4 = fEvent->GetHit(h4);
                        Int_t gch4=fHit4.GetChannel();
                        if (gch4==1344 ) ++count_hodoV_1350;
                        if (gch4==1345 ) ++count_hodoV_1350;
                        if (gch4==1346 ) ++count_hodoV_1350;
                        if (gch4==1347 ) ++count_hodoV_1350;
                        if (gch4==1348 ) ++count_hodoV_1350;
                        if (gch4==1349 ) ++count_hodoV_1350;
                        if (gch4==1351 ) ++count_hodoV_1350;
                        if (gch4==1352 ) ++count_hodoV_1350;
                        if (gch4==1353 ) ++count_hodoV_1350;
                        if (gch4==1354 ) ++count_hodoV_1350;
                        if (gch4==1355 ) ++count_hodoV_1350;
                        if (gch4==1356 ) ++count_hodoV_1350;
                        if (gch4==1357 ) ++count_hodoV_1350;
                        if (gch4==1358 ) ++count_hodoV_1350;
                        if (gch4==1359 ) ++count_hodoV_1350;
                        if (count_hodoV_1350==0) {
                            if (0<=(gch4-1344) && (gch4-1344)<16)  hodoV_select->Fill(gch4-1344);
                        }
                    }
                }
                if (gch==1351 && bad_bool_hodoV==false && bad_bool_hodoH== false) { // Vertical
                    hodo_bool_1351=true;
                    for(Int_t h2=0; h2<fEvent->GetHitSize(); h2++) {
                        fHit2 = fEvent->GetHit(h2);
                        Int_t gch2=fHit2.GetChannel();
                        if (gch2==1344 ) ++count_hodoV_1351;
                        if (gch2==1345 ) ++count_hodoV_1351;
                        if (gch2==1346 ) ++count_hodoV_1351;
                        if (gch2==1347 ) ++count_hodoV_1351;
                        if (gch2==1348 ) ++count_hodoV_1351;
                        if (gch2==1349 ) ++count_hodoV_1351;
                        if (gch2==1350 ) ++count_hodoV_1351;
                        if (gch2==1352 ) ++count_hodoV_1351;
                        if (gch2==1353 ) ++count_hodoV_1351;
                        if (gch2==1354 ) ++count_hodoV_1351;
                        if (gch2==1355 ) ++count_hodoV_1351;
                        if (gch2==1356 ) ++count_hodoV_1351;
                        if (gch2==1357 ) ++count_hodoV_1351;
                        if (gch2==1358 ) ++count_hodoV_1351;
                        if (gch2==1359 ) ++count_hodoV_1351;
                        if (count_hodoV_1351==0) {
                            if (0<=(gch2-1344) && (gch2-1344)<16)  hodoV_select->Fill(gch2-1344); // fill vertical fiber
                        }
                    }

                }
                if (gch==1352 && bad_bool_hodoV==false && bad_bool_hodoH== false) {  // Vertical
                    hodo_bool_1352=true;
                    for(Int_t h5=0; h5<fEvent->GetHitSize(); h5++) {
                        fHit5 = fEvent->GetHit(h5);
                        Int_t gch5=fHit5.GetChannel();
                        if (gch5==1344 ) ++count_hodoV_1352;
                        if (gch5==1345 ) ++count_hodoV_1352;
                        if (gch5==1346 ) ++count_hodoV_1352;
                        if (gch5==1347 ) ++count_hodoV_1352;
                        if (gch5==1348 ) ++count_hodoV_1352;
                        if (gch5==1349 ) ++count_hodoV_1352;
                        if (gch5==1350 ) ++count_hodoV_1352;
                        if (gch5==1351 ) ++count_hodoV_1352;
                        if (gch5==1353 ) ++count_hodoV_1352;
                        if (gch5==1354 ) ++count_hodoV_1352;
                        if (gch5==1355 ) ++count_hodoV_1352;
                        if (gch5==1356 ) ++count_hodoV_1352;
                        if (gch5==1357 ) ++count_hodoV_1352;
                        if (gch5==1358 ) ++count_hodoV_1352;
                        if (gch5==1359 ) ++count_hodoV_1352;
                        if (count_hodoV_1352==0) {
                            if (0<=(gch5-1344) && (gch5-1344)<16)  hodoV_select->Fill(gch5-1344); // fill vertical fiber
                        }
                    }
                }
                /////////////////////////////////////////////////////
                if (gch==1368 && bad_bool_hodoV==false && bad_bool_hodoH== false) { // Horizental
                    hodo_bool_1368=true;
                    for(Int_t h6=0; h6<fEvent->GetHitSize(); h6++) {
                        fHit6 = fEvent->GetHit(h6);
                        Int_t gch6=fHit6.GetChannel();
                        if (gch6==1360 ) ++count_hodoH_1368;
                        if (gch6==1361 ) ++count_hodoH_1368;
                        if (gch6==1362 ) ++count_hodoH_1368;
                        if (gch6==1363 ) ++count_hodoH_1368;
                        if (gch6==1364 ) ++count_hodoH_1368;
                        if (gch6==1365 ) ++count_hodoH_1368;
                        if (gch6==1366 ) ++count_hodoH_1368;
                        if (gch6==1367 ) ++count_hodoH_1368;
                        if (gch6==1369 ) ++count_hodoH_1368;
                        if (gch6==1370 ) ++count_hodoH_1368;
                        if (gch6==1371 ) ++count_hodoH_1368;
                        if (gch6==1372 ) ++count_hodoH_1368;
                        if (gch6==1373 ) ++count_hodoH_1368;
                        if (gch6==1374 ) ++count_hodoH_1368;
                        if (count_hodoH_1368==0) {
                            if (16<=(gch6-1344) && (gch6-1344)<32) hodoH_select->Fill(gch6-1344-16);
                        }
                    }
                }
                if (gch==1369 && bad_bool_hodoV==false && bad_bool_hodoH== false) { // Horizental
                    hodo_bool_1369=true;
                    for(Int_t h3=0; h3<fEvent->GetHitSize(); h3++) {
                        fHit3 = fEvent->GetHit(h3);
                        Int_t gch3=fHit3.GetChannel();
                        if (gch3==1360 ) ++count_hodoH_1369;
                        if (gch3==1361 ) ++count_hodoH_1369;
                        if (gch3==1362 ) ++count_hodoH_1369;
                        if (gch3==1363 ) ++count_hodoH_1369;
                        if (gch3==1364 ) ++count_hodoH_1369;
                        if (gch3==1365 ) ++count_hodoH_1369;
                        if (gch3==1366 ) ++count_hodoH_1369;
                        if (gch3==1367 ) ++count_hodoH_1369;
                        if (gch3==1368 ) ++count_hodoH_1369;
                        if (gch3==1370 ) ++count_hodoH_1369;
                        if (gch3==1371 ) ++count_hodoH_1369;
                        if (gch3==1372 ) ++count_hodoH_1369;
                        if (gch3==1373 ) ++count_hodoH_1369;
                        if (gch3==1374 ) ++count_hodoH_1369;
                        if (count_hodoH_1369==0) {
                            if (16<=(gch3-1344) && (gch3-1344)<32) hodoH_select->Fill(gch3-1344-16);
                        }
                    }

                }
                if (gch==1370 && bad_bool_hodoV==false && bad_bool_hodoH== false) { // Horizental
                    hodo_bool_1370=true;
                    for(Int_t h7=0; h7<fEvent->GetHitSize(); h7++) {
                        fHit7 = fEvent->GetHit(h7);
                        Int_t gch7=fHit7.GetChannel();
                        if (gch7==1360 ) ++count_hodoH_1370;
                        if (gch7==1361 ) ++count_hodoH_1370;
                        if (gch7==1362 ) ++count_hodoH_1370;
                        if (gch7==1363 ) ++count_hodoH_1370;
                        if (gch7==1364 ) ++count_hodoH_1370;
                        if (gch7==1365 ) ++count_hodoH_1370;
                        if (gch7==1366 ) ++count_hodoH_1370;
                        if (gch7==1367 ) ++count_hodoH_1370;
                        if (gch7==1368 ) ++count_hodoH_1370;
                        if (gch7==1369 ) ++count_hodoH_1370;
                        if (gch7==1371 ) ++count_hodoH_1370;
                        if (gch7==1372 ) ++count_hodoH_1370;
                        if (gch7==1373 ) ++count_hodoH_1370;
                        if (gch7==1374 ) ++count_hodoH_1370;
                        if (count_hodoH_1370==0) {
                            if (16<=(gch7-1344) && (gch7-1344)<32) hodoH_select->Fill(gch7-1344-16);
                        }
                    }
                }

            }// end of hit loop
            countmulti_tof1->Fill(count_tof1);
            countmulti_tof2->Fill(count_tof2);
            countmulti_trg1->Fill(count_trg1);
            countmulti_trg2->Fill(count_trg2);
            countmulti_trgmzH->Fill(count_trgmzH);
            countmulti_trgmzV->Fill(count_trgmzV);
            ///////
            countmulti_hodoV->Fill(count_hodoV_1351);
            countmulti_hodoH->Fill(count_hodoH_1369);
            // hodo all
            {
                for(int h=0; h<16; h++)
                    for(int v=0; v<16; v++)
                        hodoF->Fill(v,(15-h),hodoH->GetBinContent(h)+hodoV->GetBinContent(v));
            } // end Fill

            // hodo Vertical
            {
                for(int h=0; h<16; h++)
                    for(int v=0; v<16; v++)
                        hodo_multi_withmedVfiber->Fill(v,(15-h),hodoV_select->GetBinContent(v));
            }
            // hodo Horizental
            {
                for(int h=0; h<16; h++)
                    for(int v=0; v<16; v++)
                        hodo_multi_withmedHfiber->Fill(v,(15-h),hodoH_select->GetBinContent(h));
            }


            //////////////////////////////////////////////////////////////////////////////////////////
            delta_tof2tof1=LE_tof2-LE_tof1;

            //htof1_le->Fill(LE_tof1);
            //htof1_tot->Fill(TOT_tof1);
            htof1_le_tot->Fill(LE_tof1,TOT_tof1);
            //htof2_le->Fill(LE_tof2);
            //htof2_tot->Fill(TOT_tof2);
            htof2_le_tot->Fill(LE_tof2,TOT_tof2);
            //htrg1_le->Fill(LE_trg1);
            //htrg1_tot->Fill(TOT_trg1);
            htrg1_le_tot->Fill(LE_trg1, TOT_trg1);
            //htrg2_le->Fill(LE_trg2);
            //htrg2_tot->Fill(TOT_trg2);
            htrg2_le_tot->Fill(LE_trg2, TOT_trg2);
            //htrgmzV_le->Fill(LE_trgmzV);
            //htrgmzV_tot->Fill(TOT_trgmzV);
            htrgmzV_le_tot->Fill(LE_trgmzV, TOT_trgmzV);
            //htrgmzH_le->Fill(LE_trgmzH);
            //htrgmzH_tot->Fill(TOT_trgmzH);
            htrgmzH_le_tot->Fill(LE_trgmzH, TOT_trgmzH);

            hdelta_tof2tof1->Fill(delta_tof2tof1);

            if(/*hodo_bool_1351==true  &&  hodo_bool_1369==true &&
                count_hodoV_1351==0 && count_hodoH_1369==0 &&*/  // tight cut


                bad_bool_hodoV==false && bad_bool_hodoH== false &&
                (hodo_bool_1350==true || hodo_bool_1351==true|| hodo_bool_1352==true) && (hodo_bool_1368==true || hodo_bool_1369==true || hodo_bool_1370==true ) &&
                count_hodoV_1350 ==0 && count_hodoV_1351==0 && count_hodoV_1352==0 && count_hodoH_1368==0 && count_hodoH_1369==0 && count_hodoH_1370==0 &&




                tof1_bool == true && tof2_bool == true &&
                trg1_bool==true && trg2_bool==true && trgmzH_bool==true && trgmzV_bool==true &&

                LE_trg1>-144 && LE_trg1 < -132 &&
                TOT_trg1>118.4 && TOT_trg1<119.8 &&

                LE_tof1>-282 && LE_tof1<-270 &&
                TOT_tof1>45 && TOT_tof1<50 &&

                LE_tof2>-250 && LE_tof2 < -238 &&
                TOT_tof2>44 && TOT_tof2<50 &&

                LE_trg2> -125 && LE_trg2<-110 &&
                LE_trg1> -143 && LE_trg1<-130 &&

                multi_bool ==false ) { // TOF1, TOF2, TRG1

                // All good candidates
                //hdelta_tof2tof1->Fill(delta_tof2tof1);
                //goto goon;

                if (fMethod==2 && delta_tof2tof1 > 32.4 && delta_tof2tof1< 32.9 && recoP==1 && recoPi==0) goto goon; //proton Candidate
                if (fMethod==2 && delta_tof2tof1 > 31.5 && delta_tof2tof1< 32.0 && recoP==0 && recoPi==1 ) goto goon; //pi Candidate
                //std::cout<<"No Problem  methos =  " <<fMethod<<std::endl;
                if (fMethod==3 && fEvent->GetType()==0 && ((delta_tof2tof1 > 32.4 && delta_tof2tof1< 32.9) || (delta_tof2tof1 > 31.5 && delta_tof2tof1< 32.0)) ) goto goon; //with tof cuts
                //if (fMethod==3 && fEvent->GetType()==0 ) goto goon; //with tof cuts
            }
            continue; // outside the event loop
        }// if data condition

goon:
        std::cout<<"No Problem  goon  " <<std::endl;
        hdelta_tof2tof1_isproton->Fill(delta_tof2tof1);
        nHits_dac_syscut_p->Fill(nHits);
        // hodo after cuts
        {
            for(int h=0; h<16; h++)
                for(int v=0; v<16; v++)
                    hodo_afterCut->Fill(v,(15-h),hodoH_select->GetBinContent(h)+hodoV_select->GetBinContent(v));
        } // end FillF()

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

        Int_t nGoodPhotons=0;
        Int_t nGoodPhotons_TimeCutOnly=0;
        Int_t nCandidat=0;
        Int_t nGoodPhotons_true_sim=0;
        //std::cout<<"######### Hit loop will start"<<std::endl;
        Int_t counter_hit_loop=0;
        for(Int_t h=0; h<nHits; h++) {
            fHit = fEvent->GetHit(h);
            hitTime = fHit.GetLeadTime();
            //std::cout<<"hitTime  "<<hitTime <<std::endl;
            hit_b4_phs->Fill(hitTime);
            if(fEvent->GetType()!=0) hitTime+=fRand.Gaus(0,test1); // time resol. in case it was not simulated 200 ps 250 switch of smearing
            photonEnergy= fHit.GetMomentum().Mag()*1E6;
            if(fEvent->GetType()==0) {// shift hit time for data only
                if(prtangle==20) hitTime=hitTime-0.189155;
                if(prtangle==30) hitTime=hitTime+0.030525;
                if(prtangle==40) hitTime=hitTime+0.196482;
                if(prtangle==50) hitTime=hitTime+0.149906;
                if(prtangle==60) hitTime=hitTime+0.081199;
                if(prtangle==70) hitTime=hitTime-0.102829;
                if(prtangle==80) hitTime=hitTime-0.228011;
                if(prtangle==90) hitTime=hitTime+0.038770;
                if(prtangle==100) hitTime=hitTime-0.182761;
                if(prtangle==110) hitTime=hitTime-0.189065;
                if(prtangle==120) hitTime=hitTime+0.050843;
                if(prtangle==130) hitTime=hitTime+0.222726;
                if(prtangle==140) hitTime=hitTime+0.235055;
                if(prtangle==150) hitTime=hitTime+0.019390;
            }
            //======================================== dynamic cuts for sim and data
            if(fEvent->GetType()==1 || fEvent->GetType()==0) {
                Double_t cut1(0);
                {   //time cuts
                    if(prtangle>19.0 && prtangle<21.0) {
                        if(hitTime<9.0 || hitTime>25 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle>24.0 && prtangle<26.0 ) {
                        if(hitTime<9.5 || hitTime>30 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle>29.0 && prtangle<31.0 ) {
                        if(hitTime<9.5 || hitTime>35 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle>34.0 && prtangle<36.0 ) {
                        if(hitTime<9.5 || hitTime>35 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle>39.0 && prtangle<41.0 ) {
                        if(hitTime<10 || hitTime>40 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle>44.0 && prtangle<46.0 ) {
                        if(hitTime<10 || hitTime>43 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle>49.0 && prtangle<51.0 ) {
                        if(hitTime<10 || hitTime>50 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle>54.0 && prtangle<56.0 ) {
                        if(hitTime<10 || hitTime>50 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle>59.0 && prtangle<61.0 ) {
                        if(hitTime<10.5 || hitTime>50 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle>64.0 && prtangle<66.0 ) {
                        if(hitTime<11 || hitTime>50 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle>69.0 && prtangle<71.0 ) {
                        if(hitTime<11 || hitTime>50 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle>74.0 && prtangle<76.0 ) {
                        if(hitTime<11 || hitTime>50 ) continue;
                        reflected = kTRUE;
                    }
                    else if(prtangle>79.0 && prtangle<81.0 ) {
                        if(hitTime<11 || hitTime>50 ) continue;
                        reflected = kTRUE;
                    }
                    else if (prtangle>84.0 && prtangle<86.0 ) {
                        if(hitTime< 3.0 || hitTime>50.0 ) continue;
                        //if(hitTime< 13.5 && hitTime>10.0 ) continue;
                        if(hitTime<=12)   reflected = kFALSE;
                        if(hitTime>12) reflected = kTRUE;
                    }
                    else if (prtangle>89.0 && prtangle<91.0) {
                        if(hitTime< 1.0 || hitTime>50.0 ) continue;
                        //if(hitTime< 17.0 || hitTime>30.0 ) continue;
                        //if(hitTime< 1.0 || hitTime>4.0 ) continue;


                        //if(hitTime< 3.0 || hitTime>50.0 ) continue;
                        //if(hitTime< 14.0 && hitTime>10.0 ) continue;
                        if(hitTime<=12)   reflected = kFALSE;
                        if(hitTime>12) reflected = kTRUE;
                    }
                    else if(prtangle>94) {
                        //if(hitTime<3 || hitTime>25) continue;
                        reflected = kFALSE;
                    }
                }
            }
            //======================================== End of dynamic cuts
            // Double_t radiatorL = (radiator==2)? 1224.9 : 1250; //plate : bar
            Double_t radiatorL = (radiator==2)? 1224.9 : 1200; //plate : bar // changed
            if( fEvent->GetType()==1) {
                lenz = radiatorL/2.-fHit.GetPosition().Z();
            } else {
                lenz = fHit.GetPosition().Z();
            }// end of else

            if(fVerbose==3) {
                TVector3 cd = fHit.GetMomentum();
                fHist5->Fill(cd.Theta()*TMath::Sin(cd.Phi()),cd.Theta()*TMath::Cos(cd.Phi()));
            }
            //std::cout<<"###hitTime2  "<<hitTime <<std::endl;
            Int_t pixid=fHit.GetPixelId()-1;
            Int_t mcpid=fHit.GetMcpId();
            if(reflected) lenz = 2*radiatorL - lenz;
            Int_t ch = map_mpc[mcpid][pixid];
            //if(prt_isBadChannel(ch)) continue;
            Int_t sensorId = 100*mcpid+fHit.GetPixelId();
            Int_t sensorId_phs = 100*mcpid+pixid;
            //std::cout<< "sensorId = "<< sensorId << "      GetPixelId = "<<fHit.GetPixelId()  <<"    mcpid= " <<mcpid <<std::endl;
            //if(sensorId==1) continue;
            Bool_t isGoodHit(false);
            Bool_t isGoodHit_true_sim(false);
            Bool_t isGoodHit_TimeCutOnly(false);
            Bool_t isCandidat(false);
            //  if(radiator==2) isGoodHit=true;
            Int_t size =fLutNode[sensorId]->Entries();
            Double_t min_time = 1000;
            Double_t min_diff_time_tangle = 0;

            Int_t photon_ambiguity_counter_wo=0;
            Int_t photon_ambiguity_counter_wt=0;
            Int_t photon_ambiguity_counter_wtc=0;
            // start of lut loop



            // top2
            //if((pixid!=0 && mcpid!=0)||(pixid!=63 && mcpid!=11)) continue;

            //cout<< "@@@@@@@@@@@@@@@@@@  fill = "<<pixid%8<<"	"<<pixid/8<<endl;
            //if ( mcpid >= 11 && pixid> 50 ) continue;
            //            prt_hdigi[mcpid]->Fill(pixid%8, pixid/8);
            momInBar_phs=TVector3(0,0,-1);
            Double_t ref_point =9.8;
            //cout<<"sensorId="<<sensorId<<"   "<<"mcpid="<<mcpid<<"  "<<"pixID="<<fHit.GetPixelId()<<endl;
            Int_t size_phs =fLutNode_phs[sensorId_phs]->Entries();
            //if(pixid==0 && mcpid==0)cout<< "@@@@@@@@@@@@@@@@@@  size_phs first= "<<size_phs<<endl;
            //if(pixid==63 && mcpid==11)cout<< "@@@@@@@@@@@@@@@@@@  size_phs last= "<<size_phs<<endl;
            //if (size_phs == 0) cout<<"sensorId="<<sensorId<<"   "<<"mcpid="<<mcpid<<"  "<<"pixID="<<fHit.GetPixelId()<<endl;
            //if (size_phs == 0) continue;
            //std::cout<< "sensorId = "<< sensorId << "      GetPixelId = "<<fHit.GetPixelId()  <<"    mcpid= " <<mcpid <<std::endl;
            Int_t phs_solution_counter=0;
            for(Int_t i=0; i<size_phs; i++) { // size_phs  5000
                weight = 1; //fLutNode[sensorId]->GetWeight(i);
                dird   = fLutNode_phs[sensorId_phs]->GetEntry(i).Unit();
                time_phs = fLutNode_phs[sensorId_phs]->GetTime(i);//
                if (dird.X() > 0 ) continue; // for 20 deg
                if (time_phs<ref_point) continue;

                //if (dird.X() > 0 && hitTime> ref_point) continue;
                //if (dird.X() < 0 && hitTime< ref_point) continue;

                //prt_hdigi[mcpid]->Fill(pixid%8, pixid/8);
                pos_x = fLutNode_phs[sensorId_phs]->GetHitPos(i).X();
                pos_y = fLutNode_phs[sensorId_phs]->GetHitPos(i).Y();
                //if(fabs(pos_x)>1) continue;
                //if(fabs(pos_y)>3) continue;
                //hist_phs_xy->Fill(pos_x,pos_y);
                //cout<< "@@@@@@@@@@@@@@@@@@  time_phs= "<<time_phs<<" dirx= " <<dird.X()<<endl;
                tangle = momInBar_phs.Angle(dird);
                //if(tangle > TMath::PiOver2()) tangle = TMath::Pi()-tangle;
                if(fabs(tangle-fAngleP)<0.04) {
                    fHist0->Fill(hitTime-time_phs);
                    fHist0i->Fill(time_phs);
                    fHist0i_bg->Fill(hitTime);
                    //std::cout<< "##########################  diff =  "<< hitTime - time_phs<<std::endl;
                }
                fcut->SetParameters(0.9,8.5,0.5);
                fcut2->SetParameters(3.5,9.5,0.5);
                //if (tangle > fcut->Eval(time_phs)) continue;
                //if (tangle < fcut2->Eval(time_phs)) continue;
                hist_time_angle_all->Fill(time_phs,tangle);
                //hist_time_angle_all_correlation->Fill(hitTime - time_phs, tangle);
                if(fabs(hitTime - time_phs)>test2 ) continue;
                hist_time_angle_all_tcut->Fill(time_phs,tangle);
                //std::cout<< "##########################  tangle "<< tangle<<std::endl;
                fHist_correction->Fill(tangle);
                // prt_hdigi[mcpid]->Fill(pixid%8, pixid/8);
                if(gPDF ==1) fHistCh_2D[ch]->Fill(time_phs,tangle); // good after time cut
                //if(true) {
                if(method_type == 3 && tangle>0.6 && tangle<1.0 && gPDF ==2) { // 0.6 1.0
                    // use graphs
                    // use histograms
                    Int_t testt1 = fHistCh_2D_read_p[ch]->GetEntries();
                    Int_t testt2 = fHistCh_2D_read_pi[ch]->GetEntries();
                    Int_t xkp = fHistCh_2D_read_p[ch]->GetXaxis()->FindBin(hitTime);
                    Int_t ykp = fHistCh_2D_read_p[ch]->GetYaxis()->FindBin(tangle);
                    Int_t xkpi = fHistCh_2D_read_pi[ch]->GetXaxis()->FindBin(hitTime);
                    Int_t ykpi = fHistCh_2D_read_pi[ch]->GetYaxis()->FindBin(tangle);
                    //std::cout<< "##########################  testt1 "<< testt1<<"  xkp="<<xkp<<"  "<< xkp<<"content "<< fHistCh_2D_read_p[ch]->GetBinContent(xkp,ykp)<<std::endl;
                    //std::cout<< "##########################  testt2 "<< testt2<<"  xkpi="<<xkpi<<" "<< xkpi<<" content "<< fHistCh_2D_read_pi[ch]->GetBinContent(xkpi,ykpi)<<std::endl;

                    //sum1 += TMath::Log(fHistCh_2D_read_p[ch]->GetBinContent(xkp,ykp));
                    //sum2 += TMath::Log(fHistCh_2D_read_pi[ch]->GetBinContent(xkpi,ykpi));
                    Double_t content_pi = fHistCh_2D_read_pi[ch]->GetBinContent(xkpi,ykpi);
                    Double_t content_p  = fHistCh_2D_read_p[ch]->GetBinContent(xkp,ykp);

                    if(content_p !=0 )sum1 += TMath::Log(content_p);
                    if(content_pi!=0 )sum2 += TMath::Log(content_pi);

                    //std::cout<< "##########################  sum= "<< sum1<<"	"<<sum2<<"  "<<std::endl;
                    //std::cout<<"No Problem  separation  " <<kp<<" "<<kp<<std::endl;
                }



                // standard
                if(method_type == 3 && tangle>0.6 && tangle<1.0 && gPDF ==3) {
                    // use standared
                    sum1 += TMath::Log(gF1->Eval(tangle)+noise);
                    sum2 += TMath::Log(gF2->Eval(tangle)+noise);
                }



                ++phs_solution_counter;
                //std::cout<< "##########################  phs_solution_counter "<< phs_solution_counter<<std::endl;
            }// end of phs lut loop
            hist_phs_solution_number->Fill(phs_solution_counter);

            histo_photon_ambiguity_wo->Fill(photon_ambiguity_counter_wo,weight);
            histo_photon_ambiguity_wt->Fill(photon_ambiguity_counter_wt,weight);
            histo_photon_ambiguity_wtc->Fill(photon_ambiguity_counter_wtc,weight);
            /*
             // closest approach method
             //fHist_copy->Fill(min_diff_time_tangle ,weight);
             if(gPDF ==1) fHistCh[ch]->Fill(min_diff_time_tangle ,weight);
             if(method_type == 3 && min_diff_time_tangle>0.65 && min_diff_time_tangle<0.95 && isGoodHit_TimeCutOnly && isGoodHit) {
             sum1 += TMath::Log(fHistCh_graph_p[ch]->Eval(min_diff_time_tangle));
             sum2 += TMath::Log(fHistCh_graph_pi[ch]->Eval(min_diff_time_tangle));
             solution_number_approach_selection++;
             //sum1 += TMath::Log(gF1->Eval(tangle)+noise);
             //sum2 += TMath::Log(gF2->Eval(tangle)+noise);
             //std::cout<<"No Problem  separation  " <<sum1<<" "<<sum2<<std::endl;
             }
             */
            //std::cout<<"No Problem  lut loop  " <<std::endl;
            if(isGoodHit) {
                nsHits++;
                //prt_hdigi[mcpid]->Fill(pixid%8, pixid/8);
                nGoodPhotons++;
            }
            if(isGoodHit_true_sim)nGoodPhotons_true_sim++;
            if(isGoodHit_TimeCutOnly) nGoodPhotons_TimeCutOnly++;
            if(isCandidat) nCandidat++;

            ++counter_hit_loop;
            //std::cout<<" b4 hit loop  " <<counter_hit_loop<<std::endl;
        }// end of hit loop
        //std::cout<<" after hit loop  " <<std::endl;
        fnHits_p_good->Fill(nGoodPhotons);
        fnHits_p->Fill(nGoodPhotons_TimeCutOnly);
        fnHits->Fill(nCandidat);
        fnHits_true_sim->Fill(nGoodPhotons_true_sim);
        //std::cout<<"@@@@@@@@@@@"<<" nGoodPhotons= "<<nGoodPhotons<<std::endl;
        // for(Int_t j=0; j<prt_nmcp; j++){
        //   for(Int_t i=0; i<65; i++){
        // 	mcpdata[j][i]=0;
        // 	cluster[j][i]=0;
        //   }
        // }

        Double_t sum = sum1-sum2;
        if(sum!=0) {
            if(tofPid==2212) hLnDiffP->Fill(sum);
            if(tofPid==211) hLnDiffPi->Fill(sum);
            //std::cout<<"@@@@@@@@@@@@@ tofPid "<<tofPid <<std::endl;
            likelihood=sum;
            std::cout<<"@@@@@@@@@@@@@ tofPid "<<tofPid << "    likelihood="<< likelihood<<std::endl;
        }


        std::cout<<"No Problem  likelihood  " <<std::endl;
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
                //test1 = fTest;
                distPid = FindPdg(momentum,cangle);
                nph = nsHits/(Double_t)ninfit;
                //std::cout<<"@@@@@@@@@@@@@ nph "<<nph <<std::endl;
                spr = spr*1000;
                //std::cout<<"@@@@@@@@@@@@@ spr "<<spr <<std::endl;
                trr = spr/sqrt(nph);
                //std::cout<<"@@@@@@@@@@@@@ trr "<<trr <<std::endl;
                theta = fEvent->GetAngle(); // here
                // here
                //std::cout<<"@@@@@@@@@@@@@ theta "<<theta <<std::endl;
                par3 = fEvent->GetTest1();
                //std::cout<<"@@@@@@@@@@@@@ par3 "<<par3 <<std::endl;
                tree.Fill(); // alarma
                //std::cout<<"no problem 3 1"<<std::endl;
            }
            //std::cout<<"no problem 3 2"<<std::endl;
            ResetHists();
            nsHits=0;
            //std::cout<<"no problem 3 3"<<std::endl;
        }
        //std::cout<<"no problem 3 4"<<std::endl;
        if(++nsEvents>=end) break;
        //std::cout<<"no problem 3 5"<<std::endl;

        ++counter_event_loop;
    }// end of event loop
    //std::cout<<"no problem 3 6"<<std::endl;
    /*
        prt_drawDigi("m,p,v\n",2017,0,0);
        prt_cdigi->SetName(Form("hp_dataProtonS332_%d_%2.1f",(Int_t)prt_theta,prt_phi));
        prt_canvasAdd(prt_cdigi);
        prt_cdigi_palette->Draw();
        prt_canvasSave(1,0);
    */

    if(fMethod==2) {
        FindPeak(cangle,spr, prtangle);
        nph = nsHits/(Double_t)nsEvents;
        std::cout<<"@@@@@@@@@@@@@  nsEvents ="<<nsEvents <<std::endl;
        std::cout<<"@@@@@@@@@@@@@  nsHits ="<<nsHits <<std::endl;
        std::cout<<"@@@@@@@@@@@@@  nph ="<<nph <<std::endl;
        spr = spr*1000;
        trr = spr/sqrt(nph);
        theta = fEvent->GetAngle(); // here
        // here
        //theta = 60.0; // here
        par3 = fEvent->GetTest1();
        if(fVerbose) std::cout<<Form("SPR=%2.2F N=%2.2f",spr,nph)<<std::endl;
        tree.Fill(); // alarma 2
    } else {


        TString nid = Form("_%2.0f", prtangle);

        if(fVerbose<2) gROOT->SetBatch(1);
        prt_canvasAdd("r_lhood"+nid,800,400);
        prt_normalize(hLnDiffP,hLnDiffPi);
        hLnDiffP->SetLineColor(2);
        TF1 *ff;
        Double_t m1(1),m2(1),s1(1),s2(1);
        std::cout<<"no problem 3 7"<<std::endl;
        if(hLnDiffP->GetEntries()>10) {
            hLnDiffP->Fit("gaus","S");
            ff = hLnDiffP->GetFunction("gaus");
            m1=ff->GetParameter(1);
            s1=ff->GetParameter(2);
        }
        std::cout<<"no problem 3 8"<<std::endl;
        if(hLnDiffPi->GetEntries()>10) {
            hLnDiffPi->Fit("gaus","S");
            ff = hLnDiffPi->GetFunction("gaus");
            m2=ff->GetParameter(1);
            s2=ff->GetParameter(2);
        }
        std::cout<<" hLnDiffP->GetEntries()  "<< hLnDiffP->GetEntries() <<std::endl;
        std::cout<<" hLnDiffPi->GetEntries()  "<< hLnDiffPi->GetEntries() <<std::endl;
        separation = (fabs(m2-m1))/(0.5*(s1+s2));
        std::cout<<"################ separation "<< separation <<std::endl;
        gStyle->SetOptFit(0);
        gStyle->SetOptStat(0);
        //hLnDiffP->SetName(Form("s_%2.2f",separation));
        hLnDiffP->Draw();
        hLnDiffPi->SetLineColor(4);
        hLnDiffPi->Draw("same");



        ///////////////////////////////////////////////////////////////////////
        /*
        prt_canvasAdd("r_delta"+nid,800,400);
        hdelta_tof2tof1->SetTitle(Form("theta %3.1f", prtangle));
        hdelta_tof2tof1->SetStats(0);
        hdelta_tof2tof1->SetLineColor(kRed);
        hdelta_tof2tof1_isproton->SetLineColor(kBlue);
        hdelta_tof2tof1->GetXaxis()->SetTitle("LE TOF2 - TOF1 [ns]");
        hdelta_tof2tof1->GetYaxis()->SetTitle("count[#]");
        hdelta_tof2tof1->Draw();
        hdelta_tof2tof1_isproton->Draw("same");
        */
        ///////////////////////////////////////////////////////////////////////


        //prt_waitPrimitive("r_test_hist"+nid);// wait here
        prt_canvasSave(2,0);
        //prt_waitPrimitive("r_lhood","w");
        if(fVerbose) gROOT->SetBatch(0);
        tree.Fill(); // alarma
    }
    //std::cout<<"solution_number_approach_selection= "<<solution_number_approach_selection<<std::endl;
    //std::cout<<"solution_number= "<<solution_number<<std::endl;
    tree.Write();
    std::cout<<"wrtitten "<<std::endl;

    //    for(Int_t i=0; i<960; i++) {
    //        hit_time_ch_ccut[i]->Write();
    //        hit_time_phs_ch_ccut[i]->Write();
    //        hit_time_diff_ch_ccut[i]->Write();
    //
    //        hit_time_ch[i]->Write();
    //        hit_time_phs_ch[i]->Write();
    //        hit_time_diff_ch[i]->Write();
    //    }

    fHist_correction->Write();
    fHist0->Write();
    fHist0i->Write();
    fHist0i_bg->Write();
    hit_b4_phs->Write();

    hist_time_angle_all->Write();
    hist_time_angle_all_tcut->Write();
    //hist_time_angle_all_correlation->Write();
    hist_phs_solution_number->Write();
    //hist_phs_xy->Write();
    //hdelta_tof2tof1->Write();
    //hdelta_tof2tof1_isproton->Write();
    //hist_time_angle_all_cut->Write();
    //hist_time_angle_all_ok->Write();
    /*
        for(Int_t i=0; i<960; i++) {
              hist_time_angle_ch[i]->Write();
    	  hist_time_angle_ch2[i]->Write();
              hist_time_angle_ch3[i]->Write();
         }
    */

    if(gPDF ==1) {
        for(Int_t i=0; i<800; i++) {
            fHistCh_2D[i]->Write();
        }
    }

    if(gPDF ==2) {

        hLnDiffP->Write();
        hLnDiffPi->Write();
    }



    //    lut_pix_pos_xy->Write();
    //    hist_ambiguity->Write();
    //    fHist->Write();
    //    fHist_copy->Write();
    //
    //
    //    fHist_same_path->Write();
    //    fHist_bg->Write();
    //    fnHits->Write();
    //    fnHits_true_sim->Write();
    //    fnHits_p->Write();
    //    fnHits_p_good->Write();
    //
    //    fHist1->Write();
    //    fHist2->Write();
    //    nHits_dac_syscut_p->Write();
    //    nHits_dac->Write();
    //    hdelta_tof2tof1->Write();
    //
    //    fHist_same_path_wotc->Write();
    //
    //    histo_photon_ambiguity_wo->Write();
    //    histo_photon_ambiguity_wt->Write();
    //    histo_photon_ambiguity_wtc->Write();
    //


    //    for(Int_t i=0; i<prt_nmcp; i++) {
    //        fHistMcp[i]->Write();
    //    }
    //    for(Int_t i=0; i<prt_nmcp; i++) {
    //        fHistMcp_same_path[i]->Write();
    //    }
    file.Write();
    std::cout<<"no problem 3 13"<<std::endl;
    ///if(fVerbose) ResetHists();
}// end of the function

Int_t g_num =0;
Bool_t PrtLutReco::FindPeak(Double_t& cangle, Double_t& spr, Double_t a, Int_t tofpdg) {
    cangle=0;
    spr=0;
    //  gStyle->SetCanvasPreferGL(kTRUE);
    if(fHist->GetEntries()>20 || fHistPi->GetEntries()>20) {
        gROOT->SetBatch(1);
        Int_t nfound = fSpect->Search(fHist,1,"",0.9); //0.6
        if(nfound>0) cangle = fSpect->GetPositionX()[0];
        else cangle =  fHist->GetXaxis()->GetBinCenter(fHist->GetMaximumBin());
        cangle =  fHist->GetXaxis()->GetBinCenter(fHist->GetMaximumBin());
        if(cangle>0.85) cangle=0.82;
        fFit->SetParameters(100,cangle,0.010);
        fFit->SetParNames("p0","#theta_{c}","#sigma_{c}","p3","p4");
        fFit->SetParLimits(0,0.1,1E6);
        fFit->SetParLimits(1,cangle-0.04,cangle+0.04); // changed 0.04
        fFit->SetParLimits(2,0.005,0.018); // width 7-10 // changed 0.014
        //fFit->FixParameter(2,0.01);
        //fFit->FixParameter(3,0);
        //fFit->FixParameter(4,0);
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
            // Cherenkov correction
            if (false) {
                fFit->SetParLimits(2,0.004,0.008); // width 7-10
                for(Int_t i=0; i<prt_nmcp; i++) {
                    prt_canvasAdd(Form("r_tangle_%d",i),800,400);
                    fHistMcp[i]->Fit("fgaus","lq","",fAngleP-0.03,fAngleP+0.03);
                    std::cout<<"if(mcpid=="<< i<<") tangle += "<<fAngleP-fFit->GetParameter(1)<<";" <<std::endl;
                    fHistMcp[i]->Draw();
                    prt_waitPrimitive("r_tangle_"+i);// wait here
                }
            }
            /*
             fFit->SetParLimits(2,0.004,0.008); // width 7-10
             for(Int_t i=0; i<prt_nmcp; i++){
             prt_canvasAdd(Form("r_tangle_%d",i),800,400);
             fHistMcp[i]->Fit("fgaus","lq","",fAnglePi-0.03,fAnglePi+0.03);
             std::cout<<"if(mcpid=="<< i<<") tangle += "<<fAnglePi-fFit->GetParameter(1)<<";" <<std::endl;
             fHistMcp[i]->Draw();
             }
             */
            // for(Int_t i=0; i<960; i++){
            // 	prt_canvasAdd(Form("r_tangle_ch_%d",i),800,400);
            // 	fHistCh[i]->Fit("fgaus","lq","",fAngleP-0.03,fAngleP+0.03);
            // 	std::cout<<"if(ch=="<< i<<") tangle += "<<fAngleP-fFit->GetParameter(1)<<";" <<std::endl;
            // 	fHistCh[i]->Draw();
            // }
            //      TString name = Form("r_tangle_%3.1f",PrtManager::Instance()->GetTest3());
            TString nid = Form("_%2.0f",a);

            //~
            //~ for(Int_t pix=495; pix<500; pix++) {
            //~ prt_canvasAdd(Form("r_pix_%d",pix),800,400);
            //~ fHistCh_graph_p[pix]->Draw();
            //~ prt_canvasGet(Form("r_pix_%d",pix))->Update();
            //~ TLine *lin_ch_p_v = new TLine(0,0,0,1000);
            //~ lin_ch_p_v->SetX1(fAngleP);
            //~ lin_ch_p_v->SetX2(fAngleP);
            //~ lin_ch_p_v->SetY1(gPad->GetUymin());
            //~ lin_ch_p_v->SetY2(gPad->GetUymax());
            //~ lin_ch_p_v->SetLineColor(kRed);
            //~ lin_ch_p_v->Draw();
            //~
            //~ }


            /*
             ///////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_tangle"+nid,800,400);
             //fHist->SetTitle(Form("theta %3.1f , TOF PID = %d", a, tofpdg));
             fHist->SetTitle(Form("theta %3.1f , proton", a));
             fHist->SetMinimum(0);
             //fHist->Scale(1/fHist->GetMaximum());
             //prt_normalize(fHist,fHistPi);
             prt_normalize(fHist,fHist_same_path);

             fHistPi->SetLineColor(2);
             fHist_same_path->SetLineColor(2);
             fHist->Draw();
             fHist_same_path->Draw("same");
             fHist_bg->Draw("same");
             fHistPi->Draw("same");
             // gF1->Draw("same");
             // gF2->Draw("same");
             fHisti->SetLineColor(kRed+2);
             if(fHisti->GetEntries()>5) fHisti->Draw("same");
             prt_canvasGet("r_tangle"+nid)->Update();
             TLine *line = new TLine(0,0,0,1000);
             line->SetX1(fAngleP);
             line->SetX2(fAngleP);
             line->SetY1(gPad->GetUymin());
             line->SetY2(gPad->GetUymax());
             line->SetLineColor(kRed);
             line->Draw();
             TLine *line2 = new TLine(0,0,0,1000);
             line2->SetX1(fAnglePi);
             line2->SetX2(fAnglePi);
             line2->SetY1(gPad->GetUymin());
             line2->SetY2(gPad->GetUymax());
             line2->SetLineColor(kBlue);
             line2->Draw();
             std::cout<<"fAnglePi "<< fAnglePi<<std::endl;
             std::cout<<"fAngleP "<< fAngleP<<std::endl;
             //prt_waitPrimitive("r_tangle"+nid);// wait here

             //////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_fnHits"+nid,800,400);
             fnHits->SetTitle(Form("Theta %3.1f", a));
             //fnHits->SetStats(0);
             fnHits_p->SetLineColor(kGreen);
             fnHits_p_good->SetLineColor(kRed);
             fnHits->Draw();
             fnHits_p->Draw("same");
             fnHits_p_good->Draw("same");
             //std::cout<<"@@@@@@@@@@@@@  fnHits_p_good mean  ="<< fnHits_p_good->GetMean() <<std::endl;
             ///////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_delta"+nid,800,400);
             hdelta_tof2tof1->SetTitle(Form("theta %3.1f", a));
             hdelta_tof2tof1->SetStats(0);
             hdelta_tof2tof1->SetLineColor(kRed);
             hdelta_tof2tof1_isproton->SetLineColor(kBlue);
             hdelta_tof2tof1->GetXaxis()->SetTitle("LE TOF2 - TOF1 [ns]");
             hdelta_tof2tof1->GetYaxis()->SetTitle("count[#]");
             hdelta_tof2tof1->Draw();
             hdelta_tof2tof1_isproton->Draw("same");
             ///////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_htof1_le_tot"+nid,800,400);
             htof1_le_tot->SetTitle(Form("theta %3.1f", a));
             htof1_le_tot->SetStats(0);
             htof1_le_tot->Draw("colz");
             htof1_le_tot->SetXTitle("LE [ns]");
             htof1_le_tot->SetYTitle("TOT [ns]");
             ///////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_htof2_le_tot"+nid,800,400);
             htof2_le_tot->SetTitle(Form("theta %3.1f", a));
             htof2_le_tot->SetStats(0);
             htof2_le_tot->Draw("colz");
             htof2_le_tot->SetXTitle("LE [ns]");
             htof2_le_tot->SetYTitle("TOT [ns]");
             ///////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_htrg1_le_tot"+nid,800,400);
             htrg1_le_tot->SetTitle(Form("theta %3.1f", a));
             htrg1_le_tot->SetStats(0);
             htrg1_le_tot->Draw("colz");
             htrg1_le_tot->SetXTitle("LE [ns]");
             htrg1_le_tot->SetYTitle("TOT [ns]");
             ///////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_htrg2_le_tot"+nid,800,400);
             htrg2_le_tot->SetTitle(Form("theta %3.1f", a));
             htrg2_le_tot->SetStats(0);
             htrg2_le_tot->Draw("colz");
             htrg2_le_tot->SetXTitle("LE [ns]");
             htrg2_le_tot->SetYTitle("TOT [ns]");
             ///////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_htrgmzH_le_tot"+nid,800,400);
             htrgmzH_le_tot->SetTitle(Form("theta %3.1f", a));
             htrgmzH_le_tot->SetStats(0);
             htrgmzH_le_tot->Draw("colz");
             htrgmzH_le_tot->SetXTitle("LE [ns]");
             htrgmzH_le_tot->SetYTitle("TOT [ns]");
             ///////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_htrgmzV_le_tot"+nid,800,400);
             htrgmzV_le_tot->SetTitle(Form("theta %3.1f", a));
             htrgmzV_le_tot->SetStats(0);
             htrgmzV_le_tot->Draw("colz");
             htrgmzV_le_tot->SetXTitle("LE [ns]");
             htrgmzV_le_tot->SetYTitle("TOT [ns]");
             /////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_time"+nid,800,400);
             fHist1->SetTitle(Form("theta %3.1f", a));
             fHist1->SetLineColor(2);
             fHist1->Draw();
             fHist2->Draw("same");
             ////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_diff"+nid,800,400);
             fHist0->SetTitle(Form("theta %3.1f", a));
             fHist0->Draw();
             fHist0i->SetLineColor(kRed+2);
             if(fHist0i->GetEntries()>5)  fHist0i->Draw("same");
             //////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_measured_time"+nid,800,400);
             fHist1->SetTitle(Form("theta %3.1f", a));
             fHist1->SetLineColor(2);
             fHist1->Draw();
             //fHist2->Draw("same");
             /////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_cm"+nid,800,400);
             fHist3->SetTitle(Form("theta %3.1f", a));
             fHist3->Draw("colz");
             ////////////////////////////////////////////////////////////////////
             prt_drawDigi("m,p,v\n", 2017);//2
             //TPaveText* tit;
             //tit = new TPaveText(17.0,5,25.,5,"NB");
             //tit->SetFillColor(0);
             //tit->AddText(Form("theta = %f", a));
             prt_canvasAdd(prt_cdigi);
             prt_cdigi->SetName("r_hp"+nid);
             prt_cdigi->SetTitle(Form("theta %3.1f", a));
             //prt_cdigi->Draw();
             //tit->Draw();
             std::cout<<"no problem 4"<<std::endl;
             ///////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_hodo"+nid,800,400);
             hodoF->SetTitle(Form("theta %3.1f", a));
             hodoF->SetStats(0);
             hodoF->SetTitle("hodo");
             hodoF->GetXaxis()->SetTitle("[mm]");
             hodoF->GetYaxis()->SetTitle("[mm]");
             hodoF->Draw("colz");
             ///////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_hodo_afterCut"+nid,800,400);
             hodo_afterCut->SetTitle(Form("theta %3.1f", a));
             hodo_afterCut->SetStats(0);
             hodo_afterCut->SetTitle("hodo");
             hodo_afterCut->GetXaxis()->SetTitle("[mm]");
             hodo_afterCut->GetYaxis()->SetTitle("[mm]");
             hodo_afterCut->Draw("colz");
             //////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_hodo_V_pos"+nid,800,500);
             hodo_multi_withmedVfiber->SetTitle(Form("theta %3.1f", a));
             hodo_multi_withmedVfiber->SetStats(0);
             hodo_multi_withmedVfiber->SetTitle("hodo V");
             hodo_multi_withmedVfiber->GetXaxis()->SetTitle("[mm]");
             hodo_multi_withmedVfiber->GetYaxis()->SetTitle("[mm]");
             hodo_multi_withmedVfiber->Draw("colz");
             //////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_hodo_H_pos"+nid,800,500);
             hodo_multi_withmedHfiber->SetTitle(Form("theta %3.1f", a));
             hodo_multi_withmedHfiber->SetStats(0);
             hodo_multi_withmedHfiber->SetTitle("hodo H");
             hodo_multi_withmedHfiber->GetXaxis()->SetTitle("[mm]");
             hodo_multi_withmedHfiber->GetYaxis()->SetTitle("[mm]");
             hodo_multi_withmedHfiber->Draw("colz");
             //////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_multi"+nid,800,400);
             countmulti_trg1->SetTitle(Form("theta %3.1f", a));
             countmulti_trg1->SetStats(0);
             countmulti_tof1->SetLineColor(kViolet);
             countmulti_tof2->SetLineColor(kGreen);
             countmulti_trg1->SetLineColor(kBlue);
             countmulti_trg2->SetLineColor(kCyan);
             countmulti_trgmzV->SetLineColor(kRed);
             countmulti_trgmzH->SetLineColor(kMagenta);
             countmulti_trg1->GetXaxis()->SetTitle("Multiplicity");
             countmulti_trg1->GetYaxis()->SetTitle("count [#]");
             countmulti_trg1->Draw();
             countmulti_tof2->Draw("same");
             countmulti_trg2->Draw("same");
             countmulti_tof1->Draw("same");
             countmulti_trgmzV->Draw("same");
             countmulti_trgmzH->Draw("same");
             TLegend *leg_multi = new TLegend(0.5,0.6,0.8,0.85);
             leg_multi->SetFillColor(0);
             leg_multi->SetFillStyle(0);
             leg_multi->SetBorderSize(0);
             leg_multi->SetFillStyle(0);
             leg_multi->AddEntry(countmulti_trg1,"Trigger 1","lp");
             leg_multi->AddEntry(countmulti_trg2,"Trigger 2","lp");
             leg_multi->AddEntry(countmulti_trgmzV,"Trigger MZ V","lp");
             leg_multi->AddEntry(countmulti_trgmzH,"Trigger MZ H","lp");
             leg_multi->AddEntry(countmulti_tof1,"TOF 1","lp");
             leg_multi->AddEntry(countmulti_tof2,"TOF 2","lp");
             leg_multi->Draw();
             ///////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_hodo_H_number"+nid,800,500);
             countmulti_hodoH->SetTitle(Form("theta %3.1f", a));
             countmulti_hodoH->SetStats(0);
             countmulti_hodoH->SetLineColor(kBlue);
             countmulti_hodoH->GetXaxis()->SetTitle("number of spikes[#]");
             countmulti_hodoH->GetYaxis()->SetTitle("count [#]");
             countmulti_hodoH->Draw();
             TLegend *leg_countmulti_hodoH = new TLegend(0.261905, 0.386076 , 0.907268,  0.995781 );
             leg_countmulti_hodoH->SetFillColor(0);
             leg_countmulti_hodoH->SetFillStyle(0);
             leg_countmulti_hodoH->SetBorderSize(0);
             leg_countmulti_hodoH->SetFillStyle(0);
             leg_countmulti_hodoH->AddEntry(countmulti_hodoH,"number of 'spikes' associated with selected H fiber","lp");
             leg_countmulti_hodoH->Draw();
             ///////////////////////////////////////////////////////////////////////
             prt_canvasAdd("r_hodo_V_number"+nid,800,500);
             countmulti_hodoV->SetTitle(Form("theta %3.1f", a));
             countmulti_hodoV->SetStats(0);
             countmulti_hodoV->SetLineColor(kViolet);
             countmulti_hodoV->GetXaxis()->SetTitle("number of spikes[#]");
             countmulti_hodoV->GetYaxis()->SetTitle("count [#]");
             countmulti_hodoV->Draw();
             TLegend *leg_countmulti_hodoV = new TLegend(0.261905, 0.386076 , 0.907268,  0.995781 );
             leg_countmulti_hodoV->SetFillColor(0);
             leg_countmulti_hodoV->SetFillStyle(0);
             leg_countmulti_hodoV->SetBorderSize(0);
             leg_countmulti_hodoV->SetFillStyle(0);
             leg_countmulti_hodoV->AddEntry(countmulti_hodoV,"number of 'spikes' associated with selected V fiber","lp");
             leg_countmulti_hodoV->Draw();
             */
            ///////////////////////////////////////////////////////////////////////
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
            //////////////////////////////////////////////////////////////////////
            //std::cout<<"no problem 4.1"<<std::endl;
            /*
             TCanvas * can1 = new TCanvas("can1","can1",0,0,800,400);
             fHist1-> Draw();
             fHist2-> Draw("Same");
             can1->Print("time.png");
             TCanvas * can2 = new TCanvas("can2","can2",0,0,800,400);
             fHist->Draw();
             line->Draw();
             line2->Draw();
             can2->Print("chere.png");
             */
            /*
             TCanvas * can3 = new TCanvas("can3","can3",0,0,800,400);
             hdelta_tof2tof1->Draw();
             hdelta_tof2tof1->SetLineColor(kRed+2);
             hdelta_tof2tof1_isproton->Draw("same");
             can3->Print("delta.png");
             */
            /*
             TCanvas * can3 = new TCanvas("can3","can3",0,0,800,400);
             hdelta_tof2tof1->Draw();
             hdelta_tof2tof1->SetLineColor(kRed+2);
             hdelta_tof2tof1_isproton->Draw("same");
             */
            /*
             prt_canvasAdd("r_alpha"+nid,800,400);
             falpha->SetTitle(Form("alpha - theta %3.1f", a));
             falpha->Draw();
             falphai->SetLineColor(kRed+2);
             if(falphai->GetEntries()>5)  falphai->Draw("same");
             prt_canvasAdd("r_photonEnergy"+nid,800,400);
             fHistPhotonEnergy->SetTitle(Form("PhotonEnergy - Theta %3.1f", a));
             fHistPhotonEnergy->Draw();
             */
            //std::cout<<"no problem 5"<<std::endl;

            /*
             TCanvas * can6 = new TCanvas("can6","can6",0,0,800,400);
             //hHodo->Draw();
             //hHodo->SetLineColor(kRed+2);
             //hHodo->SetMaximum(3100);
             //hHodo->SetMinimum(1000);
             hodoF->SetStats(0);
             hodoF->SetTitle("hodo");
             hodoF->GetXaxis()->SetTitle("[mm]");
             hodoF->GetYaxis()->SetTitle("[mm]");
             hodoF->Draw("colz");
             //can6->Print("Hodo.png");
             */

            /*
             TCanvas * can4 = new TCanvas("can4","can4",0,0,800,400);
             fHist0->Draw();
             can4->Print("diff.png");
             */
            /*
             TCanvas * can5 = new TCanvas("can5","can5",0,0,800,400);
             fHist->Draw();
             can5->Print("tangle.png");
             */
            /*
             TCanvas * can7 = new TCanvas("LE","LE",0,0,800,400);
             htrg1_le->SetLineColor(kGreen);
             htof1_le->SetLineColor(kOrange);
             htof2_le->SetLineColor(kViolet);
             htrg2_le->SetLineColor(kRed);
             htrgmzH_le->SetLineColor(kCyan);
             htrgmzV_le->SetLineColor(kMagenta);
             htrg1_le->Draw();
             htof1_le->Draw("same");
             htof2_le->Draw("same");
             htrg2_le->Draw("same");
             htrgmzH_le->Draw("same");
             htrgmzV_le->Draw("same");
             can7->Print("le.png");
             */

            /*
             TCanvas * can8 = new TCanvas("TOT","TOT",0,0,800,400);
             htrg1_tot->SetLineColor(kGreen);
             htof1_tot->SetLineColor(kOrange);
             htof2_tot->SetLineColor(kViolet);
             htrg2_tot->SetLineColor(kRed);
             htrgmzH_tot->SetLineColor(kCyan);
             htrgmzV_tot->SetLineColor(kMagenta);
             htrg1_tot->Draw();
             htof1_tot->Draw("same");
             htof2_tot->Draw("same");
             htrg2_tot->Draw("same");
             htrgmzH_tot->Draw("same");
             htrgmzV_tot->Draw("same");
             //can8->Print("tot.png");
             */
            /*
             htof1_tot->Fill(TOT_tof1);
             htof2_le->Fill(LE_tof2);
             htof2_tot->Fill(TOT_tof2);
             htrg1_le->Fill(LE_trg1);
             htrg1_tot->Fill(TOT_trg1);
             htrg2_le->Fill(LE_trg2);
             htrg2_tot->Fill(TOT_trg2);
             htrgmzV_le->Fill(LE_trgmzV);
             htrgmzV_tot->Fill(TOT_trgmzV);
             htrgmzH_le->Fill(LE_trgmzH);
             htrgmzH_tot->Fill(TOT_trgmzH);
             */
            std::cout<<"no problem 5.1"<<std::endl;
            //prt_waitPrimitive("r_cm"+nid);// wait here
            prt_canvasSave(2,0);
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
                //leg->SetFillColor(0);
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
                //c2->Print("example.pdf");
                //c2->Print(Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/prtdirc/build/spr/tcorr_%3.1f.png", a));
                //c2->Print(Form("/lustre/nyx/panda/aali/prtdrc_2017/final_2017/prtdirc/build/spr/tcorr_%3.1f.root", a));
                //c2->Print(Form("spr/tcorr_%3.1f.png", a));
                std::cout<<"no problem 5.2"<<std::endl;
                c2->Modified();
                c2->Update();
                //c2->WaitPrimitive("");
                std::cout<<"no problem 6"<<std::endl;
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
    //falpha->Reset();
    //falphai->Reset();
    fHist1->Reset();
    fHist2->Reset();
    fHist3->Reset();
    fHist4->Reset();

    fHist->Reset();
    fnHits->Reset();
    fnHits_p->Reset();
    fnHits_p_good->Reset();


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
            //prt_waitPrimitive("ff"); // wait
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
            //prt_waitPrimitive("ff"); // wait
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












