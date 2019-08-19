// -----------------------------------------
// PrtLutReco.h
//
// Created on: 13.07.2013
// Author: R.Dzhygadlo at gsi.de
// -----------------------------------------
// Class for reconstruction in DIRC using look-up table method

#ifndef PrtLutReco_h
#define PrtLutReco_h 1

#include "PrtEvent.h"
#include "PrtHit.h"

#include "TString.h"
#include "TFile.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TRandom.h"
#include "TGraph.h"
#include "TMultiGraph.h"
class PrtLutReco {

    public:

        // Standard constructors
        PrtLutReco(TString infile, TString lutfile, Int_t verbose=0);

        // Destructor
        ~PrtLutReco();
        void Run(Int_t start=0, Int_t end=0);

    private:
        Bool_t FindPeak(Double_t& cherenkovreco, Double_t& spr,Double_t a, Int_t tofpid=0);
        Int_t FindPdg(Double_t mom, Double_t cangle);
        void FitRing(Double_t& x0, Double_t& y0, Double_t& theta);
        Double_t fillLnDiffPPi(Double_t cangle, Int_t tofPid, Double_t mom);
        Double_t fillLnDiffPPi2(Double_t cangle, Int_t tofPid, Double_t mom);
        void ResetHists();
        Int_t fDetectorID;
        Double_t fBboxNum,fPipehAngle,fDphi,fBarPhi;
        TRandom fRand;
        Int_t fMethod;

        TClonesArray *fLut;
        TClonesArray *fTrackInfoArray;

        TFile *fFile;
        TTree *fTree;
        TChain *fChain;


        PrtEvent *fEvent, *phs_event;

    
        PrtHit fHit;
        PrtHit fHit2,fHit3, fHit4, fHit5, fHit6, fHit7;
 

        // Verbosity level
        Int_t fVerbose;
        Int_t nevents;
        TString fInputFile;
        TH1F *fHist, *fHist_copy, *fHist_correction, *fHist_same_path, *fHist_bg, *fHist_same_path_wotc;
        TH1F *fHistPi;
        TH1F *fHisti;
        TF1 *fFit;
        TSpectrum *fSpect;
        Double_t fAngleP;
        Double_t fAnglePi;
        Double_t fTest;
    
    TString phs;
    Int_t prt_entries_phs ;



TMultiGraph *mg;


    TGraph *gr_pi,*gr_p;
    
    
    //TChain *chain_phs[12][64];
    //Int_t entries_phs[12][64];
    //PrtHit hit_phs[12][64];
    TFile *fFile_phs;
    TTree *fTree_phs;  
    TClonesArray *fLut_phs;


    Int_t mcpid_phs, pixid_phs;
    
   

};

#endif

