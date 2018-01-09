#include "TFile.h"
#include "TTree.h"
#include "TList.h"
#include "TCanvas.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TArrayD.h"
#include "TMath.h"
#include "TDirectory.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include <iostream>
#include <unistd.h>
#include <complex>
#include <cmath>


static const int ncentbins = 11;
static const double centBins[]={0,5,10,15,20,25,30,35,40,50,60,70};
TH1D * centbins;
#include "/home/sanders/PbPb_2015/CMSSW_7_5_8_patch7/src/RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
using namespace hi;

//----------------------------------
// Tree Variables:
//
double centval;
int noff;
double vtx;
double epang[hi::NumEPNames];
Double_t qx[hi::NumEPNames];
Double_t qy[hi::NumEPNames];
Double_t q[hi::NumEPNames];
Double_t epmult[hi::NumEPNames];
unsigned int  runno_;
Double_t rescor[hi::NumEPNames];
Double_t rescorErr[hi::NumEPNames];
Double_t sumw[hi::NumEPNames];
TH2D * qxtrk2_;
TH2D * qytrk2_;
TH2D * qxtrk3_;
TH2D * qytrk3_;
TH2D * qxtrk4_;
TH2D * qytrk4_;
TH2D * qxtrk5_;
TH2D * qytrk5_;
TH2D * qxtrk6_;
TH2D * qytrk6_;
TH2D * qxtrk7_;
TH2D * qytrk7_;
TH2D * qcnt_;
TH2D * avpt_;

void ReadTree();
void QCheck() {
    
  centbins = new TH1D("centbins","centbins",ncentbins,centBins);
  ReadTree();
}

void ReadTree(){ 
  centbins->Reset();
  TFile * tfin;
  int NumEvnts = 0;
  int nbins = ncentbins;
  int filecnt=0;
  NumEvnts = 0;
  filecnt = 0;
  TString inFile="../../VNAnalysis/vnanal.root";
  tfin    = new TFile(inFile.Data(),"read");
  TH2D * qxy = new TH2D("qxy","qxy",50,-1,1,50,-1,1);
  TTree * tree = (TTree * ) tfin->Get("vnanalyzer/tree");
  tree->SetBranchAddress("Cent",       &centval);
  tree->SetBranchAddress("NtrkOff",    &noff);
  tree->SetBranchAddress("Vtx",        &vtx);
  tree->SetBranchAddress("epang",      &epang);
  tree->SetBranchAddress("sumw", &sumw);
  tree->SetBranchAddress("qx",         &qx);
  tree->SetBranchAddress("qy",         &qy);
  tree->SetBranchAddress("q",          &q);
  tree->SetBranchAddress("mult",       &epmult);
  tree->SetBranchAddress("Run",        &runno_);
  tree->SetBranchAddress("Rescor",     &rescor);
  tree->SetBranchAddress("RescorErr",  &rescorErr);
  tree->SetBranchAddress("qxtrk2",      &qxtrk2_);
  tree->SetBranchAddress("qytrk2",      &qytrk2_);
  tree->SetBranchAddress("qxtrk3",      &qxtrk3_);
  tree->SetBranchAddress("qytrk3",      &qytrk3_);
  tree->SetBranchAddress("qxtrk4",      &qxtrk4_);
  tree->SetBranchAddress("qytrk4",      &qytrk4_);
  tree->SetBranchAddress("qxtrk5",      &qxtrk5_);
  tree->SetBranchAddress("qytrk5",      &qytrk5_);
  tree->SetBranchAddress("qxtrk6",      &qxtrk6_);
  tree->SetBranchAddress("qytrk6",      &qytrk6_);
  tree->SetBranchAddress("qxtrk7",      &qxtrk7_);
  tree->SetBranchAddress("qytrk7",      &qytrk7_);
  tree->SetBranchAddress("qcnt",       &qcnt_);
  tree->SetBranchAddress("avpt",        &avpt_);
  cout<<"tree created"<<endl;
  for(int ievent = 0; ievent<tree->GetEntries(); ievent++) {
    tree->GetEntry(ievent);
    if(fabs(vtx)>15.) continue;
    int bin = -1; 
    bin =  centbins->FindBin(centval)-1;
    if(bin>=ncentbins) continue;
    if(bin<0) continue;
    centbins->Fill(centval);
    if(bin!=5) continue;
    for(int i = 1; i<qxtrk2_->GetNbinsX(); i++) {
      if(qxtrk2_->GetXaxis()->GetBinCenter(i)<0.3
 || qxtrk2_->GetXaxis()->GetBinCenter(i)>3.0)continue;
      for(int j = 5; j<=8; j++) {
	double x = qxtrk2_->GetBinContent(i,j);
	double y = qytrk2_->GetBinContent(i,j);
	double mag = qcnt_->GetBinContent(i,j);
	if(mag<=0) continue;
  
	qxy->Fill(x/mag,y/mag);
      }
    }
    ++NumEvnts; 
  }
  qxy->Draw("colz");
  cout<<"NumEvnts: "<<NumEvnts<<endl;
  
}
