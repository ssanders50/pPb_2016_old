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
#include <iostream>
#include <unistd.h>
#include <complex>
#include <cmath>

#include "HiEvtPlaneList.h"
using namespace hi;

string rpnames[hi::NumEPNames];
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
TH2D * qxtrk_;
TH2D * qytrk_;
TH2D * qxtrk3_;
TH2D * qytrk3_;
TH2D * qcnt_;
TH2D * avpt_;

Int_t NumEvents[40];
Int_t TotNumEvents; 
TString KeyNames[40];
int NumKeys;

string reac_;
//----------------------------------

void MergeFiles(string inlist="inlist.dat", string tag="h_", int num=10) {
  int evcnt=0;
  FILE *  flist;
  flist = fopen(inlist.data(),"r");
  char buf[200];
  FILE * fout = fopen("merge.sh","w");
  int ncount = 0;
  int m = 0;
  string line;
  while(fgets(buf,200,flist)!=NULL) {
    buf[strlen(buf)-1]=0;
    string inFile=buf;
    FILE *ftest = fopen(inFile.data(),"r");
    if(ftest==NULL) continue;
    fclose(ftest);
    TFile * tfin    = new TFile(inFile.data(),"read");
    if(tfin->IsZombie())                 {
      cout<<"ZOMBIE:    " <<inFile.data()<<endl; 
      string remove = "rm "+inFile;
      system(remove.data());
      continue;
    }
    if(tfin->TestBit(TFile::kRecovered)) {cout<<"RECOVERED: " <<inFile.data()<<endl; continue;}
    tfin->ResetErrno();
    tfin->Close();
    if(fmod(m,num)==0) {
      if(m!=0) {
	cout<<line<<endl;
	system(line.data());
      }
      line = "hadd -f tmp/"+ tag + to_string(ncount++)+".root "+inFile;
    } else {
      line += " "+inFile;
    }
    m++;

  }
  if(m>0) system(line.data());
}
