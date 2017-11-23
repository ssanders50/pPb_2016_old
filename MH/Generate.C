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


TRandom3 * ran;

typedef complex<double> comp;

int epord_ = 2.;
bool trkoff = true;

static const double MaxCent = 70;
static const int MaxEvents = -1;
int MaxFiles = 50;
static const int ntrkbins = 15;
static const double trkBins[]={0,20,30,40,50,60,80,100,120,150,185,220,260,300,350,500};
static const int ncentbins = 11;
static const double centBins[]={0,5,10,15,20,25,30,35,40,50,60,70};
static const int nanals = 20;
enum AnalType {
  N2,       N3,     N4,      N5,       N6,     N7,      N42  , N523,  N523A, N63,  N62, N723, N723A, D24, D26, D34, D2232, D2232A, D2432, D2432A
};
string AnalNames[]={
  "N2",    "N3",   "N4",    "N5",    "N6",    "N7",     "N42"  ,"N523","N523A","N63","N62","N723","N723A", "D24","D26","D34","D2232","D2232A","D2432","D2432A"
};
int ANAL;
int epa;
int epb;
int epc;
int ep3a;
int ep3b;
int ep3c;

TH1D * trkbins;
TH1D * centbins;
TH2D * ptav[ncentbins];
TH2D * ptcnt[ncentbins];
TH2D * badcnt[ncentbins];
TH2D * qA[ncentbins][11];
TH2D * qB[ncentbins][11];
TH1D * qres;
TH1D * qBA[ncentbins][11];
TH1D * qCA[ncentbins][11];
TH1D * qCB[ncentbins][11];
TH1D * qBAcnt[ncentbins][11];
TH1D * qCAcnt[ncentbins][11];
TH1D * qCBcnt[ncentbins][11];
TH2D * qxav1[ncentbins];
TH2D * qyav1[ncentbins];
TH2D * qxav2[ncentbins];
TH2D * qyav2[ncentbins];
TH2D * qxav3[ncentbins];
TH2D * qyav3[ncentbins];
TH2D * qxav4[ncentbins];
TH2D * qyav4[ncentbins];
TH2D * qxav5[ncentbins];
TH2D * qyav5[ncentbins];
TH2D * qxav6[ncentbins];
TH2D * qyav6[ncentbins];
TH2D * qxav7[ncentbins];
TH2D * qyav7[ncentbins];
TH2D * qxycnt[ncentbins];
TH2D * wnA[ncentbins][11];
TH2D * wnB[ncentbins][11];
#include "HiEvtPlaneList.h"
using namespace hi;
Bool_t ispPb ;
struct qvec {
  TH2D * qA[ncentbins][11];
  TH2D * qB[ncentbins][11];
  TH2D * wnA[ncentbins][11];
  TH2D * wnB[ncentbins][11];
  TH1D * qBA[ncentbins][11];
  TH1D * qCA[ncentbins][11];
  TH1D * qCB[ncentbins][11];
  TH1D * qBAcnt[ncentbins][11];
  TH1D * qCAcnt[ncentbins][11];
  TH1D * qCBcnt[ncentbins][11];
} qanal[nanals];

#include "src/Fill.h"
int GetMidIndx(Int_t epord, TString midn) {
  int mid = -1;
  string match = midn.Data();
  for(int i = 0; i<hi::NumEPNames; i++) {
    string ct = hi::EPNames[i];
    if(epord==3) ct=hi::EPNames[i];
    if(match==ct) mid = i;
  }
  if(mid<0) {cout<<midn.Data()<< "not found"<<endl; return -1;}
  return mid;
}


Int_t flipENUM(int ein){
  if(ein<0) {
    cout<<"flipENUM called with ein = "<<ein<<"  EXPECT TO CRASH!"<<endl;
  }
  if(!ispPb) return ein;
  string ename = hi::EPNames[ein];
  if(ename.find("mid") != std::string::npos) return ein;
  if(ename.find("m") != std::string::npos) {
    ename.replace(ename.find("m"),1,"p");
  } else if (ename.find("p") != std::string::npos){
    ename.replace(ename.find("p"),1,"m");
  }
  //int epord = 2;
  //if(ename.find("3")!=std::string::npos) epord = 3;
  return GetMidIndx(epord_,ename);
}

#include "src/GetEventInfo.h"

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
TH2D * qxtrk1_;
TH2D * qytrk1_;
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
TH2D * qw_;
TH2D * avpt_;

Int_t NumEvents[40];
Int_t TotNumEvents; 
TString KeyNames[40];
int NumKeys;

string reac_;
//----------------------------------

void ReadTree(GetEventInfo * info,  string inlist);
int GetMidIndx(Int_t epord, TString midn);
void  GetNumEvents( string inlist);
void Generate(string anal="", TString reac="", string inlist="", int FileLimit=50) {
  MaxFiles = FileLimit;
  cout<<"ANAL: "<<anal<<endl;
  cout<<"REAC: "<<reac<<endl;
  FILE *  flist;
  flist = fopen(inlist.data(),"r");
  TString mid2n;
  TString mid3n;
#include "src/PbPbSetup.h"
  //#include "src/pPbSetup.h"
  ran = new TRandom3();
  reac_ = reac.Data();
  
  //
  //Locate information about data structure
  //
  GetEventInfo * info=0;
  char buf[120];
  while(fgets(buf,120,flist)!=NULL) {
    buf[strlen(buf)-1]=0;
    TString inFile = buf;
    FILE *ftest = fopen(inFile.Data(),"r");
    if(ftest==NULL) continue;
    fclose(ftest);
    TFile * tf    = new TFile(inFile.Data(),"read");
    if(tf->IsZombie())                 {cout<<"ZOMBIE:    " <<inFile.Data()<<endl; continue;}
    if(tf->TestBit(TFile::kRecovered)) {cout<<"RECOVERED: " <<inFile.Data()<<endl; continue;}
    tf->Close();
    cout<<"Found: "<<inFile.Data()<<endl;
    info = new GetEventInfo(inFile);
    if(info->status == 0) {cout<<inFile.Data()<<" not found or has error"<<endl; continue;}
    if(info->isPbPb2015pixel() || info->isMC()) trkoff = false;
    cout<<"setup with epord_: "<<epord_<<" ANAL: "<<AnalNames[ANAL].data()<<endl;
    int nbins = 0;
    if(trkoff) {
      nbins = ntrkbins;
    } else {
      nbins = ncentbins;
    }
    
    trkbins = new TH1D("trkbins","trkbins",ntrkbins,trkBins);
    trkbins->Sumw2();
    centbins = new TH1D("centbins","centbins",ncentbins,centBins);
    centbins->Sumw2();
    trkbins->SetDirectory(0);
    centbins->SetDirectory(0);

    qxy = (TH2D *) info->getTemplate()->Clone("qxy");
    qxxy = (TH2D *) info->getTemplate()->Clone("qxxy");
    qxyy = (TH2D *) info->getTemplate()->Clone("qxyy");
    qx2y3 = (TH2D *) info->getTemplate()->Clone("qx2y3");
    qy2x3 = (TH2D *) info->getTemplate()->Clone("qy2x3");
    qx2x2x3= (TH2D *) info->getTemplate()->Clone("qx2x2x3");
    qx2x3y2= (TH2D *) info->getTemplate()->Clone("qx2x3y2");
    qx3y2y2= (TH2D *) info->getTemplate()->Clone("qx3y2y2");
    qx2x2y3= (TH2D *) info->getTemplate()->Clone("qx2x2y3");
    qx2y2y3= (TH2D *) info->getTemplate()->Clone("qx2y2y3");
    qy2y2y3= (TH2D *) info->getTemplate()->Clone("qy2y2y3");    
    qcnt3  = (TH2D *) info->getTemplate()->Clone("qcnt3");
    for(int i = 0; i<ncentbins; i++) {
      ptav[i] = (TH2D *) info->getTemplate()->Clone(Form("ptav_%d",i));
      ptav[i]->SetDirectory(0);
      ptcnt[i] = (TH2D *) info->getTemplate()->Clone(Form("ptcnt_%d",i));
      ptcnt[i]->SetDirectory(0);
      badcnt[i] = (TH2D *) info->getTemplate()->Clone(Form("badcnt_%d",i));
      badcnt[i]->SetDirectory(0);
      qxav1[i] = (TH2D *) info->getTemplate()->Clone(Form("qxav1_%d",i));
      qyav1[i] = (TH2D *) info->getTemplate()->Clone(Form("qxav1_%d",i));
      qxav2[i] = (TH2D *) info->getTemplate()->Clone(Form("qxav2_%d",i));
      qyav2[i] = (TH2D *) info->getTemplate()->Clone(Form("qxav2_%d",i));
      qxav3[i] = (TH2D *) info->getTemplate()->Clone(Form("qxav3_%d",i));
      qyav3[i] = (TH2D *) info->getTemplate()->Clone(Form("qxav3_%d",i));
      qxav4[i] = (TH2D *) info->getTemplate()->Clone(Form("qxav4_%d",i));
      qyav4[i] = (TH2D *) info->getTemplate()->Clone(Form("qxav4_%d",i));
      qxav5[i] = (TH2D *) info->getTemplate()->Clone(Form("qxav5_%d",i));
      qyav5[i] = (TH2D *) info->getTemplate()->Clone(Form("qxav5_%d",i));
      qxav6[i] = (TH2D *) info->getTemplate()->Clone(Form("qxav6_%d",i));
      qyav6[i] = (TH2D *) info->getTemplate()->Clone(Form("qxav6_%d",i));
      qxav7[i] = (TH2D *) info->getTemplate()->Clone(Form("qxav7_%d",i));
      qyav7[i] = (TH2D *) info->getTemplate()->Clone(Form("qxav7_%d",i));
      qxycnt[i] = (TH2D *) info->getTemplate()->Clone(Form("qxcnt_%d",i));
    }

    for(int ian = 0; ian<nanals; ian++) {
      for(int i = 0; i<ncentbins; i++) {
	
	qanal[ian].qA[i][0] = (TH2D *) info->getTemplate()->Clone(Form("qA_%d",i));
	qanal[ian].qA[i][0]->SetDirectory(0);
	qanal[ian].qB[i][0] = (TH2D *) info->getTemplate()->Clone(Form("qB_%d",i));
	qanal[ian].qB[i][0]->SetDirectory(0);
	qanal[ian].wnA[i][0] = (TH2D *) info->getTemplate()->Clone(Form("wnA_%d",i));
	qanal[ian].wnA[i][0]->SetDirectory(0);
	qanal[ian].wnB[i][0] = (TH2D *) info->getTemplate()->Clone(Form("wnB_%d",i));
	qanal[ian].wnB[i][0]->SetDirectory(0);
	
	for(int j = 1; j<=10; j++) {
	  
	  qanal[ian].qA[i][j] = (TH2D *) info->getTemplate()->Clone(Form("qA_%d_%d",i,j)); 
	  qanal[ian].qA[i][j]->SetDirectory(0);
	  qanal[ian].qB[i][j] = (TH2D *) info->getTemplate()->Clone(Form("qB_%d_%d",i,j));
	  qanal[ian].qB[i][j]->SetDirectory(0);
	  qanal[ian].wnA[i][j] = (TH2D *) info->getTemplate()->Clone(Form("wnA_%d_%d",i,j));
	  qanal[ian].wnA[i][j]->SetDirectory(0);
	  qanal[ian].wnB[i][j] = (TH2D *) info->getTemplate()->Clone(Form("wnB_%d_%d",i,j));
	  qanal[ian].wnB[i][j]->SetDirectory(0);
	  
	}
	
	
	qanal[ian].qBA[i][0] = new TH1D(Form("qBA_%d",i),Form("qBA_%d",i),1,0,1);
	qanal[ian].qBA[i][0]->SetDirectory(0);
	qanal[ian].qCA[i][0] = new TH1D(Form("qCA_%d",i),Form("qCA_%d",i),1,0,1);
	qanal[ian].qCA[i][0]->SetDirectory(0);
	qanal[ian].qCB[i][0] = new TH1D(Form("qCB_%d",i),Form("qCB_%d",i),1,0,1);
	qanal[ian].qCB[i][0]->SetDirectory(0);
	qanal[ian].qBAcnt[i][0] = new TH1D(Form("qBAcnt_%d",i),Form("qBAcnt_%d",i),1,0,1);
	qanal[ian].qBAcnt[i][0]->SetDirectory(0);
	qanal[ian].qCAcnt[i][0] = new TH1D(Form("qCAcnt_%d",i),Form("qCAcnt_%d",i),1,0,1);
	qanal[ian].qCAcnt[i][0]->SetDirectory(0);
	qanal[ian].qCBcnt[i][0] = new TH1D(Form("qCBcnt_%d",i),Form("qCBcnt_%d",i),1,0,1);
	qanal[ian].qCBcnt[i][0]->SetDirectory(0);
	
	for(int j = 1; j<=10; j++) {
	  
	  
	  qanal[ian].qBA[i][j] = new TH1D(Form("qBA_%d_%d",i,j),Form("qBA_%d_%d",i,j),1,0,1);
	  qanal[ian].qCA[i][j] = new TH1D(Form("qCA_%d_%d",i,j),Form("qCA_%d_%d",i,j),1,0,1);      
	  qanal[ian].qCB[i][j] = new TH1D(Form("qCB_%d_%d",i,j),Form("qCB_%d_%d",i,j),1,0,1);
	  qanal[ian].qBA[i][j]->SetDirectory(0);
	  qanal[ian].qCA[i][j]->SetDirectory(0);
	  qanal[ian].qCB[i][j]->SetDirectory(0);
	  
	  qanal[ian].qBAcnt[i][j] = new TH1D(Form("qBAcnt_%d_%d",i,j),Form("qBAcnt_%d_%d",i,j),1,0,1);
	  qanal[ian].qCAcnt[i][j] = new TH1D(Form("qCAcnt_%d_%d",i,j),Form("qCAcnt_%d_%d",i,j),1,0,1);      
	  qanal[ian].qCBcnt[i][j] = new TH1D(Form("qCBcnt_%d_%d",i,j),Form("qCBcnt_%d_%d",i,j),1,0,1);
	  qanal[ian].qBAcnt[i][j]->SetDirectory(0);
	  qanal[ian].qCAcnt[i][j]->SetDirectory(0);
	  qanal[ian].qCBcnt[i][j]->SetDirectory(0);
	  
	}
	
      }
    }
    fclose(flist);
    break;
  }
  GetNumEvents(inlist);
  for(int i = 0; i<NumKeys; i++) cout<<i<<"\t"<<KeyNames[i].Data()<<"\t"<<NumEvents[i]<<endl;
  ReadTree(info, inlist);
}

void ReadTree(GetEventInfo * info, string inlist){ 
  TStopwatch * sw = new TStopwatch();
  trkbins->Reset();
  centbins->Reset();
  TFile * tfin;
  int NumEvnts = 0;
  int nbins = ncentbins;
  if(trkoff) nbins=ntrkbins;
  int filecnt=0;
  int NEvt = TotNumEvents*0.7;
  NumEvnts = 0;
  filecnt = 0;
  FILE * flist = fopen(inlist.data(),"r");
  char buf[120];
  sw->Start();
  TDirectory * foutdirs[ncentbins];
  //  string outFile = "/rfs/sanders/MH/"+AnalNames[ANAL];
  string outFile = "results/"+AnalNames[ANAL];
  string outFext = outFile+".root";
  //cout<<outFile<<endl;
  TFile * fout =new TFile(outFext.data(),"RECREATE");
  for(int i=0;i<nbins;i++) {
    foutdirs[i] = fout->mkdir(Form("%d_%d",(int)centBins[i],(int)centBins[i+1]));
  }
  for(int i=0;i<nbins;i++) {
    ptav[i]->Reset();
    ptcnt[i]->Reset();
    badcnt[i]->Reset();
  }
  while(fgets(buf,120,flist)!=NULL) {
      if(filecnt>=MaxFiles) continue;
      buf[strlen(buf)-1]=0;
      TString inFile=buf;
      ispPb = kFALSE;
      if(inFile.Contains("pPb")) ispPb = kTRUE;
      string pPbTag = "Pbp";
      if(ispPb) pPbTag = "pPb";
      if(reac_.find("PbPb")!=std::string::npos) pPbTag="PbPb";
      FILE *ftest = fopen(inFile.Data(),"r");
      if(ftest==NULL) continue;
      fclose(ftest);
      string inf = inFile.Data();
      tfin    = new TFile(inFile.Data(),"read");
      if(tfin->IsZombie())                 {cout<<"ZOMBIE:    " <<inFile.Data()<<endl; continue;}
      if(tfin->TestBit(TFile::kRecovered)) {cout<<"RECOVERED: " <<inFile.Data()<<endl; continue;}
      tfin->ResetErrno();
      ++filecnt;
      for(int i = 0; i<info->getNumKeys(); i++) {
	TTree * tree = (TTree * ) tfin->Get(Form("%s/tree",info->getKeyName(i).Data()));
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
	tree->SetBranchAddress("qxtrk1",      &qxtrk1_);
	tree->SetBranchAddress("qytrk1",      &qytrk1_);
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
	for(int ievent = 0; ievent<tree->GetEntries(); ievent++) {
	  if(MaxEvents>0&&NumEvnts>=MaxEvents) {
	    cout<<"MaxEvents: "<<MaxEvents<<endl;
	    cout<<"NumEvnts: "<<NumEvnts<<endl;
	    break;
	  }
	  tree->GetEntry(ievent);
	  if(fabs(vtx)>15.) continue;
	  if(centval>MaxCent) continue;
	  int bin = -1; 
	  if(trkoff) {
	    bin =  trkbins->FindBin(noff)-1;
	    if(bin>=ntrkbins) continue;
	    if(bin<0) continue;
	  } else {
	    bin =  centbins->FindBin(centval)-1;
	    if(bin>=ncentbins) continue;
	    if(bin<0) continue;
	  }
	  if((int)fmod( NumEvnts, NEvt/50) == 0 ) {
	    sw->Continue();
	    double elapse = sw->RealTime();
	    if(NumEvnts>100) {
	      string foutcopy = "cp "+ outFext +" "+ outFile+ "_"+to_string((int) (100*(NumEvnts/(double)NEvt)+0.5))+"_"+to_string(NumEvnts)+".root";
	      cout<<(int) (100*(NumEvnts/(double)NEvt)+0.5)<<"  Elapsed: "<<elapse<<"\t"<<"  Time per event: "<<elapse/(double)NumEvnts<< " : "<<foutcopy<<endl;
	      fout->cd();
	      trkbins->Write();
	      centbins->Write();
	      for(int k = 0; k<nbins; k++) {
	  	foutdirs[k]->cd();
	  	ptav[k]->Write("ptav");
	  	ptcnt[k]->Write("ptcnt");
	  	badcnt[k]->Write("badcnt");
		qxav1[k]->Write("qxav1");
		qyav1[k]->Write("qyav1");
		qxav2[k]->Write("qxav2");
		qyav2[k]->Write("qyav2");
		qxav3[k]->Write("qxav3");
		qyav3[k]->Write("qyav3");
		qxav4[k]->Write("qxav4");
		qyav4[k]->Write("qyav4");
		qxav5[k]->Write("qxav5");
		qyav5[k]->Write("qyav5");
		qxav6[k]->Write("qxav6");
		qyav6[k]->Write("qyav6");
		qxav7[k]->Write("qxav7");
		qyav7[k]->Write("qyav7");
		qxycnt[k]->Write("qxycnt");

		for(int ian = 0; ian<nanals; ian++) {
		  TDirectory * analdir = foutdirs[k]->mkdir(AnalNames[ian].data());
		  analdir->cd();
		  qanal[ian].qA[k][0]->Write("qA");
		  qanal[ian].qB[k][0]->Write("qB");
		  qanal[ian].qBA[k][0]->Write("qBA");
		  qanal[ian].qCA[k][0]->Write("qCA");
		  qanal[ian].qCB[k][0]->Write("qCB");
		  qanal[ian].qBAcnt[k][0]->Write("qBAcnt");
		  qanal[ian].qCAcnt[k][0]->Write("qCAcnt");
		  qanal[ian].qCBcnt[k][0]->Write("qCBcnt");
		  qanal[ian].wnA[k][0]->Write("wnA");
		  qanal[ian].wnB[k][0]->Write("wnB");
		  TDirectory * errdir = analdir->mkdir("SubEvents");
		  errdir->cd();
		  for(int j = 1; j<=10; j++) {
		    qanal[ian].qA[k][j]->Write(Form("qA_%d",j));
		    qanal[ian].qB[k][j]->Write(Form("qB_%d",j));
		    qanal[ian].qBA[k][j]->Write(Form("qBA_%d",j));
		    qanal[ian].qCA[k][j]->Write(Form("qCA_%d",j));
		    qanal[ian].qCB[k][j]->Write(Form("qCB_%d",j));
		    qanal[ian].qBAcnt[k][j]->Write(Form("qBAcnt_%d",j));
		    qanal[ian].qCAcnt[k][j]->Write(Form("qCAcnt_%d",j));
		    qanal[ian].qCBcnt[k][j]->Write(Form("qCBcnt_%d",j));
		    qanal[ian].wnA[k][j]->Write(Form("wnA_%d",j));
		    qanal[ian].wnB[k][j]->Write(Form("wnB_%d",j));
		  }
		}
	      }
	      //cout<<foutcopy<<endl;
	      fout->Close();
	      system(foutcopy.data());
	      fout =new TFile(outFext.data(),"RECREATE");
	      for(int i=0;i<nbins;i++) {
	  	foutdirs[i] = fout->mkdir(Form("%d_%d",(int)centBins[i],(int)centBins[i+1]));
	      }
	    }
	  }

	  if(epang[HFm2]<-5 || epang[HFp2]<-5 || epang[trackmid2]<-5) {
	    badcnt[bin]->Add(qcnt_);
	    continue;
	  }
	  trkbins->Fill(noff);
	  centbins->Fill(centval);
	  ptav[bin]->Add(avpt_);
	  ptcnt[bin]->Add(qcnt_);
	  qxav1[bin]->Add(qxtrk1_);
	  qyav1[bin]->Add(qytrk1_);
	  qxav2[bin]->Add(qxtrk2_);
	  qyav2[bin]->Add(qytrk2_);
	  qxav3[bin]->Add(qxtrk3_);
	  qyav3[bin]->Add(qytrk3_);
	  qxav4[bin]->Add(qxtrk4_);
	  qyav4[bin]->Add(qytrk4_);
	  qxav5[bin]->Add(qxtrk5_);
	  qyav5[bin]->Add(qytrk5_);
	  qxav6[bin]->Add(qxtrk6_);
	  qyav6[bin]->Add(qytrk6_);
	  qxav7[bin]->Add(qxtrk7_);
	  qyav7[bin]->Add(qytrk7_);
	  qxycnt[bin]->Add(qcnt_);

	  int evtchar = centval;
	  if(trkoff) evtchar = noff;
	  for(int ian = 0; ian<nanals; ian++) {
	    if(ian==N2) Fill_N( N2, bin, qxtrk2_, qytrk2_, qcnt_, qx[HFp2], qy[HFp2], qx[HFm2], qy[HFm2], qx[trackmid2], qy[trackmid2], sumw[HFp2], sumw[HFm2], sumw[trackmid2]);
	    if(ian==N3) Fill_N( N3,  bin, qxtrk3_, qytrk3_, qcnt_, qx[HFp3], qy[HFp3], qx[HFm3], qy[HFm3], qx[trackmid3], qy[trackmid3], sumw[HFp3], sumw[HFm3], sumw[trackmid3]);
	    if(ian==N4) Fill_N( N4,  bin, qxtrk4_, qytrk4_, qcnt_, qx[HFp4], qy[HFp4], qx[HFm4], qy[HFm4], qx[trackmid4], qy[trackmid4], sumw[HFp4], sumw[HFm4], sumw[trackmid4]);
	    if(ian==N5) Fill_N( N5,  bin, qxtrk5_, qytrk5_, qcnt_, qx[HFp5], qy[HFp5], qx[HFm5], qy[HFm5], qx[trackmid5], qy[trackmid5], sumw[HFp5], sumw[HFm5], sumw[trackmid5]);
	    if(ian==N6) Fill_N( N6,  bin, qxtrk6_, qytrk6_, qcnt_, qx[HFp6], qy[HFp6], qx[HFm6], qy[HFm6], qx[trackmid6], qy[trackmid6], sumw[HFp6], sumw[HFm6], sumw[trackmid6]);
	    if(ian==N7) Fill_N( N7,  bin, qxtrk7_, qytrk7_, qcnt_, qx[HFp7], qy[HFp7], qx[HFm7], qy[HFm7], qx[trackmid7], qy[trackmid7], sumw[HFp7], sumw[HFm7], sumw[trackmid7]);
	    if(ian==N523)  Fill_N523(N523, bin, qxtrk5_, qytrk5_, qcnt_, qx, qy, sumw);
	    if(ian==N523A) Fill_N523A(N523A,bin, qxtrk5_, qytrk5_, qcnt_, qx, qy, sumw);
	    if(ian==N63)   Fill_N63(N63, bin, qxtrk6_, qytrk6_, qcnt_, qx, qy, sumw);
	    if(ian==N62)   Fill_N62(N62, bin, qxtrk6_, qytrk6_, qcnt_, qx, qy, sumw);
	    if(ian==N723)  Fill_N723(N723,bin, qxtrk7_, qytrk7_, qcnt_, qx, qy, sumw);
	    if(ian==N723A) Fill_N723A(N723A,bin, qxtrk7_, qytrk7_, qcnt_, qx, qy, sumw);
	    if(ian==D24)   Fill_D24(D24,bin, qxtrk2_, qytrk2_, qcnt_, qx, qy, sumw);
	    if(ian==D26)     Fill_D26(D26,bin, qxtrk2_, qytrk2_, qcnt_, qx, qy, sumw);
	    if(ian==D34)   Fill_D34(D34,bin, qxtrk3_, qytrk3_, qcnt_, qx, qy, sumw);
	    if(ian==D2232) Fill_D2232(D2232,bin,  qxtrk2_, qytrk2_, qxtrk3_, qytrk3_, qcnt_, qx, qy, sumw);
	    if(ian==D2232A)Fill_D2232A(D2232A,bin, qxtrk2_, qytrk2_, qxtrk3_, qytrk3_, qcnt_, qx, qy, sumw);
	    if(ian==D2432) Fill_D2432(D2432,bin,  qxtrk2_, qytrk2_, qxtrk3_, qytrk3_, qcnt_, qx, qy, sumw);
	    if(ian==D2432A)Fill_D2432A(D2432A,bin, qxtrk2_, qytrk2_, qxtrk3_, qytrk3_, qcnt_, qx, qy, sumw);
	  }
	  ++NumEvnts;
	}
	tfin->Close();
      }
  }
  fout->cd();
  trkbins->Write();
  centbins->Write();
  for(int k = 0; k<nbins; k++) {
    foutdirs[k]->cd();
    ptav[k]->Write();
    ptcnt[k]->Write();
    badcnt[k]->Write();
    qxav1[k]->Write();
    qyav1[k]->Write();
    qxav2[k]->Write();
    qyav2[k]->Write();
    qxav3[k]->Write();
    qyav3[k]->Write();
    qxav4[k]->Write();
    qyav4[k]->Write();
    qxav5[k]->Write();
    qyav5[k]->Write();
    qxav6[k]->Write();
    qyav6[k]->Write();
    qxav7[k]->Write();
    qyav7[k]->Write();
    qxycnt[k]->Write();
    for(int ian = 0; ian<nanals; ian++) {
      TDirectory * analdir = foutdirs[k]->mkdir(AnalNames[ian].data());
      analdir->cd();
      qanal[ian].qA[k][0]->Write("qA");
      qanal[ian].qB[k][0]->Write("qB");
      qanal[ian].qBA[k][0]->Write("qBA");
      qanal[ian].qCA[k][0]->Write("qCA");
      qanal[ian].qCB[k][0]->Write("qCB");
      qanal[ian].qBAcnt[k][0]->Write("qBAcnt");
      qanal[ian].qCAcnt[k][0]->Write("qCAcnt");
      qanal[ian].qCBcnt[k][0]->Write("qCBcnt");
      qanal[ian].wnA[k][0]->Write("wnA");
      qanal[ian].wnB[k][0]->Write("wnB");
      TDirectory * errdir = analdir->mkdir("SubEvents");
      errdir->cd();
      for(int j = 1; j<=10; j++) {
	qanal[ian].qA[k][j]->Write(Form("qA_%d",j));
	qanal[ian].qB[k][j]->Write(Form("qB_%d",j));
	qanal[ian].qBA[k][j]->Write(Form("qBA_%d",j));
	qanal[ian].qCA[k][j]->Write(Form("qCA_%d",j));
	qanal[ian].qCB[k][j]->Write(Form("qCB_%d",j));
	qanal[ian].qBAcnt[k][j]->Write(Form("qBAcnt_%d",j));
	qanal[ian].qCAcnt[k][j]->Write(Form("qCAcnt_%d",j));
	qanal[ian].qCBcnt[k][j]->Write(Form("qCBcnt_%d",j));
	qanal[ian].wnA[k][j]->Write(Form("wnA_%d",j));
	qanal[ian].wnB[k][j]->Write(Form("wnB_%d",j));
      }
    }
  }
  fout->Close();
  
}
void GetNumEvents(string inlist){ 
  //Loop over datasets
  cout<<"GetNumEvents"<<endl;

  int ntrig = 1;
  TotNumEvents = 0;
  for(int i = 0; i<40; i++) NumEvents[i] = 0;
  TString tr = "";
  TFile * tfin;
  int filecnt = 0;
  char buf[120];
  FILE * flist = fopen(inlist.data(),"r");
  while(fgets(buf,120,flist)!=NULL) {
    if(filecnt>=MaxFiles) continue;
    buf[strlen(buf)-1]=0;
    TString inFile=buf;
    FILE *ftest = fopen(inFile.Data(),"r");
    if(ftest==NULL) continue;
    fclose(ftest);
    TList * flist;
    int nkeys = 0;
    tfin    = new TFile(inFile.Data(),"read");
    if(tfin->IsZombie())                 {cout<<"ZOMBIE:    " <<inFile.Data()<<endl; continue;}
    if(tfin->TestBit(TFile::kRecovered)) {cout<<"RECOVERED: " <<inFile.Data()<<endl; continue;}
    tfin->ResetErrno();
    flist = tfin->GetListOfKeys();
    ++filecnt;
    while(flist->At(nkeys) != flist->Last()) ++nkeys;
    ++nkeys;
    for(int i = 0; i<nkeys; i++) {
      TTree * tree = (TTree *) tfin->Get(Form("%s/tree",flist->At(i)->GetName()));
      if(NumKeys==0) KeyNames[i] = flist->At(i)->GetName();
      ++NumEvents[i] += tree->GetEntries();
      ++TotNumEvents += tree->GetEntries();
    }
    if(NumKeys==0) NumKeys = nkeys;
    tfin->Close();
    for(int i = 0; i<NumKeys; i++) {
      cout<<filecnt<<"\t"<<inFile.Data()<<"\t"<<tr.Data()<<"\t"<<KeyNames[i].Data()<<"\t NumEvents: "<<NumEvents[i]<<endl;
    }
  }
  fclose(flist);
  cout<<"Total number of events: "<<TotNumEvents<<endl;
  return;
}
