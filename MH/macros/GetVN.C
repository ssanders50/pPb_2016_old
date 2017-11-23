#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include <iostream>
static const int ncentbins = 11;
static const int centBins[]={0,5,10,15,20,25,30,35,40,50,60,70};
static const double centRefBins[]={0,5,10,15,20,25,30,35,40,50,60,70};
static const int cbins = 14;
static const int cmin[]={0, 5,10,15,20,25,30,35,40,50,60,  0,20, 60};
static const int cmax[]={5,10,15,20,25,30,35,40,50,60,70, 20,60,100};
double EtaMin = -0.8;
double EtaMax = 0.8;
string FigDir = "";
string FigSubDir = "";
enum AnalType {
    N1MCm22SUB3,       N1MCm18SUB3,        N1MCm14SUB3,      N1MCm10SUB3, 
    N1MCm06SUB3,       N1MCm02SUB3,        N1MCp22SUB3,      N1MCp18SUB3, 
    N1MCp14SUB3,       N1MCp10SUB3,        N1MCp06SUB3,       N1MCp02SUB3,
         N1SUB2,            N1SUB3,                 
      N112ASUB2,         N112ASUB3,          N123ASUB2,        N123ASUB3,            
         N2SUB2,            N2SUB3,             N3SUB2,           N3SUB3,     
         N4SUB2,            N4SUB3,            N42SUB2,          N42SUB3,          
       N42ASUB2,          N42ASUB3,             N5SUB2,           N5SUB3,             
         N6SUB2,            N6SUB3,             N7SUB2,           N7SUB3,          
       N523SUB2,          N523SUB3,          N523ASUB2,        N523ASUB3,
       N723SUB2,          N723SUB3,          N723ASUB2,        N723ASUB3,   
        N62SUB2,           N62SUB3,           N62ASUB2,         N62ASUB3,
        N63SUB2,       	   N63SUB3,           N63ASUB2,         N63ASUB3,         
        D24SUB2,           D24SUB3,           D24ASUB2,         D24ASUB3,          
	D34SUB2,           D34SUB3,           D34ASUB2,         D34ASUB3,          
      D2232SUB2,        D2232SUB3,        
      D2432SUB2,         D2432SUB3,         D2232ASUB2,       D2232ASUB3,        
     D2432ASUB2,        D2432ASUB3,            D26SUB2,          D26SUB3,  
       D26ASUB2,          D26ASUB3,            
	   CHI4,             CHI4A,               CHI5,            CHI5A,  
	  CHI62,            CHI62A,              CHI63,           CHI63A, 
	   CHI7,             CHI7A,
          N2EFF,           N2NOEFF,            N723EFF,        N723NOEFF,         
       D2432EFF,        D2432NOEFF,            CHI7EFF,        CHI7NOEFF,            
         N42EFF,              LAST
};
string AnalNames[]={
  "N1MCm22SUB3",     "N1MCm18SUB3",      "N1MCm14SUB3",    "N1MCm10SUB3",
  "N1MCm06SUB3",     "N1MCm02SUB3",      "N1MCp22SUB3",    "N1MCp18SUB3", 
  "N1MCp14SUB3",     "N1MCp10SUB3",      "N1MCp06SUB3",    "N1MCp02SUB3",
       "N1SUB2",          "N1SUB3",              
    "N112ASUB2",       "N112ASUB3",        "N123ASUB2",      "N123ASUB3",          
       "N2SUB2",          "N2SUB3",           "N3SUB2",         "N3SUB3", 
       "N4SUB2",          "N4SUB3",          "N42SUB2",        "N42SUB3",        
     "N42ASUB2",        "N42ASUB3",           "N5SUB2",         "N5SUB3",          
       "N6SUB2",          "N6SUB3",           "N7SUB2",         "N7SUB3",      
     "N523SUB2",        "N523SUB3",        "N523ASUB2",      "N523ASUB3",      
     "N723SUB2",        "N723SUB3",        "N723ASUB2",      "N723ASUB3",       
      "N62SUB2",         "N62SUB3",         "N62ASUB2",       "N62ASUB3",      
      "N63SUB2",         "N63SUB3",         "N63ASUB2",       "N63ASUB3",        
      "D24SUB2",         "D24SUB3",         "D24ASUB2",       "D24ASUB3",        
      "D34SUB2",         "D34SUB3",         "D34ASUB2",       "D34ASUB3",        
     "D2232SUB2",      "D2232SUB3",      
    "D2432SUB2",       "D2432SUB3",       "D2232ASUB2",     "D2232ASUB3",    
   "D2432ASUB2",      "D2432ASUB3",          "D26SUB2",        "D26SUB3",  
     "D26ASUB2",        "D26ASUB3",     
         "CHI4",           "CHI4A",             "CHI5",          "CHI5A",        
        "CHI62",          "CHI62A",            "CHI63",         "CHI63A",           
         "CHI7",           "CHI7A", 
        "N2EFF",         "N2NOEFF",          "N723EFF",      "N723NOEFF",        
     "D2432EFF",      "D2432NOEFF",          "CHI7EFF",      "CHI7NOEFF",         
       "N42EFF",             "LAST"
};
string ytitle[]={
"v_{1}\{#Psi_{1,MC} (-2.4<#eta<-2.0)\}", "v_{1}\{#Psi_{1,MC} (-2.0<#eta<-1.6)\}", "v_{1}\{#Psi_{1,MC} (-1.6<#eta<-1.2)\}", "v_{1}\{#Psi_{1,MC} (-1.2<#eta<-0.8)\}",
"v_{1}\{#Psi_{1,MC} (-0.8<#eta<-0.4)\}",  "v_{1}\{#Psi_{1,MC} (-0.4<#eta<0.0)\}",   "v_{1}\{#Psi_{1,MC} (2.0<#eta<2.4)\}",   "v_{1}\{#Psi_{1,MC} (1.6<#eta<2.0)\}",
  "v_{1}\{#Psi_{1,MC} (1.2<#eta<1.6)\}",   "v_{1}\{#Psi_{1,MC} (0.8<#eta<1.2)\}",   "v_{1}\{#Psi_{1,MC} (0.4<#eta<0.8)\}",   "v_{1}\{#Psi_{1,MC} (0.0<#eta<0.4)\}",
            "v_{1}\{#Psi_{1}\}",             "v_{1}\{#Psi_{1}\}",             
     "v_{1}\{#Psi{1},#Psi{2}\}",     "v_{1}\{#Psi{1},#Psi{2}\}",        "v_{1}\{#Psi{2},#Psi{3}\}",         "v_{1}\{#Psi{2},#Psi{3}\}",                        
                        "v_{2}",                        "v_{2}",                           "v_{3}",                            "v_{3}",   
                        "v_{4}",                        "v_{4}",               "v_{4}\{#Psi_{2}\}",                "v_{4}\{#Psi_{2}\}", 
 "v_{4}\{#Psi_{2A},#Psi_{2B}\}", "v_{4}\{#Psi_{2A},#Psi_{2B}\}",                           "v_{5}",                            "v_{5}",                       
                        "v_{6}",                        "v_{6}",                           "v_{7}",                            "v_{7}",       
    "v_{5}\{#Psi_{2},#Psi_{3}\}",  "v_{5}\{#Psi_{2},#Psi_{3}\}",    "v_{5}\{#Psi_{2A},#Psi_{3B}\}",     "v_{5}\{#Psi_{2A},#Psi_{3B}\}",    
    "v_{7}\{#Psi_{2},#Psi_{3}\}",    "v_{7}\{#Psi_{2},#Psi_{3}\}",  "v_{7}\{#Psi_{2A},#Psi_{3B}\}",     "v_{7}\{#Psi_{2A},#Psi_{3B}\}",          
            "v_{6}\{#Psi_{2}\}",            "v_{6}\{#Psi_{2}\}", "v_{6}\{#Psi_{2A}^{2},#Psi_{2B}\}",  "v_{6}\{#Psi_{2A}^{2}#Psi_{2B}\}",          
            "v_{6}\{#Psi_{3}\}",            "v_{6}\{#Psi_{3}\}",    "v_{6}\{#Psi_{3A},#Psi_{3B}\}",      "v_{6}\{#Psi_{3A},#Psi_{3B}\}",            
	              "D24SUB2",                      "D24SUB3",                        "D24ASUB2",                          "D24ASUB3",      
                      "D34SUB2",                      "D34SUB3",                        "D34ASUB2",                          "D34ASUB3",            
                    "D2232SUB2",                         "D2232SUB3",            
                    "D2432SUB2",                    "D2432SUB3",                      "D2232ASUB2",                        "D2232ASUB3",  
                   "D2432ASUB2",                   "D2432ASUB3",                         "D26SUB2",                           "D26SUB3",           
                     "D26ASUB2",                     "D26ASUB3",     
                     "#chi_{4}",                    "#chi_{4A}",                        "#chi_{5}",                         "#chi_{5A}",            
                    "#chi_{62}",                   "#chi_{62A}",                       "#chi_{63}",                        "#chi_{63A}",                     
                     "#chi_{7}",                    "#chi_{7A}",           
                   "v_{2}(eff)",                 "v_{2}(noeff)",                    "v_{723}(eff)",                    "v_{723}(noeff)",               
                     "D2432EFF",                   "D2432NOEFF",                   "#chi_{7}(eff)",                   "#chi_{7}(moeff)",       
                  "v_{42}(eff)",                     "LAST"
};
FILE * outint;
string soutint;
string sspec;
FILE * outspec;
TFile * fin;
double fw = 1.3;
int ANAL = -1;
double fakescale = 1.;
bool isNominal = true;
bool isWide = false;
bool isNarrow = false;
bool isLoose = false;
bool isTight = false;
bool rew = false;
string stag;
TH1D * centbins;
TH1I * centRef;

double FakeAndEff(int cent, double pt, double &eff) {
  double cb[14]={0,5,10,15,20,25,30,35,40,50,60,70,80,100};
  double cbe[6]={0,5,10,30,50,100};
  TFile * f=0;
  TFile * e=0;
  if(isNominal) {
    f = new TFile("EffAndFake/FakeRatesPixelPbPb_nominal.root");
    e = new TFile("EffAndFake/EffCorrectionsPixelPbPb_nominal.root");
  }
  if(isWide) { 
    f = new TFile("EffAndFake/FakeRatesPixelPbPb_wide.root");
    e = new TFile("EffAndFake/EffCorrectionsPixelPbPb_wide.root");
  }
  if(isNarrow) {
    f = new TFile("EffAndFake/FakeRatesPixelPbPb_narrow.root");
    e = new TFile("EffAndFake/EffCorrectionsPixelPbPb_narrow.root");
  }
  if(isLoose){ 
    f = new TFile("EffAndFake/FakeRatesPixelPbPb_loose.root");
    e = new TFile("EffAndFake/EffCorrectionsPixelPbPb_loose.root");
  }
  if(isTight){ 
    f = new TFile("EffAndFake/FakeRatesPixelPbPb_tight.root");
    e = new TFile("EffAndFake/EffCorrectionsPixelPbPb_tight.root");
  }
  TH1I * cen = new TH1I("cen","cen",13,cb);
  TH1I * cene = new TH1I("cene","cene",5,cbe);
  int ib = cen->FindBin(cent)-1;
  int ibe = cene->FindBin(cent)-1;
  string rc = "hfak_"+to_string((int)cb[ib])+"_"+to_string((int)cb[ib+1]);
  TH2D * hf = (TH2D *) f->Get(rc.data());
  string re = "Eff_"+to_string((int)cbe[ibe])+"_"+to_string((int)cbe[ibe+1]);
  TH2D * he = (TH2D *) e->Get(re.data());
  int ptbin = hf->GetYaxis()->FindBin(pt);
  int etabinmin = hf->GetXaxis()->FindBin(-0.8);
  int etabinmax = hf->GetXaxis()->FindBin(0.79);
  double val = 0;
  eff = 0;
  for(int i = etabinmin; i<=etabinmax; i++) {
    val += hf->GetBinContent(i,ptbin);
    eff += he->GetBinContent(i,ptbin);
  }
  val /=(double)(etabinmax-etabinmin+1);
  eff /=(double)(etabinmax-etabinmin+1);
  cen->Delete();
  cene->Delete();
  f->Close();
  e->Close();
  return val ;
}



TH2D * ptcntEff(TH2D * ptcnt, double cent) {
  double cbe[6]={0,5,10,30,50,100};
  TFile * e=0;
  if(isNominal) {
    e = new TFile("EffAndFake/EffCorrectionsPixelPbPb_nominal.root");
  }
  if(isWide) { 
    e = new TFile("EffAndFake/EffCorrectionsPixelPbPb_wide.root");
  }
  if(isNarrow) {
    e = new TFile("EffAndFake/EffCorrectionsPixelPbPb_narrow.root");
  }
  if(isLoose){ 
    e = new TFile("EffAndFake/EffCorrectionsPixelPbPb_loose.root");
  }
  if(isTight){ 
    cout<<"use tight cut efficiency"<<endl;
    e = new TFile("EffAndFake/EffCorrectionsPixelPbPb_tight.root");
  }
  TH1I * cene = new TH1I("cene","cene",5,cbe);
  int ibe = cene->FindBin(cent)-1;
  string re = "Eff_"+to_string((int)cbe[ibe])+"_"+to_string((int)cbe[ibe+1]);
  TH2D * he = (TH2D *) e->Get(re.data());
  TH2D * hsEff = (TH2D *) ptcnt->Clone("ptcntEff");
  hsEff->Reset();
  hsEff->SetDirectory(0);
  bool skipeff = false;
  if(AnalNames[ANAL]==AnalNames[N2EFF] || AnalNames[ANAL]==AnalNames[N723EFF] || AnalNames[ANAL]==AnalNames[D2432EFF]
     || AnalNames[ANAL]==AnalNames[CHI7EFF] || AnalNames[ANAL]==AnalNames[N42EFF]) skipeff=true;
  for(int i = 1; i<=ptcnt->GetNbinsX(); i++) {
    for(int j = 1; j<=ptcnt->GetNbinsY(); j++) {
      double pt = ptcnt->GetXaxis()->GetBinCenter(i);
      double eta = ptcnt->GetYaxis()->GetBinCenter(j);
      int ptbin = he->GetYaxis()->FindBin(pt);
      int etabin = he->GetXaxis()->FindBin(eta);
      double eff =  he->GetBinContent(etabin,ptbin);
     
      if(eff<=0 || skipeff) eff = 1;
      hsEff->SetBinContent(i,j,ptcnt->GetBinContent(i,j)/eff);
    }
  }
  cene->Delete();
  e->Close();
  return hsEff ;
}
string rootFile;
TGraphErrors * GetVNPt(int replay, int bin, double etamin, double etamax, TGraphErrors * &gA, TGraphErrors * &gB, TGraphErrors * &gspec, double *  resA, double * resB, double & vint, double & vinte, bool nonorm=false){
  cout<<"========================== "<<AnalNames[replay]<<"  with bin,etamin,etamax: "<<bin<<"\t"<<etamin<<"\t"<<etamax<<endl;
  TH1D * hspec = 0;
  TH1D * xpt=0;
  TH1D * sp=0;
  TH1D * vn = 0;
  TH1D * vnA = 0;
  TH1D * vnB = 0;
  TH1D * vn2 = 0;
  TH1D * vnA2 = 0;
  TH1D * vnB2 = 0;
  TH1D * qA1=0;
  TH1D * qB1=0;
  TH1D * wA1=0;
  TH1D * wB1 = 0;
  TH2D * qAe[10];
  TH2D * qBe[10];
  TH2D * wnAe[10];
  TH2D * wnBe[10];
  TH1D * qAe1[10]={0};
  TH1D * qBe1[10]={0};
  TH1D * wAe1[10]={0};
  TH1D * wBe1[10]={0};
  TH1D * vnAe=0;
  TH1D * vnBe=0;
  TH1D * vne=0;
  TH2D * ptav=0;
  TH2D * ptcnt=0;
  TH2D * badcnt=0;
  TH2D * qA=0;
  TH2D * qB=0;
  TH2D * wnA=0;
  TH2D * wnB=0;
  int jmin = centRef->FindBin(cmin[bin])-1;
  int jmax = centRef->FindBin(cmax[bin]-0.01)-1;
  cout<<"jmin,jmax: "<<jmin<<"\t"<<jmax<<endl;
  ANAL = replay;
  string strip = AnalNames[replay];
  bool sub2 = false;
  TString subtest = AnalNames[replay];
  if(subtest.Contains("SUB2")) sub2 = true;
  if(sub2) {
    cout<<"Use 2 subevent weighting"<<endl;
  strip = strip.substr(0,strip.find("SUB2"));
  } else {
    cout<<"Use 3 subevent weighting"<<endl;
  strip = strip.substr(0,strip.find("SUB3"));
  }
  double qBA = 0;
  double qCA = 0;
  double qCB = 0;
  double qBAcnt = 0;
  double qCAcnt = 0;
  double qCBcnt = 0;
  double qBAe[10]={0};
  double qCAe[10]={0};
  double qCBe[10]={0};
  double qBAecnt[10]={0};
  double qCAecnt[10]={0};
  double qCBecnt[10]={0};
  double centcnt = 0;
  for(int j = jmin; j<=jmax; j++) {
    string crange = to_string(cmin[j])+"_"+to_string(cmax[j]);
    cout<<"crange: "<<crange<<endl;
    if(j==jmin) {
      ptav = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/ptav",crange.data()));
      ptcnt = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/ptcnt",crange.data()));
      badcnt = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/badcnt",crange.data()));
      centbins = (TH1D * ) fin->Get("vnanalyzer/centres");
      centcnt+=centbins->GetBinContent(j);
      qA = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qA",crange.data(),strip.data()));
      qB = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qB",crange.data(),strip.data()));
      wnA = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/wnA",crange.data(),strip.data()));
      wnB = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/wnB",crange.data(),strip.data()));
      qBA += ((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qBA",crange.data(),strip.data())))->GetBinContent(1);
      qBAcnt+=((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qBAcnt",crange.data(),strip.data())))->GetBinContent(1);
      qCA += ((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qCA",crange.data(),strip.data())))->GetBinContent(1);
      qCAcnt+=((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qCAcnt",crange.data(),strip.data())))->GetBinContent(1); 
      qCB += ((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qCB",crange.data(),strip.data())))->GetBinContent(1);
      qCBcnt+=((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qCBcnt",crange.data(),strip.data())))->GetBinContent(1);
      
      for(int i = 0; i<10; i++) {
	qAe[i] = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qA_%d",crange.data(),strip.data(),i+1));
	qBe[i] = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qB_%d",crange.data(),strip.data(),i+1));
	wnAe[i] = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/wnA_%d",crange.data(),strip.data(),i+1));
	wnBe[i] = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/wnB_%d",crange.data(),strip.data(),i+1));
	qBAe[i] += ((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qBA_%d",crange.data(),strip.data(),i+1)))->GetBinContent(1);
	qBAecnt[i]+=((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qBAcnt_%d",crange.data(),strip.data(),i+1)))->GetBinContent(1);
	qCAe[i] += ((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qCA_%d",crange.data(),strip.data(),i+1)))->GetBinContent(1);
	qCAecnt[i]+=((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qCAcnt_%d",crange.data(),strip.data(),i+1)))->GetBinContent(1);
	qCBe[i] += ((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qCB_%d",crange.data(),strip.data(),i+1)))->GetBinContent(1);
	qCBecnt[i]+=((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qCBcnt_%d",crange.data(),strip.data(),i+1)))->GetBinContent(1);
      }
    } else {
      centcnt+=centbins->GetBinContent(j);
      ptav->Add( (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/ptav",crange.data())));
      ptcnt->Add( (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/ptcnt",crange.data())));
      badcnt->Add( (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/badcnt",crange.data())));
      //centbins = (TH1D * ) fin->Get("vnanalyzer/centres");
      qA->Add( (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qA",crange.data(),strip.data())));
      qB->Add( (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qB",crange.data(),strip.data())));
      wnA->Add( (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/wnA",crange.data(),strip.data())));
      wnB->Add( (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/wnB",crange.data(),strip.data())));
      qBA += ((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qBA",crange.data(),strip.data())))->GetBinContent(1);
      qBAcnt+=((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qBAcnt",crange.data(),strip.data())))->GetBinContent(1);
      qCA += ((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qCA",crange.data(),strip.data())))->GetBinContent(1);
      qCAcnt+=((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qCAcnt",crange.data(),strip.data())))->GetBinContent(1); 
      qCB += ((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qCB",crange.data(),strip.data())))->GetBinContent(1);
      qCBcnt+=((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qCBcnt",crange.data(),strip.data())))->GetBinContent(1);
      
      for(int i = 0; i<10; i++) {
	qAe[i]->Add( (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qA_%d",crange.data(),strip.data(),i+1)));
	qBe[i]->Add( (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qB_%d",crange.data(),strip.data(),i+1)));
	wnAe[i]->Add( (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/wnA_%d",crange.data(),strip.data(),i+1)));
	wnBe[i]->Add( (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/wnB_%d",crange.data(),strip.data(),i+1)));
	qBAe[i] += ((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qBA_%d",crange.data(),strip.data(),i+1)))->GetBinContent(1);
	qBAecnt[i]+=((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qBAcnt_%d",crange.data(),strip.data(),i+1)))->GetBinContent(1);
	qCAe[i] += ((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qCA_%d",crange.data(),strip.data(),i+1)))->GetBinContent(1);
	qCAecnt[i]+=((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qCAcnt_%d",crange.data(),strip.data(),i+1)))->GetBinContent(1);
	qCBe[i] += ((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qCB_%d",crange.data(),strip.data(),i+1)))->GetBinContent(1);
	qCBecnt[i]+=((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qCBcnt_%d",crange.data(),strip.data(),i+1)))->GetBinContent(1);
      }
    }
  }
  qBA/=qBAcnt;
  qCA/=qCAcnt;
  qCB/=qCBcnt;
  for(int i = 0; i<10; i++) {
  qBAe[i]/=qBAecnt[i];
  qCAe[i]/=qCAecnt[i];
  qCBe[i]/=qCBecnt[i];
    
  }
  ptav->Divide(ptcnt);
  int ietamin1=0;
  int ietamax1=0;
  int ietamin2=0;
  int ietamax2=0;
  double sign = 1.;
  if(replay==N112ASUB2 || replay==N112ASUB3) sign=-1.;
  if(etamin*etamax<0) {
    ietamin1 = qA->GetYaxis()->FindBin(etamin);
    ietamax1 = qA->GetYaxis()->FindBin(-0.001);
    ietamin2 = qA->GetYaxis()->FindBin(0.);
    ietamax2 = qA->GetYaxis()->FindBin(etamax);
    qA1 = (TH1D *) qA->ProjectionX("qA1",ietamin1,ietamax1);
    qB1 = (TH1D *) qB->ProjectionX("qB1",ietamin2,ietamax2);
    wA1 = (TH1D *) wnA->ProjectionX("wA1",ietamin1,ietamax1);
    wB1 = (TH1D *) wnB->ProjectionX("wB1",ietamin2,ietamax2);
  } else {
    ietamin1 = qA->GetYaxis()->FindBin(etamin);
    ietamax1 = qA->GetYaxis()->FindBin(etamax);
    qA1 = (TH1D *) qA->ProjectionX("qA1",ietamin1,ietamax1);
    qB1 = (TH1D *) qB->ProjectionX("qB1",ietamin1,ietamax1);
    wA1 = (TH1D *) wnA->ProjectionX("wA1",ietamin1,ietamax1);
    wB1 = (TH1D *) wnB->ProjectionX("wB1",ietamin1,ietamax1);
  }
  qA->Divide(wnA);
  qB->Divide(wnB);
  qA1->Divide(wA1);
  qB1->Divide(wB1);
  if(!nonorm) {
    qBA=fabs(qBA);
    qCA=fabs(qCA);
    qCB=fabs(qCB);
    if(sub2) {
      qA->Scale(1./sqrt(qBA));
      qB->Scale(1./sqrt(qBA));
      qA1->Scale(1./sqrt(qBA));
      qB1->Scale(1./sqrt(qBA));
      resA[0]= sqrt(qBA);
      resB[0]= sqrt(qBA);
    } else {
      qA->Scale(1./sqrt(qBA*qCA/qCB));
      qB->Scale(1./sqrt(qBA*qCB/qCA));
      qA1->Scale(1./sqrt(qBA*qCA/qCB));
      qB1->Scale(1./sqrt(qBA*qCB/qCA));
      resA[0]= sqrt(qBA*qCA/qCB);
      resB[0]= sqrt(qBA*qCB/qCA);
    }
  }
  for(int i = 0; i<10; i++) {
    
    qAe1[i] = (TH1D *) qAe[i]->ProjectionX(Form("qAe1_%d",i),ietamin1,ietamax1);
    qBe1[i] = (TH1D *) qBe[i]->ProjectionX(Form("qBe1_%d",i),ietamin2,ietamax2);
    wAe1[i] = (TH1D *) wnAe[i]->ProjectionX(Form("wA1_%d",i),ietamin1,ietamax1);
    wBe1[i] = (TH1D *) wnBe[i]->ProjectionX(Form("wB1_%d",i),ietamin2,ietamax2);
    
    qAe[i]->Divide(wnAe[i]);
    qBe[i]->Divide(wnBe[i]);
    qAe1[i]->Divide(wAe1[i]);
    qBe1[i]->Divide(wBe1[i]);
    if(!nonorm) {
      qBAe[i]=fabs(qBAe[i]);
      qCAe[i]=fabs(qCAe[i]);
      qCBe[i]=fabs(qCBe[i]);
      if(sub2) {
	qAe[i]->Scale(1./sqrt(qBAe[i]));
	qBe[i]->Scale(1./sqrt(qBAe[i]));  
	resA[i+1]= sqrt(qBAe[i]);
	resB[i+1]= sqrt(qBAe[i]);
      } else {
	qAe[i]->Scale(1./sqrt(qBAe[i]*qCAe[i]/qCBe[i]));
	qBe[i]->Scale(1./sqrt(qBAe[i]*qCBe[i]/qCAe[i]));
	resA[i+1]= sqrt(qBAe[i]*qCAe[i]/qCBe[i]);
	resB[i+1]= sqrt(qBAe[i]*qCBe[i]/qCAe[i]);
      }
    }
  }
  if(etamin*etamax<0) {
    xpt = (TH1D *) ptav->ProjectionX("xpt",ietamin1,ietamax2);
    double c = (cmin[bin]+cmax[bin])/2.;
    TH2D * hsEff = ptcntEff(ptcnt,c);
  
    //hsEff->Draw("colz");

    sp = (TH1D *) hsEff->ProjectionX("sp",ietamin1,ietamax2);
    double ebinsA = ietamax1-ietamin1+1 ;
    double ebinsB = ietamax2-ietamin2+1;
    xpt->Scale(1./(ebinsA+ebinsB));
    sp->Scale(1./(ebinsA+ebinsB));
    if(!rew) {
      vnA = (TH1D *) qA->ProjectionX("vnA",ietamin1,ietamax1);
      vnB = (TH1D *) qB->ProjectionX("vnB",ietamin2,ietamax2);
      vn = (TH1D *) vnA->Clone("vn");
      vn->Add(vnB,sign);
      vn->Scale(sign);
      vnA->Scale(1./ebinsA);
      vnB->Scale(1./ebinsB);
      vn->Scale(1./(ebinsA+ebinsB));
    } else {
      vnA = (TH1D *) qA1->Clone("vnA");
      vnB = (TH1D *) qB1->Clone("vnB");
      vn = (TH1D *) vnA->Clone("vn");
      vn->Add(vnB,sign);
      vn->Scale(sign);
      vn->Scale(0.5);
    }
    double vnm[50] ={0};
    double vnAm[50] = {0};
    double vnBm[50] = {0};
    double vn2[50] ={0};
    double vnA2[50] = {0};
    double vnB2[50] = {0};
    for(int i = 0; i<50; i++) {
      vnm[i] = 0;
      vnAm[i] = 0;
      vnBm[i] = 0;
      vn2[i] = 0;
      vnA2[i] = 0;
      vnB2[i] = 0;
    }
    for(int i = 0; i<10; i++) {
      if(!rew) {
	vnAe = (TH1D *) qAe[i]->ProjectionX(Form("vnA%d",i),ietamin1,ietamax1);
	vnBe = (TH1D *) qBe[i]->ProjectionX(Form("vnB%d",i),ietamin2,ietamax2);
	vne = (TH1D *) vnAe->Clone(Form("vn%d",i));
	vne->Add(vnBe,sign);
	vn->Scale(sign);
	vnAe->Scale(1./ebinsA);
	vnBe->Scale(1./ebinsB);
	vne->Scale(1./(ebinsA+ebinsB));
      } else {
	vnAe = (TH1D *) qAe1[i]->Clone(Form("vnA%d",i));
	vnBe = (TH1D *) qBe1[i]->Clone(Form("vnB%d",i));
	vne = (TH1D *) vnAe->Clone(Form("vn%d",i));
	vne->Add(vnBe,sign);
	vne->Scale(sign);
	vne->Scale(0.5);
      }
      
      for(int j = 0; j<vne->GetNbinsX(); j++) {
	vnm[j]+= vne->GetBinContent(j+1);
	vnAm[j]+= vnAe->GetBinContent(j+1);
	vnBm[j]+= vnBe->GetBinContent(j+1);
	vn2[j] += pow(vne->GetBinContent(j+1),2);
	vnA2[j]+= pow(vnAe->GetBinContent(j+1),2);
	vnB2[j]+= pow(vnBe->GetBinContent(j+1),2);
      }
    }
    for(int j = 0; j<vn->GetNbinsX(); j++) {
      vnm[j]/=10.;
      vnAm[j]/=10.;
      vnBm[j]/=10.;
      vn2[j]/=10.;
      vnA2[j]/=10.;
      vnB2[j]/=10.;
      vn->SetBinError(j+1, sqrt((1./9.)*( vn2[j] - pow(vnm[j], 2))));
      vnA->SetBinError(j+1,sqrt((1./9.)*(vnA2[j] - pow(vnAm[j],2))));
      vnB->SetBinError(j+1,sqrt((1./9.)*(vnB2[j] - pow(vnBm[j],2))));
    }
  } else {
    xpt = (TH1D *) ptav->ProjectionX("xpt",ietamin1,ietamax1);
    double c = (cmin[bin]+cmax[bin])/2.;
    TH2D * hsEff = ptcntEff(ptcnt,c);
  
    //hsEff->Draw("colz");

    sp = (TH1D *) hsEff->ProjectionX("sp",ietamin1,ietamax1);
    double ebinsA = ietamax1-ietamin1+1 ;
    xpt->Scale(1./ebinsA);
    sp->Scale(1./ebinsA);
    if(!rew) {
      vnA = (TH1D *) qA->ProjectionX("vnA",ietamin1,ietamax1);
      vnB = (TH1D *) qB->ProjectionX("vnB",ietamin1,ietamax1);
      vn = (TH1D *) vnA->Clone("vn");
      vn->Add(vnB,sign);
      vn->Scale(sign);
      vnA->Scale(1./ebinsA);
      vnB->Scale(1./ebinsA);
      vn->Scale(1./(ebinsA+ebinsA));
    } else {
      vnA = (TH1D *) qA1->Clone("vnA");
      vnB = (TH1D *) qB1->Clone("vnB");
      vn = (TH1D *) vnA->Clone("vn");
      vn->Add(vnB,sign);
      vn->Scale(sign);
      vn->Scale(0.5);
    }
    double vnm[50] ={0};
    double vnAm[50] = {0};
    double vnBm[50] = {0};
    double vn2[50] ={0};
    double vnA2[50] = {0};
    double vnB2[50] = {0};
    for(int i = 0; i<50; i++) {
      vnm[i] = 0;
      vnAm[i] = 0;
      vnBm[i] = 0;
      vn2[i] = 0;
      vnA2[i] = 0;
      vnB2[i] = 0;
    }
    for(int i = 0; i<10; i++) {
      if(!rew) {
	vnAe = (TH1D *) qAe[i]->ProjectionX(Form("vnA%d",i),ietamin1,ietamax1);
	vnBe = (TH1D *) qBe[i]->ProjectionX(Form("vnB%d",i),ietamin1,ietamax1);
	vne = (TH1D *) vnAe->Clone(Form("vn%d",i));
	vne->Add(vnBe,sign);
	vne->Scale(sign);
	vnAe->Scale(1./ebinsA);
	vnBe->Scale(1./ebinsA);
	vne->Scale(1./(ebinsA+ebinsA));
      } else {
	vnAe = (TH1D *) qAe1[i]->Clone(Form("vnA%d",i));
	vnBe = (TH1D *) qBe1[i]->Clone(Form("vnB%d",i));
	vne = (TH1D *) vnAe->Clone(Form("vn%d",i));
	vne->Add(vnBe,sign);
	vne->Scale(sign);
	vne->Scale(0.5);
      }
      
      for(int j = 0; j<vne->GetNbinsX(); j++) {
	vnm[j]+= vne->GetBinContent(j+1);
	vnAm[j]+= vnAe->GetBinContent(j+1);
	vnBm[j]+= vnBe->GetBinContent(j+1);
	vn2[j] += pow(vne->GetBinContent(j+1),2);
	vnA2[j]+= pow(vnAe->GetBinContent(j+1),2);
	vnB2[j]+= pow(vnBe->GetBinContent(j+1),2);
      }
    }
    for(int j = 0; j<vn->GetNbinsX(); j++) {
      vnm[j]/=10.;
      vnAm[j]/=10.;
      vnBm[j]/=10.;
      vn2[j]/=10.;
      vnA2[j]/=10.;
      vnB2[j]/=10.;
      vn->SetBinError(j+1, sqrt((1./9.)*( vn2[j] - pow(vnm[j], 2))));
      vnA->SetBinError(j+1,sqrt((1./9.)*(vnA2[j] - pow(vnAm[j],2))));
      vnB->SetBinError(j+1,sqrt((1./9.)*(vnB2[j] - pow(vnBm[j],2))));
    }







  }
  TH1D * yld;
  if(ietamax2>0) {
  yld = ptcnt->ProjectionX("yld",ietamin1,ietamax2);
  } else {
  yld = ptcnt->ProjectionX("yld",ietamin1,ietamax1);
  }
  double x[20];
  double y[20];
  double yA[20];
  double yB[20];
  double ex[20];
  double ey[20];
  double eyA[20];
  double eyB[20];
  double xspec[20];
  double yspec[20];
  double exspec[20];
  double eyspec[20];
  for(int i = 0; i<20; i++) {
    x[i] = 0;
    y[i] = 0;
    yA[i] = 0;
    yB[i] = 0;
    ex[i] = 0;
    ey[i] = 0;
    eyA[i] = 0;
    eyB[i] = 0;
    xspec[i] = 0;
    yspec[i] = 0;
    exspec[i] = 0;
    eyspec[i] = 0;
  }
  int npt = 0;
  double wvn = 0;
  double wvne = 0;
  double w = 0;
  for(int i = 1; i<=xpt->GetNbinsX(); i++) {
    double pt = xpt->GetBinContent(i);
    if(sp->GetBinContent(i)<10) break;
    if(pt>=0.3 && pt < 12.) {
      x[npt]=pt;
      y[npt]=vn->GetBinContent(i);
      yA[npt]=vnA->GetBinContent(i);
      yB[npt]=vnB->GetBinContent(i);
      ex[npt] = 0;
      ey[npt] = vn->GetBinError(i);
      eyA[npt] = vnA->GetBinError(i);
      eyB[npt] = vnB->GetBinError(i);
      xspec[npt] = pt;
      yspec[npt] = sp->GetBinContent(i)/sp->GetBinWidth(i)/(etamax-etamin)/centcnt;
      exspec[npt] = 0;
      eyspec[npt] = 0;
      if(sp->GetBinContent(i)>1) eyspec[npt] = sqrt(sp->GetBinContent(i))/sp->GetBinWidth(i)/(etamax-etamin)/centcnt;
      if(pt<3.) {
	double eff = 0;
	double cent = (cmin[bin] + cmax[bin])/2.;
	
	double fake = FakeAndEff(cent,pt,eff);
	wvn += y[npt]*yld->GetBinContent(i)/eff;
	wvne+= ey[npt]*yld->GetBinContent(i)/eff;
	w   += yld->GetBinContent(i)/eff;;
      }
      ++npt;
    }
  }
  vint = wvn/w;
  vinte = wvne/w;
  cout<<"INTEGRAL: "<<vint<<"\t"<<vinte<<endl;
  TGraphErrors * g = new TGraphErrors(npt,x,y,ex,ey);
  g->SetMarkerStyle(20);
  g->SetMarkerColor(kBlue);
  g->SetLineColor(kBlue);
  g->SetLineWidth(2);
  gspec = new TGraphErrors(npt,xspec,yspec,exspec,eyspec);
  gspec->SetMarkerStyle(20);
  gspec->SetMarkerColor(kBlue);
  gspec->SetLineColor(kBlue);
  gspec->SetLineWidth(2);
  gA = new TGraphErrors(npt,x,yA,ex,eyA);
  gA->SetMarkerStyle(28);
  gA->SetMarkerColor(kMagenta);
  gA->SetLineColor(kMagenta);
  gA->SetLineWidth(2);
  gB = new TGraphErrors(npt,x,yB,ex,eyB);
  gB->SetMarkerStyle(28);
  gB->SetMarkerColor(kCyan);
  gB->SetLineColor(kCyan);
  gB->SetLineWidth(2);
  return g;

}
TH1D * h = 0;
void GetVNCreate(int replay = N42SUB3, int bin = 0, bool NumOnly=false, bool DenomOnly=false ){
  TH1D * hspec = 0;
  FILE * ftest;
  if(isTight) {
    if((ftest = fopen(Form("%s/%s",FigDir.data(),AnalNames[replay].data()),"r"))==NULL) {
      system(Form("mkdir %s/%s",FigDir.data(),AnalNames[replay].data()));
    } else {
      fclose(ftest);
    }

    if((ftest = fopen(Form("%s/%s/data",FigDir.data(),AnalNames[replay].data()),"r"))==NULL) {
      system(Form("mkdir %s/%s/data",FigDir.data(),AnalNames[replay].data()));
    } else {
      fclose(ftest);
    }
    system(Form("touch %s/%s/data/integral.dat",FigDir.data(),AnalNames[replay].data()));
    soutint = Form("%s/%s/data/integral.dat",FigDir.data(),AnalNames[replay].data());
  }
  string cname = AnalNames[replay]+"_"+to_string(cmin[bin])+"_"+to_string(cmax[bin]);
  bool chi4 = false;
  if(replay==CHI4) chi4 = true;
  bool chi5 = false;
  if(replay==CHI5) chi5 = true;
  bool chi62 = false;
  if(replay==CHI62) chi62 = true;
  bool chi63 = false;
  if(replay==CHI63) chi63 = true;
  bool chi7 = false;
  if(replay==CHI7) chi7 = true;
  bool chi7eff = false;
  if(replay==CHI7EFF) chi7eff = true;
  bool chi7noeff = false;
  if(replay==CHI7NOEFF) chi7noeff = true;
  TCanvas * c = new TCanvas(cname.data(),cname.data(),650,500);
  
  h = new TH1D("h","h",100,0,12.);
  h->SetDirectory(0);
  h->SetMinimum(0);
  h->SetMaximum(0.6);
  TGraphErrors * hpt;
  TGraphErrors * hdenom;
  TGraphErrors * hA;
  TGraphErrors * hB;
  TGraphErrors * hAdenom;
  TGraphErrors * hBdenom;
  TGraphErrors * nwspec;
  TGraphErrors * nwspec2;
  double vint = 0;
  double vinte = 0;
  double vintdenom = 0;
  double vintedenom = 0;
  double resA[cbins];
  double resB[cbins];
  double resAdenom[cbins];
  double resBdenom[cbins];
  cout<<"In create"<<endl;
  if(replay==N112ASUB2||replay==N112ASUB3) {
    hdenom = GetVNPt(N2SUB3, bin, EtaMin,EtaMax, hAdenom, hBdenom,nwspec2, resAdenom,resBdenom, vintdenom,vintedenom, false);
    fin->Close();
    fin = new TFile(rootFile.data(),"r");
    hpt    = GetVNPt(N112ASUB3, bin, EtaMin,EtaMax, hA, hB, nwspec, resA, resB, vint,vinte, false);
    double res = (resAdenom[0]+resBdenom[0])/2.;
      for(int i = 0; i<hpt->GetN(); i++) {
         double ef = hpt->GetEY()[i]/hpt->GetY()[i];
         hpt->GetY()[i]/=res;
         hpt->GetEY()[i]/=res;
      
         ef = hA->GetEY()[i]/hA->GetY()[i];
         hA->GetY()[i]/=resBdenom[0];
         hA->GetEY()[i]/=resBdenom[0];
      
         ef = hB->GetEY()[i]/hB->GetY()[i];
         hB->GetY()[i]/=resAdenom[0];
         hB->GetEY()[i]/=resAdenom[0];
       }
   
  } else if (replay==N123ASUB2||replay==N123ASUB3){
    hdenom = GetVNPt(N3SUB3, bin, EtaMin,EtaMax, hAdenom, hBdenom,nwspec2, resAdenom,resBdenom,vintdenom,vintedenom, false);
    fin->Close();
    fin = new TFile(rootFile.data(),"r");
    hpt    = GetVNPt(N123ASUB3, bin, EtaMin,EtaMax, hA, hB, nwspec, resA, resB, vint,vinte, false);
    double res = (resAdenom[0]+resBdenom[0])/2.;
    for(int i = 0; i<hpt->GetN(); i++) {
       double ef = hpt->GetEY()[i]/hpt->GetY()[i];
       hpt->GetY()[i]/=res;
       hpt->GetEY()[i]/=res;
      
       ef = hA->GetEY()[i]/hA->GetY()[i];
       hA->GetY()[i]/=resBdenom[0];
       hA->GetEY()[i]/=resBdenom[0];
      
       ef = hB->GetEY()[i]/hB->GetY()[i];
       hB->GetY()[i]/=resAdenom[0];
       hB->GetEY()[i]/=resAdenom[0];
     }
    
  } else if (replay==N42ASUB2||replay==N42ASUB3){
    hdenom = GetVNPt(N42SUB3, bin, EtaMin,EtaMax, hAdenom, hBdenom,nwspec2, resAdenom,resBdenom,vintdenom,vintedenom, false);
    fin->Close();
    fin = new TFile(rootFile.data(),"r");
    hpt    = GetVNPt(N42ASUB3, bin, EtaMin,EtaMax, hA, hB, nwspec, resA, resB, vint,vinte, true);
    double res = (resAdenom[0]+resBdenom[0])/2.;
    for(int i = 0; i<hpt->GetN(); i++) {
       double ef = hpt->GetEY()[i]/hpt->GetY()[i];
       hpt->GetY()[i]/=res;
       hpt->GetEY()[i]/=res;
      
       ef = hA->GetEY()[i]/hA->GetY()[i];
       hA->GetY()[i]/=resBdenom[0];
       hA->GetEY()[i]/=resBdenom[0];
      
       ef = hB->GetEY()[i]/hB->GetY()[i];
       hB->GetY()[i]/=resAdenom[0];
       hB->GetEY()[i]/=resAdenom[0];
     }
  } else if (replay==D24ASUB2||replay==D24ASUB3){
    hdenom = GetVNPt(D24SUB3, bin, EtaMin,EtaMax, hAdenom, hBdenom,nwspec2, resAdenom,resBdenom,vintdenom,vintedenom, false);
    fin->Close();
    fin = new TFile(rootFile.data(),"r");
    hpt    = GetVNPt(D24ASUB3, bin, EtaMin,EtaMax, hA, hB, nwspec, resA, resB, vint,vinte, true);
    double res = (resAdenom[0]+resBdenom[0])/2.;
    for(int i = 0; i<hpt->GetN(); i++) {
       double ef = hpt->GetEY()[i]/hpt->GetY()[i];
       hpt->GetY()[i]/=res;
       hpt->GetEY()[i]/=res;
      
       ef = hA->GetEY()[i]/hA->GetY()[i];
       hA->GetY()[i]/=resBdenom[0];
       hA->GetEY()[i]/=resBdenom[0];
      
       ef = hB->GetEY()[i]/hB->GetY()[i];
       hB->GetY()[i]/=resAdenom[0];
       hB->GetEY()[i]/=resAdenom[0];
     }
    
  } else if (replay==N62ASUB2||replay==N62ASUB3){
    hdenom = GetVNPt(N62SUB3, bin, EtaMin,EtaMax, hAdenom, hBdenom,nwspec2, resAdenom,resBdenom,vintdenom,vintedenom, false);
    fin->Close();
    fin = new TFile(rootFile.data(),"r");
    hpt    = GetVNPt(N62ASUB3, bin, EtaMin,EtaMax, hA, hB, nwspec, resA, resB, vint,vinte, true);
    double res = (resAdenom[0]+resBdenom[0])/2.;
    for(int i = 0; i<hpt->GetN(); i++) {
       double ef = hpt->GetEY()[i]/hpt->GetY()[i];
       hpt->GetY()[i]/=res;
       hpt->GetEY()[i]/=res;
      
       ef = hA->GetEY()[i]/hA->GetY()[i];
       hA->GetY()[i]/=resBdenom[0];
       hA->GetEY()[i]/=resBdenom[0];
      
       ef = hB->GetEY()[i]/hB->GetY()[i];
       hB->GetY()[i]/=resAdenom[0];
       hB->GetEY()[i]/=resAdenom[0];
     }
  } else if (replay==D26ASUB2||replay==D26ASUB3){
    hdenom = GetVNPt(D26SUB3, bin, EtaMin,EtaMax, hAdenom, hBdenom,nwspec2, resAdenom,resBdenom,vintdenom,vintedenom, false);
    fin->Close();
    fin = new TFile(rootFile.data(),"r");
    hpt    = GetVNPt(D26ASUB3, bin, EtaMin,EtaMax, hA, hB, nwspec, resA, resB, vint,vinte, true);
    double res = (resAdenom[0]+resBdenom[0])/2.;
    for(int i = 0; i<hpt->GetN(); i++) {
       double ef = hpt->GetEY()[i]/hpt->GetY()[i];
       hpt->GetY()[i]/=res;
       hpt->GetEY()[i]/=res;
      
       ef = hA->GetEY()[i]/hA->GetY()[i];
       hA->GetY()[i]/=resBdenom[0];
       hA->GetEY()[i]/=resBdenom[0];
      
       ef = hB->GetEY()[i]/hB->GetY()[i];
       hB->GetY()[i]/=resAdenom[0];
       hB->GetEY()[i]/=resAdenom[0];
     }
  } else if (replay==N63ASUB3||replay==N63ASUB3){
    hdenom = GetVNPt(N63SUB3, bin, EtaMin,EtaMax, hAdenom, hBdenom,nwspec2, resAdenom,resBdenom,vintdenom,vintedenom, false);
    fin->Close();
    fin = new TFile(rootFile.data(),"r");
    hpt    = GetVNPt(N63ASUB3, bin, EtaMin,EtaMax, hA, hB, nwspec, resA, resB, vint,vinte, true);
    double res = (resAdenom[0]+resBdenom[0])/2.;
    for(int i = 0; i<hpt->GetN(); i++) {
       double ef = hpt->GetEY()[i]/hpt->GetY()[i];
       hpt->GetY()[i]/=res;
       hpt->GetEY()[i]/=res;
      
       ef = hA->GetEY()[i]/hA->GetY()[i];
       hA->GetY()[i]/=resBdenom[0];
       hA->GetEY()[i]/=resBdenom[0];
      
       ef = hB->GetEY()[i]/hB->GetY()[i];
       hB->GetY()[i]/=resAdenom[0];
       hB->GetEY()[i]/=resAdenom[0];
     }
  } else if (replay==D34ASUB2||replay==D34ASUB3){
    hdenom = GetVNPt(D34SUB3, bin, EtaMin,EtaMax, hAdenom, hBdenom,nwspec2, resAdenom,resBdenom,vintdenom,vintedenom, false);
    fin->Close();
    fin = new TFile(rootFile.data(),"r");
    hpt    = GetVNPt(D34ASUB3, bin, EtaMin,EtaMax, hA, hB, nwspec, resA, resB, vint,vinte, true);
    double res = (resAdenom[0]+resBdenom[0])/2.;
    for(int i = 0; i<hpt->GetN(); i++) {
       double ef = hpt->GetEY()[i]/hpt->GetY()[i];
       hpt->GetY()[i]/=res;
       hpt->GetEY()[i]/=res;
      
       ef = hA->GetEY()[i]/hA->GetY()[i];
       hA->GetY()[i]/=resBdenom[0];
       hA->GetEY()[i]/=resBdenom[0];
      
       ef = hB->GetEY()[i]/hB->GetY()[i];
       hB->GetY()[i]/=resAdenom[0];
       hB->GetEY()[i]/=resAdenom[0];
     }
    
    
  } else if(chi4) {
    hdenom = GetVNPt(D24SUB3, bin, EtaMin,EtaMax, hAdenom, hBdenom,nwspec2, resAdenom, resBdenom, vintdenom,vintedenom, true);
    fin->Close();
    fin = new TFile(rootFile.data(),"r");
    hpt    = GetVNPt(N42SUB3, bin, EtaMin,EtaMax, hA, hB, nwspec,resA,resB,  vint,vinte, true);
    vint/=vintdenom;
    vinte/=vintdenom;
    if(!NumOnly && !DenomOnly) { 
      for(int i = 0; i<hpt->GetN(); i++) {
    	double ef = hpt->GetEY()[i]/hpt->GetY()[i];
    	hpt->GetY()[i]/=hdenom->GetY()[i];
    	hpt->GetEY()[i]=hpt->GetY()[i]*sqrt(ef*ef+pow(hdenom->GetEY()[i]/hdenom->GetY()[i],2));
	
    	ef = hA->GetEY()[i]/hA->GetY()[i];
    	hA->GetY()[i]/=hAdenom->GetY()[i];
    	hA->GetEY()[i]=hA->GetY()[i]*sqrt(ef*ef+pow(hAdenom->GetEY()[i]/hAdenom->GetY()[i],2));
	
    	ef = hB->GetEY()[i]/hB->GetY()[i];
    	hB->GetY()[i]/=hBdenom->GetY()[i];
    	hB->GetEY()[i]=hB->GetY()[i]*sqrt(ef*ef+pow(hBdenom->GetEY()[i]/hBdenom->GetY()[i],2));
      }
    }
  } else if(replay==CHI4A) {
    hdenom = GetVNPt(D24ASUB3, bin, EtaMin,EtaMax, hAdenom, hBdenom,nwspec2, resAdenom, resBdenom, vintdenom,vintedenom, true);
    fin->Close();
    fin = new TFile(rootFile.data(),"r");
    hpt    = GetVNPt(N42ASUB3, bin, EtaMin,EtaMax, hA, hB, nwspec,resA,resB,  vint,vinte, true);
    vint/=vintdenom;
    vinte/=vintdenom;
    if(!NumOnly && !DenomOnly) { 
      for(int i = 0; i<hpt->GetN(); i++) {
    	double ef = hpt->GetEY()[i]/hpt->GetY()[i];
    	hpt->GetY()[i]/=hdenom->GetY()[i];
    	hpt->GetEY()[i]=hpt->GetY()[i]*sqrt(ef*ef+pow(hdenom->GetEY()[i]/hdenom->GetY()[i],2));
	
    	ef = hA->GetEY()[i]/hA->GetY()[i];
    	hA->GetY()[i]/=hAdenom->GetY()[i];
    	hA->GetEY()[i]=hA->GetY()[i]*sqrt(ef*ef+pow(hAdenom->GetEY()[i]/hAdenom->GetY()[i],2));
	
    	ef = hB->GetEY()[i]/hB->GetY()[i];
    	hB->GetY()[i]/=hBdenom->GetY()[i];
    	hB->GetEY()[i]=hB->GetY()[i]*sqrt(ef*ef+pow(hBdenom->GetEY()[i]/hBdenom->GetY()[i],2));
      }
    }
  } else if(chi5) {
    hdenom = GetVNPt(D2232SUB3, bin, EtaMin,EtaMax, hAdenom, hBdenom, nwspec2, resAdenom, resBdenom,vintdenom,vintedenom, true);
    fin->Close();
    fin = new TFile(rootFile.data(),"r");
    hpt = GetVNPt(N523SUB3, bin, EtaMin,EtaMax, hA, hB,nwspec,resA,resB, vint,vinte, true);
    vint/=vintdenom;
    vinte/=vintdenom;
    if(DenomOnly) {
      hA=hAdenom;
      hB=hBdenom;
      hpt=hdenom;
    }
    if(!NumOnly && !DenomOnly) {
      for(int i = 0; i<hpt->GetN(); i++) {
	double ef = hpt->GetEY()[i]/hpt->GetY()[i];
	hpt->GetY()[i]/=hdenom->GetY()[i];
	hpt->GetEY()[i]=hpt->GetY()[i]*sqrt(ef*ef+pow(hdenom->GetEY()[i]/hdenom->GetY()[i],2));
	
	ef = hA->GetEY()[i]/hA->GetY()[i];
	hA->GetY()[i]/=hAdenom->GetY()[i];
	hA->GetEY()[i]=hA->GetY()[i]*sqrt(ef*ef+pow(hAdenom->GetEY()[i]/hAdenom->GetY()[i],2));
	
      ef = hB->GetEY()[i]/hB->GetY()[i];
      hB->GetY()[i]/=hBdenom->GetY()[i];
      hB->GetEY()[i]=hB->GetY()[i]*sqrt(ef*ef+pow(hBdenom->GetEY()[i]/hBdenom->GetY()[i],2));
      }
    }
  } else if(replay==CHI5A) {
    hdenom = GetVNPt(D2232ASUB3, bin, EtaMin,EtaMax, hAdenom, hBdenom, nwspec2, resAdenom, resBdenom,vintdenom,vintedenom, true);
    fin->Close();
    fin = new TFile(rootFile.data(),"r");
    hpt = GetVNPt(N523ASUB3, bin, EtaMin,EtaMax, hA, hB,nwspec,resA,resB, vint,vinte, true);
    vint/=vintdenom;
    vinte/=vintdenom;
    if(DenomOnly) {
      hA=hAdenom;
      hB=hBdenom;
      hpt=hdenom;
    }
    if(!NumOnly && !DenomOnly) {
      for(int i = 0; i<hpt->GetN(); i++) {
	double ef = hpt->GetEY()[i]/hpt->GetY()[i];
	hpt->GetY()[i]/=hdenom->GetY()[i];
	hpt->GetEY()[i]=hpt->GetY()[i]*sqrt(ef*ef+pow(hdenom->GetEY()[i]/hdenom->GetY()[i],2));
	
	ef = hA->GetEY()[i]/hA->GetY()[i];
	hA->GetY()[i]/=hAdenom->GetY()[i];
	hA->GetEY()[i]=hA->GetY()[i]*sqrt(ef*ef+pow(hAdenom->GetEY()[i]/hAdenom->GetY()[i],2));
	
      ef = hB->GetEY()[i]/hB->GetY()[i];
      hB->GetY()[i]/=hBdenom->GetY()[i];
      hB->GetEY()[i]=hB->GetY()[i]*sqrt(ef*ef+pow(hBdenom->GetEY()[i]/hBdenom->GetY()[i],2));
      }
    }
  } else if(chi62) {
    hpt = GetVNPt(N62SUB3, bin, EtaMin,EtaMax, hA, hB,nwspec, resA,resB, vint,vinte, true);
    fin->Close();
    fin = new TFile(rootFile.data(),"r");
    hdenom = GetVNPt(D26SUB3, bin, EtaMin,EtaMax, hAdenom, hBdenom,nwspec2,resAdenom,resBdenom, vintdenom,vintedenom, true);
    vint/=vintdenom;
    vinte/=vintdenom;
    for(int i = 0; i<hpt->GetN(); i++) {
      double ef = hpt->GetEY()[i]/hpt->GetY()[i];
      hpt->GetY()[i]/=hdenom->GetY()[i];
      hpt->GetEY()[i]=hpt->GetY()[i]*sqrt(ef*ef+pow(hdenom->GetEY()[i]/hdenom->GetY()[i],2));
      
      ef = hA->GetEY()[i]/hA->GetY()[i];
      hA->GetY()[i]/=hAdenom->GetY()[i];
      hA->GetEY()[i]=hA->GetY()[i]*sqrt(ef*ef+pow(hAdenom->GetEY()[i]/hAdenom->GetY()[i],2));
      
      ef = hB->GetEY()[i]/hB->GetY()[i];
      hB->GetY()[i]/=hBdenom->GetY()[i];
      hB->GetEY()[i]=hB->GetY()[i]*sqrt(ef*ef+pow(hBdenom->GetEY()[i]/hBdenom->GetY()[i],2));
    }
  } else if(replay==CHI62A) {
    hpt = GetVNPt(N62ASUB3, bin, EtaMin,EtaMax, hA, hB,nwspec, resA,resB, vint,vinte, true);
    fin->Close();
    fin = new TFile(rootFile.data(),"r");
    hdenom = GetVNPt(D26ASUB3, bin, EtaMin,EtaMax, hAdenom, hBdenom,nwspec2,resAdenom,resBdenom, vintdenom,vintedenom, true);
    vint/=vintdenom;
    vinte/=vintdenom;
    for(int i = 0; i<hpt->GetN(); i++) {
      double ef = hpt->GetEY()[i]/hpt->GetY()[i];
      hpt->GetY()[i]/=hdenom->GetY()[i];
      hpt->GetEY()[i]=hpt->GetY()[i]*sqrt(ef*ef+pow(hdenom->GetEY()[i]/hdenom->GetY()[i],2));
      
      ef = hA->GetEY()[i]/hA->GetY()[i];
      hA->GetY()[i]/=hAdenom->GetY()[i];
      hA->GetEY()[i]=hA->GetY()[i]*sqrt(ef*ef+pow(hAdenom->GetEY()[i]/hAdenom->GetY()[i],2));
      
      ef = hB->GetEY()[i]/hB->GetY()[i];
      hB->GetY()[i]/=hBdenom->GetY()[i];
      hB->GetEY()[i]=hB->GetY()[i]*sqrt(ef*ef+pow(hBdenom->GetEY()[i]/hBdenom->GetY()[i],2));
    }
  } else if(chi63) {
    hpt = GetVNPt(N63SUB3, bin, EtaMin,EtaMax, hA, hB,nwspec, resA, resB,vint,vinte, true);
    fin->Close();
    fin = new TFile(rootFile.data(),"r");
    hdenom = GetVNPt(D34SUB3, bin, EtaMin,EtaMax, hAdenom, hBdenom,nwspec2, resAdenom,resBdenom, vintdenom,vintedenom, true);
    vint/=vintdenom;
    vinte/=vintdenom;
    for(int i = 0; i<hpt->GetN(); i++) {
      double ef = hpt->GetEY()[i]/hpt->GetY()[i];
      hpt->GetY()[i]/=hdenom->GetY()[i];
      hpt->GetEY()[i]=hpt->GetY()[i]*sqrt(ef*ef+pow(hdenom->GetEY()[i]/hdenom->GetY()[i],2));
      
      ef = hA->GetEY()[i]/hA->GetY()[i];
      hA->GetY()[i]/=hAdenom->GetY()[i];
      hA->GetEY()[i]=hA->GetY()[i]*sqrt(ef*ef+pow(hAdenom->GetEY()[i]/hAdenom->GetY()[i],2));
      
      ef = hB->GetEY()[i]/hB->GetY()[i];
      hB->GetY()[i]/=hBdenom->GetY()[i];
      hB->GetEY()[i]=hB->GetY()[i]*sqrt(ef*ef+pow(hBdenom->GetEY()[i]/hBdenom->GetY()[i],2));
    }
  } else if(replay==CHI63A) {
    hpt = GetVNPt(N63ASUB3, bin, EtaMin,EtaMax, hA, hB,nwspec, resA, resB,vint,vinte, true);
    fin->Close();
    fin = new TFile(rootFile.data(),"r");
    hdenom = GetVNPt(D34ASUB3, bin, EtaMin,EtaMax, hAdenom, hBdenom,nwspec2, resAdenom,resBdenom, vintdenom,vintedenom, true);
    vint/=vintdenom;
    vinte/=vintdenom;
    for(int i = 0; i<hpt->GetN(); i++) {
      double ef = hpt->GetEY()[i]/hpt->GetY()[i];
      hpt->GetY()[i]/=hdenom->GetY()[i];
      hpt->GetEY()[i]=hpt->GetY()[i]*sqrt(ef*ef+pow(hdenom->GetEY()[i]/hdenom->GetY()[i],2));
      
      ef = hA->GetEY()[i]/hA->GetY()[i];
      hA->GetY()[i]/=hAdenom->GetY()[i];
      hA->GetEY()[i]=hA->GetY()[i]*sqrt(ef*ef+pow(hAdenom->GetEY()[i]/hAdenom->GetY()[i],2));
      
      ef = hB->GetEY()[i]/hB->GetY()[i];
      hB->GetY()[i]/=hBdenom->GetY()[i];
      hB->GetEY()[i]=hB->GetY()[i]*sqrt(ef*ef+pow(hBdenom->GetEY()[i]/hBdenom->GetY()[i],2));
    }
  } else if(chi7) {
    hpt = GetVNPt(N723SUB3, bin, EtaMin,EtaMax, hA, hB,nwspec, resA,resB, vint,vinte,true);
    fin->Close();
    fin = new TFile(rootFile.data(),"r");
    hdenom = GetVNPt(D2432SUB3, bin, EtaMin,EtaMax, hAdenom, hBdenom,nwspec2,resAdenom,resBdenom, vintdenom,vintedenom, true);
    vint/=vintdenom;
    vinte/=vintdenom;
    for(int i = 0; i<hpt->GetN(); i++) {
      double ef = hpt->GetEY()[i]/hpt->GetY()[i];
      hpt->GetY()[i]/=hdenom->GetY()[i];
      hpt->GetEY()[i]=hpt->GetY()[i]*sqrt(ef*ef+pow(hdenom->GetEY()[i]/hdenom->GetY()[i],2));
      
      ef = hA->GetEY()[i]/hA->GetY()[i];
      hA->GetY()[i]/=hAdenom->GetY()[i];
      hA->GetEY()[i]=hA->GetY()[i]*sqrt(ef*ef+pow(hAdenom->GetEY()[i]/hAdenom->GetY()[i],2));
      
      ef = hB->GetEY()[i]/hB->GetY()[i];
      hB->GetY()[i]/=hBdenom->GetY()[i];
      hB->GetEY()[i]=hB->GetY()[i]*sqrt(ef*ef+pow(hBdenom->GetEY()[i]/hBdenom->GetY()[i],2));
    }
  } else if(replay==CHI7A) {
    hpt = GetVNPt(N723ASUB3, bin, EtaMin,EtaMax, hA, hB,nwspec, resA,resB, vint,vinte,true);
    fin->Close();
    fin = new TFile(rootFile.data(),"r");
    hdenom = GetVNPt(D2432ASUB3, bin, EtaMin,EtaMax, hAdenom, hBdenom,nwspec2,resAdenom,resBdenom, vintdenom,vintedenom, true);
    vint/=vintdenom;
    vinte/=vintdenom;
    for(int i = 0; i<hpt->GetN(); i++) {
      double ef = hpt->GetEY()[i]/hpt->GetY()[i];
      hpt->GetY()[i]/=hdenom->GetY()[i];
      hpt->GetEY()[i]=hpt->GetY()[i]*sqrt(ef*ef+pow(hdenom->GetEY()[i]/hdenom->GetY()[i],2));
      
      ef = hA->GetEY()[i]/hA->GetY()[i];
      hA->GetY()[i]/=hAdenom->GetY()[i];
      hA->GetEY()[i]=hA->GetY()[i]*sqrt(ef*ef+pow(hAdenom->GetEY()[i]/hAdenom->GetY()[i],2));
      
      ef = hB->GetEY()[i]/hB->GetY()[i];
      hB->GetY()[i]/=hBdenom->GetY()[i];
      hB->GetEY()[i]=hB->GetY()[i]*sqrt(ef*ef+pow(hBdenom->GetEY()[i]/hBdenom->GetY()[i],2));
    }
  } else if(chi7eff) {
    hpt = GetVNPt(N723EFF, bin, EtaMin,EtaMax, hA, hB, nwspec,resA, resB, vint,vinte,true);
    fin->Close();
    fin = new TFile(rootFile.data(),"r");
    hdenom = GetVNPt(D2432EFF, bin, EtaMin,EtaMax, hAdenom, hBdenom, nwspec2,resAdenom, resBdenom, vintdenom,vintedenom,true);
    vint/=vintdenom;
    vinte/=vintdenom;
    for(int i = 0; i<hpt->GetN(); i++) {
      double ef = hpt->GetEY()[i]/hpt->GetY()[i];
      hpt->GetY()[i]/=hdenom->GetY()[i];
      hpt->GetEY()[i]=hpt->GetY()[i]*sqrt(ef*ef+pow(hdenom->GetEY()[i]/hdenom->GetY()[i],2));
      
      ef = hA->GetEY()[i]/hA->GetY()[i];
      hA->GetY()[i]/=hAdenom->GetY()[i];
      hA->GetEY()[i]=hA->GetY()[i]*sqrt(ef*ef+pow(hAdenom->GetEY()[i]/hAdenom->GetY()[i],2));
      
      ef = hB->GetEY()[i]/hB->GetY()[i];
      hB->GetY()[i]/=hBdenom->GetY()[i];
      hB->GetEY()[i]=hB->GetY()[i]*sqrt(ef*ef+pow(hBdenom->GetEY()[i]/hBdenom->GetY()[i],2));
    }
  } else if(chi7noeff) {
    hpt = GetVNPt(N723NOEFF, bin, EtaMin,EtaMax, hA, hB,nwspec,resA, resB, vint,vinte,true);
    fin->Close();
    fin = new TFile(rootFile.data(),"r");
    hdenom = GetVNPt(D2432NOEFF, bin, EtaMin,EtaMax, hAdenom, hBdenom,nwspec2,resAdenom,resBdenom, vintdenom,vintedenom,true);
    vint/=vintdenom;
    vinte/=vintdenom;
    for(int i = 0; i<hpt->GetN(); i++) {
      double ef = hpt->GetEY()[i]/hpt->GetY()[i];
      hpt->GetY()[i]/=hdenom->GetY()[i];
      hpt->GetEY()[i]=hpt->GetY()[i]*sqrt(ef*ef+pow(hdenom->GetEY()[i]/hdenom->GetY()[i],2));
      
      ef = hA->GetEY()[i]/hA->GetY()[i];
      hA->GetY()[i]/=hAdenom->GetY()[i];
      hA->GetEY()[i]=hA->GetY()[i]*sqrt(ef*ef+pow(hAdenom->GetEY()[i]/hAdenom->GetY()[i],2));
      
      ef = hB->GetEY()[i]/hB->GetY()[i];
      hB->GetY()[i]/=hBdenom->GetY()[i];
      hB->GetEY()[i]=hB->GetY()[i]*sqrt(ef*ef+pow(hBdenom->GetEY()[i]/hBdenom->GetY()[i],2));
    }
  } else {
    hpt = GetVNPt(replay, bin,EtaMin,EtaMax, hA, hB,nwspec,resA,resB, vint,vinte);
  }
  outint = fopen(soutint.data(),"a+");

  fprintf(outint,"%d\t%d\t%15.10f\t%15.10f\n",cmin[bin],cmax[bin],vint,vinte);
  fclose(outint);
  double ymax = -1;
  double ymin = 1;
  double ymaxspec = 0;
  for(int n = 0; n<hpt->GetN(); n++) {
    if(hpt->GetX()[n]>4) break;
    if(hpt->GetY()[n]+hpt->GetEY()[n]>ymax) ymax = hpt->GetY()[n]+hpt->GetEY()[n];
    if(hpt->GetY()[n]-hpt->GetEY()[n]<ymin) ymin = hpt->GetY()[n]-hpt->GetEY()[n];
   if(hA->GetY()[n]+hA->GetEY()[n]>ymax) ymax = hA->GetY()[n]+hA->GetEY()[n];
    if(hA->GetY()[n]-hA->GetEY()[n]<ymin) ymin = hA->GetY()[n]-hA->GetEY()[n];
   if(hB->GetY()[n]+hB->GetEY()[n]>ymax) ymax = hB->GetY()[n]+hB->GetEY()[n];
    if(hB->GetY()[n]-hB->GetEY()[n]<ymin) ymin = hB->GetY()[n]-hB->GetEY()[n];

    if(nwspec->GetY()[n]+nwspec->GetEY()[n]>ymaxspec) ymaxspec = nwspec->GetY()[n]+nwspec->GetEY()[n];
  }

  if(ymax<0.0002) {
    ymax=0.0002;
  } else if(ymax<0.004) {
    ymax=0.001;
  } else if(ymax<0.01) {
    ymax = 0.02;
  } else if(ymax<0.02) {
    ymax = 0.03;
  } else if(ymax<0.05) {
    ymax = 0.07;
  } else if (ymax<0.1) {
    ymax = 0.2;
  } else if (ymax<0.4) {
    ymax = 0.6;
  } else if( ymax<0.8) {
    ymax = 1.0;
  } else if( ymax < 2.) {
    ymax = 2.0;
  } else {
    ymax =10;
  }

  if(ymin>0) {
    ymin=0;
  }  else if (ymin>-0.003) {
    ymin=-0.003;
  }  else if (ymin>-0.005) {
    ymin=-0.005;
  }  else if (ymin>-0.01) {
    ymin=-0.01;
  } else if(ymin>-0.02) {
    ymin = -0.02;
  } else if(ymin>-0.04) {
    ymin = -0.04;
  } else if(ymin>-0.1) {
    ymin = -0.1;
  } else if (ymin>-0.2) {
    ymin = -0.2;
  } else if (ymin>-0.4) {
    ymin = -0.4;
  } else if( ymin>-0.8) {
    ymin = -0.8;
  } else if( ymin > -2) {
    ymin = -2;
  } else {
    ymin =-1;
  }
  h->SetMinimum(ymin);
  h->SetMaximum(ymax);
  gPad->SetGrid(1,1);
  h->Draw();
  h->SetXTitle("P_{T} (GeV/c)");
  string numdenom = "";
  if(NumOnly) numdenom=" (Numerator) ";
  if(DenomOnly) numdenom=" (Denominator) ";
  string yt = ytitle[replay]+numdenom+" ("+to_string(cmin[bin])+" - "+to_string(cmax[bin])+"%)";
  h->SetYTitle(yt.data());
  hpt->Draw("p");
  string prevname = "";
  string shengquan = "";
  if(ANAL==N2SUB3 || ANAL==N2SUB2 || ANAL==N2EFF) {
    prevname = Form("data/EP_10_002_PtDists/v2_pt_ep_cen%d_%d_eta08.txt",cmin[bin],cmax[bin]);
    if(cmin[bin]<=50) shengquan = Form("data/forSteveNov11/v22_%d_%d.txt",cmin[bin],cmax[bin]);
  }
  if(ANAL==N42SUB3 || ANAL==N42SUB2 || ANAL==N42EFF || ANAL==N42ASUB3 || ANAL==N42ASUB2) {
    prevname = Form("data/hin_11_005_data/EPResults/PtDists/v4_2_%d_%d.txt",cmin[bin],cmax[bin]);
    if(cmin[bin]<=50) shengquan = Form("data/forSteveNov11/v42_%d_%d.txt",cmin[bin],cmax[bin]);
  }
  if(ANAL==N3SUB3 || ANAL==N3SUB2) prevname = Form("data/hin_11_005_data/EPResults/PtDists/v3_%d_%d.txt",cmin[bin],cmax[bin]);
  if(ANAL==N4SUB3 || ANAL==N4SUB2) prevname = Form("data/hin_11_005_data/EPResults/PtDists/v4_%d_%d.txt",cmin[bin],cmax[bin]);
  if(ANAL==N5SUB3 || ANAL==N5SUB2) prevname = Form("data/hin_11_005_data/EPResults/PtDists/v5_%d_%d.txt",cmin[bin],cmax[bin]);
  if(ANAL==N6SUB3 || ANAL==N6SUB2) prevname = Form("data/hin_11_005_data/EPResults/PtDists/v6_%d_%d.txt",cmin[bin],cmax[bin]);
  if(ANAL==N62SUB3 || ANAL==N62SUB2||ANAL==N62ASUB3 || ANAL==N62ASUB2) {
    prevname = Form("data/hin_11_005_data/EPResults/PtDists/v6_2_%d_%d.txt",cmin[bin],cmax[bin]);
    if(cmin[bin]<=50) shengquan = Form("data/forSteveNov11/v62_%d_%d.txt",cmin[bin],cmax[bin]);
   }
  if(ANAL==N63SUB3 || ANAL==N63SUB2 || ANAL==N63ASUB2 || ANAL==N63ASUB3) {
    if(cmin[bin]<=50) shengquan = Form("data/forSteveNov11/v63_%d_%d.txt",cmin[bin],cmax[bin]);
   }
  if(ANAL==N723SUB3 || ANAL==N723SUB2 ||ANAL==N723ASUB3 || ANAL==N723ASUB2 || ANAL==N723EFF || ANAL==N723NOEFF) {
    if(cmin[bin]<=50) shengquan = Form("data/forSteveNov11/v723_%d_%d.txt",cmin[bin],cmax[bin]);
  }
  if(ANAL==N523ASUB3 || ANAL==N523ASUB2 || ANAL==N523SUB3 || ANAL==N523SUB2 ) {
    if(cmin[bin]<=50) shengquan = Form("data/forSteveNov11/N523_%d_%d.txt",cmin[bin],cmax[bin]);
  }
  if(ANAL==D2232ASUB3 || ANAL==D2232ASUB2 ) {
    if(cmin[bin]<=50) shengquan = Form("data/forSteveNov11/D2232_%d_%d.txt",cmin[bin],cmax[bin]);
  }
  if(ANAL==D2432ASUB3 || ANAL==D2432ASUB2 ) {
    if(cmin[bin]<=50) shengquan = Form("data/forSteveNov11/D2432_%d_%d.txt",cmin[bin],cmax[bin]);
  }
  if(ANAL==D24SUB3 || ANAL==D24SUB2 || ANAL==D24ASUB3 || ANAL==D24ASUB2 ) {
    if(cmin[bin]<=50) shengquan = Form("data/forSteveNov11/D24_%d_%d.txt",cmin[bin],cmax[bin]);
  }
  if(ANAL==D26SUB3 || ANAL==D26SUB2 ) {
    if(cmin[bin]<=50) shengquan = Form("data/forSteveNov11/D26_%d_%d.txt",cmin[bin],cmax[bin]);
  }
  if(replay == CHI4||replay == CHI4A) {
    if(cmin[bin]<=50) shengquan = Form("data/forSteveNov11/Chi42_%d_%d.txt",cmin[bin],cmax[bin]);
    prevname = "";
  }
  if(replay == CHI5||replay == CHI5A) {
    if(cmin[bin]<=50) {
      shengquan = Form("data/forSteveNov11/Chi523_%d_%d.txt",cmin[bin],cmax[bin]);
      if(NumOnly) shengquan = Form("data/forSteveNov11/N523_%d_%d.txt",cmin[bin],cmax[bin]);
    }
    prevname = "";
 
  }
  if(replay == CHI62||replay == CHI62A) {
    if(cmin[bin]<=50) shengquan = Form("data/forSteveNov11/Chi62_%d_%d.txt",cmin[bin],cmax[bin]);
    prevname = "";
  }
  if(replay == CHI63||replay == CHI63A) {
    if(cmin[bin]<=50) shengquan = Form("data/forSteveNov11/Chi63_%d_%d.txt",cmin[bin],cmax[bin]);
    prevname = "";
  }
  if(replay == CHI7 ||replay == CHI7A ||replay == CHI7EFF) {
    if(cmin[bin]<=50) shengquan = Form("data/forSteveNov11/Chi723_%d_%d.txt",cmin[bin],cmax[bin]);
    prevname = "";
  }
  
  TLegend * leg = new TLegend(0.65,0.65,0.9,0.9);
  leg->SetTextFont(43);
  leg->SetTextSize(20);
  leg->SetFillColor(kWhite);
  leg->SetBorderSize(0);
  leg->AddEntry(hpt,AnalNames[replay].data(),"lp");
  leg->AddEntry(hA,"HF+ only","lp");
  leg->AddEntry(hB,"HF- only","lp");
  cout<<"leg formed"<<endl;
  if(prevname.length()>1 && bin<11) { 
    double x[40];
    double y[40];
    double stat[40];
    double sys[40];
    FILE * fin = fopen(prevname.data(),"r");
    char buf[80];
    int n = 0;
    while(fgets(buf,80,fin)!=NULL) {
      sscanf(buf,"%lf\t%lf\t%lf\t%lf",&x[n],&y[n],&stat[n],&sys[n]);
      ++n;
    }
    TGraphErrors * gold = new TGraphErrors(n,x,y,0,stat);
    gold->SetMarkerStyle(25);
    gold->SetMarkerColor(kRed);
    gold->SetLineColor(kRed);
    gold->Draw("p");
    leg->AddEntry(gold,"CMS Published","lp");
  }

  if(shengquan.length()>1&&bin<11) { 
    double x[40];
    double y[40];
    double stat[40];
    double sys[40];
    FILE * fin = fopen(shengquan.data(),"r");
    char buf[80];
    int n = 0;
    while(fgets(buf,80,fin)!=NULL) {
      sscanf(buf,"%lf\t%lf\t%lf",&x[n],&y[n],&stat[n]);
      x[n]+=0.05;
      ++n;
    }
    TGraphErrors * sheng = new TGraphErrors(n,x,y,0,stat);
    sheng->SetMarkerStyle(24);
    sheng->SetMarkerColor(kGreen);
    sheng->SetLineColor(kGreen);
    sheng->Draw("p");
    leg->AddEntry(sheng,"Shengquan (offset)","lp");
  }
  if(replay==N523SUB3) {
    TGraphErrors * hpt2 = GetVNPt(N523SUB3, bin,EtaMin,EtaMax, hA, hB,nwspec,resA, resB,vint,vinte);
    hpt2->SetMarkerColor(kGreen);
    hpt2->SetLineColor(kGreen);
    hpt2->Draw("p");
    leg->AddEntry(hpt2,"mixed HF","lp");
  }
  if(replay==N523ASUB3) {
    TGraphErrors * hpt2 = GetVNPt(N523ASUB3, bin,EtaMin,EtaMax, hA, hB,nwspec,resA,resB,vint,vinte);
    hpt2->SetMarkerColor(kGreen);
    hpt2->SetLineColor(kGreen);
    hpt2->Draw("p");
    leg->AddEntry(hpt2,"mixed HF","lp");
  }
  if(replay==N723SUB3) {
    //TGraphErrors * hpt2 = GetVNPt(N723ASUB3, bin,EtaMin,EtaMax, hA, hB, nwspec2,resA,resB,vint,vinte);
    //hpt2->SetMarkerColor(kGreen);
    //hpt2->SetLineColor(kGreen);
    //hpt2->Draw("p");
    //leg->AddEntry(hpt2,"mixed HF","lp");
  }
  if(replay==D2232SUB3) {
    TGraphErrors * hpt2 = GetVNPt(D2232ASUB3, bin,EtaMin,EtaMax, hA, hB,nwspec,resA,resB,vint,vinte);
    hpt2->SetMarkerColor(kGreen);
    hpt2->SetLineColor(kGreen);
    hpt2->Draw("p");
    leg->AddEntry(hpt2,"mixed HF","lp");
  }
  if(replay==D2432SUB3) {
    //TGraphErrors * hpt2 = GetVNPt(D2432ASUB3, bin,EtaMin,EtaMax,hA, hB,nwspec2,resA,resB,vint,vintey);
    //hpt2->SetMarkerColor(kGreen);
    // hpt2->SetLineColor(kGreen);
    // hpt2->Draw("p");
    // leg->AddEntry(hpt2,"mixed HF","lp");
  }
  leg->Draw();
  hA->Draw("p");
  hB->Draw("p");
  hpt->Draw("p");
  TLatex * text = new TLatex(1,0.87*ymax,AnalNames[replay].data());
  text->SetTextFont(43);
  text->SetTextSize(28);
  text->Draw();
  TLatex * t2 = new TLatex(1,0.77*ymax,Form("%d - %d%c",cmin[bin],cmax[bin],'%'));
  t2->SetTextFont(43);
  t2->SetTextSize(22);
  t2->Draw();
  TLatex * t3 = new TLatex(0.8,0.67*ymax,Form("%4.1f < #eta < %4.1f",EtaMin,EtaMax));
  t3->SetTextFont(43);
  t3->SetTextSize(22);
  t3->Draw();
  FILE * fout;
  if(isTight) {
    c->Print(Form("%s/%s/%s.pdf",FigDir.data(),AnalNames[replay].data(),cname.data()),"pdf");
    fout = fopen(Form("%s/%s/data/%s.dat",FigDir.data(),AnalNames[replay].data(),cname.data()),"w");
  } else {
    c->Print(Form("%s/%s/%s.pdf",FigDir.data(),AnalNames[replay].data(),cname.data()),"pdf");
    fout = fopen(Form("%s/%s/data/%s.dat",FigDir.data(),AnalNames[replay].data(),cname.data()),"w");
  }
  for(int i = 0; i<hpt->GetN(); i++){
    fprintf(fout,"%5.3f\t%9.7f\t%9.7f\n",hpt->GetX()[i],hpt->GetY()[i],hpt->GetEY()[i]);
  }
  fclose(fout);
  bool drawSpec = true;
  if(drawSpec){
    string c2name = "c2_"+AnalNames[replay]+"_"+to_string(cmin[bin])+"_"+to_string(cmax[bin]);
    TCanvas * c2 = new TCanvas(c2name.data(),c2name.data(),700,500);
    hspec = new TH1D("hspec","hspec",600,0,12);
    hspec->SetDirectory(0);
    hspec->SetMaximum(100*pow(10.,(double)((int) TMath::Log10(ymaxspec))));
    hspec->SetMinimum(0.00001);
    cout<<"c2: "<<c2<<endl;
    c2->cd();
    gPad->SetLogy();
    hspec->Draw();
    hspec->SetXTitle("p_{T} (GeV/c)");
    hspec->SetYTitle("1/(N_{ev}) d^{2}N/dp_{T}d#eta");
    nwspec->Draw("p");
    double ym = hspec->GetMaximum();
    TLatex * t4 = new TLatex(8,0.1*ym,AnalNames[replay].data());
    t4->SetTextFont(43);
    t4->SetTextSize(28);
    t4->Draw();
    TLatex * t6 = new TLatex(8,0.02*ym,Form("%d - %d%c",cmin[bin],cmax[bin],'%'));
    t6->SetTextFont(43);
    t6->SetTextSize(22);
    t6->Draw();
    TLatex * t7 = new TLatex(7.8,0.004*ym,Form("%4.1f < #eta < %4.1f",EtaMin,EtaMax));
    t7->SetTextFont(43);
    t7->SetTextSize(22);
    t7->Draw();
    if(isTight) {
      c2->Print(Form("%s/%s/%s.pdf",FigDir.data(),AnalNames[replay].data(),c2name.data()),"pdf");
      sspec = FigDir+"/"+AnalNames[replay]+"/data/spec_"+to_string(cmin[bin])+"_"+to_string(cmax[bin])+".dat";
    } else {
      c2->Print(Form("%s/%s/%s.pdf",FigDir.data(),AnalNames[replay].data(),c2name.data()),"pdf");
      sspec = FigDir+"/"+AnalNames[replay]+"/data/spec_"+to_string(cmin[bin])+"_"+to_string(cmax[bin])+".dat";
    }
    }
    if(isTight) {
    sspec = FigDir+"/"+AnalNames[replay]+"/data/spec_"+to_string(cmin[bin])+"_"+to_string(cmax[bin])+".dat";
    } else {
    sspec = FigDir+"/"+AnalNames[replay]+"/data/spec_"+to_string(cmin[bin])+"_"+to_string(cmax[bin])+".dat";
    }
  outspec = fopen(sspec.data(),"w");
  for(int i = 0; i<nwspec->GetN(); i++) fprintf(outspec,"%7.2f\t%12.5f\t%12.5f\n",nwspec->GetX()[i],nwspec->GetY()[i],nwspec->GetEY()[i]);

}

void GetVN(string name="N523ASUB3", string tag="useTight", double mineta = -0.8, double maxeta = 0.8, bool override = false){
  bool found = false;
  centRef = new TH1I("centRef","centRef",11,centRefBins);
  EtaMin = mineta;
  EtaMax = maxeta;
  stag = "_"+tag;
  rootFile = "";
  if(tag=="useTight" || tag=="noRecenter" || tag=="noRecenter_reweight") {
    if(tag=="useTight") stag = "";
    isTight   = true;
    isNominal = false;
    isWide    = false;
    isNarrow  = false;
    isLoose   = false;
    if(tag=="useTight") {
      rootFile = "../MH.root";
      fin = new TFile("../MH.root","read");
    }else {
      rootFile = "../MH_noRecenter.root";
      fin = new TFile("../MH_noRecenter.root","read");
      if(tag=="noRecenter_reweight") rew=true;
    }
  }
  if(tag=="useNominal") {
    isTight   = false;
    isNominal = true;
    isWide    = false;
    isNarrow  = false;
    isLoose   = false;
  }
  int en = 0;
  for(int indx = 0; indx<LAST; ++indx){
    if(AnalNames[indx]==name) {
      found = true;
      en = indx;
      break;
    }
  }
  if(!found) {
    cout<<"Failed to locate analysis type"<<endl;
    return;
  }
  FILE * ftest;
  FigDir = Form("figures%s_%04.1f_%04.1f",stag.data(),EtaMin,EtaMax);
  FigSubDir = FigDir+"/"+name.data();
  if(fopen(FigDir.data(),"r")!=NULL) {
    cout<<"Output directory "<<FigDir.data()<<" exists."<<endl;
  } else {
    system(Form("mkdir %s",FigDir.data()));
  }
  if((ftest=fopen(FigSubDir.data(),"r"))==NULL) {
    system(Form("mkdir %s",FigSubDir.data()));
  } else {
    if(override) {
      system(Form("rm -rf %s",FigSubDir.data()));
    } else {
      cout<<"Directory "<<FigSubDir.data()<<" exists.  ABORT."<<endl;
      return;
    }
    fclose(ftest);
  }

  for(int bin = 0; bin<13; bin++) {
   GetVNCreate(en,bin);
    fin->Close();
    fin = new TFile(rootFile.data(),"read");
   }
}
