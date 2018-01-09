#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TMath.h"
#include "TStopwatch.h"
#include <iostream>

#include "src/HiEvtPlaneList.h"
using namespace hi;
static const int ncentbinsNOFF = 25;
static const int centBinsNOFF[]={0, 10, 20, 30, 40, 50, 60, 70, 80, 100, 120, 135, 150, 160, 185, 210, 230, 250, 270, 300, 330, 350, 370, 390, 420, 500};
static const double centRefBinsNOFF[]={0, 10, 20, 30, 40, 50, 60, 70, 80, 100, 120, 135, 150, 160, 185, 210, 230, 250, 270, 300, 330, 350, 370, 390, 420, 500};
static const int cbinsNOFF = 25;
static const int cminNOFF[]={1, 10, 20, 30, 40, 50, 60, 70, 80, 100, 120, 135, 150, 160, 185, 210, 230, 250, 270, 300, 330, 350, 370, 390, 420};
static const int cmaxNOFF[]={10, 20, 30, 40, 50, 60, 70, 80, 100, 120, 135, 150, 160, 185, 210, 230, 250, 270, 300, 330, 350, 370, 390, 420, 500};
static const int ncentbinsCENT = 11;
static const int centBinsCENT[]={0,5,10,15,20,25,30,35,40,50,60,70};
static const double centRefBinsCENT[]={0,5,10,15,20,25,30,35,40,50,60,70};
static const int cbinsCENT = 13;
static const int cminCENT[]={0, 5,10,15,20,25,30,35,40,50,60,  0,20, 60};
static const int cmaxCENT[]={5,10,15,20,25,30,35,40,50,60,70, 20,60,100};
bool ntrkbinning = false;
int ncentbins = 0;
int centBins[50];
double centRefBins[50];
int cbins = 0;
int cmin[50];
int cmax[50];

string tag;
bool SetToA = false;
bool Decor;
double EtaMin = -0.8;
double EtaMax = 0.8;
string FigDir = "";
string FigSubDir = "";
string FigSubSubDir = "";
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
string stag;
TH1D * centbins;
TH1I * centRef;
TH1D * rcnt;
string rootFile;
double resA[50];
double resB[50];
double resAdenom[50];
double resBdenom[50];
#include "src/Types.h"
#include "src/Efficiency.C"
#include "src/GetVNPt.C"
#include "src/Harmonics.h"
#include "src/drawSpec.C"
TH1D * h = 0;
void GetVNCreate(int replay , int bin , TGraphErrors * & gint, TGraphErrors * & gintA, TGraphErrors *& gintB,  bool plotit=true, bool NumOnly=false, bool DenomOnly=false ){
  TH1D * hspec = 0;
  string cname = ANALS[replay][0]+"_"+to_string(cmin[bin])+"_"+to_string(cmax[bin])+Form("_eta_%03.1f_%03.1f",EtaMin,EtaMax);
  h = new TH1D("h","h",100,0,12.);
  h->SetDirectory(0);
  h->SetMinimum(0);

  h->SetMaximum(0.6);
  TGraphErrors * g;
  TGraphErrors * gA;
  TGraphErrors * gB;
  TGraphErrors * gspec;
  TGraphErrors * nwspec2;
  double vint = 0;
  double vinte = 0;
  double vintdenom = 0;
  double vintedenom = 0;
  double ymin = 0;
  double ymax = 1;
  if(replay==N1SUB2   || replay==N1SUB3)   g =   N1(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N1ASUB2   || replay==N1ASUB3)   g =   N1(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N1BSUB2   || replay==N1BSUB3)   g =   N1(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N1MCm22SUB3 || replay==N1MCm18SUB3 || replay==N1MCm14SUB3 || replay==N1MCm10SUB3 ||
     replay==N1MCm06SUB3 || replay==N1MCm02SUB3 || replay==N1MCp22SUB3 || replay==N1MCp18SUB3 ||
     replay==N1MCp14SUB3 || replay==N1MCp10SUB3 ||  replay==N1MCp06SUB3 || replay==N1MCp02SUB3 )
    g =   N1EVEN(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N1MCm22SUB2 || replay==N1MCm18SUB2 || replay==N1MCm14SUB2 || replay==N1MCm10SUB2 ||
     replay==N1MCm06SUB2 || replay==N1MCm02SUB2 || replay==N1MCp22SUB2 || replay==N1MCp18SUB2 ||
     
     replay==N1MCp14SUB2 || replay==N1MCp10SUB2 ||  replay==N1MCp06SUB2 || replay==N1MCp02SUB2 )
    g =   N1EVEN(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N1EVENSUB2 || replay==N1EVENSUB3) g =   N1EVEN(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N112SUB2   || replay==N112SUB3)   g =   N112(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N112ASUB2   || replay==N112ASUB3)   g =   N112(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N112BSUB2   || replay==N112BSUB3)   g =   N112(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  // if(replay==N123SUB2   || replay==N123SUB3)   g =   N123(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  // if(replay==N123ASUB2   || replay==N123ASUB3)   g =   N123(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  // if(replay==N123BSUB2   || replay==N123BSUB3)   g =   N123(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N2SUB2   || replay==N2SUB3)   g =   N2(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N3SUB2   || replay==N3SUB3)   g =   N3(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N4SUB2   || replay==N4SUB3)   g =   N4(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N5SUB2   || replay==N5SUB3)   g =   N5(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N6SUB2   || replay==N6SUB3)   g =   N6(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N7SUB2   || replay==N7SUB3)   g =   N7(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N523SUB2 || replay==N523SUB3) g = N523(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N523ASUB2 || replay==N523ASUB3) g = N523(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N42SUB2   || replay==N42SUB3)   g =   N42(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N42ASUB2   || replay==N42ASUB3)   g =   N42(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if( replay==N42BSUB3)   g =   N42(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if( replay==N42CSUB3)   g =   N42(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N62SUB2   || replay==N62SUB3)   g =   N62(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N62ASUB3)   g =   N62(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N63SUB2   || replay==N63SUB3)   g =   N63(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N63ASUB2 || replay==N63ASUB3)   g =   N63(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if( replay==N63BSUB3)   g =   N63(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if( replay==N63CSUB3)   g =   N63(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if( replay==Chi4)    g = CHI4(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if( replay==Chi4A)   g = CHI4(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if( replay==Chi5)    g = CHI5(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if( replay==Chi5A)   g = CHI5(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if( replay==Chi62)    g = CHI62(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if( replay==Chi62A)   g = CHI62(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if( replay==Chi63)    g = CHI63(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if( replay==Chi63A)   g = CHI63(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if( replay==Chi7)    g = CHI7(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if( replay==Chi7A)   g = CHI7(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N723SUB2 || replay==N723SUB3) g = N723(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);
  if(replay==N723ASUB2 || replay==N723ASUB3) g = N723(replay,bin,EtaMin,EtaMax,ymin,ymax,g,gA,gB,gspec,gint,gintA, gintB);

  double ymaxspec = 0;
  h->SetMinimum(ymin);
  h->SetMaximum(ymax);
  TCanvas * c=NULL;
  if(plotit) {
    c = new TCanvas(cname.data(),cname.data(),650,500);
    gPad->SetGrid(1,1);
    h->Draw();

    h->SetXTitle("P_{T} (GeV/c)");
  }
  string numdenom = "";
  if(NumOnly) numdenom=" (Numerator) ";
  if(DenomOnly) numdenom=" (Denominator) ";
 

  string yt = ANALS[replay][1]+numdenom+" ("+to_string(cmin[bin])+" #leq N_{trk}^{off} < "+to_string(cmax[bin])+")";
  if(!ntrkbinning) yt = ANALS[replay][1]+numdenom+" ("+to_string(cmin[bin])+" - "+to_string(cmax[bin])+"%)";

  FILE * fout;
  if(plotit) {

    h->SetYTitle(yt.data());
    g->Draw("p");  

    TLegend * leg = new TLegend(0.65,0.65,0.9,0.9);
    leg->SetTextFont(43);

    leg->SetTextSize(20);
    leg->SetFillColor(kWhite);
    leg->SetBorderSize(0);
    if(!SetToA) {
      leg->AddEntry(g,"HF both sides","lp");
    } else {
      leg->AddEntry(g,"A side only is good","lp");
    }
    leg->AddEntry(gA,"A only","lp");
    leg->AddEntry(gB,"B only","lp");
    leg->Draw();

    gA->Draw("p");
    gB->Draw("p");

    TLatex * text = new TLatex(1,0.87*ymax,ANALS[replay][0].data());
    text->SetTextFont(43);
    text->SetTextSize(28);
    text->Draw();
    TLatex * t2;
    if(ntrkbinning) {
      t2 = new TLatex(1,0.77*ymax,Form("%d #leq N_{tkr}^{off} < %d",cmin[bin],cmax[bin]));
    } else {
      t2 = new TLatex(1,0.77*ymax,Form("%d - %d%c",cmin[bin],cmax[bin],'%'));
    }
    t2->SetTextFont(43);
    t2->SetTextSize(22);
    t2->Draw();
    TLatex * t3 = new TLatex(0.8,0.67*ymax,Form("%4.1f < #eta < %4.1f",EtaMin,EtaMax));
    t3->SetTextFont(43);
    t3->SetTextSize(22);
    t3->Draw();
    c->Print(Form("%s/%s.pdf",FigSubSubDir.data(), cname.data()),"pdf");
  }
  fout = fopen(Form("%s/data/%s.dat",FigSubSubDir.data(),cname.data()),"w");
  for(int i = 0; i<g->GetN(); i++){
    fprintf(fout,"%5.3f\t%9.7f\t%9.7f\n",g->GetX()[i],g->GetY()[i],g->GetEY()[i]);
  }
  fclose(fout);

  fout = fopen(Form("%s/data/%s_A.dat",FigSubSubDir.data(),cname.data()),"w");
  for(int i = 0; i<gA->GetN(); i++){
    fprintf(fout,"%5.3f\t%9.7f\t%9.7f\n",gA->GetX()[i],gA->GetY()[i],gA->GetEY()[i]);
  }
  fclose(fout);

  fout = fopen(Form("%s/data/%s_B.dat",FigSubSubDir.data(),cname.data()),"w");
  for(int i = 0; i<gB->GetN(); i++){
    fprintf(fout,"%5.3f\t%9.7f\t%9.7f\n",gB->GetX()[i],gB->GetY()[i],gB->GetEY()[i]);
  }

  fclose(fout);

 
}
void GetVN(string rootfile = "../MH.root", string name="N2SUB3",  double mineta = -0.8, double maxeta = 0.8, bool decor = false){
  TStopwatch * timer = new TStopwatch();
  bool found = false;
  Decor = decor;
  string nlabel = name;
  if(name.find("SUB2")!=std::string::npos) Decor = false;
  if(name.find("N1MC")!=std::string::npos) Decor = false;
  if(name.find("N523")!=std::string::npos) Decor = false;
  if(name.find("N42")!=std::string::npos) Decor = false;
  if(name.find("N62")!=std::string::npos) Decor = false;
  if(name.find("N63")!=std::string::npos) Decor = false;
  if(Decor) nlabel+="_decor";
  rootFile = rootfile;
  SetTracking();
  tag = rootfile.substr(rootfile.find("/")+1,rootfile.find(".root")-rootfile.find("/")-1);
  TGraphErrors * gint[cbins];
  TGraphErrors * gintA[cbins];
  TGraphErrors * gintB[cbins];
  double x[12];
  double y[12];
  double ey[12];
  for(int i = 0; i<12; i++) x[i]=-2.2+0.4*i;
  for(int i = 0; i<cbins; i++) {
    gint[i] = new TGraphErrors(12,x,y,0,ey);
    gintA[i] = new TGraphErrors(12,x,y,0,ey);
    gintB[i] = new TGraphErrors(12,x,y,0,ey);
    gint[i]->SetMarkerStyle(20);
    gintA[i]->SetMarkerStyle(24);
    gintB[i]->SetMarkerStyle(24);
    gintA[i]->SetMarkerColor(kRed);
    gintB[i]->SetMarkerColor(kBlue);
    gintA[i]->SetLineColor(kRed);
    gintB[i]->SetLineColor(kBlue);
    gint[i]->SetName(Form("sum_%s",nlabel.data()));
    gintA[i]->SetName(Form("A_%s",nlabel.data()));
    gintB[i]->SetName(Form("B_%s",nlabel.data()));
    gint[i]->SetTitle(Form("sum_%s",nlabel.data()));
    gintA[i]->SetTitle(Form("A_%s",nlabel.data()));
    gintB[i]->SetTitle(Form("B_%s",nlabel.data()));
  }

  centRef = new TH1I(Form("centRef_%s",nlabel.data()),Form("centRef_%s",nlabel.data()),ncentbins,centRefBins);
  EtaMin = mineta;
  EtaMax = maxeta;
  stag = "_"+tag;
  if(tag=="") stag = "";
  isTight   = false;
  isNominal = true;
  isWide    = false;
  isNarrow  = false;
  isLoose   = false;

  int en = 0;
  for(int indx = 0; indx<LAST; ++indx){
    if(ANALS[indx][0]==name) {
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
  FigDir = Form("figures%s",stag.data());
  if((ftest=fopen(FigDir.data(),"r"))!=NULL) {
    //cout<<"Output directory "<<FigDir.data()<<" exists."<<endl;
    fclose(ftest);
  } else {
    system(Form("mkdir %s",FigDir.data()));
  }
  FigSubDir = FigDir+"/"+name.data();
  if((ftest=fopen(FigSubDir.data(),"r"))==NULL) {
    system(Form("mkdir %s",FigSubDir.data()));
  } else {
    cout<<"Directory "<<FigSubDir.data()<<" exists.  Will overwrite."<<endl;
    fclose(ftest);
  }
  TCanvas * ceta[cbins];
  timer->Start();
  for(int bin = 0; bin<cbins; bin++) {
    timer->Stop();
    
    int minb = rcnt->FindBin(cmin[bin]);
    int maxb = rcnt->FindBin(cmax[bin])-1;
    if(maxb<minb) maxb=minb;
    int cnt = rcnt->Integral(minb,maxb);
    cout<<timer->CpuTime()<<"\t"<<cmin[bin]<<"\t"<<cmax[bin]<<"\t"<<cnt<<endl;
    timer->ResetCpuTime();
    timer->Start();
    if(cnt<5000) continue;
    //pt distribution
    FigSubSubDir = FigSubDir+Form("/eta_%04.1f_%04.1f",EtaMin,EtaMax);
    if(Decor) FigSubSubDir+="_decor";
    if((ftest=fopen(FigSubSubDir.data(),"r"))==NULL) {
      system(Form("mkdir %s",FigSubSubDir.data()));
      system(Form("mkdir %s/data",FigSubSubDir.data()));
      system(Form("touch %s/data/integral.dat",FigSubSubDir.data()));
      soutint = Form("%s/data/integral.dat",FigSubSubDir.data());
    } else {
      soutint = Form("%s/data/integral.dat",FigSubSubDir.data());
      fclose(ftest);
    }
    GetVNCreate(en,bin,gint[bin],gintA[bin],gintB[bin]);
    //eta distribution
    string FigEtaSubDir = FigSubDir;
    FigEtaSubDir+="/EtaDistributions";
    if(Decor) FigEtaSubDir+="_decor";
    if((ftest=fopen(FigEtaSubDir.data(),"r"))==NULL) {
      system(Form("mkdir %s",FigEtaSubDir.data()));
      system(Form("mkdir %s/data",FigEtaSubDir.data()));
    } else {
      fclose(ftest);
    }
    
    ceta[bin] = new TCanvas(Form("EtaInt_%s_%d_%d",nlabel.data(),cmin[bin],cmax[bin]),Form("EtaInt_%s_%d_%d",nlabel.data(),cmin[bin],cmax[bin]),800,500);
    double xmin,xmax,ymin,ymax;
    double xminA,xmaxA,yminA,ymaxA;
    double xminB,xmaxB,yminB,ymaxB;
    string nl2 = ANALS[en][0];
    nl2+=" ("+to_string(cmin[bin])+"-"+to_string(cmax[bin]);
    if(!ntrkbinning) nl2+="%";

    nl2+=")";
    string nl3 = ANALS[en][1];
    nl3+=" ("+to_string(cmin[bin])+"-"+to_string(cmax[bin]);
    if(!ntrkbinning) nl3+="%";

    nl3+=")";
    
    TH1D * heta = new TH1D(Form("heta_%s",nl2.data()),Form("heta_%s",nl2.data()),100,-2.5,2.5);
    gint[bin]->ComputeRange(xmin,ymin,xmax,ymax);
    gintA[bin]->ComputeRange(xminA,yminA,xmaxA,ymaxA);
    gintB[bin]->ComputeRange(xminB,yminB,xmaxB,ymaxB);
    ymin = plotmin(min(ymin,yminA));
    //ymin = min(ymin,yminB);
    ymax = plotmax(max(ymax,ymaxA));
    //ymax = max(ymax,ymaxB);

    heta->SetMaximum(ymax);
    heta->SetMinimum(ymin);
    heta->SetXTitle("#eta");
    heta->SetYTitle(ANALS[en][1].data());
    heta->Draw();
    gint[bin]->Draw("p");
    gintA[bin]->Draw("p");
    gintB[bin]->Draw("p");
    TLegend * leg2 = new TLegend(0.2,0.2,0.3,0.4);
    leg2->SetTextFont(43);
    leg2->SetTextSize(20);
    leg2->SetFillColor(kWhite);
    leg2->SetBorderSize(0);
    leg2->AddEntry(gint[bin],"Adopted","lp");
    leg2->AddEntry(gintA[bin],"A side","lp");
    leg2->AddEntry(gintB[bin],"B side","lp");
    leg2->Draw();
    TLatex * tl = new TLatex( -1,0.3*(ymax-ymin)+ymin,nl3.data());
    tl->Draw();
    if(ANALS[en][2]!="") {
      string tmp = ANALS[en][2];
      if(Decor) tmp+=" (EP Decorrelation corrected)";
      TLatex * tl2 = new TLatex( -1,0.2*(ymax-ymin)+ymin,tmp.data());
      tl2->SetTextFont(43);
      tl2->SetTextSize(16);
      tl2->Draw();
    }
    ceta[bin]->Print(Form("%s/%s.pdf",FigEtaSubDir.data(),ceta[bin]->GetName()),"pdf");

    // FILE * fint = fopen(Form("%s/data/%s.dat",FigEtaSubDir.data(),ceta[bin]->GetName()),"w");
    // for(int i = 0; i<gint[bin]->GetN(); i++) {
    //   fprintf(fint,"%5.1f\t%10.6f\t%10.6f\n",gint[bin]->GetX()[i],gint[bin]->GetY()[i],gint[bin]->GetEY()[i]);
    // }
    }
}
