#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include <iostream>
static const int ncentbins = 11;
static const int centBins[]={0,5,10,15,20,25,30,35,40,50,60,70};

enum AnalType {
  N2SUB2,       N2SUB3,      N3SUB2,     N3SUB3,     N4SUB2,      N4SUB3,  
  N42SUB2,      N42SUB3,     N5SUB2,     N5SUB3,     N6SUB2,      N6SUB3,   
  N7SUB2,       N7SUB3,      N523SUB2,   N523SUB3,   N723SUB2,    N723SUB3,
  N723ASUB2,    N723ASUB3,   N62SUB2,    N62SUB3,    N63SUB2,     N63SUB3,  
  D24SUB2,      D24SUB3,     D34SUB2,    D34SUB3,    D2232SUB2,   D2232SUB3,
  D2432SUB2,    D2432SUB3,   D2232ASUB2, D2232ASUB3, D2432ASUB2,  D2432ASUB3,
  N523ASUB2,    N523ASUB3,   D26SUB2,    D26SUB3,    CHI4,        CHI5, 
  CHI63,        CHI7,        N2EFF,      N2NOEFF,    N723EFF,     N723NOEFF,
  D2432EFF,     D2432NOEFF,  CHI7EFF,    CHI7NOEFF,  N42EFF
};
string AnalNames[]={
  "N2SUB2",       "N2SUB3",      "N3SUB2",     "N3SUB3",     "N4SUB2",      "N4SUB3",  
  "N42SUB2",      "N42SUB3",     "N5SUB2",     "N5SUB3",     "N6SUB2",      "N6SUB3",   
  "N7SUB2",       "N7SUB3",      "N523SUB2",   "N523SUB3",   "N723SUB2",    "N723SUB3",
  "N723ASUB2",    "N723ASUB3",   "N62SUB2",    "N62SUB3",    "N63SUB2",     "N63SUB3",  
  "D24SUB2",      "D24SUB3",     "D34SUB2",    "D34SUB3",    "D2232SUB2",   "D2232SUB3",
  "D2432SUB2",    "D2432SUB3",   "D2232ASUB2", "D2232ASUB3", "D2432ASUB2",  "D2432ASUB3",
  "N523ASUB2",    "N523ASUB3",   "D26SUB2",    "D26SUB3",    "CHI4",        "CHI5", 
  "CHI63",        "CHI7",        "N2EFF",      "N2NOEFF",    "N723EFF",     "N723NOEFF",
  "D2432EFF",     "D2432NOEFF",  "CHI7EFF",    "CHI7NOEFF",  "N42EFF"
};
string ytitle[]={
  "v_{2}",                   "v_{2}",                   "v_{3}",               "v_{3}",                 "v_{4}",               "v_{4}",
  "v_{4}\{#Psi_{2}\}",       "v_{4}\{#Psi_{2}\}",       "v_{5}",               "v_{5}",                 "v_{6}",               "v_{6}",
  "v_{7}",                   "v_{7}",                   "v_{5}\{#Psi_{2,3}\}", "v_{5}\{#Psi_{2,3}\}",   "v_{7}\{#Psi_{2,3}\}", "v_{7}\{#Psi_{2,3}\}",
  "v_{7}\{#Psi_{2,3^{*}\}",  "v_{7}\{#Psi_{2,3^{*}}\}", "v_{6}\{#Psi_{2}\}",   "v_{6}\{#Psi_{2}\}",     "v_{6}\{#Psi_{3}\}",    "v_{6}\{#Psi_{3}\}",
  "D24SUB2",                 "D24SUB3",                 "D34SUB2",             "D34SUB3",               "D2232SUB2",            "D2232SUB3",
  "D2432SUB2",               "D2432SUB3",               "D2232ASUB2",          "D2232ASUB3",            "D2432ASUB2",           "D2432ASUB3",
  "v_{5}\{#Psi_{2,3^{*}}\}", "v_{5}\{#Psi_{2,3^{*}}\}", "D26SUB2",             "D26SUB3",               "#chi_{4}",             "#chi_{5}", 
  "#chi_{63}",               "#chi_{7}",                "v_{2}(eff)",          "v_{2}(noeff)",          "v_{723}(eff)",          "v_{723}(noeff)",
  "D2432EFF",                "D2432NOEFF",              "#chi_{7}(eff)",       "#chi_{7}(moeff)",       "v_{42}(eff)"
};
TFile * fin;
double fw = 1.3;
int ANAL = -1;
TGraphErrors * GetVNPt(string file){
  TGraphErrors * g = new TGraphErrors(file.data(),"%lg%lg%lg");
  g->SetMarkerStyle(20);
  g->SetMarkerColor(kBlue);
  g->SetLineColor(kBlue);
  return g;

}
void GetVNCreate(int replayA = N42SUB3, int replayB = N42SUB3, int bin = 0){
  string baseA = Form("figures/%s/data/",AnalNames[replayA].data());
  string nameA = baseA+AnalNames[replayA]+"_"+to_string(centBins[bin])+"_"+to_string(centBins[bin+1])+".dat";
  string baseB = Form("figures/%s/data/",AnalNames[replayB].data());
  string nameB = baseB+AnalNames[replayB]+"_"+to_string(centBins[bin])+"_"+to_string(centBins[bin+1])+".dat";
  string cname = nameA+"_by_"+nameB;
  cout<<"baseA: "<<baseA<<endl;
  cout<<"nameA: "<<nameA<<endl;
  TCanvas * c = new TCanvas(cname.data(),cname.data(),650,500);

  TH1D * h = new TH1D("h","h",100,0,12.);
  h->SetMinimum(0.9);
  h->SetMaximum(1.1);
  h->SetXTitle("p_{T} (GeV/c)");
  h->SetYTitle(Form("%s / %s",AnalNames[replayA].data(),AnalNames[replayB].data()));

  h->Draw();
  TLatex * lt = new TLatex(1.,1.07,Form("%d - %d%c",centBins[bin],centBins[bin+1],'%'));
  lt->Draw();
  TGraphErrors * hA;
  TGraphErrors * hB;
  hA = GetVNPt(nameA);
  cout<<hA->GetN()<<endl;
  cout<<"return from GetVNPt "<<nameA<<"\t"<<hA<<endl;
  hB = GetVNPt(nameB);
  for(int i = 0; i<hA->GetN(); i++) {
    hA->GetY()[i]/=hB->GetY()[i];
    hA->GetEY()[i]/=hB->GetY()[i];
  }
  hA->Draw("p");
  string outname = "figures/ratios/"+AnalNames[replayA]+"_by_"+AnalNames[replayB]+"_"+to_string(centBins[bin])+"_"+to_string(centBins[bin+1])+".pdf";
  c->Print(outname.data(),"pdf");
  cout<<"return from GetVNPt "<<nameB<<"\t"<<hB<<endl;
  double ymax = 0;
}

void Ratio(){
 for(int bin = 0; bin<11; bin++) {

  GetVNCreate(N2NOEFF,N2EFF, bin);
 GetVNCreate(N42SUB3,N42EFF, bin);
 //GetVNCreate(CHI7NOEFF,CHI7EFF, bin);
    //GetVNCreate(N3SUB3,bin);
    //GetVNCreate(N4SUB3,bin);
    //GetVNCreate(N5SUB3,bin);
    //GetVNCreate(N6SUB3,bin);
    //GetVNCreate(N7SUB3,bin);
    //GetVNCreate(N42SUB3,bin);
    //GetVNCreate(N523SUB3,bin);
    //GetVNCreate(N523ASUB3,bin);
    //GetVNCreate(N62SUB3,bin);
    //GetVNCreate(N63SUB3,bin);
    //GetVNCreate(N723SUB3,bin);
    //GetVNCreate(N723ASUB3,bin);
    //GetVNCreate(D24SUB3,bin);
    //GetVNCreate(D26SUB3,bin);
    //GetVNCreate(D34SUB3,bin);
    //GetVNCreate(D2232SUB3,bin);
    //GetVNCreate(D2232ASUB3,bin);
    //GetVNCreate(D2432SUB3,bin);
    //GetVNCreate(D2432ASUB3,bin);
    //GetVNCreate(CHI4,bin);
    //GetVNCreate(CHI5,bin);
    //GetVNCreate(CHI63,bin);
    //GetVNCreate(CHI7,bin);
    //GetVNCreate(N2EFF,bin);
    //GetVNCreate(N2NOEFF,bin);
    //GetVNCreate(N723EFF,bin);
    //GetVNCreate(N723NOEFF,bin);
    //GetVNCreate(CHI7EFF,bin);
    //GetVNCreate(CHI7NOEFF,bin);
    }
}
