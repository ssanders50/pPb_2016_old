TFile * offsets[6];
TH2D * wqx2[6][11];
TH2D * wqy2[6][11];
TH2D * wcnt[6][11];
TH1D * ox2[6][11];
TH1D * oy2[6][11];
static const int ncentbins = 11;
static const int centBins[]={0,5,10,15,20,25,30,35,40,50,60,70,100};
void RecenteringCorrections(){
  offsets[0]=new TFile("../../VNAnalysis/offset_PbPb2015_1_262799.root","read");
  offsets[1]=new TFile("../../VNAnalysis/offset_PbPb2015_262800_263230.root","read");
  offsets[2]=new TFile("../../VNAnalysis/offset_PbPb2015_263231_263359.root","read");
  offsets[3]=new TFile("../../VNAnalysis/offset_PbPb2015_263360_263379.root","read");
  offsets[4]=new TFile("../../VNAnalysis/offset_PbPb2015_263380_263614.root","read");
  offsets[5]=new TFile("../../VNAnalysis/offset_PbPb2015_263615_263757.root","read");
  for(int i = 0; i<5; i++) {
    for(int j = 0; j<ncentbins; j++) {
      wqx2[i][j]=(TH2D *) offsets[i]->Get(Form("wqxtrk_2_%d",j));
      wqy2[i][j]=(TH2D *) offsets[i]->Get(Form("wqytrk_2_%d",j));
      ox2[i][j] = wqx2[i][j]->ProjectionX(Form("ox_%d_%d_%d",i,centBins[j],centBins[j+1]),5,8);
      oy2[i][j] = wqy2[i][j]->ProjectionX(Form("oy_%d_%d_%d",i,centBins[j],centBins[j+1]),5,8);
      ox2[i][j]->Scale(0.25);
      oy2[i][j]->Scale(0.25);
      ox2[i][j]->GetXaxis()->SetRange(1,16);
      oy2[i][j]->GetXaxis()->SetRange(1,16);

    }
  }
}
