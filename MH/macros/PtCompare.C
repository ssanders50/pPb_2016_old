TH1D * hpt;
TH1D * hptRef;
TGraphErrors * gA;
TGraphErrors * gB;
TGraphErrors * gR;
void PtCompare(string repname="figures_-0.8_00.8/N523ASUB2/data/N523ASUB2_20_25.dat", string refrepname="figures_-0.8_00.8/N523ASUB3/data/N523ASUB3_20_25.dat" ){
  cout<<repname<<endl;
  cout<<refrepname<<endl;
  FILE * frep = fopen(repname.data(),"r");
  FILE * frefrep = fopen(refrepname.data(),"r");
  cout<<frep<<"\t"<<frefrep<<endl;
  char repbuf[80];
  char refrepbuf[80];
  double x[20];
  double r[20];
  double re[20];
  int n = 0;
  while( fgets(repbuf,80,frep)!=NULL) {
    fgets(refrepbuf,80,frefrep);
    double reppt, repvn, repvnerr;
    double refreppt, refrepvn, refrepvnerr;
    sscanf(repbuf,"%lg\t%lg\t%lg",&reppt,&repvn, &repvnerr);
    sscanf(refrepbuf,"%lg\t%lg\t%lg",&refreppt,&refrepvn, &refrepvnerr);
    if(fabs(reppt-refreppt)>0.05) {
      cout<<"pt mismatch: "<<reppt<<"\t"<<refreppt<<endl;
      break;
    } 
    x[n] = reppt;
    r[n] = repvn/refrepvn;
    re[n] = repvnerr*r[n]/repvn;
    cout<<n<<"\t"<<x[n]<<"\t"<<r[n]<<"\t"<<re[n]<<endl;
    ++n;
  }
  gR = new TGraphErrors(n,x,r,0,re);
  gR->SetMarkerStyle(20);
  gR->SetMarkerColor(kBlue);
  gR->SetLineColor(kBlue);
  gR->SetMarkerSize(1.4);
  TCanvas * c = new TCanvas("c","c",1000,800);
  //c->Divide(2,2);
  TH1D * hr = new TH1D("hr","hr",100, 0,12);
  hr->SetMinimum(0.9);
  hr->SetMaximum(1.1);
  hr->SetXTitle("p_{T} (GeV/c)");
  hr->SetYTitle("Ratio v_{5}(#Psi_{2A}, #Psi_{3B}) (2 subevent/3 subevent)");
  hr->GetYaxis()->SetTitleFont(43);
  hr->GetYaxis()->SetTitleSize(32);
  hr->GetXaxis()->SetTitleFont(43);
  hr->GetXaxis()->SetTitleSize(32);
  c->cd(2);
  hr->Draw();
  gPad->SetGrid(1,1);
  gR->Draw("p");
  c->Print("TwoThreeSubCompare.pdf");
}
