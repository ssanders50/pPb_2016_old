enum    TrackType {typeUndefined = 0, ppReco = 1, HIReco, Pixel};
enum TrackQuality {qualityUndefined = 0, loose = 1, normal, tight, narrow, wide};
enum TrackReaction {reacUndefined = 0, pp = 1, pPb, XeXe, PbPb};

TrackType sTrackType;
TrackQuality sTrackQuality;
TrackReaction sTrackReaction;

static const  double cb[14]={0,5,10,15,20,25,30,35,40,50,60,70,80,100};
static const  double cbe[6]={0,5,10,30,50,100};
static const  double cbeXeXe[6]={0,10,30,50,70,100};
TH1I * cen=NULL;
TH1I * cene=NULL;
TH1I * ceneXeXe=NULL;
TFile * fakeFile=NULL;
TFile * effFile=NULL;

TrackType SetTracking( ){
  //
  //Use directory path to determine reaction
  //
  system("pwd > pwd.lis");
  FILE * fpwd = fopen("pwd.lis","r");
  char buf[80];
  fgets(buf,80,fpwd);
  fclose(fpwd);
  system("rm pwd.lis");
  string spwd = buf;
  cout<<spwd<<endl;
  if(spwd.find("pp")!=std::string::npos) {
    sTrackReaction=pp;
  } else if (spwd.find("pPb")!=std::string::npos) {
    sTrackReaction = pPb;
  } else if (spwd.find("XeXe")!=std::string::npos) {
    sTrackReaction = XeXe;
  } else if (spwd.find("PbPb")!=std::string::npos) {
    sTrackReaction = PbPb;
  } else {
    sTrackReaction = reacUndefined;
  }

  TFile * fin = new TFile(rootFile.data(),"read");
  TDirectory * d = (TDirectory *) fin->Get("vnanalyzer/Conditions");
  TList * l = (TList *) d->GetListOfKeys();
  int indx = 0;
  sTrackType = typeUndefined;
  sTrackQuality = qualityUndefined;
  while(l->At(indx) != l->Last()) {
    string condition = l->At(indx++)->GetName();
    if(condition == "hiGeneralAndPixelTracks") sTrackType = Pixel;
    if(condition == "hiGeneralTracks") sTrackType = HIReco;
    if(condition == "generalTracks") sTrackType = ppReco;
    if(sTrackType == Pixel) {
      cout<<"Pixel:"<<condition<<":"<<endl;
      if(condition == "dzdzerror_0002.00") sTrackQuality = normal;
      if(condition == "dzdzerror_0005.00") sTrackQuality = loose;
      if(condition == "dzdzerror_0001.50") sTrackQuality = tight;
    } else if (sTrackType == ppReco && sTrackReaction == pPb) {
      cout<<"pPb with ppReco:"<<condition<<":"<<endl;
      if(condition == "dzdzerror_0003.00") sTrackQuality = normal;
      if(condition == "dzdzerror_0005.00") sTrackQuality = loose;
      if(condition == "dzdzerror_0002.00") sTrackQuality = tight;
    }
    if(condition == "vtx_-15.0_ 3.0") sTrackQuality = narrow;
    if(condition == "vtx_- 3.0_15.0") sTrackQuality = wide;

    cout<<condition<<endl;
  }
  if(sTrackReaction==pp || sTrackReaction==pPb) {
    ncentbins = ncentbinsNOFF;
    cbins = cbinsNOFF;
    for(int i = 0; i<=ncentbins; i++) {
      centBins[i] = centBinsNOFF[i];
      centRefBins[i] = centRefBinsNOFF[i];
    }
    for(int i = 0; i<cbins;i++) {
      cmin[i] = cminNOFF[i];
      cmax[i] = cmaxNOFF[i];
    }
    ntrkbinning = true;
  } else {
    ncentbins = ncentbinsCENT;
    cbins = cbinsCENT;
    for(int i = 0; i<=ncentbins; i++) {
      centBins[i] = centBinsCENT[i];
      centRefBins[i] = centRefBinsCENT[i];
    }
    for(int i = 0; i<cbins;i++) {
      cmin[i] = cminCENT[i];
      cmax[i] = cmaxCENT[i];
    }
  }
  if(ntrkbinning) {
    rcnt = (TH1D *) fin->Get("vnanalyzer/Noff");
  } else {
    rcnt = (TH1D *) fin->Get("vnanalyzer/cent");
  }
  rcnt->SetDirectory(0);
  fin->Close();
  if(sTrackReaction==PbPb) {
    if(sTrackType == Pixel && sTrackQuality == normal) {
      fakeFile = new TFile("EffAndFake/PbPb/FakeRatesPixelPbPb_tight.root");
      effFile = new TFile("EffAndFake/PbPb/EffCorrectionsPixelPbPb_tight.root");
    }
    if(sTrackType == Pixel && sTrackQuality == tight){ 
      fakeFile = new TFile("EffAndFake/PbPb/FakeRatesPixelPbPb_tight.root");
      effFile = new TFile("EffAndFake/PbPb/EffCorrectionsPixelPbPb_tightB.root");
    }
    if(sTrackType == Pixel && sTrackQuality == loose){ 
      fakeFile = new TFile("EffAndFake/PbPb/FakeRatesPixelPbPb_loose.root");
      effFile = new TFile("EffAndFake/PbPb/EffCorrectionsPixelPbPb_loose.root");
    }
  } else if (sTrackReaction==pPb) {
    if(sTrackQuality == normal) {
      effFile = new TFile("EffAndFake/pPb_8TeV/Hijing_8TeV_dataBS.root");
    } else if (sTrackQuality == loose) {
      effFile = new TFile("EffAndFake/pPb_8TeV/Hijing_8TeV_MB_eff_v3_loose.root");
    } else if (sTrackQuality == tight) {
      effFile = new TFile("EffAndFake/pPb_8TeV/Hijing_8TeV_MB_eff_v3_tight.root");
    } else if (sTrackQuality == narrow) {
      effFile = new TFile("EffAndFake/pPb_8TeV/Hijing_8TeV_eff_v4_narrow.root");
    } else if (sTrackQuality == wide) {
      effFile = new TFile("EffAndFake/pPb_8TeV/Hijing_8TeV_v4_wide.root");
    } 
  } else if (sTrackReaction==XeXe) {
    if(sTrackQuality == normal) {
      effFile = new TFile("EffAndFake/XeXe/XeXe_eff_table_92x_cent.root");
    }
  }
  cen = new TH1I("cen","cen",13,cb);
  cene = new TH1I("cene","cene",5,cbe);
  ceneXeXe = new TH1I("ceneXeXe","ceneXeXe",5,cbeXeXe);
  cen->SetDirectory(0);
  cene->SetDirectory(0);
  ceneXeXe->SetDirectory(0);
   
  return sTrackType;
}
double FakeAndEff(int cent, double pt, double emin, double emax, double &eff) {
  eff = 1.;
  double val = 0;
  if(sTrackType == typeUndefined) return 0.;
  int ib = cen->FindBin(cent)-1;
  int ibe = cene->FindBin(cent)-1;
  if(sTrackReaction==XeXe) ibe = ceneXeXe->FindBin(cent)-1;
  if(effFile!=NULL) {
    string re = "Eff_"+to_string((int)cbe[ibe])+"_"+to_string((int)cbe[ibe+1]);
    if(sTrackReaction==pPb) re = "rTotalEff3D_0";
    if(sTrackReaction==XeXe) re = "rTotalEff3D_"+to_string((int)cbeXeXe[ibe])+"_"+to_string((int)cbeXeXe[ibe+1]);
    TH2D * he = (TH2D *) effFile->Get(re.data());
    int ptbin = he->GetYaxis()->FindBin(pt);
    if(pt>he->GetYaxis()->GetXmax()) ptbin = he->GetYaxis()->GetLast();
    int etabinmin = he->GetXaxis()->FindBin(emin+0.001);
    int etabinmax = he->GetXaxis()->FindBin(emax-0.001);
    eff = 0;
    for(int i = etabinmin; i<=etabinmax; i++) {
      eff += he->GetBinContent(i,ptbin);
    }
    eff /=(double)(etabinmax-etabinmin+1);
    he->Delete();

  } 

  if(fakeFile!=NULL) {
    string rc = "hfak_"+to_string((int)cb[ib])+"_"+to_string((int)cb[ib+1]);
    TH2D * hf = (TH2D *) fakeFile->Get(rc.data());
    if(hf==NULL) {
      cout<<"Failed to find "<<rc<<endl;
    }
    int ptbin = hf->GetYaxis()->FindBin(pt);
    int etabinmin = hf->GetXaxis()->FindBin(emin+0.001);
    int etabinmax = hf->GetXaxis()->FindBin(emax-0.001);
    val = 0;
    for(int i = etabinmin; i<=etabinmax; i++) {
      val += hf->GetBinContent(i,ptbin);
    }
    val /=(double)(etabinmax-etabinmin+1);
    hf->Delete();
  }
  return val ;
}



TH2D * ptcntEff(TH2D * ptcnt, double cent) {
  TH2D * hsEff = (TH2D *) ptcnt->Clone("ptcntEff");
  hsEff->Reset();
  hsEff->SetDirectory(0);
  for(int i = 1; i<=ptcnt->GetNbinsX(); i++) {
    for(int j = 1; j<=ptcnt->GetNbinsY(); j++) {
      double pt = ptcnt->GetXaxis()->GetBinCenter(i);
      double etmin = ptcnt->GetYaxis()->GetBinLowEdge(j);
      double etmax = etmin+ptcnt->GetYaxis()->GetBinWidth(j)-0.001;
      double eff = 1.;
      FakeAndEff(centbins->FindBin(cent), pt, etmin, etmax, eff);
      if(eff<=0 ) eff = 1;
      hsEff->SetBinContent(i,j,ptcnt->GetBinContent(i,j)/eff);
    }
  }
  return hsEff ;
}
