TGraphErrors * GetVNPt(int replay, int bin, int epindx,  double etamin, double etamax, TGraphErrors * &gA, TGraphErrors * &gB, TGraphErrors * &gspec,  double & vint, double & vinte,    double & vintA, double & vintAe, double & vintB, double & vintBe, bool nonorm=false){

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
  TH2D * res2D = 0;
  TH2D * resw2D = 0;
  vint = 0;
  vinte = 0;
  vintA = 0;
  vintAe = 0;
  vintB = 0;
  vintBe = 0;
  fin = new TFile(rootFile.data(),"r");
  int jmin = centRef->FindBin(cmin[bin]+0.001)-1;
  int jmax = centRef->FindBin(cmax[bin]-0.001)-1;
  ANAL = replay;
  string strip = ANALS[replay][0];
  bool sub2 = false;
  TString subtest = ANALS[replay][0];
  if(subtest.Contains("SUB2")) sub2 = true;
  if(sub2) {
    strip = strip.substr(0,strip.find("SUB2"));
    Decor = false;
  } else {
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
    TDirectory * d = (TDirectory *) fin->Get("vnanalyzer/Resolutions");
    TList * l = (TList *) d->GetListOfKeys();
    string rdir = l->At(j)->GetName();
    if(cmin[j]==1) crange = "0_"+to_string(cmax[j]);
    if(j==jmin) {
      if(epindx>=0) {
	res2D  = (TH2D *) fin->Get(Form("vnanalyzer/Resolutions/%s/res%d",rdir.data(),EPOrder[epindx]));
	res2D->SetDirectory(0);
	resw2D = (TH2D *) fin->Get(Form("vnanalyzer/Resolutions/%s/resw%d",rdir.data(),EPOrder[epindx]));
	resw2D->SetDirectory(0);
      }
      ptav = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/ptav",crange.data()));
      ptav->SetDirectory(0);
      ptcnt = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/ptcnt",crange.data()));
      ptcnt->SetDirectory(0);
      badcnt = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/badcnt",crange.data()));
      badcnt->SetDirectory(0);
      centbins = (TH1D * ) fin->Get("vnanalyzer/centbins");
      centbins->SetDirectory(0);
      centcnt+=centbins->GetBinContent(j+1);
      qA = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qA",crange.data(),strip.data()));
      qA->SetDirectory(0);
      qB = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qB",crange.data(),strip.data()));
      qB->SetDirectory(0);
      wnA = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/wnA",crange.data(),strip.data()));
      wnA->SetDirectory(0);
      wnB = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/wnB",crange.data(),strip.data()));
      wnB->SetDirectory(0);
      qBA += ((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qBA",crange.data(),strip.data())))->GetBinContent(1);
      qBAcnt+=((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qBAcnt",crange.data(),strip.data())))->GetBinContent(1);
      qCA += ((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qCA",crange.data(),strip.data())))->GetBinContent(1);
      qCAcnt+=((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qCAcnt",crange.data(),strip.data())))->GetBinContent(1); 
      qCB += ((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qCB",crange.data(),strip.data())))->GetBinContent(1);
      qCBcnt+=((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/qCBcnt",crange.data(),strip.data())))->GetBinContent(1);
      for(int i = 0; i<10; i++) {
	qAe[i] = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qA_%d",crange.data(),strip.data(),i+1));
	qAe[i]->SetDirectory(0);
	qBe[i] = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qB_%d",crange.data(),strip.data(),i+1));
	qBe[i]->SetDirectory(0);
	wnAe[i] = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/wnA_%d",crange.data(),strip.data(),i+1));
	wnAe[i]->SetDirectory(0);
	wnBe[i] = (TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/wnB_%d",crange.data(),strip.data(),i+1));
	wnBe[i]->SetDirectory(0);
	qBAe[i] += ((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qBA_%d",crange.data(),strip.data(),i+1)))->GetBinContent(1);
	qBAecnt[i]+=((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qBAcnt_%d",crange.data(),strip.data(),i+1)))->GetBinContent(1);
	qCAe[i] += ((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qCA_%d",crange.data(),strip.data(),i+1)))->GetBinContent(1);
	qCAecnt[i]+=((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qCAcnt_%d",crange.data(),strip.data(),i+1)))->GetBinContent(1);
	qCBe[i] += ((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qCB_%d",crange.data(),strip.data(),i+1)))->GetBinContent(1);
	qCBecnt[i]+=((TH2D *) fin->Get(Form("vnanalyzer/Harmonics/%s/%s/SubEvents/qCBcnt_%d",crange.data(),strip.data(),i+1)))->GetBinContent(1);
      }
    } else {
      centcnt+=centbins->GetBinContent(j);
      if(epindx>=0) {
	res2D->Add( (TH2D *) fin->Get(Form("vnanalyzer/Resolutions/%s/res%d",rdir.data(),EPOrder[epindx])));
	resw2D->Add( (TH2D *) fin->Get(Form("vnanalyzer/Resolutions/%s/resw%d",rdir.data(),EPOrder[epindx])));
      }
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
  fin->Close();
  if(centcnt<100) {
    cout<<"centcnt: "<<centcnt<<endl;
    return NULL;
  }

  ptav->Divide(ptcnt);
  
  int ietamin1=0;

  int ietamax1=0;
  int ietamin2=0;
  int ietamax2=0;
  double sign = 1.;
  //if(replay==N112ASUB2 || replay==N112ASUB3) sign=-1.;
  if(etamin*etamax<0) {
    ietamin1 = qA->GetYaxis()->FindBin(etamin+0.001);

    ietamax1 = qA->GetYaxis()->FindBin(-0.001);
    ietamin2 = qA->GetYaxis()->FindBin(0.001);
    ietamax2 = qA->GetYaxis()->FindBin(etamax-0.001);
    qA1 = (TH1D *) qA->ProjectionX("qA1",ietamin1,ietamax1);
    qB1 = (TH1D *) qB->ProjectionX("qB1",ietamin2,ietamax2);
    wA1 = (TH1D *) wnA->ProjectionX("wA1",ietamin1,ietamax1);
    wB1 = (TH1D *) wnB->ProjectionX("wB1",ietamin2,ietamax2);
  } else {

    ietamin1 = qA->GetYaxis()->FindBin(etamin+0.001);
    ietamax1 = qA->GetYaxis()->FindBin(etamax-0.001);
    qA1 = (TH1D *) qA->ProjectionX("qA1",ietamin1,ietamax1);
    qB1 = (TH1D *) qB->ProjectionX("qB1",ietamin1,ietamax1);
    wA1 = (TH1D *) wnA->ProjectionX("wA1",ietamin1,ietamax1);
    wB1 = (TH1D *) wnB->ProjectionX("wB1",ietamin1,ietamax1);
  }

  qBA/=qBAcnt;
  qCA/=qCAcnt;
  qCB/=qCBcnt;
  for(int i = 0; i<10; i++) {
    qBAe[i]/=qBAecnt[i];
    qCAe[i]/=qCAecnt[i];
    qCBe[i]/=qCBecnt[i];    
  }
  if(epindx>=0 && ietamin1==ietamax1 && Decor) {
    res2D->Divide(resw2D);
    int iorder = EPOrder[epindx];
    int epmin = 0;
    int epmax = 0;
    if(iorder == 1 ) {
      epmin = HFm1;
      epmax = HFp1f;
    }else if (iorder == 2) {
      epmin = HFm2;
      epmax = HFp2f;
    }else if (iorder == 3) {
      epmin = HFm3;
      epmax = HFp3f;
    }else if (iorder == 4) {
      epmin = HFm4;
      epmax = HFp4f;
    }else if (iorder == 5) {
      epmin = HFm5;
      epmax = trackp522;
    }else if (iorder == 6) {
      epmin = HFm6;
      epmax = trackp622;
    }else if (iorder == 7) {
      epmin = HFm7;
      epmax = trackp722;
    }
    int epA = epindx;
    int epB = RCMate1[epindx];
    int epC = RCMate2[epindx];
    epC = min(epA,epB)+4+ietamin1;
    double rAB = res2D->GetBinContent(min(epA-epmin+1,epB-epmin+1),max(epA-epmin+1,epB-epmin+1));
    double rAC = res2D->GetBinContent(min(epA-epmin+1,epC-epmin+1),max(epA-epmin+1,epC-epmin+1));
    double rBC = res2D->GetBinContent(min(epB-epmin+1,epC-epmin+1),max(epB-epmin+1,epC-epmin+1));
    qBA = rAB;
    qCA = rAC;
    qCB = rBC;
    for(int i = 0; i<10; i++) {
      qBAe[i] = qBA;
      qCAe[i] = qCA;
      qCBe[i] = qCB;    
    }
    
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
  TH2D * hsEff;
  if(etamin*etamax<0) {
    xpt = (TH1D *) ptav->ProjectionX("xpt",ietamin1,ietamax2);
    double c = (cmin[bin]+cmax[bin])/2.;
    hsEff = ptcntEff(ptcnt,c);
    sp = (TH1D *) hsEff->ProjectionX("sp",ietamin1,ietamax2);
    double ebinsA = ietamax1-ietamin1+1 ;
    double ebinsB = ietamax2-ietamin2+1;
    xpt->Scale(1./(ebinsA+ebinsB));
    sp->Scale(1./(ebinsA+ebinsB));
    vnA = (TH1D *) qA->ProjectionX("vnA",ietamin1,ietamax1);
    vnA->SetDirectory(0);
    vnB = (TH1D *) qB->ProjectionX("vnB",ietamin2,ietamax2);
    vnB->SetDirectory(0);
    vn = (TH1D *) vnA->Clone("vn");
    vn->SetDirectory(0);
    vn->Add(vnB,sign);
    vn->Scale(sign);
    vnA->Scale(1./ebinsA);
    vnB->Scale(1./ebinsB);
    vn->Scale(1./(ebinsA+ebinsB));
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
      vnAe = (TH1D *) qAe[i]->ProjectionX(Form("vnA%d",i),ietamin1,ietamax1);
      vnBe = (TH1D *) qBe[i]->ProjectionX(Form("vnB%d",i),ietamin2,ietamax2);
      vne = (TH1D *) vnAe->Clone(Form("vn%d",i));
      vne->Add(vnBe,sign);
      vn->Scale(sign);
      vnAe->Scale(1./ebinsA);
      vnBe->Scale(1./ebinsB);
      vne->Scale(1./(ebinsA+ebinsB));
      
      for(int j = 0; j<vne->GetNbinsX(); j++) {
	vnm[j]+= vne->GetBinContent(j+1);
	vnAm[j]+= vnAe->GetBinContent(j+1);
	vnBm[j]+= vnBe->GetBinContent(j+1);
	vn2[j] += pow(vne->GetBinContent(j+1),2);
	vnA2[j]+= pow(vnAe->GetBinContent(j+1),2);
	vnB2[j]+= pow(vnBe->GetBinContent(j+1),2);
      }
      vnAe->Delete();
      vnBe->Delete();
      vne->Delete();

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
    hsEff = ptcntEff(ptcnt,c);
    //hsEff->Draw("colz");

    sp = (TH1D *) hsEff->ProjectionX("sp",ietamin1,ietamax1);
    double ebinsA = ietamax1-ietamin1+1 ;
    xpt->Scale(1./ebinsA);
    sp->Scale(1./ebinsA);
    vnA = (TH1D *) qA->ProjectionX("vnA",ietamin1,ietamax1);
    vnA->SetDirectory(0);
    vnB = (TH1D *) qB->ProjectionX("vnB",ietamin1,ietamax1);
    vnB->SetDirectory(0);
    vn = (TH1D *) vnA->Clone("vn");
    vn->SetDirectory(0);
    vn->Add(vnB,sign);
    vn->Scale(sign);
    vnA->Scale(1./ebinsA);

    vnB->Scale(1./ebinsA);
    vn->Scale(1./(ebinsA+ebinsA));
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
	vnAe = (TH1D *) qAe[i]->ProjectionX(Form("vnA%d",i),ietamin1,ietamax1);
	vnBe = (TH1D *) qBe[i]->ProjectionX(Form("vnB%d",i),ietamin1,ietamax1);
	vne = (TH1D *) vnAe->Clone(Form("vn%d",i));
	vne->Add(vnBe,sign);
	vne->Scale(sign);

	vnAe->Scale(1./ebinsA);
	vnBe->Scale(1./ebinsA);
	vne->Scale(1./(ebinsA+ebinsA));
      
      for(int j = 0; j<vne->GetNbinsX(); j++) {
	vnm[j]+= vne->GetBinContent(j+1);
	vnAm[j]+= vnAe->GetBinContent(j+1);
	vnBm[j]+= vnBe->GetBinContent(j+1);
	vn2[j] += pow(vne->GetBinContent(j+1),2);
	vnA2[j]+= pow(vnAe->GetBinContent(j+1),2);
	vnB2[j]+= pow(vnBe->GetBinContent(j+1),2);
      }
      vnAe->Delete();
      vnBe->Delete();
      vne->Delete();
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
  double wvnA = 0;
  double wvnAe = 0;
  double wA = 0;
  double wvnB = 0;
  double wvnBe = 0;
  double wB = 0;
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
	
	double fake = FakeAndEff(cent,pt,etamin,etamax,eff);
	wvn += y[npt]*yld->GetBinContent(i)/eff;
	wvne+= ey[npt]*yld->GetBinContent(i)/eff;
	w   += yld->GetBinContent(i)/eff;;
	wvnA += yA[npt]*yld->GetBinContent(i)/eff;
	wvnAe+= eyA[npt]*yld->GetBinContent(i)/eff;
	wA   += yld->GetBinContent(i)/eff;;
	wvnB += yB[npt]*yld->GetBinContent(i)/eff;
	wvnBe+= eyB[npt]*yld->GetBinContent(i)/eff;
	wB   += yld->GetBinContent(i)/eff;;
      }
      ++npt;
    }
  }
  yld->Delete();
  vint = wvn/w;
  vinte = wvne/w;
  vintA = wvnA/wA;
  vintAe = wvnAe/wA;
  vintB = wvnB/wB;
  vintBe = wvnBe/wB;
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

  hsEff->Delete();
  xpt->Delete();
  sp->Delete();
  vn->Delete();
  vnA->Delete();
  vnB->Delete();
  qA1->Delete();
  qB1->Delete();
  wA1->Delete();
  wB1->Delete();
  ptav->Delete();
  ptcnt->Delete();
  badcnt->Delete();
  res2D->Delete();
  resw2D->Delete();
  qA->Delete();
  qB->Delete();
  wnA->Delete();
  wnB->Delete();
  for(int i = 0; i<10; i++) {
    qAe[i]->Delete();
    qBe[i]->Delete();
    wnAe[i]->Delete();
    wnBe[i]->Delete();
    qAe1[i]->Delete();
    qBe1[i]->Delete();
  }

  return g;

}
