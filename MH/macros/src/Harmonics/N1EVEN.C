TGraphErrors * N1EVEN(int replay, int bin, double eMin, double eMax, double & ymin, double & ymax, TGraphErrors * &g,TGraphErrors * &gA, TGraphErrors * &gB, TGraphErrors * &gSpec, TGraphErrors * &gint, TGraphErrors *& gintA, TGraphErrors *& gintB ){
  fin = new TFile(rootFile.data(),"r");
  double vint = 0;
  double vinte = 0;
  double vintA = 0;
  double vintAe = 0;
  double vintB = 0;
  double vintBe = 0;
  int epindx = -1;
  TGraphErrors * gtmp;
  Decor = false;
  //
  // Start with eta distribution
  //
  if(replay==N1EVENSUB2 || replay==N1EVENSUB3) {
    int ioff = 0;
    for(int i = 0+ioff; i<12-ioff; i++) {
      double EtaMin = -2.4 + 0.4*i;
      double EtaMax = EtaMin+0.4;
      if(fabs(EtaMin)<0.001) EtaMin = 0.;
      if(fabs(EtaMax)<0.001) EtaMax = 0.;
      if(i<6) {
	int rep = N1MCp18SUB2;
	if(replay==N1EVENSUB3) rep = N1MCp18SUB3;
	gtmp = GetVNPt(rep,bin,epindx,EtaMin,EtaMax,gtmp, gtmp, gtmp, vint, vinte, vintA, vintAe, vintB, vintBe, false);
	gint->GetY()[i-ioff] = vintA;
	gint->GetEY()[i-ioff] = vintAe;
      } else {
	int rep = N1MCm18SUB2;
	if(replay==N1EVENSUB3) rep = N1MCm18SUB3;
	gtmp = GetVNPt(rep,bin,epindx,EtaMin,EtaMax,gtmp, gtmp, gtmp, vint, vinte, vintA, vintAe, vintB, vintBe, false);
	gint->GetY()[i-ioff] = vintA;
	gint->GetEY()[i-ioff] = vintAe;
      }
      gint->GetX()[i-ioff] = gint->GetX()[i];
      gintA->GetX()[i-ioff] = gintA->GetX()[i];
      gintB->GetX()[i-ioff] = gintB->GetX()[i];
      gintA->GetY()[i-ioff]=vintA;
      gintA->GetEY()[i-ioff]=vintAe;
      gintB->GetY()[i-ioff]=vintB;
      gintB->GetEY()[i-ioff]=vintBe;
    }
    gint->Set(12-2*ioff);
    gintA->Set(12-2*ioff);
    gintB->Set(12-2*ioff);
  } else {
    for(int i = 0; i<12; i++) {
      double EtaMin = -2.4 + 0.4*i;
      double EtaMax = EtaMin+0.4;
      if(fabs(EtaMin)<0.001) EtaMin = 0.;
      if(fabs(EtaMax)<0.001) EtaMax = 0.;
      gtmp = GetVNPt(replay,bin,epindx,EtaMin,EtaMax,gtmp, gtmp, gtmp, vint, vinte, vintA, vintAe, vintB, vintBe, false);
      gint->GetY()[i] = vint;
      gint->GetEY()[i] = vinte;
      gintA->GetY()[i]=vintA;
      gintA->GetEY()[i]=vintAe;
      gintB->GetY()[i]=vintB;
      gintB->GetEY()[i]=vintBe;
    }
  }
  //
  // Now do requested calculation
  //
  if(replay != N1EVENSUB3  && replay !=N1EVENSUB2){
    g = GetVNPt(replay,bin,epindx,eMin,eMax,gA, gB, gSpec, vint, vinte, vintA,vintAe,vintB,vintBe,false);
    for(int i = 0; i<g->GetN(); i++) {
      g->GetY()[i] = gA->GetY()[i];
      g->GetEY()[i] = gA->GetEY()[i];
    }
    SetToA = true;
    ymin = setYmin(g,gA,gB);
    ymax = setYmax(g,gA,gB);
    fin->Close();
    outint = fopen(soutint.data(),"a+");
    fprintf(outint,"%d\t%d\t%15.10f\t%15.10f\n",cmin[bin],cmax[bin],vint,vinte);
    fclose(outint);
  } else {
    int settype = 0;
    if(eMin>=-0.8) {
      settype = N1MCm22SUB2;
    } else if (eMax<=0.8) {
      settype = N1MCp22SUB2;
    }
    double etab = ((EtaMin+EtaMax)/2.+2.2)/0.4;
    double x[10];
    double y[10];
    double ey[10];
    
    if(fabs(fmod(etab,1.))<0.05) {
      int ie = (int)(etab+0.05);
      if (ie==0) {
       settype = N1MCp18SUB2;
      } else if (ie>=1&&ie<=10) {
	if(etab<0) {
	  settype = N1MCp22SUB2;
	} else {
	  settype = N1MCm22SUB2;
	}
      } else {
	settype = N1MCm18SUB2;
      }
  }
    g = GetVNPt(settype,bin,epindx,eMin,eMax,gA, gB, gSpec, vint, vinte,vintA,vintAe,vintB,vintBe, false);
    
    ymin = setYmin(g,gA,gB);
    ymax = setYmax(g,gA,gB);
  }
  

  return g;
}
