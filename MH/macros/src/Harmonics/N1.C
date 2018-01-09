TGraphErrors * N1(int replay, int bin, double eMin, double eMax, double & ymin, double & ymax, TGraphErrors * &g,TGraphErrors * &gA, TGraphErrors * &gB, TGraphErrors * &gSpec, TGraphErrors * &gint, TGraphErrors *& gintA, TGraphErrors *& gintB ){
  fin = new TFile(rootFile.data(),"r");
  double vint = 0;
  double vinte = 0;
  double vintA = 0;
  double vintAe = 0;
  double vintB = 0;
  double vintBe = 0;
  TGraphErrors * g2 ;
  TGraphErrors * gA2 ;
  TGraphErrors * gB2 ;
  TGraphErrors * gSpec2 ;
  TGraphErrors * gint2 ;
  TGraphErrors * gintA2 ;
  TGraphErrors * gintB2 ;
  TGraphErrors * gtmp;
  double x[12];
  double y[12];
  double ey[12];
  gint2 = new TGraphErrors(12,x,y,0,ey);
  gintA2 = new TGraphErrors(12,x,y,0,ey);
  gintB2 = new TGraphErrors(12,x,y,0,ey);
  double vint2=0;
  double vinte2=0;
  double vintA2=0;
  double vintAe2=0;
  double vintB2=0;
  double vintBe2=0;
  int epindx = -1;
  Decor = false;
  //
  // Start with eta distribution
  //
  if(replay == N1ASUB2 || replay == N1ASUB3 || replay == N1BSUB2 || replay == N1BSUB3){
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
  } else {
    for(int i = 0; i<12; i++) {
      double EtaMin = -2.4 + 0.4*i;
      double EtaMax = EtaMin+0.4;
      if(fabs(EtaMin)<0.001) EtaMin = 0.;
      if(fabs(EtaMax)<0.001) EtaMax = 0.;
      if(replay == N1SUB2) { 
	g = GetVNPt(N1ASUB2,bin,epindx,EtaMin,EtaMax,gA, gB, gSpec, vint, vinte,vintA, vintAe, vintB, vintBe, false);
      }  else {
	g = GetVNPt(N1ASUB3,bin,epindx,EtaMin,EtaMax,gA, gB, gSpec, vint, vinte, vintA, vintAe, vintB, vintBe,false);
      }
      fin->Close();
      fin = new TFile(rootFile.data(),"r");
      if(replay == N1SUB2) { 
	g2 = GetVNPt(N1BSUB2,bin,epindx,EtaMin,EtaMax,gA2, gB2, gSpec2, vint2, vinte2,vintA2, vintAe2, vintB2, vintBe2, false);
      }  else {
	g2 = GetVNPt(N1BSUB3,bin,epindx,EtaMin,EtaMax,gA2, gB2, gSpec2, vint2, vinte2,vintA2, vintAe2, vintB2, vintBe2, false);
      }
      fin->Close();
 	gint->GetY()[i] = (vintA+vintA2+vintB+vintB2)/4.;
	gint->GetEY()[i] = (vintAe+vintAe2+vintBe+vintBe2)/4.;
      gintA->GetY()[i]=(vintA+vintA2)/2.;
      gintA->GetEY()[i]=(vintAe+vintAe2)/2.;
      gintB->GetY()[i]=(vintB+vintB2)/2.;
      gintB->GetEY()[i]=(vintBe+vintBe2)/2.;
    }
  }
  //
  // Now do requested calculation
  //
  if(replay == N1ASUB2 || replay == N1ASUB3 || replay == N1BSUB2 || replay == N1BSUB3){
    g = GetVNPt(replay,bin,epindx,eMin,eMax,gA, gB, gSpec, vint, vinte,vintA, vintAe, vintB, vintBe, false);
    ymin = setYmin(g,gA,gB);
    ymax = setYmax(g,gA,gB);
    fin->Close();
    outint = fopen(soutint.data(),"a+");
    fprintf(outint,"%d\t%d\t%15.10f\t%15.10f\n",cmin[bin],cmax[bin],vint,vinte);
    fclose(outint);
    return g;
  }
  if(replay == N1SUB2) { 
    g = GetVNPt(N1ASUB2,bin,epindx,eMin,eMax,gA, gB, gSpec, vint, vinte, vintA, vintAe, vintB, vintBe,false);
  }  else {
    g = GetVNPt(N1ASUB3,bin,epindx,eMin,eMax,gA, gB, gSpec, vint, vinte,vintA, vintAe, vintB, vintBe, false);
  }
  fin->Close();
  fin = new TFile(rootFile.data(),"r");
  if(replay == N1SUB2) { 
    g2 = GetVNPt(N1BSUB2,bin,epindx,eMin,eMax,gA2, gB2, gSpec2, vint2, vinte2, vintA2, vintAe2, vintB2, vintBe2,false);
  }  else {
    g2 = GetVNPt(N1BSUB3,bin,epindx,eMin,eMax,gA2, gB2, gSpec2, vint2, vinte2,vintA2, vintAe2, vintB2, vintBe2, false);
  }
  fin->Close();
  for(int i = 0; i<g->GetN(); i++) {
    g->GetY()[i] = (g->GetY()[i]+g2->GetY()[i])/2.;
    gA->GetY()[i] = (gA->GetY()[i]+gA2->GetY()[i])/2.;
    gB->GetY()[i] = (gB->GetY()[i]+gB2->GetY()[i])/2.;
  }
  ymin = setYmin(g,gA,gB);
  ymax = setYmax(g,gA,gB);
  return g;
}
