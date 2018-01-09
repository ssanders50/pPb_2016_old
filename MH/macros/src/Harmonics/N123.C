TGraphErrors * N123(int replay, int bin, double eMin, double eMax, double & ymin, double & ymax, TGraphErrors * &g,TGraphErrors * &gA, TGraphErrors * &gB, TGraphErrors * &gSpec, TGraphErrors * &gint, TGraphErrors *& gintA, TGraphErrors *& gintB ){
  fin = new TFile(rootFile.data(),"r");
  double vint = 0;
  double vinte = 0;
  int epindx = -1;
  if(replay == N123ASUB2 || replay == N123ASUB3 || replay == N123BSUB2 || replay == N123BSUB3){
    g = GetVNPt(replay,bin,epindx,eMin,eMax,gA, gB, gSpec, gint, gintA,gintB,vint, vinte, false);
    ymin = setYmin(g,gA,gB);
    ymax = setYmax(g,gA,gB);
    fin->Close();
    outint = fopen(soutint.data(),"a+");
    fprintf(outint,"%d\t%d\t%15.10f\t%15.10f\n",cmin[bin],cmax[bin],vint,vinte);
    fclose(outint);
    return g;
  }
  if(replay == N123SUB2) { 
    g = GetVNPt(N123ASUB2,bin,epindx,eMin,eMax,gA, gB, gSpec, gint, gintA,gintB,vint, vinte, false);
  }  else {
    g = GetVNPt(N123ASUB3,bin,epindx,eMin,eMax,gA, gB, gSpec, gint, gintA,gintB,vint, vinte, false);
  }
  fin->Close();
  fin = new TFile(rootFile.data(),"r");
  TGraphErrors * g2 ;
  TGraphErrors * gA2 ;
  TGraphErrors * gB2 ;
  TGraphErrors * gSpec2 ;
  TGraphErrors * gint2 ;
  TGraphErrors * gintA2 ;
  TGraphErrors * gintB2 ;
  double x[12];
  double y[12];
  double ey[12];
  gint2 = new TGraphErrors(12,x,y,0,ey);
  gintA2 = new TGraphErrors(12,x,y,0,ey);
  gintB2 = new TGraphErrors(12,x,y,0,ey);
  double vint2=0;
  double vinte2=0;
  if(replay == N123SUB2) { 
    g2 = GetVNPt(N123BSUB2,bin,epindx,eMin,eMax,gA2, gB2, gSpec2, gint2, gintA2,gintB2,vint2, vinte2, false);
  }  else {
    g2 = GetVNPt(N123BSUB3,bin,epindx,eMin,eMax,gA2, gB2, gSpec2, gint2, gintA2,gintB2,vint2, vinte2, false);
  }
  fin->Close();
  //Only A has the correct symmetry...
  for(int i = 0; i<g->GetN(); i++) {
    gA->GetY()[i] = (gA->GetY()[i]+gA2->GetY()[i])/2.;
    gB->GetY()[i] = gA->GetY()[i];
    g->GetY()[i] = gA->GetY()[i];
    gA->GetEY()[i] = (gA->GetEY()[i]+gA2->GetEY()[i])/2.;
    gB->GetEY()[i] = gA->GetEY()[i];
    g->GetEY()[i] = gA->GetEY()[i];
  }
  ymin = setYmin(g,gA,gB);
  ymax = setYmax(g,gA,gB);
  double etab = ((EtaMin+EtaMax)/2.+2.2)/0.4;
  if(fabs(fmod(etab,1.))<0.05) {
    int ie = (int)(etab+0.05);
    
    gintA->GetY()[ie] = (gintA->GetY()[ie]+gintA2->GetY()[ie])/2.;
    gintA->GetEY()[ie] = (gintA->GetEY()[ie]+gintA2->GetEY()[ie])/2.;
    
    gintB->GetY()[ie] = gintA->GetY()[ie];
    gintB->GetEY()[ie] = gintA->GetEY()[ie];

    gint->GetY()[ie] = gintA->GetY()[ie];
    gint->GetEY()[ie] = gintA->GetEY()[ie];
}
  return g;
}
