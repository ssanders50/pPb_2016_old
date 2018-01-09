TGraphErrors * CHI7(int replay, int bin, double eMin, double eMax, double & ymin, double & ymax, TGraphErrors * &g,TGraphErrors * &gA, TGraphErrors * &gB, TGraphErrors * &gSpec,TGraphErrors * &gint, TGraphErrors *& gintA, TGraphErrors *& gintB){
  double vint = 0;
  double vinte = 0;
  double vintA = 0;
  double vintAe = 0;
  double vintB = 0;
  double vintBe = 0;

  double vint2 = 0;
  double vinte2 = 0;
  double vintA2 = 0;
  double vintAe2 = 0;
  double vintB2 = 0;
  double vintBe2 = 0;
  bool Decor = false;
  TGraphErrors * gtmp;
  TGraphErrors * gtmp2;
  TGraphErrors * g2;
  TGraphErrors * gA2;
  TGraphErrors * gB2;
  //
  // Start with eta distribution
  //
  for(int i = 0; i<12; i++) {
    double EtaMin = -2.4 + 0.4*i;
    double EtaMax = EtaMin+0.4;
    if(fabs(EtaMin)<0.001) EtaMin = 0.;
    if(fabs(EtaMax)<0.001) EtaMax = 0.;
    if(replay==Chi5) {
      gtmp  = GetVNPt(N723SUB2,bin,-1,EtaMin,EtaMax,gtmp, gtmp, gtmp, vint, vinte, vintA, vintAe, vintB, vintBe, true);
      gtmp2 = GetVNPt(D2432SUB2,bin,-1,EtaMin,EtaMax,gtmp2, gtmp2, gtmp2, vint2, vinte2, vintA2, vintAe2, vintB2, vintBe2, true);
    } else {
      gtmp  = GetVNPt(N723ASUB2,bin,-1,EtaMin,EtaMax,gtmp, gtmp, gtmp, vint, vinte, vintA, vintAe, vintB, vintBe, true);
      gtmp2 = GetVNPt(D2432ASUB2,bin,-1,EtaMin,EtaMax,gtmp2, gtmp2, gtmp2, vint2, vinte2, vintA2, vintAe2, vintB2, vintBe2, true);
    }
    fin->Close();
    fin = new TFile(rootFile.data(),"r");
    //for five subevent calculations N42B and N42C, only the A calculation has the correct resulution calculation  
    gint->GetY()[i] = vint/vint2;
    gint->GetEY()[i] = gint->GetY()[i]*sqrt(pow(vinte/vint,2)+pow(vinte2/vint2,2));;
    gintA->GetY()[i]=vintA/vintA2;
    gintA->GetEY()[i]=gintA->GetY()[i]*sqrt(pow(vintAe/vintA,2)+pow(vintAe2/vintA2,2));
    gintB->GetY()[i]=vintB/vintB2;
    gintB->GetEY()[i]=gintB->GetY()[i]*sqrt(pow(vintBe/vintB,2)+pow(vintBe2/vintB2,2));
  }
  //
  // Now do requested calculation
  //
  if(replay==Chi7) { 
    g = GetVNPt(N723SUB2,bin,-1,eMin,eMax,gA, gB, gSpec, vint, vinte,vintA, vintAe, vintB, vintBe, true);
    g2 = GetVNPt(D2432SUB2,bin,-1,EtaMin,EtaMax,gA2, gB2, gtmp2, vint2, vinte2, vintA2, vintAe2, vintB2, vintBe2, true);
  } else {
    g = GetVNPt(N723ASUB2,bin, -1,eMin,eMax,gA, gB, gSpec, vint, vinte,vintA, vintAe, vintB, vintBe, true);
    g2 = GetVNPt(D2432ASUB2,bin,-1,EtaMin,EtaMax,gA2, gB2, gtmp2, vint2, vinte2, vintA2, vintAe2, vintB2, vintBe2, true);
  }
  for(int i = 0; i<g->GetN(); i++) {
    g->GetY()[i] = g->GetY()[i]/g2->GetY()[i];
    g->GetEY()[i]=g->GetY()[i]*sqrt(pow(g->GetEY()[i]/g->GetY()[i],2)+pow(g2->GetEY()[i]/g2->GetY()[i],2));
    gA->GetY()[i] = gA->GetY()[i]/gA2->GetY()[i];
    gA->GetEY()[i]=gA->GetY()[i]*sqrt(pow(gA->GetEY()[i]/gA->GetY()[i],2)+pow(gA2->GetEY()[i]/gA2->GetY()[i],2));
    gB->GetY()[i] = gB->GetY()[i]/gB2->GetY()[i];
    gB->GetEY()[i]=gB->GetY()[i]*sqrt(pow(gB->GetEY()[i]/gB->GetY()[i],2)+pow(gB2->GetEY()[i]/gB2->GetY()[i],2));
  }
  ymin = setYmin(g,gA,gB);
  ymax = setYmax(g,gA,gB);
  fin->Close();
  outint = fopen(soutint.data(),"a+");
  fprintf(outint,"%d\t%d\t%15.10f\t%15.10f\n",cmin[bin],cmax[bin],vint,vinte);
  fclose(outint);
  return g;
}
