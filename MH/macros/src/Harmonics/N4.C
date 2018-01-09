TGraphErrors * N4(int replay, int bin, double eMin, double eMax, double & ymin, double & ymax, TGraphErrors * &g,TGraphErrors * &gA, TGraphErrors * &gB, TGraphErrors * &gSpec, TGraphErrors * &gint, TGraphErrors *& gintA, TGraphErrors *& gintB){
  double vint = 0;
  double vinte = 0;
  double vintA = 0;
  double vintAe = 0;
  double vintB = 0;
  double vintBe = 0;
  int epindx = HFp4;
  TGraphErrors * gtmp;
  //
  // Start with eta distribution
  //
  for(int i = 0; i<12; i++) {
    double EtaMin = -2.4 + 0.4*i;
    double EtaMax = EtaMin+0.4;
    if(fabs(EtaMin)<0.001) EtaMin = 0.;
    if(fabs(EtaMax)<0.001) EtaMax = 0.;
    gtmp = GetVNPt(replay,bin,epindx,EtaMin,EtaMax,gtmp, gtmp, gtmp, vint, vinte, vintA, vintAe, vintB, vintBe, false);
    if(i<6) {
      gint->GetY()[i] = vintA;
      gint->GetEY()[i] = vintAe;
    } else {
      gint->GetY()[i] = vintB;
      gint->GetEY()[i] = vintBe;
    }
    gintA->GetY()[i]=vintA;
    gintA->GetEY()[i]=vintAe;
    gintB->GetY()[i]=vintB;
    gintB->GetEY()[i]=vintBe;
  }
  //
  // Now do requested calculation
  
  g = GetVNPt(replay,bin,epindx,eMin,eMax,gA, gB, gSpec, vint, vinte, vintA, vintAe, vintB, vintBe, false);
  ymin = setYmin(g,gA,gB);
  ymax = setYmax(g,gA,gB);
  outint = fopen(soutint.data(),"a+");
  fprintf(outint,"%d\t%d\t%15.10f\t%15.10f\n",cmin[bin],cmax[bin],vint,vinte);
  fclose(outint);
  cout<<"Integral: "<<eMin<<"\t"<<eMax<<"\t"<<vint<<endl;
  return g;
}
