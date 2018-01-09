TGraphErrors * N42(int replay, int bin, double eMin, double eMax, double & ymin, double & ymax, TGraphErrors * &g,TGraphErrors * &gA, TGraphErrors * &gB, TGraphErrors * &gSpec,TGraphErrors * &gint, TGraphErrors *& gintA, TGraphErrors *& gintB){
  double vint = 0;
  double vinte = 0;
  double vintA = 0;
  double vintAe = 0;
  double vintB = 0;
  double vintBe = 0;
  bool Decor = false;
  TGraphErrors * gtmp;
  //
  // Start with eta distribution
  //
  cout<<"Enter N42 ####"<<endl;
  for(int i = 0; i<12; i++) {
    double EtaMin = -2.4 + 0.4*i;
    double EtaMax = EtaMin+0.4;
    if(fabs(EtaMin)<0.001) EtaMin = 0.;
    if(fabs(EtaMax)<0.001) EtaMax = 0.;
    gtmp = GetVNPt(replay,bin,-1,EtaMin,EtaMax,gtmp, gtmp, gtmp, vint, vinte, vintA, vintAe, vintB, vintBe, false);
    fin->Close();
    fin = new TFile(rootFile.data(),"r");
    //for five subevent calculations N42B and N42C, only the A calculation has the correct resulution calculation  
    gint->GetY()[i] = vint;
    gint->GetEY()[i] = vinte;
    gintA->GetY()[i]=vintA;
    gintA->GetEY()[i]=vintAe;
    gintB->GetY()[i]=vintB;
    gintB->GetEY()[i]=vintBe;
    if( replay==N42BSUB3 || replay==N42CSUB3) {
      for(int i = 0; i<gint->GetN(); i++) {
    	gint->GetY()[i] = gintA->GetY()[i];
    	gintB->GetY()[i] = gintA->GetY()[i];
    	gint->GetEY()[i] = gintA->GetEY()[i];
    	gintB->GetEY()[i] = gintA->GetEY()[i];
      }
    }
  }
  //
  // Now do requested calculation
  //
  g = GetVNPt(replay,bin,-1,eMin,eMax,gA, gB, gSpec, vint, vinte,vintA, vintAe, vintB, vintBe, false);
  //for five subevent calculations N42B and N42C, only the A calculation has the correct resulution calculation  
  if( replay==N42BSUB3 || replay==N42CSUB3) {
    for(int i = 0; i< g->GetN(); i++) {
      g->GetY()[i] = gA->GetY()[i];
      gB->GetY()[i] = gA->GetY()[i];
      g->GetEY()[i] = gA->GetEY()[i];
      gB->GetEY()[i] = gA->GetEY()[i];
    }
  }
  ymin = setYmin(g,gA,gB);
  ymax = setYmax(g,gA,gB);
  fin->Close();
  outint = fopen(soutint.data(),"a+");
  fprintf(outint,"%d\t%d\t%15.10f\t%15.10f\n",cmin[bin],cmax[bin],vint,vinte);
  fclose(outint);
  return g;
}
