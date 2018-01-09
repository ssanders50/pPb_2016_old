void drawSpec(){
  // string c2name = "c2_"+AnalNames[replay]+"_"+to_string(cmin[bin])+"_"+to_string(cmax[bin]);
  // TCanvas * c2 = new TCanvas(c2name.data(),c2name.data(),700,500);
  // hspec = new TH1D("hspec","hspec",600,0,12);
  // hspec->SetDirectory(0);
  // hspec->SetMaximum(100*pow(10.,(double)((int) TMath::Log10(ymaxspec))));
  // hspec->SetMinimum(0.00001);
  // c2->cd();
  // gPad->SetLogy();
  // hspec->Draw();
  // hspec->SetXTitle("p_{T} (GeV/c)");
  // hspec->SetYTitle("1/(N_{ev}) d^{2}N/dp_{T}d#eta");
  // nwspec->Draw("p");
  // double ym = hspec->GetMaximum();
  // TLatex * t4 = new TLatex(8,0.1*ym,AnalNames[replay].data());
  // t4->SetTextFont(43);
  // t4->SetTextSize(28);
  // t4->Draw();
  // TLatex * t6 = new TLatex(8,0.02*ym,Form("%d #leq N_{trk}^{off} < %d",cmin[bin],cmax[bin]));
  // t6->SetTextFont(43);
  // t6->SetTextSize(22);
  // t6->Draw();
  // TLatex * t7 = new TLatex(7.8,0.004*ym,Form("%4.1f < #eta < %4.1f",EtaMin,EtaMax));
  // t7->SetTextFont(43);
  // t7->SetTextSize(22);
  // t7->Draw();
  // if(isTight) {
  //   c2->Print(Form("%s/%s/%s.pdf",FigDir.data(),AnalNames[replay].data(),c2name.data()),"pdf");
  //   sspec = FigDir+"/"+AnalNames[replay]+"/data/spec_"+to_string(cmin[bin])+"_"+to_string(cmax[bin])+".dat";
  // } else {
  //   c2->Print(Form("%s/%s/%s.pdf",FigDir.data(),AnalNames[replay].data(),c2name.data()),"pdf");
  //   sspec = FigDir+"/"+AnalNames[replay]+"/data/spec_"+to_string(cmin[bin])+"_"+to_string(cmax[bin])+".dat";
  // }
}
