// -*- C++ -*-
//
// Package:    VNAnalyzer
// Class:      VNAnalyzer
// 


// system include files
#include <memory>
#include <algorithm>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Math/Vector3D.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/HeavyIonEvent/interface/Centrality.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h>
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidate.h>

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "CondFormats/DataRecord/interface/HeavyIonRPRcd.h"
#include "CondFormats/DataRecord/interface/HeavyIonRcd.h"
#include "CondFormats/HIObjects/interface/CentralityTable.h"
#include "CondFormats/HIObjects/interface/RPFlatParams.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"
#include "HeavyIonsAnalysis/VNAnalysis/interface/TrackEfficiency.h"

#include "TROOT.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TTree.h"
#include "TH1I.h"
#include "TF1.h"
#include "TMath.h"
#include "TRandom.h"
#include <time.h>
#include <cstdlib>
	
#include <vector>
using std::vector;
using std::rand;
using namespace std;
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneFlatten.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/HiEvtPlaneList.h"
#include "RecoHI/HiEvtPlaneAlgos/interface/LoadEPDB.h"
using namespace hi;
using namespace edm;

static const int ntrkbins = 25;
static const  double trkBins[]={0, 10, 20, 30, 40, 50, 60, 70, 80, 100, 120, 135, 150, 160, 185, 210, 230, 250, 270, 300, 330, 350, 370, 390, 420, 500};
static const int nptbins = 28;
static const float ptbins[]={0.3, 0.4, 0.5,  0.6,  0.8,  1.0,  1.25,  1.50,  2.0,
			      2.5,  3.0,  3.5,  4.0,  5.0,  6.0,  7.0, 8.0, 10., 12.0, 14.0, 16.0,  20.0, 26.0, 35.0, 45.0, 60.0, 80.0, 100., 200.};

static const int MaxTracks = 50;

static const int netabinsDefault = 12;
static const float etabinsDefault[]={-2.4, -2.0, -1.6, -1.2, -0.8, -0.4, 0, 0.4, 0.8, 1.2, 1.6, 2.0, 2.4};


//
// class declaration
//

class VNAnalyzer : public edm::EDAnalyzer {
public:
  explicit VNAnalyzer(const edm::ParameterSet&);
  ~VNAnalyzer();
      
private:
  int NtrkToBin(int ntrk){
    for(int i = 0; i<=ntrkbins; i++) {
      if(ntrk < trkBins[i]) return i;
    }
    return ntrkbins;
  }
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  bool CaloMatch(const reco::Track & track, const edm::Event & iEvent, unsigned int idx);
  // ----------member data ---------------------------
  int eporder_;


  std::string centralityVariable_;
  std::string centralityLabel_;
  std::string centralityMC_;

  edm::InputTag centralityBinTag_;
  edm::EDGetTokenT<int> centralityBinToken;
  edm::Handle<int> cbin_;
  edm::EDGetTokenT<int> tag_;

  edm::InputTag centralityTag_;
  edm::EDGetTokenT<reco::Centrality> centralityToken;
  edm::Handle<reco::Centrality> centrality_;


  edm::InputTag pfTag_;
  edm::EDGetTokenT<reco::PFCandidateCollection> pfToken_;


  edm::InputTag vertexTag_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> vertexToken;
  edm::Handle<std::vector<reco::Vertex>> vertex_;

  edm::InputTag trackTag_;
  edm::EDGetTokenT<reco::TrackCollection> trackToken;
  edm::Handle<reco::TrackCollection> trackCollection_;

  edm::InputTag inputPlanesTag_;
  edm::EDGetTokenT<reco::EvtPlaneCollection> inputPlanesToken;
  edm::Handle<reco::EvtPlaneCollection> inputPlanes_;

  edm::Service<TFileService> fs;
  TFile * frecenter;
  string offsetFileName;

  double caloCentRef_;
  double caloCentRefWidth_;
  int caloCentRefMinBin_;
  int caloCentRefMaxBin_;

  double nCentBins_;
  bool useNtrk_;

  int vs_sell;   // vertex collection size
  float vzr_sell;
  float vzErr_sell;
  TH1D * hcent;
  TH1D * hcentbins;
  TH1D * hptNtrk;
  TH1D * hptNtrkGood;
  TH1I * hNtrkRet;
  TH2D * hEff[ntrkbins];
  double centval;
  int ntrkval;
  double vtx;
  int Noff;
  double reso_;
  bool bCaloMatching_;
  int nvtx_;
  double minvz_;
  double maxvz_;
  double dzerr_;
  double chi2_;

  double dzdzerror_;
  double d0d0error_;
  double pterror_;

  Double_t epang[NumEPNames];
  Double_t eporig[NumEPNames];
  Double_t epsin[NumEPNames];
  Double_t epcos[NumEPNames];

  Double_t qx[NumEPNames];
  Double_t qy[NumEPNames];
  Double_t q[NumEPNames];
  Double_t epmult[NumEPNames];
  Double_t sumw[NumEPNames];
  Double_t sumw2[NumEPNames];
  Double_t vn[NumEPNames];

  Double_t rescor[NumEPNames];
  Double_t rescorErr[NumEPNames];
  TH1D * hPsi[NumEPNames];
  TH1D * hPsiOffset[NumEPNames];
  TH1D * hPsiFlat[NumEPNames];


  unsigned int runno_;

  TH1D * hNtrkoff;
  int nEtaBins;
  TH1I * hrun;
  string rpnames[NumEPNames];
  string effTable_;
  TTree * tree;
  TrackEfficiency *teff;
  int FlatOrder_;
  int NumFlatBins_;
  int CentBinCompression_;
  int Noffmin_;
  int Noffmax_;
  int EPOrder_;
  TH2D * qxtrk;
  TH2D * qytrk;
  TH2D * qxtrk3;
  TH2D * qytrk3;
  TH2D * qcnt;
  TH2D * avpt;
  HiEvtPlaneFlatten * flat[NumEPNames];
  bool loadDB_;
  bool useNtrkBins_; 
  bool bypassCentrality_;
  bool FirstEvent_;
  bool MB_;
  int minrun_;
  int maxrun_;
  TH2D * wqxtrkRef[7][40];
  TH2D * wqytrkRef[7][40];


  int ntrack;
  float sppt[MaxTracks];
  float spphi[MaxTracks];
  float speta[MaxTracks];

  int getNoff(const edm::Event& iEvent, const edm::EventSetup& iSetup, int bin)
  {
    int Noff = 0;
    using namespace edm;
    using namespace reco;
    qxtrk->Reset();
    qytrk->Reset();
    if(EPOrder_ == 2) {
      qxtrk3->Reset();
      qytrk3->Reset();
    }
    qcnt->Reset();
    avpt->Reset();
    
    iEvent.getByToken(vertexToken,vertex_);
    VertexCollection recoVertices = *vertex_;
    if ( recoVertices.size() > 100 ) return -1;
    sort(recoVertices.begin(), recoVertices.end(), [](const reco::Vertex &a, const reco::Vertex &b){
	if ( a.tracksSize() == b.tracksSize() ) return a.chi2() < b.chi2();
	return a.tracksSize() > b.tracksSize();
      });
   
    int primaryvtx = 0;
   
    double vz = recoVertices[primaryvtx].z();
    if (fabs(vz) < -15 || fabs(vz) > 15) {
      return -1;
    }
 
    math::XYZPoint v1( recoVertices[primaryvtx].position().x(), recoVertices[primaryvtx].position().y(), recoVertices[primaryvtx].position().z() );
    double vxError = recoVertices[primaryvtx].xError();
    double vyError = recoVertices[primaryvtx].yError();
    double vzError = recoVertices[primaryvtx].zError();
    

    iEvent.getByLabel(trackTag_,trackCollection_);

    for(TrackCollection::const_iterator itTrack = trackCollection_->begin(); itTrack != trackCollection_->end(); ++itTrack) {    
      if ( !itTrack->quality(reco::TrackBase::highPurity) ) continue;
      if ( itTrack->charge() == 0 ) continue;
      if ( fabs(itTrack->eta()) > 2.4 ) continue;
      bool bPix = false;
      int nHits = itTrack->numberOfValidHits();
      if ( itTrack->pt() < 2.4 and (nHits==3 or nHits==4 or nHits==5 or nHits==6) ) bPix = true;
      if ( not bPix ) {
	if ( nHits < 11 ) continue;
	if ( itTrack->normalizedChi2() / itTrack->hitPattern().trackerLayersWithMeasurement() > 0.15 ) {
	  continue;
	}
	if ( itTrack->ptError()/itTrack->pt() > 0.1 ) {
	  continue;
	}
	
	double d0 = -1.* itTrack->dxy(v1);
	double derror=sqrt(itTrack->dxyError()*itTrack->dxyError()+vxError*vyError);
	if ( fabs( d0/derror ) > 3.0 ) {
	  continue;
	}
	
	double dz=itTrack->dz(v1);
	double dzerror=sqrt(itTrack->dzError()*itTrack->dzError()+vzError*vzError);
	if ( fabs( dz/dzerror ) > 3.0 ) {
	  continue;
	}
      }
      hptNtrkGood->Fill(itTrack->pt());
      if(itTrack->pt()>0.4) ++Noff;
      int ipt = qxtrk->GetXaxis()->FindBin(itTrack->pt());
      int ieta = qxtrk->GetYaxis()->FindBin(itTrack->eta());
      qxtrk->Fill(itTrack->pt(), itTrack->eta(), TMath::Cos(EPOrder_*itTrack->phi()) -wqxtrkRef[EPOrder_-1][bin]->GetBinContent(ipt,ieta));
      qytrk->Fill(itTrack->pt(),itTrack->eta(), TMath::Sin(EPOrder_*itTrack->phi()) - wqytrkRef[EPOrder_-1][bin]->GetBinContent(ipt,ieta));
      if(EPOrder_ == 2) {
	qxtrk3->Fill(itTrack->pt(), itTrack->eta(), TMath::Cos(3.*itTrack->phi()) - wqxtrkRef[2][bin]->GetBinContent(ipt,ieta));
	qytrk3->Fill(itTrack->pt(), itTrack->eta(), TMath::Sin(3.*itTrack->phi()) - wqytrkRef[2][bin]->GetBinContent(ipt,ieta));
      }
      qcnt->Fill(itTrack->pt(), itTrack->eta());
      avpt->Fill(itTrack->pt(), itTrack->eta(), itTrack->pt());
      
      if( itTrack->pt() < 0.2 ) continue;
      hEff[bin]->Fill(itTrack->phi(),itTrack->eta());
    }
    return Noff;
  }


};


//
// constructors and destructor
//
VNAnalyzer::VNAnalyzer(const edm::ParameterSet& iConfig):runno_(0)
  
{
  runno_ = 0;
  loadDB_ = kTRUE;
  FirstEvent_ = kTRUE;
  for(int i = 0; i<NumEPNames; i++) {
    epang[i] = -10;
    eporig[i] = -10;
    epsin[i] = 0;
    epcos[i] = 0;
    qx[i] = 0;
    qy[i] = 0;
    q[i] = 0;
    epmult[i] = 0;
    rescor[i] = 0;
    rescorErr[i] = 0;
  }

  centralityVariable_ = iConfig.getParameter<std::string>("centralityVariable");
  if(iConfig.exists("nonDefaultGlauberModel")){
    centralityMC_ = iConfig.getParameter<std::string>("nonDefaultGlauberModel");
  }
  centralityLabel_ = centralityVariable_+centralityMC_;

  centralityBinTag_ = iConfig.getParameter<edm::InputTag>("centralityBinTag_");
  centralityBinToken = consumes<int>(centralityBinTag_);

  centralityTag_ = iConfig.getParameter<edm::InputTag>("centralityTag_");
  centralityToken = consumes<reco::Centrality>(centralityTag_);
  if(centralityToken.isUninitialized()) {
    std::cout<<"centralityToken is uninitialized."<<std::endl;
  }
  vertexTag_  = iConfig.getParameter<edm::InputTag>("vertexTag_");
  vertexToken = consumes<std::vector<reco::Vertex>>(vertexTag_);
  if(vertexToken.isUninitialized()) {
    std::cout<<"vertexToken is uninitialized."<<std::endl;
  }

  trackTag_ = iConfig.getParameter<edm::InputTag>("trackTag_");
  trackToken = consumes<reco::TrackCollection>(trackTag_);
  if(trackToken.isUninitialized()) {
    std::cout<<"trackToken is uninitialized."<<std::endl;
  }
  useNtrk_ = iConfig.getUntrackedParameter<bool>("useNtrk",false);
  if(useNtrk_) {
    NumFlatBins_ = ntrkbins;
    CentBinCompression_ = 1;
  }

  inputPlanesTag_ = iConfig.getParameter<edm::InputTag>("inputPlanesTag_");
  inputPlanesToken = consumes<reco::EvtPlaneCollection>(inputPlanesTag_);
  if(inputPlanesToken.isUninitialized()) {
    std::cout<<"inputPlanesToken is uninitialized."<<std::endl;
  }
  tag_ = consumes<int>(iConfig.getParameter<edm::InputTag>("BinLabel"));

  EPOrder_ = iConfig.getUntrackedParameter<int>("EPOrder_",2);
  FlatOrder_ = iConfig.getUntrackedParameter<int>("FlatOrder_", 9);
  NumFlatBins_ = iConfig.getUntrackedParameter<int>("NumFlatBins_",20);
  caloCentRef_ = iConfig.getUntrackedParameter<double>("caloCentRef_",80.);
  caloCentRefWidth_ = iConfig.getUntrackedParameter<double>("caloCentRefWidth_",5.);
  CentBinCompression_ = iConfig.getUntrackedParameter<int>("CentBinCompression_",5);
  Noffmin_ = iConfig.getUntrackedParameter<int>("Noffmin_", 0);
  Noffmax_ = iConfig.getUntrackedParameter<int>("Noffmax_", 50000);	
  minrun_ = iConfig.getUntrackedParameter<int>("minrun_", 0);
  maxrun_ = iConfig.getUntrackedParameter<int>("maxrun_", 50000);	
  effTable_ = iConfig.getParameter<std::string>("effTable_");
  bCaloMatching_ = iConfig.getUntrackedParameter<bool>("bCaloMaching", true);
  MB_ = iConfig.getUntrackedParameter<bool>("MB_",true);

  nvtx_ = iConfig.getUntrackedParameter<int>("nvtx_", 100);
  reso_ = iConfig.getUntrackedParameter<double>("reso", 0.2);
  if(reso_<0.01) bCaloMatching_ = false;
  if(bCaloMatching_) {
    pfTag_ = iConfig.getUntrackedParameter<edm::InputTag>("pfTag");
    pfToken_ = consumes<reco::PFCandidateCollection>(pfTag_);
  }
  dzdzerror_ = iConfig.getUntrackedParameter<double>("dzdzerror_", 3.);
  d0d0error_ = iConfig.getUntrackedParameter<double>("d0d0error_", 3.);
  pterror_ = iConfig.getUntrackedParameter<double>("pterror_",0.1);
  teff = 0;
  if(!effTable_.empty()) teff = new TrackEfficiency(effTable_.data());
  minvz_ = iConfig.getUntrackedParameter<double>("minvz_", -15.);
  maxvz_ = iConfig.getUntrackedParameter<double>("maxvz_", 15.);
  dzerr_ = iConfig.getParameter<double>("dzerr") ;
  chi2_  = iConfig.getParameter<double>("chi2") ;
  offsetFileName = iConfig.getUntrackedParameter<std::string>("offsetFile");
  frecenter = new TFile(offsetFileName.data(),"read");
  int mx = ntrkbins;
  if(!useNtrk_) {
    mx = 0;
    cout<<"need to set this up"<<endl;
    //	mx = nCentBins_;
  }
  for(int i = 0; i<mx; i++) {
    for(int j = 1; j<=7; j++){
      wqxtrkRef[j-1][i] = (TH2D *) frecenter->Get(Form("wqxtrk_%d_%d",j,i));
      wqytrkRef[j-1][i] = (TH2D *) frecenter->Get(Form("wqytrk_%d_%d",j,i));
    }
  }

  std::cout<<"==============================================="<<std::endl;
  std::cout<<"centralityBinTag_           "<<centralityBinTag_.encode()<<std::endl;
  std::cout<<"centralityTag_              "<<centralityTag_.encode()<<std::endl;
  std::cout<<"vertexTag_                  "<<vertexTag_.encode()<<std::endl;
  std::cout<<"trackTag_                   "<<trackTag_.encode()<<std::endl;
  std::cout<<"inputPlanesTag_             "<<inputPlanesTag_.encode()<<std::endl;
  std::cout<<"EPOrder_                    "<<EPOrder_<<std::endl;
  std::cout<<"FlatOrder_                  "<<FlatOrder_<<std::endl;
  std::cout<<"NumFlatBins_                "<<NumFlatBins_<<std::endl;
  std::cout<<"caloCentRef_                "<<caloCentRef_<<std::endl;
  std::cout<<"caloCentRefWidth_           "<<caloCentRefWidth_<<std::endl;
  std::cout<<"CentBinCompression_         "<<CentBinCompression_<<std::endl;
  std::cout<<"Noffmin_                    "<<Noffmin_<<std::endl;
  std::cout<<"Noffmax_                    "<<Noffmax_<<std::endl;
  std::cout<<"minrun_                     "<<minrun_<<std::endl;
  std::cout<<"maxrun_                     "<<maxrun_<<std::endl;
  std::cout<<"effTable_                   "<<effTable_<<std::endl;
  std::cout<<"dzerror_                    "<<dzdzerror_<<endl;
  std::cout<<"d0error_                    "<<d0d0error_<<endl;
  std::cout<<"pterror_                    "<<pterror_<<endl;
  std::cout<<"nvtx_                       "<<nvtx_<<endl;
  if(bCaloMatching_) { 
    std::cout<<"bCaloMatching_              true"<<std::endl;
    std::cout<<"reso_                     "<<reso_<<std::endl;   
  }
  std::cout<<"dzerr_                       "<<dzerr_<<std::endl;
  std::cout<<"chi2_                        "<<chi2_<<std::endl;
  std::cout<<"==============================================="<<std::endl;

  hNtrkoff = fs->make<TH1D>("Ntrkoff","Ntrkoff",1001,0,3000);
  int npt = nptbins;
  qxtrk = fs->make<TH2D>(Form("qxtrk_v%d",EPOrder_),Form("qxtrk_v%d",EPOrder_),npt,ptbins, netabinsDefault, etabinsDefault);
  qytrk = fs->make<TH2D>(Form("qytrk_v%d",EPOrder_),Form("qytrk_v%d",EPOrder_),npt,ptbins, netabinsDefault, etabinsDefault);
  if(EPOrder_ == 2) {
    qxtrk3 = fs->make<TH2D>("qxtrk_v3","qxtrk_v3",npt,ptbins, netabinsDefault, etabinsDefault);
    qytrk3 = fs->make<TH2D>("qytrk_v3","qytrk_v3",npt,ptbins, netabinsDefault, etabinsDefault);
  }
  qcnt =  fs->make<TH2D>(Form("qcnt_v%d",EPOrder_), Form("qcnt_v%d",EPOrder_),npt,ptbins, netabinsDefault, etabinsDefault);
  avpt =  fs->make<TH2D>("avpt","avpt",npt,ptbins, netabinsDefault, etabinsDefault);

  hcent = fs->make<TH1D>("cent","cent",220,-10,110);
  hcentbins = fs->make<TH1D>("centbins","centbins",201,0,200);
  hrun = fs->make<TH1I>("runs","runs",maxrun_-minrun_+1,minrun_,maxrun_);
  hptNtrk = fs->make<TH1D>("ptNtrk","ptNtrk",npt,ptbins);
  hptNtrk->SetXTitle("p_{T} (GeV/c)");
  hptNtrk->SetYTitle("Ntrks (|#eta|<1; 0-5)");
  hptNtrkGood = fs->make<TH1D>("ptNtrkGood","ptNtrkGood",npt,ptbins);
  hptNtrkGood->SetXTitle("p_{T} (GeV/c)");
  hptNtrkGood->SetYTitle("Ntrks (Good) (|#eta|<1; 0-5)");
  hNtrkRet = fs->make<TH1I>("NtrkRet","NtrkRet", 10,0,10);
  for(int i = 0; i<10; i++) {
    TString hn = Form("Eff_%d_%d",10*i,10*(i+1));
    hEff[i] = fs->make<TH2D>(hn.Data(),hn.Data(),50,-TMath::Pi(),TMath::Pi(),50,-2.4,2.4);
    hEff[i]->Sumw2();
    hEff[i]->SetXTitle("#phi (radians)");
    hEff[i]->SetYTitle("#eta");
    hEff[i]->SetOption("colz");
  }
  TString epnames = EPNames[0].data();
  epnames = epnames+"/D";
  NumFlatBins_ = ntrkbins;
  for(int i = 0; i<NumEPNames; i++) {
    if(i>0) epnames = epnames + ":" + EPNames[i].data() + "/D";
    TFileDirectory subdir = fs->mkdir(Form("%s",EPNames[i].data()));
    flat[i] = new HiEvtPlaneFlatten();
    flat[i]->init(FlatOrder_,NumFlatBins_,EPNames[i],EPOrder[i]);
    Double_t psirange = 4;
    if(EPOrder[i]==1 ) psirange = 3.5;
    if(EPOrder[i]==2 ) psirange = 2;
    if(EPOrder[i]==3 ) psirange = 1.5;
    if(EPOrder[i]==4 ) psirange = 1;
    if(EPOrder[i]==5) psirange = 0.8;
    if(EPOrder[i]==6) psirange = 0.6;

    hPsi[i] = subdir.make<TH1D>("psi","psi",800,-psirange,psirange);
    hPsi[i]->SetXTitle("#Psi");
    hPsi[i]->SetYTitle(Form("Counts (cent<80%c)",'%'));
    
    hPsiOffset[i] = subdir.make<TH1D>("psiOffset","psiOffset",800,-psirange,psirange);
    hPsiOffset[i]->SetXTitle("#Psi");
    hPsiOffset[i]->SetYTitle(Form("Counts (cent<80%c)",'%'));

    
    hPsiFlat[i] = subdir.make<TH1D>("psiFlat","psiFlat",800,-psirange,psirange);
    hPsiFlat[i]->SetXTitle("#Psi");
    hPsiFlat[i]->SetYTitle(Form("Counts (cent<80%c)",'%'));

  }
  std::cout<<"Done with flat init"<<std::endl;
  tree = fs->make<TTree>("tree","EP tree");

  tree->Branch("Cent",&centval,"cent/D");
  tree->Branch("NtrkOff",&Noff,"Noff/I");
  tree->Branch("ntrkflat",&ntrkval,"nofftrak/I");
  tree->Branch("Vtx",&vtx,"vtx/D");
  tree->Branch("epang",&epang, epnames.Data());
  tree->Branch("eporig",&eporig, epnames.Data());
  tree->Branch("qx",      &qx,       epnames.Data());
  tree->Branch("qy",      &qy,       epnames.Data());
  tree->Branch("q",       &q,       epnames.Data());
  tree->Branch("vn", &vn, epnames.Data());
  tree->Branch("mult",    &epmult,  epnames.Data());
  tree->Branch("sumw",    &sumw,  epnames.Data());
  tree->Branch("sumw2",    &sumw2,  epnames.Data());
  tree->Branch("Run",     &runno_,   "run/i");
  tree->Branch("Rescor",  &rescor,   epnames.Data());
  tree->Branch("RescorErr",  &rescorErr,   epnames.Data());
  tree->Branch(Form("qxtrk_v%d",EPOrder_),   "TH2D",  &qxtrk, 128000, 0);
  tree->Branch(Form("qytrk_v%d",EPOrder_),   "TH2D",  &qytrk, 128000, 0);
  if(EPOrder_ == 2) {
    tree->Branch("qxtrk_v3", "TH2D",  &qxtrk3, 128000, 0);
    tree->Branch("qytrk_v3", "TH2D",  &qytrk3, 128000, 0);
  }
  tree->Branch(Form("qcnt_v%d",EPOrder_),    "TH2D",  &qcnt, 128000, 0);
  tree->Branch("avpt",    "TH2D",  &avpt, 128000, 0);
}



VNAnalyzer::~VNAnalyzer()
{
  frecenter->Close();  
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
VNAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  using namespace std;
  using namespace reco;
  Bool_t newrun = kFALSE;
  if(runno_ != iEvent.id().run()) newrun = kTRUE;
  runno_ = iEvent.id().run();
  hrun->Fill(runno_);
  if(FirstEvent_ || newrun) {
    FirstEvent_ = kFALSE;
    newrun = kFALSE;
    //
    //Get Size of Centrality Table
    //
    if(!useNtrk_) {
      edm::ESHandle<CentralityTable> centDB_;
      iSetup.get<HeavyIonRcd>().get(centralityLabel_,centDB_);
      nCentBins_ = (int) centDB_->m_table.size();
      for(int i = 0; i<NumEPNames; i++) {
	flat[i]->setCaloCentRefBins(-1,-1);
	if(caloCentRef_>0) {
	  int minbin = (caloCentRef_-caloCentRefWidth_/2.)*nCentBins_/100.;
	  int maxbin = (caloCentRef_+caloCentRefWidth_/2.)*nCentBins_/100.;
	  minbin/=CentBinCompression_;
	  maxbin/=CentBinCompression_;
	  if(minbin>0 && maxbin>=minbin) {
	    if(EPDet[i]==HF || EPDet[i]==Castor) flat[i]->setCaloCentRefBins(minbin,maxbin);
	  }
	}
      }
    }
    //
    //Get flattening parameter file.  
    //
    edm::ESHandle<RPFlatParams> flatparmsDB_;
    iSetup.get<HeavyIonRPRcd>().get(flatparmsDB_);
    LoadEPDB * db = new LoadEPDB(flatparmsDB_,flat);
    if(!db->IsSuccess()) {
 	std::cout<<"Flattening db inconsistancy, will continue without: "<<std::endl;
     loadDB_ = kFALSE;
    }        
  } //First event
  
  
  // //
  // //Get Centrality
  // //

   int Noff=0;

  int bin = 0;
  if(!useNtrk_) {
    int ntrkval=0;
    if(Noffmin_>=0) {
      iEvent.getByToken(centralityToken, centrality_);
      ntrkval = centrality_->Ntracks();
      if ( (ntrkval < Noffmin_) || (ntrkval >= Noffmax_) ) {
	return;
      }
    }

   iEvent.getByToken(centralityBinToken, cbin_);
   int cbin = *cbin_;
   bin = cbin/CentBinCompression_; 
   double cscale = 100./nCentBins_;
   centval = cscale*cbin;
   
   
  } else {
    iEvent.getByToken(tag_,cbin_);
    ntrkval = *cbin_;
    hNtrkoff->Fill(ntrkval);
    bin = NtrkToBin(ntrkval)-1;
    centval = bin;
  }

  hcent->Fill(centval);
  hcentbins->Fill(bin);
  // //
  // //Get Vertex
  // //
  iEvent.getByToken(vertexToken,vertex_);
  if(!vertex_.isValid()) {
    std::cout<<"Error! Can't get vertex!"<<std::endl;
    return;
  }

  iEvent.getByToken(vertexToken, vertex_);
  VertexCollection recoVertices = *vertex_;
  if ( recoVertices.size() > 100 ) return;
  sort(recoVertices.begin(), recoVertices.end(), [](const reco::Vertex &a, const reco::Vertex &b){
      if ( a.tracksSize() == b.tracksSize() ) return a.chi2() < b.chi2();
      return a.tracksSize() > b.tracksSize();
    });
  
  int primaryvtx = 0;
  
  double vz = recoVertices[primaryvtx].z();
  if (fabs(vz) < -15 || fabs(vz) > 15) {
    return;
  }
  vtx = vz; 
  // //
  // //Get Event Planes
  // //
  iEvent.getByToken(inputPlanesToken,inputPlanes_);
  
  if(!inputPlanes_.isValid()){
     cout << "Error! Can't get hiEvtPlaneFlat product!" << endl;
     return ;
   }
  
   Int_t indx = 0;
   for(int i = 0; i<NumEPNames; i++) {
     epang[i] = -10;
     epsin[i] = 0;
     epcos[i] = 0;
     qx[i] = 0;
     qy[i] = 0;
     q[i] = 0;
     vn[i] = 0;
     epmult[i] = 0;
     sumw[i] = 0;
     sumw2[i] = 0;
   }
   for (EvtPlaneCollection::const_iterator rp = inputPlanes_->begin();rp !=inputPlanes_->end(); rp++) {
     if(indx != rp->indx() ) {
       cout<<"EP inconsistency found. Abort."<<endl;
       return;
     }
     if(rp->sumSin()!=0 || rp->sumCos()!=0) {
       epang[indx]=rp->angle();
       eporig[indx]=rp->angle(0);
       epsin[indx] = rp->sumSin();
       epcos[indx] = rp->sumCos();
       if(centval<80 && fabs(vtx)<15) {
	 hPsi[indx]->Fill(rp->angle(0));
	 hPsiOffset[indx]->Fill(rp->angle(1));
	 hPsiFlat[indx]->Fill(rp->angle(2));
       }
      
       qx[indx]  = rp->qx(); 
       qy[indx]  = rp->qy();
       q[indx]   = rp->q();
       vn[indx] = rp->vn(0);
       epmult[indx] = (double) rp->mult();
       sumw[indx] = rp->sumw();
       sumw2[indx] = rp->sumw2();
       
       rescor[indx] = flat[indx]->getCentRes1((int) centval);
       rescorErr[indx] = flat[indx]->getCentResErr1((int) centval);
     }
     ++indx; 
   }


  ntrkval = Noff;
  if ( Noff == -2 ) {
    return;
  }
  getNoff(iEvent, iSetup, bin);
  hNtrkoff->Fill(ntrkval);
  
  tree->Fill(); 
}



// ------------ method called once each job just before starting event loop  ------------
void 
VNAnalyzer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
VNAnalyzer::endJob() {
}
bool
VNAnalyzer::CaloMatch(const reco::Track & track, const edm::Event & iEvent, unsigned int idx)
{
  if ( !bCaloMatching_ ) return true;
  edm::Handle<reco::PFCandidateCollection> pfCand;
  iEvent.getByToken( pfToken_, pfCand );
  double energy = 0;
  for ( reco::PFCandidateCollection::const_iterator it = pfCand->begin(); it != pfCand->end(); ++it ) {
    reco::TrackRef trackRef = it->trackRef();
    if ( !((it->particleId() != reco::PFCandidate::h) ||
	 (it->particleId() != reco::PFCandidate::e) ||
	   (it->particleId() != reco::PFCandidate::mu) )) continue;
    if ( idx == trackRef.key() ) {
      energy = it->ecalEnergy() + it->hcalEnergy();
      break;
    }
  }
  
  if( track.pt() < 20 || ( energy/( track.pt()*TMath::CosH(track.eta() ) ) > reso_ && (energy)/(TMath::CosH(track.eta())) > (track.pt() - 80.0) )  ) return true;
  else {
    return false;
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(VNAnalyzer);

