#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "fhiclcpp/ParameterSet.h"
#include "lardata/Utilities/AssociationUtil.h"


#include "larcore/Geometry/Geometry.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "nusimdata/SimulationBase/MCFlux.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/Simulation/SimChannel.h"
#include "lardataobj/Simulation/AuxDetSimChannel.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RawData/RawDigit.h"
#include "lardataobj/RawData/raw.h"
#include "lardataobj/RawData/BeamInfo.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardata/Utilities/AssociationUtil.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "larcoreobj/SummaryData/POTSummary.h"
#include "larsim/MCCheater/BackTrackerService.h"
#include "larsim/MCCheater/ParticleInventoryService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/Wire.h"
#include "lardataobj/RecoBase/EndPoint2D.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/OpFlash.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"
#include "lardataobj/AnalysisBase/CosmicTag.h"
#include "lardataobj/AnalysisBase/FlashMatch.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"

#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TProfile.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TMath.h"
#include "TGraphErrors.h"
#include "TMinuit.h"
#include "TString.h"
#include "TTimeStamp.h"
#include "TVectorD.h"
#include "TCanvas.h"
#include "TFrame.h"
#include "TLine.h"
#include "TAxis.h"
#include "TTimeStamp.h"

#include <vector>
#include <fstream>
#include "TPaveStats.h"
#include <iostream>
#include <string>
#include "math.h"
#include "stdio.h"
#include <iterator>

const int kMaxTracks  = 30;

using namespace std;

namespace dune{

  class michelremovingpdhd : public art::EDAnalyzer {
  public:

    explicit michelremovingpdhd(fhicl::ParameterSet const& pset);
    virtual ~michelremovingpdhd();

    void beginJob();
    void endJob();
    void beginRun(const art::Run& run);
    void analyze(const art::Event& evt);
    void reset();

  private:
    TTree* fEventTree;
    bool     isData;
    Int_t    run;
    Int_t    subrun;
    Int_t    event;
    Double_t evttime;
    Int_t    year_month_date;
    Int_t    hour_min_sec;
    Int_t    cross_trks; //are the total number of good T0 tagged tracks stopping + crossing
    Int_t stopping_trks;
    Int_t    all_trks;
    Int_t unbroken_trks; //these are unbroken stopping tracks
    Float_t  trackthetaxz[kMaxTracks];
    Float_t  trackthetayz[kMaxTracks];
    Float_t trkstartx[kMaxTracks];
    Float_t trkstarty[kMaxTracks];
    Float_t trkstartz[kMaxTracks];
    Float_t trkendx[kMaxTracks];
    Float_t trkendy[kMaxTracks];
    Float_t trkendz[kMaxTracks];
    Float_t trklen[kMaxTracks];
    Int_t    TrkID[kMaxTracks];
    Float_t  trkstartcosxyz[kMaxTracks][3];
    Float_t  trkendcosxyz[kMaxTracks][3];
    Int_t    ntrkhits[kMaxTracks][3];
    Float_t  trkdqdx[kMaxTracks][3][3000];
    Float_t  trkdedx[kMaxTracks][3][3000];
    Float_t  trkresrange[kMaxTracks][3][3000];
    Float_t  trkhitx[kMaxTracks][3][3000];
    Float_t  trkhity[kMaxTracks][3][3000];
    Float_t  trkhitz[kMaxTracks][3][3000];
    Float_t  trkpitch[kMaxTracks][3][3000];
    Float_t  trkhitpeakt[kMaxTracks][3][3000];
    Int_t    trkhittpc[kMaxTracks][3][3000];
    Int_t    trkhitwire[kMaxTracks][3][3000];
    Float_t  trkhitrms[kMaxTracks][3][3000];
    Float_t  trkhitdeltat[kMaxTracks][3][3000];
    Float_t peakT_max[kMaxTracks];
    Float_t peakT_min[kMaxTracks];
    Float_t dist_min[kMaxTracks];
    Int_t adjacent_hits[kMaxTracks];
    Int_t lastwire[kMaxTracks];
    Int_t endtpc[kMaxTracks];
    Float_t lastpeakt[kMaxTracks];
    Int_t trkanodetoanode[kMaxTracks];
    float true_trkstartx[kMaxTracks];
    float true_trkstarty[kMaxTracks];
    float true_trkstartz[kMaxTracks];
    float true_trkendx[kMaxTracks];
    float true_trkendy[kMaxTracks];
    float true_trkendz[kMaxTracks];
    std::string fHitsModuleLabel;
    std::string fTrackModuleLabel;
    std::string fCalorimetryModuleLabel;
    std::string fPFParticleModuleLabel;
    std::string fT0ModuleLabel;
    bool  fSaveCaloInfo;
    bool  fSaveTrackInfo;
    bool  fRequireT0;

    // cumulative drop counters, reported in endJob to diagnose empty output
    long fNTracksSeen;
    long fNPassPFP;
    long fNPassT0;
    long fNWritten;

    // drift geometry, resolved per run, used to flag anode-to-anode crossers
    int    fDriftAxis;
    double fCathodeX;
    double fAnodeLoX;
    double fAnodeHiX;
    bool   fCathodeValid;

    protoana::ProtoDUNETruthUtils truthUtil;

  };

  //========================================================================
  michelremovingpdhd::michelremovingpdhd(fhicl::ParameterSet const& pset) :
    EDAnalyzer(pset),
    fHitsModuleLabel          (pset.get< std::string >("HitsModuleLabel","")         ),
    fTrackModuleLabel         (pset.get< std::string >("TrackModuleLabel","")        ),
    fCalorimetryModuleLabel   (pset.get< std::string >("CalorimetryModuleLabel","")  ),
    fPFParticleModuleLabel    (pset.get< std::string >("PFParticleModuleLabel","pandora") ),
    fT0ModuleLabel            (pset.get< std::string >("T0ModuleLabel","pandora")     ),
    fSaveCaloInfo             (pset.get< bool>("SaveCaloInfo",false)),
    fSaveTrackInfo            (pset.get< bool>("SaveTrackInfo",false)),
    fRequireT0                (pset.get< bool>("RequireT0",true)),
    fNTracksSeen(0), fNPassPFP(0), fNPassT0(0), fNWritten(0),
    fDriftAxis                (0),
    fCathodeX                 (0.),
    fAnodeLoX                 (0.),
    fAnodeHiX                 (0.),
    fCathodeValid             (false)
  {
    if (fSaveTrackInfo == false) fSaveCaloInfo = false;
  }

  //========================================================================
  michelremovingpdhd::~michelremovingpdhd(){
  }
  //========================================================================

  //========================================================================
  void michelremovingpdhd::beginJob(){
    std::cout<<"job begin..."<<std::endl;
    art::ServiceHandle<art::TFileService> tfs;
    fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");
    fEventTree->Branch("isData", &isData, "isData/O");
    fEventTree->Branch("event", &event,"event/I");
    fEventTree->Branch("evttime",&evttime,"evttime/D");
    fEventTree->Branch("run", &run,"run/I");
    fEventTree->Branch("subrun", &subrun,"surbrun/I");
    fEventTree->Branch("year_month_date", &year_month_date,"year_month_date/I");
    fEventTree->Branch("hour_min_sec", &hour_min_sec,"hour_min_sec/I");
    fEventTree->Branch("cross_trks",&cross_trks,"cross_trks/I");
    fEventTree->Branch("stopping_trks",&stopping_trks,"stopping_trks/I");
    fEventTree->Branch("all_trks",&all_trks,"all_trks/I");
    fEventTree->Branch("unbroken_trks",&unbroken_trks,"unbroken_trks/I");
    fEventTree->Branch("trackthetaxz",trackthetaxz,"trackthetaxz[cross_trks]/F");
    fEventTree->Branch("trackthetayz",trackthetayz,"trackthetayz[cross_trks]/F");
    fEventTree->Branch("trkstartx",trkstartx,"trkstartx[cross_trks]/F");
    fEventTree->Branch("trkstarty",trkstarty,"trkstarty[cross_trks]/F");
    fEventTree->Branch("trkstartz",trkstartz,"trkstartz[cross_trks]/F");
    fEventTree->Branch("trkendx",trkendx,"trkendx[cross_trks]/F");
    fEventTree->Branch("trkendy",trkendy,"trkendy[cross_trks]/F");
    fEventTree->Branch("trkendz",trkendz,"trkendz[cross_trks]/F");
    fEventTree->Branch("trklen",trklen,"trklen[cross_trks]/F");
    fEventTree->Branch("peakT_max",peakT_max,"peakT_max[cross_trks]/F");
    fEventTree->Branch("peakT_min",peakT_min,"peakT_min[cross_trks]/F");
    fEventTree->Branch("TrkID",TrkID,"TrkID[cross_trks]/I");
    fEventTree->Branch("trkstartcosxyz",trkstartcosxyz,"trkstartcosxyz[cross_trks][3]/F");
    fEventTree->Branch("trkendcosxyz",trkendcosxyz,"trkendcosxyz[cross_trks][3]/F");
    fEventTree->Branch("ntrkhits",ntrkhits,"ntrkhits[cross_trks][3]/I");
    fEventTree->Branch("trkdqdx",trkdqdx,"trkdqdx[cross_trks][3][3000]/F");
    fEventTree->Branch("trkdedx",trkdedx,"trkdedx[cross_trks][3][3000]/F");
    fEventTree->Branch("trkresrange",trkresrange,"trkresrange[cross_trks][3][3000]/F");
    fEventTree->Branch("trkhitx",trkhitx,"trkhitx[cross_trks][3][3000]/F");
    fEventTree->Branch("trkhity",trkhity,"trkhity[cross_trks][3][3000]/F");
    fEventTree->Branch("trkhitz",trkhitz,"trkhitz[cross_trks][3][3000]/F");
    fEventTree->Branch("trkpitch",trkpitch,"trkpitch[cross_trks][3][3000]/F");
    fEventTree->Branch("trkhitpeakt",trkhitpeakt,"trkhitpeakt[cross_trks][3][3000]/F");
    fEventTree->Branch("trkhittpc",trkhittpc,"trkhittpc[cross_trks][3][3000]/I");
    fEventTree->Branch("trkhitwire",trkhitwire,"trkhitwire[cross_trks][3][3000]/I");
    fEventTree->Branch("trkhitrms",trkhitrms,"trkhitrms[cross_trks][3][3000]/F");
    fEventTree->Branch("trkhitdeltat",trkhitdeltat,"trkhitdeltat[cross_trks][3][3000]/F");
    fEventTree->Branch("dist_min",dist_min,"dist_min[cross_trks]/F");
    fEventTree->Branch("adjacent_hits",adjacent_hits,"adjacent_hits[cross_trks]/I");
    fEventTree->Branch("lastwire",lastwire,"lastwire[cross_trks]/I");
    fEventTree->Branch("lastpeakt",lastpeakt,"lastpeakt[cross_trks]/F");
    fEventTree->Branch("endtpc",endtpc,"endtpc[cross_trks]/I");
    fEventTree->Branch("trkanodetoanode",trkanodetoanode,"trkanodetoanode[cross_trks]/I");
    fEventTree->Branch("true_trkstartx",true_trkstartx,"true_trkstartx[cross_trks]/F");
    fEventTree->Branch("true_trkstarty",true_trkstarty,"true_trkstarty[cross_trks]/F");
    fEventTree->Branch("true_trkstartz",true_trkstartz,"true_trkstartz[cross_trks]/F");
    fEventTree->Branch("true_trkendx",true_trkendx,"true_trkendx[cross_trks]/F");
    fEventTree->Branch("true_trkendy",true_trkendy,"true_trkendy[cross_trks]/F");
    fEventTree->Branch("true_trkendz",true_trkendz,"true_trkendz[cross_trks]/F");

  }

  //========================================================================
  void michelremovingpdhd::endJob(){
    // If the output tree is empty, this tells you which gate dropped everything.
    std::cout<<"[michelremovingpdhd] tracks seen="<<fNTracksSeen
             <<"  pass PFParticle="<<fNPassPFP
             <<"  pass T0="<<fNPassT0
             <<"  written(cross_trks)="<<fNWritten
             <<"   (RequireT0="<<fRequireT0
             <<", PFParticleModuleLabel="<<fPFParticleModuleLabel
             <<", T0ModuleLabel="<<fT0ModuleLabel
             <<", CalorimetryModuleLabel="<<fCalorimetryModuleLabel<<")"<<std::endl;
  }

  //========================================================================
  void michelremovingpdhd::beginRun(const art::Run&){
    mf::LogInfo("michelremovingpdhd")<<"begin run..."<<std::endl;

    // Resolve the drift axis and the anode/cathode drift coordinates from the
    // geometry so anode-to-anode crossers can be flagged without hardcoding
    // detector dimensions.
    auto const& geom = *art::ServiceHandle<geo::Geometry>();

    fDriftAxis = 0;
    for (auto const& tpc : geom.Iterate<geo::TPCGeo>()){
      auto const d = tpc.DriftDir();
      double ax=std::abs(d.X()), ay=std::abs(d.Y()), az=std::abs(d.Z());
      if (ax>=ay && ax>=az)      fDriftAxis = 0;
      else if (ay>=ax && ay>=az) fDriftAxis = 1;
      else                       fDriftAxis = 2;
      break;
    }

    auto boxMin=[&](geo::BoxBoundedGeo const& b){ return fDriftAxis==0?b.MinX():(fDriftAxis==1?b.MinY():b.MinZ()); };
    auto boxMax=[&](geo::BoxBoundedGeo const& b){ return fDriftAxis==0?b.MaxX():(fDriftAxis==1?b.MaxY():b.MaxZ()); };

    double gmin=1e12, gmax=-1e12;
    for (auto const& tpc : geom.Iterate<geo::TPCGeo>()){
      auto const box = tpc.BoundingBox();
      gmin = std::min(gmin, boxMin(box));
      gmax = std::max(gmax, boxMax(box));
    }
    double gc = 0.5*(gmin+gmax);
    fAnodeLoX = -360;
    fAnodeHiX = 360;
    double loMax=-1e12, hiMin=1e12;
    bool haveLo=false, haveHi=false;
    for (auto const& tpc : geom.Iterate<geo::TPCGeo>()){
      auto const box = tpc.BoundingBox();
      double c = 0.5*(boxMin(box)+boxMax(box));
      if (c < gc){ loMax = std::max(loMax, boxMax(box)); haveLo=true; }
      else       { hiMin = std::min(hiMin, boxMin(box)); haveHi=true; }
    }
    if (haveLo && haveHi){
      fCathodeX = 0.5*(loMax+hiMin);
      fCathodeValid = true;
    } else {
      fCathodeX = gc;
      fCathodeValid = false;
    }
  }
  //========================================================================

  //========================================================================

  //========================================================================

  void michelremovingpdhd::analyze( const art::Event& evt){
    reset();
    std::vector<art::Ptr<recob::Track> > tracklist;
    auto trackListHandle = evt.getHandle< std::vector<recob::Track> >("pandoraTrack");
    if (trackListHandle){
      art::fill_ptr_vector(tracklist, trackListHandle);
    }
    else return;

    std::vector<art::Ptr<recob::PFParticle> > pfplist;
    auto PFPListHandle = evt.getHandle< std::vector<recob::PFParticle> >(fPFParticleModuleLabel);
    if (PFPListHandle) art::fill_ptr_vector(pfplist, PFPListHandle);

    /******new lines*************************/
    std::vector<art::Ptr<recob::Hit>> hitlist;
    auto hitListHandle = evt.getHandle< std::vector<recob::Hit> >(fHitsModuleLabel) ; // to get information about the hits
    if (hitListHandle) art::fill_ptr_vector(hitlist, hitListHandle);

    art::FindManyP<recob::Hit, recob::TrackHitMeta> fmthm(trackListHandle, evt, fTrackModuleLabel); // to associate tracks and hits

    art::FindManyP<recob::Hit> fmtht(trackListHandle, evt, fTrackModuleLabel); // to associate tracks and hits
    // art::FindManyP<recob::Track, recob::TrackHitMeta> thass(hitListHandle, evt, fTrackModuleLabel); //to associate hit
    art::FindManyP<recob::Track> thass(hitListHandle, evt, fTrackModuleLabel); //to associate hit just trying

    art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
    if (!fmcal.isValid()){
      std::cout<<"art::Assns<recob::Track,anab::Calorimetry,void> with module label: "<<fCalorimetryModuleLabel<<" not found."<<std::endl;
      return;
    }

    art::FindManyP<anab::T0> trk_t0_assn_v(PFPListHandle, evt , fT0ModuleLabel);


    art::FindManyP<recob::PFParticle> pfp_trk_assn(trackListHandle,evt,fTrackModuleLabel);
    art::FindManyP<anab::T0> fmT0(trackListHandle, evt ,"pmtrack");

    // Get services
    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);

    isData = evt.isRealData();
    run = evt.run();
    subrun = evt.subRun();
    event = evt.id().event();
    art::Timestamp ts = evt.time();
    TTimeStamp tts(ts.timeHigh(), ts.timeLow());
    evttime=tts.AsDouble();

    UInt_t year=0;
    UInt_t month=0;
    UInt_t day=0;

    year_month_date=tts.GetDate(kTRUE,0,&year,&month,&day);

    UInt_t hour=0;
    UInt_t min=0;
    UInt_t sec=0;

    hour_min_sec=tts.GetTime(kTRUE,0,&hour,&min,&sec);

    cross_trks=0;
    stopping_trks=0;
    all_trks=0;
    unbroken_trks=0;
    // size_t NHits=hitlist.size();
    size_t NTracks = tracklist.size();
    for(size_t i=0; i<NTracks;++i){
      if(cross_trks>=kMaxTracks) break;
      art::Ptr<recob::Track> ptrack(trackListHandle, i);
      fNTracksSeen++;
      if(fTrackModuleLabel=="pandoraTrack"){
	std::vector<art::Ptr<recob::PFParticle>> pfps=pfp_trk_assn.at(i);
	if(!pfps.size()) continue;
	fNPassPFP++;
	// Require a Pandora PFParticle->T0 association. On samples that do not
	// carry one (e.g. some ProtoDUNE-HD reco), set RequireT0:false to keep
	// all tracks; the endJob summary reports where tracks are dropped.
	if(fRequireT0){
	  std::vector<art::Ptr<anab::T0>> t0s=trk_t0_assn_v.at(pfps[0].key());
	  if(!t0s.size()) continue;
	}
	fNPassT0++;
      }




      if(fTrackModuleLabel=="pmtrack"){
	std::vector<art::Ptr<anab::T0>> T0s=fmT0.at(i);
	if(T0s.size()==0)
	  continue;
      }
      all_trks++;
      std::vector<art::Ptr<anab::Calorimetry>> calos=fmcal.at(i);
      const recob::Track& track = *ptrack;
      //TVector3 pos, dir_start, dir_end, end;
      auto pos = track.Vertex();
      auto dir_start = track.VertexDirection();
      auto dir_end   = track.EndDirection();
      auto end = track.End();
      double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
      double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());
      float startx=pos.X();
      float starty=pos.Y();
      float startz=pos.Z();
      float endx=end.X();
      float endy=end.Y();
      float endz=end.Z();

      float startcosx=dir_start.X();
      float startcosy=dir_start.Y();
      float startcosz=dir_start.Z();
      float endcosx=dir_end.X();
      float endcosy=dir_end.Y();
      float endcosz=dir_end.Z();
      float tracklength=track.Length();
      size_t count = 0;

      // anode-to-anode flag: track starts near one anode and ends near the other
      //double sdrift = (fDriftAxis==0?pos.X():(fDriftAxis==1?pos.Y():pos.Z()));
      //double edrift = (fDriftAxis==0?end.X():(fDriftAxis==1?end.Y():end.Z()));
      //double anodeTol = 20.;
            bool crossesCathode = (startx*endx>0);
      bool startAtLoAnode = (startx>350);
      bool startAtHiAnode = (startx<-350);
      bool endAtLoAnode   = (endx>350);
      bool endAtHiAnode   = (endx<-350);
      bool piercesAnode   = startAtLoAnode||startAtHiAnode||endAtLoAnode||endAtHiAnode;
      bool anodeToAnode   = (startAtLoAnode&&endAtHiAnode)||(startAtHiAnode&&endAtLoAnode);
      /********************************************************************************************************/

      /*********************storing hits for each tracks*******************************/

      std::vector<art::Ptr<recob::Hit>> allHits=fmtht.at(i); //storing hits for ith track
      std::vector<float> hitpeakT;
      for(size_t h1=0;h1<allHits.size();h1++){
	hitpeakT.push_back(allHits[h1]->PeakTime());
	//	if(allHits[h1]->WireID().Plane==2 && allHits[h1]->PeakTime()>420 && allHits[h1]->PeakTime()<530 && allHits[h1]->WireID().TPC==9 && allHits[h1]->WireID().Wire>60 && allHits[h1]->WireID().Wire<80) cout<<"wire id and peak time from allHits this is from track hit association "<<allHits[h1]->WireID().Wire<<"  "<<allHits[h1]->PeakTime()<<endl;
      }
      float max=-999999;
      float min=-999999;
      min=*min_element(hitpeakT.begin(),hitpeakT.end());
      max=*max_element(hitpeakT.begin(),hitpeakT.end());
      hitpeakT.clear();

	float dist0=99999;
	float dist=99999;
	float peaktime=-1;
	int wireno=99999;
	int counter1=99999;
	int tpcno=-1;

      if(((std::abs(pos.X())>330 ||  pos.Y()<50 || pos.Y()>550 || pos.Z()<50 || pos.Z()>645) && !(std::abs(end.X())>330 ||  end.Y()<50 || end.Y()>550 || end.Z()<50 || end.Z()>645))||(!(std::abs(pos.X())>330 ||  pos.Y()<50 || pos.Y()>550 || pos.Z()<50 || pos.Z()>645) && (std::abs(end.X())>330 ||  end.Y()<50 || end.Y()>550 || end.Z()<50 || end.Z()>645))){
	stopping_trks++; //these are total stopping tracks, could be broken as well

	/**********************************broken tracks removal**************************************************************************************/
	for(size_t k=0;k<NTracks;++k){
	  art::Ptr<recob::Track> ptrack_k(trackListHandle, k);
	  const recob::Track& track_k = *ptrack_k;
	  // TVector3 pos_k, dir_pos_k, dir_end_k, end_k;
	  auto pos_k = track_k.Vertex();
	  auto dir_pos_k= track_k.VertexDirection();
	  auto dir_end_k   = track_k.EndDirection();
	  auto end_k = track_k.End();
	  if(k==i) continue;
	  if((std::abs(((end_k.Y()-pos_k.Y())/(end_k.Z()-pos_k.Z()))*(endz-pos_k.Z())+pos_k.Y()-endy)<30||std::abs(((end_k.Y()-pos_k.Y())/(end_k.Z()-pos_k.Z()))*(startz-pos_k.Z())+pos_k.Y()-starty)<30)&&(std::abs(endcosx*dir_pos_k.X()+endcosy*dir_pos_k.Y()+endcosz*dir_pos_k.Z())>0.97||std::abs(startcosx*dir_pos_k.X()+startcosy*dir_pos_k.Y()+startcosz*dir_pos_k.Z())>0.97||std::abs(endcosx*dir_end_k.X()+endcosy*dir_end_k.Y()+endcosz*dir_end_k.Z())>0.97||std::abs(startcosx*dir_end_k.X()+startcosy*dir_end_k.Y()+startcosz*dir_end_k.Z())>0.97)) break;
	  if((std::abs(((end_k.Y()-pos_k.Y())/(end_k.Z()-pos_k.Z()))*(endz-pos_k.Z())+pos_k.Y()-endy)<50||std::abs(((end_k.Y()-pos_k.Y())/(end_k.Z()-pos_k.Z()))*(startz-pos_k.Z())+pos_k.Y()-starty)<50)&&(std::abs(endcosx*dir_pos_k.X()+endcosy*dir_pos_k.Y()+endcosz*dir_pos_k.Z())>0.998||std::abs(startcosx*dir_pos_k.X()+startcosy*dir_pos_k.Y()+startcosz*dir_pos_k.Z())>0.998||std::abs(endcosx*dir_end_k.X()+endcosy*dir_end_k.Y()+endcosz*dir_end_k.Z())>0.998||std::abs(startcosx*dir_end_k.X()+startcosy*dir_end_k.Y()+startcosz*dir_end_k.Z())>0.998)) break;
	  count++;
	}

	if(!(count==NTracks-1)|| tracklength<100) continue;
	unbroken_trks++;
	//*******************************new stuff*****************************//

	std::vector<int> wirenos;
	std::vector<float> peakts;



	int planenum1=999;
	float xpos=-9999;
	float ypos=-9999;
	float zpos=-9999;
	if(fmthm.isValid()){
	  auto vhit=fmthm.at(i);
	  auto vmeta=fmthm.data(i);
	  for (size_t ii = 0; ii<vhit.size(); ++ii){ //loop over all meta data hit
	    bool fBadhit = false;
	    if (vmeta[ii]->Index() == static_cast<unsigned int>(std::numeric_limits<int>::max())){
	      fBadhit = true;

	      continue;
	    }
	    if (vmeta[ii]->Index()>=tracklist[i]->NumberTrajectoryPoints()){
	      throw cet::exception("Calorimetry_module.cc") << "Requested track trajectory index "<<vmeta[ii]->Index()<<" exceeds the total number of trajectory points "<<tracklist[i]->NumberTrajectoryPoints()<<" for track index "<<i<<". Something is wrong with the track reconstruction. Please contact tjyang@fnal.gov!!";
	    }
	    if (!tracklist[i]->HasValidPoint(vmeta[ii]->Index())){
	      fBadhit = true;
	      continue;
	    }
	    // TVector3 loc = tracklist[i]->LocationAtPoint(vmeta[ii]->Index());
	     auto loc = tracklist[i]->LocationAtPoint(vmeta[ii]->Index());

	    xpos=loc.X();
	    ypos=loc.Y();
	    zpos=loc.Z();
	    if (fBadhit) continue; //HY::If BAD hit, skip this hit and go next
	    if (zpos<-100) continue; //hit not on track
	    planenum1=vhit[ii]->WireID().Plane;
	    if(planenum1==2){
	      if(starty<endy) dist=sqrt(pow(xpos-startx,2)+pow(ypos-starty,2)+pow(zpos-startz,2));
	      if(starty>endy) dist=sqrt(pow(xpos-endx,2)+pow(ypos-endy,2)+pow(zpos-endz,2));
	      //if(vhit[ii]->WireID().Plane==2 && vhit[ii]->PeakTime()>420 && vhit[ii]->PeakTime()<530 && vhit[ii]->WireID().TPC==9 && vhit[ii]->WireID().Wire>60 && vhit[ii]->WireID().Wire<80){
	      //cout<<"wire_no and peak time from metadata"<<vhit[ii]->WireID().Wire<<"  "<<vhit[ii]->PeakTime()<<endl;

	      wirenos.push_back(vhit[ii]->WireID().Wire);
	      peakts.push_back(vhit[ii]->PeakTime());
	      // }

	      if(dist<dist0){
		dist0=dist;
		wireno=vhit[ii]->WireID().Wire;
		peaktime=vhit[ii]->PeakTime();
		tpcno=vhit[ii]->WireID().TPC;
	      }
	    }
	  }//loop over vhit
	}//fmthm valid
	/*************************************filling the values****************************/
	counter1=0;
	// bool test=false;
	for(size_t hitl=0;hitl<hitlist.size();hitl++){
	  auto & tracks = thass.at(hitlist[hitl].key());
	  //auto & vmeta = thass.data(hitlist[hitl].key());
	  //if (!tracks.empty()&&tracks[0].key() == ptrack.key()&&vmeta[0]->Index()!=static_cast<unsigned int>(std::numeric_limits<int>::max())) continue;
	  //if (!tracks.empty()&&tracks[0].key() == ptrack.key()&&hitlist[hitl]->WireID().Plane==2)
	  //if(!(hitlist[hitl]->WireID().Plane==2 && tracks[0].key() != ptrack.key() && hitlist[hitl]->PeakTime()>420 && hitlist[hitl]->PeakTime()<530 && hitlist[hitl]->WireID().TPC==9 && hitlist[hitl]->WireID().Wire>60 && hitlist[hitl]->WireID().Wire<80)) continue;
	  //cout<<"wire no and peak time from hit to track association...................."<<hitlist[hitl]->WireID().Wire<<"  "<<hitlist[hitl]->PeakTime()<<endl;
	  if (!tracks.empty() && tracks[0].key()!=ptrack.key() && tracklist[tracks[0].key()]->Length()>100) continue;
	  bool test=true;
	  float peakth1=hitlist[hitl]->PeakTime();
	  int wireh1=hitlist[hitl]->WireID().Wire;
	  for(size_t m=0;m<wirenos.size();m++){
	    if(wireh1==wirenos[m] && peakth1==peakts[m]){
	      test=false;
	      break;
	    }
	  }
	  if(!test) continue;

	  int planeid=hitlist[hitl]->WireID().Plane;
	  int tpcid=hitlist[hitl]->WireID().TPC;
	  if(abs(wireh1-wireno)<6 && abs(peakth1-peaktime)<50 && planeid==2 && tpcid==tpcno){
	    counter1++;
	  }
	}
	// for(size_t t=0;t<peakts.size();t++) cout<<"wire numbers and peak times in the vector "<<wirenos[t]<<"  "<<peakts[t]<<endl;

	wirenos.clear();
	peakts.clear();


        cout<<"no of hits closeby  "<<counter1<<"   "<<"event "<<event<<" TrkackID "<<track.ID()<<" startx, y, z "<<startx<<" "<<starty<<" "<<startz<<"  wireno, peakt tpcno "<<wireno<<" "<<peaktime<<" "<<tpcno<<" dist "<<dist0<<"min T, max_T"<<min<<" "<<max<<endl;
      }


      cross_trks++; // these are good tracks unbroken
      fNWritten++;
      dist_min[cross_trks-1]=dist0;
      adjacent_hits[cross_trks-1]=counter1;
      lastwire[cross_trks-1]=wireno;
      lastpeakt[cross_trks-1]=peaktime;
      endtpc[cross_trks-1]=tpcno;
      trkanodetoanode[cross_trks-1]=anodeToAnode?1:0;



      trackthetaxz[cross_trks-1]=theta_xz;
      trackthetayz[cross_trks-1]=theta_yz;
      trkstartx[cross_trks-1]=pos.X();
      trkstarty[cross_trks-1]=pos.Y();
      trkstartz[cross_trks-1]=pos.Z();
      trkendx[cross_trks-1]=end.X();
      trkendy[cross_trks-1]=end.Y();
      trkendz[cross_trks-1]=end.Z();
      trklen[cross_trks-1]=track.Length();
      TrkID[cross_trks-1]=track.ID();
      peakT_max[cross_trks-1]=max;
      peakT_min[cross_trks-1]=min;
      trkstartcosxyz[cross_trks-1][0]=dir_start.X();
      trkstartcosxyz[cross_trks-1][1]=dir_start.Y();
      trkstartcosxyz[cross_trks-1][2]=dir_start.Z();
      trkendcosxyz[cross_trks-1][0]=dir_end.X();
      trkendcosxyz[cross_trks-1][1]=dir_end.Y();
      trkendcosxyz[cross_trks-1][2]=dir_end.Z();
      for(size_t ical = 0; ical<calos.size(); ++ical){
	if(!calos[ical]) continue;
	if(!calos[ical]->PlaneID().isValid) continue;
	int planenum = calos[ical]->PlaneID().Plane;
	if(planenum<0||planenum>2) continue;
	const size_t NHits = calos[ical] -> dEdx().size();
	ntrkhits[cross_trks-1][planenum]=int(NHits);
	const std::vector<size_t>& hitIndices = calos[ical]->TpIndices();
	for(size_t iHit = 0; iHit < NHits; ++iHit){
	  if(iHit>=3000) break;
	  const auto& TrkPos = (calos[ical] -> XYZ())[iHit];
	  trkdqdx[cross_trks-1][planenum][iHit]=(calos[ical] -> dQdx())[iHit];
	  trkdedx[cross_trks-1][planenum][iHit]=(calos[ical] -> dEdx())[iHit];
	  trkresrange[cross_trks-1][planenum][iHit]=(calos[ical]->ResidualRange())[iHit];
	  trkhitx[cross_trks-1][planenum][iHit]=TrkPos.X();
	  trkhity[cross_trks-1][planenum][iHit]=TrkPos.Y();
	  trkhitz[cross_trks-1][planenum][iHit]=TrkPos.Z();
	  trkpitch[cross_trks-1][planenum][iHit]=(calos[ical]->TrkPitchVec())[iHit];
	  if(iHit<hitIndices.size() && hitIndices[iHit]<hitlist.size()){
	    const art::Ptr<recob::Hit>& hit = hitlist[hitIndices[iHit]];
	    trkhitpeakt[cross_trks-1][planenum][iHit]=hit->PeakTime();
	    trkhittpc[cross_trks-1][planenum][iHit]=hit->WireID().TPC;
	    trkhitwire[cross_trks-1][planenum][iHit]=hit->WireID().Wire;
	    trkhitrms[cross_trks-1][planenum][iHit]=hit->RMS();
	    trkhitdeltat[cross_trks-1][planenum][iHit]=hit->EndTick()-hit->StartTick();
	  }
	} // loop over iHit..
      } // loop over ical 2nd time...
      if (!isData){
        auto particles = truthUtil.GetMCParticleListFromRecoTrack(clockData, *tracklist[i], evt, fTrackModuleLabel);
        for (auto & part : particles){
          //cout<<(part.first)->PdgCode()<<" "<<part.second<<endl;
          if (!part.first) continue;
          true_trkstartx[cross_trks-1] = (part.first)->Vx();
          true_trkstarty[cross_trks-1] = (part.first)->Vy();
          true_trkstartz[cross_trks-1] = (part.first)->Vz();
          true_trkendx[cross_trks-1] = (part.first)->EndX();
          true_trkendy[cross_trks-1] = (part.first)->EndY();
          true_trkendz[cross_trks-1] = (part.first)->EndZ();
          break;
        }
      }
    } // loop over trks...

    fEventTree->Fill();
  } // end of analyze function

    /////////////////// Defintion of reset function ///////////
  void michelremovingpdhd::reset(){
    isData = false;
    run = -99999;
    subrun = -99999;
    event = -99999;
    evttime = -99999;
    cross_trks = -99999;
    //all_trks = -99999;
    //unbroken_trks=-99999;
    year_month_date=-99999;
    hour_min_sec=-99999;
    for(int i=0; i<kMaxTracks; i++){
      trackthetaxz[i]=-99999;
      trackthetayz[i]=-99999;
      trkstartx[i]=-99999;
      trkstarty[i]=-99999;
      trkstartz[i]=-99999;
      trkendx[i]=-99999;
      trkendy[i]=-99999;
      trkendz[i]=-99999;
      trklen[i]=-99999;
      TrkID[i]=-99999;

      dist_min[i]=-1;
      adjacent_hits[i]=-1;
      lastwire[i]=-1;
      lastpeakt[i]=-1;
      endtpc[i]=-1;
      trkanodetoanode[i]=-1;
      peakT_max[i]=-99999;
      peakT_min[i]=-99999;
      trkstartcosxyz[i][0]=-99999;
      trkstartcosxyz[i][1]=-99999;
      trkstartcosxyz[i][2]=-99999;
      trkendcosxyz[i][0]=-99999;
      trkendcosxyz[i][1]=-99999;
      trkendcosxyz[i][2]=-99999;
      ntrkhits[i][0] = -99999;
      ntrkhits[i][1] = -99999;
      ntrkhits[i][2] = -99999;
      true_trkstartx[i]=-99999;
      true_trkstarty[i]=-99999;
      true_trkstartz[i]=-99999;
      true_trkendx[i]=-99999;
      true_trkendy[i]=-99999;
      true_trkendz[i]=-99999;
      for(int j=0; j<3; j++){
	for(int k=0; k<3000; k++){
	  trkdqdx[i][j][k]=-99999;
	  trkdedx[i][j][k]=-99999;
	  trkresrange[i][j][k]=-99999;
	  trkhitx[i][j][k]=-99999;
	  trkhity[i][j][k]=-99999;
	  trkhitz[i][j][k]=-99999;
	  trkpitch[i][j][k]=-99999;
	  trkhitpeakt[i][j][k]=-99999;
	  trkhittpc[i][j][k]=-99999;
	  trkhitwire[i][j][k]=-99999;
	  trkhitrms[i][j][k]=-99999;
	  trkhitdeltat[i][j][k]=-99999;
	}
      }
    }
  }
  //////////////////////// End of definition ///////////////

  DEFINE_ART_MODULE(michelremovingpdhd)
}

