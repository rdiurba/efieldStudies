////////////////////////////////////////////////////////////////////////
// Class:       laseranalysispdhd
// File:        laseranalysispdhd_module.cc
//
// Laser-track analyzer for ProtoDUNE Horizontal Drift, Pandora reconstruction.
//
// PDHD copy of laseranalysispdvd (identical code and tree layout). Like that
// module it is based on the michelremovingpdvd analyzer: it writes the SAME tree/branch
// layout ("Event" tree) so laser trees can be read by the existing
// michelremoving MakeClass headers and analysis macros. The difference is the
// track selection -- the michel-specific cuts (T0 tag requirement and the
// stopping/unbroken-track matching) are removed, so EVERY reconstructed Pandora
// track is saved. The single laser track (and the cosmic tracks in the same
// event) can then be selected offline by start position / direction.
////////////////////////////////////////////////////////////////////////

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
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
#include "nusimdata/SimulationBase/MCParticle.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardata/DetectorInfoServices/DetectorClocksService.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/TrackHitMeta.h"
#include "lardataobj/RecoBase/Hit.h"
#include "larcoreobj/SimpleTypesAndConstants/geo_types.h"
#include "lardataobj/AnalysisBase/T0.h"
#include "lardataobj/RawData/RDTimeStamp.h"
#include "cetlib_except/exception.h"

#include "protoduneana/Utilities/ProtoDUNETruthUtils.h"

#include "TTree.h"
#include "TMath.h"
#include "TTimeStamp.h"

#include <vector>
#include <algorithm>
#include <string>
#include <limits>
#include <cmath>

const int kMaxTracks = 1000;

using namespace std;

namespace pdhdana {

  class laseranalysispdhd : public art::EDAnalyzer {
  public:
    explicit laseranalysispdhd(fhicl::ParameterSet const& pset);
    virtual ~laseranalysispdhd();

    void beginJob();
    void beginRun(const art::Run& run);
    void endJob();
    void analyze(const art::Event& evt);
    void reset();

  private:
    void bookBranches(TTree* t);   // book the shared branches on a given tree

    TTree* fEventTree;             // "Event"
    bool     isData;
    Int_t    run;
    Int_t    subrun;
    Int_t    event;
    Double_t evttime;
    Int_t    year_month_date;
    Int_t    hour_min_sec;
    ULong64_t triggertime;    // raw art event timestamp evt.time().value(), unconverted
    ULong64_t rdtimestamp;    // raw::RDTimeStamp::GetTimeStamp() of the first RDTimeStamp
                              // in the event, unconverted (0 if the product is absent).
                              // The laser CSV TIMESTAMP column is DTS 62.5 MHz ticks since
                              // the epoch; match on whichever of these two is in the same
                              // units (~1.08e17 for Dec 2024). CSV rows are 6250000 ticks
                              // (0.1 s) apart.
    Int_t    cross_trks;      // number of saved tracks (array length, kept for format)
    Int_t    stopping_trks;
    Int_t    all_trks;
    Int_t    unbroken_trks;
    Float_t  trackthetaxz[kMaxTracks];
    Float_t  trackthetayz[kMaxTracks];
    Float_t  trkstartx[kMaxTracks];
    Float_t  trkstarty[kMaxTracks];
    Float_t  trkstartz[kMaxTracks];
    Float_t  trkendx[kMaxTracks];
    Float_t  trkendy[kMaxTracks];
    Float_t  trkendz[kMaxTracks];
    Float_t  trklen[kMaxTracks];
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

    Float_t  peakT_max[kMaxTracks];
    Float_t  peakT_min[kMaxTracks];
    Float_t  dist_min[kMaxTracks];
    Int_t    adjacent_hits[kMaxTracks];
    Int_t    lastwire[kMaxTracks];
    Int_t    endtpc[kMaxTracks];
    Float_t  lastpeakt[kMaxTracks];
    Int_t    trkcathodecrosser[kMaxTracks];
    Int_t    trkanodepiercer[kMaxTracks];
    Int_t    trkanodetoanode[kMaxTracks];
    Int_t    trkflipped[kMaxTracks];      // 1 = start/end swapped so the track starts at the periscope
    float    true_trkstartx[kMaxTracks];
    float    true_trkstarty[kMaxTracks];
    float    true_trkstartz[kMaxTracks];
    float    true_trkendx[kMaxTracks];
    float    true_trkendy[kMaxTracks];
    float    true_trkendz[kMaxTracks];

    std::string fHitsModuleLabel;
    std::string fTrackModuleLabel;
    std::string fCalorimetryModuleLabel;
    art::InputTag fRDTimeStampTag;  // label of the std::vector<raw::RDTimeStamp> product
    double      fMinTrackLength;  // save tracks longer than this [cm] (0 = save all)
    // "bubble" selection, one bubble per laser (PDHD has two; they never fire
    // at the same time). Per event: try laser 0; if no track passes within
    // fBubbleRadius of (fLaserX[0],fLaserY[0],fLaserZ[0]), try laser 1, and so
    // on. Only tracks in the first laser's bubble that has any are saved.
    // Distance is 3D (X,Y,Z); if LaserX is left empty it is Y-Z only.
    // fBubbleRadius <= 0 disables the cut (save all tracks).
    std::vector<double> fLaserX;
    std::vector<double> fLaserY;
    std::vector<double> fLaserZ;
    double      fBubbleRadius;
    Int_t       laser_id;         // which laser's bubble was used (-1 = none / cut off)

    // PDHD active-volume bounding box, from the geometry service (beginRun).
    // Written to the tree so the laser ray can be clipped to the real detector
    // boundaries offline (see add_laser_truth.C).
    double      fActiveMin[3];
    double      fActiveMax[3];
    Float_t     det_min[3];       // {xmin,ymin,zmin} branch
    Float_t     det_max[3];       // {xmax,ymax,zmax} branch

    protoana::ProtoDUNETruthUtils truthUtil;
  };

  laseranalysispdhd::laseranalysispdhd(fhicl::ParameterSet const& pset) :
    EDAnalyzer(pset),
    fHitsModuleLabel        (pset.get< std::string >("HitsModuleLabel","")        ),
    fTrackModuleLabel       (pset.get< std::string >("TrackModuleLabel","")       ),
    fCalorimetryModuleLabel (pset.get< std::string >("CalorimetryModuleLabel","") ),
    fRDTimeStampTag         (pset.get< art::InputTag >("RDTimeStampTag", art::InputTag("tpcrawdecoder","daq")) ),
    fMinTrackLength         (pset.get< double      >("MinTrackLength", 0.0)       ),
    fLaserX                 (pset.get< std::vector<double> >("LaserX", {})     ),
    fLaserY                 (pset.get< std::vector<double> >("LaserY", {})     ),
    fLaserZ                 (pset.get< std::vector<double> >("LaserZ", {})     ),
    fBubbleRadius           (pset.get< double      >("BubbleRadius", -1.0)        )
  {
    if (fLaserY.size() != fLaserZ.size())
      throw cet::exception("laseranalysispdhd") << "LaserY and LaserZ must have the same length\n";
    if (!fLaserX.empty() && fLaserX.size() != fLaserY.size())
      throw cet::exception("laseranalysispdhd") << "LaserX must be empty or the same length as LaserY/LaserZ\n";
    for(int a=0;a<3;a++){ fActiveMin[a]=0.; fActiveMax[a]=0.; det_min[a]=0.; det_max[a]=0.; }
  }

  laseranalysispdhd::~laseranalysispdhd(){
  }

  // book the tree branches.
  void laseranalysispdhd::bookBranches(TTree* t){
    t->Branch("isData", &isData, "isData/O");
    t->Branch("event", &event,"event/I");
    t->Branch("evttime",&evttime,"evttime/D");
    t->Branch("run", &run,"run/I");
    t->Branch("subrun", &subrun,"surbrun/I");
    t->Branch("year_month_date", &year_month_date,"year_month_date/I");
    t->Branch("hour_min_sec", &hour_min_sec,"hour_min_sec/I");
    t->Branch("rdtimestamp", &rdtimestamp,"rdtimestamp/l");
    t->Branch("laser_id", &laser_id,"laser_id/I");
    t->Branch("det_min",det_min,"det_min[3]/F");   // PDHD active box {xmin,ymin,zmin}
    t->Branch("det_max",det_max,"det_max[3]/F");   // PDHD active box {xmax,ymax,zmax}
    t->Branch("cross_trks",&cross_trks,"cross_trks/I");
    t->Branch("stopping_trks",&stopping_trks,"stopping_trks/I");
    t->Branch("all_trks",&all_trks,"all_trks/I");
    t->Branch("unbroken_trks",&unbroken_trks,"unbroken_trks/I");
    t->Branch("trackthetaxz",trackthetaxz,"trackthetaxz[cross_trks]/F");
    t->Branch("trackthetayz",trackthetayz,"trackthetayz[cross_trks]/F");
    t->Branch("trkstartx",trkstartx,"trkstartx[cross_trks]/F");
    t->Branch("trkstarty",trkstarty,"trkstarty[cross_trks]/F");
    t->Branch("trkstartz",trkstartz,"trkstartz[cross_trks]/F");
    t->Branch("trkendx",trkendx,"trkendx[cross_trks]/F");
    t->Branch("trkendy",trkendy,"trkendy[cross_trks]/F");
    t->Branch("trkendz",trkendz,"trkendz[cross_trks]/F");
    t->Branch("trklen",trklen,"trklen[cross_trks]/F");
    t->Branch("peakT_max",peakT_max,"peakT_max[cross_trks]/F");
    t->Branch("peakT_min",peakT_min,"peakT_min[cross_trks]/F");
    t->Branch("TrkID",TrkID,"TrkID[cross_trks]/I");
    t->Branch("trkstartcosxyz",trkstartcosxyz,"trkstartcosxyz[cross_trks][3]/F");
    t->Branch("trkendcosxyz",trkendcosxyz,"trkendcosxyz[cross_trks][3]/F");
    t->Branch("ntrkhits",ntrkhits,"ntrkhits[cross_trks][3]/I");
    t->Branch("trkdqdx",trkdqdx,"trkdqdx[cross_trks][3][3000]/F");
    t->Branch("trkdedx",trkdedx,"trkdedx[cross_trks][3][3000]/F");
    t->Branch("trkresrange",trkresrange,"trkresrange[cross_trks][3][3000]/F");
    t->Branch("trkhitx",trkhitx,"trkhitx[cross_trks][3][3000]/F");
    t->Branch("trkhity",trkhity,"trkhity[cross_trks][3][3000]/F");
    t->Branch("trkhitz",trkhitz,"trkhitz[cross_trks][3][3000]/F");
    t->Branch("trkpitch",trkpitch,"trkpitch[cross_trks][3][3000]/F");
    t->Branch("trkhitpeakt",trkhitpeakt,"trkhitpeakt[cross_trks][3][3000]/F");
    t->Branch("trkhittpc",trkhittpc,"trkhittpc[cross_trks][3][3000]/I");
    t->Branch("trkhitwire",trkhitwire,"trkhitwire[cross_trks][3][3000]/I");
    t->Branch("trkhitrms",trkhitrms,"trkhitrms[cross_trks][3][3000]/F");
    t->Branch("trkhitdeltat",trkhitdeltat,"trkhitdeltat[cross_trks][3][3000]/F");
    t->Branch("dist_min",dist_min,"dist_min[cross_trks]/F");
    t->Branch("adjacent_hits",adjacent_hits,"adjacent_hits[cross_trks]/I");
    t->Branch("lastwire",lastwire,"lastwire[cross_trks]/I");
    t->Branch("lastpeakt",lastpeakt,"lastpeakt[cross_trks]/F");
    t->Branch("endtpc",endtpc,"endtpc[cross_trks]/I");
    t->Branch("trkcathodecrosser",trkcathodecrosser,"trkcathodecrosser[cross_trks]/I");
    t->Branch("trkanodepiercer",trkanodepiercer,"trkanodepiercer[cross_trks]/I");
    t->Branch("trkanodetoanode",trkanodetoanode,"trkanodetoanode[cross_trks]/I");
    t->Branch("trkflipped",trkflipped,"trkflipped[cross_trks]/I");
    t->Branch("true_trkstartx",true_trkstartx,"true_trkstartx[cross_trks]/F");
    t->Branch("true_trkstarty",true_trkstarty,"true_trkstarty[cross_trks]/F");
    t->Branch("true_trkstartz",true_trkstartz,"true_trkstartz[cross_trks]/F");
    t->Branch("true_trkendx",true_trkendx,"true_trkendx[cross_trks]/F");
    t->Branch("true_trkendy",true_trkendy,"true_trkendy[cross_trks]/F");
    t->Branch("true_trkendz",true_trkendz,"true_trkendz[cross_trks]/F");
  }

  void laseranalysispdhd::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;

    fEventTree = tfs->make<TTree>("Event", "Event Tree from Reco");
    bookBranches(fEventTree);
    fEventTree->Branch("triggertime",&triggertime,"triggertime/l");
  }

  void laseranalysispdhd::beginRun(const art::Run&){
    // PDHD active-volume box = union of all TPC bounding boxes.
    auto const& geom = *art::ServiceHandle<geo::Geometry>();
    fActiveMin[0]=fActiveMin[1]=fActiveMin[2]= 1e12;
    fActiveMax[0]=fActiveMax[1]=fActiveMax[2]=-1e12;
    for (auto const& tpc : geom.Iterate<geo::TPCGeo>()){
      auto const box = tpc.BoundingBox();
      fActiveMin[0]=std::min(fActiveMin[0],box.MinX()); fActiveMax[0]=std::max(fActiveMax[0],box.MaxX());
      fActiveMin[1]=std::min(fActiveMin[1],box.MinY()); fActiveMax[1]=std::max(fActiveMax[1],box.MaxY());
      fActiveMin[2]=std::min(fActiveMin[2],box.MinZ()); fActiveMax[2]=std::max(fActiveMax[2],box.MaxZ());
    }
    for(int a=0;a<3;a++){ det_min[a]=fActiveMin[a]; det_max[a]=fActiveMax[a]; }
  }

  void laseranalysispdhd::endJob(){
  }

  void laseranalysispdhd::analyze(const art::Event& evt){
    reset();

    std::vector<art::Ptr<recob::Track> > tracklist;
    auto trackListHandle = evt.getHandle< std::vector<recob::Track> >(fTrackModuleLabel);
    if (trackListHandle) art::fill_ptr_vector(tracklist, trackListHandle);
    else return;

    std::vector<art::Ptr<recob::Hit>> hitlist;
    auto hitListHandle = evt.getHandle< std::vector<recob::Hit> >(fHitsModuleLabel);
    if (hitListHandle) art::fill_ptr_vector(hitlist, hitListHandle);

    art::FindManyP<recob::Hit> fmtht(trackListHandle, evt, fTrackModuleLabel);
    art::FindManyP<anab::Calorimetry> fmcal(trackListHandle, evt, fCalorimetryModuleLabel);
    if (!fmcal.isValid()) return;

    auto const clockData = art::ServiceHandle<detinfo::DetectorClocksService>()->DataFor(evt);

    isData = evt.isRealData();
    run    = evt.run();
    subrun = evt.subRun();
    event  = evt.id().event();
    art::Timestamp ts = evt.time();
    triggertime = ts.value();                 // raw 64-bit art timestamp
    TTimeStamp tts(ts.timeHigh(), ts.timeLow());
    evttime = tts.AsDouble();

    // RDTimeStamp: take the first entry of the RDTimeStamp collection
    auto rdtsHandle = evt.getHandle< std::vector<raw::RDTimeStamp> >(fRDTimeStampTag);
    if (rdtsHandle && !rdtsHandle->empty()) rdtimestamp = rdtsHandle->front().GetTimeStamp();

    UInt_t year=0, month=0, day=0;
    year_month_date = tts.GetDate(kTRUE,0,&year,&month,&day);
    UInt_t hour=0, min=0, sec=0;
    hour_min_sec = tts.GetTime(kTRUE,0,&hour,&min,&sec);

    cross_trks    = 0;
    stopping_trks = 0;
    all_trks      = 0;
    unbroken_trks = 0;

    size_t NTracks = tracklist.size();

    // smallest distance between a track's trajectory points and laser k's bubble
    // centre: 3D if LaserX is given, else Y-Z only
    auto dminLaser = [&](const recob::Track& trk, size_t k){
      const bool useX = !fLaserX.empty();
      double dmin = 1e12;
      for (size_t p = 0; p < trk.NumberTrajectoryPoints(); ++p){
        if (!trk.HasValidPoint(p)) continue;
        const auto loc = trk.LocationAtPoint(p);
        const double dX = useX ? loc.X()-fLaserX[k] : 0.;
        const double dY = loc.Y()-fLaserY[k], dZ = loc.Z()-fLaserZ[k];
        dmin = std::min(dmin, std::sqrt(dX*dX + dY*dY + dZ*dZ));
      }
      return dmin;
    };

    // pick the laser: the first one (fcl order) with any track in its bubble.
    // If none has a track, laser_id stays -1 and no tracks are saved.
    if (fBubbleRadius > 0.){
      for (size_t k = 0; k < fLaserY.size() && laser_id < 0; ++k){
        for (size_t i = 0; i < NTracks; ++i){
          const recob::Track& trk = *tracklist[i];
          if (trk.Length() < fMinTrackLength) continue;
          if (dminLaser(trk, k) <= fBubbleRadius){ laser_id = int(k); break; }
        }
      }
    }

    for(size_t i=0; i<NTracks; ++i){
      if(cross_trks>=kMaxTracks) break;
      art::Ptr<recob::Track> ptrack(trackListHandle, i);
      all_trks++;

      const recob::Track& track = *ptrack;
      std::vector<art::Ptr<anab::Calorimetry>> calos = fmcal.at(i);

      auto pos       = track.Vertex();
      auto dir_start = track.VertexDirection();
      auto dir_end   = track.EndDirection();
      auto end       = track.End();

      // Orientation: the laser track must START at the periscope. If an
      // endpoint is within fBubbleRadius of a laser position (LaserX/Y/Z; 3D if
      // LaserX is given, else Y-Z) and that endpoint is the track END (it is
      // closer to a laser than the start is), swap start and end: positions,
      // directions (reversed), and the order of the calorimetry points below.
      auto dLaser = [&](double x, double y, double z){
        const bool useX = !fLaserX.empty();
        double dmin = 1e12;
        for (size_t k = 0; k < fLaserY.size(); ++k){
          const double dX = useX ? x-fLaserX[k] : 0.;
          const double dY = y-fLaserY[k], dZ = z-fLaserZ[k];
          dmin = std::min(dmin, std::sqrt(dX*dX + dY*dY + dZ*dZ));
        }
        return dmin;
      };
      bool flipped = false;
      if (fBubbleRadius > 0. && !fLaserY.empty()){
        const double dS = dLaser(pos.X(), pos.Y(), pos.Z());
        const double dE = dLaser(end.X(), end.Y(), end.Z());
        if (dE <= fBubbleRadius && dE < dS){
          std::swap(pos, end);
          const auto ds = dir_start;
          dir_start = -dir_end;
          dir_end   = -ds;
          flipped = true;
        }
      }
      double theta_xz = std::atan2(dir_start.X(), dir_start.Z());
      double theta_yz = std::atan2(dir_start.Y(), dir_start.Z());
      float startx = pos.X(), starty = pos.Y(), startz = pos.Z();
      float endx   = end.X(), endy   = end.Y(), endz   = end.Z();
      float tracklength = track.Length();

      // save every track above the (configurable) length. Default
      // MinTrackLength=0 keeps all tracks.
      if (tracklength < fMinTrackLength) continue;

      // YZ bubble selection around the laser chosen above. Disabled if radius<=0.
      if (fBubbleRadius > 0.){
        if (laser_id < 0) continue;                                        // no laser found
        if (dminLaser(track, laser_id) > fBubbleRadius) continue;
      }

      // geometry flags, same convention/thresholds as laseranalysispdvd
      bool crossesCathode = (startx*endx>0);
      bool startAtLoAnode = (startx>330);
      bool startAtHiAnode = (startx<-330);
      bool endAtLoAnode   = (endx>330);
      bool endAtHiAnode   = (endx<-330);
      bool piercesAnode   = startAtLoAnode||startAtHiAnode||endAtLoAnode||endAtHiAnode;
      bool anodeToAnode   = (startAtLoAnode&&endAtHiAnode)||(startAtHiAnode&&endAtLoAnode);

      // hit peak-time span for this track
      std::vector<art::Ptr<recob::Hit>> allHits = fmtht.at(i);
      float maxpt = -999999, minpt = -999999;
      if(!allHits.empty()){
        std::vector<float> hitpeakT;
        hitpeakT.reserve(allHits.size());
        for(size_t h1=0; h1<allHits.size(); h1++) hitpeakT.push_back(allHits[h1]->PeakTime());
        minpt = *min_element(hitpeakT.begin(),hitpeakT.end());
        maxpt = *max_element(hitpeakT.begin(),hitpeakT.end());
      }

      cross_trks++;
      trackthetaxz[cross_trks-1]=theta_xz;
      trackthetayz[cross_trks-1]=theta_yz;
      trkstartx[cross_trks-1]=startx;
      trkstarty[cross_trks-1]=starty;
      trkstartz[cross_trks-1]=startz;
      trkendx[cross_trks-1]=endx;
      trkendy[cross_trks-1]=endy;
      trkendz[cross_trks-1]=endz;
      trklen[cross_trks-1]=tracklength;
      TrkID[cross_trks-1]=track.ID();
      peakT_max[cross_trks-1]=maxpt;
      peakT_min[cross_trks-1]=minpt;
      trkstartcosxyz[cross_trks-1][0]=dir_start.X();
      trkstartcosxyz[cross_trks-1][1]=dir_start.Y();
      trkstartcosxyz[cross_trks-1][2]=dir_start.Z();
      trkendcosxyz[cross_trks-1][0]=dir_end.X();
      trkendcosxyz[cross_trks-1][1]=dir_end.Y();
      trkendcosxyz[cross_trks-1][2]=dir_end.Z();
      trkcathodecrosser[cross_trks-1]=crossesCathode?1:0;
      trkanodepiercer[cross_trks-1]=piercesAnode?1:0;
      trkanodetoanode[cross_trks-1]=anodeToAnode?1:0;
      trkflipped[cross_trks-1]=flipped?1:0;

      for(size_t ical = 0; ical<calos.size(); ++ical){
        if(!calos[ical]) continue;
        if(!calos[ical]->PlaneID().isValid) continue;
        int planenum = calos[ical]->PlaneID().Plane;
        if(planenum<0||planenum>2) continue;
        const size_t NHits = calos[ical]->dEdx().size();
        ntrkhits[cross_trks-1][planenum]=int(NHits);
        const std::vector<size_t>& hitIndices = calos[ical]->TpIndices();
        for(size_t iHit = 0; iHit < NHits; ++iHit){
          if(iHit>=3000) break;
          if(hitIndices[iHit]>=hitlist.size()) continue;
          const auto& TrkPos = (calos[ical]->XYZ())[iHit];
          const art::Ptr<recob::Hit>& hit = hitlist[hitIndices[iHit]];
          trkhittpc[cross_trks-1][planenum][iHit]=hit->WireID().TPC;
          trkhitwire[cross_trks-1][planenum][iHit]=hit->WireID().Wire;
          trkhitrms[cross_trks-1][planenum][iHit]=hit->RMS();
          trkhitdeltat[cross_trks-1][planenum][iHit]=hit->EndTick()-hit->StartTick();
          trkhitpeakt[cross_trks-1][planenum][iHit]=hit->PeakTime();
          trkdqdx[cross_trks-1][planenum][iHit]=(calos[ical]->dQdx())[iHit];
          trkdedx[cross_trks-1][planenum][iHit]=(calos[ical]->dEdx())[iHit];
          trkresrange[cross_trks-1][planenum][iHit]=(calos[ical]->ResidualRange())[iHit];
          trkhitx[cross_trks-1][planenum][iHit]=TrkPos.X();
          trkhity[cross_trks-1][planenum][iHit]=TrkPos.Y();
          trkhitz[cross_trks-1][planenum][iHit]=TrkPos.Z();
          trkpitch[cross_trks-1][planenum][iHit]=(calos[ical]->TrkPitchVec())[iHit];
        }
      }

      // flipped track: reverse the calorimetry points of every plane so they
      // run from the new start (periscope) to the new end
      if (flipped){
        const int it = cross_trks-1;
        for (int pl = 0; pl < 3; ++pl){
          const int n = std::min(ntrkhits[it][pl], 3000);
          if (n <= 1) continue;
          std::reverse(trkdqdx[it][pl],     trkdqdx[it][pl]+n);
          std::reverse(trkdedx[it][pl],     trkdedx[it][pl]+n);
          std::reverse(trkresrange[it][pl], trkresrange[it][pl]+n);
          std::reverse(trkhitx[it][pl],     trkhitx[it][pl]+n);
          std::reverse(trkhity[it][pl],     trkhity[it][pl]+n);
          std::reverse(trkhitz[it][pl],     trkhitz[it][pl]+n);
          std::reverse(trkpitch[it][pl],    trkpitch[it][pl]+n);
          std::reverse(trkhitpeakt[it][pl], trkhitpeakt[it][pl]+n);
          std::reverse(trkhittpc[it][pl],   trkhittpc[it][pl]+n);
          std::reverse(trkhitwire[it][pl],  trkhitwire[it][pl]+n);
          std::reverse(trkhitrms[it][pl],   trkhitrms[it][pl]+n);
          std::reverse(trkhitdeltat[it][pl],trkhitdeltat[it][pl]+n);
        }
      }

      if (!isData){
        auto particles = truthUtil.GetMCParticleListFromRecoTrack(clockData, *tracklist[i], evt, fTrackModuleLabel);
        for (auto & part : particles){
          if (!part.first) continue;
          true_trkstartx[cross_trks-1] = (part.first)->Vx();
          true_trkstarty[cross_trks-1] = (part.first)->Vy();
          true_trkstartz[cross_trks-1] = (part.first)->Vz();
          true_trkendx[cross_trks-1]   = (part.first)->EndX();
          true_trkendy[cross_trks-1]   = (part.first)->EndY();
          true_trkendz[cross_trks-1]   = (part.first)->EndZ();
          break;
        }
      }
    }

    fEventTree->Fill();
  }

  void laseranalysispdhd::reset(){
    isData = false;
    run = -99999;
    subrun = -99999;
    event = -99999;
    evttime = -99999;
    triggertime = 0;
    rdtimestamp = 0;
    laser_id = -1;
    cross_trks = -99999;
    stopping_trks = -99999;
    all_trks = -99999;
    unbroken_trks = -99999;
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
      trkcathodecrosser[i]=-1;
      trkanodepiercer[i]=-1;
      trkanodetoanode[i]=-1;
      trkflipped[i]=-1;
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

  DEFINE_ART_MODULE(laseranalysispdhd)
}
