//////////////////////////////////////////////////////////
// Reader for the michelremovingpdhd Event tree.
//////////////////////////////////////////////////////////

#ifndef data_crossing_pdhd_h
#define data_crossing_pdhd_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TString.h>

#include <fstream>
#include <string>
#include <vector>

class data_crossing_pdhd {
public :
   TTree          *fChain;
   Int_t           fCurrent;

   Int_t           event;
   Double_t        evttime;
   Int_t           run;
   Int_t           subrun;
   Int_t           year_month_date;
   Int_t           hour_min_sec;
   Int_t           cross_trks;
   Int_t           stopping_trks;
   Int_t           all_trks;
   Int_t           unbroken_trks;
   Float_t         trackthetaxz[1000];
   Float_t         trackthetayz[1000];
   Float_t         trkstartx[1000];
   Float_t         trkstarty[1000];
   Float_t         trkstartz[1000];
   Float_t         trkendx[1000];
   Float_t         trkendy[1000];
   Float_t         trkendz[1000];
   Float_t         trklen[1000];
   Float_t         peakT_max[1000];
   Float_t         peakT_min[1000];
   Int_t           TrkID[1000];
   Float_t         trkstartcosxyz[1000][3];
   Float_t         trkendcosxyz[1000][3];
   Int_t           ntrkhits[1000][3];
   Float_t         trkdqdx[1000][3][3000];
   Float_t         trkdedx[1000][3][3000];
   Float_t         trkresrange[1000][3][3000];
   Float_t         trkhitx[1000][3][3000];
   Float_t         trkhity[1000][3][3000];
   Float_t         trkhitz[1000][3][3000];
   Float_t         trkpitch[1000][3][3000];
   Float_t         trkhitpeakt[1000][3][3000];
   Int_t           trkhittpc[1000][3][3000];
   Int_t           trkhitwire[1000][3][3000];
   Float_t         trkhitrms[1000][3][3000];
   Float_t         trkhitdeltat[1000][3][3000];
   Float_t         dist_min[1000];
   Int_t           adjacent_hits[1000];
   Int_t           lastwire[1000];
   Float_t         lastpeakt[1000];
   Int_t           endtpc[1000];
   Int_t           trkanodetoanode[1000];

   TBranch        *b_event;
   TBranch        *b_evttime;
   TBranch        *b_run;
   TBranch        *b_surbrun;
   TBranch        *b_year_month_date;
   TBranch        *b_hour_min_sec;
   TBranch        *b_cross_trks;
   TBranch        *b_stopping_trks;
   TBranch        *b_all_trks;
   TBranch        *b_unbroken_trks;
   TBranch        *b_trackthetaxz;
   TBranch        *b_trackthetayz;
   TBranch        *b_trkstartx;
   TBranch        *b_trkstarty;
   TBranch        *b_trkstartz;
   TBranch        *b_trkendx;
   TBranch        *b_trkendy;
   TBranch        *b_trkendz;
   TBranch        *b_trklen;
   TBranch        *b_peakT_max;
   TBranch        *b_peakT_min;
   TBranch        *b_TrkID;
   TBranch        *b_trkstartcosxyz;
   TBranch        *b_trkendcosxyz;
   TBranch        *b_ntrkhits;
   TBranch        *b_trkdqdx;
   TBranch        *b_trkdedx;
   TBranch        *b_trkresrange;
   TBranch        *b_trkhitx;
   TBranch        *b_trkhity;
   TBranch        *b_trkhitz;
   TBranch        *b_trkpitch;
   TBranch        *b_trkhitpeakt;
   TBranch        *b_trkhittpc;
   TBranch        *b_trkhitwire;
   TBranch        *b_trkhitrms;
   TBranch        *b_trkhitdeltat;
   TBranch        *b_dist_min;
   TBranch        *b_adjacent_hits;
   TBranch        *b_lastwire;
   TBranch        *b_lastpeakt;
   TBranch        *b_endtpc;
   TBranch        *b_trkanodetoanode;

   data_crossing_pdhd(TTree *tree=0);   // default: chain the hardcoded filelist
   // read many files into a TChain: explicit list, or a text filelist path
   data_crossing_pdhd(const std::vector<std::string>& files,
                      const char* treeName="michelremoving/Event");
   data_crossing_pdhd(const TString& fileList,
                      const char* treeName="michelremoving/Event");
   virtual ~data_crossing_pdhd();
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
};

#endif

#ifdef data_crossing_pdhd_cxx

// Default inputs, hardcoded here: a text file listing one ROOT file path per
// line (blank lines and '#' comments are skipped). Edit this path -- or the
// list it points to -- to change what the no-argument constructor reads.
static const char* kDataCrossingFileList = "/exp/dune/data/users/rdiurba/michelremoving_pdhd.list";

// Append every ROOT file named in a text filelist to the chain.
static void data_crossing_pdhd_AddList(TChain* chain, const char* listPath)
{
   std::ifstream in(listPath);
   std::string line;
   while (std::getline(in, line)) {
      const size_t a = line.find_first_not_of(" \t\r\n");
      if (a == std::string::npos) continue;                 // blank line
      const size_t b = line.find_last_not_of(" \t\r\n");
      const std::string path = line.substr(a, b - a + 1);
      if (path.empty() || path[0] == '#') continue;         // comment
      chain->Add(path.c_str());
   }
}

// Default: chain every file in the hardcoded filelist (kDataCrossingFileList).
data_crossing_pdhd::data_crossing_pdhd(TTree *tree) : fChain(0)
{
   if (tree == 0) {
      TChain* chain = new TChain("michelremovingpdhd/Event");
      data_crossing_pdhd_AddList(chain, kDataCrossingFileList);
      tree = chain;
   }
   Init(tree);
}

// Build a TChain from an explicit list of ROOT files.
data_crossing_pdhd::data_crossing_pdhd(const std::vector<std::string>& files, const char* treeName) : fChain(0)
{
   TChain* chain = new TChain(treeName);
   for (const std::string& fn : files) chain->Add(fn.c_str());
   Init(chain);
}

// Build a TChain from a text filelist (one ROOT file path per line).
data_crossing_pdhd::data_crossing_pdhd(const TString& fileList, const char* treeName) : fChain(0)
{
   TChain* chain = new TChain(treeName);
   data_crossing_pdhd_AddList(chain, fileList.Data());
   Init(chain);
}

data_crossing_pdhd::~data_crossing_pdhd()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t data_crossing_pdhd::GetEntry(Long64_t entry)
{
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}

Long64_t data_crossing_pdhd::LoadTree(Long64_t entry)
{
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void data_crossing_pdhd::Init(TTree *tree)
{
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;

   fChain->SetBranchAddress("event", &event, &b_event);
   fChain->SetBranchAddress("evttime", &evttime, &b_evttime);
   fChain->SetBranchAddress("run", &run, &b_run);
   fChain->SetBranchAddress("subrun", &subrun, &b_surbrun);
   fChain->SetBranchAddress("year_month_date", &year_month_date, &b_year_month_date);
   fChain->SetBranchAddress("hour_min_sec", &hour_min_sec, &b_hour_min_sec);
   fChain->SetBranchAddress("cross_trks", &cross_trks, &b_cross_trks);
   fChain->SetBranchAddress("stopping_trks", &stopping_trks, &b_stopping_trks);
   fChain->SetBranchAddress("all_trks", &all_trks, &b_all_trks);
   fChain->SetBranchAddress("unbroken_trks", &unbroken_trks, &b_unbroken_trks);
   fChain->SetBranchAddress("trackthetaxz", trackthetaxz, &b_trackthetaxz);
   fChain->SetBranchAddress("trackthetayz", trackthetayz, &b_trackthetayz);
   fChain->SetBranchAddress("trkstartx", trkstartx, &b_trkstartx);
   fChain->SetBranchAddress("trkstarty", trkstarty, &b_trkstarty);
   fChain->SetBranchAddress("trkstartz", trkstartz, &b_trkstartz);
   fChain->SetBranchAddress("trkendx", trkendx, &b_trkendx);
   fChain->SetBranchAddress("trkendy", trkendy, &b_trkendy);
   fChain->SetBranchAddress("trkendz", trkendz, &b_trkendz);
   fChain->SetBranchAddress("trklen", trklen, &b_trklen);
   fChain->SetBranchAddress("peakT_max", peakT_max, &b_peakT_max);
   fChain->SetBranchAddress("peakT_min", peakT_min, &b_peakT_min);
   fChain->SetBranchAddress("TrkID", TrkID, &b_TrkID);
   fChain->SetBranchAddress("trkstartcosxyz", trkstartcosxyz, &b_trkstartcosxyz);
   fChain->SetBranchAddress("trkendcosxyz", trkendcosxyz, &b_trkendcosxyz);
   fChain->SetBranchAddress("ntrkhits", ntrkhits, &b_ntrkhits);
   fChain->SetBranchAddress("trkdqdx", trkdqdx, &b_trkdqdx);
   fChain->SetBranchAddress("trkdedx", trkdedx, &b_trkdedx);
   fChain->SetBranchAddress("trkresrange", trkresrange, &b_trkresrange);
   fChain->SetBranchAddress("trkhitx", trkhitx, &b_trkhitx);
   fChain->SetBranchAddress("trkhity", trkhity, &b_trkhity);
   fChain->SetBranchAddress("trkhitz", trkhitz, &b_trkhitz);
   fChain->SetBranchAddress("trkpitch", trkpitch, &b_trkpitch);
   fChain->SetBranchAddress("trkhitpeakt", trkhitpeakt, &b_trkhitpeakt);
   fChain->SetBranchAddress("trkhittpc", trkhittpc, &b_trkhittpc);
   fChain->SetBranchAddress("trkhitwire", trkhitwire, &b_trkhitwire);
   fChain->SetBranchAddress("trkhitrms", trkhitrms, &b_trkhitrms);
   fChain->SetBranchAddress("trkhitdeltat", trkhitdeltat, &b_trkhitdeltat);
   fChain->SetBranchAddress("dist_min", dist_min, &b_dist_min);
   fChain->SetBranchAddress("adjacent_hits", adjacent_hits, &b_adjacent_hits);
   fChain->SetBranchAddress("lastwire", lastwire, &b_lastwire);
   fChain->SetBranchAddress("lastpeakt", lastpeakt, &b_lastpeakt);
   fChain->SetBranchAddress("endtpc", endtpc, &b_endtpc);
   fChain->SetBranchAddress("trkanodetoanode", trkanodetoanode, &b_trkanodetoanode);
   Notify();
}

Bool_t data_crossing_pdhd::Notify()
{
   return kTRUE;
}
#endif // #ifdef data_crossing_pdhd_cxx

