#define velocity_pdhd_cxx
#include "velocity_pdhd.h"

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph.h>
#include <TF1.h>
#include <TFile.h>
#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <set>
#include <vector>

namespace {

// ---- ProtoDUNE-HD geometry / reco constants (set to your values) -----------
constexpr double kX0      = -360.0;   // -X anode [cm]
constexpr double kX1      =  360.0;   // +X anode [cm]
constexpr double kVnom    =  0.1566;  // drift velocity assumed by reco [cm/us]
constexpr double kTickUs  =  0.5;     // sampling period [us/tick]
constexpr int    kPlane   =  1;       // hit plane used everywhere (1 = induction)
constexpr double kBin     =  20.0;    // X bin size [cm]
constexpr int    kNbinX   =  36;      // bins over [-360,360]
constexpr int    kFirstPos = kNbinX / 2 + 1;  // first bin with center X > 0

// central Y-Z window in which v(X) is measured (avoids edge effects)
constexpr double kYmin =  100.0, kYmax = 500.0;
constexpr double kZmin =  100.0, kZmax = 360.0;

constexpr double kMinZspan = 30.0;    // required |z0 - z1| of the crosser
constexpr double kCathXtol = 40.0;    // each side's cathode-end reco X must be within this of X=0 [cm]
constexpr double kAnodeXmin = 340.0;  // anode-end hit must have |x| > this
constexpr double kCathXmax  = 10.0;   // cathode-end hit must have |x| < this (within 10 cm of X=0)
constexpr double kXcathP =  0.15, kXcathN = -0.15;   // +/- cathode X [cm]
constexpr double kZtpc   =  230.0;   // TPC boundary in Z [cm] (X boundary is the cathode at 0)

// A genuine anode-to-anode crosser must span the FULL anode->cathode drift on
// both sides, with its two anode-end hits arriving together (t0-aligned).
constexpr double kTanodeTol = 10.0;   // |t_min_pos - t_min_neg| tolerance [ticks]
// accepted summed drift span dT: a narrow window centred on the full-transit peak,
// as in the original analysis (~9200 ticks; originally +/-20, broadened here).
// Independent of kVnom.
// NOTE: this dT is built from hits reaching |x|~355 (anode) to |x|~5 (cathode),
// i.e. ~350 cm per side, so it peaks near ~8970, not the full-360-cm ~9200.
constexpr double kDTcenter = 9000.0;  // window centre [ticks]
constexpr double kDThalf   = 200.0;   // half-width [ticks]
constexpr double kDTlo = kDTcenter - kDThalf;
constexpr double kDThi = kDTcenter + kDThalf;

// how many per-bin (t, x_true) pol1 fit graphs to save for inspection
constexpr size_t kMaxSaveFits = 300;

// LArSoft drift velocity: v [mm/us] vs E [kV/cm] at T [K]. Exact port of
// eField_2.py _drift_raw: piecewise linear / ICARUS / ICARUS-Walkowiak blend /
// Walkowiak (no fudge factor). Same model as smooth_efield_pdhd.
double VDriftModel(double efield, double T = 87.5) {
  const double tshift = -87.203 + T;
  const double xFit = 0.0938163 - 0.0052563*tshift - 0.0001470*tshift*tshift;
  const double uFit = 5.18406 + 0.01448*tshift - 0.003497*tshift*tshift - 0.000516*tshift*tshift*tshift;
  // ICARUS parameters
  const double P1=-0.04640, P2=0.01712, P3=1.88125, P4=0.99408, P5=0.01172, P6=4.20214, T0=105.749;
  // Walkowiak parameters
  const double P1W=-0.01481, P2W=-0.0075, P3W=0.141, P4W=12.4, P5W=1.627, P6W=0.317, T0W=90.371;
  auto icarus = [&](double E){ return (P1*(T-T0)+1.0)*(P3*E*std::log(1.0+P4/E)+P5*std::pow(E,P6)) + P2*(T-T0); };
  auto walk   = [&](double E){ return (P1W*(T-T0W)+1.0)*(P3W*E*std::log(1.0+P4W/E)+P5W*std::pow(E,P6W)) + P2W*(T-T0W); };
  if (efield < xFit)       return efield * uFit;
  if (efield < 0.619)      return icarus(efield);
  if (efield < 0.699) {    const double wW = 12.5*(efield-0.619), wI = 12.5*(0.699-efield);
                           return wW*walk(efield) + wI*icarus(efield); }
  return walk(efield);
}

// invert VDriftModel: v [mm/us] -> E [kV/cm] by bisection.
double EFieldFromV(double v_mm_us) {
  double lo = 0.02, hi = 2.0;
  for (int i = 0; i < 60; ++i) {
    const double mid = 0.5 * (lo + hi);
    if (VDriftModel(mid) < v_mm_us) lo = mid; else hi = mid;
  }
  return 0.5 * (lo + hi);
}

// recombination correction relative to the nominal field (0.4889 kV/cm)
float LAr_density=1.39;
float alp=0.93;
float bet=0.212;
//float dedx=2.08;
float dedx=1.9;
bool userecom=true;

float recom_factor(float totEf){
  if (!userecom) return 1;
  float xsi=bet*dedx/(LAr_density*totEf);
  float xsi0=bet*dedx/(LAr_density*0.4889);
  float rec0=log(alp+xsi0)/xsi0;
  return (rec0*xsi)/log(alp+xsi);
}

// ---- SCE truth E field from the HD map, for comparison with the measured
//      velocity-derived field. Same map/convention as the calibration scripts.
constexpr double kE0 = 0.489;  // nominal drift field [kV/cm]
TFile* g_sce = TFile::Open("/pnfs/dune/persistent/stash/ProtoDUNE/HD/SpaceCharge/Sim_NoFluidFlow/v2/SCE_PDHD_Sim_NoFluidFlow_E500_v2.root");
//TFile *g_sce=TFile::Open("SCE_PDHD_DataDriven_surface.root");

TH3F* g_ex_pos = g_sce ? (TH3F*)g_sce->Get("Reco_ElecField_X_Pos") : nullptr;
TH3F* g_ey_pos = g_sce ? (TH3F*)g_sce->Get("Reco_ElecField_Y_Pos") : nullptr;
TH3F* g_ez_pos = g_sce ? (TH3F*)g_sce->Get("Reco_ElecField_Z_Pos") : nullptr;
TH3F* g_ex_neg = g_sce ? (TH3F*)g_sce->Get("Reco_ElecField_X_Neg") : nullptr;
TH3F* g_ey_neg = g_sce ? (TH3F*)g_sce->Get("Reco_ElecField_Y_Neg") : nullptr;
TH3F* g_ez_neg = g_sce ? (TH3F*)g_sce->Get("Reco_ElecField_Z_Neg") : nullptr;

double TruthEfield(double x, double y, double z) {
  TH3F* hx = (x >= 0) ? g_ex_pos : g_ex_neg;
  TH3F* hy = (x >= 0) ? g_ey_pos : g_ey_neg;
  TH3F* hz = (x >= 0) ? g_ez_pos : g_ez_neg;
  if (!hx || !hy || !hz) return 0;
  const double ex = kE0 + kE0 * hx->GetBinContent(hx->FindBin(x, y, z));
  const double ey =       kE0 * hy->GetBinContent(hy->FindBin(x, y, z));
  const double ez =       kE0 * hz->GetBinContent(hz->FindBin(x, y, z));
  return std::sqrt(ex*ex );//+ ey*ey + ez*ez);
}

}  // namespace

void velocity_pdhd::Loop() {
  if (fChain == 0) return;
 std::cout << "v(0.488) = " << VDriftModel(0.488, 87.5) << " mm/us" << std::endl;
          TH1F* X_plot     = new TH1F("X_plot", "", kNbinX, -360, 360);
  TH1F* med_vel    = new TH1F("med_vel", "median drift velocity;X [cm];v [mm/#mus]",
                              kNbinX, -360, 360);
  TH1F* med_Ef     = new TH1F("med_Ef", "median E field;X [cm];E [kV/cm]",
                              kNbinX, -360, 360);
  TH1F* med_vel_Zlo = new TH1F("med_vel_Zlo", "median drift velocity, Z<230;X [cm];v [mm/#mus]", kNbinX, -360, 360);
  TH1F* med_vel_Zhi = new TH1F("med_vel_Zhi", "median drift velocity, Z>230;X [cm];v [mm/#mus]", kNbinX, -360, 360);
  TH1F* med_Ef_Zlo  = new TH1F("med_Ef_Zlo",  "median E field, Z<230;X [cm];E [kV/cm]", kNbinX, -360, 360);
  TH1F* med_Ef_Zhi  = new TH1F("med_Ef_Zhi",  "median E field, Z>230;X [cm];E [kV/cm]", kNbinX, -360, 360);
  TH1F* input_Ef   = new TH1F("input_Ef", "MC-estimated E field;X [cm];E [kV/cm]",
                              kNbinX, -360, 360);
  TH1F* deltaT     = new TH1F("deltaT", "summed drift span (all selected);dT [ticks];tracks",
                              1200, 0, 12000);   // 10-tick bins
  TH1F* deltaT_cut = new TH1F("deltaT_cut", "summed drift span (accepted window);dT [ticks];tracks",
                              1200, 0, 12000);
  TH1F* hPeakTspan = new TH1F("peakT_span",
                              "peakT_max #minus peakT_min (one drift);peakT_max #minus peakT_min [ticks];tracks",
                              1000, 0, 12000);
  // where the dT-defining hits actually sit: |x| of the anode-end hits (both sides)
  // and of the cathode-end hits.
  TH1F* hAnodeX = new TH1F("anode_hit_x", "anode-end hit |x|;|x| [cm];hits", 220, 340, 362);
  TH1F* hCathX  = new TH1F("cath_hit_x",  "cathode-end hit |x|;|x| [cm];hits", 200, 0, 20);
  TH1F* hVelA2A    = new TH1F("vel_anode2anode",
                              "anode-to-anode average velocity;v [mm/#mus];tracks", 200, 0.5, 3.0);
  TH1F* hEfieldA2A = new TH1F("Efield_anode2anode",
                              "anode-to-anode E field;E [kV/cm];tracks", 200, 0.2, 1.0);
  std::vector<double> ef_a2a;    // per-track full-drift-derived E field [kV/cm]

  // velocity from JUST the longest drift: v = (one drift length) / (peakT_max-peakT_min).
  // No anode-to-anode splitting -- just the single deepest-to-shallowest hit time.
  TH1F* hVelPeakT = new TH1F("vel_peakT",
                             "longest-drift velocity (peakT);v [mm/#mus];tracks", 200, 0.5, 3.0);
  TH1F* hEfPeakT  = new TH1F("Efield_peakT",
                             "longest-drift E field (peakT);E [kV/cm];tracks", 200, 0.2, 1.0);
  std::vector<double> ef_peakT;

  // velocity from the ENDPOINT hits' actual positions: drift distance per side is
  // |x_anode| (reco X at the anode ~ true X, since SCE~0 there) down to the
  // cathode at X=0 -- NOT the fixed 360, so an angled track whose end hit sits a
  // few cm inside the anode is handled correctly.
  TH1F* hVelEnd = new TH1F("vel_endpoint",
                           "endpoint-hit velocity;v [mm/#mus];tracks", 200, 0.5, 3.0);
  TH1F* hEfEnd  = new TH1F("Efield_endpoint",
                           "endpoint-hit E field;E [kV/cm];tracks", 200, 0.2, 1.0);
  std::vector<double> ef_end;

  // same measurement, but for each drift side of the cathode separately:
  // one anode->cathode transit (360 cm) per side.
  TH1F* hVelNeg = new TH1F("vel_neg", "-X anode-to-cathode velocity;v [mm/#mus];tracks", 200, 0.5, 3.0);
  TH1F* hVelPos = new TH1F("vel_pos", "+X anode-to-cathode velocity;v [mm/#mus];tracks", 200, 0.5, 3.0);
  TH1F* hEfNeg  = new TH1F("Efield_neg", "-X anode-to-cathode E field;E [kV/cm];tracks", 200, 0.2, 1.0);
  TH1F* hEfPos  = new TH1F("Efield_pos", "+X anode-to-cathode E field;E [kV/cm];tracks", 200, 0.2, 1.0);
  std::vector<double> ef_neg, ef_pos;

  // same measurement, split further by TPC: 4 TPCs = (X<0,X>0) x (Z<kZtpc,Z>kZtpc).
  const char* tpcname[4] = {"Xneg_Zlo", "Xneg_Zhi", "Xpos_Zlo", "Xpos_Zhi"};
  std::vector<TH1F*> hVelTPC(4), hEfTPC(4);
  std::vector<std::vector<double> > ef_tpc(4);
  for (int t = 0; t < 4; ++t) {
    hVelTPC[t] = new TH1F(Form("vel_tpc_%s", tpcname[t]),
                          Form("TPC %s anode-to-cathode velocity;v [mm/#mus];tracks", tpcname[t]), 200, 0.5, 3.0);
    hEfTPC[t]  = new TH1F(Form("Efield_tpc_%s", tpcname[t]),
                          Form("TPC %s anode-to-cathode E field;E [kV/cm];tracks", tpcname[t]), 200, 0.2, 1.0);
  }

  // per-X-bin velocity distributions, one median per bin (one entry per track)
  // diagnostic: per-bin fit velocity vs the track's Z span. The Z->X map scales
  // any Z error by ~720/|z1-z0|, so short-Z-span tracks give noisy slopes; if high
  // v clusters at small |z1-z0|, raise kMinZspan.
  TH2F* hVvsDz = new TH2F("v_vs_dz", "per-bin fit v vs track Z span;|z1-z0| [cm];v [mm/#mus]",
                          46, 0, 460, 125, 0.5, 3.0);
  std::vector<TH1F*> gr_velocity(kNbinX);
  for (int i = 0; i < kNbinX; ++i) {
    gr_velocity[i] = new TH1F(Form("gr_vel_%d", i),
                              Form("drift velocity, X = %.0f cm;v [mm/#mus];entries",
                                   -360 + i * kBin + 0.5 * kBin),
                              120, 1.0, 2.2);
  }

  // Z-distortion correction from the SCE map (RecoBkwd = true - reco), split by
  // drift side (Pos = x>=0, Neg = x<0).
  //TFile* fdist = TFile::Open("/pnfs/dune/persistent/stash/ProtoDUNE/HD/SpaceCharge/Sim_NoFluidFlow/v2/SCE_PDHD_Sim_NoFluidFlow_E500_v2.root", "READ");
  TFile* fdist=TFile::Open("SCE_PDHD_DataDriven_surface.root");
  
  const bool okdist = (fdist && !fdist->IsZombie());
  TH3F* zBkwdPos = okdist ? (TH3F*)fdist->Get("RecoBkwd_Displacement_Z_Pos") : nullptr;
  TH3F* zBkwdNeg = okdist ? (TH3F*)fdist->Get("RecoBkwd_Displacement_Z_Neg") : nullptr;
  std::cout << ((zBkwdPos && zBkwdNeg)
                ? "using RecoBkwd_Displacement_Z"
                : "no RecoBkwd_Displacement_Z -> using raw Z") << std::endl;

  // expected summed anode->cathode->anode drift span for a full crosser
  const double vTest = VDriftModel(0.488, 87.5) / 10.0;   // [cm/us] at 0.488 kV/cm, 87.5 K
  const double dT_full = 2.0 * kX1 / (kTickUs * vTest);
  std::cout << "dT_full (0.488 kV/cm, 87.5 K) = " << dT_full
            << " ticks ; accepted window [" << kDTlo << ", " << kDThi << "]" << std::endl;

  // one velocity per (track, bin); the bin's median over tracks is med_vel
  std::vector<std::vector<double> > velocities(kNbinX);
  std::vector<std::vector<double> > velocities_lo(kNbinX), velocities_hi(kNbinX);  // per TPC (Z region)

  Long64_t nSel = 0, nWindow = 0;   // selected crossers, and those past the dT window
  std::vector<TGraph*> g_fits;      // a sample of the per-track per-bin fits, saved for inspection

  const Long64_t nentries = fChain->GetEntries();
  for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
    if (LoadTree(jentry) < 0) break;
    fChain->GetEntry(jentry);
    if (jentry % 10000 == 0) std::cout << jentry << "/" << nentries << std::endl;

    for (int i = 0; i < cross_trks; ++i) {
      if (!trkanodetoanode[i]) continue;
      //if (std::fabs(trkstartx[i]) < 355 || std::fabs(trkendx[i]) < 355) continue;

      // --- split hits by drift side (sign of reco X), using plane kPlane ------
      std::vector<double> t_pos, z_pos, y_pos, x_pos, t_neg, z_neg, y_neg, x_neg;
      const int nh = std::min(ntrkhits[i][kPlane], 3000);
      for (int j = 0; j < nh; ++j) {
        const double x = trkhitx[i][kPlane][j];
        const double y = trkhity[i][kPlane][j];
        const double z = trkhitz[i][kPlane][j];
        const double t = trkhitpeakt[i][kPlane][j];
        if (x == 0 && y == 0 && z == 0) continue;
        if (t <= 0) continue;                      // no valid hit time
        if (x > 0) { t_pos.push_back(t); z_pos.push_back(z); y_pos.push_back(y); x_pos.push_back(x); }
        else       { t_neg.push_back(t); z_neg.push_back(z); y_neg.push_back(y); x_neg.push_back(x); }
      }
      if (t_pos.size() < 5 || t_neg.size() < 5) continue;
      ++nSel;
      hPeakTspan->Fill(peakT_max[i] - peakT_min[i]);   // one-drift span from the ntuple extremes

      // sort each side BY TIME (carrying z,y): neg ascending (anode->cathode),
      // pos descending, so front()/back() are the true min/max-time (anode/cathode)
      // hits and the concatenation neg+pos runs monotonically along the track.
      // (A plain end-to-end reverse fails whenever the hits are not already stored
      // in time order -- the extremes then sit in the interior and the span breaks.)
      auto sortByTime = [](std::vector<double>& t, std::vector<double>& z,
                           std::vector<double>& y, std::vector<double>& x, bool ascending) {
        std::vector<size_t> idx(t.size());
        for (size_t k = 0; k < idx.size(); ++k) idx[k] = k;
        std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b){
          return ascending ? (t[a] < t[b]) : (t[a] > t[b]); });
        std::vector<double> tt(t.size()), zz(t.size()), yy(t.size()), xx(t.size());
        for (size_t k = 0; k < idx.size(); ++k) { tt[k]=t[idx[k]]; zz[k]=z[idx[k]]; yy[k]=y[idx[k]]; xx[k]=x[idx[k]]; }
        t.swap(tt); z.swap(zz); y.swap(yy); x.swap(xx);
      };
      sortByTime(t_neg, z_neg, y_neg, x_neg, true);
      sortByTime(t_pos, z_pos, y_pos, x_pos, false);

      const double t_neg_anode = t_neg.front();    // min time on -X side (anode)
      const double t_neg_cath  = t_neg.back();      // max time on -X side (cathode)
      const double t_pos_anode = t_pos.back();      // min time on +X side (anode)
      const double t_pos_cath  = t_pos.front();     // max time on +X side (cathode)
      // reco X of those same four hits -- should be: anodes ~ +/-360, cathodes ~ 0
      const double x_neg_anode = x_neg.front();
      const double x_neg_cath  = x_neg.back();
      const double x_pos_anode = x_pos.back();
      const double x_pos_cath  = x_pos.front();
      hAnodeX->Fill(std::fabs(x_neg_anode)); hAnodeX->Fill(std::fabs(x_pos_anode));
      hCathX->Fill(std::fabs(x_neg_cath));   hCathX->Fill(std::fabs(x_pos_cath));

      // Both drift volumes share the same references: the cathode (X=0) at the
      // deepest drift, and the anode (drift=0) at the common t0. A side can be
      // missing plane-1 hits near one end (e.g. trk-44: -X shallowest hit at
      // X~-295, not -360), so take the shared cathode = max of the two max-times
      // and the shared anode = min of the two min-times. The full anode-to-anode
      // transit is then 2*(t_cath - t_anode) -> ~2*(peakT_max-peakT_min).
      const double t_cath   = std::max(t_neg_cath, t_pos_cath);   // shared cathode (X=0)
      const double t_anode  = std::min(t_neg_anode, t_pos_anode); // shared anode (t0)
      const double span_neg = t_cath - t_neg_anode;   // per-side spans (own anode, shared cathode)
      const double span_pos = t_cath - t_pos_anode;
      const double dT       = 2.0 * (t_cath - t_anode);
      deltaT->Fill(dT);

      // require the dT-defining hits to genuinely span anode(|x|~350) -> cathode(|x|<10)
      // on BOTH sides (what trkanodetoanode should guarantee). Drops tracks whose
      // anode-end hit sits well inside the anode (e.g. trk-44 with x_neg_anode ~ -295).
      if (std::fabs(x_neg_anode) < kAnodeXmin || std::fabs(x_pos_anode) < kAnodeXmin) continue;
      if (std::fabs(x_neg_cath)  > kCathXmax  || std::fabs(x_pos_cath)  > kCathXmax)  continue;

      if (std::fabs(t_pos_anode - t_neg_anode) > kTanodeTol) continue;

      // --- CHECK: the X and T of the four hits that define dT, for the first 50
      //     tracks that PASS the cuts. Anodes should be at |X|>340 (smallest T),
      //     cathodes at |X|<10 (largest T).
      if (nSel <= 50) {
        std::cout << "trk " << nSel
                  << " | -X anode (X,T)=(" << x_neg_anode << "," << t_neg_anode << ")"
                  << " cath (X,T)=("       << x_neg_cath  << "," << t_neg_cath  << ")"
                  << " | +X anode (X,T)=(" << x_pos_anode << "," << t_pos_anode << ")"
                  << " cath (X,T)=("       << x_pos_cath  << "," << t_pos_cath  << ")"
                  << " | t_anode=" << t_anode << " t_cath=" << t_cath << " dT=" << dT
                  << " 2*(peakTmax-peakTmin)=" << 2.0*(peakT_max[i]-peakT_min[i]) << std::endl;
      }

      // --- endpoint-hit velocity: use the anode hits' ACTUAL positions as the
      //     drift distance (|x_anode| ~ true X since SCE~0 at the anode) down to
      //     the cathode at X=0, over the shared-cathode time span. Corrects the
      //     angled-track case where the end hit isn't right against the anode.
      const double distEnd = std::fabs(x_neg_anode) + std::fabs(x_pos_anode);   // cm
      const double timeEnd = (t_cath - t_neg_anode) + (t_cath - t_pos_anode);   // ticks
      if (timeEnd > 0) {
        const double vend = distEnd / (timeEnd * kTickUs) * 10.0;   // mm/us
        if (vend > 0.5 && vend < 3.0) {
          hVelEnd->Fill(vend);
          const double E = EFieldFromV(vend); ef_end.push_back(E); hEfEnd->Fill(E);
        }
      }

      // --- longest-drift velocity from peakT ALONE: one drift length (kX1) over
      //     the single deepest-to-shallowest hit-time span. No anode-to-anode.
      const double dtPeak = peakT_max[i] - peakT_min[i];
      if (dtPeak > 0) {
        const double vpk = kX1 / (dtPeak * kTickUs) * 10.0;   // mm/us
        if (vpk > 0.5 && vpk < 3.0) {
          hVelPeakT->Fill(vpk);
          const double E = EFieldFromV(vpk); ef_peakT.push_back(E); hEfPeakT->Fill(E);
        }
      }

      // --- anode-to-anode average velocity from the full transit time ---------
      const double xCathNeg = kX0 + (t_neg_cath - t_neg_anode) * kTickUs * kVnom;
      const double xCathPos = kX1 - (t_pos_cath - t_pos_anode) * kTickUs * kVnom;
      if (dT > 0 &&
          std::fabs(xCathNeg) < kCathXtol &&
          std::fabs(xCathPos) < kCathXtol) {
        const double vavg = (2.0 * kX1) / (dT * kTickUs) * 10.0;   // mm/us
        if (vavg > 0.5 && vavg < 3.0) {
          hVelA2A->Fill(vavg);
          const double Ea = EFieldFromV(vavg);
          ef_a2a.push_back(Ea); hEfieldA2A->Fill(Ea);
        }
      }

      // --- per cathode side, and per TPC (Z region) ---------------------------
      auto meanOf = [](const std::vector<double>& v)->double{
        if (v.empty()) return 0.0; double s = 0; for (double x : v) s += x; return s / v.size(); };
      // span_neg / span_pos already use the shared cathode time (computed above)

      if (span_neg > 0) {
        const double vneg = (-kX0) / (span_neg * kTickUs) * 10.0;
        if (vneg > 0.5 && vneg < 3.0) {
          hVelNeg->Fill(vneg);
          const double E = EFieldFromV(vneg); ef_neg.push_back(E); hEfNeg->Fill(E);
          const int it = (meanOf(z_neg) > kZtpc) ? 1 : 0;
          hVelTPC[it]->Fill(vneg); ef_tpc[it].push_back(E); hEfTPC[it]->Fill(E);
        }
      }
      if (span_pos > 0) {
        const double vpos = (kX1) / (span_pos * kTickUs) * 10.0;
        if (vpos > 0.5 && vpos < 3.0) {
          hVelPos->Fill(vpos);
          const double E = EFieldFromV(vpos); ef_pos.push_back(E); hEfPos->Fill(E);
          const int it = 2 + ((meanOf(z_pos) > kZtpc) ? 1 : 0);
          hVelTPC[it]->Fill(vpos); ef_tpc[it].push_back(E); hEfTPC[it]->Fill(E);
        }
      }

      if (dT < kDTlo || dT > kDThi) continue;   // keep the full-transit peak, cut the tail
      ++nWindow;
      deltaT_cut->Fill(dT);                     // accepted-window distribution

      // concatenate (neg then pos) into the full anode-to-anode path.
      std::vector<double> t_all, z_all, y_all;
      std::vector<int>    side_all;                 // -1 = neg, +1 = pos
      auto append = [&](const std::vector<double>& t, const std::vector<double>& z,
                        const std::vector<double>& y, int side) {
        for (size_t k = 0; k < t.size(); ++k) {
          t_all.push_back(t[k]); z_all.push_back(z[k]);
          y_all.push_back(y[k]); side_all.push_back(side);
        }
      };
      append(t_neg, z_neg, y_neg, -1);
      append(t_pos, z_pos, y_pos, +1);
      if (t_all.size() < 10) continue;

      // drop interior hits with a repeated Z (keep the two endpoints).
      std::vector<double> tc, zc, yc;
      std::vector<int>    sc;
      const size_t N = z_all.size();
      for (size_t k = 0; k < N; ++k) {
        bool keep = (k == 0 || k == N - 1);
        if (!keep) {
          keep = true;
          for (size_t m = 0; m < N; ++m)
            if (m != k && z_all[m] == z_all[k]) { keep = false; break; }
        }
        if (keep) { tc.push_back(t_all[k]); zc.push_back(z_all[k]);
                    yc.push_back(y_all[k]); sc.push_back(side_all[k]); }
      }

      const double z0 = zc.front();
      const double z1 = zc.back();
      if (std::fabs(z0 - z1) < kMinZspan) continue;

      // --- accumulate this track's (true X, time) per X bin -------------------
      std::vector<std::vector<double> > xval(kNbinX), tval(kNbinX);              // all hits -> med_vel
      std::vector<std::vector<double> > xlo(kNbinX), tlo(kNbinX), xhi(kNbinX), thi(kNbinX);  // per TPC
      for (size_t it = 0; it < tc.size(); ++it) {
        if (it < 5 || it + 5 > tc.size()) continue;                // skip ends
        const double y = yc[it], z = zc[it], t = tc[it];
        if (!(z > kZmin && z < kZmax && y > kYmin && y < kYmax)) continue;

        const double xcal = (sc[it] > 0)
            ? kX1 - (t - t_pos_anode) * kTickUs * kVnom
            : kX0 + (t - t_neg_anode) * kTickUs * kVnom;
        if (xcal < -360 || xcal > 360) continue;

        double zcorr = z;
        TH3F* zmap = (xcal >= 0) ? zBkwdPos : zBkwdNeg;
        if (zmap) {
          double xcal2 = xcal;
          if (std::fabs(xcal) > 350) xcal2 = 349;
          const double dz = zmap->Interpolate(xcal2, y, z);
          if (std::fabs(dz) < 40.0) zcorr = z + dz;
        }

        // Z -> true X along the anode-to-anode line, anchored at the ACTUAL end hits
        // (reco X ~ true X at the anodes, where SCE ~ 0) instead of forcing them to
        // +/-360. Angled tracks' end hits sit a few cm inside the anode; stretching
        // them to +/-360 inflates xtrue, and hence the per-bin slope (v), by
        // 720/(|x_neg_anode|+|x_pos_anode|).
        const double xtrue = x_neg_anode + (x_pos_anode - x_neg_anode) * (zcorr - z0) / (z1 - z0);
        const int xbin = X_plot->FindBin(xtrue);
        if (sc[it] > 0 && xbin < kFirstPos) continue;   // +X hits -> +X bins
        if (sc[it] < 0 && xbin >= kFirstPos) continue;  // -X hits -> -X bins
        if (xbin < 1 || xbin > kNbinX) continue;

        xval[xbin - 1].push_back(xtrue);        // combined (all Z) -> med_vel
        tval[xbin - 1].push_back(t);
        std::vector<std::vector<double> >& xv = (z > kZtpc) ? xhi : xlo;   // per TPC
        std::vector<std::vector<double> >& tv = (z > kZtpc) ? thi : tlo;
        xv[xbin - 1].push_back(xtrue);
        tv[xbin - 1].push_back(t);
      }

      // --- per X bin, fit this track's (t, x_true) -> slope = drift velocity.
      //     Straight-line fit; no ordering requirement. isCombined => med_vel.
      auto fitBins = [&](std::vector<std::vector<double> >& xval_,
                         std::vector<std::vector<double> >& tval_,
                         std::vector<std::vector<double> >& vel,
                         bool isCombined){
        for (int k = 0; k < kNbinX; ++k) {
          if (xval_[k].size() <= 5) continue;
          TGraph gr((int)tval_[k].size(), &tval_[k][0], &xval_[k][0]);
          gr.Fit("pol1", "Q");
          TF1* f1 = gr.GetFunction("pol1");
          if (!f1) continue;
          const double v = std::fabs(20.0 * f1->GetParameter(1));   // mm/us
          if (v <= 0.5 || v >= 3.0) continue;                       // sanity
          vel[k].push_back(v);
          if (isCombined) {
            gr_velocity[k]->Fill(v);
            hVvsDz->Fill(std::fabs(z1 - z0), v);
            if (g_fits.size() < kMaxSaveFits) {
              TGraph* g = new TGraph((int)tval_[k].size(), &tval_[k][0], &xval_[k][0]);
              g->SetName(Form("fit_b%02d_n%03zu", k, g_fits.size()));
              g->SetTitle(Form("X=%.0f cm  v=%.3f mm/#mus;t [ticks];x_{true} [cm]",
                               -360 + k*kBin + 0.5*kBin, v));
              g->SetMarkerStyle(20); g->SetMarkerSize(0.6);
              g->Fit("pol1", "Q");
              g_fits.push_back(g);
            }
          }
        }
      };
      fitBins(xval, tval, velocities,    true);    // med_vel: one fit per track per bin
      fitBins(xlo,  tlo,  velocities_lo, false);   // per-TPC: Z<230
      fitBins(xhi,  thi,  velocities_hi, false);   // per-TPC: Z>230
    }  // track loop
  }    // event loop

  std::cout << "selected crossers (>=5 hits/side) = " << nSel
            << " ; passing dT window = " << nWindow << std::endl;
  std::cout << "tracks contributing per X bin:";
  for (int k = 0; k < kNbinX; ++k) std::cout << " " << velocities[k].size();
  std::cout << std::endl;

  TH1F* recom_x = new TH1F("recom_factor", "recombination factor from med_Ef;X [cm];recombination factor",
                           kNbinX, -360, 360);
  for (int k = 0; k < kNbinX; ++k) {
    if (velocities[k].empty()) continue;
    const double vmed = TMath::Median(velocities[k].size(), &velocities[k][0]);
    const double xc = -360 + k * kBin + 0.5 * kBin;   // bin-center X
    med_vel->SetBinContent(k + 1, vmed);              // median velocity vs X
    med_Ef->SetBinContent(k + 1, EFieldFromV(vmed));  // median E field vs X
    recom_x->SetBinContent(k + 1, recom_factor(EFieldFromV(vmed)));  // recombination factor vs X
    input_Ef->SetBinContent(k + 1, TruthEfield(xc, 0.5*(kYmin+kYmax), 0.5*(kZmin+kZmax)));
  }

  // per-TPC med_vel / med_Ef (slope-fit method): Z<kZtpc (Zlo) and Z>kZtpc (Zhi)
  for (int k = 0; k < kNbinX; ++k) {
    if (!velocities_lo[k].empty()) {
      const double v = TMath::Median(velocities_lo[k].size(), &velocities_lo[k][0]);
      med_vel_Zlo->SetBinContent(k + 1, v);
      med_Ef_Zlo->SetBinContent(k + 1, EFieldFromV(v));
    }
    if (!velocities_hi[k].empty()) {
      const double v = TMath::Median(velocities_hi[k].size(), &velocities_hi[k][0]);
      med_vel_Zhi->SetBinContent(k + 1, v);
      med_Ef_Zhi->SetBinContent(k + 1, EFieldFromV(v));
    }
  }

  // --- average measured E field, two independent ways ------------------------
  double sumE = 0; int nv = 0;
  for (int k = 0; k < kNbinX; ++k)
    if (!velocities[k].empty()) { sumE += med_Ef->GetBinContent(k + 1); ++nv; }
  const double E_profile = nv ? sumE / nv : 0.0;
  const double E_a2a = ef_a2a.empty() ? 0.0 : TMath::Median(ef_a2a.size(), &ef_a2a[0]);
  std::cout << "average E (profile, " << nv << " bins) = " << E_profile << " kV/cm" << std::endl;
  std::cout << "average E (anode-to-anode) = " << E_a2a << " kV/cm  ("
            << ef_a2a.size() << " tracks)" << std::endl;
  auto medp = [](std::vector<double>& v)->double{ return v.empty() ? 0.0 : TMath::Median(v.size(), &v[0]); };
  std::cout << "average E (peakT longest)  = " << medp(ef_peakT) << " kV/cm  ("
            << ef_peakT.size() << " tracks)" << std::endl;
  std::cout << "average E (endpoint hits)  = " << medp(ef_end) << " kV/cm  ("
            << ef_end.size() << " tracks)" << std::endl;
  auto med = [](std::vector<double>& v)->double{ return v.empty() ? 0.0 : TMath::Median(v.size(), &v[0]); };
  std::cout << "E (-X side)  = " << med(ef_neg) << " kV/cm  (" << ef_neg.size() << " tracks)" << std::endl;
  std::cout << "E (+X side)  = " << med(ef_pos) << " kV/cm  (" << ef_pos.size() << " tracks)" << std::endl;
  for (int t = 0; t < 4; ++t)
    std::cout << "E (TPC " << tpcname[t] << ") = " << med(ef_tpc[t])
              << " kV/cm  (" << ef_tpc[t].size() << " tracks)" << std::endl;
  TH1D* Efield_summary = new TH1D("Efield_summary", "measured average E field;;E [kV/cm]", 2, 0, 2);
  Efield_summary->GetXaxis()->SetBinLabel(1, "E_profile");
  Efield_summary->GetXaxis()->SetBinLabel(2, "E_anode2anode");
  Efield_summary->SetBinContent(1, E_profile);
  Efield_summary->SetBinContent(2, E_a2a);

  TFile fout("velocity_pdhd.root", "RECREATE");
  X_plot->Write();
  med_vel->Write();
  med_Ef->Write();
  recom_x->Write();
  med_vel_Zlo->Write(); med_vel_Zhi->Write();
  med_Ef_Zlo->Write();  med_Ef_Zhi->Write();
  input_Ef->Write();
  deltaT->Write();
  deltaT_cut->Write();
  hPeakTspan->Write();
  hAnodeX->Write();
  hCathX->Write();
  hVelA2A->Write();
  hEfieldA2A->Write();
  hVelPeakT->Write();
  hEfPeakT->Write();
  hVelEnd->Write();
  hEfEnd->Write();
  hVelNeg->Write();  hVelPos->Write();
  hEfNeg->Write();   hEfPos->Write();
  for (int t = 0; t < 4; ++t) { hVelTPC[t]->Write(); hEfTPC[t]->Write(); }
  Efield_summary->Write();
  for (int k = 0; k < kNbinX; ++k) gr_velocity[k]->Write();
  hVvsDz->Write();
  fout.mkdir("pol1_fits")->cd();
  for (TGraph* g : g_fits) g->Write();
  fout.cd();
  fout.Close();
  std::cout << "velocity_pdhd finished" << std::endl;
}
