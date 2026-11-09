#define velocity_pdvd_cxx
#include "velocity_pdvd.h"

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

// ---- ProtoDUNE-VD geometry / reco constants (set to your values) -----------
constexpr double kX0      = -341.56;  // -X anode [cm]
constexpr double kX1      =  341.56;  // +X anode [cm]
constexpr double kVnom    =  0.151409;  // drift velocity assumed by reco [cm/us]
constexpr double kTickUs  =  0.5;     // sampling period [us/tick]
constexpr double kBin     =  20.0;    // X bin size [cm]
constexpr int    kNbinX   =  34;      // bins over [-340,340]
constexpr int    kFirstPos = kNbinX / 2 + 1;  // first bin with center X > 0

// central Y-Z window in which v(X) is measured (avoids edge effects)
constexpr double kYmin = -300.0, kYmax = 300.0;
constexpr double kZmin =  40.0, kZmax = 260.0;

constexpr double kMinZspan = 30.0;    // required |z0 - z1| of the crosser
constexpr double kCathCoinc = 20.0;   // the two sides' near-cathode hits must coincide [cm]
constexpr double kXcathP =  0.15, kXcathN = -0.15;   // +/- cathode X [cm]

// A genuine anode-to-anode crosser must span the FULL anode->cathode drift on
// both sides, with its two anode-end hits arriving together (t0-aligned). The
// original enforced this with a tight window on the summed drift time; without
// it, tracks that stop short of an anode have their z-endpoints stretched onto
// the anode by the straight-line xtrue map, which spuriously inflates the fitted
// velocity right at the anodes. The expected full two-sided span (in ticks) is
//   dT_full = 2 * (drift length) / (tick * v_nominal),   drift length = kX1.
constexpr double kTanodeTol = 10.0;   // |t_min_pos - t_min_neg| tolerance [ticks]
constexpr double kDTloFrac  = 0.98;   // accepted summed drift span, low fraction
constexpr double kDThiFrac  = 1.01;   // accepted summed drift span, high fraction

// LArSoft drift velocity: v [mm/us] vs E [kV/cm] at T [K]. Exact port of
// eField_2.py _drift_raw: piecewise linear / ICARUS / ICARUS-Walkowiak blend /
// Walkowiak (no fudge factor). Same model as smooth_efield_pdvd.
double VDriftModel(double efield, double T = 87.4) {
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

// ---- SCE truth E field from the VD map, for comparison with the measured
//      velocity-derived field. Same map/convention as the calibration scripts.
constexpr double kE0 = 0.450;  // nominal drift field [kV/cm]
TFile* g_sce = TFile::Open("/pnfs/dune/persistent/stash/ProtoDUNE/VD/SpaceCharge/Sim_NoFluidFlow/v2/SCE_PDVD_Sim_NoFluidFlow_E500_v2.root");
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

void velocity_pdvd::Loop() {
  if (fChain == 0) return;
 std::cout << "v(0.450) = " << VDriftModel(0.450, 87.4) << " mm/us" << std::endl;
          TH1F* X_plot     = new TH1F("X_plot", "", kNbinX, -340, 340);
  TH1F* med_vel    = new TH1F("med_vel", "median drift velocity;X [cm];v [mm/#mus]",
                              kNbinX, -340, 340);
  TH1F* med_Ef     = new TH1F("med_Ef", "median E field;X [cm];E [kV/cm]",
                              kNbinX, -340, 340);
  TH1F* input_Ef   = new TH1F("input_Ef", "MC-estimated E field;X [cm];E [kV/cm]",
                              kNbinX, -340, 340);
  TH1F* deltaT     = new TH1F("deltaT", "summed drift span;dT [ticks];tracks",
                              1000, 0, 12000);
  TH1F* hVelA2A    = new TH1F("vel_anode2anode",
                              "anode-to-anode average velocity;v [mm/#mus];tracks", 200, 1.0, 3.0);
  TH1F* hEfieldA2A = new TH1F("Efield_anode2anode",
                              "anode-to-anode E field;E [kV/cm];tracks", 200, 0.2, 1.0);
  std::vector<double> ef_a2a;    // per-track full-drift-derived E field [kV/cm]

  // per-X-bin velocity distributions, one median per bin
  std::vector<TH1F*> gr_velocity(kNbinX);
  for (int i = 0; i < kNbinX; ++i) {
    gr_velocity[i] = new TH1F(Form("gr_vel_%d", i),
                              Form("drift velocity, X = %.0f cm;v [mm/#mus];entries",
                                   -340 + i * kBin + 0.5 * kBin),
                              120, 1.0, 2.2);
  }

  // Z-distortion correction from data_crossing_pdvd, using the RecoBkwd (true -
  // reco) displacement maps, split by drift side (Pos = x>=0, Neg = x<0).
  TFile* fdist = TFile::Open("data_crossing_pdvd.root", "READ");
  const bool okdist = (fdist && !fdist->IsZombie());
  TH3F* zBkwdPos = okdist ? (TH3F*)fdist->Get("RecoBkwd_Displacement_Z_Pos") : nullptr;
  TH3F* zBkwdNeg = okdist ? (TH3F*)fdist->Get("RecoBkwd_Displacement_Z_Neg") : nullptr;
  std::cout << ((zBkwdPos && zBkwdNeg)
                ? "using RecoBkwd_Displacement_Z from data_crossing_pdvd.root"
                : "no RecoBkwd_Displacement_Z -> using raw Z") << std::endl;

  // TH3::Interpolate is only valid between the FIRST and LAST bin CENTERS, not
  // the axis edges: for a +/-340 / 20 cm map the outer centers are at +/-330, so
  // a query in the outer half-bin (330 < |x| <= 340, which our +/-340 cut still
  // allows) is "out of range" and returns 0. Clamp the query into the bin-center
  // hull so the edge just takes the nearest interpolable value.
  auto safeInterp = [](TH3F* h, double x, double y, double z)->double{
    const TAxis* ax = h->GetXaxis(); const TAxis* ay = h->GetYaxis(); const TAxis* az = h->GetZaxis();
    x = std::min(std::max(x, ax->GetBinCenter(1)), ax->GetBinCenter(ax->GetNbins()));
    y = std::min(std::max(y, ay->GetBinCenter(1)), ay->GetBinCenter(ay->GetNbins()));
    z = std::min(std::max(z, az->GetBinCenter(1)), az->GetBinCenter(az->GetNbins()));
    return h->Interpolate(x, y, z);
  };

  // expected summed anode->cathode->anode drift span for a full crosser
  const double dT_full = 2.0*kX1 / (kTickUs * kVnom);

  std::vector<std::vector<double> > velocities(kNbinX);

  const Long64_t nentries = fChain->GetEntries();
  for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
    if (LoadTree(jentry) < 0) break;
    fChain->GetEntry(jentry);
    if (jentry % 10000 == 0) std::cout << jentry << "/" << nentries << std::endl;

    for (int i = 0; i < cross_trks; ++i) {
      if (!trkanodetoanode[i]) continue;

      // --- split collection-plane hits by drift side (sign of reco X) --------
      std::vector<double> t_pos, z_pos, y_pos, t_neg, z_neg, y_neg;
      const int nh = std::min(ntrkhits[i][2], 3000);
      for (int j = 0; j < nh; ++j) {
        const double x = trkhitx[i][2][j];
        const double y = trkhity[i][2][j];
        const double z = trkhitz[i][2][j];
        const double t = trkhitpeakt[i][2][j];
        if (x == 0 && y == 0 && z == 0) continue;
        if (t <= 0) continue;                      // no valid hit time
        if (x > 0) { t_pos.push_back(t); z_pos.push_back(z); y_pos.push_back(y); }
        else       { t_neg.push_back(t); z_neg.push_back(z); y_neg.push_back(y); }
      }
      if (t_pos.size() < 5 || t_neg.size() < 5) continue;

      // order neg ascending in time (anode->cathode), pos descending, so the
      // concatenation neg+pos runs monotonically along the track.
      auto reorder = [](std::vector<double>& t, std::vector<double>& z,
                        std::vector<double>& y, bool ascending) {
        const bool isAsc = t.front() < t.back();
        if (isAsc != ascending) {
          std::reverse(t.begin(), t.end());
          std::reverse(z.begin(), z.end());
          std::reverse(y.begin(), y.end());
        }
      };
      reorder(t_neg, z_neg, y_neg, true);
      reorder(t_pos, z_pos, y_pos, false);

      const double t_neg_anode = t_neg.front();    // min time on -X side (anode)
      const double t_neg_cath  = t_neg.back();      // max time on -X side (cathode)
      const double t_pos_anode = t_pos.back();      // min time on +X side (anode)
      const double t_pos_cath  = t_pos.front();     // max time on +X side (cathode)

      // require a genuine full crosser: both anode ends t0-aligned, and the
      // summed drift span consistent with a full anode->cathode->anode path.
      const double dT = (t_pos_cath - t_pos_anode) + (t_neg_cath - t_neg_anode);
      deltaT->Fill(dT);
      if (std::fabs(t_pos_anode - t_neg_anode) > kTanodeTol) continue;

      // --- actual average drift velocity from the FULL anode-to-anode transit
      //     time, independent of the reco kVnom: total drift distance is 2*kX1
      //     (cathode-to-anode on each side) and total drift time is dT. Require
      //     the two near-cathode hits to coincide so both sides really reach the
      //     cathode-crossing point.
      if (dT > 0 &&
          std::fabs(y_neg.back() - y_pos.front()) < kCathCoinc &&
          std::fabs(z_neg.back() - z_pos.front()) < kCathCoinc) {
        // Total drift distance is 2*kX1 (cathode-to-anode on each side); the
        // average speed is that distance over the total transit time dT.
        const double vavg = (2.0 * kX1) / (dT * kTickUs) * 10.0;   // mm/us
        if (vavg > 0.5 && vavg < 3.0) {
          hVelA2A->Fill(vavg);
          const double Ea = EFieldFromV(vavg);                   // also as E
          ef_a2a.push_back(Ea); hEfieldA2A->Fill(Ea);
        }
      }

      //if (dT < kDTloFrac * dT_full || dT > kDThiFrac * dT_full) continue;

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

      // drop interior hits with a repeated Z (keep the two endpoints), as the
      // straight-line xtrue map is degenerate for duplicated Z.
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

      // --- accumulate (true X, time) per X bin -------------------------------
      std::vector<std::vector<double> > xval(kNbinX), tval(kNbinX);
      for (size_t it = 0; it < tc.size(); ++it) {
        if (it < 5 || it + 5 > tc.size()) continue;                // skip ends
        const double y = yc[it], z = zc[it], t = tc[it];
        if (!(z > kZmin && z < kZmax && y > kYmin && y < kYmax)) continue;

        // reco (drift) X, from the near-anode time; used for the side cut and to
        // look up the Z distortion at the correct drift position.
        const double xcal = (sc[it] > 0)
            ? kX1 - (t - t_pos_anode) * kTickUs * kVnom
            : kX0 + (t - t_neg_anode) * kTickUs * kVnom;
        if (xcal < -340 || xcal > 340) continue;

        // undo the reconstructed Z distortion (RecoBkwd = true - reco) before
        // inferring the true X, using the side-appropriate map.
        double zcorr = z;
        TH3F* zmap = (xcal >= 0) ? zBkwdPos : zBkwdNeg;
        if (zmap) {
	   double xcal2=xcal;
	  if (abs(xcal)>330) xcal2=329;
          const double dz = zmap->Interpolate(xcal2, y, z);
          if (std::fabs(dz) < 40.0) zcorr = z + dz;
        }

        // true X from (corrected) Z along the straight anode-to-anode line
        const double xtrue = kX1 + (kX1 - kX0) * (zcorr - z1) / (z1 - z0);
        const int xbin = X_plot->FindBin(xtrue);
        if (sc[it] > 0 && xbin < kFirstPos) continue;   // +X hits -> +X bins
        if (sc[it] < 0 && xbin >= kFirstPos) continue;  // -X hits -> -X bins
        if (xbin < 1 || xbin > kNbinX) continue;

        xval[xbin - 1].push_back(xtrue);
        tval[xbin - 1].push_back(t);
      }

      // --- per X bin: slope of true-X vs time = drift velocity ---------------
      for (int k = 0; k < kNbinX; ++k) {
        if (xval[k].size() <= 5) continue;

        // clean fit: unique X, and both X and time monotonic (as the original)
        std::set<double> uniq(xval[k].begin(), xval[k].end());
        if (uniq.size() != xval[k].size()) continue;
        std::vector<double> xs = xval[k];
        std::sort(xs.begin(), xs.end());
        if (xs != xval[k] &&
            !std::equal(xs.rbegin(), xs.rend(), xval[k].begin())) continue;
        std::vector<double> ts = tval[k];
        std::sort(ts.begin(), ts.end());
        if (ts != tval[k] &&
            !std::equal(ts.rbegin(), ts.rend(), tval[k].begin())) continue;

        TGraph gr(tval[k].size(), &tval[k][0], &xval[k][0]);
        gr.Fit("pol1", "Q");
        TF1* f1 = gr.GetFunction("pol1");
        if (!f1) continue;
        const double v = std::fabs(20.0 * f1->GetParameter(1));   // mm/us
        if (v <= 0.5 || v >= 3.0) continue;                       // sanity
        velocities[k].push_back(v);
        gr_velocity[k]->Fill(v);
      }
    }
  }

  for (int k = 0; k < kNbinX; ++k) {
    if (velocities[k].empty()) continue;
    const double vmed = TMath::Median(velocities[k].size(), &velocities[k][0]);
    const double xc = -340 + k * kBin + 0.5 * kBin;   // bin-center X
    med_vel->SetBinContent(k + 1, vmed);              // median velocity vs X
    med_Ef->SetBinContent(k + 1, EFieldFromV(vmed));  // median E field vs X
    // MC-estimated E field at the center of the measurement window
    input_Ef->SetBinContent(k + 1, TruthEfield(xc, 0.5*(kYmin+kYmax), 0.5*(kZmin+kZmax)));
  }

  // --- actual average E field, two independent ways --------------------------
  //   (1) average of the per-bin median E field profile E(X);
  //   (2) median of the per-track full anode-to-anode transit E field.
  // Both are already in kV/cm, so the actual field (which is not exactly the
  // nominal 0.450 kV/cm) can be read off directly.
  double sumE = 0; int nv = 0;
  for (int k = 0; k < kNbinX; ++k)
    if (!velocities[k].empty()) { sumE += med_Ef->GetBinContent(k + 1); ++nv; }
  const double E_profile = nv ? sumE / nv : 0.0;
  const double E_a2a = ef_a2a.empty() ? 0.0
                        : TMath::Median(ef_a2a.size(), &ef_a2a[0]);
  std::cout << "average E (profile)        = " << E_profile << " kV/cm" << std::endl;
  std::cout << "average E (anode-to-anode) = " << E_a2a << " kV/cm  ("
            << ef_a2a.size() << " tracks)" << std::endl;
  TH1D* Efield_summary = new TH1D("Efield_summary", "measured average E field;;E [kV/cm]", 2, 0, 2);
  Efield_summary->GetXaxis()->SetBinLabel(1, "E_profile");
  Efield_summary->GetXaxis()->SetBinLabel(2, "E_anode2anode");
  Efield_summary->SetBinContent(1, E_profile);
  Efield_summary->SetBinContent(2, E_a2a);

  TFile fout("velocity_pdvd.root", "RECREATE");
  X_plot->Write();
  med_vel->Write();
  med_Ef->Write();
  input_Ef->Write();
  deltaT->Write();
  hVelA2A->Write();
  hEfieldA2A->Write();
  Efield_summary->Write();
  for (int k = 0; k < kNbinX; ++k) gr_velocity[k]->Write();
  fout.Close();
  std::cout << "velocity_pdvd finished" << std::endl;
}
