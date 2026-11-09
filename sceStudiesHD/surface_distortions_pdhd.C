#include <TChain.h>
#include <TFile.h>
#include <TH2.h>
#include <TMath.h>
#include <TString.h>

// Face maps store the RecoBkwd displacement (reco -> true = true - reco):
// e.g. on the downstream face kZmax - z_reco.

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace {

constexpr double kXmin = -360.0, kXmax = 360.0;
constexpr double kYmin =    0.0, kYmax = 600.0;
constexpr double kZmin =    0.0, kZmax = 460.0;
constexpr double kBin     = 20.0;
constexpr double kFaceBuf = 40.0;
constexpr double kXanodeLo = -360.0, kXanodeHi = 360.0;
constexpr int    kMinEntries = 20; // a face bin needs >= this many tracks to be measured;
                                    // sparser bins are interpolated (0 at the X anodes)
constexpr int    kSmoothPasses = 2; // gentle 3x3 box low-pass passes (like smooth_efield),
                                    // applied after the fill; not across the cathode at X=0
constexpr double kTruncLo = 2.0, kTruncHi = 1.75, kTruncTol = 1e-4; // iterative truncated mean (MicroBooNE wire modification)
constexpr double kAnodeTagX = 340.0; // trkSel=1: an endpoint with |X| above this is on the anode

// Track selection applied to every track (dE/dx-calibration cuts, without the
// fiducial-volume and cathode-crossing requirements):
// angular: remove 65 < |theta_xz| < 115 deg and 70 < |theta_yz| < 110 deg
// (trackthetaxz/yz are atan2(dx,dz), atan2(dy,dz) in radians, as in michelremoving)
	constexpr double kThXZlo = 45.0, kThXZhi = 135.0;
constexpr double kThYZlo = 45.0, kThYZhi = 135.0;
// broken tracks: two tracks with endpoints within 30 cm and < 14 deg apart -> both removed
constexpr double kBrokenDist = 30.0, kBrokenAngle = 14.0;
// early/late hits: remove tracks with any hit peak time < 250 or > 5900 ticks
constexpr double kPeakTmin = 250.0, kPeakTmax = 5900.0;
// Michel hits: remove tracks with adjacent_hits > 0 (99999 = not evaluated by
// michelremoving, i.e. not a stopping track, so not cut)
constexpr int kAdjNotEvaluated = 99999;

inline int exitFace(double px, double py, double pz,
                    double dx, double dy, double dz) {
  double bestT = 1e30; int face = -1;
  auto consider = [&](double t, int f){ if (t > 1e-6 && t < bestT) { bestT = t; face = f; } };
  if (dx < 0) consider((kXanodeLo - px) / dx, 0);
  if (dx > 0) consider((kXanodeHi - px) / dx, 1);
  if (dy < 0) consider((kYmin     - py) / dy, 2);
  if (dy > 0) consider((kYmax     - py) / dy, 3);
  if (dz < 0) consider((kZmin     - pz) / dz, 4);
  if (dz > 0) consider((kZmax     - pz) / dz, 5);
  return face;
}

}

// Anode-to-anode tracks are always skipped (they do not pierce a Y/Z face).
// trkSel: 0 = all tracks (anode and cathode crossers)
//         1 = anode piercing: an endpoint at |X| > kAnodeTagX (either side, may also cross the cathode)
//         2 = cathode only: crosses the cathode (X changes sign) and no endpoint at |X| > kAnodeTagX
// output: surface_distortions_pdhd[_anode|_cathode].root, containing
//   distZ_zlo, distZ_zhi, distY_ylo, distY_yhi
void surface_distortions_pdhd(const char* fileList = "/exp/dune/app/users/rdiurba/sceStudiesHD/v10.txt",
                              const char* treeName = "michelremoving/Event",
                              int trkSel = 0) {

  TChain* chain = new TChain(treeName);
  {
    std::ifstream in(fileList);
    std::string line;
    while (std::getline(in, line)) {
      const size_t a = line.find_first_not_of(" \t\r\n");
      if (a == std::string::npos) continue;
      const size_t b = line.find_last_not_of(" \t\r\n");
      const std::string p = line.substr(a, b - a + 1);
      if (p.empty() || p[0] == '#') continue;
      chain->Add(p.c_str());
    }
  }
  if (chain->GetNtrees() == 0) { std::cout << "no files from " << fileList << std::endl; return; }

  Int_t   cross_trks = 0;
  Int_t   trkanodetoanode[1000];
  Float_t trkstartx[1000], trkstarty[1000], trkstartz[1000];
  Float_t trkendx[1000],   trkendy[1000],   trkendz[1000];
  Float_t trackthetaxz[1000], trackthetayz[1000];
  Float_t peakT_min[1000], peakT_max[1000];
  Int_t   adjacent_hits[1000];
  Float_t trklen[1000];
  chain->SetBranchStatus("*", 0);
  const char* used[] = {"cross_trks", "trkanodetoanode", "trklen",
                        "trkstartx", "trkstarty", "trkstartz",
                        "trkendx", "trkendy", "trkendz"};
  for (const char* bn : used) chain->SetBranchStatus(bn, 1);
  chain->SetBranchAddress("cross_trks",      &cross_trks);
  chain->SetBranchAddress("trkanodetoanode", trkanodetoanode);
  chain->SetBranchAddress("trklen",trklen);
  chain->SetBranchAddress("trkstartx",       trkstartx);
  chain->SetBranchAddress("trkstarty",       trkstarty);
  chain->SetBranchAddress("trkstartz",       trkstartz);
  chain->SetBranchAddress("trkendx",         trkendx);
  chain->SetBranchAddress("trkendy",         trkendy);
  chain->SetBranchAddress("trkendz",         trkendz);
  chain->SetBranchStatus("trackthetaxz", 1);  chain->SetBranchAddress("trackthetaxz",  trackthetaxz);
  chain->SetBranchStatus("trackthetayz", 1);  chain->SetBranchAddress("trackthetayz",  trackthetayz);
  chain->SetBranchStatus("peakT_min", 1);     chain->SetBranchAddress("peakT_min",     peakT_min);
  chain->SetBranchStatus("peakT_max", 1);     chain->SetBranchAddress("peakT_max",     peakT_max);
  chain->SetBranchStatus("adjacent_hits", 1); chain->SetBranchAddress("adjacent_hits", adjacent_hits);
  long nIn = 0, nAngle = 0, nBroken = 0, nPeakT = 0, nMichel = 0, nPass = 0;
  long nNotThrough = 0;                    // tracks not entering AND exiting the detector

  const int nX = int((kXmax - kXmin) / kBin);
  const int nY = int((kYmax - kYmin) / kBin);
  const int nZ = int((kZmax - kZmin) / kBin);
  const double wX = (kXmax - kXmin) / nX;
  const double wY = (kYmax - kYmin) / nY;
  const double wZ = (kZmax - kZmin) / nZ;
  auto xbin = [&](double x){ return int((x - kXmin) / wX); };
  auto ybin = [&](double y){ return int((y - kYmin) / wY); };
  auto zbin = [&](double z){ return int((z - kZmin) / wZ); };

  std::vector<std::vector<double> > accZlo(nX*nY), accZhi(nX*nY), accYlo(nX*nZ), accYhi(nX*nZ);
  auto fillXY = [&](std::vector<std::vector<double> >& acc, int bx, int by, double v){
    if (bx < 0 || bx >= nX || by < 0 || by >= nY) return;
    acc[bx*nY + by].push_back(v);
  };
  auto fillXZ = [&](std::vector<std::vector<double> >& acc, int bx, int bz, double v){
    if (bx < 0 || bx >= nX || bz < 0 || bz >= nZ) return;
    acc[bx*nZ + bz].push_back(v);
  };

  const Long64_t nentries = chain->GetEntries();
  for (Long64_t j = 0; j < nentries && j<200000; ++j) {
    chain->GetEntry(j);
    if (j % 10000 == 0) std::cout << j << "/" << nentries << std::endl;

    const int ntrk = std::min(cross_trks, 1000);

    // broken tracks among the tracks in this event
    std::vector<bool> broken(ntrk, false);
    for (int a = 0; a < ntrk; ++a)
      for (int b = a + 1; b < ntrk; ++b) {
        const double pa[2][3] = {{trkstartx[a], trkstarty[a], trkstartz[a]}, {trkendx[a], trkendy[a], trkendz[a]}};
        const double pb[2][3] = {{trkstartx[b], trkstarty[b], trkstartz[b]}, {trkendx[b], trkendy[b], trkendz[b]}};
        double dmin = 1e30;
        for (int u = 0; u < 2; ++u)
          for (int v = 0; v < 2; ++v)
            dmin = std::min(dmin, std::sqrt(std::pow(pa[u][0]-pb[v][0], 2) + std::pow(pa[u][1]-pb[v][1], 2) + std::pow(pa[u][2]-pb[v][2], 2)));
        if (dmin >= kBrokenDist) continue;
        double da[3], db[3], na = 0, nb = 0, dot = 0;
        for (int c = 0; c < 3; ++c) { da[c] = pa[1][c]-pa[0][c]; db[c] = pb[1][c]-pb[0][c]; na += da[c]*da[c]; nb += db[c]*db[c]; dot += da[c]*db[c]; }
        if (na <= 0 || nb <= 0) continue;
        if (std::fabs(dot) / std::sqrt(na*nb) > std::cos(kBrokenAngle * TMath::DegToRad())) broken[a] = broken[b] = true;
      }

    for (int i = 0; i < ntrk; ++i) {
      //if (trkanodetoanode[i]) continue;   // anode-to-anode: never reaches a Y/Z face
      // through-going: the track must enter and exit the detector, i.e. BOTH end
      // points within kFaceBuf of a detector boundary (X anodes, Y or Z faces);
      // stopping or broken-off tracks that end inside the volume are dropped
      {
        auto nearBoundary = [&](double x, double y, double z) {
          const double d = std::min({kXanodeHi - std::fabs(x), y - kYmin, kYmax - y, z - kZmin, kZmax - z});
          return d < kFaceBuf;
        };
        if (!nearBoundary(trkstartx[i], trkstarty[i], trkstartz[i]) ||
            !nearBoundary(trkendx[i],   trkendy[i],   trkendz[i])) { ++nNotThrough; continue; }
      }
      ++nIn;
      const double thxz = std::fabs(trackthetaxz[i]) * TMath::RadToDeg();
      const double thyz = std::fabs(trackthetayz[i]) * TMath::RadToDeg();
      //if ((thxz > kThXZlo && thxz < kThXZhi) || (thyz > kThYZlo && thyz < kThYZhi)) { ++nAngle; continue; }
      if (trklen[i]<200) continue;
      if (broken[i]) { ++nBroken; continue; }
      if (peakT_min[i] < kPeakTmin || peakT_max[i] > kPeakTmax) { ++nPeakT; continue; }
      //if (adjacent_hits[i] > 0 && adjacent_hits[i] != kAdjNotEvaluated) { ++nMichel; continue; }
      ++nPass;
      const bool crossesCathode = trkstartx[i]*trkendx[i] < 0;
      const bool piercesAnode   = std::fabs(trkstartx[i]) > kAnodeTagX || std::fabs(trkendx[i]) > kAnodeTagX;
      if (trkSel == 1 && !piercesAnode) continue;
      if (trkSel == 2 && (!crossesCathode || piercesAnode)) continue;
      const double ex[2] = {trkstartx[i], trkendx[i]};
      const double ey[2] = {trkstarty[i], trkendy[i]};
      const double ez[2] = {trkstartz[i], trkendz[i]};
      for (int e = 0; e < 2; ++e) {
        const double px = ex[e], py = ey[e], pz = ez[e];
        if (px == 0 && py == 0 && pz == 0) continue;

        const double dx = ex[e] - ex[1-e];
        const double dy = ey[e] - ey[1-e];
        const double dz = ez[e] - ez[1-e];
        const int face = exitFace(px, py, pz, dx, dy, dz);
        if (face < 0) continue;

        if (px >= kXmin && px <= kXmax && py >= kYmin && py <= kYmax) {
          if (face == 5 && std::fabs(pz - kZmax) < kFaceBuf) {
            fillXY(accZhi, xbin(px), ybin(py), kZmax - pz);
          }
          if (face == 4 && std::fabs(pz - kZmin) < kFaceBuf) {
            fillXY(accZlo, xbin(px), ybin(py), kZmin - pz);
          }
        }

        if (px >= kXmin && px <= kXmax && pz >= kZmin && pz <= kZmax) {
          if (face == 3 && std::fabs(py - kYmax) < kFaceBuf) {
            fillXZ(accYhi, xbin(px), zbin(pz), kYmax - py);
          }
          if (face == 2 && std::fabs(py - kYmin) < kFaceBuf) {
            fillXZ(accYlo, xbin(px), zbin(pz), kYmin - py);
          }
        }
      }
    }
  }

  auto build = [&](const char* name, const char* title,
                   const std::vector<std::vector<double> >& acc,
                   int nA, double aMin, double wA, const char* aLab,
                   int nB, double bMin, double wB, const char* bLab) -> TH2F* {
    const double aMax = aMin + nA*wA;   // A axis is X; the anodes are its two edges
    TH2F* h = new TH2F(name, Form("%s;%s;%s", title, aLab, bLab),
                       nA, aMin, aMax, nB, bMin, bMin + nB*wB);
    h->SetDirectory(nullptr);

    // per-bin measured value (median) and statistical weight (track count).
    // Bins below kMinEntries count as unmeasured (weight 0 -> purely interpolated).
    std::vector<double> m(nA*nB, 0.0), w(nA*nB, 0.0);
    for (int ia = 0; ia < nA; ++ia)
      for (int ib = 0; ib < nB; ++ib) {
        const std::vector<double>& v = acc[ia*nB + ib];
        if ((int)v.size() >= kMinEntries) {
          std::vector<double> kept(v);
          double mean = TMath::Mean(kept.size(), &kept[0]);
          for (int it = 0; it < 1000; ++it) {
            const double med = TMath::Median(kept.size(), &kept[0]);
            const double sd = TMath::RMS(kept.size(), &kept[0]);
            std::vector<double> next;
            for (double x : kept) if (x >= med - kTruncLo * sd && x <= med + kTruncHi * sd) next.push_back(x);
            if (next.empty()) break;
            kept.swap(next);
            const double newMean = TMath::Mean(kept.size(), &kept[0]);
            const bool done = std::fabs(newMean - mean) < kTruncTol;
            mean = newMean;
            if (done) break;
          }
          m[ia*nB+ib] = mean;
          w[ia*nB+ib] = (double)kept.size();
        }
      }

    // Measured bins keep their median exactly. Unmeasured bins are filled with the
    // mean of their four neighbours (Gauss-Seidel relaxation), with the hard
    // boundary condition that the displacement is 0 at the two X anodes (no drift
    // there -> no distortion).
    std::vector<double> u(nA*nB, 0.0);
    for (int ia = 0; ia < nA; ++ia)
      for (int ib = 0; ib < nB; ++ib)
        if (w[ia*nB+ib] > 0) u[ia*nB+ib] = m[ia*nB+ib];      // seed measured bins
    for (int iter = 0; iter < 5000; ++iter)
      for (int ia = 0; ia < nA; ++ia)
        for (int ib = 0; ib < nB; ++ib) {
          double sum = 0.0; int cnt = 0;
          sum += (ia == 0)    ? 0.0 : u[(ia-1)*nB+ib]; ++cnt;  // -X (or -anode = 0)
          sum += (ia == nA-1) ? 0.0 : u[(ia+1)*nB+ib]; ++cnt;  // +X (or +anode = 0)
          if (ib > 0)    { sum += u[ia*nB+(ib-1)]; ++cnt; }    // -Y (natural)
          if (ib < nB-1) { sum += u[ia*nB+(ib+1)]; ++cnt; }    // +Y (natural)
          if (w[ia*nB+ib] > 0) continue;                       // measured: fixed at its median
          u[ia*nB+ib] = sum / cnt;
        }
    for (int ia = 0; ia < nA; ++ia)
      for (int ib = 0; ib < nB; ++ib)
        h->SetBinContent(ia+1, ib+1, u[ia*nB+ib]);

    // gentle low-pass, like smooth_efield's boxSmooth: replace each bin by the
    // mean of its 3x3 neighbours on the SAME side of the cathode (X=0), applied
    // kSmoothPasses times. Removes bin-scale scatter without smearing across the
    // central cathode feature.
    for (int p = 0; p < kSmoothPasses; ++p) {
      std::vector<double> src(nA*nB);
      for (int ia = 0; ia < nA; ++ia)
        for (int ib = 0; ib < nB; ++ib) src[ia*nB+ib] = h->GetBinContent(ia+1, ib+1);
      for (int ia = 0; ia < nA; ++ia) {
        const bool xpos = (aMin+(ia+0.5)*wA > 0);
        for (int ib = 0; ib < nB; ++ib) {
          double sum = 0.0; int cnt = 0;
          for (int da = -1; da <= 1; ++da)
            for (int db = -1; db <= 1; ++db) {
              const int ja = ia+da, jb = ib+db;
              if (ja < 0 || ja >= nA || jb < 0 || jb >= nB) continue;
              if ((aMin+(ja+0.5)*wA > 0) != xpos) continue;   // don't smooth across the cathode
              sum += src[ja*nB+jb]; ++cnt;
            }
          if (cnt > 0) h->SetBinContent(ia+1, ib+1, sum / cnt);
        }
      }
    }
    return h;
  };

  std::vector<TH2F*> maps;
  maps.push_back(build("distZ_zlo", "#DeltaZ upstream face",   accZlo, nX, kXmin, wX, "X [cm]", nY, kYmin, wY, "Y [cm]"));
  maps.push_back(build("distZ_zhi", "#DeltaZ downstream face", accZhi, nX, kXmin, wX, "X [cm]", nY, kYmin, wY, "Y [cm]"));
  maps.push_back(build("distY_ylo", "#DeltaY bottom face",     accYlo, nX, kXmin, wX, "X [cm]", nZ, kZmin, wZ, "Z [cm]"));
  maps.push_back(build("distY_yhi", "#DeltaY top face",        accYhi, nX, kXmin, wX, "X [cm]", nZ, kZmin, wZ, "Z [cm]"));

  std::cout << "tracks in: " << nIn << "  removed by angle: " << nAngle << "  broken: " << nBroken
            << "  early/late hits: " << nPeakT << "  Michel: " << nMichel
            << "  passed: " << nPass << " (before the anode/cathode choice and the face fill)" << std::endl;
  std::cout << "tracks not entering and exiting the detector (an end point > " << kFaceBuf
            << " cm from every boundary): " << nNotThrough << std::endl;
  std::string outName = "surface_distortions_pdhd";
  if (trkSel == 1) outName += "_anode";
  if (trkSel == 2) outName += "_cathode";
  outName += ".root";
  TFile fout(outName.c_str(), "RECREATE");
  for (TH2F* h : maps) h->Write();
  fout.Close();
  std::cout << "surface_distortions_pdhd finished: " << outName << std::endl;
}
