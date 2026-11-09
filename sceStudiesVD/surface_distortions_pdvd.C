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

constexpr double kXmin = -341.56, kXmax = 341.56;
constexpr double kYmin = -340.0, kYmax = 340.0;
constexpr double kZmin =    0.0, kZmax = 300.0;
constexpr double kBin     = 20.0;
constexpr double kFaceBuf = 30.0;
constexpr double kXanodeLo = -341.56, kXanodeHi = 341.56;
constexpr int    kMinEntries = 5;   // a face bin needs >= this many tracks to be measured;
                                    // sparser bins are interpolated (0 at the X anodes)
constexpr double kAnodeTagX = 321.56; // _anode maps: an endpoint with |X| above this is on the anode
                                      // (20 cm from the anode, as 340 is for PDHD's 360)

// Track selection applied to every track (dE/dx-calibration cuts, without the
// fiducial-volume and cathode-crossing requirements):
// angular: remove 65 < |theta_xz| < 115 deg and 70 < |theta_yz| < 110 deg
// (trackthetaxz/yz are atan2(dx,dz), atan2(dy,dz) in radians, as in michelremoving)
constexpr double kThXZlo = 65.0, kThXZhi = 115.0;
constexpr double kThYZlo = 70.0, kThYZhi = 110.0;
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
// output: surface_distortions_pdvd.root, containing
//   distZ_zlo, distZ_zhi, distY_ylo, distY_yhi   nominal: all tracks, any direction
// and the same four maps for sub-samples of the nominal tracks (map stability
// test), with the name suffix
//   _anode                 anode piercing: an endpoint at |X| > kAnodeTagX (either side)
//   _cathode               cathode only: crosses the cathode, no endpoint at |X| > kAnodeTagX
//   _dXpos / _dXneg        heading +X / -X
//   _dZpos / _dZneg        heading +Z / -Z
//   _dYsteep / _dYshallow  |dy| above / below cos(45 deg) of the track length
// The direction is start->end, flipped so the track goes DOWN (dy <= 0), as
// cosmics do, so that the sign of the X and Z components means something.
void surface_distortions_pdvd(const char* fileList = "/exp/dune/app/users/rdiurba/sceStudiesVD/all.txt",
                              const char* treeName = "michelremoving/Event") {

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
  chain->SetBranchStatus("*", 0);
  const char* used[] = {"cross_trks", "trkanodetoanode",
                        "trkstartx", "trkstarty", "trkstartz",
                        "trkendx", "trkendy", "trkendz"};
  for (const char* bn : used) chain->SetBranchStatus(bn, 1);
  chain->SetBranchAddress("cross_trks",      &cross_trks);
  chain->SetBranchAddress("trkanodetoanode", trkanodetoanode);
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

  const int nX = int((kXmax - kXmin) / kBin);
  const int nY = int((kYmax - kYmin) / kBin);
  const int nZ = int((kZmax - kZmin) / kBin);
  const double wX = (kXmax - kXmin) / nX;
  const double wY = (kYmax - kYmin) / nY;
  const double wZ = (kZmax - kZmin) / nZ;
  auto xbin = [&](double x){ return int((x - kXmin) / wX); };
  auto ybin = [&](double y){ return int((y - kYmin) / wY); };
  auto zbin = [&](double z){ return int((z - kZmin) / wZ); };

  // variant 0 = nominal (all tracks, any direction), 1..6 = direction splits,
  // 7 = anode piercing, 8 = cathode only
  const int kNvar = 9;
  const char* varSuffix[kNvar] = {"", "_dXpos", "_dXneg", "_dYsteep", "_dYshallow", "_dZpos", "_dZneg",
                                  "_anode", "_cathode"};
  const char* varTitle[kNvar]  = {"", ", heading +X", ", heading -X", ", steep (|cos_{Y}|>0.71)",
                                  ", shallow (|cos_{Y}|<0.71)", ", heading +Z", ", heading -Z",
                                  ", anode piercing", ", cathode only"};
  const double kSteepCos = std::cos(45.0 * TMath::DegToRad());
  std::vector<std::vector<std::vector<double> > > accZlo(kNvar, std::vector<std::vector<double> >(nX*nY));
  std::vector<std::vector<std::vector<double> > > accZhi(kNvar, std::vector<std::vector<double> >(nX*nY));
  std::vector<std::vector<std::vector<double> > > accYlo(kNvar, std::vector<std::vector<double> >(nX*nZ));
  std::vector<std::vector<std::vector<double> > > accYhi(kNvar, std::vector<std::vector<double> >(nX*nZ));
  // fills the nominal map and the direction-split maps this track belongs to
  std::vector<int> vars;
  auto fillXY = [&](std::vector<std::vector<std::vector<double> > >& acc, int bx, int by, double v){
    if (bx < 0 || bx >= nX || by < 0 || by >= nY) return;
    for (int k : vars) acc[k][bx*nY + by].push_back(v);
  };
  auto fillXZ = [&](std::vector<std::vector<std::vector<double> > >& acc, int bx, int bz, double v){
    if (bx < 0 || bx >= nX || bz < 0 || bz >= nZ) return;
    for (int k : vars) acc[k][bx*nZ + bz].push_back(v);
  };

  const Long64_t nentries = chain->GetEntries();
  for (Long64_t j = 0; j < nentries; ++j) {
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
      if (trkanodetoanode[i]) continue;   // anode-to-anode: never reaches a Y/Z face
      ++nIn;
      const double thxz = std::fabs(trackthetaxz[i]) * TMath::RadToDeg();
      const double thyz = std::fabs(trackthetayz[i]) * TMath::RadToDeg();
      if ((thxz > kThXZlo && thxz < kThXZhi) || (thyz > kThYZlo && thyz < kThYZhi)) { ++nAngle; continue; }
      if (broken[i]) { ++nBroken; continue; }
      if (peakT_min[i] < kPeakTmin || peakT_max[i] > kPeakTmax) { ++nPeakT; continue; }
      if (adjacent_hits[i] > 0 && adjacent_hits[i] != kAdjNotEvaluated) { ++nMichel; continue; }
      ++nPass;
      const bool crossesCathode = trkstartx[i]*trkendx[i] < 0;
      const bool piercesAnode   = std::fabs(trkstartx[i]) > kAnodeTagX || std::fabs(trkendx[i]) > kAnodeTagX;
      {
        const double d[3] = {trkendx[i]-trkstartx[i], trkendy[i]-trkstarty[i], trkendz[i]-trkstartz[i]};
        const double n = std::sqrt(d[0]*d[0] + d[1]*d[1] + d[2]*d[2]);
        if (n <= 0) continue;
        const double flip = (d[1] > 0) ? -1.0 : 1.0;          // downward-going
        const double cx = flip*d[0]/n, cy = flip*d[1]/n, cz = flip*d[2]/n;
        vars.assign(1, 0);                                      // nominal: always
        vars.push_back(cx >= 0 ? 1 : 2);
        vars.push_back(std::fabs(cy) > kSteepCos ? 3 : 4);
        vars.push_back(cz >= 0 ? 5 : 6);
        if (piercesAnode) vars.push_back(7);
        else if (crossesCathode) vars.push_back(8);
      }
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
          if (face == 5 && std::fabs(pz - kZmax) < kFaceBuf) fillXY(accZhi, xbin(px), ybin(py), kZmax - pz);
          if (face == 4 && std::fabs(pz - kZmin) < kFaceBuf) fillXY(accZlo, xbin(px), ybin(py), kZmin - pz);
        }

        if (px >= kXmin && px <= kXmax && pz >= kZmin && pz <= kZmax) {
          if (face == 3 && std::fabs(py - kYmax) < kFaceBuf) fillXZ(accYhi, xbin(px), zbin(pz), kYmax - py);
          if (face == 2 && std::fabs(py - kYmin) < kFaceBuf) fillXZ(accYlo, xbin(px), zbin(pz), kYmin - py);
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

    // bins with >= kMinEntries get the median; sparser ones are marked empty
    for (int ia = 0; ia < nA; ++ia)
      for (int ib = 0; ib < nB; ++ib) {
        const std::vector<double>& v = acc[ia*nB + ib];
        h->SetBinContent(ia+1, ib+1,
          (int)v.size() < kMinEntries ? -999999 : TMath::Median(v.size(), &v[0]));
      }

    // fill the empty bins by LINEAR interpolation along X within each B-row,
    // using the measured bins as knots and the boundary condition that the
    // displacement is 0 at the two X anodes (X = aMin and X = aMax).
    for (int ib = 0; ib < nB; ++ib) {
      std::vector<double> kx, kv;                   // knots for this row, ascending in X
      kx.push_back(aMin); kv.push_back(0.0);        // -anode BC
      for (int ia = 0; ia < nA; ++ia) {
        const double val = h->GetBinContent(ia+1, ib+1);
        if (val > -1e6) { kx.push_back(aMin+(ia+0.5)*wA); kv.push_back(val); }
      }
      kx.push_back(aMax); kv.push_back(0.0);        // +anode BC

      for (int ia = 0; ia < nA; ++ia) {
        if (h->GetBinContent(ia+1, ib+1) > -1e6) continue;   // keep measured bins
        const double xq = aMin+(ia+0.5)*wA;
        size_t k = 0;                                // bracket: kx[k] <= xq < kx[k+1]
        while (k+1 < kx.size() && kx[k+1] < xq) ++k;
        const double x0 = kx[k], x1 = kx[k+1];
        const double t = (x1 > x0) ? (xq - x0) / (x1 - x0) : 0.0;
        h->SetBinContent(ia+1, ib+1, kv[k] + t * (kv[k+1] - kv[k]));
      }
    }
    return h;
  };

  std::vector<TH2F*> maps;
  for (int k = 0; k < kNvar; ++k) {
    // names/titles copied to std::string so no two Form() results are used at once
    const std::string sfx = varSuffix[k], ttl = varTitle[k];
    maps.push_back(build(("distZ_zlo" + sfx).c_str(), ("#DeltaZ upstream face" + ttl).c_str(), accZlo[k],
                         nX, kXmin, wX, "X [cm]", nY, kYmin, wY, "Y [cm]"));
    maps.push_back(build(("distZ_zhi" + sfx).c_str(), ("#DeltaZ downstream face" + ttl).c_str(), accZhi[k],
                         nX, kXmin, wX, "X [cm]", nY, kYmin, wY, "Y [cm]"));
    maps.push_back(build(("distY_ylo" + sfx).c_str(), ("distY on y=ymin face" + ttl).c_str(), accYlo[k],
                         nX, kXmin, wX, "X [cm]", nZ, kZmin, wZ, "Z [cm]"));
    maps.push_back(build(("distY_yhi" + sfx).c_str(), ("distY on y=ymax face" + ttl).c_str(), accYhi[k],
                         nX, kXmin, wX, "X [cm]", nZ, kZmin, wZ, "Z [cm]"));
  }

  std::cout << "tracks in: " << nIn << "  removed by angle: " << nAngle << "  broken: " << nBroken
            << "  early/late hits: " << nPeakT << "  Michel: " << nMichel
            << "  passed: " << nPass << " (before the face fill)" << std::endl;
  const std::string outName = "surface_distortions_pdvd.root";
  TFile fout(outName.c_str(), "RECREATE");
  for (TH2F* h : maps) h->Write();
  fout.Close();
  std::cout << "surface_distortions_pdvd finished: " << outName << std::endl;
}
