#define data_crossing_pdvd_cxx
#include "data_crossing_pdvd.h"

#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TGraph2D.h>
#include <TMath.h>
#include <TFile.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <utility>
#include <vector>

namespace {

constexpr double kX0     = -341.56;
constexpr double kX1     =  341.56;
constexpr double kVnom   =  0.156305;
constexpr double kTickUs =  0.5;
constexpr double kBin    =  20.0;
constexpr double kXcathP =  0.15;
constexpr double kXcathN = -0.15;

constexpr double kXmin = -340.0, kXmax = 340.0;
constexpr double kYmin = -340.0, kYmax = 340.0;
constexpr double kZmin =    0.0, kZmax = 300.0;

constexpr double kMaxDist   = 40.0;
constexpr double kPairDrMax = 30.0;

constexpr double kEfieldNominal = 0.45;

TFile* g_sce = TFile::Open("/pnfs/dune/persistent/stash/ProtoDUNE/VD/SpaceCharge/Sim_NoFluidFlow/v2/SCE_PDVD_Sim_NoFluidFlow_E500_v2.root");
TH3F* g_dx_pos = g_sce ? (TH3F*)g_sce->Get("RecoFwd_Displacement_X_Pos") : nullptr;
TH3F* g_dy_pos = g_sce ? (TH3F*)g_sce->Get("RecoFwd_Displacement_Y_Pos") : nullptr;
TH3F* g_dz_pos = g_sce ? (TH3F*)g_sce->Get("RecoFwd_Displacement_Z_Pos") : nullptr;
TH3F* g_dx_neg = g_sce ? (TH3F*)g_sce->Get("RecoFwd_Displacement_X_Neg") : nullptr;
TH3F* g_dy_neg = g_sce ? (TH3F*)g_sce->Get("RecoFwd_Displacement_Y_Neg") : nullptr;
TH3F* g_dz_neg = g_sce ? (TH3F*)g_sce->Get("RecoFwd_Displacement_Z_Neg") : nullptr;
TH3F* g_bx_pos = g_sce ? (TH3F*)g_sce->Get("RecoBkwd_Displacement_X_Pos") : nullptr;
TH3F* g_by_pos = g_sce ? (TH3F*)g_sce->Get("RecoBkwd_Displacement_Y_Pos") : nullptr;
TH3F* g_bz_pos = g_sce ? (TH3F*)g_sce->Get("RecoBkwd_Displacement_Z_Pos") : nullptr;
TH3F* g_bx_neg = g_sce ? (TH3F*)g_sce->Get("RecoBkwd_Displacement_X_Neg") : nullptr;
TH3F* g_by_neg = g_sce ? (TH3F*)g_sce->Get("RecoBkwd_Displacement_Y_Neg") : nullptr;
TH3F* g_bz_neg = g_sce ? (TH3F*)g_sce->Get("RecoBkwd_Displacement_Z_Neg") : nullptr;

double TruthDist(TH3F* hpos, TH3F* hneg, double x, double y, double z) {
  TH3F* h = (x >= 0) ? hpos : hneg;
  if (!h) return -999999;
  return h->GetBinContent(h->FindBin(x, y, z));
}

TFile* g_face = TFile::Open("surface_distortions_pdvd.root", "READ");
TH2F* g_distZ_lo = g_face ? (TH2F*)g_face->Get("distZ_zlo") : nullptr;
TH2F* g_distZ_hi = g_face ? (TH2F*)g_face->Get("distZ_zhi") : nullptr;
TH2F* g_distY_lo = g_face ? (TH2F*)g_face->Get("distY_ylo") : nullptr;
TH2F* g_distY_hi = g_face ? (TH2F*)g_face->Get("distY_yhi") : nullptr;

double FaceDistZ(double x, double y, double z) {
  TH2F* h = (std::fabs(z - kZmin) <= std::fabs(z - kZmax)) ? g_distZ_lo : g_distZ_hi;
  return h ? h->GetBinContent(h->FindBin(x, y)) : -999999;
}
double FaceDistY(double x, double y, double z) {
  TH2F* h = (std::fabs(y - kYmin) <= std::fabs(y - kYmax)) ? g_distY_lo : g_distY_hi;
  return h ? h->GetBinContent(h->FindBin(x, z)) : -999999;
}

double DriftVelRaw(double E, double T = 87.4) {
  constexpr double T0 = 90.371;
  constexpr double P1 = -0.01481, P2 = -0.0075, P3 = 0.141;
  constexpr double P4 = 12.4, P5 = 1.627, P6 = 0.317;
  return (P1 * (T - T0) + 1.0) * (P3 * E * std::log(1.0 + P4 / E) + P5 * std::pow(E, P6))
         + P2 * (T - T0);
}

double EfieldFromSpeed(double v) {
  double lo = 0.001, hi = 4.0;
  for (int it = 0; it < 200; ++it) {
    const double mid = 0.5*(lo+hi);
    if (DriftVelRaw(mid) < v) lo = mid; else hi = mid;
  }
  return 0.5*(lo+hi);
}

void EfieldFromDistortions(TH3F* dX, TH3F* dY, TH3F* dZ,
                           TH3F* Ex, TH3F* Ey, TH3F* Ez,
                           double wX, int nbinX, int nbinY, int nbinZ) {
  const double V0 = DriftVelRaw(kEfieldNominal);
  auto filled = [](TH3F* h, int ix, int iy, int iz) {
    return h->GetBinContent(ix, iy, iz) > -1e6;
  };
  TAxis* xax = dX->GetXaxis();
  for (int ix = 1; ix <= nbinX; ++ix) {
    const double xc = xax->GetBinCenter(ix);
    if (xc == 0) continue;

    const int jx = (xc > 0) ? ix + 1 : ix - 1;
    if (jx < 1 || jx > nbinX) continue;
    const double xj = xax->GetBinCenter(jx);
    if ((xc > 0) != (xj > 0)) continue;
    const double sDX = xj - xc;
    for (int iy = 1; iy <= nbinY; ++iy)
      for (int iz = 1; iz <= nbinZ; ++iz) {
        if (!filled(dX, ix, iy, iz) || !filled(dX, jx, iy, iz)) continue;
        if (!filled(dY, ix, iy, iz) || !filled(dY, jx, iy, iz)) continue;
        if (!filled(dZ, ix, iy, iz) || !filled(dZ, jx, iy, iz)) continue;
        const double DX = wX;
        const double Rx = sDX + (dX->GetBinContent(jx,iy,iz) - dX->GetBinContent(ix,iy,iz));
        const double Ry =       (dY->GetBinContent(jx,iy,iz) - dY->GetBinContent(ix,iy,iz));
        const double Rz =       (dZ->GetBinContent(jx,iy,iz) - dZ->GetBinContent(ix,iy,iz));
        const double Rmag = std::sqrt(Rx*Rx + Ry*Ry + Rz*Rz);
        if (Rmag < 1e-9 || DX < 1e-9) continue;
        const double vloc = (Rmag / DX) * V0;
        const double Eloc = EfieldFromSpeed(vloc);
        const double ExAbs = Eloc * Rx / Rmag;
        const double EyAbs = Eloc * Ry / Rmag;
        const double EzAbs = Eloc * Rz / Rmag;

        Ex->SetBinContent(ix, iy, iz, (std::fabs(ExAbs) - kEfieldNominal) / kEfieldNominal);
        Ey->SetBinContent(ix, iy, iz, EyAbs / kEfieldNominal);
        Ez->SetBinContent(ix, iy, iz, EzAbs / kEfieldNominal);
      }
  }
}

bool inFid(double x, double y, double z) {
  return x >= -kX1 && x <= kX1 && y >= kYmin && y <= kYmax && z >= kZmin && z <= kZmax;
}

void crossp(const double a[3], const double b[3], double o[3]) {
  o[0] = a[1]*b[2] - a[2]*b[1];
  o[1] = a[2]*b[0] - a[0]*b[2];
  o[2] = a[0]*b[1] - a[1]*b[0];
}
void shortest_int(const double A1[3], const double B1[3],
                  const double A2[3], const double B2[3], double r[7]) {
  double d1[3] = {B1[0]-A1[0], B1[1]-A1[1], B1[2]-A1[2]};
  double d2[3] = {B2[0]-A2[0], B2[1]-A2[1], B2[2]-A2[2]};
  double p2p1[3] = {A2[0]-A1[0], A2[1]-A1[1], A2[2]-A1[2]};
  double n[3];  crossp(d1, d2, n);
  double n2[3]; crossp(d2, n, n2);
  double n1[3]; crossp(d1, n, n1);
  double pv  =  (p2p1[0]*n2[0]+p2p1[1]*n2[1]+p2p1[2]*n2[2]) / (d1[0]*n2[0]+d1[1]*n2[1]+d1[2]*n2[2]);
  double pv1 = -(p2p1[0]*n1[0]+p2p1[1]*n1[1]+p2p1[2]*n1[2]) / (d2[0]*n1[0]+d2[1]*n1[1]+d2[2]*n1[2]);
  r[0]=A1[0]+pv*d1[0]; r[1]=A1[1]+pv*d1[1]; r[2]=A1[2]+pv*d1[2];
  r[3]=A2[0]+pv1*d2[0]; r[4]=A2[1]+pv1*d2[1]; r[5]=A2[2]+pv1*d2[2];
  const double dx=r[0]-r[3], dy=r[1]-r[4], dz=r[2]-r[5];
  r[6]=dx*dx+dy*dy+dz*dz;
}

struct Crosser {
  std::vector<double> t_neg, z_neg, y_neg;
  std::vector<double> t_pos, z_pos, y_pos;
};

}

void data_crossing_pdvd::Loop() {
  if (fChain == 0) return;

  const int nbinX = int((kXmax - kXmin) / kBin);
  const int nbinY = int((kYmax - kYmin) / kBin);
  const int nbinZ = int((kZmax - kZmin) / kBin);

  const double wX = (kXmax - kXmin) / nbinX;
  const double wY = (kYmax - kYmin) / nbinY;
  const double wZ = (kZmax - kZmin) / nbinZ;
  auto xbin = [&](double x) { return int((x - kXmin) / wX); };
  auto ybin = [&](double y) { return int((y - kYmin) / wY); };
  auto zbin = [&](double z) { return int((z - kZmin) / wZ); };

  TH2F* Ydist_pos_fd = new TH2F("Ydist_pos_fd", "Y dist +X (true bin);Z [cm];Y [cm]", nbinZ, kZmin, kZmax, nbinY, kYmin, kYmax);
  TH2F* Ydist_neg_fd = new TH2F("Ydist_neg_fd", "Y dist -X (true bin);Z [cm];Y [cm]", nbinZ, kZmin, kZmax, nbinY, kYmin, kYmax);
  TH2F* Zdist_pos_fd = new TH2F("Zdist_pos_fd", "Z dist +X (true bin);Z [cm];Y [cm]", nbinZ, kZmin, kZmax, nbinY, kYmin, kYmax);
  TH2F* Zdist_neg_fd = new TH2F("Zdist_neg_fd", "Z dist -X (true bin);Z [cm];Y [cm]", nbinZ, kZmin, kZmax, nbinY, kYmin, kYmax);
  TH2F* Xdist_pos_fd = new TH2F("Xdist_pos_fd", "X dist +X (true bin);Z [cm];Y [cm]", nbinZ, kZmin, kZmax, nbinY, kYmin, kYmax);
  TH2F* Xdist_neg_fd = new TH2F("Xdist_neg_fd", "X dist -X (true bin);Z [cm];Y [cm]", nbinZ, kZmin, kZmax, nbinY, kYmin, kYmax);
  TH2F* Ydist_pos_bd = new TH2F("Ydist_pos_bd", "Y dist +X (reco bin);Z [cm];Y [cm]", nbinZ, kZmin, kZmax, nbinY, kYmin, kYmax);
  TH2F* Ydist_neg_bd = new TH2F("Ydist_neg_bd", "Y dist -X (reco bin);Z [cm];Y [cm]", nbinZ, kZmin, kZmax, nbinY, kYmin, kYmax);
  TH2F* Zdist_pos_bd = new TH2F("Zdist_pos_bd", "Z dist +X (reco bin);Z [cm];Y [cm]", nbinZ, kZmin, kZmax, nbinY, kYmin, kYmax);
  TH2F* Zdist_neg_bd = new TH2F("Zdist_neg_bd", "Z dist -X (reco bin);Z [cm];Y [cm]", nbinZ, kZmin, kZmax, nbinY, kYmin, kYmax);
  TH2F* Xdist_pos_bd = new TH2F("Xdist_pos_bd", "X dist +X (reco bin);Z [cm];Y [cm]", nbinZ, kZmin, kZmax, nbinY, kYmin, kYmax);
  TH2F* Xdist_neg_bd = new TH2F("Xdist_neg_bd", "X dist -X (reco bin);Z [cm];Y [cm]", nbinZ, kZmin, kZmax, nbinY, kYmin, kYmax);

  TH3F* XYZ_Xdist = new TH3F("XYZ_Xdist", "3D X dist;X [cm];Y [cm];Z [cm]", nbinX, kXmin, kXmax, nbinY, kYmin, kYmax, nbinZ, kZmin, kZmax);
  TH3F* XYZ_Ydist = new TH3F("XYZ_Ydist", "3D Y dist;X [cm];Y [cm];Z [cm]", nbinX, kXmin, kXmax, nbinY, kYmin, kYmax, nbinZ, kZmin, kZmax);
  TH3F* XYZ_Zdist = new TH3F("XYZ_Zdist", "3D Z dist;X [cm];Y [cm];Z [cm]", nbinX, kXmin, kXmax, nbinY, kYmin, kYmax, nbinZ, kZmin, kZmax);

  TH2F* Ydist_pos_truth = new TH2F("Ydist_pos_truth", "SCE Y dist +X;Z [cm];Y [cm]", nbinZ, kZmin, kZmax, nbinY, kYmin, kYmax);
  TH2F* Ydist_neg_truth = new TH2F("Ydist_neg_truth", "SCE Y dist -X;Z [cm];Y [cm]", nbinZ, kZmin, kZmax, nbinY, kYmin, kYmax);
  TH2F* Zdist_pos_truth = new TH2F("Zdist_pos_truth", "SCE Z dist +X;Z [cm];Y [cm]", nbinZ, kZmin, kZmax, nbinY, kYmin, kYmax);
  TH2F* Zdist_neg_truth = new TH2F("Zdist_neg_truth", "SCE Z dist -X;Z [cm];Y [cm]", nbinZ, kZmin, kZmax, nbinY, kYmin, kYmax);
  TH2F* Xdist_pos_truth = new TH2F("Xdist_pos_truth", "SCE X dist +X;Z [cm];Y [cm]", nbinZ, kZmin, kZmax, nbinY, kYmin, kYmax);
  TH2F* Xdist_neg_truth = new TH2F("Xdist_neg_truth", "SCE X dist -X;Z [cm];Y [cm]", nbinZ, kZmin, kZmax, nbinY, kYmin, kYmax);
  TH3F* XYZ_Xdist_truth = new TH3F("XYZ_Xdist_truth", "SCE 3D X dist;X [cm];Y [cm];Z [cm]", nbinX, kXmin, kXmax, nbinY, kYmin, kYmax, nbinZ, kZmin, kZmax);
  TH3F* XYZ_Ydist_truth = new TH3F("XYZ_Ydist_truth", "SCE 3D Y dist;X [cm];Y [cm];Z [cm]", nbinX, kXmin, kXmax, nbinY, kYmin, kYmax, nbinZ, kZmin, kZmax);
  TH3F* XYZ_Zdist_truth = new TH3F("XYZ_Zdist_truth", "SCE 3D Z dist;X [cm];Y [cm];Z [cm]", nbinX, kXmin, kXmax, nbinY, kYmin, kYmax, nbinZ, kZmin, kZmax);

  TH2F *XY_coverage=new TH2F("XY_coverage","XY coverage;Y coordinate[cm];X [cm]",nbinY,kYmin,kYmax,nbinX,kXmin,kXmax);
  TH2F *XZ_coverage=new TH2F("XZ_coverage","XZ coverage;Z coordinate[cm];X [cm]",nbinZ,kZmin,kZmax,nbinX,kXmin,kXmax);
  TH2F *YZ_coverage=new TH2F("YZ_coverage","YZ coverage;Z [cm];Y [cm]",nbinZ,kZmin,kZmax,nbinY,kYmin,kYmax);

  auto mk3 = [&](const char* n, const char* ti){
    return new TH3F(n, Form("%s;X [cm];Y [cm];Z [cm]", ti),
                    nbinX, kXmin, kXmax, nbinY, kYmin, kYmax, nbinZ, kZmin, kZmax); };
  TH3F* XYZ_Occupancy = mk3("XYZ_Occupancy", "crossing-pair count per voxel");
  TH3F* XYZ_Xstd      = mk3("XYZ_Xstd", "X-distortion STD per voxel [cm]");
  TH3F* XYZ_Ystd      = mk3("XYZ_Ystd", "Y-distortion STD per voxel [cm]");
  TH3F* XYZ_Zstd      = mk3("XYZ_Zstd", "Z-distortion STD per voxel [cm]");

  TH1F* hThetaXZ = new TH1F("trk_thetaxz", "anode-to-anode #theta_{xz};#theta_{xz} [deg];tracks", 180, -180, 180);
  TH1F* hThetaYZ = new TH1F("trk_thetayz", "anode-to-anode #theta_{yz};#theta_{yz} [deg];tracks", 180, -180, 180);
  TH1F* hTheta   = new TH1F("trk_theta",   "anode-to-anode polar #theta;#theta [deg];tracks", 180, 0, 180);
  TH1F* hPhi     = new TH1F("trk_phi",     "anode-to-anode azimuth #phi;#phi [deg];tracks", 180, -180, 180);
  TH1F* hT0proxy = new TH1F("trk_t0proxy", "anode-end peak time (T0 reference);peakT_min [ticks];tracks", 200, 0, 6000);
  TH1F* hDCA     = new TH1F("pair_dca", "crossing-pair distance of closest approach;DCA [cm];pairs", 200, 0, 5);

  auto make2 = [](int a, int b){ return std::vector<std::vector<double> >(a*b); };
  std::vector<std::vector<double> > accYp_fd = make2(nbinZ, nbinY), accYn_fd = make2(nbinZ, nbinY);
  std::vector<std::vector<double> > accZp_fd = make2(nbinZ, nbinY), accZn_fd = make2(nbinZ, nbinY);
  std::vector<std::vector<double> > accXp_fd = make2(nbinZ, nbinY), accXn_fd = make2(nbinZ, nbinY);
  std::vector<std::vector<double> > accYp_bd = make2(nbinZ, nbinY), accYn_bd = make2(nbinZ, nbinY);
  std::vector<std::vector<double> > accZp_bd = make2(nbinZ, nbinY), accZn_bd = make2(nbinZ, nbinY);
  std::vector<std::vector<double> > accXp_bd = make2(nbinZ, nbinY), accXn_bd = make2(nbinZ, nbinY);
  std::vector<std::vector<double> > accX3(nbinX*nbinY*nbinZ), accY3(nbinX*nbinY*nbinZ), accZ3(nbinX*nbinY*nbinZ);

  auto fillCell = [&](std::vector<std::vector<double> >& acc, int zb, int yb, double v) {
    if (zb < 0 || zb >= nbinZ || yb < 0 || yb >= nbinY) return;
    acc[zb*nbinY + yb].push_back(v);
  };

  std::vector<Crosser> crossers;

  const Long64_t nentries = fChain->GetEntries();
  for (Long64_t jentry = 0; jentry < nentries; ++jentry) {
    if (LoadTree(jentry) < 0) break;
    fChain->GetEntry(jentry);
    if (jentry % 10000 == 0) std::cout << jentry << "/" << nentries << std::endl;

    for (int i = 0; i < cross_trks && i < 1000; ++i) {
      if (!trkanodetoanode[i]) continue;

      {
        const double cx = trkstartcosxyz[i][0], cy = trkstartcosxyz[i][1], cz = trkstartcosxyz[i][2];
        if (cx != 0 || cy != 0 || cz != 0) {
          hThetaXZ->Fill(std::atan2(cx, cz) * TMath::RadToDeg());
          hThetaYZ->Fill(std::atan2(cy, cz) * TMath::RadToDeg());
          hTheta->Fill(std::acos(std::max(-1.0, std::min(1.0, cz))) * TMath::RadToDeg());
          hPhi->Fill(std::atan2(cy, cx) * TMath::RadToDeg());
        }
        hT0proxy->Fill(peakT_min[i]);
      }

      Crosser c;
      const int nh = std::min(ntrkhits[i][2], 3000);
      for (int j = 0; j < nh; ++j) {
        const double x = trkhitx[i][2][j];
        const double y = trkhity[i][2][j];
        const double z = trkhitz[i][2][j];
        const double t = trkhitpeakt[i][2][j];
        if (x == 0 && y == 0 && z == 0) continue;
        if (t <= 0) continue;
        if (x > 0) { c.t_pos.push_back(t); c.z_pos.push_back(z); c.y_pos.push_back(y); }
        else       { c.t_neg.push_back(t); c.z_neg.push_back(z); c.y_neg.push_back(y); }
      }
      if (c.t_pos.size() < 5 || c.t_neg.size() < 5) continue;

      auto reorder = [](std::vector<double>& t, std::vector<double>& z,
                        std::vector<double>& y, bool ascending) {
        if ((t.front() < t.back()) != ascending) {
          std::reverse(t.begin(), t.end());
          std::reverse(z.begin(), z.end());
          std::reverse(y.begin(), y.end());
        }
      };
      reorder(c.t_neg, c.z_neg, c.y_neg, true);
      reorder(c.t_pos, c.z_pos, c.y_pos, false);

      const double z0 = c.z_neg.front(), z1 = c.z_pos.back();
      const double y0 = c.y_neg.front(), y1 = c.y_pos.back();

      const double z_true_pos = z0 + (z1 - z0) * (kXcathP - kX0) / (kX1 - kX0);
      const double z_true_neg = z0 + (z1 - z0) * (kXcathN - kX0) / (kX1 - kX0);
      const double y_true_pos = y0 + (y1 - y0) * (kXcathP - kX0) / (kX1 - kX0);
      const double y_true_neg = y0 + (y1 - y0) * (kXcathN - kX0) / (kX1 - kX0);

      const double y_meas_pos = c.y_pos.front(), z_meas_pos = c.z_pos.front();
      const double y_meas_neg = c.y_neg.back(),  z_meas_neg = c.z_neg.back();

      const double x_meas_pos = kX1 - (c.t_pos.front() - c.t_pos.back()) * kTickUs * kVnom;
      const double x_meas_neg = kX0 + (c.t_neg.back()  - c.t_neg.front()) * kTickUs * kVnom;
      const double dX_pos = x_meas_pos - kXcathP, dX_neg = x_meas_neg - kXcathN;
      const double dY_pos = y_meas_pos - y_true_pos, dY_neg = y_meas_neg - y_true_neg;
      const double dZ_pos = z_meas_pos - z_true_pos, dZ_neg = z_meas_neg - z_true_neg;

      if (inFid(kXcathP, y_true_pos, z_true_pos) &&
          std::fabs(dY_pos) < kMaxDist && std::fabs(dZ_pos) < kMaxDist && std::fabs(dX_pos) < kMaxDist) {
        const int zb = zbin(z_true_pos), yb = ybin(y_true_pos);
        fillCell(accYp_fd, zb, yb, dY_pos); fillCell(accZp_fd, zb, yb, dZ_pos); fillCell(accXp_fd, zb, yb, dX_pos);
      }
      if (inFid(kXcathN, y_true_neg, z_true_neg) &&
          std::fabs(dY_neg) < kMaxDist && std::fabs(dZ_neg) < kMaxDist && std::fabs(dX_neg) < kMaxDist) {
        const int zb = zbin(z_true_neg), yb = ybin(y_true_neg);
        fillCell(accYn_fd, zb, yb, dY_neg); fillCell(accZn_fd, zb, yb, dZ_neg); fillCell(accXn_fd, zb, yb, dX_neg);
      }

      if (inFid(kXcathP, y_meas_pos, z_meas_pos) &&
          std::fabs(dY_pos) < kMaxDist && std::fabs(dZ_pos) < kMaxDist && std::fabs(dX_pos) < kMaxDist) {
        const int zb = zbin(z_meas_pos), yb = ybin(y_meas_pos);
        fillCell(accYp_bd, zb, yb, -dY_pos); fillCell(accZp_bd, zb, yb, -dZ_pos); fillCell(accXp_bd, zb, yb, -dX_pos);
      }
      if (inFid(kXcathN, y_meas_neg, z_meas_neg) &&
          std::fabs(dY_neg) < kMaxDist && std::fabs(dZ_neg) < kMaxDist && std::fabs(dX_neg) < kMaxDist) {
        const int zb = zbin(z_meas_neg), yb = ybin(y_meas_neg);
        fillCell(accYn_bd, zb, yb, -dY_neg); fillCell(accZn_bd, zb, yb, -dZ_neg); fillCell(accXn_bd, zb, yb, -dX_neg);
      }

      crossers.push_back(std::move(c));
    }
  }

  auto fill2 = [&](TH2F* h, const std::vector<std::vector<double> >& acc) {
    for (int iz = 0; iz < nbinZ; ++iz)
      for (int iy = 0; iy < nbinY; ++iy) {
        const std::vector<double>& v = acc[iz*nbinY + iy];
        h->SetBinContent(iz+1, iy+1, v.empty() ? -999999 : TMath::Median(v.size(), &v[0]));
      }
  };
  fill2(Ydist_pos_fd, accYp_fd); fill2(Ydist_neg_fd, accYn_fd);
  fill2(Zdist_pos_fd, accZp_fd); fill2(Zdist_neg_fd, accZn_fd);
  fill2(Xdist_pos_fd, accXp_fd); fill2(Xdist_neg_fd, accXn_fd);
  fill2(Ydist_pos_bd, accYp_bd); fill2(Ydist_neg_bd, accYn_bd);
  fill2(Zdist_pos_bd, accZp_bd); fill2(Zdist_neg_bd, accZn_bd);
  fill2(Xdist_pos_bd, accXp_bd); fill2(Xdist_neg_bd, accXn_bd);

  auto patchMap = [&](TH2F* h, TH3F* hpos, TH3F* hneg, double xcath) {
    TGraph2D g;
    for (int iz = 0; iz < nbinZ; ++iz)
      for (int iy = 0; iy < nbinY; ++iy) {
        const double zc = kZmin + (iz + 0.5) * wZ;
        const double yc = kYmin + (iy + 0.5) * wY;
        const double v  = h->GetBinContent(iz+1, iy+1);
        const bool edge = (iz == 0 || iz == nbinZ-1 || iy == 0 || iy == nbinY-1);
        if (v > -1e6) {
          g.SetPoint(g.GetN(), zc, yc, v);
        } else if (edge) {
          const double s = TruthDist(hpos, hneg, xcath, yc, zc);
          h->SetBinContent(iz+1, iy+1, s);
          g.SetPoint(g.GetN(), zc, yc, s);
        }
      }
    for (int iz = 1; iz < nbinZ-1; ++iz)
      for (int iy = 1; iy < nbinY-1; ++iy) {
        if (h->GetBinContent(iz+1, iy+1) > -1e6) continue;
        const double zc = kZmin + (iz + 0.5) * wZ;
        const double yc = kYmin + (iy + 0.5) * wY;
        double v = (g.GetN() >= 3) ? g.Interpolate(zc, yc) : 0.0;
        if (v == 0.0) v = TruthDist(hpos, hneg, xcath, yc, zc);
        h->SetBinContent(iz+1, iy+1, v);
      }
  };
  patchMap(Xdist_pos_fd, g_dx_pos, g_dx_neg, kXcathP); patchMap(Xdist_neg_fd, g_dx_pos, g_dx_neg, kXcathN);
  patchMap(Ydist_pos_fd, g_dy_pos, g_dy_neg, kXcathP); patchMap(Ydist_neg_fd, g_dy_pos, g_dy_neg, kXcathN);
  patchMap(Zdist_pos_fd, g_dz_pos, g_dz_neg, kXcathP); patchMap(Zdist_neg_fd, g_dz_pos, g_dz_neg, kXcathN);
  patchMap(Xdist_pos_bd, g_bx_pos, g_bx_neg, kXcathP); patchMap(Xdist_neg_bd, g_bx_pos, g_bx_neg, kXcathN);
  patchMap(Ydist_pos_bd, g_by_pos, g_by_neg, kXcathP); patchMap(Ydist_neg_bd, g_by_pos, g_by_neg, kXcathN);
  patchMap(Zdist_pos_bd, g_bz_pos, g_bz_neg, kXcathP); patchMap(Zdist_neg_bd, g_bz_pos, g_bz_neg, kXcathN);

  for (int iz = 0; iz < nbinZ; ++iz)
    for (int iy = 0; iy < nbinY; ++iy) {
      const double zc = kZmin + (iz + 0.5) * wZ;
      const double yc = kYmin + (iy + 0.5) * wY;
      Xdist_pos_truth->SetBinContent(iz+1, iy+1, TruthDist(g_dx_pos, g_dx_neg, kXcathP, yc, zc));
      Ydist_pos_truth->SetBinContent(iz+1, iy+1, TruthDist(g_dy_pos, g_dy_neg, kXcathP, yc, zc));
      Zdist_pos_truth->SetBinContent(iz+1, iy+1, TruthDist(g_dz_pos, g_dz_neg, kXcathP, yc, zc));
      Xdist_neg_truth->SetBinContent(iz+1, iy+1, TruthDist(g_dx_pos, g_dx_neg, kXcathN, yc, zc));
      Ydist_neg_truth->SetBinContent(iz+1, iy+1, TruthDist(g_dy_pos, g_dy_neg, kXcathN, yc, zc));
      Zdist_neg_truth->SetBinContent(iz+1, iy+1, TruthDist(g_dz_pos, g_dz_neg, kXcathN, yc, zc));
    }

  const int ntr = crossers.size();
  std::cout << "crossers: " << ntr << std::endl;

  struct CrosserGeom {
    std::vector<double> X, Y, Z;
    double A[3], B[3];
  };
  std::vector<CrosserGeom> geo(ntr);
  for (int i = 0; i < ntr; ++i) {
    const Crosser& c = crossers[i];
    CrosserGeom& g = geo[i];
    const size_t nh = c.t_neg.size() + c.t_pos.size();
    g.X.reserve(nh); g.Y.reserve(nh); g.Z.reserve(nh);
    for (size_t k = 0; k < c.t_neg.size(); ++k) {
      g.X.push_back(kX0 + (c.t_neg[k] - c.t_neg.front()) * kTickUs * kVnom);
      g.Y.push_back(c.y_neg[k]); g.Z.push_back(c.z_neg[k]);
    }
    for (size_t k = 0; k < c.t_pos.size(); ++k) {
      g.X.push_back(kX1 - (c.t_pos[k] - c.t_pos.back()) * kTickUs * kVnom);
      g.Y.push_back(c.y_pos[k]); g.Z.push_back(c.z_pos[k]);
    }
    g.A[0] = kX0; g.A[1] = c.y_neg.front(); g.A[2] = c.z_neg.front();
    g.B[0] = kX1; g.B[1] = c.y_pos.back();  g.B[2] = c.z_pos.back();
  }
  std::vector<Crosser>().swap(crossers);

  for (int t = 0; t < ntr; ++t) {
    if (t % 10000 == 0) std::cout << "bulk " << t << "/" << ntr << std::endl;
    const CrosserGeom& gt = geo[t];
    const std::vector<double>& Xt = gt.X;
    const std::vector<double>& Yt = gt.Y;
    const std::vector<double>& Zt = gt.Z;

    for (int u = t+1; u < ntr; ++u) {
      const CrosserGeom& gu = geo[u];
      double tru[7];
      shortest_int(gt.A, gt.B, gu.A, gu.B, tru);

      if (std::fabs(tru[0]) > kXmax || std::fabs(tru[3]) > kXmax) continue;
      if (!inFid(tru[0], tru[1], tru[2])) continue;
      hDCA->Fill(std::sqrt(std::max(0.0, tru[6])));
      if (tru[6] > 1.0) continue;

      const std::vector<double>& Xu = gu.X;
      const std::vector<double>& Yu = gu.Y;
      const std::vector<double>& Zu = gu.Z;
      double bestSq = 1e18, mid[3] = {0,0,0};
      for (size_t a = 0; a < Xt.size(); ++a) {
        if (std::fabs(Xt[a] - tru[0]) > 100) continue;
        const double xa = Xt[a], ya = Yt[a], za = Zt[a];
        for (size_t b = 0; b < Xu.size(); ++b) {
          const double dx = xa - Xu[b];
          if (std::fabs(dx) > 5.0) continue;
          const double dy = ya - Yu[b], dz = za - Zu[b];
          const double ddSq = dx*dx + dy*dy + dz*dz;
          if (ddSq < bestSq) { bestSq = ddSq; mid[0]=0.5*(xa+Xu[b]); mid[1]=0.5*(ya+Yu[b]); mid[2]=0.5*(za+Zu[b]); }
        }
      }
      if (bestSq >= 4.0) continue;

      const double truMid[3] = {0.5*(tru[0]+tru[3]), 0.5*(tru[1]+tru[4]), 0.5*(tru[2]+tru[5])};
      if (!inFid(truMid[0], truMid[1], truMid[2]) || !inFid(mid[0], mid[1], mid[2])) continue;
      const int ix = xbin(mid[0]), iy = ybin(mid[1]), iz = zbin(mid[2]);
      if (ix < 0 || ix >= nbinX || iy < 0 || iy >= nbinY || iz < 0 || iz >= nbinZ) continue;
      const int idx = (ix*nbinY + iy)*nbinZ + iz;
      if (std::fabs(truMid[0]-mid[0]) < kMaxDist) accX3[idx].push_back(truMid[0]-mid[0]);
      if (std::fabs(truMid[1]-mid[1]) < kMaxDist) accY3[idx].push_back(truMid[1]-mid[1]);
      if (std::fabs(truMid[2]-mid[2]) < kMaxDist) accZ3[idx].push_back(truMid[2]-mid[2]);
      XY_coverage->Fill(mid[1],mid[0]);
      YZ_coverage->Fill(mid[2],mid[1]);
      XZ_coverage->Fill(mid[2],mid[0]);
    }
  }

  auto stdOf = [](const std::vector<double>& v) -> double {
    if (v.size() < 2) return 0.0;
    double m = 0; for (double x : v) m += x; m /= v.size();
    double s = 0; for (double x : v) s += (x - m) * (x - m);
    return std::sqrt(s / (v.size() - 1));
  };
  for (int ix = 0; ix < nbinX; ++ix)
    for (int iy = 0; iy < nbinY; ++iy)
      for (int iz = 0; iz < nbinZ; ++iz) {
        const int idx = (ix*nbinY + iy)*nbinZ + iz;
        XYZ_Xdist->SetBinContent(ix+1, iy+1, iz+1, accX3[idx].size()<1 ? -99999999 : TMath::Median(accX3[idx].size(), &accX3[idx][0]));
        XYZ_Ydist->SetBinContent(ix+1, iy+1, iz+1, accY3[idx].size()<1 ? -99999999 : TMath::Median(accY3[idx].size(), &accY3[idx][0]));
        XYZ_Zdist->SetBinContent(ix+1, iy+1, iz+1, accZ3[idx].size()<1 ? -99999999 : TMath::Median(accZ3[idx].size(), &accZ3[idx][0]));

        XYZ_Occupancy->SetBinContent(ix+1, iy+1, iz+1, (double)accX3[idx].size());
        XYZ_Xstd->SetBinContent(ix+1, iy+1, iz+1, stdOf(accX3[idx]));
        XYZ_Ystd->SetBinContent(ix+1, iy+1, iz+1, stdOf(accY3[idx]));
        XYZ_Zstd->SetBinContent(ix+1, iy+1, iz+1, stdOf(accZ3[idx]));
        const double xc = kXmin + (ix + 0.5) * wX;
        const double yc = kYmin + (iy + 0.5) * wY;
        const double zc = kZmin + (iz + 0.5) * wZ;
        XYZ_Xdist_truth->SetBinContent(ix+1, iy+1, iz+1, TruthDist(g_dx_pos, g_dx_neg, xc, yc, zc));
        XYZ_Ydist_truth->SetBinContent(ix+1, iy+1, iz+1, TruthDist(g_dy_pos, g_dy_neg, xc, yc, zc));
        XYZ_Zdist_truth->SetBinContent(ix+1, iy+1, iz+1, TruthDist(g_dz_pos, g_dz_neg, xc, yc, zc));
      }

  auto ringOf = [&](int iy, int iz){
    return std::min(std::min(iy, nbinY-1-iy), std::min(iz, nbinZ-1-iz));
  };
  auto fillBoundary = [&](TH3F* h, int comp) {

    std::vector<double> AX, AY, AZ, AV;
    if (comp == 0) {
      for (int ix = 1; ix <= nbinX; ++ix) {
        const double xa = h->GetXaxis()->GetBinCenter(ix);
        for (int iy = 0; iy < nbinY; ++iy)
          for (int iz = 0; iz < nbinZ; ++iz) {
            if (ringOf(iy, iz) < 2) continue;
            const double v = h->GetBinContent(ix, iy+1, iz+1);
            if (v > -1e6) {
              AX.push_back(xa);
              AY.push_back(h->GetYaxis()->GetBinCenter(iy+1));
              AZ.push_back(h->GetZaxis()->GetBinCenter(iz+1));
              AV.push_back(v);
            }
          }
      }
    }
    auto idw3 = [&](double xq, double yq, double zq)->double{
      double sw = 0.0, sv = 0.0;
      for (size_t k = 0; k < AV.size(); ++k) {
        const double dx = AX[k]-xq, dy = AY[k]-yq, dz = AZ[k]-zq;
        const double d2 = dx*dx + dy*dy + dz*dz;
        if (d2 < 1e-9) return AV[k];
        const double w = 1.0/(d2*d2);
        sw += w; sv += w*AV[k];
      }
      return sw > 0.0 ? sv/sw : 0.0;
    };
    for (int ix = 1; ix <= nbinX; ++ix) {
      const double xc = h->GetXaxis()->GetBinCenter(ix);
      auto isPrimary = [&](int iy, int iz){
        return (comp == 2 && (iz == 0 || iz == nbinZ-1)) ||
               (comp == 1 && (iy == 0 || iy == nbinY-1));
      };
      auto surfVal = [&](int iy, int iz)->double{
        const double yc = h->GetYaxis()->GetBinCenter(iy+1);
        const double zc = h->GetZaxis()->GetBinCenter(iz+1);
        const double s = (comp == 2) ? FaceDistZ(xc, yc, zc) : FaceDistY(xc, yc, zc);
        return (s > -1e6) ? s : -1e9;
      };

      TGraph2D g;
      std::vector<double> ay, az, av;
      for (int iy = 0; iy < nbinY; ++iy)
        for (int iz = 0; iz < nbinZ; ++iz) {
          const double yc = h->GetYaxis()->GetBinCenter(iy+1);
          const double zc = h->GetZaxis()->GetBinCenter(iz+1);
          double val = -1e30;
          if (ringOf(iy, iz) == 0 && isPrimary(iy, iz)) {
            const double s = surfVal(iy, iz);
            if (s > -1e8) val = s;
          } else if (ringOf(iy, iz) >= 2) {
            const double v = h->GetBinContent(ix, iy+1, iz+1);
            if (v > -1e6) val = v;
          }
          if (val > -1e8) { g.SetPoint(g.GetN(), yc, zc, val); ay.push_back(yc); az.push_back(zc); av.push_back(val); }
        }

      auto idw = [&](double yc, double zc)->double{
        std::vector<std::pair<double,double> > d;
        d.reserve(av.size());
        for (size_t k = 0; k < av.size(); ++k)
          d.push_back(std::make_pair((ay[k]-yc)*(ay[k]-yc) + (az[k]-zc)*(az[k]-zc), av[k]));
        std::sort(d.begin(), d.end());
        const int K = std::min<int>(6, (int)d.size());
        if (d[0].first < 1e-9) return d[0].second;
        double sw = 0.0, sv = 0.0;
        for (int k = 0; k < K; ++k) {
          const double w = 1.0 / (d[k].first * d[k].first);
          sw += w; sv += w * d[k].second;
        }
        return sw > 0.0 ? sv / sw : d[0].second;
      };

      auto nearD2 = [&](double yc, double zc)->double{
        double best = 1e30;
        for (size_t k = 0; k < av.size(); ++k) {
          const double d = (ay[k]-yc)*(ay[k]-yc) + (az[k]-zc)*(az[k]-zc);
          if (d < best) best = d;
        }
        return best;
      };
      const double kHull2 = (2.0*kBin) * (2.0*kBin);

      for (int iy = 0; iy < nbinY; ++iy)
        for (int iz = 0; iz < nbinZ; ++iz) {
          const int r = ringOf(iy, iz);
          if (r == 0 && isPrimary(iy, iz)) {
            const double s = surfVal(iy, iz);
            if (s > -1e8) { h->SetBinContent(ix, iy+1, iz+1, s); continue; }
          }
          if (r >= 2 && h->GetBinContent(ix, iy+1, iz+1) > -1e6) continue;

          if ((comp == 0 ? AV.empty() : av.empty())) continue;
          const double yc = h->GetYaxis()->GetBinCenter(iy+1);
          const double zc = h->GetZaxis()->GetBinCenter(iz+1);

          double v;
          if (comp == 0) {
            v = idw3(xc, yc, zc);
          } else {
            v = -1e30;
            if (g.GetN() >= 3 && nearD2(yc, zc) <= kHull2) v = g.Interpolate(yc, zc);
            if (v <= -1e29 || v == 0.0) v = idw(yc, zc);
          }
          h->SetBinContent(ix, iy+1, iz+1, v);
        }
    }
  };
  fillBoundary(XYZ_Xdist, 0);
  fillBoundary(XYZ_Ydist, 1);
  fillBoundary(XYZ_Zdist, 2);

  auto makeReco = [&](const char* name, TH3F* myE) -> TH3F* {
    TH3F* out = (TH3F*)myE->Clone(name);
    out->SetDirectory(nullptr);
    TH3F* sce = g_sce ? (TH3F*)g_sce->Get(name) : nullptr;
    out->SetTitle(sce ? sce->GetTitle() : name);
    return out;
  };

  TH3F* Bkwd_X_Pos = makeReco("RecoBkwd_Displacement_X_Pos", XYZ_Xdist);
  TH3F* Bkwd_X_Neg = makeReco("RecoBkwd_Displacement_X_Neg", XYZ_Xdist);
  TH3F* Bkwd_Y_Pos = makeReco("RecoBkwd_Displacement_Y_Pos", XYZ_Ydist);
  TH3F* Bkwd_Y_Neg = makeReco("RecoBkwd_Displacement_Y_Neg", XYZ_Ydist);
  TH3F* Bkwd_Z_Pos = makeReco("RecoBkwd_Displacement_Z_Pos", XYZ_Zdist);
  TH3F* Bkwd_Z_Neg = makeReco("RecoBkwd_Displacement_Z_Neg", XYZ_Zdist);

  TFile fout("data_crossing_pdvd.root", "RECREATE");

  Ydist_pos_fd->Write(); Ydist_neg_fd->Write(); Zdist_pos_fd->Write(); Zdist_neg_fd->Write();
  Xdist_pos_fd->Write(); Xdist_neg_fd->Write();
  Ydist_pos_bd->Write(); Ydist_neg_bd->Write(); Zdist_pos_bd->Write(); Zdist_neg_bd->Write();
  Xdist_pos_bd->Write(); Xdist_neg_bd->Write();
  XYZ_Xdist->Write(); XYZ_Ydist->Write(); XYZ_Zdist->Write();

  Bkwd_X_Pos->Write(); Bkwd_X_Neg->Write();
  Bkwd_Y_Pos->Write(); Bkwd_Y_Neg->Write();
  Bkwd_Z_Pos->Write(); Bkwd_Z_Neg->Write();
  XZ_coverage->Write(); YZ_coverage->Write(); XY_coverage->Write();

  XYZ_Occupancy->Write();
  XYZ_Xstd->Write(); XYZ_Ystd->Write(); XYZ_Zstd->Write();
  hThetaXZ->Write(); hThetaYZ->Write(); hTheta->Write(); hPhi->Write();
  hT0proxy->Write(); hDCA->Write();

  Ydist_pos_truth->Write(); Ydist_neg_truth->Write();
  Zdist_pos_truth->Write(); Zdist_neg_truth->Write();
  Xdist_pos_truth->Write(); Xdist_neg_truth->Write();
  XYZ_Xdist_truth->Write(); XYZ_Ydist_truth->Write(); XYZ_Zdist_truth->Write();
  fout.Close();
  std::cout << "data_crossing_pdvd finished" << std::endl;
}
