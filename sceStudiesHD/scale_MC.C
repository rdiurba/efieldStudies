#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TKey.h>
#include <TList.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace sce_bulk {

constexpr int    kSmoothXPasses = 2;  // 1-2-1 passes along X on the scaled RecoBkwd dX map (E_x comes from its X slope)
constexpr int    kSlopeMaxPasses = 50;  // max rounds of the X->Y->Z slope check against the MC
constexpr double kRegSim =  1.0;   // [cm] ratio regularization scale (see above): smooth -> 1 as |MC| -> 0

// ---- E field from RecoBkwd displacements: same as data_crossing_pdhd.C ----
constexpr double kEfieldNominal = 0.488;

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

void EfieldFromDistortions(TH3* dX, TH3* dY, TH3* dZ,
                           TH3* Ex, TH3* Ey, TH3* Ez,
                           double wX, int nbinX, int nbinY, int nbinZ) {
  const double V0 = DriftVelRaw(kEfieldNominal);
  auto filled = [](TH3* h, int ix, int iy, int iz) {
    return h->GetBinContent(ix, iy, iz) > -1e6;
  };
  TAxis* xax = dX->GetXaxis();
  for (int ix = 1; ix <= nbinX; ++ix) {
    const double xc = xax->GetBinCenter(ix);
    if (xc == 0) continue;

    const int step = (xc > 0) ? 1 : -1;
    int ia = ix, ib = ix + step;
    if (ib < 1 || ib > nbinX) { ia = ix - step; ib = ix; }   // outermost bin: use it and its inner neighbour
    if (ia < 1 || ia > nbinX) continue;
    const double xa = xax->GetBinCenter(ia), xb = xax->GetBinCenter(ib);
    if ((xa > 0) != (xb > 0)) continue;
    const double sDX = xb - xa;
    for (int iy = 1; iy <= nbinY; ++iy)
      for (int iz = 1; iz <= nbinZ; ++iz) {
        if (!filled(dX, ia, iy, iz) || !filled(dX, ib, iy, iz)) continue;
        if (!filled(dY, ia, iy, iz) || !filled(dY, ib, iy, iz)) continue;
        if (!filled(dZ, ia, iy, iz) || !filled(dZ, ib, iy, iz)) continue;
        const double DX = wX;
        const double Rx = sDX + (dX->GetBinContent(ib,iy,iz) - dX->GetBinContent(ia,iy,iz));
        const double Ry =       (dY->GetBinContent(ib,iy,iz) - dY->GetBinContent(ia,iy,iz));
        const double Rz =       (dZ->GetBinContent(ib,iy,iz) - dZ->GetBinContent(ia,iy,iz));
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

}  // namespace sce_bulk

void scale_MC(
    const char* surffile = "surface_distortions_pdhd.root",
    const char* outfile  = "SCE_PDHD_DataDriven_surface.root",
    const char* simfile  = "/pnfs/dune/persistent/stash/ProtoDUNE/HD/SpaceCharge/Sim_NoFluidFlow/v2/SCE_PDHD_Sim_NoFluidFlow_E500_v2.root",
    bool        mcOnly   = false)   // test: all scale factors = 1 (pure MC), E field still recomputed
{
  using namespace sce_bulk;

  TFile* fs = TFile::Open(surffile, "READ");
  if (!fs || fs->IsZombie()) { std::cout << "could not open " << surffile << std::endl; return; }
  TH2* zlo = dynamic_cast<TH2*>(fs->Get("distZ_zlo"));
  TH2* zhi = dynamic_cast<TH2*>(fs->Get("distZ_zhi"));
  TH2* ylo = dynamic_cast<TH2*>(fs->Get("distY_ylo"));
  TH2* yhi = dynamic_cast<TH2*>(fs->Get("distY_yhi"));
  if (!zlo || !zhi || !ylo || !yhi) { std::cout << "missing distZ_zlo/zhi or distY_ylo/yhi in " << surffile << std::endl; return; }

  TFile* fsim = TFile::Open(simfile, "READ");
  if (!fsim || fsim->IsZombie()) { std::cout << "could not open " << simfile << std::endl; return; }
  TH3* sx_p = dynamic_cast<TH3*>(fsim->Get("RecoBkwd_Displacement_X_Pos"));
  TH3* sx_n = dynamic_cast<TH3*>(fsim->Get("RecoBkwd_Displacement_X_Neg"));
  TH3* sy_p = dynamic_cast<TH3*>(fsim->Get("RecoBkwd_Displacement_Y_Pos"));
  TH3* sy_n = dynamic_cast<TH3*>(fsim->Get("RecoBkwd_Displacement_Y_Neg"));
  TH3* sz_p = dynamic_cast<TH3*>(fsim->Get("RecoBkwd_Displacement_Z_Pos"));
  TH3* sz_n = dynamic_cast<TH3*>(fsim->Get("RecoBkwd_Displacement_Z_Neg"));
  if (!sx_p || !sx_n || !sy_p || !sy_n || !sz_p || !sz_n) { std::cout << "missing RecoBkwd_Displacement_* in " << simfile << std::endl; return; }

  // side index: 0 = Neg (x < 0), 1 = Pos (x >= 0). The X, Y and Z maps of a side
  // must share one binning: the scale factors are made on exactly those bins.
  TH3* SX[2] = {sx_n, sx_p};
  TH3* SY[2] = {sy_n, sy_p};
  TH3* SZ[2] = {sz_n, sz_p};
  const char* sideName[2] = {"Neg", "Pos"};
  for (int s = 0; s < 2; ++s)
    for (TH3* h : {SX[s], SY[s]})
      if (h->GetNbinsX() != SZ[s]->GetNbinsX() || h->GetNbinsY() != SZ[s]->GetNbinsY() || h->GetNbinsZ() != SZ[s]->GetNbinsZ()) {
        std::cout << "RecoBkwd X/Y/Z " << sideName[s] << " maps have different binnings" << std::endl; return;
      }

  // surface map value at (a,b), interpolated; X is kept on one side of the
  // cathode (X = 0 is a bin edge of the surface maps) so the two drift volumes
  // are never mixed
  auto face2 = [](TH2* h, double a, double b, int s) {
    TAxis* xa = h->GetXaxis(); TAxis* ya = h->GetYaxis();
    const int c = xa->FindFixBin(0.0);
    const double lo = s ? xa->GetBinCenter(c) : xa->GetBinCenter(1);
    const double hi = s ? xa->GetBinCenter(xa->GetNbins()) : xa->GetBinCenter(c - 1);
    a = std::min(std::max(a, lo), hi);
    b = std::min(std::max(b, ya->GetBinCenter(1)), ya->GetBinCenter(ya->GetNbins()));
    return h->Interpolate(a, b);
  };
  // data / MC ratio on a face bin (1 = MC unchanged where the MC displacement is ~0)
  auto ratio = [&](double meas, double sim) {
    if (mcOnly) return 1.0;
    return (meas * sim + kRegSim * kRegSim) / (sim * sim + kRegSim * kRegSim);
  };

  // face ratio maps on the MC binning, one set per side: measured surface
  // distortion / MC in the first row of MC bins at that face
  auto mk2 = [&](const char* name, const char* title, TAxis* a, TAxis* b, const char* lab) {
    TH2F* h = new TH2F(name, Form("%s;X [cm];%s", title, lab), a->GetNbins(), a->GetXmin(), a->GetXmax(), b->GetNbins(), b->GetXmin(), b->GetXmax());
    if (a->GetXbins()->GetSize()) h->GetXaxis()->Set(a->GetNbins(), a->GetXbins()->GetArray());
    if (b->GetXbins()->GetSize()) h->GetYaxis()->Set(b->GetNbins(), b->GetXbins()->GetArray());
    h->SetDirectory(nullptr);
    return h;
  };
  TH2F* rZlo[2]; TH2F* rZhi[2]; TH2F* rYlo[2]; TH2F* rYhi[2];
  for (int s = 0; s < 2; ++s) {
    TAxis* ax = SZ[s]->GetXaxis(); TAxis* ay = SZ[s]->GetYaxis(); TAxis* az = SZ[s]->GetZaxis();
    rZlo[s] = mk2(Form("ratioZ_zlo_%s", sideName[s]), "data/MC #DeltaZ, upstream face",   ax, ay, "Y [cm]");
    rZhi[s] = mk2(Form("ratioZ_zhi_%s", sideName[s]), "data/MC #DeltaZ, downstream face", ax, ay, "Y [cm]");
    rYlo[s] = mk2(Form("ratioY_ylo_%s", sideName[s]), "data/MC #DeltaY, bottom face",     ax, az, "Z [cm]");
    rYhi[s] = mk2(Form("ratioY_yhi_%s", sideName[s]), "data/MC #DeltaY, top face",        ax, az, "Z [cm]");
    for (int i = 1; i <= ax->GetNbins(); ++i) {
      const double x = ax->GetBinCenter(i);
      for (int j = 1; j <= ay->GetNbins(); ++j) {
        const double y = ay->GetBinCenter(j);
        const double mlo = face2(zlo, x, y, s), mhi = face2(zhi, x, y, s);
        rZlo[s]->SetBinContent(i, j, ratio(mlo, SZ[s]->GetBinContent(i, j, 1)));
        rZhi[s]->SetBinContent(i, j, ratio(mhi, SZ[s]->GetBinContent(i, j, az->GetNbins())));
      }
      for (int k = 1; k <= az->GetNbins(); ++k) {
        const double z = az->GetBinCenter(k);
        const double mlo = face2(ylo, x, z, s), mhi = face2(yhi, x, z, s);
        rYlo[s]->SetBinContent(i, k, ratio(mlo, SY[s]->GetBinContent(i, 1, k)));
        rYhi[s]->SetBinContent(i, k, ratio(mhi, SY[s]->GetBinContent(i, ay->GetNbins(), k)));
      }
    }
  }

  // scale factor on MC bin (i,j,k): Z maps use rZ (linear in z between the first
  // and last Z bin centres), Y maps use rY (linear in y between the first and last
  // Y bin centres), X maps use rX
  auto frac = [](TAxis* a, int i) {
    const int n = a->GetNbins();
    return n > 1 ? (a->GetBinCenter(i) - a->GetBinCenter(1)) / (a->GetBinCenter(n) - a->GetBinCenter(1)) : 0.0;
  };
  auto scaleZ = [&](int s, int i, int j, int k) {
    const double f = frac(SZ[s]->GetZaxis(), k);
    return (1.0 - f) * rZlo[s]->GetBinContent(i, j) + f * rZhi[s]->GetBinContent(i, j);
  };
  auto scaleY = [&](int s, int i, int j, int k) {
    const double f = frac(SZ[s]->GetYaxis(), j);
    return (1.0 - f) * rYlo[s]->GetBinContent(i, k) + f * rYhi[s]->GetBinContent(i, k);
  };
  // X: bin-by-bin excess/deficit of the transverse (Y,Z) distortion at this voxel,
  //   rX = |(rY dY_MC, rZ dZ_MC)| / |(dY_MC, dZ_MC)|
  // with rY, rZ the Y/Z ratios above and dY_MC, dZ_MC the MC RecoBkwd displacements
  // in the same bin; the mean of rY and rZ where the MC transverse displacement is ~0.
  auto scaleX = [&](int s, int i, int j, int k) {
    if (mcOnly) return 1.0;
    const double sy = SY[s]->GetBinContent(i, j, k), sz = SZ[s]->GetBinContent(i, j, k);
    const double ry = scaleY(s, i, j, k), rz = scaleZ(s, i, j, k);
    const double mc = std::sqrt(sy*sy + sz*sz);
    const double scaled = std::sqrt((ry*sy)*(ry*sy) + (rz*sz)*(rz*sz));   // = rX * mc
    return (scaled * mc + kRegSim * kRegSim * 0.5 * (ry + rz)) / (mc * mc + kRegSim * kRegSim);
  };

  // ---- output: the MC file with every RecoBkwd displacement map scaled ----------
  // Every key of simfile is copied. RecoBkwd_Displacement_{X,Y,Z}_{Pos,Neg} are
  // scaled bin by bin on their own binning; everything else is copied unchanged.
  TFile fout(outfile, "RECREATE");
  int nScaled = 0, nCopied = 0;
  std::map<std::string, TH3*> bkwd;                          // scaled RecoBkwd maps by name
  TIter next(fsim->GetListOfKeys());
  while (TKey* key = (TKey*)next()) {
    TObject* obj = key->ReadObj();
    TString name = obj->GetName();
    char comp = 0;
    if (name.Contains("Displacement_X")) comp = 'X';
    else if (name.Contains("Displacement_Y")) comp = 'Y';
    else if (name.Contains("Displacement_Z")) comp = 'Z';
    TH3* h3 = dynamic_cast<TH3*>(obj);
    fout.cd();
    if (h3 && comp && name.BeginsWith("RecoBkwd")) {
      const int s = name.EndsWith("_Pos") ? 1 : 0;
      TH3* out = (TH3*)h3->Clone(name);
      out->SetDirectory(&fout);
      TAxis* ax = out->GetXaxis(); TAxis* ay = out->GetYaxis(); TAxis* az = out->GetZaxis();
      // anode bin: the X bin farthest from the cathode. No drift there, so all
      // three displacement components are set to 0 in it.
      const int iAnode = (std::fabs(ax->GetBinCenter(ax->GetNbins())) >= std::fabs(ax->GetBinCenter(1)))
                         ? ax->GetNbins() : 1;
      for (int i = 1; i <= ax->GetNbins(); ++i)
        for (int j = 1; j <= ay->GetNbins(); ++j)
          for (int k = 1; k <= az->GetNbins(); ++k) {
            const double sc = (comp == 'X') ? scaleX(s, i, j, k) : (comp == 'Y') ? scaleY(s, i, j, k) : scaleZ(s, i, j, k);
            out->SetBinContent(i, j, k, (i == iAnode) ? 0.0 : sc * h3->GetBinContent(i, j, k));
          }
      // dX only: smooth along X (1-2-1, kSmoothXPasses times) so its slope, i.e.
      // E_x, goes smoothly from the cathode to the anode. Each Pos/Neg map is one
      // side of the cathode, so nothing is mixed across it; the end bins use the
      // neighbours they have.
      if (comp == 'X') {
        const int nx = ax->GetNbins(), ny = ay->GetNbins(), nz = az->GetNbins();
        std::vector<double> row(nx + 2);
        for (int p = 0; p < kSmoothXPasses; ++p)
          for (int j = 1; j <= ny; ++j)
            for (int k = 1; k <= nz; ++k) {
              for (int i = 1; i <= nx; ++i) row[i] = out->GetBinContent(i, j, k);
              for (int i = 1; i <= nx; ++i) {
                double sum = 2.0 * row[i], w = 2.0;
                if (i > 1)  { sum += row[i-1]; w += 1.0; }
                if (i < nx) { sum += row[i+1]; w += 1.0; }
                out->SetBinContent(i, j, k, sum / w);
              }
            }
      }
      // Slope check against the MC, along each coordinate (X, then Y, then Z),
      // repeated until no bin changes (at most kSlopeMaxPasses rounds). Where the
      // MC is locally monotonic but the scaled map's slope points the other way,
      // the bin is replaced by the average of its two neighbours along that
      // coordinate. Edge bins (one neighbour) take that neighbour's value. The
      // anode bin stays 0.
      {
        const int nb[3] = {ax->GetNbins(), ay->GetNbins(), az->GetNbins()};
        auto sgn = [](double v) { return (v > 0) - (v < 0); };
        int changed = 1;
        for (int round = 0; round < kSlopeMaxPasses && changed; ++round) {
          changed = 0;
          for (int axis = 0; axis < 3; ++axis) {
            TH3* src = (TH3*)out->Clone("slopefix_src");   // neighbours from before this pass
            src->SetDirectory(nullptr);
            for (int i = 1; i <= nb[0]; ++i)
              for (int j = 1; j <= nb[1]; ++j)
                for (int k = 1; k <= nb[2]; ++k) {
                  if (i == iAnode) continue;
                  const int c[3] = {i, j, k};
                  const bool hasLo = c[axis] > 1, hasHi = c[axis] < nb[axis];
                  int lo[3] = {i, j, k}, hi[3] = {i, j, k};
                  lo[axis] -= 1; hi[axis] += 1;
                  const double m1 = h3->GetBinContent(i, j, k);
                  const double d1 = src->GetBinContent(i, j, k);
                  const int sL = hasLo ? sgn(m1 - h3->GetBinContent(lo[0], lo[1], lo[2])) : 0;
                  const int sR = hasHi ? sgn(h3->GetBinContent(hi[0], hi[1], hi[2]) - m1) : 0;
                  const double d0 = hasLo ? src->GetBinContent(lo[0], lo[1], lo[2]) : 0.0;
                  const double d2 = hasHi ? src->GetBinContent(hi[0], hi[1], hi[2]) : 0.0;
                  bool wrong = false; double v = d1;
                  if (hasLo && hasHi) {
                    if (sL == 0 || sL != sR) continue;                // MC not monotonic here
                    wrong = sgn(d1 - d0) == -sL || sgn(d2 - d1) == -sR;
                    v = 0.5 * (d0 + d2);
                  } else if (hasLo) {
                    if (sL == 0) continue;
                    wrong = sgn(d1 - d0) == -sL;
                    v = d0;
                  } else if (hasHi) {
                    if (sR == 0) continue;
                    wrong = sgn(d2 - d1) == -sR;
                    v = d2;
                  }
                  if (wrong && v != d1) { out->SetBinContent(i, j, k, v); ++changed; }
                }
            delete src;
          }
        }
      }
      out->Write(name, TObject::kOverwrite);
      bkwd[std::string(name.Data())] = out;
      ++nScaled;
    } else {
      obj->Write(name, TObject::kOverwrite);
      ++nCopied;
    }
  }
  std::cout << "scaled " << nScaled << " RecoBkwd displacement maps, copied " << nCopied
            << " other objects unchanged" << std::endl;

  // ---- E field from the scaled RecoBkwd maps (data_crossing_pdhd.C method) ----
  // Written as Reco_ElecField_{X,Y,Z}_{Pos,Neg} in the MC file's convention
  // (X: (|E_x| - E0)/E0, Y/Z: E/E0), replacing the copied MC E-field maps.
  for (const char* side : {"Pos", "Neg"}) {
    TH3* bx = bkwd[Form("RecoBkwd_Displacement_X_%s", side)];
    TH3* by = bkwd[Form("RecoBkwd_Displacement_Y_%s", side)];
    TH3* bz = bkwd[Form("RecoBkwd_Displacement_Z_%s", side)];
    if (!bx || !by || !bz) { std::cout << "no scaled RecoBkwd X/Y/Z " << side << " maps: E field not made" << std::endl; continue; }
    TH3* e[3];
    const char* comp3[3] = {"X", "Y", "Z"};
    for (int a = 0; a < 3; ++a) {
      e[a] = (TH3*)bx->Clone(Form("Reco_ElecField_%s_%s", comp3[a], side));
      e[a]->Reset();
      e[a]->SetTitle(Form("Reco_ElecField_%s_%s (from surface-scaled RecoBkwd)", comp3[a], side));
      e[a]->SetDirectory(&fout);
    }
    EfieldFromDistortions(bx, by, bz, e[0], e[1], e[2],
                          bx->GetXaxis()->GetBinWidth(1),
                          bx->GetNbinsX(), bx->GetNbinsY(), bx->GetNbinsZ());
    for (int a = 0; a < 3; ++a) e[a]->Write(e[a]->GetName(), TObject::kOverwrite);
    std::cout << "E field " << side << ": Reco_ElecField_{X,Y,Z}_" << side << std::endl;
  }

  fout.cd();
  for (int s = 0; s < 2; ++s) { rZlo[s]->Write(); rZhi[s]->Write(); rYlo[s]->Write(); rYhi[s]->Write(); }
  fout.Close();
  fsim->Close();
  fs->Close();
  std::cout << "wrote " << outfile << std::endl;
}

