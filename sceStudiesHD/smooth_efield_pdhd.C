#include <TFile.h>
#include <TH1.h>
#include <TH3.h>
#include <TMath.h>
#include <TString.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace {

// ---- ProtoDUNE-HD reco constants (match data_crossing_pdhd) -----------------
constexpr double kEfieldNominal = 0.488;    // nominal drift field [kV/cm]
constexpr double kPeakFactor    = 1.5;    // despike: |v| must beat each neighbour by this
constexpr int    kSmoothPasses  = 0;      // gentle 3x3x3 low-pass passes on the displacement before
                                          // the field (the SCE map is physically smooth; this removes
                                          // the bin-scale scatter the derivative amplifies). 0 disables.
constexpr int    kWallRings     = 3;      // 3rd-ring fix: interior starts at this Y/Z ring; the
                                          // intermediate rings 1..kWallRings-1 are interpolated from
                                          // the surface (ring 0) and the interior (0 disables)

// Track-validation histograms + per-voxel diagnostics written by
// data_crossing_pdhd, carried through into the smoothed output unchanged.
const char* const kValidationNames[] = {
  "trk_thetaxz", "trk_thetayz", "trk_theta", "trk_phi",
  "trk_t0proxy", "pair_dca",
  "XYZ_Occupancy", "XYZ_Xstd", "XYZ_Ystd", "XYZ_Zstd"
};

// SCE MC input map, only for the output names/titles.
TFile* g_sce = TFile::Open("/pnfs/dune/persistent/stash/ProtoDUNE/HD/SpaceCharge/Sim_NoFluidFlow/v2/SCE_PDHD_Sim_NoFluidFlow_E500_v2.root");

// LArSoft drift velocity: v [mm/us] vs E [kV/cm] at T [K]. Exact port of
// eField_2.py _drift_raw: piecewise linear / ICARUS / ICARUS-Walkowiak blend /
// Walkowiak (no fudge factor applied, as in the script).
double DriftVelRaw(double efield, double T = 87.4) {
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

// invert DriftVelRaw: speed [mm/us] -> E [kV/cm] by bisection (monotonic).
double EfieldFromSpeed(double v) {
  double lo = 0.001, hi = 4.0;
  for (int it = 0; it < 200; ++it) {
    const double mid = 0.5*(lo+hi);
    if (DriftVelRaw(mid) < v) lo = mid; else hi = mid;
  }
  return 0.5*(lo+hi);
}

// Reconstruct the relative E-field distortion from the 3D displacement maps by
// the field-line step, taken from uB paper
// Stored relative to nominal: Ex map = (|Ex|-E0)/E0 , Ey = Ey/E0 , Ez = Ez/E0.
void EfieldFromDistortions(TH3F* dX, TH3F* dY, TH3F* dZ,
                           TH3F* Ex, TH3F* Ey, TH3F* Ez,
                           double /*wX*/, int nbinX, int nbinY, int nbinZ) {
  const double V0 = DriftVelRaw(kEfieldNominal);
  auto bc     = [](TH3F* h, int ix, int iy, int iz) { return h->GetBinContent(ix, iy, iz); };
  auto filled = [](TH3F* h, int ix, int iy, int iz) { return h->GetBinContent(ix, iy, iz) > -1e6; };
  TAxis* xax = dX->GetXaxis();
  for (int ix = 1; ix <= nbinX; ++ix) {
    const double xc = xax->GetBinCenter(ix);
    if (xc == 0) continue;                             // cathode bin: field undefined
    const int jx = (xc > 0) ? ix + 1 : ix - 1;         // next bin toward the anode
    if (jx < 1 || jx > nbinX) continue;
    const double xj = xax->GetBinCenter(jx);
    if ((xc > 0) != (xj > 0)) continue;                // never cross the cathode
    const double sDX = xj - xc;                        // signed reco step; DX = |sDX| = bin width
    const double DX  = std::fabs(sDX);
    for (int iy = 1; iy <= nbinY; ++iy)
      for (int iz = 1; iz <= nbinZ; ++iz) {
        if (!filled(dX,ix,iy,iz) || !filled(dX,jx,iy,iz)) continue;
        if (!filled(dY,ix,iy,iz) || !filled(dY,jx,iy,iz)) continue;
        if (!filled(dZ,ix,iy,iz) || !filled(dZ,jx,iy,iz)) continue;
        const double Rx = sDX + (bc(dX,jx,iy,iz) - bc(dX,ix,iy,iz));
        const double Ry =       (bc(dY,jx,iy,iz) - bc(dY,ix,iy,iz));
        const double Rz =       (bc(dZ,jx,iy,iz) - bc(dZ,ix,iy,iz));
        const double Rmag = std::sqrt(Rx*Rx + Ry*Ry + Rz*Rz);
        if (Rmag < 1e-9 || DX < 1e-9) continue;
        const double Eloc = EfieldFromSpeed((Rmag / DX) * V0);   // Eq 6.2 then invert v(E)
        Ex->SetBinContent(ix, iy, iz, (std::fabs(Eloc*Rx/Rmag) - kEfieldNominal) / kEfieldNominal);
        Ey->SetBinContent(ix, iy, iz, (Eloc*Ry/Rmag) / kEfieldNominal);
        Ez->SetBinContent(ix, iy, iz, (Eloc*Rz/Rmag) / kEfieldNominal);
      }
  }
}

// Third-ring fix: the INTERMEDIATE Y/Z wall rings of the reconstructed field are
// unreliable -- the displacement there is low-stat and the single-bin difference
// amplifies it into the near-wall "third ring" bands. Interpolate rings
// 1..kWallRings-1 from BOTH the surface-anchored wall column (ring 0, which
// carries the surface distortion) AND the interior (ring >= kWallRings), by a
// same-side 3D inverse-distance blend. ring = min bin-distance to the Y/Z
// perimeter (0 = wall column). Ring 0 (surface) and the interior are kept as-is,
// so the wall constraint is preserved, not overwritten. This touches ONLY the
// Y/Z walls: the cathode (X=0) and anodes are X boundaries at large ring, left
// exactly as computed.
void smoothWallRings(TH3F* E) {
  if (kWallRings <= 1) return;
  const int nbinX = E->GetNbinsX(), nbinY = E->GetNbinsY(), nbinZ = E->GetNbinsZ();
  TAxis* xax = E->GetXaxis(); TAxis* yax = E->GetYaxis(); TAxis* zax = E->GetZaxis();
  auto ringOf = [&](int iy, int iz){
    return std::min(std::min(iy-1, nbinY-iy), std::min(iz-1, nbinZ-iz));
  };
  // anchors: ring 0 (surface) + interior (ring >= kWallRings)
  std::vector<double> PX, PY, PZ, PV; std::vector<int> PS;
  for (int ix = 1; ix <= nbinX; ++ix) {
    const double xc = xax->GetBinCenter(ix); const int side = (xc > 0);
    for (int iy = 1; iy <= nbinY; ++iy)
      for (int iz = 1; iz <= nbinZ; ++iz) {
        const int r = ringOf(iy, iz);
        if (r != 0 && r < kWallRings) continue;       // skip the intermediate rings
        const double v = E->GetBinContent(ix, iy, iz);
        if (v > -1e6) { PX.push_back(xc); PY.push_back(yax->GetBinCenter(iy));
                        PZ.push_back(zax->GetBinCenter(iz)); PV.push_back(v); PS.push_back(side); }
      }
  }
  if (PV.empty()) return;
  for (int ix = 1; ix <= nbinX; ++ix) {
    const double xc = xax->GetBinCenter(ix); const int side = (xc > 0);
    for (int iy = 1; iy <= nbinY; ++iy)
      for (int iz = 1; iz <= nbinZ; ++iz) {
        const int r = ringOf(iy, iz);
        if (r == 0 || r >= kWallRings) continue;      // keep surface (ring 0) and interior
        const double yc = yax->GetBinCenter(iy), zc = zax->GetBinCenter(iz);
        double sw = 0.0, sv = 0.0;
        for (size_t k = 0; k < PV.size(); ++k) {
          if (PS[k] != side) continue;                // stay on this drift side
          const double dx = PX[k]-xc, dy = PY[k]-yc, dz = PZ[k]-zc;
          const double d2 = dx*dx + dy*dy + dz*dz;
          if (d2 < 1e-9) { sw = 1.0; sv = PV[k]; break; }
          const double w = 1.0/(d2*d2);               // ~1/dist^4: nearest dominates
          sw += w; sv += w*PV[k];
        }
        if (sw > 0.0) E->SetBinContent(ix, iy, iz, sv/sw);
      }
  }
}

// Despike: replace a voxel by the neighbour median only if it is a local
// magnitude outlier -- |v0| exceeds every 3x3x3 neighbour (same side of the
// cathode) by more than kPeakFactor. Same as data_crossing_pdhd's medianSmooth.
void medianSmooth(TH3F* h) {
  const int nbinX = h->GetNbinsX(), nbinY = h->GetNbinsY(), nbinZ = h->GetNbinsZ();
  TH3F* src = (TH3F*)h->Clone(Form("%s_presmooth", h->GetName()));
  src->SetDirectory(nullptr);
  TAxis* xax = h->GetXaxis();
  std::vector<double> nb;
  for (int ix = 1; ix <= nbinX; ++ix) {
    const bool xpos = (xax->GetBinCenter(ix) > 0);
    for (int iy = 1; iy <= nbinY; ++iy)
      for (int iz = 1; iz <= nbinZ; ++iz) {
        const double v0 = src->GetBinContent(ix, iy, iz);
        nb.clear();
        for (int dx = -1; dx <= 1; ++dx)
          for (int dy = -1; dy <= 1; ++dy)
            for (int dz = -1; dz <= 1; ++dz) {
              const int jx = ix+dx, jy = iy+dy, jz = iz+dz;
              if (jx < 1 || jx > nbinX || jy < 1 || jy > nbinY || jz < 1 || jz > nbinZ) continue;
              if ((xax->GetBinCenter(jx) > 0) != xpos) continue;   // don't smooth across the cathode
              const double v = src->GetBinContent(jx, jy, jz);
              if (v > -1e6) nb.push_back(v);
            }
        if (nb.size() < 5) continue;
        bool outlier = true;
        for (double v : nb)
          if (std::fabs(v0) <= kPeakFactor * std::fabs(v)) { outlier = false; break; }
        if (outlier)
          h->SetBinContent(ix, iy, iz, TMath::Median(nb.size(), &nb[0]));
      }
  }
  delete src;
}

// Gentle low-pass: replace each voxel by the MEAN of its filled 3x3x3 neighbours
// on the same side of the cathode, applied `passes` times. The SCE displacement
// is physically smooth, so this removes the bin-scale statistical scatter (which
// the field derivative amplifies into speckle) while preserving the large-scale
// shape and the mean (normalization). passes=0 leaves the map untouched.
void boxSmooth(TH3F* h, int passes) {
  if (passes <= 0) return;
  const int nbinX = h->GetNbinsX(), nbinY = h->GetNbinsY(), nbinZ = h->GetNbinsZ();
  TAxis* xax = h->GetXaxis();
  for (int p = 0; p < passes; ++p) {
    TH3F* src = (TH3F*)h->Clone(Form("%s_prebox", h->GetName()));
    src->SetDirectory(nullptr);
    for (int ix = 1; ix <= nbinX; ++ix) {
      const bool xpos = (xax->GetBinCenter(ix) > 0);
      for (int iy = 1; iy <= nbinY; ++iy)
        for (int iz = 1; iz <= nbinZ; ++iz) {
          if (src->GetBinContent(ix, iy, iz) <= -1e6) continue;
          double sum = 0.0; int n = 0;
          for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
              for (int dz = -1; dz <= 1; ++dz) {
                const int jx = ix+dx, jy = iy+dy, jz = iz+dz;
                if (jx < 1 || jx > nbinX || jy < 1 || jy > nbinY || jz < 1 || jz > nbinZ) continue;
                if ((xax->GetBinCenter(jx) > 0) != xpos) continue;   // don't smooth across the cathode
                const double v = src->GetBinContent(jx, jy, jz);
                if (v > -1e6) { sum += v; ++n; }
              }
          if (n > 0) h->SetBinContent(ix, iy, iz, sum / n);
        }
    }
    delete src;
  }
}

// Clone a map onto the SCE MC input map's name and title (binning stays ours).
TH3F* makeReco(const char* name, TH3F* my) {
  TH3F* out = (TH3F*)my->Clone(name);
  out->SetDirectory(nullptr);
  TH3F* sce = g_sce ? (TH3F*)g_sce->Get(name) : nullptr;
  out->SetTitle(sce ? sce->GetTitle() : name);
  return out;
}

}  // namespace

void smooth_efield_pdhd(const char* infile  = "data_crossing_pdhd.root",
                        const char* outfile = "efield_pdhd.root") {
  TFile* f = TFile::Open(infile);
  if (f == nullptr || f->IsZombie()) { std::cout << "could not open " << infile << std::endl; return; }

  // the filled, un-smoothed backward displacement maps from data_crossing_pdhd
  TH3F* dX = dynamic_cast<TH3F*>(f->Get("XYZ_Xdist"));
  TH3F* dY = dynamic_cast<TH3F*>(f->Get("XYZ_Ydist"));
  TH3F* dZ = dynamic_cast<TH3F*>(f->Get("XYZ_Zdist"));
  if (dX == nullptr || dY == nullptr || dZ == nullptr) {
    std::cout << "missing XYZ_Xdist/Ydist/Zdist in " << infile << std::endl;
    f->Close();
    return;
  }
  dX->SetDirectory(nullptr); dY->SetDirectory(nullptr); dZ->SetDirectory(nullptr);

  // grab the track-validation histograms + per-voxel diagnostics before closing
  // the input, detached, so they can be written into the output unchanged.
  std::vector<TH1*> validation;
  for (const char* n : kValidationNames) {
    TH1* h = dynamic_cast<TH1*>(f->Get(n));
    if (h == nullptr) continue;
    TH1* c = static_cast<TH1*>(h->Clone(n));
    c->SetDirectory(nullptr);
    validation.push_back(c);
  }
  f->Close();

  const int nbinX = dX->GetNbinsX(), nbinY = dX->GetNbinsY(), nbinZ = dX->GetNbinsZ();
  const double wX = dX->GetXaxis()->GetBinWidth(1);

  // --- despike outliers, then gentle low-pass (both iterate without the 4 h
  //     crosser pass). kSmoothPasses=0 recovers the un-low-passed field.
  medianSmooth(dX);
  medianSmooth(dY);
  medianSmooth(dZ);
  boxSmooth(dX, kSmoothPasses);
  boxSmooth(dY, kSmoothPasses);
  boxSmooth(dZ, kSmoothPasses);

  // --- recompute the E field from the smoothed displacement ------------------
  TH3F* Ex = (TH3F*)dX->Clone("Efield_X"); Ex->SetDirectory(nullptr); Ex->Reset();
  TH3F* Ey = (TH3F*)dY->Clone("Efield_Y"); Ey->SetDirectory(nullptr); Ey->Reset();
  TH3F* Ez = (TH3F*)dZ->Clone("Efield_Z"); Ez->SetDirectory(nullptr); Ez->Reset();
  EfieldFromDistortions(dX, dY, dZ, Ex, Ey, Ez, wX, nbinX, nbinY, nbinZ);

  // --- third-ring fix: extend the interior field over the outer Y/Z wall rings
  smoothWallRings(Ex);
  smoothWallRings(Ey);
  smoothWallRings(Ez);

  // --- write smoothed backward displacement + recomputed E field -------------
  TH3F* Bkwd_X_Pos = makeReco("RecoBkwd_Displacement_X_Pos", dX);
  TH3F* Bkwd_X_Neg = makeReco("RecoBkwd_Displacement_X_Neg", dX);
  TH3F* Bkwd_Y_Pos = makeReco("RecoBkwd_Displacement_Y_Pos", dY);
  TH3F* Bkwd_Y_Neg = makeReco("RecoBkwd_Displacement_Y_Neg", dY);
  TH3F* Bkwd_Z_Pos = makeReco("RecoBkwd_Displacement_Z_Pos", dZ);
  TH3F* Bkwd_Z_Neg = makeReco("RecoBkwd_Displacement_Z_Neg", dZ);
  TH3F* Reco_Ex_Pos = makeReco("Reco_ElecField_X_Pos", Ex);
  TH3F* Reco_Ex_Neg = makeReco("Reco_ElecField_X_Neg", Ex);
  TH3F* Reco_Ey_Pos = makeReco("Reco_ElecField_Y_Pos", Ey);
  TH3F* Reco_Ey_Neg = makeReco("Reco_ElecField_Y_Neg", Ey);
  TH3F* Reco_Ez_Pos = makeReco("Reco_ElecField_Z_Pos", Ez);
  TH3F* Reco_Ez_Neg = makeReco("Reco_ElecField_Z_Neg", Ez);

  TFile fout(outfile, "RECREATE");
  dX->Write("XYZ_Xdist"); dY->Write("XYZ_Ydist"); dZ->Write("XYZ_Zdist");
  Bkwd_X_Pos->Write(); Bkwd_X_Neg->Write();
  Bkwd_Y_Pos->Write(); Bkwd_Y_Neg->Write();
  Bkwd_Z_Pos->Write(); Bkwd_Z_Neg->Write();
  Reco_Ex_Pos->Write(); Reco_Ex_Neg->Write();
  Reco_Ey_Pos->Write(); Reco_Ey_Neg->Write();
  Reco_Ez_Pos->Write(); Reco_Ez_Neg->Write();
  for (TH1* h : validation) h->Write();   // track-validation histos + diagnostics passthrough
  fout.Close();
  std::cout << "smooth_efield_pdhd wrote " << outfile << std::endl;
}
