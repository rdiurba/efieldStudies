#include <TCanvas.h>
#include <TFile.h>
#include <TH3.h>
#include <TH2.h>
#include <TH1.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TString.h>
#include <TROOT.h>
#include <TColor.h>
#include <TLatex.h>
#include <TLine.h>

#include <algorithm>
#include <vector>
#include <iostream>

// CET-CBTD1 (colorcet "diverging_tritanopic_cwr_75_98_c20"): a cyan-white-red
// diverging map built to stay legible under tritanopia (blue-yellow CVD) -- it
// uses cyan and red, the two ends a tritanope can still separate, with white at
// the center. Control points sampled along the colorcet map (ROOT can't load
// colorcet directly), spanning cyan -> white -> salmon-red.
//   For signed maps set a SYMMETRIC z-range so white lands on 0.
static void setCBTD1Palette(int nCont = 255) {
  const int n = 9;
  double stops[n] = {0.000, 0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1.000};
  double red[n]   = {0.000, 0.400, 0.620, 0.800, 0.960, 0.980, 0.990, 0.980, 0.960};
  double green[n] = {0.770, 0.830, 0.880, 0.920, 0.960, 0.850, 0.720, 0.600, 0.490};
  double blue[n]  = {0.880, 0.900, 0.930, 0.940, 0.950, 0.810, 0.660, 0.540, 0.430};
  TColor::CreateGradientColorTable(n, stops, red, green, blue, nCont);
  gStyle->SetNumberContours(nCont);
}

// load the shared ProtoDUNE style so title sizes, fonts, palette and margins all
// come from the style guide.  The style puts the title at Y=0.98; shift it in X so
// it clears the DUNE watermark drawn at (0.20, 0.94) on the line below it.
// The CET-CBTD1 palette is installed AFTER ForceStyle so it overrides the style's
// palette -- every COLZ TH2 drawn through these plotters then uses it.
static void useProtoDUNEStyle() {
  static bool loaded = false;
  if (!loaded) { gROOT->LoadMacro("protoDUNEStyle.C"); loaded = true; }
  gROOT->SetStyle("protoDUNEStyle");
  gROOT->ForceStyle();
  gStyle->SetTitleX(0.35);
  setCBTD1Palette();
}

// the DUNE:ProtoDUNE-<det> watermark, drawn after the histogram at the standard
// (0.20, 0.94) NDC spot -- one line below the (shifted) ROOT title so both show.
static void drawWatermark(const TString& tag) {
  TLatex tL; tL.SetNDC(); tL.SetTextFont(42);
  tL.DrawLatex(0.20, 0.94, Form("#bf{DUNE:%s}", tag.Data()));
}

// clean, filesystem-friendly token for a slice position: -200 -> "m200", 0 -> "0".
static TString posTag(double p) {
  return p < 0 ? Form("m%.0f", -p) : Form("%.0f", p);
}

// write the current canvas as both .pdf and .png under one clean base name.
static void printBoth(TCanvas* c, const TString& base) {
  c->Print(base + ".pdf");
  c->Print(base + ".png");
}

// The agreed detector binning (must match data_crossing_pdhd.C): 20 cm bins over
// X[-360,360], Y[0,600], Z[0,460]. Every slice is drawn on THIS grid (and the
// source map is sampled by coordinate), so the displayed bins are always uniform
// 20 cm regardless of what binning the input histogram happens to carry.
static const double gXmin = -340, gXmax = 340;
static const double gYmin = -340, gYmax = 340;
static const double gZmin =    0, gZmax = 300;
static const double gBin  =   20;
static void axisGeom(int idx, int& n, double& lo, double& hi) {
  if      (idx == 0) { lo = gXmin; hi = gXmax; }
  else if (idx == 1) { lo = gYmin; hi = gYmax; }
  else               { lo = gZmin; hi = gZmax; }
  n = int((hi - lo) / gBin);
}

void plot_distortion_slices(const char* infile = "data_crossing_pdvd.root",
                            const char* distAxis = "Z",
                            const char* sliceAxis = "X",
                            std::vector<double> planes = {-300, -150, 0, 150, 300},
                            bool truth = false,
                            const char* detector = "VD",
                            const char* histName = "") {
  TFile* f = TFile::Open(infile);
  if (f == nullptr || f->IsZombie()) {
    std::cout << "could not open " << infile << std::endl;
    return;
  }

  const TString dA(distAxis), sA(sliceAxis);
  // histName overrides the default XYZ_<dist>dist name, so any 3D map in the
  // file can be sliced, e.g. "Reco_ElecField_Z_Neg", "RecoBkwd_Displacement_Z_Pos".
  const TString hname = (TString(histName).Length() > 0)
                          ? TString(histName)
                          : (truth ? Form("XYZ_%sdist_truth", dA.Data())
                          : Form("Bkwd_%s_Pos", dA.Data()));
  TH3* h = dynamic_cast<TH3*>(f->Get(hname));
  if (h == nullptr) { std::cout << "missing histogram: " << hname << std::endl; return; }

  // axis index: 0=X, 1=Y, 2=Z
  auto axisIndex = [](const TString& a) -> int { return a == "X" ? 0 : (a == "Y" ? 1 : 2); };
  auto axisPtr   = [&](int idx) -> TAxis* {
    return idx == 0 ? h->GetXaxis() : (idx == 1 ? h->GetYaxis() : h->GetZaxis());
  };
  const char* axisNm[3] = {"X", "Y", "Z"};

  const int sIdx = axisIndex(sA);
  // horizontal / vertical axes for the 2D slice -- X is always drawn on the
  // horizontal axis when it is one of the two plotted coordinates.
  int hIdx, vIdx;
  if      (sIdx == 0) { hIdx = 2; vIdx = 1; }   // slice X -> Z(x) vs Y(y)
  else if (sIdx == 1) { hIdx = 0; vIdx = 2; }   // slice Y -> X(x) vs Z(y)
  else                { hIdx = 0; vIdx = 1; }   // slice Z -> X(x) vs Y(y)

  TAxis* sAxis = axisPtr(sIdx);
  TAxis* hAxis = axisPtr(hIdx);
  TAxis* vAxis = axisPtr(vIdx);

  const TString tag = (TString(detector) == "SP") ? "ProtoDUNE-SP"
                    : (TString(detector) == "VD") ? "ProtoDUNE-VD"
                                                  : "ProtoDUNE-HD";

  useProtoDUNEStyle();

  static int ncall = 0; ++ncall;
  TCanvas* c = new TCanvas(Form("cSlice_%d", ncall), "distortion slices", 950, 720);
  c->SetRightMargin(0.16);   // room for the COLZ palette
  c->SetLeftMargin(0.13);

  int drawn = 0;
  for (double pv : planes) {
    const int sb = sAxis->FindBin(pv);
    if (sb < 1 || sb > sAxis->GetNbins()) {
      std::cout << "slice " << sA << " = " << pv << " out of range, skipping" << std::endl;
      continue;
    }

    int hN, vN; double hlo, hhi, vlo, vhi;
    axisGeom(hIdx, hN, hlo, hhi);
    axisGeom(vIdx, vN, vlo, vhi);
    TH2D* hs = new TH2D(Form("slice_%s_%g_%d", sA.Data(), pv, ncall), "",
                        hN, hlo, hhi, vN, vlo, vhi);
    for (int ih = 1; ih <= hN; ++ih) {
      for (int iv = 1; iv <= vN; ++iv) {
        double p[3];
        p[sIdx] = pv;
        p[hIdx] = hlo + (ih - 0.5) * (hhi - hlo) / hN;
        p[vIdx] = vlo + (iv - 0.5) * (vhi - vlo) / vN;
        const double val = h->GetBinContent(h->FindBin(p[0], p[1], p[2]));
        if (val < -1e6) continue;           // sentinel (empty bin) -> leave blank
        hs->SetBinContent(ih, iv, val);
      }
    }

    const bool named = (TString(histName).Length() > 0);
    hs->SetTitle(Form("%s  |  %s = %.0f cm%s;%s [cm];%s [cm]",
                      named ? hname.Data() : Form("#Delta%s distortion", dA.Data()),
                      sA.Data(), sAxis->GetBinCenter(sb), truth ? " (SCE truth)" : "",
                      axisNm[hIdx], axisNm[vIdx]));
    hs->GetXaxis()->CenterTitle();
    hs->GetYaxis()->CenterTitle();
    hs->GetZaxis()->SetTitle(named ? "value" : Form("#Delta%s [cm]", dA.Data()));
    hs->GetZaxis()->CenterTitle();
    hs->GetZaxis()->SetRangeUser(-30,30);
    c->cd();
    hs->Draw("COLZ");
    drawWatermark(tag);
    // one clean-named file per slice (pdf + png), no overwrite
    const TString base = named
             ? Form("slice_%s_%s%s_%s", hname.Data(), sA.Data(), posTag(pv).Data(),detector)
             : Form("distortion_%sdist%s_%s%s_%s", dA.Data(), truth ? "_truth" : "",
                    sA.Data(), posTag(pv).Data(),detector);
    printBoth(c, base);
    ++drawn;
    delete hs;
  }

  std::cout << "wrote " << drawn << " slices (pdf+png) for " << infile << std::endl;
  f->Close();
}

// -----------------------------------------------------------------------------
// Generic COLZ slicer for a single named 3D map (STD, occupancy, ...), drawn the
// same way as the distortion slices but with an AUTO z-range (these maps are not
// signed distortions, so the fixed [-30,30] window does not apply) and its own
// title / z-axis label / output base name.
// -----------------------------------------------------------------------------
static void plot_map_slices(const char* infile, const TString& hname,
                            const TString& titleBase, const TString& zTitle,
                            const TString& fileBase, const char* sliceAxis,
                            std::vector<double> planes, const char* detector) {
  TFile* f = TFile::Open(infile);
  if (f == nullptr || f->IsZombie()) { std::cout << "could not open " << infile << std::endl; return; }

  TH3* h = dynamic_cast<TH3*>(f->Get(hname));
  if (h == nullptr) { std::cout << "missing histogram: " << hname << std::endl; f->Close(); return; }

  auto axisIndex = [](const TString& a) -> int { return a == "X" ? 0 : (a == "Y" ? 1 : 2); };
  auto axisPtr   = [&](int idx) -> TAxis* {
    return idx == 0 ? h->GetXaxis() : (idx == 1 ? h->GetYaxis() : h->GetZaxis());
  };
  const char* axisNm[3] = {"X", "Y", "Z"};

  const TString sA(sliceAxis);
  const int sIdx = axisIndex(sA);
  int hIdx, vIdx;
  if      (sIdx == 0) { hIdx = 2; vIdx = 1; }   // slice X -> Z(x) vs Y(y)
  else if (sIdx == 1) { hIdx = 0; vIdx = 2; }   // slice Y -> X(x) vs Z(y)
  else                { hIdx = 0; vIdx = 1; }   // slice Z -> X(x) vs Y(y)

  TAxis* sAxis = axisPtr(sIdx);

  const TString tag = (TString(detector) == "SP") ? "ProtoDUNE-SP"
                    : (TString(detector) == "VD") ? "ProtoDUNE-VD"
                                                  : "ProtoDUNE-HD";
  useProtoDUNEStyle();

  static int ncall = 0; ++ncall;
  TCanvas* c = new TCanvas(Form("cDiag_%d", ncall), "diagnostic slices", 950, 720);
  c->SetRightMargin(0.16);
  c->SetLeftMargin(0.13);

  int drawn = 0;
  for (double pv : planes) {
    const int sb = sAxis->FindBin(pv);
    if (sb < 1 || sb > sAxis->GetNbins()) {
      std::cout << "slice " << sA << " = " << pv << " out of range, skipping" << std::endl;
      continue;
    }

    int hN, vN; double hlo, hhi, vlo, vhi;
    axisGeom(hIdx, hN, hlo, hhi);
    axisGeom(vIdx, vN, vlo, vhi);
    TH2D* hs = new TH2D(Form("diagslice_%s_%g_%d", sA.Data(), pv, ncall), "",
                        hN, hlo, hhi, vN, vlo, vhi);
    for (int ih = 1; ih <= hN; ++ih) {
      for (int iv = 1; iv <= vN; ++iv) {
        double p[3];
        p[sIdx] = pv;
        p[hIdx] = hlo + (ih - 0.5) * (hhi - hlo) / hN;
        p[vIdx] = vlo + (iv - 0.5) * (vhi - vlo) / vN;
        const double val = h->GetBinContent(h->FindBin(p[0], p[1], p[2]));
        if (val < -1e6) continue;           // sentinel (empty bin) -> leave blank
        hs->SetBinContent(ih, iv, val);
      }
    }

    hs->SetTitle(Form("%s  |  %s = %.0f cm;%s [cm];%s [cm]",
                      titleBase.Data(), sA.Data(), sAxis->GetBinCenter(sb),
                      axisNm[hIdx], axisNm[vIdx]));
    hs->GetXaxis()->CenterTitle();
    hs->GetYaxis()->CenterTitle();
    hs->GetZaxis()->SetTitle(zTitle);
    hs->GetZaxis()->CenterTitle();
    c->cd();
    hs->Draw("COLZ");
    drawWatermark(tag);
    printBoth(c, Form("%s_%s%s_%s", fileBase.Data(), sA.Data(), posTag(pv).Data(), detector));
    ++drawn;
    delete hs;
  }

  std::cout << "wrote " << drawn << " slices (pdf+png) for " << hname << std::endl;
  f->Close();
}

// -----------------------------------------------------------------------------
// Companion: median velocity and E-field vs X from velocity_pdvd.root.
//   velocity.png : median drift velocity vs X (med_vel)
//   efield.png   : median E field vs X, Data (med_Ef) vs MC-estimated (input_Ef)
//
//   root -l -b -q 'plot_distortions_pdhd.C+' -e 'plot_velocity("velocity_pdvd.root","HD")'
// -----------------------------------------------------------------------------
void plot_velocity(const char* velfile = "velocity_pdvd.root", const char* detector = "VD") {
  TFile* f = TFile::Open(velfile);
  if (f == nullptr || f->IsZombie()) {
    std::cout << "could not open " << velfile << std::endl;
    return;
  }

  TH1* medVel = dynamic_cast<TH1*>(f->Get("med_vel"));     // median v(X) [mm/us]
  TH1* medEf  = dynamic_cast<TH1*>(f->Get("med_Ef"));      // Data E(X) [kV/cm]
  TH1* inEf   = dynamic_cast<TH1*>(f->Get("input_Ef"));    // MC-estimated E(X) [kV/cm]
  if (medVel == nullptr || medEf == nullptr || inEf == nullptr) {
    std::cout << "missing med_vel / med_Ef / input_Ef in " << velfile << std::endl;
    f->Close();
    return;
  }

  const TString tag = (TString(detector) == "SP") ? "ProtoDUNE-SP"
                    : (TString(detector) == "VD") ? "ProtoDUNE-VD"
                                                  : "ProtoDUNE-HD";
  useProtoDUNEStyle();

  // y-range over the filled bins only (empty bins are exactly zero)
  auto filledRange = [](std::initializer_list<TH1*> hs, double& lo, double& hi) {
    lo = 1e30; hi = -1e30;
    for (TH1* h : hs)
      for (int b = 1; b <= h->GetNbinsX(); ++b) {
        const double y = h->GetBinContent(b);
        if (y == 0) continue;
        lo = std::min(lo, y); hi = std::max(hi, y);
      }
  };

  static int ncall = 0; ++ncall;
  TCanvas* c = new TCanvas(Form("cVel_%d", ncall), "velocity / E field", 950, 720);
  c->SetLeftMargin(0.13);
  c->SetRightMargin(0.025);

  // --- velocity.png: median drift velocity vs X -----------------------------
  double wlo, whi; filledRange({medVel}, wlo, whi);
  const double wpad = 0.1 * (whi - wlo);
  medVel->SetTitle("Median Drift Velocity;X [cm];v [mm/#mus]");
  medVel->SetLineColor(kRed);
  medVel->GetXaxis()->CenterTitle();
  medVel->GetYaxis()->CenterTitle();
  medVel->GetYaxis()->SetRangeUser(1.2, 1.7);  // clips the empty (0) bins
  medVel->Draw("HIST");
  drawWatermark(tag);
  c->Print("velocity.png");

  // --- efield.png: Data vs MC-estimated E field vs X ------------------------
  double elo, ehi; filledRange({medEf, inEf}, elo, ehi);
  const double epad = 0.1 * (ehi - elo);
  medEf->SetTitle("Median E field;X [cm];E [kV/cm]");
  medEf->SetLineColor(kRed);      // Data
  medEf->GetXaxis()->CenterTitle();
  medEf->GetYaxis()->SetRangeUser(0.3,0.7);
  medEf->GetYaxis()->CenterTitle();
  //medEf->GetYaxis()->SetRangeUser(elo - epad, ehi + epad);
  inEf->SetLineColor(kBlue);     // MC estimated
  inEf->SetLineStyle(2);              // dashed to separate from Data

  medEf->Draw("HIST");
  inEf->Draw("HIST SAME");

  TLegend* leg = new TLegend(0.55, 0.75, 0.94, 0.90);
  leg->AddEntry(medEf, "Data",         "l");
  leg->AddEntry(inEf,  "MC Estimated", "l");
  leg->Draw();
  drawWatermark(tag);
  c->Print("efield.png");

  std::cout << "wrote velocity.png and efield.png" << std::endl;
  f->Close();
}

// -----------------------------------------------------------------------------
// 2D slices of the reconstructed E-field maps (Reco_ElecField_<comp>_{Pos,Neg}),
// analogous to plot_distortion_slices but combining the two drift sides by sign
// of X. The maps store the relative field distortion (SCE convention).
//   comp      : "X" | "Y" | "Z"   (field component)
//   sliceAxis : "X" | "Y" | "Z"   (coordinate held fixed)
//   planes    : positions [cm] along sliceAxis to slice at
// X is always drawn on the horizontal axis when it is one of the two plotted axes.
//
//   root -l -b -q 'plot_distortions_pdhd.C+' -e 'plot_efield_slices("data_crossing_pdvd.root","Z","X",{-200,0,200})'
// -----------------------------------------------------------------------------
void plot_efield_slices(const char* infile = "efield_pdvd.root",
                        const char* comp = "Z",
                        const char* sliceAxis = "X",
                        std::vector<double> planes = {-300, -150, 0, 150, 300},
                        const char* detector = "VD") {
  TFile* f = TFile::Open(infile);
  if (f == nullptr || f->IsZombie()) {
    std::cout << "could not open " << infile << std::endl;
    return;
  }

  TH3* hPos = dynamic_cast<TH3*>(f->Get(Form("Reco_ElecField_%s_Pos", comp)));
  TH3* hNeg = dynamic_cast<TH3*>(f->Get(Form("Reco_ElecField_%s_Neg", comp)));
  if (hPos == nullptr || hNeg == nullptr) {
    std::cout << "missing Reco_ElecField_" << comp << "_{Pos,Neg}" << std::endl;
    return;
  }

  // draw on the agreed 20 cm grid (same as the displacement slices); efval picks
  // the side by sign(x) and samples the map by coordinate.
  struct Ax { int n; double lo, hi; };
  Ax ax[3];
  axisGeom(0, ax[0].n, ax[0].lo, ax[0].hi);
  axisGeom(1, ax[1].n, ax[1].lo, ax[1].hi);
  axisGeom(2, ax[2].n, ax[2].lo, ax[2].hi);
  auto center = [&](int a, int i) { return ax[a].lo + (i - 0.5) * (ax[a].hi - ax[a].lo) / ax[a].n; };
  const char* axisNm[3] = {"X", "Y", "Z"};
  auto axisIndex = [](const TString& a) -> int { return a == "X" ? 0 : (a == "Y" ? 1 : 2); };

  const int sIdx = axisIndex(sliceAxis);
  int hIdx, vIdx;
  if      (sIdx == 0) { hIdx = 2; vIdx = 1; }   // slice X -> Z(x) vs Y(y)
  else if (sIdx == 1) { hIdx = 0; vIdx = 2; }   // slice Y -> X(x) vs Z(y)
  else                { hIdx = 0; vIdx = 1; }   // slice Z -> X(x) vs Y(y)

  // field value at a 3D point, from the side-appropriate map
  auto efval = [&](double x, double y, double z) -> double {
    TH3* m = (x >= 0) ? hPos : hNeg;
    return m->GetBinContent(m->FindBin(x, y, z));
  };

  const TString tag = (TString(detector) == "SP") ? "ProtoDUNE-SP"
                    : (TString(detector) == "VD") ? "ProtoDUNE-VD"
                                                  : "ProtoDUNE-HD";
  useProtoDUNEStyle();

  // maps hold the relative field distortion (SCE convention), plotted as stored
  // (dimensionless): (E_X - E_0)/E_0 for the drift component, E_Y/E_0 and
  // E_Z/E_0 for the transverse ones.
  const bool isX = (TString(comp) == "X");

  static int ncall = 0; ++ncall;
  TCanvas* c = new TCanvas(Form("cEfield_%d", ncall), "E field slices", 950, 720);
  c->SetRightMargin(0.16);   // room for the COLZ palette
  c->SetLeftMargin(0.13);

  int drawn = 0;
  for (double pv : planes) {
    TH2D* hs = new TH2D(Form("efslice_%s_%g_%d", sliceAxis, pv, ncall), "",
                        ax[hIdx].n, ax[hIdx].lo, ax[hIdx].hi,
                        ax[vIdx].n, ax[vIdx].lo, ax[vIdx].hi);
    for (int ih = 1; ih <= ax[hIdx].n; ++ih)
      for (int iv = 1; iv <= ax[vIdx].n; ++iv) {
        double p[3];
        p[sIdx] = pv; p[hIdx] = center(hIdx, ih); p[vIdx] = center(vIdx, iv);
        const double val = efval(p[0], p[1], p[2]);
        if (val < -1e6) continue;                 // sentinel (empty) -> blank
        hs->SetBinContent(ih, iv, val);           // relative distortion, as stored
      }
    hs->SetTitle(Form("Reco E_{%s}  |  %s = %.0f cm;%s [cm];%s [cm]",
                      comp, sliceAxis, pv, axisNm[hIdx], axisNm[vIdx]));
    hs->GetZaxis()->SetRangeUser(-0.5,0.5);
    hs->GetXaxis()->CenterTitle();
    hs->GetYaxis()->CenterTitle();
    hs->GetZaxis()->SetTitle(isX ? "(E_{X} #minus E_{0})/E_{0}"
                                 : Form("E_{%s}/E_{0}", comp));
    hs->GetZaxis()->CenterTitle();
    c->cd();
    hs->Draw("COLZ");
    drawWatermark(tag);
    printBoth(c, Form("efield_%s_%s%s_%s", comp, sliceAxis, posTag(pv).Data(),detector));
    ++drawn;
    delete hs;
  }
  std::cout << "wrote " << drawn << " E-field slices (pdf+png) for " << infile << std::endl;
  f->Close();
}
void plot_trackcuts(const char* infile = "data_crossing_pdvd.root", const char* detector = "VD") {
  TFile* f = TFile::Open(infile);
  if (f == nullptr || f->IsZombie()) { std::cout << "could not open " << infile << std::endl; return; }
  const TString tag = (TString(detector) == "SP") ? "ProtoDUNE-SP"
                    : (TString(detector) == "VD") ? "ProtoDUNE-VD" : "ProtoDUNE-HD";
  useProtoDUNEStyle();
  const char* names[] = {"trk_thetaxz", "trk_thetayz", "trk_theta", "trk_phi", "trk_t0proxy", "pair_dca"};
  static int nc = 0; int drawn = 0;
  for (const char* nm : names) {
    TH1* h = dynamic_cast<TH1*>(f->Get(nm));
    if (h == nullptr) { std::cout << "missing " << nm << std::endl; continue; }
    TCanvas* c = new TCanvas(Form("cTC_%d", ++nc), "track cuts", 900, 700);
    c->SetLeftMargin(0.13); c->SetRightMargin(0.05);
    h->SetLineColor(kAzure + 2);
    h->GetXaxis()->CenterTitle(); h->GetYaxis()->CenterTitle();
    h->Draw("HIST");
    drawWatermark(tag);
    printBoth(c, nm);
    ++drawn; delete c;
  }
  std::cout << "wrote " << drawn << " track-selection TH1s" << std::endl;
  f->Close();
}

// Occupancy TH2 slices (XYZ_Occupancy), drawn like the distortion slices.
void plot_occupancy(const char* infile = "data_crossing_pdvd.root",
                    const char* sliceAxis = "X",
                    std::vector<double> planes = {-300, -150, 0, 150, 300},
                    const char* detector = "VD") {
  plot_map_slices(infile, "XYZ_Occupancy", "Voxel occupancy", "entries / voxel",
                  "occupancy", sliceAxis, planes, detector);
}

// Per-voxel STD TH2 slices: comp = "X"|"Y"|"Z" -> XYZ_<comp>std, drawn like the
// distortion slices.
void plot_std_slices(const char* infile = "data_crossing_pdvd.root",
                     const char* comp = "X",
                     const char* sliceAxis = "X",
                     std::vector<double> planes = {-300, -150, 0, 150, 300},
                     const char* detector = "VD") {
  const TString cA(comp);
  plot_map_slices(infile, Form("XYZ_%sstd", cA.Data()),
                  Form("#Delta%s STD", cA.Data()),
                  Form("#sigma(#Delta%s) [cm]", cA.Data()),
                  Form("std_%sstd", cA.Data()), sliceAxis, planes, detector);
}

void plot_std(const char* infile = "data_crossing_pdvd.root", const char* detector = "VD") {
  TFile* f = TFile::Open(infile);
  if (f == nullptr || f->IsZombie()) { std::cout << "could not open " << infile << std::endl; return; }
  TH3* sx = dynamic_cast<TH3*>(f->Get("XYZ_Xstd"));
  TH3* sy = dynamic_cast<TH3*>(f->Get("XYZ_Ystd"));
  TH3* sz = dynamic_cast<TH3*>(f->Get("XYZ_Zstd"));
  if (sx == nullptr || sy == nullptr || sz == nullptr) {
    std::cout << "missing XYZ_{X,Y,Z}std in " << infile << std::endl; f->Close(); return;
  }
  const TString tag = (TString(detector) == "SP") ? "ProtoDUNE-SP"
                    : (TString(detector) == "VD") ? "ProtoDUNE-VD" : "ProtoDUNE-HD";
  useProtoDUNEStyle();
  TH1F* hx = new TH1F("std_all_X", "#Delta STD;STD [cm];voxels", 100, 0, 25);
  TH1F* hy = (TH1F*)hx->Clone("std_all_Y");
  TH1F* hz = (TH1F*)hx->Clone("std_all_Z");
  auto fillStd = [](TH3* s, TH1* h) {
    for (int ix = 1; ix <= s->GetNbinsX(); ++ix)
      for (int iy = 1; iy <= s->GetNbinsY(); ++iy)
        for (int iz = 1; iz <= s->GetNbinsZ(); ++iz) {
          const double v = s->GetBinContent(ix, iy, iz);
          if (v > 0) h->Fill(v);
        }
  };
  fillStd(sx, hx); fillStd(sy, hy); fillStd(sz, hz);
  hx->SetLineColor(kRed + 1); hy->SetLineColor(kAzure + 2); hz->SetLineColor(kGreen + 2);
  hx->GetXaxis()->CenterTitle(); hx->GetYaxis()->CenterTitle();
  const double ymax = std::max({hx->GetMaximum(), hy->GetMaximum(), hz->GetMaximum()});
  hx->SetMaximum(1.15 * ymax);
  static int nc = 0;
  TCanvas* c = new TCanvas(Form("cStd_%d", ++nc), "STD across bins", 900, 700);
  c->SetLeftMargin(0.13); c->SetRightMargin(0.05);
  hx->Draw("HIST"); hy->Draw("HIST SAME"); hz->Draw("HIST SAME");
  TLegend* leg = new TLegend(0.62, 0.72, 0.94, 0.90);
  leg->AddEntry(hx, "#DeltaX STD", "l");
  leg->AddEntry(hy, "#DeltaY STD", "l");
  leg->AddEntry(hz, "#DeltaZ STD", "l");
  leg->Draw();
  drawWatermark(tag);
  printBoth(c, "std_overall");
  std::cout << "wrote std_overall (pdf+png)" << std::endl;
  f->Close();
}

// -----------------------------------------------------------------------------
// Driver: every distortion + E-field slice (all components x all slice axes)
// plus the velocity/E-field summary, using per-slice-axis plane presets.
//   root -l -b -q 'plot_distortions_pdhd.C+' -e 'plot_distortions()'
// -----------------------------------------------------------------------------
// defined below (after plot_distortions); declared here so this driver can call them
void plot_velocity_split(const char* velfile = "velocity_pdvd.root",
                         const char* detector = "VD", double eNominal = 0.488);
void plot_med_vel_tpc(const char* velfile = "velocity_pdvd.root",
                      const char* detector = "VD");
void plot_surface_distortions(const char* infile = "surface_distortions_pdvd.root",
                              const char* detector = "VD");

void plot_distortions(const char* infile   = "efield_pdvd.root",
                      const char* velfile  = "velocity_pdvd.root",
                      const char* detector = "VD") {
  const char* c[3] = {"X", "Y", "Z"};
  std::vector<double> planes[3] = { {-300,-150,0,150,300},   // slice X
                                    {-300,-150,0,150,300},   // slice Y
                                    { 50,150,250} };         // slice Z
  for (int i = 0; i < 3; ++i)        // distortion / field component
    for (int s = 0; s < 3; ++s) {    // axis held fixed
      plot_distortion_slices(infile, c[i], c[s], planes[s], false, detector);
      plot_efield_slices(infile, c[i], c[s], planes[s], detector);
      plot_std_slices(infile, c[i], c[s], planes[s], detector);   // per-voxel STD TH2s
    }
  for (int s = 0; s < 3; ++s)        // occupancy TH2s (single map, per slice axis)
    plot_occupancy(infile, c[s], planes[s], detector);
  plot_velocity(velfile, detector);
  plot_velocity_split(velfile, detector);   // PDHD anode-to-anode / per-side / per-TPC
  plot_med_vel_tpc(velfile, detector);      // PDHD med_vel(X) split by TPC
  plot_surface_distortions("surface_distortions_pdvd.root", detector);  // detector-face maps
  plot_trackcuts(infile, detector);
    plot_std(infile, detector);

}
// -----------------------------------------------------------------------------
// PDHD-only: the three-level anode-to-cathode drift-velocity / E-field
// measurements made in velocity_pdhd.C --
//   (1) full anode-to-anode        vel_anode2anode / Efield_anode2anode
//   (2) per cathode side           vel_neg,vel_pos / Efield_neg,Efield_pos
//   (3) per TPC (Z<230 / Z>230)    vel_tpc_*        / Efield_tpc_*
// Draws the velocity distributions overlaid (side + a2a), the four TPC
// distributions overlaid, and an E-field summary point per measurement against
// the nominal field.
//   root -l -b -q 'plot_distortions_pdhd.C+' -e 'plot_velocity_split("velocity_pdvd.root","HD",0.488)'
// -----------------------------------------------------------------------------
void plot_velocity_split(const char* velfile, const char* detector, double eNominal) {
  TFile* f = TFile::Open(velfile);
  if (f == nullptr || f->IsZombie()) { std::cout << "could not open " << velfile << std::endl; return; }

  const TString tag = (TString(detector)=="SP") ? "ProtoDUNE-SP"
                    : (TString(detector)=="VD") ? "ProtoDUNE-VD" : "ProtoDUNE-HD";
  useProtoDUNEStyle();
  auto get = [&](const char* n)->TH1*{ return dynamic_cast<TH1*>(f->Get(n)); };

  static int nc = 0; ++nc;

  // --- (1)+(2) anode-to-anode + each cathode side ----------------------------
  TH1* vA = get("vel_anode2anode");
  TH1* vN = get("vel_neg");
  TH1* vP = get("vel_pos");
  if (vA && vN && vP) {
    TCanvas* c1 = new TCanvas(Form("cVsplit_%d", nc), "velocity split", 950, 720);
    c1->SetLeftMargin(0.13); c1->SetRightMargin(0.025);
    vA->SetLineColor(kBlack); vN->SetLineColor(kBlue); vP->SetLineColor(kRed);
    vA->SetLineWidth(2); vN->SetLineWidth(2); vP->SetLineWidth(2);
    vA->SetTitle("Drift velocity: anode-to-anode vs each side;v [mm/#mus];tracks");
    vA->GetXaxis()->CenterTitle(); vA->GetYaxis()->CenterTitle();
    vA->GetXaxis()->SetRangeUser(1.3, 1.9);
    const double m = std::max({vA->GetMaximum(), vN->GetMaximum(), vP->GetMaximum()});
    vA->SetMaximum(1.25 * m);
    vA->Draw("HIST"); vN->Draw("HIST SAME"); vP->Draw("HIST SAME");
    TLegend* l = new TLegend(0.52, 0.70, 0.94, 0.90);
    l->AddEntry(vA, Form("anode-to-anode (#mu=%.3f)", vA->GetMean()), "l");
    l->AddEntry(vN, Form("#minusX side (#mu=%.3f)",    vN->GetMean()), "l");
    l->AddEntry(vP, Form("+X side (#mu=%.3f)",         vP->GetMean()), "l");
    l->Draw();
    drawWatermark(tag);
    c1->Print("velocity_split.png");
  }

  // --- (3) per-TPC velocity distributions ------------------------------------
  const char* tn[4]  = {"Xneg_Zlo", "Xneg_Zhi", "Xpos_Zlo", "Xpos_Zhi"};
  const int   col[4] = {kBlue, kAzure+7, kRed, kOrange+7};
  TH1* vt[4]; bool okt = true;
  for (int t = 0; t < 4; ++t) { vt[t] = get(Form("vel_tpc_%s", tn[t])); if (!vt[t]) okt = false; }
  if (okt) {
    TCanvas* c2 = new TCanvas(Form("cVtpc_%d", nc), "velocity per TPC", 950, 720);
    c2->SetLeftMargin(0.13); c2->SetRightMargin(0.025);
    double m = 0; for (int t = 0; t < 4; ++t) m = std::max(m, vt[t]->GetMaximum());
    vt[0]->SetTitle("Anode-to-cathode drift velocity per TPC;v [mm/#mus];tracks");
    vt[0]->GetXaxis()->CenterTitle(); vt[0]->GetYaxis()->CenterTitle();
    vt[0]->GetXaxis()->SetRangeUser(1.3, 1.9); vt[0]->SetMaximum(1.25 * m);
    TLegend* l = new TLegend(0.52, 0.66, 0.94, 0.90);
    for (int t = 0; t < 4; ++t) {
      vt[t]->SetLineColor(col[t]); vt[t]->SetLineWidth(2);
      vt[t]->Draw(t == 0 ? "HIST" : "HIST SAME");
      l->AddEntry(vt[t], Form("%s (#mu=%.3f)", tn[t], vt[t]->GetMean()), "l");
    }
    l->Draw();
    drawWatermark(tag);
    c2->Print("velocity_tpc.png");
  }

  // --- E-field summary across all measurements vs the nominal field ----------
  const char* en[7] = {"Efield_anode2anode", "Efield_neg", "Efield_pos",
                       "Efield_tpc_Xneg_Zlo", "Efield_tpc_Xneg_Zhi",
                       "Efield_tpc_Xpos_Zlo", "Efield_tpc_Xpos_Zhi"};
  const char* el[7] = {"a2a", "#minusX", "+X", "Xn_Zlo", "Xn_Zhi", "Xp_Zlo", "Xp_Zhi"};
  TH1F* sum = new TH1F("Efield_measurements", "measured E field;;E [kV/cm]", 7, 0, 7);
  sum->SetStats(0);
  for (int i = 0; i < 7; ++i) {
    TH1* h = get(en[i]);
    sum->GetXaxis()->SetBinLabel(i + 1, el[i]);
    if (h && h->GetEntries() > 0) { sum->SetBinContent(i + 1, h->GetMean()); sum->SetBinError(i + 1, h->GetMeanError()); }
  }
  TCanvas* c3 = new TCanvas(Form("cEsum_%d", nc), "E summary", 950, 720);
  c3->SetLeftMargin(0.13); c3->SetRightMargin(0.025);
  sum->SetMarkerStyle(20); sum->SetMarkerSize(1.3); sum->SetLineColor(kBlack);
  sum->GetYaxis()->CenterTitle();
  sum->GetYaxis()->SetRangeUser(std::min(0.40, eNominal - 0.05), std::max(0.60, eNominal + 0.05));
  sum->Draw("E1");
  TLine* ln = new TLine(0, eNominal, 7, eNominal);
  ln->SetLineColor(kGreen + 2); ln->SetLineStyle(2); ln->SetLineWidth(2); ln->Draw();
  TLatex tl; tl.SetNDC(); tl.SetTextFont(42); tl.SetTextColor(kGreen + 2);
  tl.DrawLatex(0.50, 0.30, Form("nominal = %.3f kV/cm", eNominal));
  drawWatermark(tag);
  c3->Print("efield_measurements.png");

  std::cout << "wrote velocity_split.png, velocity_tpc.png, efield_measurements.png" << std::endl;
  f->Close();
}

// -----------------------------------------------------------------------------
// PDHD-only: median drift velocity vs X, split by TPC (Z<230 / Z>230). Each
// curve's X<0 half is the -X-side TPC and its X>0 half the +X-side TPC at that Z,
// so the two curves show all four TPCs.
//   root -l -b -q 'plot_distortions_pdhd.C+' -e 'plot_med_vel_tpc("velocity_pdvd.root","HD")'
// -----------------------------------------------------------------------------
void plot_med_vel_tpc(const char* velfile, const char* detector) {
  TFile* f = TFile::Open(velfile);
  if (f == nullptr || f->IsZombie()) { std::cout << "could not open " << velfile << std::endl; return; }
  TH1* lo = dynamic_cast<TH1*>(f->Get("med_vel_Zlo"));
  TH1* hi = dynamic_cast<TH1*>(f->Get("med_vel_Zhi"));
  if (lo == nullptr || hi == nullptr) { std::cout << "missing med_vel_Zlo/Zhi in " << velfile << std::endl; f->Close(); return; }

  const TString tag = (TString(detector)=="SP") ? "ProtoDUNE-SP"
                    : (TString(detector)=="VD") ? "ProtoDUNE-VD" : "ProtoDUNE-HD";
  useProtoDUNEStyle();
  static int nc = 0; ++nc;
  TCanvas* c = new TCanvas(Form("cMedVelTPC_%d", nc), "med_vel per TPC", 950, 720);
  c->SetLeftMargin(0.13); c->SetRightMargin(0.025);

  lo->SetLineColor(kBlue); lo->SetMarkerColor(kBlue); lo->SetMarkerStyle(20); lo->SetLineWidth(2);
  hi->SetLineColor(kRed);  hi->SetMarkerColor(kRed);  hi->SetMarkerStyle(21); hi->SetLineWidth(2);
  lo->SetTitle("Median drift velocity per TPC;X [cm];v [mm/#mus]");
  lo->GetXaxis()->CenterTitle(); lo->GetYaxis()->CenterTitle();
  lo->GetYaxis()->SetRangeUser(1.2, 1.7);            // clips the empty (0) bins
  lo->Draw("P"); hi->Draw("P SAME");

  TLine* cath = new TLine(0, 1.4, 0, 1.7);           // cathode at X=0
  cath->SetLineStyle(2); cath->SetLineColor(kGray+2); cath->Draw();
  TLegend* l = new TLegend(0.40, 0.76, 0.72, 0.90);
  l->AddEntry(lo, "Z < 230 cm", "lp");
  l->AddEntry(hi, "Z > 230 cm", "lp");
  l->Draw();
  TLatex t; t.SetNDC(); t.SetTextFont(42); t.SetTextSize(0.03);
  t.DrawLatex(0.15, 0.18, "X<0: #minusX side    X>0: +X side   (4 TPCs)");
  drawWatermark(tag);
  c->Print("med_vel_tpc.png");
  std::cout << "wrote med_vel_tpc.png" << std::endl;
  f->Close();
}

// -----------------------------------------------------------------------------
// Detector-face surface distortions from surface_distortions_*.C:
//   distZ_zlo / distZ_zhi : DeltaZ on the upstream (Z=0) / downstream face (X-Y plane)
//   distY_ylo / distY_yhi : DeltaY on the y=ymin / y=ymax face (X-Z plane)
// Each face is one-signed; the COLZ range is chosen from the data (0..25 for a
// positive face, -25..0 for a negative one), read BEFORE SetRangeUser so the
// clamp doesn't feed back into the decision. Palette comes from useProtoDUNEStyle.
//   root -l -b -q 'plot_distortions_pdvd.C+' -e 'plot_surface_distortions("surface_distortions_pdvd.root","VD")'
// -----------------------------------------------------------------------------
void plot_surface_distortions(const char* infile, const char* detector) {
  TFile* f = TFile::Open(infile);
  if (f == nullptr || f->IsZombie()) { std::cout << "could not open " << infile << std::endl; return; }

  const TString tag = (TString(detector) == "SP") ? "ProtoDUNE-SP"
                    : (TString(detector) == "VD") ? "ProtoDUNE-VD" : "ProtoDUNE-HD";
  useProtoDUNEStyle();
  const char* titles[4] = {"#DeltaZ Upstream Face", "#DeltaZ Downstream Face",
                           "#DeltaY -Y Face",  "#DeltaY +Y Face"};
  int i=0;
    const char* zlab[4]   = {"#DeltaZ [cm]", "#DeltaZ [cm]", "#DeltaY [cm]", "#DeltaY [cm]"};
  static int nc = 0;
  int drawn = 0;
  for (const char* hn : {"distZ_zlo", "distZ_zhi", "distY_ylo", "distY_yhi"}) {
    TH1* hd = dynamic_cast<TH1*>(f->Get(hn));
    if (hd == nullptr) { std::cout << "missing " << hn << std::endl; continue; }

    TCanvas* c = new TCanvas(Form("cSurf_%d", ++nc), "surface distortions", 950, 720);
    c->SetLeftMargin(0.13); c->SetRightMargin(0.16);   // room for the COLZ palette
    hd->SetTitle(titles[i]);
        hd->GetZaxis()->SetTitle(zlab[i]);
    hd->GetXaxis()->CenterTitle();
    hd->GetYaxis()->CenterTitle();
    hd->GetZaxis()->CenterTitle();

    const double dmin = hd->GetMinimum();              // read data extremes BEFORE SetRangeUser
    const double dmax = hd->GetMaximum();
    //if (std::fabs(dmin) > std::fabs(dmax)) hd->GetZaxis()->SetRangeUser(-25, 0);
    
                                       hd->GetZaxis()->SetRangeUser(-25, 25);

    hd->Draw("COLZ");
    drawWatermark(tag);
    printBoth(c, Form("surfdist_%s_%s", hn, detector));
    ++drawn;
    delete c;
    i++;
  }

  std::cout << "wrote " << drawn << " surface-distortion faces (pdf+png)" << std::endl;
  f->Close();
}

void plot_distortions_MC(
                         const char* scefile  = "/pnfs/dune/persistent/stash/ProtoDUNE/VD/SpaceCharge/Sim_NoFluidFlow/v2/SCE_PDVD_Sim_NoFluidFlow_E500_v2.root",
                         const char* detector = "VD") {
  const char* c[3] = {"X", "Y", "Z"};
  std::vector<double> planes[3] = { {-300,-150,0,150,300},   // slice X
                                    {-300,-150,0,150,300},   // slice Y
                                    { 50,150,250} };         // slice Z
  for (int i = 0; i < 3; ++i)        // distortion / field component
    for (int s = 0; s < 3; ++s) {    // axis held fixed
      plot_distortion_slices(scefile, c[i], c[s], planes[s], false, "PDVDMC");
      plot_efield_slices(scefile, c[i], c[s], planes[s], "PDVDMC");
    }
}
