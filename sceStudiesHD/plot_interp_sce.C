// plot_interp_sce.C
//
// Plots of the INTERPOLATED (surface-scaled, data-driven) SCE file made by
// interpolate_sce_bulk_pdhd.C. Every title starts with "Interpolated" and every
// file name starts with "interp_" + the quantity, e.g.
//   interp_dX_Xm150.png   Delta X,  slice at X = -150 cm
//   interp_EZ_Y300.png    E_z,      slice at Y = 300 cm
//   interp_ratio_Zlo.png  data/MC ratio on the upstream (Z low) face
// (pdf as well as png). Slices use the same 20 cm grid, style, palette and
// watermark as plot_distortions.C; Pos/Neg maps are picked by the sign of X.
//
//   root -l -b -q 'plot_interp_sce.C("SCE_PDHD_DataDriven_surface.root","HD")'

#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TH2.h>
#include <TH3.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>

#include <iostream>
#include <vector>

namespace interp_plot {

const double kXmin = -360, kXmax = 360, kYmin = 0, kYmax = 600, kZmin = 0, kZmax = 460, kBin = 20;

// CET-CBTD1 cyan-white-red diverging palette (same as plot_distortions.C)
void setPalette() {
  const int n = 9;
  double stops[n] = {0.000, 0.125, 0.250, 0.375, 0.500, 0.625, 0.750, 0.875, 1.000};
  double red[n]   = {0.000, 0.400, 0.620, 0.800, 0.960, 0.980, 0.990, 0.980, 0.960};
  double green[n] = {0.770, 0.830, 0.880, 0.920, 0.960, 0.850, 0.720, 0.600, 0.490};
  double blue[n]  = {0.880, 0.900, 0.930, 0.940, 0.950, 0.810, 0.660, 0.540, 0.430};
  TColor::CreateGradientColorTable(n, stops, red, green, blue, 255);
  gStyle->SetNumberContours(255);
}

void style() {
  static bool loaded = false;
  if (!loaded) { gROOT->LoadMacro("protoDUNEStyle.C"); loaded = true; }
  gROOT->SetStyle("protoDUNEStyle");
  gROOT->ForceStyle();
  gStyle->SetTitleX(0.35);
  setPalette();
}

TString detTag(const char* d) {
  const TString s(d);
  return s == "SP" ? "ProtoDUNE-SP" : (s == "VD" ? "ProtoDUNE-VD" : "ProtoDUNE-HD");
}

void watermark(const TString& tag) {
  TLatex t; t.SetNDC(); t.SetTextFont(42);
  t.DrawLatex(0.20, 0.94, Form("#bf{DUNE:%s}", tag.Data()));
}

TString pos(double p) { return p < 0 ? Form("m%.0f", -p) : Form("%.0f", p); }

void save(TCanvas* c, const TString& base) { c->Print(base + ".png"); c->Print(base + ".pdf"); }

void geom(int a, int& n, double& lo, double& hi) {
  if (a == 0) { lo = kXmin; hi = kXmax; } else if (a == 1) { lo = kYmin; hi = kYmax; } else { lo = kZmin; hi = kZmax; }
  n = int((hi - lo) / kBin);
}

TCanvas* canvas(const char* name) {
  static int nc = 0;
  TCanvas* c = new TCanvas(Form("%s_%d", name, ++nc), name, 950, 720);
  c->SetLeftMargin(0.13); c->SetRightMargin(0.16);
  return c;
}

// slices of a Pos/Neg pair of 3D maps. file names: interp_<qty>_<axis><pos>
void slices(TFile* f, const char* posName, const char* negName,
            const TString& qty, const TString& title, const TString& ztitle, double zmax,
            const char* det) {
  TH3* hp = dynamic_cast<TH3*>(f->Get(posName));
  TH3* hn = dynamic_cast<TH3*>(f->Get(negName));
  if (!hp || !hn) { std::cout << "missing " << posName << " / " << negName << std::endl; return; }

  const char* an[3] = {"X", "Y", "Z"};
  const std::vector<double> planes[3] = {{-300, -150, 0, 150, 300},
                                         {100, 200, 300, 400, 500},
                                         {50, 150, 250, 350, 450}};
  TCanvas* c = canvas("cInterp");
  for (int s = 0; s < 3; ++s) {
    const int ih = (s == 0) ? 2 : 0, iv = (s == 1) ? 2 : 1;   // X horizontal whenever shown
    int nh, nv; double hlo, hhi, vlo, vhi;
    geom(ih, nh, hlo, hhi); geom(iv, nv, vlo, vhi);
    for (double pv : planes[s]) {
      TH2D* h = new TH2D("h_interp_slice", "", nh, hlo, hhi, nv, vlo, vhi);
      for (int a = 1; a <= nh; ++a)
        for (int b = 1; b <= nv; ++b) {
          double p[3];
          p[s] = pv; p[ih] = hlo + (a - 0.5) * kBin; p[iv] = vlo + (b - 0.5) * kBin;
          TH3* m = (p[0] >= 0) ? hp : hn;
          h->SetBinContent(a, b, m->GetBinContent(m->FindBin(p[0], p[1], p[2])));
        }
      h->SetTitle(Form("Interpolated %s  |  %s = %.0f cm;%s [cm];%s [cm];%s",
                       title.Data(), an[s], pv, an[ih], an[iv], ztitle.Data()));
      h->GetXaxis()->CenterTitle(); h->GetYaxis()->CenterTitle(); h->GetZaxis()->CenterTitle();
      h->GetZaxis()->SetRangeUser(-zmax, zmax);
      c->cd(); h->Draw("COLZ"); watermark(detTag(det));
      save(c, Form("interp_%s_%s%s", qty.Data(), an[s], pos(pv).Data()));
      delete h;
    }
  }
  delete c;
}

// data/MC face-ratio maps. file names: interp_ratio_<face>
void ratios(TFile* f, const char* det) {
  const char* names[4] = {"ratioZ_zlo", "ratioZ_zhi", "ratioY_ylo", "ratioY_yhi"};
  const char* faces[4] = {"Zlo", "Zhi", "Ylo", "Yhi"};
  const char* titles[4] = {"data/MC #DeltaZ, upstream face", "data/MC #DeltaZ, downstream face",
                           "data/MC #DeltaY, bottom face", "data/MC #DeltaY, top face"};
  TCanvas* c = canvas("cInterpRatio");
  for (int i = 0; i < 4; ++i) {
    TH2* h = dynamic_cast<TH2*>(f->Get(names[i]));
    if (!h) { std::cout << "missing " << names[i] << std::endl; continue; }
    h->SetTitle(Form("Interpolated %s", titles[i]));
    h->GetZaxis()->SetTitle("data / MC");
    h->GetXaxis()->CenterTitle(); h->GetYaxis()->CenterTitle(); h->GetZaxis()->CenterTitle();
    h->GetZaxis()->SetRangeUser(0, 2);        // white at 1 = data equals MC
    c->cd(); h->Draw("COLZ"); watermark(detTag(det));
    save(c, Form("interp_ratio_%s", faces[i]));
  }
  delete c;
}

}  // namespace interp_plot

void plot_interp_sce(const char* infile = "SCE_PDHD_DataDriven_surface.root",
                     const char* detector = "HD") {
  using namespace interp_plot;
  TFile* f = TFile::Open(infile);
  if (!f || f->IsZombie()) { std::cout << "could not open " << infile << std::endl; return; }
  style();

  const char* c[3] = {"X", "Y", "Z"};
  for (int i = 0; i < 3; ++i) {
    // RecoBkwd displacement: interp_d<comp>_...
    slices(f, Form("RecoBkwd_Displacement_%s_Pos", c[i]), Form("RecoBkwd_Displacement_%s_Neg", c[i]),
           Form("d%s", c[i]), Form("#Delta%s (RecoBkwd)", c[i]), Form("#Delta%s [cm]", c[i]), 30, detector);
    // E field from the interpolated displacements: interp_E<comp>_...
    slices(f, Form("Reco_ElecField_%s_Pos", c[i]), Form("Reco_ElecField_%s_Neg", c[i]),
           Form("E%s", c[i]), Form("E_{%s}", c[i]),
           i == 0 ? "(E_{X} #minus E_{0})/E_{0}" : Form("E_{%s}/E_{0}", c[i]), 0.5, detector);
  }
  ratios(f, detector);   // interp_ratio_<face>

  f->Close();
  std::cout << "wrote interp_d{X,Y,Z}_*, interp_E{X,Y,Z}_*, interp_ratio_* (png+pdf)" << std::endl;
}
