#include <TFile.h>
#include <TTree.h>
#include <TBranch.h>
#include <TLeaf.h>
#include <TObjArray.h>
#include <TDirectory.h>
#include <TSystem.h>
#include <TString.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <sstream>
#include <string>
#include <cstring>
#include <vector>

namespace laser_truth {


bool gVerbose = false;

const double kP1HD[3] = {-245.3, 583.8, 352.5};
const double kFallbackDist = 60.0;

struct LaserRow { unsigned long long ts; double pos[3]; double dir[3]; };

struct CsvFile {
  std::string path;
  unsigned long long tFirst = 0, tLast = 0;
  bool loaded = false;
  std::vector<LaserRow> rows;
  long nMatched = 0;
};

std::string norm(const std::string& s) {
  std::string o;
  for (char c : s) if (std::isalnum((unsigned char)c)) o += (char)std::tolower((unsigned char)c);
  if (o.size() > 2 && o.compare(o.size() - 2, 2, "cm") == 0) o.erase(o.size() - 2);
  return o;
}

std::vector<std::string> split(const std::string& line) {
  std::vector<std::string> out;
  std::string f;
  std::stringstream ss(line);
  while (std::getline(ss, f, ';')) out.push_back(f);
  return out;
}

bool blank(const std::string& s) { return s.find_first_not_of(" \t\r\n") == std::string::npos; }

bool rowTs(const std::string& line, unsigned long long& ts) {
  const std::vector<std::string> v = split(line);
  if (v.size() < 4) return false;
  const size_t p = v[0].find_first_not_of(" \t\r");
  if (p == std::string::npos || !std::isdigit((unsigned char)v[0][p])) return false;
  try { ts = std::stoull(v[0].substr(p)); } catch (...) { return false; }
  return true;
}

bool scanRange(CsvFile& c) {
  std::ifstream in(c.path, std::ios::binary);
  if (!in) return false;
  std::string line;
  if (!std::getline(in, line)) return false;
  bool gotFirst = false;
  while (std::getline(in, line)) if (rowTs(line, c.tFirst)) { gotFirst = true; break; }
  if (!gotFirst) return false;
  in.clear();
  in.seekg(0, std::ios::end);
  const std::streamoff size = in.tellg();
  for (std::streamoff back = 65536; ; back *= 4) {
    const std::streamoff n = std::min<std::streamoff>(size, back);
    in.clear();
    in.seekg(size - n);
    const std::string tail((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
    std::vector<std::string> lines;
    std::stringstream ts(tail);
    while (std::getline(ts, line)) lines.push_back(line);
    const size_t firstFull = (n < size) ? 1 : 0;
    for (size_t k = lines.size(); k-- > firstFull; )
      if (rowTs(lines[k], c.tLast)) return true;
    if (n == size) return false;
  }
}

bool loadRows(CsvFile& c) {
  if (c.loaded) return !c.rows.empty();
  c.loaded = true;
  std::ifstream in(c.path);
  std::string line;
  if (!std::getline(in, line)) return false;
  const std::vector<std::string> h = split(line);
  std::map<std::string, int> col;
  for (size_t i = 0; i < h.size(); ++i) col[norm(h[i])] = (int)i;
  auto find = [&](const std::vector<std::string>& names) -> int {
    for (const auto& n : names) { auto it = col.find(n); if (it != col.end()) return it->second; }
    return -1;
  };
  const char* ax[3] = {"x", "y", "z"};
  int ip[3], id[3];
  const int it = find({"timestamp"});
  int imax = it;
  bool ok = (it >= 0);
  for (int a = 0; a < 3; ++a) {
    ip[a] = find({std::string("laserpos") + ax[a] + "larsoft"});
    id[a] = find({std::string("laserdirlarsoft") + ax[a], std::string("laserdir") + ax[a]});
    ok = ok && ip[a] >= 0 && id[a] >= 0;
    imax = std::max(imax, std::max(ip[a], id[a]));
  }
  if (!ok) {
    if (gVerbose) std::cout << "  " << c.path << ": missing TIMESTAMP / LaserPos*LArSoft / LaserDir* column, skipped" << std::endl;
    return false;
  }
  while (std::getline(in, line)) {
    if (blank(line)) continue;
    unsigned long long tsCheck;
    if (!rowTs(line, tsCheck)) continue;
    const std::vector<std::string> v = split(line);
    if ((int)v.size() <= imax) continue;
    LaserRow r;
    try {
      r.ts = std::stoull(v[it]);
      for (int a = 0; a < 3; ++a) { r.pos[a] = std::stod(v[ip[a]]); r.dir[a] = std::stod(v[id[a]]); }
    } catch (...) { continue; }
    c.rows.push_back(r);
  }
  std::sort(c.rows.begin(), c.rows.end(),
            [](const LaserRow& a, const LaserRow& b) { return a.ts < b.ts; });
  if (gVerbose) std::cout << "  loaded " << c.rows.size() << " rows from " << c.path << std::endl;
  return !c.rows.empty();
}

unsigned long long absdiff(unsigned long long a, unsigned long long b) { return a > b ? a - b : b - a; }

const LaserRow* nearest(const CsvFile& c, unsigned long long ts) {
  if (c.rows.empty()) return nullptr;
  auto itr = std::lower_bound(c.rows.begin(), c.rows.end(), ts,
                              [](const LaserRow& r, unsigned long long t) { return r.ts < t; });
  if (itr == c.rows.end()) return &c.rows.back();
  if (itr == c.rows.begin()) return &*itr;
  auto prev = itr - 1;
  return (absdiff(prev->ts, ts) <= absdiff(itr->ts, ts)) ? &*prev : &*itr;
}

struct Box { double lo[3], hi[3]; };
const Box kBoxHD = {{-360.0, 7.1, 0.46}, {360.0, 607.0, 460.46}};
const Box kBoxVD = {{-341.56, -340.0, 0.0}, {341.56, 340.0, 300.0}};

void exitPoint(const double p[3], const double d[3], const Box& b, Float_t out[3]) {
  double tmax = 1e300;
  for (int a = 0; a < 3; ++a) {
    if (std::fabs(d[a]) < 1e-12) continue;
    tmax = std::min(tmax, std::max((b.lo[a] - p[a]) / d[a], (b.hi[a] - p[a]) / d[a]));
  }
  for (int a = 0; a < 3; ++a) out[a] = p[a] + tmax * d[a];
}

void entryPoint(const double p[3], const double d[3], const Box& b, Float_t out[3]) {
  double tmin = 0.0;
  for (int a = 0; a < 3; ++a) {
    if (std::fabs(d[a]) < 1e-12) continue;
    tmin = std::max(tmin, std::min((b.lo[a] - p[a]) / d[a], (b.hi[a] - p[a]) / d[a]));
  }
  for (int a = 0; a < 3; ++a) out[a] = p[a] + tmin * d[a];
}

double lineDist(const double p[3], const double d[3], const double q[3]) {
  const double v[3] = {q[0]-p[0], q[1]-p[1], q[2]-p[2]};
  const double c[3] = {v[1]*d[2]-v[2]*d[1], v[2]*d[0]-v[0]*d[2], v[0]*d[1]-v[1]*d[0]};
  const double nd = std::sqrt(d[0]*d[0]+d[1]*d[1]+d[2]*d[2]);
  return nd > 0 ? std::sqrt(c[0]*c[0]+c[1]*c[1]+c[2]*c[2]) / nd : 1e30;
}

std::vector<std::string> listCsv(const std::string& dir) {
  std::vector<std::string> out;
  void* dp = gSystem->OpenDirectory(dir.c_str());
  if (!dp) return out;
  while (const char* e = gSystem->GetDirEntry(dp)) {
    const std::string n(e);
    if (n.size() > 4 && n.compare(n.size() - 4, 4, ".csv") == 0) out.push_back(dir + "/" + n);
  }
  gSystem->FreeDirectory(dp);
  std::sort(out.begin(), out.end());
  return out;
}

std::vector<std::string> listSubdirs(const std::string& dir, const std::string& prefix) {
  std::vector<std::string> out;
  void* dp = gSystem->OpenDirectory(dir.c_str());
  if (!dp) return out;
  while (const char* e = gSystem->GetDirEntry(dp)) {
    const std::string n(e);
    if (n.compare(0, prefix.size(), prefix) == 0) out.push_back(dir + "/" + n);
  }
  gSystem->FreeDirectory(dp);
  std::sort(out.begin(), out.end());
  return out;
}

}

void add_laser_truth(const char* infile,
                     const char* detector = "HD")
{
  const int         periscope   = -1;
  const char*       csvBase     = "/exp/dune/data/users/rdiurba";
  const Long64_t    maxDtTicks  = 0;
  const bool        verbose     = false;
  const int         cutRun      = 31437;
  const int         minEvent    = 6710;
  const double      maxAnodeDeg = 30.0;
  const double      minCos      = 0.95;
  const double      maxStartDist = 25.0;
  const TString tmpName = TString(infile) + ".tmp";
  const char* outfile = tmpName.Data();
  using namespace laser_truth;
  gVerbose = verbose;
  const TString det(detector);
  const bool isVD = det.EqualTo("VD", TString::kIgnoreCase);
  const TString treeDir = isVD ? "laseranalysispdvd" : "laseranalysispdhd";
  const Box& box = isVD ? kBoxVD : kBoxHD;

  std::string base(csvBase);
  while (base.size() > 1 && base.back() == '/') base.pop_back();
  const std::string detDir = base + "/laserFiles" + (isVD ? "VD" : "HD");
  std::vector<std::string> dirs;
  if (periscope >= 0) dirs.push_back(detDir + "/periscope" + std::to_string(periscope));
  else {
    dirs.push_back(detDir);
    for (const auto& d : listSubdirs(detDir, "periscope")) dirs.push_back(d);
  }
  std::vector<CsvFile> csvs;
  for (const auto& d : dirs)
    for (const auto& p : listCsv(d)) {
      CsvFile c; c.path = p;
      if (scanRange(c)) csvs.push_back(c);
      else if (gVerbose) std::cout << "  could not read timestamps from " << p << ", skipped" << std::endl;
    }
  if (csvs.empty()) { std::cout << "no laser CSVs found under " << detDir << std::endl; return; }
  if (gVerbose) std::cout << csvs.size() << " laser CSV(s):" << std::endl;
  if (gVerbose) for (const auto& c : csvs)
    std::cout << "  " << c.path << "  TIMESTAMP " << c.tFirst << " .. " << c.tLast << std::endl;

  TFile* fin = TFile::Open(infile, "READ");
  if (!fin || fin->IsZombie()) { std::cout << "could not open " << infile << std::endl; return; }
  TTree* t = dynamic_cast<TTree*>(fin->Get(treeDir + "/Event"));
  if (!t) { std::cout << "no " << treeDir << "/Event in " << infile << std::endl; fin->Close(); return; }

  const int maxTracks = std::max(1, (int)t->GetMaximum("cross_trks"));
  struct Buf { std::vector<char> data; size_t stride = 0; size_t lenType = 0; bool perTrack = false; };
  std::map<std::string, Buf> bufs;
  TObjArray* brs = t->GetListOfBranches();
  for (int ib = 0; ib < brs->GetEntries(); ++ib) {
    TBranch* br = (TBranch*)brs->At(ib);
    TLeaf* lf = (TLeaf*)br->GetListOfLeaves()->At(0);
    if (!lf) continue;
    Buf& b = bufs[br->GetName()];
    b.perTrack = (lf->GetLeafCount() != nullptr);
    b.stride   = (size_t)lf->GetLenType() * lf->GetLenStatic();
    b.lenType  = (size_t)lf->GetLenType();
    b.data.assign(b.stride * (b.perTrack ? maxTracks : 1), 0);
    t->SetBranchAddress(br->GetName(), (void*)b.data.data());
  }
  const char* need[] = {"triggertime","cross_trks","trkstartx","trkstarty","trkstartz","trkendx","trkendy","trkendz"};
  for (const char* nm : need)
    if (!bufs.count(nm)) { std::cout << "branch " << nm << " missing in " << infile << std::endl; fin->Close(); return; }

  const ULong64_t& triggertime = *reinterpret_cast<ULong64_t*>(bufs["triggertime"].data.data());
  const ULong64_t* rdtsBranch = bufs.count("rdtimestamp")
      ? reinterpret_cast<ULong64_t*>(bufs["rdtimestamp"].data.data()) : nullptr;
  ULong64_t rdtimestamp = 0;
  Int_t&     cross_trks  = *reinterpret_cast<Int_t*>(bufs["cross_trks"].data.data());
  const Float_t* sx = reinterpret_cast<Float_t*>(bufs["trkstartx"].data.data());
  const Float_t* sy = reinterpret_cast<Float_t*>(bufs["trkstarty"].data.data());
  const Float_t* sz = reinterpret_cast<Float_t*>(bufs["trkstartz"].data.data());
  const Float_t* ex = reinterpret_cast<Float_t*>(bufs["trkendx"].data.data());
  const Float_t* ey = reinterpret_cast<Float_t*>(bufs["trkendy"].data.data());
  const Float_t* ez = reinterpret_cast<Float_t*>(bufs["trkendz"].data.data());
  TBranch* bCross = t->GetBranch("cross_trks");
  const Int_t* runNum = bufs.count("run")   ? reinterpret_cast<Int_t*>(bufs["run"].data.data())   : nullptr;
  const Int_t* evtNum = bufs.count("event") ? reinterpret_cast<Int_t*>(bufs["event"].data.data()) : nullptr;

  TFile* fout = TFile::Open(outfile, "RECREATE");
  if (!fout || fout->IsZombie()) { std::cout << "could not create " << outfile << std::endl; fin->Close(); return; }
  TDirectory* outDir = fout->mkdir(treeDir);
  outDir->cd();
  TTree* out = t->CloneTree(0);
  Float_t start[3], end[3], dir[3], thetaxz, thetayz, thetaxy, trkDist;
  Int_t timed;
  out->Branch("laser_startx", &start[0], "laser_startx/F");
  out->Branch("laser_starty", &start[1], "laser_starty/F");
  out->Branch("laser_startz", &start[2], "laser_startz/F");
  out->Branch("laser_endx",   &end[0],   "laser_endx/F");
  out->Branch("laser_endy",   &end[1],   "laser_endy/F");
  out->Branch("laser_endz",   &end[2],   "laser_endz/F");
  out->Branch("laser_dirx",   &dir[0],   "laser_dirx/F");
  out->Branch("laser_diry",   &dir[1],   "laser_diry/F");
  out->Branch("laser_dirz",   &dir[2],   "laser_dirz/F");
  out->Branch("laser_thetaxz", &thetaxz, "laser_thetaxz/F");
  out->Branch("laser_thetayz", &thetayz, "laser_thetayz/F");
  out->Branch("laser_thetaxy", &thetaxy, "laser_thetaxy/F");
  out->Branch("laser_trk_dist", &trkDist, "laser_trk_dist/F");
  out->Branch("laser_timed",    &timed,   "laser_timed/I");

  long nMatch = 0, nNoRow = 0, nNoTrk = 0, nFallback = 0, nEvCut = 0, nNoGood = 0;
  const Long64_t n = t->GetEntries();
  for (Long64_t i = 0; i < n; ++i) {
    bCross->GetEntry(i);
    if (cross_trks <= 0 || cross_trks > maxTracks) { ++nNoTrk; continue; }
    t->GetEntry(i);
    if (runNum && evtNum && *runNum == cutRun && *evtNum <= minEvent) { ++nEvCut; continue; }

    const bool dirCut = runNum && *runNum == cutRun;
    const double sinMax = std::sin(maxAnodeDeg * M_PI / 180.0);
    auto okDirPts = [&](const double s3[3], const double e3[3]) {
      if (!dirCut) return true;
      if (!(e3[2] < s3[2])) return false;
      const double dx = e3[0]-s3[0], dy = e3[1]-s3[1], dz = e3[2]-s3[2];
      const double len = std::sqrt(dx*dx + dy*dy + dz*dz);
      if (len <= 0) return false;
      const double towardsAnode = (s3[0] >= 0 ? dx : -dx) / len;
      return towardsAnode <= sinMax;
    };
    auto okDir = [&](int j) {
      const double s3[3] = {sx[j], sy[j], sz[j]}, e3[3] = {ex[j], ey[j], ez[j]};
      return okDirPts(s3, e3);
    };
    rdtimestamp = (rdtsBranch && *rdtsBranch != 0)
        ? *rdtsBranch
        : (triggertime >> 32) * 62500000ULL + ((triggertime & 0xffffffffULL) * 62500000ULL) / 1000000000ULL;

    const LaserRow* best = nullptr;
    LaserRow fb;
    CsvFile* bestCsv = nullptr;
    unsigned long long bestDt = 0;
    if (rdtimestamp != 0)
      for (auto& c : csvs) {
        if (rdtimestamp + (unsigned long long)maxDtTicks < c.tFirst) continue;
        if (rdtimestamp > c.tLast + (unsigned long long)maxDtTicks) continue;
        if (!loadRows(c)) continue;
        const LaserRow* r = nearest(c, rdtimestamp);
        if (!r) continue;
        const unsigned long long dt = absdiff(r->ts, rdtimestamp);
        if (!best || dt < bestDt) { best = r; bestCsv = &c; bestDt = dt; }
      }
    if (!best || bestDt > (unsigned long long)maxDtTicks) {
      if (gVerbose && nNoRow < 10) {
        bool covered = false;
        for (const auto& c : csvs) if (rdtimestamp >= c.tFirst && rdtimestamp <= c.tLast) covered = true;
        std::cout << "  no match: entry " << i << " rdtimestamp " << rdtimestamp
                  << (covered ? "  (inside a CSV range; nearest row dt = " + std::to_string(best ? (long long)bestDt : -1LL) + " ticks)"
                              : "  (not inside any CSV time range)") << std::endl;
      }

      int jFb = -1; double dFb = 1e30;
      if (!isVD)
        for (int j = 0; j < cross_trks; ++j) {
          if (!okDir(j)) continue;
          const double dx = sx[j]-kP1HD[0], dy = sy[j]-kP1HD[1], dz = sz[j]-kP1HD[2];
          const double d = std::sqrt(dx*dx + dy*dy + dz*dz);
          if (d < dFb) { dFb = d; jFb = j; }
        }
      if (jFb < 0 || dFb > kFallbackDist) { ++nNoRow; continue; }
      fb.ts = 0;
      const double te[3] = {ex[jFb], ey[jFb], ez[jFb]};
      for (int a = 0; a < 3; ++a) { fb.pos[a] = kP1HD[a]; fb.dir[a] = te[a] - kP1HD[a]; }
      best = &fb; bestCsv = nullptr;
      timed = 0; ++nFallback;
    } else timed = 1;

    entryPoint(best->pos, best->dir, box, start);
    exitPoint(best->pos, best->dir, box, end);
    const double ls[3] = {start[0], start[1], start[2]};
    const double ld[3] = {end[0]-start[0], end[1]-start[1], end[2]-start[2]};
    const double nd = std::sqrt(best->dir[0]*best->dir[0] + best->dir[1]*best->dir[1] + best->dir[2]*best->dir[2]);
    for (int a = 0; a < 3; ++a) dir[a] = nd > 0 ? best->dir[a] / nd : -99999;
    thetaxz = std::atan2(best->dir[0], best->dir[2]);
    thetayz = std::atan2(best->dir[1], best->dir[2]);
    thetaxy = std::atan2(best->dir[1], best->dir[0]);

    int jBest = -1; double dBest = 1e30; bool flipBest = false;
    for (int j = 0; j < cross_trks; ++j) {
      if (timed == 1) {
        const double tx = ex[j]-sx[j], ty = ey[j]-sy[j], tz = ez[j]-sz[j];
        const double tl = std::sqrt(tx*tx + ty*ty + tz*tz);
        if (tl <= 0 || nd <= 0) continue;
        const double c = (tx*best->dir[0] + ty*best->dir[1] + tz*best->dir[2]) / (tl * nd);
        if (std::fabs(c) < minCos) continue;
        const bool flip = (c < 0);
        const double a3[3] = {sx[j], sy[j], sz[j]}, b3[3] = {ex[j], ey[j], ez[j]};
        if (!(flip ? okDirPts(b3, a3) : okDirPts(a3, b3))) continue;
        if (lineDist(ls, ld, flip ? b3 : a3) > maxStartDist) continue;
        const double d = std::max(lineDist(best->pos, best->dir, a3), lineDist(best->pos, best->dir, b3));
        if (d < dBest) { dBest = d; jBest = j; flipBest = flip; }
      } else {
        if (!okDir(j)) continue;
        const double dx = sx[j]-best->pos[0], dy = sy[j]-best->pos[1], dz = sz[j]-best->pos[2];
        const double d = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (d < dBest) { dBest = d; jBest = j; flipBest = false; }
      }
    }
    if (jBest < 0) { if (timed == 1) ++nNoGood; else ++nNoTrk; continue; }

    for (auto& kv : bufs) {
      Buf& b = kv.second;
      if (!b.perTrack || jBest == 0) continue;
      std::memmove(b.data.data(), b.data.data() + (size_t)jBest * b.stride, b.stride);
    }
    cross_trks = 1;

    if (flipBest) {
      auto F = [&](const char* nm) -> Float_t* { return bufs.count(nm) ? reinterpret_cast<Float_t*>(bufs[nm].data.data()) : nullptr; };
      const char* se[3][2] = {{"trkstartx","trkendx"}, {"trkstarty","trkendy"}, {"trkstartz","trkendz"}};
      for (auto& p : se) { Float_t* a = F(p[0]); Float_t* b = F(p[1]); if (a && b) std::swap(a[0], b[0]); }
      Float_t* cs = F("trkstartcosxyz"); Float_t* ce = F("trkendcosxyz");
      if (cs && ce) for (int a = 0; a < 3; ++a) { const Float_t t0 = cs[a]; cs[a] = -ce[a]; ce[a] = -t0; }
      if (cs) {
        if (Float_t* txz = F("trackthetaxz")) txz[0] = std::atan2(cs[0], cs[2]);
        if (Float_t* tyz = F("trackthetayz")) tyz[0] = std::atan2(cs[1], cs[2]);
      }

      const Int_t* nh = bufs.count("ntrkhits") ? reinterpret_cast<Int_t*>(bufs["ntrkhits"].data.data()) : nullptr;
      for (auto& kv : bufs) {
        Buf& b = kv.second;
        if (!b.perTrack || !nh || b.lenType == 0 || b.stride != b.lenType * 3 * 3000) continue;
        for (int pl = 0; pl < 3; ++pl) {
          const int nn = std::min<int>(nh[pl], 3000);
          if (nn <= 1) continue;
          char* base = b.data.data() + (size_t)pl * 3000 * b.lenType;
          for (int lo = 0, hi = nn - 1; lo < hi; ++lo, --hi)
            std::swap_ranges(base + lo * b.lenType, base + (lo + 1) * b.lenType, base + hi * b.lenType);
        }
      }
      if (bufs.count("trkflipped")) { Int_t* fl = reinterpret_cast<Int_t*>(bufs["trkflipped"].data.data()); fl[0] = (fl[0] == 1) ? 0 : 1; }
    }
    trkDist = dBest;
    ++nMatch; if (bestCsv) ++bestCsv->nMatched;
    out->Fill();
  }

  outDir->cd();
  out->Write();
  fout->Close();
  fin->Close();
  if (gSystem->Rename(outfile, infile) != 0) { std::cout << "could not overwrite " << infile << std::endl; return; }

  std::cout << "events: " << n << "  written (matched): " << nMatch
            << " (of which P1 position fallback, no time match: " << nFallback << ")"
            << "  | not written: no CSV row " << (maxDtTicks > 0 ? "within " + std::to_string(maxDtTicks) + " ticks" : std::string("with the same timestamp")) << ": " << nNoRow
            << ", run " << cutRun << " event <= " << minEvent << ": " << nEvCut
            << ", CSV match but no track with |cos| >= " << minCos << " and start within " << maxStartDist << " cm of the laser line: " << nNoGood
            << ", no saved track: " << nNoTrk << std::endl;
  if (gVerbose) for (const auto& c : csvs) if (c.nMatched) std::cout << "  " << c.nMatched << " from " << c.path << std::endl;
  if (gVerbose) std::cout << "wrote " << infile << ":" << treeDir << "/Event" << std::endl;
}
