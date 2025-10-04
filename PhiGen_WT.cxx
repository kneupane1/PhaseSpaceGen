// PhiGen_truth.cxx  (no detector smearing)
#include <cmath>
#include <fstream>
#include <iostream>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGenPhaseSpace.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TStyle.h"

// ===== masses (GeV) =====
constexpr double MASS_E = 0.000511;
constexpr double MASS_P = 0.938272;
constexpr double MASS_PHI = 1.019461;
constexpr double WIDTH_PHI = 0.00426;  // natural width ~4.26 MeV
constexpr double MASS_K = 0.493677;

// ===== small helpers =====
namespace physics {
inline double nu_from_W_Q2(double W, double Q2) { return (W * W + Q2 - MASS_P * MASS_P) / (2.0 * MASS_P); }
inline double Q2_calc(const TLorentzVector& k, const TLorentzVector& kprime) {
  TLorentzVector q = k - kprime;
  return -q.M2();
}
inline double W_calc(const TLorentzVector& k, const TLorentzVector& kprime) {
  TLorentzVector q = k - kprime;
  TLorentzVector had = q + TLorentzVector(0, 0, 0, MASS_P);
  return had.M();
}
}  // namespace physics

// sample phi mass from a narrow Breit–Wigner (accept–reject in ±50 MeV)
double sample_phi_mass(double m0, double gamma) {
  double x, y;
  do {
    x = gRandom->Uniform(m0 - 0.05, m0 + 0.05);
    if (x < 2 * MASS_K) continue;
    y = gRandom->Uniform(0.0, 1.0);
  } while (y > (gamma / 2.0) / ((x - m0) * (x - m0) + (gamma * gamma / 4.0)));
  return x;
}

int main(int argc, char* argv[]) {
  // ---------- settings ----------
  int Nevents = 5000;    // number of ACCEPTED events
  double Ebeam = 7.546;  // GeV
  double Q2_min = 1.0, Q2_max = 8.0;
  double W_min = 1.6, W_max = 3.0;

  double frac_bkg = 0.45;     // fraction of non-resonant KK
  double cont_lambda = 25.0;  // slope for continuum
  double t_slope_b = 3.0;     // GeV^-2, exp(b t)

  TRandom3 rng(0);
  gRandom = &rng;

  // initial state (lab)
  TLorentzVector target(0, 0, 0, MASS_P);
  TLorentzVector beam(0, 0, std::sqrt(Ebeam * Ebeam - MASS_E * MASS_E), Ebeam);

  std::ofstream out("phi_gen_truth.lund");

  // histos
  TH1D* hW = new TH1D("W", "W distribution;W (GeV);Counts", 200, 1.2, 3.2);
  TH1D* hQ2 = new TH1D("Q2", "Q^{2} distribution;Q^{2} (GeV^{2});Counts", 200, 0.0, 10.0);
  TH2D* hWQ2 = new TH2D("WvsQ2", "W vs Q^{2};W (GeV);Q^{2} (GeV^{2})", 160, 1.2, 3.2, 160, 0.0, 10.0);
  TH1D* hMkk = new TH1D("Mkk", "K^{+}K^{-} invariant mass;M(K^{+}K^{-}) [GeV];Counts", 240, 0.96, 1.12);

  // flux weight ~ 1/(W*Q2)
  auto fluxWeight = [](double W, double Q2) { return 1.0 / (W * Q2); };
  const double wmax = fluxWeight(W_min, Q2_min);

  int acc = 0, trials = 0;

  while (acc < Nevents) {
    ++trials;

    // sample (W,Q2) uniformly, accept–reject with flux
    double W = gRandom->Uniform(W_min, W_max);
    double Q2 = gRandom->Uniform(Q2_min, Q2_max);
    double w = fluxWeight(W, Q2);
    if (gRandom->Uniform() > w / wmax) continue;

    // electron kinematics
    double nu = physics::nu_from_W_Q2(W, Q2);
    double Eprime = Ebeam - nu;
    if (Eprime <= MASS_E) continue;

    double costh = 1.0 - Q2 / (2.0 * Ebeam * Eprime);
    if (costh < -1.0 || costh > 1.0) continue;
    double theta = std::acos(costh);
    double phi_e = gRandom->Uniform(0, 2 * TMath::Pi());

    double pe = std::sqrt(std::max(0.0, Eprime * Eprime - MASS_E * MASS_E));
    TLorentzVector eprime(pe * std::sin(theta) * std::cos(phi_e), pe * std::sin(theta) * std::sin(phi_e),
                          pe * std::cos(theta), Eprime);

    // hadronic system X (mass=W)
    TLorentzVector X = (beam + TLorentzVector(0, 0, 0, MASS_P)) - eprime;

    bool isBkg = (gRandom->Uniform() < frac_bkg);
    TLorentzVector pprime, Kplus, Kminus, phi4;
    bool ok = true;

    if (!isBkg) {
      // Resonant φ
      double mphi = sample_phi_mass(MASS_PHI, WIDTH_PHI);
      if (W < MASS_P + mphi) continue;

      TGenPhaseSpace two;
      Double_t m2[2] = {MASS_P, mphi};
      ok = two.SetDecay(X, 2, m2) && two.Generate();
      if (!ok) continue;

      pprime = *two.GetDecay(0);
      phi4 = *two.GetDecay(1);

      TGenPhaseSpace phidec;
      Double_t km[2] = {MASS_K, MASS_K};
      ok = phidec.SetDecay(phi4, 2, km) && phidec.Generate();
      if (!ok) continue;

      Kplus = *phidec.GetDecay(0);
      Kminus = *phidec.GetDecay(1);
    } else {
      // Continuum KK
      double M_lo = 0.99, M_hi = 1.12;
      double Mkk, y;
      const double ymax = 1.0;
      do {
        Mkk = gRandom->Uniform(M_lo, M_hi);
        y = gRandom->Uniform(0.0, 1.0);
      } while (y > std::exp(-cont_lambda * (Mkk - M_lo)) / ymax);

      if (W < MASS_P + Mkk) continue;

      TGenPhaseSpace two;
      Double_t m2[2] = {MASS_P, Mkk};
      ok = two.SetDecay(X, 2, m2) && two.Generate();
      if (!ok) continue;

      pprime = *two.GetDecay(0);
      TLorentzVector KK = *two.GetDecay(1);

      TGenPhaseSpace kkdec;
      Double_t km[2] = {MASS_K, MASS_K};
      ok = kkdec.SetDecay(KK, 2, km) && kkdec.Generate();
      if (!ok) continue;

      Kplus = *kkdec.GetDecay(0);
      Kminus = *kkdec.GetDecay(1);
      phi4 = KK;
    }

    // apply t-slope weight
    TLorentzVector q4 = beam - eprime;
    double t = (q4 - phi4).M2();
    double wt = std::exp(t_slope_b * t);
    if (gRandom->Uniform() > wt) continue;

    // invariant mass
    double Mkk = (Kplus + Kminus).M();

    hW->Fill(W);
    hQ2->Fill(Q2);
    hWQ2->Fill(W, Q2);
    hMkk->Fill(Mkk);

    out << "Event " << acc + 1 << "\n";
    out << "e' " << eprime.Px() << " " << eprime.Py() << " " << eprime.Pz() << " " << eprime.E() << "\n";
    out << "p' " << pprime.Px() << " " << pprime.Py() << " " << pprime.Pz() << " " << pprime.E() << "\n";
    out << "K+ " << Kplus.Px() << " " << Kplus.Py() << " " << Kplus.Pz() << " " << Kplus.E() << "\n";
    out << "K- " << Kminus.Px() << " " << Kminus.Py() << " " << Kminus.Pz() << " " << Kminus.E() << "\n";

    ++acc;
  }

  out.close();
  std::cout << "Generated " << acc << " / " << trials << " events accepted.\n";

  // ----------- save ROOT & PNGs -----------
  gStyle->SetOptStat(0);

  // Fit Mkk with Exp + BreitWigner
  TF1* fTot = new TF1("fTot", "[0]*exp(-[1]*x) + [2]*TMath::BreitWigner(x,[3],[4])", 0.99, 1.10);
  fTot->SetParNames("BkgNorm", "Lambda", "SigAmp", "Mu", "Gamma");
  fTot->SetParameters(500, 7.0, 1000, 1.0195, 0.0042);

  hMkk->Fit(fTot, "RQ");

  // save ROOT
  TFile fout("phi_gen_truth.root", "RECREATE");
  hW->Write();
  hQ2->Write();
  hWQ2->Write();
  hMkk->Write();
  fTot->Write();
  fout.Close();

  // save PNGs
  {
    TCanvas c;
    hW->Draw();
    c.SaveAs("W.png");
  }
  {
    TCanvas c;
    hQ2->Draw();
    c.SaveAs("Q2.png");
  }
  {
    TCanvas c;
    hWQ2->Draw("COLZ");
    c.SaveAs("WvsQ2.png");
  }
  {
    TCanvas c("cM", "Mkk", 900, 700);
    hMkk->Draw("E1");
    fTot->SetLineColor(kRed);
    fTot->Draw("SAME");
    c.SaveAs("Mkk_fit.png");
  }

  return 0;
}
