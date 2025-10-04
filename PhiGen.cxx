#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include "TCanvas.h"
#include "TFile.h"
#include "TGenPhaseSpace.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TLorentzVector.h"
#include "TRandom3.h"

// masses (GeV)
constexpr double MASS_E = 0.000511;
constexpr double MASS_P = 0.938272;
constexpr double MASS_PHI = 1.019461;
constexpr double WIDTH_PHI = 0.00426;
constexpr double MASS_K = 0.493677;

// physics helper functions
namespace physics {
double W_calc(const TLorentzVector &beam, const TLorentzVector &scattered) {
  TLorentzVector target(0, 0, 0, MASS_P);
  TLorentzVector q = beam - scattered;  // virtual photon
  TLorentzVector hadronic = q + target;
  return hadronic.M();
}
double Q2_calc(const TLorentzVector &beam, const TLorentzVector &scattered) {
  TLorentzVector q = beam - scattered;
  return -q.M2();  // Q^2 = - (q^2)
}
}  // namespace physics

// sample phi mass from Breit-Wigner
double breitWigner(double m0, double gamma) {
  double x, y;
  do {
    x = gRandom->Uniform(m0 - 0.05, m0 + 0.05);  // 50 MeV window
    y = gRandom->Uniform(0, 1);
  } while (y > (gamma / 2.0) / ((x - m0) * (x - m0) + (gamma * gamma / 4.0)));
  return x;
}

int main(int argc, char *argv[]) {
  int Nevents = 50000;  // number of events
  double Ebeam = 7.5;   // beam energy (GeV)
  double q2_min = 1.0, q2_max = 10.0;

  // initial states: e + p
  TLorentzVector target(0, 0, 0, MASS_P);
  TLorentzVector beam(0, 0, sqrt(Ebeam * Ebeam - MASS_E * MASS_E), Ebeam);
  TLorentzVector cms = beam + target;

  std::ofstream out("phi_gen.lund");

  // histograms
  TH1D *W_hist = new TH1D("W", "W distribution;W (GeV);Counts", 200, 1.0, 5.0);
  TH1D *Q2_hist = new TH1D("Q2", "Q^{2} distribution;Q^{2} (GeV^{2});Counts", 200, 0.0, 10.0);
  TH2D *WvsQ2 = new TH2D("WvsQ2", "W vs Q^{2};W (GeV);Q^{2} (GeV^{2})", 200, 1.0, 5.0, 200, 0.0, 10.0);

  TGenPhaseSpace event;
  int acc = 0;

  while (acc < Nevents) {
    // draw phi mass, ensure it can decay to KK
    double mphi;
    do {
      mphi = breitWigner(MASS_PHI, WIDTH_PHI);
    } while (mphi < 2 * MASS_K);

    // 3-body: e' + p' + phi(mphi)
    Double_t masses3[3] = {MASS_E, MASS_P, mphi};
    if (!event.SetDecay(cms, 3, masses3)) continue;
    if (!event.Generate()) continue;

    TLorentzVector *eprime = event.GetDecay(0);
    TLorentzVector *pprime = event.GetDecay(1);
    TLorentzVector phi = *event.GetDecay(2);

    // phi -> K+K-
    TGenPhaseSpace decay;
    Double_t kmasses[2] = {MASS_K, MASS_K};
    if (!decay.SetDecay(phi, 2, kmasses)) continue;
    if (!decay.Generate()) continue;

    TLorentzVector Kplus = *decay.GetDecay(0);
    TLorentzVector Kminus = *decay.GetDecay(1);

    // compute W and Q2
    double W = physics::W_calc(beam, *eprime);
    double Q2 = physics::Q2_calc(beam, *eprime);

    if (Q2 > q2_min && Q2 < q2_max) {
      W_hist->Fill(W);
      Q2_hist->Fill(Q2);
      WvsQ2->Fill(W, Q2);
    }

    // write in LUND-like format (simplified)
    out << "Event " << acc + 1 << "\n";
    out << "e' " << eprime->Px() << " " << eprime->Py() << " " << eprime->Pz() << " " << eprime->E() << "\n";
    out << "p' " << pprime->Px() << " " << pprime->Py() << " " << pprime->Pz() << " " << pprime->E() << "\n";
    out << "K+ " << Kplus.Px() << " " << Kplus.Py() << " " << Kplus.Pz() << " " << Kplus.E() << "\n";
    out << "K- " << Kminus.Px() << " " << Kminus.Py() << " " << Kminus.Pz() << " " << Kminus.E() << "\n";

    acc++;
  }

  out.close();
  std::cout << "Generated " << acc << " phi events\n";

  // save ROOT file
  TFile fout("phi_gen.root", "RECREATE");
  W_hist->Write();
  Q2_hist->Write();
  WvsQ2->Write();
  fout.Close();

  // save PNG plots
  TCanvas c1;
  W_hist->Draw();
  c1.SaveAs("W.png");

  TCanvas c2;
  Q2_hist->Draw();
  c2.SaveAs("Q2.png");

  TCanvas c3;
  WvsQ2->Draw("COLZ");
  c3.SaveAs("WvsQ2.png");

  return 0;
}
