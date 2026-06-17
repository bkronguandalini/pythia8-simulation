#include "Pythia8/Pythia.h"
#include "TH1F.h"
#include "TFile.h"
#include "TLorentzVector.h"

using namespace Pythia8;

int main() {
  Pythia pythia;
  pythia.readString("Beams:idA = 2212"); //próton
  pythia.readString("Beams:idB = 2212"); //próton
  pythia.readString("Beams:eCM = 13000."); //energia de centro de massa (13TeV)
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on"); //ativa produção de Z
  pythia.readString("23:onMode = off"); // Z decays
  pythia.readString("23:onIfAny = 13"); // apenas Z → μ⁺μ⁻ (13 é o codigo do PDG)
  pythia.init(); 

  // ROOT file and histogram
  TFile *outfile = new TFile("z_mumu.root", "RECREATE");
  TH1F *h_mll = new TH1F("h_mll", "Mass of mu+mu-;m_{#mu#mu} [GeV];Events", 100, 60., 120.); //Histograma da massa invariante dos múons

  for (int iEvent = 0; iEvent < 10000; ++iEvent) {
    if (!pythia.next()) continue;

    std::vector<TLorentzVector> muons; //coleta múons do evento
    for (int i = 0; i < pythia.event.size(); ++i) {
      if (pythia.event[i].id() == 13 || pythia.event[i].id() == -13) { //seleciona múons (13) e antimúons (-13)
        TLorentzVector mu(pythia.event[i].px(), pythia.event[i].py(),
                          pythia.event[i].pz(), pythia.event[i].e());
        muons.push_back(mu); //armazena o 4-vetor do múon
      }
    }

	// encontra 2 múons e calcula a massa invariante		
    if (muons.size() == 2) {
      TLorentzVector dimuon = muons[0] + muons[1]; // soma os 4-vetores
      h_mll->Fill(dimuon.M()); // preenche o histograma com a massa invariante
    }
  }

  pythia.stat(); // imprime estatísticas da simulção
  h_mll->Write();
  outfile->Close(); 

  return 0;
}
