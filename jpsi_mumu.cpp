#include "Pythia8/Pythia.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TFile.h"

using namespace Pythia8;

class JpsiDecay : public DecayHandler {
private:
  int times;               // contador de decaimentos do J/Psi
  ParticleData* pdtPtr;    // ponteiro para tabela de dados de partículas
  Rndm* rndmPtr;           // ponteiro para gerador de números aleatórios

public:
  JpsiDecay(ParticleData* pdtPtrIn, Rndm* rndmPtrIn) {
    times = 0;
    pdtPtr = pdtPtrIn; 
    rndmPtr = rndmPtrIn; 
  }

  bool decay(vector<int>& idProd, vector<double>& mProd,
             vector<Vec4>& pProd, int iDec, const Event& event);
};

bool JpsiDecay::decay(vector<int>& idProd, vector<double>& mProd,
                      vector<Vec4>& pProd, int iDec, const Event& event) {
  
  // decaimento: mumu
  idProd.push_back(-13);  // mu+ PDG (-13)
  idProd.push_back(13);   // mu- PDG (13)

  // massas dos múons
  double mMuon = pdtPtr->m0(13);
  mProd.push_back(mMuon);
  mProd.push_back(mMuon);

  double eMuon = 0.5 * mProd[0];  // energia de cada muon
  double pAbsMuon = sqrt(eMuon * eMuon - mMuon * mMuon);  // momento absoluto

  double cosTheta = 2. * rndmPtr->flat() - 1.; // cos(θ) [-1,1]
  double sinTheta = sqrt(max(0., 1. - cosTheta*cosTheta));
  double phi = 2. * M_PI * rndmPtr->flat(); // phi [0,2π]

  // Componentes do momento
  double pxMuon = pAbsMuon * sinTheta * cos(phi);
  double pyMuon = pAbsMuon * sinTheta * sin(phi);
  double pzMuon = pAbsMuon * cosTheta;

  // Quadrivetores 
  Vec4 pMuPlus(pxMuon, pyMuon, pzMuon, eMuon);
  Vec4 pMuMinus(-pxMuon, -pyMuon, -pzMuon, eMuon);
  
  pMuPlus.bst(pProd[0]); // aplicar a transformação de Lorentz
  pMuMinus.bst(pProd[0]); // usar o momento do J/Psi
  pProd.push_back(pMuPlus);
  pProd.push_back(pMuMinus);

  // mensagem de depuração (os primeiros 10 decaimentos)
  if (times++ < 10) {
    int iMother = event[iDec].mother1();
    int idMother = event[iMother].id();
    cout << "\n J/psi decay #" << times << " performed, J/psi in line " << iDec
         << ", mother id = " << idMother << "\n";
    cout << " Mu+ pT = " << pMuPlus.pT() << " GeV, eta = " << pMuPlus.eta() << "\n";
    cout << " Mu- pT = " << pMuMinus.pT() << " GeV, eta = " << pMuMinus.eta() << "\n";
  }

  return true;
}

//==========================================================================

int main() {
  // configura a simulação
  int nEvent = 10000;      // numero de eventos
  int nList = 2;           // numero de eventos para listar
  int nAbort = 5;          // numero máximo de erros

  // Inicialização do PYTHIA
  Pythia pythia;
  
  // configura a produção de J/Psi
  pythia.readString("Charmonium:all = on");
  pythia.readString("Beams:eCM = 7000.");  // energia em TeV
  pythia.readString("PhaseSpace:pTHatMin = 0.5");

  // configuração do decaimento
  DecayHandlerPtr handleDecays = make_shared<JpsiDecay>(&pythia.particleData, 
                                                       &pythia.rndm);
  vector<int> handledParticles;
  handledParticles.push_back(443);  // J/Psi PDG (443)
  pythia.setDecayPtr(handleDecays, handledParticles);

  // configuração de saida
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Inicia
  if (!pythia.init()) return 1;

  // cria os histogramas em ROOT
  TFile* outFile = new TFile("jpsi_mumu.root", "RECREATE");
  
  TH1D* pT_Jpsi = new TH1D("pT_Jpsi", "Transverse momentum of J/#psi; p_{T} (GeV); Counts", 
                             100, 0, 50);
  TH1D* eta_Jpsi = new TH1D("eta_Jpsi", "Pseudorapidity of J/#psi; #eta; Counts", 
                              100, -5, 5);
  TH1D* pT_mu = new TH1D("pT_mu", "Transverse momentum of muons; p_{T} (GeV); Counts", 
                           100, 0, 25);
  TH1D* eta_mu = new TH1D("eta_mu", "Pseudorapidity of muons; #eta; Counts", 
                            100, -5, 5);
  TH1D* mass_mumu = new TH1D("mass_mumu", "Invariant mass of #mu^{+}#mu^{-}; m (GeV); Counts", 
                                 100, 2.8, 3.4);

  // Loop dos eventos
  int iList = 0;
  int iAbort = 0;
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << "Event generation aborted after too many failures." << endl;
      break;
    }

    // verifica decaimentos externos
    bool externalDecay = false;
    for (int i = 0; i < pythia.event.size(); ++i) {
      if (pythia.event[i].statusAbs() == 93 || pythia.event[i].statusAbs() == 94) {
        externalDecay = true;
        break;
      }
    }

    // Lista alguns eventos
    if (externalDecay && ++iList <= nList) {
      pythia.process.list();
      pythia.event.list();
    }

    // preenche os histogramas
    for (int i = 0; i < pythia.event.size(); ++i) {
      if (pythia.event[i].id() == 443) {  // J/ψ
        pT_Jpsi->Fill(pythia.event[i].pT());
        eta_Jpsi->Fill(pythia.event[i].eta());
      }
      if (abs(pythia.event[i].id()) == 13) {  // mu+ ou mu-
        pT_mu->Fill(pythia.event[i].pT());
        eta_mu->Fill(pythia.event[i].eta());
      }
    }

    // calcula a massa invariante dos pares de muons
    for (int i = 0; i < pythia.event.size(); ++i) {
      if (pythia.event[i].id() == -13) {  // mu+
        for (int j = i+1; j < pythia.event.size(); ++j) {
          if (pythia.event[j].id() == 13) {  // mu-
            Vec4 pMuPlus = pythia.event[i].p();
            Vec4 pMuMinus = pythia.event[j].p();
            double mass = (pMuPlus + pMuMinus).mCalc();
            mass_mumu->Fill(mass);
          }
        }
      }
    }
  }

  // Finaliza
  pythia.stat();
  
  // salva os histogramas
  outFile->Write();
  outFile->Close();

  // Plota os histogramas
  TCanvas* c1 = new TCanvas("c1", "J/psi mumu", 1200, 800);
  c1->Divide(3,2);
  
  c1->cd(1); pT_Jpsi->Draw();
  c1->cd(2); eta_Jpsi->Draw();
  c1->cd(3); pT_mu->Draw();
  c1->cd(4); eta_mu->Draw();
  c1->cd(5); mass_mumu->Draw();
  
  c1->SaveAs("jpsi_mumu.pdf");

  return 0;
}
