// Desenvolvido por: Bruno Kron Guandalini
// Descrição: Simulação de produção de W⁺/W⁻ com decaimento em μ⁺ν_μ / μ⁻ν̄_μ usando Pythia8 e ROOT.
// Objetivo: Estudar as propriedades do W⁺/W⁻ e seu decaimento em μ⁺ν_μ / μ⁻ν̄_μ a 13 TeV.

#include "Pythia8/Pythia.h" // Biblioteca principal do Pythia8
#include "TH1D.h"           // Histograma unidimensional do ROOT
#include "TCanvas.h"        // Canvas para desenho do histograma
#include "TFile.h"          // Arquivo para salvar os dados

using namespace Pythia8;

int main() {
  // ================================================================
  // 1. Configuração do Pythia e do processo de simulação
  // ================================================================
  Pythia pythia;

  // Configurações básicas 
  pythia.readString("Beams:idA = 2212"); // Declara que é o feixe A é próton (PDG 2212)  
  pythia.readString("Beams:idB = 2212"); // Declara que é o feixe B é próton (PDG 2212)
  pythia.readString("Beams:eCM = 13000."); // Energia total no centro de massa (13 TeV)

  // Produção de W⁺/W⁻ (e Z⁰) via Drell-Yan
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");

  // Força o decaimento leptônico do W → μ ν
  pythia.readString("24:onMode = off");        // Desliga todos os canais do W⁺ (PDG 24)
  pythia.readString("24:onIfAny = 13 14");     // ativa W⁺ → μ⁺ ν_μ (PDGs 13 e 14)
  pythia.readString("-24:onMode = off");       // Desliga todos os canais do W⁻ (PDG -24)
  pythia.readString("-24:onIfAny = -13 -14");  // ativa W⁻ → μ⁻ ν̄_μ (PDGs -13 e -14)

  // Configuração da saída do Pythia (para evitar que o código fique muito difícil de ler e compreender por conta de muitas mensagens)
  pythia.readString("Next:numberShowInfo = 0"); // Desativa a exibição de informações gerais sobre o processo
  pythia.readString("Next:numberShowProcess = 0"); // Desativa a exibição de detalhes do processo de geração
  pythia.readString("Next:numberShowEvent = 0"); // Desativa a exibição de detalhes de cada evento gerado

  // Inicializa o Pythia
  if (!pythia.init()) {
    cerr << "Erro na inicialização do Pythia." << endl;
    return 1;
  }

  // ================================================================
  // 2. Configuração dos histogramas do ROOT para análise dos eventos
  // ================================================================
 
  // Cria um arquivo ROOT para salvar os histogramas
  TFile outFile("w_munu.root", "RECREATE");

  // Massa Transversa do μ + ν
  TH1D h_mass_W("h_mass_W", "Transverse Mass of #mu + #nu; m_{#mu#nu} (GeV/c^{2}); Events", 100, 0., 150.);

  // Momento Transverso do W 
  TH1D h_pt_W("h_pt_W", "Transverse Momentum of W^{#pm}; p_{T} (GeV/c); Events", 100, 0., 200.);

  // Pseudorapidez do W
  TH1D h_eta_W("h_eta_W", "Pseudorapidity of W^{#pm}; #eta; Events", 100, -5., 5.);

  // Momento transverso do múon
  TH1D h_pt_mu("h_pt_mu", "Transverse Momentum of #mu; p_{T} (GeV/c); Events", 100, 0., 100.);

  // Momento transverso do neutrino (Energia faltante transversal)
  TH1D h_pt_nu("h_pt_nu", "Transverse Momentum of #nu; p_{T}^{#nu} (GeV/c); Events", 100, 0., 100.);

  // Pseudorapidez do múon
  TH1D h_eta_mu("h_eta_mu", "Pseudorapidity of #mu; #eta; Events", 100, -5., 5.);

  // Pseudorapidez do neutrino
  TH1D h_eta_nu("h_eta_nu", "Pseudorapidity of #nu; #eta; Events", 100, -5., 5.);

  // ================================================================
  // 3. Loop de eventos
  // ================================================================
  const int nEvent = 10000;  // Número de eventos a serem gerados (altere conforme necessário)
  int nAbort = 5; // Número máximo de falhas consecutivas antes de abortar a geração
  int iAbort = 0; // Contador de falhas consecutivas

  // Loop principal de geração de eventos
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cerr << "Geração abortada após muitas falhas." << endl;
      break;
    }

    // Procura W⁺ ou W⁻ no evento
    for (int i = 0; i < pythia.event.size(); ++i) {
      int pid = pythia.event[i].id(); // Identificador PDG da partícula

      // Identifica as filhas do W⁺ (PDG 24) ou W⁻ (PDG -24)
      if (abs(pid) == 24) {
        int d1 = pythia.event[i].daughter1(); // Índice da primeira filha
        int d2 = pythia.event[i].daughter2(); // Índice da segunda filha
 
        // Garante que as filhas existam
        if (d1 <= 0 || d2 <= 0) continue;

        int id1 = pythia.event[d1].id();
        int id2 = pythia.event[d2].id();

        // Identifica qual filha é o muon (id = 13) e qual é o neutrino (id = 14)
        int iMuon = -1, iNu = -1;
        if (abs(id1) == 13 && abs(id2) == 14) {
          iMuon = d1; iNu = d2;
        } else if (abs(id1) == 14 && abs(id2) == 13) {
          iMuon = d2; iNu = d1;
        }
        if (iMuon < 0 || iNu < 0) continue; // não encontrou o padrão esperado (μ + ν)

        // Quadrivetores
        Vec4 pMu = pythia.event[iMuon].p(); // Quadrivetor do muon 
        Vec4 pNu = pythia.event[iNu].p(); // Quadrivetor do neutrino

        // ---- Preenche histogramas básicos ----
        double pt_mu = pMu.pT();
        double eta_mu = pMu.eta();
        double pt_nu = pNu.pT();
        double eta_nu = pNu.eta();
        double dphi = abs(pMu.phi() - pNu.phi());
        if (dphi > M_PI) dphi = 2*M_PI - dphi;  // garante Δφ está entre [0, π]

        h_pt_mu.Fill(pt_mu); // pT do múon
        h_eta_mu.Fill(eta_mu); // η do múon
        h_pt_nu.Fill(pt_nu); // pT do neutrino (Energia faltante transversal)
        h_eta_nu.Fill(eta_nu); // η do neutrino

        // ---- Massa transversa ----
        double mass_W = sqrt(2. * pt_mu * pt_nu * (1. - cos(dphi))); // mass_W = sqrt(2 * pT_μ * pT_ν * (1 - cos(Δφ)))
        h_mass_W.Fill(mass_W); // Massa transversa do W

        // ---- W reconstruído (μ + ν) ----
        Vec4 pW = pMu + pNu; // Quadrivetor do W reconstruído a partir do muon e do neutrino
        h_pt_W.Fill(pW.pT()); // pT do W reconstruído
        h_eta_W.Fill(pW.eta()); // η do W reconstruído
      }
    }
  }

  // ================================================================
  // 4. Finalização e geração de PDF
  // ================================================================
  pythia.stat();
  outFile.Write();
  outFile.Close();

  TCanvas c1("c1", "W^{#pm} #rightarrow #mu + #nu_{#mu}", 1200, 800);
  c1.Divide(3, 3);

  // Cores
  h_mass_W.SetLineColor(kBlue);
  h_pt_W.SetLineColor(kBlue);
  h_eta_W.SetLineColor(kBlue);
  h_pt_mu.SetLineColor(kRed);
  h_pt_nu.SetLineColor(kGreen);
  h_eta_mu.SetLineColor(kRed);
  h_eta_nu.SetLineColor(kGreen);

  // ---- Desenha os histogramas ----

  // Massa transversa do W
  c1.cd(1); 
  h_mass_W.Draw(); 

  // Momento transverso (pT) do W (escala log)
  c1.cd(2); 
  gPad->SetLogy(1); 
  h_pt_W.Draw();    // pT do W (log)

  // Pseudorapidez (η) do W
  c1.cd(3); 
  h_eta_W.Draw();      

  // Momento transverso (pT) do muon (escala log)
  c1.cd(4); 
  gPad->SetLogy(1); 
  h_pt_mu.Draw();  

  // Momento transverso (pT) do neutrino (escala log)
  c1.cd(5);  
  gPad->SetLogy(1); 
  h_pt_nu.Draw();     
  
  // Pseudorapidez (η) do muon
  c1.cd(6);
  h_eta_mu.Draw();      

  // Pseudorapidez (η) do neutrino
  c1.cd(7); 
  h_eta_nu.Draw();      

  c1.SaveAs("w_munu.pdf");

  cout << "\nArquivo ROOT: w_munu.root" << endl;
  cout << "PDF gerado: w_munu.pdf" << endl;

  return 0;
}
