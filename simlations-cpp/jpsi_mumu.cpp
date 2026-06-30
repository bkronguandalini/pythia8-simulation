// Desenvolvido por: Bruno Kron Guandalini
// Descrição: Simulação de produção de J/ψ com decaimento em μ⁺μ⁻ usando Pythia8 e ROOT.
// Objetivo: Estudar as propriedades do J/ψ e seu decaimento em μ⁺μ⁻ a 13 TeV.

#include "Pythia8/Pythia.h" // Biblioteca principal do Pythia8
#include "TH1D.h" // Histograma unidimensional do ROOT
#include "TCanvas.h" // Canvas para desenho do histograma
#include "TFile.h" // Arquivo para salvar os dados

using namespace Pythia8;

int main() {
  // ================================================================
  // 1. Configuração do Pythia e do processo de simulação
  // ================================================================
  Pythia pythia;

  // Configurações iniciais
  pythia.readString("Beams:idA = 2212");      // Definição do próton (PDG 2212)
  pythia.readString("Beams:idB = 2212");      // Definição do próton (PDG 2212)
  pythia.readString("Beams:eCM = 13000.");    // Energia do centro de massa (13 TeV)

  // Processo: produção de J/ψ via fusão glúon-glúon
  pythia.readString("Charmonium:gg2ccbar(3S1)[3S1(1)]g = on");
  pythia.readString("PhaseSpace:pTHatMin = 1.0"); // pT mínimo do processo duro

  // Força o decaimento do J/ψ apenas em μ⁺μ⁻
  pythia.readString("443:onMode = off");      // desliga todos os canais de decaimento do J/ψ
  pythia.readString("443:onIfAny = 13");      // ativa apenas canal com muon (PDG 13)

  // Configuração da saída do Pythia (reduz mensagens na tela)
  pythia.readString("Next:numberShowInfo = 0");
  pythia.readString("Next:numberShowProcess = 0");
  pythia.readString("Next:numberShowEvent = 0");

  // Inicializa o gerador
  if (!pythia.init()) {
    cerr << "Erro na inicialização do Pythia." << endl;
    return 1;
  }

  // ================================================================
  // 2. Configuração dos histogramas e dos arquivos de saída
  // ================================================================

  TFile outFile("jpsi_mumu.root", "RECREATE");

  // Histograma da massa invariante dos pares μ⁺μ⁻
  TH1D h_mass_mumu("h_mass_mumu", "Invariant Mass of #mu^{+}#mu^{-}; m_{#mu^{+}#mu^{-}} (GeV/c^{2}); Eventos", 80, 3.0, 3.2);

  // Histograma do momento transverso do J/ψ
  TH1D h_pt_Jpsi("h_pt_Jpsi", "Transverse Momentum of J/#psi; p_{T} (GeV/c); Eventos", 100, 0., 50.);

  // Histograma do pseudorapidez do J/ψ
  TH1D h_eta_Jpsi("h_eta_Jpsi", "Pseudorapidity of J/#psi; #eta; Eventos", 100, -5., 5.);

  // Histograma do momento transverso dos muons
  TH1D h_pt_muPlus("h_pt_muPlus", "Transverse Momentum of #mu^{+}; p_{T} (GeV/c); Eventos", 100, 0., 25.);
  TH1D h_pt_muMinus("h_pt_muMinus", "Transverse Momentum of #mu^{-}; p_{T} (GeV/c); Eventos", 100, 0., 25.);

  // Histograma do pseudorapidez dos muons
  TH1D h_eta_muPlus("h_eta_muPlus", "Pseudorapidity of #mu^{+}; #eta; Eventos", 100, -5., 5.);
  TH1D h_eta_muMinus("h_eta_muMinus", "Pseudorapidity of #mu^{-}; #eta; Eventos", 100, -5., 5.);

  // ================================================================
  // 3. Loop de geração de eventos e preenchimento dos histogramas
  // ================================================================

  const int nEvent = 10000; // Número de eventos a serem gerados (altere conforme necessário)
  int nAbort = 5;  // Número máximo de falhas permitidas antes de abortar a geração
  int iAbort = 0;  // Contador de falhas

  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cerr << "Geração abortada após muitas falhas." << endl;
      break;
    }
    iAbort = 0; // reseta o contador a cada evento bem-sucedido

    // Preenche os histogramas com os dados do evento
    for (int i = 0; i < pythia.event.size(); ++i) {
      int pid = pythia.event[i].id();

      // Preenche os histogramas de pT e η para J/ψ e μ⁺, μ⁻
      if (pid == 443 && !pythia.event[i].isFinal()) {
        h_pt_Jpsi.Fill(pythia.event[i].pT());
        h_eta_Jpsi.Fill(pythia.event[i].eta());
      }
      // Preenche os histogramas de pT e η para μ⁺
      else if (pid == -13 && pythia.event[i].isFinal()) {
        h_pt_muPlus.Fill(pythia.event[i].pT());
        h_eta_muPlus.Fill(pythia.event[i].eta());
      }
      // Preenche os histogramas de pT e η para μ⁻
      else if (pid == 13 && pythia.event[i].isFinal()) {
        h_pt_muMinus.Fill(pythia.event[i].pT());
        h_eta_muMinus.Fill(pythia.event[i].eta());
      }
    }

    // Procura o J/ψ (PDG = 443) que decaiu e suas filhas (μ⁺ e μ⁻)
    for (int i = 0; i < pythia.event.size(); ++i) {
      if (pythia.event[i].id() == 443) {  // J/ψ encontrado
        int d1 = pythia.event[i].daughter1(); // Índice da primeira filha
        int d2 = pythia.event[i].daughter2(); // Índice da segunda filha

        if (d1 > 0 && d2 > 0 &&
            abs(pythia.event[d1].id()) == 13 &&
            abs(pythia.event[d2].id()) == 13) {
          Vec4 pSum = pythia.event[d1].p() + pythia.event[d2].p();
          double massa = pSum.mCalc(); // Calcula a massa do par μ⁺μ⁻
          h_mass_mumu.Fill(massa);     // Preenche o histograma da massa invariante
        }
      }
    }
  }

  // ---------- Finalização ----------
  pythia.stat();   // Exibe estatísticas da simulação
  outFile.Write(); // Escreve os histogramas no arquivo ROOT
  outFile.Close(); // Fecha o arquivo

  // ==================================================================
  // 4. PDF multipágina
  // ==================================================================

  // Cores
  h_mass_mumu.SetLineColor(kBlue);
  h_pt_Jpsi.SetLineColor(kBlue);
  h_eta_Jpsi.SetLineColor(kBlue);
  h_pt_muPlus.SetLineColor(kRed);
  h_pt_muMinus.SetLineColor(kGreen+2);
  h_eta_muPlus.SetLineColor(kRed);
  h_eta_muMinus.SetLineColor(kGreen+2);

  // Estrutura auxiliar para agrupar histograma + flag de escala log
  struct Plot {
      TH1D* hist;
      bool logy;
  };

  Plot h_logs[] = {
      {&h_mass_mumu, false},
      {&h_pt_Jpsi,   true},
      {&h_eta_Jpsi,  false},
      {&h_pt_muPlus, true},
      {&h_pt_muMinus,true},
      {&h_eta_muPlus,false},
      {&h_eta_muMinus,false}
  };
  int nPlots = sizeof(h_logs) / sizeof(h_logs[0]);

  // Canvas único e PDF
  TCanvas c("c", "J/#psi #rightarrow #mu^{+}#mu^{-}", 1200, 800);
  TString pdf_jpsi = "jpsi_mumu.pdf";
  c.Print(pdf_jpsi + "[");   // inicia o PDF multipágina

  for (int i = 0; i < nPlots; ++i) {
      c.Clear();
      c.SetLogy(h_logs[i].logy ? 1 : 0);
      h_logs[i].hist->Draw("HIST");
      c.Update();
      c.Print(pdf_jpsi);
  }

  c.Print(pdf_jpsi + "]");   // finaliza o PDF

  cout << "\nArquivo ROOT: jpsi_mumu.root" << endl;
  cout << "PDF gerado: jpsi_mumu.pdf" << endl;

  return 0;
}
