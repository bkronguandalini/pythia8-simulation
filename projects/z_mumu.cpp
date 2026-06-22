// Desenvolvido por: Bruno Kron Guandalini
// Descrição: Simulação de produção de Z⁰ com decaimento em μ⁺μ⁻ usando Pythia8 e ROOT.
// Objetivo: Estudar a massa invariante do Z⁰ a 13 TeV.

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

  // Processo: produção de Z⁰ via aniquilação de quark-antiquark
  pythia.readString("WeakSingleBoson:ffbar2gmZ = on");

  // Força o decaimento do Z⁰ apenas em μ⁺μ⁻
  pythia.readString("23:onMode = off");       // desliga todos os canais de decaimento do Z⁰
  pythia.readString("23:onIfAny = 13");       // ativa apenas canal com muon (PDG 13)

  // Configuração da saída do Pythia (para evitar que o código fique muito difícil de ler e compreender por conta de muitas mensagens)
  pythia.readString("Next:numberShowInfo = 0"); // Desativa a exibição de informações gerais sobre o processo
  pythia.readString("Next:numberShowProcess = 0"); // Desativa a exibição de detalhes do processo de geração
  pythia.readString("Next:numberShowEvent = 0"); // Desativa a exibição de detalhes de cada evento gerado

  // Inicializa o gerador
  if (!pythia.init()) {
    cerr << "Erro na inicialização do Pythia." << endl;
    return 1;
  }

  // ================================================================
  // 2. Configuração dos histogramas e dos arquivos de saída
  // ================================================================

  TFile outFile("z_mumu.root", "RECREATE");

  // Histograma da massa invariante dos pares μ⁺μ⁻
  TH1D h_mass_mumu("h_mumu", "Invariant Mass of #mu^{+}#mu^{-}; m_{#mu^{+}#mu^{-}} (GeV/c^{2}); Eventos", 100, 60., 120.);

  // Histograma do momento transverso do Z⁰
  TH1D h_pt_Z("h_pt_Z", "Transverse Momentum of Z^{0}; p_{T} (GeV/c); Eventos", 100, 0., 200.);

  // Histograma do pseudorapidez do Z⁰
  TH1D h_eta_Z("h_eta_Z", "Pseudorapidity of Z^{0}; #eta; Eventos", 100, -5., 5.);

  // Histograma do momento transverso dos muons
  TH1D h_pt_muPlus("h_pt_muPlus", "Transverse Momentum of #mu^{+}; p_{T} (GeV/c); Eventos", 100, 0., 100.);
  TH1D h_pt_muMinus("h_pt_muMinus", "Transverse Momentum of #mu^{-}; p_{T} (GeV/c); Eventos", 100, 0., 100.);

  // Histograma do pseudorapidez dos muons
  TH1D h_eta_muPlus("h_eta_muPlus", "Pseudorapidity of #mu^{+}; #eta; Eventos", 100, -5., 5.);
  TH1D h_eta_muMinus("h_eta_muMinus", "Pseudorapidity of #mu^{-}; #eta; Eventos", 100, -5., 5.);

  // ================================================================
  // 3. Loop de geração de eventos e preenchimento dos histogramas
  // ================================================================

  const int nEvent = 10000; // Número de eventos a serem gerados (altere conforme necessário)
  int nAbort = 5;  // Número máximo de falhas permitidas antes de abortar a geração
  int iAbort = 0; // Contador de falhas

  // Loop de geração de eventos
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cerr << "Geração abortada após muitas falhas." << endl;
      break;
    }

  // Preenche os histogramas com os dados do evento
  for (int i = 0; i < pythia.event.size(); ++i) {
      int pid = pythia.event[i].id();

  // Preenche os histogramas de pT e η para Z° e μ⁺, μ⁻
    if (pid == 23) {
        h_pt_Z.Fill(pythia.event[i].pT());
        h_eta_Z.Fill(pythia.event[i].eta());
      }
      // Preenche os histogramas de pT e η para μ⁺
      else if (pid == -13) { 
        h_pt_muPlus.Fill(pythia.event[i].pT());
        h_eta_muPlus.Fill(pythia.event[i].eta());
      }
      // Preenche os histogramas de pT e η para μ⁻
      else if (pid == 13) {   
        h_pt_muMinus.Fill(pythia.event[i].pT());
        h_eta_muMinus.Fill(pythia.event[i].eta());
      }
    }

    // Procura o Z⁰ (PDG = 23) que decaiu e suas filhas (μ⁺ e μ⁻)
    for (int i = 0; i < pythia.event.size(); ++i) {
      if (pythia.event[i].id() == 23) {  // Z⁰ encontrado
        int d1 = pythia.event[i].daughter1(); // Índice da primeira filha
        int d2 = pythia.event[i].daughter2(); // Índice da segunda filha

        if (d1 > 0 && d2 > 0 &&
            abs(pythia.event[d1].id()) == 13 &&
            abs(pythia.event[d2].id()) == 13) {
          Vec4 pSum = pythia.event[d1].p() + pythia.event[d2].p();
          double massa = pSum.mCalc(); // Calcula (soma) a massa do par μ⁺μ⁻
          h_mass_mumu.Fill(massa); // Preenche o histograma da massa invariante
        }
      }
    }
  }

  // ---------- Finalização ----------
  pythia.stat();   // Exibe estatísticas da simulação
  outFile.Write(); // Escreve os histogramas no arquivo ROOT
  outFile.Close(); // Fecha o arquivo

  // ==================================================================
  // 4. Desenha os histogramas usando o ROOT e salva em um arquivo PDF
  // ==================================================================

  TCanvas c1("c1", "Z^{0} #rightarrow #mu^{+}#mu^{-}", 1200, 800);
  c1.Divide(3, 3);

  // Cores para diferenciar os histogramas de cada partícula (constantes do ROOT) 
  h_mass_mumu.SetLineColor(kBlue); 
  h_pt_Z.SetLineColor(kBlue);
  h_eta_Z.SetLineColor(kBlue);
  h_pt_muPlus.SetLineColor(kRed);
  h_pt_muMinus.SetLineColor(kGreen+2);
  h_eta_muPlus.SetLineColor(kRed);
  h_eta_muMinus.SetLineColor(kGreen+2);

  // Desenha os histogramas em cada subcanvas

  // Massa invariante
  c1.cd(1); h_mass_mumu.Draw();

  // Momento Transverso (pT) do Z° (escala log)
  c1.cd(2); gPad->SetLogy(1);
  h_pt_Z.Draw();

  // Pseudorapidez (η) do Z°
  c1.cd(3); h_eta_Z.Draw();

  // Momento Transverso (pT) do μ⁺ (escala log)
  c1.cd(4); gPad->SetLogy(1);
  h_pt_muPlus.Draw();

  // Momento Transverso (pT) do μ⁻ (escala log)
  c1.cd(5); gPad->SetLogy(1);
  h_pt_muMinus.Draw();

  // Pseudorapidez (η) do μ⁺ 
  c1.cd(6); h_eta_muPlus.Draw();

  // Pseudorapidez (η) do μ⁻
  c1.cd(7); h_eta_muMinus.Draw();

  // Salva o canvas com os histogramas em um arquivo PDF
  c1.SaveAs("z_mumu.pdf");

  cout << "\nArquivo ROOT: z_mumu.root" << endl;
  cout << "PDF gerado: z_mumu.pdf" << endl;

  return 0;
}
