// Desenvolvido por: Bruno Kron Guandalini
// Descrição: Simulação de produção de J/ψ com decaimento em μ⁺μ⁻ usando o Pythia8 e ROOT.
// Objetivo: Estudar as propriedades do J/ψ e seu decaimento em μ⁺μ⁻ a 13 TeV.

#include "Pythia8/Pythia.h" // inclui a classe principal do Pythia
#include "TH1D.h" // histogramas 1D do ROOT
#include "TCanvas.h" // para desenhar os histogramas
#include "TFile.h" // para salvar os histogramas em um arquivo ROOT

using namespace Pythia8;

// ================================================================
// Classe para manipular o decaimento do J/ψ
// ================================================================
class JpsiDecay : public DecayHandler {
private:
  int nJpsiTotal;          // contador de decaimentos do J/ψ
  ParticleData* pdtPtr;    // acesso à tabela de partículas
  Rndm* rndmPtr;           // gerador de números aleatórios do Pythia

// ================================================================
// Implementação do decaimento J/ψ → μ⁺ μ⁻
// ================================================================
public:
  // Construtor: guarda os ponteiros necessários
  JpsiDecay(ParticleData* pdtPtrIn, Rndm* rndmPtrIn) {
    nJpsiTotal = 0; // contador de decaimentos do J/ψ, começando em zero
    pdtPtr = pdtPtrIn; // ponteiro para a tabela de partículas
    rndmPtr = rndmPtrIn; // ponteiro para o gerador de números aleatórios
  }

  // Método chamado pelo Pythia sempre que um J/ψ precisa decair
  bool decay(vector<int>& idProd, vector<double>& mProd,
             vector<Vec4>& pProd, int iDec, const Event& evento);
};

// ----------------------------------------------------------------
bool JpsiDecay::decay(vector<int>& idProd, vector<double>& mProd,
                      vector<Vec4>& pProd, int iDec, const Event& evento) {

// ================================================================
// Implementação do decaimento J/ψ → μ⁺ μ⁻
// ================================================================
  // 1. Definir as partículas filhas μ⁺ e μ⁻ (PDG IDs -13 e 13)
  idProd.push_back(-13);   // muon mais (μ⁺)
  idProd.push_back(13);   // muon menos (μ⁻)

  // 2. Massa dos muons, obtida da tabela de partículas (as duas têm a mesma massa)
  double mMuon = pdtPtr->m0(13); // massa do muon (PDG ID 13)
  mProd.push_back(mMuon); // massa do μ⁺
  mProd.push_back(mMuon); // massa do μ⁻

  // 3. Massa do J/ψ (obtida do evento)
  double mJpsi = evento[iDec].m();

  // 4. Cinemática no referencial de repouso do J/ψ
  double eMuon = 0.5 * mJpsi;                     // energia de cada muon
  double pMuon  = sqrt(eMuon*eMuon - mMuon*mMuon); // módulo do momento do muon

  // 5. Ângulos isotrópicos para a direção do momento do μ⁺ 
  double cosTheta = 2. * rndmPtr->flat() - 1.; // cos(θ) uniformemente distribuído em [-1, 1]
  double sinTheta = sqrt(max(0., 1. - cosTheta*cosTheta)); // sin(θ) correspondente
  double phi      = 2. * M_PI * rndmPtr->flat(); // φ uniformemente distribuído em [0, 2π]

  // 6. Construção dos quadrivetores (px, py, pz, E) para μ⁺ e μ⁻
  double px = pMuon * sinTheta * cos(phi); // componente x do momento
  double py = pMuon * sinTheta * sin(phi); // componente y do momento
  double pz = pMuon * cosTheta; // componente z do momento

  // Criação dos quadrivetores para μ⁺ e μ⁻
  Vec4 pMuPlus ( px,  py,  pz, eMuon);  // momento do μ⁺
  Vec4 pMuMinus(-px, -py, -pz, eMuon);  // momento do μ⁻ (sentido oposto)

  // Adiciona os quadrivetores à lista de produções
  pProd.push_back(pMuPlus); // adiciona o μ⁺ à lista de produções
  pProd.push_back(pMuMinus); // adiciona o μ⁻ à lista de produções

  // 7. Mostra os primeiros 10 decaimentos
  if (nJpsiTotal++ < 10) {
    int iJpsi   = evento[iDec].mother1(); // índice J/ψ (particula mãe)
    int idJpsi  = (iJpsi > 0) ? evento[iJpsi].id() : 0; // PDG ID J/ψ (particula mãe)
    cout << "\n Decaimento J/psi #" << nJpsiTotal 
         << " na linha " << iDec
         << ", mae id = " << idJpsi << "\n";

    // Mostra os quadrivetores no referencial do laboratório
    Vec4 pJpsi = evento[iDec].p(); // quadrivetor da partícula mãe (J/ψ)
    Vec4 muPlusLab  = pMuPlus;   muPlusLab.bst(pJpsi); // referencial do laboratório
    Vec4 muMinusLab = pMuMinus;  muMinusLab.bst(pJpsi); // referencial do laboratório
    cout << " μ⁺ pT = " << muPlusLab.pT()
         << " GeV, η = " << muPlusLab.eta() << "\n";
    cout << " μ⁻ pT = " << muMinusLab.pT()
         << " GeV, η = " << muMinusLab.eta() << "\n";
  }

  return true;  // decaimento bem-sucedido
}

// ================================================================
// Programa principal
// ================================================================
int main() {
  // ---------- Parâmetros da simulação ----------
  const int nEvent = 10000;   // número de eventos (altere conforme necessário)
  const int nAbort = 5;       // falhas consecutivas permitidas

  // Inicialização do Pythia (responsável por gerenciar toda a simulação)
  Pythia pythia;
  
  /* CONFIGURAÇÃO DA FÍSICA: 
  A lógica aqui é preparar o Pythia para gerar especificamente o processo de produção de J/ψ seguido pelo decaimento em μ⁺μ⁻ */
  pythia.readString("Charmonium:all = on");     // Ativa a produção de charmonium
  pythia.readString("Beams:eCM = 13000.");        // Energia de 13 TeV no centro de massa (colisão padrão: pp)
  pythia.readString("PhaseSpace:pTHatMin = 0.5"); // Este parâmetro atua apenas na geração do processo duro (pT mínimo), não substituindo cortes cinemáticos aplicados na análise.

  // Instala o decaimento customizado J/ψ → μ⁺ μ⁻
  DecayHandlerPtr handleDecays =
      make_shared<JpsiDecay>(&pythia.particleData, &pythia.rndm);
  vector<int> particulasManipuladas;
  particulasManipuladas.push_back(443);  // PDG ID do J/ψ
  pythia.setDecayPtr(handleDecays, particulasManipuladas);

    // Configuração da saída do Pythia (para evitar que o código fique muito difícil de ler e compreender por conta de muitas mensagens)
  pythia.readString("Next:numberShowInfo = 0"); // Desativa a exibição de informações gerais sobre o processo
  pythia.readString("Next:numberShowProcess = 0"); // Desativa a exibição de detalhes do processo de geração
  pythia.readString("Next:numberShowEvent = 0"); // Desativa a exibição de detalhes de cada evento gerado

  // Inicia o Pythia e verifica se houve erro na inicialização
  if (!pythia.init()) {
    cout << "Erro na inicialização do Pythia!" << endl;
    return 1;
  }

  cout << "Pythia inicializado com sucesso!" << endl;

  /* LÓGICA DOS HISTOGRAMAS:
    Os histogramas servem para estudar os diferentes aspectos físicos do J/ψ e seu decaimento em μ⁺μ⁻. 
    Cada histograma apresenta:
     - Nome interno para o ROOT (sem espaços)
     - Título para display (com formatação LaTeX com pequenas alterações) 
     (substítuir \ por # e sempre que adicionar ^ colocar { } em seguida, isso no ROOT)
     - Rótulos dos eixos (x e y)
     - Número de bins e range */

  // Cria um arquivo ROOT para armazenar os histogramas
  TFile outFile("jpsi_mumu.root", "RECREATE");

  // Histograma da massa invariante μ⁺ μ⁻
  TH1D h_mass_mumu("mass_mumu","Invariant Mass of #mu^{+}#mu^{-}; m_{#mu^{+}#mu^{-}} (GeV/c^{2}); Events", 80, 3.0, 3.2);

  // Histograma do momento transverso do J/ψ
  TH1D h_pt_Jpsi("pt_Jpsi","Transverse Momentum of J/#psi; p_{T} (GeV); Events", 100, 0, 50);
  
  // Histograma da pseudorapidez do J/ψ
  TH1D h_eta_Jpsi("eta_Jpsi","Pseudorapidity of J/#psi; #eta; Events", 100, -5, 5);
  
  // Histograma do momento transverso dos muons
  TH1D h_pt_muPlus("pt_muPlus","Transverse Momentum of #mu^{+}; p_{T} (GeV); Events", 100, 0, 25);
  TH1D h_pt_muMinus("pt_muMinus","Transverse Momentum of #mu^{-}; p_{T} (GeV); Events", 100, 0, 25);
  
  // Histograma da pseudorapidez dos muons
  TH1D h_eta_muPlus("eta_muPlus","Pseudorapidity of #mu^{+}; #eta; Events", 100, -5, 5);
  TH1D h_eta_muMinus("eta_muMinus","Pseudorapidity of #mu^{-}; #eta; Events", 100, -5, 5);
  
  // ================================================================
  // Loop de eventos
  // ================================================================
  int iAbort = 0; // contador de falhas consecutivas na geração de eventos
    for (int iEvent = 0; iEvent < nEvent; ++iEvent) {
      if (!pythia.next()) {
        if (++iAbort < nAbort) continue;
        cerr << "Geracao abortada apos muitas falhas." << endl;
        break;
      }

    // Preenche os histogramas com os dados do evento
    for (int i = 0; i < pythia.event.size(); ++i) {
      int pid = pythia.event[i].id();

      // Preenche os histogramas de pT e η para J/ψ e μ⁺, μ⁻
      if (pid == 443) {
        h_pt_Jpsi.Fill(pythia.event[i].pT());
        h_eta_Jpsi.Fill(pythia.event[i].eta());
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

    // Massa invariante μ⁺ μ⁻ apenas do decaimento do J/ψ
    for (int i = 0; i < pythia.event.size(); ++i) {
      if (pythia.event[i].id() == 443 && pythia.event[i].isFinal()) {
        int d1 = pythia.event[i].daughter1(); // índice da primeira filha (μ⁺)
        int d2 = pythia.event[i].daughter2(); // índice da segunda filha (μ⁻)

        // Verifica se as filhas são realmente μ⁺ e μ⁻ antes de calcular a massa invariante
        if (d1 > 0 && d2 > 0 &&
            abs(pythia.event[d1].id()) == 13 &&
            abs(pythia.event[d2].id()) == 13) {
          Vec4 pSum = pythia.event[d1].p() + pythia.event[d2].p(); // soma dos quadrivetores das duas filhas
          double massa = pSum.mCalc(); // calcula (soma) a massa invariante das duas filhas
          h_mass_mumu.Fill(massa); // preenche o histograma da massa invariante
        }
      }
    }
  } // fim do loop de eventos

  // Finalização
  pythia.stat(); // estatísticas do Pythia
  outFile.Write(); // escreve os histogramas no arquivo ROOT
  outFile.Close(); // fecha o arquivo ROOT

  // ===============================================================
  // Desenha os histogramas usando o ROOT e salva em um arquivo PDF
  // ===============================================================

  TCanvas c1("c1", "J/psi #rightarrow #mu^{+}#mu^{-}", 1200, 800);
  c1.Divide(3, 3);

  // Cores para diferenciar os histogramas de cada partícula (constantes do ROOT) 
  h_mass_mumu.SetLineColor(kBlue); 
  h_pt_Jpsi.SetLineColor(kBlue);
  h_eta_Jpsi.SetLineColor(kBlue);
  h_pt_muPlus.SetLineColor(kRed);
  h_pt_muMinus.SetLineColor(kGreen+2);
  h_eta_muPlus.SetLineColor(kRed);
  h_eta_muMinus.SetLineColor(kGreen+2);

  // Desenha os histogramas em cada subcanvas

  // Massa invariante
  c1.cd(1); h_mass_mumu.Draw();

  // Momento Transverso (pT) do J/ψ (escala log)
  c1.cd(2); gPad->SetLogy(1);
  h_pt_Jpsi.Draw();

  // Pseudorapidez (η) do J/ψ (linear)
  c1.cd(3); h_eta_Jpsi.Draw();

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
  c1.SaveAs("jpsi_mumu.pdf");

  cout << "\nArquivo ROOT: jpsi_mumu.root" << endl;
  cout << "PDF gerado: jpsi_mumu.pdf" << endl;

  return 0;
}
