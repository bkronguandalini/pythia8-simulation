"""
Desenvolvido por: Bruno Kron Guandalini
Descrição: Simulação de produção de Z⁰ com decaimento em μ⁺μ⁻ usando Pythia8 e ROOT.
Objetivo: Estudar as propriedades do Z⁰ e seu decaimento em μ⁺μ⁻ a 13 TeV.
"""

import pythia8 # Importa a biblioteca Pythia8
import ROOT # Importa a biblioteca ROOT

def main():
    # ================================================================
    # 1. Configuração do Pythia
    # ================================================================
    
    ROOT.gROOT.SetBatch(ROOT.kTRUE) # Ativa modo batch (sem janelas gráficas)

    pythia = pythia8.Pythia()

    # Definições do feixe 
    pythia.readString("Beams:idA = 2212")      # Define o feixe A (próton PDG 2212)
    pythia.readString("Beams:idB = 2212")      # Define o feixe B (próton PDG 2212)
    pythia.readString("Beams:eCM = 13000.")    # Define a Energia do centro de massa (13 TeV)

    # Processo: produção de Z⁰ via Drell-Yan
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on")

    # Força o decaimento do Z⁰ em μ⁺μ⁻
    pythia.readString("23:onMode = off")       # desliga todos os canais do Z⁰ (PDG 23) 
    pythia.readString("23:onIfAny = 13")       # ativa apenas canal com muon (PDG 13)

    # Reduz mensagens na tela
    pythia.readString("Next:numberShowInfo = 0")
    pythia.readString("Next:numberShowProcess = 0")
    pythia.readString("Next:numberShowEvent = 0")

    if not pythia.init():
        print("Erro na inicialização do Pythia.")
        return 1

    # ================================================================
    # 2. Histogramas e arquivo ROOT
    # ================================================================
    outFile = ROOT.TFile("z_mumu.root", "RECREATE")

    # Massa invariante do Z° (μ⁺μ⁻)
    h_mass_mumu = ROOT.TH1D("h_mass_mumu","Invariant Mass of #mu^{+}#mu^{-}; m_{#mu^{+}#mu^{-}} (GeV/c^{2}); Events", 100, 60.0, 120.0)

    # Momento Transversal do Z°
    h_pt_Z = ROOT.TH1D("h_pt_Z","Transverse Momentum of Z^{0}; p_{T} (GeV/c); Events", 100, 0.0, 200.0)
    
    # Pseudorapidez do Z°
    h_eta_Z = ROOT.TH1D("h_eta_Z","Pseudorapidity of Z^{0}; #eta; Events", 100, -5.0, 5.0)

    # Momento Transversal dos muons (μ⁺ e μ⁻)
    h_pt_muPlus = ROOT.TH1D("h_pt_muPlus", "Transverse Momentum of #mu^{+}; p_{T} (GeV/c); Events", 100, 0.0, 100.0)
    h_pt_muMinus = ROOT.TH1D("h_pt_muMinus", "Transverse Momentum of #mu^{-}; p_{T} (GeV/c); Events", 100, 0.0, 100.0)
    
    # Pseudorapidez dos muons (μ⁺ e μ⁻)
    h_eta_muPlus = ROOT.TH1D("h_eta_muPlus", "Pseudorapidity of #mu^{+}; #eta; Events", 100, -5.0, 5.0)
    h_eta_muMinus = ROOT.TH1D("h_eta_muMinus", "Pseudorapidity of #mu^{-}; #eta; Events", 100, -5.0, 5.0)

    # Cores
    h_mass_mumu.SetLineColor(ROOT.kBlue)
    h_pt_Z.SetLineColor(ROOT.kBlue)
    h_eta_Z.SetLineColor(ROOT.kBlue)
    h_pt_muPlus.SetLineColor(ROOT.kRed)
    h_pt_muMinus.SetLineColor(ROOT.kGreen)
    h_eta_muPlus.SetLineColor(ROOT.kRed)
    h_eta_muMinus.SetLineColor(ROOT.kGreen)

    # ===============================================================
    # 3. Loop de eventos
    # ===============================================================
    nEvent = 10000 # Número de eventos a serem gerados (ajuste conforme necessário)
    nAbort = 5 # Número máximo de falhas consecutivas antes de abortar a geração
    iAbort = 0 # Contador de falhas consecutivas

    # Loop de geração de eventos abortados após muitas falhas
    for iEvent in range(nEvent):
        if not pythia.next():
            iAbort += 1
            if iAbort < nAbort:
                continue
            print("Geração abortada após muitas falhas.")
            break

        # Preenche histogramas de Momento Transversal (pT) e Pseudorapidez (η) das partículas
        for i in range(pythia.event.size()):
            pid = pythia.event[i].id() # Identificador PDG da partícula (pid = particle ID)

            if pid == 23: # Z⁰ (PDG 23)
                h_pt_Z.Fill(pythia.event[i].pT()) # Preenche histograma de pT do Z⁰
                h_eta_Z.Fill(pythia.event[i].eta()) # Preenche histograma de η do Z⁰
            elif pid == -13: # μ⁺ (PDG -13)
                h_pt_muPlus.Fill(pythia.event[i].pT()) # Preenche histograma de pT do μ⁺
                h_eta_muPlus.Fill(pythia.event[i].eta()) # Preenche histograma de η do μ⁺
            elif pid == 13: # μ⁻ (PDG 13)
                h_pt_muMinus.Fill(pythia.event[i].pT()) # Preenche histograma de pT do μ⁻
                h_eta_muMinus.Fill(pythia.event[i].eta()) # Preenche histograma de η do μ⁻

        # Massa invariante do par μ⁺μ⁻ (filhas diretas do Z)
        for i in range(pythia.event.size()):
            if pythia.event[i].id() == 23:
                d1 = pythia.event[i].daughter1()  # Índice da primeira filha
                d2 = pythia.event[i].daughter2()  # Índice da segunda filha

                # Verifica se as filhas são realmente múons
                if (d1 > 0 and d2 > 0
                    and abs(pythia.event[d1].id()) == 13
                    and abs(pythia.event[d2].id()) == 13):
                    pSum = pythia.event[d1].p() + pythia.event[d2].p()
                    massa = pSum.mCalc()
                    h_mass_mumu.Fill(massa)

    pythia.stat() # Estatísticas da geração de eventos

    # =======================================================
    # 4. PDF MULTIPÁGINA
    # =======================================================
    
    # Lista de histogramas com flag para escala log (no eixo Y)
    h_logs = [
        (h_mass_mumu, False),
        (h_pt_Z, True),
        (h_eta_Z, False),
        (h_pt_muPlus, True),
        (h_pt_muMinus, True),
        (h_eta_muPlus, False),
        (h_eta_muMinus, False),
    ]

    # Cria um canvas para desenhar cada histograma em uma página
    c = ROOT.TCanvas("c", "Z^{0} #rightarrow #mu^{+}#mu^{-}", 1200, 800)
    pdf_z = "z_mumu.pdf"

    # Inicia o PDF multipágina
    c.Print(pdf_z + "[")  

    for hist, logy in h_logs:
        c.Clear()  # Limpa o canvas
        c.SetLogy(1 if logy else 0)  # Escala log ou linear
        hist.Draw("HIST")            # Desenha o histograma
        c.Update()                   # Atualiza o canvas 
        c.Print(pdf_z)               # Salva a página


    c.Print(pdf_z + "]")  # Finaliza o PDF

    outFile.Write() # Salva todos os histogramas no arquivo
    outFile.Close() # Fecha o arquivo ROOT

    print("\nArquivo ROOT: z_mumu.root") # Informa o usuário sobre o arquivo ROOT gerado
    print("PDF gerado: z_mumu.pdf") # Informa o usuário sobre o PDF gerado

    return 0 # Retorna 0 para indicar sucesso na execução

# Executa a função main() se o script for executado diretamente
if __name__ == "__main__":
    main()
