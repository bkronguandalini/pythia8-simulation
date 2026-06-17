# PYTHIA8 Tutorial

> **Autor:** Bruno Kron Guandalini  

Tutorial passo a passo para instalação do PYTHIA8 com suporte ao ROOT em sistemas Linux, e execução de simulações de eventos em física de altas energias.

---

## Sumário

1. [Introdução](#1-introdução)
2. [Pré-requisitos](#2-pré-requisitos)
3. [Passo 1 — Instalação do Miniconda](#3-passo-1--instalação-do-miniconda)
4. [Passo 2 — Criação do ambiente Conda com ROOT](#4-passo-2--criação-do-ambiente-conda-com-root)
5. [Passo 3 — Instalação do PYTHIA8](#5-passo-3--instalação-do-pythia8)
6. [Passo 4 — Configuração das variáveis de ambiente](#6-passo-4--configuração-das-variáveis-de-ambiente)
7. [Passo 5 — Rodando exemplo: Simulação Z → µ⁺µ⁻](#7-passo-5--rodando-exemplo-simulação-z--µµ)
8. [Passo 6 — Criação de atalhos](#8-passo-6--criação-de-atalhos)
9. [Solução de problemas](#9-solução-de-problemas)
10. [Referências](#10-referências)

---

## 1. Introdução

O PYTHIA8 é um programa para simulações de eventos em física de altas energias, amplamente utilizado para gerar colisões hadrônicas e processos de decaimento. Este tutorial descreve passo a passo a instalação do PYTHIA8 com suporte ao ROOT em sistemas Linux.

---

## 2. Pré-requisitos

Antes de instalar o PYTHIA8, é necessário ter o **Conda** e o **ROOT** instalados. Caso já possua esses sistemas configurados, pule para o [Passo 3](#5-passo-3--instalação-do-pythia8).

---

## 3. Passo 1 — Instalação do Miniconda

O Miniconda é uma distribuição mínima do Anaconda que facilita o gerenciamento de ambientes Python.

```bash
# Download do instalador do Miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Execução do instalador
bash Miniconda3-latest-Linux-x86_64.sh
```

> **Observação:** Após a instalação, feche e reabra o terminal para que as alterações no `PATH` tenham efeito.

---

## 4. Passo 2 — Criação do ambiente Conda com ROOT

Crie um ambiente Conda dedicado para o ROOT:

```bash
# Criação do ambiente com ROOT
conda create -n root_env root -c conda-forge

# Ativação do ambiente
conda activate root_env
```

- `conda create -n root_env root` → cria um novo ambiente chamado `root_env` com o ROOT instalado
- `-c conda-forge` → especifica o canal conda-forge para download
- `conda activate root_env` → ativa o ambiente recém-criado

---

## 5. Passo 3 — Instalação do PYTHIA8

### 5.1 Download e extração

```bash
# Download do PYTHIA8
wget https://pythia.org/download/pythia83/pythia8315.tgz

# Extração do arquivo
tar xvfz pythia8315.tgz

# Entrar no diretório
cd pythia8315
```

### 5.2 Configuração com suporte ao ROOT

```bash
./configure --with-root=$CONDA_PREFIX --prefix=$HOME/pythia8
```

- `--with-root=$CONDA_PREFIX` → habilita o suporte ao ROOT usando o caminho do ambiente Conda
- `--prefix=$HOME/pythia8` → define o diretório de instalação do PYTHIA8

### 5.3 Compilação e instalação

```bash
# Compilação usando múltiplos núcleos
make -j$(nproc)

# Instalação
make install
```

- `make -j$(nproc)` → compila usando todos os núcleos disponíveis
- `make install` → instala os arquivos no diretório especificado

---

## 6. Passo 4 — Configuração das variáveis de ambiente

Adicione as seguintes linhas ao `~/.bashrc`:

```bash
export PATH=$HOME/pythia8/bin:$PATH
export LD_LIBRARY_PATH=$HOME/pythia8/lib:$LD_LIBRARY_PATH
export PYTHIA8=$HOME/pythia8
export PYTHIA8DATA=$HOME/pythia8/share/Pythia8/xmldoc
```

| Variável | Descrição |
|----------|-----------|
| `PATH` | Permite executar programas do PYTHIA8 de qualquer diretório |
| `LD_LIBRARY_PATH` | Informa ao sistema onde encontrar as bibliotecas |
| `PYTHIA8` | Aponta para o diretório de instalação |
| `PYTHIA8DATA` | Localização dos arquivos de dados do PYTHIA8 |

Recarregue o terminal:

```bash
source ~/.bashrc
```

---

## 7. Passo 5 — Rodando exemplo: Simulação Z → µ⁺µ⁻

**Objetivo:** Criar, compilar e executar uma simulação completa de produção do bóson Z.

### 7.1 Criar o arquivo de simulação

```bash
nano z_mumu_simulation.cc
```

Cole o seguinte código:

```cpp
#include "Pythia8/Pythia.h"
#include "TH1F.h"
#include "TFile.h"
#include "TLorentzVector.h"

using namespace Pythia8;

int main() {
    Pythia pythia;
    pythia.readString("Beams:idA = 2212");               // próton
    pythia.readString("Beams:idB = 2212");               // próton
    pythia.readString("Beams:eCM = 13000.");             // energia de centro de massa (13 TeV)
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on"); // ativa produção de Z
    pythia.readString("23:onMode = off");                // desativa todos os decaimentos do Z
    pythia.readString("23:onIfAny = 13");                // apenas Z -> mu+mu- (13 é o código PDG)
    pythia.init();

    TFile* outfile = new TFile("z_mumu.root", "RECREATE");
    TH1F* h_mll = new TH1F("h_mll", "Mass of mu+mu-; m_{#mu#mu} [GeV]; Events", 100, 60., 120.);

    for (int iEvent = 0; iEvent < 10000; ++iEvent) {
        if (!pythia.next()) continue;

        std::vector<TLorentzVector> muons;
        for (int i = 0; i < pythia.event.size(); ++i) {
            if (pythia.event[i].id() == 13 || pythia.event[i].id() == -13) {
                TLorentzVector mu(pythia.event[i].px(), pythia.event[i].py(),
                                  pythia.event[i].pz(), pythia.event[i].e());
                muons.push_back(mu);
            }
        }

        if (muons.size() == 2) {
            TLorentzVector dimuon = muons[0] + muons[1];
            h_mll->Fill(dimuon.M());
        }
    }

    pythia.stat();
    h_mll->Write();
    outfile->Close();

    return 0;
}
```

Para salvar no editor `nano`:

1. Pressione `Ctrl+O` para salvar
2. Pressione `Enter` para confirmar o nome do arquivo
3. Pressione `Ctrl+X` para sair

### 7.2 Compilar e executar

```bash
# Compilar
g++ z_mumu_simulation.cc -o z_mumu_simulation \
    -I$PYTHIA8/include \
    -L$PYTHIA8/lib -lpythia8 \
    `root-config --cflags --glibs`

# Executar
./z_mumu_simulation
```

### 7.3 Analisar os resultados

```bash
# Abrir o arquivo ROOT gerado
root -l z_mumu.root

# No ambiente ROOT, visualize o histograma:
h_mll->Draw()
```

### 7.4 Entendendo o código

A simulação produz o bóson Z em colisões próton-próton a 13 TeV, forçando o decaimento Z → µ⁺µ⁻:

| Parâmetro | Descrição |
|-----------|-----------|
| `Beams:idA = 2212` e `Beams:idB = 2212` | Definem prótons nos feixes |
| `Beams:eCM = 13000.` | Energia no centro de massa de 13 TeV |
| `WeakSingleBoson:ffbar2gmZ = on` | Ativa produção do bóson Z |
| `23:onMode = off` | Desativa todos os decaimentos padrão do Z |
| `23:onIfAny = 13` | Força Z → µ⁺µ⁻ (código PDG do múon: 13) |

O código reconstrói a massa invariante do par µ⁺µ⁻ e preenche um histograma. O resultado esperado é um pico na massa do Z (~91 GeV).

---

## 8. Passo 6 — Criação de atalhos

### 8.1 Editar o arquivo de configuração

```bash
nano ~/.bashrc
```

### 8.2 Adicionar os atalhos

```bash
# Atalho para ativar root_env e entrar no pythia8315
alias pythia='conda activate root_env && cd ~/pythia8315'

# Atalho para compilar exemplos comuns
alias compilar_z_mumu='g++ z_mumu_simulation.cc -o z_mumu_simulation \
    -I$PYTHIA8/include -I$ROOTSYS/include \
    -L$PYTHIA8/lib -L$ROOTSYS/lib \
    -Wl,-rpath,$PYTHIA8/lib -Wl,-rpath,$ROOTSYS/lib \
    -lpythia8 $(root-config --glibs) \
    -std=c++20 -O2'
```

### 8.3 Recarregar as configurações

```bash
source ~/.bashrc
```

### 8.4 Usar o atalho

```bash
pythia
```

**Resultado esperado:**
```
(root_env) SEU-USUÁRIO@SUA-MÁQUINA:~/pythia8315$
```

> `SEU-USUÁRIO` e `SUA-MÁQUINA` serão automaticamente substituídos pelos nomes do seu sistema.

---

## 9. Solução de problemas

| Problema | Solução |
|----------|---------|
| Erro de biblioteca não encontrada | Verifique se as variáveis de ambiente estão configuradas corretamente |
| Problemas de compilação | Certifique-se de que o ambiente Conda está ativado (`conda activate root_env`) |
| Arquivos de dados não encontrados | Confirme se a variável `PYTHIA8DATA` aponta para o diretório correto |
| Erro no ROOT | Verifique se todas as bibliotecas do ROOT estão disponíveis no ambiente |

---

## 10. Referências

- [Site oficial do PYTHIA](https://pythia.org)
- [Documentação do PYTHIA8](https://pythia.org/latest-manual/)
- [Guia sobre o PYTHIA8 (PDF)](https://pythia.org/download/pdf/pythia8300.pdf)
- [Códigos PDG para partículas](https://pdg.lbl.gov)
