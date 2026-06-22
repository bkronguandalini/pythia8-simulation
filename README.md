# PYTHIA8 Tutorial

> **Author:** Bruno Kron Guandalini  

Step-by-step tutorial for installing PYTHIA8 with ROOT support on Linux systems, and running event simulations in high-energy physics.

---

## Table of Contents

1. [Introduction](#1-introduction)
2. [Prerequisites](#2-prerequisites)
3. [Step 1 — Installing Miniconda](#3-step-1--installing-miniconda)
4. [Step 2 — Creating the Conda environment with ROOT](#4-step-2--creating-the-conda-environment-with-root)
5. [Step 3 — Installing PYTHIA8](#5-step-3--installing-pythia8)
6. [Step 4 — Setting up environment variables](#6-step-4--setting-up-environment-variables)
7. [Step 5 — Running an example: Z → µ⁺µ⁻ simulation](#7-step-5--running-an-example-z--µµ-simulation)
8. [Step 6 — Creating shortcuts (Optional)](#8-step-6--creating-shortcuts-optional)
9. [Troubleshooting](#9-troubleshooting)
10. [References](#10-references)

---

## 1. Introduction

PYTHIA8 is a program for event simulations in high-energy physics, widely used to generate hadronic collisions and decay processes. This tutorial describes step by step the installation of PYTHIA8 with ROOT support on Linux systems.

---

## 2. Prerequisites

Before installing PYTHIA8, you need to have **Conda** and **ROOT** installed. If you already have these systems set up, skip to [Step 3](#5-step-3--installing-pythia8).

---

## 3. Step 1 — Installing Miniconda

Miniconda is a minimal distribution of Anaconda that facilitates Python environment management.

```bash
# Download the Miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run the installer
bash Miniconda3-latest-Linux-x86_64.sh
```

> **NOTE:** After installation, close and reopen the terminal for the `PATH` changes to take effect.

---

## 4. Step 2 — Creating the Conda environment with ROOT

Create a dedicated Conda environment for ROOT:

```bash
# Create the environment with ROOT
conda create -n root_env root -c conda-forge

# Activate the environment
conda activate root_env
```

- `conda create -n root_env root` → creates a new environment called `root_env` with ROOT installed
- `-c conda-forge` → specifies the conda-forge channel for downloading
- `conda activate root_env` → activates the newly created environment

---

## 5. Step 3 — Installing PYTHIA8

### 5.1 Download and extraction

```bash
# Download PYTHIA8
wget https://pythia.org/download/pythia83/pythia8315.tgz

# Extract the file
tar xvfz pythia8315.tgz

# Enter the directory
cd pythia8315
```

### 5.2 Configuration with ROOT support

```bash
./configure --with-root=$CONDA_PREFIX --prefix=$HOME/pythia8
```

- `--with-root=$CONDA_PREFIX` → enables ROOT support using the Conda environment path
- `--prefix=$HOME/pythia8` → sets the PYTHIA8 installation directory

### 5.3 Compilation and installation

```bash
# Compile using multiple cores
make -j$(nproc)

# Install
make install
```

- `make -j$(nproc)` → compiles using all available cores
- `make install` → installs the files to the specified directory

---

## 6. Step 4 — Setting up environment variables

Add the following lines to `~/.bashrc`:

```bash
export PATH=$HOME/pythia8/bin:$PATH
export LD_LIBRARY_PATH=$HOME/pythia8/lib:$LD_LIBRARY_PATH
export PYTHIA8=$HOME/pythia8
export PYTHIA8DATA=$HOME/pythia8/share/Pythia8/xmldoc
```

| Variable | Description |
|----------|-----------|
| `PATH` | Allows running PYTHIA8 programs from any directory |
| `LD_LIBRARY_PATH` | Tells the system where to find the libraries |
| `PYTHIA8` | Points to the installation directory |
| `PYTHIA8DATA` | Location of PYTHIA8 data files |

Reload the terminal:

```bash
source ~/.bashrc
```

---

## 7. Step 5 — Running an example: Z → µ⁺µ⁻ simulation

**Objective:** Create, compile, and run a complete simulation of Z boson production.

### 7.1 Create the simulation file

```bash
nano test.cpp
```
Tip: The name `test` can be replaced by any name you prefer. For example, you can use `nano z_mumu_simulation.cpp`.

Paste the following code:

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

To save in the `nano` editor:

1. Press `Ctrl+O` to salve
2. Press `Enter` to confirm the filename
3. Press `Ctrl+X` to exit

### 7.2 Compile and run

```bash
# Compile
g++ test.cpp -o test \
    -I$PYTHIA8/include \
    -L$PYTHIA8/lib -lpythia8 \
    `root-config --cflags --glibs`

# Run
./test
```

Tip: Just remember to adjust the compilation command according to the name you choose in [7.1](#7.1-create-the-simulation-file) (replace `test.cpp` and `-o test`).

### 7.3 Analyze the results

```bash
# Open the generated ROOT file
root -l z_mumu.root

# In the ROOT environment, visualize the histogram:
h_mll->Draw()
```

### 7.4 Understanding the code

The simulation produces Z bosons in proton-proton collisions at 13 TeV, forcing the decay Z → µ⁺µ⁻:

| Parâmetro | Descrição |
|-----------|-----------|
| `Beams:idA = 2212` and `Beams:idB = 2212` | Set protons in the beams |
| `Beams:eCM = 13000.` | ECenter-of-mass energy of 13 TeV |
| `WeakSingleBoson:ffbar2gmZ = on` | Enables Z boson production |
| `23:onMode = off` | Disables all standard Z decays |
| `23:onIfAny = 13` | Force Z → µ⁺µ⁻ (PDG code for muon: 13) |

The code reconstructs the invariant mass of the µ⁺µ⁻ pair and fills a histogram. The expected result is a peak at the Z mass (~91 GeV).

---

## 8. Step 6 — Creating shortcuts (Optional)

### 8.1 Edit the configuration file

```bash
nano ~/.bashrc
```

### 8.2 Add the shortcuts

```bash
# Shortcut to activate root_env and enter pythia8315
alias pythia='conda activate root_env && cd ~/pythia8315'

# Shortcut to compile common examples
alias compile_test='g++ test.cpp -o test \
    -I$PYTHIA8/include -I$ROOTSYS/include \
    -L$PYTHIA8/lib -L$ROOTSYS/lib \
    -Wl,-rpath,$PYTHIA8/lib -Wl,-rpath,$ROOTSYS/lib \
    -lpythia8 $(root-config --glibs) \
    -std=c++20 -O2'
```

Tip: Just remember to adjust the compilation command according to the name you choose in [7.1](#7.1-create-the-simulation-file) (replace `test.cpp` and `-o test`).

### 8.3 Reload the configurations

```bash
source ~/.bashrc
```

### 8.4 se the shortcut

```bash
pythia
```

**Expected result:**
```
(root_env) YOUR-USER@YOUR-MACHINE:~/pythia8315$
```

> `YOUR-USER` and `YOUR-MACHINE` will be automatically replaced with your system's names.

---

## 9. Troubleshooting

| Problem | Solution |
|----------|---------|
| Library not found error |	Verify that the environment variables are set correctly |
| Compilation issues | Make sure the Conda environment is activated (`conda activate root_env`) |
| Data files not found | Confirm that the `PYTHIA8DATA` variable points to the correct directory |
| ROOT error | Check that all ROOT libraries are available in the environment |

---

## 10. References

- [Miniconda Documetation](https://www.anaconda.com/docs/getting-started/miniconda/main)
- [PYTHIA Official Website](https://pythia.org)
- [PYTHIA8 Documetation](https://pythia.org/latest-manual/)
- [PYTHIA8 Guide (PDF)](https://pythia.org/download/pdf/pythia8300.pdf)
- [PDG Particle Codes](https://pdg.lbl.gov)
