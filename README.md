# PYTHIA8 Tutorial

> **Author:** Bruno Kron Guandalini

[![Linux Mint](https://img.shields.io/badge/Linux-Mint-87CF3E?logo=linuxmint&logoColor=white)](https://linuxmint.com/)
[![ROOT](https://img.shields.io/badge/ROOT-6-blue)](https://root.cern/)
[![PYTHIA8](https://img.shields.io/badge/PYTHIA-8-orange)](https://www.pythia.org/)
[![C++20](https://img.shields.io/badge/C%2B%2B-20-red)](https://cppreference.com/cpp/20)
[![Python](https://img.shields.io/badge/Python-3.10-yellow)](https://www.python.org/)

Step-by-step tutorial for installing PYTHIA8 with ROOT support on Linux systems and running event simulations in high-energy physics.

> **Tested on Linux Mint.** The same instructions should work on Ubuntu and other Debian-based distributions, perhaps with minor modifications.
---

# Table of Contents

1. [Introduction](#1-introduction)
2. [Prerequisites](#2-prerequisites)
3. [C++ Interface](#3-c-interface)
   - [3.1 Installing Miniconda](#31-installing-miniconda)
   - [3.2 Creating the ROOT environment](#32-creating-the-conda-environment-with-root)
   - [3.3 Installing PYTHIA8](#33-installing-pythia8)
      - [3.3.1 Download and extraction](#331-download-and-extraction)
      - [3.3.2 Configuration with ROOT support](#332-configuration-with-root-support)
      - [3.3.3 Compilation and installation](#333-compilation-and-installation)
   - [3.4 Setting up environment variables](#34-setting-up-environment-variables)
   - [3.5 Running an example: Z → µ⁺µ⁻ simulation](#35-running-an-example-z--µµ-simulation)
      - [3.5.1 Create the simulation file](#351-create-the-simulation-file)
      - [3.5.2 Compile and run](#352-compile-and-run)
      - [3.5.3 Analyze the results](#353-analyze-the-results)
      - [3.5.4 Understanding the code](#354-understanding-the-code)
   - [3.6 Creating shortcuts (Optional)](#36-creating-shortcuts-optional)
      - [3.6.1 Edit the configuration file](#361-edit-the-configuration-file)
      - [3.6.2 Add the shortcuts](#362-add-the-shortcuts)
      - [3.6.3 Reload the configurations](#363-reload-the-configurations)
      - [3.6.4 Use the shortcut](#364-use-the-shortcut)
4. [Python Interface](#4-python-interface)
   - [4.1 Creating a clean Conda environment](#41-creating-a-clean-conda-environment)
   - [4.2 Installing PYTHIA8 and ROOT via conda-forge](#42-installing-pythia8-and-root-via-conda-forge)
   - [4.3 Running a script](#43-running-a-script)
   - [4.4 Quick library check (Optional)](#44-quick-library-check-optional)
5. [Troubleshooting](#5-troubleshooting)
   - [5.1 C++](#51-c)
   - [5.2 Python](#52-python)
6. [References](#6-references)
---

## 1. Introduction

PYTHIA8 is a program for event simulations in high-energy physics, widely used to generate hadronic collisions and decay processes. This tutorial describes step by step the installation of PYTHIA8 with ROOT support on Linux systems.

---

## 2. Prerequisites

Before installing PYTHIA8, you need to have **Conda** and **ROOT** installed. If you already have these systems set up, skip to [Installing PYTHIA8](#33-installing-pythia8).

---

## 3. C++ Interface

This section describes how to compile and run simulations using C++ and ROOT.

---

### 3.1 Installing Miniconda

Miniconda is a minimal distribution of Anaconda that facilitates Python environment management.

```bash
# Download the Miniconda installer
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Run the installer
bash Miniconda3-latest-Linux-x86_64.sh
```

> **NOTE:** After installation, close and reopen the terminal for the `PATH` changes to take effect.

---

### 3.2 Creating the Conda environment with ROOT

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

### 3.3 Installing PYTHIA8

#### 3.3.1 Download and extraction

```bash
# Download PYTHIA8
wget https://pythia.org/download/pythia83/pythia8315.tgz

# Extract the file
tar xvfz pythia8315.tgz

# Enter the directory
cd pythia8315
```

#### 3.3.2 Configuration with ROOT support

```bash
./configure --with-root=$CONDA_PREFIX --prefix=$HOME/pythia8
```

- `--with-root=$CONDA_PREFIX` → enables ROOT support using the Conda environment path
- `--prefix=$HOME/pythia8` → sets the PYTHIA8 installation directory

#### 3.3.3 Compilation and installation

```bash
# Compile using multiple cores
make -j$(nproc)

# Install
make install
```

- `make -j$(nproc)` → compiles using all available cores
- `make install` → installs the files to the specified directory

---

### 3.4 Setting up environment variables

Add the following lines to `~/.bashrc`:

```bash
export PATH=$HOME/pythia8/bin:$PATH
export LD_LIBRARY_PATH=$HOME/pythia8/lib:$LD_LIBRARY_PATH
export PYTHIA8=$HOME/pythia8
export PYTHIA8DATA=$HOME/pythia8/share/Pythia8/xmldoc
```

| Variable | Description |
|------------|------------|
| `PATH` | Allows running PYTHIA8 programs from any directory |
| `LD_LIBRARY_PATH` | Tells the system where to find the libraries |
| `PYTHIA8` | Points to the installation directory |
| `PYTHIA8DATA` | Location of the XML files |

Reload the terminal:

```bash
source ~/.bashrc
```

---

### 3.5 Running an example: Z → µ⁺µ⁻ simulation

**Objective:** Create, compile, and run a complete simulation of Z boson production.

#### 3.5.1 Create the simulation file

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
    pythia.readString("Beams:idA = 2212");               // proton
    pythia.readString("Beams:idB = 2212");               // proton
    pythia.readString("Beams:eCM = 13000.");             // center-of-mass energy (13 TeV)
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on"); // enable Z production
    pythia.readString("23:onMode = off");                // disable all Z decays
    pythia.readString("23:onIfAny = 13");                // only Z -> mu+mu- (PDG Code: 13)
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

---

#### 3.5.2 Compile and run

```bash
# Compile
g++ test.cpp -o test \
    -I$PYTHIA8/include \
    -L$PYTHIA8/lib -lpythia8 \
    `root-config --cflags --glibs`

# Run
./test
```

Tip: Just remember to adjust the compilation command according to the name you choose in [3.5.1](#351-create-the-simulation-file) (replace `test.cpp` and `-o test`).

#### 3.5.3 Analyze the results

```bash
# Open the generated ROOT file
root -l z_mumu.root

# In the ROOT environment, visualize the histogram:
h_mll->Draw()
```

#### 3.5.4 Understanding the code

The simulation produces Z bosons in proton-proton collisions at 13 TeV, forcing the decay Z → µ⁺µ⁻:

| Parameter | Description |
|-----------|-----------|
| `Beams:idA = 2212` and `Beams:idB = 2212` | Set protons in the beams |
| `Beams:eCM = 13000.` | Center-of-mass energy of 13 TeV |
| `WeakSingleBoson:ffbar2gmZ = on` | Enables Z boson production |
| `23:onMode = off` | Disables all standard Z decays |
| `23:onIfAny = 13` | Force Z → µ⁺µ⁻ (PDG code for muon: 13) |

The code reconstructs the invariant mass of the µ⁺µ⁻ pair and fills a histogram. The expected result is a peak at the Z mass (~91 GeV).

---

### 3.6 Creating shortcuts (Optional)

#### 3.6.1 Edit the configuration file

```bash
nano ~/.bashrc
```

#### 3.6.2 Add the shortcuts

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

Tip: Just remember to adjust the compilation command according to the name you choose in [3.5.1](#351-create-the-simulation-file) (replace `test.cpp` and `-o test`).

#### 3.6.3 Reload the configurations

```bash
source ~/.bashrc
```

#### 3.6.4 Use the shortcut

```bash
pythia
```

**Expected result:**
```
(root_env) YOUR-USER@YOUR-MACHINE:~/pythia8315$
```

> `YOUR-USER` and `YOUR-MACHINE` will be automatically replaced with your system's names.

---

## 4. Python Interface

This interface is the simplest and safest way to get started. Everything is installed automatically via `conda-forge` — no manual compilation needed. Recommended for anyone who prefers writing simulations in Python.

> **Important:** Do not try to compile Pythia8 or ROOT manually when using this approach. Using conda-forge avoids library conflicts and `undefined symbol` errors.

---

### 4.1 Creating a clean Conda environment

Create a separate environment to avoid conflicts with other installations:

```bash
conda create -n python_env python=3.10 -c conda-forge
conda activate python_env
```

> **Tip:** Use Python 3.9 or 3.10 — these are stable and well-supported versions.

**Expected result:**
```
(python_env) YOUR-USER@YOUR-MACHINE:~$
```

> `YOUR-USER` and `YOUR-MACHINE` will be automatically replaced with your system's names.
---

### 4.2 Installing PYTHIA8 and ROOT via conda-forge

```bash
conda install pythia8 root -c conda-forge
```

This single command installs:
- The Pythia8 and ROOT C++ libraries
- The Python bindings for both (`pythia8` and `ROOT`)

---

### 4.3 Running a script

With the environment activated (`conda activate python_env`), simply run:

```bash
python test.py
```

If everything is set up correctly, the program will run and generate the output files (`.root`, `.pdf`, etc.).

---

### 4.4 Quick library check (Optional)

Before running your full script, verify that both libraries are working:

```bash
python -c "import pythia8; print('Pythia OK')"

python -c "import ROOT; print('ROOT OK')"
```

---

## 5. Troubleshooting

### 5.1 C++

| Problem | Solution |
|----------|---------|
| Library not found error |	Verify that the environment variables are set correctly |
| Compilation issues | Make sure the Conda environment is activated (`conda activate root_env`) |
| Data files not found | Confirm that the `PYTHIA8DATA` variable points to the correct directory |
| ROOT error | Check that all ROOT libraries are available in the environment |

### 5.2 Python

| Problem | Cause | Solution |
|---------|-------|----------|
| `NameError: name 'null' is not defined` | Tried to run a `.ipynb` (JSON) file with Python | Extract the code with `nbconvert` or copy it manually to a `.py` file |
| `ImportError: undefined symbol: _Py_HashPointer` | Python is loading a manually compiled ROOT built with a different Python version | Deactivate `$ROOTSYS` and `$PYTHONPATH` pointing to manual installations. Use only the ROOT from conda-forge |
| `ModuleNotFoundError: No module named 'pythia8'` | Package not installed or wrong environment active | Activate the correct environment (`conda activate python_env`) and install via conda-forge |
| Compilation error with `g++` | Trying to compile Python code with a C++ compiler | Python is interpreted: use `python script.py`, not `g++` |

---

## 6. References

- [Miniconda Documetation](https://www.anaconda.com/docs/getting-started/miniconda/main)
- [PYTHIA Official Website](https://pythia.org)
- [PYTHIA8 Documetation](https://pythia.org/latest-manual/)
- [PYTHIA8 Guide (PDF)](https://pythia.org/download/pdf/pythia8300.pdf)
- [PDG Particle Codes](https://pdg.lbl.gov)
