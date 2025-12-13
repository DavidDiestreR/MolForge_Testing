# MolForge Testing — **CPU** (Conda) · WSL + Windows (PC) o Linux (Lab)

Aquest repositori està configurat per executar **MolForge exclusivament en CPU** utilitzant entorns **Conda**.  
L'entorn oficial de MolForge (`MolForge_env`) està basat en el `environment.yml` original del projecte i està pensat **només per a inferència en CPU**.  
Les tasques de **preprocessat** (com convertir SMILES → fingerprints) es realitzen amb l'entorn auxiliar `molforge-tools`.

---

## 📁 Estructura

```
MolForge_Testing/
├─ envs/
│  ├─ molforge/environment.yml            # environment oficial de MolForge (amb la instal·lació del paquet afegida)
│  └─ tools/environment.yml               # environment personal per fer data science
├─ data/
│  ├─ SMILES/                             # entrades amb SMILES
│  │  └─ raw/
│  │     ├─ CoCoGraph                     # datasets del repo de CoCoGraph (descarregar a part en "📄")
│  │     ├─ MolForge                      # datasets del repo de MolForge
│  │     └─ PubChem                       # dataset de PubChem
│  ├─ MolForge_input/                     # fingerprints (input per a MolForge)
│  ├─ MolForge_output/                    # resultats de MolForge
│  └─ analysis_output/
├─ notebooks/
│  ├─ preprocessing/
│  │  ├─ noise.ipynb
│  │  ├─ raw_to_SMILES.ipynb
│  │  └─ SMILES_to_MolForge_input.ipynb
│  ├─ MolForge/
│  │  ├─ data/sp/                         # vocabulari (importat del repo de MolForge)
│  │  ├─ saved_models/                    # checkpoints del repo de MolForge (descarregar a part en "🧠")
│  │  └─ fps_to_smiles_MolForge.ipynb
│  └─ analysis/
│     └─ output_analysis.ipynb
├─ src/
│  ├─ smiles_to_fp.py                     # funcions per passar SMILES a fingerprints amb el format desitjat
│  └─ fingerprints.py                     # codi que transforma objectes mol de rdkit a SMILES (importat del repo de MolForge)
└─ .gitignore
```

---

## 🧩 Requisits

1) **WSL2 + Ubuntu 22.04** instal·lats a Windows (en cas d'estar en PC Windows, veure guia a baix en "🐧").  
2) **Conda/Miniconda** instal·lat (veure guia d’instal·lació a baix en "📦").  
3) **Environment oficial de MolForge**: (NO CAL FER-HO perquè ja està importat en aquest repositori)
   - Copia el `environment.yml` del repo oficial de MolForge a `envs/molforge/environment.yml`.
   - Afegeix a la secció `- pip:` la instal·lació del paquet:
     ```yaml
     - "MolForge @ git+https://github.com/knu-lcbc/MolForge.git"
     ```
4) **Checkpoints del repositori de MolForge** descarregats (veure guia a baix en "🧠")

---

## 🐧 Instal·lar Ubuntu (WSL) per primera vegada

**En PowerShell (Administrador):**
```powershell
wsl --install -d Ubuntu-22.04
wsl --update
wsl -l -v      # ha de mostrar Ubuntu amb VERSION 2
```

**Obrir Ubuntu ja dins del projecte:**
```powershell
# PC (WSL):
wsl -d Ubuntu --cd /mnt/d/MolForge_Testing
```

---

## 📦 Instal·lar Conda a Ubuntu (WSL)

En la **terminal d'Ubuntu**:
```bash
# instal·lació ràpida (no interactiva)
curl -fsSL -o /tmp/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash /tmp/miniconda.sh -b -p "$HOME/miniconda3"
"$HOME/miniconda3/bin/conda" init bash
exec bash

# comprovar
conda --version
```

---

## 📦 Instal·lar Conda a Ubuntu (lab) — s'ha de fer una vegada a cada PC

Descarregar **[Miniconda](https://repo.anaconda.com/archive/Anaconda3-2025.06-0-Linux-x86_64.sh)**

Des de downloads:
```bash
# donar-li permisos d'execució
chmod +x Anaconda3-2025.06-0-Linux-x86_64.sh

# executar
./Anaconda3-2025.06-0-Linux-x86_64.sh

# comprobar
conda
```

---

## 📦 Instal·lar Conda a Windows

Instal·lar **[Miniconda](https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe)** (Selecciona "Just Me" i deixa la ruta per defecte).

Un cop instal·lat MiniConda per primera vegada obrir AnacondaPrompt i escriure "conda init bash" (D'aquesta manera es podrà fer servir comta des de Git Bash)

```bash
conda init bash
```

---

## 🧠 Checkpoints i carpeta `saved_models`

`saved_models/` està **buida** perquè els **checkpoints** de MolForge son massa pesats per a versionar-los en GitHub.

Per utilitzar els checkpoints del model ja entrenat, descarrega'ls del repo oficial de MolForge i col·loca'ls a `saved_models/` (mantenint els noms d'arxiu):

- [**top-performing**](https://drive.google.com/uc?id=1zl6HBdwYsnA4JcnOi1o6OmcrRDB5iySK) — models amb millor rendiment.
- [**all the other models**](https://drive.google.com/uc?id=1jCtbc9lMacCyiZ3iZFEtFgOfOQYtWEuD) — la resta de models disponibles.

Després de descarregar, descomprimeix i verifica que el checkpoint que utilitzaràs existeix a `saved_models/`.

---

## 📄 raw smiles de CoCoGraph

El directori `SMILES/raw/CoCoGraph/` està **buit** en aquest repositori perquè alguns dels datasets importats de CoCoGraph són massa pesants per versionar-los a GitHub.

Per tenir-los al projecte, descarrega’ls del repo de CoCoGraph i col·loca’ls a `CoCoGraph/` (mantenint els noms dels fitxers):

- [**novel_molecules.txt**](https://github.com/manurubo/CoCoGraph/raw/refs/heads/main/Data/generated_database/novel_molecules.csv?download=) — molècules novel generades per CoCoGraph.
- [**molecules_lt70atoms_annotated.txt**](https://github.com/manurubo/CoCoGraph/raw/refs/heads/main/Data/molecules_lt70atoms_annotated.csv?download=) — algunes molècules d’entrenament de CoCoGraph.

---

## 🛠️ Entorns (Conda)

### A) Entorn **MolForge**

> Farem servir el `environment.yml` oficial.

```bash
# PC (WSL)
wsl -d Ubuntu --cd /mnt/d/MolForge_Testing   # PC
cd /mnt/d/MolForge_Testing   # PC (ja amb Ubuntu obert)

# Laboratori (Linux)
cd /export/home/ddiestre/MolForge_Testing

conda env create -f envs/molforge/environment.yml -n MolForge_env
conda activate MolForge_env
```

### B) Entorn **molforge-tools** (RDKit + pandas)

El pots crear i fer-lo servir **a Windows** o **a Ubuntu** (mateix YAML).

**Windows (PowerShell):**
```powershell
cd D:\MolForge_Testing
conda env create -f envs\tools\environment.yml
conda activate molforge-tools
```

**Ubuntu (WSL o Lab):**
```bash
# PC (WSL)
wsl -d Ubuntu --cd /mnt/d/MolForge_Testing   # PC
cd /mnt/d/MolForge_Testing   # PC (ja amb Ubuntu obert)

# Laboratori (Linux)
cd /export/home/ddiestre/MolForge_Testing

conda env create -f envs/tools/environment.yml
conda activate molforge-tools
```

---

## ✅ Comprovacions d’instal·lació

**MolForge:**
```bash
conda activate MolForge_env
python - << 'PY'
import torch
print("torch:", torch.__version__)
print("build CUDA:", torch.version.cuda)
print("cuda available?:", torch.cuda.is_available())
print("selected device:", "cuda" if torch.cuda.is_available() else "cpu")
from MolForge import main as _mf
print("MolForge import OK (MolForge)")
PY
```

**molforge-tools:**
```bash
conda activate molforge-tools
python - << 'PY'
import pandas as pd
print("pandas:", pd.__version__)
try:
    from rdkit import Chem
    print("rdkit MolFromSmiles test:", Chem.MolFromSmiles("CCO") is not None)
except Exception as e:
    print("RDKit import FAILED:", e)
PY
```

---

## 🔄 Actualitzar entorns

- **tools (Conda):**
  ```bash
  conda env update -f envs/tools/environment.yml --prune
  ```
- **MolForge_env (Conda):** si cambia el YAML oficial, el més net és recrear:
  ```bash
  conda env remove -n MolForge_env
  conda env create -f envs/molforge/environment.yml -n MolForge_env
  ```

---

## 📓 Executar notebooks amb Jupyter

Activa primer el teu entorn `MolForge_env`/`molforge-tools` y executa:

```bash
# En MolForge_env (línia per línia)
conda install -y ipykernel
conda install -y jupyterlab
conda install -y ipywidgets
python -m ipykernel install --user --name MolForge_env --display-name "Python (MolForge_env)"

# En molforge-tools (els paquets ja venen instal·lats al yml)
python -m ipykernel install --user --name molforge-tools --display-name "Python (molforge-tools)"
```

Una cop instal·lat l'entorn, executa:
```bash
jupyter lab --no-browser --ip=0.0.0.0
```

> Obre la URL amb token que imprimeix Jupyter al teu navegador de Windows. Dins de JupyterLab, seleccioneu el kernell **Python (MolForge_env WSL)** per executar els notebooks amb aquest entorn.
---

### 🚨 Neteja de kernels

```bash
jupyter kernelspec list
jupyter kernelspec uninstall MolForge_env -y
jupyter kernelspec uninstall molforge-tools -y
```

---

## 🔁 Flux de treball

Executar **en aquest ordre exactament**:

```
1) raw_to_SMILES.ipynb
2) SMILES_to_MolForge_input.ipynb
3) noise.ipynb   (opcional)
4) fps_to_smiles_MolForge.ipynb
5) output_analysis.ipynb
```

### 1️⃣ raw_to_SMILES.ipynb — RAW → dataset d'SMILES amb el format desitjat  
### 2️⃣ SMILES_to_MolForge_input.ipynb — SMILES → fingerprints  
### 3️⃣ noise.ipynb — Afegeix soroll als fingerprints (opcional)  
### 4️⃣ fps_to_smiles_MolForge.ipynb — Fingerprints → SMILES predits  
### 5️⃣ output_analysis.ipynb — Avaluació i mètriques  

---

