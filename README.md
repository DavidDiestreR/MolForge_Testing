# MolForge Testing — **CPU-only** (Conda) · WSL + Windows

Este repo está preparado para **inferir en CPU**, sin usar GPU/CUDA, respetando el `environment.yml` **oficial** de MolForge.  
La idea: ejecutar **MolForge** desde **Ubuntu (WSL)** con el entorno `MolForge_env` creado a partir del YAML oficial, y usar el entorno de utilidades `molforge-tools` con **Conda** (puede ser en Windows o en Ubuntu).

---

## 📁 Estructura

```
MolForge_Testing/
├─ envs/
│  ├─ molforge/environment.yml      # environment oficial de MolForge
│  └─ tools/environment.yml         # RDKit + pandas (ligero)
├─ data/
│  ├─ SMILES/                       # entradas con SMILES
│  ├─ MolForge_input/               # fingerprints generados (input para MolForge)
│  └─ MolForge_output/              # resultados de MolForge
├─ scripts/
│  ├─ smiles_to_fps.py              # convierte SMILES → fingerprints (CPU)
│  └─ run_molforge.py               # ejecuta MolForge (CPU) fila a fila y guarda resultados
├─ notebooks/
│  ├─ 01_smiles_to_fps.ipynb        # guía paso a paso (CPU)
│  └─ 02_run_molforge_cpu.ipynb     # guía paso a paso (CPU)
└─ .gitignore
```

---

## 🧩 Requisitos

1) **WSL2 + Ubuntu 22.04** instalados en Windows (ver guía abajo).  
2) **Conda/Miniconda** instalado dentro de **Ubuntu**.  
3) **Environment oficial de MolForge**:
   - Copia el `environment.yml` **del repo oficial de MolForge** a `envs/molforge/environment.yml`.
   - **Añade en la sección `- pip:` la instalación del paquete**, por ejemplo:
     ```yaml
     - "MolForge @ git+https://github.com/knu-lcbc/MolForge.git"
     ```
   - Este proyecto asume **CPU** (no CUDA).

---

## 🐧 Instalar Ubuntu (WSL) por primera vez

**En PowerShell (Administrador):**
```powershell
wsl --install -d Ubuntu-22.04
wsl --update
wsl -l -v      # debe mostrar Ubuntu con VERSION 2
```

**Rutas:**
- Windows/PowerShell → `D:\MolForge_Testing`  
- Ubuntu/WSL → `/mnt/d/MolForge_Testing`

Abrir Ubuntu ya dentro del proyecto:
```powershell
wsl -d Ubuntu --cd /mnt/d/MolForge_Testing
```

---

## 📦 Instalar Conda en Ubuntu (WSL)

En la **terminal de Ubuntu**:
```bash
# instalación rápida (no interactiva)
curl -fsSL -o /tmp/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash /tmp/miniconda.sh -b -p "$HOME/miniconda3"
"$HOME/miniconda3/bin/conda" init bash
exec bash

# comprobar
conda --version
```

---

## 🛠️ Entornos

### A) Entorno **MolForge** (Ubuntu/WSL, CPU)

> Usaremos el `environment.yml` oficial **intacto**, solo añadiendo la línea `pip` para instalar MolForge (ver Requisitos).

```bash
# dentro de Ubuntu, en la carpeta del proyecto
cd /mnt/d/MolForge_Testing

conda env create -f envs/molforge/environment.yml -n MolForge_env
conda activate MolForge_env

# Forzar modo CPU en esta sesión
export CUDA_VISIBLE_DEVICES=-1
```

### B) Entorno **molforge-tools** (RDKit + pandas)

Puedes crearlo y usarlo **en Windows** o **en Ubuntu** (mismo YAML).

**Windows (PowerShell):**
```powershell
cd D:\MolForge_Testing
conda env create -f envs\tools\environment.yml
conda activate molforge-tools
```

**Ubuntu (WSL):**
```bash
conda env create -f envs/tools/environment.yml
conda activate molforge-tools
```

**Actualizar el entorno de tools cuando cambies el YAML:**
```bash
conda env update -f envs/tools/environment.yml --prune
```

---

## 🔁 Flujo de trabajo (solo CPU)

### 1) SMILES → Fingerprints (RDKit)

**Opción Notebook (recomendada la primera vez):**  
Abre `notebooks/01_smiles_to_fps.ipynb` y sigue las celdas.  
Entrada: CSV/Parquet con columna `smiles` (opcional `id`).  
Salida: fichero en `data/MolForge_input/` con columnas `id`, `smiles`, `fp_0000...`.

**Opción Script (rápido/automatizable):**
```bash
conda activate molforge-tools
python scripts/smiles_to_fps.py   --input data/SMILES/molecules.csv   --smiles-col smiles   --fp morgan --radius 2 --nBits 2048   --output data/MolForge_input/morgan_2048.parquet
```

### 2) Fingerprints → MolForge (CPU)

**Opción Notebook:**  
Abre `notebooks/02_run_molforge_cpu.ipynb` y define:
- `fps_path` → fichero de `data/MolForge_input/`
- `checkpoint_path` → ruta a tu `.pth`
- `fp_name` → p. ej. `ECFP4`
- `model_type` → `smiles` (o `selfies`)
- `decode` → `greedy` (o `beam` si tu repo lo soporta)

**Opción Script (CPU):**
```bash
conda activate MolForge_env
export CUDA_VISIBLE_DEVICES=-1

python scripts/run_molforge.py   --fps data/MolForge_Input/morgan_2048.parquet   --checkpoint /ruta/a/tu/checkpoint.pth   --fp-name ECFP4   --model-type smiles   --decode greedy   --out data/MolForge_output/molforge_outputs.parquet
```

---

## ✅ Comprobaciones de instalación

**MolForge_env (Ubuntu/WSL):**
```bash
conda activate MolForge_env
export CUDA_VISIBLE_DEVICES=-1
python - << 'PY'
import torch
print("cuda available?:", torch.cuda.is_available())  # esperado False
from MolForge import main as _mf
print("MolForge import OK (MolForge)")
PY
```

**molforge-tools (Windows o Ubuntu):**
```bash
conda activate molforge-tools
python - << 'PY'
import pandas as pd
print("pandas:", pd.__version__)
from rdkit import Chem
print("rdkit MolFromSmiles test:", Chem.MolFromSmiles("CCO") is not None)
PY
```

---

## 📂 Gestión de datos

- `data/SMILES/` → entradas con SMILES (`.csv`/`.parquet`).  
- `data/MolForge_input/` → fingerprints generados (input a MolForge).  
- `data/MolForge_output/` → resultados de MolForge (SMILES generados, logs).  

**Sugerencia de nombres:** incluye el tipo/size del FP (`morgan_2048.parquet`) y fecha/modo en la salida (`molforge_outputs_YYYYMMDD.parquet`).

---

## 🔄 Actualizar entornos

- **tools (Windows o Ubuntu):**
  ```bash
  conda env update -f envs/tools/environment.yml --prune
  ```
- **MolForge_env (Ubuntu):** si cambia el YAML oficial, lo más limpio es recrear:
  ```bash
  conda env remove -n MolForge_env
  conda env create -f envs/molforge/environment.yml -n MolForge_env
  ```
