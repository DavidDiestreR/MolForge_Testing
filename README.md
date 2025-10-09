# MolForge Testing — **CPU-only** (Conda - en PC **o** venv - en lab) · WSL + Windows (PC) o Linux (Lab)

Este repo está preparado para **inferir en CPU**, sin usar GPU/CUDA, respetando el `environment.yml` **oficial** de MolForge.  
La idea: ejecutar **MolForge** desde **Ubuntu (WSL)** con el entorno `MolForge_env` creado a partir del YAML oficial, y usar el entorno de utilidades `molforge-tools`.

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
│  ├─ sp/                           # vocabulari (importat del repo de MolForge)
│  └─ MolForge_output/              # resultados de MolForge
├─ saved_models/                    # Checkpoints del repo de MolForge (descarregar a banda)
├─ scripts/
│  ├─ smiles_to_fps.py              # convierte SMILES → fingerprints (CPU)
│  └─ run_molforge.py               # ejecuta MolForge (CPU) fila a fila y guarda resultats
├─ notebooks/
│  ├─ 01_smiles_to_fps.ipynb
│  └─ 02_run_molforge_cpu.ipynb
└─ .gitignore
```

---

## 🧩 Requisitos

1) **WSL2 + Ubuntu 22.04** instalados en Windows (ver guía abajo en caso de estar en PC Windows).  
2) **Conda/Miniconda** instalado dentro de **Ubuntu** (si usas la opción Conda).  
3) **Environment oficial de MolForge**: (No hace falta hacerlo porque ya está importado en este repositorio)
   - Copia el `environment.yml` **del repo oficial de MolForge** a `envs/molforge/environment.yml`.
   - **Añade en la sección `- pip:` la instalación del paquete**, por ejemplo:
     ```yaml
     - "MolForge @ git+https://github.com/knu-lcbc/MolForge.git"
     ```
   - Este proyecto asume **CPU** (no CUDA).

> **Alternativa sin Conda:** puedes usar **`venv`** (ver más abajo) tanto para el entorno de MolForge como para las utilidades.

---

## 🐧 Instalar Ubuntu (WSL) por primera vez

**En PowerShell (Administrador):**
```powershell
wsl --install -d Ubuntu-22.04
wsl --update
wsl -l -v      # debe mostrar Ubuntu con VERSION 2
```

**Rutas del proyecto según máquina (adaptar a la ruta personal de cada uno):**  
- **PC (WSL):** `D:\MolForge_Testing` (Windows) ↔ `/mnt/d/MolForge_Testing` (Ubuntu)  
- **Laboratorio (Linux):** `/export/home/ddiestre/MolForge_Testing`

Abrir Ubuntu ya dentro del proyecto:
```powershell
# PC (WSL):
wsl -d Ubuntu --cd /mnt/d/MolForge_Testing
```
```bash
# Laboratorio (Linux):
cd /export/home/ddiestre/MolForge_Testing
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

## 🛠️ Entornos (Conda)

### A) Entorno **MolForge** (Ubuntu/WSL, CPU)

> Usaremos el `environment.yml` oficial **intacto**, solo añadiendo la línea `pip` para instalar MolForge (ver Requisitos).

```bash
# PC (WSL) — o adapta a la ruta del Lab si trabajas allí
cd /mnt/d/MolForge_Testing          # PC (WSL)
# cd /export/home/ddiestre/MolForge_Testing   # Laboratorio (Linux)

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

**Ubuntu (WSL o Lab):**
```bash
# PC (WSL):
cd /mnt/d/MolForge_Testing
# o Laboratorio:
# cd /export/home/ddiestre/MolForge_Testing

conda env create -f envs/tools/environment.yml
conda activate molforge-tools
```

**Actualizar el entorno de tools cuando cambies el YAML:**
```bash
conda env update -f envs/tools/environment.yml --prune
```

---

## 🛠️ Entornos (venv) — Alternativa sin Conda

> Útil si prefieres `pip` puro o no puedes usar Conda.  
> **Nota RDKit:** en `pip` suele usarse el wheel `rdkit-pypi` (no oficial conda-forge). Si RDKit te da problemas, usa la opción Conda.

### A) Entorno **MolForge** con venv (CPU)

```bash
# PC (WSL):
cd /mnt/d/MolForge_Testing
# Laboratorio (Linux):
# cd /export/home/ddiestre/MolForge_Testing

# Crear y activar el venv (nombre sugerido: .venv_mf)
python3 -m venv .venv_mf
source .venv_mf/bin/activate

# Actualizar pip y herramientas básicas
python -m pip install --upgrade pip setuptools wheel

# Instalar PyTorch CPU (ajusta si tu sistema requiere índice específico)
pip install --index-url https://download.pytorch.org/whl/cpu torch torchvision torchaudio

# Instalar MolForge (desde el repo oficial)
pip install "MolForge @ git+https://github.com/knu-lcbc/MolForge.git"

# Dependencias comunes que suelen requerirse
pip install pandas pyarrow tqdm sentencepiece selfies rich

# Forzar modo CPU
export CUDA_VISIBLE_DEVICES=-1
```

### B) Entorno **molforge-tools** con venv

```bash
# PC (WSL):
cd /mnt/d/MolForge_Testing
# Laboratorio (Linux):
# cd /export/home/ddiestre/MolForge_Testing

# Crear y activar el venv (nombre sugerido: .venv_tools)
python3 -m venv .venv_tools
source .venv_tools/bin/activate

python -m pip install --upgrade pip setuptools wheel

# Paquetes mínimos de utilidades
pip install pandas pyarrow tqdm "rdkit-pypi>=2022.9.5"
```

> Para cambiar de un entorno a otro: `deactivate` y vuelve a `conda activate ...` o `source .venv_.../bin/activate` según corresponda.

---

## 🔁 Flujo de trabajo (solo CPU)

### 1) SMILES → Fingerprints (RDKit)

**Opción Notebook (recomendada la primera vez):**  
Abre `notebooks/01_smiles_to_fps.ipynb` y sigue las celdas.  
Entrada: CSV/Parquet con columna `smiles` (opcional `id`).  
Salida: fichero en `data/MolForge_input/` con columnas `id`, `smiles`, `fp_0000...`.

**Opción Script (rápido/automatizable):**
```bash
# Activa el entorno de tools (Conda o venv)
# conda activate molforge-tools
# source .venv_tools/bin/activate

python scripts/smiles_to_fps.py \
  --input data/SMILES/molecules.csv \
  --smiles-col smiles \
  --fp morgan --radius 2 --nBits 2048 \
  --output data/MolForge_input/morgan_2048.parquet
```

### 2) Fingerprints → MolForge (CPU)

**Opción Notebook:**  
Abre `notebooks/02_run_molforge_cpu.ipynb` y define:
- `fps_path` → fichero de `data/MolForge_input/`
- `checkpoint_path` → ruta a tu `.pth` (local, **no versionado**)
- `fp_name` → p. ej. `ECFP4`
- `model_type` → `smiles` (o `selfies`)
- `decode` → `greedy` (o `beam` si tu repo lo soporta)

**Opción Script (CPU):**
```bash
# Activa el entorno de MolForge (Conda o venv)
# conda activate MolForge_env
# source .venv_mf/bin/activate
export CUDA_VISIBLE_DEVICES=-1

python scripts/run_molforge.py \
  --fps data/MolForge_input/morgan_2048.parquet \
  --checkpoint /ruta/a/tu/checkpoint.pth \
  --fp-name ECFP4 \
  --model-type smiles \
  --decode greedy \
  --out data/MolForge_output/molforge_outputs.parquet
```

---

## ✅ Comprobaciones de instalación

**MolForge (Conda o venv):**
```bash
# conda activate MolForge_env   # o: source .venv_mf/bin/activate
export CUDA_VISIBLE_DEVICES=-1
python - << 'PY'
import torch
print("cuda available?:", torch.cuda.is_available())  # esperado False
from MolForge import main as _mf
print("MolForge import OK (MolForge)")
PY
```

**molforge-tools (Conda o venv):**
```bash
# conda activate molforge-tools  # o: source .venv_tools/bin/activate
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

## 📂 Gestión de datos

- `data/SMILES/` → entradas con SMILES (`.csv`/`.parquet`).  
- `data/MolForge_input/` → fingerprints generados (input a MolForge).  
- `data/MolForge_output/` → resultados de MolForge (SMILES generados, logs).  

**Sugerencia de nombres:** incluye el tipo/size del FP (`morgan_2048.parquet`) y fecha/modo en la salida (`molforge_outputs_YYYYMMDD.parquet`).

---

## 🔄 Actualizar entornos

- **tools (Conda):**
  ```bash
  conda env update -f envs/tools/environment.yml --prune
  ```
- **MolForge_env (Conda):** si cambia el YAML oficial, lo más limpio es recrear:
  ```bash
  conda env remove -n MolForge_env
  conda env create -f envs/molforge/environment.yml -n MolForge_env
  ```
- **venv (MolForge o tools):** para actualizar, vuelve a activar el venv y usa `pip install -U paquete` o reejecuta la lista de `pip install` correspondiente.

---

## 📝 Notas y buenas prácticas

- **`saved_models/`** existe en el repo pero su contenido (pesos) **no se versiona**; deja un `.gitkeep` como marcador.
- Forzar CPU: `export CUDA_VISIBLE_DEVICES=-1` antes de ejecutar.
- Si RDKit con `pip` falla, usa la opción **Conda** para las herramientas.
- Mantén separados los entornos de **MolForge** y de **tools** para evitar conflictos de dependencias.
