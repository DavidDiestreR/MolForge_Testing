# MolForge Testing — **CPU** (Conda) · WSL + Windows (PC) o Linux (Lab)

Este repositorio está configurado para ejecutar **MolForge exclusivamente en CPU** usando entornos **Conda**.  
El entorno oficial de MolForge (`MolForge_env`) está basado en el `environment.yml` original del proyecto y está pensado **solo para inferencia en CPU**.  
Las tareas de **preprocesado** (como convertir SMILES → fingerprints) se realizan con el entorno auxiliar `molforge-tools`. 

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
│  ├─ smiles_to_fps.py              # convierte SMILES → fingerprints
│  └─ run_molforge.py               # ejecuta MolForge fila a fila y guarda resultats
├─ notebooks/
│  ├─ 01_smiles_to_fps.ipynb
│  └─ 02_run_molforge_cpu.ipynb
└─ .gitignore
```

---

## 🧩 Requisitos

1) **WSL2 + Ubuntu 22.04** instalados en Windows (ver guía abajo en caso de estar en PC Windows).  
2) **Conda/Miniconda** instalado.  
3) **Environment oficial de MolForge**: (No hace falta hacerlo porque ya está importado en este repositorio)
   - Copia el `environment.yml` **del repo oficial de MolForge** a `envs/molforge/environment.yml`.
   - **Añade en la sección `- pip:` la instalación del paquete**:
     ```yaml
     - "MolForge @ git+https://github.com/knu-lcbc/MolForge.git"
     ```
---

## 🐧 Instalar Ubuntu (WSL) por primera vez

**En PowerShell (Administrador):**
```powershell
wsl --install -d Ubuntu-22.04
wsl --update
wsl -l -v      # debe mostrar Ubuntu con VERSION 2
```

**Abrir Ubuntu ya dentro del proyecto:**
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

## 📦 Instalar Conda en Windows

Instalar **[Miniconda](https://repo.anaconda.com/miniconda/Miniconda3-latest-Windows-x86_64.exe)** (Selecciona "Just Me" y deja la ruta por defecto).

Una vez instalado MiniConda por primera vez abrir AnacondaPrompt y escribir "conda init bash" (De esta forma se podrá usar conda desde Git Bash)

```bash
conda init bash
```

---

## 🛠️ Entornos (Conda)

### A) Entorno **MolForge**

> Usaremos el `environment.yml` oficial (incluye PyTorch y, si se desea, soporte CUDA).

```bash
# PC (WSL) — o adapta a la ruta del Lab si trabajas allí
cd /mnt/d/MolForge_Testing          # PC (WSL)
# cd /export/home/ddiestre/MolForge_Testing   # Laboratorio (Linux)

conda env create -f envs/molforge/environment.yml -n MolForge_env
conda activate MolForge_env
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

## 🧠 Checkpoints y carpeta `saved_models`

La carpeta `saved_models/` está **vacía** en este repositorio porque los **checkpoints** de MolForge son demasiado pesados para versionarlos en GitHub.

Para usar los checkpoints del modelo ya entrenado, descárgalos del repo oficial de MolForge y colócalos en `saved_models/` (manteniendo los nombres de archivo):

- [**top-performing**](https://drive.google.com/uc?id=1zl6HBdwYsnA4JcnOi1o6OmcrRDB5iySK) — modelos con mejor rendimiento (recomendado).
- [**all the other models**](https://drive.google.com/uc?id=1jCtbc9lMacCyiZ3iZFEtFgOfOQYtWEuD) — el resto de modelos disponibles.

> Tras descargar, descomprime y verifica que el checkpoint que usarás existe en `saved_models/`.

---

## 🔁 Flujo de trabajo

### 1) SMILES → Fingerprints (RDKit)

**Opción Notebook (recomendada la primera vez):**  
Abre `notebooks/01_smiles_to_fps.ipynb` y sigue las celdas.  
Entrada: CSV/Parquet con columna `smiles` (opcional `id`).  
Salida: fichero en `data/MolForge_input/` con columnas `id`, `smiles`, `fp_0000...`.

**Opción Script (rápido/automatizable):**
```bash
conda activate molforge-tools

python scripts/smiles_to_fps.py \
  --input data/SMILES/molecules.csv \
  --smiles-col smiles \
  --fp morgan --radius 2 --nBits 2048 \
  --output data/MolForge_input/morgan_2048.parquet
```

### 2) Fingerprints → MolForge

**Opción Notebook:**  
Abre `notebooks/02_run_molforge_cpu.ipynb` y define:
- `fps_path` → fichero de `data/MolForge_input/`
- `checkpoint_path` → ruta a tu `.pth` (local, **no versionado**)
- `fp_name` → p. ej. `ECFP4`
- `model_type` → `smiles` (o `selfies`)
- `decode` → `greedy` (o `beam` si tu repo lo soporta)

**Opción Script:**
```bash
conda activate MolForge_env

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

**MolForge:**
```bash
conda activate MolForge_env
python - << 'PY'
import torch
print("torch:", torch.__version__)
print("build CUDA:", torch.version.cuda)            # None => build CPU-only
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

## 📂 Gestión de datos

- `data/SMILES/` → entradas con SMILES (`.csv`/`.parquet`).  
- `data/MolForge_Input/` → fingerprints generados (input a MolForge).  
- `data/MolForge_output/` → resultados de MolForge (SMILES generados, logs).  

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

---

## 📝 Notas y buenas prácticas

- **`saved_models/`** existe en el repo pero su contenido (pesos) **no se versiona**; deja un `.gitkeep` como marcador.
- Mantén separados los entornos de **MolForge** y de **tools** para evitar conflictos de dependencias.

---

## 📓 Ejecutar notebooks con el entorno Conda

Activa primero tu entorno `MolForge_env`/`molforge-tools` y ejecuta:

```bash
# En MolForge_env (línea por línea)
conda install -y ipykernel
conda install -y jupyterlab
conda install -y ipywidgets
python -m ipykernel install --user --name MolForge_env --display-name "Python (MolForge_env)"

# En molforge-tools (los paquetes ya van instalados en el yml)
python -m ipykernel install --user --name molforge-tools --display-name "Python (molforge-tools)"
```

Una vez instalado el entorno, ejecuta:
```bash
jupyter lab --no-browser --ip=0.0.0.0
```

> Abre la URL con token que imprime Jupyter en tu navegador de Windows. Dentro de JupyterLab, selecciona el kernel **Python (MolForge_env WSL)** para ejecutar los notebooks con ese entorno.

---

### 🚨 Limpieza si eliminas entornos
Quitar kernel “zombi” después de borrar un entorno:
```bash
jupyter kernelspec list
jupyter kernelspec uninstall MolForge_env -y
jupyter kernelspec uninstall molforge-tools -y
```

