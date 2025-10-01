# MolForge Testing — **CPU-only** (Conda) · WSL + Windows

Este repo está preparado para **inferir en CPU**, sin usar GPU/CUDA, respetando el `environment.yml` **original** de MolForge.
La idea: ejecutar **MolForge** desde **Ubuntu (WSL)** con el entorno `MolForge_env` creado a partir del YAML oficial, y usar el entorno
de utilidades (`molforge-tools`) con **Conda** (puede ser en Windows o en Ubuntu).

---

## Estructura

```
MolForge_Testing/
├─ envs/
│  ├─ molforge/environment.yml      # environment oficial de MolForge (intacto)
│  └─ tools/environment.yml         # RDKit + pandas (ligero)
├─ data/
│  ├─ SMILES/                       # entradas con SMILES
│  ├─ MolForge_input/               # fingerprints generados (input para MolForge)
│  └─ MolForge_output/              # resultados de MolForge
├─ scripts/
│  ├─ smiles_to_fps.py              # convierte SMILES → fingerprints (CPU)
│  └─ run_molforge.py               # ejecuta MolForge (CPU) fila a fila y guarda resultados
├─ notebooks/
│  ├─ 01_smiles_to_fps.ipynb        # versión notebook con explicación pasito a pasito
│  └─ 02_run_molforge_cpu.ipynb     # idem para inferencia con MolForge (CPU)
└─ .gitignore
```

---

## 1) Instalar Ubuntu (WSL) por primera vez (solo si aún no lo tienes)

### 1.1 Instalar WSL y Ubuntu
Abre **PowerShell (Administrador)** y ejecuta:
```powershell
wsl --install -d Ubuntu-22.04
```
Reinicia si te lo pide. Abre **Ubuntu** desde el menú Inicio y crea usuario/contraseña.

### 1.2 (Opcional) Actualizar WSL
```powershell
wsl --update
```

### 1.3 Confirmar que estás en WSL2
```powershell
wsl -l -v
```
Debe poner: `Ubuntu ... VERSION 2`.

### 1.4 Diferencias de rutas
- En **Windows/PowerShell**: `D:\MolForge_Testing`
- En **Ubuntu/WSL**: `/mnt/d/MolForge_Testing`

Para abrir Ubuntu ya dentro del proyecto:
```powershell
wsl -d Ubuntu --cd /mnt/d/MolForge_Testing
```

---

## 2) Instalar **Conda en Ubuntu** (si aún no la tienes)

En la terminal de **Ubuntu**:
```bash
# instalación no interactiva rápida
curl -fsSL -o /tmp/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash /tmp/miniconda.sh -b -p "$HOME/miniconda3"
"$HOME/miniconda3/bin/conda" init bash
exec bash

# comprobar
conda --version
```

---

## 3) Entornos

### 3.1 Entorno de **MolForge** (Ubuntu/WSL, CPU)
> Usaremos el `environment.yml` **oficial** de MolForge **sin modificarlo**.

```bash
# dentro de Ubuntu, en la carpeta del proyecto
cd /mnt/d/MolForge_Testing

conda env create -f envs/molforge/environment.yml -n MolForge_env
conda activate MolForge_env

# Forzar CPU (sin GPU) en esta sesión — el script y el notebook también lo hacen
export CUDA_VISIBLE_DEVICES=-1

# verificación rápida
python - << 'PY'
import torch
print("cuda available?:", torch.cuda.is_available())  # debería imprimir False
print("device:", "cpu")
PY
```

### 3.2 Entorno de **utilidades (RDKit)** — `molforge-tools`
Puedes crearlo y usarlo **en Windows** (Miniconda para Windows) o **en Ubuntu**. El YAML es el mismo.

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

---

## 4) Flujo de trabajo (CPU)

### 4.1 SMILES → Fingerprints (RDKit)
**Opción A: Notebook (recomendado para entender el proceso)**
- Abre `notebooks/01_smiles_to_fps.ipynb` y sigue las celdas.
- Parámetros: `--fp` (`morgan|maccs|rdkit|atompair|tt`), `radius`, `nBits`, etc.
- Entrada: `data/SMILES/*.csv` o `.parquet` con columna `smiles` (y opcional `id`).
- Salida: `data/MolForge_input/*.parquet` (o `.csv`).

**Opción B: Script CLI**
```bash
conda activate molforge-tools
python scripts/smiles_to_fps.py   --input data/SMILES/molecules.csv   --smiles-col smiles   --fp morgan --radius 2 --nBits 2048   --output data/MolForge_input/morgan_2048.parquet
```

### 4.2 Fingerprints → MolForge (CPU)
**Opción A: Notebook**
- Abre `notebooks/02_run_molforge_cpu.ipynb`, indica:
  - `fps_path`: fichero en `data/MolForge_input/`
  - `checkpoint_path`: ruta al `.pth` del modelo entrenado
  - `fp_name`: p. ej. `ECFP4`
  - `model_type`: `smiles` (o `selfies` si aplica)
  - `decode`: `greedy` o `beam` (si el repo lo soporta)
- Ejecuta las celdas y generará `data/MolForge_output/molforge_outputs.parquet`.

**Opción B: Script CLI (CPU, sin selector GPU)**
```bash
conda activate MolForge_env
export CUDA_VISIBLE_DEVICES=-1

python scripts/run_molforge.py   --fps data/MolForge_input/morgan_2048.parquet   --checkpoint /ruta/a/tu/checkpoint.pth   --fp-name ECFP4   --model-type smiles   --decode greedy   --out data/MolForge_output/molforge_outputs.parquet
```

---

## 5) Consejos de uso (WSL vs Windows)

- Trabajar en `/mnt/d/...` desde Ubuntu es cómodo, pero **algo más lento** en I/O que trabajar dentro de `~/` (Linux). Para repos grandes, puedes clonar también en `~/proyectos/`.
- Los notebooks puedes abrirlos con **VSCode + Remote WSL** o con **Jupyter** dentro de Ubuntu:
  ```bash
  pip install jupyterlab
  jupyter lab
  ```
  y abres `http://127.0.0.1:8888` desde Windows.

---

## 6) Problemas comunes

- **`ModuleNotFoundError: MolForge/molforge`** → Asegúrate de haber instalado MolForge según su repo (el `environment.yml` oficial suele añadir la instalación por `pip`).  
- **MolForge lento** → Es normal en CPU. Para lotes grandes, considera paralelizar o ajustar batch si el repo lo permite.  
- **Rutas** → En Windows `D:\MolForge_Testing\...`, en Ubuntu `/mnt/d/MolForge_Testing/...`.
