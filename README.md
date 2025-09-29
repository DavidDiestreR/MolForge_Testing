# MolForge Minimal Pipeline — **Conda**, CUDA/GPU requerido

Este repo mínimo te permite:
1) Convertir **SMILES → fingerprints** con RDKit (entorno `molforge-tools`).
2) Ejecutar **MolForge** desde esos fingerprints y guardar la **salida** (entorno `MolForge_env` con **PyTorch CUDA**).

> Este proyecto **asume GPU** con CUDA: no incluye rutas de instalación ni ejemplos para CPU‑only.

---

## 📁 Estructura de carpetas

```
MolForge_Testing/
├─ envs/
│  ├─ molforge/environment.yml      # environment oficial de MolForge + instalador de MolForge
│  └─ tools/environment.yml         # RDKit + pandas (ligero, sin torch)
├─ data/
│  ├─ SMILES/                       # entradas con SMILES
│  ├─ MolForge_input/               # fingerprints generados (input para MolForge)
│  └─ MolForge_output/              # resultados de MolForge
├─ scripts/
│  ├─ analysis/.gitkeep
│  ├─ molforge/run_molforge.py      # ejecuta MolForge y guarda resultados
│  └─ tools/smiles_to_fps.py        # convierte SMILES → fingerprints         
└─ .gitignore
```

> **Windows (primera vez con Conda)**: ejecuta `conda init powershell` (o `conda init bash` si usas Git Bash), cierra y reabre la terminal.


---

## 🧩 Requisitos

- **Conda/Miniconda** instalado.
- **GPU NVIDIA con drivers adecuados** y soporte CUDA.
- **Checkpoints** de MolForge (archivo `.pth`) compatibles con tu fingerprint/representación.
- El `environment.yml` **oficial** del repo de MolForge copiado en `envs/molforge/environment.yml` y, dentro de `- pip:`, añade la línea para instalar el paquete, por ejemplo:
  ```yaml
  - "MolForge @ git+https://github.com/knu-lcbc/MolForge.git"
  ```
- Asegúrate de que ese `environment.yml` instala **PyTorch con CUDA** (no la build de CPU).

---

## ⚙️ Entornos con **Conda** (crear, activar, actualizar)

### 1) Entorno de herramientas → `molforge-tools` (RDKit/pandas)

Crear:
```bash
conda env create -f envs/tools/environment.yml
```

Activar / Desactivar:
```bash
conda activate molforge-tools
conda deactivate
```

Actualizar:
```bash
conda env update -f envs/tools/environment.yml --prune
```

### 2) Entorno de MolForge → `MolForge_env` (PyTorch **CUDA**)

Crear:
```bash
conda env create -f envs/molforge/environment.yml
```

Activar / Desactivar:
```bash
conda activate MolForge_env
conda deactivate
```

Actualizar:
```bash
conda env update -f envs/molforge/environment.yml --prune
```

**Verificación (CUDA obligatoria):**
```bash
conda activate MolForge_env
python - << 'PY'
import torch
assert torch.cuda.is_available(), "CUDA no disponible. Revisa tu instalación de PyTorch/CUDA y drivers NVIDIA."
print("CUDA OK | versión CUDA build:", torch.version.cuda, "| device 0:", torch.cuda.get_device_name(0))
try:
    from MolForge import main as _mf
    print("MolForge import OK")
except Exception:
    from molforge import main as _mf
    print("molforge import OK")
PY
```

---

## 🔁 Flujo de trabajo

### 0) Prepara un fichero con SMILES
Crea un CSV con columna **`smiles`** (y opcionalmente `id`), por ejemplo `data/SMILES/molecules.csv`:
```text
id,smiles
mol1,CCO
mol2,c1ccccc1
```

### A) **SMILES → Fingerprints** (entorno `molforge-tools`)
```bash
conda activate molforge-tools

python scripts/smiles_to_fps.py   --input data/SMILES/molecules.csv   --smiles-col smiles   --fp morgan --radius 2 --nBits 2048   --output data/MolForge_input/morgan_2048.parquet
```
- `--fp`: `morgan`, `maccs`, `rdkit`, `atompair`, `tt`.
- Para **morgan**: `radius=2` ≈ **ECFP4**.
- Salida: `id`, `smiles`, `fp_0000...` en `data/MolForge_input/` (Parquet o CSV).

### B) **Fingerprints → MolForge** (entorno `MolForge_env` con CUDA)
```bash
conda activate MolForge_env

python scripts/run_molforge.py   --fps data/MolForge_input/morgan_2048.parquet   --checkpoint /ruta/a/tu/checkpoint.pth   --fp-name ECFP4   --model-type smiles   --decode greedy   --out data/MolForge_output/molforge_outputs.parquet
```
- `--fp-name` debe **coincidir** con el modelo (p. ej., `ECFP4`).
- `--model-type`: `smiles` o `selfies` según el checkpoint.
- Salida: `data/MolForge_output/molforge_outputs.parquet` con `id`, `molforge_smiles`, `raw_stdout`.

---

## 🗃️ Gestión de datos

- `data/SMILES/` → entradas con SMILES.  
- `data/MolForge_input/` → fingerprints generados (input a MolForge).  
- `data/MolForge_output/` → resultados de MolForge.  

**Sugerencia:** incluye tipo/size del FP en el nombre (`morgan_2048.parquet`) y fecha/modo en la salida (`molforge_outputs_YYYYMMDD.parquet`).

---

## 🔄 Sincronizar entornos entre PCs

Tras `git pull`, actualiza los entornos para reflejar cambios en YAML:
```bash
conda env update -f envs/tools/environment.yml --prune
conda env update -f envs/molforge/environment.yml --prune
```

---

## 🧯 Problemas comunes

- **`AssertionError: CUDA no disponible`** → Estás usando la build de CPU o faltan drivers CUDA. Reinstala/actualiza el entorno de MolForge con la build CUDA correcta y comprueba `nvidia-smi`.  
- **`ModuleNotFoundError: MolForge/molforge`** → Revisa que el `pip` del YAML de MolForge incluya la línea para instalar el paquete desde Git y vuelve a actualizar el entorno.
