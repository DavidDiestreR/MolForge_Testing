# MolForge Minimal Pipeline (Data layout personalizado)

Este repo mínimo te permite:
1) Convertir **SMILES → fingerprints** con RDKit.
2) Ejecutar **MolForge** desde esos fingerprints y guardar la **salida**.
3) Mantener **dos entornos** separados: `MolForge_env` (para MolForge) y `molforge-tools` (para RDKit/pandas).

---

## 📁 Estructura de carpetas

```
molforge-minimal/
├─ envs/
│  ├─ molforge/environment.yml      # usa el del repo oficial de MolForge (añade su línea pip)
│  └─ tools/environment.yml         # RDKit + pandas
├─ data/
│  ├─ SMILES/                       # aquí pones tus archivos de entrada con SMILES
│  ├─ MolForge_input/               # aquí guardaremos los fingerprints (input de MolForge)
│  └─ MolForge_output/              # aquí se guardará la salida de MolForge
├─ scripts/
│  ├─ smiles_to_fps.py              # convierte SMILES → fingerprints (Morgan/MACCS/RDK/AtomPair/TT)
│  └─ run_molforge.py               # ejecuta MolForge fila a fila y guarda resultados
└─ .gitignore
```

> **Nota:** He ajustado los comandos y rutas para usar exactamente estas tres carpetas dentro de `data/`:
> `SMILES/`, `MolForge_input/` y `MolForge_output/`.

---

## 🧩 Requisitos previos

- Tener **Conda o Mamba** instalado. (Recomiendo `mamba` por rapidez).
- Conseguir el/los **checkpoints** del modelo de MolForge que vayas a usar (archivo `.pth`). Guárdalos donde prefieras; en los ejemplos uso `data/models/...` pero puedes usar otra ruta.
- El **environment.yml de MolForge** (del repo oficial). En ese archivo, **añade** la línea `pip` que instala MolForge desde Git.

---

## ⚙️ Entornos: crear, activar, desactivar y actualizar

### 1) Entorno de herramientas (RDKit/pandas) → `molforge-tools`

**Crear (una sola vez):**
```bash
mamba env create -f envs/tools/environment.yml
```

**Activar:**
```bash
mamba activate molforge-tools
```

**Desactivar:**
```bash
mamba deactivate
```

**Actualizar (si editas `envs/tools/environment.yml`):**
```bash
mamba env update -f envs/tools/environment.yml --prune
```

### 2) Entorno de MolForge → `MolForge_env`

1. Copia el `environment.yml` **oficial** del repo de MolForge dentro de `envs/molforge/environment.yml`.
2. Dentro de ese `environment.yml`, en la sección `- pip:`, añade la línea para instalar MolForge desde Git, por ejemplo:
   ```yaml
   - "MolForge @ git+https://github.com/knu-lcbc/MolForge.git"
   ```

**Crear:**
```bash
mamba env create -f envs/molforge/environment.yml
```

**Activar / Desactivar:**
```bash
mamba activate MolForge_env
mamba deactivate
```

**Actualizar:**
```bash
mamba env update -f envs/molforge/environment.yml --prune
```

**Comprobar que está bien instalado:**
```bash
mamba activate MolForge_env
python -c "import torch, MolForge; print('torch cuda disponible?:', torch.cuda.is_available())"
```

---

## 🔁 Flujo de trabajo paso a paso (para tontos y con ejemplos)

### 0) Prepara un fichero con SMILES
- Crea un CSV con una columna llamada **`smiles`** (y opcionalmente `id`), por ejemplo:
  ```text
  id,smiles
  mol1,CCO
  mol2,c1ccccc1
  ```
- Guárdalo como: `data/SMILES/molecules.csv`

### A) **SMILES → Fingerprints** con RDKit (entorno `molforge-tools`)
Activa el entorno de tools y ejecuta:

```bash
mamba activate molforge-tools

python scripts/smiles_to_fps.py   --input data/SMILES/molecules.csv   --smiles-col smiles   --fp morgan --radius 2 --nBits 2048   --output data/MolForge_input/morgan_2048.parquet
```

- `--fp` admite: `morgan`, `maccs`, `rdkit`, `atompair`, `tt`.
- Para **morgan**: `radius=2` ≈ **ECFP4** (esto suele ser lo que los modelos llaman ECFP4).
- **Salida**: un archivo con columnas `id`, `smiles` y `fp_0000 ... fp_2047` (0/1).

### B) **Fingerprints → MolForge** (entorno `MolForge_env`)
Activa el entorno de MolForge y ejecuta:

```bash
mamba activate MolForge_env

python scripts/run_molforge.py   --fps data/MolForge_input/morgan_2048.parquet   --checkpoint /ruta/a/tu/checkpoint.pth   --fp-name ECFP4   --model-type smiles   --decode greedy   --out data/MolForge_output/molforge_outputs.parquet
```

- `--fp-name` debe coincidir con lo que espera tu modelo (p. ej., `ECFP4`).
- `--model-type` suele ser `smiles` (también soporta `selfies` si tu modelo lo usa).
- `--decode`: `greedy` o `beam` (si el repositorio soporta beam).

**¿Qué guarda?** Un fichero con columnas:
- `id`: copiado del input (o índice si no había `id`),  
- `molforge_smiles`: resultado devuelto por MolForge para ese fingerprint,  
- `raw_stdout`: el texto completo que imprimió MolForge (por si quieres auditar).

---

## 🗃️ Gestión de archivos de datos (qué va en cada carpeta)

- `data/SMILES/` → **tus entradas crudas** con SMILES (CSV/Parquet).  
- `data/MolForge_input/` → **fingerprints** que genera `smiles_to_fps.py`.  
- `data/MolForge_output/` → **resultados** de `run_molforge.py` (SMILES generados por MolForge).  

> Sugerencia de nombres: incluye el tipo de FP y tamaño, p.ej. `morgan_2048.parquet`, y el modo de decodificación/salida, p.ej. `molforge_outputs.parquet`.

---

## 🧯 Problemas comunes y soluciones rápidas

- **`ModuleNotFoundError: MolForge`**  
  → No está instalado en `MolForge_env`. Revisa el `environment.yml` de MolForge y la línea pip. Luego:
  ```bash
  mamba env update -f envs/molforge/environment.yml --prune
  ```

- **`torch.cuda.is_available()` devuelve `False`**  
  → Estás en CPU o tu PyTorch/CUDA no coincide con tu GPU/driver. Comprueba que el `environment.yml` de MolForge trae la build correcta de PyTorch para tu GPU.

- **El script de MolForge va lento**  
  → Este ejemplo lanza MolForge **una fila por vez** por simplicidad. Si necesitas velocidad, puedes adaptar `run_molforge.py` para agrupar entradas (si el repositorio de MolForge lo permite) o paralelizar.

---

## 🧪 Comandos rápidos (copy-paste)

```bash
# Crear entornos una vez
mamba env create -f envs/tools/environment.yml
mamba env create -f envs/molforge/environment.yml

# Generar FP
mamba activate molforge-tools
python scripts/smiles_to_fps.py --input data/SMILES/molecules.csv --smiles-col smiles --fp morgan --radius 2 --nBits 2048 --output data/MolForge_input/morgan_2048.parquet

# Ejecutar MolForge
mamba activate MolForge_env
python scripts/run_molforge.py --fps data/MolForge_input/morgan_2048.parquet --checkpoint /ruta/model.pth --fp-name ECFP4 --model-type smiles --decode greedy --out data/MolForge_output/molforge_outputs.parquet
```

¡Listo!
