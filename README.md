# Proyecto de Evaluación con MolForge

Este repositorio contiene código para evaluar moléculas generadas con **[MolForge](https://github.com/knu-lcbc/MolForge)** usando **RDKit** y distintos fingerprints.  
Aquí encontrarás instrucciones claras para preparar el entorno y ejecutar el código, tanto en **CPU** como en **GPU (NVIDIA)**.

---

## ⚙️ Requisitos previos

- Tener instalado **[Conda](https://docs.conda.io/en/latest/miniconda.html)** (se recomienda Miniconda o Mambaforge).
- Una vez instalado MiniConda por primera vez abrir AnacondaPrompt y escribir "conda init bash"
- Conexión a internet para instalar los paquetes.
- (Opcional) GPU NVIDIA con drivers actualizados, si deseas acelerar el cálculo.

---

## 🖥️ Preparación del entorno

En este repositorio incluimos un archivo `environment.yml` **general** con todas las dependencias comunes.  
La única diferencia entre CPU y GPU está en cómo se instala **PyTorch**.

### 1. Crear entorno base (común a todos)
```bash
# Clonar este repositorio
git clone https://github.com/DavidDiestreR/MolForge_Testing.git
cd MolForge_Testing

# Crear el entorno base a partir de environment.yml
conda env create -f environment.yml -n env
conda activate env

# Comprovar que estamos en el entorno
python -c "import sys; print(sys.executable)"

# Eliminar un paquete del entorno
conda remove <nombre_paquete>

# eliminar el entorno ('env')
conda remove -n env --all

```

Este paso instala:
- Python 3.8
- RDKit
- Pandas, Numpy, Matplotlib, Seaborn
- TQDM, Rich
- SELFIES, SentencePiece, Gdown
- JupyterLab + ipykernel
- MolForge (instalado directamente desde GitHub)

---

### 2. Instalar PyTorch según tu máquina

#### 🔹 Opción CPU (portátiles o PCs sin GPU NVIDIA)
```bash
conda install -y pytorch -c pytorch -c conda-forge
```

#### 🔹 Opción GPU (PC con NVIDIA, recomendado en torre)
```bash
python -m pip install --upgrade pip
python -m pip install --index-url https://download.pytorch.org/whl/cu121 torch torchvision torchaudio
```

> ⚠️ Importante: no instales ambas variantes a la vez; usa **solo una** según tu hardware.

---

## ✅ Verificación

Comprueba que todo funciona:

```bash
python - << 'PY'
import torch, molforge, rdkit, pandas, matplotlib, seaborn, selfies, tqdm, sentencepiece, gdown, rich
print("env: OK")
print("torch:", torch.__version__, "| cuda avail:", torch.cuda.is_available())
if torch.cuda.is_available():
    print("device:", torch.cuda.get_device_name(0))
PY
```

- En CPU debería mostrar `cuda avail: False`.  
- En GPU debería mostrar `cuda avail: True`.

---

## 🚀 Uso del proyecto

### Opción A: fijar CPU (simple y universal)
Si quieres que tu código funcione en cualquier PC sin preocuparte de GPUs:
```python
device = "cpu"
```

Esto garantiza compatibilidad máxima (ejecutará todo en CPU).  

---

### Opción B: autodetección CPU/GPU (más flexible)
Si prefieres que el código use GPU cuando esté disponible:
```python
import torch
device = "cuda" if torch.cuda.is_available() else "cpu"
print("Using device:", device)
```

Esto selecciona GPU automáticamente si existe; si no, usará CPU.  

Puedes además permitir un **override manual** con una variable de entorno:
```bash
DEVICE=cpu python mi_script.py
```

---

## 🔄 Mantener el entorno actualizado entre varios PCs

Es habitual que instales un paquete nuevo en tu portátil y quieras tenerlo también en tu PC de torre (o viceversa).  
Para ello:

1. Instala el paquete en tu entorno actual:
   ```bash
   conda install -c conda-forge scikit-learn
   ```
   (o `pip install paquete` si es vía pip).

2. Edita tu `environment.yml` y añade el nuevo paquete en la sección `dependencies`.

3. Haz commit y push:
   ```bash
   git add environment.yml
   git commit -m "add scikit-learn"
   git push
   ```

4. En el otro PC, actualiza el entorno:
   ```bash
   git pull
   conda env update -n env -f environment.yml
   ```

👉 Esto instalará **solo los paquetes nuevos o actualizados**.  
**No reinstalará PyTorch** porque en este `environment.yml` no está incluido: lo instalamos a mano (CPU o GPU) en cada máquina.

---

## 📂 Estructura del repositorio

```
mi-proyecto/
├─ src/                    # código propio (funciones, evaluaciones, utils…)
├─ notebooks/              # experimentos Jupyter interactivos
├─ environment.yml         # entorno general (común a CPU y GPU)
├─ README.md               # este documento
└─ .gitignore
```

---

## 🔄 Buenas prácticas

- Mantén **un solo `environment.yml` general** en el repo.  
- Instala PyTorch aparte, según CPU o GPU, para no forzarlo en el YAML.  
- Si en el futuro añades paquetes propios, actualiza el `environment.yml` y haz commit.  
- En otros PCs, simplemente `git pull` + `conda env update` para sincronizar.

---
