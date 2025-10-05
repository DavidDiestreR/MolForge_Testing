#!/usr/bin/env python
"""
Ejecuta MolForge en **CPU** por cada fila de un fichero de fingerprints y guarda los SMILES predichos.

Entradas soportadas en --fps:
  • Parquet/CSV/TSV/TXT con una columna 'indices' que contiene "i j k ..." (espacios).
  • Parquet/CSV/TSV/TXT de 1 sola columna (sin cabecera): se interpreta como 'indices'.
  • Parquet/CSV/TSV con columnas 'fp_0000'..'fp_N' (0/1): se convierten a 'indices'.

Uso (ejemplos):
  conda activate MolForge_env
  export CUDA_VISIBLE_DEVICES=-1
  python scripts/run_molforge.py \
    --fps data/MolForge_Input/test_1.csv \
    --checkpoint ECFP4_smiles_checkpoint.pth \
    --fp-name ECFP4 \
    --model-type smiles \
    --decode greedy \
    --out data/MolForge_output/molforge_outputs_demo.csv
"""

import argparse
import os
import sys
import io
from contextlib import redirect_stdout

import pandas as pd

# Forzar CPU explícitamente (sin selector GPU/CPU)
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

# Import MolForge (ajusta si tu entrypoint es distinto)
from MolForge.predict import main as molforge_main


# ---------- Utilidades de E/S ----------

def _read_any_table(path: str) -> pd.DataFrame:
    """Lee .parquet/.csv/.tsv/.txt con separador inferido para .txt."""
    ext = os.path.splitext(path)[1].lower()
    if ext == ".parquet":
        return pd.read_parquet(path)
    if ext == ".tsv":
        return pd.read_csv(path, sep="\t")
    if ext == ".csv":
        return pd.read_csv(path, sep=",")
    if ext == ".txt":
        # engine='python' evita el ParserWarning cuando sep=None
        return pd.read_csv(path, sep=None, engine="python", header="infer")
    # Si no hay extensión, intenta CSV con inferencia
    return pd.read_csv(path, sep=None, engine="python", header="infer")


def _active_indices(row_bits) -> str:
    """Convierte una fila de bits 0/1 (fp_0000..fp_N) en 'i j k'."""
    return " ".join(str(i) for i, v in enumerate(row_bits) if int(v) == 1)


def read_fps_as_indices_df(path: str) -> pd.DataFrame:
    """
    Normaliza el fichero de entrada a un DataFrame con:
      - 'id' (opcional; si no existe, se autogenera como '0','1',...)
      - 'indices' (OBLIGATORIA): cadena con 'i j k ...'
    """
    df = _read_any_table(path)

    # Caso 1: ya viene una columna 'indices'
    if "indices" in df.columns:
        out = pd.DataFrame()
        out["id"] = df["id"].astype(str) if "id" in df.columns else df.index.astype(str)
        out["indices"] = df["indices"].astype(str).fillna("")
        return out

    # Caso 2: un solo campo sin cabecera (o con cabecera distinta): tómalo como 'indices'
    if df.shape[1] == 1:
        col = df.columns[0]
        out = pd.DataFrame()
        out["id"] = df.index.astype(str)
        out["indices"] = df[col].astype(str).fillna("")
        return out

    # Caso 3: formato denso con fp_0000..fp_N -> convertir a 'indices'
    bit_cols = [c for c in df.columns if c.startswith("fp_")]
    if bit_cols:
        bit_cols = sorted(bit_cols, key=lambda c: int(c.split("_")[1]))
        out = pd.DataFrame()
        out["id"] = df["id"].astype(str) if "id" in df.columns else df.index.astype(str)
        out["indices"] = [
            _active_indices(df.loc[i, bit_cols]) for i in range(len(df))
        ]
        return out

    raise SystemExit(
        "Formato no soportado: aporta una columna 'indices', un fichero de una sola columna, "
        "o columnas 'fp_0000..fp_N' con 0/1."
    )


# ---------- Ejecución MolForge ----------

def extract_result_from_stdout(stdout_str: str):
    """Extrae el texto tras 'Result:' (sin espacios); None si no aparece."""
    for line in stdout_str.splitlines():
        if line.strip().startswith("Result:"):
            return line.split("Result:", 1)[1].strip().replace(" ", "")
    return None


def run_for_fp(indices_str: str, fp_name: str, model_type: str, checkpoint: str, decode: str):
    """Ejecuta MolForge (simulando CLI), captura stdout y devuelve (resultado, stdout)."""
    argv_backup = sys.argv
    sys.argv = [
        "predict.py",
        f"--fp={fp_name}",
        f"--model_type={model_type}",
        f"--input={indices_str}",
        f"--checkpoint={checkpoint}",
        f"--decode={decode}",
    ]
    buf = io.StringIO()
    try:
        with redirect_stdout(buf):
            molforge_main()
    finally:
        sys.argv = argv_backup

    out = buf.getvalue()
    smiles = extract_result_from_stdout(out)
    return smiles, out


# ---------- CLI ----------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fps", required=True,
                    help=("Fichero Parquet/CSV/TSV/TXT con:\n"
                          "  - columna 'indices' (\"i j k ...\"), o\n"
                          "  - una sola columna (se toma como 'indices'), o\n"
                          "  - columnas 'fp_0000..fp_N' con 0/1 (se convierten)."))
    ap.add_argument("--checkpoint", required=True, help="Ruta al checkpoint de MolForge (.pth)")
    ap.add_argument("--fp-name", default="ECFP4", help="Nombre del fingerprint para MolForge (ej. ECFP4)")
    ap.add_argument("--model-type", default="smiles", choices=["smiles", "selfies"], help="Representación destino")
    ap.add_argument("--decode", default="greedy", choices=["greedy", "beam"], help="Algoritmo de decodificación")
    ap.add_argument("--out", required=True, help="Salida Parquet/CSV/TSV/TXT")
    args = ap.parse_args()

    df_idx = read_fps_as_indices_df(args.fps)  # -> columnas: id, indices

    results = []
    for i, (idx, indices_str) in enumerate(zip(df_idx["id"], df_idx["indices"]), start=1):
        if not str(indices_str).strip():
            results.append((idx, None, ""))
            print(f"[{i}/{len(df_idx)}] id={idx} -> (sin índices)")
            continue

        smiles, raw = run_for_fp(str(indices_str), args.fp_name, args.model_type, args.checkpoint, args.decode)
        results.append((idx, smiles, raw))
        print(f"[{i}/{len(df_idx)}] id={idx} -> {smiles}")

    out_df = pd.DataFrame(results, columns=["id", "molforge_smiles", "raw_stdout"])

    # Guardar
    out_ext = os.path.splitext(args.out)[1].lower()
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    if out_ext == ".parquet":
        out_df.to_parquet(args.out, index=False)
    elif out_ext in (".csv", ".tsv", ".txt"):
        sep = "\t" if out_ext == ".tsv" else ","
        out_df.to_csv(args.out, index=False, sep=sep)
    else:
        raise ValueError(f"Extensión de salida no soportada: {out_ext}")

    print(f"[OK] Guardado: {args.out} | filas={len(out_df)}")


if __name__ == "__main__":
    main()
