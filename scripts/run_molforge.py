#!/usr/bin/env python
"""
Ejecuta MolForge en **CPU** por cada fila de un fichero de fingerprints y guarda los SMILES predichos.

Uso:
  conda activate MolForge_env
  export CUDA_VISIBLE_DEVICES=-1
  python scripts/run_molforge.py \
    --fps data/MolForge_input/morgan_2048.parquet \
    --checkpoint /ruta/model.pth \
    --fp-name ECFP4 \
    --model-type smiles \
    --decode greedy \
    --out data/MolForge_output/molforge_outputs.parquet
"""

import argparse
import os
import sys
import io
import pandas as pd
from contextlib import redirect_stdout

# Forzar CPU explícitamente (sin selector GPU/CPU)
os.environ["CUDA_VISIBLE_DEVICES"] = "-1"

# Import robusto por si el módulo se llama MolForge o molforge
try:
    from MolForge import main as molforge_main
except Exception:
    from molforge import main as molforge_main


def read_table(path: str):
    ext = os.path.splitext(path)[1].lower()
    if ext == ".parquet":
        return pd.read_parquet(path)
    elif ext in (".csv", ".tsv", ".txt"):
        sep = "\t" if ext == ".tsv" else None
        return pd.read_csv(path, sep=sep)
    else:
        raise ValueError(f"Extensión no soportada: {ext}")


def active_indices(row_bits):
    on = [str(i) for i, v in enumerate(row_bits) if int(v) == 1]
    return " ".join(on)


def extract_result_from_stdout(stdout_str: str):
    result = None
    for line in stdout_str.splitlines():
        if line.strip().startswith("Result:"):
            tokens = line.split("Result:", 1)[1].strip()
            result = tokens.replace(" ", "")
            break
    return result


def run_for_fp(indices_str: str, fp_name: str, model_type: str, checkpoint: str, decode: str):
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fps", required=True, help="Fichero de fingerprints (Parquet/CSV) con columnas fp_0000...")
    ap.add_argument("--checkpoint", required=True, help="Ruta al checkpoint de MolForge (.pth)")
    ap.add_argument("--fp-name", default="ECFP4", help="Nombre del fingerprint para MolForge (ej. ECFP4)")
    ap.add_argument("--model-type", default="smiles", choices=["smiles", "selfies"], help="Representación destino")
    ap.add_argument("--decode", default="greedy", choices=["greedy", "beam"], help="Algoritmo de decodificación")
    ap.add_argument("--out", required=True, help="Salida Parquet/CSV")
    args = ap.parse_args()

    df = read_table(args.fps)
    bit_cols = [c for c in df.columns if c.startswith("fp_")]
    if not bit_cols:
        raise SystemExit("No se encontraron columnas 'fp_XXXX' en el fichero de fingerprints.")
    bit_cols = sorted(bit_cols, key=lambda c: int(c.split("_")[1]))

    ids = df["id"].astype(str).tolist() if "id" in df.columns else [str(i) for i in range(len(df))]

    results = []
    for i, idx in enumerate(ids):
        row_bits = df.loc[df.index[i], bit_cols]
        indices_str = active_indices(row_bits)
        if not indices_str:
            results.append((idx, None, ""))
            continue

        smiles, raw = run_for_fp(indices_str, args.fp_name, args.model_type, args.checkpoint, args.decode)
        results.append((idx, smiles, raw))
        print(f"[{i+1}/{len(ids)}] id={idx} -> {smiles}")

    out_df = pd.DataFrame(results, columns=["id", "molforge_smiles", "raw_stdout"])
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
