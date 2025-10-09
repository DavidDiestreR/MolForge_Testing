#!/usr/bin/env python
"""
SMILES → fingerprints (CPU, RDKit).

Uso:
  conda activate molforge-tools
  python scripts/smiles_to_fps.py \
    --input data/SMILES/molecules.csv \
    --smiles-col smiles \
    --fp morgan --radius 2 --nBits 2048 \
    --output data/MolForge_input/morgan_2048.parquet
"""

import argparse
import os
import pandas as pd
import numpy as np

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, MACCSkeys, rdMolDescriptors as rdDesc


def to_numpy(bitvect):
    n_bits = bitvect.GetNumBits()
    arr = np.zeros((n_bits,), dtype=np.uint8)
    DataStructs.ConvertToNumpyArray(bitvect, arr)
    return arr


def fp_from_mol(mol, kind: str, nBits: int, radius: int):
    kind = kind.lower()
    if kind == "morgan":
        bv = AllChem.GetMorganFingerprintAsBitVect(mol, radius=radius, nBits=nBits)
    elif kind == "maccs":
        bv = MACCSkeys.GenMACCSKeys(mol)  # 167 bits
    elif kind == "rdkit":
        bv = Chem.RDKFingerprint(mol, fpSize=nBits)
    elif kind == "atompair":
        bv = rdDesc.GetHashedAtomPairFingerprintAsBitVect(mol, nBits=nBits)
    elif kind == "tt":
        bv = rdDesc.GetHashedTopologicalTorsionFingerprintAsBitVect(mol, nBits=nBits)
    else:
        raise ValueError(f"Fingerprint desconocido: {kind}")
    return to_numpy(bv)


def read_table(path: str, smiles_col: str):
    ext = os.path.splitext(path)[1].lower()
    if ext == ".parquet":
        df = pd.read_parquet(path)
    elif ext in (".csv", ".tsv", ".txt"):
        sep = "\t" if ext == ".tsv" else None  # autodetect CSV, tab en .tsv
        df = pd.read_csv(path, sep=sep)
    else:
        raise ValueError(f"Extensión no soportada: {ext}")
    if smiles_col not in df.columns:
        raise ValueError(f"No existe la columna SMILES: {smiles_col}")
    return df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="CSV/Parquet con columna de SMILES")
    ap.add_argument("--smiles-col", default="smiles", help="Nombre de la columna SMILES")
    ap.add_argument("--id-col", default=None, help="Columna con IDs (opcional)")
    ap.add_argument("--fp", required=True, choices=["morgan", "maccs", "rdkit", "atompair", "tt"])
    ap.add_argument("--radius", type=int, default=2, help="Solo para morgan (e.g., 2≈ECFP4)")
    ap.add_argument("--nBits", type=int, default=2048, help="Tamaño hash (ignorado en MACCS)")
    ap.add_argument("--output", required=True, help="Salida (.parquet o .csv)")
    args = ap.parse_args()

    df_in = read_table(args.input, args.smiles_col)

    # id por defecto: índice
    if args.id_col and args.id_col in df_in.columns:
        ids = df_in[args.id_col].astype(str).tolist()
    else:
        ids = df_in.index.astype(str).tolist()

    smiles_list = df_in[args.smiles_col].astype(str).tolist()

    fps_list = []
    keep_ids = []
    keep_smiles = []

    for _id, smi in zip(ids, smiles_list):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue  # minimal: saltar inválidos
        arr = fp_from_mol(mol, args.fp, args.nBits, args.radius)
        fps_list.append(arr)
        keep_ids.append(_id)
        keep_smiles.append(smi)

    if not fps_list:
        raise SystemExit("No se generó ningún fingerprint (¿SMILES inválidos?)")

    # DataFrame con columnas fp_0000..
    nbits = len(fps_list[0])
    col_names = [f"fp_{i:04d}" for i in range(nbits)]
    df_fp = pd.DataFrame(np.vstack(fps_list), columns=col_names)
    df_fp.insert(0, "smiles", keep_smiles)
    df_fp.insert(0, "id", keep_ids)

    out_ext = os.path.splitext(args.output)[1].lower()
    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    if out_ext == ".parquet":
        df_fp.to_parquet(args.output, index=False)
    elif out_ext in (".csv", ".tsv", ".txt"):
        sep = "\t" if out_ext == ".tsv" else ","
        df_fp.to_csv(args.output, index=False, sep=sep)
    else:
        raise ValueError(f"Extensión de salida no soportada: {out_ext}")

    print(f"[OK] Guardado: {args.output} | filas={len(df_fp)} | bits={nbits}")


if __name__ == "__main__":
    main()
