"""
Mòdul per convertir SMILES en fingerprints de RDKit.

Fingerprints suportats (paper de MolForge):
- MACCS
- Avalon
- HashAP
- RDK4
- RDK4-L (RDK4_L, RDK4L)
- TT
- HashTT
- AEs
- ECFP0, ECFP2, ECFP4
- FCFP2, FCFP4
- ECFP2*, ECFP4*  (versions explícites en bits)

Comportament especial:
- Si l'entrada és np.nan -> retorna np.nan
- Si el SMILES no és vàlid -> retorna el string "InvalidSMILE"
"""

from rdkit import Chem
from rdkit.Chem import MACCSkeys, rdMolDescriptors, AllChem
import numpy as np
import src.fingerprints as fingerprints


def smiles_to_mol(smiles: str):
    """
    Converteix un SMILES en un objecte Mol de RDKit.

    - Si 'smiles' és np.nan -> retorna np.nan
    - Si no es pot parsejar el SMILES -> retorna np.nan
    """
    if smiles is np.nan:
        return np.nan

    # Ens assegurem que tractem un string
    smiles_str = str(smiles)

    # Si hi ha espais (p.ex. SMILES tokenitzats), els traiem
    smiles_str = smiles_str.replace(" ", "")

    mol = Chem.MolFromSmiles(smiles_str, sanitize=True) # sanitize retorna none si el mol no és vàlid
    if mol is None:
        return np.nan # retornem un np.nan si el mol no és vàlid
    return mol


def mol_to_fingerprint(mol, fp_type: str, n_bits: int = 2048, return_bits = True):
    """
    Converteix un Mol en el fingerprint especificat, seguint la implementació de MolForge.

    Paràmetres
    ----------
    mol : rdkit.Chem.Mol o np.nan
        Molècula RDKit (assumim que és vàlida). Si és np.nan -> retorna np.nan.
    fp_type : str
        Nom del fingerprint (veure llista al docstring del mòdul).
    n_bits : int
        Mida del fingerprint per als tipus "hashed" (ECFP*, FCFP, HashAP, HashTT, RDK4, RDK4-L).
        Per MACCS i Avalon s'ignora (166 i 512).

    Retorn
    ------
    Un objecte de fingerprint de RDKit (ExplicitBitVect o vector espars),
    np.nan si mol és np.nan,
    o bé llença ValueError si el tipus de fingerprint no està suportat.
    """
    if mol is np.nan:
        return np.nan

    key = fp_type

    # Fem servir les funcions definides a fingerprints.py
    # Mapegem el nom normalitzat (key) amb la funció corresponent

    if key == "MACCS":
        # A fingerprints.py es diu MAACS (mateix concepte)
        return fingerprints.MAACS(mol, return_bits=return_bits)

    elif key == "Avalon":
        return fingerprints.Avalon(mol, return_bits=return_bits)

    # --- Path-based fingerprints ---
    elif key == "RDK4":
        return fingerprints.RDK4(mol, return_bits=return_bits)

    elif key == "RDK4-L":
        # A fingerprints.py és RDK4_L
        return fingerprints.RDK4_L(mol, return_bits=return_bits)

    elif key == "HashAP":
        return fingerprints.HashAP(mol, return_bits=return_bits)

    # --- Topological torsions ---
    elif key == "TT":
        return fingerprints.TT(mol, return_bits=return_bits)

    elif key == "HashTT":
        return fingerprints.HashTT(mol, return_bits=return_bits)

    # --- Circular env (Morgan / AEs / ECFP / FCFP) ---
    elif key == "AEs":
        return fingerprints.AEs(mol, return_bits=return_bits)

    elif key == "ECFP0":
        return fingerprints.ECFP0(mol, return_bits=return_bits)

    elif key == "ECFP2":
        return fingerprints.ECFP2(mol, return_bits=return_bits)

    elif key == "ECFP4":
        return fingerprints.ECFP4(mol, return_bits=return_bits)

    # --- ECFP2* / ECFP4* = versions "bit vector" (equivalent a ECFP2_ / ECFP4_) ---
    elif key == "ECFP2*":
        # ECFP2* = versió explícita en bits; a fingerprints.py és ECFP2_
        return fingerprints.ECFP2_(mol, return_bits=return_bits)

    elif key == "ECFP4*":
        # ECFP4* = versió explícita en bits; a fingerprints.py és ECFP4_
        return fingerprints.ECFP4_(mol, return_bits=return_bits)

    # --- FCFP (feature-based Morgan) ---
    elif key == "FCFP2":
        return fingerprints.FCFP2(mol, return_bits=return_bits)

    elif key == "FCFP4":
        return fingerprints.FCFP4(mol, return_bits=return_bits)

    else:
        raise ValueError(
            f"Tipus de fingerprint no suportat: '{fp_type}'. "
            "Comprova que coincideix amb algun dels utilitzats al paper "
            "(MACCS, Avalon, HashAP, RDK4, RDK4-L, TT, HashTT, "
            "AEs, ECFP0/2/4, ECFP2*/4*, FCFP2/4)."
        )


def smiles_to_fingerprint(smiles, fp_type: str = "ECFP4", n_bits: int = 2048, return_bits=True):
    """
    Converteix un SMILES directament en fingerprint, gestionant np.nan i SMILES invàlids.

    - Si 'smiles' és np.nan -> retorna np.nan
    - Si el SMILES no és vàlid -> retorna el string "InvalidSMILE"
    - Si és vàlid -> retorna el fingerprint corresponent

    Aquest és el punt d'entrada recomanat si treballes amb un DataFrame.
    """
    # Cas especial: valor "missing"
    if smiles is np.nan:
        return np.nan

    mol = smiles_to_mol(smiles)

    if mol is np.nan:
        # SMILES rebutjat per RDKit
        return "InvalidSMILE"

    return mol_to_fingerprint(mol, fp_type=fp_type, n_bits=n_bits, return_bits=return_bits)


def get_supported_fingerprints():
    """
    Retorna la llista de noms de fingerprints suportats (per conveniència).
    """
    return [
        "MACCS",
        "Avalon",
        "RDK4",
        "RDK4-L",
        "HashAP",
        "TT",
        "HashTT",
        "ECFP0",
        "ECFP2",
        "ECFP4",
        "FCFP2",
        "FCFP4",
        "AEs",
        "ECFP2*",
        "ECFP4*",
    ]