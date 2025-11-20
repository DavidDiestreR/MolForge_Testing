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
- Si l'entrada és None -> retorna None
- Si el SMILES no és vàlid -> retorna el string "InvalidSMILE"
"""

from rdkit import Chem
from rdkit.Chem import MACCSkeys, rdMolDescriptors, rdFingerprintGenerator
from rdkit.Avalon import pyAvalonTools


def _normalize_fp_type(fp_type: str) -> str:
    """
    Normalitza el nom del fingerprint per facilitar la comparació.

    Exemples:
    - "ecfp4"     -> "ECFP4"
    - "RDK4-L"    -> "RDK4L"
    - "RDK4_L"    -> "RDK4L"
    - "  hashap " -> "HASHAP"
    - "ecfp4*"    -> "ECFP4*"
    """
    t = fp_type.strip().upper()
    t = t.replace(" ", "")
    # Unifiquem guions i guions baixos
    t = t.replace("-", "").replace("_", "")
    return t


def smiles_to_mol(smiles: str):
    """
    Converteix un SMILES en un objecte Mol de RDKit.

    - Si 'smiles' és None -> retorna None
    - Si no es pot parsejar el SMILES -> retorna None
    """
    if smiles is None:
        return None

    # Ens assegurem que tractem un string
    smiles_str = str(smiles)

    # Si hi ha espais (p.ex. SMILES tokenitzats), els traiem
    smiles_str = smiles_str.replace(" ", "")

    mol = Chem.MolFromSmiles(smiles_str)
    return mol


def mol_to_fingerprint(mol, fp_type: str, n_bits: int = 2048):
    """
    Converteix un Mol en el fingerprint especificat.

    Paràmetres
    ----------
    mol : rdkit.Chem.Mol o None
        Molècula RDKit (assumim que és vàlida). Si és None -> retorna None.
    fp_type : str
        Nom del fingerprint (veure llista al docstring del mòdul).
    n_bits : int
        Mida del fingerprint per als tipus "hashed" (ECFP*, FCFP, HashAP, HashTT, RDK4, RDK4-L).
        Per MACCS i Avalon s'ignora i es fa servir la mida estàndard (166 i 512).

    Retorn
    ------
    Un objecte de fingerprint de RDKit (ExplicitBitVect o vector espars),
    None si mol és None,
    o bé llença ValueError si el tipus de fingerprint no està suportat.
    """
    if mol is None:
        # Cas general: si l'entrada ja és None, propaguem None
        return None

    key = _normalize_fp_type(fp_type)

    # --- Predefined substructures ---
    if key == "MACCS":
        # MACCS: 166 bits basat en SMARTS
        return MACCSkeys.GenMACCSKeys(mol)

    # --- Paths and feature classes (Avalon) ---
    elif key == "AVALON":
        # Al paper fan servir 512 bits per Avalon
        return pyAvalonTools.GetAvalonFP(mol, nBits=512)

    # --- Path-based fingerprints (RDKit + atom pairs) ---
    elif key == "RDK4":
        # RDKit fingerprint: subgràfics ramificats i lineals de longitud 2–4
        return Chem.RDKFingerprint(
            mol,
            fpSize=n_bits,
            minPath=2,
            maxPath=4,
            branchedPaths=True,
        )

    elif key == "RDK4L":
        # Variant sense ramificació (només camins lineals)
        return Chem.RDKFingerprint(
            mol,
            fpSize=n_bits,
            minPath=2,
            maxPath=4,
            branchedPaths=False,
        )

    elif key == "HASHAP":
        # Atom pair fingerprint hashed a n_bits, amb distància 1–6 (paper)
        return rdMolDescriptors.GetHashedAtomPairFingerprintAsBitVect(
            mol,
            nBits=n_bits,
            minLength=1,
            maxLength=6,
        )

    # --- 4-atom paths (Topological torsion) ---
    elif key == "TT":
        # Topological torsion "sparse" (vector de comptes)
        # Equivalent al que anomenen TT al paper.
        return rdMolDescriptors.GetTopologicalTorsionFingerprint(mol)

    elif key == "HASHTT":
        # Versió hashed a n_bits (HashTT al paper)
        return rdMolDescriptors.GetHashedTopologicalTorsionFingerprintAsBitVect(
            mol,
            nBits=n_bits,
        )

    # --- Circular fingerprints (Morgan / AEs / FCFP / ECFP) ---
    elif key == "AES":
        # AEs: Morgan radius 1 en forma esparsa (sense hashing a n_bits)
        # Aquest és el "sparse ECFP" que fan servir com a atomic environments.
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=1)
        return gen.GetSparseCountFingerprint(mol)

    elif key == "ECFP0":
        # ECFP0: Morgan radius 0 hashed a n_bits (bit vector)
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=0, fpSize=n_bits)
        return gen.GetFingerprint(mol)

    # --- ECFP2 / ECFP4 ES PARS (COUNTED) ---
    elif key == "ECFP2":
        # ECFP2 (sense estrella): forma esparsa / counted (no bit vector explícit)
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=1)
        return gen.GetSparseCountFingerprint(mol)

    elif key == "ECFP4":
        # ECFP4 (sense estrella): forma esparsa / counted
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2)
        return gen.GetSparseCountFingerprint(mol)

    # --- ECFP2* / ECFP4* EXPLÍCITS EN BITS ---
    elif key == "ECFP2*":
        # ECFP2*: versió explícita en bits (hashing a n_bits)
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=1, fpSize=n_bits)
        return gen.GetFingerprint(mol)

    elif key == "ECFP4*":
        # ECFP4*: versió explícita en bits (hashing a n_bits)
        gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=n_bits)
        return gen.GetFingerprint(mol)

    elif key == "FCFP2":
        # FCFP2: Morgan radius 1 amb features farmacofòriques (bit vector)
        gen = rdFingerprintGenerator.GetMorganGenerator(
            radius=1,
            fpSize=n_bits,
            useFeatures=True,
        )
        return gen.GetFingerprint(mol)

    elif key == "FCFP4":
        # FCFP4: Morgan radius 2 amb features farmacofòriques (bit vector)
        gen = rdFingerprintGenerator.GetMorganGenerator(
            radius=2,
            fpSize=n_bits,
            useFeatures=True,
        )
        return gen.GetFingerprint(mol)

    # --- Cap cas coincideix: fingerprint no suportat ---
    else:
        raise ValueError(
            f"Tipus de fingerprint no suportat: '{fp_type}'. "
            "Comprova que coincideix amb algun dels utilitzats al paper "
            "(MACCS, Avalon, HashAP, RDK4, RDK4-L, TT, HashTT, "
            "AEs, ECFP0/2/4, ECFP2*/4*, FCFP2/4)."
        )


def smiles_to_fingerprint(smiles, fp_type: str = "ECFP4", n_bits: int = 2048):
    """
    Converteix un SMILES directament en fingerprint, gestionant None i SMILES invàlids.

    - Si 'smiles' és None -> retorna None
    - Si el SMILES no és vàlid -> retorna el string "InvalidSMILE"
    - Si és vàlid -> retorna el fingerprint corresponent

    Aquest és el punt d'entrada recomanat si treballes amb un DataFrame.
    """
    # Cas especial: valor "missing"
    if smiles is None:
        return None

    mol = smiles_to_mol(smiles)

    if mol is None:
        # SMILES rebutjat per RDKit
        return "InvalidSMILE"

    return mol_to_fingerprint(mol, fp_type=fp_type, n_bits=n_bits)


# Funció auxiliar: fingerprint RDKit -> string d'índexs activats
def fp_to_indices_str(fp):
    # Propaguem els casos especials
    if fp is None:
        return None
    if fp == "InvalidSMILE":
        return "InvalidSMILE"

    # ECFP4 (sense estrella) ve com ULongSparseIntVect (sparse counted)
    # Altres bitvectors tenen GetOnBits()
    if hasattr(fp, "GetNonzeroElements"):
        # dict {index: count}; ens quedem només amb els índexs
        elems = fp.GetNonzeroElements()
        indices = sorted(elems.keys())
    elif hasattr(fp, "GetOnBits"):
        indices = list(fp.GetOnBits())
    else:
        raise TypeError(f"No sé convertir fingerprint de tipus {type(fp)}")

    return " ".join(str(i) for i in indices)


def get_supported_fingerprints():
    """
    Retorna la llista de noms de fingerprints suportats (per conveniència).
    """
    return [
        "MACCS",
        "Avalon",
        "HashAP",
        "RDK4",
        "RDK4-L",
        "TT",
        "HashTT",
        "AEs",
        "ECFP0",
        "ECFP2",
        "ECFP4",
        "ECFP2*",
        "ECFP4*",
        "FCFP2",
        "FCFP4",
    ]