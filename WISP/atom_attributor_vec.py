from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import PeriodicTable
import numpy as np
from standardizer.io import else_none
np.seterr(divide='ignore', invalid='ignore')
import pandas as pd

from WISP.ml_helper import *

# Organic subset as defined by RDKit
ORGANIC_ATOM_SYMBOLS = [
    'H', 'B', 'C', 'N', 'O', 'F', 'Si', 'P', 'S', 'Cl',
    'Br', 'I'
]


import io
from contextlib import redirect_stderr, redirect_stdout
class _Silencer:
    """
    A useful tool for silencing stdout and stderr.
    Usage:
    >>> with _Silencer() as s:
    ...         print("kasldjf")

    >>> print("I catched:",s.out.getvalue())
    I catched: kasldjf
    <BLANKLINE>

    Note that nothing was printed and that we can later
    access the stdout via the out field. Similarly,
    stderr will be redirected to the err field.
    """

    def __init__(self):
        self.out = io.StringIO()
        self.err = io.StringIO()

    def __enter__(self):
        from rdkit import RDLogger

        RDLogger.DisableLog("rdApp.*")
        self.rs = redirect_stdout(self.out)
        self.re = redirect_stderr(self.err)
        self.rs.__enter__()
        self.re.__enter__()
        return self

    def __exit__(self, exctype, excinst, exctb):
        from rdkit import RDLogger

        RDLogger.EnableLog("rdApp.*")
        self.rs.__exit__(exctype, excinst, exctb)
        self.re.__exit__(exctype, excinst, exctb)

def mutate_atoms(smiles, mutation_subset=None):
    """
    Generate single‐atom point mutations of a molecule by substituting each atom
    in turn with elements from a specified set.

    For each atom index in the input SMILES:
      • The atom’s original symbol is recorded.
      • Each symbol in `mutation_subset` (or the global ORGANIC_ATOM_SYMBOLS)—
        except the atom’s current symbol—is tried as a replacement.
      • The molecule is sanitized after the swap; if valid, the new SMILES is yielded.
      • On failure, a tuple with replacement_symbol=None and new_smiles='' is yielded.

    Parameters:
        smiles (str): A SMILES string of the molecule to mutate.
        mutation_subset (iterable of str, optional): Atom symbols to attempt as replacements.
            Defaults to the global ORGANIC_ATOM_SYMBOLS.

    Yields:
        tuple:
            (atom_idx (int),
             original_symbol (str),
             replacement_symbol (str or None),
             new_smiles (str))
    """
    global ORGANIC_ATOM_SYMBOLS
    if mutation_subset is None:
        mutation_subset = ORGANIC_ATOM_SYMBOLS
    mol = Chem.MolFromSmiles(smiles, sanitize=False) 
    Chem.SanitizeMol(mol) # to keep the explicit hydrogens
    if mol is None or mol.GetNumAtoms() == 0:
        return mol # nothing to mutate

    mol = Chem.RWMol(mol)

    for atom_idx in range(mol.GetNumAtoms()):
        original_atom = mol.GetAtomWithIdx(atom_idx)
        original_symbol = original_atom.GetSymbol()

        for replacement_symbol in ORGANIC_ATOM_SYMBOLS:
            if replacement_symbol == original_symbol:
                continue

            try:
                with _Silencer() as _:
                    new_mol = Chem.RWMol(mol)
                    atom = new_mol.GetAtomWithIdx(atom_idx)

                    new_atomic_num = Chem.GetPeriodicTable().GetAtomicNumber(replacement_symbol)
                    atom.SetAtomicNum(new_atomic_num)

                    Chem.SanitizeMol(new_mol) # ?
                    new_smiles = Chem.MolToSmiles(new_mol, isomericSmiles=True)

                yield (atom_idx, original_symbol, replacement_symbol, new_smiles)
            except Exception:
                # Skip invalid mutations
                yield (atom_idx, original_symbol, None, '')
                continue


def predictor_on_smiles(smiles, featureMETHOD, model):
    """
    Compute a model prediction for a single molecule/mutation given its SMILES string.

    Parameters:
        smiles (str): A SMILES representation of the molecule to predict on.
        feature_method (callable): A function that maps a SMILES string to feature array.
        model: A fitted estimator with a .predict(X) method, where X is the output of get_features().

    Returns:
        np.ndarray: The prediction output (often a 1-element array) from model.predict().
    """
    data_smiles = {'SMILES': [smiles]}
    df_smiles = pd.DataFrame(data_smiles)
    df_smiles['Feature'] = list(featureMETHOD(df_smiles['SMILES'].tolist()))
    prep_features_smiles = get_features(df_smiles, ['Feature'])
    prediction = model.predict(prep_features_smiles)
    return prediction

def attribute_atoms(smiles_list: "list[str]", model, featureMETHOD, chunk_size=256,) -> np.array:
    rslt = []
    chnks = list(range(0,len(smiles_list),chunk_size))
    for i in chnks:
        print(f"processing chunk: {i} / {len(chnks)}")
        chnk = smiles_list[i:i+chunk_size]
        if len(chnk):
            rslt = rslt.extend(_attribute_atoms_chunk(chnk,model,featureMETHOD))

    assert len(rslt) == len(smiles_list)
    return rslt

def _attribute_atoms_chunk(smiles_list: "list[str]", model, featureMETHOD) -> np.array:
    df_muts = [] # smiles_org, atom_idx, smiles_mut, feature, y_org, y_mut, y_diff

    for smiles in smiles_list:
        for atom_idx, _, _, mutated in mutate_atoms(smiles):
            df_muts.append({
                "smiles_org": smiles,
                "smiles_mut": mutated,
                "atom_idx": atom_idx,
            })

    df_muts = pd.DataFrame(df_muts)
    df_muts = df_muts[df_muts["smiles_mut"] != ""]
    
    df_pred = pd.DataFrame({"smiles": list(smiles_list)+df_muts.smiles_mut.tolist(),})
    df_pred['Feature'] = list(featureMETHOD(df_pred["smiles"]))
    prep_features = get_features(df_pred, ['Feature'])
    prediction = model.predict(prep_features)

    df_pred["y_pred"] = list(prediction)
    smi_to_pred = {row["smiles"]: row["y_pred"] for _,row in df_pred.iterrows()}

    df_muts["y_org"] = df_muts.smiles_org.map(smi_to_pred)
    df_muts["y_mut"] = df_muts.smiles_mut.map(smi_to_pred)
    df_muts["y_diff"] = df_muts["y_org"] - df_muts["y_mut"]

    attributions_all = []
    for smiles in smiles_list:
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        Chem.SanitizeMol(mol)  # to keep the explicit hydrogens

        g_smi = df_muts[df_muts.smiles_org == smiles]

        attributions = np.zeros(mol.GetNumAtoms())
        for atm_idx in g_smi["atom_idx"].unique():
            g_smi_atom = g_smi[g_smi["atom_idx"] == atm_idx]
            attributions[atm_idx] += g_smi_atom["y_diff"].mean()

        attributions_all.append(attributions)

    return attributions_all

if __name__ == "__main__":
    smiles_input = "CCO"  # Ethanol
    print(f"Original SMILES: {smiles_input}")
    print("Generated Mutations:")

    for atom_idx, old_sym, new_sym, mutated in mutate_atoms(smiles_input):
        print(f"Atom {atom_idx} ({old_sym} → {new_sym}): {mutated}")
