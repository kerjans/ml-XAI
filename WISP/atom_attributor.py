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

                    Chem.SanitizeMol(new_mol)
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
    if "mordred" in featureMETHOD.__name__.lower():
        # First attempt to run in batch mode if possible
        df_smiles['Feature'] = list(featureMETHOD(df_smiles['SMILES'].tolist()))
    else:
        # Feature function does not support batch mode, run one by one
        df_smiles['Feature'] = df_smiles['SMILES'].apply(featureMETHOD)

    prep_features_smiles = get_features(df_smiles, ['Feature'])
    prediction = model.predict(prep_features_smiles)
    return prediction


def attribute_atoms(smiles: str, model, featureMETHOD) -> np.array:
    """
    Compute per-atom attribution scores by single-atom substitutions.

    For each atom in the input SMILES:
      1. Generate all valid single-atom replacements via `mutate_atoms()`.
      2. Predict the original molecule’s output.
      3. For each atom, mutate it in turn, predict the mutated molecules,
         and compute the average Δprediction relative to the original.
      4. Return a list of attribution scores, one per atom index.

    Parameters:
        smiles (str): SMILES string of the molecule.
        model: A fitted estimator with `.predict(X)` and compatible with `predictor_on_smiles`.
        feature_method (callable): Function mapping a SMILES to its feature vector(s).

    Returns:
        np.array: An array of attribution scores (floats), one per atom in the molecule.

    """
    mutated_dict = {}

    for atom_idx, _, _, mutated in mutate_atoms(smiles):
        if atom_idx not in mutated_dict:
            mutated_dict[atom_idx] = []
        if mutated:
            mutated_dict[atom_idx].append(mutated)

    y_org = predictor_on_smiles(smiles, featureMETHOD, model)
    attributions = []
    for index in mutated_dict:
        mutated_df = pd.DataFrame(mutated_dict[index])
        if mutated_df.empty:#no mutations were generated
            attributions.append(np.nan)
        else:

            if "mordred" in featureMETHOD.__name__.lower():
                mutated_df['Feature'] = list(featureMETHOD(mutated_df[0]))
            else:
                mutated_df['Feature'] = mutated_df[0].apply(else_none(featureMETHOD))

            # It can happen that during morgan fingerprint calculation,
            # we find that the mutation is somehow an invalid molecule.
            # In this case, we find out only here.
            mutated_df = mutated_df[~mutated_df["Feature"].isna()]

            if mutated_df.empty:#no mutations were generated
                attributions.append(np.nan)
            else:
                prep_features_mutat = get_features(mutated_df, ['Feature'])
                prediction = model.predict(prep_features_mutat)
                y_diff = y_org - prediction
                attributions.append(y_diff.mean())

    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    Chem.SanitizeMol(mol)  # to keep the explicit hydrogens
    assert len(attributions) == mol.GetNumAtoms()
    return attributions

if __name__ == "__main__":
    smiles_input = "CCO"  # Ethanol
    print(f"Original SMILES: {smiles_input}")
    print("Generated Mutations:")

    for atom_idx, old_sym, new_sym, mutated in mutate_atoms(smiles_input):
        print(f"Atom {atom_idx} ({old_sym} → {new_sym}): {mutated}")
