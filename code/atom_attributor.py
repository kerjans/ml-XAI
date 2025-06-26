from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdchem import PeriodicTable
import numpy as np
import pandas as pd

from xai_selfies.ml_helper import *

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
    global ORGANIC_ATOM_SYMBOLS
    if mutation_subset is None:
        mutation_subset = ORGANIC_ATOM_SYMBOLS
    mol = Chem.MolFromSmiles(smiles, sanitize=False) # to keep the explicit hydrogens
    Chem.SanitizeMol(mol)  # to keep the explicit hydrogens
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
    data_smiles = {'SMILES': [smiles]}
    df_smiles = pd.DataFrame(data_smiles)
    df_smiles['Feature'] = df_smiles['SMILES'].apply(featureMETHOD)
    prep_features_smiles = get_features(df_smiles, ['Feature'])
    prediction = model.predict(prep_features_smiles)
    return prediction


def attribute_atoms(smiles: str, model, featureMETHOD) -> np.array:
    mutated_dict = {}

    for atom_idx, old_sym, new_sym, mutated in mutate_atoms(smiles):
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
            mutated_df['Feature'] = mutated_df[0].apply(featureMETHOD)
            prep_features_mutat = get_features(mutated_df, ['Feature'])
            prediction = model.predict(prep_features_mutat)
            y_diff = y_org - prediction
            attributions.append(y_diff.mean())

    mol = Chem.MolFromSmiles(smiles, sanitize=False) # to keep the explicit hydrogens
    Chem.SanitizeMol(mol)  # to keep the explicit hydrogens
    assert len(attributions) == mol.GetNumAtoms()
    return attributions

if __name__ == "__main__":
    smiles_input = "CCO"  # Ethanol
    print(f"Original SMILES: {smiles_input}")
    print("Generated Mutations:")

    for atom_idx, old_sym, new_sym, mutated in mutate_atoms(smiles_input):
        print(f"Atom {atom_idx} ({old_sym} → {new_sym}): {mutated}")
