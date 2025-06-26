import copy
import random
from rdkit import Chem
import numpy as np
import pandas as pd
from xai_selfies.ml_helper import *

def mutations_dropout(mol,dropout=0.25,seed=None,max_iterations=50,):
    """
    
    >>> smi = "CCCO"
    >>> list(mutations_dropout(Chem.MolFromSmiles(smi)))
    [('C.O', [1, 2]), ('CCCO', []), ('C.CO', [1]), ('CC', [0, 3]), ('CO', [0, 1]), ('CC', [2, 3]), ('C', [0, 2, 3]), ('CCO', [0]), ('CCC', [3]), ('CC.O', [2]), ('C.C', [1, 3]), ('C.O', [0, 2])]

    """
    random.seed(seed)
    visited = []
    for _ in range(max_iterations):
        ais_deleted = []
        mut = copy.deepcopy(mol)
        for ai in range(mol.GetNumAtoms()):
            if random.random() < dropout:
                ais_deleted.append(ai)
                mut.GetAtomWithIdx(ai).SetAtomicNum(0)
        mut = Chem.DeleteSubstructs(mut, Chem.MolFromSmarts("[#0]"))
        mut = Chem.MolToSmiles(mut)
        try:
            test = Chem.MolFromSmiles(mut, sanitize=False)
            Chem.SanitizeMol(test)
            if ais_deleted not in visited:
                visited.append(ais_deleted)
                yield mut,ais_deleted
        except Exception:
            yield '',ais_deleted 

def predictor_on_smiles(smiles, featureMETHOD, model):
    data_smiles = {'SMILES': [smiles]}
    df_smiles = pd.DataFrame(data_smiles)
    df_smiles['Feature'] = df_smiles['SMILES'].apply(featureMETHOD)
    prep_features_smiles = get_features(df_smiles, ['Feature'])
    prediction = model.predict(prep_features_smiles)
    return prediction

def attribute_atoms_dropout(smiles: str, model, featureMETHOD) -> np.array:
    mutated_dict = {}

    mol = Chem.MolFromSmiles(smiles, sanitize=False) # to keep the explicit hydrogens
    Chem.SanitizeMol(mol)  # to keep the explicit hydrogens

    for mutated, atom_idx_list in mutations_dropout(mol):
        for atom_idx in atom_idx_list: 
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

    assert len(attributions) == mol.GetNumAtoms()
    return attributions