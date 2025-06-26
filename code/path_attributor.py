from rdkit.Chem import rdmolops
from rdkit import Chem
import numpy as np
import pandas as pd
from xai_selfies.ml_helper import *

def find_all_paths(mol, start_idx, end_idx, max_depth=10):
    paths = []

    def dfs(current_path, visited):
        current_atom = current_path[-1]
        if current_atom == end_idx:
            paths.append(list(current_path))
            return
        if len(current_path) > max_depth:
            return
        for neighbor in mol.GetAtomWithIdx(current_atom).GetNeighbors():
            neighbor_idx = neighbor.GetIdx()
            if neighbor_idx not in visited and neighbor.GetAtomicNum() != 0:
                dfs(current_path + [neighbor_idx], visited | {neighbor_idx})

    dfs([start_idx], {start_idx})
    return paths


def remove_atom_at_index(mol: Chem.Mol, atm_idx: int) -> Chem.Mol:
    """
    Removes the atom at the given index atm_idx.

    """
    mol = copy.deepcopy(mol)
    if isinstance(atm_idx,list) or isinstance(atm_idx,tuple):
        for ai in atm_idx:
            mol.GetAtomWithIdx(ai).SetAtomicNum(0)
    else:
        mol.GetAtomWithIdx(atm_idx).SetAtomicNum(0)
    mol = Chem.DeleteSubstructs(mol, Chem.MolFromSmarts("[#0]"))
    return mol

import copy
def mutations_along_paths(mol):
    """
    
    >>> smi = "CCCO"

    >>> list(mutations_along_paths(Chem.MolFromSmiles(smi)))
    [('CCO', [0]), ('CO', [0, 1]), ('O', [0, 1, 2]), ('', [0, 1, 2, 3]), ('CO', [1, 0]), ('C.CO', [1]), ('C.O', [1, 2]), ('C', [1, 2, 3]), ('O', [2, 1, 0]), ('C.O', [2, 1]), ('CC.O', [2]), ('CC', [2, 3]), ('', [3, 2, 1, 0]), ('C', [3, 2, 1]), ('CC', [3, 2]), ('CCC', [3])]
    >>> smi = "CCCOC1CCCCC1CCCBr"
    >>> list(mutations_along_paths(Chem.MolFromSmiles(smi)))
    [('CCOC1CCCCC1CCCBr', [0]), ('COC1CCCCC1CCCBr', [0, 1]), ('OC1CCCCC1CCCBr', [0, 1, 2]), ('BrCCCC1CCCCC1', [0, 1, 2, 3]), ('CCCCCCCCBr', [0, 1, 2, 3, 4]), ('CCCBr', [0, 1, 2, 3, 4, 5, 0, 1, 2, 3, 4, 9, 8, 7, 6, 5]), ('CCCBr', [0, 1, 2, 3, 4, 5, 6, 0, 1, 2, 3, 4, 9, 8, 7, 6]), ('CCCBr', [0, 1, 2, 3, 4, 5, 6, 7, 0, 1, 2, 3, 4, 9, 8, 7]), ('CCCBr', [0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 1, 2, 3, 4, 9, 8]), ('CCCBr', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 0, 1, 2, 3, 4, 9]), ('CCBr', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 0, 1, 2, 3, 4, 9, 10]), ('CBr', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 0, 1, 2, 3, 4, 9, 10, 11]), ('Br', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 0, 1, 2, 3, 4, 9, 10, 11, 12]), ('', [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 0, 1, 2, 3, 4, 9, 10, 11, 12, 13]), ('COC1CCCCC1CCCBr', [1, 0]), ('C.COC1CCCCC1CCCBr', [1]), ('C.OC1CCCCC1CCCBr', [1, 2]), ('BrCCCC1CCCCC1.C', [1, 2, 3]), ('C.CCCCCCCCBr', [1, 2, 3, 4]), ('C.CCCBr', [1, 2, 3, 4, 5, 1, 2, 3, 4, 9, 8, 7, 6, 5]), ('C.CCCBr', [1, 2, 3, 4, 5, 6, 1, 2, 3, 4, 9, 8, 7, 6]), ('C.CCCBr', [1, 2, 3, 4, 5, 6, 7, 1, 2, 3, 4, 9, 8, 7]), ('C.CCCBr', [1, 2, 3, 4, 5, 6, 7, 8, 1, 2, 3, 4, 9, 8]), ('C.CCCBr', [1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 9]), ('C.CCBr', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 1, 2, 3, 4, 9, 10]), ('C.CBr', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 1, 2, 3, 4, 9, 10, 11]), ('Br.C', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 9, 10, 11, 12]), ('C', [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 1, 2, 3, 4, 9, 10, 11, 12, 13]), ('OC1CCCCC1CCCBr', [2, 1, 0]), ('C.OC1CCCCC1CCCBr', [2, 1]), ('CC.OC1CCCCC1CCCBr', [2]), ('BrCCCC1CCCCC1.CC', [2, 3]), ('CC.CCCCCCCCBr', [2, 3, 4]), ('CC.CCCBr', [2, 3, 4, 5, 2, 3, 4, 9, 8, 7, 6, 5]), ('CC.CCCBr', [2, 3, 4, 5, 6, 2, 3, 4, 9, 8, 7, 6]), ('CC.CCCBr', [2, 3, 4, 5, 6, 7, 2, 3, 4, 9, 8, 7]), ('CC.CCCBr', [2, 3, 4, 5, 6, 7, 8, 2, 3, 4, 9, 8]), ('CC.CCCBr', [2, 3, 4, 5, 6, 7, 8, 9, 2, 3, 4, 9]), ('CC.CCBr', [2, 3, 4, 5, 6, 7, 8, 9, 10, 2, 3, 4, 9, 10]), ('CBr.CC', [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 2, 3, 4, 9, 10, 11]), ('Br.CC', [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 2, 3, 4, 9, 10, 11, 12]), ('CC', [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 2, 3, 4, 9, 10, 11, 12, 13]), ('BrCCCC1CCCCC1', [3, 2, 1, 0]), ('BrCCCC1CCCCC1.C', [3, 2, 1]), ('BrCCCC1CCCCC1.CC', [3, 2]), ('BrCCCC1CCCCC1.CCC', [3]), ('CCC.CCCCCCCCBr', [3, 4]), ('CCC.CCCBr', [3, 4, 5, 3, 4, 9, 8, 7, 6, 5]), ('CCC.CCCBr', [3, 4, 5, 6, 3, 4, 9, 8, 7, 6]), ('CCC.CCCBr', [3, 4, 5, 6, 7, 3, 4, 9, 8, 7]), ('CCC.CCCBr', [3, 4, 5, 6, 7, 8, 3, 4, 9, 8]), ('CCC.CCCBr', [3, 4, 5, 6, 7, 8, 9, 3, 4, 9]), ('CCBr.CCC', [3, 4, 5, 6, 7, 8, 9, 10, 3, 4, 9, 10]), ('CBr.CCC', [3, 4, 5, 6, 7, 8, 9, 10, 11, 3, 4, 9, 10, 11]), ('Br.CCC', [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 3, 4, 9, 10, 11, 12]), ('CCC', [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 3, 4, 9, 10, 11, 12, 13]), ('CCCCCCCCBr', [4, 3, 2, 1, 0]), ('C.CCCCCCCCBr', [4, 3, 2, 1]), ('CC.CCCCCCCCBr', [4, 3, 2]), ('CCC.CCCCCCCCBr', [4, 3]), ('CCCCCCCCBr.CCCO', [4]), ('CCCBr.CCCO', [4, 5, 4, 9, 8, 7, 6, 5]), ('CCCBr.CCCO', [4, 5, 6, 4, 9, 8, 7, 6]), ('CCCBr.CCCO', [4, 5, 6, 7, 4, 9, 8, 7]), ('CCCBr.CCCO', [4, 5, 6, 7, 8, 4, 9, 8]), ('CCCBr.CCCO', [4, 5, 6, 7, 8, 9, 4, 9]), ('CCBr.CCCO', [4, 5, 6, 7, 8, 9, 10, 4, 9, 10]), ('CBr.CCCO', [4, 5, 6, 7, 8, 9, 10, 11, 4, 9, 10, 11]), ('Br.CCCO', [4, 5, 6, 7, 8, 9, 10, 11, 12, 4, 9, 10, 11, 12]), ('CCCO', [4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 4, 9, 10, 11, 12, 13]), ('CCCBr', [5, 4, 3, 2, 1, 0, 5, 6, 7, 8, 9, 4, 3, 2, 1, 0]), ('C.CCCBr', [5, 4, 3, 2, 1, 5, 6, 7, 8, 9, 4, 3, 2, 1]), ('CC.CCCBr', [5, 4, 3, 2, 5, 6, 7, 8, 9, 4, 3, 2]), ('CCC.CCCBr', [5, 4, 3, 5, 6, 7, 8, 9, 4, 3]), ('CCCBr.CCCO', [5, 4, 5, 6, 7, 8, 9, 4]), ('CCCOCC(CCC)CCCBr', [5]), ('CCCBr.CCCO', [5, 4, 9, 8, 7, 6, 5, 6]), ('CCCBr.CCCO', [5, 4, 9, 8, 7, 5, 6, 7]), ('CCCBr.CCCO', [5, 4, 9, 8, 5, 6, 7, 8]), ('CCCBr.CCCO', [5, 4, 9, 5, 6, 7, 8, 9]), ('CCBr.CCCO', [5, 4, 9, 10, 5, 6, 7, 8, 9, 10]), ('CBr.CCCO', [5, 4, 9, 10, 11, 5, 6, 7, 8, 9, 10, 11]), ('Br.CCCO', [5, 4, 9, 10, 11, 12, 5, 6, 7, 8, 9, 10, 11, 12]), ('CCCO', [5, 4, 9, 10, 11, 12, 13, 5, 6, 7, 8, 9, 10, 11, 12, 13]), ('CCCBr', [6, 5, 4, 3, 2, 1, 0, 6, 7, 8, 9, 4, 3, 2, 1, 0]), ('C.CCCBr', [6, 5, 4, 3, 2, 1, 6, 7, 8, 9, 4, 3, 2, 1]), ('CC.CCCBr', [6, 5, 4, 3, 2, 6, 7, 8, 9, 4, 3, 2]), ('CCC.CCCBr', [6, 5, 4, 3, 6, 7, 8, 9, 4, 3]), ('CCCBr.CCCO', [6, 5, 4, 6, 7, 8, 9, 4]), ('CCCBr.CCCO', [6, 5, 6, 7, 8, 9, 4, 5]), ('CCCOC(C)C(CC)CCCBr', [6]), ('CCCBr.CCCO', [6, 5, 4, 9, 8, 7, 6, 7]), ('CCCBr.CCCO', [6, 5, 4, 9, 8, 6, 7, 8]), ('CCCBr.CCCO', [6, 5, 4, 9, 6, 7, 8, 9]), ('CCBr.CCCO', [6, 5, 4, 9, 10, 6, 7, 8, 9, 10]), ('CBr.CCCO', [6, 5, 4, 9, 10, 11, 6, 7, 8, 9, 10, 11]), ('Br.CCCO', [6, 5, 4, 9, 10, 11, 12, 6, 7, 8, 9, 10, 11, 12]), ('CCCO', [6, 5, 4, 9, 10, 11, 12, 13, 6, 7, 8, 9, 10, 11, 12, 13]), ('CCCBr', [7, 6, 5, 4, 3, 2, 1, 0, 7, 8, 9, 4, 3, 2, 1, 0]), ('C.CCCBr', [7, 6, 5, 4, 3, 2, 1, 7, 8, 9, 4, 3, 2, 1]), ('CC.CCCBr', [7, 6, 5, 4, 3, 2, 7, 8, 9, 4, 3, 2]), ('CCC.CCCBr', [7, 6, 5, 4, 3, 7, 8, 9, 4, 3]), ('CCCBr.CCCO', [7, 6, 5, 4, 7, 8, 9, 4]), ('CCCBr.CCCO', [7, 6, 5, 7, 8, 9, 4, 5]), ('CCCBr.CCCO', [7, 6, 7, 8, 9, 4, 5, 6]), ('CCCOC(CC)C(C)CCCBr', [7]), ('CCCBr.CCCO', [7, 6, 5, 4, 9, 8, 7, 8]), ('CCCBr.CCCO', [7, 6, 5, 4, 9, 7, 8, 9]), ('CCBr.CCCO', [7, 6, 5, 4, 9, 10, 7, 8, 9, 10]), ('CBr.CCCO', [7, 6, 5, 4, 9, 10, 11, 7, 8, 9, 10, 11]), ('Br.CCCO', [7, 6, 5, 4, 9, 10, 11, 12, 7, 8, 9, 10, 11, 12]), ('CCCO', [7, 6, 5, 4, 9, 10, 11, 12, 13, 7, 8, 9, 10, 11, 12, 13]), ('CCCBr', [8, 7, 6, 5, 4, 3, 2, 1, 0, 8, 9, 4, 3, 2, 1, 0]), ('C.CCCBr', [8, 7, 6, 5, 4, 3, 2, 1, 8, 9, 4, 3, 2, 1]), ('CC.CCCBr', [8, 7, 6, 5, 4, 3, 2, 8, 9, 4, 3, 2]), ('CCC.CCCBr', [8, 7, 6, 5, 4, 3, 8, 9, 4, 3]), ('CCCBr.CCCO', [8, 7, 6, 5, 4, 8, 9, 4]), ('CCCBr.CCCO', [8, 7, 6, 5, 8, 9, 4, 5]), ('CCCBr.CCCO', [8, 7, 6, 8, 9, 4, 5, 6]), ('CCCBr.CCCO', [8, 7, 8, 9, 4, 5, 6, 7]), ('CCCOC(CCC)CCCCBr', [8]), ('CCCBr.CCCO', [8, 7, 6, 5, 4, 9, 8, 9]), ('CCBr.CCCO', [8, 7, 6, 5, 4, 9, 10, 8, 9, 10]), ('CBr.CCCO', [8, 7, 6, 5, 4, 9, 10, 11, 8, 9, 10, 11]), ('Br.CCCO', [8, 7, 6, 5, 4, 9, 10, 11, 12, 8, 9, 10, 11, 12]), ('CCCO', [8, 7, 6, 5, 4, 9, 10, 11, 12, 13, 8, 9, 10, 11, 12, 13]), ('CCCBr', [9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 9, 4, 3, 2, 1, 0]), ('C.CCCBr', [9, 8, 7, 6, 5, 4, 3, 2, 1, 9, 4, 3, 2, 1]), ('CC.CCCBr', [9, 8, 7, 6, 5, 4, 3, 2, 9, 4, 3, 2]), ('CCC.CCCBr', [9, 8, 7, 6, 5, 4, 3, 9, 4, 3]), ('CCCBr.CCCO', [9, 8, 7, 6, 5, 4, 9, 4]), ('CCCBr.CCCO', [9, 8, 7, 6, 5, 9, 4, 5]), ('CCCBr.CCCO', [9, 8, 7, 6, 9, 4, 5, 6]), ('CCCBr.CCCO', [9, 8, 7, 9, 4, 5, 6, 7]), ('CCCBr.CCCO', [9, 8, 9, 4, 5, 6, 7, 8]), ('CCCBr.CCCCCOCCC', [9]), ('CCBr.CCCCCOCCC', [9, 10]), ('CBr.CCCCCOCCC', [9, 10, 11]), ('Br.CCCCCOCCC', [9, 10, 11, 12]), ('CCCCCOCCC', [9, 10, 11, 12, 13]), ('CCBr', [10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 10, 9, 4, 3, 2, 1, 0]), ('C.CCBr', [10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 10, 9, 4, 3, 2, 1]), ('CC.CCBr', [10, 9, 8, 7, 6, 5, 4, 3, 2, 10, 9, 4, 3, 2]), ('CCBr.CCC', [10, 9, 8, 7, 6, 5, 4, 3, 10, 9, 4, 3]), ('CCBr.CCCO', [10, 9, 8, 7, 6, 5, 4, 10, 9, 4]), ('CCBr.CCCO', [10, 9, 8, 7, 6, 5, 10, 9, 4, 5]), ('CCBr.CCCO', [10, 9, 8, 7, 6, 10, 9, 4, 5, 6]), ('CCBr.CCCO', [10, 9, 8, 7, 10, 9, 4, 5, 6, 7]), ('CCBr.CCCO', [10, 9, 8, 10, 9, 4, 5, 6, 7, 8]), ('CCBr.CCCCCOCCC', [10, 9]), ('CCBr.CCCOC1CCCCC1', [10]), ('CBr.CCCOC1CCCCC1', [10, 11]), ('Br.CCCOC1CCCCC1', [10, 11, 12]), ('CCCOC1CCCCC1', [10, 11, 12, 13]), ('CBr', [11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 11, 10, 9, 4, 3, 2, 1, 0]), ('C.CBr', [11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 11, 10, 9, 4, 3, 2, 1]), ('CBr.CC', [11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 11, 10, 9, 4, 3, 2]), ('CBr.CCC', [11, 10, 9, 8, 7, 6, 5, 4, 3, 11, 10, 9, 4, 3]), ('CBr.CCCO', [11, 10, 9, 8, 7, 6, 5, 4, 11, 10, 9, 4]), ('CBr.CCCO', [11, 10, 9, 8, 7, 6, 5, 11, 10, 9, 4, 5]), ('CBr.CCCO', [11, 10, 9, 8, 7, 6, 11, 10, 9, 4, 5, 6]), ('CBr.CCCO', [11, 10, 9, 8, 7, 11, 10, 9, 4, 5, 6, 7]), ('CBr.CCCO', [11, 10, 9, 8, 11, 10, 9, 4, 5, 6, 7, 8]), ('CBr.CCCCCOCCC', [11, 10, 9]), ('CBr.CCCOC1CCCCC1', [11, 10]), ('CBr.CCCOC1CCCCC1C', [11]), ('Br.CCCOC1CCCCC1C', [11, 12]), ('CCCOC1CCCCC1C', [11, 12, 13]), ('Br', [12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 12, 11, 10, 9, 4, 3, 2, 1, 0]), ('Br.C', [12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 12, 11, 10, 9, 4, 3, 2, 1]), ('Br.CC', [12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 12, 11, 10, 9, 4, 3, 2]), ('Br.CCC', [12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 12, 11, 10, 9, 4, 3]), ('Br.CCCO', [12, 11, 10, 9, 8, 7, 6, 5, 4, 12, 11, 10, 9, 4]), ('Br.CCCO', [12, 11, 10, 9, 8, 7, 6, 5, 12, 11, 10, 9, 4, 5]), ('Br.CCCO', [12, 11, 10, 9, 8, 7, 6, 12, 11, 10, 9, 4, 5, 6]), ('Br.CCCO', [12, 11, 10, 9, 8, 7, 12, 11, 10, 9, 4, 5, 6, 7]), ('Br.CCCO', [12, 11, 10, 9, 8, 12, 11, 10, 9, 4, 5, 6, 7, 8]), ('Br.CCCCCOCCC', [12, 11, 10, 9]), ('Br.CCCOC1CCCCC1', [12, 11, 10]), ('Br.CCCOC1CCCCC1C', [12, 11]), ('Br.CCCOC1CCCCC1CC', [12]), ('CCCOC1CCCCC1CC', [12, 13]), ('', [13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, 13, 12, 11, 10, 9, 4, 3, 2, 1, 0]), ('C', [13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 13, 12, 11, 10, 9, 4, 3, 2, 1]), ('CC', [13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 13, 12, 11, 10, 9, 4, 3, 2]), ('CCC', [13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 13, 12, 11, 10, 9, 4, 3]), ('CCCO', [13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 13, 12, 11, 10, 9, 4]), ('CCCO', [13, 12, 11, 10, 9, 8, 7, 6, 5, 13, 12, 11, 10, 9, 4, 5]), ('CCCO', [13, 12, 11, 10, 9, 8, 7, 6, 13, 12, 11, 10, 9, 4, 5, 6]), ('CCCO', [13, 12, 11, 10, 9, 8, 7, 13, 12, 11, 10, 9, 4, 5, 6, 7]), ('CCCO', [13, 12, 11, 10, 9, 8, 13, 12, 11, 10, 9, 4, 5, 6, 7, 8]), ('CCCCCOCCC', [13, 12, 11, 10, 9]), ('CCCOC1CCCCC1', [13, 12, 11, 10]), ('CCCOC1CCCCC1C', [13, 12, 11]), ('CCCOC1CCCCC1CC', [13, 12]), ('CCCOC1CCCCC1CCC', [13])]
    
    """
    for idx_a in range(mol.GetNumAtoms()):
        for idx_b in range(mol.GetNumAtoms()):
            paths = find_all_paths(mol,idx_a,idx_b,max_depth=120,)
            mut = copy.deepcopy(mol)
            for path in paths:
                for ai in path:
                    mut.GetAtomWithIdx(ai).SetAtomicNum(0)

            mut = Chem.DeleteSubstructs(mut, Chem.MolFromSmarts("[#0]"))
            try:############################################### sanity check
                test = Chem.MolFromSmiles(mut, sanitize=False)
                Chem.SanitizeMol(test)
                yield Chem.MolToSmiles(mut),sum(paths,[])
            except Exception:
                    yield '',sum(paths,[])

def predictor_on_smiles(smiles, featureMETHOD, model):
    data_smiles = {'SMILES': [smiles]}
    df_smiles = pd.DataFrame(data_smiles)
    df_smiles['Feature'] = df_smiles['SMILES'].apply(featureMETHOD)
    prep_features_smiles = get_features(df_smiles, ['Feature'])
    prediction = model.predict(prep_features_smiles)
    return prediction

def attribute_atoms_paths(smiles: str, model, featureMETHOD) -> np.array:
    mutated_dict = {}

    mol = Chem.MolFromSmiles(smiles, sanitize=False) # to keep the explicit hydrogens
    Chem.SanitizeMol(mol)  # to keep the explicit hydrogens

    for mutated, atom_idx_list in mutations_along_paths(mol):
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

