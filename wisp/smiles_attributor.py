"""
A SmilesAttributor as implemented in 
https://github.com/Bayer-Group/xBCF/blob/master/src/attributor.py
(see also research article: https://doi.org/10.1016/j.ailsci.2022.100047)
but adapted to work in general (independent of CDDDs).

Follows the same logic as outlined in the publication, but adapts it
to be mutating on the SMILES without encoding them as CDDDs.

Necessarily, this means that we are checking the validity of the generated 
SMILES, as we cannot rely on the CDDD model embedding it regardless
of validity.

Besides this aspect, the implementation is identical from the one
of the paper described in the reference.

Only works correctly for Python Versions 3.7 and upwards, as it contains
an ordering bug. We did not fix this bug as we wanted to keep with
the original implementation.



"""

import os
import re
import numpy as np
import pandas as pd
from random import randint
import pickle
from rdkit import Chem
from rdkit.Chem.Draw import SimilarityMaps

try:
    from cddd.inference import InferenceModel
except:
    pass


char_dict = {
    0: '</s>',
    1: '#',
    2: '%',
    3: ')',
    4: '(',
    5: '+',
    6: '-',
    7: '1',
    8: '0',
    9: '3',
    10: '2',
    11: '5',
    12: '4',
    13: '7',
    14: '6',
    15: '9',
    16: '8',
    17: ':',
    18: '=',
    19: '@',
    20: 'C',
    21: 'B',
    22: 'F',
    23: 'I',
    24: 'H',
    25: 'O',
    26: 'N',
    27: 'P',
    28: 'S',
    29: '[',
    30: ']',
    31: 'c',
    32: 'i',
    33: 'o',
    34: 'n',
    35: 'p',
    36: 's',
    37: 'Cl',
    38: 'Br',
    39: '<s>'
    }

from rdkit import Chem


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


def num_chlorines_model(smi):
    mol = Chem.MolFromSmiles(smi)
    return sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'Cl')

def predictor_smiles(smiles, featureMETHOD, model):
    data_smiles = {'SMILES': [smiles]}
    df_smiles = pd.DataFrame(data_smiles)
    df_smiles['Feature'] = df_smiles['SMILES'].apply(featureMETHOD)
    prep_features_smiles = get_features(df_smiles, ['Feature'])
    prediction = model.predict(prep_features_smiles)
    return prediction

def attribute_smiles(smiles:str,model,featureMETHOD) -> np.array:
    """
    
    >>> attribute_smiles("CC(Cl)CCCCl",num_chlorines_model)
    [-0.1, 0.0, 0.9, 0.0, 0.0, 0.0, 0.9]
    
    
    """
    per_atom_mutations = generate_mutations(smiles)

    y_org = predictor_smiles(smiles, featureMETHOD, model)
    attributions = []
    for mutations in per_atom_mutations:
        y_mut = [predictor_smiles(mutation, featureMETHOD, model) for mutation in mutations]
        y_diff = y_org - np.array(y_mut)
        attributions.append(y_diff.mean())
    
    assert len(attributions) == Chem.MolFromSmiles(smiles).GetNumAtoms()
    return attributions


def generate_mutations(sml:str, method="substitution", validate=True,) -> list[str]:
    """
    Generates all mutations for the given smiles.

    validate: Whether to only return valid smiles.

    method: Currently only substitution supported.

    >>> generate_mutations("C(=O)COC")
    [['C(=O)COC', 'B(=O)COC', 'I(=O)COC', 'N(=O)COC', 'P(=O)COC', 'S(=O)COC'], ['C(=C)COC', 'C(=B)COC', 'C(=I)COC', 'C(=O)COC', 'C(=N)COC', 'C(=P)COC', 'C(=S)COC'], ['C(=O)-OC', 'C(=O):OC', 'C(=O)COC', 'C(=O)BOC', 'C(=O)IOC', 'C(=O)OOC', 'C(=O)NOC', 'C(=O)POC', 'C(=O)SOC'], ['C(=O)C#C', 'C(=O)C-C', 'C(=O)C:C', 'C(=O)C=C', 'C(=O)CCC', 'C(=O)CBC', 'C(=O)CIC', 'C(=O)COC', 'C(=O)CNC', 'C(=O)CPC', 'C(=O)CSC'], ['C(=O)COC', 'C(=O)COB', 'C(=O)COF', 'C(=O)COI', 'C(=O)COO', 'C(=O)CON', 'C(=O)COP', 'C(=O)COS', 'C(=O)COCl', 'C(=O)COBr']]

    Testing it on two randomly taken examples from the Moleculenet Lipophilicity dataset:

    First example:
    >>> smi = "OC(=O)c1cccc(c1)N2CCC(CN3CCC(CC3)Oc4ccc(Cl)c(Cl)c4)CC2"
    >>> muts = generate_mutations(smi)

    There have to be as many mutations as there atoms, as we generate a list
    of mutations for each atom:
    >>> assert len(muts) == Chem.MolFromSmiles(smi).GetNumAtoms() 

    Check that the first mutations only affect the first atom:
    >>> print(muts[0][1])
    BC(=O)c1cccc(c1)N2CCC(CN3CCC(CC3)Oc4ccc(Cl)c(Cl)c4)CC2

    >>> print(muts[0][-1])
    BrC(=O)c1cccc(c1)N2CCC(CN3CCC(CC3)Oc4ccc(Cl)c(Cl)c4)CC2

    Check that the last mutations only affect the last atom:
    >>> print(muts[-1][0])
    OC(=O)c1cccc(c1)N2CCC(CN3CCC(CC3)Oc4ccc(Cl)c(Cl)c4)C-2

    >>> print(muts[-1][-1])
    OC(=O)c1cccc(c1)N2CCC(CN3CCC(CC3)Oc4ccc(Cl)c(Cl)c4)Cs2

    Second Example:
    >>> smi = "Cc1cnc(cn1)C(=O)NCCc2ccc(cc2)S(=O)(=O)NC(=O)NC3CCCCC3"
    >>> muts = generate_mutations(smi)
    >>> assert len(muts) == Chem.MolFromSmiles(smi).GetNumAtoms() 

    Check that the first mutations only affect the first atom:
    >>> print(muts[0][1])
    Bc1cnc(cn1)C(=O)NCCc2ccc(cc2)S(=O)(=O)NC(=O)NC3CCCCC3

    >>> print(muts[0][-1])
    Brc1cnc(cn1)C(=O)NCCc2ccc(cc2)S(=O)(=O)NC(=O)NC3CCCCC3

    And the last mutations only affect the last atom:
    >>> print(muts[-1][0])
    Cc1cnc(cn1)C(=O)NCCc2ccc(cc2)S(=O)(=O)NC(=O)NC3CCCC-3

    >>> print(muts[-1][5])
    Cc1cnc(cn1)C(=O)NCCc2ccc(cc2)S(=O)(=O)NC(=O)NC3CCCCI3

    """
    n_atoms = Chem.MolFromSmiles(sml).GetNumAtoms()
    REGEX_SML = r'Cl|Br|[#%\)\(\+\-1032547698:=@CBFIHONPS\[\]cionps]'
    atom_chars = ["Cl","Br","C","B","F","I","H","O","N","P","S","c","i","o","n","p","s"]
    char_list = re.findall(REGEX_SML, sml)
    atom_char_list = [i for i in range(len(char_list)) if char_list[i] in atom_chars]
    #matches = re.finditer(REGEX_SML, sml)
    #indices = [match.start() for match in matches]

    assert len(atom_char_list) == n_atoms
    char_vocab = list(char_dict.values())[1:-1] # this line is a bug for Python versions earlier than 3.7
    mutations_by_char = []
    if method == "substitution":
        for i, sml_char in enumerate(char_list):
            if sml_char not in atom_chars:
                continue

            sml_copy = char_list
            mutated_smls = ["".join(sml_copy[:i] + [w] + sml_copy[i+1:]) for w in char_vocab]
            if validate:
                good_muts = []
                for mut in mutated_smls:
                    try:
                        with _Silencer() as _:
                            mol = Chem.MolFromSmiles(mut)
                        if mol is not None:
                            good_muts.append(mut)
                    except:
                        pass
                mutated_smls = good_muts

            mutations_by_char.append(mutated_smls)
    else:
        assert False
    return mutations_by_char
