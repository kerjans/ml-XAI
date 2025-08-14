from rdkit import Chem
import numpy as np

class _CountCarbonsModel:

    def predict(self,smiles_list):
        rslt = []
        for smiles in smiles_list.reshape(-1):
            mol = Chem.MolFromSmiles(smiles)
            carbs = [atm
                for atm in mol.GetAtoms()
                if atm.GetSymbol() in ["C"]
            ]
            rslt.append( len(carbs))
        return np.array(rslt)


def _identity(x):
    return x



def test_atom_attributor():
    model = _CountCarbonsModel()
    from .atom_attributor import attribute_atoms
    rslt = attribute_atoms("CCOCC",model,_identity)
    dst = abs(sum(rslt) - 4)
    print("org",rslt)
    assert dst < 0.3

def test_vec_atom_attributor():
    model = _CountCarbonsModel()
    from .atom_attributor_vec import attribute_atoms
    rslt = attribute_atoms(["CCOCC",],model,_identity)
    print("vec",rslt)
    rslt = rslt[0]
    dst = abs(sum(rslt) - 4)
    assert dst < 0.3