from rdkit import Chem
from rdkit.Chem.Draw import SimilarityMaps
import numpy as np

def get_pred(fp, pred_function):
    fp = np.array([list(fp)])
    return pred_function(fp)[0]

def RDKit_attributor(smiles,fpFunction,model):
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    Chem.SanitizeMol(mol)#to keep the explicit hydrogens
    attributions = Chem.Draw.SimilarityMaps.GetAtomicWeightsForModel(mol, fpFunction,lambda x : get_pred(x, model.predict))
    return attributions