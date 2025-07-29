from rdkit import Chem
from rdkit.Chem.Draw import SimilarityMaps
import numpy as np

def get_pred(fp, pred_function):
    """
    Wrap a 1D fingerprint vector into a batch of size one, call the prediction
    function, and return the single prediction.

    Parameters:
        fp (array-like):
            A fingerprint or feature vector (e.g., list or 1D numpy array).
        pred_function (callable):
            A function or model method that accepts a 2D array of shape (n_samples, n_features)
            and returns an array-like of predictions.

    Returns:
        The first (and only) element of the prediction array for this single sample.
    """
    fp = np.array([list(fp)])
    return pred_function(fp)[0]

def RDKit_attributor(smiles, fpFunction, model):
    """
    Compute per‐atom attribution weights using RDKit’s SimilarityMaps.

    This function:
      1. Parses the SMILES into an RDKit Mol (preserving explicit hydrogens).
      2. Sanitizes the molecule.
      3. Calls RDKit’s GetAtomicWeightsForModel with:
         - mol: the RDKit molecule,
         - fpFunction: a function that generates a fingerprint feature array for a molecule,
         - a callback that wraps model.predict via get_pred.
      4. Returns the list of atomic weights produced by SimilarityMaps.

    Parameters:
        smiles (str): The input molecule as a SMILES string.
        fpFunction (callable): Function that takes an RDKit Mol and returns a fingerprint (iterable).
        model: A fitted model with a .predict(X) method accepting a 2D array.

    Returns:
        list of float: Attribution weight for each atom in the molecule.
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    Chem.SanitizeMol(mol)#to keep the explicit hydrogens
    attributions = Chem.Draw.SimilarityMaps.GetAtomicWeightsForModel(mol, fpFunction,lambda x: get_pred(x, model.predict))
    return attributions