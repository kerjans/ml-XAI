from WISP.ml_helper import *

from rdkit import Chem
from rdkit.Chem import AllChem

import shap

def calculate_atom_weights(mol, coefficients_total, bit_info):
    """
    Compute per-atom weight contributions from fingerprint bits and distribute
    each bit’s importance to the central atom and its neighbors up to radius 2.

    For each bit index `b` in `bit_info`:
      • bit_info[b] is a list of (atom_idx, radius) pairs.
      • coefficients_total[b] is the coefficient to add.
      • radius = 0 → add only to atom_idx.
      • radius = 1 → add to atom_idx + its immediate neighbors.
      • radius = 2 → add to atom_idx + its neighbors + second‐order neighbors
        (excluding the original atom).

    Parameters:
        mol (rdkit.Chem.Mol):
            An RDKit molecule whose atoms will be weighted.
        coefficients_total (dict[int, float]):
            Mapping from bit index to a importance value.
        bit_info (dict[int, list[tuple[int, int]]]):
            Mapping from bit index to a list of (atom index, radius) entries.

    Returns:
        dict[int, float]:
            A mapping from each atom index to its accumulated weight.
    """
    atom_weights = {atom.GetIdx(): 0.0 for atom in mol.GetAtoms()}

    for x in bit_info:
        for substructure_atoms in bit_info[x]:
            central = mol.GetAtomWithIdx(substructure_atoms[0])
            if substructure_atoms[1] == 0:
                atom_weights[substructure_atoms[0]] += coefficients_total[x]
            if substructure_atoms[1] == 1:
                atom_weights[substructure_atoms[0]] += coefficients_total[x]
                surr = [neighbor.GetIdx() for neighbor in central.GetNeighbors()]
                for neig_atoms in surr:
                    atom_weights[neig_atoms] += coefficients_total[x]
            if substructure_atoms[1] == 2:
                atom_weights[substructure_atoms[0]] += coefficients_total[x]
                surr = [neighbor.GetIdx() for neighbor in central.GetNeighbors()]
                for neig_atoms in surr:
                    atom_weights[neig_atoms] += coefficients_total[x]
                    neig_atoms_indx = mol.GetAtomWithIdx(neig_atoms)
                    surr_second = [neighbor_second.GetIdx() for neighbor_second in neig_atoms_indx.GetNeighbors()]
                    for sec_neighbor in surr_second:
                        if sec_neighbor != substructure_atoms[0]:
                            atom_weights[sec_neighbor] += coefficients_total[x]

    return atom_weights

def weights_morgan(smiles, coefficients_total):
    """
    Compute per‐atom weights from a SMILES string using a radius‐2 Morgan fingerprint.

    Steps:
      1. Parse the SMILES into an RDKit Mol (without automatic sanitization),
         then explicitly sanitize to preserve explicit H’s.
      2. Generate a 2048‐bit Morgan fingerprint with AdditionalOutput enabled
         to collect bit‐to‐atom mapping (bit_info).
      3. Delegate to `calculate_atom_weights` to distribute each bit’s coefficient
         to its central atom and neighbors up to radius 2.

    Parameters:
        smiles (str):
            The input molecule in SMILES format.
        coefficients_total (dict[int, float]):
            Mapping from fingerprint bit index to its coefficient.

    Returns:
        dict[int, float]:
            Atom‐indexed weights (atom idx → accumulated coefficient).

    Dependencies:
        from rdkit import Chem, AllChem
        from rdkit.Chem.AllChem import GetMorganGenerator
        calculate_atom_weights

    """
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    Chem.SanitizeMol(mol)

    generator = GetMorganGenerator(radius=2, fpSize=2048)
    ao = AllChem.AdditionalOutput()
    ao.CollectBitInfoMap()
    _ = generator.GetFingerprint(mol, additionalOutput=ao)
    bit_info = ao.GetBitInfoMap()

    atom_weights = calculate_atom_weights(mol, coefficients_total, bit_info)

    return atom_weights

def get_SHAP_Morgan_attributions(data, feature_column, smiles_column, model, explainer):
    """
    Compute per‐atom SHAP attributions for each molecule using a Morgan fingerprint explainer.

    Parameters:
        data (pd.DataFrame):
            Input DataFrame containing at least `feature_column` and `smiles_column`.
        feature_column (str):
            Name of the column whose values were used to fit the SHAP explainer.
        smiles_column (str):
            Name of the column with the SMILES strings for atom‐level mapping.
        model (Pipeline):
            A fitted sklearn Pipeline with a 'scaler' step and a downstream predictor.
        explainer (shap.Explainer):
            A SHAP explainer whose .shap_values() accepts scaled feature arrays.

    Returns:
        pd.DataFrame:
            The original DataFrame augmented with:
              - 'SHAP Attributions': a list of atom‐weight arrays for each row.
    """
    prep_data = get_features(data, [feature_column])

    shap_values = explainer.shap_values(model.named_steps['scaler'].transform(prep_data))
    print("SHAP values are calculated.")

    atom_weights_list = []

    for nr, (_, row) in enumerate(data.iterrows(), 1):
        smiles = row[smiles_column]
        shap_nr = nr - 1
        atom_weights = weights_morgan(smiles, shap_values[shap_nr])
        atom_weights_values = list(atom_weights.values())
        atom_weights_list.append(atom_weights_values)

    data['SHAP Attributions'] = atom_weights_list

    return data

def pick_shap_explainer(model, data):
    """
    Choose and construct a SHAP explainer appropriate for the fitted model in a Pipeline.

    Parameters:
        model (Pipeline): A fitted sklearn Pipeline with steps [('scaler', …), ('model', estimator)].
        background_data (pd.DataFrame): DataFrame to sample background for Kernel/Linear explainers.
        feature_columns (list of str): Names of columns used to build X for explainer (e.g. Morgan bits).

    Returns:
        shap.Explainer: A SHAP explainer instance suited to the model type.
    """
    model_type = model.named_steps['model'].__class__.__name__
    if model_type in ['GradientBoostingRegressor', 'RandomForestRegressor']:
        explainer = shap.TreeExplainer(model.named_steps['model'])
    if model_type in ['MLPRegressor', 'SVR', 'GaussianProcessRegressor'] :
        prep_data = get_features(data, ['Morgan_Fingerprint 2048Bit 2rad'])
        explainer = shap.KernelExplainer(model.predict, model.named_steps['scaler'].transform(prep_data))
    if model_type in ['BayesianRidge', 'Lasso', 'LinearRegression'] :
        prep_data = get_features(data, ['Morgan_Fingerprint 2048Bit 2rad'])
        explainer = shap.LinearExplainer(model.named_steps['model'],model.named_steps['scaler'].transform(prep_data))
    return explainer