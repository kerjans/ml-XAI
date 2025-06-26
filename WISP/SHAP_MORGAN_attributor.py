from WISP.ml_helper import *

from rdkit import Chem
from rdkit.Chem import AllChem

import shap

def calculate_atom_weights(mol, coefficients_total, bit_info):
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
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    Chem.SanitizeMol(mol)#to keep the explicit hydrogens

    bit_info = {}
    #fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048, bitInfo=bit_info)
    generator = GetMorganGenerator(radius=2, fpSize=2048)
    ao = AllChem.AdditionalOutput()
    ao.CollectBitInfoMap()
    fp = generator.GetFingerprint(mol, additionalOutput=ao)
    bit_info = ao.GetBitInfoMap()

    atom_weights = calculate_atom_weights(mol, coefficients_total, bit_info)

    return atom_weights

def get_SHAP_Morgan_attributions(data, feature_column, smiles_column, model, explainer):
    prep_data = get_features(data, [feature_column])

    shap_values = explainer.shap_values(model.named_steps['scaler'].transform(prep_data))
    print("Shap values are calculated.")

    atom_weights_list = []

    for nr, (index, row) in enumerate(data.iterrows(), 1):
        smiles = row[smiles_column]
        shap_nr = nr - 1
        atom_weights = weights_morgan(smiles, shap_values[shap_nr])
        atom_weights_values = list(atom_weights.values())
        atom_weights_list.append(atom_weights_values)

    data['SHAP Attributions'] = atom_weights_list

    return data

def pick_shap_explainer(model):
    model_type = model.named_steps['model'].__class__.__name__
    if model_type in ['GradientBoostingRegressor', 'RandomForestRegressor']:
        explainer = shap.TreeExplainer(model.named_steps['model'])
    if model_type in ['MLPRegressor', 'SVR', 'GaussianProcessRegressor'] :
        prep_data = get_features(data, ['Morgan_Fingerprint 2048Bit 2rad'])
        explainer = shap.KernelExplainer(model.predict, model.named_steps['scaler'].transform(prep_data))
    if model_type in ['BayesianRidge', 'Lasso', 'LinearRegression'] :
        prep_data = get_features(data, ['Morgan_Fingerprint 2048Bit 2rad'])
        explainer = shap.LinearExplainer(model.named_steps['model'],model.named_steps['scaler'].transform(prep_data))
    if model_type in ['ProbabilisticRandomForest']:#'RandomForestClassifier', 
        explainer = shap.TreeExplainer(model.named_steps['model'].model, model_output="probability")    #does not work
    return explainer