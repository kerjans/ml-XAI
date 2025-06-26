"""
As example, we provide a demo scikit-learn model trained on a public dataset of lipophilicity (it's the MoleculeNet
reference dataset that can be downloaded from
https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/Lipophilicity.csv
We define a simple Morgan fingerprint embedder class that is a valid embedding class for computing input features for our
model.
"""
from pathlib import Path
import pickle
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs
from sklearn.svm import SVR
import urllib.request
from xai_selfies import SOL, MODEL, DATA_DIR
from xai_selfies.main import get_attributions_df


def read_data():
    if Path.exists(SOL):
        sol_df = pd.read_csv(SOL)
    else:
        print('Downloading solubility dataset for demo purposes...')
        urllib.request.urlretrieve("https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/Lipophilicity.csv", SOL)
        sol_df = pd.read_csv(SOL)
    return sol_df


class RDKitEmbedder:
    def __init__(self, radius, folding):
        self.radius = radius
        self.dim = folding

    def encode(self, smiles):
        mols = [Chem.MolFromSmiles(smi) for smi in smiles]
        valid_mols = [mol for mol in mols if mol is not None]
        fps = [AllChem.GetMorganFingerprintAsBitVect(mol, self.radius, nBits=self.dim) for mol in valid_mols]
        X = []
        for i, fp in enumerate(fps):
            array = np.zeros((0, ), dtype=np.int8)
            DataStructs.ConvertToNumpyArray(fp, array)
            X.append(array)
        return np.array(X)


class Model:
    def __init__(self):
        self.trained_model = get_demo_model()

    def predict(self, X):
        return self.trained_model.predict(X)


def get_demo_model():
    if Path.exists(MODEL):
        with open(MODEL, 'rb') as reader:
            return pickle.load(reader)
    else:
        print('Training a simple SVM model on the demo data...')
        df = read_data()
        smiles = df['smiles'].tolist()
        embedder = RDKitEmbedder(3, 1024)
        X = embedder.encode(smiles)
        y = np.array(df['exp'].tolist())
        model = SVR()  # in 5-fold CV we get R2 of [0.5523648  0.51164119 0.4997524  0.51903977 0.51958938]
        model.fit(X, y)
        with open(MODEL, 'wb') as writer:
            pickle.dump(model, writer, protocol=pickle.HIGHEST_PROTOCOL)
        return model


if __name__ == '__main__':
    data_df = read_data()
    print(data_df)
    molids = data_df['CMPD_CHEMBLID'].tolist()[:200]
    molecules = data_df['smiles'].tolist()[:200]
    embedder = RDKitEmbedder(3, 1024)
    trained_model = Model()
    attributions_df = get_attributions_df(molecules, trained_model, embedder, ids=molids, similarity=0.)
    attributions_df['logd'] = data_df['exp'].tolist()[:200]
    attributions_df.to_csv(Path(DATA_DIR) / 'demo_cime_input_logd.csv')
    print(attributions_df)
