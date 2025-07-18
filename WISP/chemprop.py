import chemprop
import random

import pandas as pd
import matplotlib.pyplot as plt
from lightning import pytorch as pl
from pathlib import Path

from lightning.pytorch import seed_everything

from chemprop import data, featurizers, models, nn
from chemprop.nn import metrics
from chemprop.models import multi

from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error, mean_squared_error

import pdb
import numpy as np

from lightning.pytorch.callbacks import ModelCheckpoint

import torch

import os
from glob import glob

#to not get the output
import logging
logging.getLogger("lightning.pytorch").setLevel(logging.ERROR)

random.seed(6)
np.random.seed(6)
torch.manual_seed(6)
torch.use_deterministic_algorithms(True)

class SklChemprop:

    def __init__(self,problem_type:str,max_epochs:int,Smiles_Column_Name:str, Target_Column_Name:str, working_dir:str):
        # TODO: Hyperparameters should be passed here
        self.problem_type = problem_type
        self.max_epochs = max_epochs
        self.Smiles_Column_Name = Smiles_Column_Name
        self.Target_Column_Name = Target_Column_Name
        self.working_dir = working_dir

        self.mpnn = None
        self._load_latest_checkpoint()

    def _load_latest_checkpoint(self):
        checkpoint_dir = self.working_dir + 'checkpoints/'
        if not os.path.exists(checkpoint_dir):#skipps loading if there is no checkpoint dir
            return
        ckpt_files = glob(os.path.join(checkpoint_dir, 'best*.ckpt'))
        if not ckpt_files:#skipps loading if there is no checkpoint file
            return
        # Get most recent checkpoint and set as self model
        latest_ckpt = max(ckpt_files, key=os.path.getmtime)
        self.mpnn = models.MPNN.load_from_checkpoint(latest_ckpt) 

    def _df_to_loader(self,df):
        md = [data.MoleculeDatapoint.from_smi(smis, np.array(y).reshape(1,1)) for smis, y in zip(df.smiles.tolist(), df.y.tolist())]
        #ds = [data.MoleculeDataset(md[0][i], self.featurizer) for i in range(len(df))]
        ds = data.MoleculeDataset(md, self.featurizer)
        loader = data.build_dataloader(ds,shuffle=True,)
        return loader

    def fit(self,df_input:"pd.Dataframe"):
        df_input.reset_index(inplace=True)
        num_workers = 0 # number of workers for dataloader. 0 means using main process for data loading
        smiles_column = self.Smiles_Column_Name # name of the column containing SMILES strings
        target_columns = [self.Target_Column_Name] # list of names of the columns containing targets
        smis = df_input.loc[:, smiles_column].values
        ys = df_input.loc[:, target_columns].values
        
        all_data = [data.MoleculeDatapoint.from_smi(smi, y) for smi, y in zip(smis, ys)]
        
        #train validation split
        nr_val_samples = round(len(df_input) / 5) # 80/20 split
        val = df_input.sample(n=nr_val_samples, random_state=6)
        train_set = df_input.drop(val.index)

        # get the indices/data
        train_indices = [train_set.index.tolist()]
        val_indices = [val.index.tolist()]
        train_data, val_data, _ = data.split_data_by_indices(all_data, train_indices, val_indices)

        #preprocess data
        featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()

        train_dset = data.MoleculeDataset(train_data[0], featurizer)
        scaler = train_dset.normalize_targets()
        val_dset = data.MoleculeDataset(val_data[0], featurizer)
        val_dset.normalize_targets(scaler)

        train_loader = data.build_dataloader(train_dset, num_workers=num_workers, shuffle=False, seed=6)
        val_loader = data.build_dataloader(val_dset, num_workers=num_workers, shuffle=False, seed=6)

        mp = nn.BondMessagePassing()
        agg = nn.MeanAggregation()

        output_transform = nn.UnscaleTransform.from_standard_scaler(scaler)

        ffn = nn.RegressionFFN(output_transform=output_transform)

        batch_norm = True
        metric_list = [metrics.MAE()]#[nn.metrics.RMSE(), nn.metrics.R2score()]

        seed_everything(6, workers=True)
        mpnn = models.MPNN(mp, agg, ffn, batch_norm, metric_list)

        # Checkpointing
        checkpointing = ModelCheckpoint(
        self.working_dir + "checkpoints",  # Directory where model checkpoints will be saved
        "best-{epoch}-{val_loss:.2f}",  # Filename format for checkpoints, including epoch and validation loss
        "val_loss",  # Metric used to select the best checkpoint (based on validation loss)
        mode="min",  # Save the checkpoint with the lowest validation loss (minimization objective)
        save_last=True,  # Always save the most recent checkpoint, even if it's not the best
        )

        trainer = pl.Trainer(
            logger=False,
            enable_checkpointing=True, # Use `True` if you want to save model checkpoints. The checkpoints will be saved in the `checkpoints` folder.
            enable_progress_bar=False,
            accelerator="cpu",
            devices=1,
            max_epochs=self.max_epochs, # number of epochs to train for
            callbacks=[checkpointing], # Use the configured checkpoint callback
            deterministic=True,
        )

        trainer.fit(mpnn, train_loader, val_loader)
        self._load_latest_checkpoint()

    
    def predict(self, smiles):#df_input:"pd.Dataframe"
        num_workers = 0 # number of workers for dataloader. 0 means using main process for data loading
        #smiles_column = self.Smiles_Column_Name # name of the column containing SMILES strings
        #smis = df_input.loc[:, smiles_column].values
        smis = [smi for sublist in smiles for smi in sublist]
        
        all_data = [data.MoleculeDatapoint.from_smi(smi) for smi in smis]

        #preprocess data
        featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()
        test_dset = data.MoleculeDataset(all_data, featurizer)
        test_loader = data.build_dataloader(test_dset, num_workers=num_workers, shuffle=False, seed=6)

        #loader etc from previous
        with torch.inference_mode():
            trainer = pl.Trainer(
                logger=False,
                enable_progress_bar=False,
                accelerator="cpu",
                devices=1,
                deterministic=True
            )
            test_preds = trainer.predict(self.mpnn, test_loader)
        test_preds = np.concatenate(test_preds, axis=0)

        return test_preds
    
def get_features(data, CLOUMS):
    """
    Preprocesses the features according to theyr type.

    Keyword arguments:
    -- data: Table where the features are stored.
    -- CLOUMS: Colum where the feature is stored.

    Returns:
    -- X: Preprocessed feature.
    """
    features = []

    for i in CLOUMS:
        if type(data[i].values[0]) is np.ndarray:
            x = np.stack(data[i].values)
        else:
            x = data[i].values.reshape(-1, 1)
        features.append(x)

    X = np.concatenate(features, axis=1)
    return X
    
def train_GNN(train, Smiles_Column_Name, Target_Column_Name, working_dir):
    torch.manual_seed(6)
    model_GNN = SklChemprop(problem_type="regression", max_epochs=50, Smiles_Column_Name=Smiles_Column_Name, Target_Column_Name=Target_Column_Name, working_dir=working_dir)
    model_GNN.fit(train)
    
    prep_smiles = get_features(train, [Smiles_Column_Name])

    predictions = model_GNN.predict(prep_smiles)
    target = train[Target_Column_Name].values

    r2 = np.corrcoef(target.flatten(), predictions.flatten())[0,1]**2
    RMSE_z = np.sqrt(np.mean((target.flatten() - predictions.flatten())**2))
    RMSE_n = np.sqrt(np.mean((target.flatten() - np.mean(target.flatten()))**2))
    R2 = 1 - RMSE_z**2/RMSE_n**2
    MAE = mean_absolute_error(target.flatten(), predictions.flatten())
    mse = mean_squared_error(target.flatten(), predictions.flatten())
    rmse = np.sqrt(mse)

    return r2, R2, MAE, rmse

def identity(x):
    return x