import random
import numpy as np
np.seterr(divide='ignore', invalid='ignore')

import pandas as pd
import matplotlib.pyplot as plt

import torch
from lightning import pytorch as pl
from lightning.pytorch import seed_everything
from lightning.pytorch.callbacks import ModelCheckpoint

from pathlib import Path

import chemprop
from chemprop import data, featurizers, models, nn
from chemprop.nn import metrics
from chemprop.models import multi

from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_absolute_error, mean_squared_error

import os
from glob import glob

import logging
logging.getLogger("lightning.pytorch").setLevel(logging.ERROR)

random.seed(6)
np.random.seed(6)
torch.manual_seed(6)
torch.use_deterministic_algorithms(True)

class SklChemprop:
    """
    Wrapper around a Chemprop MPNN model for regression tasks
    using PyTorch Lightning. Provides methods to train from a pandas DataFrame of
    SMILES and targets, load the best checkpoint, and predict on new data.

    Attributes:
        problem_type (str): "regression" or "classification".
        max_epochs (int): Maximum number of training epochs.
        Smiles_Column_Name (str): Name of the SMILES column in the DataFrame.
        Target_Column_Name (str): Name of the target column in the DataFrame.
        working_dir (str): Directory where checkpoints and logs are stored.
        mpnn (Chemprop model): The loaded or newly trained MPNN model.
    """
    def __init__(self,problem_type:str,max_epochs:int,Smiles_Column_Name:str, Target_Column_Name:str, working_dir:str):
        self.problem_type = problem_type
        self.max_epochs = max_epochs
        self.Smiles_Column_Name = Smiles_Column_Name
        self.Target_Column_Name = Target_Column_Name
        self.working_dir = working_dir

        self.mpnn = None
        self._load_latest_checkpoint()

    def _load_latest_checkpoint(self):
        """
        Search the `working_dir/checkpoints/` folder for the most recently modified
        'best*.ckpt' file and load it into `self.mpnn`. If none exists, does nothing.
        """
        checkpoint_dir = self.working_dir + 'checkpoints/'
        if not os.path.exists(checkpoint_dir):#skips loading if there is no checkpoint dir
            return
        ckpt_files = glob(os.path.join(checkpoint_dir, 'best*.ckpt'))
        if not ckpt_files:#skips loading if there is no checkpoint file
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
        """
        Train the MPNN on the given DataFrame.

        Steps:
          1. Reset index and split into train/validation (80/20).
          2. Featurize molecules and normalize targets.
          3. Instantiate the MPNN, set up a PyTorch Lightning Trainer with checkpointing.
          4. Train for up to `self.max_epochs` epochs.
          5. Reload the best checkpoint into `self.mpnn`.

        Parameters:
            df_input (pd.DataFrame): Contains SMILES and target columns.
        """
        df_input.reset_index(inplace=True)
        num_workers = 0 # number of workers for dataloader. 0 means using main process for data loading
        smiles_column = self.Smiles_Column_Name 
        target_columns = [self.Target_Column_Name] 
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

    
    def predict(self, smiles):
        """
        Predict target values for a batch of SMILES strings.

        Parameters:
            smiles (list of str or nested list): SMILES to predict on.

        Returns:
            np.ndarray: Concatenated array of predictions from the loaded MPNN.
        """
        num_workers = 0 
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
    
def get_features(data, CLOUMNS):
    """
    Preprocesses the features according to theyr type.

    Keyword arguments:
    -- data: Table where the features are stored.
    -- CLOUMNS: Column where the feature is stored.

    Returns:
    -- X: Preprocessed feature.
    """
    features = []

    for i in CLOUMNS:
        if type(data[i].values[0]) is np.ndarray:
            x = np.stack(data[i].values)
        else:
            x = data[i].values.reshape(-1, 1)
        features.append(x)

    X = np.concatenate(features, axis=1)
    return X
    
def train_GNN(train, Smiles_Column_Name, Target_Column_Name, working_dir):
    """
    Train a Chemprop graph‐neural network on a regression task and evaluate on the same split.

    This function:
      1. Sets a fixed random seed for reproducibility.
      2. Instantiates SklChemprop with the given SMILES and target column names,
         training for up to 50 epochs.
      3. Fits the model on the `train` DataFrame.
      4. Uses `get_features(train, [Smiles_Column_Name])` as input to predict on `train`.
      5. Computes:
         - r2: squared Pearson correlation between true and predicted values,
         - R2: alternate R² via 1 − (σ_residual²/σ_total²),
         - MAE: mean absolute error,
         - rmse: root mean squared error.
      6. Returns these four metrics.

    Parameters:
        train (pd.DataFrame):
            DataFrame containing columns `Smiles_Column_Name` and `Target_Column_Name`.
        Smiles_Column_Name (str):
            Name of the column with SMILES strings in `train`.
        Target_Column_Name (str):
            Name of the numeric target column in `train`.
        working_dir (str):
            Directory where Chemprop checkpoints and logs will be written.

    Returns:
        tuple:
            r2 (float): Squared Pearson correlation coefficient.
            R2 (float): Coefficient of determination via MSE ratio.
            MAE (float): Mean absolute error.
            rmse (float): Root mean squared error.
    """
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