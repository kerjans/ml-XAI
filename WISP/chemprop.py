import chemprop

import pandas as pd
import matplotlib.pyplot as plt
from lightning import pytorch as pl
from pathlib import Path

from chemprop import data, featurizers, models, nn
from chemprop.nn import metrics
from chemprop.models import multi

from sklearn.preprocessing import StandardScaler

import pdb
import numpy as np

from lightning.pytorch.callbacks import ModelCheckpoint

import torch

import os
from glob import glob


class SklChemprop:

    def __init__(self,problem_type:str,max_epochs:int,Smiles_Column_Name:str, Target_Column_Name:str):
        # TODO: Hyperparameters should be passed here
        self.problem_type = problem_type
        self.max_epochs = max_epochs
        self.Smiles_Column_Name = Smiles_Column_Name
        self.Target_Column_Name = Target_Column_Name

    def _df_to_loader(self,df):
        md = [data.MoleculeDatapoint.from_smi(smis, np.array(y).reshape(1,1)) for smis, y in zip(df.smiles.tolist(), df.y.tolist())]
        #ds = [data.MoleculeDataset(md[0][i], self.featurizer) for i in range(len(df))]
        ds = data.MoleculeDataset(md, self.featurizer)
        loader = data.build_dataloader(ds,shuffle=True,)
        return loader

    def fit(self,df_input:"pd.Dataframe"):
        num_workers = 0 # number of workers for dataloader. 0 means using main process for data loading
        smiles_column = self.Smiles_Column_Name # name of the column containing SMILES strings
        target_columns = [self.Target_Column_Name] # list of names of the columns containing targets
        smis = df_input.loc[:, smiles_column].values
        ys = df_input.loc[:, target_columns].values
        all_data = [data.MoleculeDatapoint.from_smi(smi, y) for smi, y in zip(smis, ys)]
        
        mols = [d.mol for d in all_data]  # RDkit Mol objects are use for structure based splits
        train_indices, val_indices, test_indices = data.make_split_indices(mols, "random", (0.8, 0.1, 0.1))  # unpack the tuple into three separate lists
        train_data, val_data, self.test_data = data.split_data_by_indices(
            all_data, train_indices, val_indices, test_indices
        )

        featurizer = featurizers.SimpleMoleculeMolGraphFeaturizer()

        train_dset = data.MoleculeDataset(train_data[0], featurizer)
        scaler = train_dset.normalize_targets()
        val_dset = data.MoleculeDataset(val_data[0], featurizer)
        val_dset.normalize_targets(scaler)
        test_dset = data.MoleculeDataset(self.test_data[0], featurizer)

        train_loader = data.build_dataloader(train_dset, num_workers=num_workers)
        val_loader = data.build_dataloader(val_dset, num_workers=num_workers, shuffle=False)
        self.test_loader = data.build_dataloader(test_dset, num_workers=num_workers, shuffle=False)

        mp = nn.BondMessagePassing()
        agg = nn.MeanAggregation()

        output_transform = nn.UnscaleTransform.from_standard_scaler(scaler)

        ffn = nn.RegressionFFN(output_transform=output_transform)

        batch_norm = True
        metric_list = [metrics.R2Score(), metrics.MAE()]#[nn.metrics.RMSE(), nn.metrics.MAE()]
        mpnn = models.MPNN(mp, agg, ffn, batch_norm, metric_list)

        # Checkpointing
        checkpointing = ModelCheckpoint(
        "/Users/kerrinjanssen/Nextcloud/PhD/Bayer/XAI/test-GNN/checkpoints",  # Directory where model checkpoints will be saved
        "best-{epoch}-{val_loss:.2f}",  # Filename format for checkpoints, including epoch and validation loss
        "val_loss",  # Metric used to select the best checkpoint (based on validation loss)
        mode="min",  # Save the checkpoint with the lowest validation loss (minimization objective)
        save_last=True,  # Always save the most recent checkpoint, even if it's not the best
        )

        trainer = pl.Trainer(
            logger=False,
            enable_checkpointing=True, # Use `True` if you want to save model checkpoints. The checkpoints will be saved in the `checkpoints` folder.
            enable_progress_bar=True,
            accelerator="auto",
            devices=1,
            max_epochs=20, # number of epochs to train for
            callbacks=[checkpointing], # Use the configured checkpoint callback
        )


        trainer.fit(mpnn, train_loader, val_loader)

    def predict(self,):
    
        # Get most recent checkpoint
        checkpoint_dir = '/Users/kerrinjanssen/Nextcloud/PhD/Bayer/XAI/test-GNN/checkpoints/'
        ckpt_files = glob(os.path.join(checkpoint_dir, '*.ckpt'))
        latest_ckpt = max(ckpt_files, key=os.path.getmtime)
        mpnn = models.MPNN.load_from_checkpoint(latest_ckpt)

        #loader etc from previous
        with torch.inference_mode():
            trainer = pl.Trainer(
                logger=None,
                enable_progress_bar=True,
                accelerator="cpu",
                devices=1
            )
            test_preds = trainer.predict(mpnn, self.test_loader)
        test_preds = np.concatenate(test_preds, axis=0)

        return test_preds, self.test_data