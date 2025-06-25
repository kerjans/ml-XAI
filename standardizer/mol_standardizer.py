"""
A module for the standardization of molecules.

First, instantiate a standardizer instance:
>>> std = Standardizer(max_num_atoms=100,max_num_tautomers=10,include_stereoinfo=False,)

Then, call it on a given smiles to obtain a tuple of standardized smiles and error:
>>> std("CCCCCC[C@@](O)(F)C")
('CCCCCCC(C)(O)F', None)

In this case, no error was generated, hence the second entry of the tuple is None.
We can simply apply the standardization also on a dataframe:
>>> df_alkanes = pd.DataFrame({"smiles":["C(C)"*i for i in range(10)]})
>>> df_alkanes["smiles_std"] = df_alkanes["smiles"].apply(lambda smi: std(smi)[0]) 
>>> df_alkanes                          # doctest:+NORMALIZE_WHITESPACE
                                     smiles                        smiles_std
    0
    1                                  C(C)                                CC
    2                              C(C)C(C)                              CCCC
    3                          C(C)C(C)C(C)                          CCC(C)CC
    4                      C(C)C(C)C(C)C(C)                      CCC(C)C(C)CC
    5                  C(C)C(C)C(C)C(C)C(C)                  CCC(C)C(C)C(C)CC
    6              C(C)C(C)C(C)C(C)C(C)C(C)              CCC(C)C(C)C(C)C(C)CC
    7          C(C)C(C)C(C)C(C)C(C)C(C)C(C)          CCC(C)C(C)C(C)C(C)C(C)CC
    8      C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)      CCC(C)C(C)C(C)C(C)C(C)C(C)CC
    9  C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)C(C)  CCC(C)C(C)C(C)C(C)C(C)C(C)C(C)CC



"""
import sys


from math import nan
from pathlib import Path
import tempfile
from typing import Any, Tuple, Union
import uuid
import pandas as pd
from rdkit import rdBase
from rdkit.Chem.MolStandardize import rdMolStandardize
#from rdkit.Chem.MolStandardize.tautomer import TautomerTransform
from rdkit import Chem

from standardizer.config import _DIR_DATA
from standardizer.io import mol_to_image

rdBase.DisableLog("rdApp.*")
from rdkit.Chem import MolToSmiles, MolFromSmiles
from rdkit.Chem.AllChem import GetHashedMorganFingerprint, GetMorganFingerprint

import PIL
from PIL import ImageDraw
from PIL.Image import Image

# from .helper import *

# from .config import ConfigDict, SecretDict
import copy

import math
import random
import re
import numpy as np
from typing import List


class SmilesParseException(Exception):
    def __init__(
        self,
        smiles: str,
    ):
        super().__init__(f"Could not parse smiles: '{smiles}'")
        self.smiles = smiles


class TooManyAtomsException(Exception):
    def __init__(
        self,
        num_atoms_found: int,
        num_atoms_max: int,
        smiles: str,
    ):
        super().__init__(
            "number of atoms {0} exceeds limit of {1} for smiles {2}".format(
                num_atoms_found, num_atoms_max, smiles
            ),
        )
        self.num_atoms_found = num_atoms_found
        self.num_atoms_max = num_atoms_found
        self.smiles = smiles


class Standardizer(object):
    def __init__(
        self,
        max_num_atoms: int,
        max_num_tautomers: int,
        include_stereoinfo: bool,
        verbosity=0,
        canonicalize_tautomers: bool = None,
        keep_largest_fragment: bool = None,
        remove_defined_ions: bool = None,
        ions_sdf_file: Path = None,
        sanitize_mol: bool = None,
        disconnect_metals: bool = None,
        reionize: bool = None,
        normalize: bool = None,
    ):
        """
        Standardizes the given smiles.

        :param max_num_atoms: Maximum number of atoms until an error is raised.
        :param max_num_tautomers: The number of tautomers considered within the tautomer canonicalization.
        :param sanitize_mol: whether to apply RDKit sanitization.


        Examples:
        ----------
        >>> std = Standardizer(max_num_atoms=100,max_num_tautomers=10,include_stereoinfo=False,)
        >>> std("CCCCCC[C@@](O)(F)C")
        ('CCCCCCC(C)(O)F', None)
        >>> std = Standardizer(max_num_atoms=100,max_num_tautomers=10,include_stereoinfo=True,)
        >>> std("CCCCCC[C@@](O)(F)C")
        ('CCCCCC[C@@](C)(O)F', None)
        >>> long_mol = "C" * 101
        >>> std = Standardizer(max_num_atoms=100,max_num_tautomers=10,include_stereoinfo=False,)
        >>> std(long_mol)
        (None, TooManyAtomsException('number of atoms 101 exceeds limit of 100 for smiles CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC'))
        >>> type(std(long_mol)[-1])
        <class 'smal.mol_standardizer.TooManyAtomsException'>
        >>> std("this-is-not-a-proper-smiles")
        (None, SmilesParseException("Could not parse smiles: 'this-is-not-a-proper-smiles'"))
        >>> std(MolFromSmiles("CCCOCCCP(C)(CCC)"))
        ('CCCOCCCP(C)CCC', None)
        """
        if canonicalize_tautomers is None:
            canonicalize_tautomers = False
        if keep_largest_fragment is None:
            keep_largest_fragment = True
        if remove_defined_ions is None:
            remove_defined_ions = False
        if sanitize_mol is None:
            sanitize_mol = True
        if disconnect_metals is None:
            disconnect_metals = True
        if reionize is None:
            reionize = True
        if normalize is None:
            normalize = True
        self.sanitize_mol = sanitize_mol
        self.disconnect_metals = disconnect_metals
        self.reionize = reionize
        self.normalize = normalize
        self.max_num_atoms = max_num_atoms
        self.max_num_tautomers = max_num_tautomers
        self.include_stereoinfo = include_stereoinfo
        self.verbosity = verbosity

        self.canonicalize_tautomers = canonicalize_tautomers

        ## Load new tautomer enumarator/canonicalizer
        self.tautomerizer = rdMolStandardize.TautomerEnumerator()
        self.tautomerizer.SetMaxTautomers(self.max_num_tautomers)
        self.tautomerizer.SetRemoveSp3Stereo(
            False
        )  # Keep stereo information of keto/enol tautomerization

        self.remove_defined_ions = remove_defined_ions
        self.keep_largest_fragment = keep_largest_fragment
        if ions_sdf_file is None:
            ions_sdf_file = _DIR_DATA / "PIx5-trivials.sdf"
        self.ions_sdf_file = Path(ions_sdf_file)
        if self.remove_defined_ions:
            assert (
                ions_sdf_file.exists()
            ), "when removing defined ions then ions_sdf_file must be specified"

            self.ions = {
                Chem.CanonSmiles(
                    Chem.MolToSmiles(mol), useChiral=self.include_stereoinfo
                )
                for mol in Chem.SDMolSupplier(str(self.ions_sdf_file))
                if mol is not None
            }

    def __call__(self, smiles: Union[str, Chem.Mol]) -> Any:
        return self.calculate_single(smiles)

    # functions enabling pickle
    def __getstate__(self):
        return (
            self.max_num_atoms,
            self.max_num_tautomers,
            self.include_stereoinfo,
            self.verbosity,
        )

    def __setstate__(self, state):
        self.__init__(*state)

    @classmethod
    def from_param_dict(cls, method_param_dict, verbosity=0):
        return cls(**method_param_dict, verbosity=verbosity)

    def _my_standardizer(self, mol: Chem.Mol) -> Chem.Mol:
        mol = copy.deepcopy(mol)
        if self.sanitize_mol:
            Chem.SanitizeMol(mol)
        mol = Chem.RemoveHs(mol)

        if self.disconnect_metals:
            disconnector = rdMolStandardize.MetalDisconnector()
            mol = disconnector.Disconnect(mol)
        mol = Chem.MolToSmiles(mol)

        if self.remove_defined_ions:
            mol = Chem.MolFromSmiles(
                ".".join(
                    [
                        part
                        for part in mol.split(".")
                        if Chem.CanonSmiles(part) not in self.ions
                    ]
                )
            )

        if self.keep_largest_fragment:
            mol = rdMolStandardize.ChargeParent(Chem.MolFromSmiles(mol))

        if self.normalize:
            normalizer = rdMolStandardize.Normalizer()
            mol = normalizer.normalize(mol)
        if self.reionize:
            reionizer = rdMolStandardize.Reionizer()
            mol = reionizer.reionize(mol)
        Chem.AssignStereochemistry(mol, force=True, cleanIt=True)
        # TODO: Check this removes symmetric stereocenters
        return mol

    def _tautomerize(self, mol: Chem.Mol) -> Chem.Mol:
        return self.tautomerizer.Canonicalize(mol)

    @staticmethod
    def _isotope_parent(mol: Chem.Mol) -> Chem.Mol:
        mol = copy.deepcopy(mol)
        # Replace isotopes with common weight
        for atom in mol.GetAtoms():
            atom.SetIsotope(0)
        return mol

    def calculate_single(self, smiles) -> Tuple:
        if isinstance(smiles, str):
            try:
                mol = MolFromSmiles(
                    smiles
                )  # Read SMILES and convert it to RDKit mol object.
            except (TypeError, ValueError, AttributeError) as e:
                return None, e
        else:
            mol = smiles
        # Check, if the input SMILES has been converted into a mol object.
        if mol is None:
            return None, SmilesParseException(smiles)
        # check size of the molecule based on the non-hydrogen atom count.
        if mol.GetNumAtoms() >= self.max_num_atoms:
            return (
                None,
                TooManyAtomsException(mol.GetNumAtoms(), self.max_num_atoms, smiles),
            )
        try:
            assert bool(self.keep_largest_fragment) != bool(
                self.remove_defined_ions
            ), "set either or"
            if self.keep_largest_fragment:
                mol = rdMolStandardize.ChargeParent(mol)

            mol = self._isotope_parent(mol)
            if self.include_stereoinfo is False:
                Chem.RemoveStereochemistry(mol)
            if self.canonicalize_tautomers:
                mol = self._tautomerize(mol)
            mol_clean_tmp = self._my_standardizer(mol)
            smi_clean_tmp = MolToSmiles(
                mol_clean_tmp
            )  # convert mol object back to SMILES
            ## Double check if standardized SMILES is a valid mol object
            mol_clean = MolFromSmiles(smi_clean_tmp)
            smi_clean = MolToSmiles(mol_clean)
        except (TypeError, ValueError, AttributeError) as e:
            return None, e
        return smi_clean, None

