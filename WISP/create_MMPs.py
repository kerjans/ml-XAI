import pandas as pd
# This here is needed to suppress annoying warnings:
pd.options.mode.chained_assignment = None  # default='warn'
import subprocess
from rdkit import Chem

def count_atoms_in_smarts(smarts):
    """
    Parse a SMARTS pattern and return the number of atoms it contains.

    Parameters:
        smarts (str): A SMARTS string defining a substructure query.

    Returns:
        int: The total number of atoms in the parsed SMARTS pattern.

    """
    mol = Chem.MolFromSmarts(smarts)        
    num_atoms = mol.GetNumAtoms()
    return num_atoms

def create_MMP_database(smi_data_dir, output_dir,data_to_add_information_from, added_columns):
    """
    Generate a matched‐molecular‐pair (MMP) database from SMILES, filter by smallest transformation,
    and append specified columns from an external DataFrame.

    This function:
      1. Calls the `mmpdb` CLI to fragment and index the input SMILES file.
      2. Loads the resulting index CSV into a DataFrame with columns:
         ['smiles_1', 'smiles_2', 'ID_1', 'ID_2', 'transformation', 'constant'].
      3. Computes the atom count of the constant substructure and keeps the MMP
         with the largest constant (i.e. smallest variable part) per ID pair.
      4. Maps each specified column in `added_columns` from `data_to_add_information_from`
         (using its 'ID' column) to new columns `<col>_1` and `<col>_2` in the MMP DataFrame.

    Parameters:
        smi_data_dir (str):
            Path to the input SMILES file (TSV or CSV) for MMPdb fragmentation.
        output_dir (str):
            Directory where intermediate files ('fragments.fragdb', 'Index.csv') are written.
        data_to_add_information_from (pd.DataFrame):
            DataFrame containing an 'ID' column and the columns listed in `added_columns`.
        added_columns (list of str):
            Names of columns in `data_to_add_information_from` to attach to each MMP entry.

    Returns:
        pd.DataFrame:
            Filtered MMP database with columns:
              - 'smiles_1', 'smiles_2', 'ID_1', 'ID_2', 'transformation', 'constant'
              - 'constant_atom_count'
              - For each col in `added_columns`, '<col>_1' and '<col>_2'.
    """
    subprocess.run(["mmpdb", "fragment", 
                "--num-cuts", "1",
                smi_data_dir, "-o", output_dir + "fragments.fragdb"])
    
    subprocess.run(["mmpdb", "index", 
                "--max-variable-ratio", "0.2", 
                output_dir + "fragments.fragdb", 
                "-o", output_dir + "Index.csv",
                "--out", "csv"])
    
    column_headers = []
    column_headers.append("smiles_1")
    column_headers.append("smiles_2")
    column_headers.append("ID_1")
    column_headers.append("ID_2")
    column_headers.append("transformation")
    column_headers.append("constant")

    MMP_data = pd.read_csv(output_dir + "Index.csv", header=None, names=column_headers, sep='\t')

    MMP_data['constant_atom_count'] = MMP_data['constant'].apply(count_atoms_in_smarts)
    MMP_data = MMP_data.sort_values(by=['ID_1', 'ID_2', 'constant_atom_count'], ascending=[True, True, False])
    data_MMPs_filtered = MMP_data.drop_duplicates(subset=['ID_1', 'ID_2'], keep='first')

    for i in added_columns:
        data_MMPs_filtered[i + '_1'] = data_MMPs_filtered['ID_1'].map(data_to_add_information_from.set_index('ID')[i])
        data_MMPs_filtered[i + '_2'] = data_MMPs_filtered['ID_2'].map(data_to_add_information_from.set_index('ID')[i])

    return data_MMPs_filtered