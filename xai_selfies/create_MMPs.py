import pandas as pd
import subprocess
from rdkit import Chem

def count_atoms_in_smarts(smarts):
    mol = Chem.MolFromSmarts(smarts)        
    num_atoms = mol.GetNumAtoms()
    return num_atoms

def create_MMP_database(smi_data_dir, output_dir,data_to_add_information_from, added_colums):
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

    #remove duplicate MMPs: only keep the one where the constant part has the most atoms (the smallest transformation)
    MMP_data['constant_atom_count'] = MMP_data['constant'].apply(count_atoms_in_smarts)
    MMP_data = MMP_data.sort_values(by=['ID_1', 'ID_2', 'constant_atom_count'], ascending=[True, True, False])
    data_MMPs_filtered = MMP_data.drop_duplicates(subset=['ID_1', 'ID_2'], keep='first')


    #add target data
    for i in added_colums:
        data_MMPs_filtered[i + '_1'] = data_MMPs_filtered['ID_1'].map(data_to_add_information_from.set_index('ID')[i])
        data_MMPs_filtered[i + '_2'] = data_MMPs_filtered['ID_2'].map(data_to_add_information_from.set_index('ID')[i])

    return data_MMPs_filtered

def get_explicit_smiles(row, ed_or_pro):
    transformation = row['transformation']
    prod = transformation.split(">>")[ed_or_pro]
    mol_prod = Chem.MolFromSmiles(prod)
    mol_prod = Chem.AddHs(mol_prod)#to add hydrogens in the explicit part
    
    constant = row['constant']
    mol_constant = Chem.MolFromSmiles(constant)
    
    combo = Chem.CombineMols(mol_prod, mol_constant)
    
    dummy_atoms = [atom.GetIdx() for atom in combo.GetAtoms() if atom.GetAtomicNum() == 0]

    # Create a bond between the dummy atoms' neighbors
    neighbors = [list(combo.GetAtomWithIdx(idx).GetNeighbors())[0].GetIdx() for idx in dummy_atoms]

    editable = Chem.EditableMol(combo)
    editable.AddBond(neighbors[0], neighbors[1], order=Chem.BondType.SINGLE)
    
    # Remove dummy atoms
    for idx in sorted(dummy_atoms, reverse=True):  # reverse to avoid reindexing issues
        editable.RemoveAtom(idx)
    
    # Finalize molecule
    mol = editable.GetMol()
    Chem.SanitizeMol(mol)
    return Chem.MolToSmiles(mol)