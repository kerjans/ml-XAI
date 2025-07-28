from rdkit import Chem
import re

def get_unmatched_atom_indices_from_mol(mol, constant_smarts):

    smarts = re.sub(r'\[\*\:\d+\]', '', constant_smarts)#delete dummy atom completely
    
    fragments = smarts.split('.')

    match_atoms = set()

    for fragment in fragments:
    
        constant_mol = Chem.MolFromSmarts(fragment)
        matches = mol.GetSubstructMatches(constant_mol)
        
        if len(matches) >= 2:
            
            best_match = None
            for match in matches:
                match_set = set(match)
                external_bonds = 0

                for atom_idx in match:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    for neighbor in atom.GetNeighbors():
                        neighbor_idx = neighbor.GetIdx()
                        if neighbor_idx not in match_set:
                            external_bonds += 1

                if external_bonds == 1:
                    best_match = match
                    break
            match_atoms.update(best_match)
            
        elif matches:
            # Fall back to the first match if only one is found
            match_atoms.update(matches[0])

    all_atoms = set(range(mol.GetNumAtoms()))
    unmatched = all_atoms - match_atoms
    return unmatched

def get_unmatched_atom_indices_fragments(smiles1, smiles2, constant_smarts):
    mol1 = Chem.MolFromSmiles(smiles1, sanitize=False)
    mol2 = Chem.MolFromSmiles(smiles2, sanitize=False)
    Chem.SanitizeMol(mol2)#to keep the explicit hydrogens
    Chem.SanitizeMol(mol1)#to keep the explicit hydrogens

    if mol1 is None or mol2 is None:
        raise ValueError("Invalid SMILES input")

    # Get unmatched atoms for mol1 and mol2
    unmatched_mol1 = get_unmatched_atom_indices_from_mol(mol1, constant_smarts)
    unmatched_mol2 = get_unmatched_atom_indices_from_mol(mol2, constant_smarts)

    return unmatched_mol1, unmatched_mol2

def get_neighbors(smiles: str, atom_indices: list[int]) -> dict[int, list[int]]:
    
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    Chem.SanitizeMol(mol)#to keep the explicit hydrogens
    atom_set = set(atom_indices)
    unmatched_and_neighbors = set()

    for idx in atom_indices:
        unmatched_and_neighbors.add(idx)
        atom = mol.GetAtomWithIdx(idx)
        neighbor_indices = [
            nbr.GetIdx() for nbr in atom.GetNeighbors()
            if nbr.GetIdx() not in atom_set
        ]
        for index in neighbor_indices:
            unmatched_and_neighbors.add(index)

    return unmatched_and_neighbors

def get_unselected_atom_indices(smiles, selected_indices):
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    Chem.SanitizeMol(mol)#to keep the explicit hydrogens
    all_indices = list(range(mol.GetNumAtoms()))
    unselected_indices = [i for i in all_indices if i not in selected_indices]
    return unselected_indices