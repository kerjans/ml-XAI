from rdkit import Chem
import re

def get_unmatched_atom_indices_from_mol(mol, constant_smarts):
    """
    Identify atoms in a molecule that are not part of a “constant” substructure.

    This function:
      1. Removes any dummy‐atom labels of the form “[*:n]” from the SMARTS string.
      2. Splits the cleaned SMARTS on “.” to handle multiple fragments.
      3. For each fragment:
         a. Finds all substructure matches in `mol`.
         b. If multiple matches exist, selects the one where exactly one bond
            connects the matched atoms to the rest of the molecule.
         c. Otherwise falls back to the first match.
      4. Collects the union of all matched atom indices.
      5. Returns all atom indices in `mol` that were not matched.

    Parameters:
        mol (rdkit.Chem.Mol):
            The RDKit molecule to search.
        constant_smarts (str):
            A SMARTS pattern for the constant fragment; may include dummy‐atom
            tags (e.g. “[*:1]”), which are stripped prior to matching.

    Returns:
        set of int:
            The indices of atoms in `mol` not included in any matched fragment.
    """
    smarts = re.sub(r'\[\*\:\d+\]', '', constant_smarts)
    
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
            match_atoms.update(matches[0])

    all_atoms = set(range(mol.GetNumAtoms()))
    unmatched = all_atoms - match_atoms
    return unmatched

def get_unmatched_atom_indices_fragments(smiles1, smiles2, constant_smarts):
    """
    For a pair of SMILES strings, identify which atoms are *not* part of the
    shared “constant” fragment defined by a SMARTS pattern.

    Steps:
      1. Parse each SMILES into an RDKit Mol without auto‐sanitization.
      2. Explicitly sanitize both molecules to preserve explicit hydrogens.
      3. Use `get_unmatched_atom_indices_from_mol` on each Mol to find the set
         of atom indices not matched by the constant SMARTS.
      4. Return a tuple of two sets: (unmatched_indices_mol1, unmatched_indices_mol2).

    Parameters:
        smiles1 (str): SMILES of the first molecule.
        smiles2 (str): SMILES of the second molecule.
        constant_smarts (str): SMARTS pattern for the constant fragment shared by both.

    Returns:
        (set[int], set[int]):
            A tuple where the first element is the set of atom indices in mol1
            not in the constant fragment, and the second is the same for mol2.

    """
    mol1 = Chem.MolFromSmiles(smiles1, sanitize=False)
    mol2 = Chem.MolFromSmiles(smiles2, sanitize=False)
    Chem.SanitizeMol(mol2)
    Chem.SanitizeMol(mol1)

    if mol1 is None or mol2 is None:
        raise ValueError("Invalid SMILES input")

    unmatched_mol1 = get_unmatched_atom_indices_from_mol(mol1, constant_smarts)
    unmatched_mol2 = get_unmatched_atom_indices_from_mol(mol2, constant_smarts)

    return unmatched_mol1, unmatched_mol2

def get_neighbors(smiles: str, atom_indices: list[int]) -> dict[int, list[int]]:
    """
    Given a SMILES string and a list of atom indices, return the set containing
    each specified index plus the indices of its immediate neighbors (excluding
    any that were already listed).

    Parameters:
        smiles (str): A SMILES representation of the molecule.
        atom_indices (list[int]): Atom indices to include and whose neighbors to add.

    Returns:
        set[int]: The union of `atom_indices` and all atom indices directly bonded
                  to those atoms (excluding any in the original list).
    """
    
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
    """
    Identify atoms *not* in a given list of indices for a molecule.

    Parameters:
        smiles (str):
            A SMILES string describing the molecule.
        selected_indices (iterable of int):
            Atom indices to exclude from the result.

    Returns:
        list of int:
            All atom indices in the molecule that are not in `selected_indices`.

    Notes:
        - Parses without sanitization, then explicitly sanitizes to retain
          explicit hydrogens.
        - If SMILES parsing fails (mol is None), calling SanitizeMol will error.
    """
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    Chem.SanitizeMol(mol)
    all_indices = list(range(mol.GetNumAtoms()))
    unselected_indices = [i for i in all_indices if i not in selected_indices]
    return unselected_indices