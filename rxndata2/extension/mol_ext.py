from rdkit import Chem


def remove_atom_mapping(mol):
    [atom.SetAtomMapNum(0) for atom in mol.GetAtoms()]
    smiles = Chem.MolToSmiles(mol)
    new_mol = Chem.MolFromSmiles(smiles)
    # new_mol = Chem.RemoveHs(mol)  # Could the hydrogen atoms be labeled? For some reason, H atoms cannot be removed.
    return new_mol
