from typing import Optional

from rdkit import Chem
from rdkit.Chem import AllChem

from enzsrp.extension.mol_ext import remove_atom_mapping


def remove_atom_mapping_and_h_from_reaction_smiles(reaction_smiles: str) -> Optional[str]:
    try:
        reactants_smiles, _, products_smiles = reaction_smiles.split('>')
        rxn = AllChem.ReactionFromSmarts(f"{reactants_smiles}>>{products_smiles}")
        reactants = [remove_atom_mapping(mol) for mol in rxn.GetReactants()]
        products = [remove_atom_mapping(mol) for mol in rxn.GetProducts()]
        reactants_smiles = '.'.join([Chem.MolToSmiles(mol) for mol in reactants])
        products_smiles = '.'.join([Chem.MolToSmiles(mol) for mol in products])
        return f"{reactants_smiles}>>{products_smiles}"
    except Exception as e:
        print(f"Error processing reaction SMILES: {e}")
        return None
