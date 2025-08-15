import re
from typing import List, Tuple

from rdkit import Chem
from rdkit.Chem import AllChem

TITLE_KEY = 'title'


# noinspection PyArgumentList
def get_all_mol_names_in_rxn(rx: AllChem.ChemicalReaction):
    assert len(rx.GetAgents()) == 0
    return set(
        [mol.GetProp(TITLE_KEY) for mol in rx.GetReactants()] + [mol.GetProp(TITLE_KEY) for mol in rx.GetProducts()])


def _parse_rxn_block_manually(rxn_block: str, title_key: str) -> Tuple[List[Chem.Mol], List[Chem.Mol]]:
    """
    Parse each mol block individually to extract the title from the mol block.
    """
    split_by_mol_symbol = rxn_block.split("$MOL")
    header_block = split_by_mol_symbol[0]
    mol_blocks = split_by_mol_symbol[1:]

    n_reactant_product = [int(n) for n in header_block.strip().split('\n')[-1].strip().split()]
    if len(n_reactant_product) != 2:
        raise ValueError("Unexpected format")
    n_reactants, n_products = n_reactant_product

    reactants = []
    products = []
    for block in mol_blocks:
        block = block.lstrip('\n')
        title = block.strip().split('\n')[0].strip()
        mol_text = block.split("M  END")[0] + "M  END\n"
        mol = Chem.MolFromMolBlock(mol_text, removeHs=False, sanitize=True)
        if mol is None:
            raise ValueError("Failed to parse mol block")
        mol.SetProp(title_key, title)
        if len(reactants) < n_reactants:
            reactants.append(mol)
        else:
            products.append(mol)
    return reactants, products


def reaction_from_rxn_block_with_cleaning(rxn_block: str) -> AllChem.ChemicalReaction:
    try:
        return AllChem.ReactionFromRxnBlock(rxn_block, sanitize=True)
    except ValueError:
        before = r"M  END(\n+)(\s+)\$MOL"
        after = r"M  END\n$MOL"
        result = re.sub(before, after, rxn_block, flags=re.MULTILINE)
        return AllChem.ReactionFromRxnBlock(result, sanitize=True)


# NOTE: Considered using HasSubstructMatch, but it did not yield sufficient matches.
# Attach titles parsed manually from the reaction block to each Mol object parsed by RDKit as properties using SetProp.
# noinspection PyArgumentList
def reaction_from_rxn_block_with_mol_title(rxn_block: str):
    reaction = reaction_from_rxn_block_with_cleaning(rxn_block)
    reactants, products = _parse_rxn_block_manually(rxn_block, TITLE_KEY)
    for my_reactant, rd_reactant in zip(reactants, reaction.GetReactants()):
        assert Chem.MolToSmiles(my_reactant, True) == Chem.MolToSmiles(rd_reactant), "unexpected"
        # Chem.AssignStereochemistry(my_reactant, cleanIt=True, force=True)
        # my_reactant = Chem.RemoveHs(my_reactant)
        rd_reactant.SetProp(TITLE_KEY, my_reactant.GetProp(TITLE_KEY))
    for my_product, rd_product in zip(products, reaction.GetProducts()):
        assert Chem.MolToSmiles(my_product, True) == Chem.MolToSmiles(rd_product), "unexpected"
        # Chem.AssignStereochemistry(my_product, cleanIt=True, force=True)
        # my_product = Chem.RemoveHs(my_product)
        rd_product.SetProp(TITLE_KEY, my_product.GetProp(TITLE_KEY))
    return reaction
