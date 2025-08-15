from unittest import TestCase

from rdkit.Chem import AllChem


# noinspection SpellCheckingInspection
class TestReactionFromSmarts(TestCase):
    """
    Check how AllChem.ReactionFromSmarts works
    """

    def test_reaction_from_smarts_no_reactants(self):
        reaction_smiles = '>>CCOC'
        try:
            rxn = AllChem.ReactionFromSmarts(reaction_smiles)
            self.assertIsNotNone(rxn)
        except Exception as e:
            self.fail(f"exception raised {e}")

    def test_reaction_from_smarts_no_products(self):
        reaction_smiles = 'CCOC>>'
        try:
            rxn = AllChem.ReactionFromSmarts(reaction_smiles)
            self.assertIsNotNone(rxn)
        except Exception as e:
            self.fail(f"exception raised {e}")

    def test_reaction_from_smarts_no_mols(self):
        reaction_smiles = '>>'
        try:
            rxn = AllChem.ReactionFromSmarts(reaction_smiles)
            self.assertIsNotNone(rxn)
        except Exception as e:
            self.fail(f"exception raised {e}")

    def test_reaction_from_invalid_format(self):
        reaction_smiles = 'CCOC>'
        with self.assertRaises(ValueError):
            AllChem.ReactionFromSmarts(reaction_smiles)
