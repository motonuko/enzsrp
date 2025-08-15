import unittest

from rdkit.Chem import AllChem

from rxndata2.extension.rdkit_warning import RDKitWarningInterceptor, AmbiguousStereoChemistryWarningException
from tests.test_utils.test_default_path import TestDefaultPath


class TestRDKitWarningInterceptor(unittest.TestCase):
    def setUp(self):
        self.interceptor = RDKitWarningInterceptor()

    def test_ambiguous_stereochemistry_warning(self):
        with self.assertRaises(AmbiguousStereoChemistryWarningException):
            with self.interceptor:
                AllChem.ReactionFromRxnFile(TestDefaultPath().test_data.joinpath('71500.rxn').as_posix())

    def test_ambiguous_stereochemistry_warning2(self):
        interceptor2 = RDKitWarningInterceptor()
        with self.assertRaises(RuntimeError):
            with self.interceptor:
                with interceptor2:
                    pass


if __name__ == "__main__":
    unittest.main()
