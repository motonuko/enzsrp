from unittest import TestCase

from rdkit.Chem import AllChem

from enzsrp.data.datasource.original.rhea.original_rhea_datasource import OriginalRheaDatasource
from enzsrp.domain.entity.reaction_direciton import RheaDirectionName
from tests.test_utils.test_default_path import TestDefaultPath


class TestOriginalRheaDatasource(TestCase):

    def setUp(self):
        self.datasource = OriginalRheaDatasource(
            rhea_rxn_dir=TestDefaultPath().original_rhea_rxn_dir,
            rhea2metacyc=TestDefaultPath().original_rhea2metacyc_file,
            rhea_directions=TestDefaultPath().original_rhea_directions_file)

    def test_get_single_reaction(self):
        rxn = self.datasource.get_single_reaction('10001')
        self.assertIsInstance(rxn, AllChem.ChemicalReaction)

    def test_map_master_id_to_direction_id(self):
        lr_id = self.datasource.map_master_id_to_direction_id('10000', RheaDirectionName.RHEA_ID_LR)
        self.assertEqual('10001', lr_id)

        rl_id = self.datasource.map_master_id_to_direction_id('10000', RheaDirectionName.RHEA_ID_RL)
        self.assertEqual('10002', rl_id)

        bi_id = self.datasource.map_master_id_to_direction_id('10000', RheaDirectionName.RHEA_ID_BI)
        self.assertEqual('10003', bi_id)

    def test_map_master_id_to_direction_id_error(self):
        with self.assertRaises(ValueError):
            self.datasource.map_master_id_to_direction_id('10001', RheaDirectionName.RHEA_ID_LR)

    def test_map_rhea_id_to_metacyc_id(self):
        result = self.datasource.map_rhea_id_to_metacyc_id('10000')
        # noinspection SpellCheckingInspection
        self.assertEqual({'PENTANAMIDASE-RXN'}, result)

    def test_map_rhea_id_to_metacyc_id_error(self):
        result = self.datasource.map_rhea_id_to_metacyc_id('10001')
        self.assertEqual(0, len(result))
