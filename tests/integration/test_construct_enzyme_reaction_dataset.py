import tempfile
from pathlib import Path
from unittest import TestCase

from enzsrp.presentation.build_enzyme_reaction_dataset import EnzymeReactionDatasetBuilder
from tests.test_utils.test_default_path import TestDefaultPath


class TestEnzymeReactionDatasetConstructor(TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_dir_path = Path(self.temp_dir.name)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_construct_full_directed_dataset_with_metacyc(self):
        assert TestDefaultPath().original_metacyc_reactions_dat is not None
        dataset_constructor = EnzymeReactionDatasetBuilder(
            uniprot_entries_json=TestDefaultPath().test_data.joinpath(
                'uniprotkb_A0A0E3KBH3_OR_A0A072ULZ1_OR_A_2024_09_24.json'),
            uniprot_isoform_uniparc_mapping_json=TestDefaultPath().test_data.joinpath('output',
                                                                                      "idmapping_2024_11_13_isoform_uniparc.json"),
            rhea_rxn_dir=TestDefaultPath().original_rhea_rxn_dir,
            rhea_2_metacyc_file=TestDefaultPath().original_rhea2metacyc_file,
            rhea_directions_file=TestDefaultPath().original_rhea_directions_file,
            mcsa_data_dir=None,
            metacyc_reactions_dat=TestDefaultPath().original_metacyc_reactions_dat,
            output_path=self.temp_dir_path / 'sample.csv',
            use_undefined_direction_rxn=False,
            allow_non_exp_evidence=False
        )
        df = dataset_constructor.construct_full_directed_dataset()
        self.assertEqual(7, len(df))

    def test_construct_full_directed_dataset_no_metacyc(self):
        dataset_constructor = EnzymeReactionDatasetBuilder(
            uniprot_entries_json=TestDefaultPath().test_data.joinpath(
                'uniprotkb_A0A0E3KBH3_OR_A0A072ULZ1_OR_A_2024_09_24.json'),
            uniprot_isoform_uniparc_mapping_json=TestDefaultPath().test_data.joinpath('output',
                                                                                      "idmapping_2024_09_21_isoform_uniparc.json"),
            rhea_rxn_dir=TestDefaultPath().original_rhea_rxn_dir,
            rhea_2_metacyc_file=TestDefaultPath().original_rhea2metacyc_file,
            rhea_directions_file=TestDefaultPath().original_rhea_directions_file,
            mcsa_data_dir=None,
            metacyc_reactions_dat=None,
            output_path=self.temp_dir_path / 'sample.csv',
            use_undefined_direction_rxn=False,
            allow_non_exp_evidence=False
        )
        df = dataset_constructor.construct_full_directed_dataset()
        self.assertEqual(6, len(df))
