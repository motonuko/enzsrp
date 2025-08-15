from pathlib import Path
from unittest import TestCase
from unittest.mock import patch

from enzsrp.data.datasource.original.uniprot.uniprot_entity import RheaReference, RheaID
from enzsrp.domain.entity.reaction_direciton import MetaCycDirection, ReactionDirection
from enzsrp.presentation.build_enzyme_reaction_dataset import EnzymeReactionDatasetBuilder
from tests.test_utils.test_default_path import TestDefaultPath


class TestEnzymeReactionDatasetConstructor(TestCase):

    def setUp(self):
        assert TestDefaultPath().original_metacyc_reactions_dat is not None
        self.builder = EnzymeReactionDatasetBuilder(
            uniprot_entries_json=TestDefaultPath().test_data.joinpath(
                'uniprotkb_accession_A0A059TC02_OR_A0A05_2024_09_21.json'),
            uniprot_isoform_uniparc_mapping_json=TestDefaultPath().output_dir.joinpath(
                "idmapping_2024_11_13_isoform_uniparc.json"),
            rhea_rxn_dir=TestDefaultPath().original_rhea_rxn_dir,
            rhea_2_metacyc_file=TestDefaultPath().original_rhea2metacyc_file,
            rhea_directions_file=TestDefaultPath().original_rhea_directions_file,
            mcsa_data_dir=None,
            metacyc_reactions_dat=TestDefaultPath().original_metacyc_reactions_dat,
            output_path=Path("./exp_activities.csv"),
            use_undefined_direction_rxn=False,
            allow_non_exp_evidence=False,
        )

        self.rhea_ref = RheaReference(database='Rhea', db_id=RheaID('RHEA:10001'))

    @patch(
        'enzsrp.data.datasource.original.rhea.original_rhea_datasource.OriginalRheaDatasource.map_rhea_id_to_metacyc_id')
    @patch(
        'enzsrp.data.datasource.original.metacyc.parse_reactions_dat.MetaCycDirectionMapper.metacyc_id_to_direction')
    def test_get_metacyc_directions_1(self, metacyc_rhea_id_to_cyc_id_func, rhea_id_to_metacyc_func):
        rhea_id_to_metacyc_func.return_value = ['cyc_id_1', 'cyc_id_2']
        metacyc_rhea_id_to_cyc_id_func.return_value = MetaCycDirection.RIGHT_TO_LEFT
        result = self.builder.get_metacyc_directions(self.rhea_ref)
        self.assertEqual([ReactionDirection.RIGHT_TO_LEFT], result)

    @patch(
        'enzsrp.data.datasource.original.rhea.original_rhea_datasource.OriginalRheaDatasource.map_rhea_id_to_metacyc_id')
    @patch(
        'enzsrp.data.datasource.original.metacyc.parse_reactions_dat.MetaCycDirectionMapper.metacyc_id_to_direction')
    def test_get_metacyc_directions_2(self, metacyc_rhea_id_to_cyc_id_func, rhea_id_to_metacyc_func):
        rhea_id_to_metacyc_func.return_value = ['cyc_id_1', 'cyc_id_2']

        def return_direction_mock(metacyc_id):
            if metacyc_id == 'cyc_id_1':
                return MetaCycDirection.RIGHT_TO_LEFT
            return MetaCycDirection.LEFT_TO_RIGHT

        metacyc_rhea_id_to_cyc_id_func.side_effect = return_direction_mock
        result = self.builder.get_metacyc_directions(self.rhea_ref)
        self.assertEqual([ReactionDirection.LEFT_TO_RIGHT, ReactionDirection.RIGHT_TO_LEFT], result)

    @patch(
        'enzsrp.data.datasource.original.rhea.original_rhea_datasource.OriginalRheaDatasource.map_rhea_id_to_metacyc_id')
    @patch(
        'enzsrp.data.datasource.original.metacyc.parse_reactions_dat.MetaCycDirectionMapper.metacyc_id_to_direction')
    def test_get_metacyc_directions_3(self, metacyc_rhea_id_to_cyc_id_func, rhea_id_to_metacyc_func):
        rhea_id_to_metacyc_func.return_value = ['cyc_id_1', 'cyc_id_2']

        def return_direction_mock(metacyc_id):
            if metacyc_id == 'cyc_id_1':
                return MetaCycDirection.LEFT_TO_RIGHT
            return MetaCycDirection.RIGHT_TO_LEFT

        metacyc_rhea_id_to_cyc_id_func.side_effect = return_direction_mock
        result = self.builder.get_metacyc_directions(self.rhea_ref)
        self.assertEqual([ReactionDirection.LEFT_TO_RIGHT, ReactionDirection.RIGHT_TO_LEFT], result)

    @patch(
        'enzsrp.data.datasource.original.rhea.original_rhea_datasource.OriginalRheaDatasource.map_rhea_id_to_metacyc_id')
    @patch(
        'enzsrp.data.datasource.original.metacyc.parse_reactions_dat.MetaCycDirectionMapper.metacyc_id_to_direction')
    def test_get_metacyc_directions_4(self, metacyc_rhea_id_to_cyc_id_func, rhea_id_to_metacyc_func):
        rhea_id_to_metacyc_func.return_value = ['cyc_id_1', 'cyc_id_2']
        metacyc_rhea_id_to_cyc_id_func.returns = None
        result = self.builder.get_metacyc_directions(self.rhea_ref)
        self.assertEqual([], result)
