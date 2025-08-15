from unittest import TestCase

from enzsrp.data.datasource.original.uniprot.original_uniprot_data_source import OriginalUniprotDataSource
from tests.test_utils.test_default_path import TestDefaultPath


class TestOriginalUniprotDataSource(TestCase):
    def test_stream_entries_with_catalytic_activity(self):
        file_path = TestDefaultPath().test_data.joinpath("uniprotkb_accession_A1L3X0_OR_accession_2024_07_27.json")
        source = OriginalUniprotDataSource(file_path)
        loaded = []
        for entry in source.stream_entries_with_catalytic_activity():
            loaded.append(entry)
        self.assertEqual(len(loaded), 2)
        self.assertEqual(all(element is not None for element in loaded), True)

    def test_get_all_isoform_ids(self):
        file_path = TestDefaultPath().test_data.joinpath("uniprotkb_accession_F1MAB7_OR_O14975_2024_09_21.json")
        source = OriginalUniprotDataSource(file_path)
        result = source.get_all_isoform_ids()
        self.assertEqual(result, {'F1MAB7-2', 'F1MAB7-1', 'O14975-1', 'O14975-2', 'F1MAB7-3'})
