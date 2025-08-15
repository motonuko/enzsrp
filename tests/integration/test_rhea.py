import os
from unittest import TestCase

from tqdm import tqdm

from rxndata2.data.datasource.original.rhea.original_rhea_datasource import OriginalRheaDatasource
from tests.test_utils.test_default_path import TestDefaultPath


class TestOriginalRheaDatasource(TestCase):

    def setUp(self):
        self.datasource = OriginalRheaDatasource(
            rhea_rxn_dir=TestDefaultPath().original_rhea_rxn_dir,
            rhea2metacyc=TestDefaultPath().original_rhea2metacyc_file,
            rhea_directions=TestDefaultPath().original_rhea_directions_file)

    def test_get_single_reaction_for_all_rxn_files(self):
        file_list = []
        for root, dirs, files in os.walk(TestDefaultPath().original_rhea_rxn_dir):
            for file in files:
                file_list.append(file)
        file_list = [name.split('.')[0] for name in file_list]
        for na in tqdm(file_list):
            self.datasource.get_single_reaction(na)
        print(self.datasource.warned_rhea_ids)
