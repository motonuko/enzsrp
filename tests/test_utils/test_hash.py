import tempfile
from pathlib import Path
from unittest import TestCase

from tests.test_utils.hash import calculate_file_hash


class TestEnzymeReactionDatasetConstructor(TestCase):
    def setUp(self):
        self.temp_dir = tempfile.TemporaryDirectory()
        self.temp_dir_path = Path(self.temp_dir.name)

    def tearDown(self):
        self.temp_dir.cleanup()

    def test_hash(self):
        path = self.temp_dir_path.joinpath("hash_test.txt")
        content = "text text text text"
        with open(path, 'w') as f:
            f.write(content)
        first_hash = calculate_file_hash(path)
        with open(path, 'w') as f:
            f.write(content)
        second_hash = calculate_file_hash(path)
        self.assertEqual(first_hash, second_hash)
