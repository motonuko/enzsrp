import unittest

from rxndata2.domain.entity.evidence_and_clusion_ontology import ECO


class TestUniqueEnumValues(unittest.TestCase):
    def test_unique_enum_values(self):
        values = [enum.value for enum in ECO]
        unique_values = set(values)
        self.assertEqual(len(values), len(unique_values), "Enum values are not unique")


if __name__ == '__main__':
    unittest.main()
