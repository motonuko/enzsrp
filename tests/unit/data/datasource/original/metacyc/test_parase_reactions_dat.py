import unittest

# noinspection PyProtectedMember
from enzsrp.data.datasource.original.metacyc.parse_reactions_dat import MetaCycDirectionMapper, \
    _extract_rhea_id_from_line
from enzsrp.domain.entity.reaction_direciton import MetaCycDirection
from tests.test_utils.test_default_path import TestDefaultPath


class MyTestCase(unittest.TestCase):
    def test_mapper_success(self):
        self.metacyc_mapper = MetaCycDirectionMapper(TestDefaultPath().original_metacyc_reactions_dat)

        result = self.metacyc_mapper.metacyc_id_to_direction('RXN-22349')
        self.assertEqual(MetaCycDirection.REVERSIBLE, result)

        # noinspection SpellCheckingInspection
        result = self.metacyc_mapper.metacyc_id_to_direction('PENTANAMIDASE-RXN')
        self.assertEqual(MetaCycDirection.PHYSIOL_LEFT_TO_RIGHT, result)

        # noinspection SpellCheckingInspection
        result = self.metacyc_mapper.metacyc_id_to_direction('DIMETHYLALLYLCISTRANSFERASE-RXN')
        self.assertEqual(MetaCycDirection.LEFT_TO_RIGHT, result)

        result = self.metacyc_mapper.metacyc_id_to_direction('RXN-8790')
        self.assertEqual(MetaCycDirection.RIGHT_TO_LEFT, result)

        result = self.metacyc_mapper.metacyc_id_to_direction("SAMPLE_ID_XXXXXXXXXXXXXXXXXXXXXXXXXXXX")
        self.assertEqual(None, result)

    def test_extract_rhea_id_from_line(self):
        sample_line = 'DBLINKS - (RHEA "99999" NIL |xxxxxxx| 9999999999 NIL NIL)'
        output = _extract_rhea_id_from_line(sample_line)
        self.assertEqual(output, '99999')


if __name__ == '__main__':
    unittest.main()
