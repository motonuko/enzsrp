import json
import unittest
from datetime import datetime

from enzsrp.data.datasource.original.uniprot.uniprot_entity import EntryWithCatalyticActivity, LocationModifier, \
    Isoform, RheaReference, RheaID
from enzsrp.domain.entity.evidence_and_clusion_ontology import ECO
from enzsrp.domain.entity.reaction_direciton import ReactionDirection
from tests.test_utils.test_default_path import TestDefaultPath


class MyTestCase(unittest.TestCase):
    def test_1(self):
        file_path = TestDefaultPath().test_data.joinpath('single_entry', 'O00115.json')
        with open(file_path, 'r') as file:
            data = json.load(file)
        entry = EntryWithCatalyticActivity.from_dict(data)

        # top level properties
        self.assertEqual('O00115', entry.primary_accession)
        self.assertEqual(['B2RD06', 'B7Z4K6', 'O43910'], entry.secondary_accessions)
        self.assertEqual('DNS2A_HUMAN', entry.uniprot_kb_id)
        self.assertEqual('1: Evidence at protein level', entry.protein_existence)
        self.assertEqual(1, len(entry.catalytic_activities))
        self.assertEqual(0, len(entry.binding_sites))
        self.assertEqual(1, len(entry.active_sites))
        self.assertNotEqual(0, len(entry.features))
        self.assertEqual(360, len(entry.sequence))
        self.assertEqual('MIPLLLAALL', entry.sequence[:10])
        self.assertEqual('ARKPSRAYKI', entry.sequence[-10:])
        self.assertEqual(datetime.strptime('1998-07-15', "%Y-%m-%d"), entry.last_sequence_update_date)
        self.assertEqual(0, len(entry.cofactor_data_list))
        self.assertNotEqual(None, entry.alternative_products_data)

        activity = entry.catalytic_activities[0]
        self.assertEqual(None, activity.molecule)
        self.assertEqual(0, len(activity.phy_reactions))

        reaction = activity.reaction
        self.assertEqual('3.1.22.1', reaction.ec_number)
        self.assertEqual(2, len(reaction.evidences))
        self.assertEqual(0, len(reaction.cross_references))
        self.assertEqual(
            'Endonucleolytic cleavage to nucleoside 3\'-phosphates and 3\'-phosphooligonucleotide end-products.',
            reaction.name)
        self.assertEqual(False, reaction.has_rhea_rxn_reference)
        self.assertEqual(None, reaction.rhea_main_reference)
        self.assertEqual(True, reaction.has_exp_evidence)

        self.assertEqual(ECO.EXPERIMENTAL, reaction.evidences[0].type)
        self.assertEqual(ECO.EXPERIMENTAL, reaction.evidences[1].type)

        active_site = entry.active_sites[0]
        self.assertEqual("Active site", active_site.type)
        self.assertEqual(0, len(active_site.featureCrossReference))
        self.assertEqual(1, len(active_site.evidences))
        self.assertEqual("", active_site.description)

        self.assertEqual(295, active_site.location.start_value)
        self.assertEqual(LocationModifier.EXACT, active_site.location.start_modifier)
        self.assertEqual(295, active_site.location.end_value)
        self.assertEqual(LocationModifier.EXACT, active_site.location.end_modifier)
        self.assertEqual('295', active_site.location.string_notation)

        self.assertEqual(ECO.CURATOR_INFERENCE_EVIDENCE, active_site.evidences[0].type)

        alternative_products = entry.alternative_products_data
        self.assertEqual(2, len(alternative_products.isoforms))
        self.assertEqual(1, len(alternative_products.search_isoform_by_isoform_name('1')))
        self.assertEqual(1, len(alternative_products.search_isoform_by_isoform_name('2')))

        isoform1: Isoform = alternative_products.isoforms[0]
        self.assertEqual('1', isoform1.name)
        self.assertEqual(['O00115-1'], isoform1.isoformIds)
        self.assertEqual('Displayed', isoform1.isoformSequenceStatus)
        self.assertNotEqual(None, isoform1.sequenceIds)

        isoform2: Isoform = alternative_products.isoforms[1]
        self.assertEqual('Described', isoform2.isoformSequenceStatus)

    def test_2(self):
        file_path = TestDefaultPath().test_data.joinpath('single_entry', 'A4Q9F3.json')
        with open(file_path, 'r') as file:
            data = json.load(file)

        entry = EntryWithCatalyticActivity.from_dict(data)
        self.assertEqual(1, len(entry.catalytic_activities))

        activity = entry.catalytic_activities[0]
        self.assertEqual(1, len(activity.phy_reactions))
        self.assertEqual('Protein polyglycylase TTLL10', activity.molecule)

        reaction = activity.reaction
        self.assertEqual(True, reaction.has_rhea_rxn_reference)
        self.assertNotEqual(None, reaction.rhea_main_reference)

        main_ref = reaction.rhea_main_reference
        self.assertEqual(RheaReference(database='Rhea', db_id=RheaID('RHEA:67184')), main_ref)
        self.assertEqual('67184', main_ref.db_id.id_intstr)

        phy_reaction = activity.phy_reactions[0]
        self.assertEqual(ReactionDirection.LEFT_TO_RIGHT, phy_reaction.direction)
        self.assertEqual(3, len(phy_reaction.evidences))
        self.assertEqual(ECO.CURATOR_INFERENCE_EVIDENCE, phy_reaction.evidences[0].type)
        self.assertEqual(False, phy_reaction.has_exp_evidence)
        self.assertEqual(0, len(phy_reaction.cross_references))

        self.assertEqual(entry.sequence,
                         entry.get_new_ptm_sequence(activity))  # ptm seq is same as primary seq in this entry

    def test_3(self):
        file_path = TestDefaultPath().test_data.joinpath('single_entry', 'P0DKX7.json')
        with open(file_path, 'r') as file:
            data = json.load(file)

        entry = EntryWithCatalyticActivity.from_dict(data)
        self.assertEqual(1, len(entry.catalytic_activities))
        activity = entry.catalytic_activities[0]
        ptm_seq = entry.get_new_ptm_sequence(activity)
        expected = ("MQQSHQAGYANAADRESGIPAAVLDGIKAVAKEKNATLMFRLVNPHSTSLIAEGVATKGLGVHAKSSDWGLQAGYIPVNPNLSKLFGRAPE"
                    "VIARADNDVNSSLAHGHTAVDLTLSKERLDYLRQAGLVTGMADGVVASNHAGYEQFEFRVKETSDGRYAVQYRRKGGDDFEAVKVIGNAAG"
                    "IPLTADIDMFAIMPHLSNFRDSARSSVTSGDSVTDYLARTRRAASEATGGLDRERIDLLWKIARAGARSAVGTEARRQFRYDGDMNIGVIT"
                    "DFELEVRNALNRRAHAVGAQDVVQHGTEQNNPFPEADEK")
        self.assertEqual(expected, ptm_seq)
