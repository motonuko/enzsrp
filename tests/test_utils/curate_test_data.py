from collections import UserDict

from enzsrp.data.datasource.original.uniprot.original_uniprot_data_source import OriginalUniprotDataSource
from tests.test_utils.test_default_path import TestDefaultPath


class WriteOnceDict(UserDict):
    def __setitem__(self, key, value):
        if key in self.data:
            pass  # ignore
        super().__setitem__(key, value)


def curate_test_data():
    path = TestDefaultPath().original_uniprot_reviewed_catalytic_activity_json_file
    test_accession_ids = WriteOnceDict()
    for entry in OriginalUniprotDataSource(path).stream_entries_with_catalytic_activity():
        if len(entry.catalytic_activities) == 2:
            if all([len(a.phy_reactions) == 1 and a.phy_reactions[0].has_exp_evidence for a in
                    entry.catalytic_activities]):
                test_accession_ids[
                    'two activity, all activity is single direction with experimental evidence'] = entry.primary_accession
        if len(entry.catalytic_activities) == 1:
            activity = entry.catalytic_activities[0]
            if len(activity.phy_reactions) == 2 and all([r.has_exp_evidence for r in activity.phy_reactions]):
                test_accession_ids[
                    'bidirectional, all physiological reactions have experimental evidences'] = entry.primary_accession
            if len(activity.phy_reactions) == 1 and activity.phy_reactions[0].has_exp_evidence:
                test_accession_ids[
                    'single direction, physiological reaction has experimental evidence'] = entry.primary_accession
            if len(activity.phy_reactions) == 1 and not activity.phy_reactions[
                0].has_exp_evidence and activity.reaction.has_exp_evidence:
                test_accession_ids[
                    'single direction, physiological reaction does not have experimental evidence but reaction does'] = entry.primary_accession
            if len(activity.phy_reactions) == 1 and not activity.phy_reactions[
                0].has_exp_evidence and not activity.reaction.has_exp_evidence:
                test_accession_ids[
                    'single direction, physiological reaction and reaction don\'t have experimental evidence'] = entry.primary_accession
            if len(activity.phy_reactions) == 0:
                test_accession_ids['no direction (no  physiological reaction)'] = entry.primary_accession
        if len(test_accession_ids.keys()) == 6:
            break
    for description, accession_id in test_accession_ids.items():
        print(f"{accession_id}: {description}")

    print('Query')
    print(' OR '.join(test_accession_ids.values()))
    # expected final reaction count = 6 (2 + 2 + 1 + 1 + 0 + 0)
    # 'A0A0E3KBH3' hits


if __name__ == '__main__':
    curate_test_data()
