from pathlib import Path
from typing import Generator

import ijson
from tqdm import tqdm

from enzsrp.data.datasource.original.uniprot.uniprot_entity import EntryWithCatalyticActivity


class OriginalUniprotDataSource:

    def __init__(self, uniprot_entries_with_catalytic_activity_json: Path):
        self._uniprot_json_path = uniprot_entries_with_catalytic_activity_json

    def stream_entries_with_catalytic_activity(self) -> Generator[EntryWithCatalyticActivity, None, None]:
        with open(self._uniprot_json_path, 'rb') as file:
            for entry in tqdm(ijson.items(file, 'results.item'),
                              desc='Processing UniProt entries (250K+ entries)', unit=' entries'):
                yield EntryWithCatalyticActivity.from_dict(entry)

    def get_all_isoform_ids(self):
        result = set()
        for entry in self.stream_entries_with_catalytic_activity():
            for activity in entry.catalytic_activities:
                if activity.molecule is not None:
                    if entry.alternative_products_data is not None and len(
                            entry.alternative_products_data.isoforms) >= 1:
                        isoform_ids = [isoform_id for isoform in entry.alternative_products_data.isoforms for isoform_id
                                       in isoform.isoformIds]
                        result.update(isoform_ids)
        return result
