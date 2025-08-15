import json
from pathlib import Path
from typing import Optional, List, Dict


class IsoformIdMappingDataSource:

    def __init__(self, path: Path):
        with open(path, 'r') as file:
            self.data = json.load(file)
            assert len([item['from'] for item in self.data["results"]]) == len(set(
                [item['from'] for item in self.data["results"]])), "'from' is not unique"
            self.isoform_id_to_sequence = {item["from"]: item["to"]["sequence"]["value"] for item in
                                           self.data["results"]}
            self.isoform_id_to_active_isoform_id = self._construct_isoform_id_map()

    def map_isoform_to_seq(self, isoform_id: str) -> Optional[str]:
        return self.isoform_id_to_sequence.get(isoform_id, None)

    def bulk_map_isoform_id_to_seq(self, isoform_ids: List[str]) -> Dict:
        mapped = {isoform_id: self.isoform_id_to_sequence.get(isoform_id, None) for isoform_id in isoform_ids}
        mapped = {isoform_id: seq for (isoform_id, seq) in mapped.values() if seq is not None}
        return mapped

    def map_isoform_id_to_active_isoform_id(self, isoform_id: str) -> Optional[str]:
        return self.isoform_id_to_active_isoform_id.get(isoform_id, None)

    def bulk_map_isoform_id_to_active_isoform_id(self, isoform_ids: List[str]) -> Dict[str, str]:
        mapped = {isoform_id: self.isoform_id_to_active_isoform_id.get(isoform_id, None) for isoform_id in isoform_ids}
        return {isoform_id: active_id for isoform_id, active_id in mapped.items() if active_id is not None}

    # NOTE: Multiple IDs are linked to each entry's isoform, so we remap them to single active id.
    # NOTE: If multiple active IDs exist for a single key, the original key is used for consistency.
    def _construct_isoform_id_map(self) -> dict[str, str]:
        uniprot_id_map = {}
        for result in self.data["results"]:
            key_isoform_id = result["from"]
            to = [item["to"] for item in self.data["results"] if item["from"] == key_isoform_id]
            cross_refs = to[0]["uniParcCrossReferences"]
            active_cross_refs = [ref for ref in cross_refs if ref["active"] == True]
            uniprot_swiss = [ref for ref in active_cross_refs if
                             ref["database"] == "UniProtKB/Swiss-Prot protein isoforms"]
            # NOTE: if non-swiss entries exists, it might result in multiple hits.
            if len(uniprot_swiss) > 1:
                uniprot_swiss = [u for u in uniprot_swiss if u['id'] == key_isoform_id]
            uniprot_id_map[key_isoform_id] = uniprot_swiss[0]["id"]
        return uniprot_id_map
