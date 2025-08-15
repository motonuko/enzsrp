import json
import os
import warnings
from collections import defaultdict
from dataclasses import dataclass, asdict
from typing import List, Optional

from rxndata2.data.datasource.original.mcsa.mcsa_entity import ResidueSequence, McsaEntry
from rxndata2.data.datasource.original.uniprot.uniprot_entity import EntryWithCatalyticActivity
from rxndata2.domain.entity.amino_acid import three_to_one


# Intended for use by other classes
@dataclass(frozen=True)
class MCSAResidueSequence:
    mcsa_id: str
    uniprot_id: Optional[str]  # NOTE: Some residue is not linked to UniProt
    code: str  # e.g. 'Asp'
    is_reference: bool
    resid: int  # e.g. 7

    @classmethod
    def from_dict(cls, residue_sequence: ResidueSequence, mcsa_id: str):
        return cls(mcsa_id=mcsa_id, **asdict(residue_sequence))


class OriginalMCSADataSource:

    def __init__(self, mcsa_data_dir):
        self.mcsa_data_dir = mcsa_data_dir
        self.uniprot_id_residues_map = self._load_uniprot_id_residues_map()

    def _load_json_files(self):
        json_data = []
        for root, dirs, files in os.walk(self.mcsa_data_dir):
            for file in files:
                if file.endswith('.json'):
                    filepath = os.path.join(root, file)
                    with open(filepath, 'r', encoding='utf-8') as f:
                        data = json.load(f)
                        json_data.append(data)
        entries = [j for jo in json_data for j in jo["results"]]
        entries = sorted(entries, key=lambda x: x["mcsa_id"])
        return entries

    # def load_entries(self):
    #     entries = self._load_json_files()
    #     entries = [McsaEntry.from_dict(entry) for entry in entries]
    #     return entries

    def _load_uniprot_id_residues_map(self) -> dict[str, List[MCSAResidueSequence]]:
        entries = self._load_json_files()
        entries = [McsaEntry.from_dict(entry) for entry in entries]
        result = defaultdict(list)
        for entry in entries:
            for residue in entry.residues:
                # NOTE: entry.reference_uniprot_id may contain multiple UniProt IDs (?),
                # so the residue's UniProt ID is used instead
                if residue.residue_sequences.uniprot_id is not None:
                    result[residue.residue_sequences.uniprot_id].append(
                        MCSAResidueSequence.from_dict(residue.residue_sequences, entry.mcsa_id))
        return result

    def get_residues_by_entry(self, entry: EntryWithCatalyticActivity):
        result = set()
        accessions = [entry.primary_accession, *entry.secondary_accessions]
        if entry.last_sequence_update_date.year >= 2018:
            warnings.warn("MCS-A mapping might be legacy")
        for accession in accessions:
            residues = self.uniprot_id_residues_map.get(accession, [])
            try:
                for residue in residues:
                    if len(entry.sequence) < residue.resid or entry.sequence[residue.resid - 1] != three_to_one(
                            residue.code):
                        raise ValueError(
                            f"residue does not match to Uniprot Entry Sequence, MCS-A ID: {residue.mcsa_id}")
                result.update(residues)
            except ValueError as e:
                warnings.warn(str(e))
                return set()  # If even one residue doesn't match, discard all. Data might be old.
        return result
