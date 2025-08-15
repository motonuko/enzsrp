from pathlib import Path
from typing import Optional

import pandas as pd

from rxndata2.domain.entity.reaction_direciton import MetaCycDirection
from rxndata2.utils.list_utils import get_first_from_single_list

UNIQUE_ID_TAG = "UNIQUE-ID"
RHEA_TAG = "DBLINKS - (RHEA"
REACTION_DIRECTION_TAG = "REACTION-DIRECTION"


def _extract_rhea_id_from_line(text: str) -> str:
    start_marker = 'DBLINKS - (RHEA "'
    end_marker = '"'

    start_index = text.find(start_marker)
    end_index = text.find(end_marker, start_index + len(start_marker))

    if start_index != -1 and end_index != -1:
        extracted_text = text[start_index + len(start_marker):end_index]
        return str(extracted_text)
    else:
        raise RuntimeError()


def _extract_value_from_line(text: str) -> str:
    return text.replace("\n", "").strip().split(" - ")[1]


def _extract_reaction_direction_from_line(text: str) -> MetaCycDirection:
    direction_text = _extract_value_from_line(text)
    return MetaCycDirection.from_text(direction_text)


def _extract_rxn_direction_from_one_rxn(lines: list[str]):
    rhea_lines = [hog for hog in lines if RHEA_TAG in hog]
    reaction_direction_lines = [hog for hog in lines if REACTION_DIRECTION_TAG in hog]
    unique_id_lines = [hog for hog in lines if UNIQUE_ID_TAG in hog]
    if len(reaction_direction_lines) < 1 or len(unique_id_lines) < 1:
        return None, None, None, None
    unique_id = get_first_from_single_list(unique_id_lines)
    unique_id = _extract_value_from_line(unique_id)
    rhea_ids = [_extract_rhea_id_from_line(line) for line in rhea_lines]
    reaction_direction = get_first_from_single_list(reaction_direction_lines)
    reaction_direction = _extract_reaction_direction_from_line(reaction_direction)
    if len(rhea_ids) > 2:
        raise RuntimeError("unexpected. should be <= 2")
    rhea_ids += [None] * (2 - len(rhea_ids))
    return unique_id, rhea_ids[0], rhea_ids[1], reaction_direction


def _parse_single_rxn_format_file(file_path: Path) -> list[dict]:
    with open(file_path, 'r', encoding='unicode_escape') as file:
        lines = file.readlines()

    loading_mol_lines = []
    result = []
    for line in lines:
        if line.startswith("UNIQUE-ID - "):
            loading_mol_lines = []

        loading_mol_lines.append(line)

        if line.startswith("//"):
            unique_id, rhea, rhea2, direction = _extract_rxn_direction_from_one_rxn(loading_mol_lines)
            if direction is None:
                continue
            result.append(
                {"unique_id_rxn": unique_id, "rhea_id1": rhea, "rhea_id2": rhea2, "rxn_direction": direction.value})
    return result


class MetaCycDirectionMapper:

    def __init__(self, reaction_dat: Path):
        parsed = _parse_single_rxn_format_file(reaction_dat)
        self.df_metacyc_id_to_direction = pd.DataFrame(parsed)
        pass

    # If direction is not defined, set to None
    def metacyc_id_to_direction(self, metacyc_id: str) -> Optional[MetaCycDirection]:
        result = self.df_metacyc_id_to_direction[self.df_metacyc_id_to_direction["unique_id_rxn"] == metacyc_id]
        if len(result) == 1:
            direction_str = result["rxn_direction"].iloc[0]
            return MetaCycDirection.from_text(direction_str)
        elif len(result) >= 1:
            raise ValueError("Unexpected. MetaCyc ID should be unique(?)")
        else:
            return None
