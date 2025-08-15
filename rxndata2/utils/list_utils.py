from typing import List


def get_first_from_single_list(target_list: List, additional_message=""):
    if len(target_list) == 1:
        return target_list[0]
    else:
        raise ValueError(f"Invalid list length: len=={len(target_list)}, {additional_message}")


def get_first_from_single_set(target_set: set, additional_message=""):
    if len(target_set) == 1:
        return target_set.pop()
    else:
        raise ValueError(f"Invalid list length: len=={len(target_set)}, {additional_message}")
