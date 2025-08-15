from pathlib import Path

import pandas as pd

from tests.test_utils.test_default_path import TestDefaultPath


def check_output_file_content(file1: Path, file2: Path):
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)

    ignore_cols = []
    cols_to_compare = [c for c in df1.columns if c not in ignore_cols]

    # Expected to be False because 'rxn_evidence' and 'phy_rxn_evidence' have randomness.
    assert not df1[cols_to_compare].equals(df2[cols_to_compare])

    # Expected to be True if 'rxn_evidence' and 'phy_rxn_evidence' are removed.
    ignore_cols = ["rxn_evidence", "phy_rxn_evidence"]
    cols_to_compare = [c for c in df1.columns if c not in ignore_cols]
    assert df1[cols_to_compare].equals(df2[cols_to_compare])


if __name__ == '__main__':
    # NOTE: Change the paths to where you placed original file and generated file.
    path1 = TestDefaultPath().output_dir / 'enzsrp_full.csv'
    path2 = TestDefaultPath().output_dir / 'enzsrp_full_final.csv'
    check_output_file_content(path1, path2)
