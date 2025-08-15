import os
from pathlib import Path

import pandas as pd
from dotenv import load_dotenv

from rxndata2.utils import env_var_names

load_dotenv()


def check_data_identity():
    path = Path(os.getenv(env_var_names.output_dir))

    # Final dataset for the paper
    enzsrp_final = pd.read_csv(path / 'enzsrp_final.csv')
    enzsrp_full_final = pd.read_csv(path / 'enzsrp_full_final.csv')

    # After fixing '' column randomness
    enzsrp = pd.read_csv(path / 'enzsrp.csv')
    enzsrp_full = pd.read_csv(path / 'enzsrp_full.csv')

    exclude_cols = ['rxn_evidence', 'phy_rxn_evidence']
    # exclude_cols = []
    compare_cols = [col for col in enzsrp_final.columns if col not in exclude_cols]

    assert enzsrp[compare_cols].equals(enzsrp_final[compare_cols])
    assert enzsrp_full[compare_cols].equals(enzsrp_full_final[compare_cols])

    print('Except for the excluded columns, the contents of the DataFrames are identical.')


if __name__ == '__main__':
    check_data_identity()
