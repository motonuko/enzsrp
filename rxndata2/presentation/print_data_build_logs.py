import os
from pathlib import Path

import pandas as pd
from dotenv import load_dotenv

from rxndata2.utils import env_var_names

load_dotenv()


def main():
    output_dir = Path(os.getenv(env_var_names.output_dir))
    skipped_activities = output_dir / 'skipped_activities_enzsrp_full.csv'

    df_skipped_activities = pd.read_csv(skipped_activities)
    print(df_skipped_activities['reason'].value_counts())
    print(f"total: {len(df_skipped_activities)}")
    print()

    skipped_directions = output_dir / 'skipped_directions_enzsrp_full.csv'
    df_skipped_directions = pd.read_csv(skipped_directions)
    print(df_skipped_directions['reason'].value_counts())
    print(f"unique RheaIDs failing to map: {len(df_skipped_directions['rhea_id'].unique())}")
    print(f"total: {len(df_skipped_directions)}")


if __name__ == '__main__':
    main()
