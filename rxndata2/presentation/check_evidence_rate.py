import os
from pathlib import Path

import pandas as pd
from dotenv import load_dotenv
from matplotlib import pyplot as plt

load_dotenv()

from rxndata2.utils import env_var_names


def check_evidence_ratio(data_path: Path):
    df = pd.read_csv(data_path)
    column = 'rxn_evidence'
    # 'phy_rxn_evidence'
    value_counts_sorted = df['rxn_evidence'].value_counts().sort_values(ascending=False)

    plt.figure(figsize=(8, 8))
    plt.pie(value_counts_sorted, labels=value_counts_sorted.index, autopct='%1.1f%%', startangle=90, counterclock=False,
            wedgeprops={'edgecolor': 'black'})
    plt.title('Count of rxn_evidence Types ')
    plt.show()


if __name__ == '__main__':
    output_dir = Path(os.getenv(env_var_names.output_dir))
    check_evidence_ratio(output_dir / 'enzsrp.csv')
    check_evidence_ratio(output_dir / 'enzsrp_full.csv')
