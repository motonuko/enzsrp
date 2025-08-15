import os
from pathlib import Path

import click
import pandas as pd
import tqdm
from PIL import ImageDraw
from rdkit.Chem import AllChem, Draw

from rxndata2.utils import env_var_names


@click.command()
@click.option('--enz-srp-file', type=click.Path(exists=True),
              default=Path(os.getenv(env_var_names.output_dir)) / 'enzsrp.csv'
              if os.getenv(env_var_names.output_dir) is not None else None,
              help='Specify the file path for the protein data downloaded from UniProt')
@click.option('--output-file', type=click.Path(),
              default=Path(os.getenv(env_var_names.output_dir)) / 'reactions.pdf'
              if os.getenv(env_var_names.output_dir) is not None else None,
              help='Specify the file path for the downloaded isoform mapping file')
def draw_reactions(enz_srp_file: Path, output_file: Path):
    dfx = pd.read_csv(enz_srp_file, dtype={"rhea_id": str})
    rxns = set(dfx["rhea_id"] + "|" + dfx["rxn"])
    rxns = list(rxns)
    images = []
    for rhea_rxn in tqdm.tqdm(rxns):
        rxn = rhea_rxn.split("|")[1]
        reaction = AllChem.ReactionFromSmarts(rxn)
        img = Draw.ReactionToImage(reaction)
        draw = ImageDraw.Draw(img)
        draw.text((10, 10), rhea_rxn, (0, 0, 0))
        images.append(img)
    images[0].save(output_file, save_all=True, append_images=images[1:])

# if __name__ == '__main__':
#     print(os.getenv(env_var_names.output_dir))
#     output_dir = Path(os.getenv(env_var_names.output_dir))
#     draw_reactions(output_dir.joinpath('enzsrp.csv'),
#                    output_dir.joinpath('reactions.pdf'))
