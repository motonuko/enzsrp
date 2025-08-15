import click

from dotenv import load_dotenv

# NOTE: To use environment variables as click.option default values, .env file should be loaded before importing modules
load_dotenv()

from enzsrp.presentation.download_isoform_id_uniparc_mapping import \
    download_isoform_id_uniparc_mapping

from enzsrp.presentation.build_enzyme_reaction_dataset import build_enzyme_reaction_dataset


@click.group()
def create():
    pass


create.add_command(download_isoform_id_uniparc_mapping)  # type: ignore
create.add_command(build_enzyme_reaction_dataset)  # type: ignore

if __name__ == "__main__":
    create()
