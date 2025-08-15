import os
from datetime import datetime
from pathlib import Path
from typing import Optional

import click

from enzsrp.data.datasource.original.uniprot.original_uniprot_data_source import OriginalUniprotDataSource
from enzsrp.data.datasource.remote.uniprot_id_mapping_ext import submit_id_mapping_task_and_download_file
from enzsrp.utils import env_var_names


def _download_isoform_id_uniparc_mapping(uniprot_entries_json: Path, output_dir: Path):
    datasource = OriginalUniprotDataSource(uniprot_entries_json)
    isoform_ids = datasource.get_all_isoform_ids()
    formatted_date = datetime.today().strftime('%Y_%m_%d')
    output_file_path = output_dir.joinpath(f"idmapping_{formatted_date}_isoform_uniparc.json")
    output_file_path.parent.mkdir(parents=True, exist_ok=True)
    submit_id_mapping_task_and_download_file(from_db='UniProtKB_AC-ID', to_db='UniParc', ids=isoform_ids,
                                             output_file_path=output_file_path)


@click.command()
@click.option('--uniprot-entries-json', type=click.Path(exists=True),
              default=os.getenv(env_var_names.original_uniprot_reviewed_catalytic_activity_json_file),
              help='Specify the file path for the protein data downloaded from UniProt'
                   'The protein isoforms with catalytic activity in this file will be mapped to their UniProt IDs.')
@click.option('--output-dir', type=click.Path(),
              default=os.getenv(env_var_names.output_dir),
              help='Specify the directory where the output file will be saved.')
def download_isoform_id_uniparc_mapping(uniprot_entries_json: Optional[str], output_dir: Optional[str]):
    assert uniprot_entries_json is not None, \
        (f"Please define {env_var_names.original_uniprot_reviewed_catalytic_activity_json_file} in .env file "
         f"or pass --output-dir argument")
    assert output_dir is not None, f"Please define {env_var_names.output_dir} in .env file or pass --output-dir argument"
    _download_isoform_id_uniparc_mapping(uniprot_entries_json=Path(uniprot_entries_json), output_dir=Path(output_dir))
