import json
from typing import List

from rxndata2.data.datasource.remote.uniprot_id_mapping import submit_id_mapping, check_id_mapping_results_ready, \
    get_id_mapping_results_link, get_id_mapping_results_search


def submit_id_mapping_task_and_download_file(from_db: str, to_db: str, ids: List[str], output_file_path):
    job_id = submit_id_mapping(from_db, to_db, ids)
    print(f"jobid: {job_id}")
    if check_id_mapping_results_ready(job_id):
        link = get_id_mapping_results_link(job_id)
        results = get_id_mapping_results_search(link)
        with open(output_file_path, 'w') as f:
            json.dump(results, f, ensure_ascii=False, indent=4)
    print(f"you can also check result in: https://www.uniprot.org/id-mapping/uniparc/{job_id}/overview")
