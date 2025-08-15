import json
from pathlib import Path
from unittest import TestCase
from unittest.mock import patch, Mock, mock_open

from rxndata2.data.datasource.remote.uniprot_id_mapping import API_URL
from rxndata2.data.datasource.remote.uniprot_id_mapping_ext import submit_id_mapping_task_and_download_file


class UniprotIdMappingExtTest(TestCase):
    # NOTE: the order of patches and arguments is reversed
    @patch('builtins.open', new_callable=mock_open)
    @patch('rxndata2.data.datasource.remote.uniprot_id_mapping.requests.post')
    @patch('rxndata2.data.datasource.remote.uniprot_id_mapping.requests.Session.get')
    @patch('rxndata2.data.datasource.remote.uniprot_id_mapping_ext.get_id_mapping_results_search')
    def test_submit_id_mapping_task_and_download_file(self, mock_get_id_mapping_results_search, mock_session_get,
                                                      mock_requests_post, mock_open_func):
        mock_job_id = 'mock_job_id'
        sample_url = 'https://example.com/'
        result_content = {'x': 'y'}
        expected_json = json.dumps(result_content, ensure_ascii=False, indent=4)
        write_path = Path('sample.json')

        # Define a response function
        # noinspection PyUnusedLocal
        def mock_post(url, *args, **kwargs):
            if url == f"{API_URL}/idmapping/run":
                mock_response = Mock()
                mock_response.status_code = 200
                mock_response.json.return_value = {'jobId': mock_job_id}
                return mock_response
            raise ValueError('unexpected')

        # noinspection PyUnusedLocal
        def mock_session_get_response(url, *args, **kwargs):
            if url == f"{API_URL}/idmapping/status/{mock_job_id}":
                mock_response = Mock()
                mock_response.status_code = 200
                mock_response.json.return_value = {'results': 'X'}
                return mock_response
            elif url == f"{API_URL}/idmapping/details/{mock_job_id}":
                mock_response = Mock()
                mock_response.status_code = 200
                mock_response.json.return_value = {'redirectURL': sample_url}
                return mock_response
            elif sample_url in url:
                mock_response = Mock()
                mock_response.status_code = 200
                mock_response.json.return_value = {'redirectURL': sample_url}
                return mock_response

            raise ValueError('unexpected')

        def mock_get_id_mapping_results_search_method(url):
            self.assertEqual(url, sample_url)
            return result_content

        mock_requests_post.side_effect = mock_post
        mock_session_get.side_effect = mock_session_get_response
        mock_get_id_mapping_results_search.side_effect = mock_get_id_mapping_results_search_method

        submit_id_mapping_task_and_download_file(from_db='UniProtKB_AC-ID', to_db='UniParc', ids=['SAMPLE_ISOFORM'],
                                                 output_file_path=write_path)
        written_content = ''.join(call.args[0] for call in mock_open_func().write.call_args_list)
        self.assertEqual(written_content, expected_json)
