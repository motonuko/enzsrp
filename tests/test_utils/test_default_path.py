import os
import warnings
from functools import cached_property
from pathlib import Path
from typing import Optional

from dotenv import load_dotenv

from rxndata2.utils import env_var_names

load_dotenv()


def safe_load_env_var_path(env_variable: str) -> Optional[Path]:
    value: Optional[str] = os.getenv(env_variable)
    if value is not None:
        return Path(value)
    else:
        warnings.warn(f"Warning: {env_variable} is not defined. Returning None.")
        return None


class TestDefaultPath:
    _instance = None

    # singleton
    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super().__new__(cls, *args, **kwargs)
        return cls._instance

    @cached_property
    def original_uniprot_reviewed_catalytic_activity_json_file(self) -> Optional[Path]:
        return safe_load_env_var_path(env_var_names.original_uniprot_reviewed_catalytic_activity_json_file)

    @cached_property
    def original_rhea_directions_file(self) -> Optional[Path]:
        return safe_load_env_var_path(env_var_names.original_rhea_directions_file)

    @cached_property
    def original_rhea2metacyc_file(self) -> Optional[Path]:
        return safe_load_env_var_path(env_var_names.original_rhea2metacyc_file)

    @cached_property
    def original_rhea_rxn_dir(self) -> Optional[Path]:
        return safe_load_env_var_path(env_var_names.original_rhea_rxn_dir)

    @cached_property
    def original_metacyc_reactions_dat(self) -> Optional[Path]:
        return safe_load_env_var_path(env_var_names.original_metacyc_reactions_dat_file)

    @cached_property
    def output_dir(self) -> Optional[Path]:
        return safe_load_env_var_path(env_var_names.output_dir)

    @cached_property
    def test_module_dir(self) -> Optional[Path]:
        return safe_load_env_var_path("ENZSRP_TEST_DIR")

    @cached_property
    def test_data(self) -> Optional[Path]:
        test_dir = self.test_module_dir
        if test_dir is None:
            return None
        return test_dir / 'data'
