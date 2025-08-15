import warnings
from pathlib import Path
from typing import Optional, Set

import pandas as pd
from rdkit.Chem import AllChem

from rxndata2.domain.entity.reaction_direciton import RheaDirectionName
from rxndata2.extension.parse_reaction import reaction_from_rxn_block_with_mol_title
from rxndata2.extension.rdkit_warning import RDKitWarningInterceptor, AmbiguousStereoChemistryWarningException, \
    Tagged3DBut2DMarkersFoundException


class OriginalRheaDatasource:

    def __init__(self, rhea_rxn_dir: Path, rhea2metacyc: Path, rhea_directions: Path):
        self.single_reaction_cache = {}
        self.rhea_rxn_dir = rhea_rxn_dir
        self.rhea_metacyc_df = pd.read_csv(rhea2metacyc, sep="\t", dtype={'RHEA_ID': str, 'MASTER_ID': str})
        self.rhea_id_map = pd.read_csv(rhea_directions, sep="\t", dtype=str)
        self.master_ids = self.rhea_id_map['RHEA_ID_MASTER'].values
        self.warned_rhea_ids = set()

    def map_master_id_to_direction_id(self, rhea_id: str, direction: RheaDirectionName):
        selected_row = self.rhea_id_map[self.rhea_id_map[RheaDirectionName.RHEA_ID_MASTER.value] == rhea_id]
        if len(selected_row) == 1:
            return selected_row[direction.value].iloc[0]
        elif len(selected_row) > 1:
            raise ValueError(f"Unexpected multiple hit while mapping rhea MASTER_ID {rhea_id}")
        else:
            raise ValueError(
                f"rhea ID: {rhea_id} is not listed in rhea MASTER_ID column in id table. It might be other ID (e.g. RHEA_ID_LR)")

    def get_single_reaction(self, rhea_id: str) -> Optional[AllChem.ChemicalReaction]:
        try:
            path = self.rhea_rxn_dir.joinpath(f"{rhea_id}.rxn")
            with open(path, 'r') as file:
                rxn_block = file.read()
                if rhea_id == '55541':
                    rxn_block = rxn_block.replace('93123', '93 123')
            interceptor = RDKitWarningInterceptor()
            # NOTE: You can ignore warnings by comment out following code, otherwise reactions with these warning won't be used.
            interceptor.set_ignore_warning_exceptions(
                [AmbiguousStereoChemistryWarningException, Tagged3DBut2DMarkersFoundException])
            with interceptor:
                result = reaction_from_rxn_block_with_mol_title(rxn_block)
            if rhea_id == '61542':
                print(result)
            return result
        except FileNotFoundError as e:
            # NOTE: If a molecule in ChEBI does not have a provided structure, the rxn file is likely missing
            # (the Rhea page also does not have a download button).
            warnings.warn(f"File not found Rhea ID: {rhea_id}. This Rhea ID won't be used. {e}")
            return None
        except ValueError as e:
            warnings.warn(f"Value error at rhea_id: {rhea_id}, {e}")
            self.warned_rhea_ids.add(rhea_id)
            return None
        except AmbiguousStereoChemistryWarningException as e:
            # NOTE: Some rxn files contain ambiguous stereo chemistry.
            # {'61542', '10929', '71497', '21306', '77973', '61538', '26431', '14790', '61537', '77972', '71504',
            # '13046', '61541', '61546', '26432', '11661', '71501', '71505', '10930', '13047', '11662', '71500',
            # '61545', '55541', '21305', '71496', '14791'}
            warnings.warn(f"rhea: {rhea_id}, message: {e}")
            self.warned_rhea_ids.add(rhea_id)
            return None
        # except NotRemovingHydrogenAtomWithoutNeighbors as e:
        #     return None
        # except NotRemovingHydrogenAtomWithDummyAtomNeighbors as e:
        #     return None
        except Exception as e:
            warnings.warn(f"rhea: {rhea_id}, message: {e}")
            self.warned_rhea_ids.add(rhea_id)
            return None

    def get_single_reaction_with_cache(self, rhea_id: str) -> Optional[AllChem.ChemicalReaction]:
        assert isinstance(rhea_id, str), "rhea_id should be `str` type"
        if rhea_id in self.single_reaction_cache:
            return self.single_reaction_cache[rhea_id]
        result = self.get_single_reaction(rhea_id)
        self.single_reaction_cache[rhea_id] = result
        return result

    # There may be multiple hits. If there is no reference from Rhea to MetaCyc, an empty set is returned.
    def map_rhea_id_to_metacyc_id(self, rhea_id: str) -> Set[str]:
        if rhea_id not in self.master_ids:
            warnings.warn('MASTER_ID is expected')
            return set()
        metacyc_ids = self.rhea_metacyc_df[self.rhea_metacyc_df['MASTER_ID'] == rhea_id]["ID"].values
        return set([metacyc_id.strip() for metacyc_id in metacyc_ids])
