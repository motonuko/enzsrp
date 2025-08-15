import os
import warnings
from enum import Enum
from pathlib import Path
from typing import Optional, List, Set

import click
import pandas as pd
from rdkit.Chem import AllChem

from enzsrp.data.datasource.original.mcsa.original_mcsa_datasource import MCSAResidueSequence, OriginalMCSADataSource
from enzsrp.data.datasource.original.metacyc.parse_reactions_dat import MetaCycDirectionMapper
from enzsrp.data.datasource.original.rhea.original_rhea_datasource import OriginalRheaDatasource
from enzsrp.data.datasource.original.uniprot.id_mapping_data_source import IsoformIdMappingDataSource
from enzsrp.data.datasource.original.uniprot.original_uniprot_data_source import OriginalUniprotDataSource
from enzsrp.data.datasource.original.uniprot.uniprot_entity import BindingSite, EntryWithCatalyticActivity, \
    RheaReference, CatalyticActivity
from enzsrp.domain.entity.evidence_and_clusion_ontology import ECO
from enzsrp.domain.entity.reaction_direciton import ReactionDirection
from enzsrp.extension.parse_reaction import get_all_mol_names_in_rxn
from enzsrp.utils import env_var_names
from enzsrp.utils.list_utils import get_first_from_single_set


class DirectionSource(Enum):
    UNIPROT = "UniProt_Physiological_Reaction"
    METACYC = "MetaCyc"
    FORCE_L2R = "Force_left-to-right"


def data_formatter(data_id: str, primary_accession: str,
                   isoform_id: Optional[str],
                   ptm_key: Optional[str],
                   sequence: str,
                   ec_number: Optional[str],
                   has_physiological_rxn: bool,
                   rxn_has_exp_evidence: bool,
                   physiological_rxn_has_exp_evidence: Optional[bool],
                   rhea_master_id: str,
                   rhea_id: str,
                   rhea_rxn: AllChem.ChemicalReaction,
                   binding_site_all: List[BindingSite],
                   binding_site_exp: List[BindingSite],
                   mcsa_residues: Optional[Set[MCSAResidueSequence]],
                   direction_source: DirectionSource,
                   rxn_ecos: Optional[List[ECO]],
                   phy_rxn_ecos: Optional[List[ECO]]):
    rhea_rxn = AllChem.ReactionToSmiles(rhea_rxn)
    binding_site_all_sorted = sorted([site.location.string_notation for site in binding_site_all])
    binding_site_exp_sorted = sorted([site.location.string_notation for site in binding_site_exp])
    sorted_mcsa_residues = sorted(
        [str(residue.resid) for residue in mcsa_residues]) if mcsa_residues is not None else None
    return {
        "data_id": data_id,
        "primary_accession": primary_accession,
        "isoform_id": isoform_id,
        "ptm_key": ptm_key,
        "sequence": sequence,
        "ec_number": ec_number,
        "has_physiological_rxn": has_physiological_rxn,
        "physiological_rxn_has_exp_evidence": physiological_rxn_has_exp_evidence,
        "rxn_has_exp_evidence": rxn_has_exp_evidence,
        "rhea_master_id": rhea_master_id,
        "rhea_id": rhea_id,
        "rxn": rhea_rxn,
        "binding_site_all": "|".join(binding_site_all_sorted),
        "binding_site_exp": "|".join(binding_site_exp_sorted),
        "mcsa_residues": "|".join(sorted_mcsa_residues) if sorted_mcsa_residues is not None else None,
        "direction_source": direction_source.value,
        "rxn_evidence": '|'.join(sorted(list(set([eco.value for eco in rxn_ecos])))) if rxn_ecos else None,
        "phy_rxn_evidence": '|'.join(sorted(list(set([eco.value for eco in phy_rxn_ecos])))) if phy_rxn_ecos else None,
    }


def get_binding_site(binding_sites: List[BindingSite], rxn: AllChem.ChemicalReaction, only_exp: bool):
    if only_exp:
        return [site for site in binding_sites if
                any([evidence.type.value == ECO.EXPERIMENTAL.value for evidence in
                     site.evidences]) and site.ligand.ligand_id is not None and site.
                ligand.ligand_id.id in get_all_mol_names_in_rxn(rxn)]
    else:
        return [site for site in binding_sites if
                site.ligand.ligand_id is not None and site.ligand.ligand_id.id in get_all_mol_names_in_rxn(rxn)]


class ActivityDiscardReason(Enum):
    NO_EXP_EVIDENCE = 'No experimental evidence'
    NO_RHEA_REFERENCE = 'Rhea ID to RXN mapping failed'
    NO_DIRECTION = 'Direction could not be defined'


class DirectionDiscardReason(Enum):
    FAILED_RHEA_ID_TO_RXN_MAPPING = 'Failed to map rhea id to reaction'


class EnzymeReactionDatasetBuilder:
    ALLOWED_DIRECTION_SETS = [
        {ReactionDirection.LEFT_TO_RIGHT},
        {ReactionDirection.RIGHT_TO_LEFT},
        {ReactionDirection.LEFT_TO_RIGHT, ReactionDirection.RIGHT_TO_LEFT}
    ]

    def __init__(self, uniprot_entries_json: Path, uniprot_isoform_uniparc_mapping_json: Path, rhea_rxn_dir: Path,
                 rhea_2_metacyc_file: Path, rhea_directions_file: Path, mcsa_data_dir: Optional[Path],
                 metacyc_reactions_dat: Optional[Path], output_path: Path, use_undefined_direction_rxn: bool,
                 allow_non_exp_evidence: bool):
        self.source = OriginalUniprotDataSource(uniprot_entries_json)
        self.rhea_source = OriginalRheaDatasource(rhea_rxn_dir=rhea_rxn_dir, rhea2metacyc=rhea_2_metacyc_file,
                                                  rhea_directions=rhea_directions_file)
        self.isoform_source = IsoformIdMappingDataSource(uniprot_isoform_uniparc_mapping_json)

        self.mcsa_source = OriginalMCSADataSource(mcsa_data_dir) if mcsa_data_dir is not None else None
        self.metacyc_mapper = MetaCycDirectionMapper(
            metacyc_reactions_dat) if metacyc_reactions_dat is not None else None

        self.save_path = output_path
        self.use_undefined_direction_rxn = use_undefined_direction_rxn
        self.allow_non_exp_evidence = allow_non_exp_evidence

    def _mapping_isoform(self, entry: EntryWithCatalyticActivity, isoform_ids):
        assert isoform_ids is not None
        active_isoform_ids = self.isoform_source.bulk_map_isoform_id_to_active_isoform_id(isoform_ids)
        if len(set(active_isoform_ids.values())) == 0:
            warnings.warn('Unexpected. Are the download dates for UniProt data and UniParc data the same?')
            return None, None
        isoform_active_id = get_first_from_single_set(set(active_isoform_ids.values()))
        if isoform_active_id not in isoform_ids:
            warnings.warn(f"active_id is not present in the original isoform_ids. {isoform_ids} {isoform_active_id}")
        isoform_seq = self.isoform_source.map_isoform_to_seq(isoform_active_id)
        if len(entry.alternative_products_data.isoforms) == 1:
            if isoform_seq != entry.sequence:
                warnings.warn('Unexpected. If only one isoform exits, isoform sequence should be '
                              'equal to entry primary sequence.')
                return None, None
        return isoform_active_id, isoform_seq

    def get_isoform_sequence(self, entry, activity):
        assert activity in entry.catalytic_activities
        assert entry.alternative_products_data is not None
        isoform_ids = entry.alternative_products_data.search_isoform_by_isoform_name(activity.isoform_molecule_name)
        isoform_main_id, isoform_or_ptm_seq = self._mapping_isoform(entry, isoform_ids)
        return isoform_main_id, isoform_or_ptm_seq

    @staticmethod
    def get_post_translational_mod_seq(entry, activity: CatalyticActivity):
        assert activity in entry.catalytic_activities
        ptm_seq = entry.get_new_ptm_sequence(activity)
        return activity.molecule, ptm_seq

    def get_metacyc_directions(self, ref: RheaReference) -> List[ReactionDirection]:
        metacyc_ids = self.rhea_source.map_rhea_id_to_metacyc_id(ref.db_id.id_intstr)
        cyc_directions = [self.metacyc_mapper.metacyc_id_to_direction(metacyc_id) for metacyc_id in
                          metacyc_ids]
        cyc_directions = [cyc_direction for cyc_direction in cyc_directions if cyc_direction is not None]
        directions = set([direction for cyc_direction in cyc_directions for direction in
                          cyc_direction.to_reaction_directions()])
        directions = sorted(directions, key=lambda x: x.value)
        assert len(directions) == 0 or set(directions) in self.ALLOWED_DIRECTION_SETS
        return directions

    def construct_full_directed_dataset(self) -> pd.DataFrame:
        data = []
        skipped_activities = []
        skipped_directions = []
        for entry in self.source.stream_entries_with_catalytic_activity():
            mcsa_residues = None
            if self.mcsa_source is not None:
                mcsa_residues = self.mcsa_source.get_residues_by_entry(entry)
            for activity in entry.catalytic_activities:
                if self.allow_non_exp_evidence:
                    verified_phy_rxns = [phy_rxn for phy_rxn in activity.phy_reactions]
                else:
                    verified_phy_rxns = [phy_rxn for phy_rxn in activity.phy_reactions if
                                         phy_rxn.has_exp_evidence or activity.reaction.has_exp_evidence]
                direction_source = None
                if len(verified_phy_rxns) >= 1:
                    direction_source = DirectionSource.UNIPROT
                    for phy_rxn in verified_phy_rxns:
                        # Isoform or PTM processing
                        isoform_main_id, ptm_key, isoform_or_ptm_seq = None, None, None
                        if activity.has_isoform_molecule:
                            isoform_main_id, isoform_or_ptm_seq = self.get_isoform_sequence(entry, activity)
                        elif activity.has_non_isoform_molecule:
                            ptm_key, isoform_or_ptm_seq = self.get_post_translational_mod_seq(entry, activity)

                        rhea_id = self.rhea_source.map_master_id_to_direction_id(
                            activity.reaction.rhea_main_reference.db_id.id_intstr,
                            phy_rxn.direction.rhea_diction_name)
                        if not (rxn := self.rhea_source.get_single_reaction_with_cache(rhea_id)):
                            skipped_directions.append({
                                'accession': entry.primary_accession,
                                'reaction': activity.reaction.name,
                                'rhea_id': rhea_id,
                                'direction': phy_rxn.direction.short_name,
                                'reason': DirectionDiscardReason.FAILED_RHEA_ID_TO_RXN_MAPPING.value})
                            continue
                        data.append(data_formatter(
                            data_id=self.construct_unique_id(entry.primary_accession, isoform_main_id,
                                                             activity.reaction.rhea_main_reference, phy_rxn.direction),
                            primary_accession=entry.primary_accession,
                            isoform_id=isoform_main_id,
                            ptm_key=ptm_key,
                            sequence=isoform_or_ptm_seq if isoform_or_ptm_seq else entry.sequence,
                            ec_number=activity.reaction.ec_number,
                            has_physiological_rxn=True,
                            physiological_rxn_has_exp_evidence=phy_rxn.has_exp_evidence,
                            rxn_has_exp_evidence=activity.reaction.has_exp_evidence,
                            rhea_master_id=activity.reaction.rhea_main_reference.db_id.id_intstr,
                            rhea_id=rhea_id,
                            rhea_rxn=rxn,
                            binding_site_all=get_binding_site(entry.binding_sites, rxn,
                                                              only_exp=False),
                            binding_site_exp=get_binding_site(entry.binding_sites, rxn,
                                                              only_exp=True),
                            mcsa_residues=mcsa_residues,
                            direction_source=direction_source,
                            rxn_ecos=[evidence.type for evidence in activity.reaction.evidences],
                            phy_rxn_ecos=[evidence.type for evidence in activity.reaction.evidences]))
                else:
                    # Direction will be defined by MetaCyc data
                    if not activity.reaction.has_exp_evidence and not self.allow_non_exp_evidence:
                        # no exp evidence
                        skipped_activities.append({
                            'accession': entry.primary_accession,
                            'reaction': activity.reaction.name,
                            'reason': ActivityDiscardReason.NO_EXP_EVIDENCE.value})
                        continue
                    if not activity.reaction.has_rhea_rxn_reference:
                        # no Rhea reference
                        skipped_activities.append({
                            'accession': entry.primary_accession,
                            'reaction': activity.reaction.name,
                            'reason': ActivityDiscardReason.NO_RHEA_REFERENCE.value})
                        continue
                    # Isoform or PTM processing
                    isoform_main_id, ptm_key, isoform_or_ptm_seq = None, None, None
                    if activity.has_isoform_molecule:
                        isoform_main_id, isoform_or_ptm_seq = self.get_isoform_sequence(entry, activity)
                    elif activity.has_non_isoform_molecule:
                        ptm_key, isoform_or_ptm_seq = self.get_post_translational_mod_seq(entry, activity)

                    cyc_directions = []
                    if self.metacyc_mapper is not None:
                        ref = activity.reaction.rhea_main_reference
                        cyc_directions = self.get_metacyc_directions(ref)
                        if len(cyc_directions) != 0:
                            direction_source = DirectionSource.METACYC
                    directions = cyc_directions
                    if self.use_undefined_direction_rxn and len(directions) == 0:
                        directions.append(ReactionDirection.LEFT_TO_RIGHT)
                        direction_source = DirectionSource.FORCE_L2R
                    if len(directions) == 0:
                        # Reaction direction could not be defined
                        skipped_activities.append({
                            'accession': entry.primary_accession,
                            'reaction': activity.reaction.name,
                            'reason': ActivityDiscardReason.NO_DIRECTION.value})
                        continue
                    for direction in directions:
                        rhea_id = self.rhea_source.map_master_id_to_direction_id(
                            activity.reaction.rhea_main_reference.db_id.id_intstr,
                            direction.rhea_diction_name)
                        if not (rxn := self.rhea_source.get_single_reaction_with_cache(rhea_id)):
                            skipped_directions.append({
                                'accession': entry.primary_accession,
                                'reaction': activity.reaction.name,
                                'rhea_id': rhea_id,
                                'direction': direction.short_name,
                                'reason': DirectionDiscardReason.FAILED_RHEA_ID_TO_RXN_MAPPING.value})
                            continue
                        data.append(data_formatter(
                            data_id=self.construct_unique_id(entry.primary_accession,
                                                             isoform_main_id,
                                                             activity.reaction.rhea_main_reference,
                                                             direction),
                            primary_accession=entry.primary_accession,
                            isoform_id=isoform_main_id,
                            ptm_key=ptm_key,
                            sequence=isoform_or_ptm_seq if isoform_or_ptm_seq else entry.sequence,
                            ec_number=activity.reaction.ec_number,
                            has_physiological_rxn=False,
                            physiological_rxn_has_exp_evidence=None,
                            rxn_has_exp_evidence=activity.reaction.has_exp_evidence,
                            rhea_master_id=activity.reaction.rhea_main_reference.db_id.id_intstr,
                            rhea_id=rhea_id,
                            rhea_rxn=rxn,
                            binding_site_all=get_binding_site(entry.binding_sites, rxn,
                                                              only_exp=False),
                            binding_site_exp=get_binding_site(entry.binding_sites, rxn,
                                                              only_exp=True),
                            mcsa_residues=mcsa_residues,
                            direction_source=direction_source,
                            rxn_ecos=[evidence.type for evidence in activity.reaction.evidences],
                            phy_rxn_ecos=None))

        df = pd.DataFrame(data)
        assert df.duplicated().sum() == 0, "Unexpected behavior: duplicated rows has found"
        self.save_path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(self.save_path, index=False)

        df1 = pd.DataFrame(skipped_activities)
        df1.to_csv(self.save_path.parent.joinpath(f"skipped_activities_{self.save_path.stem}.csv"), index=False)
        df2 = pd.DataFrame(skipped_directions)
        df2.to_csv(self.save_path.parent.joinpath(f"skipped_directions_{self.save_path.stem}.csv"), index=False)
        return df

    @staticmethod
    def construct_unique_id(primary_accession: str, isoform_name: Optional[str], rhea_master_ref: RheaReference,
                            direction: ReactionDirection) -> str:
        isoform_text = isoform_name if isoform_name else "#"
        return f"{primary_accession}_{isoform_text}_{rhea_master_ref.db_id.id_intstr}_{direction.short_name}"

    # def print_statistics(self):
    #     df = pd.read_csv(self.save_path)
    #     nan_count = df['binding_site_all'].isna().sum()
    #     not_nan_count = df['binding_site_all'].notna().sum()
    #     print("binding site all NaN:", nan_count)
    #     print("binding site all not NaN:", not_nan_count)
    #     nan_count = df['binding_site_exp'].isna().sum()
    #     not_nan_count = df['binding_site_exp'].notna().sum()
    #     print("binding site exp NaN:", nan_count)
    #     print("binding site exp not NaN:", not_nan_count)


@click.command()
@click.option('--uniprot-entries-json', type=click.Path(exists=True),
              default=os.getenv(env_var_names.original_uniprot_reviewed_catalytic_activity_json_file),
              help='Specify the file path for the protein data downloaded from UniProt')
@click.option('--uniprot-isoform-uniparc-mapping-json', type=click.Path(exists=True),
              default=os.getenv(env_var_names.isoform_uniparc_id_mapping_file),
              help='Specify the file path for the downloaded isoform mapping file')
@click.option('--rhea-rxn-dir', type=click.Path(exists=True),
              default=os.getenv(env_var_names.original_rhea_rxn_dir),
              help='Specify the dir path that contains Rhea reactions rxn files')
@click.option('--rhea-2-metacyc-file', type=click.Path(exists=True),
              default=os.getenv(env_var_names.original_rhea2metacyc_file),
              help='Specify rhea2metacyc.tsv file path downloaded from Rhea')
@click.option('--rhea-directions-file', type=click.Path(exists=True),
              default=os.getenv(env_var_names.original_rhea_directions_file),
              help='Specify rhea-directions.tsv file path downloaded from Rhea')
@click.option('--mcsa-data-dir', type=click.Path(exists=True),
              default=os.getenv(env_var_names.original_mcsa_data_dir),
              help='Specify the M-CSA data dir downloaded from M-CSA database')
@click.option('--metacyc-reactions-dat', type=click.Path(exists=True),
              default=os.getenv(env_var_names.original_metacyc_reactions_dat_file),
              help='Specify reactions.dat file path downloaded from MetaCyc. If this path is not specified, '
                   'Some reactions won\'t be appeared in output file due to the lack of reaction direction information'
                   'You can set `--use-undefined-direction-rxn option to force to use all undefined direction reaction`')
@click.option('--output-dir', type=click.Path(),
              default=Path(os.getenv(env_var_names.output_dir)) if os.getenv(
                  env_var_names.output_dir) is not None else None,
              help='Specify output file path')
@click.option('--use-undefined-direction-rxn', is_flag=True,  # NOTE: not tested
              help='All undefined direction reactions will be treated as forward-direction reaction. This option WILL '
                   'PRODUCE MISLABELED REACTIONS, but the usage of undefined direction reactions possibly useful in some'
                   ' cases because of its data size')
@click.option('--allow-non-exp-evidence', is_flag=True,
              help='')
def build_enzyme_reaction_dataset(
        uniprot_entries_json: str,
        uniprot_isoform_uniparc_mapping_json: str,
        rhea_rxn_dir: str,
        rhea_2_metacyc_file: str,
        rhea_directions_file: str,
        mcsa_data_dir: str,
        metacyc_reactions_dat: str,
        output_dir: str,
        use_undefined_direction_rxn: bool,
        allow_non_exp_evidence,
):
    assert use_undefined_direction_rxn == False, use_undefined_direction_rxn
    print(metacyc_reactions_dat)
    output_path = Path(output_dir) / 'enzsrp_full.csv' if allow_non_exp_evidence else Path(output_dir) / 'enzsrp.csv'
    dataset_constructor = EnzymeReactionDatasetBuilder(
        uniprot_entries_json=Path(uniprot_entries_json),
        uniprot_isoform_uniparc_mapping_json=Path(uniprot_isoform_uniparc_mapping_json),
        rhea_rxn_dir=Path(rhea_rxn_dir),
        rhea_2_metacyc_file=Path(rhea_2_metacyc_file),
        rhea_directions_file=Path(rhea_directions_file),
        mcsa_data_dir=Path(mcsa_data_dir) if mcsa_data_dir is not None else None,
        metacyc_reactions_dat=Path(metacyc_reactions_dat) if metacyc_reactions_dat is not None else None,
        output_path=Path(output_path),
        use_undefined_direction_rxn=use_undefined_direction_rxn,
        allow_non_exp_evidence=allow_non_exp_evidence
    )
    dataset_constructor.construct_full_directed_dataset()
