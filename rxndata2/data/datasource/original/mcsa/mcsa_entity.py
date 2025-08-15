from dataclasses import dataclass
from typing import List, Optional


@dataclass
class Role:
    group_function: str
    function_type: str
    function: str
    emo: str

    @classmethod
    def from_dict(cls, data):
        return cls(**data)


@dataclass
class ResidueChain:
    chain_name: str
    pdb_id: str
    assembly_chain_name: str
    assembly: str
    code: str
    resid: int
    auth_resid: int
    is_reference: bool
    domain_name: str
    domain_cath_id: str

    @classmethod
    def from_dict(cls, data):
        return cls(**data)


@dataclass
class ResidueSequence:
    uniprot_id: Optional[str]  # NOTE: Some residue is not linked to UniProt
    code: str  # e.g. 'Asp'
    is_reference: bool
    resid: int  # e.g. 7

    @classmethod
    def from_dict(cls, data):
        uniprot_id = data.pop('uniprot_id', "")
        uniprot_id = uniprot_id.strip()
        if uniprot_id != "":
            assert uniprot_id.isalnum(), print(uniprot_id)  # Check if multiple IDs are not concatenated with a comma
        if uniprot_id == "":
            uniprot_id = None
        return cls(uniprot_id=uniprot_id, **data)


@dataclass
class Residue:
    mcsa_id: str
    roles_summary: str
    function_location_abv: str
    main_annotation: str
    ptm: str
    roles: List[Role]
    residue_chains: List[ResidueChain]  # Likely associated with PDB; can sometimes be 0.
    residue_sequences: ResidueSequence  # Likely associated with UniProt; will never be 0.

    @classmethod
    def from_dict(cls, data):
        assert len(data['residue_sequences']) == 1
        assert len(data['residue_chains']) <= 1
        roles_data = data.pop('roles', [])
        residue_chains_data = data.pop('residue_chains', [])
        residue_sequences_data = data.pop('residue_sequences', [])

        roles = [Role.from_dict(role_data) for role_data in roles_data]
        residue_chains = [ResidueChain.from_dict(chain_data) for chain_data in residue_chains_data]
        residue_sequences = ResidueSequence.from_dict(residue_sequences_data[0])
        return cls(roles=roles, residue_chains=residue_chains, residue_sequences=residue_sequences, **data)


@dataclass
class Protein:
    sequences: List[str]

    @classmethod
    def from_dict(cls, data):
        return cls(**data)


@dataclass
class McsaEntry:
    mcsa_id: str
    enzyme_name: str
    is_reference_uniprot_id: bool
    reference_uniprot_id: str
    url: str
    description: str
    protein: Protein
    all_ecs: List[str]
    residues: List[Residue]
    reaction: dict

    @classmethod
    def from_dict(cls, data):
        assert data["is_reference_uniprot_id"]
        assert data["reference_uniprot_id"] is not None and len(data["reference_uniprot_id"]) > 0
        protein_data = data.pop('protein', {})
        protein = Protein.from_dict(protein_data)
        residues_data = data.pop('residues', {})
        residues_data = [Residue.from_dict(residue) for residue in residues_data]
        return cls(protein=protein, residues=residues_data, **data)
