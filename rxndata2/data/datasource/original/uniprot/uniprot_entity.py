import re
import warnings
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
from typing import Optional, List, Union

from rxndata2.domain.entity.evidence_and_clusion_ontology import ECO
from rxndata2.domain.entity.reaction_direciton import ReactionDirection


@dataclass(frozen=True)
class Evidence:
    type: ECO
    evd_id: Optional[str]
    source: Optional[str]

    @classmethod
    def from_dict(cls, data: dict):
        eco_type = ECO.get_eco_enum(data["evidenceCode"])
        evd_id = data.get("id", None)
        source = data.get("source", None)
        return cls(eco_type, evd_id, source)

    def to_dict(self) -> dict:
        return {
            "type": self.type.value,
            "evd_id": self.evd_id,
            "source": self.source if self.source else None,
        }


@dataclass(frozen=True)
class RheaID:
    id: str
    pattern = r'^RHEA:\d+$'

    def __post_init__(self):
        if not re.match(self.pattern, self.id):
            raise ValueError("Invalid Rhea ID format.")

    @property
    def id_intstr(self) -> str:
        match = re.search(r'^RHEA:(\d+)$', self.id)
        if match:
            return match.group(1)
        else:
            raise ValueError("Rhea id has invalid format")


@dataclass(frozen=True)
class RheaCompID:
    id: str
    pattern = r'^RHEA-COMP:\d+$'

    def __post_init__(self):
        if not re.match(self.pattern, self.id):
            raise ValueError("Invalid Rhea ID format.")

    @property
    def id_intstr(self) -> str:
        match = re.search(r'^RHEA-COMP:(\d+)$', self.id)
        if match:
            return match.group(1)
        else:
            raise ValueError("Rhea id has invalid format")


def map_to_rhea_id_or_rhea_comp_id(text: str) -> Union[RheaID, RheaCompID]:
    try:
        return RheaID(text)
    except ValueError as e:
        try:
            return RheaCompID(text)
        except ValueError as e2:
            raise ValueError(str(e) + str(e2))


@dataclass(frozen=True)
class ChEBIID:
    """
    id should be something like this -> ChEBI:CHEBI:57540
    """
    id: str

    def __post_init__(self):
        pattern = r'^CHEBI:\d+$'
        if not re.match(pattern, self.id):
            raise ValueError(f"Invalid ChEBI ID format. {self.id}")

    @property
    def short_id(self) -> str:
        match = re.search(r'^CHEBI:(\d+)$', self.id)
        if match:
            return match.group(1)
        else:
            raise ValueError("ChEBI id has invalid format")


# Parent class
@dataclass(frozen=True)
class CrossReference:
    database: str
    db_id: str


@dataclass(frozen=True)
class RheaReference(CrossReference):
    database = 'Rhea'
    db_id: Union[RheaID, RheaCompID]

    @classmethod
    def from_dict(cls, data: dict) -> "RheaReference":
        assert data["database"] == cls.database, print(data)
        return cls(cls.database, db_id=map_to_rhea_id_or_rhea_comp_id(data["id"]))

    @property
    def is_comp(self):
        return isinstance(self.db_id, RheaCompID)


@dataclass(frozen=True)
class ChEBIReference(CrossReference):
    database = 'ChEBI'
    db_id: ChEBIID

    @classmethod
    def from_dict(cls, data: dict) -> "ChEBIReference":
        assert data["database"] == cls.database, data
        return cls(cls.database, db_id=ChEBIID(data["id"]))


@dataclass(frozen=True)
class Reaction:
    ec_number: Optional[str]
    evidences: List[Evidence]  # if evidences is None, it will be empty list
    cross_references: List[CrossReference]  # only rhea or chebi
    name: str

    @classmethod
    def from_dict(cls, data: dict):
        evidences = data.get("evidences", [])
        evidences = [Evidence.from_dict(evidence) for evidence in evidences]
        # NOTE: Information for RHEA-COMP and CHEBI is also defined in the Rhea DB rxn files
        # and can be supplemented from there.
        if "reactionCrossReferences" not in data:
            return cls(ec_number=data.get("ecNumber", None), evidences=evidences, cross_references=[],
                       name=data["name"])

        cross_references = []
        for ref in data["reactionCrossReferences"]:
            if ref["database"] == RheaReference.database:
                cross_references.append(RheaReference.from_dict(ref))
            elif ref["database"] == ChEBIReference.database:
                cross_references.append(ChEBIReference.from_dict(ref))
            else:
                raise RuntimeError(ref)
        return cls(ec_number=data.get("ecNumber", None), evidences=evidences, cross_references=cross_references,
                   name=data["name"])

    @property
    def has_rhea_rxn_reference(self) -> bool:
        return any(ref.database == "Rhea" for ref in self.cross_references)

    @property
    def rhea_main_reference(self) -> Optional[RheaReference]:
        al = [ref for ref in self.cross_references if isinstance(ref, RheaReference) and not ref.is_comp]
        if len(al) == 1:
            return al[0]
        elif len(al) == 0:
            return None
        else:
            raise ValueError("unexpected")

    @property
    def has_exp_evidence(self) -> bool:
        return any([evidence.type.value == ECO.EXPERIMENTAL.value for evidence in self.evidences])


@dataclass(frozen=True)
class PhysiologicalReaction:
    direction: ReactionDirection
    evidences: List[Evidence]
    cross_references: List[CrossReference]

    @classmethod
    def from_dict(cls, data: dict):

        evidences = data.get("evidences", [])
        evidences = [Evidence.from_dict(evidence) for evidence in evidences]

        if "reactionCrossReferences" not in data:
            return cls(direction=ReactionDirection.from_string(data["directionType"]), evidences=evidences,
                       cross_references=[])

        cross_references = []
        for ref in data["reactionCrossReferences"]:
            if ref["database"] == RheaReference.database:
                cross_references.append(RheaReference.from_dict(ref))
            elif ref["database"] == ChEBIReference.database:
                cross_references.append(ChEBIReference.from_dict(ref))
            else:
                raise RuntimeError(ref)
        return cls(direction=ReactionDirection.from_string(data["directionType"]), evidences=evidences,
                   cross_references=cross_references)

    @property
    def has_exp_evidence(self) -> bool:
        return any([evidence.type.value == ECO.EXPERIMENTAL.value for evidence in self.evidences])


@dataclass(frozen=True)
class CatalyticActivity:
    molecule: Optional[str]  # may include references to isoforms, e.g., https://rest.uniprot.org/uniprotkb/A0PJK1.json
    reaction: Reaction
    phy_reactions: List[PhysiologicalReaction]

    def __post_init__(self):
        assert len(self.phy_reactions) <= 2
        assert len(set([rxn.direction for rxn in self.phy_reactions])) == len(self.phy_reactions)

    @classmethod
    def from_dict(cls, data: dict):
        molecule = data.get("molecule")
        reaction = Reaction.from_dict(data["reaction"])
        phy_reaction = []
        if "physiologicalReactions" in data:
            phy_reaction = [PhysiologicalReaction.from_dict(phy) for phy in data["physiologicalReactions"]]
        return cls(molecule, reaction, phy_reaction)

    @property
    def isoform_molecule_name(self) -> Optional[str]:
        if self.molecule is None or "Isoform" not in self.molecule:
            return None
        return self.molecule.replace("Isoform ", "").strip()

    @property
    def has_isoform_molecule(self) -> Optional[str]:
        return self.molecule is not None and "Isoform" in self.molecule

    @property
    def has_non_isoform_molecule(self) -> Optional[str]:
        return self.molecule is not None and "Isoform" not in self.molecule


class LocationModifier(Enum):
    EXACT = "EXACT"
    OUTSIDE = "OUTSIDE"
    UNSURE = "UNSURE"
    UNKNOWN = "UNKNOWN"

    @classmethod
    def from_string(cls, text):
        for item in cls:
            if item.value == text:
                return item
        raise ValueError(f"Invalid string: {text}")


@dataclass(frozen=True)
class Location:
    start_value: int  # starts form '1'
    start_modifier: LocationModifier
    end_value: int
    end_modifier: LocationModifier

    @classmethod
    def from_dict(cls, data: dict) -> "Location":
        return cls(
            start_value=data["start"]["value"],
            start_modifier=LocationModifier.from_string(data["start"]["modifier"]),
            end_value=data["end"]["value"],
            end_modifier=LocationModifier.from_string(data["end"]["modifier"]),
        )

    @property
    def string_notation(self) -> str:
        if self.start_modifier.value == "OUTSIDE" or self.end_modifier.value == "OUTSIDE":
            warnings.warn("This entry has OUTSIDE location annotation")
        if self.start_value == self.end_value:
            return str(self.start_value)
        else:
            return f"{self.start_value}-{self.end_value}"


@dataclass(frozen=True)
class Ligand:
    name: str
    ligand_id: Optional[ChEBIID]
    label: Optional[str]
    note: Optional[str]

    @classmethod
    def from_dict(cls, data: dict) -> "Ligand":
        assert not (set(data.keys()) - {"name", "id", "label", "note"}), data
        if data.get("id") is None:
            assert data["name"] == "substrate"  # id is None only when the name is "substrate"

        # The format of Ligand ID is like: ChEBI:CHEBI:nnnnn
        return cls(
            name=data["name"],
            ligand_id=ChEBIID(data.get("id").split(":", 1)[1]) if data.get("id") else None,
            label=data.get("label"),
            note=data.get("note"))

    @property
    def ligand_id_simple(self):
        return self.ligand_id.short_id if self.ligand_id else None


# For cases where heme has iron, e.g., A0A075TMP8
@dataclass(frozen=True)
class LigandPart:
    name: str
    ligand_id: Optional[str]
    label: Optional[str]
    note: Optional[str]

    @classmethod
    def from_dict(cls, data: dict) -> "LigandPart":
        assert not (set(data.keys()) - {"name", "id", "label", "note"}), data
        return cls(
            name=data["name"],
            ligand_id=data.get("id"),
            label=data.get("label"),
            note=data.get("note"))

    @property
    def ligand_id_simple(self):
        return self.ligand_id.split(':', 1)[-1].strip() if self.ligand_id else None


@dataclass(frozen=True)
class BindingSite:
    type: str
    location: Location
    featureCrossReference: List[ChEBIReference]
    description: str
    evidences: List[Evidence]
    ligand: Ligand
    ligandPart: Optional[LigandPart]

    @classmethod
    def from_dict(cls, data: dict) -> "BindingSite":
        assert not (set(data.keys()) - {"type", "location", "featureCrossReferences", "evidences", "ligand",
                                        "description", "ligandPart"}), data
        return cls(
            type=data["type"],
            location=Location.from_dict(data["location"]),
            featureCrossReference=[ChEBIReference.from_dict(f) for f in
                                   data["featureCrossReferences"]] if data.get("featureCrossReferences") else [],
            description=data["description"],
            evidences=[Evidence.from_dict(f) for f in data["evidences"]] if data.get("evidences") else [],
            ligand=Ligand.from_dict(data["ligand"]),
            ligandPart=LigandPart.from_dict(data.get("ligandPart")) if data.get("ligandPart") else None,
        )

    @property
    def has_exp_evidence(self) -> bool:
        return any([evidence.type.value == ECO.EXPERIMENTAL.value for evidence in self.evidences])


@dataclass(frozen=True)
class ActiveSite:
    type: str
    location: Location
    featureCrossReference: List[ChEBIReference]
    description: str
    evidences: List[Evidence]

    @classmethod
    def from_dict(cls, data: dict) -> "ActiveSite":
        location = Location.from_dict(data.pop("location"))
        feature_cross_reference = [ChEBIReference.from_dict(ref) for ref in data.pop("featureCrossReferences", [])]
        evidences = [Evidence.from_dict(ref) for ref in data.pop("evidences", [])]
        return cls(
            location=location,
            featureCrossReference=feature_cross_reference,
            evidences=evidences,
            **data
        )

    @property
    def has_exp_evidence(self) -> bool:
        return any([evidence.type.value == ECO.EXPERIMENTAL.value for evidence in self.evidences])


@dataclass(frozen=True)
class Site:
    type: str
    location: Location
    featureCrossReference: List[ChEBIReference]
    description: str
    evidences: List[Evidence]

    @classmethod
    def from_dict(cls, data: dict) -> "Site":
        location = Location.from_dict(data.pop("location"))
        feature_cross_reference = [ChEBIReference.from_dict(ref) for ref in data.pop("featureCrossReferences", [])]
        evidences = [Evidence.from_dict(ref) for ref in data.pop("evidences", [])]
        return cls(
            location=location,
            featureCrossReference=feature_cross_reference,
            evidences=evidences,
            **data
        )

    @property
    def has_exp_evidence(self) -> bool:
        return any([evidence.type.value == ECO.EXPERIMENTAL.value for evidence in self.evidences])


@dataclass(frozen=True)
class Feature:
    type: str
    location: Optional[Location]
    description: str

    @classmethod
    def from_dict(cls, data: dict) -> "Feature":
        location = data.get("location")
        if location is not None:
            location = Location.from_dict(location)
        return cls(
            type=data['type'],
            location=location,
            description=data['description']
        )


@dataclass(frozen=True)
class Cofactor:
    name: str
    evidences: List[Evidence]
    cofactorCrossReference: ChEBIReference

    @classmethod
    def from_dict(cls, data):
        return cls(
            name=data['name'],
            evidences=[Evidence.from_dict(e) for e in data.get('evidences', [])],
            cofactorCrossReference=ChEBIReference.from_dict(data['cofactorCrossReference'])
        )


@dataclass(frozen=True)
class NoteText:
    evidences: List[Evidence]
    value: str

    @classmethod
    def from_dict(cls, data):
        return cls(
            evidences=[Evidence.from_dict(e) for e in data.get('evidences', [])],
            value=data['value']
        )


@dataclass(frozen=True)
class CofactorData:
    cofactors: List[Cofactor]
    note_texts: List[NoteText]

    @classmethod
    def from_dict(cls, data):
        return cls(
            cofactors=[Cofactor.from_dict(cofactor) for cofactor in data.get('cofactors', [])],
            note_texts=[NoteText.from_dict(text) for text in data.get('note', {}).get('texts', [])]
        )


@dataclass(frozen=True)
class Isoform:
    name: str
    isoformIds: List[str]
    isoformSequenceStatus: str  # Probably "displayed" indicates the canonical (main) sequence.
    sequenceIds: Optional[List[str]] = field(default_factory=list)  # possibly non-optional

    @staticmethod
    def from_dict(data: dict) -> 'Isoform':
        return Isoform(
            name=data['name']['value'],
            isoformIds=data['isoformIds'],
            isoformSequenceStatus=data['isoformSequenceStatus'],
            sequenceIds=data.get('sequenceIds', [])
        )


@dataclass(frozen=True)
class AlternativeProductsData:
    commentType: str
    events: List[str]
    isoforms: List[Isoform]  # possibly len == 1

    @staticmethod
    def from_dict(data: dict) -> 'AlternativeProductsData':
        isoforms = [Isoform.from_dict(isoform) for isoform in data["isoforms"]]
        return AlternativeProductsData(
            commentType=data["commentType"],
            events=data["events"],
            isoforms=isoforms
        )

    # iso_name is like: activity.isoform_molecule_name
    # returns multiple isoform ids
    def search_isoform_by_isoform_name(self, query_isoform_name: str) -> Optional[list[str]]:
        hit_isoforms = [isoform.isoformIds for isoform in self.isoforms if isoform.name == query_isoform_name]
        assert len(hit_isoforms) <= 1, hit_isoforms
        return hit_isoforms[0] if len(hit_isoforms) == 1 else None


@dataclass
class EntryWithCatalyticActivity:
    primary_accession: str
    secondary_accessions: List[str]  # NOTE: if it is not defined, it will be empty list
    uniprot_kb_id: str
    protein_existence: str
    catalytic_activities: List[CatalyticActivity]
    binding_sites: List[BindingSite]
    active_sites: List[ActiveSite]
    sites: List[Site]
    features: List[Feature]
    sequence: str
    last_sequence_update_date: datetime
    cofactor_data_list: List[CofactorData]
    alternative_products_data: Optional[AlternativeProductsData]

    def __post_init__(self):
        pass

    @classmethod
    def from_dict(cls, data: dict):
        activities = [comment for comment in data["comments"] if comment["commentType"] == "CATALYTIC ACTIVITY"]
        activities = [CatalyticActivity.from_dict(activity) for activity in activities]
        assert len(activities) != 0

        binding_sites = [feature for feature in data["features"] if feature["type"] == "Binding site"]
        binding_sites = [BindingSite.from_dict(site) for site in binding_sites]

        active_sites = [feature for feature in data["features"] if feature["type"] == "Active site"]
        active_sites = [ActiveSite.from_dict(site) for site in active_sites]

        sites = [feature for feature in data["features"] if feature["type"] == "Site"]
        sites = [Site.from_dict(site) for site in sites]

        features = [feature for feature in data["features"]]
        features = [Feature.from_dict(feature) for feature in features]

        cofactor_data = [CofactorData.from_dict(comment) for comment in data["comments"] if
                         comment["commentType"] == "COFACTOR"]

        alternative = [AlternativeProductsData.from_dict(comment) for comment in data["comments"] if
                       comment["commentType"] == "ALTERNATIVE PRODUCTS"]
        assert len(alternative) <= 1
        alternative = alternative[0] if alternative else None

        # ASSERT
        for activity in activities:
            if activity.molecule is None:
                continue
            if activity.has_isoform_molecule:
                continue
            match_features = [feature for feature in features if feature.description == activity.molecule]
            if len(match_features) > 0:
                types = [f.type for f in match_features]
                assert len(types) == 1 and types[0] == 'Chain'
                continue
            raise ValueError()

        for activity in activities:
            if activity.has_isoform_molecule:
                assert alternative is not None, activity.molecule

        return cls(primary_accession=data["primaryAccession"],
                   secondary_accessions=data.get("secondaryAccessions", []),
                   uniprot_kb_id=data["uniProtkbId"],
                   protein_existence=data["proteinExistence"],
                   catalytic_activities=activities,
                   binding_sites=binding_sites,
                   active_sites=active_sites,
                   sites=sites,
                   features=features,
                   sequence=data["sequence"]["value"],
                   last_sequence_update_date=datetime.strptime(data["entryAudit"]["lastSequenceUpdateDate"],
                                                               "%Y-%m-%d"),
                   cofactor_data_list=cofactor_data,
                   alternative_products_data=alternative)

    def get_new_ptm_sequence(self, target_activity: CatalyticActivity) -> Optional[str]:
        assert target_activity in self.catalytic_activities
        assert target_activity.has_non_isoform_molecule
        assert target_activity.molecule in [f.description for f in self.features]
        location = [f.location for f in self.features if f.description == target_activity.molecule][0]
        if location.start_modifier == LocationModifier.EXACT and location.end_modifier == LocationModifier.EXACT:
            return self.sequence[location.start_value - 1:location.end_value]
        return None
