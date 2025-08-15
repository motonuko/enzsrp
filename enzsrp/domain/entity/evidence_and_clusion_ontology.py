from enum import Enum


# Evidence and Conclusion Ontology
# see: https://www.uniprot.org/help/evidences
class ECO(Enum):
    # Evidence used only in manual assertions
    EXPERIMENTAL = 'ECO:0000269'
    NON_TRACEABLE_AUTHOR_STATEMENT = 'ECO:0000303'  # manually curated based on scientific literature without experimental evidence
    CURATOR_INFERENCE_EVIDENCE = 'ECO:0000305'
    SEQUENCE_SIMILARITY = 'ECO:0000250'

    # Evidence used in manual and automatic assertions
    SEQUENCE_MODEL_MANUAL = 'ECO:0000255'
    SEQUENCE_MODEL_AUTO = 'ECO:0000256'
    SEQUENCE_MODEL_AUTO2 = 'ECO:0000259'
    IMPORTED_INFO_MANUAL = 'ECO:0000312'
    IMPORTED_INFO_AUTO = 'ECO:0000313'
    COMBINATORIAL_MANUAL = 'ECO:0007744'
    COMBINATORIAL_AUTO = 'ECO:0007829'
    DEEP_LEARNING_AUTO = 'ECO:0008006'

    # GO annotations (It seems they are not used for rxn and binding site)
    # GO_EXP = 'ECO:0000269'  # same as EXPERIMENTAL
    GO_IDA = 'ECO:0000314'
    GO_IPI = 'ECO:0000353'
    GO_IMP = 'ECO:0000315'
    GO_IGI = 'ECO:0000316'
    GO_IEP = 'ECO:0000270'
    GO_HTP = 'ECO:0006056'
    GO_HDA = 'ECO:0007005'
    GO_HMP = 'ECO:0007001'
    GO_HGI = 'ECO:0007003'
    GO_HEP = 'ECO:0007007'
    GO_IBA = 'ECO:0000318'
    GO_IBD = 'ECO:0000319'
    GO_IKR = 'ECO:0000320'
    GO_IRD = 'ECO:0000321'

    @classmethod
    def get_eco_enum(cls, value):
        for eco_enum in ECO:
            if eco_enum.value == value:
                return eco_enum
        raise ValueError("No matched item")
