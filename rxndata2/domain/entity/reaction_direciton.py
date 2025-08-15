from enum import Enum
from typing import Tuple


class RheaDirectionName(Enum):
    RHEA_ID_MASTER = "RHEA_ID_MASTER"
    RHEA_ID_LR = "RHEA_ID_LR"
    RHEA_ID_RL = "RHEA_ID_RL"
    RHEA_ID_BI = "RHEA_ID_BI"

    @classmethod
    def from_text(cls, text):
        for direction in cls:
            if direction.value == text:
                return direction
        raise ValueError(f"No matching enum value for {text}")


class ReactionDirection(Enum):
    LEFT_TO_RIGHT = "left-to-right"
    RIGHT_TO_LEFT = "right-to-left"

    @classmethod
    def from_string(cls, direction_str):
        for direction in cls:
            if direction.value == direction_str:
                return direction
        raise ValueError(f"Invalid direction string: {direction_str}")

    @property
    def rhea_diction_name(self):
        if self == ReactionDirection.LEFT_TO_RIGHT:
            return RheaDirectionName.RHEA_ID_LR
        elif self == ReactionDirection.RIGHT_TO_LEFT:
            return RheaDirectionName.RHEA_ID_RL

    @property
    def short_name(self):
        if self == ReactionDirection.LEFT_TO_RIGHT:
            return "l2r"
        elif self == ReactionDirection.RIGHT_TO_LEFT:
            return "r2l"


# see: https://biocyc.org/PGDBConceptsGuide.shtml#TAG:__tex2page_sec_4.2
class MetaCycDirection(Enum):
    REVERSIBLE = "REVERSIBLE"
    PHYSIOL_LEFT_TO_RIGHT = "PHYSIOL-LEFT-TO-RIGHT"
    PHYSIOL_RIGHT_TO_LEFT = "PHYSIOL-RIGHT-TO-LEFT"  # does not exist
    IRREVERSIBLE_LEFT_TO_RIGHT = "IRREVERSIBLE-LEFT-TO-RIGHT"  # does not exist
    IRREVERSIBLE_RIGHT_TO_LEFT = "IRREVERSIBLE-RIGHT-TO-LEFT"  # does not exist
    LEFT_TO_RIGHT = "LEFT-TO-RIGHT"
    RIGHT_TO_LEFT = "RIGHT-TO-LEFT"

    @classmethod
    def from_text(cls, text):
        for direction in cls:
            if direction.value == text:
                return direction
        raise ValueError(f"No matching enum value for {text}")

    def to_reaction_directions(self) -> Tuple[ReactionDirection, ...]:
        if self.value == MetaCycDirection.REVERSIBLE.value:
            return ReactionDirection.LEFT_TO_RIGHT, ReactionDirection.RIGHT_TO_LEFT

        elif (self.value == MetaCycDirection.LEFT_TO_RIGHT.value or
              self.value == MetaCycDirection.IRREVERSIBLE_LEFT_TO_RIGHT.value or
              self.value == MetaCycDirection.PHYSIOL_LEFT_TO_RIGHT.value):
            return ReactionDirection.LEFT_TO_RIGHT,

        elif (self.value == MetaCycDirection.RIGHT_TO_LEFT.value or
              self.value == MetaCycDirection.IRREVERSIBLE_RIGHT_TO_LEFT.value or
              self.value == MetaCycDirection.PHYSIOL_RIGHT_TO_LEFT.value):
            return ReactionDirection.RIGHT_TO_LEFT,
        else:
            raise ValueError("unexpected")
