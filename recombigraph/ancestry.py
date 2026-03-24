from dataclasses import dataclass, replace
from typing import Optional

@dataclass
class Segment:
    left: float
    right: float
    parent_homolog_id: Optional[int]
    founder_homolog_id: int
    
    def copy(self, **kwargs):
        return replace(self, **kwargs)

@dataclass
class Homolog:
    homolog_id: int
    chromosome: str
    individual_id: str
    time: int
    length: float
    segments: list[Segment]

    def to_slot(self, slot_id):
        return Slot(
            slot_id=slot_id,
            homolog_id=self.homolog_id,
            chromosome=self.chromosome,
            length=self.length,
            segments=[
                seg.copy(parent_homolog_id=self.homolog_id)
                for seg in self.segments
            ],
        )

@dataclass
class Slot:
    """temp object for chromatids"""
    slot_id: int
    homolog_id: int
    chromosome: str
    length: float
    segments: list[Segment]

    def copy(self, **kwargs):
        return replace(self, **kwargs)