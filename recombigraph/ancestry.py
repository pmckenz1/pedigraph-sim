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
