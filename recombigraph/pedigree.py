from dataclasses import dataclass
from typing import Optional

@dataclass
class PedigreeRecord:
    name: str
    parent1: Optional[str]
    parent2: Optional[str]

class Pedigree:
    def __init__(self, records):
        self.records = records
