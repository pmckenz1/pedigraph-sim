from dataclasses import dataclass

@dataclass
class ChromosomeSpec:
    name: str
    length_cm: float
    centromere_cm: float

@dataclass
class GenomeSpec:
    chromosomes: list[ChromosomeSpec]
    ploidy: int = 2
    map_function: str = "HALDANE"
