from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Any

@dataclass(frozen=True)
class LocalForest:
    chromosome: str
    left: float
    right: float
    edges: frozenset[tuple[int, int]]
    sample_homolog_ids: tuple[int, ...]

    def nodes(self) -> tuple[int, ...]:
        """All node IDs present in this local forest."""
        out = set(self.sample_homolog_ids)
        for parent, child in self.edges:
            out.add(parent)
            out.add(child)
        return tuple(sorted(out))

    def roots(self) -> tuple[int, ...]:
        """Nodes with no parent within this local forest."""
        nodes = set(self.nodes())
        children = {child for _, child in self.edges}
        return tuple(sorted(nodes - children))

    def children_map(self) -> dict[int, tuple[int, ...]]:
        """Parent -> children mapping."""
        out: dict[int, list[int]] = {}
        for parent, child in self.edges:
            out.setdefault(parent, []).append(child)
        return {k: tuple(sorted(v)) for k, v in out.items()}

    def parent_map(self) -> dict[int, int]:
        """Child -> parent mapping."""
        return {child: parent for parent, child in self.edges}

    def span(self) -> float:
        return self.right - self.left

    def is_empty(self) -> bool:
        return self.right <= self.left


@dataclass(frozen=True)
class LocalForestSequence:
    chromosome: str
    forests: tuple[LocalForest, ...]

    def __iter__(self):
        return iter(self.forests)

    def __len__(self) -> int:
        return len(self.forests)

    def __getitem__(self, idx):
        return self.forests[idx]

    @property
    def left(self) -> float:
        return self.forests[0].left if self.forests else 0.0

    @property
    def right(self) -> float:
        return self.forests[-1].right if self.forests else 0.0

    def breakpoints(self) -> tuple[float, ...]:
        if not self.forests:
            return tuple()
        vals = [self.forests[0].left]
        vals.extend(f.right for f in self.forests)
        return tuple(vals)

    def sample_homolog_ids(self) -> tuple[int, ...]:
        if not self.forests:
            return tuple()
        return self.forests[0].sample_homolog_ids

    def nodes(self) -> tuple[int, ...]:
        out = set()
        for forest in self.forests:
            out.update(forest.nodes())
        return tuple(sorted(out))

def _build_homolog_lookup(result) -> dict[int, Any]:
    out = {}
    for individual in result.individuals.values():
        for homologs in individual.homologs_by_chromosome.values():
            for homolog in homologs:
                out[homolog.homolog_id] = homolog
    return out


def _segment_covering(homolog, position: float):
    for seg in homolog.segments:
        if seg.left <= position < seg.right:
            return seg
    return None


def _local_forest_at_position(result, chromosome: str, sample_homolog_ids, position: float):
    """
    Return a set of (parent_homolog_id, child_homolog_id) edges describing
    the local forest at one genomic position.
    """
    homolog_lookup = _build_homolog_lookup(result)
    seen = set()
    edges = set()

    def visit(hid: int):
        if hid in seen:
            return
        seen.add(hid)

        h = homolog_lookup[hid]
        if h.chromosome != chromosome:
            return

        seg = _segment_covering(h, position)
        if seg is None:
            return

        parent_id = seg.parent_homolog_id
        if parent_id is not None:
            edges.add((parent_id, hid))
            visit(parent_id)

    for hid in sample_homolog_ids:
        visit(hid)

    return edges


def _ancestral_breakpoints(result, chromosome: str, sample_homolog_ids):
    """
    Collect all breakpoint coordinates encountered anywhere in the ancestry
    of the sampled homologs on this chromosome.
    """
    homolog_lookup = _build_homolog_lookup(result)
    seen = set()
    breaks = set()

    def visit(hid: int):
        if hid in seen:
            return
        seen.add(hid)

        h = homolog_lookup[hid]
        if h.chromosome != chromosome:
            return

        for seg in h.segments:
            breaks.add(seg.left)
            breaks.add(seg.right)
            if seg.parent_homolog_id is not None:
                visit(seg.parent_homolog_id)

    for hid in sample_homolog_ids:
        visit(hid)

    return sorted(breaks)


def _local_forests(result, chromosome: str, sample_homolog_ids):
    """
    Yield (left, right, edges) across all ancestry-defined intervals.
    """
    homolog_lookup = _build_homolog_lookup(result)

    chrom_lengths = {
        h.length for h in homolog_lookup.values() if h.chromosome == chromosome
    }
    if not chrom_lengths:
        raise ValueError(f"Unknown chromosome: {chromosome!r}")
    if len(chrom_lengths) != 1:
        raise ValueError(f"Inconsistent lengths found for chromosome: {chromosome!r}")
    chrom_length = chrom_lengths.pop()

    breaks = set(_ancestral_breakpoints(result, chromosome, sample_homolog_ids))
    breaks.add(0.0)
    breaks.add(chrom_length)
    breaks = sorted(breaks)

    out = []
    for left, right in zip(breaks[:-1], breaks[1:]):
        if right <= left:
            continue
        mid = 0.5 * (left + right)
        edges = _local_forest_at_position(result, chromosome, sample_homolog_ids, mid)
        out.append((left, right, edges))

    return out