from .ancestry import Homolog

def get_segments_extent(
    h: Homolog,
    start_loc: float,
    stop_loc: float,
    update_parent_id: bool = True,
) -> list[Segment]:
    """Return copies of segments from homolog h overlapping [start_loc, stop_loc)."""
    if start_loc >= stop_loc:
        raise ValueError("start_loc must be less than stop_loc")

    segs = []
    for seg in h.segments:
        if seg.left < stop_loc and start_loc < seg.right:
            trimmed_seg = seg.copy()
            trimmed_seg.left = max(trimmed_seg.left, start_loc)
            trimmed_seg.right = min(trimmed_seg.right, stop_loc)
            if update_parent_id:
                trimmed_seg.parent_homolog_id = h.homolog_id
            segs.append(trimmed_seg)

    return segs


def breakpoints_to_intervals(breakpoints: list[float], length: float) -> list[tuple[float, float]]:
    """Convert breakpoints into half-open intervals [left, right)."""
    bps = sorted(set(bp for bp in breakpoints if 0 < bp < length))
    edges = [0.0] + bps + [length]
    return [(edges[i], edges[i + 1]) for i in range(len(edges) - 1)]


def merge_adjacent_segments(
    segments: list[Segment],
) -> list[Segment]:
    if not segments:
        return []

    merged = [segments[0].copy()]

    for seg in segments[1:]:
        last = merged[-1]
        if (
            abs(last.right - seg.left) < 1e-12
            and last.parent_homolog_id == seg.parent_homolog_id
            and last.founder_homolog_id == seg.founder_homolog_id
        ):
            last.right = seg.right
        else:
            merged.append(seg.copy())

    return merged


def recombine_two_homologs(
    h0: Homolog,
    h1: Homolog,
    breakpoints: list[float],
    start_phase: int = 0,
    merge_adjacent: bool = True,
    homolog_id: Optional[int] = None,
    individual_id: Optional[str] = None,
    time: Optional[int] = None,
) -> Homolog:
    """Deterministically recombine two homologs and return a new recombinant Homolog."""
    if h0.length != h1.length:
        raise ValueError("homologs should be same length")
    if h0.chromosome != h1.chromosome:
        raise ValueError("homologs should be from the same chromosome")
    if start_phase not in (0, 1):
        raise ValueError("start_phase must be 0 or 1")

    intervals = breakpoints_to_intervals(breakpoints, h0.length)
    new_segments = []
    phase = start_phase

    for start, stop in intervals:
        if phase == 0:
            new_segments.extend(get_segments_extent(h0, start, stop))
        else:
            new_segments.extend(get_segments_extent(h1, start, stop))
        phase = 1 - phase

    if merge_adjacent:
        new_segments = merge_adjacent_segments(new_segments)

    return Homolog(
        homolog_id=homolog_id,
        chromosome=h0.chromosome,
        individual_id=individual_id,
        length=h0.length,
        time=time,
        segments=new_segments,
    )