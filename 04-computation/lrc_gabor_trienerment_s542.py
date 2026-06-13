#!/usr/bin/env python3
"""Probe Gabor-cell t(r)ienerments for small LRC sector languages.

The experiment follows HYP-2027's suggestion that the right complex
vertices may be time-frequency cells rather than runners or arcs alone.
For a fixed denominator n and a sector occupancy vector c, use vertices

    (sector k, harmonic m), 0 <= k < n, 1 <= m < n.

The local Gabor atom is represented by the exact pair

    (amplitude c_k, phase residue m*k mod n).

Tournament Analysis declaration:
  * vertices: Gabor cells (sector, harmonic), not runners or arcs;
  * pairwise observable: the amplitude/phase key of c_k zeta_n^(mk);
  * switch/gauge: larger amplitude wins, equal amplitude is compared by
    cyclic phase majority, and equal or antipodal phase is a t(r)ienerment tie;
  * tie Hamiltonian path: lexicographic (phase, sector, harmonic), used as
    the canonical path inside tie classes;
  * fingerprints: tie counts/rates, ternary Krawtchouk B1 axis, score
    histograms, strict directed 3-cycles, SCC counts, and Hamiltonian path
    counts when the vertex count is at most 12.

This is intentionally a first quotient. It preserves the LRC target predicate
through the underlying sector occupancy vector (the good arc is empty at
sectors 0 and n-1), but it destroys runner identities and wall-crossing order.
The challenged assumption is that useful tournament vertices must be runners
or arcs; the probe tests whether Gabor cells expose a ternary uncertainty axis.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import cos, gcd, pi, sin
from statistics import mean
from typing import Iterable


CONFIGS = [(3, 10), (4, 10), (5, 8), (6, 7), (7, 7)]
HP_VERTEX_LIMIT = 12
EPS = 1.0e-9


@dataclass(frozen=True)
class Fingerprint:
    n: int
    vertex_count: int
    edge_count: int
    tie_count: int
    tie_rate_ppm: int
    ternary_b1: int
    score_histogram: tuple[tuple[int, int], ...]
    strict_3cycles: int
    tied_triples: int
    scc_count: int
    largest_scc: int
    hamiltonian_paths: int | None


@dataclass(frozen=True)
class VectorRecord:
    counts: tuple[int, ...]
    good: bool
    sector_support: int
    nonzero_harmonic_support: int
    full_fourier_support: int
    non_dc_gabor_product: int
    uncertainty_product: int
    fingerprint: Fingerprint


def primitive_speed_sets(n: int, max_speed: int) -> Iterable[tuple[int, ...]]:
    """Yield primitive speed sets of size n-1 in [1, max_speed]."""
    for speeds in combinations(range(1, max_speed + 1), n - 1):
        if gcd_many(speeds) == 1:
            yield speeds


def gcd_many(values: Iterable[int]) -> int:
    g = 0
    for value in values:
        g = gcd(g, value)
    return g


def open_cell_midpoints(n: int, speeds: tuple[int, ...]) -> list[Fraction]:
    """Return one midpoint from each open sector-cell of the runner arrangement."""
    walls = {Fraction(0, 1), Fraction(1, 1)}
    for speed in speeds:
        for sector in range(n):
            walls.add(Fraction(sector, n * speed))
    ordered = sorted(walls)
    midpoints: list[Fraction] = []
    for left, right in zip(ordered, ordered[1:]):
        if left < right:
            midpoints.append((left + right) / 2)
    return midpoints


def sector_occupancy(n: int, speeds: tuple[int, ...], time: Fraction) -> tuple[int, ...]:
    counts = [0] * n
    for speed in speeds:
        position = (speed * time) % 1
        sector = int(position * n)
        if sector == n:
            sector = 0
        counts[sector] += 1
    return tuple(counts)


def is_lrc_good(counts: tuple[int, ...]) -> bool:
    return counts[0] == 0 and counts[-1] == 0


def dft_harmonic_support(counts: tuple[int, ...]) -> int:
    n = len(counts)
    support = 0
    for harmonic in range(1, n):
        real = 0.0
        imag = 0.0
        for sector, count in enumerate(counts):
            angle = -2.0 * pi * harmonic * sector / n
            real += count * cos(angle)
            imag += count * sin(angle)
        if real * real + imag * imag > EPS:
            support += 1
    return support


def gabor_vertices(n: int) -> tuple[tuple[int, int], ...]:
    return tuple((sector, harmonic) for sector in range(n) for harmonic in range(1, n))


def atom_key(counts: tuple[int, ...], vertex: tuple[int, int]) -> tuple[int, int, int, int]:
    """Return the exact comparison key for a Gabor atom.

    The first two entries are the mathematical observable. The final entries
    provide a stable tie Hamiltonian path inside equal or antipodal phase
    classes.
    """
    sector, harmonic = vertex
    return counts[sector], (sector * harmonic) % len(counts), sector, harmonic


def compare_vertices(counts: tuple[int, ...], left: tuple[int, int], right: tuple[int, int]) -> int:
    """Return 1 if left beats right, -1 if right beats left, and 0 for a tie."""
    left_amp, left_phase, _, _ = atom_key(counts, left)
    right_amp, right_phase, _, _ = atom_key(counts, right)
    if left_amp != right_amp:
        return 1 if left_amp > right_amp else -1
    return cyclic_phase_compare(len(counts), left_phase, right_phase)


def cyclic_phase_compare(n: int, left_phase: int, right_phase: int) -> int:
    """Compare two phase residues by the cyclic tournament gauge.

    For odd n this is the regular cyclic tournament on residues. For even n,
    antipodal residues are tied, producing a genuine t(r)ienerment.
    """
    diff = (right_phase - left_phase) % n
    if diff == 0:
        return 0
    if n % 2 == 0 and diff == n // 2:
        return 0
    return 1 if diff < n / 2 else -1


def gabor_fingerprint(counts: tuple[int, ...]) -> Fingerprint:
    n = len(counts)
    vertices = gabor_vertices(n)
    vertex_count = len(vertices)
    edge_count = vertex_count * (vertex_count - 1) // 2
    adjacency = [0] * vertex_count
    strict_orientation: dict[tuple[int, int], int] = {}
    tie_count = 0

    for i, j in combinations(range(vertex_count), 2):
        cmp_value = compare_vertices(counts, vertices[i], vertices[j])
        if cmp_value > 0:
            adjacency[i] |= 1 << j
            strict_orientation[(i, j)] = 1
        elif cmp_value < 0:
            adjacency[j] |= 1 << i
            strict_orientation[(i, j)] = -1
        else:
            tie_count += 1
            adjacency[i] |= 1 << j
            adjacency[j] |= 1 << i
            strict_orientation[(i, j)] = 0

    outdegrees = [adjacency[i].bit_count() for i in range(vertex_count)]
    score_histogram = tuple(sorted(Counter(outdegrees).items()))
    strict_3cycles, tied_triples = triple_counts(vertex_count, strict_orientation)
    scc_sizes = strongly_connected_component_sizes(adjacency, vertex_count)
    hp_count = None
    if vertex_count <= HP_VERTEX_LIMIT:
        hp_count = hamiltonian_path_count(adjacency, vertex_count)

    tie_rate_ppm = round(1_000_000 * tie_count / edge_count) if edge_count else 0
    ternary_b1 = 2 * edge_count - 3 * tie_count
    return Fingerprint(
        n=n,
        vertex_count=vertex_count,
        edge_count=edge_count,
        tie_count=tie_count,
        tie_rate_ppm=tie_rate_ppm,
        ternary_b1=ternary_b1,
        score_histogram=score_histogram,
        strict_3cycles=strict_3cycles,
        tied_triples=tied_triples,
        scc_count=len(scc_sizes),
        largest_scc=max(scc_sizes) if scc_sizes else 0,
        hamiltonian_paths=hp_count,
    )


def triple_counts(vertex_count: int, strict_orientation: dict[tuple[int, int], int]) -> tuple[int, int]:
    strict_cycles = 0
    tied_triples = 0
    for a, b, c in combinations(range(vertex_count), 3):
        ab = strict_orientation[(a, b)]
        ac = strict_orientation[(a, c)]
        bc = strict_orientation[(b, c)]
        if ab == 0 or ac == 0 or bc == 0:
            tied_triples += 1
            continue
        # In the natural order a < b < c, the two cyclic orientations are:
        # a->b, b->c, c->a and b->a, c->b, a->c.
        if (ab, ac, bc) in ((1, -1, 1), (-1, 1, -1)):
            strict_cycles += 1
    return strict_cycles, tied_triples


def strongly_connected_component_sizes(adjacency: list[int], vertex_count: int) -> list[int]:
    reverse = [0] * vertex_count
    for source, mask in enumerate(adjacency):
        remaining = mask
        while remaining:
            low_bit = remaining & -remaining
            target = low_bit.bit_length() - 1
            reverse[target] |= 1 << source
            remaining ^= low_bit

    seen = 0
    order: list[int] = []

    def dfs_forward(start: int) -> None:
        nonlocal seen
        stack = [(start, False)]
        while stack:
            node, expanded = stack.pop()
            if expanded:
                order.append(node)
                continue
            if seen & (1 << node):
                continue
            seen |= 1 << node
            stack.append((node, True))
            remaining = adjacency[node] & ~seen
            while remaining:
                low_bit = remaining & -remaining
                stack.append((low_bit.bit_length() - 1, False))
                remaining ^= low_bit

    for node in range(vertex_count):
        if not (seen & (1 << node)):
            dfs_forward(node)

    assigned = 0
    sizes: list[int] = []
    for start in reversed(order):
        if assigned & (1 << start):
            continue
        size = 0
        stack = [start]
        assigned |= 1 << start
        while stack:
            node = stack.pop()
            size += 1
            remaining = reverse[node] & ~assigned
            while remaining:
                low_bit = remaining & -remaining
                target = low_bit.bit_length() - 1
                assigned |= low_bit
                stack.append(target)
                remaining ^= low_bit
        sizes.append(size)
    return sorted(sizes, reverse=True)


def hamiltonian_path_count(adjacency: list[int], vertex_count: int) -> int:
    full_mask = (1 << vertex_count) - 1
    dp = [[0] * vertex_count for _ in range(1 << vertex_count)]
    for vertex in range(vertex_count):
        dp[1 << vertex][vertex] = 1
    for mask in range(1 << vertex_count):
        for end in range(vertex_count):
            count = dp[mask][end]
            if not count:
                continue
            available = adjacency[end] & ~mask
            while available:
                low_bit = available & -available
                nxt = low_bit.bit_length() - 1
                dp[mask | low_bit][nxt] += count
                available ^= low_bit
    return sum(dp[full_mask])


def vector_record(counts: tuple[int, ...], cache: dict[tuple[int, ...], Fingerprint]) -> VectorRecord:
    if counts not in cache:
        cache[counts] = gabor_fingerprint(counts)
    sector_support = sum(1 for count in counts if count)
    nonzero_harmonic_support = dft_harmonic_support(counts)
    full_fourier_support = nonzero_harmonic_support + 1
    return VectorRecord(
        counts=counts,
        good=is_lrc_good(counts),
        sector_support=sector_support,
        nonzero_harmonic_support=nonzero_harmonic_support,
        full_fourier_support=full_fourier_support,
        non_dc_gabor_product=sector_support * nonzero_harmonic_support,
        uncertainty_product=sector_support * full_fourier_support,
        fingerprint=cache[counts],
    )


def scan_config(n: int, max_speed: int) -> str:
    fingerprint_cache: dict[tuple[int, ...], Fingerprint] = {}
    vector_records: dict[tuple[int, ...], VectorRecord] = {}
    cell_counter: Counter[tuple[int, ...]] = Counter()
    speed_set_count = 0
    cell_count = 0

    for speeds in primitive_speed_sets(n, max_speed):
        speed_set_count += 1
        for midpoint in open_cell_midpoints(n, speeds):
            counts = sector_occupancy(n, speeds, midpoint)
            cell_counter[counts] += 1
            cell_count += 1
            if counts not in vector_records:
                vector_records[counts] = vector_record(counts, fingerprint_cache)

    records = list(vector_records.values())
    good_records = [record for record in records if record.good]
    bad_records = [record for record in records if not record.good]
    good_cells = sum(cell_counter[record.counts] for record in good_records)
    bad_cells = cell_count - good_cells
    all_uncertainty_ok = all(record.uncertainty_product >= n for record in records)

    lines = [
        f"=== n={n}, max_speed={max_speed} ===",
        f"primitive speed sets: {speed_set_count}",
        f"open cells scanned: {cell_count}",
        f"distinct sector vectors: {len(records)}",
        f"LRC-good cells/vectors: {good_cells}/{len(good_records)}",
        f"LRC-bad cells/vectors: {bad_cells}/{len(bad_records)}",
        f"full Fourier support uncertainty S*(F_nonzero+1) >= n: {all_uncertainty_ok}",
    ]
    lines.extend(summary_block("good vectors", good_records, cell_counter))
    lines.extend(summary_block("bad vectors", bad_records, cell_counter))
    lines.extend(fingerprint_separation(records))
    lines.extend(example_block("minimum uncertainty good", good_records, cell_counter, prefer_min=True))
    lines.extend(example_block("maximum tie-rate bad", bad_records, cell_counter, prefer_min=False))
    return "\n".join(lines)


def summary_block(label: str, records: list[VectorRecord], cell_counter: Counter[tuple[int, ...]]) -> list[str]:
    if not records:
        return [f"{label}: none"]
    weighted_count = sum(cell_counter[record.counts] for record in records)
    return [
        f"{label}:",
        f"  cell mass: {weighted_count}",
        f"  sector support avg/range: {mean_field(records, 'sector_support')}",
        f"  nonzero harmonic support avg/range: {mean_field(records, 'nonzero_harmonic_support')}",
        f"  full Fourier support avg/range: {mean_field(records, 'full_fourier_support')}",
        f"  non-DC Gabor S*F avg/range: {mean_field(records, 'non_dc_gabor_product')}",
        f"  full uncertainty S*F avg/range: {mean_field(records, 'uncertainty_product')}",
        f"  tie count avg/range: {mean_fp(records, 'tie_count')}",
        f"  tie rate ppm avg/range: {mean_fp(records, 'tie_rate_ppm')}",
        f"  ternary B1 avg/range: {mean_fp(records, 'ternary_b1')}",
        f"  strict 3-cycles avg/range: {mean_fp(records, 'strict_3cycles')}",
        f"  SCC count avg/range: {mean_fp(records, 'scc_count')}",
        f"  HP counts distribution: {hp_distribution(records)}",
    ]


def mean_field(records: list[VectorRecord], field: str) -> str:
    values = [getattr(record, field) for record in records]
    return format_avg_range(values)


def mean_fp(records: list[VectorRecord], field: str) -> str:
    values = [getattr(record.fingerprint, field) for record in records]
    return format_avg_range(values)


def format_avg_range(values: list[int]) -> str:
    return f"{mean(values):.3f} / {min(values)}..{max(values)}"


def hp_distribution(records: list[VectorRecord]) -> str:
    values = [record.fingerprint.hamiltonian_paths for record in records]
    if any(value is None for value in values):
        return "not computed for this vertex count"
    return str(dict(sorted(Counter(values).items())))


def fingerprint_key(record: VectorRecord) -> tuple[object, ...]:
    fp = record.fingerprint
    return (
        record.sector_support,
        record.nonzero_harmonic_support,
        record.full_fourier_support,
        record.non_dc_gabor_product,
        record.uncertainty_product,
        fp.tie_count,
        fp.ternary_b1,
        fp.strict_3cycles,
        fp.tied_triples,
        fp.scc_count,
        fp.largest_scc,
        fp.score_histogram,
        fp.hamiltonian_paths,
    )


def fingerprint_separation(records: list[VectorRecord]) -> list[str]:
    classes: dict[tuple[object, ...], list[VectorRecord]] = defaultdict(list)
    for record in records:
        classes[fingerprint_key(record)].append(record)
    pure_good = 0
    pure_bad = 0
    mixed = 0
    mixed_examples: list[list[VectorRecord]] = []
    for members in classes.values():
        outcomes = {member.good for member in members}
        if outcomes == {True}:
            pure_good += 1
        elif outcomes == {False}:
            pure_bad += 1
        else:
            mixed += 1
            if len(mixed_examples) < 2:
                mixed_examples.append(members)
    lines = [
        "Gabor-trienerment fingerprint separation:",
        f"  classes: {len(classes)}",
        f"  pure good / pure bad / mixed: {pure_good} / {pure_bad} / {mixed}",
    ]
    for index, members in enumerate(mixed_examples, 1):
        trimmed = sorted(member.counts for member in members)[:4]
        lines.append(f"  mixed example {index}: {trimmed}")
    return lines


def example_block(
    label: str,
    records: list[VectorRecord],
    cell_counter: Counter[tuple[int, ...]],
    prefer_min: bool,
) -> list[str]:
    if not records:
        return []
    if prefer_min:
        key = lambda record: (
            record.uncertainty_product,
            record.fingerprint.tie_rate_ppm,
            record.counts,
        )
    else:
        key = lambda record: (
            -record.fingerprint.tie_rate_ppm,
            record.uncertainty_product,
            record.counts,
        )
    record = sorted(records, key=key)[0]
    fp = record.fingerprint
    return [
        f"{label} example:",
        f"  counts={record.counts}, cell_mass={cell_counter[record.counts]}",
        f"  S={record.sector_support}, F_nonzero={record.nonzero_harmonic_support}, "
        f"F_full={record.full_fourier_support}, S*F_full={record.uncertainty_product}",
        f"  ties={fp.tie_count}/{fp.edge_count}, B1={fp.ternary_b1}, strict_3cycles={fp.strict_3cycles}",
        f"  SCCs={fp.scc_count}, largest_SCC={fp.largest_scc}, HP={fp.hamiltonian_paths}",
    ]


def main() -> None:
    print("LRC Gabor t(r)ienerment scan, S542")
    print("Gauge: larger sector amplitude wins; equal amplitude uses cyclic phase; equal/antipodal phase ties.")
    print("Target predicate: sectors 0 and n-1 are empty.")
    print()
    for config in CONFIGS:
        print(scan_config(*config))
        print()


if __name__ == "__main__":
    main()
