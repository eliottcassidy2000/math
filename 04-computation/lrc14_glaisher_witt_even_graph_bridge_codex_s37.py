"""HYP-2660/T905: Glaisher/Witt/even-graph bridge for LRC14.

This scout turns the user's equinumerosity prompt into exact finite checks:

* Euler/Glaisher: distinct partitions and odd partitions have equal counts.
* Witt ghosts: log prod(1+x^n) has coefficients sigma_odd(m)/m.
* Tournament/even-graph: tournaments on n vertices equal even graphs on n+1
  vertices by parity closure at a root.
* LRC14: one-hole AP-window cores are encoded by Glaisher dyadic tower defects.

The goal is not to prove LRC(14), but to identify a proof quotient that keeps
the dyadic/CRT address exposed instead of collapsing immediately to raw speeds.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, product
from math import gcd


TARGET = Fraction(1, 14)
MAX_PARTITION_N = 40


def divisors(n: int) -> list[int]:
    return [d for d in range(1, n + 1) if n % d == 0]


def sigma_odd(n: int) -> int:
    return sum(d for d in divisors(n) if d % 2 == 1)


def distinct_partition_counts(limit: int) -> list[int]:
    counts = [0] * (limit + 1)
    counts[0] = 1
    for part in range(1, limit + 1):
        for total in range(limit, part - 1, -1):
            counts[total] += counts[total - part]
    return counts


def odd_partition_counts(limit: int) -> list[int]:
    counts = [0] * (limit + 1)
    counts[0] = 1
    for part in range(1, limit + 1, 2):
        for total in range(part, limit + 1):
            counts[total] += counts[total - part]
    return counts


def enumerate_odd_partitions(n: int, max_part: int | None = None) -> list[tuple[int, ...]]:
    if max_part is None or max_part > n:
        max_part = n if n % 2 else n - 1
    if n == 0:
        return [()]
    rows: list[tuple[int, ...]] = []
    for part in range(min(max_part, n), 0, -2):
        for rest in enumerate_odd_partitions(n - part, part):
            rows.append((part,) + rest)
    return rows


def glaisher_to_distinct(odd_parts: tuple[int, ...]) -> tuple[int, ...]:
    multiplicities = Counter(odd_parts)
    distinct: list[int] = []
    for odd, mult in sorted(multiplicities.items()):
        bit = 0
        while mult:
            if mult & 1:
                distinct.append(odd * (2**bit))
            mult >>= 1
            bit += 1
    return tuple(sorted(distinct, reverse=True))


def distinct_to_odd(distinct_parts: tuple[int, ...]) -> tuple[int, ...]:
    odd_parts: list[int] = []
    for part in distinct_parts:
        level, odd = v2_odd(part)
        odd_parts.extend([odd] * (2**level))
    return tuple(sorted(odd_parts, reverse=True))


def witt_log_coeff(m: int) -> Fraction:
    """Coefficient of x^m in log prod_{n>=1}(1+x^n)."""
    coeff = Fraction(0, 1)
    for n in divisors(m):
        k = m // n
        coeff += Fraction(1 if k % 2 else -1, k)
    return coeff


def v2_odd(n: int) -> tuple[int, int]:
    level = 0
    while n % 2 == 0:
        n //= 2
        level += 1
    return level, n


def frac_mod1(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def circle_dist(x: Fraction) -> Fraction:
    x = frac_mod1(x)
    return min(x, 1 - x)


def forbidden_intervals_for_speed(speed: int) -> list[tuple[Fraction, Fraction]]:
    intervals: list[tuple[Fraction, Fraction]] = []
    denom = 14 * speed
    for k in range(speed):
        start = Fraction(14 * k - 1, denom)
        end = Fraction(14 * k + 1, denom)
        if start < 0:
            intervals.append((Fraction(0), end))
            intervals.append((1 + start, Fraction(1)))
        elif end > 1:
            intervals.append((start, Fraction(1)))
            intervals.append((Fraction(0), end - 1))
        else:
            intervals.append((start, end))
    return intervals


def is_safe(t: Fraction, core: tuple[int, ...]) -> bool:
    return all(circle_dist(speed * t) > TARGET for speed in core)


def safe_components(core: tuple[int, ...]) -> list[tuple[Fraction, Fraction]]:
    points = {Fraction(0), Fraction(1)}
    for speed in core:
        for a, b in forbidden_intervals_for_speed(speed):
            points.add(a)
            points.add(b)
    ordered = sorted(points)
    components: list[tuple[Fraction, Fraction]] = []
    for a, b in zip(ordered, ordered[1:]):
        if a == b:
            continue
        mid = (a + b) / 2
        if is_safe(mid, core):
            components.append((a, b))
    return components


def safe_measure(core: tuple[int, ...]) -> Fraction:
    return sum((b - a for a, b in safe_components(core)), Fraction(0))


def tower_map(core: tuple[int, ...]) -> dict[int, list[int]]:
    towers: dict[int, list[int]] = defaultdict(list)
    for speed in sorted(core):
        level, odd = v2_odd(speed)
        towers[odd].append(level)
    return dict(sorted((odd, levels) for odd, levels in towers.items()))


def single_deletion_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    full = tuple(range(1, 14))
    for deleted in full:
        core = tuple(x for x in full if x != deleted)
        level, odd = v2_odd(deleted)
        towers = tower_map(core)
        rows.append(
            {
                "deleted": deleted,
                "measure": safe_measure(core),
                "level": level,
                "odd": odd,
                "odd_deleted": deleted % 2 == 1,
                "tower_after": tuple(towers.get(odd, [])),
                "lower_present": deleted % 2 == 0 and deleted // 2 in core,
                "upper_present": 2 * deleted in core,
                "core": core,
            }
        )
    return sorted(rows, key=lambda row: (row["measure"], row["deleted"]))


def tournament_to_even_graph(n: int, bits: tuple[int, ...]) -> set[tuple[int, int]]:
    """Map a tournament on n vertices to an even graph on vertices 0..n.

    Internal edge bit 1 records the anti-order/backward arc on {0,...,n-1}.
    Root edges incident to n are forced by parity.
    """
    root = n
    graph: set[tuple[int, int]] = set()
    degrees = [0] * (n + 1)
    idx = 0
    for i, j in combinations(range(n), 2):
        if bits[idx]:
            graph.add((i, j))
            degrees[i] ^= 1
            degrees[j] ^= 1
        idx += 1
    for i in range(n):
        if degrees[i]:
            graph.add((i, root))
            degrees[i] ^= 1
            degrees[root] ^= 1
    assert all(degree == 0 for degree in degrees)
    return graph


def even_graph_to_tournament_bits(n: int, graph: set[tuple[int, int]]) -> tuple[int, ...]:
    bits: list[int] = []
    for i, j in combinations(range(n), 2):
        bits.append(1 if (i, j) in graph or (j, i) in graph else 0)
    return tuple(bits)


def verify_tournament_even_graph(max_n: int = 6) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for n in range(1, max_n + 1):
        edge_count = n * (n - 1) // 2
        seen: set[tuple[tuple[int, int], ...]] = set()
        ok = True
        for bits in product((0, 1), repeat=edge_count):
            graph = tournament_to_even_graph(n, bits)
            if even_graph_to_tournament_bits(n, graph) != bits:
                ok = False
            seen.add(tuple(sorted(graph)))
        rows.append(
            {
                "n": n,
                "tournament_count": 2**edge_count,
                "even_graphs_on_n_plus_1": len(seen),
                "formula_even": 2 ** ((n + 1) * n // 2 - (n + 1) + 1),
                "switching_classes_tournaments_n": 2 ** ((n - 1) * (n - 2) // 2),
                "even_graphs_on_n": 2 ** (n * (n - 1) // 2 - n + 1) if n else 1,
                "invertible": ok,
            }
        )
    return rows


def transitive_tournament_fingerprint(vertices: list[str]) -> dict[str, object]:
    n = len(vertices)
    scores = {vertex: n - i - 1 for i, vertex in enumerate(vertices)}
    score_hist = Counter(scores.values())
    return {
        "vertices": vertices,
        "hamiltonian_path_count": 1,
        "score_histogram": dict(sorted(score_hist.items())),
        "directed_3_cycles": 0,
        "scc_sizes": [1] * n,
    }


def fmt_fraction(x: Fraction) -> str:
    return f"{x} = {float(x):.10f}"


def main() -> None:
    print("HYP-2660/T905 LRC14 Glaisher/Witt/even-graph bridge scout")
    print()

    distinct_counts = distinct_partition_counts(MAX_PARTITION_N)
    odd_counts = odd_partition_counts(MAX_PARTITION_N)
    mismatches = [
        (n, distinct_counts[n], odd_counts[n])
        for n in range(MAX_PARTITION_N + 1)
        if distinct_counts[n] != odd_counts[n]
    ]
    print("Euler/Glaisher coefficient check")
    print(f"  range: n<= {MAX_PARTITION_N}")
    print(f"  mismatches: {mismatches}")
    print("  selected counts p_distinct(n)=p_odd(n):")
    for n in [7, 14, 21, 28, 35, 40]:
        print(f"    n={n:2d}: {distinct_counts[n]}")
    print("  Glaisher bijection samples for n=14:")
    for odd_partition in enumerate_odd_partitions(14)[:8]:
        distinct = glaisher_to_distinct(odd_partition)
        inverse = distinct_to_odd(distinct)
        print(f"    {odd_partition!s:28s} -> {distinct!s:22s} -> {inverse}")
    print()

    print("Witt/log ghost check for f=prod(1+x^n)=g=prod_odd(1-x^odd)^-1")
    ghost_mismatches = []
    for m in range(1, MAX_PARTITION_N + 1):
        coeff = witt_log_coeff(m)
        ghost = m * coeff
        if ghost != sigma_odd(m):
            ghost_mismatches.append((m, ghost, sigma_odd(m)))
    print(f"  ghost mismatches m<= {MAX_PARTITION_N}: {ghost_mismatches}")
    print("  first ghost components m*[x^m] log f = sigma_odd(m):")
    for m in range(1, 15):
        print(f"    m={m:2d}: coeff={witt_log_coeff(m)!s:>7s}, ghost={sigma_odd(m):2d}")
    print()

    print("LRC14 single-deletion Glaisher tower ledger")
    odd_core = tuple(range(1, 14, 2))
    print(f"  odd skeleton core {odd_core}: safe measure {fmt_fraction(safe_measure(odd_core))}")
    print("  rows sorted by fixed-observer safe measure at target 1/14:")
    print(
        "    del  measure                 odd level  tower_after  lower upper  odd_deleted"
    )
    for row in single_deletion_rows():
        print(
            f"    {row['deleted']:>3d}  {str(row['measure']):>18s}  "
            f"{row['odd']:>3d} {row['level']:>5d}  "
            f"{str(row['tower_after']):>11s}  "
            f"{str(row['lower_present'])[0]:>5s} {str(row['upper_present'])[0]:>5s}  "
            f"{str(row['odd_deleted'])[0]:>11s}"
        )
    print(
        "  reading: every even deletion leaves the odd skeleton intact; drop 6 deletes "
        "the level-1 bit in odd tower 3 while both 3 and 12 survive."
    )
    print()

    print("Tournament/even-graph equinumerosity via parity closure")
    for row in verify_tournament_even_graph(6):
        print(
            f"  n={row['n']}: tournaments={row['tournament_count']}, "
            f"even_graphs_on_n+1={row['even_graphs_on_n_plus_1']}, "
            f"formula={row['formula_even']}, invertible={row['invertible']}, "
            f"switching_classes(T_n)=even_graphs_on_n={row['even_graphs_on_n']}"
        )
    print(
        "  bijection: internal anti-order bits choose a graph on n labeled vertices; "
        "root edges are forced uniquely so every degree is even."
    )
    print()

    print("Tournament Analysis")
    vertices = [
        "glaisher_tower_word",
        "crt_residue_cell",
        "endpoint_owner_ledger",
        "even_graph_parity_closure",
        "witt_ghost_divisor_sum",
        "laurent_wall_polynomial",
        "raw_sumset_excess",
        "raw_speed",
    ]
    fingerprint = transitive_tournament_fingerprint(vertices)
    print(f"  vertices: {', '.join(vertices)}")
    print(
        "  pairwise observable: which quotient preserves the LRC14 lower-bound "
        "predicate while distinguishing the one-hole AP-window rows."
    )
    print(
        "  switch/gauge: choose the parity/ghost closure before scalarizing; "
        "raw speed and raw excess are terminal quotients."
    )
    print("  tie Hamiltonian path:")
    print("    " + " > ".join(vertices))
    print(f"  score histogram: {fingerprint['score_histogram']}")
    print(f"  directed 3-cycles: {fingerprint['directed_3_cycles']}")
    print(f"  SCC sizes: {fingerprint['scc_sizes']}")
    print()

    print("Assumption challenge")
    print(
        "  alternate vertices considered: runners, gaps, speed towers, fixed circle "
        "sections, section boundaries, wall-crossing events, residues, cover arcs, "
        "Fourier/Witt ghost modes, matroid circuits, even-graph parity bits, and "
        "proof obligations."
    )
    print(
        "  preserved predicate: the fixed-observer safe set G_C and its addressed "
        "component/tower ownership. Destroyed information: exact component endpoints "
        "once one passes to partition counts, raw tournament counts, or scalar excess."
    )
    print(
        "  challenged assumption: the vertices should be runners or tournament arcs. "
        "The live quotient is closer to a parity-closed tower word: Glaisher closure "
        "for dyadic speeds, even-graph closure for tournaments, and Witt ghost closure "
        "for Euler products."
    )


if __name__ == "__main__":
    main()
