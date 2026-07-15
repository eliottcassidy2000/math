#!/usr/bin/env python3
"""Exact recursive comb closure for THM-810's quartic s=3 interface.

This file is intentionally self-contained.  It studies every packet

    3 * ([12] \ R)  union  {u_r : r in R},

where ``R`` is a multiplicative coset of ``{1,5,8,12}``,
``u_r == 3*r (mod 13)``, and ``u_r mod 3`` agrees on negative pairs.

The infinite search is reduced recursively by an interval-comb discrepancy
bound.  All interval endpoints and every comparison are exact Fractions.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations, permutations, product
from math import floor
from typing import Iterable, Sequence


P = 13
DELTA = F(1, P)
DUTY = 2 * DELTA
DISCREPANCY = DUTY * (1 - DUTY)
H = frozenset((1, 5, 8, 12))
BASE = tuple(range(1, P))


Interval = tuple[F, F]


def intersect_interval_unions(
    left: Sequence[Interval], right: Sequence[Interval]
) -> tuple[Interval, ...]:
    out: list[Interval] = []
    i = j = 0
    while i < len(left) and j < len(right):
        lo = max(left[i][0], right[j][0])
        hi = min(left[i][1], right[j][1])
        if lo < hi:
            out.append((lo, hi))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return tuple(out)


def safe_bands(speed: int) -> tuple[Interval, ...]:
    return tuple(
        (F(P * k + 1, P * speed), F(P * (k + 1) - 1, P * speed))
        for k in range(speed)
    )


def strict_safe_components(speeds: Iterable[int]) -> tuple[Interval, ...]:
    current: tuple[Interval, ...] = ((F(0), F(1)),)
    for speed in sorted(speeds):
        current = intersect_interval_unions(current, safe_bands(speed))
        if not current:
            break
    return current


def subtract_danger(components: Sequence[Interval], speed: int) -> tuple[Interval, ...]:
    """Return components intersected with ``||speed*t||>1/13``."""

    return intersect_interval_unions(components, safe_bands(speed))


def measure(components: Sequence[Interval]) -> F:
    return sum((hi - lo for lo, hi in components), F(0))


def cosets() -> tuple[tuple[int, ...], ...]:
    rows: list[tuple[int, ...]] = []
    for a in BASE:
        row = tuple(sorted({a * h % P for h in H}))
        if row not in rows:
            rows.append(row)
    assert len(rows) == 3
    return tuple(rows)


def negative_pairs(labels: Sequence[int]) -> tuple[tuple[int, int], ...]:
    unused = set(labels)
    pairs: list[tuple[int, int]] = []
    while unused:
        r = min(unused)
        pair = tuple(sorted((r, -r % P)))
        assert len(pair) == 2 and set(pair) <= unused
        pairs.append(pair)
        unused -= set(pair)
    assert len(pairs) == 2
    return tuple(pairs)


def crt_base(label: int, parity: int) -> int:
    values = [
        u
        for u in range(1, 3 * P)
        if u % P == 3 * label % P and u % 3 == parity
    ]
    assert len(values) == 1
    return values[0]


def least_remaining_bound(
    components: Sequence[Interval], remaining_count: int
) -> int:
    """Necessary upper bound for the least of ``remaining_count`` combs.

    For an interval union E with K components and measure L,

        |E intersect D_w| <= (2/13)|E| + (22/169)K/w.

    If m danger combs cover E and all their speeds are at least v, summing
    this inequality gives

        v <= 22*m*K / (13*(13-2*m)*|E|).
    """

    assert components and 1 <= remaining_count <= 4
    length = measure(components)
    bound = F(
        22 * remaining_count * len(components),
        13 * (13 - 2 * remaining_count),
    ) / length
    return floor(bound)


class Audit:
    def __init__(self) -> None:
        self.nodes = Counter()
        self.bound_max = Counter()
        self.bound_min = Counter()
        self.bound_sum = Counter()
        self.dead_no_candidate = Counter()
        self.tight_prefixes: list[tuple[object, ...]] = []
        self.terminal_rows: list[
            tuple[tuple[object, ...], tuple[tuple[int, int], ...], tuple[Interval, ...]]
        ] = []
        self.digest = sha256()

    def record_bound(self, depth: int, bound: int) -> None:
        self.bound_max[depth] = max(self.bound_max[depth], bound)
        if depth not in self.bound_min:
            self.bound_min[depth] = bound
        else:
            self.bound_min[depth] = min(self.bound_min[depth], bound)
        self.bound_sum[depth] += bound


def recurse(
    *,
    components: tuple[Interval, ...],
    remaining: tuple[int, ...],
    bases: dict[int, int],
    last_speed: int,
    chosen: tuple[tuple[int, int], ...],
    context: tuple[object, ...],
    audit: Audit,
) -> None:
    depth = len(chosen)
    audit.nodes[depth] += 1
    if not components:
        audit.tight_prefixes.append((*context, chosen))
        return
    if not remaining:
        audit.terminal_rows.append((context, chosen, components))
        return

    bound = least_remaining_bound(components, len(remaining))
    audit.record_bound(depth, bound)
    candidates: list[tuple[int, int]] = []
    for label in remaining:
        speed = bases[label]
        if speed <= last_speed:
            speed += 39 * ((last_speed - speed) // 39 + 1)
        while speed <= bound:
            candidates.append((speed, label))
            speed += 39
    candidates.sort()
    if not candidates:
        audit.dead_no_candidate[depth] += 1
        return

    for speed, label in candidates:
        new_components = subtract_danger(components, speed)
        record = (*context, chosen, label, speed, len(new_components), measure(new_components))
        audit.digest.update((repr(record) + "\n").encode())
        recurse(
            components=new_components,
            remaining=tuple(r for r in remaining if r != label),
            bases=bases,
            last_speed=speed,
            chosen=(*chosen, (label, speed)),
            context=context,
            audit=audit,
        )


def tournament_fingerprint(edge: Sequence[Sequence[bool]]) -> tuple[object, ...]:
    """Standard exact fingerprint for a four-vertex tournament."""

    n = len(edge)
    scores = tuple(sum(row) for row in edge)
    score_histogram = tuple(sorted(Counter(scores).items()))
    triangles = sum(
        (edge[i][j] and edge[j][k] and edge[k][i])
        or (edge[j][i] and edge[k][j] and edge[i][k])
        for i, j, k in combinations(range(n), 3)
    )
    reach = [[i == j or edge[i][j] for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            for j in range(n):
                reach[i][j] = reach[i][j] or (reach[i][k] and reach[k][j])
    unseen = set(range(n))
    scc: list[int] = []
    while unseen:
        i = min(unseen)
        component = {j for j in unseen if reach[i][j] and reach[j][i]}
        scc.append(len(component))
        unseen -= component
    hamiltonian_paths = sum(
        all(edge[path[i]][path[i + 1]] for i in range(n - 1))
        for path in permutations(range(n))
    )
    return score_histogram, triangles, tuple(sorted(scc, reverse=True)), hamiltonian_paths


def comparison_edge(values: Sequence[F]) -> tuple[tuple[bool, ...], ...]:
    """Orient larger value first, using index order as the tie Hamilton path."""

    n = len(values)
    edge = [[False] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        if values[i] >= values[j]:
            edge[i][j] = True
        else:
            edge[j][i] = True
    return tuple(tuple(row) for row in edge)


def terminal_tournament(
    context: tuple[object, ...], chosen: tuple[tuple[int, int], ...]
) -> tuple[object, ...]:
    """Raw-versus-conditional marginal-erosion tournament telemetry.

    Vertices are the four chosen exception combs, ordered by speed.  The raw
    pair observable compares their individual erosion of the original core
    safe set.  The switched gauge first removes the other two combs, then
    compares the two conditional marginal erosions.  Index order is the fixed
    tie Hamilton path in both gauges.
    """

    labels = context[0]
    assert isinstance(labels, tuple)
    core = tuple(3 * r for r in BASE if r not in labels)
    root = strict_safe_components(core)
    speeds = tuple(speed for _, speed in chosen)
    assert tuple(sorted(speeds)) == speeds and len(speeds) == 4

    root_measure = measure(root)
    raw_values = tuple(root_measure - measure(subtract_danger(root, w)) for w in speeds)
    raw_edge = comparison_edge(raw_values)

    switched = [[False] * 4 for _ in range(4)]
    for i, j in combinations(range(4), 2):
        conditioned = root
        for k, speed in enumerate(speeds):
            if k not in (i, j):
                conditioned = subtract_danger(conditioned, speed)
        conditioned_measure = measure(conditioned)
        value_i = conditioned_measure - measure(subtract_danger(conditioned, speeds[i]))
        value_j = conditioned_measure - measure(subtract_danger(conditioned, speeds[j]))
        if value_i >= value_j:
            switched[i][j] = True
        else:
            switched[j][i] = True
    switched_edge = tuple(tuple(row) for row in switched)
    flips = sum(
        raw_edge[i][j] != switched_edge[i][j] for i, j in combinations(range(4), 2)
    )
    return (
        tournament_fingerprint(raw_edge),
        tournament_fingerprint(switched_edge),
        flips,
    )


def main() -> None:
    audit = Audit()
    root_rows: list[tuple[object, ...]] = []
    for labels in cosets():
        core_labels = tuple(r for r in BASE if r not in labels)
        core = tuple(3 * r for r in core_labels)
        components = strict_safe_components(core)
        pairs = negative_pairs(labels)
        root_rows.append(
            (
                labels,
                core_labels,
                len(components),
                measure(components),
                least_remaining_bound(components, 4),
                pairs,
            )
        )
        for pair_bits in product((1, 2), repeat=2):
            parity = {
                r: pair_bits[index]
                for index, pair in enumerate(pairs)
                for r in pair
            }
            bases = {r: crt_base(r, parity[r]) for r in labels}
            context = (labels, pair_bits, tuple(sorted(bases.items())))
            recurse(
                components=components,
                remaining=labels,
                bases=bases,
                last_speed=0,
                chosen=(),
                context=context,
                audit=audit,
            )

    expected_roots = [
        ((1, 5, 8, 12), 36, F(8893, 45045), 246),
        ((2, 3, 10, 11), 36, F(772, 4095), 258),
        ((4, 6, 7, 9), 42, F(443, 2145), 275),
    ]
    assert [(row[0], row[2], row[3], row[4]) for row in root_rows] == expected_roots
    assert audit.nodes == Counter({0: 12, 1: 324, 2: 3046, 3: 4493, 4: 34})
    assert sum(audit.nodes.values()) == 7_909
    assert audit.bound_min == Counter({0: 246, 1: 130, 2: 71, 3: 32})
    assert audit.bound_max == Counter({0: 275, 1: 398, 2: 395, 3: 265})
    assert audit.dead_no_candidate == Counter({2: 1_062, 3: 4_459})
    assert not audit.tight_prefixes
    assert len(audit.terminal_rows) == 34

    terminal_min = min(
        audit.terminal_rows,
        key=lambda row: measure(row[2]),
    )
    terminal_min_measure = measure(terminal_min[2])
    terminal_min_witness = (terminal_min[0], terminal_min[1])
    assert terminal_min_measure == F(41_681, 720_720)
    tournament_rows = [
        terminal_tournament(context, chosen)
        for context, chosen, _ in audit.terminal_rows
    ]
    transitive_fingerprint = (
        ((0, 1), (1, 1), (2, 1), (3, 1)),
        0,
        (1, 1, 1, 1),
        1,
    )
    assert all(
        raw == switched == transitive_fingerprint
        for raw, switched, _ in tournament_rows
    )
    edge_flip_histogram = Counter(flips for _, _, flips in tournament_rows)
    assert edge_flip_histogram == Counter({0: 15, 1: 14, 2: 4, 3: 1})
    search_digest = audit.digest.hexdigest()
    assert search_digest == "c603151d5d6606f0492bb6c484dc976fb57befa258d0ee823b9f24e16d7f0a15"

    lines = [
        "quartic_s3 recursive comb closure",
        "threshold=1/13 duty=2/13 discrepancy=22/169",
    ]
    for row in root_rows:
        lines.append(
            f"root.labels={row[0]} core={row[1]} components={row[2]} "
            f"measure={row[3]} least_exception_bound={row[4]} negative_pairs={row[5]}"
        )
    lines.extend(
        [
            f"nodes_by_depth={dict(sorted(audit.nodes.items()))}",
            f"bound_min_by_depth={dict(sorted(audit.bound_min.items()))}",
            f"bound_max_by_depth={dict(sorted(audit.bound_max.items()))}",
            f"dead_no_candidate_by_depth={dict(sorted(audit.dead_no_candidate.items()))}",
            f"tight_prefixes={len(audit.tight_prefixes)}",
            f"terminal_uncovered_rows={len(audit.terminal_rows)}",
            f"terminal_min_residual_measure={terminal_min_measure}",
            f"terminal_min_witness={terminal_min_witness}",
            f"tournament.raw_fingerprint={transitive_fingerprint}",
            f"tournament.conditional_fingerprint={transitive_fingerprint}",
            f"tournament.edge_flip_histogram={dict(sorted(edge_flip_histogram.items()))}",
            f"search_digest={search_digest}",
        ]
    )
    if audit.tight_prefixes:
        lines.extend(f"tight={row}" for row in audit.tight_prefixes[:20])
    output_digest = sha256(("\n".join(lines) + "\n").encode()).hexdigest()
    lines.append(f"sha256={output_digest}")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
