#!/usr/bin/env python3
"""
Event-media tournament probes for LRC.

codex-2026-06-01-S539

This script challenges the default assumption that tournament vertices should
be runners.  It studies tournaments whose vertices are media through which the
LRC clock changes:

  * holes: empty sectors of the n-sector circle;
  * sectors with hole-survival scores;
  * gates: boundaries between adjacent sectors, ranked by next crossing time.

The point is not to solve LRC.  It is to ask which tournament encodings keep
enough information for "observer loneliness" to be class-local, and which lose
the observer anchoring.

Tournament Analysis declaration
-------------------------------
Pairwise observables:
  sector hole survival time; boundary gate next-crossing priority.

Switch/gauge:
  bare isomorphism, observer-anchored coloring, and target-visible coloring.

Tie Hamiltonian path:
  cyclic sector/gate order.

Fingerprints:
  class counts and good-only/bad-only/mixed fiber counts over exact open
  sector-crossing cells.
"""

from __future__ import annotations

from collections import defaultdict
from fractions import Fraction
from itertools import combinations, permutations
from math import gcd
from functools import reduce


ZERO = Fraction(0)
ONE = Fraction(1)


def frac(x: Fraction) -> Fraction:
    return x - Fraction(x.numerator // x.denominator)


def primitive_speed_sets(n: int, max_speed: int):
    for speeds in combinations(range(1, max_speed + 1), n - 1):
        if reduce(gcd, speeds) == 1:
            yield speeds


def sector_walls(speeds: tuple[int, ...], n: int) -> list[Fraction]:
    walls = {ZERO, ONE}
    for v in speeds:
        for k in range(n * v):
            walls.add(Fraction(k, n * v))
    return sorted(walls)


def open_cell_midpoints(speeds: tuple[int, ...], n: int) -> list[Fraction]:
    walls = sector_walls(speeds, n)
    return [(a + b) / 2 for a, b in zip(walls, walls[1:]) if a < b]


def occupancy(speeds: tuple[int, ...], n: int, t: Fraction) -> tuple[int, ...]:
    counts = [0] * n
    for v in speeds:
        counts[int(n * frac(Fraction(v) * t)) % n] += 1
    return tuple(counts)


def is_good(counts: tuple[int, ...]) -> bool:
    return counts[0] == 0 and counts[-1] == 0


def next_gate_times(speeds: tuple[int, ...], n: int, t: Fraction) -> tuple[Fraction, ...]:
    """Return time until the next crossing of each boundary gate b/n."""
    out = []
    for gate in range(n):
        target = Fraction(gate, n)
        best = None
        for v in speeds:
            phase = frac(Fraction(v) * t)
            delta_phase = frac(target - phase)
            if delta_phase == ZERO:
                delta_phase = ONE
            delta = delta_phase / v
            if best is None or delta < best:
                best = delta
        out.append(best if best is not None else ONE)
    return tuple(out)


def sector_survival_times(
    speeds: tuple[int, ...], n: int, t: Fraction
) -> tuple[Fraction, ...]:
    """A sector changes when either adjacent boundary gate is crossed."""
    gate_times = next_gate_times(speeds, n, t)
    return tuple(min(gate_times[k], gate_times[(k + 1) % n]) for k in range(n))


def orient_by_scores(scores: list[tuple], n: int) -> tuple[tuple[int, ...], ...]:
    adj = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if scores[i] > scores[j]:
                adj[i][j], adj[j][i] = 1, 0
            elif scores[j] > scores[i]:
                adj[i][j], adj[j][i] = 0, 1
            else:
                # Cyclic tie path.  This is the built-in Hamiltonian tie path.
                forward = 1 <= (j - i) % n <= (n - 1) // 2
                adj[i][j], adj[j][i] = (1, 0) if forward else (0, 1)
    return tuple(tuple(row) for row in adj)


def canonical_colored(
    adj: tuple[tuple[int, ...], ...], colors: tuple
) -> tuple:
    n = len(adj)
    best = None
    for p in permutations(range(n)):
        c = tuple(colors[p[i]] for i in range(n))
        bits = tuple(adj[p[i]][p[j]] for i in range(n) for j in range(n) if i != j)
        key = (c, bits)
        if best is None or key < best:
            best = key
    return best


def hole_only_class(counts: tuple[int, ...], survival: tuple[Fraction, ...], anchored: bool):
    holes = [i for i, c in enumerate(counts) if c == 0]
    h = len(holes)
    if h == 0:
        return ("none",)
    scores = [(survival[i],) for i in holes]
    adj = orient_by_scores(scores, h)
    if anchored:
        # Distance to observer clasp; keeps "which hole is near 0" but not full position.
        colors = tuple(min(i, (len(counts) - i) % len(counts)) for i in holes)
    else:
        colors = tuple(0 for _ in holes)
    return (h, canonical_colored(adj, colors))


def sector_hole_survival_class(
    counts: tuple[int, ...], survival: tuple[Fraction, ...], anchored: bool
):
    n = len(counts)
    scores = [
        (1 if counts[k] == 0 else 0, survival[k], -counts[k])
        for k in range(n)
    ]
    adj = orient_by_scores(scores, n)
    if anchored:
        colors = tuple(
            (
                1 if k == 0 else 2 if k == n - 1 else 0,
                1 if counts[k] == 0 else 0,
            )
            for k in range(n)
        )
    else:
        colors = tuple(0 for _ in range(n))
    return canonical_colored(adj, colors)


def gate_priority_class(
    counts: tuple[int, ...], gate_times: tuple[Fraction, ...], anchored: bool
):
    n = len(counts)
    # Earlier crossing means higher priority: score is negative next time.
    scores = [(-gate_times[g],) for g in range(n)]
    adj = orient_by_scores(scores, n)
    if anchored:
        # Gate g is between sector g-1 and sector g.  Color by whether its two
        # adjacent cells are empty and whether the gate is one of the observer gates.
        colors = tuple(
            (
                1 if g == 0 else 2 if g == n - 1 else 0,
                1 if counts[(g - 1) % n] == 0 else 0,
                1 if counts[g] == 0 else 0,
            )
            for g in range(n)
        )
    else:
        colors = tuple(0 for _ in range(n))
    return canonical_colored(adj, colors)


MAPPINGS = [
    ("hole_only_bare", lambda c, s, g: hole_only_class(c, s, False)),
    ("hole_only_anchored", lambda c, s, g: hole_only_class(c, s, True)),
    ("sector_survival_bare", lambda c, s, g: sector_hole_survival_class(c, s, False)),
    ("sector_survival_anchored", lambda c, s, g: sector_hole_survival_class(c, s, True)),
    ("gate_priority_bare", lambda c, s, g: gate_priority_class(c, g, False)),
    ("gate_priority_anchored", lambda c, s, g: gate_priority_class(c, g, True)),
]


def classify(n: int, max_speed: int):
    fibers = {name: defaultdict(set) for name, _ in MAPPINGS}
    states = 0
    good_states = 0
    sets = 0
    certified = {name: set() for name, _ in MAPPINGS}

    for speeds in primitive_speed_sets(n, max_speed):
        sets += 1
        set_seen_good = {name: False for name, _ in MAPPINGS}
        for t in open_cell_midpoints(speeds, n):
            counts = occupancy(speeds, n, t)
            survival = sector_survival_times(speeds, n, t)
            gates = next_gate_times(speeds, n, t)
            good = is_good(counts)
            states += 1
            good_states += int(good)
            mark = "good" if good else "bad"
            for name, fn in MAPPINGS:
                cls = fn(counts, survival, gates)
                fibers[name][cls].add(mark)
                if good:
                    set_seen_good[name] = True
        for name, ok in set_seen_good.items():
            if ok:
                certified[name].add(speeds)

    rows = []
    for name, classmap in fibers.items():
        good_only = sum(1 for v in classmap.values() if v == {"good"})
        bad_only = sum(1 for v in classmap.values() if v == {"bad"})
        mixed = sum(1 for v in classmap.values() if v == {"good", "bad"})
        rows.append((name, len(classmap), good_only, bad_only, mixed, len(certified[name])))
    return sets, states, good_states, rows


def main() -> None:
    print("Event-media tournaments for LRC -- codex S539")
    print("=" * 78)
    print("Open cells are exact sector-crossing cells; boundary compactification is excluded.")
    print("Good means observer-adjacent sectors 0 and n-1 are both empty.")
    print()

    for n, max_speed in [(4, 10), (5, 8), (6, 7)]:
        sets, states, good_states, rows = classify(n, max_speed)
        print(f"n={n}: max_speed={max_speed}, speed_sets={sets}, open states={states}, good states={good_states}")
        for name, classes, good_only, bad_only, mixed, certified in rows:
            print(
                f"  {name:25s} classes={classes:5d} "
                f"good_only={good_only:4d} bad_only={bad_only:5d} "
                f"mixed={mixed:4d} certified_sets={certified:4d}/{sets}"
            )
        print()

    print("Interpretation")
    print("-" * 78)
    print("Bare hole and gate tournaments are beautifully small but usually mixed:")
    print("  they forget where the observer clasp is.  This is the main hidden assumption.")
    print("Anchoring the observer cells/gates restores purity at the cost of label tax.")
    print("Hole-only tournaments are the most radical: vertices are vacancies, not")
    print("  runners or sectors.  They suggest an exclusion-process version of LRC.")
    print("Gate-priority tournaments are kinetic data structures: edges change exactly")
    print("  when boundary-crossing certificates fail.  LRC becomes a question about")
    print("  whether the two observer gates can jointly lose incoming pressure.")


if __name__ == "__main__":
    main()
