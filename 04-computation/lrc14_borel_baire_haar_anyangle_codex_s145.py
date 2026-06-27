#!/usr/bin/env python3
"""S145: Borel/Baire/Haar and any-angle planning carriers for LRC14.

This scout imports the prompt's Borel sets, Baire sets, Haar theorem, and
any-angle path-planning algorithms into the current POKE/LRC14 proof route.

Mathematical reading:
  * For a finite runner set S, the strict safe set
        {t in T : ||v t|| > 1/14 for all v in S}
    is an open Borel/Baire subset of the circle.
  * The closed threshold-safe set with >= is closed Borel/Baire.
  * Haar's theorem supplies the canonical normalized translation-invariant
    measure on the compact circle, and on any finite-dimensional orbit-closure
    subtorus cut out by integer relations.
  * In compact metric tori, Baire sets and Borel sets agree for these finite
    arc events; the useful distinction is category/interior versus Haar mass.

POKE reading:
  AP and Goddyn-Wong are boundary-only rows at threshold 1/14: strict Haar
  mass zero, no Baire interior, but nonempty closed threshold support.  A proof
  route that sees only Haar mass loses the exact tight locus.  The sixth
  any-angle planner suggested here, "Haar-Baire Wave*", therefore propagates
  interval fronts labelled by (Haar mass, Baire interior, boundary support).
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from itertools import permutations
from math import gcd


DELTA = Fraction(1, 14)


def frac_part(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def circular_distance_to_integer(x: Fraction) -> Fraction:
    y = frac_part(x)
    return min(y, 1 - y)


def split_circle_interval(a: Fraction, b: Fraction) -> list[tuple[Fraction, Fraction]]:
    """Return pieces of [a,b] modulo 1 as intervals in [0,1].

    The interval length is always < 1 in this script because delta < 1/2.
    Boundaries are immaterial for Haar measure; we keep closed-style endpoints
    so exact rational union arithmetic remains simple.
    """
    while a < 0:
        a += 1
        b += 1
    while a >= 1:
        a -= 1
        b -= 1
    if b <= 1:
        return [(a, b)]
    return [(a, Fraction(1)), (Fraction(0), b - 1)]


def unsafe_intervals_for_speed(v: int, delta: Fraction = DELTA) -> list[tuple[Fraction, Fraction]]:
    radius = delta / v
    out: list[tuple[Fraction, Fraction]] = []
    for m in range(v):
        center = Fraction(m, v)
        out.extend(split_circle_interval(center - radius, center + radius))
    return out


def merge_intervals(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged: list[list[Fraction]] = [[intervals[0][0], intervals[0][1]]]
    for a, b in intervals[1:]:
        if a <= merged[-1][1]:
            if b > merged[-1][1]:
                merged[-1][1] = b
        else:
            merged.append([a, b])
    return [(a, b) for a, b in merged]


def complement_intervals(covered: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    if not covered:
        return [(Fraction(0), Fraction(1))]
    out: list[tuple[Fraction, Fraction]] = []
    cursor = Fraction(0)
    for a, b in covered:
        if cursor < a:
            out.append((cursor, a))
        cursor = max(cursor, b)
    if cursor < 1:
        out.append((cursor, Fraction(1)))
    return out


def safe_open_components(speeds: tuple[int, ...], delta: Fraction = DELTA) -> list[tuple[Fraction, Fraction]]:
    intervals: list[tuple[Fraction, Fraction]] = []
    for v in speeds:
        intervals.extend(unsafe_intervals_for_speed(v, delta))
    return complement_intervals(merge_intervals(intervals))


def interval_measure(intervals: list[tuple[Fraction, Fraction]]) -> Fraction:
    return sum((b - a for a, b in intervals), start=Fraction(0))


def threshold_safe_points(speeds: tuple[int, ...], delta: Fraction = DELTA) -> tuple[Fraction, ...]:
    candidates = {Fraction(0)}
    for v in speeds:
        for m in range(v):
            candidates.add(frac_part((Fraction(m) - delta) / v))
            candidates.add(frac_part((Fraction(m) + delta) / v))
    good = []
    for t in sorted(candidates):
        if all(circular_distance_to_integer(Fraction(v) * t) >= delta for v in speeds):
            good.append(t)
    return tuple(good)


def max_gap_at_time(speeds: tuple[int, ...], t: Fraction) -> Fraction:
    return min(circular_distance_to_integer(Fraction(v) * t) for v in speeds)


@dataclass(frozen=True)
class Row:
    name: str
    speeds: tuple[int, ...]
    note: str


ROWS = [
    Row("AP", tuple(range(1, 14)), "tight AP boundary row"),
    Row("GW 12->24", tuple(list(range(1, 12)) + [13, 24]), "tight Goddyn-Wong boundary row"),
    Row("near 12->36", tuple(list(range(1, 12)) + [13, 36]), "first unit-excess K33 near-miss"),
    Row("petal 10->20", tuple(list(range(1, 10)) + [11, 12, 13, 20]), "unit-visible C27 petal"),
    Row("petal 13->26", tuple(list(range(1, 13)) + [26]), "unit-visible C27 petal"),
    Row(
        "two-swap 10,12->20,24",
        tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 20, 24]),
        "S138 petal plus GW splice",
    ),
    Row(
        "two-swap 10,12->20,36",
        tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 13, 20, 36]),
        "S138 petal plus K33 splice",
    ),
]


@dataclass(frozen=True)
class Planner:
    name: str
    vector: tuple[int, int, int, int, int, int, int]
    note: str

    @property
    def total(self) -> int:
        return sum(self.vector)


PLANNERS = [
    Planner("Haar-Baire Wave*", (5, 5, 5, 4, 5, 5, 5), "new proof planner: Borel fronts with Haar/category/boundary labels"),
    Planner("ANYA interval-taut", (5, 5, 4, 1, 5, 2, 5), "interval nodes and taut-path pruning"),
    Planner("CWave primitive-front", (4, 4, 5, 2, 4, 2, 4), "circular/line wavefront primitives"),
    Planner("AP Theta* angle-prop", (3, 3, 4, 2, 4, 1, 3), "angle-propagated line-of-sight shortcut"),
    Planner("Theta* line-of-sight", (3, 2, 4, 2, 3, 1, 3), "parent-to-successor line-of-sight shortcut"),
    Planner("Field D* interpolation", (2, 3, 4, 5, 3, 1, 3), "dynamic interpolation over weighted cells"),
    Planner("Lazy Theta*", (2, 2, 4, 4, 3, 1, 3), "delayed visibility checks"),
    Planner("Block A*", (2, 4, 2, 2, 2, 1, 4), "local database of all small block paths"),
]


def tournament_fingerprint(planners: list[Planner]) -> dict[str, object]:
    names = [p.name for p in planners]
    n = len(planners)
    adj = [[False] * n for _ in range(n)]
    tie_order = {name: i for i, name in enumerate(names)}
    for i, left in enumerate(planners):
        for j, right in enumerate(planners):
            if i == j:
                continue
            wins = sum(a > b for a, b in zip(left.vector, right.vector))
            losses = sum(a < b for a, b in zip(left.vector, right.vector))
            if wins == losses:
                adj[i][j] = tie_order[left.name] < tie_order[right.name]
            else:
                adj[i][j] = wins > losses
    scores = [sum(row) for row in adj]

    c3 = 0
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(b + 1, n):
                cyc = (
                    adj[a][b] and adj[b][c] and adj[c][a]
                ) or (
                    adj[a][c] and adj[c][b] and adj[b][a]
                )
                if cyc:
                    c3 += 1

    # Hamiltonian path count, small n brute force.
    ham = 0
    for perm in permutations(range(n)):
        if all(adj[perm[i]][perm[i + 1]] for i in range(n - 1)):
            ham += 1

    # SCCs by reachability.
    reach = [[adj[i][j] or i == j for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                for j in range(n):
                    reach[i][j] = reach[i][j] or reach[k][j]
    seen = set()
    sccs: list[tuple[str, ...]] = []
    for i in range(n):
        if i in seen:
            continue
        comp = tuple(j for j in range(n) if reach[i][j] and reach[j][i])
        seen.update(comp)
        sccs.append(tuple(names[j] for j in comp))

    return {
        "scores": dict(zip(names, scores)),
        "score_histogram": {s: scores.count(s) for s in sorted(set(scores))},
        "directed_3cycles": c3,
        "sccs": sccs,
        "hamiltonian_paths": ham,
        "ranking": [name for _score, name in sorted(zip(scores, names), reverse=True)],
    }


def print_rows() -> None:
    print("[1] Borel/Baire/Haar audit at threshold delta=1/14")
    print(
        f"  {'row':28s} {'strict Haar':>13s} {'components':>10s} "
        f"{'closed pts':>10s} {'best boundary':>14s} {'category readout'}"
    )
    for row in ROWS:
        comps = safe_open_components(row.speeds)
        measure = interval_measure(comps)
        points = threshold_safe_points(row.speeds)
        best = max((max_gap_at_time(row.speeds, t) for t in points), default=Fraction(0))
        if measure > 0:
            cat = "nonempty open, nonmeager"
        elif points:
            cat = "boundary-only, meager"
        else:
            cat = "covered at threshold"
        print(
            f"  {row.name:28s} {str(measure):>13s} {len(comps):10d} "
            f"{len(points):10d} {str(best):>14s} {cat}"
        )
    print("  readout:")
    print("    AP and GW are not positive-Haar objects at threshold; they are")
    print("    closed boundary-support objects at the six denominator-14 unit times.")
    print("    Near-miss and petal rows have positive strict Haar mass, hence a")
    print("    Baire-open witness interval rather than only a boundary witness.")
    print()


def print_sixth_planner() -> None:
    print("[2] Sixth any-angle carrier: Haar-Baire Wave*")
    print("  Classical any-angle planners propagate geometry through grids:")
    print("    Field D* interpolates; Theta* shortcuts by line of sight;")
    print("    Block A* caches local paths; ANYA uses interval/taut nodes;")
    print("    CWave propagates primitive wavefront arcs and lines.")
    print("  Proposed proof planner:")
    print("    Haar-Baire Wave* propagates Borel interval fronts on the circle/subtorus,")
    print("    labelled by (strict Haar mass, Baire interior, closed boundary support).")
    print("  LRC use:")
    print("    line of sight  -> no unsafe arc blocks a witness interval;")
    print("    taut path      -> every heading change is a cover-arc boundary event;")
    print("    wavefront      -> exact denominator combs and wall crossings;")
    print("    Haar label     -> canonical invariant mass on the orbit-closure torus;")
    print("    Baire label    -> open/nonmeager interval versus meager boundary support.")
    print()


def print_planner_tournament() -> None:
    print("[3] Tournament Analysis on planner/proof carriers")
    criteria = (
        "exactness",
        "interval_nodes",
        "continuous_geometry",
        "dynamic_update",
        "lrc_transfer",
        "haar_baire_label",
        "anti_scalar_guard",
    )
    print("  Criteria:", ", ".join(criteria))
    for planner in PLANNERS:
        print(f"  {planner.name:24s} vector={planner.vector} total={planner.total:2d}  {planner.note}")
    fp = tournament_fingerprint(PLANNERS)
    print("  fingerprint:")
    print(f"    score_histogram={fp['score_histogram']}")
    print(f"    directed_3cycles={fp['directed_3cycles']}")
    print(f"    sccs={fp['sccs']}")
    print(f"    hamiltonian_paths={fp['hamiltonian_paths']}")
    print("    ranking=" + " > ".join(fp["ranking"]))
    print()


def print_proof_targets() -> None:
    print("[4] Proof targets generated")
    targets = [
        (
            "boundary-support lemma",
            "Classify threshold-safe but strict-Haar-zero rows by exact denominator "
            "boundary points; AP/GW should be the only q=14 survivors after current reductions.",
        ),
        (
            "Haar subtorus lift",
            "For relation-lattice residuals, push the unsafe arc cover to the orbit "
            "closure and integrate by normalized Haar measure there, not by an "
            "independent-block product measure.",
        ),
        (
            "Baire-open strictness",
            "A covering-strictness proof can aim to show every non-AP/GW reduced "
            "residual has a nonempty open safe interval; this is stronger than "
            "positive measure in finite arc arrangements but keeps category visible.",
        ),
        (
            "any-angle finite frontier",
            "Use ANYA/CWave style interval fronts as nodes: expand only wall-crossing "
            "intervals where a cover arc becomes tight, with lazy exact-M checks.",
        ),
        (
            "state-lift bridge",
            "If Haar-Baire Wave* cannot produce an open interval, the remaining "
            "boundary-only packet should carry enough labels to feed HYP-2908/THM-572.",
        ),
    ]
    for i, (name, body) in enumerate(targets, start=1):
        print(f"  {i}. {name}: {body}")
    print()


def main() -> None:
    print("=" * 78)
    print("S145 Borel/Baire/Haar and any-angle proof carriers for LRC14")
    print("=" * 78)
    print_rows()
    print_sixth_planner()
    print_planner_tournament()
    print_proof_targets()
    print("Status: proof-interface only; no LRC14 proof claimed.")


if __name__ == "__main__":
    main()
