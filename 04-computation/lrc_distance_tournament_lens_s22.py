#!/usr/bin/env python3
"""
lrc_distance_tournament_lens_s22.py

oracle-2026-05-31 S22

Probe a tournament-style metric refinement of the Lonely Runner Conjecture.

The usual LRC functional keeps only

    min_i ||v_i t||.

This script keeps the two circular neighbors of the stationary runner, the
nearest/second-nearest digraph on all runners, and the circular tournament
obtained from the runner order at selected witness times.  The goal is to make
the "two-neighbor" data explicit for the active n=14,15,16 frontier rows.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from itertools import combinations
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S356 = SourceFileLoader(
    "lonely_runner_residue_probe_s356",
    str(ROOT / "04-computation" / "lonely_runner_residue_probe_s356.py"),
).load_module()
S360 = SourceFileLoader(
    "lonely_runner_endpoint_protection_s360",
    str(ROOT / "04-computation" / "lonely_runner_endpoint_protection_s360.py"),
).load_module()

ONE = Fraction(1, 1)


@dataclass(frozen=True)
class Sample:
    label: str
    speeds: tuple[int, ...]


@dataclass(frozen=True)
class MetricTournamentRow:
    label: str
    n: int
    speeds: tuple[int, ...]
    classification: str
    threshold: Fraction
    time_kind: str
    time: Fraction
    arch_gap: Fraction
    debt: int
    left_gap: Fraction
    left_labels: tuple[int, ...]
    right_gap: Fraction
    right_labels: tuple[int, ...]
    stationary_nearest: Fraction
    stationary_nearest_labels: tuple[int, ...]
    stationary_second: Fraction
    stationary_second_labels: tuple[int, ...]
    stationary_in_nearest: int
    stationary_in_two_nearest: int
    mutual_nearest_pairs: int
    position_classes: int
    collision_pairs: int
    antipodal_pairs: int
    tie_pairs: int
    outdegree_stationary: int
    score_histogram: tuple[tuple[int, int], ...]
    cyclic_triples: int
    cyclic_triples_through_stationary: int


def fmt_frac(value: Fraction | None) -> str:
    return S356.fmt_frac(value)


def fmt_signed(value: Fraction) -> str:
    if value < 0:
        return "-" + fmt_frac(-value)
    return "+" + fmt_frac(value)


def fmt_labels(labels: tuple[int, ...], limit: int = 4) -> str:
    if not labels:
        return "-"
    shown = ",".join(str(label) for label in labels[:limit])
    if len(labels) > limit:
        shown += ",..."
    return shown


def fmt_hist(hist: tuple[tuple[int, int], ...]) -> str:
    return ",".join(f"{key}:{count}" for key, count in hist) or "-"


def circle(value: Fraction) -> Fraction:
    return value % ONE


def clockwise_delta(start: Fraction, end: Fraction) -> Fraction:
    return (end - start) % ONE


def circular_distance(a: Fraction, b: Fraction) -> Fraction:
    delta = clockwise_delta(a, b)
    return min(delta, ONE - delta)


def primitive_gcd(values: tuple[int, ...]) -> int:
    out = 0
    for value in values:
        out = gcd(out, value)
    return out


def ladder(n: int, scale: int, skip: int) -> tuple[int, ...]:
    speeds = tuple(sorted({1} | {scale * q for q in range(1, n) if q != skip}))
    if len(speeds) != n - 1 or primitive_gcd(speeds) != 1:
        raise ValueError(f"bad ladder n={n}, scale={scale}, skip={skip}: {speeds}")
    return speeds


def witness_time(speeds: tuple[int, ...]) -> tuple[str, Fraction]:
    report = S356.report("metric-tournament", list(speeds))
    if report.witness is not None:
        return "gap-mid", report.witness
    if report.boundary_witness is not None:
        return "boundary", report.boundary_witness
    raise ValueError(f"no LRC witness found for {speeds}")


def shell_labels(
    distances: tuple[tuple[Fraction, int], ...], shell_index: int
) -> tuple[Fraction, tuple[int, ...]]:
    shells: list[tuple[Fraction, list[int]]] = []
    for distance, label in distances:
        if not shells or shells[-1][0] != distance:
            shells.append((distance, [label]))
        else:
            shells[-1][1].append(label)
    if shell_index >= len(shells):
        return Fraction(0), ()
    distance, labels = shells[shell_index]
    return distance, tuple(sorted(labels))


def stationary_bracket(
    positions: tuple[Fraction, ...], labels: tuple[int, ...]
) -> tuple[Fraction, tuple[int, ...], Fraction, tuple[int, ...]]:
    right_deltas = tuple((clockwise_delta(Fraction(0), positions[i]), labels[i]) for i in range(1, len(labels)))
    left_deltas = tuple((clockwise_delta(positions[i], Fraction(0)), labels[i]) for i in range(1, len(labels)))

    right_gap = min(delta for delta, _label in right_deltas)
    left_gap = min(delta for delta, _label in left_deltas)
    right_labels = tuple(sorted(label for delta, label in right_deltas if delta == right_gap))
    left_labels = tuple(sorted(label for delta, label in left_deltas if delta == left_gap))
    return left_gap, left_labels, right_gap, right_labels


def lex_winner(i: int, j: int, positions: tuple[Fraction, ...]) -> int:
    delta = clockwise_delta(positions[i], positions[j])
    if delta == 0 or delta == Fraction(1, 2):
        return i if i < j else j
    return i if delta > Fraction(1, 2) else j


def tournament_stats(
    positions: tuple[Fraction, ...]
) -> tuple[int, tuple[tuple[int, int], ...], int, int]:
    outdegrees = [0] * len(positions)
    for i, j in combinations(range(len(positions)), 2):
        outdegrees[lex_winner(i, j, positions)] += 1

    cyclic = 0
    cyclic_through_stationary = 0
    for triple in combinations(range(len(positions)), 3):
        local = Counter()
        for i, j in combinations(triple, 2):
            local[lex_winner(i, j, positions)] += 1
        if sorted(local.values()) == [1, 1, 1]:
            cyclic += 1
            if 0 in triple:
                cyclic_through_stationary += 1

    return (
        outdegrees[0],
        tuple(sorted(Counter(outdegrees).items())),
        cyclic,
        cyclic_through_stationary,
    )


def neighbor_stats(
    positions: tuple[Fraction, ...], labels: tuple[int, ...]
) -> tuple[
    Fraction,
    tuple[int, ...],
    Fraction,
    tuple[int, ...],
    int,
    int,
    int,
]:
    ranked_by_vertex: list[tuple[tuple[Fraction, int, int], ...]] = []
    for i in range(len(labels)):
        ranked = tuple(
            sorted(
                (circular_distance(positions[i], positions[j]), labels[j], j)
                for j in range(len(labels))
                if i != j
            )
        )
        ranked_by_vertex.append(ranked)

    stationary_distances = tuple((distance, label) for distance, label, _idx in ranked_by_vertex[0])
    nearest_distance, nearest_labels = shell_labels(stationary_distances, 0)
    second_distance, second_labels = shell_labels(stationary_distances, 1)

    nearest_targets = tuple(ranked[0][2] for ranked in ranked_by_vertex)
    two_targets = tuple(tuple(item[2] for item in ranked[:2]) for ranked in ranked_by_vertex)

    stationary_in_nearest = sum(1 for idx in range(1, len(labels)) if nearest_targets[idx] == 0)
    stationary_in_two = sum(1 for idx in range(1, len(labels)) if 0 in two_targets[idx])
    mutual = sum(
        1
        for i, j in combinations(range(len(labels)), 2)
        if nearest_targets[i] == j and nearest_targets[j] == i
    )

    return (
        nearest_distance,
        nearest_labels,
        second_distance,
        second_labels,
        stationary_in_nearest,
        stationary_in_two,
        mutual,
    )


def degeneracy_stats(positions: tuple[Fraction, ...]) -> tuple[int, int, int, int]:
    multiplicities = Counter(positions)
    collision_pairs = sum(count * (count - 1) // 2 for count in multiplicities.values())
    antipodal_pairs = sum(
        1
        for i, j in combinations(range(len(positions)), 2)
        if circular_distance(positions[i], positions[j]) == Fraction(1, 2)
    )
    return len(multiplicities), collision_pairs, antipodal_pairs, collision_pairs + antipodal_pairs


def analyze(sample: Sample) -> MetricTournamentRow:
    summary = S360.summarize(list(sample.speeds))
    time_kind, time = witness_time(summary.speeds)
    labels = (0,) + summary.speeds
    positions = tuple(circle(Fraction(speed) * time) for speed in labels)
    (
        left_gap,
        left_labels,
        right_gap,
        right_labels,
    ) = stationary_bracket(positions, labels)
    (
        nearest_distance,
        nearest_labels,
        second_distance,
        second_labels,
        stationary_in_nearest,
        stationary_in_two,
        mutual,
    ) = neighbor_stats(positions, labels)
    position_classes, collision_pairs, antipodal_pairs, tie_pairs = degeneracy_stats(positions)
    out0, score_hist, cyclic, cyclic0 = tournament_stats(positions)

    return MetricTournamentRow(
        label=sample.label,
        n=len(summary.speeds) + 1,
        speeds=summary.speeds,
        classification=summary.classification,
        threshold=summary.threshold,
        time_kind=time_kind,
        time=time,
        arch_gap=summary.max_gap / summary.threshold,
        debt=summary.unprotected_count,
        left_gap=left_gap,
        left_labels=left_labels,
        right_gap=right_gap,
        right_labels=right_labels,
        stationary_nearest=nearest_distance,
        stationary_nearest_labels=nearest_labels,
        stationary_second=second_distance,
        stationary_second_labels=second_labels,
        stationary_in_nearest=stationary_in_nearest,
        stationary_in_two_nearest=stationary_in_two,
        mutual_nearest_pairs=mutual,
        position_classes=position_classes,
        collision_pairs=collision_pairs,
        antipodal_pairs=antipodal_pairs,
        tie_pairs=tie_pairs,
        outdegree_stationary=out0,
        score_histogram=score_hist,
        cyclic_triples=cyclic,
        cyclic_triples_through_stationary=cyclic0,
    )


def samples() -> tuple[Sample, ...]:
    return (
        Sample("n14 initial", tuple(range(1, 14))),
        Sample("n14 d=7 ladder", ladder(14, 7, 6)),
        Sample("n14 d=14 ladder", ladder(14, 14, 6)),
        Sample("n14 window completion", (1, 2, 3, 4, 5, 14, 15, 20, 23, 24, 25, 26, 99)),
        Sample("n15 initial", tuple(range(1, 15))),
        Sample("n15 d=3 ladder", ladder(15, 3, 6)),
        Sample("n15 d=5 ladder", ladder(15, 5, 6)),
        Sample("n15 d=15 ladder", ladder(15, 15, 6)),
        Sample("n15 window completion", (1, 2, 3, 4, 5, 6, 7, 13, 15, 16, 19, 22, 28, 180)),
        Sample("n16 initial", tuple(range(1, 16))),
        Sample("n16 d=2 ladder", ladder(16, 2, 14)),
        Sample("n16 d=4 ladder", ladder(16, 4, 14)),
        Sample("n16 d=8 ladder", ladder(16, 8, 14)),
        Sample("n16 d=16 ladder", ladder(16, 16, 14)),
        Sample("n16 window completion", (1, 2, 3, 4, 5, 6, 7, 8, 9, 16, 18, 22, 26, 30, 84)),
    )


def print_bracket_rows(rows: tuple[MetricTournamentRow, ...]) -> None:
    print("A. STATIONARY TWO-NEIGHBOR BRACKET")
    print("=" * 108)
    print(
        "  label                    class        t-kind   t          gap/th debt "
        "left        right       margin"
    )
    for row in rows:
        margin = min(row.left_gap, row.right_gap) - row.threshold
        left = f"{fmt_frac(row.left_gap)}[{fmt_labels(row.left_labels, 3)}]"
        right = f"{fmt_frac(row.right_gap)}[{fmt_labels(row.right_labels, 3)}]"
        print(
            f"  {row.label:<24} {row.classification:<12} {row.time_kind:<8} "
            f"{fmt_frac(row.time):<10} {fmt_frac(row.arch_gap):>6} "
            f"{row.debt:>4} {left:<11} {right:<11} {fmt_signed(margin):>7}"
        )
    print()


def print_neighbor_rows(rows: tuple[MetricTournamentRow, ...]) -> None:
    print("B. NEAREST AND SECOND-NEAREST DIGRAPH DATA")
    print("=" * 108)
    print(
        "  label                    NN0         2nd0        inNN0 in2NN0 mutual "
        "pos/ties/coll/anti"
    )
    for row in rows:
        nearest = f"{fmt_frac(row.stationary_nearest)}[{fmt_labels(row.stationary_nearest_labels, 3)}]"
        second = f"{fmt_frac(row.stationary_second)}[{fmt_labels(row.stationary_second_labels, 3)}]"
        degeneracy = (
            f"{row.position_classes}/{row.tie_pairs}/"
            f"{row.collision_pairs}/{row.antipodal_pairs}"
        )
        print(
            f"  {row.label:<24} {nearest:<11} {second:<11} "
            f"{row.stationary_in_nearest:>5} {row.stationary_in_two_nearest:>6} "
            f"{row.mutual_nearest_pairs:>6} {degeneracy}"
        )
    print()


def print_tournament_rows(rows: tuple[MetricTournamentRow, ...]) -> None:
    print("C. LEX-COMPLETED CIRCULAR TOURNAMENT SHADOW")
    print("=" * 108)
    print("  label                    out0 score_hist          cycles cycles_through_0")
    for row in rows:
        print(
            f"  {row.label:<24} {row.outdegree_stationary:>4} "
            f"{fmt_hist(row.score_histogram):<19} "
            f"{row.cyclic_triples:>6} {row.cyclic_triples_through_stationary:>16}"
        )
    print()


def print_synthesis(rows: tuple[MetricTournamentRow, ...]) -> None:
    safe_rows = [row for row in rows if min(row.left_gap, row.right_gap) >= row.threshold]
    ladder_margin_groups: dict[tuple[int, Fraction], list[str]] = {}
    for row in rows:
        if "ladder" not in row.label:
            continue
        margin = min(row.left_gap, row.right_gap) - row.threshold
        ladder_margin_groups.setdefault((row.n, margin), []).append(row.label)
    repeated_margins = [
        (n, margin, labels)
        for (n, margin), labels in sorted(ladder_margin_groups.items())
        if len(labels) > 1
    ]
    no_stationary_nearest = [
        row.label
        for row in rows
        if row.classification == "positive_gap" and row.stationary_in_nearest == 0
    ]
    print("D. SYNTHESIS")
    print("=" * 108)
    print(
        "  1. The LRC condition is exactly a two-sided bracket condition around "
        "the stationary vertex:"
    )
    print("     left_gap(t) >= 1/n and right_gap(t) >= 1/n.")
    print(
        f"     All {len(safe_rows)} sampled witness rows satisfy this exact local "
        "condition."
    )
    print(
        "  2. Quotient ladders preserve more local metric data than max-gap "
        "alone records.  Repeated ladder bracket margins:"
    )
    if repeated_margins:
        for n, margin, labels in repeated_margins:
            print(f"     n={n}, margin={fmt_signed(margin)}: {', '.join(labels)}")
    else:
        print("     -")
    print(
        "  3. The global nearest-neighbor digraph is extra data, not a substitute "
        "for the bracket."
    )
    print(
        "     Positive-gap rows where no moving runner has the stationary runner "
        "as its nearest neighbor:"
    )
    print(f"     {', '.join(no_stationary_nearest) if no_stationary_nearest else '-'}")
    print(
        "  4. Tournament structure enters through the circular order flow.  As t "
        "moves, adjacent swaps and half-turn ties change the lex-completed "
        "circular tournament; a counterexample would need every such order "
        "state to keep at least one of the two stationary bracket gaps below "
        "1/n."
    )
    print()


def main() -> None:
    rows = tuple(analyze(sample) for sample in samples())
    print("LRC DISTANCE-TOURNAMENT TWO-NEIGHBOR LENS (oracle-2026-05-31-S22)")
    print("=" * 108)
    print(
        "Rows include initial boundary systems, quotient ladders, and S21 "
        "coarse-row completions for n=14,15,16."
    )
    print(
        "pos/ties/coll/anti means number of position classes, lex tie pairs, "
        "collision pairs, and antipodal pairs."
    )
    print()
    print_bracket_rows(rows)
    print_neighbor_rows(rows)
    print_tournament_rows(rows)
    print_synthesis(rows)


if __name__ == "__main__":
    main()
