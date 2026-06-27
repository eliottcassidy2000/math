#!/usr/bin/env python3
"""Poincare-style symmetry ledger for LRC14 proof packets.

This is a synthesis scout, not a proof.  It reads the LRC row as a collection
of worldlines x=v*t on the time/phase cylinder and asks which transformations
preserve the observer-relative loneliness predicate, which preserve only an
observer-coupled predicate, and which are useful only as residual/guardrail
language.

Tournament Analysis is over symmetry/proof carriers, not runners.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import gcd
from typing import Callable


DELTA = Fraction(1, 14)


def normalize_row(row: list[int]) -> list[int]:
    """Primitive absolute-speed normal form, preserving multiplicities."""
    vals = [abs(v) for v in row if v != 0]
    if not vals:
        return []
    g = 0
    for v in vals:
        g = gcd(g, v)
    return sorted(v // g for v in vals)


def danger_intervals_for_speed(v: int, delta: Fraction = DELTA) -> list[tuple[Fraction, Fraction]]:
    """Return danger intervals for ||v*t|| < delta on [0,1)."""
    v = abs(v)
    if v == 0:
        return [(Fraction(0), Fraction(1))]
    intervals: list[tuple[Fraction, Fraction]] = []
    for m in range(v):
        start = (Fraction(m, v) - Fraction(delta, v)) % 1
        end = (Fraction(m, v) + Fraction(delta, v)) % 1
        if start < end:
            intervals.append((start, end))
        else:
            intervals.append((Fraction(0), end))
            intervals.append((start, Fraction(1)))
    return intervals


def union_measure(intervals: list[tuple[Fraction, Fraction]]) -> Fraction:
    if not intervals:
        return Fraction(0)
    intervals = sorted(intervals)
    total = Fraction(0)
    cur_a, cur_b = intervals[0]
    for a, b in intervals[1:]:
        if a <= cur_b:
            if b > cur_b:
                cur_b = b
        else:
            total += cur_b - cur_a
            cur_a, cur_b = a, b
    total += cur_b - cur_a
    return total


def safe_measure(row: list[int], delta: Fraction = DELTA) -> Fraction:
    intervals: list[tuple[Fraction, Fraction]] = []
    for v in row:
        intervals.extend(danger_intervals_for_speed(v, delta))
    return Fraction(1) - union_measure(intervals)


def fmt_frac(x: Fraction) -> str:
    if x.denominator == 1:
        return str(x.numerator)
    return f"{x.numerator}/{x.denominator}"


@dataclass(frozen=True)
class RowTransform:
    name: str
    transform: Callable[[list[int]], list[int]]
    interpretation: str
    required_label: str


@dataclass(frozen=True)
class Carrier:
    name: str
    predicate: int
    observer_anchor: int
    integer_lattice: int
    exact_scale: int
    packet_labels: int
    boost_compat: int
    anti_scalar: int
    proof_maturity: int
    destroyed: str

    @property
    def score(self) -> tuple[int, ...]:
        return (
            self.predicate,
            self.observer_anchor,
            self.integer_lattice,
            self.exact_scale,
            self.packet_labels,
            self.boost_compat,
            self.anti_scalar,
            self.proof_maturity,
        )


def tournament_fingerprint(carriers: list[Carrier]) -> dict[str, object]:
    n = len(carriers)
    outdeg = [0] * n
    edges: dict[tuple[int, int], int] = {}
    for i in range(n):
        for j in range(i + 1, n):
            si, sj = carriers[i].score, carriers[j].score
            wi = sum(a > b for a, b in zip(si, sj))
            wj = sum(b > a for a, b in zip(si, sj))
            winner = i if wi >= wj else j
            loser = j if winner == i else i
            edges[(winner, loser)] = 1
            outdeg[winner] += 1

    c3 = 0
    for i in range(n):
        for j in range(i + 1, n):
            for k in range(j + 1, n):
                wins = {
                    (a, b)
                    for a, b in edges
                    if a in (i, j, k) and b in (i, j, k)
                }
                if (
                    ((i, j) in wins and (j, k) in wins and (k, i) in wins)
                    or ((i, k) in wins and (k, j) in wins and (j, i) in wins)
                ):
                    c3 += 1

    order = sorted(range(n), key=lambda idx: carriers[idx].score, reverse=True)
    return {
        "score_hist": sorted((d, outdeg.count(d)) for d in set(outdeg)),
        "directed_3cycles": c3,
        "hamiltonian_path": " > ".join(carriers[i].name for i in order),
    }


def main() -> None:
    named_rows = {
        "AP13": list(range(1, 14)),
        "GW13": list(range(1, 12)) + [13, 24],
        "K33_12_to_36": list(range(1, 12)) + [13, 36],
        "q27_two_block_probe": [7, 14, 21, 28, 35, 42, 49, 56, 63, 70, 77, 91, 260],
    }

    transforms = [
        RowTransform(
            "identity",
            lambda row: list(row),
            "baseline worldline packet",
            "none",
        ),
        RowTransform(
            "global_reflection_v_to_minus_v",
            lambda row: [-v for v in row],
            "space reflection / time reversal; exact because ||-v*t||=||v*t||",
            "orientation label optional for LRC, required for winding tournaments",
        ),
        RowTransform(
            "independent_runner_parity_flips",
            lambda row: [v if i % 2 == 0 else -v for i, v in enumerate(row)],
            "runner-by-runner counterclockwise move; exact for observer distance",
            "sign cut is invisible to LRC but visible to pairwise order",
        ),
        RowTransform(
            "time_dilation_speed_scale_5",
            lambda row: [5 * v for v in row],
            "integer time dilation; measure-preserving under t -> 5t",
            "gcd / primitive scale label",
        ),
        RowTransform(
            "stationary_velocity_translation_plus_5",
            lambda row: [v + 5 for v in row],
            "Galilean/Poincare boost after forgetting observer velocity",
            "observer speed label is destroyed",
        ),
        RowTransform(
            "observer_coupled_boost_plus_5_then_recenter",
            lambda row: list(row),
            "boost all runners and observer by +5, then subtract observer speed",
            "observer-coupled frame label retained",
        ),
    ]

    print("=== LRC14 Poincare-style symmetry ledger S172 ===")
    print(f"threshold delta={fmt_frac(DELTA)}")
    print()
    print("Worldline model")
    print("  runner i: x_i(t)=v_i*t on the time/phase cylinder")
    print("  danger tube: ||x_i(t)-x_observer(t)|| < delta")
    print("  safe time: a horizontal time slice misses every danger tube")
    print("  Poincare analogy: transformations are useful only after saying")
    print("    whether the observer worldline, integer lattice, and tube width are retained.")
    print()

    print("Exact safe-measure stress tests")
    for row_name, row in named_rows.items():
        base = safe_measure(row)
        print(f"  row={row_name:<20s} base_safe={fmt_frac(base):>10s}")
        for tr in transforms[1:]:
            moved = tr.transform(row)
            val = safe_measure(moved)
            same = "same" if val == base else "CHANGED"
            primitive = normalize_row(moved)
            preview = ",".join(map(str, primitive[:8]))
            if len(primitive) > 8:
                preview += ",..."
            print(
                f"    {tr.name:<42s} safe={fmt_frac(val):>10s} {same:<7s}"
                f" primitive=[{preview}] label={tr.required_label}"
            )
        print()

    carriers = [
        Carrier(
            "observer_coupled_worldline_tube_groupoid",
            5,
            5,
            4,
            5,
            5,
            5,
            5,
            4,
            "none if observer worldline and tube metric are retained",
        ),
        Carrier(
            "individual_sign_flip_parity_kernel",
            5,
            5,
            5,
            5,
            3,
            3,
            5,
            5,
            "pairwise winding orientation and tournament order",
        ),
        Carrier(
            "integer_time_dilation_primitive_scale",
            5,
            4,
            5,
            5,
            4,
            4,
            4,
            5,
            "raw magnitude before primitive/gcd normalization",
        ),
        Carrier(
            "anchored_metric_winding_tournament",
            4,
            5,
            4,
            4,
            5,
            3,
            4,
            4,
            "continuous tube widths unless threshold gaps are attached",
        ),
        Carrier(
            "stationary_velocity_translation",
            2,
            1,
            5,
            3,
            2,
            2,
            5,
            2,
            "observer velocity; not a standard LRC automorphism",
        ),
        Carrier(
            "bare_winding_iso_class",
            2,
            1,
            3,
            2,
            2,
            1,
            2,
            3,
            "metric gap widths, observer placement, exact Farey scale",
        ),
        Carrier(
            "lorentz_velocity_addition_shadow",
            1,
            2,
            1,
            1,
            2,
            4,
            4,
            1,
            "integer-speed lattice and fixed Euclidean circle metric",
        ),
        Carrier(
            "raw_speed_scalar",
            0,
            0,
            1,
            0,
            0,
            0,
            0,
            0,
            "everything packet-level",
        ),
    ]
    fp = tournament_fingerprint(carriers)
    print("Tournament Analysis over Poincare/symmetry proof carriers")
    print(f"  score_hist={fp['score_hist']}")
    print(f"  directed_3cycles={fp['directed_3cycles']}")
    print(f"  hamiltonian_path={fp['hamiltonian_path']}")
    print()

    print("Carrier guardrails")
    for c in sorted(carriers, key=lambda x: x.score, reverse=True):
        print(f"  {c.name:<42s} destroyed={c.destroyed}")
    print()

    print("Proof-target readout")
    print("  True automorphisms for the standard anchored LRC predicate:")
    print("    runner permutation, independent sign flips, global reflection, integer dilation/primitive scaling.")
    print("  Observer-coupled boosts are exact only in the enlarged worldline groupoid:")
    print("    keep observer velocity, then recenter by relative speeds.")
    print("  Stationary velocity translation is the useful failure mode:")
    print("    it measures exactly why observer labels cannot be quotiented away.")
    print("  Lorentz/Poincare language should enter as a cone/tube ledger, not a speed scalar:")
    print("    if a boost tilts the danger tube, the metric, observer, and lattice labels become load-bearing.")


if __name__ == "__main__":
    main()
