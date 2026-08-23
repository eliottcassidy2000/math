#!/usr/bin/env python3
"""Exact safe-mass/component audit for the 58 conditional THM-3878 rows.

This is a new route, separate from the one-interval/deletion certificates.
It compares the exact pair-obstruction mass and component topology with:

* the THM-780 phase-cell floor, strengthened by the THM-789 symmetric
  difference-set move for an eleven-speed 1/12-deep pack;
* the best possible scalar Frechet coupling bound;
* inversion symmetry and the two reflected Lipschitz arcs; and
* component-slot counts at the legal boundary scale t=U.

All interval operations use Fraction arithmetic.  No canonical companion is
imported.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
import json
import sys


sys.stdout.reconfigure(newline="\n")

CHECKS = 0


def require(condition: bool, message: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(message)


Interval = tuple[Fraction, Fraction]


SCALE1 = (
    (1, 3), (1, 4), (1, 9), (1, 10),
    (2, 3), (2, 9), (2, 15), (2, 21), (2, 23),
    (3, 7), (3, 8), (3, 14), (3, 17), (3, 19), (3, 20),
    (3, 22), (3, 26), (3, 31), (3, 38),
    (4, 7), (4, 13), (4, 19), (4, 21), (4, 25), (4, 37),
    (4, 43), (4, 49), (4, 51),
    (5, 6), (5, 12), (5, 17), (5, 18), (5, 24), (5, 29),
    (5, 36), (5, 39), (5, 41), (5, 42), (5, 48), (5, 53),
    (5, 54), (5, 63),
    (6, 11), (6, 17), (6, 19), (6, 23), (6, 41), (6, 47),
    (6, 53), (6, 65),
    (7, 10), (7, 13), (7, 15), (7, 22),
    (8, 9), (8, 21), (9, 11),
)


def raw_danger(frequency: int) -> list[Interval]:
    """The open 1/14-danger comb, clipped to [0,1]."""
    radius = Fraction(1, 14 * frequency)
    answer = []
    for k in range(frequency + 1):
        centre = Fraction(k, frequency)
        left = max(Fraction(0), centre - radius)
        right = min(Fraction(1), centre + radius)
        if left < right:
            answer.append((left, right))
    return answer


def merge_open(intervals: list[Interval]) -> list[Interval]:
    """Union open intervals; equality of endpoints does not join them."""
    merged: list[list[Fraction]] = []
    for left, right in sorted(intervals):
        if not merged or left >= merged[-1][1]:
            merged.append([left, right])
        elif right > merged[-1][1]:
            merged[-1][1] = right
    return [(left, right) for left, right in merged]


def circle_component_lengths(intervals: list[Interval]) -> tuple[Fraction, ...]:
    merged = merge_open(intervals)
    require(bool(merged), "nonempty open union")
    lengths = [right - left for left, right in merged]
    if len(merged) >= 2 and merged[0][0] == 0 and merged[-1][1] == 1:
        wrap = merged[0][1] + 1 - merged[-1][0]
        lengths = [wrap] + [right - left for left, right in merged[1:-1]]
    return tuple(sorted(lengths))


def union_metrics(frequencies: tuple[int, ...]) -> tuple[Fraction, tuple[Fraction, ...]]:
    merged = merge_open(sum((raw_danger(w) for w in frequencies), []))
    measure = sum((right - left for left, right in merged), Fraction(0))
    return measure, circle_component_lengths(merged)


def complement_component_lengths(frequencies: tuple[int, ...]) -> tuple[Fraction, ...]:
    merged = merge_open(sum((raw_danger(w) for w in frequencies), []))
    gaps: list[Fraction] = []
    for (_, right), (next_left, _) in zip(merged, merged[1:]):
        # Equality leaves one safe wall point between two open danger teeth;
        # it is a genuine zero-length connected component of the closed safe
        # set and must be retained for component-count claims.
        gaps.append(next_left - right)
    wrap = merged[0][0] + 1 - merged[-1][1]
    # When wrap=0, the circle point 0 belongs to the danger comb and is not a
    # safe singleton.  Positive wrap is a genuine safe arc.
    if wrap > 0:
        gaps.append(wrap)
    return tuple(sorted(gaps))


def fmt(value: Fraction) -> str:
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    require(len(SCALE1) == 57 and len(set(SCALE1)) == 57, "57 scale-one rows")

    # THM-780 with d=11, beta=1/12, alpha=1/14 has M=84.  The
    # symmetric-difference upgrade from THM-789 doubles the heavy-cell mass.
    phase_cell_floor = Fraction(1, 84**11)
    symmetric_packet_floor = 2 * phase_cell_floor

    # Exact scalar hostile: the genuine eleven-speed AP pack.
    ap_measure, _ = union_metrics(tuple(range(1, 12)))
    ap_safe_measure = 1 - ap_measure
    ap_components = complement_component_lengths(tuple(range(1, 12)))
    require(len(ap_components) == 18, "AP11 closed safe component count")
    ap_positive_count = sum(length > 0 for length in ap_components)
    require(ap_positive_count == 14,
            "AP11 positive-length component count")
    require(max(ap_components) == Fraction(1, 77), "AP11 longest component")

    rows = []
    for p, q in SCALE1:
        obstruction_measure, components = union_metrics((p, q))
        beta = max(components)
        pair_intersection = Fraction(2, 7) - obstruction_measure

        # These rows survived exactly because one connected image arc of
        # length 1/42 can fit strictly into the obstruction.
        require(beta > Fraction(1, 42), f"survivor width {(p, q)}")

        # Any universal scalar safe-mass lower bound is <= the actual AP11
        # value.  Even that actual value is below every obstruction mass.
        require(ap_safe_measure < obstruction_measure,
                f"AP scalar hostile {(p, q)}")
        require(symmetric_packet_floor < obstruction_measure,
                f"phase-cell floor comparison {(p, q)}")

        # At the legal boundary t=U, choosing the AP hostile U=11 gives at
        # least as many positive-length obstruction slots as its 14 positive
        # safe arcs.  The four isolated closed-safe walls are a separate
        # pointwise endpoint channel and carry no mass.  This is only a
        # count-compatibility test, not an actual containment claim.
        require(11 * len(components) >= ap_positive_count,
                f"component slot compatibility {(p, q)}")

        rows.append({
            "s": 1,
            "p": p,
            "q": q,
            "measure": fmt(obstruction_measure),
            "components": len(components),
            "beta": fmt(beta),
            "pair_intersection": fmt(pair_intersection),
            "frechet_from_thm780": "0",
        })

    # The exact scale-two quotient obstruction from THM-3878.
    scale2_measure = Fraction(4, 63)
    scale2_components = 2
    scale2_beta = Fraction(2, 63)
    require(ap_safe_measure < scale2_measure, "AP scalar hostile scale two")
    require(symmetric_packet_floor < scale2_measure, "phase-cell scale two")
    require(11 * scale2_components >= ap_positive_count,
            "scale-two component slot compatibility")
    require(scale2_beta > Fraction(1, 42), "scale-two width survives")
    rows.append({
        "s": 2,
        "p": 1,
        "q": 9,
        "measure": fmt(scale2_measure),
        "components": scale2_components,
        "beta": fmt(scale2_beta),
        "pair_intersection": "quotient_not_raw_pair_moment",
        "frechet_from_thm780": "0",
    })

    scale1_measures = [Fraction(row["measure"]) for row in rows[:-1]]
    scale1_betas = [Fraction(row["beta"]) for row in rows[:-1]]
    scale1_counts = [row["components"] for row in rows[:-1]]
    min_measure = min(scale1_measures)
    min_measure_pairs = tuple(
        (row["p"], row["q"]) for row in rows[:-1]
        if Fraction(row["measure"]) == min_measure
    )
    min_beta = min(scale1_betas)
    min_beta_pairs = tuple(
        (row["p"], row["q"]) for row in rows[:-1]
        if Fraction(row["beta"]) == min_beta
    )

    row_blob = json.dumps(rows, sort_keys=True, separators=(",", ":")).encode("ascii")
    semantic = {
        "scope": "THM-3878 58 certificate survivors under t>=U",
        "safe_packet": "THM780 d11 beta1/12 alpha1/14 plus THM789 A-A",
        "conclusion": "zero closures from scalar measure, inversion symmetry, or component count",
        "missing": "mixed incidence/correlation of G(u) with pullback pair danger",
    }
    semantic_blob = json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")

    print("THM3878_SAFE_MASS_MOMENT_COMPONENT_AUDIT_20260823")
    print("scope=t>=U;scale1=57;scale2=(2,1,9);LRC14=OPEN")
    print(f"thm780_phase_cell_floor={fmt(phase_cell_floor)}")
    print(f"symmetric_return_packet_floor={fmt(symmetric_packet_floor)}")
    print(f"ap11_safe_measure={fmt(ap_safe_measure)};closed_components={len(ap_components)};positive_arcs={sum(x>0 for x in ap_components)};max_component={fmt(max(ap_components))}")
    print("thm1042_count_correction=table_20_false;canonical_script_positive_arcs_14;closed_topological_components_18")
    print(f"scale1_min_obstruction_measure={fmt(min_measure)};pairs={min_measure_pairs}")
    print(f"scale1_min_beta={fmt(min_beta)};pairs={min_beta_pairs}")
    print(f"scale1_component_count_range={min(scale1_counts)}..{max(scale1_counts)}")
    print(f"scale2_obstruction=measure:{fmt(scale2_measure)},components:2,beta:{fmt(scale2_beta)}")
    print("measure_frechet_closures=0;symmetry_closures=0;component_count_closures=0")
    print("moment_identity=mu_safe_after_pair=M0-M1+M2;mixed_M1_M2_not_supplied")
    print(f"row_census_sha256={sha256(row_blob).hexdigest()}")
    print(f"semantic_sha256={sha256(semantic_blob).hexdigest()}")
    print(f"checks={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
