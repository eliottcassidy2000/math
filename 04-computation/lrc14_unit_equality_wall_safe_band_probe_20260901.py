#!/usr/bin/env python3
"""Exact audit of the unit equality-wall safe-band lemma for d <= 7.

For a row d*C union T with |T|=d and every tail a unit modulo d, the
ordinary sheet-capacity sum is exactly one.  The boundary repair is to move a
lower-dimensional core phase into the closed quotient band on which one
selected tail is safe on *all* d lifts.  The full band, rather than one
endpoint lattice, gives the cofinal condition

    max(T) >= (14-d) * max(C).

All arithmetic below is rational.  This is a certificate audit, not a proof
of LRC(14) outside the displayed typed rows.
"""

from __future__ import annotations

import ast
from fractions import Fraction
from itertools import combinations
from math import gcd
from pathlib import Path


DELTA = Fraction(1, 14)
CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"CHECK_FAILED {label}")


def circle_norm(x: Fraction) -> Fraction:
    x %= 1
    return min(x, 1 - x)


def fraction_universe(max_denominator: int) -> tuple[Fraction, ...]:
    return tuple(
        sorted(
            {
                Fraction(numerator, denominator)
                for denominator in range(1, max_denominator + 1)
                for numerator in range(denominator)
            }
        )
    )


def safe_band(d: int) -> tuple[Fraction, Fraction]:
    return Fraction(d, 14), 1 - Fraction(d, 14)


def selected_lift_clearances(d: int, quotient_residue: Fraction) -> tuple[Fraction, ...]:
    """Clearances of a d-unit tail, after permuting its d lift labels."""
    z = quotient_residue % 1
    return tuple(circle_norm((z + label) / d) for label in range(d))


def bad_labels(d: int, speed: int, quotient_phase: Fraction) -> tuple[int, ...]:
    return tuple(
        label
        for label in range(d)
        if circle_norm(Fraction(speed) * (quotient_phase + label) / d) < DELTA
    )


def nearest_selected_safe_phase(
    quotient_phase: Fraction, d: int, selected_tail: int
) -> Fraction:
    """Nearest point whose selected-tail residue lies in the closed safe band."""
    y = quotient_phase % 1
    lo, hi = safe_band(d)
    z = (selected_tail * y) % 1
    if lo <= z <= hi:
        return y

    endpoints = (
        (Fraction(n) + endpoint) / selected_tail
        for n in range(selected_tail)
        for endpoint in (lo, hi)
    )
    return min(endpoints, key=lambda x: (circle_norm(x - y), x))


def affine_line_on_cell(speed: int, midpoint: Fraction) -> tuple[int, int]:
    value = speed * midpoint
    integer_part = value.numerator // value.denominator
    fractional_part = value - integer_part
    if fractional_part < Fraction(1, 2):
        return speed, -integer_part
    return -speed, integer_part + 1


def exact_core_maximin(core: tuple[int, ...]) -> tuple[Fraction, Fraction]:
    """Return (maximum clearance, a maximizing phase) by exact wall atomization."""
    walls = {Fraction(0), Fraction(1)}
    for speed in core:
        walls.update(Fraction(k, 2 * speed) for k in range(2 * speed + 1))
    ordered_walls = sorted(walls)
    candidates = set(ordered_walls)

    for left, right in zip(ordered_walls, ordered_walls[1:]):
        midpoint = (left + right) / 2
        lines = [affine_line_on_cell(speed, midpoint) for speed in core]
        for i, (slope_a, intercept_a) in enumerate(lines):
            for slope_b, intercept_b in lines[i + 1 :]:
                if slope_a == slope_b:
                    continue
                crossing = Fraction(intercept_b - intercept_a, slope_a - slope_b)
                if left < crossing < right:
                    candidates.add(crossing)

    return max(
        (min(circle_norm(speed * phase) for speed in core), phase)
        for phase in candidates
    )


def least_unit_at_least(lower_bound: int, d: int) -> int:
    value = lower_bound
    while gcd(value, d) != 1:
        value += 1
    return value


def deterministic_tail_set(d: int, selected_tail: int) -> tuple[int, ...]:
    tails = [selected_tail]
    candidate = 1
    while len(tails) < d:
        if gcd(candidate, d) == 1 and candidate not in tails:
            tails.append(candidate)
        candidate += 1
    return tuple(tails)


def audit_constructed_row(
    d: int,
    core: tuple[int, ...],
    tails: tuple[int, ...],
    selected_tail: int,
    core_phase: Fraction,
    core_margin: Fraction,
    label: str,
) -> tuple[Fraction, tuple[int, ...]]:
    require(2 <= d <= 7, f"{label}:d-range")
    require(len(core) == 13 - d, f"{label}:core-cardinality")
    require(len(tails) == d and len(set(tails)) == d, f"{label}:tail-cardinality")
    require(all(gcd(tail, d) == 1 for tail in tails), f"{label}:unit-typing")
    require(selected_tail in tails, f"{label}:selected-present")

    radius = core_margin - DELTA
    require(radius >= Fraction(d, 14 * (14 - d)), f"{label}:lower-lrc-margin")
    quotient_phase = nearest_selected_safe_phase(core_phase, d, selected_tail)
    movement = circle_norm(quotient_phase - core_phase)
    require(
        movement <= Fraction(d, 14 * selected_tail),
        f"{label}:safe-band-distance",
    )
    require(
        max(core) * movement <= radius,
        f"{label}:lipschitz-budget",
    )

    for speed in core:
        require(
            circle_norm(speed * quotient_phase) >= DELTA,
            f"{label}:core-after-move:{speed}",
        )

    selected_bad = bad_labels(d, selected_tail, quotient_phase)
    require(not selected_bad, f"{label}:selected-tail-zero-mask")

    burned: set[int] = set()
    for tail in tails:
        tail_bad = bad_labels(d, tail, quotient_phase)
        if tail != selected_tail:
            require(len(tail_bad) <= 1, f"{label}:unit-singleton:{tail}")
        burned.update(tail_bad)
    require(len(burned) <= d - 1, f"{label}:boundary-capacity-drop")

    free_labels = tuple(label_index for label_index in range(d) if label_index not in burned)
    require(bool(free_labels), f"{label}:free-label")
    row = tuple(d * speed for speed in core) + tails
    for free_label in free_labels:
        physical_phase = (quotient_phase + free_label) / d
        require(
            min(circle_norm(speed * physical_phase) for speed in row) >= DELTA,
            f"{label}:literal-witness:{free_label}",
        )
    return quotient_phase, free_labels


def main() -> None:
    source_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(source_tree)), "no-assert-nodes")
    require(
        not any(
            isinstance(node, ast.Constant) and isinstance(node.value, float)
            for node in ast.walk(source_tree)
        ),
        "no-float-literals",
    )

    # Exact symbolic safe-band iff and translated-grid capacity.
    symbolic_phases = fraction_universe(84)
    symbolic_band_checks = 0
    translated_grid_checks = 0
    endpoint_checks = 0
    epsilon = Fraction(1, 1008)
    for d in range(2, 8):
        lo, hi = safe_band(d)
        for z in symbolic_phases:
            actual = min(selected_lift_clearances(d, z)) >= DELTA
            expected = lo <= z <= hi
            require(actual == expected, f"band-iff:d={d}:z={z}")
            symbolic_band_checks += 1

            danger_count = sum(
                circle_norm(z + Fraction(label, d)) < DELTA for label in range(d)
            )
            require(danger_count <= 1, f"translated-grid:d={d}:z={z}")
            translated_grid_checks += 1

        require(min(selected_lift_clearances(d, lo)) == DELTA, f"left-endpoint:d={d}")
        require(min(selected_lift_clearances(d, hi)) == DELTA, f"right-endpoint:d={d}")
        require(
            min(selected_lift_clearances(d, (lo - epsilon) % 1)) < DELTA,
            f"left-outside-strict:d={d}",
        )
        require(
            min(selected_lift_clearances(d, (hi + epsilon) % 1)) < DELTA,
            f"right-outside-strict:d={d}",
        )
        endpoint_checks += 4
        if d < 7:
            require(
                min(selected_lift_clearances(d, lo + epsilon)) >= DELTA,
                f"left-inside-closed:d={d}",
            )
            require(
                min(selected_lift_clearances(d, hi - epsilon)) >= DELTA,
                f"right-inside-closed:d={d}",
            )
            endpoint_checks += 2
        else:
            require(lo == hi == Fraction(1, 2), "d7-singleton-band")
            endpoint_checks += 1

    # The safe-band preimage is never farther than d/(14w).  The y=0
    # controls attain equality, so the metric constant itself is sharp.
    preimage_phases = fraction_universe(24)
    preimage_checks = 0
    preimage_equality_checks = 0
    for d in range(2, 8):
        lo, hi = safe_band(d)
        for selected_tail in range(1, 36):
            if gcd(selected_tail, d) != 1:
                continue
            equality_phase = nearest_selected_safe_phase(Fraction(0), d, selected_tail)
            require(
                circle_norm(equality_phase) == Fraction(d, 14 * selected_tail),
                f"preimage-equality:d={d}:w={selected_tail}",
            )
            preimage_equality_checks += 1
            for phase in preimage_phases:
                moved = nearest_selected_safe_phase(phase, d, selected_tail)
                residue = (selected_tail * moved) % 1
                require(lo <= residue <= hi, f"preimage-safe:d={d}:w={selected_tail}:y={phase}")
                require(
                    circle_norm(moved - phase) <= Fraction(d, 14 * selected_tail),
                    f"preimage-bound:d={d}:w={selected_tail}:y={phase}",
                )
                preimage_checks += 2

    # Complete deterministic core atlas: every C subset {1,...,11} of the
    # cardinality 13-d.  This is exactly 1,024 cores over d=2,...,7.
    core_counts: dict[int, int] = {}
    minimum_excess: dict[int, Fraction] = {}
    constructed_rows = 0
    for d in range(2, 8):
        lower_margin = Fraction(1, 14 - d)
        count = 0
        excesses: list[Fraction] = []
        for core in combinations(range(1, 12), 13 - d):
            core_margin, core_phase = exact_core_maximin(core)
            require(core_margin >= lower_margin, f"core-lrc:d={d}:C={core}")
            excesses.append(core_margin - lower_margin)
            selected_tail = least_unit_at_least((14 - d) * max(core), d)
            tails = deterministic_tail_set(d, selected_tail)
            audit_constructed_row(
                d,
                core,
                tails,
                selected_tail,
                core_phase,
                core_margin,
                f"atlas:d={d}:C={core}",
            )
            count += 1
            constructed_rows += 1
        core_counts[d] = count
        minimum_excess[d] = min(excesses)

    # d=5 equality-wall positive control.  The raw THM-2060/2064 capacity is
    # 5*(1/5)=1, so it cannot certify.  Moving the AP8 core phase tears the
    # selected tail's mask to zero; the other four tails occupy labels 0..3.
    d5_core = tuple(range(1, 9))
    d5_tails = (72, 1, 9, 7, 8)
    d5_selected = 72
    d5_core_phase = Fraction(1, 9)
    d5_core_margin = Fraction(1, 9)
    d5_phase, d5_free = audit_constructed_row(
        5,
        d5_core,
        d5_tails,
        d5_selected,
        d5_core_phase,
        d5_core_margin,
        "d5-positive-control",
    )
    d5_masks = tuple((tail, bad_labels(5, tail, d5_phase)) for tail in d5_tails)
    d5_capacity = sum(Fraction((5 + 6) // 7, 5) for _ in d5_tails)
    require(d5_capacity == 1, "d5-raw-capacity-equality")
    require(d5_phase == Fraction(107, 1008), "d5-phase")
    require(d5_masks == ((72, ()), (1, (0,)), (9, (1,)), (7, (2,)), (8, (3,))), "d5-masks")
    require(d5_free == (4,), "d5-unique-free-label")
    d5_witness = (d5_phase + d5_free[0]) / 5
    d5_row = tuple(5 * speed for speed in d5_core) + d5_tails
    d5_clearance = min(circle_norm(speed * d5_witness) for speed in d5_row)
    require(d5_witness == Fraction(4139, 5040), "d5-witness")
    require(d5_clearance == Fraction(107, 1008), "d5-clearance")

    print("LRC14_UNIT_EQUALITY_WALL_SAFE_BAND_EXACT_AUDIT")
    print("STATUS=FINITE-EXACT_CERTIFICATE_AUDIT_NOT_LRC14")
    print("RANGE=d=2,3,4,5,6,7; |C|=13-d; |T|=d; gcd(t,d)=1 for every tail")
    print("SAFE_BAND={z: selected tail is safe on all d lifts}=[d/14,1-d/14]")
    print("PREIMAGE_DISTANCE<=d/(14w)")
    print("COFINAL_CONDITION=max(T)>=(14-d)*max(C)")
    print("ACTUAL_MARGIN_CONDITION=w>=d*R/[14*(M(C)-1/14)]")
    print(
        "THRESHOLDS="
        + ",".join(f"d{d}:{14-d}R" for d in range(2, 8))
    )
    print(f"SYMBOLIC_PHASES={len(symbolic_phases)} BAND_CHECKS={symbolic_band_checks}")
    print(f"TRANSLATED_GRID_CHECKS={translated_grid_checks} ENDPOINT_CHECKS={endpoint_checks}")
    print(
        f"PREIMAGE_PHASES={len(preimage_phases)} PREIMAGE_CHECKS={preimage_checks} "
        f"PREIMAGE_EQUALITY_CHECKS={preimage_equality_checks}"
    )
    for d in range(2, 8):
        print(
            f"CORE_ATLAS d={d} cores={core_counts[d]} "
            f"min_excess_over_1/{14-d}={minimum_excess[d]}"
        )
    print(f"CORE_ATLAS_TOTAL={constructed_rows}")
    print("D5_NEW_CONSEQUENCE=max(T)>=9*max(C)")
    print(f"D5_CONTROL C={d5_core} T={d5_tails} R=8 selected=72 threshold=72")
    print(
        f"D5_CONTROL core_phase={d5_core_phase} moved_phase={d5_phase} "
        f"movement={circle_norm(d5_phase-d5_core_phase)}"
    )
    print(f"D5_CONTROL raw_capacity={d5_capacity} masks={d5_masks} free_labels={d5_free}")
    print(f"D5_CONTROL witness={d5_witness} clearance={d5_clearance}")
    print(f"CHECKS={CHECKS} RESULT=PASS")


if __name__ == "__main__":
    main()
