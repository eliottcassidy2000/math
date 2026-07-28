#!/usr/bin/env python3
"""Exact two-target repair of the canonical E3/present incompatibility.

The later packet scouts use the one-target section t=0 of the lawful
THM-2365/2407 target action.  That section retains the unshifted c3-safe
factor and is therefore disjoint from the exclusive source E3, which
requires c3-danger.  This companion restores the second target coordinate

    eta = e_c2 - e_q1,        lambda = e_c3 - e_q2

and exhausts all (ell,s,t) in F_7 x F_13 x F_13.  It then uses the exact
prefix correlation identity to impose the actual ordinary-D^6 terminal fork
Q_(3,{1,2}), rather than stopping at source compatibility.

All sets are finite unions of rational open intervals on the canonical T
grid.  There are no floating-point calculations and no Python assertions.
"""

from bisect import bisect_right
from fractions import Fraction
from pathlib import Path
import hashlib
import importlib.util


ROOT = Path(__file__).resolve().parents[1]
P = 13
R = P**6
CARRIER_PATH = (
    ROOT / "04-computation"
    / "lrc14_replica_dichotomy_typed_row_opus_20260727.py"
)
CARRIER_SHA256 = (
    "6ba64a68a9fd008d2e06949b1f1cf75012f1f4e734f75f55ce0af58ae20ad7b9"
)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def load_present_module():
    require(
        hashlib.sha256(CARRIER_PATH.read_bytes()).hexdigest()
        == CARRIER_SHA256,
        "audited canonical present-packet dependency changed",
    )
    spec = importlib.util.spec_from_file_location(
        "two_target_present_base", CARRIER_PATH
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def intersect_sorted(left, right):
    out = []
    i = j = 0
    while i < len(left) and j < len(right):
        a = max(left[i][0], right[j][0])
        b = min(left[i][1], right[j][1])
        if a < b:
            out.append((a, b))
        if left[i][1] <= right[j][1]:
            i += 1
        else:
            j += 1
    return out


def interval_mass(intervals):
    return sum(right - left for left, right in intervals)


def common_safe(module, intervals):
    intervals = module.subtract_comb(
        intervals, module.W[module.GUARD], 91, -13, 13
    )
    for index in module.UNIT_IDX:
        intervals = module.subtract_comb(
            intervals, module.W[index], 182, -13, 13
        )
    return intervals


def exclusive_source(module, source):
    blockers = (module.C1, module.C2, module.C3)
    intervals = module.make_comb(blockers[source - 1], 182, -13, 13)
    intervals = common_safe(module, intervals)
    for index, blocker in enumerate(blockers, start=1):
        if index != source:
            intervals = module.subtract_comb(
                intervals, blocker, 182, -13, 13
            )
    return intervals


def deepest_fork(module):
    """Q_(3,{1,2}): c1,c2 dangerous and c3 safe."""
    intervals = module.make_comb(module.C1, 182, -13, 13)
    intervals = module.intersect_comb(
        intervals, module.C2, 182, -13, 13
    )
    intervals = common_safe(module, intervals)
    return module.subtract_comb(intervals, module.C3, 182, -13, 13)


def source_present_section(module, source_intervals, source_clock, s, t,
                           clock_comb):
    """Intersect E_j with the full U_(s,t)d_(1,ell) present section.

    Common guard and q3,q4,q5 safeties are already present on E_j.  Only the
    two lawful graft pairs need to be shifted here.
    """
    intervals = intersect_sorted(source_intervals, clock_comb[source_clock])
    intervals = module.subtract_comb(
        intervals, module.W[1], 182, -14 * s - 13, -14 * s + 13
    )
    intervals = module.subtract_comb(
        intervals, module.W[2], 182, -14 * t - 13, -14 * t + 13
    )
    intervals = module.subtract_comb(
        intervals, module.C2, 182, 14 * s - 13, 14 * s + 13
    )
    intervals = module.subtract_comb(
        intervals, module.C3, 182, 14 * t - 13, 14 * t + 13
    )
    return intervals


def dilated_overlap(module, intervals, target_prefix, target_mass):
    length = interval_mass(intervals)
    if not length:
        return Fraction(0)
    acc_r, acc_prefix = module.sweep_acc(
        intervals, R % module.T, *target_prefix
    )
    return module.IR_from_acc(
        R, length, target_mass, acc_r, acc_prefix
    )


def first_open_witness(module, source_present, target):
    """Return one rational x in source_present with {R*x} in target."""
    target_starts = [left for left, _ in target]
    target_ends = [right for _, right in target]
    for left, right in source_present:
        lower = left * R
        upper = right * R
        first_turn = lower // module.T - 1
        last_turn = upper // module.T + 1
        for turn in range(first_turn, last_turn + 1):
            local_low = lower - turn * module.T
            local_high = upper - turn * module.T
            index = bisect_right(target_ends, local_low)
            if index >= len(target):
                continue
            a = max(lower, turn * module.T + target_starts[index])
            b = min(upper, turn * module.T + target_ends[index])
            if a < b:
                return Fraction(a + b, 2 * R * module.T)
    raise RuntimeError("positive overlap had no open rational witness")


def strict_member(value, intervals, denominator):
    coordinate = value * denominator
    starts = [left for left, _ in intervals]
    index = bisect_right(starts, coordinate) - 1
    return (
        index >= 0
        and Fraction(intervals[index][0]) < coordinate
        < Fraction(intervals[index][1])
    )


def primitive_fourier_coordinates(values, frequency):
    """Power-basis coordinates of sum_t values[t] zeta^(frequency*t).

    For prime P, reduce zeta^(P-1) by
    1+zeta+...+zeta^(P-1)=0.  A nonzero returned vector is an exact
    certificate that the corresponding primitive Fourier coefficient does
    not vanish.
    """
    require(len(values) == P and 0 < frequency < P,
            "invalid primitive Fourier request")
    coordinates = [Fraction(0) for _ in range(P - 1)]
    for index, value in enumerate(values):
        exponent = frequency * index % P
        if exponent == P - 1:
            for slot in range(P - 1):
                coordinates[slot] -= value
        else:
            coordinates[exponent] += value
    return tuple(coordinates)


def main():
    module = load_present_module()
    require(module.T == 297836897838480, "canonical present grid changed")

    sources = tuple(exclusive_source(module, j) for j in (1, 2, 3))
    fork = deepest_fork(module)
    fork_prefix = module.make_prefix(fork)
    fork_mass = interval_mass(fork)
    clock_comb = tuple(
        module.make_comb(module.C1, 182, 26 * ell - 13, 26 * ell + 13)
        for ell in range(7)
    )

    # The inherited t=0 section has a sharp global compatibility matrix.
    one_target_positive = []
    for source, source_intervals in enumerate(sources, start=1):
        positive = []
        for ell in range(7):
            for s in range(P):
                intervals = source_present_section(
                    module, source_intervals, ell, s, 0, clock_comb
                )
                if interval_mass(intervals):
                    positive.append((ell, s))
        one_target_positive.append(tuple(positive))

    require(
        one_target_positive[0] == tuple((0, s) for s in range(P)),
        "E1/one-target compatibility law changed",
    )
    require(
        one_target_positive[1]
        == tuple((ell, s) for ell in range(1, 7) for s in range(1, P)),
        "E2/one-target compatibility law changed",
    )
    require(
        not one_target_positive[2],
        "the t=0 present section unexpectedly met E3",
    )

    source_rows = []
    semantic_rows = []
    witness_intervals = None
    witness_key = None
    for ell in range(7):
        for s in range(P):
            for t in range(P):
                intervals = source_present_section(
                    module, sources[2], ell, s, t, clock_comb
                )
                mass = interval_mass(intervals)
                if not mass:
                    continue
                source_rows.append((ell, s, t, mass, len(intervals)))
                semantic_mass = dilated_overlap(
                    module, intervals, fork_prefix, fork_mass
                )
                if semantic_mass:
                    semantic_rows.append((ell, s, t, semantic_mass))
                    if witness_intervals is None:
                        witness_intervals = intervals
                        witness_key = (ell, s, t)

    expected_keys = tuple(
        (ell, s, t)
        for ell in range(1, 7)
        for s in range(P)
        for t in range(1, P)
    )
    require(
        tuple(row[:3] for row in source_rows) == expected_keys,
        "full two-target E3 source-attachment law changed",
    )
    require(
        tuple(row[:3] for row in semantic_rows) == expected_keys,
        "some full two-target E3 section lost the D^6 terminal fork",
    )

    semantic_lookup = {
        (ell, s, t): mass for ell, s, t, mass in semantic_rows
    }
    marked_moving_count = 0
    for ell in range(1, 7):
        for s in range(P):
            profile = tuple(
                semantic_lookup.get((ell, s, t), Fraction(0))
                for t in range(P)
            )
            require(profile[0] == 0 and all(value > 0 for value in profile[1:]),
                    "marked two-target semantic profile lost its sharp zero")
            for frequency in range(1, P):
                coordinates = primitive_fourier_coordinates(
                    profile, frequency
                )
                require(any(coordinates),
                        "marked deepest-target character vanished")
                marked_moving_count += 1
    require(marked_moving_count == 6 * P * (P - 1),
            "marked moving-character count changed")

    aggregate_target_profile = tuple(
        sum(
            (
                semantic_lookup.get((ell, s, t), Fraction(0))
                for ell in range(7)
                for s in range(P)
            ),
            Fraction(0),
        )
        for t in range(P)
    )
    require(
        aggregate_target_profile[0] == 0
        and all(value > 0 for value in aggregate_target_profile[1:]),
        "aggregated deepest-target profile lost its sharp zero",
    )
    require(
        all(
            aggregate_target_profile[t]
            == aggregate_target_profile[(-t) % P]
            for t in range(P)
        ),
        "aggregated deepest-target profile lost inversion symmetry",
    )
    aggregate_moving_count = sum(
        any(primitive_fourier_coordinates(
            aggregate_target_profile, frequency
        ))
        for frequency in range(1, P)
    )
    require(aggregate_moving_count == P - 1,
            "aggregated deepest-target character vanished")

    min_source = min(source_rows, key=lambda row: row[3])
    max_source = max(source_rows, key=lambda row: row[3])
    min_semantic = min(semantic_rows, key=lambda row: row[3])
    max_semantic = max(semantic_rows, key=lambda row: row[3])
    witness = first_open_witness(module, witness_intervals, fork)
    require(
        strict_member(witness, witness_intervals, module.T),
        "witness left the full present/source section",
    )
    endpoint = (R * witness) % 1
    require(
        strict_member(endpoint, fork, module.T),
        "witness endpoint left the deepest fork",
    )

    # The reason for the old zero is definition-level, not statistical:
    # every t=0 section explicitly subtracts the unshifted c3 danger comb.
    c3_danger = module.make_comb(module.C3, 182, -13, 13)
    inherited_c3_overlap = sum(
        interval_mass(intersect_sorted(
            module.build_F(ell, s), c3_danger
        ))
        for ell in range(7)
        for s in range(P)
    )
    require(inherited_c3_overlap == 0, "t=0 c3-safe law changed")

    print("LRC14 TWO-TARGET PRESENT / SEMANTIC ATTACHMENT AUDIT")
    print(f"row={module.W} p={P} R={R} T={module.T}")
    print(f"canonical_present_dependency_sha256={CARRIER_SHA256}")
    print(
        "one_target_t0_source_section_counts="
        f"E1:{len(one_target_positive[0])},"
        f"E2:{len(one_target_positive[1])},"
        f"E3:{len(one_target_positive[2])}"
    )
    print(
        "one_target_index_laws="
        "E1 iff ell=0; E2 iff ell!=0 and s!=0; E3 never"
    )
    print(
        "uniform_t0_mechanism=every F_(ell,s,0) contains c3-safe; "
        "exclusive E3 contains c3-danger; strict intersection empty"
    )
    print(
        f"two_target_E3_source_sections={len(source_rows)} "
        "index_law=ell!=0,t!=0,all_s"
    )
    print(
        f"two_target_E3_source_mass_min={min_source} "
        f"max={max_source}"
    )
    print(
        f"two_target_E3_to_D6_fork_sections={len(semantic_rows)} "
        "index_law=ell!=0,t!=0,all_s"
    )
    print(
        f"two_target_semantic_mass_min={min_semantic} "
        f"max={max_semantic}"
    )
    print(
        f"marked_(ell,s)_deep_target_moving_characters="
        f"{marked_moving_count}/{6 * P * (P - 1)}"
    )
    print(f"aggregated_deep_target_profile={aggregate_target_profile}")
    print(
        f"aggregated_deep_target_moving_characters="
        f"{aggregate_moving_count}/{P - 1}"
    )
    print(
        f"strict_witness_key={witness_key} x={witness} "
        f"D6x={endpoint}"
    )
    print(
        "scope=canonical typed row and full lawful two-target present sheet;"
        "positive semantic mass and deepest-target coefficients;"
        "no paired left relation residue;no affine-cycle attachment;"
        "no row exclusion;no LRC14 conclusion"
    )
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
