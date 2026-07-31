#!/usr/bin/env python3
"""Exact extended-address determinant-sector probe on one gain-two wing cell.

The selected physical cell is (s,t,e)=(0,10,3), the D-multiplier corner of
the all-81 wing chamber.  Its source and target wings have a common
first-failure label (shifted E3) and its raw target coefficient is twice its
raw source coefficient.  We retain the source and target carrier harmonics
separately, reconstruct the full

    (left endpoint, source harmonic, right endpoint, target harmonic)

pullback of size 13^6, and test one fixed (q,Delta,k,l) sector.

The computation is a characteristic-zero nonvanishing certificate through
the first exact cyclotomic specialization already certified in THM-2625.  It
does not make the selected source and target carriers translates, identify a
physical transition, or exclude an LRC row.
"""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parents[2]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


PINNED = {
    "lrc14_extended_carrier_endpoint_lib.py":
        "4b3f9f195b1634e1e84a1bc8bccb878a1c8c44aec13f24d197f92547c9e36c57",
    "lrc14_root_zero_full_target_semantic_clutch_20260728.py":
        "208f71020efa19fa47f66d2da061ab03fa7bc87beeb077b4008c069f499736d8",
}
for name, expected in PINNED.items():
    actual = sha256(
        (COMP / name).read_bytes().replace(b"\r\n", b"\n")
    ).hexdigest()
    require(actual == expected, f"pinned dependency changed: {name}")

import lrc14_extended_carrier_endpoint_lib as base
import lrc14_root_zero_full_target_semantic_clutch_20260728 as physical


P = 13
LABEL = (0, 10)
PHYSICAL_CLOCK = 3
Q_TARGET = (1, 0)
DELTA_TARGET = 1
K_TARGET = 0
L_TARGET = 0


def weighted_intersection(pieces, intervals):
    return tuple(
        physical.relative.private.old.intersect_weighted_union(
            pieces, intervals
        )
    )


def weighted_mass(pieces):
    return sum((right - left) * weight for left, right, weight in pieces)


def weighted_endpoint_sum(pieces, frequency, embedding):
    prime = embedding[0]
    total = 0
    for left, right, weight in pieces:
        value = base.endpoint.endpoint_sum(
            ((left, right),), frequency, (embedding,)
        )[0]
        total = (total + weight * value) % prime
    return total


def weighted_x_sweep(pieces, terminal, q_starts, frequency, embedding, tabs):
    prime = embedding[0]
    total = 0
    overlap = 0
    for left, right, weight in pieces:
        value, local_overlap = base.endpoint.x_sweep(
            ((left, right),), terminal, q_starts, frequency,
            (embedding,), tabs,
        )
        total = (total + weight * value[0]) % prime
        overlap += weight * local_overlap
    return total, overlap


def build_selected_wing():
    module, _prefixes, _, _, rails, present, _starts = (
        physical.relative.lift.m.core.build_carrier_data()
    )
    pair_prefixes = physical.relative.private.build_pair_prefixes(module)
    _raw_source, _raw_target, _rail_pairs, overlap_details = (
        physical.overlap.overlap_vectors(
            module, pair_prefixes, rails, present, rail_index=8
        )
    )
    full_module = physical.target.load_present_module()
    e3 = physical.target.exclusive_source(full_module, 3)
    fork = physical.target.deepest_fork(full_module)
    clock_combs = tuple(
        full_module.make_comb(
            full_module.C1, 182, 26 * clock - 13, 26 * clock + 13
        )
        for clock in range(7)
    )
    q_pairs = physical.q_restricted_pair_prefixes(
        full_module, pair_prefixes, fork
    )
    s, t = LABEL
    source_base, target_base = overlap_details[PHYSICAL_CLOCK]
    section = physical.target.source_present_section(
        full_module, e3, PHYSICAL_CLOCK, s, t, clock_combs
    )
    one_source = weighted_intersection(source_base, section)
    one_target = weighted_intersection(target_base, section)
    common_source = physical.common_physical_cut(
        source_base, section, -1
    )
    common_target = physical.common_physical_cut(
        target_base, section, +1
    )
    source_wing = physical.subtract_weighted(one_source, common_source)
    target_wing = physical.subtract_weighted(one_target, common_target)

    universe = ((0, full_module.T),)
    factors = (
        e3,
        clock_combs[PHYSICAL_CLOCK],
        full_module.subtract_comb(
            universe, full_module.W[1], 182,
            -14 * s - 13, -14 * s + 13,
        ),
        full_module.subtract_comb(
            universe, full_module.W[2], 182,
            -14 * t - 13, -14 * t + 13,
        ),
        full_module.subtract_comb(
            universe, full_module.C2, 182,
            14 * s - 13, 14 * s + 13,
        ),
        full_module.subtract_comb(
            universe, full_module.C3, 182,
            14 * t - 13, 14 * t + 13,
        ),
    )

    def first_failures(one, direction):
        remaining = one
        origins = []
        for factor in factors:
            shifted = physical.overlap.shift_union(
                factor, direction * physical.SHIFT
            )
            origins.append(
                weighted_intersection(
                    remaining, physical.overlap.complement(shifted)
                )
            )
            remaining = weighted_intersection(remaining, shifted)
        return tuple(origins), remaining

    source_origins, source_remainder = first_failures(one_source, -1)
    target_origins, target_remainder = first_failures(one_target, +1)
    require(source_remainder == common_source, "source origin remainder")
    require(target_remainder == common_target, "target origin remainder")
    require(
        source_origins[0] == source_wing
        and target_origins[0] == target_wing
        and all(not pieces for pieces in source_origins[1:])
        and all(not pieces for pieces in target_origins[1:]),
        "selected gain-two cell is not a pure shifted-E3 failure",
    )

    pair = q_pairs[PHYSICAL_CLOCK]
    source_value = physical.relative.private.delayed_carry_pair(
        source_wing, pair, {}
    )[12][1]
    target_value = physical.relative.private.delayed_carry_pair(
        target_wing, pair, {}
    )[6][1]
    require(source_value != 0 and target_value == 2 * source_value,
            "raw gain-two coefficient changed")
    require(
        physical.overlap.shift_weighted(
            source_wing, physical.SHIFT
        ) != target_wing,
        "selected wing unexpectedly became a chart translate",
    )
    require(
        weighted_mass(source_wing) != weighted_mass(target_wing),
        "mass hostile unexpectedly disappeared",
    )
    return (
        source_wing,
        target_wing,
        source_value,
        target_value,
        tuple(weighted_mass(origin) for origin in source_origins),
        tuple(weighted_mass(origin) for origin in target_origins),
    )


def evaluate_left(ell, present, carrier, terminal, q_starts, tabs, embedding):
    restricted = weighted_intersection(
        carrier, present
    )
    value, _overlap = weighted_x_sweep(
        restricted, terminal, q_starts, base.X0, embedding, tabs
    )
    prime, root = embedding
    zeta = pow(root, base.NRED // P, prime)
    return value * pow(zeta, ell[base.DEEP], prime) % prime


def evaluate_right(ell, present, carrier, embedding):
    restricted = weighted_intersection(
        carrier, present
    )
    return weighted_endpoint_sum(
        restricted, -base.Y0, embedding
    )


def main():
    (
        source_wing,
        target_wing,
        source_value,
        target_value,
        source_origin_masses,
        target_origin_masses,
    ) = build_selected_wing()
    terminal = tuple(base.endpoint.build_set(base.PAT_Q12, base.ZERO))
    q_starts = [left for left, _right in terminal]
    embedding = base.endpoint.MODS[0]
    prime, root = embedding
    tabs = base.endpoint.make_tabs(
        terminal, base.X0, (embedding,)
    )
    zeta = pow(root, base.NRED // P, prime)
    powers = tuple(pow(zeta, exponent, prime) for exponent in range(P))
    present_sets = base.present_cache()
    source_shifts = tuple(
        physical.overlap.shift_weighted(
            source_wing, -harmonic * (base.T // P)
        )
        for harmonic in range(P)
    )
    target_shifts = tuple(
        physical.overlap.shift_weighted(
            target_wing, harmonic * (base.T // P)
        )
        for harmonic in range(P)
    )

    # Dual banks on the separately allocated 13^3 left and right quotients.
    left_dual = {}
    right_dual = {}
    for address, ell in base.REPS.items():
        present = present_sets[address]
        for harmonic in range(P):
            left_dual[address, harmonic] = evaluate_left(
                ell, present, source_shifts[harmonic],
                terminal, q_starts, tabs, embedding
            )
            right_dual[address, harmonic] = evaluate_right(
                ell, present, target_shifts[harmonic], embedding
            )

    # Representative-gauge hostile: ell -> ell+W must be accompanied by
    # source harmonic +1 and target harmonic -1.
    for address, harmonic in (
        ((0, 0), 0), ((1, 0), 2), ((0, 1), 7), ((3, 7), 12)
    ):
        ell = base.REPS[address]
        shifted_ell = tuple(
            (ell[index] + base.WMOD[index]) % P
            for index in range(9)
        )
        shifted_present = tuple(
            base.endpoint.build_set(base.PAT_E3, shifted_ell)
        )
        shifted_source = physical.overlap.shift_weighted(
            source_wing, -(harmonic + 1) * (base.T // P)
        )
        shifted_target = physical.overlap.shift_weighted(
            target_wing, (harmonic - 1) * (base.T // P)
        )
        require(
            evaluate_left(
                shifted_ell, shifted_present, shifted_source, terminal,
                q_starts, tabs, embedding,
            ) == left_dual[address, harmonic],
            "left representative gauge failed",
        )
        require(
            evaluate_right(
                shifted_ell, shifted_present, shifted_target, embedding,
            ) == right_dual[address, harmonic],
            "right representative gauge failed",
        )

    def left_primal(point, harmonic):
        total = 0
        for address in base.KEYS:
            pairing = (
                address[0] * point[0] + address[1] * point[1]
            ) % P
            for dual_harmonic in range(P):
                exponent = -pairing - dual_harmonic * harmonic
                total += (
                    left_dual[address, dual_harmonic]
                    * powers[exponent % P]
                )
        return total % prime

    def right_primal(point, harmonic):
        total = 0
        for address in base.KEYS:
            pairing = (
                address[0] * point[0] + address[1] * point[1]
            ) % P
            for dual_harmonic in range(P):
                exponent = pairing - dual_harmonic * harmonic
                total += (
                    right_dual[address, dual_harmonic]
                    * powers[exponent % P]
                )
        return total % prime

    # One actual nondegenerate determinant sector in the extended pullback.
    edges = []
    sector = 0
    q0, q1 = Q_TARGET
    for right0 in range(P):
        # det((q0,q1),R)=q0*R1-q1*R0=1.
        if q0:
            right1 = (
                DELTA_TARGET + q1 * right0
            ) * pow(q0, -1, P) % P
            rights = ((right0, right1),)
        else:
            rights = tuple(
                (right0, right1)
                for right1 in range(P)
                if (q0 * right1 - q1 * right0) % P == DELTA_TARGET
            )
        for right in rights:
            left = ((right[0] + q0) % P, (right[1] + q1) % P)
            left_value = left_primal(left, K_TARGET)
            right_value = right_primal(right, L_TARGET)
            coefficient = left_value * right_value % prime
            edges.append((left, right, left_value, right_value, coefficient))
            sector = (sector + coefficient) % prime
    require(len(edges) == P, "nondegenerate sector edge count")
    require(
        all(
            (
                (left[0] - right[0]) % P,
                (left[1] - right[1]) % P,
            ) == Q_TARGET
            and (left[0] * right[1] - left[1] * right[0]) % P
            == DELTA_TARGET
            for left, right, *_values in edges
        ),
        "endpoint sector typing",
    )
    require(
        sum(bool(value) for value in left_dual.values()) == 641
        and sum(bool(value) for value in right_dual.values()) == 510,
        "extended dual support census changed",
    )
    require(
        all(row[-1] for row in edges)
        and sector == 64_970_359_711_488_738,
        "selected determinant-sector certificate changed",
    )

    print("extended-address selected-wing determinant-sector probe")
    print(f"dependency_sha256={tuple(PINNED.items())}")
    print(f"selected_label_(s,t,e)={(*LABEL, PHYSICAL_CLOCK)}")
    print("first_failure_order=(shifted_E3,shifted_clock,shifted_q1,"
          "shifted_q2,shifted_c2,shifted_c3)")
    print(f"source_first_failure_masses={source_origin_masses}")
    print(f"target_first_failure_masses={target_origin_masses}")
    print(f"raw_selected_coefficients=(source={source_value},"
          f"target={target_value},target/source=2)")
    print(f"weighted_piece_counts=(source={len(source_wing)},"
          f"target={len(target_wing)})")
    print(f"weighted_masses=(source={weighted_mass(source_wing)},"
          f"target={weighted_mass(target_wing)})")
    print("chart_translate_equal=False")
    print(f"dual_support=(left={sum(bool(v) for v in left_dual.values())}/2197,"
          f"right={sum(bool(v) for v in right_dual.values())}/2197)")
    print("representative_gauge=(ell,a,b)~(ell+W,a+1,b-1):PASS")
    print("pullback_coordinates=(L,k,R,l); "
          "equivalently=(q=L-R,k,l,R)=G_full_x_V")
    print(f"pullback_size=13^6={P**6}")
    print(f"tested_sector=(q={Q_TARGET},Delta={DELTA_TARGET},"
          f"k={K_TARGET},l={L_TARGET})")
    print(f"tested_sector_edge_support={sum(bool(row[-1]) for row in edges)}/13")
    print(f"tested_sector_value_mod_{prime}={sector}")
    print(f"tested_sector_nonzero={sector != 0}")
    print(
        "scope=one exact weighted factor-labelled wing cell and one "
        "extended endpoint sector; no translate/action/current positivity/"
        "row exclusion/LRC14 conclusion"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
