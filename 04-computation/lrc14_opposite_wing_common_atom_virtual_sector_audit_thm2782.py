#!/usr/bin/env python3
"""Independent THM-2782 opposite-wing and virtual-sector audit.

This control does not import either discovery companion.  It reconstructs
the aligned one-sided supports A and B, identifies the common carrier and
opposite wing cofibres as

    A cap B, A \ B, B \ A,

on all 81 labels and seven clocks.  It then independently rebuilds the
selected D=(0,10,3) first-failure decomposition and the extended endpoint
sector in both certified THM-2625 finite-field embeddings.

The representative gauge is checked on every 169*13 dual cell of both
endpoint banks.  A fixed-carrier hostile is also required to fail.
"""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path
import ast
import sys


ROOT = Path(__file__).resolve().parents[1]
COMP = ROOT / "04-computation"
sys.path.insert(0, str(COMP))


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


PINNED = {
    "lrc14_root_zero_full_target_semantic_clutch_20260728.py":
        "208f71020efa19fa47f66d2da061ab03fa7bc87beeb077b4008c069f499736d8",
    "lrc14_natural_single_sheet_root_zero_target_bank_20260728.py":
        "591241e912c452b1985c7fc700884183a6e50440a579e1d853123d276cd416a8",
    "lrc14_extended_carrier_endpoint_lib.py":
        "4b3f9f195b1634e1e84a1bc8bccb878a1c8c44aec13f24d197f92547c9e36c57",
}
for name, expected in PINNED.items():
    actual = sha256(
        (COMP / name).read_bytes().replace(b"\r\n", b"\n")
    ).hexdigest()
    require(actual == expected, f"pinned dependency changed: {name}")

import lrc14_root_zero_full_target_semantic_clutch_20260728 as physical
import lrc14_natural_single_sheet_root_zero_target_bank_20260728 as natural
import lrc14_extended_carrier_endpoint_lib as endpoint_base


P = 13
COMMON_S = (0, 1, 2, 3, 8, 9, 10, 11, 12)
COMMON_T = (3, 4, 5, 6, 7, 8, 9, 10, 11)
SELECTED = (0, 10, 3)
Q_TARGET = (1, 0)
DELTA_TARGET = 1
K_TARGET = 0
L_TARGET = 0
PERIOD = physical.T
SHIFT = physical.SHIFT


def merge_intervals(intervals):
    result = []
    for left, right in sorted(intervals):
        if left >= right:
            continue
        if result and left <= result[-1][1]:
            result[-1] = (result[-1][0], max(result[-1][1], right))
        else:
            result.append((left, right))
    return tuple(result)


def normalize_weighted(pieces):
    result = []
    for left, right, weight in sorted(pieces):
        if left >= right or not weight:
            continue
        if result and left < result[-1][1]:
            require(
                weight == result[-1][2],
                "overlapping weighted pieces have unequal weights",
            )
            result[-1] = (
                result[-1][0], max(result[-1][1], right), weight
            )
        elif (
            result
            and left == result[-1][1]
            and weight == result[-1][2]
        ):
            result[-1] = (result[-1][0], right, weight)
        else:
            result.append((left, right, weight))
    return tuple(result)


def support_of(pieces):
    return merge_intervals(
        (left, right) for left, right, weight in pieces if weight
    )


def complement(intervals):
    intervals = merge_intervals(intervals)
    result = []
    cursor = 0
    for left, right in intervals:
        if cursor < left:
            result.append((cursor, left))
        cursor = max(cursor, right)
    if cursor < PERIOD:
        result.append((cursor, PERIOD))
    return tuple(result)


def intersect_intervals(first, second):
    first = merge_intervals(first)
    second = merge_intervals(second)
    result = []
    i = j = 0
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            result.append((left, right))
        if first[i][1] <= second[j][1]:
            i += 1
        else:
            j += 1
    return tuple(result)


def intersect_weighted(pieces, intervals):
    """Retain the weights of ``pieces`` on an unweighted interval union."""
    pieces = normalize_weighted(pieces)
    intervals = merge_intervals(intervals)
    result = []
    i = j = 0
    while i < len(pieces) and j < len(intervals):
        left = max(pieces[i][0], intervals[j][0])
        right = min(pieces[i][1], intervals[j][1])
        if left < right:
            result.append((left, right, pieces[i][2]))
        if pieces[i][1] <= intervals[j][1]:
            i += 1
        else:
            j += 1
    return normalize_weighted(result)


def intersect_weighted_sorted(pieces, intervals):
    """Fast typed scan when both inputs are already sorted and disjoint."""
    result = []
    i = j = 0
    while i < len(pieces) and j < len(intervals):
        left = max(pieces[i][0], intervals[j][0])
        right = min(pieces[i][1], intervals[j][1])
        if left < right:
            result.append((left, right, pieces[i][2]))
        if pieces[i][1] <= intervals[j][1]:
            i += 1
        else:
            j += 1
    return tuple(result)


def shift_intervals(intervals, amount):
    result = []
    for left, right in intervals:
        length = right - left
        start = (left + amount) % PERIOD
        stop = start + length
        if stop <= PERIOD:
            result.append((start, stop))
        else:
            result.extend(((start, PERIOD), (0, stop - PERIOD)))
    return merge_intervals(result)


def shift_weighted(pieces, amount):
    result = []
    for left, right, weight in pieces:
        length = right - left
        start = (left + amount) % PERIOD
        stop = start + length
        if stop <= PERIOD:
            result.append((start, stop, weight))
        else:
            result.extend(
                ((start, PERIOD, weight), (0, stop - PERIOD, weight))
            )
    return normalize_weighted(result)


def disjoint_union(first, second):
    require(
        not intersect_intervals(support_of(first), support_of(second)),
        "purported weighted union is not disjoint",
    )
    return normalize_weighted(tuple(first) + tuple(second))


def weighted_mass(pieces):
    return sum(
        (right - left) * weight for left, right, weight in pieces
    )


def build_physical_geometry():
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
    return (
        module,
        rails,
        present,
        overlap_details,
        full_module,
        e3,
        clock_combs,
        q_pairs,
    )


def build_natural_e1_profiles():
    """Reconstruct the natural sidecar instead of trusting its printed claim."""
    module, _prefixes, _, _, rails, present, _starts = (
        natural.relative.lift.m.core.build_carrier_data()
    )
    target_pullback = shift_weighted(rails[8][3], -SHIFT)
    rail_source = normalize_weighted(rails[8][3])
    rail_common = intersect_intervals(
        support_of(rail_source), support_of(target_pullback)
    )
    require(rail_common, "natural rail-eight overlap vanished")
    # The canonical overlap helper retains equal source and target rail
    # weights.  Reconstruct its two-weight object only to pass the proved
    # restriction primitive.
    rail_pairs = physical.overlap.intersect_weighted_profiles(
        rails[8][3], target_pullback
    )
    require(
        rail_pairs and all(row[2] == row[3] for row in rail_pairs),
        "natural rail weights differ on common support",
    )
    source_base = {}
    target_base = {}
    for ell in range(7):
        source_base[ell], target_base[ell] = (
            natural.clutch.restrict_to_relative_overlap(
                module, present, rail_pairs, ell
            )
        )
    e3 = natural.two.exclusive_source(module, 3)
    clock_combs = tuple(
        module.make_comb(
            module.C1, 182, 26 * ell - 13, 26 * ell + 13
        )
        for ell in range(7)
    )
    profiles = {}
    for s in COMMON_S:
        for t in COMMON_T:
            section = tuple(
                natural.two.source_present_section(
                    module, e3, 1, s, t, clock_combs
                )
            )
            profiles[s, t] = (
                intersect_weighted(source_base[1], section),
                intersect_weighted(target_base[1], section),
            )
    return profiles


def audit_all_cells():
    (
        _module,
        _rails,
        _present,
        overlap_details,
        full_module,
        e3,
        clock_combs,
        q_pairs,
    ) = build_physical_geometry()

    pulled_bases = {}
    for clock, (source_base, target_base) in enumerate(overlap_details):
        pulled = shift_weighted(target_base, -SHIFT)
        require(
            normalize_weighted(source_base) == pulled,
            f"rail base translation failed at clock {clock}",
        )
        pulled_bases[clock] = pulled

    decompositions = 0
    empty_wing_pairs = 0
    source_nonempty = 0
    target_nonempty = 0
    e1_profiles = {}
    selected_data = None

    for s in COMMON_S:
        for t in COMMON_T:
            for clock, (source_base, target_base) in enumerate(
                overlap_details
            ):
                section = tuple(
                    physical.target.source_present_section(
                        full_module, e3, clock, s, t, clock_combs
                    )
                )
                shifted_section = shift_intervals(section, -SHIFT)

                # A is the source one-sided carrier.  B is the target
                # one-sided carrier pulled back to the same source chart.
                a = intersect_weighted(source_base, section)
                b_direct = intersect_weighted(
                    pulled_bases[clock], shifted_section
                )
                target_one = intersect_weighted(target_base, section)
                b = shift_weighted(target_one, -SHIFT)
                require(
                    b == b_direct,
                    "pullback did not commute with one-sided restriction",
                )

                common_a = intersect_weighted(a, support_of(b))
                common_b = intersect_weighted(b, support_of(a))
                a_wing = intersect_weighted(a, complement(support_of(b)))
                b_wing = intersect_weighted(b, complement(support_of(a)))

                # Independent A∩B/A\B/B\A reconstruction must coincide with
                # the discovery-facing canonical constructors.
                canonical_common_source = normalize_weighted(
                    physical.common_physical_cut(
                        source_base, section, -1
                    )
                )
                canonical_common_target_pullback = shift_weighted(
                    physical.common_physical_cut(
                        target_base, section, +1
                    ),
                    -SHIFT,
                )
                canonical_source_wing = normalize_weighted(
                    physical.subtract_weighted(
                        a, canonical_common_source
                    )
                )
                canonical_target_wing_pullback = shift_weighted(
                    physical.subtract_weighted(
                        target_one,
                        physical.common_physical_cut(
                            target_base, section, +1
                        ),
                    ),
                    -SHIFT,
                )
                require(
                    common_a == canonical_common_source
                    and common_b == canonical_common_target_pullback,
                    "common carrier is not the aligned A-intersection-B",
                )
                require(
                    a_wing == canonical_source_wing
                    and b_wing == canonical_target_wing_pullback,
                    "wing is not the aligned opposite exclusive piece",
                )
                require(
                    disjoint_union(common_a, a_wing) == a
                    and disjoint_union(common_b, b_wing) == b,
                    "one-sided carrier decomposition failed",
                )
                require(
                    not intersect_intervals(
                        support_of(a_wing), support_of(b_wing)
                    ),
                    "opposite wing cofibres gained common support",
                )
                decompositions += 1
                empty_wing_pairs += 1
                source_nonempty += bool(a_wing)
                target_nonempty += bool(b_wing)

                if clock == 1:
                    e1_profiles[s, t] = (a, target_one)

                if (s, t, clock) == SELECTED:
                    selected_data = {
                        "full_module": full_module,
                        "e3": e3,
                        "clock_combs": clock_combs,
                        "q_pair": q_pairs[clock],
                        "section": section,
                        "source_one": a,
                        "target_one": target_one,
                        "source_wing": a_wing,
                        "target_wing": shift_weighted(b_wing, SHIFT),
                    }

    require(
        decompositions == empty_wing_pairs == 81 * 7,
        "all-cell wing census changed",
    )
    require(
        source_nonempty == target_nonempty == 193,
        "nonempty wing census changed",
    )
    require(selected_data is not None, "selected D cell was not built")

    natural_profiles = build_natural_e1_profiles()
    require(
        all(
            e1_profiles[key] == natural_profiles[key]
            for key in natural_profiles
        ),
        "natural fixed-e=1 sheet is not the one-sided family slice",
    )
    return (
        selected_data,
        decompositions,
        source_nonempty,
        target_nonempty,
        len(natural_profiles),
    )


def audit_selected_d(selected):
    full_module = selected["full_module"]
    s, t, clock = SELECTED
    universe = ((0, full_module.T),)
    factors = (
        selected["e3"],
        selected["clock_combs"][clock],
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
    reconstructed = universe
    for factor in factors:
        reconstructed = intersect_intervals(reconstructed, factor)
    require(
        tuple(reconstructed) == tuple(selected["section"]),
        "selected section factorization changed",
    )

    def first_failures(one, direction):
        remaining = normalize_weighted(one)
        origins = []
        for factor in factors:
            shifted = shift_intervals(factor, direction * SHIFT)
            origins.append(
                intersect_weighted(
                    remaining, complement(shifted)
                )
            )
            remaining = intersect_weighted(remaining, shifted)
        return tuple(origins), remaining

    source_origins, source_remainder = first_failures(
        selected["source_one"], -1
    )
    target_origins, target_remainder = first_failures(
        selected["target_one"], +1
    )
    require(
        source_origins[0] == selected["source_wing"]
        and target_origins[0] == selected["target_wing"]
        and all(not origin for origin in source_origins[1:])
        and all(not origin for origin in target_origins[1:]),
        "D cell is not pure at the shifted-E3 first failure",
    )
    require(
        source_remainder
        == normalize_weighted(
            physical.common_physical_cut(
                selected["source_one"], selected["section"], -1
            )
        ),
        "source first-failure remainder changed",
    )
    # The call above receives an already-sectioned source profile; imposing
    # the section once more is idempotent.  Check the target analog directly.
    require(
        target_remainder
        == normalize_weighted(
            physical.common_physical_cut(
                selected["target_one"], selected["section"], +1
            )
        ),
        "target first-failure remainder changed",
    )

    pair = selected["q_pair"]

    def coefficient(pieces, carry):
        return physical.relative.private.delayed_carry_pair(
            pieces, pair, {}
        )[carry][1]

    source_value = coefficient(selected["source_wing"], 12)
    target_value = coefficient(selected["target_wing"], 6)
    source_masses = tuple(weighted_mass(row) for row in source_origins)
    target_masses = tuple(weighted_mass(row) for row in target_origins)
    # First-failure ownership depends on the declared factor order.  Retain
    # the simultaneous failure incidences as the hostile control against
    # misreading "E3-first" as "fails E3 and no other factor".
    source_failure_masses = tuple(
        weighted_mass(
            intersect_weighted(
                selected["source_wing"],
                complement(shift_intervals(factor, -SHIFT)),
            )
        )
        for factor in factors
    )
    target_failure_masses = tuple(
        weighted_mass(
            intersect_weighted(
                selected["target_wing"],
                complement(shift_intervals(factor, +SHIFT)),
            )
        )
        for factor in factors
    )
    require(
        source_value == 2_853_968_755_527_296_447_040
        and target_value == 5_707_937_511_054_592_894_080
        and target_value == 2 * source_value,
        "D-cell gain-two coefficient changed",
    )
    require(
        source_masses
        == (2_188_067_024_426_754_240, 0, 0, 0, 0, 0)
        and target_masses
        == (1_458_711_349_617_836_160, 0, 0, 0, 0, 0),
        "D-cell first-failure masses changed",
    )
    require(
        len(selected["source_wing"]) == 3
        and len(selected["target_wing"]) == 2
        and shift_weighted(selected["source_wing"], SHIFT)
        != selected["target_wing"]
        and weighted_mass(selected["source_wing"])
        != weighted_mass(selected["target_wing"]),
        "D-cell geometric hostile changed",
    )
    return (
        source_value,
        target_value,
        source_masses,
        target_masses,
        source_failure_masses,
        target_failure_masses,
    )


def weighted_endpoint_sum(pieces, frequency, embedding):
    prime = embedding[0]
    total = 0
    for left, right, weight in pieces:
        value = endpoint_base.endpoint.endpoint_sum(
            ((left, right),), frequency, (embedding,)
        )[0]
        total = (total + weight * value) % prime
    return total


def weighted_x_sweep(
    pieces, terminal, q_starts, frequency, embedding, tabs
):
    prime = embedding[0]
    total = 0
    for left, right, weight in pieces:
        values, _overlap = endpoint_base.endpoint.x_sweep(
            ((left, right),),
            terminal,
            q_starts,
            frequency,
            (embedding,),
            tabs,
        )
        total = (total + weight * values[0]) % prime
    return total


_ENDPOINT_PRESENT_SETS = None
_ENDPOINT_SHIFTED_PRESENT_SETS = None


def endpoint_present_sets():
    global _ENDPOINT_PRESENT_SETS
    if _ENDPOINT_PRESENT_SETS is None:
        _ENDPOINT_PRESENT_SETS = {
            address: tuple(
                endpoint_base.endpoint.build_set(
                    endpoint_base.PAT_E3, ell
                )
            )
            for address, ell in endpoint_base.REPS.items()
        }
    return _ENDPOINT_PRESENT_SETS


def endpoint_shifted_present_sets(unit):
    global _ENDPOINT_SHIFTED_PRESENT_SETS
    if _ENDPOINT_SHIFTED_PRESENT_SETS is None:
        shifted = {}
        for address, ell in endpoint_base.REPS.items():
            shifted_ell = tuple(
                (ell[index] + endpoint_base.WMOD[index]) % P
                for index in range(len(ell))
            )
            shifted[address] = tuple(
                endpoint_base.endpoint.build_set(
                    endpoint_base.PAT_E3, shifted_ell
                )
            )
            require(
                merge_intervals(shifted[address])
                == shift_intervals(
                    endpoint_present_sets()[address], -unit
                ),
                "ell-to-ell+W present-set covariance failed",
            )
        _ENDPOINT_SHIFTED_PRESENT_SETS = shifted
    return _ENDPOINT_SHIFTED_PRESENT_SETS


def endpoint_field_audit(source_wing, target_wing, embedding, gauge=False):
    terminal = tuple(
        endpoint_base.endpoint.build_set(
            endpoint_base.PAT_Q12, endpoint_base.ZERO
        )
    )
    q_starts = tuple(left for left, _right in terminal)
    prime, root = embedding
    require(
        pow(root, endpoint_base.NRED, prime) == 1
        and all(
            pow(
                root,
                endpoint_base.NRED // factor,
                prime,
            )
            != 1
            for factor in endpoint_base.endpoint.NN_PRIMES
        ),
        "embedding root lost exact cyclotomic order",
    )
    tabs = endpoint_base.endpoint.make_tabs(
        terminal, endpoint_base.X0, (embedding,)
    )
    zeta = pow(root, endpoint_base.NRED // P, prime)
    powers = tuple(pow(zeta, exponent, prime) for exponent in range(P))
    unit = endpoint_base.T // P
    source_shifts = tuple(
        shift_weighted(source_wing, -harmonic * unit)
        for harmonic in range(P)
    )
    target_shifts = tuple(
        shift_weighted(target_wing, harmonic * unit)
        for harmonic in range(P)
    )
    present_sets = endpoint_present_sets()

    def evaluate_left(ell, present, carrier):
        restricted = intersect_weighted_sorted(carrier, present)
        value = weighted_x_sweep(
            restricted,
            terminal,
            q_starts,
            endpoint_base.X0,
            embedding,
            tabs,
        )
        return (
            value
            * pow(zeta, ell[endpoint_base.DEEP], prime)
            % prime
        )

    def evaluate_right(_ell, present, carrier):
        restricted = intersect_weighted_sorted(carrier, present)
        return weighted_endpoint_sum(
            restricted, -endpoint_base.Y0, embedding
        )

    left_dual = {}
    right_dual = {}
    for address, ell in endpoint_base.REPS.items():
        present = present_sets[address]
        for harmonic in range(P):
            left_dual[address, harmonic] = evaluate_left(
                ell, present, source_shifts[harmonic]
            )
            right_dual[address, harmonic] = evaluate_right(
                ell, present, target_shifts[harmonic]
            )

    gauge_checks = 0
    hostile_left = None
    hostile_right = None
    if gauge:
        shifted_present_sets = endpoint_shifted_present_sets(unit)
        for address, ell in endpoint_base.REPS.items():
            shifted_ell = tuple(
                (ell[index] + endpoint_base.WMOD[index]) % P
                for index in range(len(ell))
            )
            shifted_present = shifted_present_sets[address]
            for harmonic in range(P):
                left_shifted = evaluate_left(
                    shifted_ell,
                    shifted_present,
                    shift_weighted(
                        source_wing, -(harmonic + 1) * unit
                    ),
                )
                right_shifted = evaluate_right(
                    shifted_ell,
                    shifted_present,
                    shift_weighted(
                        target_wing, (harmonic - 1) * unit
                    ),
                )
                require(
                    left_shifted == left_dual[address, harmonic],
                    "exhaustive left representative gauge failed",
                )
                require(
                    right_shifted == right_dual[address, harmonic],
                    "exhaustive right representative gauge failed",
                )
                gauge_checks += 1

                if hostile_left is None:
                    fixed_left = evaluate_left(
                        shifted_ell,
                        shifted_present,
                        source_shifts[harmonic],
                    )
                    if fixed_left != left_dual[address, harmonic]:
                        hostile_left = (
                            address,
                            harmonic,
                            left_dual[address, harmonic],
                            fixed_left,
                        )
                if hostile_right is None:
                    fixed_right = evaluate_right(
                        shifted_ell,
                        shifted_present,
                        target_shifts[harmonic],
                    )
                    if fixed_right != right_dual[address, harmonic]:
                        hostile_right = (
                            address,
                            harmonic,
                            right_dual[address, harmonic],
                            fixed_right,
                        )
        require(
            gauge_checks == P**3
            and hostile_left is not None
            and hostile_right is not None,
            "gauge census or fixed-carrier hostile failed",
        )

    # Independent loop ordering: collapse the carrier harmonic first at
    # k=l=0, then take the two-dimensional endpoint transforms.
    collapsed_left = {
        address: sum(
            left_dual[address, harmonic]
            for harmonic in range(P)
        )
        % prime
        for address in endpoint_base.KEYS
    }
    collapsed_right = {
        address: sum(
            right_dual[address, harmonic]
            for harmonic in range(P)
        )
        % prime
        for address in endpoint_base.KEYS
    }

    def left_primal(point):
        return sum(
            collapsed_left[address]
            * powers[
                -(
                    address[0] * point[0]
                    + address[1] * point[1]
                )
                % P
            ]
            for address in endpoint_base.KEYS
        ) % prime

    def right_primal(point):
        return sum(
            collapsed_right[address]
            * powers[
                (
                    address[0] * point[0]
                    + address[1] * point[1]
                )
                % P
            ]
            for address in endpoint_base.KEYS
        ) % prime

    edges = []
    sector = 0
    q0, q1 = Q_TARGET
    for right0 in range(P):
        require(q0 != 0, "selected sector solver expects q0 nonzero")
        right1 = (
            (DELTA_TARGET + q1 * right0)
            * pow(q0, -1, P)
        ) % P
        right = (right0, right1)
        left = ((right0 + q0) % P, (right1 + q1) % P)
        left_value = left_primal(left)
        right_value = right_primal(right)
        coefficient = left_value * right_value % prime
        require(
            (
                (left[0] - right[0]) % P,
                (left[1] - right[1]) % P,
            )
            == Q_TARGET
            and (
                left[0] * right[1] - left[1] * right[0]
            )
            % P
            == DELTA_TARGET,
            "endpoint edge escaped the typed determinant sector",
        )
        edges.append(
            (left, right, left_value, right_value, coefficient)
        )
        sector = (sector + coefficient) % prime

    support = (
        sum(bool(value) for value in left_dual.values()),
        sum(bool(value) for value in right_dual.values()),
    )
    require(
        len(edges) == P
        and all(row[-1] for row in edges)
        and sector != 0,
        "selected determinant sector vanished",
    )
    return {
        "prime": prime,
        "support": support,
        "sector": sector,
        "edge_support": sum(bool(row[-1]) for row in edges),
        "gauge_checks": gauge_checks,
        "hostile_left": hostile_left,
        "hostile_right": hostile_right,
    }


def main():
    (
        selected,
        decompositions,
        source_nonempty,
        target_nonempty,
        natural_matches,
    ) = audit_all_cells()
    (
        source_value,
        target_value,
        source_masses,
        target_masses,
        source_failure_masses,
        target_failure_masses,
    ) = audit_selected_d(selected)

    field_rows = tuple(
        endpoint_field_audit(
            selected["source_wing"],
            selected["target_wing"],
            embedding,
            gauge=(index == 0),
        )
        for index, embedding in enumerate(endpoint_base.endpoint.MODS)
    )
    require(
        field_rows[0]["support"] == (641, 510),
        "first-field dual support census changed",
    )
    require(
        field_rows[0]["sector"] == 64_970_359_711_488_738,
        "first-field sector value changed",
    )

    source = Path(__file__).read_text()
    tree = ast.parse(source)
    require(
        not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
        "optimized-mode validity gate contains assert",
    )

    print("THM-2782 OPPOSITE-WING COMMON-ATOM / VIRTUAL-SECTOR AUDIT")
    print(f"dependency_sha256={tuple(PINNED.items())}")
    print(
        "aligned_object_typing="
        "A=source_one_sided;"
        "B=pullback(target_one_sided);"
        "common=A_intersect_B;"
        "wings=(A_minus_B,B_minus_A)"
    )
    print(
        f"all_81x7_cells=decompositions:{decompositions}/567;"
        f"opposite_wing_intersections_empty:{decompositions}/567"
    )
    print(
        f"nonempty_wing_supports="
        f"(source:{source_nonempty}/567,target:{target_nonempty}/567)"
    )
    print(
        f"natural_e1_one_sided_exact_matches={natural_matches}/81;"
        "gain11_is_inherited_from_the_separate_canonical_sidecar"
    )
    print(
        f"selected_D={SELECTED};"
        "first_failure_order="
        "(shifted_E3,shifted_clock,shifted_q1,shifted_q2,"
        "shifted_c2,shifted_c3)"
    )
    print(
        f"selected_D_coefficients=(source:{source_value},"
        f"target:{target_value},gain:{target_value // source_value})"
    )
    print(
        f"selected_D_first_failure_masses="
        f"(source:{source_masses},target:{target_masses})"
    )
    print(
        f"selected_D_simultaneous_failure_masses="
        f"(source:{source_failure_masses},"
        f"target:{target_failure_masses})"
    )
    print(
        f"selected_D_geometry=(pieces:"
        f"{len(selected['source_wing'])}/{len(selected['target_wing'])},"
        f"masses:{weighted_mass(selected['source_wing'])}/"
        f"{weighted_mass(selected['target_wing'])},"
        "chart_translate:False)"
    )
    for index, row in enumerate(field_rows, 1):
        print(
            f"field_{index}=prime:{row['prime']};"
            f"dual_support:{row['support']};"
            f"sector_q{Q_TARGET}_Delta{DELTA_TARGET}_k{K_TARGET}_l{L_TARGET}:"
            f"{row['sector']};edge_support:{row['edge_support']}/13"
        )
    first = field_rows[0]
    print(
        f"representative_gauge_exhaustive="
        f"left_and_right:{first['gauge_checks']}/2197"
    )
    print(
        f"fixed_carrier_gauge_hostiles="
        f"left:{first['hostile_left']};right:{first['hostile_right']}"
    )
    print(
        "CHARACTERISTIC_ZERO_CERTIFICATE="
        "NONZERO_IN_TWO_CERTIFIED_EXACT_ORDER_SPECIALIZATIONS"
    )
    print(
        "STRONGEST_CONCLUSION="
        "the_opposite_physical_wing_cofibres_cannot_supply_one_common_"
        "ancestry_(1,1)_allocation_atom;the_independent_endpoint_product_"
        "is_only_a_virtual_pullback_sector"
    )
    print(
        "OPEN_SIDECAR="
        "one_atom_dependent_one-sided/common_four-corner_lift_with_"
        "allocation_factor_origins_endpoint_origin_clock_semantic_owner_"
        "and_root-deck_intertwiner"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
