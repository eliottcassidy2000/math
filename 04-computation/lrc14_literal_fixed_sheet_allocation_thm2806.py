#!/usr/bin/env python3
"""Exact same-sheet common-carrier endpoint-square audit.

This companion keeps one literal Boolean THM-2584 ancestry triple through
the rail-eight common physical carrier, all factor masks, the separate
carrier-harmonic/endpoint transforms, and only then forms four endpoint
products.  It distinguishes four statements:

* after summing the two carrier-twist variables, the fixed-sheet central
  bare/present coefficient square has nonzero Segre--Hadamard coordinate D3;
* before that carrier DFT, the sole fourfold-supported twist point has raw
  Mobius face zero, while the nonzero central D3 is the sum over the 12-by-12
  complement where only the fully bare state survives;
* translating the sheet through its carrier orbit gives a separate
  four-sample endpoint parallelogram;
* the endpoint determinant labels on the chosen parallelogram have mixed
  face det(s,t)=1.

The last two facts do not identify the endpoint translations with the
bare/present carrier allocation toggles.  The literal present state keeps
the outer sheet selector fixed and multiplies by each carrier twist; only
twist zero survives.  Thus its intrinsic endpoint action is scalar with
zero step.  The script separately searches whether the moving-orbit array
can nevertheless be a translated/gauged bare array and finds an exact no-go.
"""

from __future__ import annotations

from bisect import bisect_right
from hashlib import sha256
from pathlib import Path
from struct import pack
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

import lrc14_extended_carrier_endpoint_lib as endpoint_base
import lrc14_root_zero_full_target_semantic_clutch_20260728 as physical


P = 13
CELL = (0, 4, 1)
COMMON_S = (0, 1, 2, 3, 8, 9, 10, 11, 12)
COMMON_T = tuple(range(3, 12))
ATOM_INTERVAL = (142004992589460, 142005019034340)
ATOM_WEIGHT = 27581135604
ANCESTRY_LABEL = (59162, 26, 56658)
ANCESTRY_ALLOWED = (138281416853580, 159555049051860)
SEMANTIC_VALUE = 103478815440
SOURCE_STATE = (7, 1, 12, 12)
TARGET_STATE = (8, 0, 6, 1)
DELAYED_HALF = 1
LEFT_ORIGIN = (0, 0)
RIGHT_ORIGIN = (0, 0)
SOURCE_STEP = (0, 1)
TARGET_STEP = (12, 0)

EXPECTED_FIELDS = {
    352341050142921841: {
        "factors": (
            345283638862540358,
            180096700909697779,
            167108689435423878,
            144423867160388279,
        ),
        "products": (
            347197813355232446,
            3438130309897238,
            320064096166488418,
            146419401152414063,
        ),
        "hadamard": (
            112437340698188483,
            236493496489149044,
            165063327916487722,
            170114988031260853,
        ),
    },
    956354278959359281: {
        "factors": (
            77115584218272234,
            620226426191367084,
            950936471790018101,
            445878683356068526,
        ),
        "products": (
            788257626617883404,
            251268809443191686,
            161379668316310808,
            381695617330433150,
        ),
        "hadamard": (
            626247442748459767,
            496451150414331132,
            316672868160569376,
            757304766188814060,
        ),
    },
}


def weighted_intersection(pieces, intervals):
    return tuple(
        physical.relative.private.old.intersect_weighted_union(
            pieces, intervals
        )
    )


def indexed_weighted_intersection(pieces, intervals, starts):
    """Restrict a tiny carrier without rebuilding a huge start array."""
    out = []
    for left, right, weight in pieces:
        index = max(0, bisect_right(starts, left) - 1)
        while index < len(intervals) and intervals[index][0] < right:
            a = max(left, intervals[index][0])
            b = min(right, intervals[index][1])
            if a < b:
                out.append((a, b, weight))
            index += 1
    return tuple(out)


def support_of(pieces):
    return physical.overlap.merge_intervals(
        (left, right) for left, right, weight in pieces if weight
    )


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
    return tuple(out)


def contained(interval, intervals):
    return intersect_sorted((interval,), intervals) == (interval,)


def add(point, step):
    return tuple((point[index] + step[index]) % P for index in range(2))


def det(left, right):
    return (left[0] * right[1] - left[1] * right[0]) % P


def bank_digest(bank, keys):
    return sha256(
        ",".join(str(bank[key]) for key in keys).encode()
    ).hexdigest()


def build_geometry():
    (
        module, _prefixes, _whole, _masses, rails, present, _starts
    ) = physical.relative.lift.m.core.build_carrier_data()
    pair_prefixes = physical.relative.private.build_pair_prefixes(module)
    _sv, _tv, rail_pairs, details = physical.overlap.overlap_vectors(
        module, pair_prefixes, rails, present, rail_index=8
    )
    full_module = physical.target.load_present_module()
    e3 = physical.target.exclusive_source(full_module, 3)
    fork = physical.target.deepest_fork(full_module)
    clocks = tuple(
        full_module.make_comb(
            full_module.C1, 182, 26 * clock - 13, 26 * clock + 13
        )
        for clock in range(7)
    )
    q_pairs = physical.q_restricted_pair_prefixes(
        full_module, pair_prefixes, fork
    )
    require(
        tuple(rails[8][:3]) == (1, 4, 12),
        "rail-eight physical metadata changed",
    )
    require(
        all(source_weight == target_weight
            for _a, _b, source_weight, target_weight in rail_pairs),
        "rail-eight aligned profile lost equal weights",
    )
    return module, full_module, details, e3, clocks, q_pairs


def section_factors(full_module, e3, clocks, s, t, clock):
    universe = ((0, full_module.T),)
    return (
        e3,
        clocks[clock],
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


def common_cell(details, full_module, e3, clocks, s, t, clock):
    source_base, target_base = details[clock]
    section = physical.target.source_present_section(
        full_module, e3, clock, s, t, clocks
    )
    source = weighted_intersection(source_base, section)
    target = weighted_intersection(target_base, section)
    target_pullback = physical.overlap.shift_weighted(
        target, -physical.SHIFT
    )
    aligned = physical.overlap.intersect_weighted_profiles(
        source, target_pullback
    )
    require(
        all(source_weight == target_weight
            for _a, _b, source_weight, target_weight in aligned),
        f"aligned common weights differ at {(s, t, clock)}",
    )
    common = tuple(
        (left, right, source_weight)
        for left, right, source_weight, _target_weight in aligned
    )
    canonical_source = physical.common_physical_cut(
        source_base, section, -1
    )
    canonical_target = physical.common_physical_cut(
        target_base, section, +1
    )
    require(common == canonical_source,
            f"source common constructor mismatch at {(s, t, clock)}")
    require(
        physical.overlap.shift_weighted(
            common, physical.SHIFT
        ) == canonical_target,
        f"target common constructor mismatch at {(s, t, clock)}",
    )
    source_wing = physical.subtract_weighted(source, common)
    target_wing_pullback = physical.subtract_weighted(
        target_pullback, common
    )
    require(
        not intersect_sorted(
            support_of(source_wing), support_of(target_wing_pullback)
        ),
        f"opposite wings intersect at {(s, t, clock)}",
    )
    return (
        source, target, target_pullback, common,
        source_wing, target_wing_pullback, section,
    )


def exhaustive_wing_control(details, full_module, e3, clocks):
    source_nonempty = target_nonempty = 0
    common_nonempty_by_clock = [0] * 7
    selected = None
    cells = 0
    for s in COMMON_S:
        for t in COMMON_T:
            for clock in range(7):
                row = common_cell(
                    details, full_module, e3, clocks, s, t, clock
                )
                source_nonempty += bool(row[4])
                target_nonempty += bool(row[5])
                common_nonempty_by_clock[clock] += bool(row[3])
                cells += 1
                if (s, t, clock) == CELL:
                    selected = row
    require(
        (cells, source_nonempty, target_nonempty) == (567, 193, 193),
        "all-567 opposite-wing census changed",
    )
    require(
        tuple(common_nonempty_by_clock) == (0, 81, 56, 56, 0, 0, 0),
        "per-clock common-cell census changed",
    )
    require(selected is not None, "selected common cell was not visited")
    return (
        selected,
        (cells, source_nonempty, target_nonempty),
        tuple(common_nonempty_by_clock),
    )


def selected_atom_audit(selected, full_module, e3, clocks, q_pairs):
    (
        _source, _target, _target_pullback, common,
        source_wing, target_wing_pullback, _section,
    ) = selected
    expected_piece = (*ATOM_INTERVAL, ATOM_WEIGHT)
    require(
        expected_piece in common,
        "selected weighted common piece disappeared",
    )
    source_atom = ((*ATOM_INTERVAL, 1),)
    target_atom = physical.overlap.shift_weighted(
        source_atom, physical.SHIFT
    )
    require(
        contained(ATOM_INTERVAL, support_of(selected[0]))
        and contained(target_atom[0][:2], support_of(selected[1])),
        "atom selector is not contained in the zero-twist carrier",
    )
    require(
        not intersect_sorted(support_of(source_atom), support_of(source_wing))
        and not intersect_sorted(
            support_of(source_atom), support_of(target_wing_pullback)
        ),
        "selected common atom was reselected from an exclusive wing",
    )

    s, t, clock = CELL
    factors = section_factors(
        full_module, e3, clocks, s, t, clock
    )
    source_interval = ATOM_INTERVAL
    target_interval = tuple(
        endpoint + physical.SHIFT for endpoint in ATOM_INTERVAL
    )
    source_native = tuple(
        contained(source_interval, factor) for factor in factors
    )
    source_adjacent = tuple(
        contained(
            source_interval,
            physical.overlap.shift_union(factor, -physical.SHIFT),
        )
        for factor in factors
    )
    target_native = tuple(
        contained(target_interval, factor) for factor in factors
    )
    target_adjacent = tuple(
        contained(
            target_interval,
            physical.overlap.shift_union(factor, physical.SHIFT),
        )
        for factor in factors
    )
    require(
        all(source_native + source_adjacent
            + target_native + target_adjacent),
        "selected atom lost a simultaneous factor",
    )

    source_semantic = physical.relative.private.delayed_carry_pair(
        source_atom, q_pairs[clock], {}
    )
    target_semantic = physical.relative.private.delayed_carry_pair(
        target_atom, q_pairs[clock], {}
    )
    require(
        source_semantic[12] == (0, SEMANTIC_VALUE)
        and target_semantic[6] == (0, SEMANTIC_VALUE),
        "selected unit sheet lost its delayed semantic coefficient",
    )
    unit = full_module.T // P

    def allocation_mask(atom_interval, carrier, sign):
        result = []
        for harmonic in range(P):
            translated = physical.overlap.shift_weighted(
                carrier, sign * harmonic * unit
            )
            result.append(
                intersect_sorted(
                    (atom_interval,), support_of(translated)
                )
            )
        return tuple(result)

    source_allocation_mask = allocation_mask(
        ATOM_INTERVAL, selected[0], -1
    )
    target_allocation_mask = allocation_mask(
        target_atom[0][:2], selected[1], +1
    )
    require(
        source_allocation_mask
        == ((ATOM_INTERVAL,),) + ((),) * (P - 1)
        and target_allocation_mask
        == ((target_atom[0][:2],),) + ((),) * (P - 1),
        "fixed atom meets a nonzero carrier twist",
    )
    return (
        source_atom,
        target_atom,
        (source_native, source_adjacent, target_native, target_adjacent),
        (source_semantic[12], target_semantic[6]),
        (source_allocation_mask, target_allocation_mask),
    )


def interval_hit(numerator, denominator, intervals, starts):
    floor_value = numerator // denominator
    index = bisect_right(starts, floor_value) - 1
    if (
        index >= 0
        and intervals[index][0] * denominator <= numerator
        < intervals[index][1] * denominator
    ):
        return index
    return None


def ancestry_sheet_audit():
    core = physical.relative.lift.m.core
    old = core.old
    e_intervals = tuple(
        old.base.build_set(old.base.PAT_E, old.base.ZELL)
    )
    q_intervals = tuple(
        old.base.build_set(old.host.PAT_QB, old.base.ZELL)
    )
    e_starts = tuple(left for left, _right in e_intervals)
    q_starts = tuple(left for left, _right in q_intervals)
    period = old.T
    depth = old.rail.DEPTH
    packet_clock = old.rail.RPACK
    require(
        (period, depth, packet_clock)
        == (297836897838480, 13**5, 13**2),
        "Boolean ancestry scales changed",
    )
    midpoint = sum(ATOM_INTERVAL) // 2
    target_midpoint = midpoint + physical.SHIFT
    displacement = period // P

    source_u_count = target_u_count = u_mismatches = 0
    source_u_digest = sha256()
    target_u_digest = sha256()
    first_common_u = None
    for ancestry in range(depth):
        source_x = midpoint + ancestry * period
        target_x = target_midpoint + ancestry * period
        source_q = interval_hit(
            source_x, depth, q_intervals, q_starts
        ) is not None
        target_q = interval_hit(
            target_x, depth, q_intervals, q_starts
        ) is not None
        if not (source_q or target_q):
            continue
        for packet in range(packet_clock):
            source_active = (
                source_q
                and interval_hit(
                    source_x + packet * depth * period,
                    depth * packet_clock,
                    e_intervals,
                    e_starts,
                ) is not None
            )
            target_active = (
                target_q
                and interval_hit(
                    target_x + packet * depth * period,
                    depth * packet_clock,
                    e_intervals,
                    e_starts,
                ) is not None
            )
            encoded = pack(">II", ancestry, packet)
            if source_active:
                source_u_count += 1
                source_u_digest.update(encoded)
            if target_active:
                target_u_count += 1
                target_u_digest.update(encoded)
            if source_active != target_active:
                u_mismatches += 1
            if source_active and target_active and first_common_u is None:
                first_common_u = (ancestry, packet)

    source_v_count = target_v_count = v_mismatches = 0
    source_v_digest = sha256()
    target_v_digest = sha256()
    first_common_v = None
    for source_sheet in range(depth):
        source_value = midpoint - displacement + source_sheet * period
        target_value = (
            target_midpoint - displacement + source_sheet * period
        )
        source_active = interval_hit(
            source_value, depth, e_intervals, e_starts
        ) is not None
        target_active = interval_hit(
            target_value, depth, e_intervals, e_starts
        ) is not None
        encoded = pack(">I", source_sheet)
        if source_active:
            source_v_count += 1
            source_v_digest.update(encoded)
        if target_active:
            target_v_count += 1
            target_v_digest.update(encoded)
        if source_active != target_active:
            v_mismatches += 1
        if source_active and target_active and first_common_v is None:
            first_common_v = source_sheet

    require(
        (
            source_u_count, target_u_count, u_mismatches,
            source_v_count, target_v_count, v_mismatches,
        ) == (966606, 966606, 0, 28534, 28534, 0),
        "source/target Boolean ancestry label sets differ",
    )
    require(
        source_u_digest.digest() == target_u_digest.digest()
        and source_v_digest.digest() == target_v_digest.digest(),
        "Boolean ancestry label-set digests differ",
    )
    require(
        source_u_count * source_v_count == ATOM_WEIGHT,
        "Boolean ancestry factorization does not recover rail weight",
    )
    require(
        (*first_common_u, first_common_v) == ANCESTRY_LABEL,
        "first common Boolean ancestry label changed",
    )

    ancestry, packet, source_sheet = ANCESTRY_LABEL

    def containing(numerator, denominator, intervals, starts):
        index = interval_hit(
            numerator, denominator, intervals, starts
        )
        require(index is not None, "declared ancestry sheet is inactive")
        return intervals[index]

    constraints = []
    for label, offset in (
        ("source", 0), ("target", physical.SHIFT)
    ):
        q_interval = containing(
            midpoint + offset + ancestry * period,
            depth, q_intervals, q_starts,
        )
        e_u_interval = containing(
            midpoint + offset + ancestry * period
            + packet * depth * period,
            depth * packet_clock, e_intervals, e_starts,
        )
        e_v_interval = containing(
            midpoint + offset - displacement + source_sheet * period,
            depth, e_intervals, e_starts,
        )
        constraints.extend((
            (
                f"{label}_U_Q",
                q_interval[0] * depth - ancestry * period - offset,
                q_interval[1] * depth - ancestry * period - offset,
            ),
            (
                f"{label}_U_E",
                e_u_interval[0] * depth * packet_clock
                - ancestry * period - packet * depth * period - offset,
                e_u_interval[1] * depth * packet_clock
                - ancestry * period - packet * depth * period - offset,
            ),
            (
                f"{label}_V_E",
                e_v_interval[0] * depth + displacement
                - source_sheet * period - offset,
                e_v_interval[1] * depth + displacement
                - source_sheet * period - offset,
            ),
        ))
    allowed = (
        max(left for _label, left, _right in constraints),
        min(right for _label, _left, right in constraints),
    )
    require(
        allowed == ANCESTRY_ALLOWED
        and allowed[0] <= ATOM_INTERVAL[0]
        < ATOM_INTERVAL[1] <= allowed[1],
        "concrete ancestry label does not cover the whole physical atom",
    )
    return {
        "u_counts": (source_u_count, target_u_count),
        "v_counts": (source_v_count, target_v_count),
        "u_digest": source_u_digest.hexdigest(),
        "v_digest": source_v_digest.hexdigest(),
        "constraints": tuple(constraints),
        "allowed": allowed,
    }


def endpoint_sum(pieces, frequency, embedding):
    prime = embedding[0]
    total = 0
    for left, right, weight in pieces:
        value = endpoint_base.endpoint.endpoint_sum(
            ((left, right),), frequency, (embedding,)
        )[0]
        total = (total + weight * value) % prime
    return total


def x_sweep(pieces, terminal, q_starts, embedding, tabs):
    prime = embedding[0]
    total = 0
    overlap = 0
    for left, right, weight in pieces:
        value, local_overlap = endpoint_base.endpoint.x_sweep(
            ((left, right),), terminal, q_starts, endpoint_base.X0,
            (embedding,), tabs,
        )
        total = (total + weight * value[0]) % prime
        overlap += weight * local_overlap
    return total, overlap


def translation_covariances(bare, carried, zeta, prime):
    """Find carried(x)=c*zeta^(u.x)*bare(x+s) on the whole endpoint plane."""
    answers = []
    points = endpoint_base.KEYS
    for shift in points:
        for character in points:
            scalar = None
            valid = True
            for point in points:
                shifted = add(point, shift)
                base_value = bare[shifted]
                phase = pow(
                    zeta,
                    (character[0] * point[0]
                     + character[1] * point[1]) % P,
                    prime,
                )
                denominator = base_value * phase % prime
                value = carried[point]
                if denominator == 0:
                    if value:
                        valid = False
                        break
                    continue
                ratio = value * pow(denominator, -1, prime) % prime
                if scalar is None:
                    scalar = ratio
                elif ratio != scalar:
                    valid = False
                    break
            if valid:
                answers.append((shift, character, scalar or 0))
    return tuple(answers)


def endpoint_field_audit(
    source_atom, target_atom, present_sets, embedding
):
    terminal = tuple(
        endpoint_base.endpoint.build_set(
            endpoint_base.PAT_Q12, endpoint_base.ZERO
        )
    )
    q_starts = tuple(left for left, _right in terminal)
    tabs = endpoint_base.endpoint.make_tabs(
        terminal, endpoint_base.X0, (embedding,)
    )
    prime, root = embedding
    zeta = pow(root, endpoint_base.NRED // P, prime)
    powers = tuple(pow(zeta, exponent, prime) for exponent in range(P))
    unit = endpoint_base.T // P
    source_shifts = tuple(
        physical.overlap.shift_weighted(source_atom, -harmonic * unit)
        for harmonic in range(P)
    )
    target_shifts = tuple(
        physical.overlap.shift_weighted(target_atom, harmonic * unit)
        for harmonic in range(P)
    )
    left_cache = {}
    right_cache = {}

    def left_value(ell, present, starts, carrier):
        restricted = indexed_weighted_intersection(
            carrier, present, starts
        )
        if restricted not in left_cache:
            left_cache[restricted] = x_sweep(
                restricted, terminal, q_starts, embedding, tabs
            )[0]
        phase = pow(zeta, ell[endpoint_base.DEEP], prime)
        return left_cache[restricted] * phase % prime

    def right_value(_ell, present, starts, carrier):
        restricted = indexed_weighted_intersection(
            carrier, present, starts
        )
        if restricted not in right_cache:
            right_cache[restricted] = endpoint_sum(
                restricted, -endpoint_base.Y0, embedding
            )
        return right_cache[restricted]

    left_dual = {}
    right_dual = {}
    gauge_left_failures = gauge_right_failures = 0
    fixed_sheet_left_gauge_rows = fixed_sheet_right_gauge_rows = 0
    fixed_left_mismatches = fixed_right_mismatches = 0
    first_fixed_left = first_fixed_right = None
    for address, ell in endpoint_base.REPS.items():
        present = present_sets[address]
        present_starts = tuple(left for left, _right in present)
        shifted_ell = tuple(
            (ell[index] + endpoint_base.WMOD[index]) % P
            for index in range(9)
        )
        shifted_present = physical.overlap.shift_union(
            present, -unit
        )
        shifted_starts = tuple(
            left for left, _right in shifted_present
        )
        for harmonic in range(P):
            key = (address, harmonic)
            left_dual[key] = left_value(
                ell, present, present_starts, source_shifts[harmonic]
            )
            right_dual[key] = right_value(
                ell, present, present_starts, target_shifts[harmonic]
            )
            shifted_left = left_value(
                shifted_ell, shifted_present, shifted_starts,
                source_shifts[(harmonic + 1) % P],
            )
            shifted_right = right_value(
                shifted_ell, shifted_present, shifted_starts,
                target_shifts[(harmonic - 1) % P],
            )
            gauge_left_failures += (
                shifted_left != left_dual[key]
            )
            gauge_right_failures += (
                shifted_right != right_dual[key]
            )
            if harmonic == 0:
                fixed_sheet_left_gauge_rows += (
                    shifted_left == left_dual[key]
                )
                fixed_sheet_right_gauge_rows += (
                    shifted_right == right_dual[key]
                )
            fixed_left = left_value(
                shifted_ell, shifted_present, shifted_starts,
                source_shifts[harmonic],
            )
            fixed_right = right_value(
                shifted_ell, shifted_present, shifted_starts,
                target_shifts[harmonic],
            )
            if fixed_left != left_dual[key]:
                fixed_left_mismatches += 1
                if first_fixed_left is None:
                    first_fixed_left = (
                        key, left_dual[key], fixed_left
                    )
            if fixed_right != right_dual[key]:
                fixed_right_mismatches += 1
                if first_fixed_right is None:
                    first_fixed_right = (
                        key, right_dual[key], fixed_right
                    )
    require(
        gauge_left_failures == gauge_right_failures == 0,
        "representative gauge failed",
    )
    require(
        fixed_sheet_left_gauge_rows == fixed_sheet_right_gauge_rows
        == P**2,
        "fixed-sheet h=0 representative gauge failed",
    )
    require(
        fixed_left_mismatches and fixed_right_mismatches,
        "fixed-carrier representative hostile vanished",
    )
    require(
        sum(bool(value) for value in left_dual.values()) == 522
        and sum(bool(value) for value in right_dual.values()) == 522,
        "atom-specific dual support census changed",
    )

    def inverse_bank(dual, sign, mode):
        result = {}
        for point in endpoint_base.KEYS:
            total = 0
            for address in endpoint_base.KEYS:
                pairing = (
                    address[0] * point[0]
                    + address[1] * point[1]
                ) % P
                if mode == "bare":
                    total += (
                        P * dual[address, 0]
                        * powers[(sign * pairing) % P]
                    )
                elif mode == "fixed_present":
                    total += (
                        dual[address, 0]
                        * powers[(sign * pairing) % P]
                    )
                elif mode == "orbit_present":
                    total += sum(
                        dual[address, harmonic]
                        * powers[(sign * pairing) % P]
                        for harmonic in range(P)
                    )
                else:
                    raise RuntimeError("unknown inverse-bank mode")
            result[point] = total % prime
        return result

    left_orbit = inverse_bank(left_dual, -1, "orbit_present")
    right_orbit = inverse_bank(right_dual, +1, "orbit_present")
    left_present = inverse_bank(left_dual, -1, "fixed_present")
    right_present = inverse_bank(right_dual, +1, "fixed_present")
    left_bare = inverse_bank(left_dual, -1, "bare")
    right_bare = inverse_bank(right_dual, +1, "bare")
    require(
        all(left_orbit.values()) and all(right_orbit.values())
        and all(left_present.values()) and all(right_present.values())
        and all(left_bare.values()) and all(right_bare.values()),
        "endpoint plane lost support",
    )
    require(
        all(
            left_bare[point] == P * left_present[point] % prime
            and right_bare[point] == P * right_present[point] % prime
            for point in endpoint_base.KEYS
        ),
        "fixed-sheet bare/present scalar law failed",
    )

    left_origin_1 = add(LEFT_ORIGIN, SOURCE_STEP)
    right_origin_1 = add(RIGHT_ORIGIN, TARGET_STEP)
    factors = (
        left_orbit[LEFT_ORIGIN],
        left_orbit[left_origin_1],
        right_orbit[RIGHT_ORIGIN],
        right_orbit[right_origin_1],
    )
    products = (
        factors[0] * factors[2] % prime,
        factors[1] * factors[2] % prime,
        factors[0] * factors[3] % prime,
        factors[1] * factors[3] % prime,
    )
    hadamard = (
        sum(products) % prime,
        (products[0] + products[1] - products[2] - products[3]) % prime,
        (products[0] - products[1] + products[2] - products[3]) % prime,
        (products[0] - products[1] - products[2] + products[3]) % prime,
    )
    require(
        factors == EXPECTED_FIELDS[prime]["factors"]
        and products == EXPECTED_FIELDS[prime]["products"]
        and hadamard == EXPECTED_FIELDS[prime]["hadamard"],
        "stored four-corner field certificate changed",
    )
    require(
        all(factors) and all(products) and all(hadamard),
        "four-corner nonvanishing certificate failed",
    )
    require(
        (
            hadamard[0] * hadamard[3]
            - hadamard[1] * hadamard[2]
        ) % prime == 0,
        "Segre--Hadamard Pluecker identity failed",
    )
    require(
        hadamard[3]
        == (factors[0] - factors[1])
        * (factors[2] - factors[3]) % prime,
        "D3 is not the coefficient mixed face",
    )

    left_literal_covariances = translation_covariances(
        left_bare, left_present, zeta, prime
    )
    right_literal_covariances = translation_covariances(
        right_bare, right_present, zeta, prime
    )
    left_orbit_covariances = translation_covariances(
        left_bare, left_orbit, zeta, prime
    )
    right_orbit_covariances = translation_covariances(
        right_bare, right_orbit, zeta, prime
    )
    fixed_flag_factors = (
        left_bare[LEFT_ORIGIN],
        left_present[LEFT_ORIGIN],
        right_bare[RIGHT_ORIGIN],
        right_present[RIGHT_ORIGIN],
    )
    translated_flag_factors = (
        left_bare[LEFT_ORIGIN],
        left_present[left_origin_1],
        right_bare[RIGHT_ORIGIN],
        right_present[right_origin_1],
    )

    def flag_square(flag_factors):
        flag_products = (
            flag_factors[0] * flag_factors[2] % prime,
            flag_factors[1] * flag_factors[2] % prime,
            flag_factors[0] * flag_factors[3] % prime,
            flag_factors[1] * flag_factors[3] % prime,
        )
        flag_hadamard = (
            sum(flag_products) % prime,
            (
                flag_products[0] + flag_products[1]
                - flag_products[2] - flag_products[3]
            ) % prime,
            (
                flag_products[0] - flag_products[1]
                + flag_products[2] - flag_products[3]
            ) % prime,
            (
                flag_products[0] - flag_products[1]
                - flag_products[2] + flag_products[3]
            ) % prime,
        )
        require(
            flag_hadamard[3]
            == (flag_factors[0] - flag_factors[1])
            * (flag_factors[2] - flag_factors[3]) % prime,
            "flag-square D3 factorization failed",
        )
        require(
            all(flag_factors)
            and all(flag_products)
            and all(flag_hadamard),
            "flag-square all-nonzero gate failed",
        )
        require(
            (
                flag_hadamard[0] * flag_hadamard[3]
                - flag_hadamard[1] * flag_hadamard[2]
            ) % prime == 0,
            "flag-square Pluecker identity failed",
        )
        return flag_products, flag_hadamard

    fixed_flag_products, fixed_flag_hadamard = flag_square(
        fixed_flag_factors
    )
    translated_flag_products, translated_flag_hadamard = flag_square(
        translated_flag_factors
    )
    base_product = (
        left_present[LEFT_ORIGIN]
        * right_present[RIGHT_ORIGIN] % prime
    )
    primal_vectors = {}
    primal_mobius = {}
    for source_twist in range(P):
        for target_twist in range(P):
            vector = (
                base_product,
                base_product if source_twist == 0 else 0,
                base_product if target_twist == 0 else 0,
                base_product
                if source_twist == 0 and target_twist == 0
                else 0,
            )
            primal_vectors[(source_twist, target_twist)] = vector
            primal_mobius[(source_twist, target_twist)] = (
                vector[0] - vector[1] - vector[2] + vector[3]
            ) % prime
    primal_support_counts = tuple(
        sum(vector[index] != 0 for vector in primal_vectors.values())
        for index in range(4)
    )
    primal_all_four = tuple(
        point for point, vector in primal_vectors.items() if all(vector)
    )
    primal_mobius_support = tuple(
        point for point, value in primal_mobius.items() if value
    )
    require(
        primal_support_counts == (P**2, P, P, 1)
        and primal_all_four == ((0, 0),)
        and primal_vectors[(0, 0)] == (base_product,) * 4
        and primal_mobius[(0, 0)] == 0
        and primal_mobius_support
        == tuple(
            (source_twist, target_twist)
            for source_twist in range(1, P)
            for target_twist in range(1, P)
        )
        and all(
            primal_mobius[point] == base_product
            for point in primal_mobius_support
        )
        and sum(primal_mobius.values()) % prime
        == fixed_flag_hadamard[3],
        "raw carrier-twist Mobius census changed",
    )
    require(
        fixed_flag_factors == (
            P * left_present[LEFT_ORIGIN] % prime,
            left_present[LEFT_ORIGIN],
            P * right_present[RIGHT_ORIGIN] % prime,
            right_present[RIGHT_ORIGIN],
        )
        and fixed_flag_products == tuple(
            multiplier * base_product % prime
            for multiplier in (P**2, P, P, 1)
        )
        and fixed_flag_hadamard == tuple(
            multiplier * base_product % prime
            for multiplier in (
                (P + 1)**2, P**2 - 1, P**2 - 1, (P - 1)**2
            )
        ),
        "universal fixed-sheet allocation law failed",
    )
    inverse_p = pow(P, -1, prime)
    require(
        ((0, 0), (0, 0), inverse_p) in left_literal_covariances
        and ((0, 0), (0, 0), inverse_p) in right_literal_covariances,
        "intrinsic zero-step scalar covariance missing",
    )
    require(
        not left_orbit_covariances and not right_orbit_covariances,
        "moving-orbit bank unexpectedly became a bare allocation translate",
    )
    return {
        "prime": prime,
        "factors": factors,
        "products": products,
        "hadamard": hadamard,
        "pluecker_residual": (
            hadamard[0] * hadamard[3]
            - hadamard[1] * hadamard[2]
        ) % prime,
        "left_dual_digest": bank_digest(
            left_dual,
            tuple(
                (address, harmonic)
                for address in endpoint_base.KEYS
                for harmonic in range(P)
            ),
        ),
        "right_dual_digest": bank_digest(
            right_dual,
            tuple(
                (address, harmonic)
                for address in endpoint_base.KEYS
                for harmonic in range(P)
            ),
        ),
        "fixed_left_mismatches": fixed_left_mismatches,
        "fixed_right_mismatches": fixed_right_mismatches,
        "first_fixed_left": first_fixed_left,
        "first_fixed_right": first_fixed_right,
        "left_bare_support": sum(bool(v) for v in left_bare.values()),
        "right_bare_support": sum(bool(v) for v in right_bare.values()),
        "fixed_sheet_left_gauge_rows": fixed_sheet_left_gauge_rows,
        "fixed_sheet_right_gauge_rows": fixed_sheet_right_gauge_rows,
        "left_literal_covariances": left_literal_covariances,
        "right_literal_covariances": right_literal_covariances,
        "left_orbit_covariances": left_orbit_covariances,
        "right_orbit_covariances": right_orbit_covariances,
        "fixed_flag_factors": fixed_flag_factors,
        "fixed_flag_products": fixed_flag_products,
        "fixed_flag_hadamard": fixed_flag_hadamard,
        "base_product": base_product,
        "primal_support_counts": primal_support_counts,
        "primal_all_four": primal_all_four,
        "primal_mobius_support_count": len(primal_mobius_support),
        "primal_mobius_sum": sum(primal_mobius.values()) % prime,
        "translated_flag_factors": translated_flag_factors,
        "translated_flag_products": translated_flag_products,
        "translated_flag_hadamard": translated_flag_hadamard,
    }


def endpoint_audit(source_atom, target_atom):
    present_sets = endpoint_base.present_cache()
    unit = endpoint_base.T // P
    # The constructor identity is proved factorwise by the common shift
    # ell_i -> ell_i+W_i.  Replay representative addresses directly; the
    # ensuing gauge evaluation is exhaustive over all 169*13 cells.
    direct_covariance_controls = (
        (0, 0), (1, 0), (0, 1), (3, 7),
    )
    for address in direct_covariance_controls:
        ell = endpoint_base.REPS[address]
        shifted_ell = tuple(
            (ell[index] + endpoint_base.WMOD[index]) % P
            for index in range(9)
        )
        direct = tuple(
            endpoint_base.endpoint.build_set(
                endpoint_base.PAT_E3, shifted_ell
            )
        )
        translated = physical.overlap.shift_union(
            present_sets[address], -unit
        )
        require(
            direct == translated,
            f"present-set representative covariance failed at {address}",
        )
    return tuple(
        endpoint_field_audit(
            source_atom, target_atom, present_sets, embedding,
        )
        for embedding in endpoint_base.endpoint.MODS
    )


def determinant_square():
    left_1 = add(LEFT_ORIGIN, SOURCE_STEP)
    right_1 = add(RIGHT_ORIGIN, TARGET_STEP)
    endpoint_pairs = (
        (LEFT_ORIGIN, RIGHT_ORIGIN),
        (left_1, RIGHT_ORIGIN),
        (LEFT_ORIGIN, right_1),
        (left_1, right_1),
    )
    q_values = tuple(
        (
            (left[0] - right[0]) % P,
            (left[1] - right[1]) % P,
        )
        for left, right in endpoint_pairs
    )
    determinant_values = tuple(
        det(left, right) for left, right in endpoint_pairs
    )
    mixed = (
        determinant_values[0] - determinant_values[1]
        - determinant_values[2] + determinant_values[3]
    ) % P
    require(
        q_values == ((0, 0), (0, 1), (1, 0), (1, 1))
        and determinant_values == (0, 0, 0, 1)
        and mixed == det(SOURCE_STEP, TARGET_STEP) == 1,
        "endpoint determinant parallelogram changed",
    )
    return endpoint_pairs, q_values, determinant_values, mixed


def main():
    (
        _module, full_module, details, e3, clocks, q_pairs
    ) = build_geometry()
    selected, wing_census, common_clock_census = exhaustive_wing_control(
        details, full_module, e3, clocks
    )
    (
        source_atom, target_atom, factor_signature, semantic_values,
        allocation_masks,
    ) = selected_atom_audit(
        selected, full_module, e3, clocks, q_pairs
    )
    ancestry = ancestry_sheet_audit()
    endpoint_pairs, q_values, determinant_values, determinant_mixed = (
        determinant_square()
    )
    fields = endpoint_audit(source_atom, target_atom)
    source_mask_bits = tuple(bool(row) for row in allocation_masks[0])
    target_mask_bits = tuple(bool(row) for row in allocation_masks[1])
    source_gauge_bits = tuple(
        source_mask_bits[(index - 1) % P] for index in range(P)
    )
    target_gauge_bits = tuple(
        target_mask_bits[(index + 1) % P] for index in range(P)
    )
    require(
        source_gauge_bits == (False, True) + (False,) * (P - 2)
        and target_gauge_bits
        == (False,) * (P - 1) + (True,)
        and sum(source_mask_bits) == sum(source_gauge_bits) == 1
        and sum(target_mask_bits) == sum(target_gauge_bits) == 1,
        "fixed-sheet representative mask transport failed",
    )

    print("same-sheet common-carrier endpoint-square exact audit")
    print(f"dependency_sha256={tuple(PINNED.items())}")
    print(f"physical_cell_(s,t,clock)={CELL}")
    print(f"aligned_objects=A_source,B_target_pullback,C=A_intersect_B")
    print(f"opposite_wing_census=(cells,source_nonempty,target_nonempty)"
          f"={wing_census};opposite_wings_disjoint=567/567")
    print(f"nonempty_common_cells_by_clock_0_to_6={common_clock_census};"
          "fixed_rail_eight_common_supports_only_3_of_7_clocks")
    print(f"common_weighted_piece={(*ATOM_INTERVAL, ATOM_WEIGHT)}")
    print(f"unit_sheet_interval_weight1={(*ATOM_INTERVAL, 1)}")
    print(f"unit_sheet_target_translate={target_atom[0]}")
    print("atom_selector_idempotence=(1_xi*C_source_0=1_xi,"
          "1_Txi*C_target_0=1_Txi):PASS")
    print(f"fixed_xi_carrier_twist_intersections="
          f"(source={allocation_masks[0]},target={allocation_masks[1]})")
    print(f"fixed_xi_carrier_twist_support_bits="
          f"(source={source_mask_bits},target={target_mask_bits})")
    print(f"fixed_xi_representative_gauge_mask_transport="
          f"(source_delta0_to_delta1={source_gauge_bits},"
          f"target_delta0_to_delta12={target_gauge_bits},"
          "k0_sums=(1,1))")
    print("factor_order=(E3,clock,q1,q2,c2,c3)")
    print(
        "complete_simultaneous_factor_signature="
        f"(source_native={factor_signature[0]},"
        f"source_adjacent={factor_signature[1]},"
        f"target_native={factor_signature[2]},"
        f"target_adjacent={factor_signature[3]})"
    )
    print("first_failure_labels=(source=common_after_all_factors,"
          "target=common_after_all_factors)")
    print(f"semantic_states=(source={SOURCE_STATE},target={TARGET_STATE},"
          f"delayed_half={DELAYED_HALF})")
    print(f"unit_delayed_values=(source_carry12={semantic_values[0]},"
          f"target_carry6={semantic_values[1]})")
    print(f"boolean_ancestry_factorization="
          f"{ancestry['u_counts'][0]}*{ancestry['v_counts'][0]}"
          f"={ATOM_WEIGHT}")
    print(f"boolean_label_sets=(U_source_target={ancestry['u_counts']},"
          f"V_source_target={ancestry['v_counts']},mismatches=(0,0))")
    print(f"boolean_label_digests=(U={ancestry['u_digest']},"
          f"V={ancestry['v_digest']})")
    print(f"concrete_boolean_label_(a,b,eprime)={ANCESTRY_LABEL}")
    print(f"concrete_label_inverse_branch_constraints="
          f"{ancestry['constraints']}")
    print(f"concrete_label_joint_source_x_range={ancestry['allowed']}")
    print(f"selected_interval_slacks_inside_joint_range="
          f"{(ATOM_INTERVAL[0]-ancestry['allowed'][0], ancestry['allowed'][1]-ATOM_INTERVAL[1])}")
    print("aggregate_forgets=(a,b,eprime);factor/clock/semantic/carry/"
          "endpoint masks depend on projected x and are constant on the "
          "selected interval; the concrete triple is retained separately")
    print("fixed_origin_full_pullback_address="
          "(r=0,k=0,l=0,L=(0,0),R=(0,0),q=0,Delta=0);"
          "r=0_is_the_distinguished_e0_section")
    print("representative_gauge=(ell,a,b)~(ell+W,a+1,b-1):"
          "PASS_2197_left_and_2197_right_in_each_field")
    for field in fields:
        print(f"field_prime={field['prime']}")
        print(f"dual_support=(left=522/2197,right=522/2197);"
              f"dual_digests=({field['left_dual_digest']},"
              f"{field['right_dual_digest']})")
        print(f"fixed_carrier_gauge_hostile=(left_mismatches="
              f"{field['fixed_left_mismatches']},right_mismatches="
              f"{field['fixed_right_mismatches']},first_left="
              f"{field['first_fixed_left']},first_right="
              f"{field['first_fixed_right']})")
        print(f"fixed_sheet_h0_gauge_rows="
              f"(source_h0_to_h1={field['fixed_sheet_left_gauge_rows']}/169,"
              f"target_h0_to_h12="
              f"{field['fixed_sheet_right_gauge_rows']}/169)")
        print(f"bare_endpoint_support=(left={field['left_bare_support']}/169,"
              f"right={field['right_bare_support']}/169)")
        print(f"literal_fixed_xi_bare_to_present_covariances="
              f"(left={field['left_literal_covariances']},"
              f"right={field['right_literal_covariances']})")
        print(f"moving_xi_orbit_bare_covariances="
              f"(left={field['left_orbit_covariances']},"
              f"right={field['right_orbit_covariances']})")
        print(f"moving_orbit_endpoint_samples_(P_L,P_LplusS,Q_R,Q_RplusT)="
              f"{field['factors']}")
        print(f"moving_orbit_sample_products_(v00,v10,v01,v11)="
              f"{field['products']}")
        print(f"moving_orbit_sample_hadamard_(D0,D1,D2,D3)="
              f"{field['hadamard']}")
        print(f"moving_orbit_sample_mixed_face_D3={field['hadamard'][3]};"
              f"nonzero={field['hadamard'][3] != 0}")
        print(f"pluecker_D0D3_minus_D1D2="
              f"{field['pluecker_residual']} (rank_one_identity)")
        print("atom_local_allocation_constructor=bare keeps fixed xi "
              "constant over the omitted carrier-twist variable; present "
              "keeps outer selector xi and multiplies by C_h, so only h=0 "
              "survives; therefore Pbare=13*Ppresent and "
              "Qbare=13*Qpresent at every endpoint")
        print(f"fixed_origin_allocation_factors_(Pbare,Ppresent,"
              f"Qbare,Qpresent)={field['fixed_flag_factors']}")
        print(f"fixed_origin_allocation_products="
              f"{field['fixed_flag_products']}")
        print(f"fixed_origin_allocation_hadamard_(D0,D1,D2,D3)="
              f"{field['fixed_flag_hadamard']}")
        print(f"fixed_origin_allocation_D3="
              f"{field['fixed_flag_hadamard'][3]};"
              f"nonzero={field['fixed_flag_hadamard'][3] != 0};"
              "status=central_post_carrier_DFT_coefficient")
        print("fixed_origin_universal_allocation_profile="
              "base_product*(169,13,13,1);"
              "hadamard=base_product*(196,168,168,144)")
        print("primal_carrier_twist_support_counts_(B,P,Q,H)="
              f"{field['primal_support_counts']};"
              f"sole_fourfold_point={field['primal_all_four']}")
        print("sole_primal_fourfold_vector_(B,P,Q,H)="
              f"({field['base_product']},{field['base_product']},"
              f"{field['base_product']},{field['base_product']});"
              "raw_mobius=0")
        print("raw_primal_mobius_support="
              f"12x12_bare_only_complement;count="
              f"{field['primal_mobius_support_count']};"
              f"value={field['base_product']};"
              f"central_sum={field['primal_mobius_sum']}")
        print("intrinsic_allocation_endpoint_steps="
              "(source=(0,0),target=(0,0));"
              "scalar_covariances=(1/13,1/13);"
              "determinant_mixed_face=0")
        print(f"translated_allocation_candidate_factors_(Pbare00,Ppresent10,"
              f"Qbare00,Qpresent01)={field['translated_flag_factors']}")
        print(f"translated_allocation_candidate_products="
              f"{field['translated_flag_products']}")
        print(f"translated_allocation_candidate_hadamard_(D0,D1,D2,D3)="
              f"{field['translated_flag_hadamard']}")
        print(f"translated_allocation_candidate_D3="
              f"{field['translated_flag_hadamard'][3]};"
              f"nonzero={field['translated_flag_hadamard'][3] != 0};"
              "status=extrinsic_nonzero_endpoint_sample_not_the_fixed_xi_"
              "allocation_toggle")
    print(f"endpoint_pairs={endpoint_pairs}")
    print(f"endpoint_target_differences_q={q_values}")
    print(f"endpoint_determinant_values={determinant_values}")
    print(f"endpoint_determinant_mixed_face=det(source_step,target_step)="
          f"{determinant_mixed}")
    print("scope=one literal common Boolean ancestry label gives a fixed-xi "
          "central post-carrier-DFT allocation coefficient square with "
          "nonzero D3 in two exact-order fields.  The sole primal twist "
          "point with all four allocation states has raw Mobius face zero; "
          "the central D3 is the sum of the 12x12 bare-only complement, so "
          "this is not a pre-Fourier common-atom mixed face.  Its intrinsic "
          "allocation action is scalar with zero endpoint steps, hence zero "
          "determinant face.  Moving "
          "xi through the carrier orbit gives a separate endpoint "
          "parallelogram with nonzero sample D3 and determinant-label face, "
          "but exhaustive empty orbit/bare covariance proves it is not the "
          "fixed-xi allocation.  This fixed rail-eight bank has common cells "
          "on only 3 of 7 clocks.  No seven-clock Cech/root-deck map, row "
          "exclusion, THM-2772 common-atom gate, or LRC14 conclusion is "
          "inferred.")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
