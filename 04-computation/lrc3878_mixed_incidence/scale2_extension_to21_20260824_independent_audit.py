#!/usr/bin/env python3
"""No-import exact audit of THM-4002's max(E)=21 extension.

The finite engine is deliberately rebuilt here.  It does not import either
scale-two production script.  For each new body it constructs the positive-
length closed safe arcs, subtracts the exact open pullback of the two-lift
obstruction, proves that the escaped set has positive rational measure, and
checks an explicit rational midpoint in physical coordinates.

Touching open danger components are merged modulo isolated boundary points.
This only discards valid equality witnesses and is therefore conservative.
Literal endpoint behavior is audited separately from the positive-length
enumeration.
"""

from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from math import comb as binomial, isqrt
import json
import sys


sys.stdout.reconfigure(newline="\n")

RADIUS = Q(1, 14)
LEFT_COMPONENT = (Q(2, 21), Q(8, 63))
RIGHT_COMPONENT = (Q(55, 63), Q(19, 21))
OBSTRUCTION = (LEFT_COMPONENT, RIGHT_COMPONENT)
OBSTRUCTION_MASS = Q(4, 63)

OLD_BODY_COUNT = 167_960
OLD_FINITE_CELL_COUNT = 574_963
EXPECTED_NEW_BODY_COUNT = 184_756
EXPECTED_NEW_FINITE_CELL_COUNT = 781_184
EXPECTED_COMBINED_BODY_COUNT = 352_716
EXPECTED_COMBINED_FINITE_CELL_COUNT = 1_356_147
EXPECTED_MAX_TAIL = 80
EXPECTED_MAX_TAIL_BODY = (1, 3, 4, 11, 13, 15, 16, 17, 18, 19, 21)
EXPECTED_SEMANTIC_SHA256 = "3735c3267d5100d7005ce58952c3a4896baa9f10cf9b31005e9a5e736eec112f"

CHECKS = 0


def require(condition, label):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(label)


def fmt(value):
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def merge_mod_points(intervals):
    """Merge positive-length intervals, conservatively filling touch points."""
    merged = []
    for left, right in sorted(intervals):
        if left >= right:
            continue
        if not merged or left > merged[-1][1]:
            merged.append([left, right])
        elif right > merged[-1][1]:
            merged[-1][1] = right
    return tuple((left, right) for left, right in merged)


def speed_danger_comb(speed):
    """Positive-length model of {y: ||speed*y|| < 1/14} on [0,1]."""
    half_width = Q(1, 14 * speed)
    intervals = []
    for residue in range(speed):
        center = Q(residue, speed)
        left, right = center - half_width, center + half_width
        if left < 0:
            intervals.append((Q(0), right))
            intervals.append((left + 1, Q(1)))
        elif right > 1:
            intervals.append((left, Q(1)))
            intervals.append((Q(0), right - 1))
        else:
            intervals.append((left, right))
    return merge_mod_points(intervals)


def positive_complement(intervals):
    answer = []
    position = Q(0)
    for left, right in merge_mod_points(intervals):
        if position < left:
            answer.append((position, left))
        position = max(position, right)
    if position < 1:
        answer.append((position, Q(1)))
    return tuple(answer)


def body_safe_arcs(body, danger_combs):
    return positive_complement(
        interval for speed in body for interval in danger_combs[speed]
    )


def pullback_obstruction(multiplier):
    return tuple(
        ((left + cell) / multiplier, (right + cell) / multiplier)
        for cell in range(multiplier)
        for left, right in OBSTRUCTION
    )


def subtract_open(closed_intervals, open_intervals):
    """Positive-length pieces of a closed union outside an open union."""
    open_intervals = tuple(sorted(open_intervals))
    answer = []
    global_index = 0
    for left, right in closed_intervals:
        position = left
        while (global_index < len(open_intervals)
               and open_intervals[global_index][1] <= position):
            global_index += 1
        index = global_index
        while index < len(open_intervals) and open_intervals[index][0] < right:
            bad_left, bad_right = open_intervals[index]
            if position < bad_left:
                answer.append((position, min(bad_left, right)))
            position = max(position, bad_right)
            if position >= right:
                break
            index += 1
        if position < right:
            answer.append((position, right))
    return tuple((left, right) for left, right in answer if left < right)


def clearance(speed, phase):
    residue = (speed * phase) % 1
    return min(residue, 1 - residue)


def unit_clearance(value):
    residue = value % 1
    return min(residue, 1 - residue)


def strict_bv_tail(mass, component_count):
    wall = (Q(component_count * component_count) * OBSTRUCTION_MASS
            / (3 * mass * mass * (1 - OBSTRUCTION_MASS)))
    candidate = isqrt(wall.numerator // wall.denominator)
    while Q(candidate * candidate) <= wall:
        candidate += 1
    return candidate, wall


def open_intersection(first, second):
    intersections = []
    i = j = 0
    first, second = tuple(first), tuple(second)
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            intersections.append((left, right))
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return merge_mod_points(intersections)


def affine_danger(multiplier, parity):
    """z with ||(multiplier*z+parity)/2|| < 1/14, as open arcs."""
    intervals = []
    for integer in range(-2, multiplier + 3):
        left = Q(2 * integer - parity, multiplier) - Q(1, 7 * multiplier)
        right = Q(2 * integer - parity, multiplier) + Q(1, 7 * multiplier)
        left, right = max(Q(0), left), min(Q(1), right)
        if left < right:
            intervals.append((left, right))
    return merge_mod_points(intervals)


def region_contains_open(region, point):
    return any(left < point < right for left, right in region)


def obstruction_contains_open(point):
    return region_contains_open(OBSTRUCTION, point)


def parity_lift_is_safe(z, parity):
    return (unit_clearance((z + parity) / 2) >= RADIUS
            and unit_clearance((9 * z + parity) / 2) >= RADIUS)


def both_parity_lifts_are_blocked(z):
    return not parity_lift_is_safe(z, 0) and not parity_lift_is_safe(z, 1)


def exact_obstruction_controls():
    """Rebuild C from the two physical lifts and probe hostile conventions."""
    blocked_by_parity = []
    for parity in (0, 1):
        t_bad = affine_danger(1, parity)
        nine_t_bad = affine_danger(9, parity)
        blocked_by_parity.append(merge_mod_points(t_bad + nine_t_bad))

    rebuilt = open_intersection(blocked_by_parity[0], blocked_by_parity[1])
    require(rebuilt == OBSTRUCTION, "independent physical reconstruction of C")
    require(RIGHT_COMPONENT == (1 - LEFT_COMPONENT[1], 1 - LEFT_COMPONENT[0]),
            "reflection of obstruction components")
    require(sum((right - left for left, right in rebuilt), Q(0))
            == OBSTRUCTION_MASS, "reconstructed obstruction mass")

    critical = {Q(0), Q(1)}
    for region in (*blocked_by_parity, rebuilt):
        for left, right in region:
            critical.update((left, right))
    ordered = sorted(critical)
    samples = set(ordered)
    samples.update((left + right) / 2 for left, right in zip(ordered, ordered[1:])
                   if left < right)
    for z in sorted(samples):
        require(both_parity_lifts_are_blocked(z) == obstruction_contains_open(z),
                f"physical/C partition equivalence at {fmt(z)}")

    endpoints = tuple(value for interval in OBSTRUCTION for value in interval)
    require(all(not both_parity_lifts_are_blocked(z) for z in endpoints),
            "literal obstruction endpoints are safe equality witnesses")
    require(all(both_parity_lifts_are_blocked((left + right) / 2)
                for left, right in OBSTRUCTION),
            "component interiors block both lifts")

    # Hostile 1: replacing the open components by closed components falsely
    # deletes all four equality witnesses.
    closed_endpoint_false_positives = sum(
        not both_parity_lifts_are_blocked(z) for z in endpoints
    )
    require(closed_endpoint_false_positives == 4,
            "closed-endpoint mutation detected four false positives")

    # Hostile 2: for even t the two physical lifts do not switch parity.
    # The corresponding one-parity blocked regions are strictly larger than C.
    parity_witnesses = []
    for region in blocked_by_parity:
        cuts = sorted({Q(0), Q(1),
                       *(value for interval in region for value in interval),
                       *(value for interval in rebuilt for value in interval)})
        witness = None
        for left, right in zip(cuts, cuts[1:]):
            mid = (left + right) / 2
            if region_contains_open(region, mid) != obstruction_contains_open(mid):
                witness = mid
                break
        require(witness is not None, "even-parity mutation has a witness")
        parity_witnesses.append(witness)

    # Hostile 3: omitting the reflected component halves the obstruction and
    # misses a point which the direct physical predicate declares blocked.
    reflected_midpoint = sum(RIGHT_COMPONENT, Q(0)) / 2
    require(both_parity_lifts_are_blocked(reflected_midpoint),
            "reflected-component omission mutation detected")
    require(not region_contains_open((LEFT_COMPONENT,), reflected_midpoint),
            "one-sided mutation really omits reflected midpoint")

    return {
        "rebuilt": [[fmt(a), fmt(b)] for a, b in rebuilt],
        "partition_samples": len(samples),
        "closed_endpoint_false_positives": closed_endpoint_false_positives,
        "even_parity_witnesses": [fmt(z) for z in parity_witnesses],
        "reflected_omission_witness": fmt(reflected_midpoint),
    }


def physical_lifts(y, odd_multiplier):
    valid = []
    for lift_index in (0, 1):
        x = (y + lift_index) / 2
        if (clearance(odd_multiplier, x) >= RADIUS
                and clearance(9 * odd_multiplier, x) >= RADIUS):
            valid.append(lift_index)
    return tuple(valid)


def audit_new_stratum():
    danger_combs = {speed: speed_danger_comb(speed) for speed in range(1, 22)}
    pullbacks = {odd: pullback_obstruction(odd) for odd in range(21, 150, 2)}
    transcript = sha256()

    body_count = 0
    finite_cell_count = 0
    maximum_tail = 0
    maximum_tail_body = None
    maximum_tail_wall = None
    minimum_escape = None

    for lower_ten in combinations(range(1, 21), 10):
        body = lower_ten + (21,)
        body_count += 1
        safe = body_safe_arcs(body, danger_combs)
        mass = sum((right - left for left, right in safe), Q(0))
        require(mass > 0, f"positive safe mass {body}")
        tail, wall = strict_bv_tail(mass, len(safe))
        require(Q((tail - 1) ** 2) <= wall < Q(tail ** 2),
                f"minimal strict BV tail {body}")
        if (tail, body) > (maximum_tail, maximum_tail_body or ()):
            maximum_tail = tail
            maximum_tail_body = body
            maximum_tail_wall = wall

        transcript.update(
            ("B|" + repr(body) + "|" + fmt(mass) + "|"
             + str(len(safe)) + "|" + str(tail) + "\n").encode("ascii")
        )

        for odd_multiplier in range(21, tail, 2):
            escaped = subtract_open(safe, pullbacks[odd_multiplier])
            escaped_mass = sum((right - left for left, right in escaped), Q(0))
            require(escaped_mass > 0,
                    f"positive exact escape {(body, odd_multiplier)}")

            y = sum(escaped[0], Q(0)) / 2
            require(all(clearance(speed, y) >= RADIUS for speed in body),
                    f"body-safe exact midpoint {(body, odd_multiplier)}")
            z = (odd_multiplier * y) % 1
            require(not obstruction_contains_open(z),
                    f"midpoint outside open quotient obstruction {(body, odd_multiplier)}")

            lifts = physical_lifts(y, odd_multiplier)
            require(bool(lifts), f"explicit physical lift {(body, odd_multiplier)}")
            require(bool(lifts) == (not both_parity_lifts_are_blocked(z)),
                    f"quotient-to-physical equivalence {(body, odd_multiplier)}")
            lift_index = lifts[0]
            x = (y + lift_index) / 2
            require(all(clearance(2 * speed, x) >= RADIUS for speed in body),
                    f"physical body clearance {(body, odd_multiplier, lift_index)}")
            require(clearance(odd_multiplier, x) >= RADIUS,
                    f"physical t clearance {(body, odd_multiplier, lift_index)}")
            require(clearance(9 * odd_multiplier, x) >= RADIUS,
                    f"physical 9t clearance {(body, odd_multiplier, lift_index)}")

            record = (escaped_mass, body, odd_multiplier, y, lifts)
            if minimum_escape is None or record < minimum_escape:
                minimum_escape = record
            transcript.update(
                ("T|" + str(odd_multiplier) + "|" + fmt(escaped_mass) + "|"
                 + fmt(y) + "|" + repr(lifts) + "\n").encode("ascii")
            )
            finite_cell_count += 1

    require(body_count == EXPECTED_NEW_BODY_COUNT, "new max-21 body count")
    require(body_count == binomial(20, 10), "max-21 combinatorial body count")
    require(finite_cell_count == EXPECTED_NEW_FINITE_CELL_COUNT,
            "new max-21 finite cell count")
    require(maximum_tail == EXPECTED_MAX_TAIL, "largest strict BV tail")
    require(maximum_tail_body == EXPECTED_MAX_TAIL_BODY,
            "lexicographically final extremal tail body")
    require(maximum_tail < max(pullbacks) + 2, "pullback cache covers finite range")

    # The combined finite-cell number is an exact disjoint-stratum sum with
    # the already independently audited [1,20] count.  This audit does not
    # rerun that old stratum.
    require(OLD_BODY_COUNT + body_count == EXPECTED_COMBINED_BODY_COUNT,
            "combined body count")
    require(binomial(21, 11) == EXPECTED_COMBINED_BODY_COUNT,
            "combined combinatorial body count")
    require(OLD_FINITE_CELL_COUNT + finite_cell_count
            == EXPECTED_COMBINED_FINITE_CELL_COUNT, "combined finite cell count")

    escape_mass, escape_body, escape_t, escape_y, escape_lifts = minimum_escape
    return {
        "new_bodies": body_count,
        "new_finite_cells": finite_cell_count,
        "combined_bodies": OLD_BODY_COUNT + body_count,
        "combined_finite_cells": OLD_FINITE_CELL_COUNT + finite_cell_count,
        "max_tail": maximum_tail,
        "max_tail_body": list(maximum_tail_body),
        "max_tail_wall": fmt(maximum_tail_wall),
        "minimum_positive_escape": {
            "measure": fmt(escape_mass),
            "body": list(escape_body),
            "t": escape_t,
            "y": fmt(escape_y),
            "lifts": list(escape_lifts),
        },
        "enumeration_transcript_sha256": transcript.hexdigest(),
    }


def main():
    controls = exact_obstruction_controls()
    census = audit_new_stratum()
    semantic_record = {"controls": controls, "census": census}
    semantic_digest = sha256(
        json.dumps(semantic_record, sort_keys=True,
                   separators=(",", ":")).encode("ascii")
    ).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_FROZEN":
        require(semantic_digest == EXPECTED_SEMANTIC_SHA256, "semantic freeze")

    minimum = census["minimum_positive_escape"]
    print("LRC14_SCALE2_EXTENSION_TO21_NO_IMPORT_EXACT_AUDIT_20260824")
    print("scope=new_stratum_E_size_11_subset_1..21_with_maxE_21;odd_t>=21")
    print(f"new_bodies={census['new_bodies']};new_finite_cells={census['new_finite_cells']};failures=0")
    print(f"combined_bodies={census['combined_bodies']};combined_finite_cells={census['combined_finite_cells']}")
    print("combined_cell_basis=already_audited_height20_count_plus_independent_new_stratum")
    print(f"max_tail={census['max_tail']};max_tail_body={tuple(census['max_tail_body'])}")
    print(f"max_tail_wall={census['max_tail_wall']};control=(79^2<=wall<80^2)")
    print("minimum_positive_escape=" + repr((
        minimum["measure"], tuple(minimum["body"]), minimum["t"],
        minimum["y"], tuple(minimum["lifts"]))))
    print("finite_method=positive_rational_measure_of_closed_body_safe_arcs_minus_open_pullback")
    print("physical_check=explicit_midpoint_and_one_selected_valid_lift_checked_at_all_13_speeds")
    print("obstruction_rebuild=" + repr(tuple(tuple(pair) for pair in controls["rebuilt"])))
    print(f"obstruction_partition_samples={controls['partition_samples']}")
    print("hostiles=closed_endpoints_rejected_4;even_parity_mutation_witnessed;reflected_component_omission_witnessed")
    print("enumeration_transcript_sha256=" + census["enumeration_transcript_sha256"])
    print("semantic_sha256=" + semantic_digest)
    print(f"checks={CHECKS}")
    print("limitations=fixed_scale_two_(1,9)_row_only;odd_t>=21;old_stratum_not_rerun;no_arbitrary_body_or_LRC14_claim")
    print("RESULT=PASS;LRC14=OPEN")


if __name__ == "__main__":
    main()
