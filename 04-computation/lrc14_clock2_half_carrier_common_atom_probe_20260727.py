#!/usr/bin/env python3
"""Exact common-x probe joining the THM-2625 clock to THM-2635 h=3.

This deliberately stops one level short of claiming an endpoint-current
cospan.  It places the literal clock-two packet

    F_2 = E_0 intersect {x : 169 x in Q_a}

and THM-2635's selected epsilon=1, q=h=3, r=9 half-carrier on the same
circle coordinate.  Both are lifted to G=169*T_DEN.  The delayed Q_b digit-3
condition is then integrated exactly at clock 13^6 on that common grid.

Positive output is a genuine common physical atom (with the inherited
nonnegative rail multiplicity), not a product of separately integrated
marginals.  It does *not* by itself prove that THM-2625's Fourier endpoint
current restricts nontrivially, nor that the construction is equivariant under
the endpoint-allocation quotient; those are separate typed gates.
"""

from bisect import bisect_right
from collections import Counter
from fractions import Fraction
from math import gcd

import lrc14_canonical_endpoint_current_thm2625 as endpoint
import lrc14_cross_time_target_future_diagonal_thm2616 as cross
import lrc14_old_wall_successor_sector_thm2630 as wall


P = 13
Q7 = 7
T = endpoint.T_DEN
R2 = endpoint.RDIL
R6 = cross.old.R
G = R2 * T


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


require(T == cross.T, "THM-2625/2635 base grids differ")
require(G == endpoint.NN, "clock-two common grid drift")


def clock_two_packet():
    """Materialize E_0 intersect (169.)^-1 Q_a on the G-grid."""
    e_intervals = endpoint.build_set(endpoint.PAT_E, endpoint.ZERO_ELL)
    q_intervals = endpoint.build_set(endpoint.PAT_QA, endpoint.ZERO_ELL)
    q_starts = [left for left, _ in q_intervals]
    count_q = len(q_intervals)
    packet = []
    for left, right in e_intervals:
        lifted_left = R2 * left
        source_left = lifted_left % T
        span = R2 * (right - left)
        require(span < T, "clock-two E interval needs multiple source wraps")
        source_right = source_left + span
        lift_offset = lifted_left - source_left
        index = bisect_right(q_starts, source_left) - 1
        offset = 0
        if index < 0:
            index = count_q - 1
            offset = -T
        while True:
            qa0, qb0 = q_intervals[index]
            qa = qa0 + offset
            if qa >= source_right:
                break
            qb = qb0 + offset
            if qb > source_left:
                low = max(source_left, qa)
                high = min(source_right, qb)
                if low < high:
                    packet.append((lift_offset + low, lift_offset + high))
            index += 1
            if index == count_q:
                index = 0
                offset += T
    require(all(0 <= a < b <= G for a, b in packet), "clock-two packet range")
    require(all(packet[i - 1][1] <= packet[i][0]
                for i in range(1, len(packet))), "clock-two packet order")
    length = sum(b - a for a, b in packet)
    require(Fraction(length, G) == Fraction(21376087, 17907461390),
            "clock-two packet measure anchor changed")
    return packet


def intersect_weighted(pieces, intervals):
    """Intersect sorted weighted pieces and an interval union on one grid."""
    starts = [left for left, _ in intervals]
    out = []
    for left, right, weight in pieces:
        index = max(0, bisect_right(starts, left) - 1)
        while index < len(intervals) and intervals[index][0] < right:
            a = max(left, intervals[index][0])
            b = min(right, intervals[index][1])
            if a < b:
                if out and out[-1][1] == a and out[-1][2] == weight:
                    out[-1] = (out[-1][0], b, weight)
                else:
                    out.append((a, b, weight))
            index += 1
    return out


def intersect_weighted_periodic(pieces, intervals, period):
    """Intersect weighted pieces with the periodized interval union."""
    starts = [left for left, _ in intervals]
    output = []
    for left, right, weight in pieces:
        first_period = left // period
        last_period = (right - 1) // period
        for branch in range(first_period, last_period + 1):
            offset = branch * period
            source_left = max(0, left - offset)
            source_right = min(period, right - offset)
            index = max(0, bisect_right(starts, source_left) - 1)
            while index < len(intervals) and intervals[index][0] < source_right:
                a = max(source_left, intervals[index][0]) + offset
                b = min(source_right, intervals[index][1]) + offset
                if a < b:
                    if output and output[-1][1] == a and output[-1][2] == weight:
                        output[-1] = (output[-1][0], b, weight)
                    else:
                        output.append((a, b, weight))
                index += 1
    return output


def scaled_prefix(prefix):
    starts, lengths, _ = prefix
    intervals = [(R2 * start, R2 * (start + length))
                 for start, length in zip(starts, lengths)]
    scaled_starts = [left for left, _ in intervals]
    scaled_lengths = [right - left for left, right in intervals]
    cumulative = [0]
    for length in scaled_lengths:
        cumulative.append(cumulative[-1] + length)
    return scaled_starts, scaled_lengths, cumulative


def phi_at(x, starts, lengths, cumulative):
    index = bisect_right(starts, x) - 1
    if index < 0:
        return 0
    return cumulative[index] + min(x - starts[index], lengths[index])


def delayed_numerator(pieces, prefix):
    """Numerator over R6*G of integral weight(x)1_Q(R6*x) dx."""
    starts, lengths, cumulative = prefix
    q_length = cumulative[-1]
    weighted_length = 0
    residue_accumulator = 0
    phi_accumulator = 0
    reduced_clock = R6 % G
    for left, right, weight in pieces:
        rleft = left * reduced_clock % G
        rright = right * reduced_clock % G
        weighted_length += weight * (right - left)
        residue_accumulator += weight * (rright - rleft)
        phi_accumulator += weight * (
            phi_at(rright, starts, lengths, cumulative)
            - phi_at(rleft, starts, lengths, cumulative)
        )
    floor_numerator = R6 * weighted_length - residue_accumulator
    require(floor_numerator % G == 0, "delayed wrap count is not integral")
    value = q_length * (floor_numerator // G) + phi_accumulator
    require(value >= 0, "negative common-atom numerator")
    return value


def delayed_vector(pieces, prefixes):
    """All thirteen delayed-digit numerators in one common-grid sweep."""
    weighted_length = 0
    residue_accumulator = 0
    phi_accumulators = [0] * P
    reduced_clock = R6 % G
    for left, right, weight in pieces:
        rleft = left * reduced_clock % G
        rright = right * reduced_clock % G
        weighted_length += weight * (right - left)
        residue_accumulator += weight * (rright - rleft)
        for digit, prefix in enumerate(prefixes):
            phi_accumulators[digit] += weight * (
                phi_at(rright, *prefix) - phi_at(rleft, *prefix)
            )
    floor_numerator = R6 * weighted_length - residue_accumulator
    require(floor_numerator % G == 0, "delayed vector wrap count is not integral")
    quotient = floor_numerator // G
    result = tuple(
        prefix[2][-1] * quotient + phi_accumulators[digit]
        for digit, prefix in enumerate(prefixes)
    )
    require(min(result) >= 0, "negative common-atom digit numerator")
    return result


def main():
    packet = clock_two_packet()
    packet_starts = [left for left, _ in packet]
    e_zero = endpoint.build_set(endpoint.PAT_E, endpoint.ZERO_ELL)
    q_zero = endpoint.build_set(endpoint.PAT_QA, endpoint.ZERO_ELL)
    e_lifted = [(R2 * left, R2 * right) for left, right in e_zero]
    (module, prefixes, _whole_prefixes, _digit_masses, rails,
     present, present_starts) = cross.build_carrier_data()

    # THM-2635's unique uniform chronological unit closure.
    q = h = 3
    probe = 9
    epsilon = 1
    shift = (-q) % P
    require((epsilon, h, 0, 4, probe) == (1, 3, 0, 4, 9),
            "chronological closure tuple drift")
    digit_prefixes = [
        [scaled_prefix(prefixes[ell5][digit]) for digit in range(P)]
        for ell5 in range(Q7)
    ]

    rows = {}
    baseline_rows = {}
    common_piece_count = 0
    positive_entries = 0
    common_content = 0
    baseline_content = 0
    prefuture_weighted_length = 0
    prefuture_positive_entries = 0
    digit_totals = [0] * P
    cell_digit_support = {}
    rail_stage_positive = Counter()
    present_packet_positive = 0
    e_only_half_positive = 0
    q_only_half_positive = 0
    e_only_h3_positive = 0
    q_only_h3_positive = 0
    e_only_h3_total = 0
    q_only_h3_total = 0
    for s, ell4, root, pieces in rails:
        theta = (root - 12) % P
        if theta != wall.selected_theta(s, ell4):
            continue
        lifted_rail = [(R2 * left, R2 * right, weight)
                       for left, right, weight in pieces]
        if intersect_weighted(lifted_rail, packet):
            rail_stage_positive["packet"] += 1
        if intersect_weighted(lifted_rail, e_lifted):
            rail_stage_positive["E"] += 1
        if intersect_weighted_periodic(lifted_rail, q_zero, T):
            rail_stage_positive["Q"] += 1
        values = []
        baseline_values = []
        row_digit_support = set()
        for ell5 in range(Q7):
            overlap = cross.old.intersect_weighted_union(
                pieces, present[ell5, shift], present_starts[ell5, shift]
            )
            lifted_overlap = [(R2 * left, R2 * right, weight)
                              for left, right, weight in overlap]
            present_packet_positive += int(bool(
                intersect_weighted(lifted_overlap, packet)
            ))
            half = cross.old.intersect_weighted_comb(
                overlap, module.C3, 182,
                14 * probe - 13, 14 * probe,
            )
            baseline = cross.old.delayed_weighted_numerator(
                half, prefixes[ell5][h]
            )
            lifted = [(R2 * left, R2 * right, weight)
                      for left, right, weight in half]
            common = intersect_weighted(lifted, packet)
            e_only = intersect_weighted(lifted, e_lifted)
            q_only = intersect_weighted_periodic(lifted, q_zero, T)
            both_route = intersect_weighted_periodic(e_only, q_zero, T)
            require(common == both_route,
                    "materialized packet and sequential E/Q routes disagree")
            e_only_half_positive += int(bool(e_only))
            q_only_half_positive += int(bool(q_only))
            # Re-run through the public intersection route as a small hostile
            # control; starts are supplied to avoid rebuilding them repeatedly.
            require(common == cross.old.intersect_weighted_union(
                lifted, packet, packet_starts
            ), "common-grid intersection routes disagree")
            common_weighted_length = sum(
                weight * (right - left) for left, right, weight in common
            )
            prefuture_weighted_length += common_weighted_length
            prefuture_positive_entries += int(common_weighted_length > 0)
            digit_values = delayed_vector(common, digit_prefixes[ell5])
            e_only_h3 = delayed_numerator(e_only, digit_prefixes[ell5][h])
            q_only_h3 = delayed_numerator(q_only, digit_prefixes[ell5][h])
            e_only_h3_positive += int(e_only_h3 > 0)
            q_only_h3_positive += int(q_only_h3 > 0)
            e_only_h3_total += e_only_h3
            q_only_h3_total += q_only_h3
            require(
                digit_values[h]
                == delayed_numerator(common, digit_prefixes[ell5][h]),
                "scalar/vector common delayed routes disagree",
            )
            value = digit_values[h]
            for digit, digit_value in enumerate(digit_values):
                digit_totals[digit] += digit_value
                if digit_value:
                    row_digit_support.add(digit)
            require(value <= R2 * baseline,
                    "clock-two restriction exceeds unrestricted carrier")
            values.append(value)
            baseline_values.append(baseline)
            common_piece_count += len(common)
            positive_entries += int(value > 0)
            common_content = gcd(common_content, value)
            baseline_content = gcd(baseline_content, baseline)
        key = (s, ell4)
        require(key not in rows, "duplicate selected rail")
        rows[key] = tuple(values)
        baseline_rows[key] = tuple(baseline_values)
        cell_digit_support[key] = tuple(sorted(row_digit_support))

    require(len(rows) == 84, "canonical selected-rail census changed")
    positive_cell_counts = Counter(sum(value > 0 for value in row)
                                   for row in rows.values())
    empty_cells = tuple(key for key, row in sorted(rows.items())
                        if not any(value > 0 for value in row))
    zero_entries = tuple(
        (key, ell5) for key, row in sorted(rows.items())
        for ell5, value in enumerate(row) if value == 0
    )

    # A deliberately labelled diagnostic: primitivity is taken only over this
    # h=3 common slice, not the full THM-2635 carrier.
    slice_units = []
    if common_content and all(value % common_content == 0
                              for row in rows.values() for value in row):
        for key, row in sorted(rows.items()):
            normalized = tuple(
                (value // common_content) * pow(probe, -1, P) % P
                for value in row
            )
            reduced = tuple((normalized[i] - normalized[-1]) % P
                            for i in range(Q7 - 1))
            if cross.old.sat.multiplication_determinant_7(reduced):
                slice_units.append(key)

    total_numerator = sum(sum(row) for row in rows.values())
    digit_support_hist = Counter(cell_digit_support.values())
    common_future_digits = tuple(sorted(
        set.intersection(*(set(values) for values in cell_digit_support.values()))
    ))

    # Local hostile atlas: retain q=3 and the entire clock-two packet, but
    # vary the physical half-tooth.  This identifies whether emptiness comes
    # from the clock itself or from the closure's specific (epsilon,r)=(1,9).
    half_atlas = {}
    half_traces = {}
    for epsilon_probe in ((epsilon, candidate_probe)
                          for epsilon in range(2)
                          for candidate_probe in range(1, P)):
        candidate_epsilon, candidate_probe = epsilon_probe
        positive_before = 0
        positive_after = 0
        total_after = 0
        trace = []
        for s, ell4, root, pieces in rails:
            theta = (root - 12) % P
            if theta != wall.selected_theta(s, ell4):
                continue
            for ell5 in range(Q7):
                overlap = cross.old.intersect_weighted_union(
                    pieces, present[ell5, shift], present_starts[ell5, shift]
                )
                if candidate_epsilon == 0:
                    low, high = 14 * candidate_probe, 14 * candidate_probe + 13
                else:
                    low, high = 14 * candidate_probe - 13, 14 * candidate_probe
                half = cross.old.intersect_weighted_comb(
                    overlap, module.C3, 182, low, high
                )
                lifted = [(R2 * left, R2 * right, weight)
                          for left, right, weight in half]
                common = intersect_weighted(lifted, packet)
                positive_before += int(bool(common))
                value = delayed_numerator(common, digit_prefixes[ell5][h])
                trace.append((tuple(common), value))
                positive_after += int(value > 0)
                total_after += value
        half_atlas[epsilon_probe] = (positive_before, positive_after, total_after)
        half_traces[epsilon_probe] = tuple(trace)

    twin_piece_equal_count = sum(
        left[0] == right[0]
        for left, right in zip(half_traces[(0, 11)], half_traces[(1, 12)])
    )
    twin_value_equal_count = sum(
        left[1] == right[1]
        for left, right in zip(half_traces[(0, 11)], half_traces[(1, 12)])
    )

    # The only other abstract chronological closure from THM-2635 is h=7,
    # epsilon=1, r=5.  It is not a uniform unit section, but it is the cheapest
    # reroute if the h=3 common atom is structurally forbidden.
    alternate_closure = []
    alt_q, alt_probe, alt_epsilon = 7, 5, 1
    alt_shift = (-alt_q) % P
    for s, ell4, root, pieces in rails:
        theta = (root - 12) % P
        if theta != wall.selected_theta(s, ell4):
            continue
        for ell5 in range(Q7):
            overlap = cross.old.intersect_weighted_union(
                pieces, present[ell5, alt_shift], present_starts[ell5, alt_shift]
            )
            half = cross.old.intersect_weighted_comb(
                overlap, module.C3, 182,
                14 * alt_probe - 13, 14 * alt_probe,
            )
            lifted = [(R2 * left, R2 * right, weight)
                      for left, right, weight in half]
            common = intersect_weighted(lifted, packet)
            value = delayed_numerator(common, digit_prefixes[ell5][alt_q])
            alternate_closure.append((bool(common), value))
    print("THM-2625 clock-two / THM-2635 h=3 common-x probe")
    print(f"T={T} G=169*T={G} delayed_clock={R6}")
    print(f"clock_two_packet_intervals={len(packet)}")
    print(f"selected_rows={len(rows)} common_piece_count={common_piece_count}")
    print(
        "rail_stage_positive(E,Q,EmeetQ)/84="
        f"({rail_stage_positive['E']},{rail_stage_positive['Q']},"
        f"{rail_stage_positive['packet']})"
    )
    print(
        f"present_packet_positive={present_packet_positive}/588 "
        f"half_stage_positive(E,Q,EmeetQ)/588="
        f"({e_only_half_positive},{q_only_half_positive},"
        f"{prefuture_positive_entries})"
    )
    print(
        f"h3_relaxed_positive(E_only,Q_only)/588="
        f"({e_only_h3_positive},{q_only_h3_positive}) "
        f"h3_relaxed_totals=({e_only_h3_total},{q_only_h3_total})"
    )
    print(
        f"prefuture_positive_entries={prefuture_positive_entries}/588 "
        f"prefuture_weighted_length={prefuture_weighted_length}"
    )
    print(
        f"positive_entries={positive_entries}/588 "
        f"positive_cell_count_hist={tuple(sorted(positive_cell_counts.items()))}"
    )
    print(f"empty_cells={empty_cells}")
    print(f"zero_entry_count={len(zero_entries)} zero_entries={zero_entries}")
    print(
        f"h3_slice_common_content={common_content} "
        f"baseline_h3_content={baseline_content} "
        f"slice_unit_cells={len(slice_units)}/84"
    )
    print(f"common_total_numerator={total_numerator}")
    print(
        f"off_diagonal_digit_totals={tuple(digit_totals)} "
        f"global_digit_support={tuple(i for i, value in enumerate(digit_totals) if value)}"
    )
    print(
        f"cell_digit_support_hist={tuple(sorted(digit_support_hist.items()))} "
        f"common_future_digits={common_future_digits}"
    )
    print(
        "q3_half_atlas(epsilon,probe:before,after,total)="
        f"{tuple((key, value) for key, value in sorted(half_atlas.items()))}"
    )
    print(
        "neighbor_tooth_equal_traces(piece,value)/588="
        f"({twin_piece_equal_count},{twin_value_equal_count})"
    )
    print(
        "alternate_h7_closure(before,after,total)="
        f"({sum(before for before, _ in alternate_closure)},"
        f"{sum(value > 0 for _, value in alternate_closure)},"
        f"{sum(value for _, value in alternate_closure)})"
    )
    print("typed verdict: the fixed zero-address clock-two / h=3 diagonal is empty")
    print("non-verdict: endpoint Fourier survival and quotient equivariance remain open")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
