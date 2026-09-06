#!/usr/bin/env python3
"""Exact verifier for the sparse-interval reduction of THM-4409.

The audit checks three identities/inequalities independently of the max-flow
implementation:

* each pair-component family and each third-sheet family is 6-separated;
* edgewise loads min(|I|,|J|) obey all vertex capacities;
* the gap from that meet envelope to physical overlap is the sum of the two
  directed-exclusive hinge losses ``min(|I minus J|,|J minus I|)``.
"""

from fractions import Fraction as Q
from itertools import combinations, permutations
from math import gcd
import sys

import lrc14_third_sheet_component_network_thm4409 as core


class AuditError(RuntimeError):
    pass


CHECKS = 0


def need(ok, data):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise AuditError(data)


def length(interval):
    return interval[1] - interval[0]


def fmt(q):
    return f"{q.numerator}/{q.denominator}"


def contacts(left, right):
    result = []
    i = j = 0
    while i < len(left) and j < len(right):
        left_right = left[i][1]
        right_right = right[j][1]
        lo = max(left[i][0], right[j][0])
        hi = min(left_right, right_right)
        if lo < hi:
            result.append((i, j, hi - lo))
        if left_right <= right_right:
            i += 1
        if right_right <= left_right:
            j += 1
    return result


def overlap_length(first, second):
    return max(Q(0), min(first[1], second[1]) - max(first[0], second[0]))


def pair_components_with_parents(first, second):
    result = []
    i = j = 0
    while i < len(first) and j < len(second):
        first_right, second_right = first[i][1], second[j][1]
        lo = max(first[i][0], second[j][0])
        hi = min(first_right, second_right)
        if lo < hi:
            result.append(((lo, hi), i, j))
        if first_right <= second_right:
            i += 1
        if second_right <= first_right:
            j += 1
    return result


def audit_family(pieces, scale, tag):
    for index, item in enumerate(pieces):
        need(length(item) <= scale, ("length", tag, index, item, scale))
        if index:
            gap = item[0] - pieces[index - 1][1]
            need(gap >= 6 * scale, ("gap", tag, index, gap, scale))


def audit_graph(left, right, tag):
    edge_data = contacts(left, right)
    left_load = [Q(0)] * len(left)
    right_load = [Q(0)] * len(right)
    envelope = exact = hinge = Q(0)
    max_left_degree = max_right_degree = 0
    left_degree = [0] * len(left)
    right_degree = [0] * len(right)
    for i, j, overlap in edge_data:
        li, lj = length(left[i]), length(right[j])
        load = min(li, lj)
        exclusive_left = li - overlap
        exclusive_right = lj - overlap
        local_hinge = min(exclusive_left, exclusive_right)
        need(load - overlap == local_hinge,
             ("local-hinge", tag, i, j, li, lj, overlap, local_hinge))
        left_load[i] += load
        right_load[j] += load
        left_degree[i] += 1
        right_degree[j] += 1
        envelope += load
        exact += overlap
        hinge += local_hinge
    for i, load in enumerate(left_load):
        need(load <= length(left[i]), ("left-capacity", tag, i, load, length(left[i])))
    for j, load in enumerate(right_load):
        need(load <= length(right[j]), ("right-capacity", tag, j, load, length(right[j])))
    need(envelope - exact == hinge, ("global-hinge", tag, envelope, exact, hinge))
    if left_degree:
        max_left_degree = max(left_degree)
    if right_degree:
        max_right_degree = max(right_degree)
    return envelope, exact, hinge, len(edge_data), max_left_degree, max_right_degree


def triple_pair(w, pair):
    first, second = pair
    third = 3 - first - second
    envelope = exact = hinge = Q(0)
    edges = 0
    max_left_degree = max_right_degree = 0
    pair_scale = min(Q(1, 7 * w[first]), Q(1, 7 * w[second]))
    third_scale = Q(1, 7 * w[third])
    for assignment in permutations((0, 1, 2)):
        first_pieces = core.sheet_intervals(w[first], assignment[first])
        second_pieces = core.sheet_intervals(w[second], assignment[second])
        parented = pair_components_with_parents(first_pieces, second_pieces)
        pair_pieces = tuple(item[0] for item in parented)
        need(pair_pieces == core.intersect(first_pieces, second_pieces),
             ("parented-pair", w, pair, assignment))
        third_pieces = core.sheet_intervals(w[third], assignment[third])
        audit_family(pair_pieces, pair_scale, (w, pair, assignment, "pair"))
        audit_family(third_pieces, third_scale, (w, pair, assignment, "third"))
        data = audit_graph(pair_pieces, third_pieces, (w, pair, assignment))
        for pair_index, third_index, triple_overlap in contacts(pair_pieces, third_pieces):
            pair_piece, first_index, second_index = parented[pair_index]
            cross_first = overlap_length(first_pieces[first_index], third_pieces[third_index])
            cross_second = overlap_length(second_pieces[second_index], third_pieces[third_index])
            need(
                triple_overlap == min(length(pair_piece), cross_first, cross_second),
                (
                    "quantitative-helly",
                    w,
                    pair,
                    assignment,
                    pair_index,
                    third_index,
                    triple_overlap,
                    length(pair_piece),
                    cross_first,
                    cross_second,
                ),
            )
        envelope += data[0]
        exact += data[1]
        hinge += data[2]
        edges += data[3]
        max_left_degree = max(max_left_degree, data[4])
        max_right_degree = max(max_right_degree, data[5])
    need(envelope - exact == hinge, ("pair-hinge", w, pair))
    return envelope, exact, hinge, edges, max_left_degree, max_right_degree


def sharp_separation_controls():
    # For every kappa<2, choose epsilon small enough that one interval meets
    # two kappa-separated unit intervals but has capacity below the two
    # edgewise unit minima.  The kappa=3/2 rational instance is explicit.
    left = ((Q(9, 10), Q(13, 5)),)  # length 17/10
    right = ((Q(0), Q(1)), (Q(5, 2), Q(7, 2)))  # unit lengths, gap 3/2
    edge_data = contacts(left, right)
    need(len(edge_data) == 2, ("hostile-contact", edge_data))
    edge_sum = sum((min(length(left[i]), length(right[j])) for i, j, _ in edge_data), Q(0))
    need(edge_sum == 2 and edge_sum > length(left[0]),
         ("hostile-competition", edge_sum, length(left[0])))

    # At separation exactly 2, positive contact with two unit intervals forces
    # the spanning interval to have length strictly above 2.
    control_left = ((Q(9, 10), Q(31, 10)),)  # length 11/5
    control_right = ((Q(0), Q(1)), (Q(3), Q(4)))
    control_edges = contacts(control_left, control_right)
    need(len(control_edges) == 2, ("threshold-contact", control_edges))
    control_sum = sum(
        (min(length(control_left[i]), length(control_right[j])) for i, j, _ in control_edges),
        Q(0),
    )
    need(control_sum == 2 < length(control_left[0]),
         ("threshold-feasible", control_sum, length(control_left[0])))
    return edge_sum, length(left[0]), control_sum, length(control_left[0])


def main(height=79):
    sys.stdout.reconfigure(newline="\n")
    controls = sharp_separation_controls()
    values = tuple(v for v in range(1, height + 1, 2) if v % 3)
    triples = tuple(
        w for w in combinations(values, 3)
        if gcd(gcd(w[0], w[1]), w[2]) == 1
    )
    pair_choices = ((0, 1), (0, 2), (1, 2))
    max_left_degree = max_right_degree = 0
    total_edges = crossing_pairs = exact_pairs = 0
    named = {}
    pair_successes = {pair: 0 for pair in pair_choices}
    best_pair_counts = {pair: 0 for pair in pair_choices}
    selected_exact = 0
    strict_physical = equality_physical = above_physical = 0
    target_rows = []
    for w in triples:
        rows = []
        for pair in pair_choices:
            data = triple_pair(w, pair)
            envelope, exact, hinge, edges, dl, dr = data
            rows.append((envelope, pair, data))
            if envelope <= core.TARGET:
                pair_successes[pair] += 1
            total_edges += edges
            crossing_pairs += hinge > 0
            exact_pairs += hinge == 0
            max_left_degree = max(max_left_degree, dl)
            max_right_degree = max(max_right_degree, dr)
            need(envelope >= exact, ("upper", w, pair, envelope, exact))
            if w in ((1, 5, 11), (1, 19, 79)):
                named[w, pair] = data
        need(len({row[2][1] for row in rows}) == 1, ("pair-independent-physical", w, rows))
        physical = rows[0][2][1]
        best = min(rows, key=lambda row: (row[0], row[1]))
        best_pair_counts[best[1]] += 1
        if best[0] == physical:
            selected_exact += 1
        if physical < core.TARGET:
            strict_physical += 1
        elif physical == core.TARGET:
            equality_physical += 1
        else:
            above_physical += 1
        if best[0] == core.TARGET:
            target_rows.append((w, best[1]))

    if height == 79:
        need(pair_successes == {(0, 1): 2818, (0, 2): 2855, (1, 2): 2859},
             ("thm4409-pair-successes", pair_successes))
        need(best_pair_counts == {(0, 1): 400, (0, 2): 533, (1, 2): 1977},
             ("thm4409-best-pairs", best_pair_counts))
        need(selected_exact == 1747, ("thm4409-selected-exact", selected_exact))
        need((strict_physical, equality_physical, above_physical) == (2909, 1, 0),
             ("thm4409-physical-trichotomy", strict_physical, equality_physical, above_physical))
        need(target_rows == [((1, 5, 11), (0, 1))], ("thm4409-target-row", target_rows))

    print("LRC14 SIX-SEPARATED CONTACT CAPACITY COLLAPSE -- THM-4414")
    print("status=PASS PROVED_ANALYTICALLY+VERIFIED_EXACT; LRC14=OPEN")
    print("height", height, "values", len(values), "triples", len(triples), "pair_rows", 3 * len(triples))
    print("checks", CHECKS, "total_contact_edges", total_edges)
    print("pair_rows_exact/crossing", exact_pairs, crossing_pairs)
    print("maximum_left/right_degree", max_left_degree, max_right_degree)
    print("pair_successes", sorted(pair_successes.items()))
    print("best_pair_counts", sorted(best_pair_counts.items()))
    print("selected_envelope_exact", selected_exact)
    print("physical_strict/equality/above", strict_physical, equality_physical, above_physical)
    print("target_rows", target_rows)
    print("separation_3/2_hostile edge_sum/capacity", fmt(controls[0]), fmt(controls[1]))
    print("separation_2_control edge_sum/capacity", fmt(controls[2]), fmt(controls[3]))
    for key in sorted(named):
        envelope, exact, hinge, edges, dl, dr = named[key]
        print(
            "named", key[0], key[1],
            "envelope", fmt(envelope), "exact", fmt(exact), "hinge", fmt(hinge),
            "edges", edges, "degrees", dl, dr,
        )
    print("verdict=PASS")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 79)
