#!/usr/bin/env python3
"""Exact Fraction certificate for THM-4170.

The certificate builds the common wall arrangement of the THM-4156 pool,
classifies every open wall cell by the optional labels which fail there, and
reconstructs all 2,925 three-deletion banks.  It then verifies an explicit
eight-edge matching whose members are all eventually above the THM-4150 Haar
threshold.  No floating point arithmetic is used.
"""

from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import comb, floor
import sys


sys.stdout.reconfigure(newline="\n")


THRESHOLD = F(4, 63)
DENSITY = F(6, 7)
OSCILLATION = F(6, 49)
ANCHORS = (120, 126, 143)
POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
    120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252,
    264, 286, 290,
)
OPTIONAL = tuple(value for value in POOL if value not in ANCHORS)
INDEX = {value: index for index, value in enumerate(OPTIONAL)}
BAD = 1 << 31
MATCHING = (
    (8, 84, 252),
    (15, 88, 95),
    (40, 85, 170),
    (63, 176, 193),
    (10, 145, 290),
    (42, 132, 264),
    (16, 168, 286),
    (80, 190, 240),
)
EXPECTED_ROWS = {
    (8, 84, 252): (
        148, F(53623011809, 701583054480),
        F(4961689987, 2455540690680),
        F(44500410884160, 4961689987), 8969,
    ),
    (15, 88, 95): (
        178, F(140987531863, 1824115941648),
        F(17603497445, 6384405795768),
        F(139153987548576, 17603497445), 7905,
    ),
    (40, 85, 170): (
        166, F(4641131947, 59611632080),
        F(6087298409, 1877766410520),
        F(38168476426080, 6087298409), 6271,
    ),
    (63, 176, 193): (
        164, F(23261449, 302928780),
        F(7400521, 3180752190),
        F(63874697040, 7400521), 8632,
    ),
    (10, 145, 290): (
        162, F(12058235549, 157251374280),
        F(1229956807, 550379809980),
        F(10917738271440, 1229956807), 8877,
    ),
    (42, 132, 264): (
        156, F(354200884867, 4560289854120),
        F(49204909241, 15961014489420),
        F(304887950246880, 49204909241), 6197,
    ),
    (16, 168, 286): (
        168, F(4772688913, 62044759920),
        F(1591026937, 651469979160),
        F(13401668142720, 1591026937), 8424,
    ),
    (80, 190, 240): (
        164, F(13415915189, 175395763620),
        F(1270909207, 613885172670),
        F(12327816528720, 1270909207), 9700,
    ),
}


def require(condition, label):
    if not condition:
        raise AssertionError(label)


def ceil_fraction(value):
    return -((-value.numerator) // value.denominator)


def safe_at(point, speed):
    residue = (speed * point.numerator) % point.denominator
    return 14 * residue >= point.denominator and 14 * residue <= 13 * point.denominator


def mask_of(labels):
    answer = 0
    for label in labels:
        answer |= 1 << INDEX[label]
    return answer


def labels_of(mask):
    return tuple(label for index, label in enumerate(OPTIONAL) if mask >> index & 1)


def global_cells():
    walls = {F(0), F(1)}
    for speed in POOL:
        for tooth in range(speed):
            walls.add(F(14 * tooth + 1, 14 * speed))
            walls.add(F(14 * tooth + 13, 14 * speed))
    walls = tuple(sorted(walls))
    failures = []
    lengths = defaultdict(F)
    histogram = Counter()
    anchor_set = frozenset(ANCHORS)
    for left, right in zip(walls, walls[1:]):
        midpoint = (left + right) / 2
        failed = tuple(speed for speed in POOL if not safe_at(midpoint, speed))
        if any(speed in anchor_set for speed in failed) or len(failed) > 3:
            failure = BAD
            histogram[4] += 1
        else:
            failure = mask_of(failed)
            histogram[len(failed)] += 1
            lengths[failure] += right - left
        failures.append(failure)
    return walls, tuple(failures), dict(lengths), histogram


def bank_summary(edge, failures, lengths):
    measure = sum((length for failure, length in lengths.items()
                   if failure & ~edge == 0), F(0))
    components = 0
    previous = False
    for failure in failures:
        included = failure & ~edge == 0
        if included and not previous:
            components += 1
        previous = included
    surplus = DENSITY * measure - THRESHOLD
    bound = OSCILLATION * components / surplus if surplus > 0 else None
    activation = ceil_fraction(bound) if bound is not None else None
    return components, measure, surplus, bound, activation


def safe_prefix(point, speed):
    """Integral of 1_[1/14,13/14]({speed*x}) from zero to point."""
    phase = speed * point
    whole = floor(phase)
    residue = phase - whole
    partial = max(F(0), min(residue, F(13, 14)) - F(1, 14))
    return (F(6, 7) * whole + partial) / speed


def comb_failure_buckets(walls, failures, speed):
    buckets = defaultdict(F)
    previous = safe_prefix(walls[0], speed)
    for index, failure in enumerate(failures):
        current = safe_prefix(walls[index + 1], speed)
        if failure != BAD:
            buckets[failure] += current - previous
        previous = current
    return dict(buckets)


def direct_bank_intervals(edge, walls, failures):
    answer = []
    for index, failure in enumerate(failures):
        if failure & ~edge != 0:
            continue
        left, right = walls[index], walls[index + 1]
        if answer and answer[-1][1] == left:
            answer[-1] = (answer[-1][0], right)
        else:
            answer.append((left, right))
    return tuple(answer)


def direct_comb_mass(intervals, speed):
    return sum((safe_prefix(right, speed) - safe_prefix(left, speed)
                for left, right in intervals), F(0))


def exact_hypergraph(triples, buckets):
    edges = []
    equalities = []
    measures = {}
    for edge in triples:
        measure = sum((length for failure, length in buckets.items()
                       if failure & ~edge == 0), F(0))
        measures[edge] = measure
        if measure >= THRESHOLD:
            edges.append(edge)
        if measure == THRESHOLD:
            equalities.append(edge)
    return tuple(edges), tuple(equalities), measures


def main():
    print("LRC14_TRIPLE_DELETION_MATCHING_EVENTUAL_HAAR_THM4170_20260826")
    print(f"P={POOL};A={ANCHORS};O={OPTIONAL}")
    walls, failures, lengths, histogram = global_cells()
    require(len(walls) == 7134, "wall count")
    require(histogram == Counter({0: 150, 1: 328, 2: 518, 3: 678, 4: 5459}),
            "failure-size histogram")
    triples = tuple(mask_of(labels) for labels in combinations(OPTIONAL, 3))
    require(len(triples) == comb(27, 3) == 2925, "triple universe")

    summaries = {edge: bank_summary(edge, failures, lengths) for edge in triples}
    strict = tuple(edge for edge in triples if summaries[edge][2] > 0)
    equalities = tuple(edge for edge in triples if summaries[edge][2] == 0)
    require(len(strict) == 1335, "strict limit edge count")
    require(not equalities, "limit equality audit")

    matching_masks = tuple(mask_of(row) for row in MATCHING)
    used = 0
    for edge in matching_masks:
        require(edge & used == 0, "matching disjointness")
        used |= edge
        require(edge in strict, "matching edge strictness")
        require(summaries[edge] == EXPECTED_ROWS[labels_of(edge)],
                "matching exact row")
    require(used.bit_count() == 24, "matching covers 24 distinct labels")
    require(max(summaries[edge][4] for edge in matching_masks) == 9700,
            "matching activation cutoff")

    activation_edges = tuple(edge for edge in strict if summaries[edge][4] <= 9700)
    previous_edges = tuple(edge for edge in strict if summaries[edge][4] <= 9699)
    previous_levels = sorted({summaries[edge][4] for edge in strict
                              if summaries[edge][4] < 9700})
    require(previous_levels[-1] == 9687, "previous activation event")
    require(len(activation_edges) == 666, "activation edge count")
    require(len(previous_edges) < len(activation_edges), "activation event is nonempty")

    semantic = "".join(
        f"{labels_of(edge)}:{summaries[edge][0]}:{summaries[edge][1]}:"
        f"{summaries[edge][2]}:{summaries[edge][4]};"
        for edge in triples
    )
    print(f"walls={len(walls)};cells={len(failures)};cell_hist={tuple(sorted(histogram.items()))}")
    print(f"triple_banks={len(triples)};strict_limit_edges={len(strict)};limit_equalities={len(equalities)}")
    print(f"matching={MATCHING};matching_labels={used.bit_count()};body_size=7;"
          f"pigeonhole_slack={len(MATCHING)-7};bodies_per_q={comb(27,7)}")
    print(f"activation_edges_q9700={len(activation_edges)};previous_event={previous_levels[-1]};"
          f"activation_edges_previous={len(previous_edges)}")
    for labels, edge in zip(MATCHING, matching_masks):
        components, measure, surplus, bound, activation = summaries[edge]
        print(f"matching_edge={labels};components={components};bank_mass={measure};"
              f"limit_surplus={surplus};B={bound};ceil_B={activation}")

    for speed in (9699, 9700):
        buckets = comb_failure_buckets(walls, failures, speed)
        edges, exact_equalities, measures = exact_hypergraph(triples, buckets)
        direct_rows = []
        for labels, edge in zip(MATCHING, matching_masks):
            intervals = direct_bank_intervals(edge, walls, failures)
            direct = direct_comb_mass(intervals, speed)
            require(direct == measures[edge], "bucket/direct exact control")
            direct_rows.append((labels, direct, direct - THRESHOLD,
                                direct >= THRESHOLD))
        require(not exact_equalities, "boundary equality audit")
        require(len(edges) == {9699: 1351, 9700: 1356}[speed],
                "boundary hyperedge count")
        if speed == 9700:
            require(all(row[3] for row in direct_rows), "q9700 matching activation")
        print(f"q={speed};hyperedges={len(edges)};equalities={exact_equalities};"
              f"matching_rows={tuple(direct_rows)}")

    print(f"semantic_sha256={sha256(semantic.encode()).hexdigest()}")
    print("quantifier=every_7_subset_hits_at_most_7_pairwise_disjoint_repair_triples")
    print("cutoff_scope=9700_is_the_fixed_matching_discrepancy_cutoff_not_an_actual_tail_minimum")
    print("tail=all_integers_q>=9700;finite_q_below_9700_not_classified")
    print("PASS")


if __name__ == "__main__":
    main()
