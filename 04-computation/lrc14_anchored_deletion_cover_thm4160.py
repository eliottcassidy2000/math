#!/usr/bin/env python3
"""Exact primary certificate for THM-4160's anchored deletion cover.

The finite universe is every positive integer newcomer q outside the fixed
30-label pool.  The bounded-discrepancy estimate proved in THM-4160 reduces
the exact singleton, pair, and triple candidate scans to q <= 4389, 4683,
and 4919 respectively.  All geometry below uses fractions.Fraction.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from math import comb, lcm


Q = Fraction
DELTA = Q(1, 14)
THRESHOLD = Q(4, 63)
SAFE_DENSITY = Q(6, 7)
PRIMITIVE_MIN = Q(-3, 49)
PRIMITIVE_MAX = Q(3, 49)
PRIMITIVE_OSCILLATION = PRIMITIVE_MAX - PRIMITIVE_MIN

ANCHORS = (120, 126, 143)
POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63,
    80, 84, 85, 88, 95, 120, 126, 132, 143, 145,
    168, 170, 176, 190, 193, 240, 252, 264, 286, 290,
)
OPTIONAL = tuple(value for value in POOL if value not in ANCHORS)

EXPECTED_SINGLETONS = (5, 66, 182, 298, 336, 340, 380, 386, 528, 572)
EXPECTED_PAIR_CANDIDATES = EXPECTED_SINGLETONS + (235,)
EXPECTED_TRIPLE_CANDIDATES = EXPECTED_PAIR_CANDIDATES + (355,)
EXPECTED_REPAIR_SETS = {
    5: (85, 88, 95, 145, 168, 176, 193, 240, 252, 264, 290),
    66: (42, 80, 85, 88, 95, 145, 168, 170, 176, 193, 240, 252, 264, 286, 290),
    182: (85, 95, 145, 168, 176, 193, 240, 252, 290),
    298: (85, 88, 95, 145, 168, 176, 193, 240, 252, 264, 290),
    336: (85, 88, 95, 145, 168, 176, 193, 252, 264),
    340: (85, 88, 95, 145, 168, 240, 252, 264),
    380: (8, 16, 85, 88, 95, 132, 145, 168, 170, 176, 193, 240, 252, 264, 286, 290),
    386: (8, 16, 42, 85, 88, 95, 132, 145, 168, 170, 176, 190, 193, 240, 252, 264, 286, 290),
    528: (8, 85, 88, 95, 145, 168, 170, 193, 240, 252, 286, 290),
    572: (85, 88, 95, 145, 168, 176, 193, 240, 252, 264, 290),
    235: (85, 88, 145, 168, 193, 240, 252),
    355: (85, 88, 95, 193, 240, 252),
}


def require(predicate: bool, label: object) -> None:
    if not predicate:
        raise RuntimeError(f"requirement failed: {label}")


def qfloor(value: Q) -> int:
    return value.numerator // value.denominator


def measure(intervals: tuple[tuple[Q, Q], ...]) -> Q:
    return sum((right - left for left, right in intervals), Q(0))


def intersect_one(
    intervals: tuple[tuple[Q, Q], ...], speed: int
) -> tuple[tuple[Q, Q], ...]:
    """Intersect a closed interval union with ||speed*x|| >= 1/14."""
    result: list[tuple[Q, Q]] = []
    for left, right in intervals:
        lo = max(0, qfloor(speed * left) - 1)
        hi = min(speed - 1, qfloor(speed * right) + 1)
        for tooth in range(lo, hi + 1):
            a = max(left, (Q(tooth) + DELTA) / speed)
            b = min(right, (Q(tooth + 1) - DELTA) / speed)
            if a <= b:
                if result and result[-1][1] == a:
                    result[-1] = (result[-1][0], b)
                else:
                    result.append((a, b))
    return tuple(result)


def components(speeds: tuple[int, ...]) -> tuple[tuple[Q, Q], ...]:
    result = ((Q(0), Q(1)),)
    for speed in speeds:
        result = intersect_one(result, speed)
    return result


def safe_prefix_on_unit_tooth(value: Q) -> Q:
    if value <= DELTA:
        return Q(0)
    if value <= 1 - DELTA:
        return value - DELTA
    return SAFE_DENSITY


def safe_cumulative(speed: int, value: Q) -> Q:
    scaled = speed * value
    whole = qfloor(scaled)
    fractional = scaled - whole
    return (SAFE_DENSITY * whole + safe_prefix_on_unit_tooth(fractional)) / speed


def intersect_measure(intervals: tuple[tuple[Q, Q], ...], speed: int) -> Q:
    return sum(
        (safe_cumulative(speed, right) - safe_cumulative(speed, left)
         for left, right in intervals),
        Q(0),
    )


class PrimitiveEvaluator:
    """Integer-coordinate form of the exact safe-comb primitive identity."""

    def __init__(self, intervals: tuple[tuple[Q, Q], ...]) -> None:
        common = 1
        for left, right in intervals:
            common = lcm(common, left.denominator, right.denominator)
        self.common = common
        self.mass_numerator = sum(
            right.numerator * (common // right.denominator)
            - left.numerator * (common // left.denominator)
            for left, right in intervals
        )
        endpoints = []
        for left, right in intervals:
            endpoints.append((-1, left.numerator % left.denominator,
                              left.denominator, common // left.denominator))
            endpoints.append((+1, right.numerator % right.denominator,
                              right.denominator, common // right.denominator))
        self.endpoints = tuple(endpoints)

    @staticmethod
    def centered_error(residue: int, denominator: int) -> int:
        phase = 14 * residue
        if phase <= denominator:
            return -12 * residue
        if phase >= 13 * denominator:
            return 12 * (denominator - residue)
        return 2 * residue - denominator

    def scaled_measure(self, speed: int) -> int:
        error = 0
        for sign, numerator, denominator, weight in self.endpoints:
            residue = (speed * numerator) % denominator
            error += sign * self.centered_error(residue, denominator) * weight
        return 12 * speed * self.mass_numerator + error

    def comparison(self, speed: int) -> int:
        # Measure denominator is 14*speed*common.  Comparing with 4/63
        # cancels a factor seven and gives 9*N versus 8*speed*common.
        return 9 * self.scaled_measure(speed) - 8 * speed * self.common

    def exact_measure(self, speed: int) -> Q:
        return Q(self.scaled_measure(speed), 14 * speed * self.common)


def repair_row(
    evaluators: dict[int, PrimitiveEvaluator],
    bounds: dict[int, Q],
    newcomer: int,
) -> tuple[tuple[int, Q], ...]:
    repaired = []
    for removed in OPTIONAL:
        if newcomer > qfloor(bounds[removed]):
            continue
        comparison = evaluators[removed].comparison(newcomer)
        if comparison >= 0:
            repaired.append((removed, evaluators[removed].exact_measure(newcomer)))
    return tuple(repaired)


def exact_multi_repair_rows(
    leave_one: dict[int, tuple[tuple[Q, Q], ...]],
    candidates: tuple[int, ...],
    size: int,
    necessary_masks: dict[tuple[int, ...], frozenset[int]],
) -> dict[tuple[int, ...], tuple[tuple[int, Q], ...]]:
    rows: dict[tuple[int, ...], tuple[tuple[int, Q], ...]] = {}
    for values in combinations(candidates, size):
        possible = set(OPTIONAL)
        for proper_size in range(1, size):
            for subvalues in combinations(values, proper_size):
                possible.intersection_update(necessary_masks[subvalues])
        repaired = []
        for removed in sorted(possible):
            intervals = leave_one[removed]
            for newcomer in values:
                intervals = intersect_one(intervals, newcomer)
            mass = measure(intervals)
            if mass >= THRESHOLD:
                repaired.append((removed, mass))
        rows[values] = tuple(repaired)
    return rows


def fraction_text(value: Q) -> str:
    return f"{value.numerator}/{value.denominator}"


def main() -> None:
    require(len(POOL) == 30 and len(OPTIONAL) == 27, "fixed universes")
    require(set(ANCHORS).isdisjoint(OPTIONAL), "anchor/optional split")
    require(PRIMITIVE_OSCILLATION == Q(6, 49), "primitive oscillation")

    base = components(POOL)
    base_mass = measure(base)
    require((len(base), base_mass) ==
            (150, Q(298133356159, 4560289854120)), "THM-4156 base control")

    leave_one = {
        removed: components(tuple(value for value in POOL if value != removed))
        for removed in OPTIONAL
    }
    bounds: dict[int, Q] = {}
    bound_rows = []
    for removed in OPTIONAL:
        mass = measure(leave_one[removed])
        component_count = len(leave_one[removed])
        denominator = THRESHOLD - SAFE_DENSITY * mass
        require(denominator > 0, ("positive discrepancy denominator", removed))
        bound = PRIMITIVE_OSCILLATION * component_count / denominator
        bounds[removed] = bound
        bound_rows.append((bound, removed, mass, component_count, denominator))
    bound_rows.sort(key=lambda row: (row[0], row[1]), reverse=True)

    expected_rank_owners = (
        252, 85, 88, 95, 168, 145, 193, 240, 264, 176, 290, 286,
        8, 170, 16, 42, 80, 132, 190, 63, 10, 15, 40, 84, 60, 30, 20,
    )
    require(tuple(row[1] for row in bound_rows) == expected_rank_owners,
            "complete discrepancy-bound ordering")
    require(bounds[20] == bounds[30] == bounds[60], "only bound tie")
    require(len({row[0] for row in bound_rows}) == 25, "no other bound tie")

    arity_cutoffs = {
        arity: qfloor(bound_rows[8 - arity][0]) for arity in range(1, 9)
    }
    require(arity_cutoffs == {
        1: 4389, 2: 4683, 3: 4919, 4: 4940,
        5: 5417, 6: 5640, 7: 5671, 8: 27816,
    }, "rank cutoffs")

    evaluators = {
        removed: PrimitiveEvaluator(leave_one[removed]) for removed in OPTIONAL
    }

    histograms = {1: Counter(), 2: Counter(), 3: Counter()}
    masks: dict[tuple[int, ...], frozenset[int]] = {}
    masses: dict[int, tuple[tuple[int, Q], ...]] = {}
    equality_count = 0
    for newcomer in range(1, arity_cutoffs[3] + 1):
        if newcomer in POOL:
            continue
        repaired = repair_row(evaluators, bounds, newcomer)
        masses[newcomer] = repaired
        masks[(newcomer,)] = frozenset(removed for removed, _ in repaired)
        equality_count += sum(mass == THRESHOLD for _, mass in repaired)
        for arity in (1, 2, 3):
            if newcomer <= arity_cutoffs[arity]:
                histograms[arity][len(repaired)] += 1

    singleton_candidates = tuple(
        newcomer for newcomer in range(1, arity_cutoffs[1] + 1)
        if newcomer not in POOL and len(masses[newcomer]) >= 8
    )
    pair_candidates = tuple(
        newcomer for newcomer in range(1, arity_cutoffs[2] + 1)
        if newcomer not in POOL and len(masses[newcomer]) >= 7
    )
    triple_candidates = tuple(
        newcomer for newcomer in range(1, arity_cutoffs[3] + 1)
        if newcomer not in POOL and len(masses[newcomer]) >= 6
    )
    require(singleton_candidates == EXPECTED_SINGLETONS, "singleton classification")
    require(pair_candidates == tuple(sorted(EXPECTED_PAIR_CANDIDATES)),
            "pair-member classification")
    require(triple_candidates == tuple(sorted(EXPECTED_TRIPLE_CANDIDATES)),
            "triple-member classification")
    require(histograms[1] == Counter({
        0: 3790, 1: 537, 2: 8, 3: 4, 4: 5, 5: 3, 6: 1, 7: 1,
        8: 1, 9: 2, 11: 3, 12: 1, 15: 1, 16: 1, 18: 1,
    }), "singleton-universe histogram")
    require(histograms[2] == Counter({
        0: 4081, 1: 540, 2: 8, 3: 4, 4: 5, 5: 3, 6: 1, 7: 1,
        8: 1, 9: 2, 11: 3, 12: 1, 15: 1, 16: 1, 18: 1,
    }), "pair-member-universe histogram")
    require(histograms[3] == Counter({
        0: 4317, 1: 540, 2: 8, 3: 4, 4: 5, 5: 3, 6: 1, 7: 1,
        8: 1, 9: 2, 11: 3, 12: 1, 15: 1, 16: 1, 18: 1,
    }), "triple-member-universe histogram")
    for newcomer, expected in EXPECTED_REPAIR_SETS.items():
        require(tuple(removed for removed, _ in masses[newcomer]) == expected,
                ("repair set", newcomer))
    require(equality_count == 0, "no qualifying equality in q<=4919 scan")

    for newcomer, removed in ((386, 190), (340, 264), (235, 193), (4390, 240)):
        primitive_mass = evaluators[removed].exact_measure(newcomer)
        cumulative_mass = intersect_measure(leave_one[removed], newcomer)
        require(primitive_mass == cumulative_mass,
                ("primitive/cumulative control", newcomer, removed))

    # A full-pool addition would repair every deletion, so the singleton
    # cutoff 4389 is already a global cutoff for this stricter test.
    base_evaluator = PrimitiveEvaluator(base)
    full_pool_additions = tuple(
        newcomer for newcomer in range(1, arity_cutoffs[1] + 1)
        if newcomer not in POOL and base_evaluator.comparison(newcomer) >= 0
    )
    require(not full_pool_additions, "no full-pool addition, equality included")

    pair_rows = exact_multi_repair_rows(
        leave_one, pair_candidates, 2, masks
    )
    for values, repaired in pair_rows.items():
        masks[values] = frozenset(removed for removed, _ in repaired)
    pair_histogram = Counter(len(repaired) for repaired in pair_rows.values())
    require(pair_histogram == Counter({1: 28, 0: 27}), "pair repair histogram")
    require(not any(
        mass == THRESHOLD for repaired in pair_rows.values() for _, mass in repaired
    ), "no pair equality")
    require(all(
        not repaired or tuple(removed for removed, _ in repaired) == (252,)
        for repaired in pair_rows.values()
    ), "unique pair repair owner")
    require(not any(len(repaired) >= 7 for repaired in pair_rows.values()),
            "no forced pair")

    # Add the singleton masks for all triple candidates, then exact-audit all
    # 220 triples.  Pair masks are necessary filters only; surviving measures
    # are recomputed from the leave-one interval union.
    all_triple_pairs = exact_multi_repair_rows(
        leave_one, triple_candidates, 2, masks
    )
    for values, repaired in all_triple_pairs.items():
        masks[values] = frozenset(removed for removed, _ in repaired)
    triple_rows = exact_multi_repair_rows(
        leave_one, triple_candidates, 3, masks
    )
    triple_histogram = Counter(len(repaired) for repaired in triple_rows.values())
    require(triple_histogram == Counter({0: 220}), "triple repair histogram")
    require(not any(len(repaired) >= 6 for repaired in triple_rows.values()),
            "no forced triple")

    strongest = max(
        EXPECTED_SINGLETONS,
        key=lambda newcomer: (len(masses[newcomer]), -newcomer),
    )
    require(strongest == 386 and len(masses[strongest]) == 18,
            "strongest singleton")
    strongest_min = min(masses[strongest], key=lambda row: row[1])
    require(strongest_min ==
            (190, Q(144877112923, 2280144927060)), "strongest minimum repair")

    boundary = 4390
    boundary_repairs = repair_row(evaluators, bounds, boundary)
    require(not boundary_repairs, "q=4390 direct hostile")

    nonzero_pairs = tuple(
        values for values, repaired in sorted(pair_rows.items()) if repaired
    )
    semantic = {
        "anchors": ANCHORS,
        "pool": POOL,
        "base": (len(base), fraction_text(base_mass)),
        "bounds": [
            (removed, fraction_text(bound), qfloor(bound), fraction_text(mass),
             component_count, fraction_text(denominator))
            for bound, removed, mass, component_count, denominator in bound_rows
        ],
        "cutoffs": arity_cutoffs,
        "histograms": {
            arity: sorted(histograms[arity].items()) for arity in (1, 2, 3)
        },
        "repair_sets": {
            newcomer: tuple(removed for removed, _ in masses[newcomer])
            for newcomer in EXPECTED_REPAIR_SETS
        },
        "pair_histogram": sorted(pair_histogram.items()),
        "nonzero_pairs": nonzero_pairs,
        "triple_histogram": sorted(triple_histogram.items()),
        "family": (comb(27, 8), 10 * comb(27, 7)),
    }

    print("LRC14_ANCHORED_DELETION_COVER_THM4160_PRIMARY_20260826")
    print(f"anchors={ANCHORS};pool_size={len(POOL)};optional_size={len(OPTIONAL)}")
    print(f"base_components={len(base)};base_mass={base_mass}")
    print("safe_comb_mean=6/7;primitive_range=[-3/49,3/49];oscillation=6/49")
    print("bound_rows_descending=(removed,B,floor,mass,components,denominator)")
    for bound, removed, mass, component_count, denominator in bound_rows:
        print(f"  {removed},{bound},{qfloor(bound)},{mass},{component_count},{denominator}")
    print(f"arity_cutoffs={arity_cutoffs}")
    print("candidate_universes="
          f"singleton:1<=q<=4389,q_notin_P,count={4389-len(POOL)};"
          f"pair_member:1<=q<=4683,q_notin_P,count={4683-len(POOL)};"
          f"triple_member:1<=q<=4919,q_notin_P,count={4919-len(POOL)}")
    print(f"repair_histogram_singleton_universe={tuple(sorted(histograms[1].items()))}")
    print(f"repair_histogram_pair_member_universe={tuple(sorted(histograms[2].items()))}")
    print(f"repair_histogram_triple_member_universe={tuple(sorted(histograms[3].items()))}")
    print(f"qualifying_repair_equalities_q_le_4919={equality_count}")
    print(f"singleton_candidates={singleton_candidates}")
    for newcomer in singleton_candidates:
        repaired = masses[newcomer]
        minimum = min(repaired, key=lambda row: row[1])
        maximum = max(repaired, key=lambda row: row[1])
        print(
            f"singleton={newcomer};repair_count={len(repaired)};"
            f"repair_set={tuple(removed for removed, _ in repaired)};"
            f"min_repair={minimum};min_margin={minimum[1]-THRESHOLD};"
            f"max_repair={maximum};equalities=0"
        )
    print(f"strongest_singleton={strongest};strongest_repair_count=18;"
          f"minimum_repair={strongest_min};minimum_margin={strongest_min[1]-THRESHOLD}")
    print(f"hostile_q=235;repair_count={len(masses[235])};"
          f"repair_set={tuple(removed for removed, _ in masses[235])}")
    print("analytic_boundary_q=4390;at_most_seven_repairs_by_eighth_bound;"
          f"direct_repair_count={len(boundary_repairs)}")
    print(f"full_pool_additions_ge_threshold={full_pool_additions};"
          "full_pool_equalities=0")
    print(f"pair_member_candidates={pair_candidates};pair_histogram={tuple(sorted(pair_histogram.items()))}")
    print(f"nonzero_pairs={nonzero_pairs};all_nonzero_pair_repair_owner=252;"
          "admitted_pairs_at_threshold_7=()")
    print(f"triple_member_candidates={triple_candidates};"
          f"triple_histogram={tuple(sorted(triple_histogram.items()))};"
          "admitted_triples_at_threshold_6=()")
    print("pigeonhole_rule=for |Q|=t, |K|=8-t and |D_Q|>8-t imply D_Q\\K nonempty")
    print("singleton_omission_minima=" + str(tuple(
        (newcomer, len(masses[newcomer]) - 7) for newcomer in singleton_candidates
    )))
    print(f"newcomers_disjoint_pool={set(singleton_candidates).isdisjoint(POOL)};"
          "common_dilation_preserves_distinct_body_labels=True;"
          "doubled_body_labels_even_and_odd_tails_disjoint=True")
    print(f"old_family={comb(27,8)};new_per_singleton={comb(27,7)};"
          f"new_family={len(singleton_candidates)*comb(27,7)};"
          f"total_family={comb(27,8)+len(singleton_candidates)*comb(27,7)}")
    print("semantic_sha256=" + sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest())


if __name__ == "__main__":
    main()
