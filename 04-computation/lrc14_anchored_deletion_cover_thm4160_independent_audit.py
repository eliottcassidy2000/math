#!/usr/bin/env python3
"""Independent exact referee for the THM-4160 deletion-cover hierarchy.

This deliberately does not import the candidate implementation.  It builds the
27 leave-one safe sets from one global rational wall arrangement, evaluates
comb intersections with a centered primitive in repair-specific integer
coordinates, and then finds repairwise pair graphs and their triangles.
"""

from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import floor, lcm


DELTA = F(1, 14)
DENSITY = F(6, 7)
OSCILLATION = F(6, 49)
THRESHOLD = F(4, 63)
ANCHORS = (120, 126, 143)
POOL = (
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
    120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252,
    264, 286, 290,
)
OPTIONAL = tuple(v for v in POOL if v not in ANCHORS)
POOL_SET = frozenset(POOL)
MAX_Q = 27_816


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def safe_at(x, speed):
    """Exact closed safe predicate ||speed*x|| >= 1/14."""
    residue = (speed * x.numerator) % x.denominator
    return 14 * residue >= x.denominator and 14 * residue <= 13 * x.denominator


def merge_cell(bank, left, right):
    if left == right:
        return
    if bank and bank[-1][1] == left:
        bank[-1] = (bank[-1][0], right)
    else:
        bank.append((left, right))


def global_wall_leave_one_banks():
    """Build all U_r from a common wall-cell classification, not intersections."""
    walls = {F(0), F(1)}
    for speed in POOL:
        for tooth in range(speed):
            walls.add(F(14 * tooth + 1, 14 * speed))
            walls.add(F(14 * tooth + 13, 14 * speed))
    walls = sorted(walls)
    banks = {removed: [] for removed in OPTIONAL}
    base = []
    for left, right in zip(walls, walls[1:]):
        midpoint = (left + right) / 2
        failures = tuple(speed for speed in POOL if not safe_at(midpoint, speed))
        if not failures:
            merge_cell(base, left, right)
            for removed in OPTIONAL:
                merge_cell(banks[removed], left, right)
        elif len(failures) == 1 and failures[0] in banks:
            merge_cell(banks[failures[0]], left, right)
    return tuple(walls), tuple(base), {r: tuple(bank) for r, bank in banks.items()}


def mass(bank):
    return sum((right - left for left, right in bank), F(0))


def intersect_comb(bank, speed):
    """Exact interval clipping, used only after the independent wall build."""
    answer = []
    for left, right in bank:
        low = max(0, floor(speed * left) - 1)
        high = min(speed - 1, floor(speed * right) + 1)
        for tooth in range(low, high + 1):
            a = max(left, F(14 * tooth + 1, 14 * speed))
            b = min(right, F(14 * tooth + 13, 14 * speed))
            if a <= b:
                merge_cell(answer, a, b)
    return tuple(answer)


class PrimitiveEvaluator:
    """Exact centered-primitive evaluator for mu(bank intersect G_q)."""

    def __init__(self, bank):
        self.bank = tuple(bank)
        self.components = len(bank)
        common = 1
        for left, right in bank:
            common = lcm(common, left.denominator, right.denominator)
        self.common = common
        self.mass_numerator = sum(
            right.numerator * (common // right.denominator)
            - left.numerator * (common // left.denominator)
            for left, right in bank
        )
        endpoints = []
        for left, right in bank:
            endpoints.append((-1, left.numerator % left.denominator,
                              left.denominator, common // left.denominator))
            endpoints.append((+1, right.numerator % right.denominator,
                              right.denominator, common // right.denominator))
        self.endpoints = tuple(endpoints)

    @staticmethod
    def centered_error(residue, denominator):
        phase = 14 * residue
        if phase <= denominator:
            return -12 * residue
        if phase >= 13 * denominator:
            return 12 * (denominator - residue)
        return 2 * residue - denominator

    def comparison(self, q):
        error = 0
        for sign, numerator, denominator, weight in self.endpoints:
            residue = (q * numerator) % denominator
            error += sign * self.centered_error(residue, denominator) * weight
        scaled_measure = 12 * q * self.mass_numerator + error
        return 9 * scaled_measure - 8 * q * self.common

    def exact_measure(self, q):
        error = 0
        for sign, numerator, denominator, weight in self.endpoints:
            residue = (q * numerator) % denominator
            error += sign * self.centered_error(residue, denominator) * weight
        scaled_measure = 12 * q * self.mass_numerator + error
        return F(scaled_measure, 14 * q * self.common)

    def discrepancy_bound(self):
        deficit = THRESHOLD - DENSITY * F(self.mass_numerator, self.common)
        require(deficit > 0, "finite discrepancy denominator must be positive")
        return OSCILLATION * self.components / deficit


def direct_measure(bank, speed):
    return mass(intersect_comb(bank, speed))


def classify_singletons(evaluators, bounds):
    masks = [0] * (MAX_Q + 1)
    equality_rows = []
    for bit, removed in enumerate(OPTIONAL):
        evaluator = evaluators[removed]
        cutoff = floor(bounds[removed])
        for q in range(1, cutoff + 1):
            if q in POOL_SET:
                continue
            comparison = evaluator.comparison(q)
            if comparison >= 0:
                masks[q] |= 1 << bit
            if comparison == 0:
                equality_rows.append((q, removed))
    return masks, tuple(equality_rows)


def pair_graphs(banks, masks):
    candidates_by_repair = {}
    edges_by_repair = {}
    equality_rows = []
    tested = 0
    discrepancy_pruned = 0
    for bit, removed in enumerate(OPTIONAL):
        candidates = tuple(q for q in range(1, MAX_Q + 1) if masks[q] >> bit & 1)
        candidates_by_repair[removed] = candidates
        edges = []
        base = banks[removed]
        for index, q1 in enumerate(candidates):
            if index + 1 == len(candidates):
                break
            first_bank = intersect_comb(base, q1)
            first_eval = PrimitiveEvaluator(first_bank)
            cutoff = floor(first_eval.discrepancy_bound())
            for q2 in candidates[index + 1:]:
                if q2 > cutoff:
                    discrepancy_pruned += 1
                    continue
                tested += 1
                comparison = first_eval.comparison(q2)
                if comparison >= 0:
                    edges.append((q1, q2))
                if comparison == 0:
                    equality_rows.append((removed, q1, q2))
        edges_by_repair[removed] = tuple(edges)
    return (candidates_by_repair, edges_by_repair, tuple(equality_rows),
            tested, discrepancy_pruned)


def triangles_and_exact_triples(banks, edges_by_repair):
    triangles_by_repair = {}
    surviving = []
    equality_rows = []
    for removed in OPTIONAL:
        edges = frozenset(edges_by_repair[removed])
        vertices = sorted({q for edge in edges for q in edge})
        triangles = []
        for a, b, c in combinations(vertices, 3):
            if (a, b) in edges and (a, c) in edges and (b, c) in edges:
                triangles.append((a, b, c))
        triangles_by_repair[removed] = tuple(triangles)
        for a, b, c in triangles:
            bank_ab = intersect_comb(intersect_comb(banks[removed], a), b)
            comparison = PrimitiveEvaluator(bank_ab).comparison(c)
            if comparison >= 0:
                surviving.append((removed, a, b, c))
            if comparison == 0:
                equality_rows.append((removed, a, b, c))
    return triangles_by_repair, tuple(surviving), tuple(equality_rows)


def mask_labels(mask):
    return tuple(removed for bit, removed in enumerate(OPTIONAL) if mask >> bit & 1)


def main():
    print("LRC14_HIGHER_ARITY_INDEPENDENT_FRACTION_REFEREE_20260825")
    print(f"universe=q positive, q notin P, q<=B_(rank 1);P={POOL};A={ANCHORS};O={OPTIONAL}")
    walls, base, banks = global_wall_leave_one_banks()
    require(len(walls) == 7_134, "THM-4156 wall count changed")
    require(len(base) == 150, "THM-4156 component count changed")
    require(mass(base) == F(298133356159, 4560289854120), "THM-4156 mass changed")
    evaluators = {r: PrimitiveEvaluator(banks[r]) for r in OPTIONAL}
    bounds = {r: evaluators[r].discrepancy_bound() for r in OPTIONAL}
    ranked = sorted(OPTIONAL, key=lambda r: bounds[r], reverse=True)
    expected_ranked = (252, 85, 88, 95, 168, 145, 193, 240, 264, 176, 290,
                       286, 8, 170, 16, 42, 80, 132, 190, 63, 10, 15, 40, 84,
                       20, 30, 60)
    require(tuple(ranked) == expected_ranked, "discrepancy ranking changed")
    print(f"walls={len(walls)};base_components={len(base)};base_mass={mass(base)}")
    for rank, removed in enumerate(ranked, 1):
        print(f"bound_rank={rank};repair={removed};B={bounds[removed]};floor={floor(bounds[removed])};"
              f"components={len(banks[removed])};mass={mass(banks[removed])}")

    cutoffs = {}
    for size in range(1, 9):
        needed = 9 - size
        owner = ranked[needed - 1]
        cutoffs[size] = floor(bounds[owner])
        print(f"rank_cutoff=size{size};needed={needed};owner={owner};"
              f"B={bounds[owner]};floor={cutoffs[size]}")
    require(cutoffs == {1: 4389, 2: 4683, 3: 4919, 4: 4940,
                        5: 5417, 6: 5640, 7: 5671, 8: 27816},
            "arity cutoffs changed")

    masks, singleton_equalities = classify_singletons(evaluators, bounds)
    rows = tuple((q, masks[q]) for q in range(1, MAX_Q + 1) if masks[q])
    histogram = Counter(mask.bit_count() for q, mask in rows)
    zero_count = MAX_Q - len(POOL) - len(rows)
    histogram[0] = zero_count
    semantic = "".join(f"{q}:{','.join(map(str, mask_labels(mask)))};" for q, mask in rows)
    print(f"singleton_positive_rows={len(rows)};hist={tuple(sorted(histogram.items()))};"
          f"equality_rows={singleton_equalities}")
    print(f"singleton_mask_semantic_sha256={sha256(semantic.encode()).hexdigest()}")
    multi_rows = tuple((q, mask_labels(mask)) for q, mask in rows if mask.bit_count() >= 2)
    print(f"singleton_multi_rows={multi_rows}")
    qualifying = tuple((q, mask_labels(masks[q])) for q in range(1, cutoffs[1] + 1)
                       if q not in POOL_SET and masks[q].bit_count() >= 8)
    print(f"singleton_qualifying={qualifying}")
    require(tuple(q for q, _ in qualifying) ==
            (5, 66, 182, 298, 336, 340, 380, 386, 528, 572),
            "singleton qualifying list changed")
    require(not singleton_equalities, "unexpected singleton equality")
    require(len(rows) == 573, "positive singleton-mask row count changed")
    require(histogram == Counter({0: 27213, 1: 541, 2: 8, 3: 4, 4: 5,
                                  5: 3, 6: 1, 7: 1, 8: 1, 9: 2, 11: 3,
                                  12: 1, 15: 1, 16: 1, 18: 1}),
            "singleton histogram changed")

    # Primitive-vs-literal controls on both sides of the threshold.
    controls = ((386, 190), (340, 264), (235, 193), (27_816, 252),
                (27_817, 252), (4390, 240))
    for q, removed in controls:
        primitive = evaluators[removed].exact_measure(q)
        literal = direct_measure(banks[removed], q)
        require(primitive == literal, f"primitive/literal disagreement at {(q, removed)}")
        print(f"control=q{q};repair={removed};measure={literal};"
              f"comparison={'EQ' if literal == THRESHOLD else 'PASS' if literal > THRESHOLD else 'FAIL'}")
    require(evaluators[252].comparison(27_817) < 0, "post-cutoff hostile unexpectedly repairs")
    require(bounds[252] < 27_817, "largest analytic bound does not close the infinite tail")

    (candidates_by_repair, edges_by_repair, pair_equalities,
     pair_tests, pair_prunes) = pair_graphs(banks, masks)
    pair_masks = defaultdict(int)
    for bit, removed in enumerate(OPTIONAL):
        for pair in edges_by_repair[removed]:
            pair_masks[pair] |= 1 << bit
    pair_hist = Counter(mask.bit_count() for mask in pair_masks.values())
    max_pair_repairs = max(pair_hist, default=0)
    edge_counts = tuple((r, len(candidates_by_repair[r]), len(edges_by_repair[r]))
                        for r in OPTIONAL)
    pair_semantic = "".join(
        f"{a},{b}:{','.join(map(str, mask_labels(mask)))};"
        for (a, b), mask in sorted(pair_masks.items())
    )
    print(f"repairwise_candidate_edge_counts={edge_counts}")
    print(f"pair_exact_tests={pair_tests};pair_discrepancy_prunes={pair_prunes};"
          f"surviving_pairs={len(pair_masks)};pair_repair_hist={tuple(sorted(pair_hist.items()))};"
          f"max_pair_repairs={max_pair_repairs};equality_rows={pair_equalities}")
    print(f"pair_semantic_sha256={sha256(pair_semantic.encode()).hexdigest()}")
    print(f"pair_rows={tuple((pair, mask_labels(mask)) for pair, mask in sorted(pair_masks.items()))}")
    require(not pair_equalities, "unexpected pair equality")
    require(len(pair_masks) == 45 and max_pair_repairs == 1,
            "pair maximum or survivor count changed")
    require(set(mask_labels(mask) for mask in pair_masks.values()) == {(252,)},
            "a pair repairs outside r=252")

    triangles, triple_survivors, triple_equalities = triangles_and_exact_triples(
        banks, edges_by_repair)
    triangle_counts = tuple((r, len(triangles[r])) for r in OPTIONAL)
    print(f"repairwise_triangle_counts={triangle_counts};"
          f"triple_survivors={triple_survivors};equality_rows={triple_equalities}")
    require(not triple_survivors and not triple_equalities,
            "a repair survives three newcomers or lands on equality")

    pair_witness = (5, 66)
    witness_bank = intersect_comb(intersect_comb(banks[252], pair_witness[0]),
                                  pair_witness[1])
    witness_measure = mass(witness_bank)
    require(witness_measure >= THRESHOLD, "declared pair witness failed")
    print(f"pair_max_witness={pair_witness};repair=252;measure={witness_measure};"
          f"margin={witness_measure - THRESHOLD}")
    print("deletion_cover_maxima=size1:18@q386;size2:1@pair(5,66),repair252;"
          "sizes3to8:0")
    print("gate_audit=size1 exactly ten rows reach 8;size2 max1<7;"
          "size3 max0<6;size4 max0<5;size5 max0<4;size6 max0<3;"
          "size7 max0<2;size8 max0<1")
    print("infinite_tail=q>=27817 has zero singleton repairs by B_252<27817;"
          "monotonicity gives D(Q) subset intersection_q D({q})")
    print("PASS")


if __name__ == "__main__":
    main()
