#!/usr/bin/env python3
"""Exact scratch verifier for the signed-(1,4,1) scale-three triple comb.

This is deliberately self-contained.  It compares the determinant/quadrature
formula with two literal physical-sheet implementations and certifies the
finite endpoint list left by the analytic large-product bound.
"""

from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
from math import gcd


RADIUS = Fraction(1, 14)
R = Fraction(3, 14)
SHARP = Fraction(12, 301)


def check(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def mod_one(x):
    return x - x.numerator // x.denominator


def circle_distance(x):
    x = mod_one(x)
    return min(x, 1 - x)


@lru_cache(maxsize=None)
def danger_intervals(w, sheet):
    radius = Fraction(1, 14 * w)
    pieces = []
    for k in range(w):
        center = mod_one(Fraction(k, w) - Fraction(sheet, 3))
        left, right = center - radius, center + radius
        if left < 0:
            pieces.extend(((Fraction(0), right), (1 + left, Fraction(1))))
        elif right > 1:
            pieces.extend(((left, Fraction(1)), (Fraction(0), right - 1)))
        else:
            pieces.append((left, right))
    pieces.sort()
    merged = []
    for left, right in pieces:
        if merged and left <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], right))
        else:
            merged.append((left, right))
    return tuple(merged)


def intersect(first, second):
    out = []
    i = j = 0
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            out.append((left, right))
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return tuple(out)


def interval_measure(intervals):
    return sum((right - left for left, right in intervals), Fraction(0))


def failure_measure_direct(triple):
    a, b, c = triple
    pos = intersect(intersect(danger_intervals(a, 0), danger_intervals(b, 1)),
                    danger_intervals(c, 2))
    neg = intersect(intersect(danger_intervals(a, 0), danger_intervals(b, 2)),
                    danger_intervals(c, 1))
    return 3 * (interval_measure(pos) + interval_measure(neg))


def owner(w, x):
    hits = [j for j in range(3)
            if circle_distance(w * (x + Fraction(j, 3))) < RADIUS]
    check(len(hits) <= 1, ("multiple owners", w, x, hits))
    return hits[0] if hits else None


@lru_cache(maxsize=None)
def third_period_walls(w):
    period = Fraction(1, 3)
    ans = {Fraction(0), period}
    for j in range(3):
        for k in range(w):
            for sign in (-1, 1):
                x = Fraction(k, w) + sign * Fraction(1, 14 * w) - Fraction(j, 3)
                q = x / period
                ans.add(x - (q.numerator // q.denominator) * period)
    return frozenset(ans)


def failure_measure_cells(triple):
    cuts = sorted(set().union(*(third_period_walls(w) for w in triple)))
    total = Fraction(0)
    for left, right in zip(cuts, cuts[1:]):
        if left == right:
            continue
        owners = [owner(w, (left + right) / 2) for w in triple]
        if None not in owners and len(set(owners)) == 3:
            total += right - left
    return 3 * total


def nearest_integer_nonnegative(x):
    return (x + Fraction(1, 2)).numerator // (x + Fraction(1, 2)).denominator


def check_owner_determinant_cells(p, q, b, branch):
    """Check the unquotiented three-lift owner criterion cell by cell."""
    triple = tuple(sorted((p, b, q)))
    cuts = sorted(set().union(*(third_period_walls(w) for w in triple)))
    count = 0
    for left, right in zip(cuts, cuts[1:]):
        if left == right:
            continue
        x = (left + right) / 2
        y = 3 * x
        n_p = nearest_integer_nonnegative(p * y)
        n_q = nearest_integer_nonnegative(q * y)
        n_b = nearest_integer_nonnegative(b * y)
        e_p = p * y - n_p
        e_q = q * y - n_q
        e_b = b * y - n_b
        endpoint_eligible = abs(e_p) < R and abs(e_q) < R
        K = q * n_p - p * n_q
        determinant_criterion = endpoint_eligible and K % 4 == 0 and (K // 4) % 3 != 0
        owners = [owner(w, x) for w in triple]
        physical_failure = None not in owners and len(set(owners)) == 3
        check(physical_failure == determinant_criterion,
              ("owner/determinant mismatch", p, q, b, branch, x, owners, K,
               e_p, e_q, e_b))
        if physical_failure:
            if branch == "mean":
                defect = 4 * n_b - n_p - n_q
            else:
                defect = 4 * n_b - n_q + n_p
            check(defect == 0, ("nonzero middle defect", triple, x, defect))
        count += 1
    return count


def positive_part(x):
    return max(Fraction(0), x)


def quadrature_error(t):
    """E(t)=sum_{k>=1,3\nmid k}(t-k)_+-t^2/3, in closed form."""
    rho = t - 3 * (t.numerator // (3 * t.denominator))
    if rho <= 1:
        return -rho * rho / 3
    if rho <= 2:
        return rho - 1 - rho * rho / 3
    return 2 * rho - 3 - rho * rho / 3


def hinge_parameters(p, q):
    return Fraction(3 * (q - p), 56), Fraction(3 * (q + p), 56)


def component_length(p, q, k):
    A, B = hinge_parameters(p, q)
    return Fraction(4, p * q) * (positive_part(B - k) - positive_part(A - k))


def failure_measure_sum(p, q):
    _, B = hinge_parameters(p, q)
    last = B.numerator // B.denominator
    return 2 * sum((component_length(p, q, k)
                    for k in range(1, last + 1) if k % 3), Fraction(0))


def failure_measure_formula(p, q):
    A, B = hinge_parameters(p, q)
    return Fraction(3, 98) + Fraction(8, p * q) * (
        quadrature_error(B) - quadrature_error(A)
    )


def endpoint_orientation(p, q):
    """The unique signed-(1,4,1) branch for coefficient-one p<q, if any."""
    hits = []
    if (p + q) % 4 == 0:
        b = (p + q) // 4
        if b > 0 and b % 2 and b % 3 and b not in (p, q):
            hits.append(("mean", b))
    if (q - p) % 4 == 0:
        b = (q - p) // 4
        if b > 0 and b % 2 and b % 3 and b not in (p, q):
            hits.append(("difference", b))
    check(len(hits) <= 1, ("multiple (1,4,1) orientations", p, q, hits))
    return hits[0] if hits else None


def has_signed_relation(triple, coefficient):
    for middle in triple:
        endpoints = tuple(x for x in triple if x != middle)
        if coefficient * middle in (sum(endpoints), abs(endpoints[0] - endpoints[1])):
            return True
    return False


def primitive_endpoint_pairs(limit, product_limit=None):
    for p in range(1, limit + 1, 2):
        if p % 3 == 0:
            continue
        for q in range(p + 2, limit + 1, 2):
            if q % 3 == 0 or gcd(p, q) != 1:
                continue
            if product_limit is not None and p * q >= product_limit:
                continue
            oriented = endpoint_orientation(p, q)
            if oriented:
                branch, b = oriented
                yield p, q, b, branch


def main():
    formula_checks = 0
    cell_checks = 0
    owner_cell_checks = 0
    overlap_121 = []
    for p, q, b, branch in primitive_endpoint_pairs(249):
        triple = tuple(sorted((p, b, q)))
        direct = failure_measure_direct(triple)
        summed = failure_measure_sum(p, q)
        formula = failure_measure_formula(p, q)
        check(direct == summed == formula,
              ("formula mismatch", p, q, b, branch, direct, summed, formula))
        check(has_signed_relation(triple, 4), ("missing relation", triple))
        if has_signed_relation(triple, 2):
            overlap_121.append(triple)
        for permuted in (triple, (triple[1], triple[0], triple[2]),
                         (triple[2], triple[1], triple[0])):
            check(failure_measure_direct(permuted) == direct,
                  ("permutation mismatch", triple, permuted))
        if q <= 61:
            check(failure_measure_cells(triple) == direct,
                  ("cell mismatch", triple, direct))
            owner_cell_checks += check_owner_determinant_cells(p, q, b, branch)
            cell_checks += 1
        formula_checks += 1
    check(not overlap_121, ("overlap with (1,2,1)", overlap_121))

    small = []
    for p, q, b, branch in primitive_endpoint_pairs(288, product_limit=289):
        value = failure_measure_formula(p, q)
        triple = tuple(sorted((p, b, q)))
        check(failure_measure_direct(triple) == value,
              ("small direct mismatch", triple, value))
        small.append((value, p, q, b, branch, triple))
    small.sort(reverse=True)
    check(len(small) == 31, ("small-pair count", len(small)))
    check(small[0] == (SHARP, 1, 43, 11, "mean", (1, 11, 43)),
          ("wrong maximum", small[0]))
    check(all(value < SHARP for value, *_ in small[1:]),
          ("nonunique maximum", small[:5]))

    coarse_289 = Fraction(3, 98) + Fraction(8, 3 * 289)
    check(coarse_289 < SHARP, ("large-product bound", coarse_289, SHARP))
    coarse_gap = SHARP - coarse_289
    check(coarse_gap == Fraction(101, 3653538), ("coarse gap", coarse_gap))

    # Common dilation is a physical Haar conjugacy, not an LRC time evolution.
    for triple in ((1, 11, 43), (5, 13, 47), (5, 11, 49), (1, 5, 23)):
        base = failure_measure_direct(triple)
        for d in (5, 7, 11):
            if d % 3:
                check(failure_measure_direct(tuple(d * x for x in triple)) == base,
                      ("dilation mismatch", triple, d))

    # Same coefficient-one endpoints but different physical middle speeds:
    # endpoint data alone do not determine the comb outside the relation.
    middle_hostile = {
        (1, 5, 43): failure_measure_direct((1, 5, 43)),
        (1, 11, 43): failure_measure_direct((1, 11, 43)),
    }
    check(middle_hostile == {(1, 5, 43): Fraction(6, 301),
                             (1, 11, 43): Fraction(12, 301)},
          ("middle hostile", middle_hostile))

    # Lexicographically first primitive positive comb outside both relation sectors.
    first_outside = None
    pool = [w for w in range(1, 24, 2) if w % 3]
    for triple in combinations(pool, 3):
        if gcd(gcd(triple[0], triple[1]), triple[2]) != 1:
            continue
        if has_signed_relation(triple, 2) or has_signed_relation(triple, 4):
            continue
        value = failure_measure_direct(triple)
        if value:
            first_outside = (triple, value)
            break
    check(first_outside == ((1, 5, 23), Fraction(2, 161)),
          ("first outside hostile", first_outside))

    residual_rows = 0
    residual_max = (Fraction(-1), None)
    pool_199 = [w for w in range(1, 200, 2) if w % 3]
    for triple in combinations(pool_199, 3):
        if gcd(gcd(triple[0], triple[1]), triple[2]) != 1:
            continue
        if has_signed_relation(triple, 2) or has_signed_relation(triple, 4):
            continue
        value = failure_measure_direct(triple)
        residual_rows += 1
        residual_max = max(residual_max, (value, triple))
    check(residual_max == (Fraction(12, 371), (1, 11, 53)),
          ("height-199 residual maximum", residual_max))

    print("SIGNED-(1,4,1) SCALE-THREE EXACT AUDIT")
    print(f"formula/direct/permutation checks through endpoint 249: {formula_checks}")
    print(f"independent wall-cell checks through q=61: {cell_checks}")
    print(f"literal owner/determinant cells through q=61: {owner_cell_checks}")
    print("signed-(1,2,1) overlap in audited unit sector: 0")
    print(f"large-product bound at pq=289: {coarse_289}")
    print(f"gap to 12/301: {coarse_gap}")
    print(f"small endpoint pairs pq<289: {len(small)}")
    print("SMALL PAIRS, DESCENDING EXACT MEASURE")
    for value, p, q, b, branch, triple in small:
        print(f"p={p} q={q} b={b} branch={branch} triple={triple} measure={value}")
    print(f"unique maximum: triple={small[0][-1]} measure={small[0][0]}")
    print(f"same-endpoint hostile: {middle_hostile}")
    print(f"first positive comb outside (1,2,1)/(1,4,1): {first_outside}")
    print(f"height-199 residual after both sectors: rows={residual_rows} max={residual_max}")
    print("PASS")


if __name__ == "__main__":
    main()
