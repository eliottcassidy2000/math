#!/usr/bin/env python3
"""Exact scout for the scale-three, three-tail owner-permutation comb.

For tail speeds prime to 3, each tail can spoil at most one of the three
physical lifts x, x+1/3, x+2/3.  This script computes the exact Haar measure
of the phases where three labelled tails spoil all three lifts.
"""

from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import gcd


RADIUS = Fraction(1, 14)


def mod_one(x: Fraction) -> Fraction:
    return x - (x.numerator // x.denominator)


def mod_period(x: Fraction, period: Fraction) -> Fraction:
    quotient = x / period
    return x - (quotient.numerator // quotient.denominator) * period


def circle_distance(x: Fraction) -> Fraction:
    x = mod_one(x)
    return min(x, 1 - x)


def owner(w: int, x: Fraction):
    hits = [j for j in range(3) if circle_distance(w * (x + Fraction(j, 3))) < RADIUS]
    assert len(hits) <= 1, (w, x, hits)
    return hits[0] if hits else None


@lru_cache(maxsize=None)
def walls(w: int):
    period = Fraction(1, 3)
    ans = {Fraction(0), period}
    # Work on one 1/3-period of the unlabelled three-sheet failure set.
    for j in range(3):
        for k in range(w):
            for sign in (-1, 1):
                x = Fraction(k, w) + sign * Fraction(1, 14 * w) - Fraction(j, 3)
                ans.add(mod_period(x, period))
    return frozenset(ans)


@lru_cache(maxsize=None)
def danger_intervals(w: int, sheet: int):
    """Half-open representatives; endpoints have measure zero."""
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


def intersect_intervals(first, second):
    answer = []
    i = j = 0
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            answer.append((left, right))
        if first[i][1] < second[j][1]:
            i += 1
        else:
            j += 1
    return tuple(answer)


def interval_measure(intervals):
    return sum((right - left for left, right in intervals), Fraction(0))


@lru_cache(maxsize=None)
def oriented_pair_intervals(a: int, b: int, b_sheet: int):
    """Cache the repeated first intersection in the height census."""
    return intersect_intervals(danger_intervals(a, 0), danger_intervals(b, b_sheet))


def has_signed_121_relation(triple):
    for doubled in triple:
        others = list(triple)
        others.remove(doubled)
        if 2 * doubled in (others[0] + others[1], abs(others[0] - others[1])):
            return True
    return False


def signed_121_endpoints(triple):
    """Return (p,q,b,branch) with b carrying coefficient 2, or None."""
    hits = []
    for b in triple:
        endpoints = sorted(w for w in triple if w != b)
        p, q = endpoints
        if 2 * b == p + q:
            hits.append((p, q, b, "mean"))
        if 2 * b == q - p:
            hits.append((p, q, b, "difference"))
    if not hits:
        return None
    if len(hits) != 1:
        raise RuntimeError(("nonunique signed-(1,2,1) orientation", triple, hits))
    return hits[0]


def positive_part(x: Fraction) -> Fraction:
    return max(Fraction(0), x)


def resonant_hinge_parameters(p: int, q: int):
    return Fraction(3 * (q - p), 28), Fraction(3 * (q + p), 28)


def resonant_component_length(p: int, q: int, k: int) -> Fraction:
    """Length of the y-circle component with endpoint determinant 2k."""
    A, B = resonant_hinge_parameters(p, q)
    return Fraction(2, p * q) * (
        positive_part(B - k) - positive_part(A - k)
    )


def resonant_failure_measure_formula(triple) -> Fraction:
    oriented = signed_121_endpoints(triple)
    if oriented is None:
        raise RuntimeError(("nonresonant formula request", triple))
    p, q, _, _ = oriented
    _, B = resonant_hinge_parameters(p, q)
    last = B.numerator // B.denominator
    return 2 * sum(
        (
            resonant_component_length(p, q, k)
            for k in range(1, last + 1)
            if k % 3
        ),
        Fraction(0),
    )


def two_of_three_quadrature_error(t: Fraction) -> Fraction:
    """Periodic error sum_(3 does not divide k)(t-k)_+ - t^2/3."""
    r = t - 3 * (t.numerator // (3 * t.denominator))
    if r < 1:
        return -r * r / 3
    if r < 2:
        return r - 1 - r * r / 3
    return 2 * r - 3 - r * r / 3


def resonant_failure_measure_quadrature(triple) -> Fraction:
    p, q, _, _ = signed_121_endpoints(triple)
    A, B = resonant_hinge_parameters(p, q)
    return Fraction(3, 49) + Fraction(4, p * q) * (
        two_of_three_quadrature_error(B) - two_of_three_quadrature_error(A)
    )


def small_resonant_endpoint_cases():
    answer = []
    for p in range(1, 81, 2):
        if p % 3 == 0:
            continue
        for q in range(p + 2, 81, 2):
            if p * q >= 81 or q % 3 == 0 or gcd(p, q) != 1:
                continue
            if q % 12 == p % 12:
                b = (p + q) // 2
            elif q % 12 == (-p) % 12:
                b = (q - p) // 2
            else:
                continue
            answer.append((p, q, b, tuple(sorted((p, b, q)))))
    return tuple(answer)


def failure_measure(triple):
    a, b, c = triple
    positive = intersect_intervals(
        oriented_pair_intervals(a, b, 1),
        danger_intervals(c, 2),
    )
    negative = intersect_intervals(
        oriented_pair_intervals(a, b, 2),
        danger_intervals(c, 1),
    )
    # The other four permutations are the two displayed orientations after a
    # common cyclic sheet shift, hence have equal Haar measure.
    return 3 * (interval_measure(positive) + interval_measure(negative)), 3 * (
        len(positive) + len(negative)
    )


def failure_measure_by_cells(triple):
    """Independent wall-cell implementation on one 1/3-period."""
    cuts = sorted(set().union(*(walls(w) for w in triple)))
    total = Fraction(0)
    cells = 0
    for left, right in zip(cuts, cuts[1:]):
        if left == right:
            continue
        own = [owner(w, (left + right) / 2) for w in triple]
        if None not in own and len(set(own)) == 3:
            total += right - left
            cells += 1
    return 3 * total, 3 * cells


def main():
    controls = [(1, 10, 11), (1, 11, 13), (1, 3, 9), (1, 5, 7)]
    print("CONTROLS")
    for triple in controls:
        if all(w % 3 for w in triple):
            measure, cells = failure_measure(triple)
            if (measure, cells) != failure_measure_by_cells(triple):
                raise RuntimeError(("control implementation mismatch", triple))
            print(f"{triple}: measure={measure} ~ {float(measure):.12f}; cells={cells}")
        else:
            print(f"{triple}: outside full-order-three regime")

    audit_pool = [w for w in range(1, 32, 2) if w % 3]
    audit_count = 0
    for triple in combinations(audit_pool, 3):
        if gcd(gcd(triple[0], triple[1]), triple[2]) != 1:
            continue
        direct = failure_measure(triple)
        independent = failure_measure_by_cells(triple)
        if direct != independent:
            raise RuntimeError(("implementation mismatch", triple, direct, independent))
        for permuted in ((triple[1], triple[0], triple[2]), (triple[2], triple[1], triple[0])):
            if failure_measure(permuted)[0] != direct[0]:
                raise RuntimeError(("label permutation mismatch", triple, permuted))
        dilated = tuple(5 * w for w in triple)
        if failure_measure(dilated)[0] != direct[0]:
            raise RuntimeError(("common dilation mismatch", triple, dilated))
        audit_count += 1
    print(f"INDEPENDENT cell/interval/permutation/dilation checks={audit_count}: PASS")

    pool = [w for w in range(1, 200, 2) if w % 3]
    records = []
    for triple in combinations(pool, 3):
        if gcd(gcd(*triple[:2]), triple[2]) != 1:
            continue
        measure, cells = failure_measure(triple)
        records.append((measure, triple, cells))
    records.sort(reverse=True)
    print(f"ENUM primitive odd nonmultiples of 3 through {pool[-1]}: {len(records)} triples")
    for measure, triple, cells in records[:30]:
        print(f"{triple}: measure={measure} ~ {float(measure):.12f}; cells={cells}")
    values = sorted({m for m, _, _ in records}, reverse=True)
    print("TOP DISTINCT MEASURES")
    for value in values[:20]:
        count = sum(m == value for m, _, _ in records)
        print(f"{value} ~ {float(value):.12f}: count={count}")

    resonant = [row for row in records if has_signed_121_relation(row[1])]
    nonresonant = [row for row in records if not has_signed_121_relation(row[1])]
    print("SIGNED-(1,2,1) RESONANCE SPLIT")
    print(f"resonant={len(resonant)} max={resonant[0][1]} measure={resonant[0][0]}")
    print(
        f"nonresonant={len(nonresonant)} max={nonresonant[0][1]} "
        f"measure={nonresonant[0][0]}"
    )
    if not all(has_signed_121_relation(triple) for _, triple, _ in records[:30]):
        raise RuntimeError("top-thirty resonance signal failed")
    print("top_30_all_signed_121_resonant=PASS")

    sharp_bound = Fraction(6, 77)
    equality_cases = []
    for measure, triple, _ in resonant:
        oriented = signed_121_endpoints(triple)
        p, q, b, branch = oriented
        if q % 12 not in (p % 12, (-p) % 12):
            raise RuntimeError(("endpoint congruence", triple, oriented))
        expected_branch = "mean" if q % 12 == p % 12 else "difference"
        if branch != expected_branch:
            raise RuntimeError(("endpoint branch", triple, oriented))
        if gcd(p, q) != 1:
            raise RuntimeError(("primitive endpoint gcd", triple, oriented))
        formula = resonant_failure_measure_formula(triple)
        quadrature = resonant_failure_measure_quadrature(triple)
        if measure != formula or formula != quadrature:
            raise RuntimeError(
                ("resonant formula mismatch", triple, measure, formula, quadrature)
            )
        coarse_bound = Fraction(3, 49) + Fraction(4, 3 * p * q)
        if measure > coarse_bound or measure > sharp_bound:
            raise RuntimeError(("resonant bound", triple, measure, coarse_bound))
        if measure == sharp_bound:
            equality_cases.append((triple, oriented))
    if equality_cases != [((1, 5, 11), (1, 11, 5, "difference"))]:
        raise RuntimeError(("resonant equality classification", equality_cases))
    print(
        "RESONANT ENDPOINT/QUADRATURE THEOREM "
        f"checks={len(resonant)} bound={sharp_bound} equality={equality_cases[0][0]}: PASS"
    )

    small_cases = small_resonant_endpoint_cases()
    expected_pairs = (
        (1, 11),
        (1, 13),
        (1, 23),
        (1, 25),
        (1, 35),
        (1, 37),
        (1, 47),
        (1, 49),
        (1, 59),
        (1, 61),
        (1, 71),
        (1, 73),
        (5, 7),
    )
    if tuple((row[0], row[1]) for row in small_cases) != expected_pairs:
        raise RuntimeError(("small endpoint case list", small_cases))
    print("RESONANT ENDPOINT CASES pq<81")
    for p, q, b, triple in small_cases:
        measure = resonant_failure_measure_formula(triple)
        print(f"p={p} q={q} b={b} triple={triple} measure={measure}")
    if max((resonant_failure_measure_formula(row[3]), row) for row in small_cases)[1][
        :2
    ] != (1, 11):
        raise RuntimeError("small endpoint maximum is not (1,11)")
    print("pq_ge_81_bound=3/49+4/(3pq)<6/77; pq_lt_81_cases=13: PASS")


if __name__ == "__main__":
    main()
