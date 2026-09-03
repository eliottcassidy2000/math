#!/usr/bin/env python3
"""Exact sharp THM-4386 audit for the non-unit coefficient rays at l1 <= 14.

This script verifies the finite part of the proof after the monotone sampling
bound.  It is optimization-safe and compares the raw-carrier series against
an independent literal physical interval comb on every retained triple.
"""

from collections import defaultdict
from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import permutations, product
from math import gcd
from pathlib import Path


HERE = Path(__file__).resolve().parent
SPEC = spec_from_file_location(
    "coverage_referee", HERE / "lrc14_relation_incidence_referee_thm4386.py"
)
BASE = module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)

R = F(3, 14)
DATA = {
    (2, 5, 5): {
        "integral": F(81, 2450),
        "maximum": F(564, 20405),
        "maximizers": {(5, 53, 55)},
        "threshold": F(204050, 1333),
        "height": 153,
        "count": 208,
    },
    (2, 5, 7): {
        "integral": F(9, 343),
        "maximum": F(12, 539),
        "maximizers": {(7, 17, 77), (7, 53, 77)},
        "threshold": F(539, 3),
        "height": 179,
        "count": 442,
    },
    (4, 5, 5): {
        "integral": F(36, 1225),
        "maximum": F(444, 18179),
        "maximizers": {(5, 49, 53)},
        "threshold": F(64925, 366),
        "height": 177,
        "count": 252,
    },
}


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def canonical(v):
    d = gcd(gcd(abs(v[0]), abs(v[1])), abs(v[2]))
    check(d > 0, ("zero-vector", v))
    v = tuple(x // d for x in v)
    for x in v:
        if x:
            return v if x > 0 else tuple(-y for y in v)
    raise RuntimeError(("zero-after-normalization", v))


def coefficient_vectors(pattern):
    ans = set()
    for magnitudes in set(permutations(pattern)):
        for signs in product((-1, 1), repeat=3):
            v = canonical(tuple(a*b for a,b in zip(magnitudes, signs)))
            if min(v) < 0 < max(v):
                ans.add(v)
    return tuple(sorted(ans))


def strict_floor(x):
    return (x.numerator - 1) // x.denominator


def carrier_length(w, raw):
    terms = [2 * R / speed for speed in w]
    for i, j, k in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
        terms.append(R/w[i] + R/w[j] - F(abs(raw[k]), w[i]*w[j]))
    return max(F(0), min(terms))


def ray_measure(w, relation):
    cap = min(
        strict_floor(R * (w[j] + w[k]) / abs(relation[i]))
        for i, j, k in ((0, 1, 2), (1, 0, 2), (2, 0, 1))
    )
    return 2 * sum(
        (carrier_length(w, tuple(k*x for x in relation))
         for k in range(1, cap + 1) if k % 3),
        F(0),
    )


def box_slice_integral(pattern):
    # Shift the boxes [-a r,a r] to [0,2 a r] and use the elementary
    # three-box convolution formula at the central sum.
    total = sum(pattern)
    alternating = F(0)
    for mask in range(8):
        residual = total - 2 * sum(pattern[i] for i in range(3) if mask >> i & 1)
        if residual > 0:
            alternating += (-1) ** mask.bit_count() * residual * residual
    return R * R * alternating / (2 * pattern[0] * pattern[1] * pattern[2])


def enumerate_sector(pattern, height):
    units = tuple(w for w in range(1, height + 1, 2) if w % 3)
    presentations = defaultdict(set)
    for relation in coefficient_vectors(pattern):
        for a in units:
            for b in units:
                if a >= b:
                    continue
                numerator = -(relation[0] * a + relation[1] * b)
                if numerator % relation[2]:
                    continue
                c = numerator // relation[2]
                if c <= b or c > height or c % 2 == 0 or c % 3 == 0:
                    continue
                w = (a, b, c)
                if gcd(gcd(a, b), c) != 1:
                    continue
                check(sum(x*y for x,y in zip(relation,w)) == 0,
                      ("relation", pattern, relation, w))
                presentations[w].add(relation)
    return presentations


def main():
    total_triples = 0
    total_component_checks = 0
    summaries = []
    for pattern, expected in DATA.items():
        integral = box_slice_integral(pattern)
        check(integral == expected["integral"],
              ("slice-integral", pattern, integral, expected["integral"]))
        baseline = F(2, 3) * integral
        gap = expected["maximum"] - baseline
        check(gap > 0, ("nonpositive-gap", pattern, baseline, expected["maximum"]))
        threshold = F(6, 7) / gap
        check(threshold == expected["threshold"],
              ("height-threshold", pattern, threshold, expected["threshold"]))
        check(strict_floor(threshold) == expected["height"],
              ("threshold-floor", pattern, strict_floor(threshold), expected["height"]))

        presentations = enumerate_sector(pattern, expected["height"])
        check(len(presentations) == expected["count"],
              ("finite-count", pattern, len(presentations), expected["count"]))
        rows = []
        for w, relations in presentations.items():
            values = {ray_measure(w, relation) for relation in relations}
            check(len(values) == 1, ("presentation-values", pattern, w, relations, values))
            value = next(iter(values))
            direct, _ = BASE.physical_measure(w)
            check(value == direct, ("direct-comb", pattern, w, relations, value, direct))
            components = BASE.physical_component_data(w)
            for row in components:
                n = row[3]
                vanishing = [sum(a*q for a,q in zip(relation,n)) for relation in relations]
                check(vanishing == [0] * len(relations),
                      ("short-defect", pattern, w, n, relations, vanishing))
                total_component_checks += 1
            if len(relations) > 1:
                check(value == 0, ("independent-short-relations", pattern, w, relations, value))
            rows.append((value, w, tuple(sorted(relations))))
        rows.sort(reverse=True)
        peak = rows[0][0]
        maximizers = {row[1] for row in rows if row[0] == peak}
        check(peak == expected["maximum"],
              ("maximum", pattern, peak, expected["maximum"]))
        check(maximizers == expected["maximizers"],
              ("maximizers", pattern, maximizers, expected["maximizers"]))
        total_triples += len(rows)
        summaries.append((pattern, integral, baseline, gap, threshold, len(rows), peak,
                          tuple(sorted(maximizers)), tuple(rows[:5])))

    print("MONOTONE_BOUND mu <= (2/3) integral + 6/(7 w_max)")
    print(f"FINITE_TOTAL triples={total_triples} base_components_checked={total_component_checks}")
    for summary in summaries:
        print("SECTOR", summary)
    print("PASS")


if __name__ == "__main__":
    main()
