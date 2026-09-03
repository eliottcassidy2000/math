#!/usr/bin/env python3
"""Exact audit of the THM-4386 raw cross-product carrier formula.

The map n mod Z*w -> w cross n is audited against direct physical combs.
Every check remains active under optimized Python.
"""

from fractions import Fraction as F
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import gcd
from pathlib import Path


HERE = Path(__file__).resolve().parent
SPEC = spec_from_file_location(
    "coverage_referee", HERE / "lrc14_relation_incidence_referee_thm4386.py"
)
BASE = module_from_spec(SPEC)
SPEC.loader.exec_module(BASE)

R = F(3, 14)


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def raw_cross(w, n):
    return (
        w[1] * n[2] - w[2] * n[1],
        w[2] * n[0] - w[0] * n[2],
        w[0] * n[1] - w[1] * n[0],
    )


def carrier_length(w, carrier):
    terms = [F(2) * R / speed for speed in w]
    for i, j, k in ((0, 1, 2), (0, 2, 1), (1, 2, 0)):
        terms.append(
            R / w[i] + R / w[j]
            - F(abs(carrier[k]), w[i] * w[j])
        )
    return max(F(0), min(terms))


def strict_integer_bound(x):
    check(x > 0, ("nonpositive-bound", x))
    return (x.numerator - 1) // x.denominator


def lattice_formula(w):
    bounds = (
        strict_integer_bound(R * (w[1] + w[2])),
        strict_integer_bound(R * (w[0] + w[2])),
        strict_integer_bound(R * (w[0] + w[1])),
    )
    rows = {}
    for c0 in range(-bounds[0], bounds[0] + 1):
        for c1 in range(-bounds[1], bounds[1] + 1):
            numerator = -(w[0] * c0 + w[1] * c1)
            if numerator % w[2]:
                continue
            c2 = numerator // w[2]
            if abs(c2) > bounds[2]:
                continue
            carrier = (c0, c1, c2)
            if any(c % 3 == 0 for c in carrier):
                continue
            length = carrier_length(w, carrier)
            if length > 0:
                rows[carrier] = length
    return rows


def direct_components(w):
    positive, negative = BASE.physical_intervals(w)
    rows = {}
    for intervals in (positive, negative):
        for lo, hi in intervals:
            x = (lo + hi) / 2
            y = BASE.mod_one(3 * x)
            n = tuple(BASE.nearest_integer(speed * y) for speed in w)
            carrier = raw_cross(w, n)
            check(carrier not in rows,
                  ("duplicate-raw-carrier", w, carrier, rows.get(carrier), (lo, hi)))
            rows[carrier] = 3 * (hi - lo)
    return rows


def cross_snf_invariants(w):
    matrix = (
        (0, -w[2], w[1]),
        (w[2], 0, -w[0]),
        (-w[1], w[0], 0),
    )
    entries = [abs(x) for row in matrix for x in row]
    d1 = 0
    for x in entries:
        d1 = gcd(d1, x)
    minors = []
    for rows in combinations(range(3), 2):
        for cols in combinations(range(3), 2):
            minors.append(abs(
                matrix[rows[0]][cols[0]] * matrix[rows[1]][cols[1]]
                - matrix[rows[0]][cols[1]] * matrix[rows[1]][cols[0]]
            ))
    determinantal_two = 0
    for x in minors:
        determinantal_two = gcd(determinantal_two, x)
    check(d1 > 0 and determinantal_two % d1 == 0,
          ("snf-divisibility", w, d1, determinantal_two))
    return d1, determinantal_two // d1, 0


def audit_one(w):
    check(gcd(gcd(w[0], w[1]), w[2]) == 1, ("not-primitive", w))
    check(cross_snf_invariants(w) == (1, 1, 0), ("cross-snf", w, cross_snf_invariants(w)))
    direct = direct_components(w)
    lattice = lattice_formula(w)
    check(direct == lattice, ("carrier-formula", w, direct, lattice))
    direct_mu, _ = BASE.physical_measure(w)
    check(sum(lattice.values(), F(0)) == direct_mu,
          ("carrier-sum", w, sum(lattice.values(), F(0)), direct_mu))
    for carrier in direct:
        check(sum(a * b for a, b in zip(carrier, w)) == 0,
              ("orthogonal", w, carrier))
        check(all(c % 3 for c in carrier), ("owner-gate", w, carrier))
    return len(direct)


def main():
    tested = 0
    components = 0
    for w in BASE.primitive_triples(79):
        components += audit_one(w)
        tested += 1
    high_controls = (
        (1, 11, 121), (1, 11, 175), (5, 53, 55),
        (35, 55, 77), (67, 131, 199), (191, 193, 199),
    )
    for w in high_controls:
        components += audit_one(w)
        tested += 1
    print(f"RAW_CARRIER_FORMULA exact_triples={tested} base_components={components}")
    print("SNF diag(1,1,0) on every control")
    for w in ((1, 5, 11), (1, 11, 25), (1, 11, 175), (35, 55, 77)):
        print(f"CONTROL {w} rows={lattice_formula(w)} sum={sum(lattice_formula(w).values(),F(0))}")
    print("PASS")


if __name__ == "__main__":
    main()
