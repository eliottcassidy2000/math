#!/usr/bin/env python3
"""Exact two-tournament and addressed KG(7,3) controls; stdlib only.
No assert statements: all verification remains enabled under python -O.
"""
import json
from fractions import Fraction
from itertools import combinations, permutations
from math import comb


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def multiply(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def identity(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def characteristic_by_newton(s):
    n = len(s)
    power = identity(n)
    traces = [n]
    for k in range(1, n+1):
        power = multiply(power, s)
        traces.append(sum(power[i][i] for i in range(n)))
    coefficients = [1]
    for k in range(1, n+1):
        numerator = -sum(coefficients[k-i]*traces[i] for i in range(1, k+1))
        require(numerator % k == 0, 'Newton coefficient not integral')
        coefficients.append(numerator // k)
    return coefficients, traces[1:]


def tournament(differences, expected_h, expected_aut, expected_coefficients):
    n = 7
    a = [[int((j-i) % n in differences) for j in range(n)] for i in range(n)]
    require(all(a[i][i] == 0 for i in range(n)), 'diagonal')
    require(all(a[i][j]+a[j][i] == 1 for i in range(n) for j in range(n) if i != j), 'not a tournament')
    require(all(sum(row) == 3 for row in a), 'not regular')
    per_start = [0]*n
    automorphisms = []
    for p in permutations(range(n)):
        if all(a[p[i]][p[i+1]] for i in range(n-1)):
            per_start[p[0]] += 1
        if all(a[i][j] == a[p[i]][p[j]] for i in range(n) for j in range(n)):
            automorphisms.append(p)
    require(sum(per_start) == expected_h, 'Hamiltonian path total')
    require(per_start == [expected_h//n]*n, 'per-start counts')
    require(len(automorphisms) == expected_aut, 'automorphism order')
    affine = {tuple((k*i+b) % n for i in range(n)) for k in differences for b in range(n)
              if all(a[i][j] == a[(k*i+b) % n][(k*j+b) % n] for i in range(n) for j in range(n))}
    require(set(automorphisms) == affine, 'unexpected nonaffine automorphism')
    reflection_anti = all(a[(-i) % n][(-j) % n] == a[j][i] for i in range(n) for j in range(n))
    require(reflection_anti, 'reflection is not anti')
    require(tuple((-i) % n for i in range(n)) not in automorphisms, 'reflection incorrectly an automorphism')
    s = [[a[i][j]-a[j][i] for j in range(n)] for i in range(n)]
    coefficients, traces = characteristic_by_newton(s)
    require(coefficients == expected_coefficients, 'characteristic polynomial')
    return dict(differences=sorted(differences), Hamiltonian_paths=expected_h,
                per_start=per_start, automorphism_order=len(automorphisms),
                reflection_anti=reflection_anti, skew_coefficients_descending=coefficients,
                skew_power_traces_1_to_7=traces)


def kneser():
    faces = [frozenset(t) for t in combinations(range(7), 3)]
    n = len(faces)
    require(n == 35, 'face count')
    k = [[int(not s & r) for r in faces] for s in faces]
    inverse = [[Fraction((-1)**len(s & r), 4*comb(3, len(s & r))) for r in faces] for s in faces]
    require(multiply(k, inverse) == identity(n), 'K times claimed inverse')
    require(multiply(inverse, k) == identity(n), 'claimed inverse times K')
    # Addressed four-set incidence is exactly K after complementing row addresses.
    universe = frozenset(range(7))
    require(all(k[i][j] == int(faces[j] <= universe-faces[i]) for i in range(n) for j in range(n)), 'four-face incidence')
    # Deterministic nonconstant rational field and its reconstructed values.
    field = [[Fraction(i*i-3*i+1, i+1)] for i in range(n)]
    require(multiply(inverse, multiply(k, field)) == field, 'field round trip')
    return dict(dimension=n, row_sum=sorted(set(map(sum,k))),
                inverse_by_intersection=[str(Fraction((-1)**t,4*comb(3,t))) for t in range(4)],
                checked_product_entries_each_direction=n*n,
                addressed_additive_field_round_trip=True)


def main():
    results = {
        'status':'FINITE-EXACT',
        'universe':'Two explicit labelled order-seven tournaments, all 7! vertex permutations each; all 35 addressed triples',
        'R7':tournament({1,2,3},175,7,[1,0,21,0,35,0,7,0]),
        'Paley7':tournament({1,2,4},189,21,[1,0,21,0,147,0,343,0]),
        'KG_7_3':kneser(),
        'scope':'No full order-seven tournament census; no unlabelled deck reconstruction; no Collatz-denominator equivalence.'
    }
    print(json.dumps(results, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
