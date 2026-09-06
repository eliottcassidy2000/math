#!/usr/bin/env python3
"""Independent projective Hermite audit using literal symbolic jets and Smith.

Requires SymPy; no producer imports. The all-height proof is analytic.
"""
from itertools import combinations
from math import gcd, lcm, prod
from sympy import Matrix, Symbol, ZZ, expand
from sympy.matrices.normalforms import smith_normal_form


def check(ok, why):
    if not ok:
        raise RuntimeError(why)


def bracket(v, w):
    return v[0]*w[1]-v[1]*w[0]


def tangent(v):
    # Deliberately use bounded literal Bezout search rather than Euclid.
    bound = max(abs(z) for z in v) + 1
    for a in range(-bound, bound+1):
        for b in range(-bound, bound+1):
            if bracket(v, (a, b)) == 1:
                return a, b
    raise RuntimeError('Bezout search bound')


def matrix(vectors):
    t = Symbol('t')
    degree = 2*len(vectors)-1
    rows = []
    for v in vectors:
        w = tangent(v)
        expansions = [expand((v[0]+t*w[0])**k *
                             (v[1]+t*w[1])**(degree-k))
                      for k in range(degree+1)]
        rows.extend(([f.coeff(t, 0) for f in expansions],
                     [f.coeff(t, 1) for f in expansions]))
    return Matrix(rows)


def formula(vectors):
    answer = 1
    for i, v in enumerate(vectors):
        others = [z for j, z in enumerate(vectors) if i != j]
        w = tangent(v)
        d = prod(bracket(v, z) for z in others)
        c = 2*sum(bracket(w, z)*prod(bracket(v, y) for k, y in enumerate(others)
                                     if k != j)
                  for j, z in enumerate(others))
        answer = lcm(answer, abs(d)**3//gcd(abs(d), c))
    return answer


def main():
    directions = ((1, 0), (0, 1), (1, 1), (-1, 1), (1, 2), (2, 1))
    bank = [tuple(c) for n in range(1, 5) for c in combinations(directions, n)]
    bank += [((1, 0), (1, 8), (1, 24)),
             ((1, 0), (0, 1), (1, 1), (2, 1), (3, 1)),
             ((1, 0), (1, 3), (1, 9), (2, 1)),
             ((1, 0), (0, 1), (1, 1), (2, 1), (3, 1), (4, 1))]
    for vectors in bank:
        a = matrix(vectors)
        smith = smith_normal_form(a, domain=ZZ)
        diagonal = tuple(abs(int(smith[i, i])) for i in range(a.rows))
        check(diagonal[-1] == formula(vectors), 'largest Smith factor')
        check(abs(a.det()) == prod(abs(bracket(v, w))**4
                                  for v, w in combinations(vectors, 2)),
              'bracket determinant')
        inverse_denominator = lcm(*(int(x.q) for x in a.inv()))
        check(inverse_denominator == diagonal[-1], 'inverse denominator')
    print('PROJECTIVE_HERMITE_INDEPENDENT_SYMPY_AUDIT')
    print(f'configurations={len(bank)} node_counts=1..6')
    print('literal_symbolic_jets_vs_full_integer_Smith=PASS')
    print('bracket_determinant_and_exact_inverse_denominator=PASS')
    print('controls=infinity;all_projective_residues_p2_p3;deep_infinity;negative_slopes')
    print('PASS')


if __name__ == '__main__':
    main()
