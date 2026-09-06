#!/usr/bin/env python3
"""Exact controls for generic dual testers and the actual integral Student map.

Run from the repository root. No imported research producer, floating point,
randomness, symbolic assumptions, or assertions (including under python -O).
"""
from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd, lcm, prod
import json

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ

GATES = 0


def check(predicate, message):
    global GATES
    GATES += 1
    if not predicate:
        raise RuntimeError(message)


def path_gcd(weights):
    """Size-indexed gcds of products of nonadjacent edges of one path."""
    cap = (len(weights) + 1) // 2
    previous_two = [1] + [0] * cap
    previous_one = previous_two[:]
    for weight in weights:
        current = [1] + [gcd(previous_one[k], weight * previous_two[k - 1])
                         for k in range(1, cap + 1)]
        previous_two, previous_one = previous_one, current
    return previous_one


def path_weights(m):
    even = []
    for k in range((m + 1) // 2):
        if k:
            even.append(12 * k)
        even.append(2 * m - 2 * k)
    odd = []
    for k in range(m // 2):
        odd.extend([6 * (2 * k + 1), 2 * m - 2 * k - 1])
    return even, odd


def student_smith(m):
    left, right = [path_gcd(w) for w in path_weights(m)]
    divisors = []
    for k in range(m + 1):
        divisors.append(reduce(gcd, (left[a] * right[k - a]
                                    for a in range(len(left))
                                    if 0 <= k - a < len(right)), 0))
    return divisors, [divisors[k] // divisors[k - 1]
                      for k in range(1, m + 1)]


def actual_matrix(m):
    """Build directly by differentiating monomials, not by the path formula."""
    x = sp.Symbol('x')
    return sp.Matrix(m + 1, m, lambda i, j: sp.Poly(
        (x*x + 6) * sp.diff(x**j, x) - 2*m*x*x**j, x).nth(i))


def student_covector(m):
    row = [Fraction(0) for _ in range(m + 1)]
    row[0] = Fraction(1)
    for k in range(1, m // 2 + 1):
        row[2*k] = row[2*k-2] * Fraction(6*(2*k-1), 2*m-2*k+1)
    denominator = lcm(*(q.denominator for q in row))
    integer = [int(q * denominator) for q in row]
    content = reduce(gcd, integer)
    return [a // content for a in integer]


def solve_response(m, rhs):
    """Unique rational inverse on its image, retaining the final raw equation."""
    theta = [Fraction(0) for _ in range(m + 2)]
    for k in range(m, 0, -1):
        theta[k-1] = (Fraction(rhs[k]) - 6*(k+1)*theta[k+1]) / (k-1-2*m)
    if Fraction(rhs[0]) != 6*theta[1]:
        return None
    return theta[:m]


def generic_dual_checks():
    # One fully symbolic source and two right-hand sides; independent bordered
    # determinants check the coefficient statement, not just numerical equality.
    n, r, bank = 3, 2, 2
    A = sp.Matrix(n, r, sp.symbols('a0:6'))
    B = sp.Matrix(n, bank, sp.symbols('b0:6'))
    variables = sp.symbols('u0:6')
    U = sp.Matrix(r, n, variables)
    G = U*A
    S = G.det()*B - A*G.adjugate()*U*B
    expansion = sp.zeros(n, bank)
    bordered_seen = set()
    for I in combinations(range(n), r):
        AI, BI = A[list(I), :], B[list(I), :]
        residue = AI.det()*B - A*AI.adjugate()*BI
        expansion += U[:, list(I)].det()*residue
        for j in range(n):
            for alpha in range(bank):
                if j in I:
                    check(sp.expand(residue[j, alpha]) == 0, 'pivot residue')
                else:
                    block = A[list(I) + [j], :].row_join(B[list(I) + [j], alpha])
                    check(sp.expand(residue[j, alpha] - block.det()) == 0,
                          'bordered determinant')
                    bordered_seen.add(sp.expand(block.det()))
    for value in S - expansion:
        check(sp.expand(value) == 0, 'generic Cauchy-Binet residue identity')
    coefficient_seen = set()
    for value in S:
        polynomial = sp.Poly(sp.expand(value), *variables)
        for monomial, coefficient in polynomial.terms():
            if coefficient != 0:
                coefficient_seen.add(sp.expand(coefficient))
                check(sum(monomial) == r, 'dual degree')
                check(any(sp.expand(coefficient - sign*minor) == 0
                          for sign in [-1, 1] for minor in bordered_seen),
                      'coefficient is signed mixed minor')
    for minor in bordered_seen:
        check(any(sp.expand(minor - sign*coefficient) == 0
                  for sign in [-1, 1] for coefficient in coefficient_seen),
              'all mixed minors retained')
    # Arbitrary rings follow from the integer polynomial identity. These tests
    # exercise additional sizes/ranks, including rank drops, with exact integers.
    for n, r, bank in [(4, 2, 3), (4, 3, 2), (5, 3, 2)]:
        for seed in range(7):
            A = sp.Matrix(n, r, lambda i, j: ((i+1)*(j+2)+seed*(i-j)) % 7-3)
            B = sp.Matrix(n, bank, lambda i, j: (i*i+2*j+seed) % 9-4)
            U = sp.Matrix(r, n, lambda i, j: (3*i+j*j+seed) % 5-2)
            G = U*A
            residue = G.det()*B - A*G.adjugate()*U*B
            expected = sp.zeros(n, bank)
            for I in combinations(range(n), r):
                AI, BI = A[list(I), :], B[list(I), :]
                expected += U[:, list(I)].det()*(AI.det()*B-A*AI.adjugate()*BI)
            check(residue == expected, 'integer specialization')
    # A full rank complex source can have zero bilinear Gram matrix and zero
    # normal-equation measurement on an incompatible target.
    A = sp.Matrix([1, sp.I, 0])
    B = sp.Matrix([0, 0, 1])
    check(A.rank() == 1 and (A.T*A)[0] == 0 and (A.T*B)[0] == 0,
          'isotropic Gram hostile')
    check(A.row_join(B).rank() == 2, 'Gram hostile actual incompatibility')
    A, B = sp.Matrix([2, 0]), sp.Matrix([1, 0])
    check(A.row_join(B).det() == 0 and Fraction(1, 2).denominator == 2,
          'raw mixed minors do not imply integral solvability')
    print('DUAL_TESTER symbolic mixed-minor identity and complex/integral hostiles PASS')


def main():
    generic_dual_checks()
    inherited = {5: [21, 0, 14, 0, 36, 0],
                 6: [77, 0, 42, 0, 84, 0, 360],
                 7: [143, 0, 66, 0, 108, 0, 360, 0],
                 8: [715, 0, 286, 0, 396, 0, 1080, 0, 5040]}
    for m, covector in inherited.items():
        check(student_covector(m) == covector, 'inherited THM4308 Student row')
    for m in range(2, 49):
        A = actual_matrix(m)
        divisors, smith = student_smith(m)
        D = smith_normal_form(A, domain=ZZ)
        independent = [abs(int(D[i, i])) for i in range(m)]
        check(smith == independent, 'SymPy Smith mismatch m=' + str(m))
        covector = student_covector(m)
        check(sp.Matrix([covector])*A == sp.zeros(1, m), 'Student annihilator')
        check(reduce(gcd, covector) == 1, 'primitive Student annihilator')
        for j in range(m):
            rhs = list(A[:, j])
            solution = solve_response(m, rhs)
            check(solution == [int(i == j) for i in range(m)], 'raw inverse column')
        if m <= 6:
            # All square minors, without using sparsity or a matching formula.
            for k in range(m + 1):
                direct = 0
                for rows in combinations(range(m + 1), k):
                    for cols in combinations(range(m), k):
                        direct = gcd(direct, int(A[list(rows), list(cols)].det()))
                check(abs(direct) == divisors[k], 'exhaustive minors')
        if m in [2, 5, 8, 14, 24, 48]:
            print('SMITH', json.dumps({'m': m, 'factors': smith,
                                       'torsion_order': divisors[-1],
                                       'uniform_denominator': smith[-1]},
                                      separators=(',', ':')))
    for m in range(2, 257):
        divisors, smith = student_smith(m)
        check(len(smith) == m and all(v > 0 for v in divisors), 'full rank DP')
        check(all(smith[i+1] % smith[i] == 0 for i in range(m-1)), 'Smith divisibility')
        check(prod(smith) == divisors[-1], 'Smith product')
        # An explicit maximal minor of each parity block bounds the torsion
        # order, and therefore proves the stated torsion-prime upper bound.
        r, s = (m+1)//2, m//2
        selected_minor = prod(2*m-2*k for k in range(r))*prod(6*(2*k+1) for k in range(s))
        check(selected_minor % divisors[-1] == 0, 'maximal-minor upper bound')
    A = actual_matrix(2)
    bank = [[0, 1, 0], [2, 0, -1]]
    solutions = [solve_response(2, b) for b in bank]
    check(solutions == [[Fraction(-1, 4), Fraction(0)],
                        [Fraction(0), Fraction(1, 3)]], 'm2 actual bank hostile')
    denominator = lcm(*(v.denominator for column in solutions for v in column))
    check(denominator == 12, 'joint denominator is sharp')
    check(solve_response(2, [1, 0, 0]) is None, 'raw Student incompatibility')
    print('ACTUAL_BANK m=2 rhs=(x,2-x^2) solutions=(-1/4,x/3) common_denominator=12')
    _, row14 = student_smith(14)
    check(row14[-1] == 720720, 'row14 arithmetic precision')
    row14_denominators = []
    for j in range(1, 15, 2):
        solution = solve_response(14, [int(k == j) for k in range(15)])
        check(solution is not None, 'row14 odd monomial compatibility')
        row14_denominators.append(lcm(*(v.denominator for v in solution)))
    check(row14_denominators == [28, 182, 2184, 4004, 20020, 18018, 16016],
          'row14 individual precision')
    check(lcm(*row14_denominators) == row14[-1], 'row14 bank attains full exponent')
    print('ROW14_BANK powers=(1,3,5,7,9,11,13) denominators=' +
          json.dumps(row14_denominators, separators=(',', ':')) +
          ' common_denominator=720720')
    localized_exponent = row14[-1]
    for prime in [2, 3, 7]:
        while localized_exponent % prime == 0:
            localized_exponent //= prime
    check(localized_exponent == 715, 'row14 localization hostile')
    print('LOCALIZED_ROW14 invert=(2,3,7) residual_exponent=715=5*11*13')
    print('UNIVERSE exact Smith m=2..48; every minor m=2..6; DP structural checks m=2..256')
    print('PASS gates=' + str(GATES))


if __name__ == '__main__':
    main()
