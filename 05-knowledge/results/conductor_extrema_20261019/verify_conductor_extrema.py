#!/usr/bin/env python3
"""Independent finite-exact controls for the uniform rank-three conductor extrema.
Run with python3 -B; also safe under -O: all checks use require, never assert.
Standard library only. All arithmetic is integral or Fraction; no random choices.
This script corroborates the all-r proofs, rather than replacing them.
"""
from collections import Counter
from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import comb, gcd
import json
from pathlib import Path


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def rank_fraction(matrix):
    """Ordinary rational Gaussian elimination."""
    a = [[Fraction(x) for x in row] for row in matrix]
    if not a:
        return 0
    pivot_row = 0
    for col in range(len(a[0])):
        p = next((i for i in range(pivot_row, len(a)) if a[i][col]), None)
        if p is None:
            continue
        a[pivot_row], a[p] = a[p], a[pivot_row]
        v = a[pivot_row][col]
        a[pivot_row] = [x/v for x in a[pivot_row]]
        for i in range(pivot_row + 1, len(a)):
            v = a[i][col]
            if v:
                a[i] = [x-v*y for x, y in zip(a[i], a[pivot_row])]
        pivot_row += 1
        if pivot_row == len(a):
            break
    return pivot_row


def rank_integer(matrix):
    """Fraction-free elimination, with primitive integer row reduction."""
    a = [list(row) for row in matrix]
    if not a:
        return 0
    pivot_row = 0
    for col in range(len(a[0])):
        p = next((i for i in range(pivot_row, len(a)) if a[i][col]), None)
        if p is None:
            continue
        a[pivot_row], a[p] = a[p], a[pivot_row]
        pivot = a[pivot_row][col]
        for i in range(pivot_row + 1, len(a)):
            v = a[i][col]
            if v:
                a[i] = [pivot*x-v*y for x, y in zip(a[i], a[pivot_row])]
                common = reduce(gcd, a[i], 0)
                if common:
                    a[i] = [x//common for x in a[i]]
        pivot_row += 1
        if pivot_row == len(a):
            break
    return pivot_row


def determinant3(a, b, c):
    return (a[0]*(b[1]*c[2]-b[2]*c[1])
            -a[1]*(b[0]*c[2]-b[2]*c[0])
            +a[2]*(b[0]*c[1]-b[1]*c[0]))


def uniform(points):
    return all(determinant3(*triple) != 0 for triple in combinations(points, 3))


def evaluate(points, q):
    monomials = [(q-a-b, a, b) for a in range(q+1) for b in range(q+1-a)]
    return [[x**c*y**a*z**b for c, a, b in monomials] for x, y, z in points]


def threshold(r):
    s = 0
    while comb(s+2, 2) < r:
        s += 1
    return s


def lower(r):
    s = threshold(r)
    return (s-1)*r - comb(s+2, 3) + 1


def upper(r):
    return (r-2)**2 // 4


def max_separator(r):
    return r // 2  # ceil((r-1)/2)


def rank_checked(matrix, label):
    a = rank_integer(matrix)
    b = rank_fraction(matrix)
    require(a == b, f'rank algorithms disagree: {label}')
    return a


def first_separator_degrees(points):
    """e_i is in the column image iff deleting row i lowers rank by one."""
    r = len(points)
    answer = [None]*r
    for q in range(max_separator(r)+1):
        matrix = evaluate(points, q)
        h = rank_integer(matrix)
        for i in range(r):
            if answer[i] is None:
                h_deleted = rank_integer(matrix[:i] + matrix[i+1:])
                if h_deleted == h-1:
                    answer[i] = q
        if all(x is not None for x in answer):
            return answer
    raise RuntimeError('separator upper bound failed')


def check_formulas():
    for r in range(3, 1001):
        s, c = threshold(r), max_separator(r)
        require(s <= c, 'thresholds reversed')
        require(lower(r) == sum(max(0, r-comb(q+2, 2)) for q in range(1, r)),
                'lower formula mismatch')
        require(upper(r) == sum(max(0, r-(2*q+1)) for q in range(1, r)),
                'upper formula mismatch')
        require(lower(r) <= upper(r), 'bounds reversed')
    return {'r_values': [3, 1000], 'cases': 998}


def check_explicit_families():
    results = []
    matrices = 0
    deleted_rows = 0
    for r in range(3, 19):
        s, c = threshold(r), max_separator(r)
        for name, points, expected_rank, d, defect in (
            ('conic', [(1, t, t*t) for t in range(1, r+1)],
             lambda q: min(r, 2*q+1), c, upper(r)),
            ('power_curve', [(1, t, t**r) for t in range(1, r+1)],
             lambda q: min(r, comb(q+2, 2)), s, lower(r))):
            require(uniform(points), f'nonuniform {name} r={r}')
            h = []
            for q in range(c+2):
                observed = rank_checked(evaluate(points, q), (name, r, q))
                require(observed == expected_rank(q), f'Hilbert rank {name} r={r} q={q}')
                h.append(observed)
                matrices += 1
            require(sum(r-x for x in h[1:]) == defect, f'defect {name} r={r}')
            before = evaluate(points, d-1)
            full_rank_before = rank_checked(before, (name, r, 'before'))
            for i in range(r):
                deleted = before[:i] + before[i+1:]
                require(rank_checked(deleted, (name, r, 'deleted', i)) == full_rank_before,
                        f'premature separator {name} r={r} i={i}')
                deleted_rows += 1
            require(rank_checked(evaluate(points, d), (name, r, 'separator')) == r,
                    f'missing endpoint separators {name} r={r}')
            if name == 'power_curve':
                for q in range(r):
                    exponents = [a+r*b for a in range(q+1) for b in range(q+1-a)]
                    require(len(set(exponents)) == comb(q+2, 2), 'exponent collision below r')
                collision = [a+r*b for a in range(r+1) for b in range(r+1-a)]
                require(len(set(collision)) < len(collision), 'missing q=r collision boundary')
            results.append({'r': r, 'family': name, 'hilbert_ranks': h,
                            'seminormal_defect': defect, 'normalization_defect': defect+r-1,
                            'all_branch_conductors': d})
    return {'r_values': [3, 18], 'families': 32, 'double_checked_full_matrices': matrices,
            'double_checked_row_deletions': deleted_rows, 'results': results}


def check_minimum_conductor_caveat():
    r = 7
    points = [(1, t, t*t) for t in range(r-1)] + [(1, r-1, (r-1)**2+1)]
    require(uniform(points), 'D7 not uniform')
    h = [rank_checked(evaluate(points, q), ('D7', q)) for q in range(5)]
    d = first_separator_degrees(points)
    require(h == [1, 3, 6, 7, 7], 'D7 Hilbert ranks')
    require(sum(r-x for x in h[1:]) == lower(r) == 5, 'D7 not minimum defect')
    require(d == [3]*6+[2], 'D7 branchwise conductor caveat')
    return {'points': points, 'hilbert_ranks': h, 'seminormal_defect': 5, 'conductors': d}


def check_finite_universe():
    # Ten conic points plus two rational points not on that conic.
    # Enumerate every subset of size >=3; retain exactly the U(3,r) realizations.
    pool = [(1, t, t*t) for t in range(1, 11)] + [(1, 11, 122), (1, 12, 146)]
    counts = Counter()
    max_counts = Counter()
    min_counts = Counter()
    quadratic_classification_cases = 0
    for r in range(3, len(pool)+1):
        for indices in combinations(range(len(pool)), r):
            points = [pool[i] for i in indices]
            if not uniform(points):
                continue
            counts[r] += 1
            c, s = max_separator(r), threshold(r)
            ranks = [rank_integer(evaluate(points, q)) for q in range(c+1)]
            delta = sum(r-h for h in ranks[1:])
            require(lower(r) <= delta <= upper(r), f'universal bounds {indices}')
            require(all(min(r, 2*q+1) <= h <= min(r, comb(q+2, 2))
                        for q, h in enumerate(ranks)), f'rank sandwich {indices}')
            exact_min_condition = ranks[s-1] == comb(s+1, 2) and ranks[s] == r
            require((delta == lower(r)) == exact_min_condition, f'minimum equality {indices}')
            if delta == lower(r):
                min_counts[r] += 1
            if delta == upper(r):
                max_counts[r] += 1
            if r >= 6:
                quadratic_classification_cases += 1
                require((delta == upper(r)) == (ranks[2] == 5),
                        f'conic maximum equivalence {indices}')
            else:
                require(delta == r-3, f'small-r forced defect {indices}')
    return {'pool': pool, 'uniform_subsets_by_size': dict(sorted(counts.items())),
            'uniform_subsets': sum(counts.values()),
            'minimum_attainers_by_size': dict(sorted(min_counts.items())),
            'maximum_attainers_by_size': dict(sorted(max_counts.items())),
            'quadratic_maximum_classification_cases': quadratic_classification_cases}


def main():
    result = {'status': 'PASS', 'arithmetic': 'integer and exact Fraction only',
              'formula_controls': check_formulas(),
              'explicit_families': check_explicit_families(),
              'minimum_does_not_force_equal_branch_conductors': check_minimum_conductor_caveat(),
              'finite_universe': check_finite_universe()}
    rendered = json.dumps(result, indent=2, sort_keys=True)
    print(rendered)


if __name__ == '__main__':
    main()
