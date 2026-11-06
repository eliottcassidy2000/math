#!/usr/bin/env python3
"""Independent rational RREF eigenray census, with no SymPy dependency."""
from fractions import Fraction
from itertools import product
from math import gcd, lcm
from functools import reduce
import json

def need(b, s):
    if not b:
        raise RuntimeError(s)

I = [[int(i == j) for j in range(4)] for i in range(4)]
M = {'O': [[1,0,1,0],[1,2,1,2],[0,0,0,0],[0,0,1,1]],
     'Z': [[2,1,1,1],[0,1,0,0],[0,0,1,1],[0,0,0,0]],
     'J': [[0,1,0,1],[0,0,0,0],[1,0,1,1],[0,0,1,0]]}

def mm(a, b):
    return [[sum(a[i][k] * b[k][j] for k in range(4)) for j in range(4)] for i in range(4)]

def kernel(a):
    r = [[Fraction(x) for x in row] for row in a]
    pivots = []
    for col in range(4):
        pivot = next((j for j in range(len(pivots), 4) if r[j][col]), None)
        if pivot is None:
            continue
        i = len(pivots)
        r[i], r[pivot] = r[pivot], r[i]
        scale = r[i][col]
        r[i] = [x / scale for x in r[i]]
        for j in range(4):
            if i != j:
                scale = r[j][col]
                r[j] = [x - scale * y for x, y in zip(r[j], r[i])]
        pivots.append(col)
    free = [c for c in range(4) if c not in pivots]
    result = []
    for col in free:
        v = [Fraction(0)] * 4
        v[col] = 1
        for i, pivot in enumerate(pivots):
            v[pivot] = -r[i][col]
        result.append(v)
    return result

def main():
    F = [2, 1, 3, 3]
    growth = {}
    for name, a in M.items():
        v = [2 * sum(F[i] * a[i][j] for i in range(4)) - 3 * F[j] for j in range(4)]
        need(min(v) >= 0, 'growth')
        growth[name] = v
    report = {}
    for q in (3, 5, 7, 9):
        bound = 0
        while 3 ** (bound + 1) <= q * 2 ** (bound + 1):
            bound += 1
        rays = {}
        checked = 0
        nonzero_kernels = 0
        for n in range(1, bound + 1):
            for word in product('OZJ', repeat=n):
                checked += 1
                a = I
                for c in word:
                    a = mm(M[c], a)
                basis = kernel([[a[i][j] - q * I[i][j] for j in range(4)] for i in range(4)])
                if not basis:
                    continue
                nonzero_kernels += 1
                need(len(basis) == 1, 'multidimensional kernel')
                v = basis[0]
                den = lcm(*(x.denominator for x in v))
                v = tuple(int(x * den) for x in v)
                if max(v) <= 0:
                    v = tuple(-x for x in v)
                if min(v) < 0:
                    continue
                common = reduce(gcd, v)
                v = tuple(x // common for x in v)
                matrix = tuple(tuple(row) for row in a)
                rays.setdefault(v, {}).setdefault(matrix, []).append(''.join(word))
        maximum = max((len(matrices) for matrices in rays.values()), default=0)
        need(maximum < q, 'possible digit grammar')
        report[q] = {'length_bound': bound, 'words_checked': checked,
                     'nonzero_kernels': nonzero_kernels,
                     'max_distinct_matrices_per_ray': maximum,
                     'rays': [{'ray': list(v), 'matrix_word_classes': list(matrices.values())}
                              for v, matrices in sorted(rays.items())]}
    print(json.dumps({'status': 'PASS', 'growth_surplus_coefficients': growth, 'bases': report}, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
