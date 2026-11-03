#!/usr/bin/env python3
"""Independent exact-Q S169 collision-bank audit; imports no producer.

Constructs wedge constraints by alternating tensors, rather than copying
the producer's six-term coordinate formula. All checks survive -O.
"""
from fractions import Fraction as Q
from itertools import combinations, product
import json

C = ((0, 1, 0), (1, 1, 1), (1, 0, 0))
D = ((1, 0, 1), (0, 2, 0), (0, 0, 0))
V = ((1, 2, 3), (2, 3, 6), (0, 1, 2))
BASIS = ('', 'C', 'D', 'CC', 'CD', 'DC', 'DD', 'CCD', 'CDD')
BANK = ('CDDD', 'DD', 'CDCD', 'CCD', 'CDC', 'DDCC',
        'DCD', 'DCDC', 'CDD', 'CCC', 'CCCC', 'C')

def check(condition, message):
    if not condition:
        raise RuntimeError(message)

def mm(a, b):
    return [[sum(x*y for x, y in zip(row, col)) for col in zip(*b)]
            for row in a]

def flat(a):
    return [x for row in a for x in row]

def word(w):
    a = [[int(i == j) for j in range(3)] for i in range(3)]
    for c in w:
        a = mm(C if c == 'C' else D, a)
    return a

def rref(a):
    a = [[Q(x) for x in row] for row in a if any(row)]
    if not a:
        return [], []
    r = 0
    pivots = []
    for j in range(len(a[0])):
        k = next((k for k in range(r, len(a)) if a[k][j]), None)
        if k is None:
            continue
        a[k], a[r] = a[r], a[k]
        p = a[r][j]
        a[r] = [x/p for x in a[r]]
        for k in range(len(a)):
            if k != r and a[k][j]:
                p = a[k][j]
                a[k] = [x-p*y for x, y in zip(a[k], a[r])]
        pivots.append(j)
        r += 1
        if r == len(a):
            break
    return a, pivots

def rank(a):
    return len(rref(a)[1])

def inverse(a):
    n = len(a)
    rr, pp = rref([row + [int(i == j) for j in range(n)]
                   for i, row in enumerate(a)])
    check(pp[:n] == list(range(n)), 'invertible coordinate basis')
    return [row[n:] for row in rr[:n]]

def normalized_constraints(t):
    """Columns are c_ij, representing A(e_i) = e_i wedge c_ij e_j."""
    n = len(t)
    pairs = [(i, j) for i in range(n) for j in range(n) if i != j]
    triples = list(combinations(range(n), 3))
    ix = {abc: r for r, abc in enumerate(triples)}
    out = [[Q(0) for _ in pairs] for _ in triples]
    for col, (i, j) in enumerate(pairs):
        for k in range(n):
            if k == i or k == j:
                continue
            abc = (k, i, j)
            inversions = sum(abc[a] > abc[b] for a in range(3)
                             for b in range(a+1, 3))
            out[ix[tuple(sorted(abc))]][col] += (-1)**inversions*t[k]*t[i]
    return out

def coords(basis, vectors):
    inv = inverse(list(map(list, zip(*basis))))
    return [flat(mm(inv, [[x] for x in v])) for v in vectors]

def main():
    vectors = [flat(mm(word(w), V)) for w in BANK]
    bcoords = coords(vectors[:9], vectors)
    positive = [row for t in bcoords[9:] for row in normalized_constraints(t)]
    positive_rank = rank(positive)
    check(positive_rank == 63, 'exact-Q positive bank constraint rank')

    relations = [[-x for x in bcoords[9+k]] +
                 [Q(j == k) for j in range(3)] for k in range(3)]
    period = []
    for ell in relations:
        check(all(sum(ell[r]*vectors[r][i] for r in range(12)) == 0
                  for i in range(9)), 'exact relation')
        for i, j in combinations(range(9), 2):
            row = [Q(0)]*108
            for r, t in enumerate(vectors):
                row[9*r+j] = ell[r]*t[i]
                row[9*r+i] = -ell[r]*t[j]
            period.append(row)
    period_rank = rank(period)
    check(period_rank == 87, 'exact-Q original period rank')

    radius_words = [''.join(w) for k in range(4) for w in product('CD', repeat=k)]
    radius_coords = coords([flat(word(w)) for w in BASIS],
                           [flat(word(w)) for w in radius_words])
    radius_rows = [row for t in radius_coords for row in normalized_constraints(t)]
    radius_rank = rank(radius_rows)
    check(radius_rank == 60, 'exact-Q radius-three constraint rank')

    pairs = [(i, j) for i in range(9) for j in range(9) if i != j]
    index = {v: k for k, v in enumerate(pairs)}
    h = [[0]*72 for _ in range(3)]
    h[0][index[(4, 5)]] = h[0][index[(5, 4)]] = 1
    for j, z in enumerate((-1, -1, -1, 2, 1, 1, 1, 0, 0)):
        if j != 7:
            h[1][index[(7, j)]] = z
    for j, z in enumerate((1, 1, 1, 0, -2, -2, -1, 0, 0)):
        if j != 8:
            h[2][index[(8, j)]] = z
    check(all(sum(x*y for x, y in zip(row, hostile)) == 0
              for row in radius_rows for hostile in h), 'exact hostile identities')
    common = [[int(j == k) for i, j in pairs] for k in range(9)]
    kernel_rank = rank(common+h)
    check(kernel_rank == 12, 'exact-Q independent common and hostile classes')
    radius_two_rank = rank([flat(word(w)) for w in radius_words if len(w) <= 2])
    check(radius_two_rank == 7, 'exact-Q radius-two span')

    generic = []
    for n in range(3, 10):
        extra = [[1]*n, list(range(1, n+1)), [i*i for i in range(1, n+1)]]
        rk = rank([r for t in extra for r in normalized_constraints(t)])
        check(rk == n*(n-2), 'exact-Q generic n+3 bank')
        generic.append({'dimension': n, 'rank_Q': rk})
    print(json.dumps({'status': 'PASS', 'positive_constraint_rank_Q': positive_rank,
                      'direct_period_rank_Q': period_rank,
                      'radius_three_constraint_rank_Q': radius_rank,
                      'radius_three_kernel_rank_Q': kernel_rank,
                      'radius_three_blind_quotient_dimension_Q': 3,
                      'radius_two_span_Q': radius_two_rank,
                      'generic_banks': generic}, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
