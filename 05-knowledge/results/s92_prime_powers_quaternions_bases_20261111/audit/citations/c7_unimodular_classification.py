#!/usr/bin/env python3
"""Auditor 'citations' - finite check of the classification used in THM-5134 (F) and THM-5135 (2):
over Z_p, p odd, a unimodular quadratic form is determined by its rank and its discriminant modulo squares.

What the theorems use is the consequence for class counts.  Checked here by brute force:
for EVERY symmetric matrix A modulo q = p^k with det A a unit, the vector of class counts
    #{x in (Z/q)^n : x^T A x = c},  c in Z/q,
depends only on (n, Legendre symbol of det A).   Cases: q = 3, 5, 7 (n = 2, 3), q = 3 (n = 4), q = 9, 25 (n = 2),
q = 9 (n = 3).
Also at q = 2, 4, 8 the analogous statement FAILS for integral forms x^T A x (A symmetric, odd determinant):
several different count vectors occur for the same rank and determinant class - the prime 2 needs its own
classification (type, i.e. even / odd, and more).
"""
import itertools
import sys
import numpy as np
from numba import njit


@njit(cache=False)
def class_counts(A, n, q):
    cnt = np.zeros(q, dtype=np.int64)
    x = np.zeros(n, dtype=np.int64)
    total = q ** n
    for idx in range(total):
        t = idx
        for i in range(n):
            x[i] = t % q
            t //= q
        s = 0
        for i in range(n):
            for j in range(n):
                s += x[i] * A[i, j] * x[j]
        cnt[s % q] += 1
    return cnt


def det_mod(A, q):
    n = A.shape[0]
    # exact integer determinant by permutation expansion (n <= 4)
    tot = 0
    for perm in itertools.permutations(range(n)):
        sgn = 1
        for i in range(n):
            for j in range(i + 1, n):
                if perm[i] > perm[j]:
                    sgn = -sgn
        prod = 1
        for i in range(n):
            prod *= int(A[i, perm[i]])
        tot += sgn * prod
    return tot % q


def legendre(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def run(q, p, n):
    idxs = [(i, j) for i in range(n) for j in range(i, n)]
    seen = {}
    nmat = 0
    for vals in itertools.product(range(q), repeat=len(idxs)):
        A = np.zeros((n, n), dtype=np.int64)
        for (i, j), v in zip(idxs, vals):
            A[i, j] = v
            A[j, i] = v
        d = det_mod(A, q)
        if d % p == 0:
            continue
        nmat += 1
        key = legendre(d, p) if p != 2 else d % 8 if q >= 8 else d % q
        vec = tuple(int(v) for v in class_counts(A, n, q))
        seen.setdefault(key, set()).add(vec)
    return nmat, seen


def main():
    ok = True
    for q, p, n in ((3, 3, 2), (3, 3, 3), (3, 3, 4), (5, 5, 2), (5, 5, 3), (7, 7, 2), (7, 7, 3), (9, 3, 2), (25, 5, 2),
                    (9, 3, 3)):
        nmat, seen = run(q, p, n)
        distinct = {k: len(v) for k, v in seen.items()}
        good = all(v == 1 for v in distinct.values()) and len(seen) == 2
        # the two classes differ
        vecs = [next(iter(v)) for v in seen.values()]
        differ = len(vecs) == 2 and vecs[0] != vecs[1]
        print("q = %2d, rank %d: %7d unimodular symmetric matrices; count vectors per discriminant class: %s; "
              "two classes with different counts: %s" % (q, n, nmat, distinct, differ))
        ok &= good and differ
        sys.stdout.flush()
    print("odd p: class counts depend only on (rank, discriminant class): %s" % ok)
    for q, n in ((2, 2), (4, 2), (8, 2), (2, 3), (4, 3), (2, 4)):
        nmat, seen = run(q, 2, n)
        distinct = {k: len(v) for k, v in seen.items()}
        print("q = %d, rank %d (prime 2, forms x^T A x): %6d matrices of odd determinant; distinct count vectors per "
              "determinant residue: %s" % (q, n, nmat, distinct))
    print("ALL CHECKS PASSED" if ok else "SOME CHECK FAILED")


if __name__ == "__main__":
    main()
