#!/usr/bin/env python3
"""Verifier (lane digits): own exact linear algebra (integer, fraction-free).  Independent of the lane's affine.py."""
from math import gcd


def rowmat(r, M):
    d = len(M[0])
    return [sum(r[k] * M[k][j] for k in range(len(r)) if r[k]) for j in range(d)]


def matmul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum(A[i][k] * B[k][j] for k in range(m)) for j in range(p)] for i in range(n)]


def matvec(M, v):
    return [sum(M[i][j] * v[j] for j in range(len(v))) for i in range(len(M))]


def dot(r, s):
    return sum(x * y for x, y in zip(r, s))


def primitive(v):
    g = 0
    for x in v:
        g = gcd(g, x)
    if g > 1:
        v = [x // g for x in v]
    return v


class Space:
    """row space over Q, integer echelon basis (pivot = first nonzero entry, pivots distinct)."""

    def __init__(self, d):
        self.d = d
        self.rows = {}        # pivot -> row

    def reduce(self, v):
        v = list(v)
        for p in sorted(self.rows):
            if v[p] != 0:
                row = self.rows[p]
                q = row[p]
                c = v[p]
                v = [q * x - c * y for x, y in zip(v, row)]
                v = primitive(v)
        return v

    def add(self, v):
        v = self.reduce(v)
        p = next((i for i, x in enumerate(v) if x != 0), None)
        if p is None:
            return False
        self.rows[p] = primitive(v)
        return True

    def contains(self, v):
        return not any(self.reduce(v))

    def dim(self):
        return len(self.rows)

    def kills(self, s):
        """True iff every row vanishes on the column vector s."""
        return all(dot(r, s) == 0 for r in self.rows.values())


def krylov(seeds, mats, d):
    V = Space(d)
    queue = []
    for r in seeds:
        if V.add(r):
            queue.append(list(r))
    while queue and V.dim() < d:
        r = queue.pop()
        for G in mats:
            r2 = rowmat(r, G)
            if V.add(r2):
                queue.append(r2)
    return V


def class_space(L, mats, mults):
    d = len(L)
    lam = []
    for G, a in zip(mats, mults):
        LG = rowmat(L, G)
        lam.append([LG[j] - a * L[j] for j in range(d)])
    seeds = []
    for l in lam:
        for G in mats:
            lg = rowmat(l, G)
            seeds.append([lg[j] - l[j] for j in range(d)])
    return krylov(seeds, mats, d), lam


def charpoly_int(M):
    """det(xI - M) for an integer matrix by exact interpolation-free expansion (Leibniz for d <= 4, else Berkowitz)."""
    d = len(M)
    # Berkowitz (division-free)
    C = [1]
    for k in range(d):
        # leading principal submatrix A_k = M[:k+1, :k+1]; a = M[k][k]; R = M[k][:k]; Ccol = M[:k][k]; A = M[:k][:k]
        a = M[k][k]
        R = [M[k][j] for j in range(k)]
        Cc = [M[i][k] for i in range(k)]
        A = [[M[i][j] for j in range(k)] for i in range(k)]
        # Toeplitz column: [1, -a, -R Cc, -R A Cc, -R A^2 Cc, ...]
        col = [1, -a]
        v = Cc
        for t in range(k):
            col.append(-sum(R[i] * v[i] for i in range(k)))
            v = [sum(A[i][j] * v[j] for j in range(k)) for i in range(k)]
        newC = [0] * (len(C) + 1)
        for i in range(len(newC)):
            s = 0
            for j in range(len(C)):
                if 0 <= i - j < len(col):
                    s += col[i - j] * C[j]
            newC[i] = s
        C = newC
    # C[0] x^d + C[1] x^(d-1) + ... + C[d]
    return C


def polyval_desc(c, x):
    v = 0
    for t in c:
        v = v * x + t
    return v


def int_roots_exact(c):
    """ALL integer roots of the integer polynomial with descending coefficients c (exact: rational root theorem)."""
    c = list(c)
    roots = set()
    while len(c) > 1 and c[-1] == 0:
        roots.add(0)
        c.pop()
    if len(c) <= 1:
        return roots
    n = abs(c[-1])
    # divisors of n by trial division
    divs = [1]
    m = n
    p = 2
    while p * p <= m:
        if m % p == 0:
            e = 0
            while m % p == 0:
                m //= p
                e += 1
            divs = [x * p ** k for x in divs for k in range(e + 1)]
        p += 1 if p == 2 else 2
    if m > 1:
        divs = [x * q for x in divs for q in (1, m)]
    for x in divs:
        for t in (x, -x):
            if polyval_desc(c, t) == 0:
                roots.add(t)
    return roots


def int_eigs(M):
    return sorted(int_roots_exact(charpoly_int(M)))


if __name__ == "__main__":
    import random
    rng = random.Random(1)
    # self-test of charpoly against sympy on random matrices
    import sympy
    ok = True
    for _ in range(40):
        d = rng.choice([3, 4, 5, 6])
        M = [[rng.randint(-5, 5) for _ in range(d)] for _ in range(d)]
        x = sympy.symbols("x")
        cp = sympy.Matrix(M).charpoly(x).all_coeffs()
        if [int(t) for t in cp] != charpoly_int(M):
            ok = False
    print("charpoly_int = sympy charpoly on 40 random integer matrices:", ok)
    print("int roots of (x-3)(x+2)(x^2+1):", sorted(int_roots_exact([1, -1, -5, -1, -6])))
