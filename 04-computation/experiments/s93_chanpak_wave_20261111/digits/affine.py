#!/usr/bin/env python3
"""Lane "digits": exact linear algebra for affine digit systems inside a linear representation.

Setting: states s in Q^d, operations G_1..G_k (integer d x d matrices), readout L (row vector), H = L s.
A class K of states closed under the G_i on which  H(G_i s) = a_i H(s) + b_i  (a_i, b_i constants of the class).

Facts used (proved in PROOFS.txt, P3):
  (F1) if H is unbounded on K then a_i is an eigenvalue of G_i;
  (F2) with lambda_i = L G_i - a_i L, the set
           C = { s : lambda_i (G_w - I) s = 0 for every word w and every i }
       is a G-invariant linear subspace, equal to the annihilator of the smallest row space V containing the rows
       lambda_i (G_j - I) and closed under right multiplication by every G_j;  a state s generates an exact affine
       class (its orbit) iff s lies in C, and then b_i = lambda_i s.
So for given operations and multipliers the classification is one Krylov closure.
"""
from fractions import Fraction

import numpy as np


def matmul(A, B):
    n, m, p = len(A), len(B), len(B[0])
    return [[sum(A[i][k] * B[k][j] for k in range(m)) for j in range(p)] for i in range(n)]


def rowmat(r, M):
    d = len(M[0])
    return [sum(r[k] * M[k][j] for k in range(len(r))) for j in range(d)]


def matvec(M, v):
    return [sum(M[i][j] * v[j] for j in range(len(v))) for i in range(len(M))]


def det_int(M):
    """Bareiss, exact for integer matrices."""
    n = len(M)
    A = [list(r) for r in M]
    sign = 1
    prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            p = next((r for r in range(k + 1, n) if A[r][k] != 0), None)
            if p is None:
                return 0
            A[k], A[p] = A[p], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


def charpoly(M):
    """coefficients c[0..d] of det(xI - M) = sum c[k] x^k (Faddeev-LeVerrier, exact)."""
    d = len(M)
    c = [0] * (d + 1)
    c[d] = 1
    Mk = [[0] * d for _ in range(d)]
    for k in range(1, d + 1):
        # Mk = M * Mk + c[d-k+1] I
        for i in range(d):
            Mk[i][i] += c[d - k + 1]
        Mk = matmul(M, Mk)
        tr = sum(Mk[i][i] for i in range(d))
        assert tr % k == 0
        c[d - k] = -tr // k
    return c


def polyval(c, x):
    v = 0
    for k in range(len(c) - 1, -1, -1):
        v = v * x + c[k]
    return v


def int_eigs(M, lo=-8, hi=64):
    """integer eigenvalues: every integer in [lo, hi] tested exactly, others through float roots + exact check."""
    c = charpoly(M)
    out = set()
    for a in range(lo, hi + 1):
        if polyval(c, a) == 0:
            out.add(a)
    try:
        roots = np.roots([float(x) for x in reversed(c)])
        for z in roots:
            if abs(z.imag) < 1e-6 * max(1.0, abs(z)):
                a = int(round(z.real))
                for t in (a - 1, a, a + 1):
                    if polyval(c, t) == 0:
                        out.add(t)
    except Exception:
        pass
    return sorted(out)


def charpoly_q(M):
    """coefficients of det(xI - M) for a rational matrix (Faddeev-LeVerrier over Q)."""
    d = len(M)
    Mq = [[Fraction(x) for x in row] for row in M]
    c = [Fraction(0)] * (d + 1)
    c[d] = Fraction(1)
    Mk = [[Fraction(0)] * d for _ in range(d)]
    for k in range(1, d + 1):
        for i in range(d):
            Mk[i][i] += c[d - k + 1]
        Mk = matmul(Mq, Mk)
        tr = sum(Mk[i][i] for i in range(d))
        c[d - k] = -tr / k
    return c


def int_eigs_q(M, lo=-8, hi=64):
    """integer eigenvalues of a rational matrix: every integer in [lo, hi] exactly, others via float roots + exact check."""
    c = charpoly_q(M)
    out = set()
    for a in range(lo, hi + 1):
        if polyval(c, a) == 0:
            out.add(a)
    try:
        roots = np.roots([float(x) for x in reversed(c)])
        for z in roots:
            if abs(z.imag) < 1e-6 * max(1.0, abs(z)):
                a = int(round(z.real))
                for t in (a - 1, a, a + 1):
                    if polyval(c, t) == 0:
                        out.add(t)
    except Exception:
        pass
    return sorted(out)


class RowSpace:
    """row space over Q in reduced echelon form."""

    def __init__(self, d):
        self.d = d
        self.rows = []      # each: list of Fractions, pivot normalised to 1
        self.piv = []

    def reduce(self, r):
        r = [Fraction(x) for x in r]
        for row, p in zip(self.rows, self.piv):
            if r[p] != 0:
                f = r[p]
                r = [x - f * y for x, y in zip(r, row)]
        return r

    def add(self, r):
        r = self.reduce(r)
        p = next((i for i, x in enumerate(r) if x != 0), None)
        if p is None:
            return False
        f = r[p]
        r = [x / f for x in r]
        for i in range(len(self.rows)):
            if self.rows[i][p] != 0:
                g = self.rows[i][p]
                self.rows[i] = [x - g * y for x, y in zip(self.rows[i], r)]
        self.rows.append(r)
        self.piv.append(p)
        return True

    def dim(self):
        return len(self.rows)

    def contains(self, r):
        return all(x == 0 for x in self.reduce(r))

    def annihilator(self):
        """basis (integer primitive vectors) of { s : rho s = 0 for every row rho }."""
        d = self.d
        free = [j for j in range(d) if j not in self.piv]
        out = []
        for f in free:
            v = [Fraction(0)] * d
            v[f] = Fraction(1)
            for row, p in zip(self.rows, self.piv):
                v[p] = -row[f]
            den = 1
            for x in v:
                den = den * x.denominator // np.gcd(den, x.denominator)
            vi = [int(x * den) for x in v]
            g = 0
            for x in vi:
                g = int(np.gcd(g, abs(x)))
            if g > 1:
                vi = [x // g for x in vi]
            out.append(vi)
        return out


def krylov_rows(seed_rows, mats, d):
    V = RowSpace(d)
    queue = []
    for r in seed_rows:
        if V.add(r):
            queue.append([Fraction(x) for x in r])
    while queue and V.dim() < d:
        r = queue.pop()
        for G in mats:
            r2 = rowmat(r, G)
            if V.add(r2):
                queue.append(r2)
    return V


def class_space(L, mats, mults):
    """V (RowSpace) and lambdas for operations mats with multipliers mults."""
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
    V = krylov_rows(seeds, mats, d)
    return V, lam


def dot(r, s):
    return sum(x * y for x, y in zip(r, s))


if __name__ == "__main__":
    # self-test on the ray step a of the one-mark system: C(a; 2) is everything with e = g (after one step)
    a = [[1, -1, 1, 0], [0, 0, 0, 1], [0, 0, 2, -1], [0, 0, 0, 1]]
    L = [1, 0, 0, 0]
    print("charpoly(a) =", charpoly(a), "integer eigenvalues", int_eigs(a))
    V, lam = class_space(L, [a], [2])
    print("dim V =", V.dim(), "rows", V.rows, "lambda", lam)
    print("C basis:", V.annihilator())
