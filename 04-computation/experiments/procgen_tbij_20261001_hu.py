#!/usr/bin/env python3
"""procgen_tbij_20261001_hu.py -- Higashitani-Ueyama (arXiv:2409.10904) switching classes in Alt_n(Z/l) and
modular Eulerian matrices: s_{l,n} = t_{l,n} for all l (Brauer's permutation lemma for finite abelian
groups).  Both sides computed by Burnside, with fixed-point counts from Smith normal forms over Z:
  |A^sigma|        = |coker [ (sigma - 1) | X_1..X_n | l*I ]|           (A = Alt_n(Z/l) / <X_v>)
  |(Ker phi)^sigma| = #{N in (Z/l)^m : (sigma - 1) N = 0, row sums of N = 0}
Coordinates: N <-> (N_ij)_{i<j}; sigma acts by (sigma N)_ij = N_{sigma^-1 i, sigma^-1 j} (signed permutation).
"""
import itertools
import math
import os
import sys
from fractions import Fraction

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import procgen_tbij_20261001_lib as L  # noqa: E402
import procgen_tbij_20261001_nat as N  # noqa: E402


def signed_perm_matrix(sigma, n):
    E = L.edge_list(n)
    idx = L.edge_index(n)
    m = len(E)
    sinv = L.inverse(sigma)
    P = [[0] * m for _ in range(m)]
    for k, (i, j) in enumerate(E):
        a, b = sinv[i], sinv[j]
        if a < b:
            P[k][idx[(a, b)]] = 1
        else:
            P[k][idx[(b, a)]] = -1
    return P


def X_vectors(n):
    E = L.edge_list(n)
    vecs = []
    for v in range(n):
        x = []
        for (i, j) in E:
            if i == v:
                x.append(-1)
            elif j == v:
                x.append(1)
            else:
                x.append(0)
        vecs.append(x)
    return vecs


def rowsum_matrix(n):
    E = L.edge_list(n)
    R = [[0] * len(E) for _ in range(n)]
    for k, (i, j) in enumerate(E):
        R[i][k] += 1     # N_ij contributes to row i
        R[j][k] -= 1     # N_ji = -N_ij contributes to row j
    return R


def invariant_factors(rows, ncols):
    if not rows:
        return []
    M = Matrix(rows)
    S = smith_normal_form(M, domain=ZZ)
    d = []
    for i in range(min(S.shape)):
        if S[i, i] != 0:
            d.append(abs(int(S[i, i])))
    return d


def coker_order_A(sigma, n, l):
    m = len(L.edge_list(n))
    P = signed_perm_matrix(sigma, n)
    cols = []
    for c in range(m):
        cols.append([P[r][c] - (1 if r == c else 0) for r in range(m)])
    cols += X_vectors(n)
    for c in range(m):
        cols.append([l if r == c else 0 for r in range(m)])
    # matrix with these columns: m x (#cols); cokernel order = product of invariant factors (full rank)
    rows = [[cols[c][r] for c in range(len(cols))] for r in range(m)]
    d = invariant_factors(rows, len(cols))
    assert len(d) == m
    return math.prod(d)


def fixed_eulerian(sigma, n, l):
    m = len(L.edge_list(n))
    P = signed_perm_matrix(sigma, n)
    B = [[P[r][c] - (1 if r == c else 0) for c in range(m)] for r in range(m)] + rowsum_matrix(n)
    d = invariant_factors(B, m)
    r = len(d)
    cnt = l ** (m - r)
    for di in d:
        cnt *= math.gcd(di, l)
    return cnt


def s_and_t(n, l):
    s = Fraction(0)
    t = Fraction(0)
    for mu in L.partitions(n):
        g = N.standard_perm(mu, n)
        z = L.zee(mu)
        s += Fraction(coker_order_A(g, n, l), z)
        t += Fraction(fixed_eulerian(g, n, l), z)
    assert s.denominator == 1 and t.denominator == 1
    return int(s), int(t)


def brute_t(n, l):
    """iso classes of modular Eulerian matrices by brute force (small n)"""
    E = L.edge_list(n)
    m = len(E)
    seen = set()
    count = 0
    perms = list(itertools.permutations(range(n)))
    for vals in itertools.product(range(l), repeat=m):
        Mx = {}
        for (i, j), v in zip(E, vals):
            Mx[(i, j)] = v
            Mx[(j, i)] = (-v) % l
        if any(sum(Mx[(i, j)] for j in range(n) if j != i) % l for i in range(n)):
            continue
        key = min(tuple(Mx[(p[i], p[j])] for (i, j) in E) for p in perms)
        if key not in seen:
            seen.add(key)
    return len(seen)


def brute_s(n, l):
    """switching classes up to iso by brute force (small n): canonical key = min over relabellings of the
    min over the switching group"""
    E = L.edge_list(n)
    m = len(E)
    Xs = X_vectors(n)
    imgs = set()
    for a in itertools.product(range(l), repeat=n):
        imgs.add(tuple(sum(a[v] * Xs[v][k] for v in range(n)) % l for k in range(m)))
    imgs = list(imgs)
    perms = list(itertools.permutations(range(n)))
    idx = L.edge_index(n)
    seen = set()
    done = set()
    for vals in itertools.product(range(l), repeat=m):
        if vals in done:
            continue
        orbit = set()
        for p in perms:
            pinv = L.inverse(p)
            # relabel: new_ij = old_{p^-1 i, p^-1 j}
            rel = []
            for (i, j) in E:
                a, b = pinv[i], pinv[j]
                rel.append(vals[idx[(a, b)]] if a < b else (-vals[idx[(b, a)]]) % l)
            for x in imgs:
                orbit.add(tuple((r + y) % l for r, y in zip(rel, x)))
        done |= orbit
        seen.add(min(orbit))
    return len(seen)
