#!/usr/bin/env python3
"""audit H: independent core routines (no imports from the author's modules).

Map: T(x) = (m_i x + r_i)/d on x = i mod d.  Pair chain u = M v + e; v's digit j, u's digit i = (M j + e) mod d;
M' = M m_i/m_j, e' = (m_i e + r_i - (m_i/m_j) M r_j)/d.  Debt x = exponent vector of M over the primes dividing the m's
(sorted increasingly) -- for the maps audited here the m's are products of distinct primes so this is the debt lattice in
its standard integer coordinates.  Root of digit j under coupling pi: v_pi(j) - v_j with v_j = exponent vector of m_j.
"""
from fractions import Fraction as Fr
import math


def factor(n):
    f = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            f[p] = f.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        f[n] = f.get(n, 0) + 1
    return f


def primes_of(m):
    return sorted({p for x in m for p in factor(x)})


def expvec(x, primes):
    """exponent vector of a positive rational x over the given primes (asserts no other primes occur)."""
    x = Fr(x)
    fn, fd = factor(x.numerator), factor(x.denominator)
    for p in list(fn) + list(fd):
        assert p in primes, (x, p)
    return tuple(fn.get(p, 0) - fd.get(p, 0) for p in primes)


def std_r(m):
    d = len(m)
    return [(-m[i] * i) % d for i in range(d)]


def lattice_rank(vecs):
    """rank over Q of a list of integer vectors (fraction-free Gaussian elimination)."""
    rows = [list(map(Fr, v)) for v in vecs if any(v)]
    if not rows:
        return 0
    n = len(rows[0])
    rk = 0
    for c in range(n):
        piv = None
        for i in range(rk, len(rows)):
            if rows[i][c] != 0:
                piv = i
                break
        if piv is None:
            continue
        rows[rk], rows[piv] = rows[piv], rows[rk]
        for i in range(len(rows)):
            if i != rk and rows[i][c] != 0:
                f = rows[i][c] / rows[rk][c]
                rows[i] = [a - f * b for a, b in zip(rows[i], rows[rk])]
        rk += 1
    return rk


def gen_subgroup(gens, mod):
    H = {1 % mod}
    todo = [1 % mod]
    while todo:
        h = todo.pop()
        for g in gens:
            y = h * g % mod
            if y not in H:
                H.add(y)
                todo.append(y)
    return H


def ratio_group(m, mod):
    gens = [(a * pow(b, -1, mod)) % mod for a in m for b in m]
    return gen_subgroup(gens, mod)


def residue(x, mod):
    """x a Fraction with denominator prime to mod -> x mod `mod`."""
    x = Fr(x)
    return (x.numerator * pow(x.denominator, -1, mod)) % mod


class Map:
    def __init__(self, d, m, r=None):
        self.d, self.m = d, list(m)
        self.r = list(r) if r is not None else std_r(m)
        for i in range(d):
            assert math.gcd(self.m[i], d) == 1
            assert (self.m[i] * i + self.r[i]) % d == 0
        self.primes = primes_of(self.m)
        self.v = [expvec(x, self.primes) for x in self.m]
        self.n = len(self.primes)
        self.rank = lattice_rank([tuple(a - b for a, b in zip(self.v[i], self.v[0])) for i in range(d)])
        self.Lam = sum(math.log(x / d) for x in self.m) / d

    def roots(self, a, b):
        d = self.d
        return [tuple(p - q for p, q in zip(self.v[(a * j + b) % d], self.v[j])) for j in range(d)]

    def Dmat(self, a, b):
        n = self.n
        D = [[0] * n for _ in range(n)]
        for xi in self.roots(a, b):
            for p in range(n):
                if xi[p]:
                    for q in range(n):
                        D[p][q] += xi[p] * xi[q]
        return D

    # exact pair-chain step on Fractions
    def step(self, M, e, j):
        d, m, r = self.d, self.m, self.r
        i = (residue(M, d) * j + residue(e, d)) % d
        N = m[i] * e + r[i] - Fr(m[i], m[j]) * M * r[j]
        e2 = N / d
        assert residue(N, d) == 0
        assert math.gcd(e2.denominator, d) == 1
        return M * Fr(m[i], m[j]), e2, i

    # memoized residue recursion A_s(h) (the recursion of THM-4611, written independently)
    def A_rec(self, s, Mres, eres, memo=None):
        if memo is None:
            memo = {}
        key = (s, Mres, eres)
        if key in memo:
            return memo[key]
        d, m, r, n = self.d, self.m, self.r, self.n
        a, b = Mres % d, eres % d
        D = self.Dmat(a, b)
        if s == 1:
            memo[key] = D
            return D
        mod, pmod = d ** s, d ** (s - 1)
        tot = [[d ** (s - 1) * D[p][q] for q in range(n)] for p in range(n)]
        for j in range(d):
            i = (a * j + b) % d
            ratio = m[i] * pow(m[j], -1, mod) % mod
            N = (m[i] * eres + r[i] - ratio * Mres * r[j]) % mod
            assert N % d == 0
            e2 = (N // d) % pmod
            M2 = ratio * Mres % pmod
            sub = self.A_rec(s - 1, M2, e2, memo)
            for p in range(n):
                for q in range(n):
                    tot[p][q] += sub[p][q]
        memo[key] = tot
        return tot


def det_int(M):
    n = len(M)
    if n == 1:
        return M[0][0]
    if n == 2:
        return M[0][0] * M[1][1] - M[0][1] * M[1][0]
    return sum((-1) ** c * M[0][c] * det_int([row[:c] + row[c + 1:] for row in M[1:]]) for c in range(n))


def adj_int(Q):
    n = len(Q)
    out = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            minor = [[Q[x][y] for y in range(n) if y != i] for x in range(n) if x != j]
            out[i][j] = (-1) ** (i + j) * det_int(minor)
    return out


def charpoly_pd(S):
    """Exact PD test for a symmetric integer matrix via the elementary symmetric functions of its eigenvalues
    (all principal-minor sums e_1..e_n > 0  <=>  PD, for symmetric matrices; Descartes' rule on a real-rooted polynomial).
    Independent of the leading-minor (Sylvester) test used by the author."""
    import itertools
    n = len(S)
    for k in range(1, n + 1):
        e = 0
        for idx in itertools.combinations(range(n), k):
            e += det_int([[S[a][b] for b in idx] for a in idx])
        if e <= 0:
            return False
    return True
