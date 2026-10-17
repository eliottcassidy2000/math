#!/usr/bin/env python3
"""S39 audit, item 1: independent recomputation of the Kloosterman moments M_k(p), k <= 8.

Method A (exact, combinatorial; no Kloosterman sums at all): count closed walks in G_p = Cay(F_(p^2), mu_(p+1)) by a
radial dynamic programme.  F_(p^2) is modelled as F_p[th]/(th^2 - th - r) (a model different from the session's
x + y sqrt(n)); N(a + b th) = a^2 + ab - r b^2.  R_k(c) = number of k-walks from 0 to a point of norm c (c = 0 means the
vertex 0 itself; norm 0 <=> point 0).  R_(k+1)(c) = sum_(u in mu) R_k(N(d_c - u)) for a fixed representative d_c of
norm c.  W_k = R_k(0) and M_k = (-1)^k (p^2 W_k - (p+1)^k)/(p+1) (THM-4663 (1), which is THM-4628 + vertex
transitivity).  The DP is run modulo three primes near 2^31 and combined by CRT (exact: W_k < (p+1)^k < product).
Method B (numerical, independent): mpmath with 40 digits, K(a) = sum_x cos(2 pi (x + a/x)/p) (K real), M_k = sum K^k,
rounded; for every prime p < 400.
Both are compared with the session's JSON for every prime in range.
Usage: python -B a1_kloosterman_moments_independent.py SESSION_JSON PMAX_EXACT PMAX_MP
"""
import sys, json, time
import numpy as np
from numba import njit
import sympy
import mpmath as mp

MODS = [2147483629, 2147483587, 2147483579]          # primes below 2^31
for m in MODS:
    assert sympy.isprime(m)


def field_model(p):
    """r with th^2 - th - r irreducible mod p: disc 1 + 4r a non-residue."""
    for r in range(1, p):
        if pow((1 + 4 * r) % p, (p - 1) // 2, p) == p - 1:
            return r
    raise RuntimeError


@njit
def build_table(p, r):
    # all elements (a, b) -> norm
    mu_a = np.empty(p + 1, dtype=np.int64); mu_b = np.empty(p + 1, dtype=np.int64)
    rep_a = -np.ones(p, dtype=np.int64); rep_b = -np.ones(p, dtype=np.int64)
    cnt = 0
    for a in range(p):
        for b in range(p):
            n = (a * a + a * b - r * b * b) % p
            if n == 1:
                mu_a[cnt] = a; mu_b[cnt] = b; cnt += 1
            if rep_a[n] < 0:
                rep_a[n] = a; rep_b[n] = b
    T = np.empty((p, p + 1), dtype=np.int64)
    for c in range(p):
        da = rep_a[c]; db = rep_b[c]
        for j in range(p + 1):
            a = (da - mu_a[j]) % p; b = (db - mu_b[j]) % p
            T[c, j] = (a * a + a * b - r * b * b) % p
    return cnt, T


@njit
def walks_mod(p, T, kmax, m):
    R = np.zeros(p, dtype=np.int64); R[0] = 1
    W = np.zeros(kmax + 1, dtype=np.int64); W[0] = 1
    for k in range(1, kmax + 1):
        Rn = np.zeros(p, dtype=np.int64)
        for c in range(p):
            s = 0
            for j in range(p + 1):
                s += R[T[c, j]]
            Rn[c] = s % m
        R = Rn
        W[k] = R[0]
    return W


def exact_W(p, kmax=8):
    r = field_model(p)
    cnt, T = build_table(p, r)
    assert cnt == p + 1
    res = [walks_mod(p, T, kmax, m) for m in MODS]
    out = []
    for k in range(kmax + 1):
        val = int(sympy.ntheory.modular.crt(MODS, [int(x[k]) for x in res])[0])
        assert 0 <= val < (p + 1) ** k + 1
        out.append(val)
    return out


def moments_from_W(p, W):
    M = []
    for k in range(1, len(W)):
        num = p * p * W[k] - (p + 1) ** k
        assert num % (p + 1) == 0, (p, k)
        M.append((-1) ** k * (num // (p + 1)))
    return M


def mp_moments(p, kmax=8):
    mp.mp.dps = 40
    cosv = [mp.cos(2 * mp.pi * j / p) for j in range(p)]
    inv = [0] + [pow(x, -1, p) for x in range(1, p)]
    M = [mp.mpf(0)] * kmax
    for a in range(1, p):
        K = mp.fsum(cosv[(x + a * inv[x]) % p] for x in range(1, p))
        pw = mp.mpf(1)
        for k in range(kmax):
            pw *= K
            M[k] += pw
    out = []
    for k in range(kmax):
        v = int(mp.nint(M[k]))
        assert abs(M[k] - v) < mp.mpf(10) ** -15, (p, k)
        out.append(v)
    return out


if __name__ == '__main__':
    J = json.load(open(sys.argv[1]))['moments']
    PE, PM = int(sys.argv[2]), int(sys.argv[3])
    t0 = time.time()
    bad, n = [], 0
    mine = {}
    for p in sympy.primerange(5, PE + 1):
        M = moments_from_W(p, exact_W(p))
        mine[p] = M
        n += 1
        if str(p) in J and J[str(p)] != M:
            bad.append(p)
    print('Method A (exact radial walk DP mod 3 primes + CRT): %d primes 5 <= p <= %d; disagreements with session JSON: %s (%.0f s)'
          % (n, PE, bad, time.time() - t0))
    print('  e.g. p = 5: %s; p = 7: %s' % (mine[5], mine[7]))
    assert not bad
    t0 = time.time()
    bad, n = [], 0
    for p in sympy.primerange(5, PM):
        Mm = mp_moments(p)
        n += 1
        if Mm != mine[p] or J[str(p)] != Mm:
            bad.append(p)
    print('Method B (mpmath 40-digit direct Kloosterman sums): %d primes 5 <= p < %d; disagreements: %s (%.0f s)'
          % (n, PM, bad, time.time() - t0))
    assert not bad
    # classical checks M_1 = 1, M_2 = p^2 - p - 1 on all exact data
    assert all(mine[p][0] == 1 and mine[p][1] == p * p - p - 1 for p in mine)
    json.dump({str(p): v for p, v in mine.items()}, open(sys.argv[4], 'w')) if len(sys.argv) > 4 else None
    print('ALL CHECKS PASSED')
