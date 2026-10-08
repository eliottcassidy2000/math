#!/usr/bin/env python3
"""Audit G: direct simulation of Matthews-Watts pair orbits (stdlib only).

T(x) = (m_i x + r_i)/d on x = i mod d.  Haar y is realised as an integer uniform mod d^(T+pad): the first T branch
digits of y (and of y + e) are then exactly Haar-distributed (d-adic Terras property), and the integer orbit IS the
d-adic orbit.  We run u_t = T^t(y + e), v_t = T^t(y) as big integers, record the first equal-time merge u_t == v_t,
and track the debt M_t = prod m_(i_s)/m_(j_s) as an exact prime-exponent vector (M_t = 1 iff the vector is 0).

Also: an exact pair chain (M, e) driven by fresh uniform digits, with Fractions, for cross-checking the update
    M' = M m_i/m_j,  e' = (m_i e + r_i - (m_i/m_j) M r_j)/d,  i = M j + e mod d,
against the direct orbits (u_t = M_t v_t + e_t exactly)."""
import random
from fractions import Fraction as Fr

def factor(n):
    f = {}; p = 2
    while p * p <= n:
        while n % p == 0: f[p] = f.get(p, 0) + 1; n //= p
        p += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

class MW:
    def __init__(self, d, m, r):
        self.d, self.m, self.r = d, list(m), list(r)
        for i in range(d):
            assert (self.m[i] * i + self.r[i]) % d == 0, (i, m, r)
            assert all(self.m[i] % p for p in factor(d)), 'm_i must be prime to d'
        primes = sorted({p for x in self.m for p in factor(x)})
        self.primes = primes
        self.E = [[factor(x).get(p, 0) for p in primes] for x in self.m]
        # step vectors for (i, j): E[i] - E[j]
        self.step = [[tuple(a - b for a, b in zip(self.E[i], self.E[j])) for j in range(d)] for i in range(d)]
        self.zero = tuple(0 for _ in primes)

    def T(self, x):
        i = x % self.d
        return (self.m[i] * x + self.r[i]) // self.d

def run_pairs(mw, e, Tmax, N, seed, checkpoints, pad=64, window_visits=True, k0=None):
    """Direct orbits.  Returns dict: merged_by[T] counts, visits stats.  k0: optional initial debt exponent multiplier
    u0 = M0*y + e with M0 = prod m^k0 (not used for translations)."""
    d = mw.d; m = mw.m; r = mw.r; step = mw.step
    rnd = random.Random(seed)
    mod = d ** (Tmax + pad)
    cps = sorted(checkpoints)
    merged_by = {c: 0 for c in cps}
    merge_times = []
    # debt-visit statistics: visits to M = 1 in (c/4, c] for chains unmerged at c/4
    vis_win = {c: [0, 0] for c in cps}      # [sum of visits, number of chains unmerged at c/4]
    any_late_visit = {c: 0 for c in cps}     # chains with a debt visit in (c/2, c] and unmerged at c/2
    unmerged_half = {c: 0 for c in cps}
    nonzero_debt_merges = 0
    for s in range(N):
        y = rnd.randrange(mod)
        u = y + e; v = y
        k = mw.zero; nz = False
        tm = None
        visits = []          # times t (1-based) with debt 0 after step t, before merge
        for t in range(1, Tmax + 1):
            i = u % d; j = v % d
            u = (m[i] * u + r[i]) // d
            v = (m[j] * v + r[j]) // d
            if i != j:
                st = step[i][j]
                if any(st):
                    k = tuple(a + b for a, b in zip(k, st))
            if u == v:
                tm = t; break
            if window_visits and not any(k):
                visits.append(t)
        if tm is not None:
            merge_times.append(tm)
            if any(k): nonzero_debt_merges += 1
        for c in cps:
            if tm is not None and tm <= c: merged_by[c] += 1
            if window_visits:
                if tm is None or tm > c // 4:
                    vis_win[c][1] += 1
                    vis_win[c][0] += sum(1 for t in visits if c // 4 < t <= c)
                if tm is None or tm > c // 2:
                    unmerged_half[c] += 1
                    if any(c // 2 < t <= c for t in visits): any_late_visit[c] += 1
    return dict(N=N, merged_by=merged_by, merge_times=merge_times, vis_win=vis_win,
                any_late_visit=any_late_visit, unmerged_half=unmerged_half, nonzero_debt_merges=nonzero_debt_merges)

def dres(x, d):
    """residue of a rational d-adic integer x mod d"""
    x = Fr(x)
    return (x.numerator * pow(x.denominator, -1, d)) % d

def pair_chain_check(mw, e0, steps, samples, seed, pad=64):
    """Cross-check the pair-chain update against direct orbits: returns number of mismatches."""
    d = mw.d; rnd = random.Random(seed); bad = 0
    for _ in range(samples):
        y = rnd.randrange(d ** (steps + pad)); u = y + e0; v = y
        M = Fr(1); e = Fr(e0)
        for t in range(steps):
            j = v % d; i_pred = (dres(M, d) * j + dres(e, d)) % d; i = u % d
            if i != i_pred: bad += 1; break
            Mn = M * Fr(mw.m[i], mw.m[j])
            e = (mw.m[i] * e + mw.r[i] - Fr(mw.m[i], mw.m[j]) * M * mw.r[j]) / d
            M = Mn
            u = mw.T(u); v = mw.T(v)
            if u != M * v + e: bad += 1; break
    return bad

def fmt_run(res, cps):
    N = res['N']; out = []
    for c in cps:
        q = res['merged_by'][c] / N
        se = (q * (1 - q) / N) ** 0.5
        vw = res['vis_win'][c]; vwin = vw[0] / vw[1] if vw[1] else float('nan')
        la = res['any_late_visit'][c] / res['unmerged_half'][c] if res['unmerged_half'][c] else float('nan')
        out.append(f"T={c}: merged {q:.4f}+-{se:.4f}, no-merge {1-q:.4f}; visits/(unmerged@T/4) in (T/4,T] {vwin:.3f}; "
                   f"P(debt visit in (T/2,T] | unmerged@T/2) {la:.3f}")
    return out
