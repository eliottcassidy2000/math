#!/usr/bin/env python3
"""
procgen_cubedist_20260925_prune.py -- smaller explicit flip sets inside Bad_k (upper bounds on delta_k),
cube-distance lane 2026-09-25 (HYP-9138).

Greedy pruning: start from R = Bad_k (class (i) by Theorem 1) and try to restore sigma = + at each r in R in a
fixed order, keeping the change whenever G_sigma still has no expanding cycle (exact test: the C engine's
potential certificate at threshold q*/r* = best lower approximation of log_3 2 with denominator <= 2^k, which is
equivalent to 'no expanding simple cycle'; the final set is re-verified by an integer potential that is checked
edge by edge in Python, and by exact Karp for k <= 15).  Several orders are tried; the best set is reported.
Also reported: the one-flip-per-necklace rules, which FAIL (a lower-bound-tight construction is not this simple).
"""
import sys
import os
import math
import random
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import procgen_cubedist_20260925_lib as L    # noqa: E402
import procgen_cubedist_20260925_bad as B    # noqa: E402

A = math.log(1.5)
Bw = math.log(2)


def classi(k, fl, thr):
    st, _ = L.certificate(k, L.mask_of(fl), thr.numerator, thr.denominator)
    return st == 'OK'


def ballot_margin(r, k):
    """min over j = 1..k of the Collatz partial log-multiplier of r (how robustly r is undecided)"""
    n, s, m = r, 0.0, 1e9
    for _ in range(k):
        if n & 1:
            s += A
            n = (3 * n + 1) // 2
        else:
            s -= Bw
            n //= 2
        m = min(m, s)
    return m


def prune(k, orders=('asc', 'desc', 'margin_lo', 'margin_hi', 'rand1', 'rand2'), log=None):
    thr = L.expanding_threshold(k)
    bad = [int(r) for r in B.bad_numpy(k)[0]]
    L.check(classi(k, bad, thr), "Bad_k not class (i)")
    best = None
    stats = {}
    for name in orders:
        order = list(bad)
        if name == 'asc':
            order.sort()
        elif name == 'desc':
            order.sort(reverse=True)
        elif name == 'margin_lo':
            order.sort(key=lambda r: (ballot_margin(r, k), r))
        elif name == 'margin_hi':
            order.sort(key=lambda r: (-ballot_margin(r, k), r))
        else:
            random.Random(1000 * k + sum(map(ord, name))).shuffle(order)
        cur = set(bad)
        for r in order:
            cur.discard(r)
            if not classi(k, cur, thr):
                cur.add(r)
        stats[name] = len(cur)
        if best is None or len(cur) < len(best):
            best = set(cur)
    return best, stats, len(bad)


def verify(k, fl):
    F = L.best_lower_approx(1 << k)
    # certificate at the exact threshold (every simple cycle has p <= 2^k)
    st, psi = L.certificate(k, L.mask_of(fl), F.numerator, F.denominator)
    L.check(st == 'OK', "pruned set has an expanding cycle")
    L.check(L.verify_certificate(k, set(fl), psi, F.numerator, F.denominator), "certificate check failed")
    d = None
    if k <= 15:
        d = L.karp_density(k, L.mask_of(fl))
        L.check(3 ** d.numerator < 2 ** d.denominator, "Karp: not class (i)")
    return d


def one_per_necklace(k, rule):
    """one flip per expanding necklace: the ballot rotation with the lowest ('low') or highest ('high') starting
    level in the periodic walk of the necklace (a necklace with more than ck ones has a ballot rotation by the
    cycle lemma).  Returns (|R|, class (i)?)."""
    K = 1 << k
    word = {r: tuple(L.parity_word(r, k)) for r in range(K)}
    inv = {w: r for r, w in word.items()}
    seen = set()
    R = set()
    for r in range(1, K, 2):
        w = word[r]
        if w in seen:
            continue
        rots = [w[i:] + w[:i] for i in range(k)]
        seen.update(rots)
        if not 3 ** sum(w) > 2 ** k:
            continue
        # exact levels: S_i = a_i log 3 - i log 2 compared through (a_i, i); ballot start i: S_(i+j) > S_i, j=1..k
        pref = [0]
        for b in w + w:
            pref.append(pref[-1] + b)
        cands = []
        for i in range(k):
            ok = all(3 ** (pref[i + j] - pref[i]) > 2 ** j for j in range(1, k + 1))
            if ok:
                cands.append(i)
        L.check(cands, "expanding necklace without a ballot rotation")
        # compare levels S_i = pref[i] log 3 - i log 2 exactly: S_i < S_i' iff 3^(pref[i]) 2^(i') < 3^(pref[i']) 2^i
        key = lambda i: (pref[i] * math.log(3) - i * math.log(2))
        i0 = min(cands, key=key) if rule == 'low' else max(cands, key=key)
        R.add(inv[rots[i0]])
    thr = L.expanding_threshold(k)
    return len(R), classi(k, R, thr)


def main(log=print, kmax=14):
    log("P1. greedy pruning inside Bad_k (upper bounds on delta_k; every set verified class (i))")
    log("   k   N_k  pruned  |Bad_k|  pruned/2^(k-1)  pruned/N_k   rho_max   orders tried -> sizes")
    rows = {}
    for k in range(3, kmax + 1):
        t0 = time.time()
        orders = ('asc', 'desc', 'margin_lo', 'margin_hi', 'rand1', 'rand2')
        best, stats, nb = prune(k, orders)
        d = verify(k, best)
        N = L.necklace_lower_bound(k)
        rows[k] = (len(best), sorted(best))
        log(f"  {k:2d} {N:5d}  {len(best):6d}  {nb:7d}   {len(best) / 2 ** (k - 1):.5f}       {len(best) / N:6.3f}    "
            f"{str(d) if d is not None else '-':>6}   {stats}  [{time.time() - t0:.0f}s]")
    ks = sorted(rows)
    # least-squares fit of log2(pruned / 2^(k-1)) against k and log2 k over k >= 8
    import numpy as np
    xs = [k for k in ks if k >= 8]
    y = np.array([math.log2(rows[k][0] / 2 ** (k - 1)) for k in xs])
    X = np.array([[1.0, k, math.log2(k)] for k in xs])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    X2 = np.array([[1.0, k] for k in xs])
    coef2, *_ = np.linalg.lstsq(X2, y, rcond=None)
    log(f"  fit log2(pruned/2^(k-1)) ~ {coef2[0]:.3f} {coef2[1]:+.4f} k   (pure exponential, k = 8..{kmax})")
    log(f"  fit log2(pruned/2^(k-1)) ~ {coef[0]:.3f} {coef[1]:+.4f} k {coef[2]:+.3f} log2 k   "
        f"(theory: slope -(1-h) = -0.0500, log-coefficient -1.5)")
    log("")
    log("P2. one flip per expanding necklace (|R| = N_k): lowest / highest ballot rotation -- class (i)?")
    for k in range(3, 15):
        nl, okl = one_per_necklace(k, 'low')
        nh, okh = one_per_necklace(k, 'high')
        L.check(nl == nh == L.necklace_lower_bound(k), "one-per-necklace size != N_k")
        log(f"  k={k:2d}: |R| = N_k = {nl:4d}; lowest: {'class (i)' if okl else 'expanding cycle left'}; "
            f"highest: {'class (i)' if okh else 'expanding cycle left'}")
    log("")
    log("P3. flipping only the long rising runs (the 2-adic neighbourhoods of -1, and of -1 and -5) fails for j >= 3:")
    log("    R_j = {r = -1 mod 2^j}, R'_j = R_j u {r = -5 mod 2^j}; an expanding cycle survives (e.g. (1^(j-1) 0)^inf)")
    for k in (8, 10, 12):
        thr = L.expanding_threshold(k)
        row = []
        for j in range(2, k + 1):
            Rj = set(r for r in range(1, 1 << k, 2) if (r + 1) % (1 << j) == 0)
            Rj2 = Rj | set(r for r in range(1, 1 << k, 2) if (r + 5) % (1 << j) == 0)
            row.append((j, classi(k, Rj, thr), classi(k, Rj2, thr)))
        good = [(j, a, b) for j, a, b in row if a or b]
        L.check(all((not a or j == 2) and (not b or j <= 4) for j, a, b in good),
                "a run-based flip set of depth >= 5 is class (i)")
        log(f"  k={k:2d}: (j, R_j class (i)?, R'_j class (i)?) true only for {good}; R_2 = R'_3 = all residues 3 mod 4")
        log("         (the all-d strategy, Haar 1/2); R'_4 = {11, 15 mod 16} is the lifted level-4 optimum (Haar 1/4)")
    return rows


if __name__ == '__main__':
    main(kmax=int(sys.argv[1]) if len(sys.argv) > 1 else 14)
