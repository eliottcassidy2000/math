#!/usr/bin/env python3
"""S91 helper (census_k2), step 4: Monte Carlo estimate of the BY-VALUE statistics of the root-change words at
scales far beyond any census (value windows [2^n, 2^(n+1)), n up to 800), to see the limits of the local factors.

Estimator (Knuth's tree-size estimator on the tail tree; unbiased for every sum over counted words of the window):
  K = 2: draw (LA, LB) uniformly in {LA + LB <= n}, uniform letters for A and B; the word A J B J (empty tail) has
         weight 2^(LA+LB); then walk a uniform random O/Z path C0 below A J B J: the prefix after k letters is the
         counted word A J B J C0 Z, value F = D - E + U, weight 2^(LA+LB+k); stop when F >= 2^(n+1) (F increases).
  K = 1: the same with one block A (LA uniform in [0, n]).
Every exactly-K-J counted word of the window is reached with probability 1/(#(LA, LB) * its weight), so weighted
sums estimate exact counts.  Real state in float64 (relative error 1e-13: window membership is affected for a
fraction 1e-13 of the words); the state mod p is carried exactly in integers, so the residue class, the confinement
events C1, S2 (p | gcd(D, U) before a root change) and T (p | gcd(x3, y3), read off the tail's action on the basis
states) are exact.  Standard errors from independent batches.

Usage:  python3 mc_byvalue.py validate            compare with the exact census on [2^23, 2^24) (K = 1 and K = 2)
        python3 mc_byvalue.py scan [nsamp_millions]   n = 30, 60, 120, 240, 480, 800; prints the tables used in the report
NUMERICAL (Monte Carlo).  Exit 1 when the validation against the exact census fails.
"""
import math
import sys
import time

import numpy as np
from numba import njit, prange

PRIMES = np.array([3, 5, 7, 11, 13, 17], dtype=np.int64)
MAXP = 17
W1 = 0.8035922284
W0 = 1 - W1
NFL = 8          # flag accumulators: 0 total weight, 1 k1, 2 l2, 3 k1&l2, 4 count of contributions, 5 sum w^2, 6 sheet0 (K=1), 7 spare


@njit(cache=True)
def one_sample(K, n0, primes, fl, cls, ccls, cat):
    """one random path; adds its weighted contributions to fl[NFL], cls[p, r], ccls[p, r], cat[p, S, T, C1, st]."""
    npr = primes.shape[0]
    X = 2.0 ** n0
    X2 = 2.0 ** (n0 + 1)
    if K == 2:
        while True:
            LA = np.random.randint(0, n0 + 1)
            LB = np.random.randint(0, n0 + 1)
            if LA + LB <= n0:
                break
    else:
        LA = np.random.randint(0, n0 + 1)
        LB = 0
    D = 1.0
    E = 0.0
    U = 2.0
    d = 1.0
    sD = np.ones(npr, np.int64)
    sE = np.zeros(npr, np.int64)
    sU = np.full(npr, 2, np.int64)
    sd = np.ones(npr, np.int64)
    bits = 0
    nbits = 0
    last = 0
    for i in range(LA):
        if nbits == 0:
            bits = np.random.randint(0, 1 << 62)
            nbits = 62
        b = bits & 1
        bits >>= 1
        nbits -= 1
        if b == 0:                       # O
            E = 0.0
            U = 2.0 * U
            for j in range(npr):
                sE[j] = 0
                sU[j] = (2 * sU[j]) % primes[j]
        else:                            # Z
            D = D - E + U
            E = d
            U = 2.0 * U - d
            for j in range(npr):
                p = primes[j]
                nd = (sD[j] - sE[j] + sU[j]) % p
                sU[j] = (2 * sU[j] - sd[j]) % p
                sE[j] = sd[j]
                sD[j] = nd
        last = b
    e1 = last                            # E1 = 1 iff A ends in Z
    C1 = np.zeros(npr, np.int64)
    S2 = np.zeros(npr, np.int64)
    for j in range(npr):
        if sD[j] == 0 and sU[j] == 0:
            C1[j] = 1
    firstB = 0
    lastB = 0
    k1 = 1
    if K == 2:
        # first root change
        nD = D - E + U
        nE = D - E
        nU = D + U
        d = D
        D = nD
        E = nE
        U = nU
        for j in range(npr):
            p = primes[j]
            a = (sD[j] - sE[j] + sU[j]) % p
            b2 = (sD[j] - sE[j]) % p
            c = (sD[j] + sU[j]) % p
            sd[j] = sD[j]
            sD[j] = a
            sE[j] = b2
            sU[j] = c
        for i in range(LB):
            if nbits == 0:
                bits = np.random.randint(0, 1 << 62)
                nbits = 62
            b = bits & 1
            bits >>= 1
            nbits -= 1
            if b == 0:
                E = 0.0
                U = 2.0 * U
                for j in range(npr):
                    sE[j] = 0
                    sU[j] = (2 * sU[j]) % primes[j]
            else:
                D = D - E + U
                E = d
                U = 2.0 * U - d
                for j in range(npr):
                    p = primes[j]
                    nd = (sD[j] - sE[j] + sU[j]) % p
                    sU[j] = (2 * sU[j] - sd[j]) % p
                    sE[j] = sd[j]
                    sD[j] = nd
            if i == 0:
                firstB = 1 + b
            lastB = 1 + b
        k1 = 1 if (e1 == 0 or firstB == 2) else 0
        for j in range(npr):
            if sD[j] == 0 and sU[j] == 0:
                S2[j] = 1
        z2 = 1 if lastB == 1 else 0
    else:
        for j in range(npr):
            S2[j] = C1[j]
            C1[j] = 0
        z2 = 1 if e1 == 0 else 0         # E = 0 before the root change
    Epre = np.empty(npr, np.int64)
    for j in range(npr):
        Epre[j] = sE[j]
    lw0 = LA + LB - n0
    # the word ending in the last J (empty tail)
    F = D - E + U
    if F >= X2:
        return
    if F >= X:
        w = 2.0 ** lw0
        l2 = z2
        fl[0] += w
        fl[4] += 1.0
        fl[5] += w * w
        if k1:
            fl[1] += w
        if l2:
            fl[2] += w
        if k1 and l2:
            fl[3] += w
        if l2 == 0:
            fl[6] += w
        for j in range(npr):
            p = primes[j]
            r = (sD[j] - sE[j] + sU[j]) % p
            cls[j, r] += w
            st = 1 if Epre[j] == 0 else 2
            cat[j, S2[j], 0, C1[j], st] += w
            if S2[j]:
                ccls[j, r] += w
    # last root change
    nD = D - E + U
    nE = D - E
    nU = D + U
    d = D
    D = nD
    E = nE
    U = nU
    for j in range(npr):
        p = primes[j]
        a = (sD[j] - sE[j] + sU[j]) % p
        b2 = (sD[j] - sE[j]) % p
        c = (sD[j] + sU[j]) % p
        sd[j] = sD[j]
        sD[j] = a
        sE[j] = b2
        sU[j] = c
    # basis states: J applied to (D, E, U) = (1, 0, 0) and to (0, 0, 1)
    aD = np.ones(npr, np.int64)
    aE = np.ones(npr, np.int64)
    aU = np.ones(npr, np.int64)
    ad = np.ones(npr, np.int64)
    bD = np.ones(npr, np.int64)
    bE = np.zeros(npr, np.int64)
    bU = np.ones(npr, np.int64)
    bd = np.zeros(npr, np.int64)
    k = 0
    firstC = 0
    while True:
        F = D - E + U
        if F >= X2:
            break
        if F >= X:
            w = 2.0 ** (lw0 + k)
            odd = 1 if (k == 0 or firstC == 2) else 0        # the tail C0 Z starts with Z
            l2 = 1 if (odd == 1 or z2 == 1) else 0
            fl[0] += w
            fl[4] += 1.0
            fl[5] += w * w
            if k1:
                fl[1] += w
            if l2:
                fl[2] += w
            if k1 and l2:
                fl[3] += w
            if l2 == 0:
                fl[6] += w
            for j in range(npr):
                p = primes[j]
                r = (sD[j] - sE[j] + sU[j]) % p
                cls[j, r] += w
                T = 0
                if (aD[j] - aE[j] + aU[j]) % p == 0 and (bD[j] - bE[j] + bU[j]) % p == 0:
                    T = 1
                if odd:
                    st = 0
                elif Epre[j] == 0:
                    st = 1
                else:
                    st = 2
                cat[j, S2[j], T, C1[j], st] += w
                if S2[j] or T:
                    ccls[j, r] += w
        if nbits == 0:
            bits = np.random.randint(0, 1 << 62)
            nbits = 62
        b = bits & 1
        bits >>= 1
        nbits -= 1
        if b == 0:
            E = 0.0
            U = 2.0 * U
            for j in range(npr):
                p = primes[j]
                sE[j] = 0
                sU[j] = (2 * sU[j]) % p
                aE[j] = 0
                aU[j] = (2 * aU[j]) % p
                bE[j] = 0
                bU[j] = (2 * bU[j]) % p
        else:
            D = D - E + U
            E = d
            U = 2.0 * U - d
            for j in range(npr):
                p = primes[j]
                nd = (sD[j] - sE[j] + sU[j]) % p
                sU[j] = (2 * sU[j] - sd[j]) % p
                sE[j] = sd[j]
                sD[j] = nd
                nd = (aD[j] - aE[j] + aU[j]) % p
                aU[j] = (2 * aU[j] - ad[j]) % p
                aE[j] = ad[j]
                aD[j] = nd
                nd = (bD[j] - bE[j] + bU[j]) % p
                bU[j] = (2 * bU[j] - bd[j]) % p
                bE[j] = bd[j]
                bD[j] = nd
        if k == 0:
            firstC = 1 + b
        k += 1


@njit(parallel=True, cache=True)
def run_batches(K, n0, nbatch, nsamp, seed, primes):
    npr = primes.shape[0]
    fl = np.zeros((nbatch, NFL))
    cls = np.zeros((nbatch, npr, MAXP))
    ccls = np.zeros((nbatch, npr, MAXP))
    cat = np.zeros((nbatch, npr, 2, 2, 2, 3))
    for bt in prange(nbatch):
        np.random.seed(seed + 7919 * bt)
        for it in range(nsamp):
            one_sample(K, n0, primes, fl[bt], cls[bt], ccls[bt], cat[bt])
    return fl, cls, ccls, cat


def ratio(num, den):
    """ratio of sums with a batch (jackknife) standard error; num, den are per-batch arrays."""
    nb = num.shape[0]
    tn, td = num.sum(), den.sum()
    est = tn / td
    jk = (tn - num) / (td - den)
    se = math.sqrt((nb - 1) / nb * ((jk - jk.mean()) ** 2).sum())
    return est, se


def summarize(K, n0, fl, cls, ccls, cat):
    """dictionary of by-value estimates (value, standard error) for exactly-K-J words on [2^n0, 2^(n0+1))."""
    out = {}
    W = fl[:, 0]
    out["ess"] = (W.sum() ** 2 / fl[:, 5].sum(), 0.0)
    out["k1"] = ratio(fl[:, 1], W)
    out["l2"] = ratio(fl[:, 2], W)
    out["k1l2"] = ratio(fl[:, 3], W)
    out["sheet0"] = ratio(fl[:, 6], W)
    for j, p in enumerate(PRIMES):
        p = int(p)
        c = cls[:, j, :p]
        out[("L0", p)] = tuple(p * v for v in ratio(c[:, 0], W))
        out[("Lm1", p)] = tuple(p * v for v in ratio(c[:, p - 1], W))
        oth = [r for r in range(p) if r not in (0, p - 1)]
        out[("Loth", p)] = tuple(p * v / len(oth) for v in ratio(c[:, oth].sum(axis=1), W))
        out[("Lnz", p)] = tuple(p * v / (p - 1) for v in ratio(c[:, 1:].sum(axis=1), W))
        num0 = c[:, 0]
        nz = c[:, 1:].sum(axis=1) / (p - 1)
        out[("rho", p)] = ratio(num0, nz)
        ct = cat[:, j]                      # [S, T, C1, st]
        out[("S", p)] = ratio(ct[:, 1].sum(axis=(1, 2, 3)), W)
        out[("T", p)] = ratio(ct[:, :, 1].sum(axis=(1, 2, 3)), W)
        out[("C1", p)] = ratio(ct[:, :, :, 1].sum(axis=(1, 2, 3)), W)
        out[("S|C1", p)] = ratio(ct[:, 1, :, 1].sum(axis=(1, 2)), ct[:, :, :, 1].sum(axis=(1, 2, 3)) + 1e-300)
        out[("S,-C1", p)] = ratio(ct[:, 1, :, 0].sum(axis=(1, 2)), W)
        out[("z|S", p)] = ratio(ct[:, 1, :, :, :2].sum(axis=(1, 2, 3)), ct[:, 1].sum(axis=(1, 2, 3)) + 1e-300)
        cc = ccls[:, j, :p]
        unc0 = c[:, 0] - cc[:, 0]
        uncnz = (c[:, 1:].sum(axis=1) - cc[:, 1:].sum(axis=1)) / (p - 1)
        out[("U0/Unz", p)] = ratio(unc0, uncnz)
        out[("conf", p)] = ratio(cc.sum(axis=1), W)
    return out


def limits_k2(p, sig):
    return 1 + (p - 1) * (sig * p - 1) / p ** 3, 1 - (sig * p - 1) / p ** 3


def validate():
    """MC at n0 = 23 against the exact census of exactly-one-J and exactly-two-J words on [2^23, 2^24)."""
    import rc_census as rc
    ok = True
    lo, hi = 23, 24
    primes = [int(p) for p in PRIMES]
    r1 = rc.census(1, 1 << lo, 1 << hi, primes, nthreads=3)
    r2 = rc.census(2, 1 << lo, 1 << hi, primes, nthreads=3)
    r0 = rc.census(0, 1 << lo, 1 << hi, primes, nthreads=3)
    nodd = (1 << hi) - (1 << lo) >> 1
    exact = {}
    for K, (ra, rb) in ((1, (r1, r0)), (2, (r2, r1))):
        totK = ra["tot"] - rb["tot"]
        words = int(totK.sum())
        ex = {"words": words}
        for j, p in enumerate(primes):
            means, n = rc.class_means(totK, 1 << lo, 1 << hi, p)
            L = means / (words / nodd)
            ex[("L0", p)] = L[0]
            ex[("Lnz", p)] = float(np.mean(L[1:]))
            c = ra["cnt"][j, K]
            g = ra["allp"][j, K].sum()
            ex[("S", p)] = c[1].sum() / g
            ex[("T", p)] = c[:, 1].sum() / g
        if K == 2:
            c2 = r2["cnt2"]
            ex["k1"] = c2[1].sum() / c2.sum()
            ex["l2"] = c2[:, 1].sum() / c2.sum()
        else:
            ex["sheet0"] = r1["allcnt"][1, 0, 0] / r1["allcnt"][1].sum()
        exact[K] = ex
    for K in (1, 2):
        t0 = time.time()
        fl, cls, ccls, cat = run_batches(K, lo, 60, 100000, 20261110 + K, PRIMES)
        s = summarize(K, lo, fl, cls, ccls, cat)
        ex = exact[K]
        nstrata = (lo + 1) * (lo + 2) / 2 if K == 2 else (lo + 1)
        est_words = fl[:, 0].sum() * 2.0 ** lo * nstrata / (60 * 100000)
        print(f"validation K = {K}, window [2^23, 2^24): exact words {ex['words']}, MC estimate {est_words:.4g} "
              f"(ratio {est_words / ex['words']:.4f}); effective sample size {s['ess'][0]:.3g}; {time.time() - t0:.1f} s")
        ok &= abs(est_words / ex["words"] - 1) < 0.01
        keys = (["k1", "l2"] if K == 2 else ["sheet0"]) + [(q, p) for p in primes for q in ("L0", "Lnz", "S", "T")]
        worst = 0.0
        for key in keys:
            v, se = s[key]
            zscore = (v - ex[key]) / se if se > 0 else 0.0
            worst = max(worst, abs(zscore))
            if key in ("k1", "l2", "sheet0") or key[1] in (3, 5, 17):
                print(f"   {str(key):>12}: exact {ex[key]:.5f}   MC {v:.5f} +- {se:.5f}   z = {zscore:+.2f}")
        print(f"   worst |z| over {len(keys)} statistics: {worst:.2f}")
        ok &= worst < 4.5
    return ok


def scan(nsamp_m):
    ok = True
    nb = 60
    ns = int(nsamp_m * 1e6 / nb)
    rows = {}
    for K in (1, 2):
        for n0 in (30, 60, 120, 240, 480, 800):
            t0 = time.time()
            fl, cls, ccls, cat = run_batches(K, n0, nb, ns, 777 + 1000 * K + n0, PRIMES)
            rows[(K, n0)] = summarize(K, n0, fl, cls, ccls, cat)
            print(f"# K = {K}, n = {n0}: {nb * ns} paths, effective sample size {rows[(K, n0)]['ess'][0]:.3g}, {time.time() - t0:.0f} s",
                  flush=True)
    ns_list = (30, 60, 120, 240, 480, 800)
    print("\nK = 1 (exactly one J): by-value shares and class factors;  limits (b) by-value: sheet 0 share w_0 = 0.19641, "
          "L(0) = 1 - 1/p^2 + w_1/p, L(-1) = 1 - 1/p^2 + w_0/p;  (a) letter-count: 1/4, 1 + (3p-4)/(4p^2), 1 - 1/p^2 + 1/(4p)")
    print(f"  {'n':>4} {'sheet0':>16} " + " ".join(f"{'L(0) p=' + str(int(p)):>16} {'L(-1) p=' + str(int(p)):>16}" for p in (3, 5, 11)))
    for n0 in ns_list:
        s = rows[(1, n0)]
        print(f"  {n0:4d} {s['sheet0'][0]:8.5f}+-{s['sheet0'][1]:.5f} " + " ".join(
            f"{s[('L0', p)][0]:8.4f}+-{s[('L0', p)][1]:.4f} {s[('Lm1', p)][0]:8.4f}+-{s[('Lm1', p)][1]:.4f}" for p in (3, 5, 11)))
    print(f"  {'(b)':>4} {W0:16.5f} " + " ".join(f"{1 - 1 / p ** 2 + W1 / p:16.4f} {1 - 1 / p ** 2 + W0 / p:16.4f}" for p in (3, 5, 11)))
    print(f"  {'(a)':>4} {0.25:16.5f} " + " ".join(f"{1 + (3 * p - 4) / (4 * p * p):16.4f} {1 - 1 / p ** 2 + 1 / (4 * p):16.4f}" for p in (3, 5, 11)))
    print("\nK = 2 (exactly two J): flags and cascade;  by-value limit: kappa = lambda = w_1 = 0.80359, P(S2) = (p-1)/p^3 + w_1/p^2;"
          "  letter-count: 3/4, 3/4, (7p-4)/(4p^3)")
    print(f"  {'n':>4} {'kappa':>16} {'lambda':>16} {'kappa&lambda':>16} " + " ".join(f"{'P(S2) p=' + str(p):>16}" for p in (3, 5, 7)))
    for n0 in ns_list:
        s = rows[(2, n0)]
        print(f"  {n0:4d} {s['k1'][0]:8.5f}+-{s['k1'][1]:.5f} {s['l2'][0]:8.5f}+-{s['l2'][1]:.5f} {s['k1l2'][0]:8.5f}+-{s['k1l2'][1]:.5f} "
              + " ".join(f"{s[('S', p)][0]:8.5f}+-{s[('S', p)][1]:.5f}" for p in (3, 5, 7)))
    print(f"  {'val':>4} {W1:16.5f} {W1:16.5f} {W1 * W1:16.5f} " + " ".join(f"{(p - 1) / p ** 3 + W1 / p ** 2:16.5f}" for p in (3, 5, 7)))
    print(f"  {'let':>4} {0.75:16.5f} {0.75:16.5f} {0.5625:16.5f} " + " ".join(f"{(7 * p - 4) / (4 * p ** 3):16.5f}" for p in (3, 5, 7)))
    print("\nK = 2 class factors of the exactly-two-J words: L(0) and rho = L(0) / mean L(r != 0);  limits: by-value "
          "L(0) = 1 + (p-1)(2 w_1 p - 1)/p^3, letter-count 1 + (3p-2)(p-1)/(2p^3)")
    for q in ("L0", "rho"):
        print(f"  {q:>4} {'n':>4} " + " ".join(f"{'p=' + str(int(p)):>15}" for p in PRIMES))
        for n0 in ns_list:
            s = rows[(2, n0)]
            print(f"  {'':>4} {n0:4d} " + " ".join(f"{s[(q, int(p))][0]:8.4f}+-{s[(q, int(p))][1]:.4f}" for p in PRIMES))
        if q == "L0":
            print(f"  {'':>4} {'val':>4} " + " ".join(f"{limits_k2(int(p), 2 * W1)[0]:15.4f}" for p in PRIMES))
            print(f"  {'':>4} {'let':>4} " + " ".join(f"{limits_k2(int(p), 1.5)[0]:15.4f}" for p in PRIMES))
        else:
            print(f"  {'':>4} {'val':>4} " + " ".join(f"{limits_k2(int(p), 2 * W1)[0] / limits_k2(int(p), 2 * W1)[1]:15.4f}" for p in PRIMES))
            print(f"  {'':>4} {'let':>4} " + " ".join(f"{limits_k2(int(p), 1.5)[0] / limits_k2(int(p), 1.5)[1]:15.4f}" for p in PRIMES))
    print("\nK = 2 confined share (S2 or T), unconfined class-0 ratio U0/Unz, and P(S2 | C1) (cascade survival, limit kappa):")
    print(f"  {'n':>4} " + " ".join(f"{'conf p=' + str(p):>15} {'U0/Unz':>15} {'S2|C1':>13}" for p in (3, 5)))
    for n0 in ns_list:
        s = rows[(2, n0)]
        print(f"  {n0:4d} " + " ".join(
            f"{s[('conf', p)][0]:8.4f}+-{s[('conf', p)][1]:.4f} {s[('U0/Unz', p)][0]:8.4f}+-{s[('U0/Unz', p)][1]:.4f} "
            f"{s[('S|C1', p)][0]:7.4f}+-{s[('S|C1', p)][1]:.4f}" for p in (3, 5)))
    np.save("mc_byvalue_scan.npy", np.array([[K, n0, float(p), rows[(K, n0)][("L0", int(p))][0], rows[(K, n0)][("L0", int(p))][1],
                                               rows[(K, n0)][("Lnz", int(p))][0], rows[(K, n0)][("Lnz", int(p))][1]]
                                              for (K, n0) in rows for p in PRIMES]))
    # decisions: at the largest n the flags must sit at w_1 (and not at 3/4)
    s = rows[(2, 800)]
    for key in ("k1", "l2"):
        v, se = s[key]
        ok &= abs(v - W1) < 0.004 + 4 * se
    s = rows[(1, 800)]
    ok &= abs(s["sheet0"][0] - W0) < 0.003 + 4 * s["sheet0"][1]
    return ok


def main():
    import numba
    numba.set_num_threads(3)
    mode = sys.argv[1] if len(sys.argv) > 1 else "validate"
    if mode == "validate":
        ok = validate()
    else:
        ok = scan(float(sys.argv[2]) if len(sys.argv) > 2 else 6.0)
    print("ALL CHECKS PASSED" if ok else "CHECK FAILED")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
