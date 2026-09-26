#!/usr/bin/env python3
"""
procgen_cubedist_20260925_bad.py -- the undecided-residue construction (HYP-9138), cube-distance lane 2026-09-25.

sigma_k := Collatz (all +) with sigma = - exactly on Bad_k, the residues r mod 2^k on which Collatz has no
descent within k steps (3^(a_j) > 2^j for j = 1..k).

Theorem 1 (proved in the note, checked here): sigma_k is in class (i); every cycle of G_(sigma_k) has odd
density <= F_k, the best lower rational approximation of log_3 2 with denominator <= k, with equality;
every integer n > n_0(k) descends below itself within k steps.  Hence delta_k <= |Bad_k| and
delta_k / 2^(k-1) <= 2 |Bad_k| / 2^k -> 0.
Theorem 2 (lower bound): delta_k >= N_k, the number of binary necklaces of length k with more than k log_3 2 ones.

Sections (every claim is a check() that raises on failure):
  B1  |Bad_k| by enumeration (numpy) = ballot DP; the block lemma, checked on the parity graph for k <= 13;
  B2  class (i) of sigma_k: exact Karp (k <= 15) and an integer potential certificate at threshold F_k,
      verified edge by edge in Python (k <= 20); rho_max = F_k (Karp for k <= 15; for every k the explicit
      first-descent word of the note survives as a cycle);
  B3  bounded lookahead: descent within k steps for n > n_0, and the finite check that every n <= n_0 reaches 1
      (so sigma_k is transitive on the positive integers) for k <= 20;
  B4  the sandwich N_k <= delta_k <= |Bad_k| and its asymptotics (exact integers, k <= 200).
"""
import sys
import os
import math
import time
from fractions import Fraction
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import procgen_cubedist_20260925_lib as L  # noqa: E402

LOG3_2 = math.log(2) / math.log(3)


def crit_table(K):
    """crit[j] = least a with 3^a > 2^j"""
    return [L.critical_ones(j) for j in range(K + 1)]


def bad_numpy(k):
    """Bad_k by vectorised simulation of Collatz on the representatives 0..2^k-1 (exact int64 for k <= 30)."""
    L.check(k <= 30, "int64 range")
    crit = crit_table(k)
    n = np.arange(1 << k, dtype=np.int64)
    a = np.zeros(1 << k, dtype=np.int64)
    alive = np.ones(1 << k, dtype=bool)
    first = np.zeros(1 << k, dtype=np.int64)       # first descent time (0 = none within k)
    for j in range(1, k + 1):
        odd = (n & 1) == 1
        a += odd
        n = np.where(odd, (3 * n + 1) >> 1, n >> 1)
        desc = a < crit[j]              # 3^a < 2^j  (3^a != 2^j)
        newly = alive & desc
        first[newly] = j
        alive &= ~desc
    return np.nonzero(alive)[0], first


def first_descent_word(F):
    """the word a_j = floor(c j) + 1 (j < d), a_d = a, for F = a/d the best lower approximation:
    all proper prefixes ballot, total contracting."""
    a, d = F.numerator, F.denominator
    word = []
    prev = 0
    for j in range(1, d + 1):
        aj = (math.floor(LOG3_2 * j) + 1) if j < d else a
        word.append(aj - prev)
        prev = aj
    return word


def block_lemma_check(k, badset):
    """For every residue s outside Bad_k, and 1 <= i < d(s): no 2-adic x = s (mod 2^k) has T^i(x) in Bad_k.
    T^i(x) mod 2^(k-i) is determined by s, and T^i(x) mod 2^k ranges over all lifts; so the check is that the
    class T^i(s) mod 2^(k-i) contains no element of Bad_k."""
    K = 1 << k
    badclasses = [None] * (k + 1)
    arr = np.array(sorted(badset), dtype=np.int64)
    for m in range(0, k + 1):
        badclasses[m] = set((arr % (1 << m)).tolist())
    viol = 0
    checked = 0
    for s in range(K):
        if s in badset:
            continue
        d = L.first_descent(s, k)
        L.check(d is not None and d <= k, "a residue outside Bad_k has no descent within k")
        n = s
        for i in range(1, d):
            n = L.step(n, 3, 1)
            m = k - i
            if (n % (1 << m)) in badclasses[m]:
                viol += 1
            checked += 1
    return viol, checked


def descent_threshold(k, badset):
    """n_0 = max over odd s outside Bad_k of c/(2^d - 3^a), where T^d(n) = (3^a n + c)/2^d on n = s mod 2^k."""
    n0 = 0
    for s in range(1, 1 << k, 2):
        if s in badset:
            continue
        a, c, n = 0, 0, s
        for j in range(1, k + 1):
            if n & 1:
                a += 1
                c = 3 * c + (1 << (j - 1))
                n = (3 * n + 1) // 2
            else:
                n //= 2
            if 3 ** a < 2 ** j:
                thr = c // (2 ** j - 3 ** a)
                n0 = max(n0, thr)
                break
    return n0


def n0_dp(K):
    """n_0(k) for all k <= K without enumerating residues: n_0(k) = max over first-descent words w of length
    d <= k (3^(a_j) > 2^j for j < d, 3^(a_d) < 2^d) of floor(c_w / (2^d - 3^(a_d))), where T^d(n) = (3^a n + c_w)/2^d.
    c_(j+1) = 3 c_j + 2^j after an odd step and c_j after an even one is increasing in c_j, so a DP over (j, a_j)
    keeping the maximal c is exact.  Returns {k: n_0(k)}."""
    cur = {0: 0}
    best = 0
    out = {}
    for j in range(0, K):
        nxt = {}
        for a, c in cur.items():
            for b in (0, 1):
                na = a + b
                nc = 3 * c + (1 << j) if b else c
                J = j + 1
                if 3 ** na > 2 ** J:
                    if nxt.get(na, -1) < nc:
                        nxt[na] = nc
                else:
                    best = max(best, nc // (2 ** J - 3 ** na))
        cur = nxt
        out[j + 1] = best
    return out


def in_bad(m, k):
    """is m mod 2^k in Bad_k (Collatz ballot for k steps)?"""
    return L.ballot(L.parity_word(m % (1 << k), k), 3)


def transitive_check_large(k, n0):
    """every n in [1, n0] reaches 1 under T_sigma_k, membership in Bad_k tested on demand"""
    reach = {1, 2}
    memo = {}
    for n in range(1, n0 + 1):
        path = []
        m = n
        onpath = set()
        while m not in reach:
            if m in onpath:
                return False, n
            onpath.add(m)
            path.append(m)
            if m & 1:
                r = m % (1 << k)
                if r not in memo:
                    memo[r] = in_bad(r, k)
                m = (3 * m + (-1 if memo[r] else 1)) // 2
            else:
                m //= 2
            if len(path) > 100000:
                return False, n
        reach.update(path)
    return True, None


def transitive_check(k, badset, n0):
    """every n in [1, n0] reaches 1 under T_sigma (sigma = - on Bad_k)."""
    K = 1 << k
    reach = {1: True, 2: True}
    for n in range(1, n0 + 1):
        path = []
        m = n
        steps = 0
        while m not in reach:
            path.append(m)
            if m & 1:
                m = (3 * m + (-1 if (m % K) in badset else 1)) // 2
            else:
                m //= 2
            steps += 1
            if steps > 100000 or m in path[-1000:]:
                return False, n
        val = reach[m]
        for p in path:
            reach[p] = val
        if not val:
            return False, n
    return True, None


def main(log=print, kmax_cert=20, kmax_karp=15, kmax_lemma=13, kmax_trans=300):
    log("B1. Bad_k: enumeration vs ballot DP; the block lemma")
    for k in range(2, kmax_cert + 1):
        bad, first = bad_numpy(k)
        L.check(len(bad) == L.ballot_count(k), f"|Bad_{k}| enumeration != DP")
        L.check(all(int(r) % 4 == 3 for r in bad), "Bad_k not inside 3 mod 4")
    log(f"  |Bad_k| enumeration = ballot DP for k = 2..{kmax_cert}; every element of Bad_k is 3 (mod 4) (k >= 2)")
    for k in range(2, kmax_lemma + 1):
        bad, _ = bad_numpy(k)
        viol, checked = block_lemma_check(k, set(int(r) for r in bad))
        L.check(viol == 0, f"block lemma fails at k={k}")
    log(f"  block lemma (no Bad_k point strictly inside a Collatz descent block), all lifts: k = 2..{kmax_lemma}: 0 violations")

    log("")
    log("B2. class (i) of sigma_k = Collatz with sigma = - on Bad_k")
    log("   k  |Bad_k|  Haar=|Bad|/2^(k-1)  F_k (bound)  Karp rho_max   potential cert (max psi)  first-descent word of F_k")
    rows = []
    for k in range(2, kmax_cert + 1):
        t0 = time.time()
        bad, _ = bad_numpy(k)
        badset = set(int(r) for r in bad)
        F = L.best_lower_approx(k)
        mask = L.mask_of(badset)
        kd = None
        if k <= kmax_karp:
            kd = L.karp_density(k, mask)
            L.check(kd == F, f"Karp rho_max != F_k at k={k}")
            L.check(3 ** kd.numerator < 2 ** kd.denominator, "not class (i)")
        st, psi = L.certificate(k, mask, F.numerator, F.denominator)
        L.check(st == 'OK', f"no potential certificate at k={k}")
        L.check(L.verify_certificate(k, badset, psi, F.numerator, F.denominator), f"certificate fails at k={k}")
        L.check(3 ** F.numerator < 2 ** F.denominator, "F_k not below log_3 2")
        w = first_descent_word(F)
        L.check(sum(w) == F.numerator and len(w) == F.denominator, "word size")
        a = 0
        for j, b in enumerate(w, 1):
            a += b
            if j < len(w):
                L.check(3 ** a > 2 ** j, "word prefix not ballot")
        L.check(3 ** a < 2 ** len(w), "word not contracting")
        rows.append((k, len(bad), F, kd, max(psi)))
        log(f"  {k:2d} {len(bad):7d}  {len(bad) / 2 ** (k - 1):.5f}            {str(F):>6}      "
            f"{str(kd) if kd is not None else '-':>6}          OK ({max(psi):6d})              "
            f"{''.join(map(str, w))}   [{time.time() - t0:.1f}s]")
    log("  every certificate: integer psi with psi(t) <= psi(s) - w(s) on all 2^(k+1) edges, w = (r-q, -q) for F_k = q/r;")
    log("  so every cycle has a/p <= F_k < log_3 2 (class (i)).  The first-descent word of F_k is a Collatz cycle of")
    log("  length <= k avoiding Bad_k, so rho_max = F_k exactly (Karp agrees for k <= %d)." % kmax_karp)

    log("")
    log("B3. bounded lookahead and transitivity: every n > n_0 descends within k steps; every n <= n_0 reaches 1")
    log("   k   n_0    all n <= n_0 reach 1")
    for k in range(2, kmax_cert + 1):
        bad, _ = bad_numpy(k)
        badset = set(int(r) for r in bad)
        n0 = descent_threshold(k, badset)
        ok, bad_n = transitive_check(k, badset, max(n0, 2))
        L.check(ok, f"sigma_{k} has an orbit not reaching 1 (n = {bad_n})")
        # direct spot check of the descent claim above n_0
        K = 1 << k
        for n in range(n0 + 1, n0 + 20001):
            m, j = n, 0
            while j < max(k, 2):
                m = m // 2 if m % 2 == 0 else (3 * m + (-1 if (m % K) in badset else 1)) // 2
                j += 1
                if m < n:
                    break
            L.check(m < n, f"n={n} does not descend within k steps at k={k}")
        log(f"  {k:2d} {n0:6d}   yes (and every n in (n_0, n_0 + 20000] descends within {max(k, 2)} steps)")

    log("  n_0(k) from the first-descent-word DP equals the enumerated value for k <= %d" % kmax_cert)
    dp = n0_dp(kmax_trans)
    for k in range(2, kmax_cert + 1):
        bad, _ = bad_numpy(k)
        L.check(dp[k] == descent_threshold(k, set(int(r) for r in bad)), f"n_0 DP mismatch at k={k}")
    jumps = []
    prev = None
    for k in range(2, kmax_trans + 1):
        ok, bad_n = transitive_check_large(k, max(dp[k], 2))
        L.check(ok, f"sigma_{k} has a positive orbit not reaching 1 (n = {bad_n})")
        if dp[k] != prev:
            jumps.append((k, dp[k]))
            prev = dp[k]
    log(f"  k = 2..{kmax_trans}: every n <= n_0(k) reaches 1 under T_sigma_k, so (with the descent above n_0) every")
    log(f"  positive orbit of T_sigma_k reaches 1.  n_0(k) changes only at (k, n_0): {jumps}")
    log("  (the jumps sit at the denominators 2, 5, 8, 27, 46, 65, ... of the best lower approximations F_k).")

    log("")
    log("B4. the sandwich N_k <= delta_k <= |Bad_k| (exact integers); h = h(log_3 2)")
    h = -(LOG3_2 * math.log2(LOG3_2) + (1 - LOG3_2) * math.log2(1 - LOG3_2))
    log(f"  h = {h:.10f}, 1 - h = {1 - h:.10f}")
    log("    k       N_k (lower)        |Bad_k| (upper)     ratio   Haar upper 2|Bad_k|/2^k   N_k k^1.5/2^(hk)  |Bad_k| k^1.5/2^(hk)")
    for k in list(range(2, 41)) + [50, 60, 80, 100, 150, 200]:
        B = L.ballot_count(k)
        N = L.necklace_lower_bound(k)
        L.check(N <= B, "N_k > |Bad_k|")
        L.check(B <= 2 ** (h * k) + 1e-9, "Chernoff bound fails")
        L.check(N * 3 * k * k >= 2 ** (h * k), "N_k >= 2^(hk)/(3k^2) fails")
        log(f"  {k:3d} {N:>22d} {B:>22d}  {B / N:6.3f}   {2 * B / 2 ** k:.3e}               "
            f"{N * k ** 1.5 / 2 ** (h * k):6.3f}           {B * k ** 1.5 / 2 ** (h * k):6.3f}")
    log("  checks: N_k <= |Bad_k| <= 2^(hk) and N_k >= 2^(hk)/(3k^2) at every listed k")
    ks = list(range(8, 201))
    ys = np.array([math.log2(2 * L.ballot_count(k) / 2 ** k) for k in ks])
    X = np.array([[1.0, k, math.log2(k)] for k in ks])
    coef, *_ = np.linalg.lstsq(X, ys, rcond=None)
    X2 = np.array([[1.0, k] for k in ks if k <= 14])
    coef2, *_ = np.linalg.lstsq(X2, np.array([math.log2(2 * L.ballot_count(k) / 2 ** k) for k in ks if k <= 14]),
                                rcond=None)
    log(f"  fit of the Haar density 2|Bad_k|/2^k of sigma_k over 8 <= k <= 200: log2 = {coef[0]:.3f} {coef[1]:+.4f} k "
        f"{coef[2]:+.3f} log2 k   (theory: -(1-h) = {-(1 - h):.4f} per k; -3/2 up to a slowly varying factor)")
    log(f"  a pure exponential fit over 8 <= k <= 14 gives slope {coef2[1]:+.4f} per k (~2^(-k/{-1 / coef2[1]:.1f})): the")
    log("  small-k decay is dominated by the k^(-3/2) factor")
    return rows


if __name__ == '__main__':
    main()
