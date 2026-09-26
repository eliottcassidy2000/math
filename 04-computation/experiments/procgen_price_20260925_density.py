#!/usr/bin/env python3
"""procgen_price_20260925 -- part 1: the exact Collatz-undecided densities rho_L (DP), and the exponent 1 - h.

Bad_L = classes n mod 2^L whose Collatz (shortcut T) parity word w_1..w_L has 3^{a_k} > 2^k for every k <= L
(a_k = number of odd steps among the first k).  Every n in such a class satisfies T^k(n) > n for k <= L; every n
in any other class with n above an explicit threshold descends within L steps (Terras).  rho_L = |Bad_L| / 2^L.

Exact DP over (k, a_k).  Outputs: |Bad_L| for L <= 40 (and the cross-check value 27,328 at L = 20), the ratios
rho_{L+1}/rho_L, the fit rho_L ~ C L^{-3/2} 2^{-(1-h)L} with h = h(log_3 2), and the constant C along L <= 4000.
Also the exact expected "compression" E[1/Rmax] over Bad_L (Rmax = the largest height 3^{a_k}/2^k on the
first L steps), the factor by which a flip placed at the highest point of each bad path is cheaper.
Runtime < 1 minute; memory < 50 MB.
"""
import math
from fractions import Fraction

LOG32 = math.log(2) / math.log(3)          # p = log_3 2
h = -(LOG32 * math.log2(LOG32) + (1 - LOG32) * math.log2(1 - LOG32))
GAP = 1 - h

def thr(k):
    """least a with 3^a > 2^k"""
    a = 0
    while 3 ** a <= 2 ** k:
        a += 1
    return a

def bad_counts(Lmax):
    """cnt[L] = |Bad_L| (exact integers), via the prefix condition a_k >= thr(k) for all k <= L."""
    T = [thr(k) for k in range(Lmax + 1)]
    dist = {0: 1}           # a -> number of words of length k that are prefix-bad so far
    out = [1]
    for k in range(1, Lmax + 1):
        nd = {}
        for a, c in dist.items():
            for up in (0, 1):
                b = a + up
                if b >= T[k]:
                    nd[b] = nd.get(b, 0) + c
        dist = nd
        out.append(sum(dist.values()))
    return out

def bad_counts_float(Lmax):
    """same, in floating point normalised by 2^L (for large L)."""
    T = None
    dist = {0: 1.0}
    out = [1.0]
    a_thr = 0
    for k in range(1, Lmax + 1):
        # threshold: least a with a*log2(3) > k
        a_thr = math.floor(k * LOG32) + 1
        nd = {}
        for a, c in dist.items():
            for up in (0, 1):
                b = a + up
                if b >= a_thr:
                    nd[b] = nd.get(b, 0.0) + 0.5 * c
        dist = nd
        out.append(sum(dist.values()))
    return out

def compression(L):
    """E[1/Rmax] and E[log2 Rmax] over Bad_L (uniform over bad words); Rmax = max_k 3^{a_k}/2^k (k <= L)."""
    # state: (a, maxexp) where maxexp = argmax pair stored as (a_j, j) maximising a_j log 3 - j log 2
    # store the max as a float key of log2 height rounded to 1e-9 together with the exact pair
    dist = {(0, (0, 0)): 1}
    T = [thr(k) for k in range(L + 1)]
    for k in range(1, L + 1):
        nd = {}
        for (a, m), c in dist.items():
            for up in (0, 1):
                b = a + up
                if b < T[k]:
                    continue
                hb = b * math.log2(3) - k
                hm = m[0] * math.log2(3) - m[1]
                mm = (b, k) if hb > hm else m
                key = (b, mm)
                nd[key] = nd.get(key, 0) + c
        dist = nd
    tot = sum(dist.values())
    inv = sum(c * 2.0 ** (-(m[0] * math.log2(3) - m[1])) for (a, m), c in dist.items()) / tot
    lg = sum(c * (m[0] * math.log2(3) - m[1]) for (a, m), c in dist.items()) / tot
    return tot, inv, lg

if __name__ == "__main__":
    print("=" * 100)
    print("PART 1. Exact Collatz-undecided densities rho_L = |Bad_L| / 2^L  (FINITE-EXACT, DP over (k, a_k))")
    print("=" * 100)
    print(f"  p = log_3 2 = {LOG32:.10f};  h = h(p) = {h:.10f} bits;  1 - h = {GAP:.10f}")
    cnt = bad_counts(64)
    print(f"  cross-check: |Bad_20| = {cnt[20]} (brackets lane / foundry value 27,328);  |Bad_4| = {cnt[4]} (n = 7,11,15 mod 16)")
    print("   L    |Bad_L|            rho_L          rho_L * L^1.5 * 2^((1-h)L)")
    for L in range(1, 41):
        r = cnt[L] / 2 ** L
        print(f"  {L:3d}  {cnt[L]:>14d}   {r:.10f}   {r * L ** 1.5 * 2 ** (GAP * L):.5f}")
    fl = bad_counts_float(4000)
    print("  large L (floating DP):   L   rho_L            -log2(rho_L)/L    C_L = rho_L L^1.5 2^((1-h)L)")
    for L in (50, 100, 200, 400, 800, 1600, 3200, 4000):
        r = fl[L]
        print(f"                        {L:5d}   {r:.6e}   {-math.log2(r) / L:.6f}         {r * L ** 1.5 * 2 ** (GAP * L):.5f}")
    print("  => rho_L = 2^{-(1-h)L + O(log L)}; the prefactor C_L stays bounded (oscillating with the Beatty")
    print("     pattern of log_3 2), consistent with rho_L = Theta(L^{-3/2} 2^{-(1-h)L}).")
    print()
    print("  Compression: a flip at the highest point of a bad path serves a progression of density 2^-L / Rmax.")
    print("   L   |Bad_L|    E[1/Rmax]    E[log2 Rmax]   rho_L * E[1/Rmax]")
    for L in (4, 5, 8, 10, 12, 16, 20, 24, 28, 32, 36, 40):
        tot, inv, lg = compression(L)
        assert tot == cnt[L]
        print(f"  {L:3d}  {tot:9d}   {inv:.5f}      {lg:.4f}         {cnt[L] / 2 ** L * inv:.6f}")
