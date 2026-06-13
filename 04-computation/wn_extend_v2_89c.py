#!/usr/bin/env python3
"""
wn_extend_v2_89c.py — Extend W(n) to n=22+ using proven DP
opus-2026-03-15-S89c

Uses the exact same DP logic as gk_verify_89c.py but optimized with
arrays instead of dict, and processes masks in order.
"""

import sys
import time
from fractions import Fraction
from math import factorial
from functools import reduce

def falling_factorial(n, k):
    return reduce(lambda a, b: a * b, range(n, n - k, -1), 1)

def compute_W_dp(n, verbose=True):
    """Compute W(n) using bitmask DP with dict. Exact copy of verified logic."""
    if verbose:
        print(f"Computing W({n})...", flush=True)
    t0 = time.time()

    full = (1 << n) - 1
    dp = {}
    for v in range(n):
        dp[((1 << v), v)] = 1

    count = 0
    for mask in range(1, full + 1):
        if verbose and mask % (1 << max(0, n-6)) == 0:
            elapsed = time.time() - t0
            pct = mask / full * 100
            print(f"  mask progress: {pct:.1f}% ({elapsed:.1f}s)", flush=True)

        for v in range(n):
            if not (mask & (1 << v)):
                continue
            cnt = dp.get((mask, v), 0)
            if cnt == 0:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if u == v - 1:
                    continue
                weight = 2 * cnt if u == v + 1 else cnt
                new_key = (mask | (1 << u), u)
                dp[new_key] = dp.get(new_key, 0) + weight

    W = sum(dp.get((full, v), 0) for v in range(n))
    elapsed = time.time() - t0
    if verbose:
        print(f"  W({n}) = {W} ({elapsed:.1f}s)")
    return W

# Known W values for verification
known_W = {
    3: 8, 4: 36, 5: 176, 6: 1080, 7: 8268, 8: 76104, 9: 822880,
    10: 10262160, 11: 145655148, 12: 2326258872, 13: 41383027792,
    14: 814116023280, 15: 17574474321708, 16: 413321769498648,
    17: 10545197020978752, 18: 290565166853498160,
}

# g_k polynomial evaluation
gk_coeffs = {
    1: (0, 0, 3, 0),
    2: (0, 3, 0, 0),
    3: (2, 0, 1, 0),
    4: (10, -33, 50, -24),
    5: (388, -2040, 3431, -1776),
    6: (69660, -380445, 653748, -342960),
    7: (19826270, -109486152, 189674605, -100014720),
    8: (7309726742, -40641958545, 70757788486, -37425556680),
    9: (3262687720240, -18232387983408, 31858349908595, -16888649645424),
}

def g_eval(k, m):
    if k in gk_coeffs:
        a, b, c, d = gk_coeffs[k]
        return Fraction(a*m**3 + b*m**2 + c*m + d, 3)
    if m == 1:
        return Fraction(1)
    if m == 2:
        return Fraction(2*k)
    return None

def extract_new_gk(n, W_n):
    """Given W(n), extract the residual after subtracting known g_k contributions.
    Returns (k_new, m_new, g_new_value) where g_{k_new}(m_new) is the extracted value.
    """
    nf = factorial(n)
    cv2 = Fraction(W_n, nf) - 1

    # Subtract all known contributions
    residual = cv2
    max_k = (n - 1) // 2
    unknown_k = None

    for k in range(1, max_k + 1):
        m = n - 2*k
        if m < 1:
            break
        g = g_eval(k, m)
        if g is not None:
            residual -= Fraction(2 * g, falling_factorial(n, 2*k))
        else:
            if unknown_k is None:
                unknown_k = k
            # Can still continue if higher k values are known via boundary conditions

    if unknown_k is not None:
        # The residual includes the unknown term: 2·g_{unknown_k}(m)/(n)_{2·unknown_k}
        m = n - 2*unknown_k
        ff = falling_factorial(n, 2*unknown_k)
        g_val = residual * ff / 2
        return unknown_k, m, g_val, residual
    else:
        return None, None, Fraction(0), residual

# First verify at n=14
print("Verification at n=14:")
W14 = compute_W_dp(14)
print(f"Known W(14) = {known_W[14]}, computed = {W14}, match = {W14 == known_W[14]}")

# Now benchmark n=16 to estimate feasibility of larger n
print(f"\nBenchmark n=16:")
t0 = time.time()
W16 = compute_W_dp(16)
t16 = time.time() - t0
print(f"Known W(16) = {known_W[16]}, match = {W16 == known_W[16]}")
print(f"Time: {t16:.1f}s")

# Estimate for n=22: time roughly scales as n² × 2^n / (16² × 2^16)
scale = (22**2 * 2**22) / (16**2 * 2**16)
print(f"\nEstimated time for n=22: {t16 * scale:.0f}s = {t16 * scale / 3600:.1f}h")
print("This is infeasible in Python. Need C implementation.")

# But we CAN try n=17, 18 and see if the timings match known values
if t16 < 60:  # only if n=16 was fast enough
    print(f"\nTrying n=17:")
    W17 = compute_W_dp(17)
    print(f"Known W(17) = {known_W[17]}, match = {W17 == known_W[17]}")

print("\nDone!")
