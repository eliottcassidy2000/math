#!/usr/bin/env python3
"""
gk_truth_v2_s112.py — Definitive test: transfer matrix g_k vs opus g_k
kind-pasteur-2026-03-15-S112

W(n) = sum_{sigma in S_n: NUD} 2^{adj1(sigma)}
CV^2 = W(n)/n! - 1

Compute W(n) via bitmask DP (feasible to n~20).
Then test both g_k families against CV^2.
"""

from fractions import Fraction
from math import factorial

def compute_W_dp(n):
    """W(n) via bitmask DP. O(2^n * n)."""
    if n == 1:
        return 1
    # dp[mask][last] = sum of weights for partial permutations using 'mask' values
    # ending with value 'last'
    dp = {}
    # Initialize: place first value
    for v in range(n):
        dp[(1 << v, v)] = 1

    for length in range(1, n):
        new_dp = {}
        for (mask, last), weight in dp.items():
            if bin(mask).count('1') != length:
                continue
            for v in range(n):
                if mask & (1 << v):
                    continue  # already used
                if v == last - 1:
                    continue  # unit descent, NUD violated
                new_weight = weight * (2 if v == last + 1 else 1)
                key = (mask | (1 << v), v)
                new_dp[key] = new_dp.get(key, 0) + new_weight
        for k, w in new_dp.items():
            dp[k] = dp.get(k, 0) + w

    full_mask = (1 << n) - 1
    return sum(dp.get((full_mask, v), 0) for v in range(n))

def falling_factorial(n, k):
    r = Fraction(1)
    for i in range(k):
        r *= (n - i)
    return r

def transfer_gk(n):
    """Compute g_k for all valid k via transfer matrix."""
    num_edges = n - 2
    k_max = (n - 1) // 2

    if k_max <= 0:
        return {}

    state = [[Fraction(1)] + [Fraction(0)] * k_max,
             [Fraction(0)] * (k_max + 1),
             [Fraction(0)] * (k_max + 1)]

    for step in range(num_edges):
        A, B, C = state
        new_A = [Fraction(0)] * (k_max + 1)
        new_B = [Fraction(0)] * (k_max + 1)
        new_C = [Fraction(0)] * (k_max + 1)

        for i in range(k_max + 1):
            new_A[i] += A[i] + C[i]
            if i > 0:
                new_B[i] += 2 * A[i-1] + C[i-1]
            new_C[i] += B[i]

        state = [new_A, new_B, new_C]

    total = [state[0][i] + state[1][i] + state[2][i] for i in range(k_max + 1)]

    result = {}
    for k in range(1, k_max + 1):
        if total[k] != 0:
            result[k] = total[k] / 2
    return result

def opus_gk(k, m):
    """Opus degree-3 g_k. Returns None if k > 9."""
    coeffs = {
        1: None,  # special
        2: None,  # special
        3: (2, 0, 1, 0),
        4: (10, -33, 50, -24),
        5: (388, -2040, 3431, -1776),
        6: (69660, -380445, 653748, -342960),
        7: (19826270, -109486152, 189674605, -100014720),
        8: (7309726742, -40641958545, 70757788486, -37425556680),
        9: (3262687720240, -18232387983408, 31858349908595, -16888649645424),
    }
    if k == 1:
        return Fraction(m)
    if k == 2:
        return Fraction(m * m)
    if k not in coeffs or coeffs[k] is None:
        return None
    a, b, c, d = coeffs[k]
    return Fraction(a * m**3 + b * m**2 + c * m + d, 3)

print("="*70)
print("W(n) and CV^2 via bitmask DP")
print("="*70)

W_vals = {}
for n in range(1, 19):
    W = compute_W_dp(n)
    W_vals[n] = W
    nfact = factorial(n)
    cv2 = Fraction(W, nfact) - 1
    print(f"n={n:2d}: W={W}, CV^2 = {cv2} = {float(cv2):.10f}")

print("\n" + "="*70)
print("COMPARISON: Transfer Matrix vs Opus")
print("Convention: sum over k with m = n-2k >= 1")
print("="*70)

for n in range(3, 19):
    nfact = factorial(n)
    cv2_true = Fraction(W_vals[n], nfact) - 1

    # Transfer matrix
    tm = transfer_gk(n)
    cv2_tm = Fraction(0)
    for k, gk in tm.items():
        m = n - 2*k
        if m < 1:
            continue
        cv2_tm += 2 * gk / falling_factorial(n, 2*k)

    # Opus
    cv2_opus = Fraction(0)
    opus_ok = True
    for k in range(1, (n-1)//2 + 1):
        m = n - 2*k
        if m < 1:
            continue
        gk = opus_gk(k, m)
        if gk is None:
            opus_ok = False
            break
        cv2_opus += 2 * gk / falling_factorial(n, 2*k)

    tm_match = cv2_tm == cv2_true
    opus_match = opus_ok and cv2_opus == cv2_true

    tm_err = float(cv2_tm - cv2_true) if not tm_match else 0
    opus_err = float(cv2_opus - cv2_true) if opus_ok and not opus_match else 0

    tm_s = "OK" if tm_match else f"err={tm_err:.2e}"
    opus_s = "OK" if opus_match else (f"err={opus_err:.2e}" if opus_ok else "incomplete")

    print(f"n={n:2d}: TM={tm_s:15s}  Opus={opus_s:15s}")

    # Show discrepant terms
    if not tm_match or (opus_ok and not opus_match):
        for k in range(1, (n-1)//2 + 1):
            m = n - 2*k
            if m < 1:
                continue
            tm_gk = tm.get(k, Fraction(0))
            op_gk = opus_gk(k, m)
            if op_gk is not None and tm_gk != op_gk:
                print(f"      k={k}, m={m}: TM g_{k}={tm_gk}, opus g_{k}={op_gk}, diff={tm_gk-op_gk}")

print("\nDone!")
