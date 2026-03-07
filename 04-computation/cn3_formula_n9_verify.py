#!/usr/bin/env python3
"""
Verify tr(c_{n-3}) = 2*(n-2)!*t_3 + const(n) at n=9.

At n=9: c_{n-3} = c_6. Coefficient of t_3 should be 2*7! = 10080.
Use circulant tournaments (scalar M) for efficiency.

Also check: does the constant follow a pattern?
  n=3: const = -0.5
  n=5: const = -30
  n=7: const = -2100

opus-2026-03-06-S26
"""

from itertools import permutations, combinations
import numpy as np
from math import factorial

def count_paths_weighted(A, verts, r_val, start=None, end=None):
    total = 0.0
    for p in permutations(verts):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        w = 1.0
        for i in range(len(p)-1):
            w *= r_val + (A[p[i]][p[i+1]] - 0.5)
        total += w
    return total

def transfer_entry_r(A, a, r_val):
    """Compute M[a,a](r) — single diagonal entry."""
    n = len(A)
    U = [v for v in range(n) if v != a]
    total = 0.0
    for k in range(len(U)+1):
        for S in combinations(U, k):
            S_set = set(S)
            R = [v for v in U if v not in S_set]
            S_verts = sorted(list(S) + [a])
            R_verts = sorted(R + [a])
            ea = count_paths_weighted(A, S_verts, r_val, end=a)
            ba = count_paths_weighted(A, R_verts, r_val, start=a)
            total += ((-1)**k) * ea * ba
    return total

def ham_path_count_dp(A):
    """Count Ham paths using DP on bitmasks."""
    n = len(A)
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for size in range(2, n+1):
        for S in range(1 << n):
            if bin(S).count('1') != size:
                continue
            for v in range(n):
                if not (S & (1 << v)):
                    continue
                S_prev = S ^ (1 << v)
                total = 0
                for u in range(n):
                    if not (S_prev & (1 << u)):
                        continue
                    if A[u][v] and (S_prev, u) in dp:
                        total += dp[(S_prev, u)]
                if total > 0:
                    dp[(S, v)] = total
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def count_3cycles(A):
    n = len(A)
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                count += (A[i][j]*A[j][k]*A[k][i] + A[i][k]*A[k][j]*A[j][i])
    return count

def circulant_tournament(n, gen_set):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j and (j-i)%n in gen_set:
                A[i][j] = 1
    return A

# =====================================================================
n = 9
print("=" * 70)
print(f"n={n}: VERIFY tr(c_{{n-3}}) = 2*(n-2)!*t_3 + const")
print("=" * 70)

# For circulant: M = (H/n)*I, so tr(M(r)) = n * M[0,0](r)
# M[0,0](r) = sum of c_k/n * r^2k (scalar coefficients)
# Need (n-1)/2 + 1 = 5 sample points for degree-4 polynomial in u = r^2

num_coeffs = (n-1)//2 + 1  # 5 coefficients: c_0, c_2, c_4, c_6, c_8
u_samples = np.linspace(0, 0.5, num_coeffs + 1)  # Extra point for robustness
r_samples = np.sqrt(u_samples)

# Generate circulant tournaments
half = list(range(1, (n+1)//2))
gen_sets = set()
for mask in range(1 << len(half)):
    gs = set()
    for k, d in enumerate(half):
        if mask & (1 << k):
            gs.add(d)
        else:
            gs.add(n - d)
    gen_sets.add(frozenset(gs))

print(f"\n  {len(gen_sets)} circulant tournaments on n={n}")
print(f"  Expected coefficient of t_3: 2*{n-2}! = 2*{factorial(n-2)} = {2*factorial(n-2)}")
print(f"  Expected c_{{n-1}} = tr(c_8) = {n}*{factorial(n-1)} = {n*factorial(n-1)}")

data = []
for gs in sorted(gen_sets)[:4]:  # Just check first 4 (for speed)
    A = circulant_tournament(n, gs)
    H = ham_path_count_dp(A)
    t3 = count_3cycles(A)

    # For scalar M: M[0,0](r) = tr(M(r))/n
    # Compute M[0,0] at sample r values
    print(f"\n  gen={sorted(gs)}: H={H}, t3={t3}", end="", flush=True)

    vals = [transfer_entry_r(A, 0, rv) for rv in r_samples]

    # Fit polynomial in u = r^2 of degree (n-1)/2 = 4
    poly = np.polyfit(u_samples, vals, num_coeffs - 1)
    # poly[0] = coefficient of u^4 = c_8/n
    # poly[1] = coefficient of u^3 = c_6/n
    # poly[2] = coefficient of u^2 = c_4/n
    # poly[3] = coefficient of u^1 = c_2/n
    # poly[4] = coefficient of u^0 = c_0/n

    tr_c8 = n * poly[0]
    tr_c6 = n * poly[1]
    tr_c4 = n * poly[2]
    tr_c2 = n * poly[3]
    tr_c0 = n * poly[4]

    # Verify reconstruction
    H_recon = tr_c0 + tr_c2/4 + tr_c4/16 + tr_c6/64 + tr_c8/256
    print(f" H_recon={H_recon:.2f}")

    data.append((sorted(gs), H, t3, tr_c0, tr_c2, tr_c4, tr_c6, tr_c8))

    print(f"    tr(c_8)={tr_c8:.1f} (expected {n*factorial(n-1)})")
    print(f"    tr(c_6)={tr_c6:.1f}")
    print(f"    tr(c_4)={tr_c4:.1f}")
    print(f"    tr(c_2)={tr_c2:.1f}")
    print(f"    tr(c_0)={tr_c0:.4f}")

# Check formula: tr(c_6) = 2*7!*t3 + const
if len(data) >= 2:
    # Use first two data points to solve 2*7!*t3 + const = tr_c6
    t3_0, c6_0 = data[0][2], data[0][5]
    t3_1, c6_1 = data[1][2], data[1][5]

    if t3_0 != t3_1:
        slope = (c6_1 - c6_0) / (t3_1 - t3_0)
        const = c6_0 - slope * t3_0
        expected_slope = 2 * factorial(n-2)

        print(f"\n  Derived: tr(c_6) = {slope:.1f}*t3 + {const:.1f}")
        print(f"  Expected slope: 2*{n-2}! = {expected_slope}")
        print(f"  Match: {abs(slope - expected_slope) < 1}")

        # Verify with other data points
        for gs, H, t3, c0, c2, c4, c6, c8 in data:
            pred = slope * t3 + const
            print(f"    t3={t3}: tr(c_6)={c6:.1f}, predicted={pred:.1f}, err={abs(c6-pred):.1f}")

# =====================================================================
# Summary of constants
# =====================================================================
print("\n" + "=" * 70)
print("CONSTANTS IN c_{n-3} FORMULA")
print("=" * 70)

consts = {3: -0.5, 5: -30, 7: -2100}
print(f"\n  n=3: tr(c_0) = 2*1!*t3 + (-0.5)")
print(f"  n=5: tr(c_2) = 2*3!*t3 + (-30)")
print(f"  n=7: tr(c_4) = 2*5!*t3 + (-2100)")
if len(data) >= 2 and t3_0 != t3_1:
    consts[9] = const
    print(f"  n=9: tr(c_6) = 2*7!*t3 + ({const:.1f})")

# Check pattern in constants
print(f"\n  Constants: {[consts[n] for n in sorted(consts)]}")
ratios = []
keys = sorted(consts.keys())
for i in range(1, len(keys)):
    r = consts[keys[i]] / consts[keys[i-1]]
    ratios.append(r)
print(f"  Ratios: {[f'{r:.2f}' for r in ratios]}")

# Check: const(n) = -(n-1)!/(n-1) * something?
for nn in sorted(consts):
    c = consts[nn]
    print(f"  n={nn}: const={c}, const/(n-1)!={c/factorial(nn-1):.6f}, const/C(n,3)={c/(nn*(nn-1)*(nn-2)/6):.6f}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
