#!/usr/bin/env python3
"""
f_poly_oeis_fast.py — Fast F(T,x) computation and OEIS sequence search.

Uses DP (Held-Karp) instead of brute-force permutation enumeration.

Author: opus-2026-03-07-S44
"""
from itertools import combinations
import math
import random
from functools import reduce
from math import gcd

def tournament_from_bits(bits, n):
    A = [[0]*n for _ in range(n)]
    pos = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> pos) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            pos += 1
    return A

def compute_F_dp(A, n):
    """Compute F(T,x) coefficients using Held-Karp DP.
    dp[mask][v][fwd] = # Hamiltonian paths ending at v using vertices in mask
                       with exactly fwd forward edges.
    """
    full = (1 << n) - 1
    # dp[mask][v] = dict mapping fwd_count -> number_of_paths
    dp = [[None]*(n) for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = {0: 1}

    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] is None:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                new_mask = mask | (1 << u)
                fwd_edge = A[v][u]  # 1 if v→u (forward), 0 if u→v
                if dp[new_mask][u] is None:
                    dp[new_mask][u] = {}
                for fwd_count, cnt in dp[mask][v].items():
                    new_fwd = fwd_count + fwd_edge
                    dp[new_mask][u][new_fwd] = dp[new_mask][u].get(new_fwd, 0) + cnt

    # Sum over all ending vertices
    F = [0] * n
    for v in range(n):
        if dp[full][v] is not None:
            for fwd_count, cnt in dp[full][v].items():
                F[fwd_count] += cnt
    return F

# ============================================================
# SPECIAL TOURNAMENT F POLYNOMIALS
# ============================================================
print("=" * 60)
print("F_k FOR SPECIAL TOURNAMENTS")
print("=" * 60)

print("\nTransitive tournament F(T,x) coefficients:")
for n in range(2, 10):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    F = compute_F_dp(A, n)
    print(f"  n={n}: F = {F}")

# ============================================================
# NUMBER OF DISTINCT F POLYNOMIALS
# ============================================================
print("\n" + "=" * 60)
print("NUMBER OF DISTINCT F(T,x) POLYNOMIALS")
print("=" * 60)

for n in range(2, 9):
    m = n*(n-1)//2
    if 1 << m > 10000000:
        # Sample
        random.seed(42)
        seen = set()
        for _ in range(50000):
            bits = random.getrandbits(m)
            A = tournament_from_bits(bits, n)
            F = compute_F_dp(A, n)
            seen.add(tuple(F))
        print(f"  n={n}: >= {len(seen)} distinct F (sampled from {1<<m} tournaments)")
    else:
        seen = set()
        for bits in range(1 << m):
            A = tournament_from_bits(bits, n)
            F = compute_F_dp(A, n)
            seen.add(tuple(F))
        print(f"  n={n}: {len(seen)} distinct F (exhaustive, {1<<m} tournaments)")

# ============================================================
# F(T, 2) FOR TRANSITIVE TOURNAMENT — OEIS SEARCH
# ============================================================
print("\n" + "=" * 60)
print("F(transitive, 2) SEQUENCE")
print("=" * 60)

seq_F2 = []
for n in range(2, 10):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            A[i][j] = 1
    F = compute_F_dp(A, n)
    F2 = sum(F[k] * 2**k for k in range(n))
    seq_F2.append(F2)

print(f"  F(transitive, 2): {seq_F2}")
# This is sum_sigma 2^{asc(sigma)} = A_n(2) where A_n is Eulerian polynomial
# A_n(2) is a known sequence!

# ============================================================
# MAX/MIN F_k AND H RANGE
# ============================================================
print("\n" + "=" * 60)
print("H RANGE AND MAX F_k")
print("=" * 60)

for n in range(2, 9):
    m = n*(n-1)//2

    if 1 << m > 1000000:
        random.seed(42)
        iterator = [random.getrandbits(m) for _ in range(10000)]
    else:
        iterator = range(1 << m)

    seen = set()
    all_F = []
    for bits in iterator:
        A = tournament_from_bits(bits, n)
        F = compute_F_dp(A, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)
        all_F.append(F)

    H_vals = sorted(set(F[n-1] for F in all_F))
    max_F1 = max(F[1] for F in all_F) if n > 1 else 0
    min_F1 = min(F[1] for F in all_F) if n > 1 else 0
    center = (n-1)//2
    max_Fc = max(F[center] for F in all_F)
    min_Fc = min(F[center] for F in all_F)

    print(f"  n={n}: {len(all_F)} distinct F, H in [{H_vals[0]}, {H_vals[-1]}], "
          f"F_1 in [{min_F1}, {max_F1}], F_{center} in [{min_Fc}, {max_Fc}]")

# ============================================================
# AVERAGE F POLYNOMIAL — COMPARE WITH n!/2^{n-1} * (1+x)^{n-1}
# ============================================================
print("\n" + "=" * 60)
print("AVERAGE F(T,x) OVER ALL TOURNAMENTS")
print("=" * 60)

for n in range(2, 8):
    m = n*(n-1)//2
    total = [0] * n
    for bits in range(1 << m):
        A = tournament_from_bits(bits, n)
        F = compute_F_dp(A, n)
        for k in range(n):
            total[k] += F[k]

    num_T = 1 << m
    avg = [total[k] / num_T for k in range(n)]
    expected = [math.factorial(n) * math.comb(n-1, k) / 2**(n-1) for k in range(n)]
    match = all(abs(avg[k] - expected[k]) < 0.01 for k in range(n))
    print(f"  n={n}: avg F_k = {[f'{a:.1f}' for a in avg]}")
    print(f"         n!/2^(n-1)*C(n-1,k) = {[f'{e:.1f}' for e in expected]}")
    print(f"         match: {match}")

# ============================================================
# NOVEL: VARIANCE OF F_k ACROSS TOURNAMENTS
# ============================================================
print("\n" + "=" * 60)
print("VARIANCE OF F_k ACROSS ALL TOURNAMENTS")
print("=" * 60)

for n in range(2, 8):
    m = n*(n-1)//2
    totals = [0] * n
    sq_totals = [0] * n
    num_T = 1 << m

    for bits in range(num_T):
        A = tournament_from_bits(bits, n)
        F = compute_F_dp(A, n)
        for k in range(n):
            totals[k] += F[k]
            sq_totals[k] += F[k]**2

    vars = [(sq_totals[k]/num_T - (totals[k]/num_T)**2) for k in range(n)]
    print(f"  n={n}: var(F_k) = {[f'{v:.1f}' for v in vars]}")

# ============================================================
# S(T) DISTRIBUTION
# ============================================================
print("\n" + "=" * 60)
print("S(T) = SIGNED HP PERMANENT DISTRIBUTION")
print("=" * 60)

for n in range(2, 9):
    m = n*(n-1)//2

    if 1 << m > 1000000:
        random.seed(42)
        iterator = [random.getrandbits(m) for _ in range(10000)]
    else:
        iterator = range(1 << m)

    seen = set()
    S_vals = []
    for bits in iterator:
        A = tournament_from_bits(bits, n)
        F = compute_F_dp(A, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)
        S = (-1)**(n-1) * sum(F[k] * (-1)**k for k in range(n))
        S_vals.append(S)

    S_set = sorted(set(S_vals))
    g = reduce(gcd, [abs(s) for s in S_set if s != 0]) if any(s != 0 for s in S_set) else 0
    print(f"  n={n}: |S| range [{min(abs(s) for s in S_set)}, {max(abs(s) for s in S_set)}], "
          f"# distinct={len(S_set)}, GCD={g}")

# ============================================================
# F(omega) DIVISIBILITY ACROSS n
# ============================================================
print("\n" + "=" * 60)
print("F(T, omega) DIVISIBILITY BY 9 — ACROSS n")
print("=" * 60)

for n in range(3, 9):
    m = n*(n-1)//2
    d = n - 1

    if 1 << m > 500000:
        random.seed(42)
        iterator = [random.getrandbits(m) for _ in range(5000)]
    else:
        iterator = range(1 << m)

    seen = set()
    f_omega_vals = []
    for bits in iterator:
        A = tournament_from_bits(bits, n)
        F = compute_F_dp(A, n)
        key = tuple(F)
        if key in seen:
            continue
        seen.add(key)

        # F(omega) when d mod 3 == 0: F(omega) = (3*S0 - n!)/2 (real)
        # When d mod 3 != 0: F(omega) is complex
        S0 = sum(F[k] for k in range(n) if k % 3 == 0)
        S1 = sum(F[k] for k in range(n) if k % 3 == 1)
        S2 = sum(F[k] for k in range(n) if k % 3 == 2)

        # 2*Re(F(omega)) = 2*S0 - S1 - S2 = 3*S0 - n!
        two_re = 3 * S0 - math.factorial(n)
        # 2*Im(F(omega))/sqrt(3) = S1 - S2
        two_im_scaled = S1 - S2

        f_omega_vals.append((two_re, two_im_scaled))

    re_vals = [v[0] for v in f_omega_vals]
    im_vals = [v[1] for v in f_omega_vals]

    re_g = reduce(gcd, [abs(v) for v in re_vals if v != 0]) if any(v != 0 for v in re_vals) else 0
    im_g = reduce(gcd, [abs(v) for v in im_vals if v != 0]) if any(v != 0 for v in im_vals) else 0

    is_real = all(v == 0 for v in im_vals)

    print(f"  n={n} (d={d}, d%3={d%3}): "
          f"2*Re GCD={re_g}, Im=0?={is_real}, "
          f"2*Im/√3 GCD={im_g if not is_real else 'N/A'}")
    if is_real:
        # F(omega) = two_re / 2
        # GCD of F(omega) = re_g / 2 (if re_g even) or re_g (need to check)
        if re_g % 2 == 0:
            print(f"         F(omega) GCD = {re_g // 2}")
        else:
            print(f"         2*F(omega) GCD = {re_g}")
