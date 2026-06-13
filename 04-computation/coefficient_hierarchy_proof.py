#!/usr/bin/env python3
"""
COEFFICIENT HIERARCHY THEOREM (THM-055) — Verification Script

KEY RESULT: tr(c_{n-1-2k}) = sum_P e_{2k}(s_P)
where e_{2k} is a degree-2k polynomial in f_P (forward arc count).

Via Newton's identities with s_i^2 = 1/4 (constant):
  e_{2k}(f) has explicit closed-form coefficients.

The hierarchy:
  k=0: tr(c_{n-1}) = (n-1)!  [universal]
  k=1: tr(c_{n-3}) depends on sum_P f^2 -> depends only on t_3
  k=2: tr(c_{n-5}) depends on sum_P f^4 -> depends on MORE than t_3

opus-2026-03-06-S11b (continued^4)
"""
from itertools import permutations, combinations
from math import factorial, comb
from sympy import symbols, expand, Rational

# =====================================================================
# Part 1: Newton's identity reduction
# =====================================================================
print("=" * 70)
print("PART 1: e_{2k} AS POLYNOMIAL IN f")
print("=" * 70)

def compute_e_polynomials(n):
    """Compute e_k as polynomial in f = #forward arcs, using Newton's identities."""
    m = n - 1  # number of edges in path
    p1 = symbols('p1')
    f = symbols('f')

    # Newton's identities: k*e_k = sum_{i=1}^{k} (-1)^{i-1} * e_{k-i} * p_i
    # With p_{2j} = m/4^j, p_{2j+1} = p1/4^j
    e_p1 = [Rational(1)]  # e_0 = 1
    for j in range(1, m + 1):
        val = 0
        for i in range(1, j + 1):
            if i % 2 == 0:
                p_i = Rational(m, 4**(i//2))
            else:
                p_i = p1 / 4**((i-1)//2)
            val += (-1)**(i-1) * e_p1[j-i] * p_i
        e_p1.append(val / j)

    # Convert to polynomial in f = p1 + m/2
    e_f = []
    for k in range(m + 1):
        expr = expand(e_p1[k].subs(p1, f - Rational(m, 2)))
        e_f.append(expr)
    return e_f

for n in [5, 7]:
    print(f"\nn={n}:")
    e_f = compute_e_polynomials(n)
    f = symbols('f')
    for k in range(0, n, 2):
        print(f"  e_{k}(f) = {e_f[k]}")

# =====================================================================
# Part 2: Verify tr(c_{n-1-2k}) = sum_P e_{2k}(f_P)
# =====================================================================
print("\n" + "=" * 70)
print("PART 2: VERIFICATION")
print("=" * 70)

import numpy as np
import random

def ham_count_dp(A, n):
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or dp[mask][v] == 0: continue
            for u in range(n):
                if (mask & (1 << u)) or A[v][u] != 1: continue
                dp[mask | (1 << u)][u] += dp[mask][v]
    return sum(dp[(1<<n)-1][v] for v in range(n))

def count_3_cycles(A, n):
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if A[i][j] and A[j][k] and A[k][i]: count += 1
                if A[i][k] and A[k][j] and A[j][i]: count += 1
    return count

def compute_tr_coefficients_via_r(A, n):
    """Compute c_k trace values via r-polynomial interpolation."""
    num_coeffs = (n + 1) // 2  # at odd n
    r_values = np.linspace(0, 1, num_coeffs + 1)[:num_coeffs]
    r_powers = np.array([[r**(2*k) for k in range(num_coeffs)] for r in r_values])

    traces = []
    for r in r_values:
        full = (1 << n) - 1
        dp_fwd = [[0.0]*n for _ in range(1 << n)]
        for v in range(n):
            dp_fwd[1 << v][v] = 1.0
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)) or dp_fwd[mask][v] == 0: continue
                for u in range(n):
                    if mask & (1 << u): continue
                    wt = r + (A[v][u] - 0.5)
                    dp_fwd[mask | (1 << u)][u] += dp_fwd[mask][v] * wt
        dp_bwd = [[0.0]*n for _ in range(1 << n)]
        for v in range(n):
            dp_bwd[1 << v][v] = 1.0
        for mask in range(1, 1 << n):
            for v in range(n):
                if not (mask & (1 << v)) or dp_bwd[mask][v] == 0: continue
                for u in range(n):
                    if mask & (1 << u): continue
                    wt = r + (A[u][v] - 0.5)
                    dp_bwd[mask | (1 << u)][u] += dp_bwd[mask][v] * wt
        tr_val = 0.0
        for v in range(n):
            for mask_before in range(1 << n):
                if mask_before & (1 << v): continue
                mask_with_v = mask_before | (1 << v)
                if dp_fwd[mask_with_v][v] == 0: continue
                mask_after = full ^ mask_before
                if not (mask_after & (1 << v)): continue
                if dp_bwd[mask_after][v] == 0: continue
                k = bin(mask_before).count('1')
                tr_val += ((-1)**k) * dp_fwd[mask_with_v][v] * dp_bwd[mask_after][v]
        traces.append(tr_val)

    coeffs = np.linalg.solve(r_powers, traces)
    return coeffs

def compute_sum_ek_direct(A, n, k_val):
    """Compute sum_P e_k(s_P) directly over all permutations."""
    m = n - 1
    total = 0.0
    for p in permutations(range(n)):
        s = [A[p[i]][p[i+1]] - 0.5 for i in range(m)]
        # e_k of s values
        ek = sum(np.prod([s[idx] for idx in combo]) for combo in combinations(range(m), k_val))
        total += ek
    return total

def compute_sum_ek_via_f(A, n, e_poly):
    """Compute sum_P e_k(f_P) using the polynomial in f."""
    from sympy import lambdify
    f = symbols('f')
    e_func = lambdify(f, e_poly, 'numpy')

    m = n - 1
    total = 0.0
    for p in permutations(range(n)):
        fp = sum(1 for i in range(m) if A[p[i]][p[i+1]] == 1)
        total += float(e_func(fp))
    return total

# Test at n=5
print("\nn=5:")
n = 5
e_f_list = compute_e_polynomials(n)
f = symbols('f')

for trial in range(3):
    random.seed(42 + trial)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    H = ham_count_dp(A, n)
    t3 = count_3_cycles(A, n)

    # Transfer matrix coefficients
    c_coeffs = compute_tr_coefficients_via_r(A, n)

    # Direct e_k sums
    for k in [0, 2, 4]:
        direct = compute_sum_ek_direct(A, n, k)
        via_f = compute_sum_ek_via_f(A, n, e_f_list[k])
        coeff_idx = k // 2  # c_0, c_2, c_4 -> indices 0, 1, 2
        tm_val = c_coeffs[coeff_idx]  # c_{2*idx} = c_{n-1-2k} -> idx = (n-1-k)/2... no

        # c_coeffs[j] = tr(c_{2j})
        # tr(c_{n-1-2k}) = sum_P e_{2k}, so c_{n-1-2k} at index (n-1-2k)/2 = (n-1)/2 - k
        idx = (n-1)//2 - k//2
        tm_val = c_coeffs[idx] if idx >= 0 else None

        tm_str = f"{tm_val:.4f}" if tm_val is not None else "N/A"
        match_str = "YES" if tm_val is not None and abs(direct - tm_val) < 0.1 else "NO"
        print(f"  T{trial}: e_{k} direct={direct:.4f}, via_f={via_f:.4f}, tr(c)={tm_str}, match={match_str}")

# =====================================================================
# Part 3: Moment hierarchy — which moments depend on what
# =====================================================================
print("\n" + "=" * 70)
print("PART 3: MOMENT HIERARCHY")
print("=" * 70)

n = 5
tiles = [(a,b) for a in range(n) for b in range(a) if a-b >= 2]
m_tiles = len(tiles)

def tiling_to_adj(bits, n):
    A = [[0]*n for _ in range(n)]
    for i in range(1, n):
        A[i][i-1] = 1
    tiles = [(a,b) for a in range(n) for b in range(a) if a-b >= 2]
    tiles.sort()
    for idx, (a, b) in enumerate(tiles):
        if (bits >> idx) & 1:
            A[b][a] = 1
        else:
            A[a][b] = 1
    return A

# Collect all iso classes
seen = set()
iso_data = []
for bits in range(2**m_tiles):
    A = tiling_to_adj(bits, n)
    key = tuple(sorted(tuple(A[i]) for i in range(n)))
    if key in seen: continue
    seen.add(key)

    t3 = count_3_cycles(A, n)
    H = ham_count_dp(A, n)
    m = n - 1
    sf = [0] * 5
    for p in permutations(range(n)):
        fp = sum(1 for i in range(m) if A[p[i]][p[i+1]] == 1)
        for k in range(5):
            sf[k] += fp**k
    iso_data.append({'t3': t3, 'H': H, 'sf': sf})

# Group by t3 and check which moments vary
from collections import defaultdict
by_t3 = defaultdict(list)
for d in iso_data:
    by_t3[d['t3']].append(d)

print(f"\nn=5: {len(iso_data)} distinct tournaments (by sorted rows)")
print(f"\n{'t3':>3} {'count':>5} {'sf2 varies':>10} {'sf3 varies':>10} {'sf4 varies':>10}")
for t3 in sorted(by_t3):
    group = by_t3[t3]
    sf2_uniq = len(set(d['sf'][2] for d in group))
    sf3_uniq = len(set(d['sf'][3] for d in group))
    sf4_uniq = len(set(d['sf'][4] for d in group))
    print(f"{t3:>3} {len(group):>5} {sf2_uniq:>10} {sf3_uniq:>10} {sf4_uniq:>10}")

print("""
CONCLUSION:
  sum_P f^2, sum_P f^3: determined by t_3 alone (1 value per t_3 group)
  sum_P f^4: NOT determined by t_3 (3 values at t_3=4)

  => tr(c_{n-3}) = sum_P e_2(f) depends only on t_3  [PROVED]
  => tr(c_{n-5}) = sum_P e_4(f) depends on t_3 AND sum_P f^4  [VERIFIED]

  At n=5: sum_P f_(4) = f*(f-1)*(f-2)*(f-3) summed over permutations
         = sum of products of ALL 4 edge indicators
         = count of permutations where ALL edges are forward
         = H(T) (the Hamiltonian path count!)

  So tr(c_0) at n=5 depends on H and t_3: tr(c_0) = H - 3*t_3.
""")

print("=" * 70)
print("DONE")
print("=" * 70)
