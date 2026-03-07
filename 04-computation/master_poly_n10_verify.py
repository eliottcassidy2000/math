#!/usr/bin/env python3
"""
Verify Universal Master Polynomial (THM-059) at n=10 (even).

At n=10, W(r) has only odd powers: r^1, r^3, r^5, r^7, r^9.
Free position count for each invariant:
  C_0(r) = F_9(r)           (free = 9)
  C_t3(r) = 2*F_7(r)        (free = 7)
  C_t5(r) = 2*F_5(r)        (free = 5)
  C_t7(r) = 2*F_3(r)        (free = 3)
  C_t9(r) = 2*F_1(r)        (free = 1)
  C_bc(r) = 4*F_5(r)        (free = 5)
  C_bc35(r) = 4*F_3(r)      (free = 3)
  C_bc37(r) = 4*F_1(r)      (free = 1)
  C_a3(r)  = 8*F_3(r)       (free = 3)

New predictions to verify:
  F_9(r) = ?*r + ?*r^3 + ?*r^5 + ?*r^7 + ?*r^9
  from b_{4,j} = [1, 85, 147, 30, 1] and the odd formula.

opus-2026-03-06-S31
"""
from itertools import combinations
from math import factorial, comb
from fractions import Fraction
import numpy as np
import random

# =====================================================================
# Compute b-triangle and master polynomials
# =====================================================================

def compute_b_triangle(max_k):
    """Compute b_{k,j} with recurrence b_{k,j} = b_{k-1,j-1} + (j+1)^2 * b_{k-1,j}."""
    b = {}
    b[(0,0)] = Fraction(1)
    for k in range(1, max_k+1):
        for j in range(k+1):
            prev_diag = b.get((k-1, j-1), Fraction(0))
            prev_above = b.get((k-1, j), Fraction(0))
            b[(k,j)] = prev_diag + (j+1)**2 * prev_above
    return b

def compute_F_even(k, b_tri):
    """Compute F_{2k}(r) as list of coefficients [r^0, r^2, r^4, ...]."""
    # F_{2k}(r) = sum_j b_{k,j} * (2j+1)! * u^j, u = r^2 - 1/4
    # Expand u^j in powers of r^2
    result = [Fraction(0)] * (k+1)
    for j in range(k+1):
        bv = b_tri[(k,j)]
        fact = Fraction(factorial(2*j+1))
        # u^j = (r^2 - 1/4)^j = sum_m C(j,m) r^{2m} (-1/4)^{j-m}
        for m in range(j+1):
            coeff = bv * fact * Fraction(comb(j,m)) * Fraction(-1,4)**(j-m)
            result[m] += coeff
    return result

def compute_F_odd(k, b_tri):
    """Compute F_{2k+1}(r) as list of coefficients [r^1, r^3, r^5, ...]."""
    # F_{2k+1}(r) = r * sum_j b_{k,j} * (2j+2)! * u^j
    result = [Fraction(0)] * (k+1)
    for j in range(k+1):
        bv = b_tri[(k,j)]
        fact = Fraction(factorial(2*j+2))
        for m in range(j+1):
            coeff = bv * fact * Fraction(comb(j,m)) * Fraction(-1,4)**(j-m)
            result[m] += coeff
    return result

b_tri = compute_b_triangle(5)

print("b-triangle:")
for k in range(6):
    row = [b_tri[(k,j)] for j in range(k+1)]
    print(f"  k={k}: {row}")

print("\nMaster polynomials:")
for idx in range(10):
    if idx % 2 == 0:
        k = idx // 2
        coeffs = compute_F_even(k, b_tri)
        powers = [f"{c}*r^{2*i}" for i, c in enumerate(coeffs) if c != 0]
    else:
        k = (idx-1) // 2
        coeffs = compute_F_odd(k, b_tri)
        powers = [f"{c}*r^{2*i+1}" for i, c in enumerate(coeffs) if c != 0]

    # Check F(1/2) = 1
    if idx % 2 == 0:
        val = sum(c * Fraction(1,2)**(2*i) for i, c in enumerate(coeffs))
    else:
        val = sum(c * Fraction(1,2)**(2*i+1) for i, c in enumerate(coeffs))

    print(f"  F_{idx}(r) = {' + '.join(powers[:4])}{'...' if len(powers)>4 else ''}  F({idx},1/2)={val}")

# =====================================================================
# Tournament computation functions
# =====================================================================

def compute_W_dp(A, n, r_val):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1.0
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            val = dp.get((mask, v), 0)
            if val == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                wt = r_val + (A[v][u] - 0.5)
                key = (mask | (1 << u), u)
                dp[key] = dp.get(key, 0) + val * wt
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def extract_W_coeffs_even(A, n):
    """Extract odd-power coefficients for even n."""
    num_coeffs = n // 2
    r_sample = [0.1 * (k+1) for k in range(num_coeffs)]
    W_vals = [compute_W_dp(A, n, r) for r in r_sample]
    V = np.array([[r**(2*k+1) for k in range(num_coeffs)] for r in r_sample])
    return np.linalg.solve(V, W_vals)

def count_t3(A, n):
    return sum(1 for a,b,c in combinations(range(n),3)
               if A[a][b]*A[b][c]*A[c][a] or A[a][c]*A[c][b]*A[b][a])

def count_hc_sub(A, n, verts):
    """Count directed Hamiltonian cycles on a subset of vertices."""
    k = len(verts)
    if k < 3: return 0
    sub = [[A[verts[i]][verts[j]] for j in range(k)] for i in range(k)]
    dp = [[0]*k for _ in range(1 << k)]
    dp[1][0] = 1
    for m in range(1, 1 << k):
        for v in range(k):
            if not (m & (1 << v)) or dp[m][v] == 0: continue
            for u in range(k):
                if m & (1 << u): continue
                if sub[v][u]: dp[m | (1 << u)][u] += dp[m][v]
    full = (1 << k) - 1
    return sum(dp[full][v] for v in range(1,k) if sub[v][0])

def count_t5(A, n):
    return sum(count_hc_sub(A, n, list(verts)) for verts in combinations(range(n), 5))

def count_t7(A, n):
    return sum(count_hc_sub(A, n, list(verts)) for verts in combinations(range(n), 7))

def count_t9(A, n):
    return sum(count_hc_sub(A, n, list(verts)) for verts in combinations(range(n), 9))

def count_bc(A, n):
    """Count pairs of vertex-disjoint directed 3-cycles."""
    cyc3 = [set(t) for t in combinations(range(n), 3)
            if A[t[0]][t[1]]*A[t[1]][t[2]]*A[t[2]][t[0]] or
               A[t[0]][t[2]]*A[t[2]][t[1]]*A[t[1]][t[0]]]
    return sum(1 for i in range(len(cyc3)) for j in range(i+1, len(cyc3))
               if cyc3[i].isdisjoint(cyc3[j]))

def count_bc35_w(A, n):
    """Weighted bc35: sum over 3-cycles C of hc(T[complement(C)] restricted to 5-subsets)."""
    total = 0
    for triple in combinations(range(n), 3):
        a, b, c = triple
        if A[a][b]*A[b][c]*A[c][a] or A[a][c]*A[c][b]*A[b][a]:
            remaining = [v for v in range(n) if v not in triple]
            for quint in combinations(remaining, 5):
                total += count_hc_sub(A, n, list(quint))
    return total

def count_bc37_w(A, n):
    """Weighted bc37: sum over 3-cycles C of hc(T[complement(C)] restricted to 7-subsets)."""
    total = 0
    for triple in combinations(range(n), 3):
        a, b, c = triple
        if A[a][b]*A[b][c]*A[c][a] or A[a][c]*A[c][b]*A[b][a]:
            remaining = [v for v in range(n) if v not in triple]
            for sept in combinations(remaining, 7):
                total += count_hc_sub(A, n, list(sept))
    return total

def count_bc55_w(A, n):
    """Pairs of vertex-disjoint 5-cycles: sum over {S1,S2} partitions into 5+5 of hc(S1)*hc(S2)."""
    total = 0
    for s1 in combinations(range(n), 5):
        s2 = tuple(v for v in range(n) if v not in s1)
        if s1[0] < s2[0]:  # avoid double counting
            hc1 = count_hc_sub(A, n, list(s1))
            if hc1 > 0:
                hc2 = count_hc_sub(A, n, list(s2))
                total += hc1 * hc2
    return total

def count_alpha3(A, n):
    """Triple of vertex-disjoint 3-cycles (weighted by cycle counts)."""
    total = 0
    for v9 in combinations(range(n), 9):
        for part in _partitions_333(v9):
            p = 1
            for triple in part:
                hc = count_hc_sub(A, n, list(triple))
                p *= hc
            total += p
    return total

def _partitions_333(vertices):
    """Generate all ways to partition 9 vertices into three triples."""
    v = list(vertices)
    results = []
    for t1 in combinations(v, 3):
        rem1 = [x for x in v if x not in t1]
        for t2 in combinations(rem1, 3):
            if t2[0] < t1[0]: continue  # avoid duplicates
            t3 = tuple(x for x in rem1 if x not in t2)
            if t3[0] < t2[0]: continue
            results.append((t1, t2, t3))
    return results

def count_H(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            c = dp.get((mask, v), 0)
            if c == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + c
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def random_tournament(n, seed):
    random.seed(seed)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    return A

# =====================================================================
# Main verification at n=10
# =====================================================================

n = 10
num_samples = 20

print(f"\n{'='*70}")
print(f"VERIFICATION AT n={n}")
print(f"{'='*70}")

# Predicted F_9 from b-triangle row k=4: [1, 85, 147, 30, 1]
F9_coeffs = compute_F_odd(4, b_tri)
print(f"\nPredicted F_9(r) = {' + '.join(f'{c}*r^{2*i+1}' for i, c in enumerate(F9_coeffs))}")
F9_at_half = sum(c * Fraction(1,2)**(2*i+1) for i, c in enumerate(F9_coeffs))
print(f"F_9(1/2) = {F9_at_half}")

# F_7 coefficients
F7_coeffs = compute_F_odd(3, b_tri)
F5_coeffs = compute_F_odd(2, b_tri)
F3_coeffs = compute_F_odd(1, b_tri)
F1_coeffs = compute_F_odd(0, b_tri)

print(f"\nExpected C_0(r) = F_9(r)")
print(f"Expected C_t3(r) = 2*F_7(r) = {[2*c for c in F7_coeffs]}")
print(f"Expected C_t5(r) = 2*F_5(r) = {[2*c for c in F5_coeffs]}")
print(f"Expected C_t7(r) = 2*F_3(r) = {[2*c for c in F3_coeffs]}")
print(f"Expected C_t9(r) = 2*F_1(r) = {[2*c for c in F1_coeffs]}")
print(f"Expected C_bc(r) = 4*F_5(r) = {[4*c for c in F5_coeffs]}")
print(f"Expected C_bc35(r) = 4*F_3(r) = {[4*c for c in F3_coeffs]}")
print(f"Expected C_bc37(r) = 4*F_1(r) = {[4*c for c in F1_coeffs]}")
print(f"Expected C_a3(r) = 8*F_3(r) = {[8*c for c in F3_coeffs]}")

print(f"\nComputing {num_samples} random tournaments at n={n}...")
data = []
for trial in range(num_samples):
    A = random_tournament(n, n*1000 + trial)
    wc = extract_W_coeffs_even(A, n)
    t3 = count_t3(A, n)
    t5 = count_t5(A, n)
    t7 = count_t7(A, n)
    t9 = count_t9(A, n)
    bc = count_bc(A, n)
    bc35 = count_bc35_w(A, n)
    bc37 = count_bc37_w(A, n)
    bc55 = count_bc55_w(A, n)
    a3 = count_alpha3(A, n)
    H = count_H(A, n)
    data.append({
        'wc': wc, 't3': t3, 't5': t5, 't7': t7, 't9': t9,
        'bc': bc, 'bc35': bc35, 'bc37': bc37, 'bc55': bc55, 'a3': a3, 'H': H
    })
    print(f"  T{trial}: t3={t3:4d} t5={t5:5d} t7={t7:6d} bc={bc:4d} bc55={bc55:5d} H={H:6d} w_1={wc[0]:.1f}")

# =====================================================================
# Regression: extract per-invariant r-polynomials
# =====================================================================

inv_names = ['const', 't3', 't5', 't7', 't9', 'bc', 'bc35', 'bc37', 'bc55', 'a3']
X = np.array([[1, d['t3'], d['t5'], d['t7'], d['t9'],
               d['bc'], d['bc35'], d['bc37'], d['bc55'], d['a3']] for d in data])

print(f"\n{'='*70}")
print("PER-COEFFICIENT REGRESSION")
print(f"{'='*70}")

per_inv = {name: [] for name in inv_names}
for k in range(n//2):
    y = np.array([d['wc'][k] for d in data])
    coeffs_k, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    y_pred = X @ coeffs_k
    err = np.max(np.abs(y - y_pred))
    power = 2*k + 1
    terms = []
    for i, name in enumerate(inv_names):
        frac = Fraction(coeffs_k[i]).limit_denominator(100000)
        if abs(coeffs_k[i]) > 0.01:
            terms.append(f"{frac}*{name}")
        per_inv[name].append(coeffs_k[i])
    print(f"  w_{power} = {' + '.join(terms)}  (err={err:.6f}, rank={rank})")

# =====================================================================
# Verify against THM-059 predictions
# =====================================================================

print(f"\n{'='*70}")
print("THM-059 PREDICTION VERIFICATION")
print(f"{'='*70}")

# Check background polynomial C_0(r) = F_9(r)
print("\nC_0(r) at n=10 vs F_9(r):")
for k in range(n//2):
    actual = per_inv['const'][k]
    predicted = float(F9_coeffs[k]) if k < len(F9_coeffs) else 0
    match = abs(actual - predicted) < 1.0
    print(f"  r^{2*k+1}: actual={actual:.2f}, predicted={predicted:.2f}, {'OK' if match else 'FAIL'}")

# Check C_t3(r) = 2*F_7(r)
print("\nC_t3(r) at n=10 vs 2*F_7(r):")
for k in range(n//2):
    actual = per_inv['t3'][k]
    predicted = float(2*F7_coeffs[k]) if k < len(F7_coeffs) else 0
    match = abs(actual - predicted) < 1.0
    print(f"  r^{2*k+1}: actual={actual:.2f}, predicted={predicted:.2f}, {'OK' if match else 'FAIL'}")

# Check C_t5(r) = 2*F_5(r)
print("\nC_t5(r) at n=10 vs 2*F_5(r):")
for k in range(min(n//2, len(F5_coeffs)+1)):
    actual = per_inv['t5'][k]
    predicted = float(2*F5_coeffs[k]) if k < len(F5_coeffs) else 0
    match = abs(actual - predicted) < 1.0
    print(f"  r^{2*k+1}: actual={actual:.2f}, predicted={predicted:.2f}, {'OK' if match else 'FAIL'}")

# Check C_bc(r) = 4*F_5(r)
print("\nC_bc(r) at n=10 vs 4*F_5(r):")
for k in range(min(n//2, len(F5_coeffs)+1)):
    actual = per_inv['bc'][k]
    predicted = float(4*F5_coeffs[k]) if k < len(F5_coeffs) else 0
    match = abs(actual - predicted) < 1.0
    print(f"  r^{2*k+1}: actual={actual:.2f}, predicted={predicted:.2f}, {'OK' if match else 'FAIL'}")

# OCF check
print(f"\n{'='*70}")
print("OCF CHECK")
print(f"{'='*70}")
for i, d in enumerate(data[:5]):
    alpha1 = d['t3'] + d['t5'] + d['t7'] + d['t9']
    alpha2 = d['bc'] + d['bc35'] + d['bc37'] + d['bc55']
    alpha3 = d['a3']
    ocf = 1 + 2*alpha1 + 4*alpha2 + 8*alpha3
    print(f"  T{i}: H={d['H']}, OCF={ocf}, match={d['H']==ocf}")

print(f"\n{'='*70}")
print("DONE")
print(f"{'='*70}")
