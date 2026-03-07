#!/usr/bin/env python3
"""
Complete w_1 formula at n=8 (even).

At n=8, w_1 is the "deepest" coefficient (analog of w_0 at odd n).
w_3 = 3024 - 480*t3 + 48*t5 + 96*bc  (exact, 0 error)
w_5 = -20160 + 1440*t3  (exact)
w_7 = 40320 = n!  (exact)

w_1 needs additional invariants. By the partition-OCF correspondence:
- Level k=3 patterns at n=8: (6,), (4,2), (2,2,2)
  - (6,): 7-vertex subtournaments -> t7 enters
  - (4,2): 5-vert + 3-vert -> bc35 type cross terms
  - (2,2,2): three 3-vert triples -> needs 9 vertices, IMPOSSIBLE at n=8

So w_1 at n=8 should depend on: t3, t5, t7, bc, bc35_w.
At n=8, bc35_w = sum over 3-cycles C of hc(T[V\C]) where V\C has 5 vertices.

opus-2026-03-06-S30
"""
from itertools import combinations
from math import factorial
from fractions import Fraction
import random
import numpy as np

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
    num_coeffs = n // 2
    r_sample = [0.1 * (k+1) for k in range(num_coeffs)]
    W_vals = [compute_W_dp(A, n, r) for r in r_sample]
    V = np.array([[r**(2*k+1) for k in range(num_coeffs)] for r in r_sample])
    return np.linalg.solve(V, W_vals)

def count_t3(A, n):
    return sum(1 for a,b,c in combinations(range(n),3)
               if A[a][b]*A[b][c]*A[c][a] or A[a][c]*A[c][b]*A[b][a])

def count_t5(A, n):
    t5 = 0
    for verts in combinations(range(n), 5):
        sub = [[A[verts[i]][verts[j]] for j in range(5)] for i in range(5)]
        dp = [[0]*5 for _ in range(1 << 5)]
        dp[1][0] = 1
        for m in range(1, 1 << 5):
            for v in range(5):
                if not (m & (1 << v)) or dp[m][v] == 0: continue
                for u in range(5):
                    if m & (1 << u): continue
                    if sub[v][u]: dp[m | (1 << u)][u] += dp[m][v]
        full = (1 << 5) - 1
        t5 += sum(dp[full][v] for v in range(1,5) if sub[v][0])
    return t5

def count_t7(A, n):
    t7 = 0
    for verts in combinations(range(n), 7):
        sub = [[A[verts[i]][verts[j]] for j in range(7)] for i in range(7)]
        dp = [[0]*7 for _ in range(1 << 7)]
        dp[1][0] = 1
        for m in range(1, 1 << 7):
            for v in range(7):
                if not (m & (1 << v)) or dp[m][v] == 0: continue
                for u in range(7):
                    if m & (1 << u): continue
                    if sub[v][u]: dp[m | (1 << u)][u] += dp[m][v]
        full = (1 << 7) - 1
        t7 += sum(dp[full][v] for v in range(1,7) if sub[v][0])
    return t7

def count_bc(A, n):
    cyc3 = [set(t) for t in combinations(range(n), 3)
            if A[t[0]][t[1]]*A[t[1]][t[2]]*A[t[2]][t[0]] or
               A[t[0]][t[2]]*A[t[2]][t[1]]*A[t[1]][t[0]]]
    return sum(1 for i in range(len(cyc3)) for j in range(i+1, len(cyc3))
               if cyc3[i].isdisjoint(cyc3[j]))

def count_bc35_w(A, n):
    """bc35_w = sum over 3-cycles C of hc(T[V\C])."""
    total = 0
    cyc3 = []
    for triple in combinations(range(n), 3):
        a, b, c = triple
        if A[a][b]*A[b][c]*A[c][a] or A[a][c]*A[c][b]*A[b][a]:
            cyc3.append(set(triple))
    for c3 in cyc3:
        remaining = sorted(v for v in range(n) if v not in c3)
        k = len(remaining)
        if k < 5: continue
        for quint in combinations(remaining, 5):
            sub = [[A[quint[i]][quint[j]] for j in range(5)] for i in range(5)]
            dp = [[0]*5 for _ in range(1 << 5)]
            dp[1][0] = 1
            for m in range(1, 1 << 5):
                for v in range(5):
                    if not (m & (1 << v)) or dp[m][v] == 0: continue
                    for u in range(5):
                        if m & (1 << u): continue
                        if sub[v][u]: dp[m | (1 << u)][u] += dp[m][v]
            full = (1 << 5) - 1
            total += sum(dp[full][v] for v in range(1,5) if sub[v][0])
    return total

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
print("=" * 70)
print("COMPLETE w_1 AT n=8")
print("=" * 70)

n = 8
num_samples = 25
data = []

print(f"\nComputing {num_samples} samples...")
for trial in range(num_samples):
    A = random_tournament(n, n*1000 + trial)
    wc = extract_W_coeffs_even(A, n)
    t3 = count_t3(A, n)
    t5 = count_t5(A, n)
    t7 = count_t7(A, n)
    bc = count_bc(A, n)
    bc35 = count_bc35_w(A, n)
    H = count_H(A, n)
    data.append({
        'wc': wc, 't3': t3, 't5': t5, 't7': t7,
        'bc': bc, 'bc35': bc35, 'H': H
    })
    if trial < 8:
        print(f"  T{trial}: t3={t3:3d} t5={t5:4d} t7={t7:5d} bc={bc:3d} bc35={bc35:5d} w1={wc[0]:.1f}")

# =====================================================================
print(f"\n{'='*70}")
print("REGRESSION: w_1 = f(t3, t5, t7, bc, bc35)")
print(f"{'='*70}")

inv_names = ['const', 't3', 't5', 't7', 'bc', 'bc35']
X = np.array([[1, d['t3'], d['t5'], d['t7'], d['bc'], d['bc35']] for d in data])
y = np.array([d['wc'][0] for d in data])

coeffs, res, rank, sv = np.linalg.lstsq(X, y, rcond=None)
y_pred = X @ coeffs
max_err = np.max(np.abs(y - y_pred))

print(f"\nw_1 = ", end="")
terms = []
for i, name in enumerate(inv_names):
    frac = Fraction(coeffs[i]).limit_denominator(10000)
    if abs(coeffs[i]) > 0.01:
        terms.append(f"{frac}*{name}")
print(" + ".join(terms))
print(f"Max error: {max_err:.6f}")

if max_err < 0.5:
    print("\nEXACT FIT!")
    print("\nExact coefficients:")
    for i, name in enumerate(inv_names):
        frac = Fraction(coeffs[i]).limit_denominator(10000)
        print(f"  {name}: {coeffs[i]:.6f} ≈ {frac}")
else:
    print(f"\nNot exact (err={max_err:.4f}). May need more invariants.")
    # Try also adding H directly
    print("\nTrying with H as additional variable:")
    X2 = np.column_stack([X, [d['H'] for d in data]])
    coeffs2, _, _, _ = np.linalg.lstsq(X2, y, rcond=None)
    y_pred2 = X2 @ coeffs2
    max_err2 = np.max(np.abs(y - y_pred2))
    inv2 = inv_names + ['H']
    print(f"  w_1 = ", end="")
    terms2 = []
    for i, name in enumerate(inv2):
        frac = Fraction(coeffs2[i]).limit_denominator(10000)
        if abs(coeffs2[i]) > 0.01:
            terms2.append(f"{frac}*{name}")
    print(" + ".join(terms2))
    print(f"  Max error: {max_err2:.6f}")

# =====================================================================
print(f"\n{'='*70}")
print("COMPLETE n=8 COEFFICIENT TABLE")
print(f"{'='*70}")

# Recompute all levels with full invariant set
for k in range(n//2):
    y = np.array([d['wc'][k] for d in data])
    coeffs_k, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    y_pred_k = X @ coeffs_k
    err_k = np.max(np.abs(y - y_pred_k))
    power = 2*k + 1
    terms = []
    for i, name in enumerate(inv_names):
        frac = Fraction(coeffs_k[i]).limit_denominator(10000)
        if abs(coeffs_k[i]) > 0.01:
            terms.append(f"{frac}*{name}")
    print(f"  w_{power} = {' + '.join(terms)}  (err={err_k:.6f})")

# =====================================================================
# Shift principle verification with correct invariants
print(f"\n{'='*70}")
print("SHIFT PRINCIPLE WITH CORRECT INVARIANTS")
print(f"{'='*70}")

# Extract per-invariant r-polynomials
per_inv = {name: [] for name in inv_names}
for k in range(n//2):
    y = np.array([d['wc'][k] for d in data])
    coeffs_k, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
    for i, name in enumerate(inv_names):
        per_inv[name].append(coeffs_k[i])

for name in inv_names:
    coeffs_list = per_inv[name]
    terms = []
    for k, c in enumerate(coeffs_list):
        if abs(c) > 0.01:
            frac = Fraction(c).limit_denominator(10000)
            terms.append(f"{frac}*r^{2*k+1}")
    if terms:
        print(f"  C_{name}(r) at n=8 = {' + '.join(terms)}")

print(f"\n  Shift principle checks:")
# C_t5(r) at n=8 should match C_t3(r) at n=6 = -8r + 48r^3
t5_poly = per_inv['t5']
print(f"    C_t5(r) at n=8: r coeff = {t5_poly[0]:.4f} (expect -8)")
print(f"                     r^3 coeff = {t5_poly[1]:.4f} (expect 48)")
print(f"                     r^5,r^7 = {t5_poly[2]:.4f}, {t5_poly[3]:.4f} (expect 0)")

# C_t7(r) at n=8 should match C_t5(r) at n=6 = 4r
t7_poly = per_inv['t7']
print(f"    C_t7(r) at n=8: r coeff = {t7_poly[0]:.4f} (expect 4)")
print(f"                     r^3,r^5,r^7 = {t7_poly[1]:.4f}, {t7_poly[2]:.4f}, {t7_poly[3]:.4f} (expect 0)")

# C_bc(r) at n=8 should match C_bc(r) at n=6 IF there's a bc shift principle
# At n=6: C_bc(r) = 8r
bc_poly = per_inv['bc']
print(f"    C_bc(r) at n=8: r coeff = {bc_poly[0]:.4f}")
print(f"                     r^3 coeff = {bc_poly[1]:.4f}")
print(f"    C_bc(r) at n=6 = 8r")

# =====================================================================
# Check W(1/2) = H with the OCF structure
print(f"\n{'='*70}")
print("OCF CHECK AT n=8")
print(f"{'='*70}")

for d in data[:5]:
    # At n=8: alpha_1 = t3+t5+t7, alpha_2 = bc+bc35, alpha_3 = 0 (need 9 verts)
    alpha_1 = d['t3'] + d['t5'] + d['t7']
    alpha_2 = d['bc'] + d['bc35']
    ocf = 1 + 2*alpha_1 + 4*alpha_2
    print(f"  H={d['H']} OCF={ocf} match={d['H']==ocf}")

print(f"\n{'='*70}")
print("DONE")
print(f"{'='*70}")
