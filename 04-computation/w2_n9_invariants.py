#!/usr/bin/env python3
"""
w_2 at n=9: what invariants appear beyond t3, t5, bc?

At n=9, w_{n-7} = w_2. This is the coefficient of r^2 in W(r).
By the position pattern analysis, sigma_6 contributes, which involves:
- Pattern (6,): 7-vertex sub-tournaments → t7 enters
- Pattern (4,2): 5-vert + 3-vert → cross terms
- Pattern (2,2,2): three 3-vert → triple cycle count (alpha_3)
- Pattern (3,3): two 4-vert → new cross terms

This script extracts the exact w_2 formula at n=9 via regression
over computed samples, identifying all necessary invariants.

opus-2026-03-06-S29
"""
from itertools import permutations, combinations
from math import factorial, comb
import random
import numpy as np

def compute_W_dp(A, n, r_val):
    """Compute W(r) via DP."""
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

def extract_W_coeffs(A, n):
    """Extract all W(r) coefficients."""
    num_coeffs = (n + 1) // 2
    r_sample = [0.05 * (k+1) for k in range(num_coeffs)]
    W_vals = [compute_W_dp(A, n, r) for r in r_sample]
    V = np.array([[r**(2*k) for k in range(num_coeffs)] for r in r_sample])
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
        hc = sum(dp[full][v] for v in range(1,5) if sub[v][0])
        t5 += hc
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
        hc = sum(dp[full][v] for v in range(1,7) if sub[v][0])
        t7 += hc
    return t7

def count_bc(A, n):
    cyc_triples = [set(t) for t in combinations(range(n), 3)
                   if A[t[0]][t[1]]*A[t[1]][t[2]]*A[t[2]][t[0]] or
                      A[t[0]][t[2]]*A[t[2]][t[1]]*A[t[1]][t[0]]]
    return sum(1 for i in range(len(cyc_triples)) for j in range(i+1, len(cyc_triples))
               if cyc_triples[i].isdisjoint(cyc_triples[j]))

def count_bc35_w(A, n):
    """bc35_w = weighted count of disjoint (3-cycle, 5-cycle) pairs."""
    total = 0
    cyc3 = []
    for triple in combinations(range(n), 3):
        a, b, c = triple
        if A[a][b]*A[b][c]*A[c][a] or A[a][c]*A[c][b]*A[b][a]:
            cyc3.append(set(triple))

    for c3 in cyc3:
        remaining = [v for v in range(n) if v not in c3]
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
            hc = sum(dp[full][v] for v in range(1,5) if sub[v][0])
            total += hc
    return total

def count_alpha3(A, n):
    """alpha_3 = # triples of pairwise disjoint directed 3-cycles."""
    cyc3 = []
    for triple in combinations(range(n), 3):
        a, b, c = triple
        if A[a][b]*A[b][c]*A[c][a] or A[a][c]*A[c][b]*A[b][a]:
            cyc3.append(set(triple))
    total = 0
    for i in range(len(cyc3)):
        for j in range(i+1, len(cyc3)):
            if not cyc3[i].isdisjoint(cyc3[j]): continue
            for k in range(j+1, len(cyc3)):
                if cyc3[k].isdisjoint(cyc3[i]) and cyc3[k].isdisjoint(cyc3[j]):
                    total += 1
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

# =====================================================================
print("=" * 70)
print("w_2 AT n=9: INVARIANT EXTRACTION")
print("=" * 70)

n = 9
num_samples = 25
data = []

print(f"\nComputing {num_samples} samples at n={n}...")
for trial in range(num_samples):
    random.seed(n*1000 + trial)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1

    w_coeffs = extract_W_coeffs(A, n)
    w2 = w_coeffs[1]  # coefficient of r^2

    t3 = count_t3(A, n)
    t5 = count_t5(A, n)
    t7 = count_t7(A, n)
    bc = count_bc(A, n)
    bc35 = count_bc35_w(A, n)
    a3 = count_alpha3(A, n)
    H = count_H(A, n)

    data.append({
        'w2': w2, 't3': t3, 't5': t5, 't7': t7,
        'bc': bc, 'bc35': bc35, 'a3': a3, 'H': H
    })
    if trial < 8:
        print(f"  T{trial}: t3={t3:3d} t5={t5:4d} t7={t7:5d} bc={bc:3d} bc35={bc35:5d} a3={a3:3d} w2={w2:.1f}")

# =====================================================================
# Regression: w2 = c0 + c1*t3 + c2*t5 + c3*t7 + c4*bc + c5*bc35 + c6*a3
# =====================================================================
print(f"\n{'='*70}")
print("REGRESSION: w2 vs invariants")
print(f"{'='*70}")

# Build matrix
inv_names = ['const', 't3', 't5', 't7', 'bc', 'bc35', 'a3']
X = np.array([[1, d['t3'], d['t5'], d['t7'], d['bc'], d['bc35'], d['a3']] for d in data])
y = np.array([d['w2'] for d in data])

# Solve least squares
coeffs, residuals, rank, sv = np.linalg.lstsq(X, y, rcond=None)
y_pred = X @ coeffs
max_err = np.max(np.abs(y - y_pred))

print(f"\nw2 = ", end="")
terms = []
for i, name in enumerate(inv_names):
    if abs(coeffs[i]) > 0.01:
        terms.append(f"{coeffs[i]:.2f}*{name}")
print(" + ".join(terms))
print(f"Max residual: {max_err:.6f}")

if max_err < 0.5:
    print("EXACT FIT! These invariants fully determine w2.")
    # Round to exact integers/fractions
    print("\nExact coefficients:")
    for i, name in enumerate(inv_names):
        from fractions import Fraction
        frac = Fraction(coeffs[i]).limit_denominator(100)
        print(f"  {name}: {coeffs[i]:.6f} ≈ {frac}")
else:
    print(f"NOT exact — residual {max_err:.4f}. May need more invariants.")
    # Try adding t3^2 or other nonlinear terms
    print("\nTrying with additional terms: t3^2, t3*t5, t5^2...")
    X2 = np.column_stack([X, [d['t3']**2 for d in data], [d['t3']*d['t5'] for d in data]])
    coeffs2, _, _, _ = np.linalg.lstsq(X2, y, rcond=None)
    y_pred2 = X2 @ coeffs2
    max_err2 = np.max(np.abs(y - y_pred2))
    print(f"  With t3^2, t3*t5: max residual = {max_err2:.6f}")

# =====================================================================
# Verify match to THM-055 n=9 formula
# =====================================================================
print(f"\n{'='*70}")
print("COMPARISON TO THM-055 n=9 FORMULA")
print(f"{'='*70}")
print("""
THM-055 (n=9) claimed:
  c_2 = 462*t_3 - 60*t_5 + 12*t_7 - 120*bc33 + 24*bc35_w + 48*a3 - 2640
""")

# Check against THM-055 formula
for d in data[:10]:
    thm055 = 462*d['t3'] - 60*d['t5'] + 12*d['t7'] - 120*d['bc'] + 24*d['bc35'] + 48*d['a3'] - 2640
    print(f"  w2={d['w2']:.1f}, THM-055={thm055}, diff={abs(d['w2']-thm055):.2f}")

# =====================================================================
# Penalty at n=9
# =====================================================================
print(f"\n{'='*70}")
print("PENALTY H - w_0 AT n=9")
print(f"{'='*70}")

for d in data[:10]:
    w_coeffs_full = None
    # Recompute w0 for these samples
    random.seed(n*1000 + data.index(d))
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    wc = extract_W_coeffs(A, n)
    w0 = wc[0]
    penalty = d['H'] - w0
    print(f"  H={d['H']:5d} w0={w0:10.2f} penalty={penalty:10.2f}  "
          f"t5={d['t5']:4d} bc={d['bc']:3d} t7={d['t7']:5d} bc35={d['bc35']:5d} a3={d['a3']:3d}")

# Fit penalty as function of invariants
print(f"\n  Fitting penalty = c0 + c1*t5 + c2*bc + c3*t7 + c4*bc35 + c5*a3:")
X_pen = np.array([[1, d['t5'], d['bc'], d['t7'], d['bc35'], d['a3']] for d in data])
penalties = []
for idx, d in enumerate(data):
    random.seed(n*1000 + idx)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    wc = extract_W_coeffs(A, n)
    penalties.append(d['H'] - wc[0])

y_pen = np.array(penalties)
coeffs_pen, _, _, _ = np.linalg.lstsq(X_pen, y_pen, rcond=None)
y_pen_pred = X_pen @ coeffs_pen
max_err_pen = np.max(np.abs(y_pen - y_pen_pred))

pen_names = ['const', 't5', 'bc', 't7', 'bc35', 'a3']
print(f"\n  penalty = ", end="")
terms_pen = []
for i, name in enumerate(pen_names):
    if abs(coeffs_pen[i]) > 0.01:
        from fractions import Fraction
        frac = Fraction(coeffs_pen[i]).limit_denominator(1000)
        terms_pen.append(f"{frac}*{name}")
print(" + ".join(terms_pen))
print(f"  Max residual: {max_err_pen:.6f}")

print(f"\n{'='*70}")
print("DONE")
print(f"{'='*70}")
