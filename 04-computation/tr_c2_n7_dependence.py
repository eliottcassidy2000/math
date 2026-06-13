#!/usr/bin/env python3
"""
What does tr(c_2) depend on at n=7?

We found: tr(c_4) = 240*t_3 - 2100 (exact, depends only on t_3).
But tr(c_2) varies even within the same t_3. What other invariant
does it depend on?

Candidates:
- Score sequence (sum of squares of scores, etc.)
- Number of 5-cycles (t_5)
- "Landau excess" sum s_i^2 - n(n-1)^2/4
- Path counts for sub-tournaments

opus-2026-03-06-S26
"""

from itertools import permutations, combinations
import numpy as np
import random

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

def transfer_trace_r(A, r_val):
    n = len(A)
    total = 0.0
    for a in range(n):
        U = [v for v in range(n) if v != a]
        val = 0.0
        for k in range(len(U)+1):
            for S in combinations(U, k):
                S_set = set(S)
                R = [v for v in U if v not in S_set]
                S_verts = sorted(list(S) + [a])
                R_verts = sorted(R + [a])
                ea = count_paths_weighted(A, S_verts, r_val, end=a)
                bb = count_paths_weighted(A, R_verts, r_val, start=b)
                val += ((-1)**k) * ea * bb
        total += val
    return total

def transfer_trace_r_v2(A, r_val):
    """Corrected: for diagonal M[a,a], both subpaths start/end at a."""
    n = len(A)
    total = 0.0
    for a in range(n):
        U = [v for v in range(n) if v != a]
        val = 0.0
        for k in range(len(U)+1):
            for S in combinations(U, k):
                S_set = set(S)
                R = [v for v in U if v not in S_set]
                S_verts = sorted(list(S) + [a])
                R_verts = sorted(R + [a])
                ea = count_paths_weighted(A, S_verts, r_val, end=a)
                ba = count_paths_weighted(A, R_verts, r_val, start=a)
                val += ((-1)**k) * ea * ba
        total += val
    return total

def ham_path_count(A):
    n = len(A)
    return sum(1 for p in permutations(range(n))
               if all(A[p[i]][p[i+1]] == 1 for i in range(n-1)))

def count_3cycles(A):
    n = len(A)
    count = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                count += (A[i][j]*A[j][k]*A[k][i] + A[i][k]*A[k][j]*A[j][i])
    return count

def count_5cycles(A):
    n = len(A)
    count = 0
    for perm in permutations(range(n), 5):
        if all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
            count += 1
    return count // 5

def score_sq_sum(A):
    n = len(A)
    scores = [sum(A[i]) for i in range(n)]
    return sum(s*s for s in scores)

# =====================================================================
n = 7
print("=" * 70)
print("n=7: tr(c_2) DEPENDENCE")
print("=" * 70)

u_samples = np.array([0.0, 0.04, 0.16, 0.36])
random.seed(17)

data = []
for trial in range(20):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    H = ham_path_count(A)
    t3 = count_3cycles(A)
    t5 = count_5cycles(A)
    sq = score_sq_sum(A)
    scores = tuple(sorted(sum(A[i]) for i in range(n)))

    trace_vals = [transfer_trace_r_v2(A, np.sqrt(u)) for u in u_samples]
    poly = np.polyfit(u_samples, trace_vals, 3)
    tr_c6 = poly[0]
    tr_c4 = poly[1]
    tr_c2 = poly[2]
    tr_c0 = poly[3]

    data.append({
        'H': H, 't3': t3, 't5': t5, 'sq': sq, 'scores': scores,
        'tr_c0': tr_c0, 'tr_c2': tr_c2, 'tr_c4': tr_c4, 'tr_c6': tr_c6
    })

# Fit tr(c_2) = a*t3 + b*t5 + c
print("\n  Fitting tr(c_2) = a*t3 + b*t5 + c:")
X = np.array([[d['t3'], d['t5'], 1] for d in data])
y = np.array([d['tr_c2'] for d in data])
coeffs, residuals, _, _ = np.linalg.lstsq(X, y, rcond=None)
print(f"    tr(c_2) = {coeffs[0]:.4f}*t3 + {coeffs[1]:.4f}*t5 + {coeffs[2]:.4f}")
if len(residuals) > 0:
    print(f"    Residual: {residuals[0]:.6f}")

# Check predictions
max_err = 0
for d in data:
    pred = coeffs[0]*d['t3'] + coeffs[1]*d['t5'] + coeffs[2]
    err = abs(d['tr_c2'] - pred)
    max_err = max(max_err, err)
print(f"    Max error: {max_err:.4f}")

# Also try: tr(c_2) = a*t3 + b*sq + c
print("\n  Fitting tr(c_2) = a*t3 + b*sum(s_i^2) + c:")
X2 = np.array([[d['t3'], d['sq'], 1] for d in data])
coeffs2, res2, _, _ = np.linalg.lstsq(X2, y, rcond=None)
print(f"    tr(c_2) = {coeffs2[0]:.4f}*t3 + {coeffs2[1]:.4f}*sq + {coeffs2[2]:.4f}")
max_err2 = max(abs(d['tr_c2'] - coeffs2[0]*d['t3'] - coeffs2[1]*d['sq'] - coeffs2[2]) for d in data)
print(f"    Max error: {max_err2:.4f}")

# Note: t3 and sq are related! t3 = C(n,3) - sum C(s_i, 2) = C(n,3) - (sq - n*(n-1)/2)/2
# So t3 = C(n,3) - sq/2 + n*(n-1)/4
# At n=7: t3 = 35 - sq/2 + 7*6/4 = 35 - sq/2 + 10.5 = 45.5 - sq/2
print(f"\n  Relation: t3 = C(7,3) - sum(C(s_i,2)) = 35 - (sq - 21)/2 = 45.5 - sq/2")
for d in data[:5]:
    pred_t3 = 45.5 - d['sq']/2
    print(f"    sq={d['sq']}, t3={d['t3']}, predicted={pred_t3:.1f}")

# Since t3 and sq are linearly related, try: tr(c_2) = a*t5 + b*sq + c
print("\n  Fitting tr(c_2) = a*t5 + b*sq + c (sq encodes t3):")
X3 = np.array([[d['t5'], d['sq'], 1] for d in data])
coeffs3, res3, _, _ = np.linalg.lstsq(X3, y, rcond=None)
print(f"    tr(c_2) = {coeffs3[0]:.4f}*t5 + {coeffs3[1]:.4f}*sq + {coeffs3[2]:.4f}")
max_err3 = max(abs(d['tr_c2'] - coeffs3[0]*d['t5'] - coeffs3[1]*d['sq'] - coeffs3[2]) for d in data)
print(f"    Max error: {max_err3:.4f}")

# Print data table
print(f"\n  {'Scores':<25s} | H  | t3 | t5 | sq | tr_c2  | pred_t3t5 | err")
print("  " + "-" * 90)
for d in sorted(data, key=lambda x: x['t3']):
    pred = coeffs[0]*d['t3'] + coeffs[1]*d['t5'] + coeffs[2]
    err = d['tr_c2'] - pred
    print(f"  {str(d['scores']):<25s} | {d['H']:3d} | {d['t3']:2d} | {d['t5']:2d} | {d['sq']:2d} | {d['tr_c2']:7.2f} | {pred:9.2f} | {err:+.2f}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
