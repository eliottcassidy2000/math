#!/usr/bin/env python3
"""
Why are tr(c_2) and tr(c_4) universal for regular tournaments?

At n=5 (regular): tr(c_2) = 12*5 - 30 = 30, tr(c_4) = 120 = 5!
At n=7 (regular): tr(c_2) = 63, tr(c_4) = 1260, tr(c_6) = 5040 = 7!

Hypothesis: For regular tournaments (score = (n-1)/2 everywhere),
the coefficients tr(c_k) for k >= 2 depend only on n, not on the
specific tournament.

Since t_3 is determined by the score sequence (t_3 = C(n,3) - sum C(s_i,2)),
and regular tournaments have s_i = (n-1)/2 for all i, we get:
t_3 = C(n,3) - n*C((n-1)/2, 2)

For n=5: t_3 = 10 - 5*1 = 5
For n=7: t_3 = 35 - 7*3 = 14
For n=9: t_3 = 84 - 9*6 = 30

This fixes tr(c_2) if tr(c_2) depends only on t_3 and n. But at n=7
ALL coefficients tr(c_2), tr(c_4) are universal — do they depend on
higher cycle counts too?

Check: for non-regular tournaments at n=7, does tr(c_2) vary with t_3?

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
                bb = count_paths_weighted(A, R_verts, r_val, start=a)
                val += ((-1)**k) * ea * bb
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

# =====================================================================
print("=" * 70)
print("n=5: NON-REGULAR tr(c_2) FORMULA")
print("=" * 70)

n = 5
u_samples = np.array([0.0, 0.04, 0.16])

# Sample various tournaments
random.seed(123)
data_n5 = []

# Some specific tournaments
specific = [
    [[0,1,1,1,1],[0,0,1,1,1],[0,0,0,1,1],[0,0,0,0,1],[0,0,0,0,0]],  # transitive
    [[0,1,1,0,0],[0,0,1,0,0],[0,0,0,1,0],[1,1,0,0,0],[1,1,1,1,0]],  # various
]

for _ in range(15):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    specific.append(A)

for A in specific:
    H = ham_path_count(A)
    t3 = count_3cycles(A)
    scores = tuple(sorted(sum(row) for row in A))

    trace_vals = [transfer_trace_r(A, np.sqrt(u)) for u in u_samples]
    poly = np.polyfit(u_samples, trace_vals, 2)
    tr_c4 = poly[0]
    tr_c2 = poly[1]
    tr_c0 = poly[2]

    # Check tr(c_2) = 12*t3 - 30
    expected_c2 = 12*t3 - 30
    match = abs(tr_c2 - expected_c2) < 0.5

    data_n5.append((scores, H, t3, tr_c0, tr_c2, tr_c4, expected_c2, match))

print(f"\n  {'Scores':<20s} | H | t3 | tr(c2) | 12t3-30 | tr(c4) | match")
print("  " + "-" * 70)
for scores, H, t3, c0, c2, c4, exp, m in sorted(data_n5):
    print(f"  {str(scores):<20s} | {H:2d} | {t3:1d} | {c2:6.1f} | {exp:7.1f} | {c4:5.0f} | {m}")

# =====================================================================
print("\n" + "=" * 70)
print("n=7: NON-REGULAR TOURNAMENTS (SAMPLED)")
print("=" * 70)

n = 7
u_samples_7 = np.array([0.0, 0.04, 0.16, 0.36])

random.seed(42)
data_n7 = []
for trial in range(8):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    H = ham_path_count(A)
    t3 = count_3cycles(A)
    scores = tuple(sorted(sum(row) for row in A))

    trace_vals = [transfer_trace_r(A, np.sqrt(u)) for u in u_samples_7]
    poly = np.polyfit(u_samples_7, trace_vals, 3)
    tr_c6 = poly[0]
    tr_c4 = poly[1]
    tr_c2 = poly[2]
    tr_c0 = poly[3]

    data_n7.append((scores, H, t3, tr_c0, tr_c2, tr_c4, tr_c6))

print(f"\n  {'Scores':<25s} | H  | t3 | tr(c0) | tr(c2) | tr(c4) | tr(c6)")
print("  " + "-" * 85)
for scores, H, t3, c0, c2, c4, c6 in sorted(data_n7):
    print(f"  {str(scores):<25s} | {H:3d} | {t3:2d} | {c0:6.2f} | {c2:7.2f} | {c4:7.1f} | {c6:7.1f}")

# Check: does tr(c_2) = f(t3) at n=7?
print("\n  tr(c_2) vs t3 at n=7:")
for scores, H, t3, c0, c2, c4, c6 in sorted(data_n7, key=lambda x: x[2]):
    print(f"    t3={t3:2d}, tr(c_2)={c2:7.2f}, scores={scores}")

# Try linear fit: tr(c_2) = a*t3 + b
X = np.array([[d[2], 1] for d in data_n7])
y = np.array([d[4] for d in data_n7])
coeffs = np.linalg.lstsq(X, y, rcond=None)[0]
print(f"\n  Best fit: tr(c_2) = {coeffs[0]:.4f}*t3 + {coeffs[1]:.4f}")
for scores, H, t3, c0, c2, c4, c6 in sorted(data_n7, key=lambda x: x[2]):
    predicted = coeffs[0]*t3 + coeffs[1]
    print(f"    t3={t3:2d}, tr(c_2)={c2:7.2f}, predicted={predicted:7.2f}, err={abs(c2-predicted):.2f}")

# Similarly for tr(c_4)
print("\n  tr(c_4) vs t3 at n=7:")
y4 = np.array([d[5] for d in data_n7])
coeffs4 = np.linalg.lstsq(X, y4, rcond=None)[0]
print(f"  Best fit: tr(c_4) = {coeffs4[0]:.4f}*t3 + {coeffs4[1]:.4f}")
for scores, H, t3, c0, c2, c4, c6 in sorted(data_n7, key=lambda x: x[2]):
    predicted = coeffs4[0]*t3 + coeffs4[1]
    print(f"    t3={t3:2d}, tr(c_4)={c4:7.1f}, predicted={predicted:7.1f}, err={abs(c4-predicted):.1f}")

print("\n" + "=" * 70)
print("THEORETICAL CHECK")
print("=" * 70)
print(f"""
  At n=5: tr(c_2) = 12*t_3 - 30
    Coefficient of t_3: 12
    Constant: -30

  Can we derive 12 from first principles?
  The coefficient c_2 comes from the r^2 term in the transfer matrix.
  In the IE formula, each edge weight is (r + s_e) where s_e = +1/2 or -1/2.
  The r^2 term collects all products with exactly 2 factors of r.

  For the TRACE (diagonal only):
  M[a,a](r) = sum over (S,R) partitions of V-a:
    (-1)^|S| * (paths ending at a in S+a) * (paths starting at a in R+a)
  The r^2 contribution picks 2 edges among the sub-paths to contribute r.

  For a regular tournament: s_i + s_j = 0 for complementary edge pairs.
  The sum over all edge-pair choices might telescope due to regularity.
""")

print("=" * 70)
print("DONE")
print("=" * 70)
