#!/usr/bin/env python3
"""
W(r) vs tr(M(r)): are they the same polynomial?

W(r) = sum over all Ham paths P of product_{edges e in P} (r + s_e)
     where s_e = A[u,v] - 1/2 in {+1/2, -1/2}

tr(M(r)) = sum_a M[a,a](r) via the IE decomposition

At r=1/2: both equal H (the Hamiltonian path count).
For general r: are they equal?

If YES, then the coefficient analysis is much simpler:
  w_{n-1} = H (number of paths)
  w_{n-3} = sum_P e_2(s_P) = sum_P sum_{i<j} s_i * s_j
  and we can derive the tr(c_{n-3}) formula from combinatorics of edge signs.

opus-2026-03-06-S27
"""

from itertools import permutations, combinations
import numpy as np

def _count_weighted(A, verts, r_val, start=None, end=None):
    total = 0.0
    for p in permutations(verts):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        w = 1.0
        for i in range(len(p)-1):
            w *= r_val + (A[p[i]][p[i+1]] - 0.5)
        total += w
    return total

def transfer_trace_at_r(A, r_val):
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
                ea = _count_weighted(A, S_verts, r_val, end=a)
                ba = _count_weighted(A, R_verts, r_val, start=a)
                val += ((-1)**k) * ea * ba
        total += val
    return total

def weighted_path_sum(A, r_val):
    """W(r) = sum over all Ham paths of product(r + s_e)."""
    n = len(A)
    total = 0.0
    for p in permutations(range(n)):
        w = 1.0
        for i in range(n-1):
            w *= r_val + (A[p[i]][p[i+1]] - 0.5)
        total += w
    return total

# =====================================================================
import random
print("=" * 70)
print("W(r) vs tr(M(r)): ARE THEY THE SAME?")
print("=" * 70)

for n in [3, 4, 5]:
    print(f"\n  n={n}:")
    random.seed(n * 100)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    for r_val in [0.0, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0]:
        W = weighted_path_sum(A, r_val)
        trM = transfer_trace_at_r(A, r_val)
        match = "✓" if abs(W - trM) < 1e-8 else "✗"
        print(f"    r={r_val:.1f}: W={W:10.4f}, tr(M)={trM:10.4f}, diff={W-trM:+.2e} {match}")

# =====================================================================
# If they ARE the same, derive the coefficient structure
print("\n" + "=" * 70)
print("W(r) COEFFICIENT ANALYSIS (assuming W = tr(M))")
print("=" * 70)

# For a single Ham path (v_0, ..., v_{n-1}) with edge signs s_i:
# w(P,r) = prod_{i=0}^{n-2} (r + s_i) = sum_k e_{n-1-k}(s) * r^k
#
# e_0 = 1 (coefficient of r^{n-1})
# e_1 = sum_i s_i (coefficient of r^{n-2})
# e_2 = sum_{i<j} s_i*s_j (coefficient of r^{n-3})
#
# Now s_i = T(v_i, v_{i+1}) - 1/2.
# Let f_i = T(v_i, v_{i+1}) in {0,1}. Then s_i = f_i - 1/2.
#
# e_1(s) = sum(f_i - 1/2) = (sum f_i) - (n-1)/2 = f - (n-1)/2
# where f = number of forward edges.
#
# Sum over all paths:
# w_{n-2} = sum_P e_1(s_P) = sum_P [f_P - (n-1)/2]
#         = (total forward edges across all paths) - (n-1)/2 * H
#
# But by symmetry (M(r) = M(-r)), w_{n-2} = 0.
# This means: sum_P f_P = (n-1)/2 * H.
# = each path has (n-1)/2 forward edges on average!

n = 5
print(f"\n  n={n}: Verify w_{{n-2}} = 0")
random.seed(55)
for trial in range(3):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    total_forward = 0
    H = 0
    for p in permutations(range(n)):
        forward = sum(1 for i in range(n-1) if A[p[i]][p[i+1]] == 1)
        w = 1 if forward == n-1 else 0  # Directed path check
        total_forward += forward  # Over ALL permutations
        H += w

    # sum_P f_P over all Ham paths (not just directed ones)
    fwd_in_ham = 0
    ham_count = 0
    for p in permutations(range(n)):
        fwd = sum(1 for i in range(n-1) if A[p[i]][p[i+1]] == 1)
        # This path exists as a weighted path; w_{n-2} = sum_P e_1(s_P)
        # where s_i = A[p[i],p[i+1]] - 1/2
        e1 = sum(A[p[i]][p[i+1]] - 0.5 for i in range(n-1))
        fwd_in_ham += e1  # This should be 0

    print(f"    Trial {trial}: sum_P e_1(s_P) = {fwd_in_ham:.4f} (should be 0)")

# =====================================================================
# Derive w_{n-3} = sum_P e_2(s_P) in terms of tournament invariants
print(f"\n  n={n}: Analyze w_{{n-3}} = sum_P e_2(s_P)")

random.seed(42)
results = []
for trial in range(20):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    # Compute e_2 sum
    e2_sum = 0.0
    for p in permutations(range(n)):
        s = [A[p[i]][p[i+1]] - 0.5 for i in range(n-1)]
        # e_2 = sum_{i<j} s_i * s_j
        e2 = sum(s[i]*s[j] for i in range(n-1) for j in range(i+1, n-1))
        e2_sum += e2

    # Count 3-cycles
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                t3 += (A[i][j]*A[j][k]*A[k][i] + A[i][k]*A[k][j]*A[j][i])

    results.append((t3, e2_sum))

# Fit: e2_sum = a*t3 + b
X = np.array([[r[0], 1] for r in results])
y = np.array([r[1] for r in results])
coeffs = np.linalg.lstsq(X, y, rcond=None)[0]
max_err = max(abs(y - X @ coeffs))
print(f"\n    sum_P e_2(s_P) = {coeffs[0]:.4f}*t3 + {coeffs[1]:.4f}")
print(f"    Max error: {max_err:.6f}")
if max_err < 0.01:
    print(f"    EXACT! This is the algebraic basis for tr(c_{{n-3}}).")

# Now verify: w_{n-3} should equal tr(c_{n-3})
# w_{n-3} = sum_P e_2(s_P)
# tr(c_{n-3}) = 12*t_3 - 30 at n=5
print(f"\n    Expected tr(c_2) = 12*t3 - 30")
print(f"    Got w_{{n-3}} = {coeffs[0]:.4f}*t3 + {coeffs[1]:.4f}")
print(f"    Ratio of slopes: {coeffs[0]/12:.4f}")

# =====================================================================
# Repeat for n=7
print(f"\n  n=7: Compute w_{{n-3}} = sum_P e_2(s_P)")
print("  (Using DP for speed — count by edge contributions)")

# For n=7, we can't enumerate all 7!=5040 permutations efficiently
# but actually we can for w_{n-3}. Let me just do it.
n = 7
random.seed(77)
results7 = []
for trial in range(5):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    print(f"    Trial {trial}...", end="", flush=True)
    e2_sum = 0.0
    for p in permutations(range(n)):
        s = [A[p[i]][p[i+1]] - 0.5 for i in range(n-1)]
        e2 = sum(s[i]*s[j] for i in range(n-1) for j in range(i+1, n-1))
        e2_sum += e2

    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                t3 += (A[i][j]*A[j][k]*A[k][i] + A[i][k]*A[k][j]*A[j][i])

    results7.append((t3, e2_sum))
    print(f" t3={t3}, w_{{n-3}}={e2_sum:.1f}")

if len(results7) >= 2:
    X7 = np.array([[r[0], 1] for r in results7])
    y7 = np.array([r[1] for r in results7])
    c7 = np.linalg.lstsq(X7, y7, rcond=None)[0]
    e7 = max(abs(y7 - X7 @ c7))
    print(f"\n    sum_P e_2(s_P) = {c7[0]:.4f}*t3 + {c7[1]:.4f}")
    print(f"    Max error: {e7:.4f}")
    print(f"    Expected: 240*t3 - 2100")
    print(f"    Ratio: slope={c7[0]/240:.4f}, const_ratio={c7[1]/-2100:.4f}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
