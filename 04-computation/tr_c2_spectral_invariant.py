#!/usr/bin/env python3
"""
What invariant determines tr(c_2) at n=7?

We know:
  - tr(c_4) = 240*t_3 - 2100 (EXACT, depends only on t_3)
  - tr(c_6) = 720 (universal)
  - tr(c_2) varies even within same score sequence
  - tr(c_0) = H - tr(c_2)/4 - tr(c_4)/16 - tr(c_6)/64

So understanding tr(c_2) would complete our understanding of all coefficients.

STRATEGY: Compute tr(c_2) for many n=7 tournaments and correlate with:
  1. Adjacency matrix eigenvalues
  2. Number of 5-cycles (t_5)
  3. Number of 7-cycles (t_7)
  4. Number of directed paths of length k (for various k)
  5. "Imbalance" invariants (sum of s_i^3, s_i^4, etc.)
  6. The skew-adjacency matrix spectrum

opus-2026-03-06-S27
"""

from itertools import permutations, combinations
import numpy as np
from math import factorial

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

def count_5cycles(A):
    """Count directed 5-cycles."""
    n = len(A)
    count = 0
    for perm in permutations(range(n), 5):
        if all(A[perm[i]][perm[(i+1)%5]] for i in range(5)):
            count += 1
    return count // 5  # Each 5-cycle counted 5 times (cyclic)

def transfer_trace_at_r(A, r_val):
    """Compute tr(M(r)) via IE formula for diagonal entries."""
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

def weighted_path_sum(A, r_val):
    """W(r) = sum over all Ham paths of product(r + s_e).
    This is NOT the same as tr(M(r)) in general!"""
    n = len(A)
    total = 0.0
    for p in permutations(range(n)):
        w = 1.0
        for i in range(n-1):
            w *= r_val + (A[p[i]][p[i+1]] - 0.5)
        total += w
    return total

def skew_adj_eigenvalues(A):
    """Eigenvalues of skew-adjacency matrix S[i,j] = A[i,j] - A[j,i]."""
    n = len(A)
    S = np.array([[A[i][j] - A[j][i] for j in range(n)] for i in range(n)], dtype=float)
    return sorted(np.linalg.eigvals(S).imag, reverse=True)

def directed_path_count(A, length):
    """Count directed paths of given length (not necessarily Hamiltonian)."""
    n = len(A)
    An = np.matrix(A, dtype=float)
    power = np.linalg.matrix_power(An, length)
    return int(np.trace(power))  # Number of closed walks

# =====================================================================
import random
n = 7
print("=" * 70)
print(f"n={n}: WHAT DETERMINES tr(c_2)?")
print("=" * 70)

# Collect data
random.seed(2026)
u_samples = np.array([0.0, 0.04, 0.16, 0.36])

data = []
count = 0
while count < 40:
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    H = ham_path_count_dp(A)
    t3 = count_3cycles(A)
    scores = tuple(sorted(sum(A[i]) for i in range(n)))

    # Compute even-r polynomial coefficients via sampling
    trace_vals = [transfer_trace_at_r(A, np.sqrt(u)) for u in u_samples]
    poly = np.polyfit(u_samples, trace_vals, 3)
    tr_c6 = poly[0]
    tr_c4 = poly[1]
    tr_c2 = poly[2]
    tr_c0 = poly[3]

    # Spectral invariants
    An = np.array(A, dtype=float)
    eigs = sorted(np.linalg.eigvals(An - An.T).imag, reverse=True)

    # Score-based invariants
    score_list = [sum(A[i]) for i in range(n)]
    sq_sum = sum(s**2 for s in score_list)
    cube_sum = sum(s**3 for s in score_list)
    fourth_sum = sum(s**4 for s in score_list)

    # Closed walk counts
    A_mat = np.array(A, dtype=float)
    cw2 = int(np.trace(A_mat @ A_mat))  # = sum of s_i (= edges)
    cw3 = int(round(np.trace(np.linalg.matrix_power(A_mat, 3))))  # Related to 3-cycles
    cw4 = int(round(np.trace(np.linalg.matrix_power(A_mat, 4))))
    cw5 = int(round(np.trace(np.linalg.matrix_power(A_mat, 5))))

    data.append({
        'H': H, 't3': t3, 'scores': scores,
        'tr_c0': tr_c0, 'tr_c2': tr_c2, 'tr_c4': tr_c4, 'tr_c6': tr_c6,
        'sq': sq_sum, 'cube': cube_sum, 'fourth': fourth_sum,
        'eigs': eigs,
        'cw2': cw2, 'cw3': cw3, 'cw4': cw4, 'cw5': cw5,
    })

    count += 1
    if count % 10 == 0:
        print(f"  Computed {count}/40...", flush=True)

# Also compute W(r) vs tr(M(r)) for a few tournaments
print("\n--- W(r) vs tr(M(r)) comparison ---")
for i in range(3):
    A = [[0]*n for _ in range(n)]
    random.seed(100 + i)
    for ii in range(n):
        for jj in range(ii+1, n):
            if random.random() < 0.5:
                A[ii][jj] = 1
            else:
                A[jj][ii] = 1

    for r_test in [0.0, 0.2, 0.5]:
        W = weighted_path_sum(A, r_test)
        trM = transfer_trace_at_r(A, r_test)
        print(f"  T{i}, r={r_test}: W(r)={W:.4f}, tr(M(r))={trM:.4f}, diff={W-trM:.6f}")

# =====================================================================
print("\n" + "=" * 70)
print("CORRELATION ANALYSIS")
print("=" * 70)

y = np.array([d['tr_c2'] for d in data])

# Test 1: tr(c_2) = a*t3 + b (score-independent part?)
X1 = np.array([[d['t3'], 1] for d in data])
c1, res1, _, _ = np.linalg.lstsq(X1, y, rcond=None)
err1 = max(abs(y - X1 @ c1))
print(f"\n  tr(c_2) ~ a*t3 + b: max_err={err1:.2f}")

# Test 2: tr(c_2) = a*t3 + b*sq + c
X2 = np.array([[d['t3'], d['sq'], 1] for d in data])
c2, res2, _, _ = np.linalg.lstsq(X2, y, rcond=None)
err2 = max(abs(y - X2 @ c2))
print(f"  tr(c_2) ~ a*t3 + b*sq + c: max_err={err2:.2f}")

# Test 3: tr(c_2) = a*cw3 + b*cw4 + c*cw5 + d
X3 = np.array([[d['cw3'], d['cw4'], d['cw5'], 1] for d in data])
c3, res3, _, _ = np.linalg.lstsq(X3, y, rcond=None)
err3 = max(abs(y - X3 @ c3))
print(f"  tr(c_2) ~ a*cw3 + b*cw4 + c*cw5 + d: max_err={err3:.2f}")

# Test 4: tr(c_2) = a*cw3 + b*cw5 + c
X4 = np.array([[d['cw3'], d['cw5'], 1] for d in data])
c4, res4, _, _ = np.linalg.lstsq(X4, y, rcond=None)
err4 = max(abs(y - X4 @ c4))
print(f"  tr(c_2) ~ a*cw3 + b*cw5 + c: max_err={err4:.2f}")

# Test 5: tr(c_2) = a*sq + b*cube + c*fourth + d
X5 = np.array([[d['sq'], d['cube'], d['fourth'], 1] for d in data])
c5, res5, _, _ = np.linalg.lstsq(X5, y, rcond=None)
err5 = max(abs(y - X5 @ c5))
print(f"  tr(c_2) ~ a*sq + b*cube + c*fourth + d: max_err={err5:.2f}")

# Test 6: Include sum of eig_i^4 (4th moment of skew-adjacency spectrum)
eig_moments = []
for d in data:
    e = d['eigs']
    eig_moments.append([sum(x**2 for x in e), sum(x**4 for x in e)])
X6 = np.array([[d['cw3'], eig_moments[i][0], eig_moments[i][1], 1]
                for i, d in enumerate(data)])
c6, res6, _, _ = np.linalg.lstsq(X6, y, rcond=None)
err6 = max(abs(y - X6 @ c6))
print(f"  tr(c_2) ~ a*cw3 + b*eig2 + c*eig4 + d: max_err={err6:.2f}")

# Test 7: ALL closed walks up to length 5
X7 = np.array([[d['cw2'], d['cw3'], d['cw4'], d['cw5'], 1] for d in data])
c7, res7, _, _ = np.linalg.lstsq(X7, y, rcond=None)
err7 = max(abs(y - X7 @ c7))
print(f"  tr(c_2) ~ cw2+cw3+cw4+cw5: max_err={err7:.2f}")
if err7 < 0.5:
    print(f"    EXACT! coeffs = {c7}")

# =====================================================================
# Check if it's related to H itself
print("\n--- Checking correlation with H ---")
X_H = np.array([[d['H'], d['t3'], 1] for d in data])
c_H, _, _, _ = np.linalg.lstsq(X_H, y, rcond=None)
err_H = max(abs(y - X_H @ c_H))
print(f"  tr(c_2) ~ a*H + b*t3 + c: max_err={err_H:.2f}")

# tr(c_0) = H - tr(c_2)/4 - tr(c_4)/16 - tr(c_6)/64, so:
# H = tr(c_0) + tr(c_2)/4 + tr(c_4)/16 + tr(c_6)/64
# => tr(c_2) = 4*H - 4*tr(c_0) - tr(c_4)/4 - tr(c_6)/16
# So if we know H AND tr(c_0), we know tr(c_2). But that's circular.

# =====================================================================
# Print detailed table for score class (3,3,3,3,3,3,3) — regular
print("\n--- REGULAR tournaments (score = (3,3,3,3,3,3,3)) ---")
reg_data = [d for d in data if d['scores'] == (3,3,3,3,3,3,3)]
print(f"  Found {len(reg_data)} regular tournaments")
for d in sorted(reg_data, key=lambda x: x['H']):
    print(f"  H={d['H']:3d} t3={d['t3']:2d} cw3={d['cw3']:3d} cw5={d['cw5']:4d} tr(c_2)={d['tr_c2']:7.1f}")

# =====================================================================
# Check: is tr(c_2) determined by cw3 AND cw5?
print("\n--- Same cw3 & cw5 but different tr(c_2)? ---")
from collections import defaultdict
cw_groups = defaultdict(list)
for d in data:
    cw_groups[(d['cw3'], d['cw5'])].append(d)
for key, group in sorted(cw_groups.items()):
    if len(group) > 1:
        trc2_vals = [d['tr_c2'] for d in group]
        if max(trc2_vals) - min(trc2_vals) > 1:
            print(f"  cw3={key[0]}, cw5={key[1]}: tr(c_2) = {[f'{v:.1f}' for v in trc2_vals]}")

print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
