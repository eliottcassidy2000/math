#!/usr/bin/env python3
"""
Test: which tr(M(r)) definition is correct?

Two definitions in the codebase:
1. transfer_trace_at_r (W_vs_trM_relationship.py):
   Both E_a and B_a use SAME tournament A → gives W(r) = sum_P prod(r+s)

2. compute_transfer_trace_bf (transfer_coeffs_general_n.py):
   Forward uses A[v][u], backward uses A[u][v] (converse!) → gives DIFFERENT polynomial

Hypothesis: Definition 1 is correct. Definition 2 has a bug (A[u][v] should be A[v][u]).

Test: compare both to the direct permutation sum at n=5 and n=7.

opus-2026-03-06-S29
"""
from itertools import permutations, combinations
import numpy as np
import random

def W_direct(A, n, r_val):
    """W(r) = sum_P prod(r + A[p_i,p_{i+1}] - 1/2)."""
    total = 0.0
    for P in permutations(range(n)):
        prod = 1.0
        for i in range(n-1):
            prod *= r_val + A[P[i]][P[i+1]] - 0.5
        total += prod
    return total

def trM_same_dir(A, n, r_val):
    """Transfer trace: both forward and backward use A[v][u]."""
    dp_fwd = {}
    for v in range(n):
        dp_fwd[(1 << v, v)] = 1.0
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            val = dp_fwd.get((mask, v), 0)
            if val == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                wt = r_val + (A[v][u] - 0.5)
                dp_fwd[(mask | (1 << u), u)] = dp_fwd.get((mask | (1 << u), u), 0) + val * wt

    dp_bwd = {}
    for v in range(n):
        dp_bwd[(1 << v, v)] = 1.0
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            val = dp_bwd.get((mask, v), 0)
            if val == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                wt = r_val + (A[v][u] - 0.5)  # SAME direction
                dp_bwd[(mask | (1 << u), u)] = dp_bwd.get((mask | (1 << u), u), 0) + val * wt

    full = (1 << n) - 1
    tr_val = 0.0
    for v in range(n):
        for mask_before in range(1 << n):
            if mask_before & (1 << v): continue
            mask_with_v = mask_before | (1 << v)
            fwd = dp_fwd.get((mask_with_v, v), 0)
            if fwd == 0: continue
            mask_after = full ^ mask_before
            if not (mask_after & (1 << v)): continue
            bwd = dp_bwd.get((mask_after, v), 0)
            if bwd == 0: continue
            k = bin(mask_before).count('1')
            tr_val += ((-1)**k) * fwd * bwd
    return tr_val

def trM_conv_dir(A, n, r_val):
    """Transfer trace: backward uses A[u][v] (converse) - the BUGGED version."""
    dp_fwd = {}
    for v in range(n):
        dp_fwd[(1 << v, v)] = 1.0
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            val = dp_fwd.get((mask, v), 0)
            if val == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                wt = r_val + (A[v][u] - 0.5)
                dp_fwd[(mask | (1 << u), u)] = dp_fwd.get((mask | (1 << u), u), 0) + val * wt

    dp_bwd = {}
    for v in range(n):
        dp_bwd[(1 << v, v)] = 1.0
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            val = dp_bwd.get((mask, v), 0)
            if val == 0: continue
            for u in range(n):
                if mask & (1 << u): continue
                wt = r_val + (A[u][v] - 0.5)  # CONVERSE direction
                dp_bwd[(mask | (1 << u), u)] = dp_bwd.get((mask | (1 << u), u), 0) + val * wt

    full = (1 << n) - 1
    tr_val = 0.0
    for v in range(n):
        for mask_before in range(1 << n):
            if mask_before & (1 << v): continue
            mask_with_v = mask_before | (1 << v)
            fwd = dp_fwd.get((mask_with_v, v), 0)
            if fwd == 0: continue
            mask_after = full ^ mask_before
            if not (mask_after & (1 << v)): continue
            bwd = dp_bwd.get((mask_after, v), 0)
            if bwd == 0: continue
            k = bin(mask_before).count('1')
            tr_val += ((-1)**k) * fwd * bwd
    return tr_val

# =====================================================================
print("=" * 70)
print("NORMALIZATION TEST: W(r) vs two tr(M(r)) definitions")
print("=" * 70)

for n in [5, 7]:
    print(f"\n  n={n}:")
    for trial in range(2):
        random.seed(n*1000 + trial)
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        print(f"\n  Trial {trial}:")
        for r in [0.0, 0.25, 0.5, 1.0]:
            w = W_direct(A, n, r)
            t1 = trM_same_dir(A, n, r)
            t2 = trM_conv_dir(A, n, r)
            print(f"    r={r:.2f}: W={w:12.4f}  trM_same={t1:12.4f}  trM_conv={t2:12.4f}")

# =====================================================================
print(f"\n{'='*70}")
print("POLYNOMIAL COEFFICIENTS")
print(f"{'='*70}")

for n in [5, 7]:
    num_coeffs = (n + 1) // 2  # even-power coefficients only
    r_sample = np.linspace(0, 1, num_coeffs + 1)[:num_coeffs]

    for trial in range(2):
        random.seed(n*1000 + trial)
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        W_vals = [W_direct(A, n, r) for r in r_sample]
        T2_vals = [trM_conv_dir(A, n, r) for r in r_sample]

        V = np.array([[r**(2*k) for k in range(num_coeffs)] for r in r_sample])
        w_coeffs = np.linalg.solve(V, W_vals)
        c_coeffs = np.linalg.solve(V, T2_vals)

        print(f"\n  n={n}, Trial {trial}:")
        for k in range(num_coeffs):
            diff = w_coeffs[k] - c_coeffs[k]
            print(f"    r^{2*k}: W={w_coeffs[k]:12.4f}  trM_conv={c_coeffs[k]:12.4f}  diff={diff:+12.4f}")

# =====================================================================
# Check what the "converse" version gives at r=1/2
print(f"\n{'='*70}")
print("H-count verification at r=1/2")
print(f"{'='*70}")

for n in [5, 7]:
    for trial in range(3):
        random.seed(n*1000 + trial)
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        # Count H directly
        H = 0
        for P in permutations(range(n)):
            if all(A[P[i]][P[i+1]] for i in range(n-1)):
                H += 1

        w = W_direct(A, n, 0.5)
        t2 = trM_conv_dir(A, n, 0.5)
        print(f"  n={n}, T{trial}: H={H}, W(1/2)={w:.1f}, trM_conv(1/2)={t2:.1f}")

print(f"\n{'='*70}")
print("DONE")
print(f"{'='*70}")
