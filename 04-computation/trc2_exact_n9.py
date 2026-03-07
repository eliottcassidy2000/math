#!/usr/bin/env python3
"""Verify the c_2 formula at n=9 with more tournaments."""
from itertools import combinations
from math import factorial, comb
import numpy as np
import random
import time

n = 9

def compute_hc_on_subset(A, verts):
    k = len(verts)
    sub = [[A[verts[i]][verts[j]] for j in range(k)] for i in range(k)]
    dp = [[0]*k for _ in range(1 << k)]
    dp[1][0] = 1
    for mask in range(1, 1 << k):
        for v in range(k):
            if not (mask & (1 << v)) or dp[mask][v] == 0: continue
            for u in range(k):
                if mask & (1 << u): continue
                if sub[v][u]:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    full = (1 << k) - 1
    return sum(dp[full][v] for v in range(1, k) if sub[v][0])

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

def compute_invariants(A, n):
    # t_k: directed k-cycle counts
    cycle_data = {}
    for k in [3, 5, 7, 9]:
        cycle_data[k] = []
        for verts in combinations(range(n), k):
            hc = compute_hc_on_subset(A, verts)
            if hc > 0:
                cycle_data[k].append((frozenset(verts), hc))

    t3 = sum(hc for (vs, hc) in cycle_data[3])
    t5 = sum(hc for (vs, hc) in cycle_data[5])
    t7 = sum(hc for (vs, hc) in cycle_data[7])
    t9 = sum(hc for (vs, hc) in cycle_data[9])

    # bc33_w: pairs of disjoint 3-cycles (hc=1 each, so unweighted)
    bc33 = 0
    sets3 = cycle_data[3]
    for i in range(len(sets3)):
        for j in range(i+1, len(sets3)):
            if sets3[i][0].isdisjoint(sets3[j][0]):
                bc33 += 1

    # bc35_w: pairs (3-cycle, 5-cycle set) disjoint, weighted by 5-cycle hc
    bc35_w = 0
    for (vs3, hc3) in cycle_data[3]:
        for (vs5, hc5) in cycle_data[5]:
            if vs3.isdisjoint(vs5):
                bc35_w += hc3 * hc5  # hc3=1 always

    # alpha_3: triples of disjoint 3-cycles
    a3 = 0
    for i in range(len(sets3)):
        for j in range(i+1, len(sets3)):
            if not sets3[i][0].isdisjoint(sets3[j][0]): continue
            for ki in range(j+1, len(sets3)):
                if sets3[ki][0].isdisjoint(sets3[i][0]) and sets3[ki][0].isdisjoint(sets3[j][0]):
                    a3 += 1

    return t3, t5, t7, t9, bc33, bc35_w, a3

def compute_transfer_coeffs(A, n, r_values):
    traces = []
    full = (1 << n) - 1
    for r in r_values:
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
    return traces

# Formula: c_2 = 462*t3 - 60*t5 + 12*t7 - 120*bc33 + 24*bc35_w + 48*a3 - 2640
r_values = [0.0, 0.25, 0.5, 0.75, 1.0]

print("Verifying c_2 = 462*t3 - 60*t5 + 12*t7 - 120*bc33 + 24*bc35_w + 48*a3 - 2640")
print("at n=9 with 30 random tournaments...")

max_err = 0
for trial in range(30):
    random.seed(trial * 47 + 13)
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1

    t3, t5, t7, t9, bc33, bc35_w, a3 = compute_invariants(A, n)

    traces = compute_transfer_coeffs(A, n, r_values)
    r_powers = np.array([[r**(2*k) for k in range(5)] for r in r_values])
    coeffs = np.linalg.solve(r_powers, traces)
    c2_actual = coeffs[1]

    c2_formula = 462*t3 - 60*t5 + 12*t7 - 120*bc33 + 24*bc35_w + 48*a3 - 2640
    err = abs(c2_actual - c2_formula)
    max_err = max(max_err, err)

    if trial < 10 or err > 1:
        print(f"  T{trial:2d}: c2={c2_actual:10.1f} formula={c2_formula:10d} err={err:.2f}")

print(f"\nMax error: {max_err:.4f}")
print(f"Formula {'VERIFIED' if max_err < 1 else 'FAILED'}")

print("\nDONE")
