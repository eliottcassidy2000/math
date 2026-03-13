#!/usr/bin/env python3
"""
regular_n7_classification.py — opus-2026-03-13-S67j

INVESTIGATING: Regular tournaments at n=7 have H ∈ {171, 175, 189}.
What distinguishes these three classes?

THM-164 predicts: degree-4 (5-cycles) and degree-6 (7-cycles) Fourier terms
create the discrimination. We compute:
1. Number of directed 3-cycles, 5-cycles, 7-cycles for each class
2. Whether Paley P_7 is special among regular tournaments
3. The automorphism group structure
"""

import numpy as np
import math
import random
import time

def ham_path_count_dp(A):
    n = A.shape[0]
    dp = np.zeros((1 << n, n), dtype=np.int64)
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if dp[mask][v] == 0 or not (mask & (1 << v)):
                continue
            for u in range(n):
                if not (mask & (1 << u)) and A[v][u] == 1:
                    dp[mask | (1 << u)][u] += dp[mask][v]
    return int(np.sum(dp[(1 << n) - 1]))

def count_directed_cycles(A, length):
    """Count directed cycles of given length in tournament A."""
    n = A.shape[0]
    if length == 3:
        # c3 = number of directed 3-cycles
        count = 0
        for i in range(n):
            for j in range(i+1, n):
                for k in range(j+1, n):
                    # Check all 2 orientations of the cycle
                    if A[i][j] and A[j][k] and A[k][i]:
                        count += 1
                    if A[i][k] and A[k][j] and A[j][i]:
                        count += 1
        return count
    elif length == 5:
        # Count 5-cycles: i->j->k->l->m->i
        count = 0
        from itertools import permutations
        for perm in permutations(range(n), length):
            valid = True
            for idx in range(length):
                if A[perm[idx]][perm[(idx+1) % length]] != 1:
                    valid = False
                    break
            if valid:
                count += 1
        return count // length  # each cycle counted 'length' times (rotations)
    elif length == 7:
        # For n=7, a 7-cycle uses all vertices = Hamiltonian cycle
        count = 0
        from itertools import permutations
        for perm in permutations(range(1, n)):
            full = (0,) + perm
            valid = True
            for idx in range(n):
                if A[full[idx]][full[(idx+1) % n]] != 1:
                    valid = False
                    break
            if valid:
                count += 1
        return count // n  # each cycle counted n times
    return 0

def score_seq(A):
    return tuple(sorted(A.sum(axis=1).astype(int)))

# =====================================================================
# Generate all regular tournaments at n=7 by sampling
# =====================================================================
print("=" * 70)
print("REGULAR TOURNAMENTS AT n=7: THREE H CLASSES")
print("=" * 70)

n = 7
edges = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(edges)

random.seed(42)

# Collect regular tournaments by sampling
print(f"  Sampling regular tournaments from 2^{m} = {2**m} total...")

regular_by_H = {}
n_samples = 100000
t0 = time.time()

for trial in range(n_samples):
    bits = random.randint(0, 2**m - 1)
    A = np.zeros((n,n), dtype=np.int8)
    for k, (i,j) in enumerate(edges):
        if (bits >> k) & 1:
            A[i][j] = 1
        else:
            A[j][i] = 1

    # Check if regular
    scores = A.sum(axis=1)
    if not np.all(scores == 3):
        continue

    h = ham_path_count_dp(A)
    if h not in regular_by_H:
        regular_by_H[h] = []
    regular_by_H[h].append((bits, A.copy()))

elapsed = time.time() - t0
total_found = sum(len(v) for v in regular_by_H.values())
print(f"  Found {total_found} regular tournaments in {elapsed:.1f}s")
print(f"  H classes: {sorted(regular_by_H.keys())}")
for h_val in sorted(regular_by_H.keys()):
    print(f"    H={h_val}: {len(regular_by_H[h_val])} found")

# =====================================================================
# Analyze cycle structure of each class
# =====================================================================
print("\n" + "=" * 70)
print("CYCLE STRUCTURE BY H CLASS")
print("=" * 70)

for h_val in sorted(regular_by_H.keys()):
    group = regular_by_H[h_val]
    # Analyze first 10 representatives
    sample_size = min(10, len(group))

    c3_vals = []
    c5_vals = []
    c7_vals = []

    for bits, A in group[:sample_size]:
        c3 = count_directed_cycles(A, 3)
        c5 = count_directed_cycles(A, 5)
        c7 = count_directed_cycles(A, 7)
        c3_vals.append(c3)
        c5_vals.append(c5)
        c7_vals.append(c7)

    print(f"\n  H={h_val} (sample of {sample_size}):")
    print(f"    c3 (3-cycles): {sorted(set(c3_vals))} (all same? {len(set(c3_vals))==1})")
    print(f"    c5 (5-cycles): {sorted(set(c5_vals))} (all same? {len(set(c5_vals))==1})")
    print(f"    c7 (7-cycles): {sorted(set(c7_vals))} (all same? {len(set(c7_vals))==1})")

    # For regular tournaments at n=7: c3 should be constant
    # c3 = n(n-1)(n-3)/24 for regular tournaments (Kendall-Babington Smith)
    expected_c3 = n * (n-1) * (n+1) // 24
    print(f"    Expected c3 for regular n=7: {expected_c3}")

    # Also compute det(I+A) and AA^T diagonal
    for bits, A in group[:1]:
        det_val = abs(np.linalg.det(np.eye(n) + A.astype(float)))
        AAT = A.astype(float) @ A.astype(float).T
        diag_AAT = np.diag(AAT)
        off_diag_var = np.var(AAT[np.triu_indices(n, 1)])
        print(f"    det(I+A) = {det_val:.2f}")
        print(f"    AA^T diagonal = {diag_AAT}")
        print(f"    AA^T off-diagonal variance = {off_diag_var:.4f}")

# =====================================================================
# Check Paley P_7 specifically
# =====================================================================
print("\n" + "=" * 70)
print("PALEY P_7 STRUCTURE")
print("=" * 70)

QR = set()
for k in range(1, 7):
    QR.add((k * k) % 7)
print(f"  QR mod 7 = {sorted(QR)}")

A_paley = np.zeros((7, 7), dtype=np.int8)
for i in range(7):
    for j in range(7):
        if i != j and ((j - i) % 7) in QR:
            A_paley[i][j] = 1

h_paley = ham_path_count_dp(A_paley)
c3_paley = count_directed_cycles(A_paley, 3)
c5_paley = count_directed_cycles(A_paley, 5)
c7_paley = count_directed_cycles(A_paley, 7)

print(f"  H(P_7) = {h_paley}")
print(f"  c3 = {c3_paley}")
print(f"  c5 = {c5_paley}")
print(f"  c7 = {c7_paley}")

# AA^T for Paley
AAT_paley = A_paley.astype(float) @ A_paley.astype(float).T
print(f"  AA^T = ")
print(AAT_paley.astype(int))

# Is AA^T circulant? (it should be, since A is circulant)
row0 = AAT_paley[0]
is_circulant = all(np.allclose(AAT_paley[i], np.roll(row0, i)) for i in range(7))
print(f"  AA^T is circulant? {is_circulant}")

off_diag = AAT_paley[np.triu_indices(7, 1)]
print(f"  Off-diagonal values: {sorted(set(off_diag.astype(int)))}")
print(f"  Off-diagonal variance: {np.var(off_diag):.6f}")

det_paley = abs(np.linalg.det(np.eye(7) + A_paley.astype(float)))
print(f"  det(I+A) = {det_paley:.2f}")
print(f"  Expected: (7+1)^4/2^7 = {8**4/128:.2f}")

# =====================================================================
# The key question: what INVARIANT distinguishes H classes?
# =====================================================================
print("\n" + "=" * 70)
print("INVARIANT ANALYSIS: WHAT DISTINGUISHES H CLASSES?")
print("=" * 70)

# Hypothesis: the directed 5-cycle count c5 and/or 7-cycle count c7
# determines the H class within regular tournaments

# Gather more data
for h_val in sorted(regular_by_H.keys()):
    group = regular_by_H[h_val]
    c5_all = []
    c7_all = []
    det_all = []

    for bits, A in group[:30]:
        c5_all.append(count_directed_cycles(A, 5))
        c7_all.append(count_directed_cycles(A, 7))
        det_all.append(abs(np.linalg.det(np.eye(n) + A.astype(float))))

    print(f"\n  H={h_val} ({len(group)} tournaments, analyzing {len(c5_all)}):")
    print(f"    c5 range: [{min(c5_all)}, {max(c5_all)}], mean={np.mean(c5_all):.1f}")
    print(f"    c7 range: [{min(c7_all)}, {max(c7_all)}], mean={np.mean(c7_all):.1f}")
    print(f"    det(I+A): [{min(det_all):.1f}, {max(det_all):.1f}], mean={np.mean(det_all):.1f}")

    # Check: is H = f(c5, c7)?
    # If c3 is constant, then H = 1 + 2*(c3 + c5 + c7 + ...)
    # For n=7: c3 constant, so H variation comes from c5, c7 variation
    if c5_all and c7_all:
        # H = constant + 2*c5 + 2*c7 (if only 3,5,7-cycles matter)
        predicted = [1 + 2*(7*8 + c5 + c7) for c5, c7 in zip(c5_all, c7_all)]
        # Wait, need the correct formula
        # OCF: H = 1 + 2*alpha_1 where alpha_1 = sum of directed cycle counts
        # Actually H = I(Ω(T), 2) = Σ_S 2^|S| where S is independent set in cycle graph
        # The simple version: H = 1 + 2*(c3 + c5 + c7) for n=7 regular?
        # Let me just check correlation
        if len(set(c5_all)) > 1:
            corr_c5 = np.corrcoef(c5_all[:len(group[:30])],
                                  [h_val]*len(c5_all))[0,1] if len(c5_all) > 1 else 0
        if len(set(c7_all)) > 1:
            corr_c7 = np.corrcoef(c7_all[:len(group[:30])],
                                  [h_val]*len(c7_all))[0,1] if len(c7_all) > 1 else 0

print("\n\nDONE — regular_n7_classification.py complete")
