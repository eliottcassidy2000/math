#!/usr/bin/env python3
"""
beta3_analysis_n7.py - Analyze beta_3 != 0 tournaments at n=7

Key question: WHY is beta_2 always 0 while beta_3 can be nonzero?
What's special about dimension 2?

Strategy:
1. Find explicit n=7 tournaments with beta_3 = 1
2. Analyze their Omega chain complex
3. Compare with beta_3 = 0 tournaments
4. Look for structural differences that explain why beta_2 = 0

Author: kind-pasteur-2026-03-08-S42
"""
import sys, os, time, random
import numpy as np
from itertools import combinations
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, build_full_boundary_matrix,
    compute_omega_basis, path_betti_numbers
)
sys.stdout = _saved

random.seed(42)


def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A


def random_tournament(n):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    return A


def count_3cycles(A, n):
    c3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[j][i] and A[i][k] and A[k][j]):
                    c3 += 1
    return c3


def full_chain_dims(A, n, max_p=None):
    """Get Omega dimensions and ranks."""
    if max_p is None:
        max_p = n - 1

    allowed = {}
    for p in range(-1, max_p + 2):
        if p < 0:
            allowed[p] = []
        else:
            allowed[p] = enumerate_allowed_paths(A, n, p)

    results = {}
    for p in range(max_p + 1):
        omega_p = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
        dim_omega = omega_p.shape[1] if omega_p.ndim == 2 else 0
        dim_A = len(allowed[p])

        # rank of d_p
        if dim_omega > 0 and len(allowed[p-1]) > 0:
            D = build_full_boundary_matrix(allowed[p], allowed[p-1])
            D_omega = D @ omega_p
            S = np.linalg.svd(D_omega, compute_uv=False)
            rank = sum(s > 1e-8 for s in S)
        else:
            rank = 0

        # im(d_{p+1})
        dim_omega_p1 = 0
        if p + 1 <= max_p:
            omega_p1 = compute_omega_basis(A, n, p+1, allowed[p+1], allowed[p])
            dim_omega_p1 = omega_p1.shape[1] if omega_p1.ndim == 2 else 0
            if dim_omega_p1 > 0 and len(allowed[p]) > 0:
                D1 = build_full_boundary_matrix(allowed[p+1], allowed[p])
                D1_omega = D1 @ omega_p1
                S1 = np.linalg.svd(D1_omega, compute_uv=False)
                im = sum(s > 1e-8 for s in S1)
            else:
                im = 0
        else:
            im = 0

        ker = dim_omega - rank
        beta = ker - im

        results[p] = {
            'dim_A': dim_A, 'dim_omega': dim_omega,
            'rank_d': rank, 'ker_d': ker, 'im_next': im,
            'beta': beta
        }

    return results


# ============================================================
# PART 1: Find beta_3 != 0 examples at n=7
# ============================================================
print("=" * 70)
print("FINDING beta_3 != 0 TOURNAMENTS AT n=7")
print("=" * 70)

n = 7
beta3_examples = []
beta3_zero_examples = []
beta1_one_examples = []
t0 = time.time()

for trial in range(5000):
    A = random_tournament(n)
    betti = path_betti_numbers(A, n, max_dim=5)
    scores = tuple(sorted([sum(row) for row in A]))
    c3 = count_3cycles(A, n)

    if len(betti) > 3 and betti[3] != 0:
        beta3_examples.append((A, betti, scores, c3))
        if len(beta3_examples) <= 3:
            print(f"  beta_3!=0 example #{len(beta3_examples)}: beta={betti}, scores={scores}, c3={c3}")
    elif betti[1] == 0 and len(beta3_zero_examples) < 3:
        beta3_zero_examples.append((A, betti, scores, c3))
    elif betti[1] == 1 and len(beta1_one_examples) < 3:
        beta1_one_examples.append((A, betti, scores, c3))

elapsed = time.time() - t0
print(f"\n  Found {len(beta3_examples)} tournaments with beta_3!=0 ({elapsed:.0f}s)")


# ============================================================
# PART 2: Analyze beta_3!=0 tournaments in detail
# ============================================================
print(f"\n{'='*70}")
print("DETAILED ANALYSIS OF beta_3!=0 TOURNAMENTS")
print("=" * 70)

for idx, (A, betti, scores, c3) in enumerate(beta3_examples[:5]):
    print(f"\n--- Example #{idx+1}: beta={betti}, scores={scores}, c3={c3} ---")
    res = full_chain_dims(A, n, max_p=5)

    for p in range(6):
        if p in res:
            r = res[p]
            print(f"  p={p}: dim(A)={r['dim_A']}, dim(Omega)={r['dim_omega']}, "
                  f"rank(d)={r['rank_d']}, ker={r['ker_d']}, im(d+1)={r['im_next']}, "
                  f"beta={r['beta']}")

    # Key ratios
    print(f"  Omega dims: {[res[p]['dim_omega'] for p in range(6)]}")
    print(f"  A dims:     {[res[p]['dim_A'] for p in range(6)]}")


# Same for beta_3=0 comparison
print(f"\n--- COMPARISON: beta_3=0 tournament ---")
for idx, (A, betti, scores, c3) in enumerate(beta3_zero_examples[:2]):
    print(f"\n  beta={betti}, scores={scores}, c3={c3}")
    res = full_chain_dims(A, n, max_p=5)
    for p in range(6):
        if p in res:
            r = res[p]
            print(f"    p={p}: dim(O)={r['dim_omega']}, rank(d)={r['rank_d']}, ker={r['ker_d']}, im(d+1)={r['im_next']}, beta={r['beta']}")
    print(f"    Omega dims: {[res[p]['dim_omega'] for p in range(6)]}")


# ============================================================
# PART 3: What tournament properties correlate with beta_3!=0?
# ============================================================
print(f"\n{'='*70}")
print("TOURNAMENT PROPERTIES vs beta_3")
print("=" * 70)

score_b3 = defaultdict(lambda: [0, 0])  # [count_b3_0, count_b3_nonzero]
c3_b3 = defaultdict(lambda: [0, 0])

for (A, betti, scores, c3) in beta3_examples:
    score_b3[scores][1] += 1
    c3_b3[c3][1] += 1

# Also sample more to get b3=0 stats
for trial in range(5000):
    A = random_tournament(n)
    betti = path_betti_numbers(A, n, max_dim=3)
    scores = tuple(sorted([sum(row) for row in A]))
    c3 = count_3cycles(A, n)

    if len(betti) > 3 and betti[3] == 0:
        score_b3[scores][0] += 1
        c3_b3[c3][0] += 1
    elif len(betti) > 3 and betti[3] != 0:
        score_b3[scores][1] += 1
        c3_b3[c3][1] += 1

print(f"\nbeta_3 by c3:")
for c3 in sorted(c3_b3.keys()):
    b0, b1 = c3_b3[c3]
    rate = 100*b1/(b0+b1) if b0+b1 > 0 else 0
    print(f"  c3={c3}: beta_3=0: {b0}, beta_3!=0: {b1} ({rate:.1f}%)")

print(f"\nbeta_3 by score (top 10 most common with beta_3!=0):")
score_items = sorted(score_b3.items(), key=lambda x: -x[1][1])
for scores, (b0, b1) in score_items[:10]:
    rate = 100*b1/(b0+b1) if b0+b1 > 0 else 0
    print(f"  scores={scores}: beta_3=0: {b0}, beta_3!=0: {b1} ({rate:.1f}%)")


# ============================================================
# PART 4: Paley T_7 detailed analysis
# ============================================================
print(f"\n{'='*70}")
print("PALEY T_7 DETAILED ANALYSIS")
print("=" * 70)

# Paley T_7: QR mod 7 = {1, 2, 4}
n = 7
qr = {1, 2, 4}
A = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(n):
        if i != j and ((j - i) % n) in qr:
            A[i][j] = 1

c3 = count_3cycles(A, n)
scores = tuple(sorted([sum(row) for row in A]))
betti = path_betti_numbers(A, n, max_dim=5)
print(f"  Paley T_7: beta={betti}, scores={scores}, c3={c3}")

res = full_chain_dims(A, n, max_p=5)
for p in range(6):
    if p in res:
        r = res[p]
        print(f"  p={p}: dim(A)={r['dim_A']}, dim(Omega)={r['dim_omega']}, "
              f"rank(d)={r['rank_d']}, ker={r['ker_d']}, im(d+1)={r['im_next']}, "
              f"beta={r['beta']}")
print(f"  Omega dims: {[res[p]['dim_omega'] for p in range(6)]}")
print(f"  A dims:     {[res[p]['dim_A'] for p in range(6)]}")


# ============================================================
# PART 5: What dimension-specific property makes beta_2 = 0?
# ============================================================
print(f"\n{'='*70}")
print("DIMENSION-SPECIFIC ANALYSIS: p=2 vs p=3")
print("=" * 70)

# For all beta_3!=0 examples, check:
# 1. rank(d_2) + rank(d_3) = dim(Omega_2)? (exactness at p=2)
# 2. rank(d_3) + rank(d_4) = dim(Omega_3)? (exactness at p=3 — should FAIL)
# 3. rank(d_4) + rank(d_5) = dim(Omega_4)?

print("\nFor beta_3 != 0 examples:")
for idx, (A, betti, scores, c3) in enumerate(beta3_examples[:5]):
    res = full_chain_dims(A, n, max_p=5)
    print(f"\n  Example #{idx+1}: beta={betti}")
    for p in range(2, 5):
        rank_dp = res[p]['rank_d']
        rank_dp1 = res[p]['im_next']  # = rank(d_{p+1} on Omega_{p+1})
        dim_Op = res[p]['dim_omega']
        exact = (rank_dp + rank_dp1 == dim_Op)
        print(f"    p={p}: rank(d_{p})={rank_dp} + rank(d_{p+1})={rank_dp1} = {rank_dp+rank_dp1}, dim(Omega_{p})={dim_Op}, exact={exact}")

print("\nFor beta_3 = 0 comparison:")
for idx, (A, betti, scores, c3) in enumerate(beta3_zero_examples[:2]):
    res = full_chain_dims(A, n, max_p=5)
    print(f"\n  beta={betti}")
    for p in range(2, 5):
        rank_dp = res[p]['rank_d']
        rank_dp1 = res[p]['im_next']
        dim_Op = res[p]['dim_omega']
        exact = (rank_dp + rank_dp1 == dim_Op)
        print(f"    p={p}: rank(d_{p})={rank_dp} + rank(d_{p+1})={rank_dp1} = {rank_dp+rank_dp1}, dim(Omega_{p})={dim_Op}, exact={exact}")


# ============================================================
# PART 6: Check if n=6 exhaustive has ANY beta_3 != 0
# ============================================================
print(f"\n{'='*70}")
print("n=6: EXHAUSTIVE beta_3 CHECK")
print("=" * 70)

n = 6
total = 1 << (n*(n-1)//2)
t0 = time.time()
beta3_nonzero = 0
beta2_nonzero = 0

for bits in range(total):
    A = build_adj(n, bits)
    betti = path_betti_numbers(A, n, max_dim=4)
    if len(betti) > 3 and betti[3] != 0:
        beta3_nonzero += 1
        if beta3_nonzero <= 3:
            scores = tuple(sorted([sum(row) for row in A]))
            print(f"  beta_3!=0: T#{bits} scores={scores} beta={betti}")
    if len(betti) > 2 and betti[2] != 0:
        beta2_nonzero += 1

    if (bits+1) % 5000 == 0:
        elapsed = time.time() - t0
        print(f"  {bits+1}/{total} ({elapsed:.0f}s) b2_nz={beta2_nonzero} b3_nz={beta3_nonzero}")

elapsed = time.time() - t0
print(f"\nn=6: ({elapsed:.0f}s)")
print(f"  beta_2 != 0: {beta2_nonzero}/{total}")
print(f"  beta_3 != 0: {beta3_nonzero}/{total}")


print("\n\nDone.")
