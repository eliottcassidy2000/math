#!/usr/bin/env python3
"""
Is H an EXACT function of det(I - A)?

At n=5: det(I - A) = prod(1 - lambda_i).
The eigenvalues have Re(lambda) = -1/2 for non-Perron, so the product
is dominated by the Perron eigenvalue.

But det is NOT a complete invariant — different tournament isomorphism classes
can have the same det but different H.

Let's check: for fixed H, how many distinct values of det(I-uA) exist?
"""
import numpy as np
from itertools import permutations, combinations
from collections import defaultdict

def ham_path_count_dp(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or (mask, v) not in dp: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

# EXHAUSTIVE at n=5
n = 5
edges = [(i,j) for i in range(n) for j in range(i+1,n)]
m = len(edges)

print(f"=== EXHAUSTIVE n={n}: H vs spectral invariants ===")

# Collect all tournaments
data = defaultdict(list)
for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(edges):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1
    H = ham_path_count_dp(A, n)
    A_np = np.array(A, dtype=float)
    eigs = sorted(np.linalg.eigvals(A_np), key=lambda x: -x.real)
    char_poly = np.round(np.real(np.poly(eigs)), 6)  # characteristic polynomial coefficients
    det_val = np.linalg.det(np.eye(n) - A_np)  # det(I - A) = char_poly(1)
    det_half = np.linalg.det(np.eye(n) - 0.5 * A_np)

    scores = tuple(sorted(sum(A[i]) for i in range(n)))

    data[H].append({
        'eigs': eigs,
        'char_poly': tuple(char_poly),
        'det1': round(det_val, 6),
        'det_half': round(det_half, 6),
        'scores': scores
    })

print(f"  {'H':>3} | {'count':>5} | {'#det(I-A)':>9} | {'#det(I-A/2)':>11} | {'#char_polys':>11} | {'#score_seqs':>11}")
print("-" * 75)
for H in sorted(data.keys()):
    items = data[H]
    det1_vals = set(d['det1'] for d in items)
    deth_vals = set(d['det_half'] for d in items)
    cp_vals = set(d['char_poly'] for d in items)
    score_vals = set(d['scores'] for d in items)
    print(f"  {H:>3} | {len(items):>5} | {len(det1_vals):>9} | {len(deth_vals):>11} | {len(cp_vals):>11} | {len(score_vals):>11}")

# Show the actual values
print(f"\n  Detail:")
for H in sorted(data.keys()):
    items = data[H]
    det1_vals = sorted(set(d['det1'] for d in items))
    deth_vals = sorted(set(d['det_half'] for d in items))
    scores = sorted(set(d['scores'] for d in items))
    print(f"  H={H}: det(I-A)={det1_vals}, det(I-A/2)={deth_vals}, scores={scores}")

# ========== Now check: does the characteristic polynomial determine H? ==========
print(f"\n\n=== Does char poly determine H? ===")
cp_to_H = defaultdict(set)
for H in data:
    for d in data[H]:
        cp_to_H[d['char_poly']].add(H)

ambiguous = {cp: Hs for cp, Hs in cp_to_H.items() if len(Hs) > 1}
print(f"  {len(cp_to_H)} distinct char polys, {len(ambiguous)} ambiguous")
for cp, Hs in ambiguous.items():
    print(f"    char_poly={cp} -> H values: {sorted(Hs)}")

# ========== Score sequence vs H ==========
print(f"\n\n=== Does score sequence determine H? ===")
score_to_H = defaultdict(set)
for H in data:
    for d in data[H]:
        score_to_H[d['scores']].add(H)

ambiguous_scores = {s: Hs for s, Hs in score_to_H.items() if len(Hs) > 1}
print(f"  {len(score_to_H)} distinct score seqs, {len(ambiguous_scores)} ambiguous")
for s, Hs in sorted(ambiguous_scores.items()):
    print(f"    scores={s} -> H values: {sorted(Hs)}")

# ========== BEAUTIFUL: At n=5, score sequence ALMOST determines H ==========
# Let's check n=7 (sampled)
print(f"\n\n=== n=7: Score sequence vs H (sampled) ===")
import random
random.seed(42)
n = 7
score_to_H7 = defaultdict(set)
for _ in range(1000):
    A = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5: A[i][j] = 1
            else: A[j][i] = 1
    H = ham_path_count_dp(A, n)
    scores = tuple(sorted(sum(A[i]) for i in range(n)))
    score_to_H7[scores].add(H)

ambiguous7 = {s: Hs for s, Hs in score_to_H7.items() if len(Hs) > 1}
print(f"  {len(score_to_H7)} distinct score seqs observed")
print(f"  {len(ambiguous7)} with multiple H values")
for s, Hs in sorted(ambiguous7.items()):
    print(f"    scores={s} -> H values: {sorted(Hs)}")
