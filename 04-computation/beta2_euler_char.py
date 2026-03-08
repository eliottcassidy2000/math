#!/usr/bin/env python3
"""
beta2_euler_char.py — Relative Euler characteristic and Omega_2 = TT span

Key questions:
1. Is dim(Omega_2) = #transitive_triples? (Yes empirically)
2. Does Omega_2 = span(TT paths)? If so, NT paths add nothing.
3. What is Sigma_v dim(Omega_p^rel) combinatorially?

Author: opus-2026-03-08-S49
"""
import sys, time
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = __import__('os').fdopen(__import__('os').open(__import__('os').devnull, __import__('os').O_WRONLY), 'w')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def dim_om(om):
    return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx): A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

def compute_omega_dims(A, n, max_p=4):
    dims = {}
    for p in range(max_p + 1):
        ap = enumerate_allowed_paths(A, n, p)
        if p == 0:
            dims[p] = n
        elif not ap:
            dims[p] = 0
        else:
            apm1 = enumerate_allowed_paths(A, n, p-1)
            om = compute_omega_basis(A, n, p, ap, apm1)
            dims[p] = dim_om(om)
    return dims

print("=" * 70)
print("OMEGA_2 STRUCTURE AND RELATIVE EULER CHAR")
print("=" * 70)

n = 5
m = n*(n-1)//2
total = 1 << m

# Part 1: dim(Omega_2) = #TT?
print("\n--- Part 1: dim(Omega_2) vs #TT ---")
match = 0
for bits in range(total):
    A = build_adj(n, bits)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap1 = enumerate_allowed_paths(A, n, 1)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    d2 = dim_om(om2)
    tt_count = sum(1 for p in ap2 if A[p[0]][p[2]])
    if d2 == tt_count:
        match += 1
print(f"  dim(Omega_2) = #TT: {match}/{total} ({'ALL' if match==total else 'NOT ALL'})")

# Part 2: Sigma_v dim(Omega_p^rel) for p=0,1,2,3
print(f"\n--- Part 2: Sigma_v dim(Omega_p^rel) ---")
sum_dims = defaultdict(Counter)
t0 = time.time()
for bits in range(total):
    A = build_adj(n, bits)
    dims_T = compute_omega_dims(A, n, max_p=4)
    sums = {p: 0 for p in range(5)}
    for v in range(n):
        others = [i for i in range(n) if i != v]
        n1 = n-1
        A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
        dims_sub = compute_omega_dims(A_sub, n1, max_p=4)
        for p in range(5):
            sums[p] += dims_T[p] - dims_sub[p]
    for p in range(5):
        sum_dims[p][sums[p]] += 1

print(f"  n={n} ({time.time()-t0:.0f}s):")
for p in range(5):
    print(f"    p={p}: Sigma_v dim(Omega_p^rel) = {dict(sorted(sum_dims[p].items()))}")

# Key: p=0 should always be n. p=1 should always be 2*C(n,2) = n(n-1).
# p=2 should always be 3*(C(n,3)-t3) (each TT counted 3 times).
# Let me verify p=2 formula.
print(f"\n--- Part 3: Verify Sigma_v dim(Omega_2^rel) = 3*(C(n,3)-t3) ---")
mismatch = 0
for bits in range(total):
    A = build_adj(n, bits)
    t3 = 0
    for i in range(n):
        for j in range(i+1,n):
            for k in range(j+1,n):
                if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
                    t3 += 1
    expected = 3*(10 - t3)  # C(5,3)=10
    dims_T = compute_omega_dims(A, n, max_p=2)
    actual = 0
    for v in range(n):
        others = [i for i in range(n) if i != v]
        n1 = n-1
        A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
        dims_sub = compute_omega_dims(A_sub, n1, max_p=2)
        actual += dims_T[2] - dims_sub[2]
    if actual != expected:
        mismatch += 1
        if mismatch <= 3:
            print(f"  bits={bits}: t3={t3}, expected={expected}, actual={actual}")

if mismatch == 0:
    print(f"  CONFIRMED: Sigma_v dim(Omega_2^rel) = 3*(C(n,3)-t3) for ALL {total}")
else:
    print(f"  {mismatch} mismatches")

# Part 4: Check Sigma_v dim(Omega_3^rel) formula
# Omega_3 = ? Let's see dim(Omega_3) distribution
print(f"\n--- Part 4: Sigma_v dim(Omega_3^rel) ---")
# Omega_3 elements are "doubly-transitive" 4-paths
# A 4-path (a,b,c,d): needs ∂(abcd) in Omega_2 = TT span
# For a single path: ∂_3(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)
# Each of these 3-element sub-paths must be in Omega_2 (or cancel).

sum3_dist = Counter()
for bits in range(total):
    A = build_adj(n, bits)
    dims_T = compute_omega_dims(A, n, max_p=3)
    actual = 0
    for v in range(n):
        others = [i for i in range(n) if i != v]
        n1 = n-1
        A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
        dims_sub = compute_omega_dims(A_sub, n1, max_p=3)
        actual += dims_T[3] - dims_sub[3]
    sum3_dist[actual] += 1

print(f"  Sigma_v dim(Omega_3^rel): {dict(sorted(sum3_dist.items()))}")

# Part 5: Relative Euler char = sum_p (-1)^p * Sigma_v dim(Omega_p^rel)
print(f"\n--- Part 5: Sigma_v chi^rel ---")
chi_dist = Counter()
for bits in range(total):
    A = build_adj(n, bits)
    dims_T = compute_omega_dims(A, n, max_p=4)
    chi = 0
    for v in range(n):
        others = [i for i in range(n) if i != v]
        n1 = n-1
        A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
        dims_sub = compute_omega_dims(A_sub, n1, max_p=4)
        for p in range(5):
            chi += (-1)**p * (dims_T[p] - dims_sub[p])
    chi_dist[chi] += 1

print(f"  Sigma_v chi^rel: {dict(sorted(chi_dist.items()))}")

# Part 6: For each (bits,v) pair, compute h2_rel and the breakdown:
# h2_rel = z2_rel - b2_rel
# z2_rel = dim(Omega_2^rel) - rk(d2^rel)
# b2_rel = rk(d3^rel)
# So h2_rel relates to dim(Omega_2^rel), dim(Omega_1^rel), dim(Omega_3^rel)

print(f"\n--- Part 6: h2_rel breakdown (z2_rel, b2_rel) ---")
breakdown = Counter()
for bits in range(total):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]
    dims_T = compute_omega_dims(A, n, max_p=3)
    
    for v in range(n):
        if scores[v] == 0 or scores[v] == n-1:
            continue
        others = [i for i in range(n) if i != v]
        n1 = n-1
        A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
        dims_sub = compute_omega_dims(A_sub, n1, max_p=3)
        d_rel = {p: dims_T[p] - dims_sub[p] for p in range(4)}
        breakdown[(d_rel[0], d_rel[1], d_rel[2], d_rel[3])] += 1

print(f"  (d0_rel, d1_rel, d2_rel, d3_rel) → count:")
for key in sorted(breakdown, key=lambda k: -breakdown[k]):
    print(f"    {key} → {breakdown[key]}")

print("\nDone.")
