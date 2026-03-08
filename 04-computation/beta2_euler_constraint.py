#!/usr/bin/env python3
"""β_2 = 0 via Euler characteristic constraints.

For a tournament T on n vertices:
χ(T) = Σ(-1)^k dim(Ω_k) = Σ(-1)^k β_k

If β_k = 0 for k ≥ 2 (as verified at n=5), then χ = β_0 - β_1 = 1 - β_1.

Question: Is χ(T) = 1 - β_1(T) for ALL tournaments?
If so, and if we can independently compute χ, then β_2 = 0 follows from:
β_2 = 0 iff χ = 1 - β_1 + β_2 - ... = 1 - β_1 (i.e., no higher contributions)

But this is circular unless we can prove the higher β_k vanish independently.

ALTERNATIVE: Maybe we can show χ = 1 - β_1 by COMPUTING χ = Σ(-1)^k dim(Ω_k)
and showing it equals 1 - β_1. This would at least constrain the sum
β_2 - β_3 + β_4 - ... = 0.

At n=7: β_3 can be 1, so β_2 - β_3 + β_4 - ... = 0 requires
β_2 = β_3 - β_4 + ... At n=7, if β_3=1 and β_k=0 for k≥4, then β_2=1.
But we verified β_2=0! So β_2 ≠ β_3 - β_4 + ... in general.
This means χ ≠ 1 - β_1 at n=7 when β_3=1.

Let me verify: what IS χ for tournaments with β_3=1 at n=7?
"""
import numpy as np
from itertools import combinations
import sys, time, random
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix,
    path_betti_numbers
)

def all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

# ===== n=5: Verify χ = 1 - β_1 =====
print("=" * 70)
print("EULER CHARACTERISTIC vs BETTI NUMBERS")
print("=" * 70)

for n in [4, 5]:
    print(f"\n--- n={n} ---")
    mismatch = 0
    total = 0
    for A in all_tournaments_gen(n):
        total += 1
        # Compute dim(Ω_k) for all k
        dims = [n]  # Ω_0 = n
        aps = {0: [[i] for i in range(n)]}
        for d in range(1, n):
            aps[d] = enumerate_allowed_paths(A, n, d)
            if d == 1:
                dims.append(len(aps[d]))
            else:
                om = compute_omega_basis(A, n, d, aps[d], aps[d-1])
                dims.append(om.shape[1] if om.ndim == 2 else 0)

        chi = sum((-1)**k * dims[k] for k in range(n))
        betti = path_betti_numbers(A, n, max_dim=n-1)
        chi_betti = sum((-1)**k * betti[k] for k in range(len(betti)))

        if chi != chi_betti:
            mismatch += 1
            print(f"  MISMATCH: χ_Ω={chi}, χ_β={chi_betti}")

        # Check: is χ = 1 - β_1?
        if chi != 1 - betti[1]:
            if n <= 5:
                t3 = sum(1 for a, b, c in combinations(range(n), 3)
                         if (A[a][b] and A[b][c] and A[c][a]) or
                            (A[b][a] and A[a][c] and A[c][b]))
                print(f"  χ={chi} ≠ 1-β_1={1-betti[1]} (β={betti}, t3={t3})")

    print(f"  χ mismatches: {mismatch}/{total}")

# ===== n=7: Sample to check χ and β =====
print(f"\n\n{'='*70}")
print("n=7: EULER CHARACTERISTIC AND BETTI NUMBERS (sample)")
print("="*70)

n = 7
random.seed(42)
chi_data = Counter()
sample_size = 200

count = 0
for A in all_tournaments_gen(n):
    count += 1
    if random.random() > sample_size / (2**21):
        continue

    # Compute dim(Ω_k) for all k
    dims = [n]
    for d in range(1, n):
        ap = enumerate_allowed_paths(A, n, d)
        if d == 1:
            dims.append(len(ap))
        else:
            ap_dm1 = enumerate_allowed_paths(A, n, d-1)
            om = compute_omega_basis(A, n, d, ap, ap_dm1)
            dims.append(om.shape[1] if om.ndim == 2 else 0)

    chi = sum((-1)**k * dims[k] for k in range(n))
    betti = path_betti_numbers(A, n, max_dim=n-1)
    chi_betti = sum((-1)**k * betti[k] for k in range(len(betti)))

    b_tuple = tuple(betti)
    chi_data[(chi, b_tuple)] += 1

    if chi_data[(chi, b_tuple)] <= 2:
        t3 = sum(1 for a, b, c in combinations(range(n), 3)
                 if (A[a][b] and A[b][c] and A[c][a]) or
                    (A[b][a] and A[a][c] and A[c][b]))
        eq_str = "=" if chi == 1 - betti[1] else "≠"
        print(f"  χ={chi}, β={betti}, Ω={dims}, t3={t3} (χ {eq_str} 1-β_1)")

    if count % 500000 == 0:
        print(f"  ... {count}", flush=True)

print(f"\nSummary:")
print(f"  (χ, β): count")
for k in sorted(chi_data):
    chi, betti = k
    eq = "=" if chi == 1 - betti[1] else "≠"
    print(f"    χ={chi} {eq} 1-β_1, β={list(betti)}: {chi_data[k]}")

print("\nDone.")
