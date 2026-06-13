#!/usr/bin/env python3
"""
beta2_simplex_deformation.py - The simplex deformation perspective

The transitive tournament T_0 has chain complex = standard simplex complex:
  Omega_p(T_0) = C(n, p+1) = dimension of p-chains on the (n-1)-simplex
  H_p(T_0) = 0 for p >= 1 (simplex is contractible)

Any tournament T is obtained from T_0 by flipping arcs. The key question:
does the "deformation" from simplex to tournament preserve H_2 = 0?

APPROACH: Track how the chain complex changes as we flip each arc.
Specifically, for a sequence of flips T_0 -> T_1 -> ... -> T,
track dim(Omega_p), rk(d_p), dim(Z_p), and beta_p at each step.

The Euler characteristic chi = sum (-1)^p dim(Omega_p) should satisfy
chi = sum (-1)^p beta_p. If beta_2 = 0, chi constrains other betas.

KEY QUESTION: Is chi constant for all tournaments? Or does it vary?

Author: kind-pasteur-2026-03-08-S41
"""
import sys, time, os
import numpy as np
from collections import Counter, defaultdict
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)

_saved = sys.stdout
sys.stdout = open(os.devnull, 'w', encoding='utf-8')
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis,
    build_full_boundary_matrix
)
sys.stdout = _saved

def build_adj(n, bits):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def compute_all_data(A, n, max_p=None):
    """Compute all chain complex data up to max_p."""
    if max_p is None:
        max_p = n - 1
    paths = {}
    omega = {}
    for p in range(max_p + 1):
        paths[p] = enumerate_allowed_paths(A, n, p)
        if p == 0:
            omega[p] = np.eye(n)
        elif len(paths[p]) > 0 and len(paths[p-1]) > 0:
            omega[p] = compute_omega_basis(A, n, p, paths[p], paths[p-1])
        else:
            omega[p] = np.zeros((max(1, len(paths[p])), 0))

    dims = {p: (omega[p].shape[1] if omega[p].ndim == 2 else 0) for p in range(max_p + 1)}
    rks = {}

    for p in range(1, max_p + 1):
        dp = dims[p]
        dp1 = dims[p-1]
        if dp > 0 and len(paths.get(p-1, [])) > 0:
            bd = build_full_boundary_matrix(paths[p], paths[p-1])
            bd_om = bd @ omega[p]
            if p >= 2 and dp1 > 0:
                coords, _, _, _ = np.linalg.lstsq(omega[p-1], bd_om, rcond=None)
                rks[p] = np.linalg.matrix_rank(coords, tol=1e-8)
            else:
                rks[p] = np.linalg.matrix_rank(bd_om, tol=1e-8)
        else:
            rks[p] = 0

    betas = {}
    for p in range(max_p):
        z_p = dims[p] - rks.get(p, 0)
        b_p = z_p - rks.get(p+1, 0)
        betas[p] = b_p

    chi_chain = sum((-1)**p * dims[p] for p in range(max_p + 1))
    chi_homology = sum((-1)**p * betas.get(p, 0) for p in range(max_p))

    return {
        'dims': dims, 'rks': rks, 'betas': betas,
        'chi_chain': chi_chain, 'chi_homology': chi_homology
    }


print("=" * 70)
print("SIMPLEX DEFORMATION PERSPECTIVE")
print("=" * 70)

# ANALYSIS 1: Euler characteristic
print("\n--- ANALYSIS 1: Euler characteristic ---")

for n in [4, 5, 6]:
    n_arcs = n*(n-1)//2
    total = 1 << n_arcs

    chi_dist = Counter()
    chi_check = 0

    if n <= 5:
        rng = range(total)
    else:
        import random
        random.seed(42)
        rng = [random.randint(0, total-1) for _ in range(3000)]

    for bits in rng:
        A = build_adj(n, bits)
        data = compute_all_data(A, n)
        chi = data['chi_chain']
        chi_dist[chi] += 1
        if data['chi_chain'] != data['chi_homology']:
            chi_check += 1

    print(f"\n  n={n}:")
    print(f"    chi (chain) distribution: {dict(sorted(chi_dist.items()))}")
    print(f"    chi_chain != chi_homology: {chi_check} cases")

    # Expected chi for the simplex:
    # chi = sum (-1)^p C(n, p+1) for p=0,...,n-1
    # = sum_{k=1}^n (-1)^{k-1} C(n,k) = 1 - sum_{k=0}^n (-1)^k C(n,k) + 1
    # = 1 - (1-1)^n + 1 = 1 (for n >= 1)
    # Wait: sum_{k=0}^n (-1)^k C(n,k) = 0 for n >= 1.
    # chi = sum_{p=0}^{n-1} (-1)^p C(n,p+1) = sum_{k=1}^n (-1)^{k-1} C(n,k)
    # = -sum_{k=1}^n (-1)^k C(n,k) = -(sum_{k=0}^n (-1)^k C(n,k) - 1) = -(0 - 1) = 1.
    print(f"    Simplex chi = 1")


# ANALYSIS 2: Full Betti number table
print(f"\n{'='*70}")
print("ANALYSIS 2: Full Betti number distributions")
print("=" * 70)

n = 5
n_arcs = n*(n-1)//2
total = 1 << n_arcs

betti_table = defaultdict(Counter)
for bits in range(total):
    A = build_adj(n, bits)
    data = compute_all_data(A, n)
    for p in range(n-1):
        betti_table[p][data['betas'].get(p, 0)] += 1

print(f"\nn={n} Betti distributions:")
for p in sorted(betti_table.keys()):
    print(f"  beta_{p}: {dict(sorted(betti_table[p].items()))}")

# Check Euler-Poincare
chi_from_betti = Counter()
for bits in range(total):
    A = build_adj(n, bits)
    data = compute_all_data(A, n)
    chi_b = sum((-1)**p * data['betas'].get(p, 0) for p in range(n-1))
    chi_from_betti[chi_b] += 1

print(f"\n  chi from betas: {dict(sorted(chi_from_betti.items()))}")


# ANALYSIS 3: dim(Omega_p) as a function of tournament
print(f"\n{'='*70}")
print("ANALYSIS 3: Omega dimensions for all tournaments at n=5")
print("=" * 70)

dim_tuples = Counter()
for bits in range(total):
    A = build_adj(n, bits)
    data = compute_all_data(A, n)
    dtuple = tuple(data['dims'][p] for p in range(n))
    dim_tuples[dtuple] += 1

print(f"Distinct dim tuples: {len(dim_tuples)}")
print(f"Distribution:")
for dtuple in sorted(dim_tuples.keys()):
    print(f"  Om={dtuple}: {dim_tuples[dtuple]} tournaments")


# ANALYSIS 4: Rank tuples
print(f"\n{'='*70}")
print("ANALYSIS 4: Rank tuples for all tournaments at n=5")
print("=" * 70)

rk_tuples = Counter()
for bits in range(total):
    A = build_adj(n, bits)
    data = compute_all_data(A, n)
    rtuple = tuple(data['rks'].get(p, 0) for p in range(1, n))
    rk_tuples[rtuple] += 1

print(f"Distinct rank tuples: {len(rk_tuples)}")
print(f"Distribution:")
for rtuple in sorted(rk_tuples.keys()):
    dims_tup = None
    # Find corresponding dim tuple
    for bits in range(total):
        A = build_adj(n, bits)
        data = compute_all_data(A, n)
        if tuple(data['rks'].get(p, 0) for p in range(1, n)) == rtuple:
            dims_tup = tuple(data['dims'][p] for p in range(n))
            betas_tup = tuple(data['betas'].get(p, 0) for p in range(n-1))
            break
    print(f"  rk={rtuple}: {rk_tuples[rtuple]} tours, dims={dims_tup}, betas={betas_tup}")


# ANALYSIS 5: Track a specific deformation path from transitive to cyclic
print(f"\n{'='*70}")
print("ANALYSIS 5: Deformation path from transitive to 3-cycle")
print("=" * 70)

# Transitive tournament: 0->1->2->3->4, 0->2, 0->3, 0->4, 1->3, 1->4, 2->4
# = all arcs go from smaller to larger index
# bits = all 1s (by our encoding: bit i corresponds to pair (a,b) with a<b, bit=1 means a->b)
all_ones = (1 << n_arcs) - 1
A_trans = build_adj(n, all_ones)
print(f"Transitive: scores={[sum(row) for row in A_trans]}")
data_trans = compute_all_data(A_trans, n)
print(f"  dims: {tuple(data_trans['dims'][p] for p in range(n))}")
print(f"  rks: {tuple(data_trans['rks'].get(p,0) for p in range(1,n))}")
print(f"  betas: {tuple(data_trans['betas'].get(p,0) for p in range(n-1))}")

# Now flip arcs one by one, tracking the chain complex
# Let's flip the arc 4->3 (reverse the "last" arc)
print("\nFlipping arcs one at a time from transitive:")
current_bits = all_ones

# List of arcs to flip to create a specific tournament (e.g., the cyclic one)
# The cyclic tournament on 5 has 0->1->2->3->4->0 and 0->2, 1->3, 2->4, 3->0, 4->1
# Let me just flip some arcs and see

arcs_to_flip = []
for i in range(n):
    for j in range(i+1, n):
        # Find which arcs to flip to go from transitive to some target
        idx = 0
        for a in range(n):
            for b in range(a+1, n):
                if a == i and b == j:
                    arcs_to_flip.append((i, j, idx))
                idx += 1

# Flip a few arcs and track
for step in range(min(5, len(arcs_to_flip))):
    i, j, idx = arcs_to_flip[-(step+1)]  # flip from the "largest" arcs first
    current_bits ^= (1 << idx)
    A = build_adj(n, current_bits)
    data = compute_all_data(A, n)
    scores = [sum(row) for row in A]
    print(f"\n  Step {step+1}: flip ({i},{j}), scores={scores}")
    print(f"    dims: {tuple(data['dims'][p] for p in range(n))}")
    print(f"    rks: {tuple(data['rks'].get(p,0) for p in range(1,n))}")
    print(f"    betas: {tuple(data['betas'].get(p,0) for p in range(n-1))}")
    print(f"    chi: {data['chi_chain']}")


# ANALYSIS 6: Is dim(Omega_2) + dim(Omega_4) constant?
# Or any other symmetric combination?
print(f"\n{'='*70}")
print("ANALYSIS 6: Symmetric combinations of Omega dimensions")
print("=" * 70)

# For the simplex: dims are C(n,1), C(n,2), ..., C(n,n)
# which satisfy C(n,k) = C(n, n-k), i.e., palindromic.
# For general tournaments, is there any palindromic-like relation?

sums = Counter()
for bits in range(total):
    A = build_adj(n, bits)
    data = compute_all_data(A, n)
    # Try: dim(Om_0) + dim(Om_4) and dim(Om_1) + dim(Om_3)
    s04 = data['dims'][0] + data['dims'].get(4, 0)
    s13 = data['dims'][1] + data['dims'].get(3, 0)
    s2 = data['dims'][2]
    sums[(s04, s13, s2)] += 1

print(f"(dim_Om0+Om4, dim_Om1+Om3, dim_Om2) distribution:")
for key in sorted(sums.keys()):
    print(f"  {key}: {sums[key]} tournaments")


# ANALYSIS 7: Dold-Kan type relationship -- is there a universal formula
# for rk(d_p) in terms of dim(Omega_q)?
print(f"\n{'='*70}")
print("ANALYSIS 7: Rank formulas")
print("=" * 70)

# We know: rk(d_p) + rk(d_{p+1}) + beta_p = dim(Omega_p)
# And beta_0 = 1, beta_2 = 0 (conjectured)
# So: rk(d_1) = dim(Om_0) - 1 = n - 1 (constant)
#     rk(d_2) + rk(d_3) = dim(Om_2)  (from beta_2 = 0)
#     rk(d_2) = dim(Om_1) - dim(Z_1) = C(n,2) - rk(d_2)... no.
#     dim(Om_1) = rk(d_1) + dim(Z_1)
#     dim(Z_1) = C(n,2) - (n-1) = constant
#     rk(d_2) = dim(Z_1) - beta_1 = C(n,2) - (n-1) - beta_1

# So the chain of equalities:
# rk(d_1) = n - 1  (constant)
# rk(d_2) = C(n-1,2) - beta_1  (depends only on beta_1)
# rk(d_3) = dim(Om_2) - rk(d_2) = dim(Om_2) - C(n-1,2) + beta_1  (from beta_2=0)
# rk(d_4) = dim(Om_3) - rk(d_3) - beta_3 = dim(Om_3) - dim(Om_2) + C(n-1,2) - beta_1 - beta_3

# From chi = 1: 1 - beta_1 + 0 - beta_3 + beta_4 - ... = 1
# So: beta_1 + beta_3 - beta_4 + ... = 0

print("Rank relationships (from beta_2 = 0 and chi):")
print(f"  rk(d_1) = n - 1 = {n-1}")
print(f"  rk(d_2) = C(n-1,2) - beta_1 = {(n-1)*(n-2)//2} - beta_1")
print(f"  rk(d_3) = dim(Om_2) - C(n-1,2) + beta_1  [from beta_2=0]")
print(f"  chi = 1 always => beta_1 + beta_3 = beta_4 + beta_5 + ...")

# Verify chi = 1 for all n=5 tournaments
chi_values = Counter()
for bits in range(total):
    A = build_adj(n, bits)
    data = compute_all_data(A, n)
    chi_values[data['chi_chain']] += 1

print(f"\n  chi distribution at n={n}: {dict(sorted(chi_values.items()))}")
if len(chi_values) == 1 and 1 in chi_values:
    print("  chi = 1 for ALL tournaments! (CONFIRMED)")
else:
    print("  chi is NOT always 1")


print("\nDone.")
