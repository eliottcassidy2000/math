#!/usr/bin/env python3
"""
Ω_2 CANCELLATION CHAINS IN TOURNAMENTS

We discovered: Ω_2 = span(transitive triples) ⊕ span(cancellation chains).

A cancellation chain is a linear combination of non-transitive 2-paths
(a,b,c) with c→a, where the non-allowed face terms (a,c) cancel.

STRUCTURE: For a fixed non-edge pair (a,c) with c→a:
  All 2-paths (a,*,c) that contribute -(a,c):
    (a,b,c) for each b with a→b and b→c and b∉{a,c}

  The non-allowed face (a,c) gets coefficient -1 from each such path.
  So if there are k such paths, the cancellation space for (a,c) is
  the (k-1)-dimensional space of chains summing to zero coefficient.

  But these chains also have ALLOWED face terms:
  ∂(a,b,c) = (b,c) - (a,c) + (a,b)
  The allowed terms (b,c) and (a,b) don't need to cancel.

  So a cancellation chain looks like:
  α₁(a,b₁,c) + α₂(a,b₂,c) + ... with α₁+α₂+... = 0
  and its boundary is:
  α₁(b₁,c) - α₁(a,c) + α₁(a,b₁) + α₂(b₂,c) - α₂(a,c) + ...
  = Σ αᵢ(bᵢ,c) + Σ αᵢ(a,bᵢ) - (Σ αᵢ)(a,c)
  = Σ αᵢ(bᵢ,c) + Σ αᵢ(a,bᵢ)  [since Σ αᵢ = 0]

This is a chain in A_1, so it's in Ω_2. ✓

KEY QUESTION: What do these elements contribute to ker(∂_2|Ω_2)?
If a cancellation chain z = Σ αᵢ(a,bᵢ,c) with Σαᵢ=0 satisfies ∂z=0:
  Σ αᵢ(bᵢ,c) + Σ αᵢ(a,bᵢ) = 0
  This requires each edge coefficient to be zero:
  For edge (bᵢ,c): coefficient αᵢ = 0 (from (bᵢ,c) appearing in ∂(a,bᵢ,c))
  Wait — could (bᵢ,c) also appear from another triple?

  Edge (bᵢ,c) appears in ∂(d,bᵢ,c) for any d with d→bᵢ and d→c.
  If d=a (which has a→bᵢ but c→a, so NOT a→c) — (a,bᵢ,c) is our path.
  If d≠a: (d,bᵢ,c) is a different path, possibly transitive.
"""
import numpy as np
from collections import Counter, defaultdict
import sys
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import enumerate_allowed_paths, compute_omega_basis

def all_tournaments_gen(n):
    edges = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(edges):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

# ===== Classify cancellation chains =====
print("=" * 70)
print("CANCELLATION CHAIN CLASSIFICATION")
print("=" * 70)

n = 5
chain_stats = Counter()

for t_idx, A in enumerate(all_tournaments_gen(n)):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)

    # For each non-edge (a,c) with c→a, count 2-paths (a,*,c)
    for a in range(n):
        for c in range(n):
            if a == c or A[a][c] == 1:
                continue
            # c→a, so (a,c) is not an edge
            # Count paths (a,b,c) with a→b, b→c
            intermediates = [b for b in range(n) if b != a and b != c
                           and A[a][b] == 1 and A[b][c] == 1]
            if len(intermediates) > 1:
                chain_stats[len(intermediates)] += 1

print("Non-edge pairs (c→a) with k intermediates (k>1):")
for k in sorted(chain_stats.keys()):
    print(f"  k={k}: {chain_stats[k]}")
    print(f"    Each gives {k-1} independent cancellation chains in Ω_2")

# ===== Do cancellation chains contribute to ker(∂_2)? =====
print(f"\n\n{'='*70}")
print("DO CANCELLATION CHAINS CONTRIBUTE TO ker(∂_2)?")
print("="*70)

# A cancellation chain z = (a,b1,c) - (a,b2,c) (simplest case, 2 paths)
# ∂z = (b1,c) - (a,c) + (a,b1) - (b2,c) + (a,c) - (a,b2)
#     = (b1,c) - (b2,c) + (a,b1) - (a,b2)  [(a,c) terms cancel]
#
# For ∂z = 0:
#   coeff of (b1,c) = 1, coeff of (b2,c) = -1, etc.
#   These are all different edges (b1≠b2), so ∂z ≠ 0 unless
#   (b1,c) and (b2,c) are the SAME edge... but b1≠b2 so they can't be.
#
# So the simplest cancellation chain NEVER has zero boundary!
# It always contributes to rank(∂_2|Ω_2) but not to ker(∂_2|Ω_2).

# Verify this:
n = 5
cancel_in_ker = 0
cancel_not_in_ker = 0

for t_idx, A in enumerate(all_tournaments_gen(n)):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a2_tuples = [tuple(p) for p in a2]
    a1_tuples = [tuple(p) for p in a1]
    a1_idx = {t: i for i, t in enumerate(a1_tuples)}

    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0
    tt_count = sum(1 for p in a2 if A[p[0]][p[2]] == 1)

    if dim_om2 <= tt_count:
        continue  # no cancellation chains

    # Build ∂_2: full A_2 → A_1
    bd2_full = np.zeros((len(a1_tuples), len(a2_tuples)))
    for j, (a, b, c) in enumerate(a2_tuples):
        bd2_full[a1_idx[(b,c)], j] += 1
        if (a,c) in a1_idx:
            bd2_full[a1_idx[(a,c)], j] -= 1
        bd2_full[a1_idx[(a,b)], j] += 1

    # ∂_2 on Ω_2
    bd2_om = bd2_full @ om2

    # Find the cancellation-chain basis vectors
    tt_set = set(tuple(p) for p in a2 if A[p[0]][p[2]] == 1)
    for col in range(dim_om2):
        vec = om2[:, col]
        terms = [(a2_tuples[j], vec[j]) for j in range(len(a2_tuples)) if abs(vec[j]) > 1e-8]
        has_ntt = any(t not in tt_set for t, _ in terms)
        if has_ntt:
            # This is (partly) a cancellation chain
            bd = bd2_om[:, col]
            bd_norm = np.linalg.norm(bd)
            if bd_norm < 1e-8:
                cancel_in_ker += 1
            else:
                cancel_not_in_ker += 1

print(f"Cancellation chain Ω_2 elements in ker(∂_2): {cancel_in_ker}")
print(f"Cancellation chain Ω_2 elements NOT in ker(∂_2): {cancel_not_in_ker}")

# ===== Does ker(∂_2|Ω_2) = ker(∂_2|TT)? =====
print(f"\n\n{'='*70}")
print("ker(∂_2|Ω_2) vs ker(∂_2|TT)")
print("="*70)

n = 5
always_match = True
for A in all_tournaments_gen(n):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a1_tuples = [tuple(p) for p in a1]
    a2_tuples = [tuple(p) for p in a2]
    a1_idx = {t: i for i, t in enumerate(a1_tuples)}

    tt = [tuple(p) for p in a2 if A[p[0]][p[2]] == 1]

    # ker(∂_2|TT): TT → A_1
    if len(tt) == 0:
        continue
    bd2_tt = np.zeros((len(a1_tuples), len(tt)))
    for j, (a, b, c) in enumerate(tt):
        bd2_tt[a1_idx[(b,c)], j] += 1
        bd2_tt[a1_idx[(a,c)], j] -= 1
        bd2_tt[a1_idx[(a,b)], j] += 1
    rank_tt = np.linalg.matrix_rank(bd2_tt, tol=1e-8)
    ker_tt = len(tt) - rank_tt

    # ker(∂_2|Ω_2): Ω_2 → A_1
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_om2 = om2.shape[1] if om2.ndim == 2 else 0

    # Build ∂_2 on A_2
    bd2_full = np.zeros((len(a1_tuples), len(a2_tuples)))
    for j, path in enumerate(a2_tuples):
        a, b, c = path
        for sign, face in [(1, (b,c)), (-1, (a,c)), (1, (a,b))]:
            if face in a1_idx:
                bd2_full[a1_idx[face], j] += sign

    bd2_om = bd2_full @ om2
    rank_om = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    ker_om = dim_om2 - rank_om

    if ker_tt != ker_om:
        always_match = False
        print(f"  MISMATCH: ker_TT={ker_tt}, ker_Ω_2={ker_om}, dim_Ω_2={dim_om2}, |TT|={len(tt)}")

if always_match:
    print(f"n={n}: ker(∂_2|TT) = ker(∂_2|Ω_2) for ALL tournaments")
    print("\nThis means: cancellation chains never contribute to ker(∂_2)!")
    print("So the β_2 = 0 question only involves transitive triples after all.")

print("\nDone.")
