#!/usr/bin/env python3
"""
beta2_rel_betti_sum.py — Sum of relative Betti numbers

From the exact sequence of the pair (T, T\\v):
β₀^rel(T,T\\v) = 1 always (source vertex contributes 1 to H₀^rel)

Wait, β₀^rel = dim H₀(T,T\\v). From LES:
H₀(T\\v) → H₀(T) → H₀(T,T\\v) → 0
H₀(T) = k (1-dim) since tournament is connected.
H₀(T\\v) = k (1-dim) since T\\v is connected for interior v.
Map i*: H₀(T\\v) → H₀(T) is surjective (induced by inclusion).
So H₀(T,T\\v) = 0 for interior v!

Actually wait: for NON-interior v (source/sink), T\\v might be disconnected.
A sink has d⁻ = n-1 (everyone beats it). T\\sink is still a tournament on n-1.
A source has d⁺ = n-1. T\\source is still a tournament on n-1.
Both are connected (every tournament is strongly connected or has a Hamiltonian path).

Hmm, actually not every tournament is strongly connected. But every tournament is
connected in the underlying undirected sense (as undirected graph). And H₀ for path
homology just counts connected components. Since T\\v is still a complete digraph
minus one vertex, it has a Hamiltonian path and is connected.

So β₀^rel = 0 for all v? Let me verify.

Also: the relative Euler char Σ_v χ^rel relates to the alternating sum of
Σ_v dim(Ω_p^rel). I computed this earlier.

If β₀^rel = 0 for all v, then:
Σ_v h₂_rel = Σ_v χ^rel + Σ_v β₁^rel - Σ_v β₃^rel + ...

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


def compute_all_betti(A, n, max_p=4):
    """Compute β_p(T) for p = 0,...,max_p."""
    betti = {}
    oms = {}
    aps = {}
    for p in range(max_p + 1):
        aps[p] = enumerate_allowed_paths(A, n, p)
        if p == 0:
            oms[p] = np.eye(n)
            betti[0] = 1  # Tournament always connected
        elif not aps[p]:
            oms[p] = np.zeros((0, 0))
            betti[p] = 0
        else:
            oms[p] = compute_omega_basis(A, n, p, aps[p], aps[p-1])
            dp = dim_om(oms[p])
            if dp == 0:
                betti[p] = 0
                continue

            # Z_p = ker(∂_p)
            bd = build_full_boundary_matrix(aps[p], aps[p-1])
            if p == 1:
                # ∂₁: Ω₁ → A₀
                d_mat = bd @ oms[1]
            else:
                d_mat = np.linalg.lstsq(oms[p-1], bd @ oms[p], rcond=None)[0]
            rk_d = np.linalg.matrix_rank(d_mat, tol=1e-8)
            zp = dp - rk_d

            # B_p = im(∂_{p+1})
            if p < max_p and aps.get(p+1) and dim_om(oms.get(p+1, np.zeros((0,0)))) > 0:
                # Need oms[p+1] which might not be computed yet
                if p+1 not in oms:
                    aps[p+1] = enumerate_allowed_paths(A, n, p+1)
                    oms[p+1] = compute_omega_basis(A, n, p+1, aps[p+1], aps[p]) if aps[p+1] else np.zeros((0,0))
                dp1 = dim_om(oms[p+1])
                if dp1 > 0:
                    bd1 = build_full_boundary_matrix(aps[p+1], aps[p])
                    d1_mat = np.linalg.lstsq(oms[p], bd1 @ oms[p+1], rcond=None)[0]
                    bp = np.linalg.matrix_rank(d1_mat, tol=1e-8)
                else:
                    bp = 0
            else:
                bp = 0

            betti[p] = zp - bp

    return betti


def compute_rel_betti(A, n, v, max_p=3):
    """Compute relative Betti numbers β_p^rel(T, T\\v) using LES.

    From LES: ... → H_p(T\\v) →^{i*} H_p(T) →^{j*} H_p(T,T\\v) →^δ H_{p-1}(T\\v) → ...

    β_p^rel = dim H_p(T,T\\v) = dim coker(i*_p) + dim ker(i*_{p-1})
    wait, this isn't right. Let me use the exact formula.

    Actually: β_p^rel can be computed directly from the relative chain complex.
    """
    others = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]

    betti_T = compute_all_betti(A, n, max_p)
    betti_sub = compute_all_betti(A_sub, n1, max_p)

    # Euler characteristic of relative complex
    # χ^rel = Σ (-1)^p dim(Ω_p^rel)
    # and χ^rel = Σ (-1)^p β_p^rel

    return betti_T, betti_sub


print("=" * 70)
print("RELATIVE BETTI NUMBER SUMS")
print("=" * 70)

n = 5
m = n*(n-1)//2
total = 1 << m

# For each tournament: compute β_p(T) and β_p(T\\v) for all v
# Then compute β_p^rel from LES

# Actually, let me just compute the key quantities directly:
# From LES at each v:
# ... → H₂(T\\v) → H₂(T) → H₂(T,T\\v) → H₁(T\\v) → H₁(T) → H₁(T,T\\v) → H₀(T\\v) → H₀(T) → H₀(T,T\\v) → 0

# H₀(T) = k, H₀(T\\v) = k (both connected for interior v)
# So: H₀(T,T\\v) has dim = dim coker(H₀(T\\v) → H₀(T)) = 0 (surjective)

# At H₁ level:
# H₁(T\\v) →^{i*} H₁(T) → H₁(T,T\\v) → H₀(T\\v) → H₀(T)
# H₀(T\\v) → H₀(T) is an isomorphism (k → k).
# So ker of this map = 0, meaning im(H₁(T,T\\v) → H₀(T\\v)) = 0.
# So H₁(T,T\\v) → H₀(T\\v) is zero.
# Hence: H₁(T) → H₁(T,T\\v) is surjective.
# β₁^rel = dim coker(i*_{H₁}) + 0 = β₁(T) - rk(i*_{H₁})... no wait.

# From exactness: H₁(T) →^{j*} H₁(T,T\\v) →^{δ} H₀(T\\v) →^{i*} H₀(T)
# ker(i*) at H₀ level = im(δ at H₁ level).
# i* is iso → ker(i*) = 0 → im(δ) = 0 → δ = 0.
# → j* at H₁ is surjective.
# So β₁^rel = dim H₁(T,T\\v) = dim coker(H₁(T\\v) → H₁(T)) = β₁(T) - rk(i*_{H₁}).

# Hmm, from exactness at H₁(T):
# H₁(T\\v) →^{i*} H₁(T) →^{j*} H₁(T,T\\v)
# im(i*) = ker(j*). So rk(j*) = β₁(T) - rk(i*).
# And j* is surjective (shown above). So β₁^rel = β₁(T) - rk(i*).
# But rk(i*) ≤ min(β₁(T\\v), β₁(T)).

# So: β₁^rel = β₁(T) - rk(i*_{H₁}) where i* is inclusion on H₁.

# Similarly at H₂ level:
# ... → H₂(T\\v) →^{i*} H₂(T) →^{j*} H₂(T,T\\v) →^{δ} H₁(T\\v) →^{i*} H₁(T) → ...
# If β₂(T\\v) = 0 (induction): i* at H₂ is zero map (from 0).
# So ker(j*) = im(i*) = 0 → j* injective.
# j* maps into H₂(T,T\\v).
# From exactness: im(j*) = ker(δ). So β₂(T) = dim ker(δ).
# δ: H₂(T,T\\v) → H₁(T\\v). rk(δ) = h₂_rel - dim ker(δ) = h₂_rel - β₂(T).
# Also: im(δ) = ker(i*_{H₁}). So rk(δ) = β₁(T\\v) - rk(i*_{H₁}).
# Combining: h₂_rel - β₂(T) = β₁(T\\v) - rk(i*_{H₁}).
# So: β₂(T) = h₂_rel - β₁(T\\v) + rk(i*_{H₁}).

# If this holds for ALL v: Σ_v β₂(T) = n·β₂(T) = Σ_v h₂_rel - Σ_v β₁(T\\v) + Σ_v rk(i*_{H₁,v})

# This is a KEY IDENTITY! Let me verify it.

print("\n--- Key identity: n·β₂ = Σ h₂_rel - Σ β₁(T\\v) + Σ rk(i*_{H₁,v}) ---")

from beta2_h2rel_sum import compute_h2_rel

count = 0
for bits in range(total):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    betti_T = compute_all_betti(A, n, 3)
    beta2_T = betti_T.get(2, 0)
    beta1_T = betti_T.get(1, 0)

    h2_sum = 0
    beta1_sub_sum = 0
    rk_i_sum = 0

    for v in range(n):
        h2r = compute_h2_rel(A, n, v)
        h2_sum += h2r

        others = [i for i in range(n) if i != v]
        n1 = n-1
        A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
        betti_sub = compute_all_betti(A_sub, n1, 2)
        beta1_sub = betti_sub.get(1, 0)
        beta1_sub_sum += beta1_sub

        # Compute rk(i*_{H₁}): rank of inclusion map H₁(T\\v) → H₁(T)
        # This requires computing the actual map on homology
        # For now, use the identity: rk(i*) = β₁(T\\v) - rk(δ)
        # And rk(δ) = h₂_rel - β₂(T) (assuming β₂(T\\v) = 0)
        rk_delta = h2r - beta2_T if beta2_T <= h2r else 0
        rk_i = beta1_sub - rk_delta

        rk_i_sum += rk_i

    # Check identity: n·β₂ = h2_sum - beta1_sub_sum + rk_i_sum
    lhs = n * beta2_T
    rhs = h2_sum - beta1_sub_sum + rk_i_sum
    if abs(lhs - rhs) > 0.1:
        count += 1
        if count <= 5:
            print(f"  bits={bits}: n·β₂={lhs}, Σh₂-Σβ₁+Σrk={rhs}")

if count == 0:
    print(f"  ✓ Identity holds for all {total} tournaments")
else:
    print(f"  ✗ {count} failures")

# Also compute: Σ_v β₁(T\\v) for each tournament
print(f"\n--- Σ_v β₁(T\\v) distribution ---")
beta1_sub_sums = Counter()
for bits in range(total):
    A = build_adj(n, bits)
    s = 0
    for v in range(n):
        others = [i for i in range(n) if i != v]
        n1 = n-1
        A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
        betti_sub = compute_all_betti(A_sub, n1, 2)
        s += betti_sub.get(1, 0)
    beta1_sub_sums[s] += 1

print(f"  Σ_v β₁(T\\v): {dict(sorted(beta1_sub_sums.items()))}")

print("\nDone.")
