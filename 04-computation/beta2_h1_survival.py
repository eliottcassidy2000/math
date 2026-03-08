#!/usr/bin/env python3
"""
beta2_h1_survival.py — H₁ survival under inclusion

From the LES: exactness at H₁(T\\v) gives ker(i*) = im(δ).
When h₂_rel=1, β₁(T\\v)=1:
  - δ injective (rk=1) ⟺ i* = 0 ⟺ H₁ class of T\\v dies in T
  - δ = 0 ⟺ i* injective ⟺ H₁ class of T\\v survives in T

So: h₂_rel contributes to β₂(T) iff H₁ survives.
And: Σ contributions = #{v : H₁(T\\v) survives under inclusion into T}.

The question is: how many vertices can have SURVIVING H₁?

Also test: does i*: H₁(T\\v) → H₁(T) being zero relate to
the 3-cycle that generates H₁(T\\v) being "filled" in T?

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


def compute_h1_data(A, n):
    """Compute β₁ and Z₁ basis, B₁ basis for tournament T."""
    ap0 = enumerate_allowed_paths(A, n, 0)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    d1 = dim_om(om1); d2 = dim_om(om2)

    # ∂₁ in Ω₁ coords: maps Ω₁ → A₀
    bd1 = build_full_boundary_matrix(ap1, ap0)
    bd1_om = bd1 @ om1  # n × d1

    # Z₁ = ker(∂₁|_{Ω₁})
    U, S, Vt = np.linalg.svd(bd1_om, full_matrices=True)
    rk = sum(s > 1e-8 for s in S)
    z1_basis = Vt[rk:, :]  # (z1_dim × d1), in Ω₁ coords

    # B₁ = im(∂₂|_{Ω₂} → Ω₁)
    if d2 > 0:
        bd2 = build_full_boundary_matrix(ap2, ap1)
        bd2_om = np.linalg.lstsq(om1, bd2 @ om2, rcond=None)[0]  # d1 × d2
        b1_rk = np.linalg.matrix_rank(bd2_om, tol=1e-8)
    else:
        bd2_om = np.zeros((d1, 0))
        b1_rk = 0

    beta1 = z1_basis.shape[0] - b1_rk

    return {
        'beta1': beta1, 'z1_basis': z1_basis, 'bd2_om': bd2_om,
        'om1': om1, 'ap1': ap1, 'd1': d1, 'b1_rk': b1_rk
    }


def check_h1_survival(A, n, v):
    """Check if the H₁ class of T\\v survives under inclusion i: T\\v → T.

    Returns: (h₂_rel, β₁_sub, survives)
    """
    others = [i for i in range(n) if i != v]
    n1 = n - 1
    A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
    remap = {i: others[i] for i in range(n1)}

    # Compute H₁ data for T and T\\v
    data_T = compute_h1_data(A, n)
    data_sub = compute_h1_data(A_sub, n1)

    beta1_T = data_T['beta1']
    beta1_sub = data_sub['beta1']

    if beta1_sub == 0:
        return (0, 0, None)  # No H₁ to survive

    # Find H₁ generator of T\\v (a Z₁ class not in B₁)
    z1_sub = data_sub['z1_basis']  # in Ω₁(T\\v) coords
    om1_sub = data_sub['om1']
    ap1_sub = data_sub['ap1']

    # Map this cycle into T's chain complex via inclusion
    # Ω₁(T\\v) → Ω₁(T) via embedding
    ap1_T = data_T['ap1']
    om1_T = data_T['om1']
    ap1_T_list = [tuple(p) for p in ap1_T]
    d1_T = data_T['d1']
    d1_sub = data_sub['d1']

    # Build inclusion matrix: Ω₁(T\\v) coords → Ω₁(T) coords
    # First: ap1_sub paths → ap1_T paths (via remap)
    embed1 = np.zeros((len(ap1_T_list), d1_sub))
    for j in range(d1_sub):
        for k, path_sub in enumerate(ap1_sub):
            path_T = tuple(remap[x] for x in path_sub)
            if path_T in ap1_T_list:
                embed1[ap1_T_list.index(path_T), j] = om1_sub[k, j]
    # Convert to Ω₁(T) coords
    incl = np.linalg.lstsq(om1_T, embed1, rcond=None)[0]  # d1_T × d1_sub

    # The H₁ generator in T\\v is a vector in z1_sub not in B₁(T\\v)
    # Map via inclusion to T: incl @ z1_sub_gen (in Ω₁(T) coords)
    # Check if it's in Z₁(T) ∩ B₁(T) (i.e., a boundary in T)
    # If it's a boundary in T: H₁ class dies
    # If it's NOT a boundary: H₁ class survives

    # The Z₁ basis of T\\v in Ω₁(T\\v) coords
    # H₁ generator = any z1 element not in B₁(T\\v)
    # For β₁=1: z1 has dim z1_dim, B₁ has dim b1_rk = z1_dim - 1
    z1_dim_sub = z1_sub.shape[0]
    b1_rk_sub = data_sub['b1_rk']

    # Find the H₁ generator: a Z₁ element not in B₁
    bd2_sub = data_sub['bd2_om']
    if bd2_sub.shape[1] > 0:
        # B₁ in Z₁ coords: project bd2_sub columns into Z₁
        b1_in_z1 = z1_sub @ bd2_sub  # z1_dim × d2
        # The H₁ generator is orthogonal to B₁ within Z₁
        U_b, S_b, _ = np.linalg.svd(b1_in_z1, full_matrices=True)
        b_rk = sum(s > 1e-8 for s in S_b)
        h1_gen = U_b[:, b_rk:]  # z1_dim × 1 (the H₁ generator in Z₁ coords)
    else:
        h1_gen = np.eye(z1_dim_sub)  # all of Z₁ is H₁

    if h1_gen.shape[1] == 0:
        return (0, beta1_sub, None)

    # Map H₁ generator to T: first to Ω₁(T\\v) coords, then to Ω₁(T)
    gen_sub = z1_sub.T @ h1_gen[:, 0]  # d1_sub vector
    gen_T = incl @ gen_sub  # d1_T vector (in Ω₁(T) coords)

    # Check: is gen_T a boundary in T?
    bd2_T = data_T['bd2_om']
    if bd2_T.shape[1] > 0:
        # Is gen_T in im(bd2_T)?
        resid = gen_T - bd2_T @ np.linalg.lstsq(bd2_T, gen_T, rcond=None)[0]
        resid_norm = np.linalg.norm(resid)
    else:
        resid_norm = np.linalg.norm(gen_T)

    # Also check: is gen_T in Z₁(T)?
    bd1_T = build_full_boundary_matrix(data_T['ap1'], enumerate_allowed_paths(A, n, 0))
    bd1_check = bd1_T @ om1_T @ gen_T
    z1_check = np.linalg.norm(bd1_check)

    survives = resid_norm > 1e-6
    is_cycle = z1_check < 1e-6

    return (1 if survives else 0, beta1_sub, survives, is_cycle)


print("=" * 70)
print("H₁ SURVIVAL UNDER INCLUSION")
print("=" * 70)

n = 5
m = n*(n-1)//2
total = 1 << m

survival_count_dist = Counter()  # #vertices with surviving H₁
survival_vs_h2sum = defaultdict(Counter)

t0 = time.time()
for bits in range(total):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    surviving = 0
    for v in range(n):
        if scores[v] == 0 or scores[v] == n-1:
            continue
        result = check_h1_survival(A, n, v)
        if result[1] > 0 and len(result) > 2 and result[2]:
            surviving += 1

    survival_count_dist[surviving] += 1

elapsed = time.time() - t0
print(f"\nn={n}: {total} tournaments in {elapsed:.0f}s")
print(f"  #(interior v with surviving H₁(T\\v)): {dict(sorted(survival_count_dist.items()))}")

# Now: compare survival count with Σ h₂_rel
# They should be EQUAL: h₂_rel contribution to β₂ iff H₁ survives
print(f"\n--- Comparing survival count with Σ h₂_rel ---")

# Re-run with h₂_rel computation
from beta2_h2rel_sum import compute_h2_rel

mismatch = 0
for bits in range(total):
    A = build_adj(n, bits)
    scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

    surviving = 0
    h2_sum = 0
    for v in range(n):
        h2r = compute_h2_rel(A, n, v)
        h2_sum += h2r

        if scores[v] == 0 or scores[v] == n-1:
            continue
        result = check_h1_survival(A, n, v)
        if result[1] > 0 and len(result) > 2 and result[2]:
            surviving += 1

    if surviving != h2_sum:
        mismatch += 1
        if mismatch <= 5:
            print(f"  MISMATCH bits={bits}: surviving={surviving}, Σ h₂_rel={h2_sum}")

if mismatch == 0:
    print(f"  ✓ PERFECT MATCH: #surviving = Σ h₂_rel for all {total} tournaments")
else:
    print(f"  ✗ {mismatch} mismatches")

print("\nDone.")
