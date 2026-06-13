#!/usr/bin/env python3
"""
beta2_z2_tt_only.py — Is Z₂ always spanned by TT (transitive triple) paths?

Key observation from multicone analysis: Z₂ basis vectors seem to consist
only of TT paths. If true universally, this is a major structural constraint.

Also test: is Ω₃ always spanned by "doubly transitive" 3-paths?
A doubly-transitive 3-path (a,b,c,d) has a→c, a→d, b→d (all transitive).

Author: opus-2026-03-08-S49
"""
import sys, time
import numpy as np
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


print("=" * 70)
print("Z₂ ⊆ span(TT)? UNIVERSAL TEST")
print("=" * 70)

for n in [4, 5, 6]:
    m = n*(n-1)//2
    total = 1 << m
    t0 = time.time()

    z2_has_nt = 0
    z2_pure_tt = 0
    z2_zero = 0
    omega2_has_nt = 0

    for bits in range(total):
        A = build_adj(n, bits)

        ap0 = enumerate_allowed_paths(A, n, 0)
        ap1 = enumerate_allowed_paths(A, n, 1)
        ap2 = enumerate_allowed_paths(A, n, 2)
        ap3 = enumerate_allowed_paths(A, n, 3)
        om1 = compute_omega_basis(A, n, 1, ap1, ap0)
        om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))

        d2 = dim_om(om2)
        if d2 == 0:
            z2_zero += 1
            continue

        # Check if Ω₂ has NT paths
        for j in range(d2):
            col = om2[:, j]
            for i, p in enumerate(ap2):
                if abs(col[i]) > 1e-8:
                    if not A[p[0]][p[2]]:  # NT: c→a, not a→c
                        omega2_has_nt += 1
                        break
            else:
                continue
            break

        # Compute Z₂
        bd2 = build_full_boundary_matrix(ap2, ap1)
        d2_mat = np.linalg.lstsq(om1, bd2 @ om2, rcond=None)[0]
        U, S, Vt = np.linalg.svd(d2_mat, full_matrices=True)
        rk = sum(s > 1e-8 for s in S)
        z2_dim = d2 - rk

        if z2_dim == 0:
            z2_zero += 1
            continue

        # Z₂ in A₂ coords
        z2_om = Vt[rk:, :]
        z2_A2 = om2 @ z2_om.T  # |A₂| × z2_dim

        # Check for NT paths in Z₂
        has_nt = False
        for j in range(z2_dim):
            z = z2_A2[:, j]
            for i, p in enumerate(ap2):
                if abs(z[i]) > 1e-8:
                    if not A[p[0]][p[2]]:  # NT
                        has_nt = True
                        break
            if has_nt:
                break

        if has_nt:
            z2_has_nt += 1
        else:
            z2_pure_tt += 1

    elapsed = time.time() - t0
    print(f"\nn={n}: {total} tournaments in {elapsed:.0f}s")
    print(f"  Z₂=0: {z2_zero}")
    print(f"  Z₂ pure TT: {z2_pure_tt}")
    print(f"  Z₂ has NT: {z2_has_nt}")
    print(f"  Ω₂ has NT component: {omega2_has_nt}")

    if z2_has_nt == 0 and z2_pure_tt > 0:
        print(f"  *** Z₂ ⊆ span(TT) CONFIRMED for all n={n} tournaments! ***")

# Now check Ω₃ and "doubly transitive" paths
print("\n" + "=" * 70)
print("Ω₃ STRUCTURE: DOUBLY TRANSITIVE?")
print("=" * 70)

n = 5
m = n*(n-1)//2
total = 1 << m

dt_count = 0
ndt_count = 0

for bits in range(total):
    A = build_adj(n, bits)

    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

    d3 = dim_om(om3)
    if d3 == 0:
        continue

    # Check each Ω₃ basis vector
    has_non_dt = False
    for j in range(d3):
        col = om3[:, j]
        for i, p in enumerate(ap3):
            if abs(col[i]) > 1e-8:
                # Check doubly transitive: a→c, a→d, b→d
                a, b, c, d = p[0], p[1], p[2], p[3]
                is_dt = A[a][c] and A[a][d] and A[b][d]
                if not is_dt:
                    has_non_dt = True
                    break
        if has_non_dt:
            break

    if has_non_dt:
        ndt_count += 1
    else:
        dt_count += 1

print(f"\nn={n}: Ω₃ doubly-transitive: {dt_count}, has non-DT: {ndt_count}")

# Now: given Z₂ is pure TT, and every TT path (a,b,c) has a→b→c and a→c,
# what are the constraints on z ∈ Z₂?
print("\n" + "=" * 70)
print("Z₂ CONSTRAINTS FOR PURE TT")
print("=" * 70)

n = 5
print(f"\nFor a TT path (a,b,c) with a→b→c, a→c:")
print(f"  ∂₂(a,b,c) = (b,c) - (a,c) + (a,b)")
print(f"  All three are edges in A₁.")
print(f"  For z = Σ α_{{abc}} (a,b,c) ∈ Z₂: Σ α_{{abc}} ∂₂(a,b,c) = 0 in Ω₁")
print(f"  This means: for each edge (i,j) ∈ A₁:")
print(f"    Σ_{{c: (i,j,c) TT}} α_{{ijc}} - Σ_{{a: (a,i,j) TT}} α_{{aij}} + Σ_{{b: (i,b,j) TT}} α_{{ibj}} = 0")

# But wait - the last term needs a→b→c format with a=i, c=j
# ∂₂(a,b,c) has edge (a,b) in position 2 (with coeff +1), (a,c) in position 1 (coeff -1), (b,c) in position 0 (coeff +1)
# So for edge (i,j):
# Contribution from (a,b,c) where (b,c)=(i,j): coeff = +α_{a,i,j}
# Contribution from (a,b,c) where (a,c)=(i,j): coeff = -α_{i,b,j}
# Contribution from (a,b,c) where (a,b)=(i,j): coeff = +α_{i,j,c}

# So: for edge (i,j), the coefficient in ∂₂(z) is:
# Σ_{a: (a,i,j) TT} α_{aij} + Σ_{c: (i,j,c) TT} α_{ijc} - Σ_{b: (i,b,j) TT} α_{ibj} = 0

# Wait, let me be more careful. ∂₂(a,b,c) = (b,c) - (a,c) + (a,b).
# For edge e=(i,j) ∈ A₁:
# (b,c)=(i,j) means b=i,c=j, path is (a,i,j): coeff = +α_{a,i,j}
# (a,c)=(i,j) means a=i,c=j, path is (i,b,j): coeff = -α_{i,b,j}
# (a,b)=(i,j) means a=i,b=j, path is (i,j,c): coeff = +α_{i,j,c}

# So: Σ_a α_{aij} - Σ_b α_{ibj} + Σ_c α_{ijc} = 0

print(f"\n  Correct constraint for edge (i,j):")
print(f"    Σ_a α_{{a,i,j}} - Σ_b α_{{i,b,j}} + Σ_c α_{{i,j,c}} = 0")
print(f"    (sum over TT paths ending at ij) - (sum over TT paths with ij at ends) + (sum over TT paths starting at ij)")

# Count: how many TT paths exist for a typical n=5 tournament?
A = build_adj(5, 0)  # transitive
tt_count = sum(1 for p in enumerate_allowed_paths(A, 5, 2) if A[p[0]][p[2]])
print(f"\n  Transitive T₅: {tt_count} TT paths out of {len(enumerate_allowed_paths(A, 5, 2))} allowed 2-paths")

A = build_adj(5, 76)  # first regular
ap2 = enumerate_allowed_paths(A, 5, 2)
tt_count = sum(1 for p in ap2 if A[p[0]][p[2]])
nt_count = sum(1 for p in ap2 if not A[p[0]][p[2]])
print(f"  Regular T₅: {tt_count} TT, {nt_count} NT out of {len(ap2)} allowed 2-paths")

print("\nDone.")
