#!/usr/bin/env python3
"""
beta2_relative_anatomy.py — What IS H₂(T, T\\v)?

Understand the relative homology H₂(T, T\\v) = ker(∂₂^rel) / im(∂₃^rel).

The relative chain complex Ω_p(T, T\\v) = Ω_p(T) / Ω_p(T\\v) consists of
"classes of Ω-chains that essentially use vertex v."

For the connecting map δ: H₂(T,T\\v) → H₁(T\\v):
- Take a relative 2-cycle z (uses v, ∂₂z lands in T\\v chains)
- δ(z) = [∂₂z̃] where z̃ is a lift of z

δ INJECTIVE means: if ∂₂z̃ is a boundary in T\\v, then z is a relative boundary.

Question: When exactly is δ NOT injective? What goes wrong at source/sink?

Author: opus-2026-03-08-S49
"""
import sys
import numpy as np
from collections import Counter
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
print("ANATOMY OF RELATIVE H₂ AND THE CONNECTING MAP δ")
print("=" * 70)

# Study the failure case: T#36, v=4 (source)
n = 5
bits = 36
A = build_adj(n, bits)
v = 4

scores = [sum(A[i][j] for j in range(n) if j != i) for i in range(n)]
print(f"\nT#{bits} scores={scores}")
print(f"v={v}, d⁺(v)={scores[v]}")

# What paths use v?
ap2 = enumerate_allowed_paths(A, n, 2)
ap3 = enumerate_allowed_paths(A, n, 3)

ap2_list = [tuple(p) for p in ap2]
ap3_list = [tuple(p) for p in ap3]

ap2_with_v = [p for p in ap2_list if v in p]
ap2_without_v = [p for p in ap2_list if v not in p]
ap3_with_v = [p for p in ap3_list if v in p]
ap3_without_v = [p for p in ap3_list if v not in p]

print(f"\n|A₂|={len(ap2_list)}: {len(ap2_with_v)} use v, {len(ap2_without_v)} don't")
print(f"|A₃|={len(ap3_list)}: {len(ap3_with_v)} use v, {len(ap3_without_v)} don't")

print(f"\nA₂ paths using v:")
for p in ap2_with_v:
    a, b, c = p
    tt = "TT" if A[a][c] else "NT"
    pos = [i for i, x in enumerate(p) if x == v][0]
    print(f"  {p} [{tt}] v at position {pos}")

print(f"\nA₃ paths using v:")
for p in ap3_with_v:
    a, b, c, d = p
    dt = "DT" if (A[a][c] and A[b][d]) else "non-DT"
    pos = [i for i, x in enumerate(p) if x == v][0]
    print(f"  {p} [{dt}] v at position {pos}")

# Since v=4 is a source (beats everyone), it can only appear at position 0.
# Because if v were at position 1, we'd need something → v, but nothing beats v.
# So ALL paths using v start with v.

# This means:
# A₂ using v: (v, *, *) with v→* and *→*
# A₃ using v: (v, *, *, *) with v→* etc.

# For the relative boundary ∂₂^rel:
# ∂₂(v, b, c) = (b, c) - (v, c) + (v, b)
# In the relative complex, (b,c) ∈ Ω₁(T\v), so it maps to 0 in the quotient.
# The relative part is: -(v, c) + (v, b) which uses v.
# Wait, the quotient Ω₁(T)/Ω₁(T\v) consists of arcs involving v.

# For a source v: arcs involving v are all v→x (since nobody beats v).
# Arcs NOT involving v form Ω₁(T\v).
# So relative Ω₁ = span of {(v, x) : v→x} = R^{n-1}.

# Relative ∂₂(v, b, c) projected to relative Ω₁:
# = -(v,c) + (v,b) (the faces involving v)

# For the connecting map:
# ∂₂(v, b, c) = (b,c) - (v,c) + (v,b) (in full Ω₁(T))
# The relative part: -(v,c) + (v,b) (in Ω₁(T)/Ω₁(T\v))
# The connecting map sends: [∂₂(lift)] mapped to Ω₁(T\v) part = (b,c)

# So δ of a relative 2-cycle sends it to [the (b,c) parts of ∂₂(v,b,c)].
# These are arcs in T\v.

# If v is a source, the paths using v are (v, *, *) and the non-v face is (* ,*).
# The connecting map sends H₂(T, T\v) → H₁(T\v) via the non-v faces.

# For the source case, the relative 2-cycles are combinations of (v, b, c) paths
# whose ∂₂ in relative Ω₁ vanishes.
# Relative ∂₂(v, b, c) in relative Ω₁ = -(v,c) + (v,b).
# Cycle condition: Σ α_{bc} (-(v,c) + (v,b)) = 0 in relative Ω₁.
# That means Σ_{c} (Σ_{b} α_{bc}) (v,c) = Σ_{b} (Σ_{c} α_{bc}) (v,b)
# In terms of arcs (v, x) for x ≠ v:
# coefficient of (v, x) = Σ_c α_{xc} - Σ_b α_{bx} = 0 for all x.

# This is a "flow conservation" condition at each x: in-flow = out-flow.

# The connecting map sends this to [non-v parts]: Σ α_{bc} (b,c) ∈ H₁(T\v).
# If all these (b,c) are boundaries in T\v, then δ = 0.

# When β₁(T\v) = 1, H₁(T\v) ≠ 0, so the image could be nonzero OR could land
# in B₁. In the source case, the image apparently lands in B₁ always.

# WHY? Because the flow conservation condition from v gives a "balanced" flow
# on T\v, which is necessarily a boundary (balanced flows on graphs are boundaries).
# Wait, is this true?

# Let me understand T\v better
others = [i for i in range(n) if i != v]
n1 = n - 1
A_sub = [[A[others[i]][others[j]] for j in range(n1)] for i in range(n1)]
scores_sub = [sum(A_sub[i][j] for j in range(n1) if j!=i) for i in range(n1)]

print(f"\nT\\v (vertices {others}):")
print(f"  scores = {scores_sub}")
for i in range(n1):
    row = [A_sub[i][j] for j in range(n1)]
    print(f"  {others[i]}: {row}")

# Check β₁(T\v)
ap0_sub = enumerate_allowed_paths(A_sub, n1, 0)
ap1_sub = enumerate_allowed_paths(A_sub, n1, 1)
ap2_sub = enumerate_allowed_paths(A_sub, n1, 2)

om1_sub = compute_omega_basis(A_sub, n1, 1, ap1_sub, ap0_sub)
om2_sub = compute_omega_basis(A_sub, n1, 2, ap2_sub, ap1_sub) if ap2_sub else np.zeros((0,0))

bd1_sub = build_full_boundary_matrix(ap1_sub, ap0_sub)
rk_d1 = np.linalg.matrix_rank(bd1_sub, tol=1e-8)
z1 = dim_om(om1_sub) - rk_d1

if dim_om(om2_sub) > 0:
    bd2_sub = build_full_boundary_matrix(ap2_sub, ap1_sub)
    coords2_sub = np.linalg.lstsq(om1_sub, bd2_sub @ om2_sub, rcond=None)[0]
    rk_d2 = np.linalg.matrix_rank(coords2_sub, tol=1e-8)
else:
    rk_d2 = 0

beta1 = z1 - rk_d2
print(f"\n  β₁(T\\v) = {beta1}")

# Find the H₁ cycle
if beta1 > 0:
    U1, S1, Vt1 = np.linalg.svd(bd1_sub @ om1_sub, full_matrices=True)
    z1_basis = Vt1[rk_d1:]  # rows

    # Find which Z₁ directions are NOT in im(∂₂)
    if rk_d2 > 0:
        proj = z1_basis @ coords2_sub
        U2, S2, _ = np.linalg.svd(proj, full_matrices=True)
        rk_proj = sum(s > 1e-8 for s in S2)
        h1_dirs = z1_basis[rk_proj:]  # in Ω₁ coords
    else:
        h1_dirs = z1_basis

    print(f"  H₁ generators (in Ω₁(T\\v) coords):")
    ap1_sub_list = [tuple(p) for p in ap1_sub]
    for k in range(h1_dirs.shape[0]):
        h1_vec = om1_sub @ h1_dirs[k]
        print(f"    Cycle {k}:")
        for j in range(len(ap1_sub_list)):
            if abs(h1_vec[j]) > 1e-8:
                a_orig = others[ap1_sub_list[j][0]]
                b_orig = others[ap1_sub_list[j][1]]
                print(f"      {h1_vec[j]:+.4f} * ({a_orig},{b_orig})")


# =============================================================
# Now contrast with an interior vertex
# =============================================================
print(f"\n{'='*70}")
print("CONTRAST: Interior vertex v=0 (d⁺=1)")
print("-" * 50)

v2 = 0
dv2 = scores[v2]
print(f"v={v2}, d⁺(v)={dv2}")

ap2_with_v2 = [p for p in ap2_list if v2 in p]
ap3_with_v2 = [p for p in ap3_list if v2 in p]

print(f"  A₂ paths using v={v2}: {len(ap2_with_v2)}")
for p in ap2_with_v2:
    a, b, c = p
    tt = "TT" if A[a][c] else "NT"
    pos = [i for i, x in enumerate(p) if x == v2][0]
    print(f"    {p} [{tt}] v at position {pos}")

print(f"\n  A₃ paths using v={v2}: {len(ap3_with_v2)}")
for p in ap3_with_v2:
    a, b, c, d = p
    dt = "DT" if (A[a][c] and A[b][d]) else "non-DT"
    pos = [i for i, x in enumerate(p) if x == v2][0]
    print(f"    {p} [{dt}] v at position {pos}")

# When v is interior:
# Paths using v can have v at position 0, 1, or 2 (for A₂).
# This means the relative chain complex has more structure.
# The connecting map sends more diverse faces to T\v.


# =============================================================
# KEY INSIGHT: When v is source, v can only be at position 0.
# =============================================================
print(f"\n{'='*70}")
print("KEY OBSERVATION: Position of v in paths")
print("-" * 50)

# For all n=5 tournaments, count position distributions
n = 5
m = n*(n-1)//2

for v_type, v_scores in [("source", [n-1]), ("sink", [0]), ("interior", list(range(1, n-1)))]:
    pos2_dist = Counter()
    pos3_dist = Counter()
    total_tours = 0

    for bits in range(1 << m):
        A = build_adj(n, bits)
        scores = [sum(A[i][j] for j in range(n) if j!=i) for i in range(n)]

        for v in range(n):
            if scores[v] not in v_scores:
                continue

            total_tours += 1
            ap2 = enumerate_allowed_paths(A, n, 2)
            ap3 = enumerate_allowed_paths(A, n, 3)

            for p in ap2:
                if v in tuple(p):
                    pos = list(p).index(v)
                    pos2_dist[pos] += 1

            for p in ap3:
                if v in tuple(p):
                    pos = list(p).index(v)
                    pos3_dist[pos] += 1

            break  # only one v per tournament for efficiency

    if total_tours > 0:
        print(f"\n  {v_type} vertex (d⁺ ∈ {v_scores}), {total_tours} (T,v) pairs:")
        print(f"    A₂ position: {dict(sorted(pos2_dist.items()))}")
        print(f"    A₃ position: {dict(sorted(pos3_dist.items()))}")


print("\nDone.")
