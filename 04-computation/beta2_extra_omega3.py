#!/usr/bin/env python3
"""
beta2_extra_omega3.py — Investigate Ω₃ elements beyond DT+cancel

At n=5: 144 tournaments have Ω₃ ⊋ DT+cancel.
Questions:
1. What are these extra elements? Are they 3-bad-face combinations?
2. Do they contribute to im(∂₃)? Or do they all go to Z₃ (kernel)?
3. Is DT+cancel already sufficient for filling Z₂ (even when Ω₃ is larger)?

NEW APPROACH: Think of the Ω₃ constraint as a LINEAR SYSTEM.
Ω₃ = ker(J₃) where J₃: R^{A₃} → R^{bad faces} is the "junk projection."
J₃ has columns indexed by A₃ paths and rows indexed by non-A₂ faces.

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


print("=" * 70)
print("EXTRA Ω₃ ELEMENTS BEYOND DT+CANCEL")
print("=" * 70)

n = 5
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

extra_cases = []

for bits in range(1 << m):
    A = build_adj(n, bits)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap0 = enumerate_allowed_paths(A, n, 0)

    if not ap3: continue

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    om3 = compute_omega_basis(A, n, 3, ap3, ap2) if ap3 else np.zeros((0,0))

    d2 = dim_om(om2)
    d3 = dim_om(om3)
    if d3 == 0 or d2 == 0: continue

    # Build DT + cancel vectors
    dt_vecs = []
    bad_groups = defaultdict(list)

    for i, p in enumerate(ap3):
        a, b, c, d = p
        if A[a][c] and A[b][d]:
            v = np.zeros(len(ap3)); v[i] = 1
            dt_vecs.append(v)
        else:
            if not A[a][c]: bad_groups[('02', a, c)].append(i)
            if not A[b][d]: bad_groups[('13', b, d)].append(i)

    cancel_vecs = []
    for key, indices in bad_groups.items():
        for j in range(1, len(indices)):
            v = np.zeros(len(ap3))
            v[indices[0]] = 1; v[indices[j]] = -1
            cancel_vecs.append(v)

    all_vecs = dt_vecs + cancel_vecs
    if not all_vecs: continue

    V = np.column_stack(all_vecs)
    coeffs = np.linalg.lstsq(V, om3, rcond=None)[0]
    residual = np.max(np.abs(om3 - V @ coeffs))

    if residual > 1e-8:
        scores = tuple(sorted(sum(A[i][j] for j in range(n) if j!=i) for i in range(n)))
        c3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
                 if A[i][j]+A[j][k]+A[k][i] in (0,3))

        # Find an extra Ω₃ element
        # It's in ker(J₃) but not in span(dt_vecs + cancel_vecs)
        # First, find it as an Ω₃ element expressed in A₃ coordinates

        # Project out the DT+cancel components from om3
        proj_residual = om3 - V @ coeffs

        # Find a nonzero column
        for col in range(d3):
            if np.max(np.abs(proj_residual[:, col])) > 1e-8:
                extra_vec = proj_residual[:, col]
                # Normalize
                extra_vec = extra_vec / np.max(np.abs(extra_vec))
                break

        # Which A₃ paths are involved?
        involved = []
        for i in range(len(ap3)):
            if abs(extra_vec[i]) > 1e-8:
                p = ap3[i]
                a, b, c, d = p
                bad02 = not A[a][c]
                bad13 = not A[b][d]
                involved.append((tuple(p), round(extra_vec[i], 4), bad02, bad13))

        # Compute ∂₃ image of this extra element
        bd3 = build_full_boundary_matrix(ap3, ap2)
        bd_image = bd3 @ extra_vec
        bd_norm = np.linalg.norm(bd_image)

        # Is this extra element in Z₃ (kernel of ∂₃)?
        in_z3 = bd_norm < 1e-8

        # Compute Z₂ and check filling
        bd2 = build_full_boundary_matrix(ap2, ap1)
        bd2_om = bd2 @ om2
        coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
        rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
        z2_dim = d2 - rk2

        # DT filling rank
        bd3_full = build_full_boundary_matrix(ap3, ap2)
        if dt_vecs:
            DT_mat = np.column_stack(dt_vecs)
            bd3_dt = bd3_full @ DT_mat
            coords_dt = np.linalg.lstsq(om2, bd3_dt, rcond=None)[0]
            rk_dt = np.linalg.matrix_rank(coords_dt, tol=1e-8)
        else:
            rk_dt = 0

        # DT+cancel filling rank
        ALL_mat = np.column_stack(all_vecs)
        bd3_all = bd3_full @ ALL_mat
        coords_all = np.linalg.lstsq(om2, bd3_all, rcond=None)[0]
        rk_all = np.linalg.matrix_rank(coords_all, tol=1e-8)

        # Full filling rank
        bd3_om = bd3_full @ om3
        coords3 = np.linalg.lstsq(om2, bd3_om, rcond=None)[0]
        rk3 = np.linalg.matrix_rank(coords3, tol=1e-8)

        extra_cases.append({
            'bits': bits, 'scores': scores, 'c3': c3,
            'd3': d3, 'n_dt': len(dt_vecs), 'n_cancel': len(cancel_vecs),
            'z2_dim': z2_dim, 'rk_dt': rk_dt, 'rk_all': rk_all, 'rk3': rk3,
            'in_z3': in_z3, 'involved': involved,
        })

print(f"\nTournaments with extra Ω₃ elements: {len(extra_cases)}")

# Summary
extra_in_z3 = sum(1 for e in extra_cases if e['in_z3'])
extra_not_z3 = sum(1 for e in extra_cases if not e['in_z3'])
print(f"  Extra elements in Z₃ (kernel of ∂₃): {extra_in_z3}")
print(f"  Extra elements NOT in Z₃: {extra_not_z3}")

# Does DT+cancel already fill Z₂?
dt_cancel_fills = sum(1 for e in extra_cases if e['rk_all'] >= e['z2_dim'])
dt_cancel_fails = sum(1 for e in extra_cases if e['rk_all'] < e['z2_dim'])
print(f"\n  DT+cancel fills Z₂: {dt_cancel_fills}/{len(extra_cases)}")
print(f"  DT+cancel FAILS to fill Z₂: {dt_cancel_fails}/{len(extra_cases)}")

# Detailed examples
print(f"\nDetailed examples:")
for e in extra_cases[:5]:
    print(f"\n  T#{e['bits']} scores={e['scores']}, c₃={e['c3']}")
    print(f"    Ω₃={e['d3']}, DT={e['n_dt']}, cancel={e['n_cancel']}")
    print(f"    Z₂={e['z2_dim']}, rk_DT={e['rk_dt']}, rk_DT+cancel={e['rk_all']}, rk(∂₃)={e['rk3']}")
    print(f"    Extra in Z₃? {e['in_z3']}")
    print(f"    Extra element paths:")
    for path, coeff, bad02, bad13 in e['involved']:
        bad_str = ""
        if bad02: bad_str += " BAD02"
        if bad13: bad_str += " BAD13"
        if not bad02 and not bad13: bad_str = " DT"
        print(f"      {coeff:+.4f} * {path}{bad_str}")

# Score distribution
score_dist = Counter(e['scores'] for e in extra_cases)
print(f"\nScore distribution of extra cases:")
for scores, count in sorted(score_dist.items()):
    print(f"  {scores}: {count}")

# KEY: Check if DT alone fills Z₂ when extras are in Z₃
print(f"\n{'='*70}")
print("KEY FINDING: DO EXTRAS MATTER FOR FILLING?")
print("=" * 70)

# If all extra elements are in Z₃, then they don't contribute to im(∂₃).
# So DT+cancel should suffice for filling.
# But if DT+cancel doesn't fill Z₂, then the extras (not in Z₃) are needed.

for e in extra_cases:
    if e['rk_all'] < e['z2_dim']:
        print(f"\n  EXTRA NEEDED: T#{e['bits']} scores={e['scores']}")
        print(f"    Z₂={e['z2_dim']}, rk_DT+cancel={e['rk_all']}, rk(∂₃)={e['rk3']}")
        print(f"    Extra in Z₃? {e['in_z3']}")

if dt_cancel_fails == 0 and extra_in_z3 == len(extra_cases):
    print("\n✓ ALL extra Ω₃ elements are in Z₃ (contribute nothing to im(∂₃))")
    print("✓ DT+cancel ALWAYS suffices for filling Z₂")
    print("\n⟹ PROOF STRATEGY: Show that DT+cancel spans a subspace")
    print("   whose ∂₃-image contains Z₂.")

# ================================================================
# PART 2: Understand the DT-filling deficit more precisely
# ================================================================
print(f"\n{'='*70}")
print("PART 2: DT-FILLING DEFICIT STRUCTURE")
print("=" * 70)

# At n=5, DT alone fills ALL Z₂ (rk_dt = z2_dim always).
# At n=6, there's a deficit for 960 tournaments.
# Let's understand what happens at n=6.

print("\nAt n=5: checking DT filling...")
n = 5
dt_ok = 0
dt_fail = 0
for bits in range(1 << (n*(n-1)//2)):
    A = build_adj(n, bits)
    ap1 = enumerate_allowed_paths(A, n, 1)
    ap2 = enumerate_allowed_paths(A, n, 2)
    ap3 = enumerate_allowed_paths(A, n, 3)
    ap0 = enumerate_allowed_paths(A, n, 0)

    om1 = compute_omega_basis(A, n, 1, ap1, ap0)
    om2 = compute_omega_basis(A, n, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    d2 = dim_om(om2)
    if d2 == 0 or not ap3: continue

    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2_dim = d2 - rk2
    if z2_dim == 0: continue

    # DT filling
    dt_idx = [i for i, p in enumerate(ap3) if A[p[0]][p[2]] and A[p[1]][p[3]]]
    bd3 = build_full_boundary_matrix(ap3, ap2)
    if dt_idx:
        bd3_dt = np.column_stack([bd3[:, i] for i in dt_idx])
        coords_dt = np.linalg.lstsq(om2, bd3_dt, rcond=None)[0]
        rk_dt = np.linalg.matrix_rank(coords_dt, tol=1e-8)
    else:
        rk_dt = 0

    if rk_dt >= z2_dim:
        dt_ok += 1
    else:
        dt_fail += 1

print(f"  DT alone fills Z₂: {dt_ok}/{dt_ok+dt_fail}")
print(f"  DT alone fails: {dt_fail}/{dt_ok+dt_fail}")

# ================================================================
# PART 3: CRUCIAL — What is the STRUCTURE of im(∂₃|_DT)?
# ================================================================
print(f"\n{'='*70}")
print("PART 3: BOUNDARY IMAGE OF DT PATHS")
print("=" * 70)

# For a DT 4-path (a,b,c,d) with a→c, b→d:
# ∂₃(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)
# All faces in A₂.
# (a,b,c) and (b,c,d) are TT.
# If a→d: (a,c,d) TT, (a,b,d) TT. All TT.
# If d→a: (a,c,d) NT (cycle a→c→d→a), (a,b,d) NT (cycle a→b→d→a).
#
# In Ω₂ coords, the NT paths are "absorbed" by the quotient.
# So ∂₃ of a DT path, viewed in Ω₂, is a combination of TT basis elements.
#
# QUESTION: What is the STRUCTURE of this combination?

print("\nBoundary of DT paths in Ω₂ coordinates:")
print("(Expressing ∂₃ of each DT path as a sum of Ω₂ basis elements)")

# Pick a specific tournament and analyze
for test_bits in [0, 1, 5, 100]:
    A = build_adj(5, test_bits)
    scores = tuple(sorted(sum(A[i][j] for j in range(5) if j!=i) for i in range(5)))

    ap1 = enumerate_allowed_paths(A, 5, 1)
    ap2 = enumerate_allowed_paths(A, 5, 2)
    ap3 = enumerate_allowed_paths(A, 5, 3)
    ap0 = enumerate_allowed_paths(A, 5, 0)

    om1 = compute_omega_basis(A, 5, 1, ap1, ap0)
    om2 = compute_omega_basis(A, 5, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    d2 = dim_om(om2)
    if d2 == 0: continue

    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2_dim = d2 - rk2

    dt_idx = [i for i, p in enumerate(ap3) if A[p[0]][p[2]] and A[p[1]][p[3]]]
    bd3 = build_full_boundary_matrix(ap3, ap2)

    print(f"\n  T#{test_bits} scores={scores}: dim(Ω₂)={d2}, dim(Z₂)={z2_dim}, |DT|={len(dt_idx)}")

    for j, idx in enumerate(dt_idx):
        path = tuple(ap3[idx])
        a, b, c, d = path
        bd_vec = bd3[:, idx]  # boundary in A₂ coords

        # Project to Ω₂ coords
        om2_coords = np.linalg.lstsq(om2, bd_vec, rcond=None)[0]

        # Which Ω₂ basis elements are involved?
        nonzero_om2 = [(k, round(om2_coords[k], 4)) for k in range(d2) if abs(om2_coords[k]) > 1e-8]

        # Classify the a↔d arc
        ad_type = "a→d" if A[a][d] else "d→a"

        print(f"    DT({a},{b},{c},{d}) [{ad_type}]: {len(nonzero_om2)} Ω₂ components, "
              f"coeffs = {[c for _, c in nonzero_om2]}")

# ================================================================
# PART 4: THE BIPARTITE GRAPH: DT paths ↔ Z₂ elements
# ================================================================
print(f"\n{'='*70}")
print("PART 4: BIPARTITE MATCHING PERSPECTIVE")
print("=" * 70)

# Think of DT paths as "fillers" and Z₂ basis elements as "targets."
# ∂₃|_DT: DT → Ω₂ gives a bipartite graph between DT and Z₂.
# β₂=0 iff the image of this map covers all of Z₂.
# By Hall's theorem, this holds iff for every subset S ⊆ Z₂,
# the neighborhood N(S) in DT satisfies |N(S)| ≥ |S|.
#
# This is a rank condition, not literally Hall's theorem (it's over R, not Z).
# But the spirit is: enough DT paths must "reach" each Z₂ direction.

print("\nFor each n=5 tournament: column rank analysis of ∂₃|_DT → Z₂")

interesting = []
for bits in range(1 << 10):
    A = build_adj(5, bits)
    ap1 = enumerate_allowed_paths(A, 5, 1)
    ap2 = enumerate_allowed_paths(A, 5, 2)
    ap3 = enumerate_allowed_paths(A, 5, 3)
    ap0 = enumerate_allowed_paths(A, 5, 0)

    om1 = compute_omega_basis(A, 5, 1, ap1, ap0)
    om2 = compute_omega_basis(A, 5, 2, ap2, ap1) if ap2 else np.zeros((0,0))
    d2 = dim_om(om2)
    if d2 == 0 or not ap3: continue

    bd2 = build_full_boundary_matrix(ap2, ap1)
    bd2_om = bd2 @ om2
    coords2 = np.linalg.lstsq(om1, bd2_om, rcond=None)[0]
    rk2 = np.linalg.matrix_rank(coords2, tol=1e-8)
    z2_dim = d2 - rk2
    if z2_dim == 0: continue

    z2_basis = np.linalg.svd(coords2, full_matrices=True)[2][rk2:].T

    dt_idx = [i for i, p in enumerate(ap3) if A[p[0]][p[2]] and A[p[1]][p[3]]]
    bd3 = build_full_boundary_matrix(ap3, ap2)

    # ∂₃|_DT in Ω₂ coords
    if dt_idx:
        bd3_dt = np.column_stack([bd3[:, i] for i in dt_idx])
        om2_coords = np.linalg.lstsq(om2, bd3_dt, rcond=None)[0]
        # Project to Z₂
        z2_proj = z2_basis.T @ om2_coords  # z2_dim × |DT|
        rk_z2_proj = np.linalg.matrix_rank(z2_proj, tol=1e-8)
    else:
        rk_z2_proj = 0

    if rk_z2_proj == z2_dim:
        pass  # DT fills Z₂
    else:
        interesting.append((bits, z2_dim, rk_z2_proj, len(dt_idx)))

if not interesting:
    print("  ✓ At n=5, DT always fills Z₂ (confirmed for first 1024 tournaments)")
else:
    print(f"  FAILURES: {len(interesting)}")
    for bits, z2, rk, n_dt in interesting[:5]:
        print(f"    bits={bits}: Z₂={z2}, rk(proj)={rk}, |DT|={n_dt}")


print("\nDone.")
