#!/usr/bin/env python3
"""
beta2_mayer_vietoris.py — Mayer-Vietoris decomposition for arc flips

When flipping u→v to v→u, the chain complex changes:
  C_p(T)  = C_p(shared) ⊕ C_p(arc_uv)
  C_p(T') = C_p(shared) ⊕ C_p(arc_vu)

where C_p(shared) = paths not using the flipped arc,
      C_p(arc_uv) = paths using u→v,
      C_p(arc_vu) = paths using v→u.

Key question: what algebraic relationship between these pieces
forces ΔZ₂ = ΔB₂?

Also investigate: is there a natural chain homotopy between C(T) and C(T')
that demonstrates β₂ invariance?

Author: opus-2026-03-08-S45
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


def classify_paths(A, n, p, u, v):
    """Split A_p(T) paths into shared/arc categories."""
    ap = enumerate_allowed_paths(A, n, p)
    shared = []
    uses_uv = []
    for path in ap:
        found = False
        for i in range(len(path)-1):
            if path[i] == u and path[i+1] == v:
                uses_uv.append(path)
                found = True
                break
        if not found:
            shared.append(path)
    return ap, shared, uses_uv


def compute_omega_and_cycles(A, n, max_p=4):
    """Compute full chain complex data."""
    ap = {}
    om = {}
    bd = {}

    for p in range(max_p + 1):
        ap[p] = enumerate_allowed_paths(A, n, p)
        if p == 0:
            om[p] = np.eye(n)
        elif ap[p]:
            om[p] = compute_omega_basis(A, n, p, ap[p], ap[p-1])
        else:
            om[p] = np.zeros((0, 0))

    for p in range(1, max_p + 1):
        if ap[p] and ap[p-1]:
            bd[p] = build_full_boundary_matrix(ap[p], ap[p-1])
        else:
            bd[p] = None

    return ap, om, bd


n = 5
pairs = [(i, j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)

print("=" * 70)
print("MAYER-VIETORIS ANALYSIS FOR ARC FLIPS")
print("=" * 70)

# Pick a specific case where ΔZ₂ ≠ 0 to study
# Find one from (a,b,c,d)=(0,1,0,2) with ΔZ₂=1
test_case = None
for bits in range(1 << m):
    A = [[0]*n for _ in range(n)]
    for idx, (i,j) in enumerate(pairs):
        if (bits >> idx) & 1: A[i][j] = 1
        else: A[j][i] = 1

    for u in range(n):
        for v in range(n):
            if u == v or not A[u][v]:
                continue
            # Compute ABCD
            sets = {'A': 0, 'B': 0, 'C': 0, 'D': 0}
            for w in range(n):
                if w == u or w == v:
                    continue
                uw = A[u][w]
                wv = A[w][v]
                if uw and wv:
                    sets['A'] += 1
                elif not uw and not wv:
                    sets['B'] += 1
                elif uw and not wv:
                    sets['C'] += 1
                else:
                    sets['D'] += 1
            if sets == {'A': 0, 'B': 1, 'C': 0, 'D': 2}:
                test_case = (bits, u, v)
                break
        if test_case:
            break
    if test_case:
        break

bits, u, v = test_case
print(f"\nTest case: T#{bits}, flip ({u}→{v}) to ({v}→{u})")
A = [[0]*n for _ in range(n)]
for idx, (i,j) in enumerate(pairs):
    if (bits >> idx) & 1: A[i][j] = 1
    else: A[j][i] = 1

for i in range(n):
    out = [j for j in range(n) if A[i][j]]
    print(f"  {i} → {out}")

# Build flipped tournament
A_flip = [row[:] for row in A]
A_flip[u][v] = 0
A_flip[v][u] = 1

# Classify paths before and after
print(f"\nPath classification (before flip):")
for p in range(5):
    all_p, shared, uses = classify_paths(A, n, p, u, v)
    print(f"  p={p}: total={len(all_p)}, shared={len(shared)}, uses_uv={len(uses)}")

print(f"\nPath classification (after flip):")
for p in range(5):
    all_p, shared, uses = classify_paths(A_flip, n, p, v, u)
    print(f"  p={p}: total={len(all_p)}, shared={len(shared)}, uses_vu={len(uses)}")

# Note: "shared" paths should be the same for both T and T'
# (paths not using either u→v or v→u)
print(f"\nShared paths (using neither u→v nor v→u):")
for p in range(5):
    _, shared_T, _ = classify_paths(A, n, p, u, v)
    _, shared_T2, _ = classify_paths(A_flip, n, p, v, u)
    # shared_T = paths not using u→v in T
    # shared_T2 = paths not using v→u in T'
    # These should be the same set!
    set_T = set(tuple(x) for x in shared_T)
    set_T2 = set(tuple(x) for x in shared_T2)
    if set_T == set_T2:
        print(f"  p={p}: |shared| = {len(set_T)} ✓")
    else:
        print(f"  p={p}: DIFFER! |T\\arc|={len(set_T)}, |T'\\arc|={len(set_T2)}")
        print(f"    In T only: {set_T - set_T2}")
        print(f"    In T' only: {set_T2 - set_T}")

# Full chain complex data
ap_T, om_T, bd_T = compute_omega_and_cycles(A, n)
ap_Tf, om_Tf, bd_Tf = compute_omega_and_cycles(A_flip, n)

def dim_om(om):
    return om.shape[1] if om.ndim == 2 and om.shape[0] > 0 else 0

print(f"\n{'='*70}")
print("CHAIN COMPLEX COMPARISON")
print("=" * 70)

for p in range(5):
    d_T = dim_om(om_T[p])
    d_Tf = dim_om(om_Tf[p])
    print(f"  Ω_{p}: T={d_T}, T'={d_Tf}, Δ={d_Tf - d_T}")

# Compute Z₂, B₂ for both
def get_z2_b2(ap, om, bd, n):
    d2 = dim_om(om[2])
    if d2 == 0:
        return 0, 0

    bd2_om = bd[2] @ om[2]
    coords = np.linalg.lstsq(om[1], bd2_om, rcond=None)[0]
    rk = np.linalg.matrix_rank(coords, tol=1e-8)
    z2 = d2 - rk

    d3 = dim_om(om[3])
    if d3 > 0:
        bd3_om = bd[3] @ om[3]
        bd3_coords = np.linalg.lstsq(om[2], bd3_om, rcond=None)[0]
        b2 = np.linalg.matrix_rank(bd3_coords, tol=1e-8)
    else:
        b2 = 0

    return z2, b2

z2_T, b2_T = get_z2_b2(ap_T, om_T, bd_T, n)
z2_Tf, b2_Tf = get_z2_b2(ap_Tf, om_Tf, bd_Tf, n)

print(f"\n  Z₂: T={z2_T}, T'={z2_Tf}, Δ={z2_Tf - z2_T}")
print(f"  B₂: T={b2_T}, T'={b2_Tf}, Δ={b2_Tf - b2_T}")
print(f"  β₂: T={z2_T-b2_T}, T'={z2_Tf-b2_Tf}")

# CRUCIAL: look at HOW the Ω bases change
# Which Ω₂ elements are "new" and which are "lost"?
print(f"\n{'='*70}")
print("Ω₂ BASIS COMPARISON (in A₂ coordinates)")
print("=" * 70)

# Express om_T[2] and om_Tf[2] in a COMMON A₂ space
# Problem: A₂(T) ≠ A₂(T') since some paths are in one but not the other
# Solution: embed both in the "union" A₂ space

ap2_T = [tuple(x) for x in ap_T[2]]
ap2_Tf = [tuple(x) for x in ap_Tf[2]]
ap2_union = sorted(set(ap2_T) | set(ap2_Tf))

print(f"  |A₂(T)| = {len(ap2_T)}, |A₂(T')| = {len(ap2_Tf)}, |A₂(union)| = {len(ap2_union)}")
print(f"  A₂(T) only: {set(ap2_T) - set(ap2_Tf)}")
print(f"  A₂(T') only: {set(ap2_Tf) - set(ap2_T)}")

# What's the structure of the "new" and "lost" allowed 2-paths?
new_2paths = set(ap2_Tf) - set(ap2_T)
lost_2paths = set(ap2_T) - set(ap2_Tf)
print(f"\n  New 2-paths (in T' but not T): {sorted(new_2paths)}")
print(f"  Lost 2-paths (in T but not T'): {sorted(lost_2paths)}")

# Check: do the new/lost paths involve v→u / u→v?
for p in sorted(new_2paths):
    uses_vu = any(p[i] == v and p[i+1] == u for i in range(len(p)-1))
    print(f"    {p}: uses v→u: {uses_vu}")
for p in sorted(lost_2paths):
    uses_uv = any(p[i] == u and p[i+1] == v for i in range(len(p)-1))
    print(f"    {p}: uses u→v: {uses_uv}")

# Similarly for 3-paths
ap3_T = set(tuple(x) for x in ap_T[3])
ap3_Tf = set(tuple(x) for x in ap_Tf[3])
new_3 = ap3_Tf - ap3_T
lost_3 = ap3_T - ap3_Tf

print(f"\n  |A₃(T)| = {len(ap3_T)}, |A₃(T')| = {len(ap3_Tf)}")
print(f"  New 3-paths: {len(new_3)}, Lost 3-paths: {len(lost_3)}")
print(f"  Δ|A₃| = {len(new_3) - len(lost_3)}")

# STRUCTURAL INSIGHT: The boundary map ∂ connects p-paths to (p-1)-paths.
# When a p-path is lost, its boundary faces may or may not be lost.
# The "boundary compatibility" between lost/gained paths may force ΔZ₂ = ΔB₂.

print(f"\n{'='*70}")
print("BOUNDARY STRUCTURE OF LOST/GAINED PATHS")
print("=" * 70)

# For each lost 2-path, compute its boundary (in A₁)
ap1_T = [tuple(x) for x in ap_T[1]]
print(f"\nLost 2-paths and their boundaries:")
for p in sorted(lost_2paths):
    # boundary of (a,b,c) = (b,c) - (a,c) + (a,b)
    a, b, c = p
    faces = [(b,c), (a,c), (a,b)]
    signs = [+1, -1, +1]
    print(f"  ∂({a},{b},{c}) = +({b},{c}) - ({a},{c}) + ({a},{b})")
    for face, sign in zip(faces, signs):
        in_T = face in set(tuple(x) for x in ap_T[1])
        in_Tf = face in set(tuple(x) for x in ap_Tf[1])
        status = "shared" if in_T and in_Tf else "lost" if in_T and not in_Tf else "new" if not in_T and in_Tf else "neither"
        print(f"    {'+' if sign > 0 else '-'}({face[0]},{face[1]}): {status}")

print(f"\nGained 2-paths and their boundaries:")
for p in sorted(new_2paths):
    a, b, c = p
    faces = [(b,c), (a,c), (a,b)]
    signs = [+1, -1, +1]
    print(f"  ∂({a},{b},{c}) = +({b},{c}) - ({a},{c}) + ({a},{b})")
    for face, sign in zip(faces, signs):
        in_T = face in set(tuple(x) for x in ap_T[1])
        in_Tf = face in set(tuple(x) for x in ap_Tf[1])
        status = "shared" if in_T and in_Tf else "lost" if in_T and not in_Tf else "new" if not in_T and in_Tf else "neither"
        print(f"    {'+' if sign > 0 else '-'}({face[0]},{face[1]}): {status}")

# Now do the same for 3-paths
print(f"\nLost 3-paths and their faces:")
for p in sorted(lost_3):
    a, b, c, d = p
    faces = [(b,c,d), (a,c,d), (a,b,d), (a,b,c)]
    signs = [+1, -1, +1, -1]
    statuses = []
    for face in faces:
        in_T = face in set(ap2_T)
        in_Tf = face in set(ap2_Tf)
        if in_T and in_Tf:
            statuses.append("shared")
        elif in_T:
            statuses.append("lost")
        elif in_Tf:
            statuses.append("new")
        else:
            statuses.append("absent")
    print(f"  ({a},{b},{c},{d}): faces [{', '.join(statuses)}]")

print(f"\nGained 3-paths and their faces:")
for p in sorted(new_3):
    a, b, c, d = p
    faces = [(b,c,d), (a,c,d), (a,b,d), (a,b,c)]
    statuses = []
    for face in faces:
        in_T = face in set(ap2_T)
        in_Tf = face in set(ap2_Tf)
        if in_T and in_Tf:
            statuses.append("shared")
        elif in_T:
            statuses.append("lost")
        elif in_Tf:
            statuses.append("new")
        else:
            statuses.append("absent")
    print(f"  ({a},{b},{c},{d}): faces [{', '.join(statuses)}]")

# KEY OBSERVATION: is there a natural BIJECTION between lost and gained
# paths that respects boundary structure?
print(f"\n{'='*70}")
print("BIJECTION SEARCH: lost ↔ gained")
print("=" * 70)

# The swap σ: u↔v maps paths using u→v to paths using v→u
# σ(a₁,...,aₖ) = (σ(a₁),...,σ(aₖ)) where σ(u)=v, σ(v)=u, σ(w)=w
def swap_uv(path, u, v):
    return tuple(v if x == u else u if x == v else x for x in path)

print(f"\nSwap σ: {u}↔{v}")
print(f"\nLost 2-paths → σ(lost):")
for p in sorted(lost_2paths):
    sp = swap_uv(p, u, v)
    in_gained = sp in new_2paths
    print(f"  {p} → {sp}: in gained? {in_gained}")

print(f"\nGained 2-paths → σ(gained):")
for p in sorted(new_2paths):
    sp = swap_uv(p, u, v)
    in_lost = sp in lost_2paths
    print(f"  {p} → {sp}: in lost? {in_lost}")

# If swap is a bijection lost ↔ gained, it preserves boundary structure
# (since σ commutes with ∂), and this might be the key to β₂ invariance.

print("\nDone.")
