#!/usr/bin/env python3
"""
beta2_omega2_relations.py — Understanding Ω₂ relations between TT and NT

KEY QUESTION: Z₂ has NT components, yet DT boundaries (TT-only) span Z₂.
How? Because within Ω₂, NT triples are linearly related to TT triples!

Specifically: Ω₂ is defined as ker(∂₁ ∘ ∂₂) = {z ∈ A₂ : ∂₁(∂₂(z)) = 0}.
Wait, that's not right. Ω₂ = ker(∂̃₂) where ∂̃₂ is the natural map.

Actually, Ω_p is defined inductively:
- Ω_0 = R^V (vertices)  
- Ω_p = {c ∈ A_p : ∂_p(c) ∈ Ω_{p-1}}

So Ω₂ = {z ∈ A₂ : ∂₂(z) ∈ Ω₁}.
And Ω₁ = {z ∈ A₁ : ∂₁(z) ∈ Ω₀} = A₁ (since Ω₀ = R^V).

For tournaments: Ω₁ = A₁ = all directed edges = R^{n(n-1)/2}.
So Ω₂ = {z ∈ A₂ : ∂₂(z) ∈ A₁}.

∂₂(a,b,c) = (b,c) - (a,c) + (a,b).
For a TT triple (a→c): all three faces (b,c), (a,c), (a,b) are in A₁. ✓
For an NT triple (c→a): face (a,c) is NOT in A₁! It's (c,a) that's in A₁.

So ∂₂(a,b,c) for NT triple has (a,c) which is NOT an edge. Thus:
- Individual NT triples are NOT in Ω₂
- But combinations like (a,b,c) - (a,b',c) where both are NT can be in Ω₂
  IF the "bad" (a,c) terms cancel

Wait: in what basis are we working? Let me reconsider.

A₂ = set of allowed 2-paths. An allowed 2-path (a,b,c) means a→b AND b→c.
The FACES of (a,b,c) are: 
  ∂₂(a,b,c) = (b,c) - (a,c) + (a,b)

Here (a,c) means the formal 1-chain e_{(a,c)} in span(A₁).
If a→c, then (a,c) ∈ A₁ and e_{(a,c)} has coefficient -1.
If c→a, then (a,c) ∉ A₁, but (c,a) ∈ A₁. The map ∂₂ gives -(a,c).

BUT in GLMY theory, (a,c) is NOT in A₁ if c→a. So the image ∂₂(a,b,c)
contains a term -(a,c) which is NOT in span(A₁).

Therefore: (a,b,c) ∈ Ω₂ iff ∂₂(a,b,c) ∈ span(A₁), which means all 
faces must be in A₁. For (a,b,c): faces are (b,c), (a,c), (a,b).
  (a,b) ∈ A₁ iff a→b ✓ (always, from allowed path)
  (b,c) ∈ A₁ iff b→c ✓ (always, from allowed path)
  (a,c) ∈ A₁ iff a→c — this is the TT condition!

SO: individual 2-paths in Ω₂ are EXACTLY the TT triples!

But wait — Ω₂ also contains LINEAR COMBINATIONS:
z = Σ aᵢ (xᵢ, yᵢ, zᵢ) ∈ Ω₂ iff ∂₂(z) ∈ span(A₁).

If z has NT triples, their "bad face" terms must cancel.

Author: opus-2026-03-08-S43
"""
import sys
import numpy as np
from collections import Counter
sys.path.insert(0, '04-computation')
sys.stdout.reconfigure(line_buffering=True)
from path_homology_v2 import (
    enumerate_allowed_paths, compute_omega_basis, build_full_boundary_matrix,
)

def all_tournaments(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(pairs)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i,j) in enumerate(pairs):
            if (mask >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
        yield A

print("=" * 70)
print("Ω₂ RELATION STRUCTURE: TT bases + NT cancellations")
print("=" * 70)

n = 5
print(f"\n--- n = {n} ---")

# For each tournament, decompose Ω₂ into:
# 1. TT triples (individually in Ω₂)
# 2. NT cancellation elements (combinations of NT triples in Ω₂)

nt_cancel_dims = Counter()
total_nt_cancel = 0
total_tournaments = 0

for A in all_tournaments(n):
    total_tournaments += 1
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    
    if not a2:
        continue
    
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_Om2 = om2.shape[1] if om2.ndim == 2 else 0
    
    tt_count = sum(1 for p in a2 if A[p[0]][p[2]] == 1)
    nt_count = len(a2) - tt_count
    
    # dim(Ω₂) = tt_count + dim(NT cancellation space)
    nt_cancel_dim = dim_Om2 - tt_count
    nt_cancel_dims[nt_cancel_dim] += 1
    if nt_cancel_dim > 0:
        total_nt_cancel += 1

print(f"  NT cancellation dimension distribution: {dict(sorted(nt_cancel_dims.items()))}")
print(f"  Tournaments with NT cancellation elements: {total_nt_cancel}/{total_tournaments}")
print(f"  dim(Ω₂) = |TT| + dim(NT_cancel)")

# Now the CRITICAL question: does Z₂ have components in NT_cancel direction?
# If so, DT boundaries can't reach those components directly.
# But the fact that DT fills Z₂ means: Z₂ ∩ NT_cancel ⊆ im(∂₃|DT) projected onto NT_cancel.
# 
# Actually this is confusing. Let me think more carefully.
#
# Ω₂ has a basis: {e_{TT_1}, ..., e_{TT_k}, c_1, ..., c_m}
# where e_{TT_i} are unit vectors for TT triples, and c_j are NT cancellation elements.
#
# A Z₂ element z = Σ αᵢ e_{TT_i} + Σ βⱼ cⱼ with ∂₂(z) = 0.
#
# DT boundary ∂₃(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c).
# For DT: all 4 faces are TT. So ∂₃(DT) ∈ span(TT).
# In Ω₂ coordinates: ∂₃(DT) has nonzero only in TT positions.
#
# If z has βⱼ ≠ 0 but DT fills z, then:
# z = DT_boundary, which is in span(TT).
# This means: βⱼ = 0 for the Z₂ elements that ARE boundaries.
#
# But we found Z₂ elements WITH NT components. These must be filled by 
# NON-DT Ω₃ elements (which CAN have NT components in their boundary).
#
# WAIT: the section script showed DT alone fills Z₂ at n=5. But the NT
# analysis shows Z₂ has NT components. This seems contradictory!
#
# Resolution: the A₂ coordinates include both TT and NT paths.
# An Ω₂ element is a LINEAR COMBINATION of A₂ paths.
# When I say "DT boundary spans Z₂", I mean in Ω₂ coordinates.
# In A₂ coordinates, the DT boundary may have nonzero NT components 
# (through the Ω₂ basis change)!
#
# Let me verify: does ∂₃(DT) in A₂ coords have NT components?

print(f"\n{'='*70}")
print("DT BOUNDARY IN A₂ COORDINATES: NT COMPONENTS?")
print("=" * 70)

n = 5
dt_has_nt_in_a2 = 0
dt_no_nt = 0

for idx, A in enumerate(all_tournaments(n)):
    if idx >= 100:
        break
    
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    
    if not a3:
        continue
    
    bd3 = build_full_boundary_matrix(a3, a2)
    
    tt_mask = np.array([A[p[0]][p[2]] == 1 for p in a2])
    
    # For each DT path
    for j, p in enumerate(a3):
        a, b, c, d = p
        if A[a][c] == 1 and A[b][d] == 1:  # DT
            boundary = bd3[:, j]  # in A₂ coords
            nt_components = boundary[~tt_mask]
            if nt_components.size > 0 and np.max(np.abs(nt_components)) > 1e-8:
                dt_has_nt_in_a2 += 1
            else:
                dt_no_nt += 1

print(f"  DT boundaries with NT components in A₂: {dt_has_nt_in_a2}")
print(f"  DT boundaries without NT components: {dt_no_nt}")

if dt_has_nt_in_a2 == 0:
    print(f"  *** DT boundaries are PURE TT in A₂! ***")
    print(f"""
  This means: ∂₃(DT path) only involves TT 2-paths.
  
  For DT path (a,b,c,d): faces are:
    (b,c,d) — TT iff b→d. YES (DT condition)
    (a,c,d) — TT iff a→d. Is a→d guaranteed? 
    (a,b,d) — TT iff a→d. Same question!
    (a,b,c) — TT iff a→c. YES (DT condition)
  
  Wait: (a,c,d) being TT needs a→d. Is that guaranteed?
  A DT path has: a→b→c→d, a→c, b→d.
  But a→d is NOT guaranteed! The edge a↔d could go either way.
  
  So DT faces CAN be NT! Let me recheck the data...
""")
else:
    print(f"  DT boundaries DO have NT components in A₂!")
    print(f"  This resolves the paradox: DT boundaries can fill NT-containing Z₂ elements")
    print(f"  because DT 4-path faces include NT 2-paths when the \"outer\" edge goes backward.")

# Let me verify: which DT faces are NT?
print(f"\n--- DT face analysis ---")
n = 5
dt_face_patterns = Counter()

for A in all_tournaments(n):
    a3 = enumerate_allowed_paths(A, n, 3)
    if not a3:
        continue
    
    for p in a3:
        a, b, c, d = p
        if not (A[a][c] == 1 and A[b][d] == 1):
            continue  # not DT
        
        # Faces: (b,c,d), (a,c,d), (a,b,d), (a,b,c)
        f0_tt = A[b][d]  # b→d for (b,c,d) to be TT — YES (DT condition)
        f1_tt = A[a][d]  # a→d for (a,c,d) to be TT
        f2_tt = A[a][d]  # a→d for (a,b,d) to be TT — same edge!
        f3_tt = A[a][c]  # a→c for (a,b,c) to be TT — YES (DT condition)
        
        pattern = (f0_tt, f1_tt, f2_tt, f3_tt)
        dt_face_patterns[pattern] += 1

print(f"  DT face TT patterns (f0,f1,f2,f3):")
for pat, count in sorted(dt_face_patterns.items()):
    print(f"    {pat}: {count}")
    # f0 and f3 should always be 1 (TT)
    # f1 and f2 depend on a→d edge

# So DT faces (a,c,d) and (a,b,d) are NT when d→a!
# This means DT boundaries DO include NT 2-paths.
# This is why DT boundaries can span Z₂ elements with NT components.

print(f"""
RESOLUTION:
  A DT path (a,b,c,d) has a→c and b→d (DT conditions).
  But the edge a↔d is free! When d→a:
    - Face (a,c,d) has c→d but d→a, so it's NT (not a→d)
    - Face (a,b,d) has a→b but d→a, so it's NT (not a→d)
  
  So ∂₃(a,b,c,d) = (b,c,d) - (a,c,d) + (a,b,d) - (a,b,c)
  can have 2 NT faces when d→a, contributing NT components.
  
  This is the KEY MECHANISM: DT boundaries reach into the NT part of A₂
  through the FREE edge a↔d, allowing them to fill Z₂ elements 
  that have NT components.
  
  But wait — are the NT faces even in A₂? (a,c,d): a→c (DT), c→d (path).
  YES, it's an allowed 2-path regardless of a↔d direction.
  
  And (a,b,d): a→b (path), b→d (DT). YES, allowed.
  Both faces are allowed 2-paths, but they're NT when d→a.
""")

# NEXT: How many NT 2-paths can DT boundaries reach?
print(f"\n--- DT boundary NT-reachability ---")
n = 5

nt_reachable_fractions = []

for idx, A in enumerate(all_tournaments(n)):
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    
    if not a3:
        continue
    
    bd3 = build_full_boundary_matrix(a3, a2)
    tt_mask = np.array([A[p[0]][p[2]] == 1 for p in a2])
    nt_indices = np.where(~tt_mask)[0]
    
    if len(nt_indices) == 0:
        continue
    
    # DT paths
    dt_indices = [j for j, p in enumerate(a3) if A[p[0]][p[2]] == 1 and A[p[1]][p[3]] == 1]
    if not dt_indices:
        continue
    
    # DT boundaries restricted to NT 2-paths
    bd_dt_nt = bd3[np.ix_(nt_indices, dt_indices)]
    
    # Rank = how many NT directions are reachable
    rk = np.linalg.matrix_rank(bd_dt_nt, tol=1e-8)
    nt_reachable_fractions.append(rk / len(nt_indices))

print(f"  NT reachability fraction: mean={np.mean(nt_reachable_fractions):.3f}, min={min(nt_reachable_fractions):.3f}, max={max(nt_reachable_fractions):.3f}")

print("\nDone.")
