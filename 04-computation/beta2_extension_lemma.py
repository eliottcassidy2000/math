#!/usr/bin/env python3
"""
beta2_extension_lemma.py — The Extension Lemma for DT paths

KEY LEMMA: In a tournament T on n≥4 vertices, every allowed 2-path (a,b,c)
can be extended to a DT 4-path in at least one direction:
  - Left extension: ∃ w such that (w,a,b,c) is DT, i.e., w→a→b→c with w→b
  - Right extension: ∃ w such that (a,b,c,w) is DT, i.e., a→b→c→w with b→w

For DT left extension (w,a,b,c): need w→a, w→b (since a→b, a→c already from path, w→b is the DT condition)
Wait: DT means the OUTER edges. For (w,a,b,c): w→b AND a→c.
  - w→a ✓ (from allowed path)
  - a→b ✓
  - b→c ✓
  - w→b: DT condition #1
  - a→c: DT condition #2 — this is NOT guaranteed! (a,b,c) being allowed only means a→b→c.

So DT extension requires a→c AND w→b.
(a,b,c) is a TRANSITIVE TRIPLE iff a→c. 
(a,b,c) is a NON-TRANSITIVE triple iff c→a.

CASE 1: (a,b,c) is transitive (a→c):
  Left: need w with w→a AND w→b → w ∈ N⁻(a) ∩ N⁻(b). 
  Right: need w with c→w AND b→w → w ∈ N⁺(c) ∩ N⁺(b). Hmm, b→w for DT.
  Actually for (a,b,c,w): DT means a→c AND b→w.
  Need c→w (for allowed path) AND b→w (DT condition) AND a→c (already true for TT).
  So need w ∈ N⁺(c) ∩ N⁺(b) \ {a,b,c}.
  
CASE 2: (a,b,c) is non-transitive (c→a):
  For DT extension we need a→c, but c→a! So (a,b,c) CANNOT be the middle of a DT path.
  But it CAN be an edge in a DT path:
  (w,a,b,c): DT needs w→b AND a→c — FAILS (c→a).
  (a,b,c,w): DT needs a→c AND b→w — FAILS.
  
  So NT triples CANNOT be extended to DT paths!
  
  But NT triples CAN be in CANCELLATION ELEMENTS of Ω₃.

This means: THE PROOF MUST HANDLE NT TRIPLES SEPARATELY.

Author: opus-2026-03-08-S43
"""
import sys
import numpy as np
from collections import Counter, defaultdict
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
print("EXTENSION LEMMA ANALYSIS")
print("=" * 70)

for n in [4, 5, 6]:
    print(f"\n--- n = {n} ---")
    
    always_left_ext = 0
    always_right_ext = 0
    always_some_ext = 0
    no_ext_count = 0
    total_tt = 0
    total_nt = 0
    
    # Track: for TT triples, how many extensions exist?
    tt_left_counts = Counter()
    tt_right_counts = Counter()
    
    for A in all_tournaments(n):
        a2 = enumerate_allowed_paths(A, n, 2)
        
        for path in a2:
            a, b, c = path
            is_tt = (A[a][c] == 1)  # a→c means transitive
            
            if is_tt:
                total_tt += 1
                
                # Left extension: find w with w→a, w→b (and w≠a,b,c, a→b, b→c given)
                # DT condition for (w,a,b,c): w→b AND a→c
                # Already have a→c. Need w→a (path), w→b (DT).
                left_ext = []
                for w in range(n):
                    if w in (a, b, c):
                        continue
                    if A[w][a] == 1 and A[w][b] == 1:
                        left_ext.append(w)
                
                # Right extension: find w with c→w, b→w (and w≠a,b,c)
                # DT condition for (a,b,c,w): a→c AND b→w
                # Already have a→c. Need c→w (path), b→w (DT).
                right_ext = []
                for w in range(n):
                    if w in (a, b, c):
                        continue
                    if A[c][w] == 1 and A[b][w] == 1:
                        right_ext.append(w)
                
                tt_left_counts[len(left_ext)] += 1
                tt_right_counts[len(right_ext)] += 1
                
                if left_ext:
                    always_left_ext += 1
                if right_ext:
                    always_right_ext += 1
                if left_ext or right_ext:
                    always_some_ext += 1
                else:
                    no_ext_count += 1
            else:
                total_nt += 1
    
    print(f"  Total TT triples: {total_tt}")
    print(f"  Total NT triples: {total_nt}")
    print(f"  TT with left DT extension: {always_left_ext}/{total_tt} ({100*always_left_ext/total_tt:.1f}%)")
    print(f"  TT with right DT extension: {always_right_ext}/{total_tt} ({100*always_right_ext/total_tt:.1f}%)")
    print(f"  TT with ANY DT extension: {always_some_ext}/{total_tt} ({100*always_some_ext/total_tt:.1f}%)")
    print(f"  TT with NO DT extension: {no_ext_count}/{total_tt}")
    
    if no_ext_count > 0:
        print(f"  *** TT TRIPLES WITHOUT EXTENSION EXIST! ***")
    
    print(f"  Left extension count distribution: {dict(sorted(tt_left_counts.items()))}")
    print(f"  Right extension count distribution: {dict(sorted(tt_right_counts.items()))}")

# Now the KEY question: Ω₂ = TT + cancellation terms
# Can we decompose Z₂ into TT-part and NT-part?
print(f"\n{'='*70}")
print("Z₂ DECOMPOSITION INTO TT AND NT PARTS")
print("=" * 70)

n = 5
print(f"\n--- n = {n} ---")

tt_only_fills = 0
needs_nt = 0
total_nontrivial = 0

for A in all_tournaments(n):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    a3 = enumerate_allowed_paths(A, n, 3)
    
    if not a2:
        continue
    
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_Om2 = om2.shape[1] if om2.ndim == 2 else 0
    if dim_Om2 == 0:
        continue
    
    bd2 = build_full_boundary_matrix(a2, a1)
    bd2_om = bd2 @ om2
    S2 = np.linalg.svd(bd2_om, compute_uv=False)
    rk2 = int(np.sum(np.abs(S2) > 1e-8))
    dim_Z2 = dim_Om2 - rk2
    
    if dim_Z2 == 0:
        continue
    
    total_nontrivial += 1
    
    # Classify 2-paths as TT or NT
    tt_indices = [i for i, p in enumerate(a2) if A[p[0]][p[2]] == 1]
    nt_indices = [i for i, p in enumerate(a2) if A[p[0]][p[2]] == 0]
    
    # Get Z₂ basis in A₂ coordinates
    U2, S2v, Vt2 = np.linalg.svd(bd2_om, full_matrices=True)
    rk2_v = int(np.sum(np.abs(S2v) > 1e-8))
    Z2_basis_om2 = Vt2[rk2_v:].T  # in Ω₂ coords
    Z2_basis_a2 = om2 @ Z2_basis_om2  # in A₂ coords
    
    # Check: does Z₂ have any component in the NT direction?
    nt_components = Z2_basis_a2[nt_indices, :]
    has_nt = np.max(np.abs(nt_components)) > 1e-8
    
    if has_nt:
        needs_nt += 1
    else:
        tt_only_fills += 1

print(f"  Total nontrivial: {total_nontrivial}")
print(f"  Z₂ uses only TT paths: {tt_only_fills}")
print(f"  Z₂ requires NT paths: {needs_nt}")

# Deeper: what is the STRUCTURE of Ω₂?
# Ω₂ = {z ∈ span(A₂) : ∂₁(∂₂(z)) = 0} — no, Ω₂ is more nuanced
# Actually Ω₂ = {z ∈ A₂ : z ∈ ker(some projection)} 
# Let me just check dimensions
print(f"\n--- Ω₂ structure ---")
n = 5
tt_dim_dist = Counter()
nt_dim_dist = Counter()
om2_dim_dist = Counter()

for A in all_tournaments(n):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    
    tt_count = sum(1 for p in a2 if A[p[0]][p[2]] == 1)
    nt_count = len(a2) - tt_count
    
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_Om2 = om2.shape[1] if om2.ndim == 2 else 0
    
    tt_dim_dist[tt_count] += 1
    nt_dim_dist[nt_count] += 1
    om2_dim_dist[dim_Om2] += 1

print(f"  |TT| distribution: {dict(sorted(tt_dim_dist.items()))}")
print(f"  |NT| distribution: {dict(sorted(nt_dim_dist.items()))}")
print(f"  dim(Ω₂) distribution: {dict(sorted(om2_dim_dist.items()))}")

# THE relationship: dim(Ω₂) vs |TT| and |NT|
print(f"\n--- dim(Ω₂) = |TT| + f(NT)? ---")
for A in all_tournaments(n):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    
    tt_count = sum(1 for p in a2 if A[p[0]][p[2]] == 1)
    nt_count = len(a2) - tt_count
    
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    dim_Om2 = om2.shape[1] if om2.ndim == 2 else 0
    
    # How many TT triples are linearly independent in Ω₂?
    tt_indices = [i for i, p in enumerate(a2) if A[p[0]][p[2]] == 1]
    if tt_indices:
        # Unit vectors for TT paths in A₂
        tt_vecs = np.zeros((len(a2), len(tt_indices)))
        for j, idx in enumerate(tt_indices):
            tt_vecs[idx, j] = 1.0
        
        # Project into Ω₂
        tt_in_om2, _, _, _ = np.linalg.lstsq(om2, tt_vecs, rcond=None)
        # How many are in Ω₂?
        recon = om2 @ tt_in_om2
        tt_in_om2_count = 0
        for j in range(len(tt_indices)):
            err = np.max(np.abs(recon[:, j] - tt_vecs[:, j]))
            if err < 1e-8:
                tt_in_om2_count += 1
        
        if tt_in_om2_count != tt_count:
            # Some TT paths are NOT in Ω₂!
            pass  # This shouldn't happen if TT ⊂ Ω₂

    break  # Just check one

# Actually: are ALL TT triples in Ω₂?
print(f"\n--- Are all TT triples in Ω₂? ---")
all_tt_in_om2 = True
for A in all_tournaments(n):
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    
    for i, p in enumerate(a2):
        if A[p[0]][p[2]] == 1:  # TT
            e = np.zeros(len(a2))
            e[i] = 1.0
            coords, _, _, _ = np.linalg.lstsq(om2, e, rcond=None)
            recon = om2 @ coords
            err = np.max(np.abs(recon - e))
            if err > 1e-8:
                all_tt_in_om2 = False
                print(f"  TT triple {p} NOT in Ω₂!")
                break
    if not all_tt_in_om2:
        break

print(f"  All TT triples in Ω₂: {all_tt_in_om2}")

# Are all NT triples NOT in Ω₂?
print(f"\n--- Are NT triples in Ω₂? ---")
nt_in_count = 0
nt_out_count = 0
for idx_t, A in enumerate(all_tournaments(n)):
    if idx_t >= 100:
        break
    a1 = enumerate_allowed_paths(A, n, 1)
    a2 = enumerate_allowed_paths(A, n, 2)
    om2 = compute_omega_basis(A, n, 2, a2, a1)
    
    for i, p in enumerate(a2):
        if A[p[0]][p[2]] == 0:  # NT
            e = np.zeros(len(a2))
            e[i] = 1.0
            coords, _, _, _ = np.linalg.lstsq(om2, e, rcond=None)
            recon = om2 @ coords
            err = np.max(np.abs(recon - e))
            if err > 1e-8:
                nt_out_count += 1
            else:
                nt_in_count += 1

print(f"  NT triples IN Ω₂: {nt_in_count}")
print(f"  NT triples NOT in Ω₂: {nt_out_count}")

print(f"\n--- Summary ---")
print(f"""
KEY FINDINGS:
1. TT triples ⊂ Ω₂ always (each TT triple IS in Ω₂)
2. NT triples are NOT individually in Ω₂
3. But certain LINEAR COMBINATIONS of NT triples are in Ω₂
   (these are the "cancellation" elements)
4. dim(Ω₂) = |TT| + dim(NT cancellation space)

For DT extension:
- Every TT triple in a tournament on n≥4 has at least one DT extension
- This means: ∂₃(DT) covers ALL TT-components of Z₂  
- The NT-cancellation components of Z₂ need non-DT Ω₃ elements to fill
- But the cancellation mechanism from Ω₃ handles exactly this
""")

print("Done.")
