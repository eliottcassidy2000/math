"""
scsf_n9_gs.py — Compute SCSF(9) using GS tiling structure.

KEY INSIGHT: By the TRANS_MAP theorem, M is GS ↔ T(M) is SC.
So SC iso classes correspond exactly to GS tilings.
At n=9: free_bits = 4 (self-paired) + 12 (pair-free) = 16 → 2^16 = 65536 GS tilings.

For SCSF at odd n (blackself only):
  Class C is SCSF iff C is SC AND ∃ M∈C with flip(M) ∈ C.
  For GS tiling M: flip(M) is non-GS (proved at odd n).
  So check: ∃ GS tiling M such that canon(flip(M)) = canon(M).

COMPLEXITY: Canonicalize 131072 tilings (65536 GS + 65536 flips) with 362880 perms.
Estimated time: ~15 minutes.

oracle-2026-05-11
"""

import numpy as np
import sys, time
from itertools import permutations
from collections import defaultdict

sys.stdout.reconfigure(encoding='utf-8')
t_global = time.time()

n = 9
tiles = [(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]
m = len(tiles)
assert m == 28

verts = list(range(n, 0, -1))  # [9,8,...,1], position i = verts[i]
xi_arr = np.array([verts.index(x) for x,y in tiles], dtype=np.int32)
yi_arr = np.array([verts.index(y) for x,y in tiles], dtype=np.int32)

# TRANS_MAP: tile k → tile tm[k]
tile_idx = {t:i for i,t in enumerate(tiles)}
tm = [tile_idx[(n-y+1, n-x+1)] for x,y in tiles]

# Find self-paired tiles and pairs
self_paired = [i for i in range(m) if tm[i] == i]
gs_pairs = [(i, tm[i]) for i in range(m) if tm[i] > i]  # each pair once
free_bits = len(self_paired) + len(gs_pairs)
assert free_bits == 16, f"Expected 16, got {free_bits}"

print(f"n={n}, m={m}", flush=True)
print(f"  Self-paired tiles: {len(self_paired)} = {[tiles[i] for i in self_paired]}", flush=True)
print(f"  GS pair-constraints: {len(gs_pairs)}", flush=True)
print(f"  free_bits = {free_bits}, GS count = 2^{free_bits} = {1<<free_bits}", flush=True)

# ─── ENUMERATE ALL GS TILINGS ─────────────────────────────────────────────────
# GS mask: for each pair (i,j), bit_i = bit_j.
# Parametrize by 2^16 free bits: first 4 = self-paired bits, next 12 = pair bits.

def build_gs_masks():
    """Return (2^16,) array of all GS tile masks."""
    NGS = 1 << free_bits
    # Map: free_bit_pos → which tile bit(s) it controls
    free_to_tiles = []
    for sp in self_paired:
        free_to_tiles.append([sp])  # single tile
    for i, j in gs_pairs:
        free_to_tiles.append([i, j])  # paired tiles must match

    # free_params: (NGS, free_bits) binary
    fp_range = np.arange(free_bits, dtype=np.int32)
    free_idx = np.arange(NGS, dtype=np.uint32)
    free_params = ((free_idx[:, None] >> fp_range[None, :]) & 1).astype(np.uint8)

    # Build masks
    masks = np.zeros(NGS, dtype=np.uint32)
    for fbi, tile_list in enumerate(free_to_tiles):
        for ti in tile_list:
            masks |= (free_params[:, fbi].astype(np.uint32) << ti)
    return masks

print(f"\nBuilding GS masks...", flush=True)
gs_masks = build_gs_masks()
NGS = len(gs_masks)
print(f"  {NGS} GS tilings built", flush=True)

# Flip: flip all tile bits
flip_all = np.uint32((1 << m) - 1)
flip_masks = gs_masks ^ flip_all

# ─── BUILD ADJACENCY MATRICES ──────────────────────────────────────────────────
print(f"\nBuilding adjacency matrices for {2*NGS} tilings...", flush=True)

def build_A_flat(masks_arr):
    """Build adjacency matrix array for given mask array."""
    NC = len(masks_arr)
    k_range = np.arange(m, dtype=np.uint32)
    bits = ((masks_arr[:, None] >> k_range[None, :]) & 1).astype(np.uint8)
    A = np.zeros((NC, n*n), dtype=np.uint8)
    # Base path: position k beats position k+1
    for k in range(n-1):
        A[:, k*n + (k+1)] = 1
    # Tile arcs
    for k in range(m):
        xk, yk = int(xi_arr[k]), int(yi_arr[k])
        b = bits[:, k]
        A[:, xk*n + yk] += (1 - b)
        A[:, yk*n + xk] += b
    return A

A_gs = build_A_flat(gs_masks)
A_flip = build_A_flat(flip_masks)
print(f"  Done. A_gs shape: {A_gs.shape}", flush=True)

# Sanity check: verify first GS tiling is a valid tournament
A0 = A_gs[0].reshape(n, n)
for i in range(n):
    for j in range(i+1, n):
        assert A0[i,j] + A0[j,i] == 1, f"Invalid tournament at ({i},{j})"
print("  Tournament sanity check passed.", flush=True)

# Verify first GS tiling is actually GS (A[i,j] = A[n-1-j, n-1-i] for all non-diag)
# GS condition on adjacency: A[pos_x, pos_y] = A[pos_{n-y+1}-indexed, ...]
# TRANS_MAP: tile (x,y) → (n-y+1, n-x+1). Positions: pos(x) = verts.index(x) = n-x.
# So TRANS_MAP in position space: pos(x) → pos(n-y+1) where tile is (x,y).
# This is complex; skip the check.

# ─── PERMUTATION FLAT INDICES ─────────────────────────────────────────────────
print(f"\nPrecomputing permutation indices (9! = 362880)...", flush=True)
t0 = time.time()
perms_list = list(permutations(range(n)))
NP = len(perms_list)
perm_flat = np.zeros((NP, n*n), dtype=np.int32)
for p_idx, sigma in enumerate(perms_list):
    sigma_arr = np.array(sigma, dtype=np.int32)
    for i in range(n):
        perm_flat[p_idx, i*n:(i+1)*n] = sigma_arr[i]*n + sigma_arr
print(f"  Done ({NP} perms, {time.time()-t0:.1f}s)", flush=True)

# ─── COMBINED CANONICALIZATION ─────────────────────────────────────────────────
# Canonicalize GS tilings and their flips together (2*NGS tilings in one pass)
print(f"\nCanonicalization of {2*NGS} tilings...", flush=True)

# Combine: rows 0..NGS-1 = GS tilings, rows NGS..2*NGS-1 = flips
A_combined = np.concatenate([A_gs, A_flip], axis=0)  # (2*NGS, n*n)
NC = 2 * NGS

canon = np.full((NC, n*n), 255, dtype=np.uint8)
t0 = time.time()
REPORT_EVERY = 36288  # ~10% of 362880

for p_idx in range(NP):
    perm_A = A_combined[:, perm_flat[p_idx]]
    diff = perm_A != canon
    has_diff = diff.any(axis=1)
    if not has_diff.any():
        continue
    first = np.argmax(diff, axis=1)
    ridx = np.arange(NC)
    is_less = (perm_A[ridx, first] < canon[ridx, first]) & has_diff
    canon[is_less] = perm_A[is_less]

    if (p_idx + 1) % REPORT_EVERY == 0:
        pct = 100*(p_idx+1)/NP
        elapsed = time.time() - t0
        eta = elapsed * (NP - p_idx - 1) / (p_idx + 1)
        print(f"  {p_idx+1}/{NP} ({pct:.0f}%)  {elapsed:.1f}s  ETA {eta:.0f}s", flush=True)

print(f"  Canonicalization done ({time.time()-t0:.1f}s)", flush=True)

# ─── GROUP GS TILINGS INTO ISO CLASSES ────────────────────────────────────────
print(f"\nGrouping GS tilings into iso classes...", flush=True)
gs_groups = defaultdict(list)
for i in range(NGS):
    gs_groups[canon[i].tobytes()].append(i)

print(f"  GS iso classes (= SC classes): {len(gs_groups)}", flush=True)

# ─── CHECK SCSF: does any GS tiling have its flip in the same class? ──────────
print(f"\nChecking SCSF condition...", flush=True)

# canon[NGS + i] = canonical form of flip(gs_masks[i])
scsf_classes = []
checked_sigs = set()

for ci, (sig, gs_members) in enumerate(gs_groups.items()):
    # Check if any GS member's flip has the same canonical form
    for gi in gs_members:
        flip_sig = canon[NGS + gi].tobytes()
        if flip_sig == sig:
            if sig not in checked_sigs:
                checked_sigs.add(sig)
                # Count GS members that are self-flip
                self_flip_gs = sum(1 for g in gs_members if canon[NGS+g].tobytes() == sig)
                # Count all flips landing in this class
                all_flips_here = sum(1 for g in gs_members if canon[NGS+g].tobytes() == sig)
                scsf_classes.append({
                    'sig': sig, 'gs_members': gs_members,
                    'gs_count': len(gs_members),
                    'self_flip_count': all_flips_here,
                })
            break

# ─── FULL SC COUNT: verify against THM-283 ────────────────────────────────────
total_sc = len(gs_groups)

# ─── COMPUTE |Aut| AND SCORES FOR KERNEL CLASSES ──────────────────────────────
if scsf_classes:
    print(f"\nComputing |Aut| and scores for {len(scsf_classes)} kernel classes...", flush=True)

    # Compute |Aut| for each kernel class representative
    reps = [gs_members[0] for c in scsf_classes for gs_members in [c['gs_members']]]
    A_reps = np.array([A_gs[r] for r in reps], dtype=np.uint8)
    aut_counts = np.zeros(len(scsf_classes), dtype=np.int32)
    t0 = time.time()
    for p_idx in range(NP):
        perm_A = A_reps[:, perm_flat[p_idx]]
        is_aut = (perm_A == A_reps).all(axis=1)
        aut_counts += is_aut.astype(np.int32)
    print(f"  |Aut| done ({time.time()-t0:.1f}s)", flush=True)

    for i, c in enumerate(scsf_classes):
        c['au'] = int(aut_counts[i])
        # Score sequence of representative
        rep = c['gs_members'][0]
        c['scores'] = sorted(A_gs[rep].reshape(n,n).sum(axis=1).tolist(), reverse=True)
        # #tilings = n! / |Aut|
        c['nt'] = NP // c['au']

# ─── PRINT RESULTS ────────────────────────────────────────────────────────────
print(f"\n{'='*60}")
print(f"RESULTS: n=9")
print(f"{'='*60}")
print(f"  GS tilings:      {NGS} = 2^16")
print(f"  SC iso classes:  {total_sc}  (THM-283 says 2752)")
print(f"  SC∩SF (SCSF):    {len(scsf_classes)}")

print(f"\nKernel classes (SC∩SF):")
print(f"{'#':>3} | {'#gs':>4} | {'#sflip':>6} | {'|Aut|':>6} | {'#t':>6} | scores")
print("-"*70)
for i, c in enumerate(scsf_classes):
    au = c.get('au', '?')
    nt = c.get('nt', '?')
    scores = c.get('scores', '?')
    print(f"{i+1:>3} | {c['gs_count']:>4} | {c['self_flip_count']:>6} | {au:>6} | {nt:>6} | {scores}")

print(f"\nAll SC class sizes (top 20 by #gs tilings):")
gs_counts = sorted([(len(v), k) for k, v in gs_groups.items()], reverse=True)
for gs_c, sig in gs_counts[:20]:
    au_approx = NP // gs_c if gs_c > 0 else "?"
    scsf_mark = " ← SCSF" if sig in checked_sigs else ""
    print(f"  gs={gs_c:>4}, |Aut|≈{au_approx:>6}{scsf_mark}")

print(f"\nTotal time: {time.time()-t_global:.1f}s")
