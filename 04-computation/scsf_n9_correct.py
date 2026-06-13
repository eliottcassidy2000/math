"""
scsf_n9_correct.py — CORRECT algorithm for SCSF(9) at odd n.

CRITICAL FIX: At odd n, SCSF classes are BLACKSELF — the self-flip tilings are
NON-GS (neither member of the flip pair is GS). The previous algorithm only
checked if flip(GS-tiling) ≅ GS-tiling, which gives 0 at all odd n (proved theorem).

CORRECT ALGORITHM:
  For each SC iso class C with GS representative M_C:
    For each permutation σ ∈ S_n:
      Compute N_σ = tiling of σ(T(M_C))   [this is a tiling in C's iso class]
      Compute T̃(N_σ)                       [flip all tile arcs]
      Check if T̃(N_σ) ≅ T(M_C)            [is T̃ in C?]
      If yes → C is SCSF.

KEY INSIGHT: N_σ is in class C (since σ(T(M_C)) ≅ T(M_C)), and if T̃(N_σ) ≅ T(M_C),
then flip(N_σ) is also in C → C is SF. Combined with C being SC → C is SCSF.

N_σ may be GS or non-GS. At odd n, the GS case gives no SCSF (proved). The
non-GS case can give SCSF (blackself).

COMPLEXITY: 2752 SC classes × 362880 perms = 1B iterations.
With numpy vectorization over all perms per class: feasible in ~5 minutes.

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

verts = list(range(n, 0, -1))
xi_arr = np.array([verts.index(x) for x,y in tiles], dtype=np.int32)
yi_arr = np.array([verts.index(y) for x,y in tiles], dtype=np.int32)
tile_idx = {t:i for i,t in enumerate(tiles)}
tm_map = [tile_idx[(n-y+1, n-x+1)] for x,y in tiles]

self_paired = [i for i in range(m) if tm_map[i] == i]
gs_pairs_idx = [(i, tm_map[i]) for i in range(m) if tm_map[i] > i]
free_bits = len(self_paired) + len(gs_pairs_idx)

print(f"n={n}, m={m}, free_bits={free_bits}", flush=True)

# ─── BUILD GS TILINGS ─────────────────────────────────────────────────────────
NGS = 1 << free_bits
free_to_tiles = [[sp] for sp in self_paired] + [[i,j] for i,j in gs_pairs_idx]
fp_range = np.arange(free_bits, dtype=np.int32)
free_idx = np.arange(NGS, dtype=np.uint32)
free_params = ((free_idx[:, None] >> fp_range[None, :]) & 1).astype(np.uint8)
gs_masks = np.zeros(NGS, dtype=np.uint32)
for fbi, tile_list in enumerate(free_to_tiles):
    for ti in tile_list:
        gs_masks |= (free_params[:, fbi].astype(np.uint32) << ti)

def build_A_batch(masks_arr):
    NC = len(masks_arr)
    k_range = np.arange(m, dtype=np.uint32)
    bits = ((masks_arr[:, None] >> k_range[None, :]) & 1).astype(np.uint8)
    A = np.zeros((NC, n, n), dtype=np.uint8)
    for k in range(n-1):
        A[:, k, k+1] = 1
    for k in range(m):
        b = bits[:, k]
        A[:, int(xi_arr[k]), int(yi_arr[k])] += (1 - b)
        A[:, int(yi_arr[k]), int(xi_arr[k])] += b
    return A

print(f"Building GS adjacency matrices ({NGS})...", flush=True)
A_gs = build_A_batch(gs_masks)   # (NGS, n, n)

# ─── PERMUTATION FLAT INDICES ─────────────────────────────────────────────────
print(f"Precomputing permutation indices (9! = 362880)...", flush=True)
perms_list = list(permutations(range(n)))
NP = len(perms_list)
# perm_flat[p, i*n+j] = sigma[i]*n + sigma[j] for permutation sigma
perm_flat = np.zeros((NP, n*n), dtype=np.int32)
for p_idx, sigma in enumerate(perms_list):
    sigma_arr = np.array(sigma, dtype=np.int32)
    for i in range(n):
        perm_flat[p_idx, i*n:(i+1)*n] = sigma_arr[i]*n + sigma_arr
print(f"  Done ({NP} perms, {time.time()-t_global:.1f}s)", flush=True)

# ─── TILE PAIR POSITIONS FOR FLIP ─────────────────────────────────────────────
# Tile arcs = pairs (i,j) with |i-j|>=2 (0-indexed positions)
# Flip: swap directions of all tile arcs
tp_i = np.array([i for i in range(n) for j in range(n) if abs(i-j) >= 2 and i < j], dtype=np.int32)
tp_j = np.array([j for i in range(n) for j in range(n) if abs(i-j) >= 2 and i < j], dtype=np.int32)
print(f"  Tile pairs for flip: {len(tp_i)} (should be {m//1} = {n*(n-1)//2 - (n-1)})", flush=True)

# Precompute: which (i,j) flattened positions to swap
tp_ij_flat = (tp_i * n + tp_j).astype(np.int32)  # positions (i,j)
tp_ji_flat = (tp_j * n + tp_i).astype(np.int32)  # positions (j,i)

# ─── CANONICALIZE GS TILINGS → SC ISO CLASSES ─────────────────────────────────
# We need canonical forms for the GS tilings to group them into SC classes.
# Use the approach: canonicalize A_gs (shape (NGS, n²)) with all NP permutations.
# Then group by canonical form.

print(f"\nCanonicalization of {NGS} GS tilings...", flush=True)
A_gs_flat = A_gs.reshape(NGS, n*n)
canon_gs = np.full((NGS, n*n), 255, dtype=np.uint8)
t0 = time.time()
for p_idx in range(NP):
    perm_A = A_gs_flat[:, perm_flat[p_idx]]
    diff = perm_A != canon_gs
    has_diff = diff.any(axis=1)
    if not has_diff.any(): continue
    first = np.argmax(diff, axis=1)
    ridx = np.arange(NGS)
    is_less = (perm_A[ridx, first] < canon_gs[ridx, first]) & has_diff
    canon_gs[is_less] = perm_A[is_less]
    if (p_idx+1) % 36288 == 0:
        pct = 100*(p_idx+1)/NP
        elapsed = time.time()-t0
        eta = elapsed*(NP-p_idx-1)/(p_idx+1)
        print(f"  {p_idx+1}/{NP} ({pct:.0f}%)  {elapsed:.1f}s  ETA {eta:.0f}s", flush=True)
print(f"  Canonicalization done ({time.time()-t0:.1f}s)", flush=True)

# Group GS tilings into SC iso classes
sc_groups = defaultdict(list)
for gi in range(NGS):
    sc_groups[canon_gs[gi].tobytes()].append(gi)

print(f"  SC iso classes: {len(sc_groups)} (expected 2752)", flush=True)

# Get one representative GS tiling per SC class
sc_reps = list(sc_groups.keys())  # canonical form bytes
sc_rep_gs_idx = [members[0] for members in sc_groups.values()]  # GS indices
sc_A_reps = A_gs[sc_rep_gs_idx]   # (NC_SC, n, n) — one per SC class
NC_SC = len(sc_reps)
print(f"  SC representatives: {NC_SC}", flush=True)

# ─── COMPUTE CERT (SCORE SEQ) FOR EACH SC CLASS ────────────────────────────────
# cert = sorted out-degree sequence of the representative
sc_certs = sc_A_reps.sum(axis=2)  # (NC_SC, n) out-degrees
sc_certs_sorted = np.sort(sc_certs, axis=1)  # (NC_SC, n) sorted
print(f"  Score certs computed", flush=True)

# ─── MAIN SCSF CHECK ─────────────────────────────────────────────────────────
# For each SC class C (with rep A_C):
#   For each σ ∈ S_n:
#     Compute A_sigma = σ(A_C) (apply permutation to adjacency)
#     Compute A_tilde = flip tile arcs of A_sigma
#     Check if out-degrees(A_tilde) sorted == out-degrees(A_C) sorted (cert check)
#     If cert passes: check full isomorphism
#   If any σ gives A_tilde ≅ A_C → C is SCSF

print(f"\nMain SCSF check ({NC_SC} SC classes × {NP} perms)...", flush=True)

scsf_classes = []
t0 = time.time()

# Process all SC classes in parallel, vectorizing over permutations
# For each class, we apply all NP permutations at once.
# Memory: NP × n² = 362880 × 81 = 29.4M uint8 = 29.4MB per class → need to process
# one class at a time to keep memory manageable.

for cls_idx in range(NC_SC):
    A_C = sc_A_reps[cls_idx]   # (n, n)
    A_C_flat = A_C.reshape(n*n)  # (n²,)
    target_cert = sc_certs_sorted[cls_idx]  # (n,) sorted out-degrees

    # Apply all NP permutations: (NP, n²)
    A_sigma_flat = A_C_flat[perm_flat]  # (NP, n²) — vectorized!

    # Flip tile arcs: for each pair (i,j) with |i-j|>=2, swap columns ij and ji
    A_tilde_flat = A_sigma_flat.copy()
    # Swap tp_ij_flat <-> tp_ji_flat columns
    tmp = A_tilde_flat[:, tp_ij_flat].copy()
    A_tilde_flat[:, tp_ij_flat] = A_tilde_flat[:, tp_ji_flat]
    A_tilde_flat[:, tp_ji_flat] = tmp

    # Cert check: sorted out-degrees of A_tilde
    A_tilde_3d = A_tilde_flat.reshape(NP, n, n)
    out_deg_tilde = A_tilde_3d.sum(axis=2)        # (NP, n)
    out_sorted_tilde = np.sort(out_deg_tilde, axis=1)  # (NP, n)
    cert_match = (out_sorted_tilde == target_cert[None, :]).all(axis=1)  # (NP,)

    n_cert_match = cert_match.sum()
    if n_cert_match == 0:
        if (cls_idx+1) % 200 == 0:
            print(f"  {cls_idx+1}/{NC_SC}  {time.time()-t0:.1f}s  scsf={len(scsf_classes)}", flush=True)
        continue

    # Full isomorphism check for cert-passing permutations
    matching_perms = np.where(cert_match)[0]
    found_scsf = False

    for p_idx in matching_perms:
        A_tilde = A_tilde_flat[p_idx].reshape(n, n)
        # Check if A_tilde ≅ A_C using stronger cert first
        # Level 2: in-degree sequence
        in_deg_tilde = A_tilde.sum(axis=0)
        in_sorted_tilde = tuple(sorted(in_deg_tilde.tolist()))
        in_sorted_C = tuple(sorted(A_C.sum(axis=0).tolist()))
        if in_sorted_tilde != in_sorted_C:
            continue
        # Level 3: 2nd order cert
        out = out_deg_tilde[p_idx]
        inp = in_deg_tilde
        out_C = A_C.sum(axis=1)
        inp_C = A_C.sum(axis=0)
        # Simple 2nd order: for each vertex, sorted (out,in) of out-neighbors
        def get_cert2(A_mat, out_v, in_v):
            vertex_sigs = []
            for v in range(n):
                out_nbrs = np.where(A_mat[v] == 1)[0]
                sig = tuple(sorted([(int(out_v[w]), int(in_v[w])) for w in out_nbrs]))
                vertex_sigs.append((int(out_v[v]), int(in_v[v]), sig))
            return tuple(sorted(vertex_sigs))
        c2_tilde = get_cert2(A_tilde, out, inp)
        c2_C = get_cert2(A_C, out_C, inp_C)
        if c2_tilde != c2_C:
            continue
        # Full brute-force isomorphism
        A_C_bytes = A_C.tobytes()
        is_iso = False
        for sigma2 in perms_list:
            # Check if sigma2(A_tilde) == A_C
            check = True
            for i in range(n):
                for j in range(n):
                    if A_tilde[sigma2[i], sigma2[j]] != A_C[i, j]:
                        check = False
                        break
                if not check: break
            if check:
                is_iso = True
                break
        if is_iso:
            found_scsf = True
            break

    if found_scsf:
        scores = tuple(sorted(A_C.sum(axis=1).tolist(), reverse=True))
        gs_count = len(sc_groups[sc_reps[cls_idx]])
        scsf_classes.append({
            'cls_idx': cls_idx,
            'scores': scores,
            'gs_count': gs_count,
        })
        print(f"  → SCSF class found! cls={cls_idx}, scores={scores}, gs={gs_count}  ({time.time()-t0:.1f}s)", flush=True)

    if (cls_idx+1) % 200 == 0:
        print(f"  {cls_idx+1}/{NC_SC}  {time.time()-t0:.1f}s  scsf={len(scsf_classes)}", flush=True)

print(f"\n{'='*60}")
print(f"SCSF(9) = {len(scsf_classes)}")
print(f"{'='*60}")
for i, c in enumerate(scsf_classes):
    print(f"  Class {i+1}: scores={c['scores']}, gs_tilings={c['gs_count']}")

print(f"\nTotal time: {time.time()-t_global:.1f}s")
