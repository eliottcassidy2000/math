"""
scsf_n9_fast.py — Fast SCSF(9) computation using certificate hashing.

KEY SPEEDUP: Instead of canonicalizing all 131072 tilings with 362880 permutations,
use a multi-level hash certificate to filter candidates, then only run full
brute-force canonicalization on the small number of matching pairs.

Certificate: degree_seq + 2nd_order_degree_seq (WL refinement, 2 rounds).
For most non-isomorphic pairs, the certificate will differ → fast reject.

Only pairs surviving the certificate check get full isomorphism testing.

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
tm = [tile_idx[(n-y+1, n-x+1)] for x,y in tiles]

self_paired = [i for i in range(m) if tm[i] == i]
gs_pairs_idx = [(i, tm[i]) for i in range(m) if tm[i] > i]
free_bits = len(self_paired) + len(gs_pairs_idx)

print(f"n={n}, m={m}, free_bits={free_bits}, GS={1<<free_bits}", flush=True)

# ─── ENUMERATE GS TILINGS ─────────────────────────────────────────────────────
NGS = 1 << free_bits
free_to_tiles = [[sp] for sp in self_paired] + [[i,j] for i,j in gs_pairs_idx]
fp_range = np.arange(free_bits, dtype=np.int32)
free_idx = np.arange(NGS, dtype=np.uint32)
free_params = ((free_idx[:, None] >> fp_range[None, :]) & 1).astype(np.uint8)
gs_masks = np.zeros(NGS, dtype=np.uint32)
for fbi, tile_list in enumerate(free_to_tiles):
    for ti in tile_list:
        gs_masks |= (free_params[:, fbi].astype(np.uint32) << ti)
flip_all = np.uint32((1 << m) - 1)
flip_masks = gs_masks ^ flip_all
print(f"  Built {NGS} GS masks", flush=True)

# ─── BUILD ADJACENCY MATRICES ──────────────────────────────────────────────────
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
    return A  # (NC, n, n)

print("Building adjacency matrices...", flush=True)
t0 = time.time()
A_gs = build_A_batch(gs_masks)    # (NGS, n, n)
A_flip = build_A_batch(flip_masks) # (NGS, n, n)
print(f"  Done ({time.time()-t0:.1f}s)", flush=True)

# ─── CERTIFICATE: WL-STYLE HASH ────────────────────────────────────────────────
def compute_cert(A_batch):
    """
    Compute a hashable certificate for each tournament in A_batch.
    Uses 3 rounds of WL refinement based on (out-deg, in-deg, neighbor-scores).
    Returns a list of frozensets/tuples usable as dict keys.
    """
    NC = len(A_batch)
    # out-degree, in-degree
    out_deg = A_batch.sum(axis=2)  # (NC, n) = row sums
    in_deg = A_batch.sum(axis=1)   # (NC, n) = col sums

    # Level 0: (out_deg_v, in_deg_v) for each vertex
    # Level 1: (out_deg_v, in_deg_v, sorted(out_deg of out-neighbors), sorted(in_deg of in-neighbors))
    certs = []
    for i in range(NC):
        out = out_deg[i]  # (n,)
        inp = in_deg[i]   # (n,)
        A = A_batch[i]    # (n,n)
        # For each vertex v, compute:
        # out-neighbor out-degs, in-neighbor in-degs
        vertex_feats = []
        for v in range(n):
            out_nbrs = np.where(A[v] == 1)[0]  # v→w
            in_nbrs = np.where(A[:, v] == 1)[0]  # w→v
            feat = (
                int(out[v]), int(inp[v]),
                tuple(sorted(out[out_nbrs].tolist())),
                tuple(sorted(inp[in_nbrs].tolist())),
                tuple(sorted(out[in_nbrs].tolist())),
            )
            vertex_feats.append(feat)
        cert = tuple(sorted(vertex_feats))
        certs.append(cert)
    return certs

print("Computing certificates for GS tilings...", flush=True)
t0 = time.time()
# Process in chunks for memory
CHUNK = 4096
cert_gs = []
cert_flip = []
for start in range(0, NGS, CHUNK):
    end = min(start + CHUNK, NGS)
    cert_gs.extend(compute_cert(A_gs[start:end]))
    cert_flip.extend(compute_cert(A_flip[start:end]))
    if (start // CHUNK) % 4 == 0:
        print(f"  {end}/{NGS} ({100*end/NGS:.0f}%) {time.time()-t0:.1f}s", flush=True)
print(f"  Certificates done ({time.time()-t0:.1f}s)", flush=True)

# ─── FAST FILTER: CERT MATCH ──────────────────────────────────────────────────
# GS tiling i might be SCSF if cert_gs[i] == cert_flip[i]
print("\nFiltering by certificate...", flush=True)
cert_match_idx = [i for i in range(NGS) if cert_gs[i] == cert_flip[i]]
print(f"  Certificate matches: {len(cert_match_idx)} out of {NGS}", flush=True)

# ─── FULL ISOMORPHISM FOR CERT-MATCHING PAIRS ─────────────────────────────────
# For each cert-matching pair (A_gs[i], A_flip[i]), check full isomorphism.
print(f"\nFull isomorphism check for {len(cert_match_idx)} pairs...", flush=True)

def is_isomorphic(A1, A2):
    """Brute-force isomorphism: check all n! permutations."""
    A1 = A1.reshape(n, n)
    A2 = A2.reshape(n, n)
    # Quick: compare score sequences first
    out1 = sorted(A1.sum(axis=1).tolist())
    out2 = sorted(A2.sum(axis=1).tolist())
    if out1 != out2: return False
    # Try all n! permutations
    for sigma in permutations(range(n)):
        # A1_sigma[sigma[i], sigma[j]] = A1[i, j]
        # Equivalently: A1 permuted by sigma as perm_A[i,j] = A1[sigma^{-1}[i], sigma^{-1}[j]]
        # Or: check if A2[sigma[i], sigma[j]] == A1[i,j] for all i,j
        ok = True
        for i in range(n):
            for j in range(n):
                if A2[sigma[i], sigma[j]] != A1[i, j]:
                    ok = False
                    break
            if not ok: break
        if ok: return True
    return False

def is_isomorphic_fast(A1, A2):
    """Faster: use color refinement first, then constrained brute-force."""
    A1 = A1.reshape(n, n)
    A2 = A2.reshape(n, n)
    # Multi-round WL coloring
    out1 = A1.sum(axis=1); in1 = A1.sum(axis=0)
    out2 = A2.sum(axis=1); in2 = A2.sum(axis=0)

    # Build vertex color: (out-deg, in-deg, sorted out-nbr-out-degs, ...)
    def get_colors(A):
        out = A.sum(axis=1); inp = A.sum(axis=0)
        colors = []
        for v in range(n):
            out_nbrs = np.where(A[v]==1)[0]
            in_nbrs = np.where(A[:,v]==1)[0]
            c = (int(out[v]), int(inp[v]),
                 tuple(sorted(out[out_nbrs])), tuple(sorted(inp[in_nbrs])))
            colors.append(c)
        return colors

    c1 = get_colors(A1); c2 = get_colors(A2)
    # Check color histogram
    from collections import Counter
    if Counter(c1) != Counter(c2): return False

    # Build constrained permutations: sigma[v] must have same color as v
    color_groups1 = defaultdict(list)
    color_groups2 = defaultdict(list)
    for v, c in enumerate(c1): color_groups1[c].append(v)
    for v, c in enumerate(c2): color_groups2[c].append(v)

    # Enumerate constrained permutations
    # For each color class, try all assignments from class2 to class1
    colors_unique = list(color_groups1.keys())
    # Build partial permutations recursively
    def try_perms(color_idx, sigma):
        if color_idx == len(colors_unique):
            # Full sigma: check all arcs
            for i in range(n):
                for j in range(n):
                    if A2[sigma[i], sigma[j]] != A1[i, j]:
                        return False
            return True
        c = colors_unique[color_idx]
        verts1 = color_groups1[c]
        verts2 = color_groups2[c]
        for perm in permutations(verts2):
            for v1, v2 in zip(verts1, perm):
                sigma[v1] = v2
            if try_perms(color_idx + 1, sigma):
                return True
        return False

    sigma = [0] * n
    return try_perms(0, sigma)

scsf_gs_indices = []
t0 = time.time()
for idx, i in enumerate(cert_match_idx):
    A1 = A_gs[i]
    A2 = A_flip[i]
    if is_isomorphic_fast(A1, A2):
        scsf_gs_indices.append(i)
    if (idx+1) % 100 == 0:
        print(f"  Checked {idx+1}/{len(cert_match_idx)}: {len(scsf_gs_indices)} SCSF so far  {time.time()-t0:.1f}s", flush=True)
print(f"  Full iso check done: {len(scsf_gs_indices)} SCSF GS tilings  ({time.time()-t0:.1f}s)", flush=True)

# ─── GROUP SCSF GS TILINGS INTO SC ISO CLASSES ────────────────────────────────
# Two GS tilings are in the same SC class iff their tournaments are isomorphic.
# Use cert-based grouping + verification.
print(f"\nGrouping {len(scsf_gs_indices)} SCSF GS tilings into SC iso classes...", flush=True)

if scsf_gs_indices:
    # Use certs to group
    scsf_gs_certs = [cert_gs[i] for i in scsf_gs_indices]
    groups = defaultdict(list)
    for k, cert in zip(scsf_gs_indices, scsf_gs_certs):
        groups[cert].append(k)

    # Verify within each cert group that they're all in the same iso class
    scsf_classes = []
    for cert, members in groups.items():
        rep = members[0]
        same_class = [rep]
        other = members[1:]
        while other:
            new_other = []
            for i in other:
                if is_isomorphic_fast(A_gs[rep], A_gs[i]):
                    same_class.append(i)
                else:
                    new_other.append(i)
            if not new_other: break
            # Split needed (multiple iso classes with same cert)
            rep2 = new_other[0]
            sub = [rep2]
            rest = []
            for i in new_other[1:]:
                if is_isomorphic_fast(A_gs[rep2], A_gs[i]):
                    sub.append(i)
                else:
                    rest.append(i)
            scsf_classes.append({'gs_members': sub, 'gs_count': len(sub)})
            other = rest
        scsf_classes.append({'gs_members': same_class, 'gs_count': len(same_class)})

    # Compute score sequences for each class
    for c in scsf_classes:
        rep = c['gs_members'][0]
        c['scores'] = tuple(sorted(A_gs[rep].sum(axis=1).tolist(), reverse=True))

    print(f"\n{'='*60}")
    print(f"SCSF(9) = {len(scsf_classes)}")
    print(f"{'='*60}")
    for i, c in enumerate(scsf_classes):
        print(f"  Class {i+1}: gs_tilings={c['gs_count']}, scores={c['scores']}")
else:
    print(f"\nSCSF(9) = 0")

print(f"\nTotal time: {time.time()-t_global:.1f}s")
