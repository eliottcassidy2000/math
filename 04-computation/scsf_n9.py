"""
scsf_n9.py — Compute SCSF(9): self-complementary AND self-flip tournament classes at n=9.

Strategy:
  1. Batch-process all 2^28 tilings to find SF-score candidates (near-regular score filter)
  2. Canonicalize survivors using all 9! permutations
  3. Check SC and SF conditions

Score conventions for n=9:
  Flip (tile-complement) maps:
    s̃(1) = 7 - s(1)    [base-path sink]
    s̃(9) = 9 - s(9)    [base-path source]
    s̃(v) = 8 - s(v)    [middle vertices v=2..8]

Near-regularity conjecture: SC∩SF classes have all scores in {3,4,5} (sum=36=9×4).

oracle-2026-05-11
"""

import numpy as np
import sys, time
from itertools import permutations
from collections import defaultdict

sys.stdout.reconfigure(encoding='utf-8')
t_global = time.time()

n = 9
# Tile arcs: (x,y) with y in 1..n-2, x in n down to y+2
tiles = [(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]
m = len(tiles)
assert m == 28, f"Expected m=28, got {m}"

# Base path: n->n-1->...->2->1 (0-indexed verts: n-1 -> n-2 -> ... -> 0)
verts = list(range(n, 0, -1))  # [9,8,7,...,1]
xi_arr = np.array([verts.index(x) for x,y in tiles], dtype=np.int32)
yi_arr = np.array([verts.index(y) for x,y in tiles], dtype=np.int32)

# TRANS_MAP for GS check: tile (x,y) -> (n-y+1, n-x+1)
tile_idx = {t:i for i,t in enumerate(tiles)}
tm = [tile_idx[(n-y+1, n-x+1)] for x,y in tiles]
gs_pairs = [(i, tm[i]) for i in range(m) if tm[i] > i]

print(f"n={n}, m={m}, 2^m={1<<m}", flush=True)
print(f"  Tile arc structure: {len(tiles)} tiles, {len(gs_pairs)} GS-pair constraints", flush=True)

# ─── SCORE COMPUTATION VIA BIT EXTRACTION ─────────────────────────────────────
# For each tile k: bit=0 means arc xi_arr[k]→yi_arr[k] (xi gets +1 from base wins)
# Actually: bit=0 → arc FROM xi_arr[k] TO yi_arr[k] (xi_arr[k] BEATS yi_arr[k])
#           bit=1 → arc FROM yi_arr[k] TO xi_arr[k] (yi_arr[k] BEATS xi_arr[k])
# Score of vertex v = (# outgoing arcs from v)
# Base path contribution: vertex at position k in base path (k=0 is vertex n, k=n-1 is vertex 1)
#   verts[k] → verts[k+1]: so verts[k] beats verts[k+1]
#   Base path score for vertex i: = 1 if i is NOT the last vertex (vertex 1, index n-1)

# Base path score: verts[0]=n beats verts[1]=n-1, ..., verts[n-2]=2 beats verts[n-1]=1
# So verts[k] gets +1 for k=0..n-2 (i.e., all except the sink vertex=1=index n-1)
base_score = np.zeros(n, dtype=np.int16)
for k in range(n-1):
    base_score[k] += 1  # verts[k] beats verts[k+1]

# ─── BATCH SCORE FILTER ───────────────────────────────────────────────────────
BATCH = 1 << 22  # 4M tilings per batch
NT = 1 << m      # 268M total

print(f"\nPhase 1: Batch score filter ({NT} tilings in batches of {BATCH})...", flush=True)

# Precompute bit extraction: for tile k, bit = (mask >> k) & 1
# score[v] = base_score[v] + sum over tiles k where yi_arr[k]==v and bit_k==1
#                           + sum over tiles k where xi_arr[k]==v and bit_k==0

# Precompute: for each vertex v, which tiles contribute and how
# tile k contributes +1 to xi_arr[k] if bit_k==0, +1 to yi_arr[k] if bit_k==1
# Equivalently: score[v] = base_score[v] + (# tiles k with xi_arr[k]==v and bit_k==0)
#                                         + (# tiles k with yi_arr[k]==v and bit_k==1)

# For batch processing, use matrix approach:
# Build (m,) arrays xi_arr, yi_arr
# For a batch of masks shape (B,), mask bits are (B, m) via (masks[:,None] >> arange(m)) & 1

# SF SCORE TEST: sorted(T scores) == sorted(T̃ scores)
# T̃ scores:
#   s̃[0] = 7 - s[0]   (vertex n = index 0 = base-path source)
#   s̃[8] = 9 - s[8]   (vertex 1 = index 8 = base-path sink)
#   s̃[v] = 8 - s[v]   (middle vertices v=1..7)
# NEAR-REGULARITY FILTER: all scores in {3,4,5}, sum=36

candidates = []
k_range = np.arange(m, dtype=np.int32)

# Precompute x_mat, y_mat for vectorized score computation:
# score = base_score + (1-bits) @ x_mat + bits @ y_mat
# where x_mat[k,v] = 1 iff xi_arr[k]==v, y_mat[k,v] = 1 iff yi_arr[k]==v
x_mat = np.zeros((m, n), dtype=np.int16)
y_mat = np.zeros((m, n), dtype=np.int16)
for k in range(m):
    x_mat[k, int(xi_arr[k])] = 1
    y_mat[k, int(yi_arr[k])] = 1

t0 = time.time()
for batch_start in range(0, NT, BATCH):
    batch_end = min(batch_start + BATCH, NT)
    B = batch_end - batch_start
    masks = np.arange(batch_start, batch_end, dtype=np.uint32)

    # Extract bits: (B, m)
    bits = ((masks[:, None] >> k_range[None, :]) & 1).astype(np.int16)

    # Vectorized score computation: (B, n)
    scores = base_score[None, :] + (1 - bits) @ x_mat + bits @ y_mat

    # Near-regularity filter: all scores in [3,5]
    near_reg = (scores >= 3).all(axis=1) & (scores <= 5).all(axis=1)

    if not near_reg.any():
        if (batch_start // BATCH) % 16 == 0:
            pct = 100*(batch_end)/NT
            print(f"  {batch_end//1000000}M/{NT//1000000}M ({pct:.0f}%) {time.time()-t0:.0f}s candidates={len(candidates)}", flush=True)
        continue

    # SF score test on near-regular candidates
    nr_masks = masks[near_reg]
    nr_scores = scores[near_reg]  # (K, 9)

    # Compute T̃ scores
    flip_scores = 8 - nr_scores.copy()  # middle vertices default: 8-s
    flip_scores[:, 0] = 7 - nr_scores[:, 0]   # source vertex (9, index 0): 7-s
    flip_scores[:, 8] = 9 - nr_scores[:, 8]   # sink vertex (1, index 8): 9-s

    # Sort both and compare
    s_sorted = np.sort(nr_scores, axis=1)
    f_sorted = np.sort(flip_scores, axis=1)
    sf_score_pass = (s_sorted == f_sorted).all(axis=1)

    if sf_score_pass.any():
        cands = nr_masks[sf_score_pass].tolist()
        candidates.extend(cands)

    if (batch_start // BATCH) % 8 == 0:
        pct = 100*(batch_end)/NT
        print(f"  {batch_end//1000000}M/{NT//1000000}M ({pct:.0f}%) {time.time()-t0:.0f}s candidates={len(candidates)}", flush=True)

print(f"\nScore filter done: {len(candidates)} candidates  ({time.time()-t0:.1f}s)", flush=True)

if len(candidates) == 0:
    print("No candidates found! Exiting.")
    sys.exit(0)

# ─── BUILD ADJACENCY MATRICES FOR CANDIDATES ──────────────────────────────────
print(f"\nPhase 2: Building adjacency matrices for {len(candidates)} candidates...", flush=True)
NC = len(candidates)
masks_cand = np.array(candidates, dtype=np.uint32)

# Build A_flat (NC, n*n) for candidates
A_flat = np.zeros((NC, n*n), dtype=np.uint8)
# Base path: verts[k] → verts[k+1] for k=0..n-2
for k in range(n-1):
    A_flat[:, verts[k]*n + verts[k+1]] = 1  # Wait, need 0-indexed positions!
    # verts[k] is the vertex label (9,8,...,1), we need 0-indexed position
    # Actually, let me use position index = n - label (so label n=9 → index 0, label 1 → index 8)
    # verts[k] is already position index (since verts=[9,8,...,1] and index 0=vert 9, etc.)
    # BUT: A_flat[mask, i*n+j] = 1 means position i beats position j
    # verts[k]=position k in the verts list → position index k
    A_flat[:, k*n + (k+1)] = 1  # position k beats position k+1

# Tile arcs
bits_cand = ((masks_cand[:, None] >> k_range[None, :]) & 1).astype(np.uint8)
for k in range(m):
    xk, yk = int(xi_arr[k]), int(yi_arr[k])
    b = bits_cand[:, k]  # (NC,)
    # bit=0 → position xk beats position yk
    # bit=1 → position yk beats position xk
    A_flat[:, xk*n + yk] += (1 - b)
    A_flat[:, yk*n + xk] += b

# Verify: A_flat should have exactly 1 in each off-diagonal (i,j) or (j,i) pair
# Quick sanity check on first candidate
A0 = A_flat[0].reshape(n, n)
for i in range(n):
    for j in range(i+1, n):
        assert A0[i,j] + A0[j,i] == 1, f"Bad tournament at ({i},{j})"
print("  Adjacency matrix sanity check passed.", flush=True)

# ─── COMPLEMENT ADJACENCY ─────────────────────────────────────────────────────
# T^op: A^op[i][j] = A[j][i] = transpose (off-diagonal)
trans_idx = np.array([j*n+i for i in range(n) for j in range(n)], dtype=np.int32)
# But diagonal is 0, so transpose is fine
comp_A_flat = A_flat[:, trans_idx]

# ─── PERMUTATION ACTION ───────────────────────────────────────────────────────
print(f"\nPhase 3: Precomputing permutation flat indices ({n}! = {np.math.factorial(n)})...", flush=True)
perms_list = list(permutations(range(n)))
NP = len(perms_list)
print(f"  {NP} permutations", flush=True)

perm_flat = np.zeros((NP, n*n), dtype=np.int32)
for p_idx, sigma in enumerate(perms_list):
    for i in range(n):
        for j in range(n):
            perm_flat[p_idx, i*n+j] = sigma[i]*n + sigma[j]

# ─── CANONICALIZATION ─────────────────────────────────────────────────────────
print(f"\nPhase 4: Canonicalizing {NC} candidates with {NP} permutations...", flush=True)

def compute_canonical_batch(src, label=""):
    """src: (NC, n²) uint8. Returns (NC, n²) uint8 canonical forms."""
    canon = np.full((NC, n*n), 255, dtype=np.uint8)
    t0 = time.time()
    for p_idx in range(NP):
        perm_A = src[:, perm_flat[p_idx]]
        diff = perm_A != canon
        has_diff = diff.any(axis=1)
        first = np.argmax(diff, axis=1)
        ridx = np.arange(NC)
        is_less = (perm_A[ridx, first] < canon[ridx, first]) & has_diff
        canon[is_less] = perm_A[is_less]
        if (p_idx+1) % 36288 == 0:  # every ~10%
            pct = 100*(p_idx+1)/NP
            print(f"  {label}: {p_idx+1}/{NP} ({pct:.0f}%)  {time.time()-t0:.1f}s", flush=True)
    return canon

t0 = time.time()
canon = compute_canonical_batch(A_flat, "canon")
print(f"  canonical done  {time.time()-t0:.1f}s", flush=True)

t0 = time.time()
comp_canon = compute_canonical_batch(comp_A_flat, "comp_canon")
print(f"  comp_canonical done  {time.time()-t0:.1f}s", flush=True)

# ─── GROUP INTO ISO CLASSES ───────────────────────────────────────────────────
print(f"\nPhase 5: Grouping into iso classes...", flush=True)
groups = defaultdict(list)
for i in range(NC):
    groups[canon[i].tobytes()].append(i)

classes = []
for sig_bytes, members in groups.items():
    rep = members[0]
    sc = np.array_equal(canon[rep], comp_canon[rep])
    # Compute scores for this rep
    score_v = A_flat[rep].reshape(n,n).sum(axis=1).tolist()  # out-degree = score
    classes.append({
        'sig': sig_bytes,
        'rep': rep,
        'members': members,
        'nt': len(members),
        'sc': sc,
        'scores': score_v,
    })

# ─── SF CHECK ─────────────────────────────────────────────────────────────────
# A class is SF if any of its members has flip (T̃) isomorphic to T.
# flip(mask) = ~mask & ((1<<m)-1)
# T̃ adjacency = flip the tile bits: bit=0→1, bit=1→0 for tile arcs, base path unchanged.
# T̃ A_flat = recompute with flipped bits = comp of tile part.

# Build T̃ A_flat for all candidates:
flip_bits = 1 - bits_cand  # (NC, m)
tilde_A_flat = np.zeros((NC, n*n), dtype=np.uint8)
for k in range(n-1):
    tilde_A_flat[:, k*n + (k+1)] = 1  # base path unchanged
for k in range(m):
    xk, yk = int(xi_arr[k]), int(yi_arr[k])
    b = flip_bits[:, k]
    tilde_A_flat[:, xk*n + yk] += (1 - b)
    tilde_A_flat[:, yk*n + xk] += b

print(f"  Computing canonical of T̃...", flush=True)
t0 = time.time()
tilde_canon = compute_canonical_batch(tilde_A_flat, "tilde")
print(f"  tilde_canonical done  {time.time()-t0:.1f}s", flush=True)

# Map canonical sig → class index
sig_to_cls = {}
for ci, c in enumerate(classes):
    sig_to_cls[c['sig']] = ci

# Check SF: class i is SF if tilde_canon[some rep] matches canon of same class
sf_by_class = [False] * len(classes)
# For each candidate, check if its tilde_canon matches its own class canon
for i in range(NC):
    cls_i = sig_to_cls.get(canon[i].tobytes())
    if cls_i is None: continue
    tilde_sig = tilde_canon[i].tobytes()
    # tilde_sig should be in sig_to_cls
    if tilde_sig in sig_to_cls:
        tilde_cls = sig_to_cls[tilde_sig]
        if tilde_cls == cls_i:
            sf_by_class[cls_i] = True

for ci, c in enumerate(classes):
    c['sf'] = sf_by_class[ci]

# ─── RESULTS ─────────────────────────────────────────────────────────────────
sc_classes = [c for c in classes if c['sc']]
sf_classes = [c for c in classes if c['sf']]
scsf_classes = [c for c in classes if c['sc'] and c['sf']]

print(f"\n{'='*60}")
print(f"RESULTS FOR n=9 (near-regular score filter)")
print(f"{'='*60}")
print(f"  Score-filter candidates: {NC}")
print(f"  Iso classes found:       {len(classes)}")
print(f"  SC classes:              {len(sc_classes)}")
print(f"  SF classes:              {len(sf_classes)}")
print(f"  SC∩SF (kernel):          {len(scsf_classes)}")
print(f"\n  SCSF(9) = {len(scsf_classes)}")
print(f"\nKernel details:")
for c in scsf_classes:
    scores_sorted = sorted(c['scores'], reverse=True)
    # GS check
    gs = sum(1 for mask_i in c['members']
             if all(((masks_cand[mask_i]>>i)&1)==((masks_cand[mask_i]>>j)&1)
                    for i,j in gs_pairs))
    print(f"  #tilings={c['nt']}, scores={scores_sorted}, gs={gs}")

print(f"\nSC class score distribution:")
score_seqs = {}
for c in sc_classes:
    s = tuple(sorted(c['scores'], reverse=True))
    score_seqs[s] = score_seqs.get(s,0) + 1
for s, cnt in sorted(score_seqs.items()):
    print(f"  {s}: {cnt} classes")

print(f"\nTotal time: {time.time()-t_global:.1f}s")
