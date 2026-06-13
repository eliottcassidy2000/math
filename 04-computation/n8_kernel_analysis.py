"""
n8_kernel_analysis.py — Detailed analysis of the 5 SC∩SF kernel classes at n=8.

Goal: Characterize all 5 kernel classes — bl value, scores, blueself/blackself status,
      automorphism group size, and the structure of the flip-pair graph.

From previous session:
  - 4 blueself classes: bl ∈ {37, 37, 39, 39}
  - 1 blackself class: bl unknown
  - Kernel class 5 is SC∩SF via NON-GS tilings (blackself)

oracle-2026-05-11
"""

import numpy as np
import sys, time
from itertools import permutations
from collections import defaultdict

sys.stdout.reconfigure(encoding='utf-8')
t_global = time.time()

n = 8
tiles = [(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]
m = len(tiles)
assert m == 21

verts = list(range(n, 0, -1))  # [8,7,...,1]
xi_arr = np.array([verts.index(x) for x,y in tiles], dtype=np.int32)
yi_arr = np.array([verts.index(y) for x,y in tiles], dtype=np.int32)

tile_idx = {t:i for i,t in enumerate(tiles)}
tm = [tile_idx[(n-y+1, n-x+1)] for x,y in tiles]
gs_pairs = [(i, tm[i]) for i in range(m) if tm[i] > i]

NT = 1 << m  # 2^21 = 2M
k_range = np.arange(m, dtype=np.int32)

print(f"n={n}, m={m}, NT={NT}", flush=True)

# ─── SCORE FILTER (SF score test, balanced {3,4,5} or {4,4,4,4,3,3,3,3}) ──────
# For n=8 SC∩SF classes: scores in {3,4} (near-regular, mean=3.5, but sum=28=C(8,2))
# Actually: sum of scores = C(8,2)=28, mean=28/8=3.5, so scores in {3,4} exactly
# SF score formula for n=8:
#   s̃(1)=s̃(pos 7) = 6 - s(pos 7)  [sink]
#   s̃(8)=s̃(pos 0) = 8 - s(pos 0)  [source]
#   s̃(v)=s̃(pos v) = 7 - s(pos v)  [middle pos 1..6]

x_mat = np.zeros((m, n), dtype=np.int16)
y_mat = np.zeros((m, n), dtype=np.int16)
for k in range(m):
    x_mat[k, int(xi_arr[k])] = 1
    y_mat[k, int(yi_arr[k])] = 1

base_score = np.zeros(n, dtype=np.int16)
for k in range(n-1):
    base_score[k] += 1  # pos k beats pos k+1

print("Phase 1: Score filter...", flush=True)
t0 = time.time()
BATCH = 1 << 20  # 1M
masks_all = np.arange(NT, dtype=np.uint32)
bits_all = ((masks_all[:, None] >> k_range[None, :]) & 1).astype(np.int16)
scores_all = base_score[None, :] + (1 - bits_all) @ x_mat + bits_all @ y_mat

# Balanced score filter: sorted scores = (4,4,4,4,3,3,3,3)
target = np.array([4,4,4,4,3,3,3,3], dtype=np.int16)
scores_sorted = np.sort(scores_all, axis=1)[:, ::-1]
balanced = (scores_sorted == target).all(axis=1)
print(f"  Balanced candidates: {balanced.sum()} ({time.time()-t0:.1f}s)", flush=True)

# SF score test on balanced candidates
bal_masks = masks_all[balanced]
bal_scores = scores_all[balanced]

flip_scores = 7 - bal_scores.copy()  # middle default
flip_scores[:, 0] = 8 - bal_scores[:, 0]   # source (pos 0)
flip_scores[:, 7] = 6 - bal_scores[:, 7]   # sink (pos 7)

s_sorted = np.sort(bal_scores, axis=1)
f_sorted = np.sort(flip_scores, axis=1)
sf_pass = (s_sorted == f_sorted).all(axis=1)
print(f"  SF-score-passing: {sf_pass.sum()}", flush=True)

cand_masks = bal_masks[sf_pass]
cand_bits = bits_all[balanced][sf_pass].astype(np.uint8)
cand_scores = bal_scores[sf_pass]
NC = len(cand_masks)

# ─── BUILD ADJACENCY MATRICES ──────────────────────────────────────────────────
print(f"\nPhase 2: Building adjacency matrices for {NC} candidates...", flush=True)
A_flat = np.zeros((NC, n*n), dtype=np.uint8)
for k in range(n-1):
    A_flat[:, k*n+(k+1)] = 1  # base path
for k in range(m):
    xk, yk = int(xi_arr[k]), int(yi_arr[k])
    b = cand_bits[:, k]
    A_flat[:, xk*n+yk] += (1-b)
    A_flat[:, yk*n+xk] += b

comp_A_flat = A_flat[:, np.array([j*n+i for i in range(n) for j in range(n)], dtype=np.int32)]

# T̃ adjacency (flip all tile arcs)
tilde_A_flat = np.zeros((NC, n*n), dtype=np.uint8)
for k in range(n-1):
    tilde_A_flat[:, k*n+(k+1)] = 1
flip_bits = 1 - cand_bits
for k in range(m):
    xk, yk = int(xi_arr[k]), int(yi_arr[k])
    b = flip_bits[:, k]
    tilde_A_flat[:, xk*n+yk] += (1-b)
    tilde_A_flat[:, yk*n+xk] += b

# ─── PERMUTATION SETUP ─────────────────────────────────────────────────────────
print(f"\nPhase 3: Building permutation indices (n!=40320)...", flush=True)
perms_list = list(permutations(range(n)))
NP = len(perms_list)
perm_flat = np.zeros((NP, n*n), dtype=np.int32)
for p_idx, sigma in enumerate(perms_list):
    for i in range(n):
        for j in range(n):
            perm_flat[p_idx, i*n+j] = sigma[i]*n + sigma[j]
print(f"  Done ({NP} perms)", flush=True)

# ─── CANONICALIZATION ─────────────────────────────────────────────────────────
def compute_canonical(src, label=""):
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
        if (p_idx+1) % 8000 == 0:
            print(f"  {label}: {p_idx+1}/{NP} ({100*(p_idx+1)/NP:.0f}%)  {time.time()-t0:.1f}s", flush=True)
    return canon

print(f"\nPhase 4: Canonicalizing...", flush=True)
t0 = time.time()
canon = compute_canonical(A_flat, "T")
print(f"  T canonical done ({time.time()-t0:.1f}s)", flush=True)
comp_canon = compute_canonical(comp_A_flat, "T^op")
print(f"  T^op canonical done ({time.time()-t0:.1f}s)", flush=True)
tilde_canon = compute_canonical(tilde_A_flat, "T~")
print(f"  T~ canonical done ({time.time()-t0:.1f}s)", flush=True)

# ─── GROUP INTO ISO CLASSES ────────────────────────────────────────────────────
print(f"\nPhase 5: Grouping...", flush=True)
groups = defaultdict(list)
for i in range(NC):
    groups[canon[i].tobytes()].append(i)

sig_to_ci = {}
classes = []
for ci, (sig, members) in enumerate(groups.items()):
    rep = members[0]
    sc = np.array_equal(canon[rep], comp_canon[rep])
    # SF: check if any member's tilde matches its own class
    scores_v = A_flat[rep].reshape(n,n).sum(axis=1).tolist()
    # GS count
    gs = sum(1 for idx in members
             if all(((cand_masks[idx]>>i)&1)==((cand_masks[idx]>>j)&1) for i,j in gs_pairs))
    classes.append({
        'ci': ci, 'sig': sig, 'rep': rep, 'members': members, 'nt': len(members),
        'sc': sc, 'scores': scores_v, 'gs': gs, 'sf': False,
    })
    sig_to_ci[sig] = ci

# SF check
for i in range(NC):
    ci = sig_to_ci[canon[i].tobytes()]
    t_sig = tilde_canon[i].tobytes()
    if t_sig in sig_to_ci and sig_to_ci[t_sig] == ci:
        classes[ci]['sf'] = True

# |Aut| computation
print(f"\nPhase 6: Computing |Aut|...", flush=True)
A_reps = np.array([A_flat[c['rep']] for c in classes], dtype=np.uint8)
aut_counts = np.zeros(len(classes), dtype=np.int32)
t0 = time.time()
for p_idx in range(NP):
    perm_A = A_reps[:, perm_flat[p_idx]]
    is_aut = (perm_A == A_reps).all(axis=1)
    aut_counts += is_aut.astype(np.int32)
for ci, c in enumerate(classes):
    c['au'] = int(aut_counts[ci])
print(f"  |Aut| done ({time.time()-t0:.1f}s)", flush=True)

# ─── KERNEL ANALYSIS ──────────────────────────────────────────────────────────
sc_cls = [c for c in classes if c['sc']]
sf_cls = [c for c in classes if c['sf']]
scsf_cls = [c for c in classes if c['sc'] and c['sf']]

print(f"\n{'='*60}")
print(f"n=8 KERNEL ANALYSIS")
print(f"{'='*60}")
print(f"Candidates: {NC}")
print(f"Iso classes: {len(classes)}")
print(f"SC:   {len(sc_cls)}")
print(f"SF:   {len(sf_cls)}")
print(f"SC∩SF: {len(scsf_cls)}")

print(f"\nKernel (SC∩SF) classes:")
print(f"{'#':>3} | {'#t':>5} | {'gs':>4} | {'|Aut|':>6} | {'bl=gs':>6} | {'BS?':>5} | scores")
print("-"*70)
for i, c in enumerate(scsf_cls):
    bl = c['gs']
    is_blue = bl > 0
    bs_type = "BLUE" if is_blue else "BLACK"
    scores_s = tuple(sorted(c['scores'], reverse=True))
    print(f"{i+1:>3} | {c['nt']:>5} | {c['gs']:>4} | {c['au']:>6} | {bl:>6} | {bs_type:>5} | {scores_s}")

print(f"\nAll SC classes:")
print(f"{'#':>3} | {'#t':>5} | {'gs':>4} | {'|Aut|':>6} | {'SF?':>4} | scores")
print("-"*60)
for i, c in enumerate(sc_cls):
    sf_mark = "SF" if c['sf'] else "--"
    scores_s = tuple(sorted(c['scores'], reverse=True))
    print(f"{i+1:>3} | {c['nt']:>5} | {c['gs']:>4} | {c['au']:>6} | {sf_mark:>4} | {scores_s}")

# ─── BLUESELF vs BLACKSELF DEEP ANALYSIS ──────────────────────────────────────
print(f"\n{'='*60}")
print(f"BLUESELF/BLACKSELF ANALYSIS")
print(f"{'='*60}")

for i, c in enumerate(scsf_cls):
    print(f"\nKernel class {i+1}: nt={c['nt']}, gs={c['gs']}, |Aut|={c['au']}")
    scores_s = tuple(sorted(c['scores'], reverse=True))
    print(f"  Scores: {scores_s}")

    if c['gs'] > 0:
        print(f"  BLUESELF: has {c['gs']} GS tilings")
        # Find GS tilings in this class and check if they are flip-paired
        gs_in_class = [(idx, cand_masks[idx]) for idx in c['members']
                       if all(((cand_masks[idx]>>pi)&1)==((cand_masks[idx]>>pj)&1)
                              for pi,pj in gs_pairs)]
        print(f"  GS tilings in class: {[hex(mask) for mask in [m for _,m in gs_in_class[:6]]]}")
        # For each GS tiling, find flip pair partner
        flip_all_bits = (1 << m) - 1
        for idx, mask in gs_in_class:
            flip_mask = (~mask) & flip_all_bits
            flip_idx = np.where(cand_masks == flip_mask)[0]
            if len(flip_idx) > 0:
                fi = flip_idx[0]
                same_class = (sig_to_ci.get(canon[fi].tobytes()) == sig_to_ci.get(canon[idx].tobytes()))
                print(f"  GS mask={hex(mask)}: flip={hex(flip_mask)}, same_class={same_class}, flip_is_GS={all(((flip_mask>>pi)&1)==((flip_mask>>pj)&1) for pi,pj in gs_pairs)}")
    else:
        print(f"  BLACKSELF: no GS tilings — SF via non-GS pairs only")
        # Find self-flip pairs in this class
        self_flip_pairs = []
        checked = set()
        flip_all = (1 << m) - 1
        for idx in c['members']:
            if idx in checked: continue
            mask = cand_masks[idx]
            flip_mask = (~mask) & flip_all
            flip_idx = np.where(cand_masks == flip_mask)[0]
            if len(flip_idx) > 0:
                fi = flip_idx[0]
                if fi not in checked:
                    same_class = (sig_to_ci.get(canon[fi].tobytes()) == sig_to_ci.get(canon[idx].tobytes()))
                    if same_class:
                        self_flip_pairs.append((mask, flip_mask))
                        checked.add(idx); checked.add(fi)
        print(f"  Non-GS self-flip pairs: {len(self_flip_pairs)}")
        for mk, fmk in self_flip_pairs[:5]:
            gs1 = all(((mk>>pi)&1)==((mk>>pj)&1) for pi,pj in gs_pairs)
            gs2 = all(((fmk>>pi)&1)==((fmk>>pj)&1) for pi,pj in gs_pairs)
            print(f"    {hex(mk)} (GS={gs1}) ↔ {hex(fmk)} (GS={gs2})")

print(f"\nTotal time: {time.time()-t_global:.1f}s")
