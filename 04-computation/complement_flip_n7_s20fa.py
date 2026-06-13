#!/usr/bin/env python3
"""
complement_flip_n7_s20fa.py — Faster complement-flip overlap at n=7
kind-pasteur-2026-03-24-S20fa

Use score-based pre-filtering and nauty-style canonicalization shortcuts
to speed up the n=7 computation.

Strategy:
1. Build ALL 2^15 = 32768 tilings as adjacency matrices
2. Compute score sequences as a cheap hash
3. Group by score sequence
4. For each tiling: compute complement, check if complement's class is
   a single-flip neighbor (using canonical forms within score groups)
"""

import sys
from math import factorial
from itertools import permutations
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  FAST COMPLEMENT-FLIP OVERLAP AT n=7")
print("  kind-pasteur-2026-03-24-S20fa")
print("=" * 80)

n = 7
TILES = []
for y in range(1, n-1):
    for x in range(n, y+1, -1):
        TILES.append((x, y))
m = len(TILES)
N = n
VERTS = list(range(n, 0, -1))

print(f"  n={n}, m={m}, 2^m={2**m}")

def bits_to_adj(bits):
    A = [[0]*N for _ in range(N)]
    for k in range(N-1): A[k][k+1] = 1
    for i in range(m):
        x, y = TILES[i]
        xi, yi = VERTS.index(x), VERTS.index(y)
        if bits[i] == 0: A[xi][yi] = 1
        else: A[yi][xi] = 1
    return A

def adj_to_tuple(A):
    return tuple(tuple(row) for row in A)

def adj_complement(A):
    return [[1-A[i][j] if i!=j else 0 for j in range(N)] for i in range(N)]

def score_seq(A):
    return tuple(sorted(sum(A[i]) for i in range(N)))

# Pre-compute permutation action on adjacency matrices
all_perms = list(permutations(range(N)))

def canonicalize(A):
    best = None
    for p in all_perms:
        s = tuple(A[p[i]][p[j]] for i in range(N) for j in range(N))
        if best is None or s < best: best = s
    return best

# Step 1: Build all tilings, group by score sequence
print("  Building tilings and score groups...")
t0 = time.time()

score_groups = defaultdict(list)  # score -> [(mask, adj_tuple)]
for mask in range(1 << m):
    bits = [(mask >> k) & 1 for k in range(m)]
    A = bits_to_adj(bits)
    sc = score_seq(A)
    score_groups[sc].append((mask, A))

print(f"  {2**m} tilings in {len(score_groups)} score groups ({time.time()-t0:.1f}s)")

# Step 2: For each score group, canonicalize (this is the expensive part)
print("  Canonicalizing within score groups...")
t1 = time.time()

canon_map = {}  # mask -> canon tuple
canon_to_masks = defaultdict(list)  # canon -> [mask]
done = 0

for sc, members in sorted(score_groups.items(), key=lambda x: -len(x[1])):
    for mask, A in members:
        cn = canonicalize(A)
        canon_map[mask] = cn
        canon_to_masks[cn].append(mask)
    done += len(members)
    if done % 5000 == 0:
        print(f"    {done}/{2**m} ({time.time()-t1:.0f}s)")

print(f"  Canonicalization done ({time.time()-t1:.0f}s)")

# Step 3: Identify classes and complement pairs
unmerged = sorted(set(canon_map.values()))
V = len(unmerged)

# Complement of each class
comp_of = {}
for cn in unmerged:
    mask = canon_to_masks[cn][0]
    bits = [(mask >> k) & 1 for k in range(m)]
    A = bits_to_adj(bits)
    Ac = adj_complement(A)
    comp_of[cn] = canonicalize(Ac)

SC = sum(1 for cn in unmerged if comp_of[cn] == cn)
comp_pairs = set()
for cn in unmerged:
    if comp_of[cn] != cn:
        comp_pairs.add((min(cn, comp_of[cn]), max(cn, comp_of[cn])))

# Merged
merged_map = {cn: min(cn, comp_of[cn]) for cn in unmerged}
V_merged = len(set(merged_map.values()))

print(f"\n  V = {V}, SC = {SC}, V_merged = {V_merged}")
print(f"  Complement pairs: {len(comp_pairs)}")

# Step 4: Find wiggly edges (unmerged)
print("\n  Finding wiggly edges...")
t2 = time.time()

wiggly_edges = defaultdict(set)  # (cn1, cn2) -> set of tile indices
for mask in range(1 << m):
    cn1 = canon_map[mask]
    for wi in range(m):
        cn2 = canon_map[mask ^ (1 << wi)]
        if cn1 < cn2:
            wiggly_edges[(cn1, cn2)].add(wi)
        elif cn2 < cn1:
            wiggly_edges[(cn2, cn1)].add(wi)

E = len(wiggly_edges)
print(f"  Wiggly edges: {E} ({time.time()-t2:.0f}s)")

# Step 5: Find overlap
overlap = set()
for pair in comp_pairs:
    if pair in wiggly_edges:
        overlap.add(pair)

print(f"\n  OVERLAP (complement AND single-flip): {len(overlap)}")
print(f"  Overlap / complement_pairs = {len(overlap)}/{len(comp_pairs)} = {len(overlap)/len(comp_pairs):.4f}")

# Step 6: Characterize overlap pairs
if overlap:
    # Compute H for overlap pairs
    def count_ham_paths(A):
        dp = [[0]*N for _ in range(1 << N)]
        for v in range(N): dp[1 << v][v] = 1
        for msk in range(1, 1 << N):
            for v in range(N):
                if not (msk & (1 << v)): continue
                if dp[msk][v] == 0: continue
                for u in range(N):
                    if msk & (1 << u): continue
                    if A[v][u]: dp[msk | (1 << u)][u] += dp[msk][v]
        return sum(dp[(1 << N) - 1])

    print(f"\n  OVERLAP PAIRS:")
    print(f"  {'H':>6} {'|Aut|':>6} {'Score':>30} {'#tiles':>7}")

    overlap_H = []
    overlap_scores = []
    overlap_tiles = []
    for pair in sorted(overlap):
        cn1, cn2 = pair
        mask1 = canon_to_masks[cn1][0]
        A1 = bits_to_adj([(mask1 >> k) & 1 for k in range(m)])
        H1 = count_ham_paths(A1)
        sc1 = score_seq(A1)
        aut1 = factorial(N) // len(canon_to_masks[cn1])
        ntiles = len(wiggly_edges[pair])
        skips = sorted(set(TILES[wi][0] - TILES[wi][1] - 1 for wi in wiggly_edges[pair]))

        print(f"  {H1:6d} {aut1:6d} {str(sc1):>30} {ntiles:7d}  skips={skips}")
        overlap_H.append(H1)
        overlap_scores.append(sc1)
        overlap_tiles.append(ntiles)

    print(f"\n  H values: {sorted(set(overlap_H))}")
    print(f"  Score sequences: {sorted(set(overlap_scores))}")
    print(f"  Tiles per pair: min={min(overlap_tiles)}, max={max(overlap_tiles)}, avg={sum(overlap_tiles)/len(overlap_tiles):.1f}")
    print(f"  Overlap nodes: {2*len(overlap)}/{V} ({2*len(overlap)/V*100:.1f}%)")

    # All overlap pairs have palindromic scores?
    palindromic = all(sc[i] + sc[N-1-i] == N-1 for sc in overlap_scores for i in range(N//2))
    print(f"  All palindromic scores? {palindromic}")

print(f"\n  Total time: {time.time()-t0:.0f}s")
print("DONE.")
