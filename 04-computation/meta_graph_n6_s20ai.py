#!/usr/bin/env python3
"""
meta_graph_n6_s20ai.py -- kind-pasteur-2026-03-22-S20ai

META-GRAPH AT n=6: 56 iso classes, weighted edges, full evaluation.

At n=5 we had 12 iso classes. At n=6: 56 (A000568).
This is large enough to see real structure but small enough to compute.

Key questions:
1. The WEIGHTED adjacency matrix W[i,j] = # (tournament, arc flip) transitions
2. Which classes are most connected? Which are isolated?
3. Does the DAG property (orient by H) still hold?
4. H, HC, independence polynomial of G_6
5. How many "chains" (paths from H_min to H_max)?
6. Width of the DAG (max antichain)

For the weight matrix: W[i,j] tells us how many ways to go from
iso class i to iso class j by flipping one arc in one tournament.
This is the TRANSITION INTENSITY between structural types.

Author: kind-pasteur-2026-03-22-S20ai
"""
import sys
import numpy as np
from math import comb
from collections import defaultdict
from itertools import permutations
import time
sys.stdout.reconfigure(line_buffering=True)

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    return sum(dp[((1 << n) - 1, v)] for v in range(n))

print("=" * 70)
print("  META-GRAPH G_6: FULL WEIGHTED EVALUATION")
print("=" * 70)

n = 6
pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
m = len(pairs)  # 15

# ================================================================
# STEP 1: Build iso classes at n=6 via canonical forms
# ================================================================
print(f"\n  Step 1: Building iso classes at n={n} ({2**m} tournaments)...")
t0 = time.time()

# Use a faster canonical form: sort by (score, then adjacency row sums, then full form)
# For n=6, full permutation check (720 perms) is feasible per tournament
# but 32768 * 720 = 23.6M checks. Let's use a hash-based approach.

def canonical_n6(A):
    """Canonical form for n=6 tournament."""
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or form < best:
            best = form
    return best

# Precompute H and canonical form for all tournaments
H_map = {}
class_map = defaultdict(list)  # canonical -> [bits]

for bits in range(2**m):
    A = np.zeros((n,n), dtype=np.int8)
    for k, (i,j) in enumerate(pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    H = count_hp(A, n)
    H_map[bits] = H
    cf = canonical_n6(A)
    class_map[cf].append(bits)

    if (bits + 1) % 2000 == 0:
        elapsed = time.time() - t0
        rate = (bits + 1) / elapsed
        eta = (2**m - bits - 1) / rate
        print(f"    {bits+1}/{2**m} ({rate:.0f}/s, ETA {eta:.0f}s)...")

# Build class list sorted by H
classes = []
for cf, members in sorted(class_map.items(), key=lambda x: (H_map[x[1][0]], len(x[1]))):
    H = H_map[members[0]]
    A = np.zeros((n,n), dtype=np.int8)
    bits0 = members[0]
    for k, (i,j) in enumerate(pairs):
        if (bits0 >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    score = tuple(sorted(A.sum(axis=1).astype(int)))
    classes.append({
        'id': len(classes), 'H': H, 'score': score,
        'size': len(members), 'members': set(members)
    })

N = len(classes)
print(f"\n  {N} iso classes found in {time.time()-t0:.1f}s")

# Build reverse lookup: bits -> class_id
bits_to_class = {}
for c in classes:
    for b in c['members']:
        bits_to_class[b] = c['id']

# ================================================================
# STEP 2: Build weighted adjacency matrix
# ================================================================
print(f"\n  Step 2: Building weighted adjacency matrix...")
t1 = time.time()

W = np.zeros((N, N), dtype=np.int64)

for bits in range(2**m):
    ci = bits_to_class[bits]
    for k in range(m):
        nb = bits ^ (1 << k)
        cj = bits_to_class[nb]
        if ci != cj:
            W[ci][cj] += 1

    if (bits + 1) % 5000 == 0:
        print(f"    {bits+1}/{2**m}...")

print(f"  Weight matrix computed in {time.time()-t1:.1f}s")

# Build unweighted adjacency
adj = (W > 0).astype(np.int8)
n_edges = sum(1 for i in range(N) for j in range(i+1, N) if adj[i][j] or adj[j][i])

# ================================================================
# STEP 3: Basic graph properties
# ================================================================
print(f"\n{'='*70}")
print(f"  G_6 BASIC PROPERTIES")
print(f"{'='*70}\n")

print(f"  Vertices: {N}")
print(f"  Edges (undir): {n_edges}")
print(f"  Density: {2*n_edges/(N*(N-1)):.4f}")

degrees = adj.sum(axis=1) + adj.sum(axis=0)
# For undirected: degree = number of distinct neighbors
deg_undir = np.array([len(set(j for j in range(N) if adj[i][j] or adj[j][i])) for i in range(N)])
print(f"  Degree sequence: min={deg_undir.min()}, max={deg_undir.max()}, mean={deg_undir.mean():.1f}")
print(f"  Degree distribution:")
for d in sorted(set(deg_undir)):
    count = sum(1 for x in deg_undir if x == d)
    print(f"    deg={d}: {count} classes")

# ================================================================
# STEP 4: H-gradient structure (DAG analysis)
# ================================================================
print(f"\n{'='*70}")
print(f"  H-GRADIENT STRUCTURE")
print(f"{'='*70}\n")

n_up = 0; n_down = 0; n_level = 0
level_edges = []
for i in range(N):
    for j in range(i+1, N):
        if adj[i][j] or adj[j][i]:
            hi, hj = classes[i]['H'], classes[j]['H']
            if hi < hj: n_up += 1
            elif hi > hj: n_down += 1
            else:
                n_level += 1
                level_edges.append((i, j, hi))

print(f"  H-increasing edges: {n_up}")
print(f"  H-decreasing edges: {n_down}")
print(f"  H-level edges: {n_level}")
print(f"  DAG? {'YES' if n_down == 0 and n_level == 0 else 'ALMOST' if n_down == 0 else 'NO'}")

if level_edges:
    print(f"\n  Level edges (same H, adjacent classes):")
    for i, j, h in level_edges[:20]:
        print(f"    class {i} -- class {j} at H={h}, sizes=({classes[i]['size']},{classes[j]['size']}), scores=({list(classes[i]['score'])},{list(classes[j]['score'])})")

# ================================================================
# STEP 5: Weight analysis (how many flips between classes?)
# ================================================================
print(f"\n{'='*70}")
print(f"  WEIGHT MATRIX ANALYSIS")
print(f"{'='*70}\n")

# Weight symmetry check
asymmetric = 0
for i in range(N):
    for j in range(i+1, N):
        if W[i][j] != W[j][i]:
            asymmetric += 1

print(f"  Weight symmetric (W[i,j] = W[j,i])? {asymmetric == 0} ({asymmetric} asymmetric pairs)")

# Total weight per class
total_weight = W.sum(axis=1) + np.diag(W)  # include self-loops? No, diag should be 0 for inter-class
total_out = W.sum(axis=1)
for c in classes[:5]:
    expected = c['size'] * m
    actual = int(total_out[c['id']])
    print(f"  Class {c['id']}: size={c['size']}, total_out_weight={actual}, expected(size*m)={expected}, match={actual==expected}")

# Heaviest edges
print(f"\n  Heaviest edges (top 15):")
edge_weights = []
for i in range(N):
    for j in range(i+1, N):
        w = W[i][j] + W[j][i]
        if w > 0:
            edge_weights.append((w, i, j))
edge_weights.sort(reverse=True)
for w, i, j in edge_weights[:15]:
    print(f"    ({i},{j}): weight={w}, H=({classes[i]['H']},{classes[j]['H']}), sizes=({classes[i]['size']},{classes[j]['size']})")

# Lightest edges
print(f"\n  Lightest nonzero edges (bottom 10):")
for w, i, j in sorted(edge_weights)[:10]:
    print(f"    ({i},{j}): weight={w}, H=({classes[i]['H']},{classes[j]['H']})")

# ================================================================
# STEP 6: H-level structure
# ================================================================
print(f"\n{'='*70}")
print(f"  H-LEVEL STRUCTURE")
print(f"{'='*70}\n")

H_levels = defaultdict(list)
for c in classes:
    H_levels[c['H']].append(c['id'])

print(f"  {'H':>4s} {'#classes':>8s} {'sizes':>30s} {'scores':>40s}")
for H_val in sorted(H_levels.keys()):
    cids = H_levels[H_val]
    sizes = [classes[i]['size'] for i in cids]
    scores = list(set(classes[i]['score'] for i in cids))
    print(f"  {H_val:>4d} {len(cids):>8d} {str(sizes)[:30]:>30s} {str(scores)[:40]:>40s}")

# ================================================================
# STEP 7: Chains from min to max
# ================================================================
print(f"\n{'='*70}")
print(f"  CHAINS FROM H_min TO H_max")
print(f"{'='*70}\n")

H_min = min(c['H'] for c in classes)
H_max = max(c['H'] for c in classes)
sources = [c['id'] for c in classes if c['H'] == H_min]
sinks = [c['id'] for c in classes if c['H'] == H_max]

print(f"  H_min = {H_min} ({len(sources)} source classes)")
print(f"  H_max = {H_max} ({len(sinks)} sink classes)")

# Count maximal chains (paths from any source to any sink in the DAG)
# Use DP: for each class, count paths from sources
# (Only if DAG -- if level edges exist, handle them)
chain_count = np.zeros(N, dtype=np.int64)
for s in sources:
    chain_count[s] = 1

# Process in H-order
for H_val in sorted(H_levels.keys()):
    for ci in H_levels[H_val]:
        if chain_count[ci] == 0: continue
        for cj in range(N):
            if adj[ci][cj] and classes[cj]['H'] > classes[ci]['H']:
                chain_count[cj] += chain_count[ci]

total_chains = sum(chain_count[s] for s in sinks)
print(f"  Maximal chains (source->sink): {total_chains}")

# Longest path (number of edges)
# BFS/DFS from sources
longest = 0
dist_from_source = np.full(N, -1, dtype=int)
for s in sources:
    dist_from_source[s] = 0
for H_val in sorted(H_levels.keys()):
    for ci in H_levels[H_val]:
        if dist_from_source[ci] < 0: continue
        for cj in range(N):
            if adj[ci][cj] and classes[cj]['H'] > classes[ci]['H']:
                dist_from_source[cj] = max(dist_from_source[cj], dist_from_source[ci] + 1)

longest = max(dist_from_source)
print(f"  Longest chain: {longest} edges")
print(f"  Number of H-levels: {len(H_levels)}")

# Width (max antichain = max classes at same H)
width = max(len(cids) for cids in H_levels.values())
width_H = [h for h, cids in H_levels.items() if len(cids) == width][0]
print(f"  Width (max antichain): {width} (at H={width_H})")

# ================================================================
# STEP 8: Compare with n=5 results
# ================================================================
print(f"\n{'='*70}")
print(f"  COMPARISON: G_5 vs G_6")
print(f"{'='*70}\n")

print(f"  {'Property':>25s} {'G_5 (n=5)':>15s} {'G_6 (n=6)':>15s} {'Ratio':>10s}")
print(f"  {'-'*25:>25s} {'-'*15:>15s} {'-'*15:>15s} {'-'*10:>10s}")

comparisons = [
    ("Vertices", 12, N, N/12),
    ("Edges", 30, n_edges, n_edges/30),
    ("Density", 0.4545, 2*n_edges/(N*(N-1)), (2*n_edges/(N*(N-1)))/0.4545),
    ("H-levels", 7, len(H_levels), len(H_levels)/7),
    ("Width", 3, width, width/3),
    ("Sources", 1, len(sources), len(sources)),
    ("Sinks", 2, len(sinks), len(sinks)/2),
    ("Chains", 99, total_chains, total_chains/99),
    ("Level edges", 1, n_level, n_level),
    ("Max degree", 7, int(deg_undir.max()), int(deg_undir.max())/7),
]

for name, v5, v6, ratio in comparisons:
    print(f"  {name:>25s} {str(v5):>15s} {str(v6):>15s} {ratio:>10.2f}")

print(f"\n  GROWTH PATTERNS:")
print(f"  Vertices: 12 -> {N} (x{N/12:.1f}) [A000568 growth]")
print(f"  Edges: 30 -> {n_edges} (x{n_edges/30:.1f})")
print(f"  Chains: 99 -> {total_chains} (x{total_chains/99:.1f})")
print(f"  Width: 3 -> {width} (x{width/3:.1f})")
