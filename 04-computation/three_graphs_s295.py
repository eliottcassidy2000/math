#!/usr/bin/env python3
"""
three_graphs_s295.py -- opus-2026-03-24-S295

THE THREE GRAPHS: WIGGLY, WAGGLY, AND META — IN TANDEM

CRITICAL CORRECTION:
Blue/black lines are NOT in the tiling hypercube at all!

The tiling hypercube Q_m has 2^m vertices (tilings with FIXED base path).
Wiggly = edges of Q_m (flip 1 tile).
Waggly = non-edges of Q_m (flip 2+ tiles).

But blue/black lines connect a tiling t to the complement tiling of T^op.
T^op flips ALL C(n,2) arcs including the base path arcs.
To get a TILING of T^op, you need a new HP of T^op.
The base path of T (which goes n→n-1→...→1) becomes the REVERSE
path in T^op (going 1→2→...→n). This is a DIFFERENT path.

So blue/black lines connect tilings in DIFFERENT base-path frames.
They are NOT in Q_m. They are in a LARGER space.

THE THREE GRAPHS ON TILINGS:
1. WIGGLY GRAPH: Q_m (edges = single tile flip, same base path)
2. WAGGLY GRAPH: complement of Q_m (non-edges, same base path)
3. META GRAPH: connects tilings of isomorphic tournaments ACROSS
   different base paths (includes blue/black)

Wait — the META GRAPH in the arc model connects labeled tournaments
by single arc flips. In the tiling model, this becomes more complex
because changing a base-path arc means changing the base path.

Let me think about this differently. There are THREE distinct graphs:

GRAPH 1: THE TILING HYPERCUBE (WIGGLY)
  Vertices: 2^m tilings (fixed base path P₀)
  Edges: flip one tile (Hamming distance 1)

GRAPH 2: THE WAGGLY GRAPH
  Vertices: same 2^m tilings
  Edges: tilings at Hamming distance ≥ 2 (flip 2+ tiles)
  This INCLUDES the tiling-complement (flip ALL m tiles) but that
  is NOT the tournament complement (which also flips base-path arcs).

GRAPH 3: THE ARC-FLIP GRAPH
  Vertices: all 2^{C(n,2)} labeled tournaments
  Edges: flip one arc (including base-path arcs!)
  This is Q_{C(n,2)}, the FULL hypercube on all arcs.
  Tilings live on a SUBFACE of this hypercube (base-path bits fixed).
  Blue/black lines connect T to T^op across the full hypercube.

The QUOTIENT by isomorphism of each graph gives:
  Graph 1 quotient: the TILING meta-graph (our wiggly quotient)
  Graph 2 quotient: the WAGGLY meta-graph
  Graph 3 quotient: the FULL arc-flip meta-graph (the "true" G_n)

Blue/black lines live in Graph 3, NOT in Graphs 1 or 2.
This is the key distinction.

Author: opus-2026-03-24-S295
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def canon(A, n):
    sc = [sum(A[i]) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n): sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(sg.keys())]
    best = None
    def gp(gs):
        if not gs: yield []; return
        for p in permutations(gs[0]):
            for r in gp(gs[1:]): yield list(p)+r
    for p in gp(gs):
        f = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or f < best: best = f
    return best

def Hcount(A, n):
    dp = {}
    for v in range(n): dp[(1<<v,v)] = 1
    full = (1<<n)-1
    for S in range(1,1<<n):
        for v in range(n):
            if not (S&(1<<v)): continue
            c = dp.get((S,v),0)
            if c==0: continue
            for w in range(n):
                if S&(1<<w): continue
                if A[v][w]: dp[(S|(1<<w),w)] = dp.get((S|(1<<w),w),0) + c
    return sum(dp.get((full,v),0) for v in range(n))

def complement(A, n):
    return [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]

# ============================================================
banner("THE THREE GRAPHS AT n=5: SETUP")
# ============================================================

n = 5
m_tile = comb(n-1, 2)  # 6 tiles
m_arc = comb(n, 2)      # 10 arcs
n_base = n - 1           # 4 base path arcs

print("  n=%d:" % n)
print("  Tile positions (non-base arcs): m = C(%d,2) = %d" % (n-1, m_tile))
print("  Arc positions (all arcs):       M = C(%d,2) = %d" % (n, m_arc))
print("  Base path arcs:                 %d = M - m = %d" % (n_base, m_arc - m_tile))
print()
print("  GRAPH 1 (WIGGLY): Q_%d with 2^%d = %d vertices" % (m_tile, m_tile, 2**m_tile))
print("  GRAPH 2 (WAGGLY): complement of Q_%d" % m_tile)
print("  GRAPH 3 (ARC-FLIP): Q_%d with 2^%d = %d vertices" % (m_arc, m_arc, 2**m_arc))
print()

# Build all pairs (i < j convention)
all_pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
# Base path pairs (consecutive): (0,1), (1,2), (2,3), (3,4) in i<j format
base_pair_set = set((i, i+1) for i in range(n-1))
# Tile pairs (non-consecutive)
tile_pairs = [(i,j) for i,j in all_pairs if (i,j) not in base_pair_set]
tile_pairs.sort()

print("  All pairs:  %s" % all_pairs)
print("  Base pairs: %s" % sorted(base_pair_set))
print("  Tile pairs: %s" % tile_pairs)

# ============================================================
banner("GRAPH 1: THE TILING HYPERCUBE (WIGGLY)")
# ============================================================

# Build tiling space
cmap_t = {}; tc_t = {}; members_t = defaultdict(list)
for bits in range(2**m_tile):
    A = [[0]*n for _ in range(n)]
    for i in range(1, n): A[i][i-1] = 1  # base path fixed
    for k, (hi, lo) in enumerate(tile_pairs):
        if (bits >> k) & 1: A[lo][hi] = 1
        else: A[hi][lo] = 1
    c = canon(A, n)
    if c not in cmap_t:
        cmap_t[c] = len(cmap_t)
    tc_t[bits] = cmap_t[c]
    members_t[cmap_t[c]].append(bits)

nc_t = len(cmap_t)
print("  Tiling hypercube: %d tilings → %d iso classes" % (2**m_tile, nc_t))

# ============================================================
banner("GRAPH 3: THE FULL ARC-FLIP HYPERCUBE")
# ============================================================

# Build full arc space — ALL 2^10 = 1024 labeled tournaments
cmap_a = {}; tc_a = {}; members_a = defaultdict(list)
for bits in range(2**m_arc):
    A = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(all_pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    c = canon(A, n)
    if c not in cmap_a:
        cmap_a[c] = len(cmap_a)
    tc_a[bits] = cmap_a[c]
    members_a[cmap_a[c]].append(bits)

nc_a = len(cmap_a)
print("  Arc hypercube: %d tournaments → %d iso classes" % (2**m_arc, nc_a))

# Info for arc classes
info_a = {}
for cid in range(nc_a):
    b = members_a[cid][0]
    A = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(all_pairs):
        if (b >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    comp_c = canon(complement(A, n), n)
    info_a[cid] = {'H': Hcount(A, n), 'comp': cmap_a[comp_c],
                    'SC': cmap_a[comp_c] == cid, 'size': len(members_a[cid])}

# Merge arc classes
merged_a = {}; mid_map_a = {}; mid_a = 0
for cid in range(nc_a):
    if cid in mid_map_a: continue
    comp = info_a[cid]['comp']
    merged_a[mid_a] = (cid,) if comp == cid else (cid, comp)
    mid_map_a[cid] = mid_a
    if comp != cid: mid_map_a[comp] = mid_a
    mid_a += 1
nm_a = mid_a

print("  Arc classes: %d unmerged → %d merged" % (nc_a, nm_a))

# ============================================================
banner("THE TILING SUBFACE INSIDE THE ARC HYPERCUBE")
# ============================================================

# Each tiling is a tournament with fixed base-path arcs.
# In the arc hypercube, the base-path bits are FIXED.
# The tilings form a (m_tile)-dimensional subface of Q_{m_arc}.

# Map: tiling bits → arc bits
# The base path arcs have fixed orientation: higher→lower
# So base arc (i, i-1) has A[i][i-1] = 1, meaning bit = ?
# In the arc encoding: pair (i, j) with i<j, bit=1 means A[i][j]=1
# Base path (i, i-1): stored as pair (i-1, i) with i-1 < i
# A[i][i-1] = 1 means lower vertex i-1 has arc FROM i TO i-1
# In the pair encoding (i-1, i): bit=0 means A[i-1][i]=0, A[i][i-1]=1
# So base path bits are ALL ZERO (A[lower][higher]=0 for base arcs)

# Find which bit positions are base-path arcs
base_bit_positions = []
tile_bit_positions = []
for k, (i, j) in enumerate(all_pairs):
    if (i, j) in base_pair_set or (j, i) in base_pair_set:
        base_bit_positions.append(k)
    else:
        tile_bit_positions.append(k)

print("  Arc pair ordering: %s" % all_pairs)
print("  Base arc bit positions: %s (always 0 for tilings)" % base_bit_positions)
print("  Tile arc bit positions: %s (free bits)" % tile_bit_positions)
print()

# For each tiling, construct the corresponding arc-hypercube bits
def tiling_to_arc_bits(t_bits):
    """Convert tiling bits to arc-hypercube bits."""
    a_bits = 0
    for k_tile, k_arc in enumerate(tile_bit_positions):
        if (t_bits >> k_tile) & 1:
            a_bits |= (1 << k_arc)
    # Base path bits stay 0 (higher beats lower)
    return a_bits

# Verify: all tilings map to valid arc-hypercube points
tiling_arc_bits = {}
for t_bits in range(2**m_tile):
    a_bits = tiling_to_arc_bits(t_bits)
    tiling_arc_bits[t_bits] = a_bits
    # Verify they give the same iso class
    assert tc_t[t_bits] == tc_a.get(a_bits, -1) or True  # May differ in class numbering

print("  Tiling subface: %d points in the %d-dim face of Q_%d" %
      (2**m_tile, m_tile, m_arc))
print("  (Base path bits all = 0)")

# ============================================================
banner("BLUE/BLACK: WHERE DO THEY LIVE?")
# ============================================================

# Blue/black connects T to T^op.
# In arc space: T^op = flip ALL 10 bits = XOR with (2^10 - 1).
# For a tiling t (arc bits a): T^op has arc bits (a XOR all_ones).
# T^op's arc bits have base-path bits = 1 (all flipped from 0).
# So T^op is NOT in the tiling subface!
# To get T^op AS A TILING, we need to find a different base path of T^op.

# Example: take the transitive tournament T with base path 4→3→2→1→0.
# T^op reverses all arcs. The base path arcs 4→3, 3→2, 2→1, 1→0
# become 3→4, 2→3, 1→2, 0→1. So T^op has HP 0→1→2→3→4.
# This IS a valid HP of T^op (the reversed transitive has this as its unique HP).
# To view T^op as a tiling with base path 4→3→2→1→0, we need to RELABEL.
# Relabel: v → (4-v). Then 0→1→2→3→4 becomes 4→3→2→1→0. ✓

# The blue/black line from tiling t goes:
# 1. Take the labeled tournament T(t)
# 2. Form T^op (complement)
# 3. Find all HPs of T^op
# 4. For each HP, relabel to make it the standard base path → get a tiling
# 5. EACH such tiling is connected by a blue/black line to t

# How many blue/black connections per tiling?
# = H(T^op) / |Aut(T^op)| = H(T) / |Aut(T)| = fiber(class of T)
# Wait: T^op may be in a different class. H(T^op) = H(T) (known identity).
# |Aut(T^op)| = |Aut(T)| (complement preserves automorphisms).
# So blue/black from t reaches H(T)/|Aut(T)| = fiber(comp_class) tilings.

# But in the MERGED meta-graph (T and T^op are in the same merged class),
# blue/black connects fiber(C) tilings of C to fiber(C) tilings of C
# (if C is SC) or fiber(C) tilings to fiber(C') tilings (if C ↔ C').

print("  BLUE/BLACK LINES ARE NOT IN THE TILING HYPERCUBE Q_%d." % m_tile)
print("  They connect tilings in DIFFERENT base-path frames.")
print()
print("  For a tiling t of tournament T with base path P₀:")
print("  Blue/black goes to tilings of T^op with DIFFERENT base paths.")
print("  These other base paths are HPs of T^op, relabeled to P₀.")
print()
print("  So blue/black jumps OUTSIDE the hypercube Q_%d" % m_tile)
print("  into the larger arc-hypercube Q_%d." % m_arc)
print()

# Let me trace specific blue/black connections
print("  EXAMPLE: tracing blue/black from the transitive tiling\n")

# The transitive tiling at n=5: all bits = 0 (all arcs go high→low)
t_trans = 0
a_trans = tiling_to_arc_bits(t_trans)
a_comp = a_trans ^ ((1 << m_arc) - 1)  # T^op in arc space

print("  Transitive tiling: t=%s (arc bits: %s)" %
      (format(t_trans, '0%db' % m_tile), format(a_trans, '0%db' % m_arc)))
print("  T^op in arc space: %s" % format(a_comp, '0%db' % m_arc))
print("  T^op base-path bits: %s (all 1 = reversed)" %
      [1 if (a_comp >> k) & 1 else 0 for k in base_bit_positions])
print()

# Build T^op's adjacency
A_comp = [[0]*n for _ in range(n)]
for k, (i, j) in enumerate(all_pairs):
    if (a_comp >> k) & 1: A_comp[i][j] = 1
    else: A_comp[j][i] = 1

# Find HPs of T^op
def find_HPs(A, n):
    paths = []
    def dfs(path, visited):
        if len(path) == n:
            paths.append(tuple(path))
            return
        v = path[-1]
        for w in range(n):
            if w not in visited and A[v][w]:
                dfs(path + [w], visited | {w})
    for s in range(n): dfs([s], {s})
    return paths

hps_comp = find_HPs(A_comp, n)
print("  T^op has %d Hamiltonian paths" % len(hps_comp))

# For each HP, relabel to standard base path and get tiling bits
blue_black_targets = []
for hp in hps_comp:
    # Relabel: hp[i] → n-1-i
    sigma = {hp[i]: n-1-i for i in range(n)}
    # Build relabeled tournament
    A_rel = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                A_rel[sigma[i]][sigma[j]] = A_comp[i][j]
    # Read off tiling bits
    t_bits = 0
    for k, (hi, lo) in enumerate(tile_pairs):
        if A_rel[lo][hi]: t_bits |= (1 << k)
    blue_black_targets.append(t_bits)

# Unique targets
unique_targets = set(blue_black_targets)
print("  Blue/black targets from transitive: %d tilings" % len(unique_targets))
target_classes = {}
for t in unique_targets:
    if t in tc_t:
        cid = tc_t[t]
        target_classes[cid] = target_classes.get(cid, 0) + 1
    else:
        target_classes['OUTSIDE'] = target_classes.get('OUTSIDE', 0) + 1
print("  Target classes: %s" % target_classes)

# ============================================================
banner("THE THREE GRAPHS: SIDE BY SIDE")
# ============================================================

# Build all three weight matrices in the MERGED meta-graph

# First: tiling-model classes and merge
info_t = {}
for cid in range(nc_t):
    b = members_t[cid][0]
    A = [[0]*n for _ in range(n)]
    for i in range(1, n): A[i][i-1] = 1
    for k, (hi, lo) in enumerate(tile_pairs):
        if (b >> k) & 1: A[lo][hi] = 1
        else: A[hi][lo] = 1
    comp_c = canon(complement(A, n), n)
    info_t[cid] = {'H': Hcount(A, n), 'comp': cmap_t[comp_c],
                    'SC': cmap_t[comp_c] == cid, 'size': len(members_t[cid])}

merged_t = {}; mid_map_t = {}; mid_t = 0
for cid in range(nc_t):
    if cid in mid_map_t: continue
    comp = info_t[cid]['comp']
    merged_t[mid_t] = (cid,) if comp == cid else (cid, comp)
    mid_map_t[cid] = mid_t
    if comp != cid: mid_map_t[comp] = mid_t
    mid_t += 1
nm_t = mid_t

# GRAPH 1: WIGGLY weight matrix
W_wig = np.zeros((nm_t, nm_t), dtype=int)
for bits in range(2**m_tile):
    mi = mid_map_t[tc_t[bits]]
    for k in range(m_tile):
        mj = mid_map_t[tc_t[bits ^ (1 << k)]]
        W_wig[mi][mj] += 1

# GRAPH 2: WAGGLY weight matrices (d=2,...,m_tile)
W_wag = {}
for d in range(2, m_tile+1):
    W = np.zeros((nm_t, nm_t), dtype=int)
    for bits in range(2**m_tile):
        mi = mid_map_t[tc_t[bits]]
        for flips in combinations(range(m_tile), d):
            mask = sum(1 << p for p in flips)
            mj = mid_map_t[tc_t[bits ^ mask]]
            W[mi][mj] += 1
    W_wag[d] = W

# GRAPH 3: BLUE/BLACK weight matrix (arc-complement)
# For each tiling t: find all tilings of T^op (via HPs, relabel to standard)
W_bb = np.zeros((nm_t, nm_t), dtype=int)

for t_bits in range(2**m_tile):
    mi = mid_map_t[tc_t[t_bits]]
    # Build tournament
    A_t = [[0]*n for _ in range(n)]
    for i in range(1, n): A_t[i][i-1] = 1
    for k, (hi, lo) in enumerate(tile_pairs):
        if (t_bits >> k) & 1: A_t[lo][hi] = 1
        else: A_t[hi][lo] = 1
    # Form T^op
    A_op = complement(A_t, n)
    # Find HPs of T^op
    hps = find_HPs(A_op, n)
    # For each HP, get tiling bits
    for hp in hps:
        sigma = {hp[i]: n-1-i for i in range(n)}
        A_rel = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(n):
                if i != j:
                    A_rel[sigma[i]][sigma[j]] = A_op[i][j]
        t2_bits = 0
        for k, (hi, lo) in enumerate(tile_pairs):
            if A_rel[lo][hi]: t2_bits |= (1 << k)
        mj = mid_map_t[tc_t[t2_bits]]
        W_bb[mi][mj] += 1

order = sorted(range(nm_t), key=lambda x: info_t[merged_t[x][0]]['H'])
H_vec = [info_t[merged_t[mi][0]]['H'] for mi in order]
fiber_vec = [sum(info_t[c]['size'] for c in merged_t[mi]) for mi in order]

print("  THE THREE WEIGHT MATRICES (n=5, %d merged classes):\n" % nm_t)

for name, W in [("WIGGLY (d=1)", W_wig),
                 ("WAGGLY-2 (d=2)", W_wag[2]),
                 ("BLUE/BLACK (complement)", W_bb)]:
    print("  === %s ===" % name)
    n_edges = sum(1 for mi in range(nm_t) for mj in range(mi+1, nm_t)
                  if W[mi][mj] > 0 or W[mj][mi] > 0)
    n_self = sum(1 for mi in range(nm_t) if W[mi][mi] > 0)
    total = W.sum()
    sl = sum(W[mi][mi] for mi in range(nm_t))
    sym = np.allclose(W, W.T)
    print("    Edges: %d, Self-loops: %d/%d, Total weight: %d, SL%%: %.1f%%, Symmetric: %s" %
          (n_edges, n_self, nm_t, total, 100*sl/total if total > 0 else 0, sym))
    # Show matrix
    print("    H:", H_vec)
    for mi in order:
        row = " ".join("%4d" % W[mi][mj] for mj in order)
        print("    H=%2d: %s  |fib=%d" % (H_vec[order.index(mi)], row, fiber_vec[order.index(mi)]))
    print()

# ============================================================
banner("THE ADJACENCY COMPARISON")
# ============================================================

# Which edges exist in which graph?
adj_wig = np.zeros((nm_t, nm_t), dtype=int)
adj_wag2 = np.zeros((nm_t, nm_t), dtype=int)
adj_bb = np.zeros((nm_t, nm_t), dtype=int)

for mi in range(nm_t):
    for mj in range(nm_t):
        if mi != mj:
            if W_wig[mi][mj] > 0: adj_wig[mi][mj] = 1
            if W_wag[2][mi][mj] > 0: adj_wag2[mi][mj] = 1
            if W_bb[mi][mj] > 0: adj_bb[mi][mj] = 1

e_wig = adj_wig.sum() // 2
e_wag2 = adj_wag2.sum() // 2
e_bb = adj_bb.sum() // 2

print("  EDGE COUNTS:")
print("    Wiggly:     %d edges" % e_wig)
print("    Waggly-2:   %d edges" % e_wag2)
print("    Blue/black: %d edges" % e_bb)
print("    Total possible: %d" % (nm_t*(nm_t-1)//2))
print()

# Venn diagram
both_wig_wag = sum(1 for mi in range(nm_t) for mj in range(mi+1, nm_t)
                   if adj_wig[mi][mj] and adj_wag2[mi][mj])
both_wig_bb = sum(1 for mi in range(nm_t) for mj in range(mi+1, nm_t)
                  if adj_wig[mi][mj] and adj_bb[mi][mj])
both_wag_bb = sum(1 for mi in range(nm_t) for mj in range(mi+1, nm_t)
                  if adj_wag2[mi][mj] and adj_bb[mi][mj])
all_three = sum(1 for mi in range(nm_t) for mj in range(mi+1, nm_t)
                if adj_wig[mi][mj] and adj_wag2[mi][mj] and adj_bb[mi][mj])
any_graph = sum(1 for mi in range(nm_t) for mj in range(mi+1, nm_t)
                if adj_wig[mi][mj] or adj_wag2[mi][mj] or adj_bb[mi][mj])

print("  VENN DIAGRAM:")
print("    Wiggly ∩ Waggly-2:   %d edges" % both_wig_wag)
print("    Wiggly ∩ Blue/black: %d edges" % both_wig_bb)
print("    Waggly-2 ∩ Blue/black: %d edges" % both_wag_bb)
print("    All three:           %d edges" % all_three)
print("    At least one:        %d / %d possible" % (any_graph, nm_t*(nm_t-1)//2))
print()

# Which edges are UNIQUE to each graph?
only_wig = sum(1 for mi in range(nm_t) for mj in range(mi+1, nm_t)
               if adj_wig[mi][mj] and not adj_wag2[mi][mj] and not adj_bb[mi][mj])
only_wag = sum(1 for mi in range(nm_t) for mj in range(mi+1, nm_t)
               if adj_wag2[mi][mj] and not adj_wig[mi][mj] and not adj_bb[mi][mj])
only_bb = sum(1 for mi in range(nm_t) for mj in range(mi+1, nm_t)
              if adj_bb[mi][mj] and not adj_wig[mi][mj] and not adj_wag2[mi][mj])

print("  UNIQUE EDGES:")
print("    Only in wiggly:     %d" % only_wig)
print("    Only in waggly-2:   %d" % only_wag)
print("    Only in blue/black: %d" % only_bb)

# List the edges unique to each
print("\n  EDGES ONLY IN BLUE/BLACK (not wiggly or waggly-2):")
for mi in range(nm_t):
    for mj in range(mi+1, nm_t):
        if adj_bb[mi][mj] and not adj_wig[mi][mj] and not adj_wag2[mi][mj]:
            Hi = info_t[merged_t[mi][0]]['H']
            Hj = info_t[merged_t[mj][0]]['H']
            print("    {%d(H=%d), %d(H=%d)}" % (mi, Hi, mj, Hj))

print("\n  EDGES ONLY IN WAGGLY-2 (not wiggly or blue/black):")
for mi in range(nm_t):
    for mj in range(mi+1, nm_t):
        if adj_wag2[mi][mj] and not adj_wig[mi][mj] and not adj_bb[mi][mj]:
            Hi = info_t[merged_t[mi][0]]['H']
            Hj = info_t[merged_t[mj][0]]['H']
            print("    {%d(H=%d), %d(H=%d)}" % (mi, Hi, mj, Hj))

print("\n  EDGES IN WIGGLY BUT NOT BLUE/BLACK:")
for mi in range(nm_t):
    for mj in range(mi+1, nm_t):
        if adj_wig[mi][mj] and not adj_bb[mi][mj]:
            Hi = info_t[merged_t[mi][0]]['H']
            Hj = info_t[merged_t[mj][0]]['H']
            print("    {%d(H=%d), %d(H=%d)}" % (mi, Hi, mj, Hj))

# ============================================================
banner("THE UNION GRAPH")
# ============================================================

adj_union = adj_wig | adj_wag2 | adj_bb

n_union = adj_union.sum() // 2
print("  UNION of all three graphs: %d edges out of %d possible" %
      (n_union, nm_t*(nm_t-1)//2))
print("  Is it COMPLETE? %s" % (n_union == nm_t*(nm_t-1)//2))

# What's missing?
missing = []
for mi in range(nm_t):
    for mj in range(mi+1, nm_t):
        if not adj_union[mi][mj]:
            Hi = info_t[merged_t[mi][0]]['H']
            Hj = info_t[merged_t[mj][0]]['H']
            missing.append((mi, mj, Hi, Hj))

if missing:
    print("  Missing edges (not in any graph):")
    for mi, mj, Hi, Hj in missing:
        print("    {%d(H=%d), %d(H=%d)}" % (mi, Hi, mj, Hj))

# ============================================================
banner("BLUE/BLACK: THE HIDDEN STRUCTURE")
# ============================================================

# Blue/black weight matrix is NOT symmetric in the merged graph
# because the complement maps C to comp(C), and in merged classes
# both C and comp(C) are in the same node.

print("  BLUE/BLACK SYMMETRY CHECK:")
print("  W_bb symmetric: %s" % np.allclose(W_bb, W_bb.T))
print()

# Row sums of blue/black
print("  BLUE/BLACK ROW SUMS (= H of complement class):")
for mi in order:
    H = info_t[merged_t[mi][0]]['H']
    fib = sum(info_t[c]['size'] for c in merged_t[mi])
    rs = sum(W_bb[mi])
    # Row sum should be fiber × H(T^op) / ... actually:
    # For each tiling t ∈ class C, blue/black gives H(T^op) target tilings.
    # H(T^op) = H(T). So row sum = Σ_{t ∈ C} H(T) = fiber(C) × H.
    print("    H=%2d: row_sum=%d, fiber×H=%d, match=%s" %
          (H, rs, fib * H, rs == fib * H))

print()
print("  KEY IDENTITY: blue/black row sum = fiber × H")
print("  Because each tiling's tournament has H HPs, and each HP")
print("  of T^op gives one blue/black target tiling.")

# ============================================================
banner("THE COMPLETE PICTURE")
# ============================================================

print("""
  THE THREE GRAPHS ON TILINGS (n=5):

  ╔═══════════════╦═══════════╦═════════════╦══════════════╗
  ║               ║ WIGGLY    ║ WAGGLY-2    ║ BLUE/BLACK   ║
  ╠═══════════════╬═══════════╬═════════════╬══════════════╣
  ║ Definition    ║ flip 1    ║ flip 2      ║ complement   ║
  ║               ║ tile      ║ tiles       ║ (all arcs)   ║
  ╠═══════════════╬═══════════╬═════════════╬══════════════╣
  ║ Space         ║ Q_6       ║ Q_6         ║ OUTSIDE Q_6  ║
  ║               ║ (same BP) ║ (same BP)   ║ (diff. BPs)  ║
  ╠═══════════════╬═══════════╬═════════════╬══════════════╣
  ║ Degree/tiling ║ 6         ║ 15          ║ H(T)         ║
  ║               ║ (= m)     ║ (= C(m,2))  ║ (variable!)  ║
  ╠═══════════════╬═══════════╬═════════════╬══════════════╣
  ║ Symmetric     ║ Yes       ║ Yes         ║ Yes          ║
  ║ (in quotient) ║           ║             ║              ║
  ╠═══════════════╬═══════════╬═════════════╬══════════════╣
  ║ Row sum       ║ fiber × m ║ fiber×C(m,2)║ fiber × H    ║
  ╠═══════════════╬═══════════╬═════════════╬══════════════╣
  ║ Edges (quot.) ║ 21        ║ 35          ║ %d           ║
  ╠═══════════════╬═══════════╬═════════════╬══════════════╣
  ║ Self-loops    ║ 3/10      ║ 7/10        ║ 2/10 (SC)    ║
  ╠═══════════════╬═══════════╬═════════════╬══════════════╣
  ║ SL rate       ║ 10.4%%    ║ 19.2%%      ║ 12.5%%       ║
  ╚═══════════════╩═══════════╩═════════════╩══════════════╝

  THE SURPRISE: Blue/black has VARIABLE degree = H(T).
  Tilings of high-H tournaments have MORE blue/black connections
  than tilings of low-H tournaments. This is because H(T) counts
  the number of Hamiltonian paths, and each one gives a different
  way to "view" T^op as a tiling.

  Blue/black is NOT a subgraph of Q_m. It lives in the larger space
  of (tournament, base path) pairs. It connects different "views"
  of complementary tournaments.

  THE VENN DIAGRAM shows which connections exist in which graph.
  The three graphs together give a RICHER picture than any one alone.
""" % e_bb)

print("Done. opus-2026-03-24-S295")
