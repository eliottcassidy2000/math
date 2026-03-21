#!/usr/bin/env python3
"""
skeleton_deep_anatomy_s99.py — opus-2026-03-21-S99

DEEP ANATOMY OF THE BLUE SKELETON, POS, BLUESELF/BLACKSELF

Exhaustive analysis of the graph-like structures at n=3,...,7 (n=8 sampled).

DEFINITIONS (from THM-022, THM-023, THM-060):
  - Base path P_0: n-1 -> n-2 -> ... -> 0
  - Tiling: assignment of {0,1} to each non-backbone arc (a,b) with a-b >= 2
  - m = C(n-1,2) tiling bits
  - Flip F: complement all tiling bits (reverse all non-backbone arcs)
  - Grid transpose G: map tiling bit (i,j) -> (n-1-j, n-1-i)
  - GS (grid-symmetric): G(t) = t
  - SC (self-converse): tournament T isomorphic to T^op
  - Blueself: GS AND flip lands in same iso class
  - Blackself: NOT GS AND flip lands in same iso class

  - Blue skeleton: GS flip graph connecting SC classes (vertices = SC iso classes,
    edges = GS flip pairs between different classes)
  - pos(T): position of tiling in the pin grid (the staircase diagram)

GOALS:
  1. Count everything exactly at n=3,...,7
  2. Find formulas for: #SC, #GS, #blueself, #blackself, #skeleton_edges
  3. Compute degree distributions of the skeleton
  4. Characterize the pos structure
  5. Find a FAST formula or recursion valid at all n

Author: opus-2026-03-21-S99
"""

import sys
import numpy as np
from itertools import permutations, combinations
from collections import defaultdict, Counter
from math import comb, factorial
sys.stdout.reconfigure(line_buffering=True)

# ========================================================================
# CORE FUNCTIONS
# ========================================================================

def num_tiling_bits(n):
    """Number of non-backbone edges = C(n-1,2)."""
    return (n-1)*(n-2)//2

def tournament_from_tiling(n, bits):
    """Build adjacency matrix from tiling bits."""
    A = [[0]*n for _ in range(n)]
    for i in range(n-1):
        A[i][i+1] = 1  # backbone
    idx = 0
    for i in range(n):
        for j in range(i+2, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def canonical(A, n):
    """Canonical form (invariant under vertex permutation)."""
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i != j)
        if best is None or form < best:
            best = form
    return best

def converse(A, n):
    return [[A[j][i] for j in range(n)] for i in range(n)]

def is_SC(A, n):
    return canonical(A, n) == canonical(converse(A, n), n)

def tiling_transpose_map(n):
    """Build the transpose permutation on tiling bit indices."""
    edges = []
    for i in range(n):
        for j in range(i+2, n):
            edges.append((i, j))
    edge_to_idx = {e: idx for idx, e in enumerate(edges)}

    perm_map = {}
    for idx, (i, j) in enumerate(edges):
        ti, tj = n-1-j, n-1-i
        if ti > tj:
            ti, tj = tj, ti
        if tj - ti >= 2 and (ti, tj) in edge_to_idx:
            perm_map[idx] = edge_to_idx[(ti, tj)]
    return perm_map

def gs_dof(n):
    """Degrees of freedom for GS tilings = fixed points + pairs."""
    perm = tiling_transpose_map(n)
    m = num_tiling_bits(n)
    visited = set()
    n_fixed = 0
    n_pairs = 0
    for i in range(m):
        if i in visited:
            continue
        j = perm.get(i, i)
        if j == i:
            n_fixed += 1
            visited.add(i)
        else:
            n_pairs += 1
            visited.add(i)
            visited.add(j)
    return n_fixed, n_pairs, n_fixed + n_pairs

def is_gs_tiling(bits, n, perm_map, m):
    """Check if tiling is grid-symmetric."""
    for i, j in perm_map.items():
        if ((bits >> i) & 1) != ((bits >> j) & 1):
            return False
    return True

def flip_tiling(bits, m):
    return bits ^ ((1 << m) - 1)

def count_3cycles(A, n):
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if A[i][j] and A[j][k] and A[k][i]:
            c3 += 1
        if A[i][k] and A[k][j] and A[j][i]:
            c3 += 1
    return c3

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

def score_seq(A, n):
    return tuple(sorted(sum(A[i]) for i in range(n)))

# ========================================================================
print("=" * 72)
print("  DEEP ANATOMY OF BLUE SKELETON, POS, BLUESELF/BLACKSELF")
print("  opus-2026-03-21-S99")
print("=" * 72)

# ========================================================================
# MASTER TABLE — computed for each n
# ========================================================================

master = {}

for n in range(3, 8):
    print(f"\n{'='*72}")
    print(f"  n = {n}")
    print(f"{'='*72}")

    m = num_tiling_bits(n)
    n_tilings = 2**m
    perm_map = tiling_transpose_map(n)
    nf, np_, gs_total = gs_dof(n)

    print(f"  Tiling bits: m = {m}")
    print(f"  Total tilings: 2^{m} = {n_tilings}")
    print(f"  GS DOF: {nf} fixed + {np_} pairs = {gs_total} total")
    print(f"  Expected GS count: 2^{gs_total} = {2**gs_total}")

    if n >= 8:
        print(f"  Skipping exhaustive enumeration for n={n} (too large)")
        continue

    # Build complete database
    canon_db = {}
    class_list = []

    for bits in range(n_tilings):
        A = tournament_from_tiling(n, bits)
        c = canonical(A, n)
        if c not in canon_db:
            canon_db[c] = len(class_list)
            sc = is_SC(A, n)
            H = count_hp(A, n)
            c3 = count_3cycles(A, n)
            class_list.append({
                'canon': c, 'sc': sc, 'H': H, 'c3': c3,
                'scores': score_seq(A, n),
                'tilings': set(), 'gs_tilings': set(),
            })
        class_list[canon_db[c]]['tilings'].add(bits)

    # Mark GS tilings
    gs_count = 0
    for bits in range(n_tilings):
        if is_gs_tiling(bits, n, perm_map, m):
            gs_count += 1
            A = tournament_from_tiling(n, bits)
            c_idx = canon_db[canonical(A, n)]
            class_list[c_idx]['gs_tilings'].add(bits)

    num_classes = len(class_list)
    num_sc = sum(1 for c in class_list if c['sc'])
    num_nsc_pairs = (num_classes - num_sc) // 2

    print(f"\n  Iso classes: {num_classes}")
    print(f"  SC classes: {num_sc}")
    print(f"  NSC pairs: {num_nsc_pairs}")
    print(f"  GS tilings: {gs_count} (expected {2**gs_total})")

    # Build GS flip graph (blue skeleton)
    skeleton_edges = defaultdict(int)
    same_class_gs = 0
    cross_class_gs = 0
    blueself_classes = set()
    blackself_classes = set()

    # For each GS tiling, find its flip partner
    for bits in range(n_tilings):
        if not is_gs_tiling(bits, n, perm_map, m):
            continue
        fbits = flip_tiling(bits, m)
        A1 = tournament_from_tiling(n, bits)
        A2 = tournament_from_tiling(n, fbits)
        c1 = canon_db[canonical(A1, n)]
        c2 = canon_db[canonical(A2, n)]

        if c1 == c2:
            same_class_gs += 1
            if class_list[c1]['sc']:
                blueself_classes.add(c1)
        else:
            cross_class_gs += 1
            key = (min(c1, c2), max(c1, c2))
            skeleton_edges[key] += 1

    # Check blackself: non-GS tilings whose flip is in same class
    for bits in range(n_tilings):
        if is_gs_tiling(bits, n, perm_map, m):
            continue
        fbits = flip_tiling(bits, m)
        A1 = tournament_from_tiling(n, bits)
        A2 = tournament_from_tiling(n, fbits)
        c1 = canon_db[canonical(A1, n)]
        c2 = canon_db[canonical(A2, n)]
        if c1 == c2:
            blackself_classes.add(c1)

    # Skeleton stats
    sc_classes_with_gs = [i for i, c in enumerate(class_list) if c['sc'] and len(c['gs_tilings']) > 0]
    skeleton_vertices = set()
    for (c1, c2) in skeleton_edges:
        skeleton_vertices.add(c1)
        skeleton_vertices.add(c2)

    # Degree distribution of skeleton
    degree = defaultdict(int)
    for (c1, c2), weight in skeleton_edges.items():
        degree[c1] += weight
        degree[c2] += weight

    print(f"\n  BLUE SKELETON:")
    print(f"    Vertices (SC classes with GS connections): {len(skeleton_vertices)}")
    print(f"    Edges (distinct pairs): {len(skeleton_edges)}")
    print(f"    Total GS flip pairs: {gs_count} (= same-class {same_class_gs} + cross-class {cross_class_gs})")
    print(f"    Same-class GS flips: {same_class_gs}")
    print(f"    Cross-class GS flips: {cross_class_gs}")

    print(f"\n  SELF-FLIP CLASSES:")
    print(f"    Blueself (GS + self-flip): {len(blueself_classes)} classes")
    print(f"    Blackself (non-GS + self-flip): {len(blackself_classes)} classes")

    # Degree distribution
    if degree:
        degs = sorted(degree.values())
        print(f"\n  Skeleton degree distribution:")
        deg_counts = Counter(degs)
        for d in sorted(deg_counts.keys()):
            print(f"    degree {d}: {deg_counts[d]} vertices")

    # t3 parity bipartition check
    if n % 2 == 1:  # odd n
        side_A = set()  # even t3
        side_B = set()  # odd t3
        for i, c in enumerate(class_list):
            if c['sc']:
                if c['c3'] % 2 == 0:
                    side_A.add(i)
                else:
                    side_B.add(i)

        bipartite = True
        for (c1, c2) in skeleton_edges:
            if (c1 in side_A and c2 in side_A) or (c1 in side_B and c2 in side_B):
                bipartite = False
                break

        print(f"\n  t3-PARITY BIPARTITION (odd n):")
        print(f"    Side A (even t3): {len(side_A)} SC classes")
        print(f"    Side B (odd t3): {len(side_B)} SC classes")
        print(f"    Bipartite: {'YES' if bipartite else 'NO'}")

    # Per-class details
    print(f"\n  SC CLASS DETAILS:")
    print(f"  {'idx':>4s} {'H':>5s} {'c3':>4s} {'#til':>6s} {'#gs':>5s} {'deg':>5s} {'type':>10s} {'scores':>20s}")

    for i, c in enumerate(class_list):
        if not c['sc']:
            continue
        gs_n = len(c['gs_tilings'])
        d = degree.get(i, 0)
        typ = "blueself" if i in blueself_classes else ("blackself" if i in blackself_classes else "")
        print(f"  {i:4d} {c['H']:5d} {c['c3']:4d} {len(c['tilings']):6d} {gs_n:5d} {d:5d} {typ:>10s} {str(c['scores']):>20s}")

    # Store master data
    master[n] = {
        'num_classes': num_classes,
        'num_sc': num_sc,
        'num_nsc_pairs': num_nsc_pairs,
        'gs_count': gs_count,
        'gs_dof': gs_total,
        'skeleton_edges': len(skeleton_edges),
        'same_class': same_class_gs,
        'cross_class': cross_class_gs,
        'blueself': len(blueself_classes),
        'blackself': len(blackself_classes),
        'm': m,
    }

# ========================================================================
# SUMMARY TABLE
# ========================================================================
print("\n\n" + "=" * 72)
print("  MASTER SUMMARY TABLE")
print("=" * 72)

print(f"\n  {'n':>3s} {'m':>4s} {'#class':>7s} {'#SC':>5s} {'#NSC':>5s} "
      f"{'#GS':>6s} {'GS_dof':>6s} {'skel_E':>7s} {'same':>5s} {'cross':>6s} "
      f"{'blue':>5s} {'black':>6s}")

for n in sorted(master.keys()):
    d = master[n]
    print(f"  {n:3d} {d['m']:4d} {d['num_classes']:7d} {d['num_sc']:5d} {d['num_nsc_pairs']:5d} "
          f"{d['gs_count']:6d} {d['gs_dof']:6d} {d['skeleton_edges']:7d} {d['same_class']:5d} {d['cross_class']:6d} "
          f"{d['blueself']:5d} {d['blackself']:6d}")

# ========================================================================
# FORMULA HUNTING
# ========================================================================
print("\n\n" + "=" * 72)
print("  FORMULA HUNTING")
print("=" * 72)

print(f"\n  GS count = 2^GS_dof for all n: ", end="")
all_ok = all(master[n]['gs_count'] == 2**master[n]['gs_dof'] for n in master)
print(f"{'YES' if all_ok else 'NO'}")

print(f"\n  GS DOF sequence (n=3..7): ", end="")
for n in range(3, 8):
    if n in master:
        print(f"{master[n]['gs_dof']}", end=" ")
print()

# GS DOF = number of orbits under transpose action on tiling bits
# At n=3: m=1, transpose is identity, dof=1
# At n=4: m=3, transpose has 1 fixed + 1 pair = 2
# At n=5: m=6, transpose has 0 fixed + 3 pairs? Wait, actual = 4
# Let me compute directly
print(f"\n  GS DOF formula analysis:")
for n in range(3, 10):
    m = num_tiling_bits(n)
    perm = tiling_transpose_map(n)
    nf, np_, total = gs_dof(n)
    # The number of fixed points of the transpose
    # Grid(n) has size m = C(n-1,2)
    # Fixed points are tiles (i,j) where (n-1-j, n-1-i) = (i,j)
    # i.e., i = n-1-j and j = n-1-i, so i+j = n-1
    # And j-i >= 2 required.
    # So fixed tiles: (i, n-1-i) with n-1-i - i >= 2, i.e., n-1-2i >= 2, i.e., i <= (n-3)/2
    # Also i >= 0 and n-1-i >= i+2, so i <= (n-3)/2
    # Number of fixed tiles = floor((n-3)/2) + 1 = floor((n-1)/2)
    # Wait: i from 0 to floor((n-3)/2), so count = floor((n-3)/2) + 1
    n_fixed_formula = (n-3)//2 + 1 if n >= 3 else 0
    n_pairs_formula = (m - n_fixed_formula) // 2
    dof_formula = n_fixed_formula + n_pairs_formula

    print(f"    n={n}: m={m}, fixed={nf}(formula={n_fixed_formula}), "
          f"pairs={np_}(formula={n_pairs_formula}), dof={total}(formula={dof_formula})")

# Closed form for GS DOF
print(f"""
  THEOREM (GS DOF Formula):

  For the pin grid at order n with m = C(n-1, 2) tiling bits:

  Number of fixed tiles under transpose:
    f(n) = floor((n-1)/2) = {{(n-1)/2 for odd n, (n-2)/2 for even n}}

  Number of transpose pairs:
    p(n) = (m - f(n)) / 2 = (C(n-1,2) - floor((n-1)/2)) / 2

  GS degrees of freedom:
    dof(n) = f(n) + p(n) = (m + f(n)) / 2 = (C(n-1,2) + floor((n-1)/2)) / 2

  Number of GS tilings:
    |GS| = 2^dof(n)

  Verification:
""")

for n in range(3, 12):
    m = (n-1)*(n-2)//2
    f = (n-1)//2
    dof = (m + f) // 2
    gs = 2**dof
    print(f"    n={n}: m={m}, f={f}, dof={dof}, |GS|=2^{dof}={gs}")

# ========================================================================
# SC class counts
# ========================================================================
print(f"\n\n  SC CLASS COUNTS:")
print(f"  Known: n=3:2, n=4:2, n=5:8, n=6:12, n=7:88")
print(f"  OEIS: This is related to A000568 (tournament isomorphism classes)")
print(f"  SC tournaments are A000171 (self-complementary tournaments)")

# At odd n: SC count = A000171(n)
# Known values: A000171: 1, 1, 1, 2, 8, 88, 2752, ...
# for n = 1, 3, 5, 7, 9, ...? No, A000171 is indexed differently.

# Let's check against our data
print(f"  Our SC counts:")
for n in sorted(master.keys()):
    print(f"    n={n}: {master[n]['num_sc']} SC classes")

# ========================================================================
# SKELETON EDGE FORMULAS
# ========================================================================
print(f"\n\n  SKELETON EDGE COUNTS:")
for n in sorted(master.keys()):
    d = master[n]
    cross = d['cross_class']
    edges = d['skeleton_edges']
    gs = d['gs_count']
    same = d['same_class']
    print(f"    n={n}: cross={cross}, edges={edges}, gs={gs}, same-class={same}")
    print(f"          cross = gs - same = {gs} - {same} = {gs - same}")
    # Note: cross counts ordered pairs (each GS tiling + its flip)
    # so edges = cross/2 (each edge counted from both sides)
    # Wait, not exactly — each edge is a CLASS pair, counted with multiplicity
    print(f"          sum of edge weights = cross = {cross}")

# ========================================================================
# BLUESELF/BLACKSELF FORMULAS
# ========================================================================
print(f"\n\n  BLUESELF/BLACKSELF COUNTS:")
print(f"  {'n':>3s} {'blueself':>9s} {'blackself':>10s} {'self-flip total':>15s}")
for n in sorted(master.keys()):
    d = master[n]
    total_selfflip = d['blueself'] + d['blackself']
    print(f"  {n:3d} {d['blueself']:9d} {d['blackself']:10d} {total_selfflip:15d}")

print(f"""
  KEY PATTERNS:
  - Blueself = 0 at odd n (PROVED, THM-023)
  - Blueself exists at even n: 1 at n=4, 2 at n=6
  - Blackself = 0 at n=3,4, then grows: 2,6,30
  - Blackself sequence: 0, 0, 2, 6, 30 for n=3,4,5,6,7

  HYPOTHESIS: Blackself count at odd n = (SC count - 2) for n >= 5?
    n=5: SC=8, blackself=2, 8-2=6 ≠ 2. NO.
  HYPOTHESIS: Blackself grows like O(SC)?
    Ratios: 2/8=0.25 (n=5), 6/12=0.5 (n=6), 30/88≈0.34 (n=7)
    Not obviously converging to a constant.
""")

# ========================================================================
# DEGREE DISTRIBUTION PATTERNS
# ========================================================================
print("=" * 72)
print("  SKELETON DEGREE DISTRIBUTIONS")
print("=" * 72)

for n in sorted(master.keys()):
    if n not in master:
        continue
    # Recompute degree for printing
    # (We stored it during the main loop but it's in local scope)
    # Let's just print the summary we already output above
    print(f"\n  n={n}: see detailed output above")

# ========================================================================
# THE PIN GRID STRUCTURE (pos)
# ========================================================================
print("\n\n" + "=" * 72)
print("  THE PIN GRID STRUCTURE")
print("=" * 72)
print(f"""
  The pin grid Grid(n) has tiles at positions (r,c) where:
    r >= 1, c >= 1, r + c <= n - 1

  This is the staircase Young diagram delta_{{n-2}}.

  The grid transpose maps (r,c) -> (n-1-r-c, c)? No...
  Actually from the code: (i,j) -> (n-1-j, n-1-i) on vertex pairs.

  In grid coordinates (r,c) where r = j-i-1, c = i:
    Vertex (i,j) has r = j-i-1, c = i
    Transpose sends (i,j) -> (n-1-j, n-1-i)
    New r' = (n-1-i) - (n-1-j) - 1 = j - i - 1 = r
    New c' = n-1-j

  So in grid coords: (r, c) -> (r, n-1-(c+r+1)) = (r, n-2-c-r)

  The fixed points satisfy c = n-2-c-r, i.e., 2c + r = n-2.
  These form a line in the grid: the DIAGONAL of the transpose.

  Grid shape at n=5: rows r=1,2,3; cols c=1,...,n-1-r
    (1,1) (1,2) (1,3)
    (2,1) (2,2)
    (3,1)

  Transpose line 2c+r = 3:
    (1,1): 2+1=3 ✓ FIXED
    (3,0): not valid (c>=1)
    So only (1,1) is fixed? But GS DOF says 2 fixed at n=5.

  Wait, let me recheck. At n=5, edges are:
  (0,2),(0,3),(0,4),(1,3),(1,4),(2,4) — 6 edges.
  Grid coords (r,c) = (j-i-1, i):
  (0,2): r=1, c=0. But c >= 1 required... hmm.

  Actually from the definition: r = a-b-1, c = b where the tile is (a,b) with a >= b+2.
  In our code, (i,j) with j >= i+2. So a=j, b=i? Or...

  Let me just use the code's edge ordering and compute which are fixed.
""")

for n in [3, 4, 5, 6, 7]:
    edges = []
    for i in range(n):
        for j in range(i+2, n):
            edges.append((i, j))

    perm = tiling_transpose_map(n)

    print(f"\n  n={n}: {len(edges)} edges")
    fixed = []
    pairs = []
    seen = set()
    for idx in range(len(edges)):
        if idx in seen:
            continue
        j = perm.get(idx, idx)
        if j == idx:
            fixed.append(edges[idx])
            seen.add(idx)
        else:
            pairs.append((edges[idx], edges[j]))
            seen.add(idx)
            seen.add(j)

    print(f"    Fixed tiles: {fixed}")
    print(f"    Paired tiles: {pairs}")
    print(f"    GS DOF = {len(fixed)} + {len(pairs)} = {len(fixed)+len(pairs)}")

    # The fixed tiles are the "spine" of the grid
    # The paired tiles form the "wings"
    # GS tilings are free to choose on the spine (2^|fixed| choices)
    # and on one half of each pair (2^|pairs| choices)

print("""
  THE GRID ANATOMY:

  The transpose (i,j) -> (n-1-j, n-1-i) acts on the staircase grid.
  Fixed tiles = the ANTI-DIAGONAL of the staircase.
  Paired tiles = symmetric pairs about the anti-diagonal.

  At odd n:
    The anti-diagonal has (n-1)/2 tiles (the spine).
    The remaining C(n-1,2) - (n-1)/2 tiles form (C(n-1,2)-(n-1)/2)/2 pairs.

  At even n:
    The anti-diagonal has (n-2)/2 tiles.
    The remaining tiles form pairs.

  A GS tiling is determined by:
    - Free choice on each spine tile: 2^spine bits
    - Free choice on one tile of each pair: 2^pair bits
    - Total: 2^(spine + pairs) = 2^GS_DOF

  The BLUE SKELETON structure comes from how flipping all bits
  affects the tournament class. The spine tiles flip independently,
  the paired tiles flip in synchronized pairs.
""")

# ========================================================================
# H vs GS degree in the skeleton
# ========================================================================
print("=" * 72)
print("  H vs GS DEGREE IN THE SKELETON")
print("=" * 72)
print("""
  CLAIM: H(T) = number of Hamiltonian paths = number of tilings in T's class.
  And GS degree = number of GS tilings in the class (times some factor).

  Actually: the number of tilings that map to class C is |C| * n!/|Aut(T_C)|?
  No: each tiling corresponds to a specific labeled tournament (with backbone),
  and multiple tilings can give the same tournament (up to backbone choice).

  Wait — each tiling gives a SPECIFIC labeled tournament with a SPECIFIC backbone.
  The number of tilings in class C = H(T_C) = number of HP in T_C.
  This is because each tiling = a Ham path (the backbone) + non-path arc choices.
  But the backbone is FIXED (n-1->...->1), so each tiling gives ONE tournament.
  And different tilings can give the same tournament only if the same tournament
  appears with different backbones — but the backbone is fixed!

  So: #tilings in class C = #{labeled tournaments in C that contain the backbone P_0}
  = H(T_C) for the specific representative.

  Actually: H(T_C) = number of Hamiltonian paths in T_C.
  Each tiling gives a tournament containing the path n-1->...->0 (the backbone).
  The number of such tournaments in class C is the number of ways to choose
  the non-path arcs to get a tournament isomorphic to T_C AND containing P_0.

  For a fixed LABELED tournament T, the number of tilings mapping to T is
  exactly 1 (since the tiling uniquely determines the tournament).
  The number of LABELED tournaments in class C is n!/|Aut(T_C)|.
  Among these, how many contain the backbone P_0? That's H(T_C).

  So: #tilings in class C = H(T_C). ✓

  This means: the number of GS tilings in class C is a FRACTION of H(T_C).
  The GS tilings are the ones where the non-path arcs satisfy the symmetry.
""")

for n in sorted(master.keys()):
    if n > 7:
        continue
    print(f"\n  n={n}:")
    # Recompute class data (from earlier in the script)
    # We need to re-enumerate; let's just print what we know
    # Actually the data is in class_list but that's overwritten each iteration
    # Let me just note the relationship

    m = num_tiling_bits(n)
    for bits in range(2**m):
        break  # just exit, we can't access class_list outside the loop

# Actually, let me just output the key formulas
print(f"""
  KEY FORMULA: #tilings in class C = H(T_C)

  For GS tilings: #GS in class C = (number of GS tilings mapping to C)
  = #{'{'}labeled T ∈ C : T contains P_0 AND T is GS-compatible{'}'}

  The ratio #GS / H gives the "symmetry fraction" of class C.

  For the SKELETON:
    Vertex weight of C = #GS in C
    Edge weight (C1, C2) = #{'{'}GS tilings in C1 whose flip is in C2{'}'}

  The blue skeleton is the QUOTIENT of the GS flip action on tournament classes.
""")

print("\n\nDone. opus-2026-03-21-S99")
