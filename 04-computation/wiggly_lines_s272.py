#!/usr/bin/env python3
"""
wiggly_lines_s272.py -- opus-2026-03-23-S272

WIGGLY LINES: FIBER CONNECTIONS THROUGH THE OVERLAP

Each tiling (labeled tournament on n vertices) decomposes as:
  T = boundary(0) + boundary(n-1) + apex(0,n-1) + overlap(interior)

The OVERLAP has C(n-2,2) cells (arcs among vertices 1..n-2).
Fixing the boundary+apex and varying the overlap gives 2^{C(n-2,2)}
tilings — a C(n-2,2)-dimensional subcube of Q_m.

WIGGLY LINES connect tilings that share the same boundary+apex
but differ in the overlap. Each tiling has 2^{C(n-2,2)} - 1
wiggly neighbors (all tilings in its overlap fiber minus itself).

QUESTIONS:
1. How do wiggly fibers project onto the meta-graph?
   (Do wiggly-connected tilings land in the SAME class or different?)
2. What is the wiggly analog of the meta-graph edge count?
3. How do wiggly fibers interact with the opaque/translucent structure?

Author: opus-2026-03-23-S272
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations
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

banner("THE WIGGLY LINE CONCEPT")

print("""
  THE RECURSIVE DECOMPOSITION:
  T on [n] = boundary_left (arcs from 0) + boundary_right (arcs to n-1)
             + apex (arc 0↔n-1) + overlap (arcs among 1..n-2)

  CELL COUNTS:
    boundary_left: n-2 cells (arcs 0→j for j=1..n-2)
    boundary_right: n-2 cells (arcs j→(n-1) for j=1..n-2)
    apex: 1 cell (arc 0↔(n-1))
    overlap: C(n-2, 2) cells (arcs among 1..n-2)
    TOTAL: 2(n-2) + 1 + C(n-2,2) = C(n,2) = m  ✓

  WIGGLY FIBER: Fix boundary+apex, vary overlap.
  Fiber size = 2^{C(n-2,2)} (one tournament per interior orientation).

  At n=3: C(1,2)=0 → fiber=1 (no interior, only boundary+apex)
  At n=4: C(2,2)=1 → fiber=2 (one interior arc, 2 choices)
  At n=5: C(3,2)=3 → fiber=8 (three interior arcs, 8 choices)
  At n=6: C(4,2)=6 → fiber=64 (six interior arcs, 64 choices)
  At n=7: C(5,2)=10 → fiber=1024

  Each fiber is a C(n-2,2)-subcube of Q_m.
  The total number of fibers = 2^{2(n-2)+1} = 2^{2n-3}
  (one per boundary+apex configuration).
""")

for n in [4, 5, 6]:
    banner("WIGGLY LINES AT n=%d" % n)
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    pair_idx = {p: k for k, p in enumerate(pairs)}

    # Classify pairs into boundary_left, boundary_right, apex, overlap
    boundary_left = [(0, j) for j in range(1, n-1)]
    boundary_right = [(j, n-1) for j in range(1, n-1)]
    apex_pair = [(0, n-1)]
    overlap_pairs = [(i, j) for i in range(1, n-1) for j in range(i+1, n-1)]

    bl_bits = [pair_idx[p] for p in boundary_left]
    br_bits = [pair_idx[p] for p in boundary_right]
    ap_bits = [pair_idx[p] for p in apex_pair]
    ov_bits = [pair_idx[p] for p in overlap_pairs]

    n_overlap = len(overlap_pairs)
    fiber_size = 2**n_overlap
    n_boundary = 2*(n-2) + 1  # boundary + apex bits
    n_fibers = 2**n_boundary

    print("  Overlap cells: %d = C(%d,2)" % (n_overlap, n-2))
    print("  Boundary+apex cells: %d = 2(%d-2)+1" % (n_boundary, n))
    print("  Fiber size: 2^%d = %d" % (n_overlap, fiber_size))
    print("  Number of fibers: 2^%d = %d" % (n_boundary, n_fibers))
    print("  Total: %d × %d = %d = 2^%d ✓\n" % (n_fibers, fiber_size, n_fibers*fiber_size, m))

    # Build iso classes
    cmap = {}; sizes = {}; tc = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        c = canon(A, n)
        if c not in cmap:
            cid = len(cmap); cmap[c] = cid; sizes[cid] = 0
        sizes[cmap[c]] += 1; tc[bits] = cmap[c]
    nc = len(cmap)

    info = {}
    for cid in range(nc):
        info[cid] = {'size': sizes[cid]}

    # Merged
    for cid in range(nc):
        A = [[0]*n for _ in range(n)]
        bits_rep = None
        for bits in range(2**m):
            if tc[bits] == cid:
                bits_rep = bits; break
        if bits_rep is None: continue
        for k,(ii,jj) in enumerate(pairs):
            if (bits_rep>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        comp_c = canon(complement(A, n), n)
        info[cid]['H'] = Hcount(A, n)
        info[cid]['comp'] = cmap[comp_c]
        info[cid]['SC'] = cmap[comp_c] == cid

    merged = {}; mid_map = {}; mid = 0
    for cid in range(nc):
        if cid in mid_map: continue
        comp = info[cid]['comp']
        merged[mid] = (cid,) if comp == cid else (cid, comp)
        mid_map[cid] = mid
        if comp != cid: mid_map[comp] = mid
        mid += 1
    nm = mid

    # ============================================================
    # WIGGLY FIBER ANALYSIS
    # ============================================================

    # Extract the boundary+apex key for each tournament
    def boundary_key(bits):
        """The boundary+apex configuration (fixing overlap)."""
        key = 0
        for i, k in enumerate(bl_bits + br_bits + ap_bits):
            if (bits >> k) & 1:
                key |= (1 << i)
        return key

    def overlap_key(bits):
        """The overlap configuration."""
        key = 0
        for i, k in enumerate(ov_bits):
            if (bits >> k) & 1:
                key |= (1 << i)
        return key

    # Group tournaments by their boundary key (= wiggly fiber)
    fibers = defaultdict(list)
    for bits in range(2**m):
        bk = boundary_key(bits)
        fibers[bk].append(bits)

    # Verify fiber sizes
    fiber_sizes = Counter(len(f) for f in fibers.values())
    print("  Fiber size distribution: %s" % dict(fiber_sizes))
    print("  All fibers have size %d: %s\n" %
          (fiber_size, all(len(f) == fiber_size for f in fibers.values())))

    # For each fiber: how many distinct merged classes does it visit?
    classes_per_fiber = []
    class_dist_per_fiber = []
    for bk, fiber_tilings in fibers.items():
        classes_in_fiber = set()
        class_counts = Counter()
        for bits in fiber_tilings:
            mi = mid_map[tc[bits]]
            classes_in_fiber.add(mi)
            class_counts[mi] += 1
        classes_per_fiber.append(len(classes_in_fiber))
        class_dist_per_fiber.append(class_counts)

    cpf_dist = Counter(classes_per_fiber)
    print("  CLASSES PER WIGGLY FIBER:")
    for k in sorted(cpf_dist.keys()):
        print("    %d classes: %d fibers" % (k, cpf_dist[k]))

    avg_classes = np.mean(classes_per_fiber)
    print("  Average classes per fiber: %.2f (out of %d merged classes)" %
          (avg_classes, nm))

    # How are the classes distributed within a typical fiber?
    print("\n  EXAMPLE FIBERS:")
    for i, (bk, fiber_tilings) in enumerate(list(fibers.items())[:5]):
        class_counts = Counter()
        for bits in fiber_tilings:
            mi = mid_map[tc[bits]]
            H = info[merged[mi][0]]['H']
            class_counts[(mi, H)] += 1
        print("    Fiber %d (bkey=%d): %s" %
              (i, bk, dict(class_counts)))

    # ============================================================
    # WIGGLY EDGES PROJECTED TO META-GRAPH
    # ============================================================

    banner("WIGGLY EDGES → META-GRAPH EDGES (n=%d)" % n)

    # A "wiggly edge" connects two tilings in the same fiber
    # (same boundary+apex, different overlap).
    # How many distinct META-GRAPH edges are "witnessed" by wiggly connections?

    wiggly_meta_edges = set()
    wiggly_self_loops = 0
    total_wiggly_pairs = 0

    for bk, fiber_tilings in fibers.items():
        # All pairs within this fiber
        for i in range(len(fiber_tilings)):
            for j in range(i+1, len(fiber_tilings)):
                bi = fiber_tilings[i]; bj = fiber_tilings[j]
                mi = mid_map[tc[bi]]; mj = mid_map[tc[bj]]
                total_wiggly_pairs += 1
                if mi == mj:
                    wiggly_self_loops += 1
                else:
                    wiggly_meta_edges.add((min(mi,mj), max(mi,mj)))

    # Compare with actual meta-graph edges
    actual_edges = set()
    for ci in range(nc):
        mi = mid_map[ci]
        for cj in range(nc):
            mj = mid_map[cj]
            if mi < mj:
                # Check if they're adjacent in meta-graph
                for bits in range(2**m):
                    if tc[bits] != ci: continue
                    for k in range(m):
                        if tc[bits ^ (1<<k)] == cj:
                            actual_edges.add((mi, mj))
                            break
                    else: continue
                    break

    E_meta = len(actual_edges)
    E_wiggly = len(wiggly_meta_edges)
    E_both = len(wiggly_meta_edges & actual_edges)

    print("  Total wiggly pairs: %d (= %d fibers × C(%d, 2))" %
          (total_wiggly_pairs, n_fibers, fiber_size))
    print("  Wiggly self-loops: %d (%.1f%%)" %
          (wiggly_self_loops, 100*wiggly_self_loops/total_wiggly_pairs if total_wiggly_pairs > 0 else 0))
    print()
    print("  Wiggly meta-edges (distinct class pairs): %d" % E_wiggly)
    print("  Actual meta-graph edges: %d" % E_meta)
    print("  Wiggly ∩ Actual: %d" % E_both)
    print("  Wiggly \\ Actual: %d (wiggly-only, NOT in meta-graph)" %
          (E_wiggly - E_both))
    print("  Actual \\ Wiggly: %d (opaque-only, NOT from wiggly)" %
          (E_meta - E_both))

    # KEY QUESTION: Do wiggly lines RECOVER the meta-graph edges?
    print("\n  Wiggly recovery rate: %d/%d = %.1f%% of meta-graph edges" %
          (E_both, E_meta, 100*E_both/E_meta if E_meta > 0 else 0))

    # KEY QUESTION: Do wiggly lines ADD new edges?
    if E_wiggly > E_both:
        print("  Wiggly ADDS %d NEW edges not in the meta-graph!" %
              (E_wiggly - E_both))
        # These new edges connect classes reachable by changing the
        # overlap (multiple interior arcs at once), not just one arc.
        print("  (These connect classes reachable via multi-flip in the overlap)")
    else:
        print("  Wiggly adds NO new edges beyond the meta-graph.")

    # ============================================================
    # THE SUM OF SQUARES IDENTITY
    # ============================================================

    banner("SUM OF SQUARES IDENTITY (n=%d)" % n)

    # From S20cv: 2^{m-1} = Σ_O fiber(O)^2
    # where O ranges over (n-2)-classes and fiber(O) = ways to extend
    # an overlap class to a full n-tournament.

    # Let me verify: how do wiggly fibers relate to the sum of squares?
    # Each wiggly fiber has a fixed boundary+apex.
    # The overlap in each fiber is a free (n-2)-tournament.
    # The ISO CLASS of the overlap determines a partition of the fiber.

    # For each fiber: the overlap classes partition the fiber_size tilings.
    # If overlap class O_i has f_i tilings in this fiber:
    # Then Σ f_i = fiber_size and the class distribution determines the wiggly structure.

    # Group by overlap iso class
    def overlap_iso_class(bits):
        """Iso class of the tournament restricted to vertices 1..n-2."""
        n2 = n - 2
        if n2 <= 1: return 0
        A2 = [[0]*n2 for _ in range(n2)]
        for (i,j) in overlap_pairs:
            k = pair_idx[(i,j)]
            if (bits >> k) & 1:
                A2[i-1][j-1] = 1
            else:
                A2[j-1][i-1] = 1
        return canon(A2, n2)

    # For one fiber: class distribution of overlaps
    example_bk = list(fibers.keys())[0]
    example_fiber = fibers[example_bk]
    overlap_classes = Counter()
    for bits in example_fiber:
        oc = overlap_iso_class(bits)
        overlap_classes[oc] += 1

    print("  Example fiber overlap class distribution:")
    for oc, count in overlap_classes.most_common():
        print("    Overlap class: size %d" % count)

    # Sum of squares check
    sum_sq = sum(c**2 for c in overlap_classes.values())
    print("  Σ (count)² = %d, should relate to 2^{m-1} = %d" %
          (sum_sq, 2**(m-1)))

    # Actually the sum of squares identity from S234 says:
    # Over ALL boundary configurations (both left and right),
    # the number of tournaments is:
    # Σ over overlap classes O: (# ways to extend O left) × (# ways to extend O right) × 2
    # = Σ fiber_left(O) × fiber_right(O) × 2

    # For a FIXED boundary, the overlap varies freely:
    # each fiber has EXACTLY one tiling per overlap configuration.
    # So within one fiber: each of the 2^{C(n-2,2)} overlap configs
    # maps to exactly one tiling, and they all have the same boundary.

    print("\n  Within each fiber: %d tilings, %d overlap configs → bijection" %
          (fiber_size, 2**n_overlap))
    print("  (One tiling per overlap, by construction)")

    # The WIGGLY GRAPH within each fiber: it's the C(n-2,2)-cube Q_{C(n-2,2)}
    # because any two overlap configs are connected by a sequence of single-arc flips.
    # But wiggly lines connect ALL pairs in the fiber, not just single-flip neighbors!
    # Wiggly = sharing boundary+apex, ANY difference in overlap (not just 1-flip).

    print("\n  THE WIGGLY GRAPH within a fiber is the COMPLETE GRAPH K_%d" % fiber_size)
    print("  (Every pair of tilings with the same boundary+apex is wiggly-connected)")
    print("  Wiggly degree: %d (= fiber_size - 1 = 2^{C(n-2,2)} - 1)" % (fiber_size - 1))

banner("SYNTHESIS: THE THREE TYPES OF TILING CONNECTIONS")

print("""
  EVERY PAIR OF TILINGS (labeled tournaments) on n vertices has one of:

  1. OPAQUE CONNECTION: differ by one arc, different iso class.
     = meta-graph edges. m neighbors per tiling.
     These form the meta-graph G_n/Z_2.

  2. BLUESELF/BLACKSELF CONNECTION: differ by one arc, SAME iso class.
     = neutral arc flips. Avg (m × arc_neutrality) neighbors per tiling.
     The translucent/self-loop structure from S270.
     The transitive class forms a permutohedron.

  3. WIGGLY CONNECTION: same boundary+apex, different overlap.
     = sharing all boundary arcs but differing in the interior.
     2^{C(n-2,2)} - 1 neighbors per tiling.
     Each wiggly fiber is a COMPLETE GRAPH on the overlap configurations.

  The three types are NOT disjoint:
  - A wiggly pair that differs in exactly 1 overlap arc is ALSO opaque or blueself.
  - A wiggly pair that differs in k overlap arcs is connected by a k-step
    path in Q_m, but directly connected by the wiggly line.

  WIGGLY CONNECTIONS ARE THE RECURSIVE FIBER:
  Fixing the boundary (how vertices 0 and n-1 connect to the interior)
  and the apex (0↔n-1 direction), the interior tournament varies freely.
  The wiggly fiber = the (n-2)-tournament space embedded in n-tournament space.

  THE META-GRAPH PROJECTION OF WIGGLY LINES:
  Some wiggly pairs land in the same class → wiggly self-loops
  Some wiggly pairs land in different classes → wiggly meta-edges

  At n=5: wiggly meta-edges = %d, actual meta-edges = %d
  Wiggly recovery: %.0f%% of meta-graph edges are wiggly-witnessed.

  THE WIGGLY STRUCTURE IS THE RECURSIVE BURNSIDE:
  Each wiggly fiber is parameterized by a (n-2)-tournament.
  The iso class distribution within a fiber depends on HOW the
  interior (n-2)-tournament interacts with the boundary.
  This IS the recursive structure we seek for computing E(G_n).
""" % (E_wiggly, E_meta, 100*E_both/E_meta if E_meta > 0 else 0))

print("Done. opus-2026-03-23-S272")
