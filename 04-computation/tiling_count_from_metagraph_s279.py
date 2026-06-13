#!/usr/bin/env python3
"""
tiling_count_from_metagraph_s279.py -- opus-2026-03-24-S279

FORCING THE TILING COUNT FROM META-GRAPH STRUCTURE

Each merged iso class has:
  - fiber = # labeled tilings = n!/|Aut| (the target we want to derive)
  - blue_count = # grid-symmetric tilings in the class
  - black_count = fiber - blue_count = # non-grid-symmetric tilings
  - wiggly_weight_to[D] = # labeled (T, tile) pairs in this class
    whose single flip lands in class D
  - wiggly_self_weight = # labeled (T, tile) pairs that stay in this class
  - total_wiggly = fiber × m (m wiggly lines per tiling)

THE KEY IDENTITY:
  total_wiggly_out + wiggly_self = fiber × m
  where total_wiggly_out = Σ_D wiggly_weight_to[D] for D ≠ self

So: fiber = (total_wiggly_out + wiggly_self) / m

But total_wiggly_out = Σ_D W[self][D] for D ≠ self in the F-matrix.
And wiggly_self = W[self][self] = self-loop weight.

So: fiber = Σ_D W[self][D] / m = (row sum of F-matrix) / m

Wait — the F-matrix row sum = size(C) × m (each labeled tournament
has m neighbors). So: fiber = F_row_sum / m = size(C) = fiber. Circular!

THE NON-TRIVIAL APPROACH: Use the RATIO of blue/black counts to
wiggly weights to determine |Aut| independently.

Author: opus-2026-03-24-S279
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

def is_grid_sym(A, n):
    """Is the tournament self-complementary under the anti-diagonal reflection?"""
    # Grid-symmetric: A[i][j] = 1-A[n-1-j][n-1-i] for all i,j
    for i in range(n):
        for j in range(n):
            if i == j: continue
            if A[i][j] != 1 - A[n-1-j][n-1-i]:
                return False
    return True

for n in [4, 5]:
    banner("DERIVING TILINGS FROM META-GRAPH STRUCTURE AT n=%d" % n)
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

    cmap = {}; sizes = {}; tc = {}; reps = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        c = canon(A, n)
        if c not in cmap:
            cid = len(cmap); cmap[c] = cid; sizes[cid] = 0
            reps[cid] = [row[:] for row in A]
        sizes[cmap[c]] += 1; tc[bits] = cmap[c]
    nc = len(cmap)

    info = {}
    for cid in range(nc):
        A = reps[cid]
        comp_c = canon(complement(A, n), n)
        info[cid] = {'H': Hcount(A, n), 'comp': cmap[comp_c],
                     'SC': cmap[comp_c] == cid, 'size': sizes[cid],
                     'aut': factorial(n)//sizes[cid]}

    merged = {}; mid_map = {}; mid = 0
    for cid in range(nc):
        if cid in mid_map: continue
        comp = info[cid]['comp']
        merged[mid] = (cid,) if comp == cid else (cid, comp)
        mid_map[cid] = mid
        if comp != cid: mid_map[comp] = mid
        mid += 1
    nm = mid

    # Compute blue/black counts per merged node
    blue_counts = {}
    for mi in range(nm):
        cids = merged[mi]
        blue = 0
        for cid in cids:
            for bits in range(2**m):
                if tc[bits] != cid: continue
                A = [[0]*n for _ in range(n)]
                for k,(ii,jj) in enumerate(pairs):
                    if (bits>>k)&1: A[ii][jj]=1
                    else: A[jj][ii]=1
                if is_grid_sym(A, n):
                    blue += 1
        blue_counts[mi] = blue

    # Wiggly weights
    W = np.zeros((nm, nm), dtype=int)
    for bits in range(2**m):
        mi = mid_map[tc[bits]]
        for k in range(m):
            mj = mid_map[tc[bits ^ (1 << k)]]
            W[mi][mj] += 1

    # Fiber sizes
    fiber = {}
    for mi in range(nm):
        fiber[mi] = sum(sizes[cid] for cid in merged[mi])

    H_vec = {mi: info[merged[mi][0]]['H'] for mi in range(nm)}
    aut_vec = {mi: info[merged[mi][0]]['aut'] for mi in range(nm)}

    print("  %-4s | %-4s | %-4s | %-6s | %-6s | %-8s | %-8s | %-8s | %-10s" %
          ("mid", "H", "|Aut|", "fiber", "blue", "SL_wt", "cross_wt", "total_wt", "fiber_from"))
    print("  " + "-"*85)

    for mi in range(nm):
        H = H_vec[mi]; aut = aut_vec[mi]; fib = fiber[mi]
        blue = blue_counts[mi]; black = fib - blue
        sl_wt = W[mi][mi]
        cross_wt = sum(W[mi][mj] for mj in range(nm) if mj != mi)
        total_wt = sl_wt + cross_wt

        # DERIVATION 1: fiber from total wiggly weight
        fiber_from_wt = total_wt // m

        # DERIVATION 2: |Aut| from blue fraction
        # blue_fraction = blue / fiber = 2^{-(n-2)} for grid-sym
        # Actually: blue = #{grid-sym tilings in class}
        # #{grid-sym tilings} = n! / |Aut_gs| where Aut_gs = automorphisms
        # preserving grid symmetry... this is complex.

        # DERIVATION 3: fiber from H and the tiling formula
        # H(T) × |Aut(T)| = (tilings × |Aut|) ... no, that's the
        # orbit-stabilizer: size(C) × |Aut(C)| = n!.
        # So fiber = n! / |Aut|.
        # Can we get |Aut| from the meta-graph?

        # DERIVATION 4: |Aut| from self-loop pattern
        # The self-loop weight = #{(T, arc) : T in class, flip arc stays in class}
        # = fiber × (avg neutral arcs per T in this class)
        # SL_wt / fiber = avg neutral arcs.
        # If |Aut| = 1: every T has the same neutral arc pattern.
        # If |Aut| > 1: automorphisms create symmetry in the neutral pattern.

        # DERIVATION 5: base_weight encodes fiber information
        # From the perfect factorization: base_weight(edge) = total_weight / m
        # And total_weight for edge {C,D} = #{(T, arc) : T∈C, T^arc∈D or T∈D, T^arc∈C}
        # The total_weight depends on BOTH fiber(C) and fiber(D).

        print("  %-4d | %-4d | %-4d | %-6d | %-6d | %-8d | %-8d | %-8d | %-10d" %
              (mi, H, aut, fib, blue, sl_wt, cross_wt, total_wt, fiber_from_wt))

    # ============================================================
    banner("THE CONSTRAINT SYSTEM AT n=%d" % n)
    # ============================================================

    # The wiggly weights give us a LINEAR SYSTEM:
    # For each edge {mi, mj}: W[mi][mj] + W[mj][mi] = total_weight(edge)
    # The total weight depends on fiber(mi) and fiber(mj).

    # From the perfect factorization:
    # W[mi][mj] = base(mi,mj) × m (with base = fiber-dependent)

    # Actually: W[mi][mj] = Σ over T in merged[mi]: #{arcs k : T^k in merged[mj]}
    # Each T has m arcs. By S_n transitivity, each arc position contributes
    # equally at the CLASS level. So W[mi][mj] = fiber(mi) × (avg arcs to mj per T).

    # Define: p(mi → mj) = probability that a random (T in mi, random arc)
    # lands in class mj. Then W[mi][mj] = fiber(mi) × m × p(mi → mj).

    # And: Σ_mj p(mi → mj) = 1 (probabilities sum to 1).
    # So: Σ_mj W[mi][mj] = fiber(mi) × m. CHECK.

    # The transition probabilities p(mi → mj) form a STOCHASTIC MATRIX.
    # This is the MARKOV CHAIN on the merged meta-graph!

    print("  TRANSITION PROBABILITIES p(mi → mj) = W[mi][mj] / (fiber(mi) × m):\n")
    P = np.zeros((nm, nm))
    for mi in range(nm):
        row_sum = sum(W[mi][mj] for mj in range(nm))
        for mj in range(nm):
            P[mi][mj] = W[mi][mj] / row_sum if row_sum > 0 else 0

    # The Markov chain P has stationary distribution proportional to fiber.
    # Because: fiber(mi) × P[mi][mj] = W[mi][mj] / m
    # And: fiber(mj) × P[mj][mi] = W[mj][mi] / m
    # Detailed balance: W[mi][mj] = W[mj][mi]? Not necessarily.
    # But: Σ_mi fiber(mi) × P[mi][mj] = (1/m) × Σ_mi W[mi][mj]
    #     = (1/m) × (column sum of W at mj) = ... not obviously fiber(mj).

    # Check: is the Markov chain reversible with stationary dist ∝ fiber?
    pi = np.array([fiber[mi] for mi in range(nm)], dtype=float)
    pi /= pi.sum()

    print("  Stationary distribution (∝ fiber):")
    for mi in range(nm):
        print("    π(%d) = %.4f (fiber=%d)" % (mi, pi[mi], fiber[mi]))

    # Check detailed balance: π(i) P(i,j) = π(j) P(j,i)
    db_ok = True
    for mi in range(nm):
        for mj in range(mi+1, nm):
            lhs = pi[mi] * P[mi][mj]
            rhs = pi[mj] * P[mj][mi]
            if abs(lhs - rhs) > 1e-10:
                db_ok = False
                print("  DB FAIL: π(%d)P(%d,%d)=%.6f ≠ π(%d)P(%d,%d)=%.6f" %
                      (mi, mi, mj, lhs, mj, mj, mi, rhs))

    if db_ok:
        print("\n  ✓ DETAILED BALANCE VERIFIED: π(i)P(i,j) = π(j)P(j,i)")
        print("    The Markov chain is REVERSIBLE with π ∝ fiber.")
    else:
        print("\n  ✗ Detailed balance FAILS.")

    # ============================================================
    banner("FORCING FIBER FROM BLUE + WIGGLY (n=%d)" % n)
    # ============================================================

    # APPROACH: The total wiggly weight row sum gives fiber × m.
    # So fiber = row_sum / m. This is trivial BUT requires the F-matrix.

    # Can we get fiber from the UNWEIGHTED meta-graph + blue counts?

    # We know:
    # (a) Each node has degree d(mi) in the meta-graph.
    # (b) Each node has blue_count(mi) blue tilings.
    # (c) Total tilings = 2^m.
    # (d) fiber(mi) = #{tilings in merged class mi}.
    # (e) Σ fiber(mi) = 2^m.
    # (f) fiber(mi) = n!/|Aut| or 2×n!/|Aut| depending on SC/non-SC.

    # The blue fraction of tilings: blue / 2^m = 2^{-(n-2)} (grid-sym fraction).
    # So Σ blue(mi) = 2^m × 2^{-(n-2)} = 2^{m-(n-2)} = 2^{C(n-1,2)}.

    total_blue = sum(blue_counts.values())
    expected_blue = 2**(m - (n-2))
    print("  Total blue tilings: %d (expected 2^{m-(n-2)} = %d, match=%s)" %
          (total_blue, expected_blue, total_blue == expected_blue))

    # CONSTRAINT from wiggly base weights:
    # base_weight(mi, mj) = total_weight(mi, mj) / m
    # And from the perfect factorization: base_weight is constant per edge.

    # The base weight relates to fibers:
    # total_weight(mi, mj) = Σ_{T in mi} #{arcs: T^arc ∈ mj} + Σ_{T in mj} #{arcs: T^arc ∈ mi}
    # = fiber(mi) × m × p(mi→mj) + fiber(mj) × m × p(mj→mi)

    # With detailed balance: fiber(mi) × p(mi→mj) = fiber(mj) × p(mj→mi)
    # So: total_weight = 2 × fiber(mi) × m × p(mi→mj)

    # And: base_weight = total_weight / m = 2 × fiber(mi) × p(mi→mj)

    # The self-loop: base_self(mi) = 2 × fiber(mi) × p(mi→mi) / m
    # Hmm, self-loop weight = fiber(mi) × m × p(mi→mi) (not doubled).

    print("\n  CAN WE DERIVE FIBER FROM THE TRANSITION MATRIX?")
    print("  The transition matrix P(mi, mj) = W[mi][mj] / (fiber(mi) × m)")
    print("  If we know P and W, we can solve: fiber(mi) = W_row_sum(mi) / m\n")

    # The point: P is what we'd measure from the meta-graph dynamics.
    # W is the RAW weight. fiber = W_row / m.

    # But measuring P requires knowing fiber (to normalize W).
    # CHICKEN-AND-EGG: you need fiber to get P, and P to get fiber.

    # HOWEVER: the UNWEIGHTED meta-graph gives us the ADJACENCY.
    # And the blue counts give us a PARTIAL fiber estimate.

    # KEY INSIGHT: For the transitive class (|Aut|=1):
    # fiber = n! (known exactly). Blue count = depends on grid-sym.
    # For a class with |Aut|=a: fiber = n!/a.

    # The RATIO of fibers between connected classes:
    # fiber(mi)/fiber(mj) = |Aut(mj)|/|Aut(mi)| = aut_ratio.
    # From the transition matrix: P(mi→mj)/P(mj→mi) = fiber(mj)/fiber(mi) = 1/aut_ratio.
    # And from W: W[mi][mj]/W[mj][mi] = fiber(mi)×P(mi→mj) / (fiber(mj)×P(mj→mi))
    #           = (fiber(mi)/fiber(mj)) × (P(mi→mj)/P(mj→mi))
    #           = aut_ratio × (1/aut_ratio) = 1.
    # So W[mi][mj] = W[mj][mi]!!! The weight matrix is SYMMETRIC.

    print("  IS W SYMMETRIC? W[mi][mj] = W[mj][mi]?")
    w_sym = True
    for mi in range(nm):
        for mj in range(mi+1, nm):
            if W[mi][mj] != W[mj][mi]:
                w_sym = False
                print("  FAIL: W[%d][%d]=%d ≠ W[%d][%d]=%d" %
                      (mi, mj, W[mi][mj], mj, mi, W[mj][mi]))
    if w_sym:
        print("  ✓ W IS SYMMETRIC: W[mi][mj] = W[mj][mi] for all pairs")

    # If W is symmetric AND detailed balance holds:
    # fiber(mi) × P(mi→mj) = fiber(mj) × P(mj→mi)
    # fiber(mi) × W[mi][mj]/(fiber(mi)×m) = fiber(mj) × W[mj][mi]/(fiber(mj)×m)
    # W[mi][mj]/m = W[mj][mi]/m. ✓ (automatically from W symmetric.)

    # So detailed balance holds trivially when W is symmetric.
    # But this means we CAN'T distinguish fiber ratios from W alone!
    # Because W[mi][mj] = W[mj][mi] regardless of the fiber ratio.

    # THE FIBER INFORMATION IS IN THE DIAGONAL:
    # W[mi][mi] = fiber(mi) × m × p(mi→mi) = fiber(mi) × (neutral arcs per T in mi)
    # And: W_row_sum(mi) = fiber(mi) × m.
    # So: p(mi→mi) = W[mi][mi] / (fiber(mi) × m) = W[mi][mi] / W_row_sum(mi).

    # The self-loop PROBABILITY is what we can measure without knowing fiber.
    # And: fiber(mi) = W_row_sum(mi) / m.

    # So: the ONLY way to get fiber from the meta-graph is via W_row_sum / m.
    # This requires the FULL WEIGHT information, not just the topology.

    print("\n  FIBER DERIVATION:")
    for mi in range(nm):
        row_sum = sum(W[mi][mj] for mj in range(nm))
        fiber_derived = row_sum // m
        fiber_actual = fiber[mi]
        H = H_vec[mi]; aut = aut_vec[mi]
        print("    Node %d (H=%d, |Aut|=%d): W_row=%d, fiber=W_row/m=%d, actual=%d, match=%s" %
              (mi, H, aut, row_sum, fiber_derived, fiber_actual, fiber_derived == fiber_actual))

    # ============================================================
    banner("THE BLUE-WIGGLY CONSTRAINT (n=%d)" % n)
    # ============================================================

    # What if we DON'T have the full weights, but only:
    # (1) The unweighted adjacency (which nodes are connected)
    # (2) The blue count per node
    # (3) The self-loop probability per node (from the neutral arc fraction)

    # Can we derive fiber?

    # From (3): p_self(mi) = neutral arcs / m. And fiber = W_row / m.
    # But W_row = fiber × m → circular again.

    # HOWEVER: The TOTAL tilings = 2^m.
    # And: Σ fiber(mi) = 2^m.
    # And: Σ blue(mi) = 2^{C(n-1,2)}.
    # These are LINEAR constraints on the fiber vector.

    # Plus: for SC classes, fiber = n!/|Aut|.
    # For non-SC pairs: fiber(mi) = 2 × n!/|Aut|.

    # The blue fraction: blue(mi)/fiber(mi) varies per class.
    # It depends on the specific automorphism structure.

    print("  BLUE FRACTION PER CLASS:")
    for mi in range(nm):
        H = H_vec[mi]; fib = fiber[mi]; blue = blue_counts[mi]
        frac = blue / fib if fib > 0 else 0
        is_sc = info[merged[mi][0]]['SC']
        print("    Node %d (H=%d, SC=%s): blue=%d/%d = %.4f" %
              (mi, H, is_sc, blue, fib, frac))

    # THE FORMULA: fiber = row_sum / m (requires weights)
    # WITHOUT weights: we need another approach.

    # CREATIVE APPROACH: H × |Aut| = H × (n!/fiber) is an invariant.
    # If we know H (computable), then fiber = n! / (H×|Aut|/H) = n!/|Aut|.
    # But we need |Aut|...

    # ORBIT-STABILIZER: fiber × |Aut| = n! (exact).
    # So: fiber = n!/|Aut|, and |Aut| = n!/fiber.
    # They're inverses. We need ONE of them to get the other.

    # FROM THE META-GRAPH: the DEGREE of a node tells us something.
    # deg(mi) = #{neighbors of mi} = #{classes reachable by 1 flip}.
    # From local=global: deg = E (every node has same neighbor set)... no,
    # deg varies! At n=5: deg ranges from 2 to 7.

    # Is there a formula: fiber = f(deg, H, blue)?

    print("\n  TESTING: fiber = f(degree, H, blue)")
    for mi in range(nm):
        H = H_vec[mi]; fib = fiber[mi]; blue = blue_counts[mi]
        deg = sum(1 for mj in range(nm) if mi != mj and W[mi][mj] > 0)
        # Try: fiber ∝ degree? fiber ∝ H? fiber ∝ blue?
        print("    Node %d: H=%d, deg=%d, blue=%d, fiber=%d, H×deg=%d, blue×m=%d" %
              (mi, H, deg, blue, fib, H*deg, blue*m))

print("\nDone. opus-2026-03-24-S279")
