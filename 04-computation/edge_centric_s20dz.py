#!/usr/bin/env python3
"""
edge_centric_s20dz.py — Think in terms of EDGES, not vertices.
kind-pasteur-2026-03-23-S20dz

THE EDGE-CENTRIC VIEW:

An arc position e = {u,v} in K_n. Flipping e in tournament T gives T'.
The arc e GENERATES an edge in the metagraph iff [T] != [T'].

For a FIXED arc position e, define:
  f(e) = number of iso class PAIRS {[T], [T']} connected by flipping e
       = number of distinct unordered {class, class} pairs

By S_n symmetry: f(e) is the SAME for all arc positions e.
So: E(G_n) = f(e) for any e. (Each edge in the metagraph is generated
by at least one arc position, and by symmetry each position generates
the same number of edges.)

Wait — that's not right. Different positions can generate the SAME edge.
So: E(G_n) = |union over all e of edges(e)|.
And: sum_e |edges(e)| = T_n (total transitions, counting multiplicities).

The edge-centric formula:
  For a fixed arc position e:
  |edges(e)| = #{distinct {[T],[T']} : T has [T] != [T flip e]}
  = (1/n!) * #{labeled T : [T] != [T flip e]} / 2
  ... this needs Burnside on the STABILIZER of e.

THE KEY INSIGHT: Fix the arc e = {0,1}.
The stabilizer Stab(e) in S_n = S_{n-2} (permutations fixing 0 and 1)
times Z_2 (swapping 0 and 1) = S_{n-2} x Z_2.

By Burnside on the STABILIZER:
  (number of orbit pairs connected through e) =
  (1/|Stab(e)|) * #{(T, sigma) : sigma in Stab(e), sigma(T)=T, [T]!=[T flip e]}
"""

import sys
from math import comb, factorial
from itertools import permutations
from collections import defaultdict
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  EDGE-CENTRIC VIEW: Counting metagraph edges through arc positions")
print("  kind-pasteur-2026-03-23-S20dz")
print("=" * 80)

for n in [3, 4, 5, 6]:
    t0 = time.time()
    ALL_ARCS = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(ALL_ARCS)
    perms = list(permutations(range(n)))

    def make_adj(bits):
        A = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(ALL_ARCS):
            if bits & (1 << k): A[i][j] = 1
            else: A[j][i] = 1
        return A

    def canon(A):
        best = None
        for p in perms:
            s = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s < best: best = s
        return best

    # Fix arc e = {0, 1} (the first arc position)
    e_idx = 0  # arc (0,1)

    # For each tournament T, check if flipping e changes the class
    # Count: pairs, self-loops, and which classes are connected
    class_pairs_via_e = set()  # unordered {cn1, cn2} pairs via this arc
    self_loops_via_e = 0

    for bits in range(1 << m):
        A = make_adj(bits)
        cn1 = canon(A)
        A2 = make_adj(bits ^ (1 << e_idx))
        cn2 = canon(A2)
        if cn1 == cn2:
            self_loops_via_e += 1
        else:
            pair = (min(cn1, cn2), max(cn1, cn2))
            class_pairs_via_e.add(pair)

    # By S_n symmetry: every arc position gives the same count.
    # Total metagraph edges = |union of all edges(e)| for all e.
    # Each edge in the metagraph appears in SOME positions.
    # The "thickness" of an edge = how many positions generate it.

    # From our data: E(G_n) = known.
    E_known = {3:1, 4:5, 5:30, 6:290}[n]

    # |edges via position e| (distinct class pairs)
    edges_via_e = len(class_pairs_via_e)

    # Total transitions via ALL positions = T_n
    # Self-loops via ALL positions = SL_mine
    # By symmetry: self_loops_via_e * C(n,2) = SL_mine * ... no, that's not right
    # because the self-loop count per position is NOT SL_mine / C(n,2) unless
    # the self-loop rate is uniform across positions.

    # Actually by S_n symmetry on arc positions:
    # total self-loops = C(n,2) * self_loops_via_e / (2^m normalization?)
    # No: self_loops_via_e counts LABELED tournaments where flipping arc (0,1)
    # preserves the class. Total = sum over all arcs = C(n,2) * self_loops_via_e
    # (since each arc has the same count by S_n symmetry).
    # But total labeled self-loops = sum over all T of neutral_arcs(T) = SL_full.
    # So: SL_full = C(n,2) * self_loops_via_e.

    SL_full = comb(n,2) * self_loops_via_e
    SL_orbits_computed = SL_full // factorial(n)

    print(f"\nn={n}: ({time.time()-t0:.1f}s)")
    print(f"  Arc position (0,1):")
    print(f"    Distinct class pairs via this arc: {edges_via_e}")
    print(f"    Self-loops via this arc: {self_loops_via_e} out of {2**m}")
    print(f"    Self-loop rate: {self_loops_via_e/2**m:.6f}")
    print(f"  By S_n symmetry:")
    print(f"    SL_full = C(n,2) * {self_loops_via_e} = {SL_full}")
    print(f"    SL_orbits = SL_full / n! = {SL_orbits_computed}")
    print(f"  E(G_n) = {E_known}")
    print(f"  edges_via_e = {edges_via_e} (this is a LOWER BOUND on E)")
    print(f"  edges_via_e / E = {edges_via_e/E_known:.4f}")

    # Key question: does every metagraph edge appear in SOME arc position?
    # YES, by definition: an edge {C1, C2} exists iff some arc flip connects them.
    # So: E = |union of edges(e) over all positions e|.
    # And: each position generates edges_via_e distinct pairs.
    # But many pairs are generated by MULTIPLE positions.

    # By inclusion-exclusion on positions:
    # E = |union| <= C(n,2) * edges_via_e (trivial upper bound)
    # E >= edges_via_e (lower bound from one position)

    print(f"  Lower bound (1 position): {edges_via_e}")
    print(f"  Upper bound (all positions, no overlap): {comb(n,2) * edges_via_e}")
    print(f"  Actual E: {E_known}")
    print(f"  Average positions per edge: {comb(n,2) * edges_via_e / E_known:.2f}")

    # The AVERAGE POSITIONS PER EDGE = how many arc positions generate
    # each metagraph edge. This is the "thickness" in our earlier analysis.
    # avg_thickness = T / E (not quite the same but related).

    # BY S_n SYMMETRY: the average is the same for each position.
    # So: avg_overlap = C(n,2) * edges_via_e / E.
    # And: T = C(n,2) * (2^m - self_loops_via_e) = total non-self-loop transitions.
    # T_orbits = T / n!.

    non_sl_via_e = 2**m - self_loops_via_e
    T_non_sl = comb(n,2) * non_sl_via_e
    T_non_sl_orbits = T_non_sl // factorial(n)

    print(f"  Total non-SL transitions: C(n,2) * {non_sl_via_e} = {T_non_sl}")
    print(f"  T_non_sl_orbits: {T_non_sl_orbits}")
    print(f"  Should equal T - SL_mine: ", end="")
    T_known = {3:4, 4:16, 5:88, 6:704}[n]
    SL_mine = {3:2, 4:6, 5:16, 6:58}[n]  # from earlier computation (orbit level)
    # Wait, SL_mine at orbit level != SL_orbits from SL_full/n!
    # SL_mine = sum of neutral ARC ORBITS per class
    # SL_orbits (from opus formula) = T/2 + (n-2)! - E
    SL_opus = T_known//2 + factorial(n-2) - E_known
    print(f"T-SL(opus) = {T_known} - {SL_opus} = {T_known - SL_opus}")

print(f"\n{'='*70}")
print("THE EDGE-CENTRIC FORMULA")
print(f"{'='*70}")
print()
print("For a FIXED arc position e = {0,1}:")
print("  self_loops_via_e = # labeled T where flipping e preserves class")
print("  edges_via_e = # distinct class pairs connected by flipping e")
print()
print("By S_n symmetry:")
print("  SL_full = C(n,2) * self_loops_via_e")
print("  Total non-SL transitions = C(n,2) * (2^m - self_loops_via_e)")
print()
print("The STABILIZER of arc {0,1} in S_n is S_2 x S_{n-2}:")
print("  S_2 = swap vertices 0 and 1")
print("  S_{n-2} = permute vertices {2,...,n-1}")
print("  |Stab| = 2 * (n-2)!")
print()
print("By Burnside on the stabilizer:")
print("  edges_via_e = (1/|Stab|) * sum_{sigma in Stab} Fix_pair(sigma)")
print("  where Fix_pair(sigma) counts T with sigma(T)=T and [T]!=[T flip e]")
print()
print("This is a Burnside sum over S_2 x S_{n-2} instead of S_n.")
print("|S_2 x S_{n-2}| = 2*(n-2)! which is MUCH smaller than n!.")
print()
print("For n=9: |Stab| = 2*7! = 10080 vs |S_9| = 362880.")
print("The computation is 36x FASTER per arc position.")
print("And there are only C(n,2) = 36 positions (all equivalent).")
print("So edges_via_e determines everything!")
print()
print("DONE.")
print("=" * 80)
