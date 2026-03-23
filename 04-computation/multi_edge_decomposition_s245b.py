#!/usr/bin/env python3
"""
multi_edge_decomposition_s245b.py -- opus-2026-03-23-S245b

COMPLETE ORBIT DECOMPOSITION: simple + multi + self-loop

edge_orbits = simple_cross + multi_cross + self_loop_orbits

where:
  simple_cross = E(G_n) = distinct class pairs connected
  multi_cross = additional orbits connecting same pair via different arcs
  self_loop_orbits = orbits where both endpoints are in the same class

We know: edge_orbits = T_n/2 + (n-2)! (the breakthrough formula)
We know: E(G_n) from direct computation through n=8
We computed: self_loop_orbits at n=3,4,5 via Burnside

So: multi_cross = edge_orbits - E - self_loop = (T/2 + (n-2)!) - E - SL_edge

This gives us ALL THREE components. The question: which has a formula?

Author: opus-2026-03-23-S245b
"""

import sys
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

banner("COMPLETE 3-WAY ORBIT DECOMPOSITION")

# Known values
T = {3:4, 4:16, 5:88, 6:704, 7:8912, 8:188288}
E = {3:1, 4:5, 5:30, 6:290, 7:4086, 8:91161}

# edge_orbits = T/2 + (n-2)!
EO = {n: T[n]//2 + factorial(n-2) for n in T}

# Compute self_loop_orbits at n=3,4,5 via direct Burnside
SL_edge = {}

for n in [3, 4, 5]:
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    pair_idx = {p: k for k, p in enumerate(pairs)}

    cmap = {}; tc = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        c = canon(A, n)
        if c not in cmap: cmap[c] = len(cmap)
        tc[bits] = cmap[c]

    fix_sl_sum = 0
    for perm_tuple in permutations(range(n)):
        perm = list(perm_tuple)
        def ap(bits):
            new_bits = 0
            for k, (i,j) in enumerate(pairs):
                pi, pj = perm[i], perm[j]
                if pi < pj:
                    ok = pair_idx[(pi,pj)]
                    if bits & (1<<ok): new_bits |= (1<<k)
                else:
                    ok = pair_idx[(pj,pi)]
                    if not (bits & (1<<ok)): new_bits |= (1<<k)
            return new_bits

        for bits in range(2**m):
            for k in range(m):
                flipped = bits ^ (1<<k)
                if flipped <= bits: continue
                if tc[bits] != tc[flipped]: continue  # only self-loops
                pb = ap(bits); pf = ap(flipped)
                if (pb == bits and pf == flipped) or (pb == flipped and pf == bits):
                    fix_sl_sum += 1

    SL_edge[n] = fix_sl_sum // factorial(n)

# Also compute self_loop_orbits at n=6 from the F matrix
# F[C][C] > 0 means class C has arcs reversing to itself.
# Self-loop edges at n=6: orbits of {T, T^e} with T ≅ T^e.
# We computed SL_edge(5) = 14 above. For n=6, we need a different approach.

# From S245: SL_A + SL_B. Case B = (n-2)! always self-loops.
# Case A self-loop orbits = SL_A.
# From n=5 data: SL_A = 8, SL_B = 6, SL_edge = 14.
# We need SL_A at n=6. This requires the Burnside computation at n=6
# which is too expensive for full permutation enumeration (720! permutations
# times 32768 tournaments × 15 arcs). Let me estimate.

# Actually, from the orbit accounting:
# SL_edge(n) = SL_A(n) + SL_B(n)
# SL_B(n) = (n-2)! (verified for n=3,4,5)
# multi_cross(n) = Case_A orbits that connect already-connected class pairs
# E(n) = Case_A orbits that connect NEW class pairs

# Total Case A = T/2.
# Case A = simple_cross + multi_cross + SL_A
# T/2 = E + multi + SL_A

# So: SL_A = T/2 - E - multi
# And: SL_edge = SL_A + (n-2)! = T/2 - E - multi + (n-2)!
# And: multi = T/2 - E - SL_A = T/2 - E - (SL_edge - (n-2)!)

# From our data at n=3,4,5:
print("  n  | edge_orb | E    | SL_edge | multi | SL_A | SL_B=(n-2)!")
print("  ---|----------|------|---------|-------|------|--------")
for n in [3, 4, 5]:
    eo = EO[n]
    e = E[n]
    sl = SL_edge[n]
    multi = eo - e - sl
    sl_b = factorial(n-2)
    sl_a = sl - sl_b
    print("  %d  | %7d  | %4d | %7d | %5d | %4d | %4d" %
          (n, eo, e, sl, multi, sl_a, sl_b))

# Let me verify: T/2 = E + multi + SL_A
for n in [3, 4, 5]:
    eo = EO[n]
    e = E[n]
    sl = SL_edge[n]
    multi = eo - e - sl
    sl_b = factorial(n-2)
    sl_a = sl - sl_b
    t_half = T[n] // 2
    check = e + multi + sl_a
    print("  n=%d: T/2=%d, E+multi+SL_A=%d+%d+%d=%d, match=%s" %
          (n, t_half, e, multi, sl_a, check, t_half == check))

banner("THE THREE SEQUENCES")

# simple_cross (= E): 1, 5, 30, 290, 4086, 91161
# multi_cross: ?, ?, ?, ...
# self_loop_edge: 2, 5, 14, ...

# From the formula:
# multi = T/2 - E - SL_A = T/2 - E - (SL_edge - (n-2)!)
# If SL_B = (n-2)! holds for all n:
# multi = T/2 - E - SL_edge + (n-2)!
# = (T/2 + (n-2)!) - E - SL_edge
# = edge_orbits - E - SL_edge
# = D_orbits - SL_edge  where D_orbits = edge_orbits - E

print("  The D_orbits sequence (= multi + SL_edge = edge_orbits - E):")
for n in range(3, 9):
    d_orb = EO[n] - E[n]
    print("    n=%d: D_orbits = %d" % (n, d_orb))

D_orb = {n: EO[n] - E[n] for n in range(3, 9)}

print("\n  D_orbits = %s" % [D_orb[n] for n in range(3, 9)])
# 2, 5, 20, 86, 490, 3703

# And SL_edge = 2, 5, 14, ?, ?, ?
# So multi = D_orbits - SL_edge = 0, 0, 6, ?, ?, ?

print("\n  Multi-edge orbit counts:")
for n in [3, 4, 5]:
    multi = D_orb[n] - SL_edge[n]
    print("    n=%d: multi = D_orbits(%d) - SL_edge(%d) = %d" %
          (n, D_orb[n], SL_edge[n], multi))

# multi: 0, 0, 6
# These are edges where TWO DIFFERENT arc-orbits connect the same class pair.

banner("WHAT ARE THE 6 MULTI-EDGES AT n=5?")

# At n=5: 6 multi-edge orbits. Let me find them.
n = 5; m = comb(n, 2)
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
pair_idx = {p: k for k, p in enumerate(pairs)}

cmap = {}; tc = {}; sizes = {}; reps = {}
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

# For each class pair (C, D) with C ≠ D:
# Count the number of edge orbits connecting them.
# An edge orbit connecting C to D = an S_n orbit on {(T, e) : T ∈ C, T^e ∈ D}
# Each such orbit maps to the undirected edge {C, D} in G_n.
# If there are k orbits → 1 edge + (k-1) multi-edges.

# From the F matrix and orbit structure, the number of arc-orbits from C to D
# is mult(C→D). Each arc-orbit is also an edge orbit (from Case A).
# So the number of Case A edge orbits connecting {C, D} = mult(C→D) + mult(D→C).
# If this sum > 1, there are multi-edges.

# Wait — each arc-orbit from C to D gives a DIRECTED edge orbit.
# Two directed orbits C→D and D→C give ONE undirected edge {C,D} (if they're paired).
# Or they might be in the same undirected orbit.

# For UNORDERED edges: {T, T^e} is the same as {T^e, T}.
# So an arc-orbit (C→D) and the corresponding "reverse" orbit (D→C) might be
# in the same unordered edge orbit, or they might be different.

# This is Case A vs Case B again.
# Case A: σ fixes both T and T^e → the arc-orbit C→D and the reverse D→C
#   are in DIFFERENT unordered edge orbits (σ fixes each separately).
# Case B: σ swaps T and T^e → pairs the arc-orbits (C→D and D→C are in the
#   SAME unordered edge orbit).

# For a cross-edge {C, D} with C ≠ D:
# The number of unordered edge orbits from Case A connecting {C, D}:
#   = (mult(C→D) + mult(D→C)) / 2 if all arc-orbits pair up
#   Actually: it's more complex. Each Case A unordered orbit needs
#   σ fixing both T∈C and T^e∈D. If mult(C→D) = mult(D→C), then
#   the pairing might work, but if they differ, some orbits are unpaired.

# This is getting complicated. Let me just count directly.

# Build all edge orbits
edge_to_orbit = {}
n_orbits = 0

for bits in range(2**m):
    for k in range(m):
        flipped = bits ^ (1<<k)
        if flipped <= bits: continue
        edge = (bits, flipped)
        if edge in edge_to_orbit: continue

        # Find canonical form of this edge under S_n
        best_edge = None
        for perm_tuple in permutations(range(n)):
            perm = list(perm_tuple)
            def ap(b):
                new_b = 0
                for kk, (i,j) in enumerate(pairs):
                    pi, pj = perm[i], perm[j]
                    if pi < pj:
                        ok = pair_idx[(pi,pj)]
                        if b & (1<<ok): new_b |= (1<<kk)
                    else:
                        ok = pair_idx[(pj,pi)]
                        if not (b & (1<<ok)): new_b |= (1<<kk)
                return new_b

            new_a = ap(bits); new_b = ap(flipped)
            e = (min(new_a, new_b), max(new_a, new_b))
            if best_edge is None or e < best_edge:
                best_edge = e

        if best_edge not in edge_to_orbit:
            edge_to_orbit[best_edge] = n_orbits
            n_orbits += 1
        edge_to_orbit[edge] = edge_to_orbit[best_edge]

# Now classify each orbit
orbit_classes = {}  # orbit_id -> (class_i, class_j) of the edge
for (bits, flipped), orbit_id in edge_to_orbit.items():
    if orbit_id not in orbit_classes:
        ci = tc[bits]; cj = tc[flipped]
        orbit_classes[orbit_id] = (min(ci,cj), max(ci,cj))

# Count orbits per class pair
pair_orbit_count = Counter()
for oid, (ci, cj) in orbit_classes.items():
    pair_orbit_count[(ci, cj)] += 1

# Multi-edge pairs: those with count > 1 (for cross) or count > 0 (for self-loops)
print("  Class pairs with multiple edge orbits (at n=5):\n")
multi_pairs = 0
for (ci, cj), count in sorted(pair_orbit_count.items()):
    if ci == cj:
        label = "SELF-LOOP"
        if count > 0:
            print("    {%d, %d}: %d orbits (%s)" % (ci, cj, count, label))
    elif count > 1:
        label = "MULTI-EDGE"
        multi_pairs += count - 1
        print("    {%d, %d}: %d orbits (%s, %d extra)" % (ci, cj, count, label, count-1))

self_loop_total = sum(c for (ci,cj), c in pair_orbit_count.items() if ci == cj)
simple_cross = sum(1 for (ci,cj), c in pair_orbit_count.items() if ci != cj and c >= 1)
multi_total = sum(c-1 for (ci,cj), c in pair_orbit_count.items() if ci != cj and c > 1)

print("\n  Summary:")
print("    Total edge orbits: %d" % n_orbits)
print("    Self-loop orbits: %d" % self_loop_total)
print("    Simple cross (= E): %d" % simple_cross)
print("    Multi-edge orbits: %d" % multi_total)
print("    Check: %d + %d + %d = %d (should be %d)" %
      (self_loop_total, simple_cross, multi_total,
       self_loop_total + simple_cross + multi_total, n_orbits))

print("\nDone. opus-2026-03-23-S245b")
