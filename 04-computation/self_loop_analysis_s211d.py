#!/usr/bin/env python3
"""
self_loop_analysis_s211d.py -- opus-2026-03-23-S211d

CORRECTING THE EDGE FORMULA: account for self-loops.

The formula E = (T-D)/2 is WRONG because it ignores self-loops.
Self-loops = arc-orbits where reversing brings you back to the same class.

The correct accounting:
T = SL + cross_total
where SL = sum_C mult(C->C) (self-loop orbits)
      cross_total = sum_{C≠D} mult(C->D)

And E = number of distinct unordered pairs {C,D} with C≠D that are connected.

The formula E = (T-D)/2 = (T - (S_2-T)/2) / 2 would require
SL + excess_cross = D_self + D_cross, which need not hold.

This script computes the EXACT relationship and finds what formula DOES work.

Author: opus-2026-03-23-S211d
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

def canonical(A, n):
    if n <= 5:
        best = None
        for perm in permutations(range(n)):
            form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j)
            if best is None or form < best: best = form
        return best
    else:
        scores = [sum(A[i]) for i in range(n)]
        score_groups = defaultdict(list)
        for v in range(n): score_groups[scores[v]].append(v)
        groups = [score_groups[s] for s in sorted(score_groups.keys())]
        best = None
        def gen_perms(grps):
            if not grps:
                yield []
                return
            for p in permutations(grps[0]):
                for rest in gen_perms(grps[1:]):
                    yield list(p) + rest
        for perm in gen_perms(groups):
            form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j)
            if best is None or form < best: best = form
        return best

def analyze(n):
    """Full analysis: T, S_2, D, E, SL, cross, excess for given n."""
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

    # Build iso classes
    cmap = {}; sizes = {}; reps = {}; tc = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        c = canonical(A, n)
        if c not in cmap:
            cid = len(cmap); cmap[c] = cid; sizes[cid] = 0
            reps[cid] = [row[:] for row in A]
        sizes[cmap[c]] += 1
        tc[bits] = cmap[c]
    nc = len(cmap)

    # Build F matrix
    F = np.zeros((nc, nc), dtype=int)
    for bits in range(2**m):
        ci = tc[bits]
        for k in range(m):
            F[ci][tc[bits^(1<<k)]] += 1

    # E from adjacency
    adj = (F > 0).astype(int)
    np.fill_diagonal(adj, 0)
    E = np.sum(adj) // 2

    # For each class: find arc-orbits, self-loops, cross-orbits
    T_total = 0
    SL_total = 0  # self-loop orbits
    cross_total = 0  # non-self-loop orbits
    D_self_total = 0
    D_cross_total = 0
    S2_total = 0
    S2_self_total = 0
    S2_cross_total = 0

    # Also track: for each undirected edge, what's the directed mult sum?
    edge_mults = {}  # (min(C,D), max(C,D)) -> (mult_CD, mult_DC)

    for cid in range(nc):
        A = reps[cid]
        # Automorphisms
        auts = []
        for perm in permutations(range(n)):
            if all(A[i][j] == A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j):
                auts.append(perm)

        # Arc orbits
        arcs = [(i,j) for i in range(n) for j in range(n) if i!=j and A[i][j]]
        seen = set(); orbits = []
        for arc in arcs:
            if arc in seen: continue
            orbit = set()
            for perm in auts:
                orbit.add((perm[arc[0]], perm[arc[1]]))
            orbits.append(orbit)
            seen.update(orbit)

        # Targets
        targets = []
        for orbit in orbits:
            e = min(orbit)
            A2 = [row[:] for row in A]
            A2[e[0]][e[1]] = 0; A2[e[1]][e[0]] = 1
            targets.append(cmap[canonical(A2, n)])

        n_orbits = len(orbits)
        T_total += n_orbits

        mult = Counter(targets)
        sl = mult.get(cid, 0)  # self-loops for this class
        SL_total += sl
        cross_total += n_orbits - sl

        # D decomposition
        for target, m_val in mult.items():
            d_contrib = comb(m_val, 2)
            s2_contrib = m_val**2
            S2_total += s2_contrib
            if target == cid:
                D_self_total += d_contrib
                S2_self_total += s2_contrib
            else:
                D_cross_total += d_contrib
                S2_cross_total += s2_contrib

                # Track edge multiplicities
                key = (min(cid, target), max(cid, target))
                if key not in edge_mults:
                    edge_mults[key] = [0, 0]
                if cid < target:
                    edge_mults[key][0] = m_val
                else:
                    edge_mults[key][1] = m_val

    D_total = D_self_total + D_cross_total

    # Cross-edge excess
    excess_cross = 0
    for (a, b), (m1, m2) in edge_mults.items():
        excess_cross += (m1 + m2 - 2)

    return {
        'n': n, 'nc': nc, 'E': E, 'T': T_total,
        'SL': SL_total, 'cross': cross_total,
        'D': D_total, 'D_self': D_self_total, 'D_cross': D_cross_total,
        'S2': S2_total, 'S2_self': S2_self_total, 'S2_cross': S2_cross_total,
        'excess_cross': excess_cross,
        'n_edges_with_mult': len(edge_mults)
    }

# ============================================================
banner("COMPLETE ANALYSIS: n=3,4,5,6")
# ============================================================

for n in [3, 4, 5, 6]:
    r = analyze(n)
    print("  n=%d (V=%d):" % (r['n'], r['nc']))
    print("    T=%d = SL(%d) + cross(%d)" % (r['T'], r['SL'], r['cross']))
    print("    E=%d (from F-matrix adjacency)" % r['E'])
    print("    D=%d = D_self(%d) + D_cross(%d)" % (r['D'], r['D_self'], r['D_cross']))
    print("    S_2=%d = S2_self(%d) + S2_cross(%d)" % (r['S2'], r['S2_self'], r['S2_cross']))
    print("    excess_cross=%d (directed orbits beyond 2 per edge)" % r['excess_cross'])
    print()

    # Check formulas
    e_naive = (r['T'] - r['D']) // 2
    e_formula_check = (3*r['T'] - r['S2']) // 4

    # The correct formula should be:
    # cross_total = 2E + excess_cross
    # T = SL + cross = SL + 2E + excess_cross
    # => E = (T - SL - excess_cross) / 2
    e_correct = (r['T'] - r['SL'] - r['excess_cross']) // 2

    # Also: excess_cross = sum over edges of (mult_forward + mult_back - 2)
    #   = cross_total - 2E = (T - SL) - 2E
    # D_cross = sum over directed (C,D) C≠D of C(mult,2)
    # D_self = sum over C of C(mult(C→C), 2)

    print("    E_naive = (T-D)/2 = %d  (WRONG?)" % e_naive)
    print("    E_correct = (T - SL - excess)/2 = %d" % e_correct)
    print("    E_formula2 = (3T-S2)/4 = %d" % e_formula_check)
    print("    actual E = %d" % r['E'])
    print("    Match naive: %s, Match correct: %s" %
          (e_naive == r['E'], e_correct == r['E']))

    # Can we express E in terms of T, SL, D_self, D_cross?
    # E = (T - SL - excess_cross) / 2
    # excess_cross = T - SL - 2E
    # D_cross = sum of C(mult,2) for cross-pairs
    # For each cross-pair: C(mult,2) = (mult^2 - mult)/2

    # sum over directed (C→D, C≠D): mult(C→D)
    #   = cross_total = T - SL
    # sum over directed (C→D, C≠D): mult(C→D)^2
    #   = S2_cross
    # sum over directed: C(mult, 2) = (S2_cross - cross) / 2 = D_cross

    # For each UNDIRECTED edge {A,B}:
    # directed orbits = mult(A→B) + mult(B→A)
    # The MINIMUM is 2 (one in each direction).
    # excess = sum_edges (mult_forward + mult_back - 2)
    #        = cross_total - 2E

    # Can we get excess from D_cross?
    # D_cross = sum_{A≠B} C(mult(A→B), 2)
    #         = sum_{edges} [C(k_AB, 2) + C(k_BA, 2)]
    # This is NOT the same as excess = sum_edges (k_AB + k_BA - 2)

    # Try: is D_self related to SL in a simple way?
    print("    SL=%d, D_self=%d, SL-2*D_self=%d" %
          (r['SL'], r['D_self'], r['SL'] - 2*r['D_self']))

    # Check: SL + excess_cross vs D_self + D_cross vs D:
    diff = r['SL'] + r['excess_cross'] - r['D_self'] - r['D_cross']
    print("    SL + excess - D = %d + %d - %d = %d (should be 0 if E=(T-D)/2 worked)" %
          (r['SL'], r['excess_cross'], r['D'], diff))
    print()

# ============================================================
banner("THE CORRECTED FORMULA")
# ============================================================

print("""
  THE FORMULA E = (T - D) / 2 IS WRONG.

  It fails at n=6 (and actually at n=3 too, masked by integer division).

  The correct decomposition:
  T = SL + 2*E + excess_cross
  where:
    SL = sum_C mult(C -> C) = self-loop orbits
    excess_cross = sum_{{C,D}} (mult(C->D) + mult(D->C) - 2) for C≠D connected
                 = cross_total - 2*E

  So: E = (T - SL - excess_cross) / 2

  Both SL and excess_cross require computing orbit targets, not just counts.
  There is no closed-form expression for E from T and D alone.
""")

# ============================================================
banner("CAN WE COMPUTE SL FROM BURNSIDE?")
# ============================================================

# SL = number of arc-orbits within each class that are self-loops.
# An arc-orbit of class C is a self-loop iff reversing any arc in the orbit
# gives back a tournament in class C.

# This is related to the "self-adjacency" of the class:
# F[C][C] = size(C) * (arcs per T going to same class).
# SL(C) = number of arc-orbits targeting C itself.

# From F: arcs_to_self_per_T(C) = F[C][C] / size(C)
# But arc-orbits_to_self ≠ arcs_to_self / |Aut|  (orbits have varying sizes!)

# However: arcs_to_self_per_T = sum over self-loop orbits of |orbit|
# And number_of_self_loop_orbits = SL(C).
# If |Aut| = 1: all orbits have size 1, so SL(C) = arcs_to_self = F[C][C]/size(C).
# If |Aut| > 1: SL(C) = (number of arcs going to self class) / avg_orbit_size_for_self_orbits.

# For the TOTAL: SL = sum_C SL(C).

# Can we compute SL from F without knowing orbit structure?
# Only if we know how arcs-to-self decompose into orbits.

# But for generic tournaments (|Aut|=1, which is almost all at large n):
# SL(C) = F[C][C] / size(C) = diagonal of the "mult" matrix.
# And F[C][C]/size(C) arcs per tournament hit the same class.

# So for |Aut|=1 classes: SL(C) = F[C][C] / size(C) (exact)
# For |Aut|>1 classes: SL(C) = (Burnside on self-arcs)

# Let's check: how many classes have |Aut| > 1?
for n in [3, 4, 5, 6]:
    m = comb(n, 2)
    pairs_n = [(i,j) for i in range(n) for j in range(i+1,n)]
    cmap_n = {}; sizes_n = {}; reps_n = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(ii,jj) in enumerate(pairs_n):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        c = canonical(A, n)
        if c not in cmap_n:
            cid = len(cmap_n); cmap_n[c] = cid; sizes_n[cid] = 0
            reps_n[cid] = [row[:] for row in A]
        sizes_n[cmap_n[c]] += 1
    nc_n = len(cmap_n)

    nontrivial = 0
    for cid in range(nc_n):
        A = reps_n[cid]
        aut_size = sum(1 for perm in permutations(range(n))
                       if all(A[i][j] == A[perm[i]][perm[j]]
                              for i in range(n) for j in range(n) if i!=j))
        if aut_size > 1:
            nontrivial += 1

    print("  n=%d: %d/%d classes have |Aut| > 1" % (n, nontrivial, nc_n))

# ============================================================
banner("THE EDGE FORMULA: WHAT WE CAN COMPUTE")
# ============================================================

# From the Burnside framework, we can compute:
# T(n) via cycle index (DONE, exact formula)
# V(n) via cycle index (DONE, A000568)
# F[C][C] at the class level (requires building the meta-graph)

# The issue: E requires the ADJACENCY pattern, not just orbit counts.
# D and S_2 capture orbit-count statistics but miss the self-loop structure.

# HOWEVER: at large n, SL → 0 as a fraction of T, and excess → 0.
# So E ~ T/2 asymptotically, and the formula E ≈ (T-D)/2 ≈ T/2.

# The precise relationship:
# T - 2E = SL + excess_cross
# So E = (T - SL - excess_cross) / 2

# At small n, SL and excess are significant.
# At large n, both become negligible.

# Let's quantify:
for n in [3, 4, 5, 6]:
    r = analyze(n)
    gap = r['T'] - 2*r['E']
    print("  n=%d: T-2E = %d = SL(%d) + excess(%d), (T-2E)/T = %.4f" %
          (n, gap, r['SL'], r['excess_cross'], gap/r['T']))

# From the known values at higher n:
T_known = {3:4, 4:16, 5:88, 6:704, 7:8912, 8:188288, 9:6847200}
E_known = {3:1, 4:5, 5:30, 6:290, 7:4086, 8:91161, 9:3380751}
print("\n  At higher n (from known T and E):")
for n in range(3, 10):
    gap = T_known[n] - 2*E_known[n]
    print("  n=%d: T-2E = %d, (T-2E)/T = %.6f" %
          (n, gap, gap/T_known[n]))

print("\n  T-2E = SL + excess_cross. This sum DECREASES relative to T.")
print("  So E ~ T/2 becomes exact at large n.")

print("\nDone. opus-2026-03-23-S211d")
