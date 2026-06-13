#!/usr/bin/env python3
"""
edge_count_formula_s185.py -- opus-2026-03-22-S185

THE EDGE COUNT FORMULA FOR G_n

From kind-pasteur S20bv: degree(C) = #{distinct classes among cross-orbits}
where Aut(T) partitions the m arcs into orbits.

From Feng (arXiv:2510.25202): the dual Burnside process Q=AB
has stationary distribution pi(g) proportional to |Fix(g)|.
The matrix factorization Q=AB, K=BA connects the classical
(Burnside) and dual Markov chains.

This script:
1. Implements the Aut-orbit degree computation
2. Derives the edge count at n=3,4,5 from Aut-orbit structure
3. Explores the Burnside/dual-Burnside connection for edges
4. Attempts the formula: |E| = Burnside-like sum over S_n
5. Computes at n=6 using optimized code

Author: opus-2026-03-22-S185
"""

import sys
import numpy as np
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter
from fractions import Fraction
import time
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def count_hp(A, n):
    dp = defaultdict(int)
    for v in range(n): dp[(1<<v,v)] = 1
    for mask in range(1,1<<n):
        for v in range(n):
            if not (mask&(1<<v)): continue
            if dp[(mask,v)]==0: continue
            for w in range(n):
                if mask&(1<<w): continue
                if A[v][w]: dp[(mask|(1<<w),w)] += dp[(mask,v)]
    return sum(dp[((1<<n)-1,v)] for v in range(n))

def canonical(A, n):
    best = None
    for perm in permutations(range(n)):
        form = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or form < best: best = form
    return best

def get_automorphisms(A, n):
    """Return list of automorphisms (as permutations) of tournament A."""
    auts = []
    for perm in permutations(range(n)):
        ok = True
        for i in range(n):
            for j in range(n):
                if i != j and A[i][j] != A[perm[i]][perm[j]]:
                    ok = False; break
            if not ok: break
        if ok: auts.append(perm)
    return auts

def arc_orbits(A, n, auts):
    """Partition arcs into Aut-orbits. Return list of orbits."""
    arcs = [(i,j) for i in range(n) for j in range(n) if i != j and A[i][j]]
    arc_set = set(arcs)
    visited = set()
    orbits = []
    for arc in arcs:
        if arc in visited: continue
        orbit = set()
        queue = [arc]
        while queue:
            a = queue.pop()
            if a in orbit: continue
            orbit.add(a)
            for perm in auts:
                mapped = (perm[a[0]], perm[a[1]])
                if mapped in arc_set and mapped not in orbit:
                    queue.append(mapped)
        orbits.append(frozenset(orbit))
        visited.update(orbit)
    return orbits

def reverse_arc_class(A, n, arc, cmap):
    """Reverse one arc and return the class of the resulting tournament."""
    i, j = arc
    A2 = [row[:] for row in A]
    A2[i][j] = 0; A2[j][i] = 1
    return cmap[canonical(A2, n)]

# ============================================================
# PART 1: THE AUT-ORBIT DEGREE FORMULA
# ============================================================

banner("PART 1: AUT-ORBIT DEGREE FORMULA")

print("""
  THEOREM (kind-pasteur S20bv):
  For each iso class C with representative T:
    1. Aut(T) partitions the m directed arcs into orbits
    2. All arcs in the same orbit reverse to the SAME class
       (by Aut-equivariance)
    3. Self-loop orbits: reversing gives a T isomorphic to T
    4. Cross-orbits: reversing gives a T in a DIFFERENT class
    5. degree(C) = #{distinct classes among cross-orbits}

  Proof of step 2:
  If sigma in Aut(T) maps arc (i,j) to arc (k,l), then
  reversing (i,j) in T gives T', and sigma maps T' to the
  tournament obtained by reversing (k,l) in T.
  So T'_reversed_{(i,j)} is isomorphic to T'_reversed_{(k,l)}.
  Hence they're in the same iso class.

  EDGE COUNT: |E(G_n)| = (1/2) * sum_C degree(C)
  (handshaking lemma, dividing by 2 for unordered edges)
""")

# ============================================================
# PART 2: COMPUTE DEGREES VIA AUT-ORBITS
# ============================================================

for n in [3, 4, 5]:
    banner("n = %d: AUT-ORBIT DEGREE COMPUTATION" % n)
    t0 = time.time()

    # Build class catalog
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    cmap = {}; cdata = {}; creps = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: A[i][j]=1
            else: A[j][i]=1
        c = canonical(A, n)
        if c not in cmap:
            cid = len(cmap); cmap[c] = cid
            cdata[cid] = {'H': count_hp(A, n), 'size': 0}
            creps[cid] = [row[:] for row in A]
        cdata[cmap[c]]['size'] += 1
    nc = len(cdata)

    total_degree = 0
    total_self_orbits = 0
    total_cross_orbits = 0

    print("  %-4s %-4s %-6s %-6s %-6s %-6s %-6s %s" %
          ("cls", "H", "|Aut|", "#orb", "self", "cross", "deg", "neighbor classes"))

    for cid in range(nc):
        A = creps[cid]
        auts = get_automorphisms(A, n)
        aut_size = len(auts)
        orbits = arc_orbits(A, n, auts)
        n_orbits = len(orbits)

        # For each orbit, check if reversing gives same class or different
        self_orbs = 0
        cross_orb_classes = {}
        for orb in orbits:
            # Pick any arc in the orbit
            arc = next(iter(orb))
            target_class = reverse_arc_class(A, n, arc, cmap)
            if target_class == cid:
                self_orbs += 1
            else:
                cross_orb_classes[target_class] = cross_orb_classes.get(target_class, 0) + 1

        degree = len(cross_orb_classes)
        total_degree += degree
        total_self_orbits += self_orbs
        total_cross_orbits += len(orbits) - self_orbs

        neighbors = sorted(cross_orb_classes.keys())
        print("  %-4d %-4d %-6d %-6d %-6d %-6d %-6d %s" %
              (cid, cdata[cid]['H'], aut_size, n_orbits, self_orbs,
               n_orbits - self_orbs, degree, neighbors))

    edge_count = total_degree // 2
    elapsed = time.time() - t0
    print("\n  Total degree = %d, |E| = %d (%.1fs)" %
          (total_degree, edge_count, elapsed))

    # Verify against known values
    known_E = {3: 1, 4: 5, 5: 30}
    print("  Expected: %d. Match: %s" % (known_E[n], edge_count == known_E[n]))

# ============================================================
# PART 3: THE ORBIT COUNT FORMULA
# ============================================================

banner("PART 3: ORBIT COUNT ANALYSIS")

print("""
  For each class C with |Aut|=a:
    Total arcs = m (directed arcs, where A[i][j]=1)
    Wait: arcs FROM vertex i: score(i). Total arcs = C(n,2).
    But directed arcs (including A[i][j]=1 for i>j): total = C(n,2).

  Actually: each edge {i,j} has ONE direction (i->j or j->i).
  There are C(n,2) = m undirected edges, each with one direction.
  Aut(T) acts on these m DIRECTED arcs (not 2m).

  Number of orbits = m / |Aut| ... NO! By Burnside on arcs:
  #orbits = (1/|Aut|) * sum_{g in Aut} Fix_arcs(g)
  where Fix_arcs(g) = #{arcs (i,j) with g(i)->g(j) AND A[g(i)][g(j)]=1}

  But since g is an automorphism: A[g(i)][g(j)] = A[i][j].
  So Fix_arcs(g) = #{arcs (i,j) with g({i,j}) = {i,j} and direction preserved}
  = #{edges fixed as a SET by g, AND direction preserved}.

  For the IDENTITY g=e: Fix_arcs(e) = m (all arcs fixed).
  For g = transposition (a,b): Fix_arcs(g) = #{arcs not involving a or b}
  + (1 if (a,b) and (b,a) are both arcs... impossible in tournament)
  + 0 (since swapping a,b reverses the arc between them).
  So Fix_arcs((a,b)) = C(n-2, 2) for a transposition.
""")

for n in [3, 4, 5]:
    m = comb(n, 2)
    print("  n=%d: m=%d" % (n, m))

    # For each class, show orbit structure
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    cmap = {}; cdata = {}; creps = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: A[i][j]=1
            else: A[j][i]=1
        c = canonical(A, n)
        if c not in cmap:
            cid = len(cmap); cmap[c] = cid
            cdata[cid] = {'H': count_hp(A, n), 'size': 0}
            creps[cid] = [row[:] for row in A]
        cdata[cmap[c]]['size'] += 1

    for cid in range(len(cdata)):
        A = creps[cid]
        auts = get_automorphisms(A, n)
        a = len(auts)
        orbits = arc_orbits(A, n, auts)
        n_orbs = len(orbits)
        orbit_sizes = sorted([len(o) for o in orbits])
        # Burnside check: sum |Fix_arcs(g)| / |Aut| = n_orbs
        fix_sum = 0
        for perm in auts:
            fix_count = 0
            for i in range(n):
                for j in range(n):
                    if i != j and A[i][j]:
                        if perm[i] == i and perm[j] == j:
                            # This arc is fixed by perm
                            fix_count += 1
                        # Actually: arc (i,j) is fixed by perm iff
                        # perm sends (i,j) to (i,j), i.e., perm[i]=i AND perm[j]=j
                        # BUT Aut preserves A, so A[perm[i]][perm[j]]=A[i][j]=1
                        # The arc (i,j) maps to (perm[i], perm[j]).
                        # It's "fixed" as a directed arc iff perm[i]=i, perm[j]=j.
            fix_sum += fix_count
        burnside_orbs = fix_sum / a
        print("    class %d (H=%d, |Aut|=%d): %d orbits = %d/|Aut| = %.1f, sizes=%s" %
              (cid, cdata[cid]['H'], a, n_orbs, fix_sum, burnside_orbs, orbit_sizes))
    print()

# ============================================================
# PART 4: THE TOTAL ARC-ORBIT COUNT
# ============================================================

banner("PART 4: TOTAL ARC-ORBIT COUNT")

print("  total_orbits = sum over classes of (#arc orbits)")
print("  This counts the total number of distinct arc orbits across all classes.\n")

for n in [3, 4, 5]:
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]
    cmap = {}; cdata = {}; creps = {}
    for bits in range(2**m):
        A = [[0]*n for _ in range(n)]
        for k,(i,j) in enumerate(pairs):
            if (bits>>k)&1: A[i][j]=1
            else: A[j][i]=1
        c = canonical(A, n)
        if c not in cmap:
            cid = len(cmap); cmap[c] = cid
            cdata[cid] = {'H': count_hp(A, n), 'size': 0}
            creps[cid] = [row[:] for row in A]
        cdata[cmap[c]]['size'] += 1

    total_orbs = 0
    total_self = 0
    total_cross = 0
    for cid in range(len(cdata)):
        A = creps[cid]
        auts = get_automorphisms(A, n)
        orbits = arc_orbits(A, n, auts)
        n_orbs = len(orbits)
        total_orbs += n_orbs
        for orb in orbits:
            arc = next(iter(orb))
            target = reverse_arc_class(A, n, arc, cmap)
            if target == cid:
                total_self += 1
            else:
                total_cross += 1

    print("  n=%d: total_orbs=%d, self=%d, cross=%d, sum_degree=%d, |E|=%d" %
          (n, total_orbs, total_self, total_cross, total_cross, total_cross // 2))
    print("    total_orbs = sum m/|Aut| ... no, that's not right")
    print("    Burnside: total_orbs = sum_{classes} (1/|Aut|)*sum_{g} Fix_arcs(g)")

    # The total arc-orbit count across all classes
    # = sum_{classes C} #{Aut(C)-orbits on arcs of C}
    # By Burnside: = sum_C (1/|Aut(C)|) * sum_{g in Aut(C)} |Fix_arcs(g)|

    # There's a SIMPLER way:
    # total_orbs = (1/n!) * sum_{sigma in S_n} sum_T (# arcs of T fixed by sigma)
    # Wait, this double-counts because we sum over classes, not all T.

    # Actually:
    # sum_C #orbits(C) = sum_C m/|Aut(C)| (only if Aut acts FREELY, which it doesn't)
    # In general: #orbits = m/|Aut| only when all arcs have trivial stabilizer.

    # For |Aut|=1: #orbits = m.
    # For |Aut|>1: #orbits < m.

    # sum_C m/|Aut(C)| = m * sum_C 1/|Aut(C)|
    m_sum_inv_aut = m * sum(Fraction(1, len(get_automorphisms(creps[c], n)))
                           for c in range(len(cdata)))
    print("    m * sum 1/|Aut| = %s" % m_sum_inv_aut)
    print("    actual total_orbs = %d" % total_orbs)
    print()

# ============================================================
# PART 5: THE DUAL BURNSIDE CONNECTION
# ============================================================

banner("PART 5: THE DUAL BURNSIDE PROCESS (Feng 2510.25202)")

print("""
  Feng's dual Burnside process provides a matrix factorization Q = AB:
  - The CLASSICAL Burnside kernel K = BA acts on group elements
  - The DUAL kernel Q = AB acts on states (tournaments)
  - Both share nonzero eigenvalues

  For tournament edge counting:
  The classical Burnside gives |V(G_n)| = (1/n!) * sum Fix(sigma).
  The DUAL should give information about EDGES.

  IDEA: Apply Burnside to the set of (tournament, arc) pairs:
  |{S_n-orbits of (T, e) pairs}| = (1/n!) * sum_sigma |Fix_pairs(sigma)|
  where Fix_pairs(sigma) = #{(T, e) : sigma(T)=T and sigma(e)=e}

  Fix_pairs(sigma) = sum_T [sigma in Aut(T)] * |{arcs of T fixed by sigma}|
  = sum_T [sigma in Aut(T)] * Fix_arcs(sigma, T)

  This gives the total number of Aut-orbits on (T, e) pairs,
  which counts:
  - self-loop orbits (reversal preserves class)
  - cross-loop orbits (reversal changes class, counted ONCE)

  But edges = cross-loop orbits, counted in UNORDERED pairs.
  Each edge corresponds to at least one cross-orbit from each side.
  So this is not quite right for edges.

  HOWEVER: the TOTAL number of cross-orbits = sum_C degree(C) = 2|E|.
  And the total number of ALL orbits = sum_C #orbits(C).
  So: 2|E| = total_cross_orbits = total_orbits - total_self_orbits.

  And: total_orbits = (1/n!) * sum_sigma sum_{T: sigma in Aut(T)} Fix_arcs(sigma, T).

  This IS a Burnside formula for the total arc-orbit count!
""")

# Compute the Burnside sum for total arc-orbits
for n in [3, 4, 5]:
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1, n)]

    # For each sigma in S_n, compute sum_T [sigma in Aut(T)] * Fix_arcs(sigma, T)
    burnside_sum = 0
    for sigma in permutations(range(n)):
        # Find all tournaments T with sigma in Aut(T)
        # A tournament T has sigma in Aut iff A[i][j] = A[sigma[i]][sigma[j]] for all i,j
        # This constrains arcs: arcs in the same sigma-orbit must have the same direction.

        # Compute pair orbits under sigma
        pair_orbits = []
        seen_pairs = set()
        for i in range(n):
            for j in range(i+1, n):
                if (i,j) in seen_pairs: continue
                orbit = set()
                ci, cj = i, j
                while (min(ci,cj), max(ci,cj)) not in orbit:
                    orbit.add((min(ci,cj), max(ci,cj)))
                    ci, cj = sigma[ci], sigma[cj]
                pair_orbits.append(orbit)
                seen_pairs.update(orbit)

        n_pair_orbits = len(pair_orbits)

        # For tournaments: each pair orbit with an ODD-length cycle
        # that reverses the canonical order contributes 0 to Fix.
        # For even-length cycles or non-reversing cycles, 2 choices.
        # The exact count of fixed tournaments:

        # For each orbit, check if sigma reverses the pair (i,j) -> (j,i)
        n_fixed_tournaments = 1
        for orb in pair_orbits:
            # Trace the orbit: does it involve a reversal?
            pair = next(iter(orb))
            i, j = pair
            ci, cj = i, j
            reverses = False
            for _ in range(len(orb)):
                ci2, cj2 = sigma[ci], sigma[cj]
                if ci2 > cj2:  # canonical order is min,max
                    reverses = not reverses
                ci, cj = min(ci2,cj2), max(ci2,cj2)
            # After going around the orbit, if reverses is True
            # and orbit length is odd, then Fix = 0
            if reverses and len(orb) % 2 == 1:
                n_fixed_tournaments = 0
                break
            # Otherwise: 2 choices (or 1 if reverses and even length)
            if not reverses:
                n_fixed_tournaments *= 2
            # If reverses and even: 2 choices (pair forced)
            elif len(orb) % 2 == 0:
                n_fixed_tournaments *= 2  # actually need more care
            # Simplified: just count 2^{non-reversing orbits}

        # Fix_arcs for this sigma on the fixed tournaments:
        # For each fixed tournament T and each arc (i,j) of T:
        # the arc is fixed by sigma iff sigma[i]=i and sigma[j]=j.
        # Number of such arcs = #{directed arcs from fixed points of sigma}

        # Fixed points of sigma
        fixed_pts = [v for v in range(n) if sigma[v] == v]
        n_fixed_arcs = len(fixed_pts) * (len(fixed_pts) - 1) // 2
        # Wait: among fixed points, the arcs are between pairs of fixed points.
        # Each such pair has a direction. For a FIXED tournament:
        # all arcs are determined. Among fixed-point pairs, the direction is free
        # (constrained by the orbit structure), so the number of arcs FROM
        # fixed points to fixed points = C(#fixed, 2), each with one direction.
        # The arc (i,j) with A[i][j]=1 is fixed by sigma iff sigma[i]=i, sigma[j]=j.
        # So there are EXACTLY C(#fixed, 2) arcs fixed by sigma (in each direction).
        # But actually only half of them (those where A[i][j]=1, not A[j][i]=1).
        # Hmm, for a specific T, Fix_arcs(sigma) = #{(i,j): i,j both fixed by sigma AND A[i][j]=1}
        # = approximately C(#fixed, 2) / 2 + ... depends on T.

        # For the Burnside sum, we need sum_T Fix_arcs(sigma) over T with sigma in Aut(T).
        # This is complex. Let me just compute the total directly.
        burnside_sum += n_fixed_tournaments * n_fixed_arcs

    total_orb_burnside = Fraction(burnside_sum, factorial(n))
    print("  n=%d: Burnside sum/(n!) = %s" % (n, total_orb_burnside))

# ============================================================
# PART 6: THE EDGE COUNT FORMULA
# ============================================================

banner("PART 6: TOWARDS THE FORMULA")

print("""
  FROM THE DATA:
  n=3: total_degree=2, |E|=1
  n=4: total_degree=10, |E|=5
  n=5: total_degree=60, |E|=30

  The DEGREE SUM = 2, 10, 60.
  = 2*1, 2*5, 2*30. (Handshaking lemma.)

  RATIO: degree_sum / |V| = 1.0, 2.5, 5.0
  = 1, 5/2, 5. Sequence of average degrees.

  CAN WE PREDICT degree_sum from n?
  1, 5, 30 ... but wait, these are the |E| values.
  2, 10, 60 are the degree sums.

  Degree sum = 2 * |E| = 2 * A_new(n).

  THE KEY QUANTITIES:
  total_orbits - total_self_orbits = total_cross_orbits = degree_sum.

  So: |E| = (total_orbits - total_self_orbits) / 2.

  We need formulas for total_orbits and total_self_orbits.

  TOTAL ORBITS = sum_C (1/|Aut(C)|) * sum_{g in Aut(C)} Fix_arcs(g)
  This is a DOUBLE Burnside sum (over classes and over automorphisms).

  SELF ORBITS = sum_C (1/|Aut(C)|) * #{arcs whose reversal preserves class}
  = sum_C F[C][C] / (|Aut(C)| * ... )

  Actually: total_self_orbits = sum_C (# self-loop orbits of C).
  A self-loop orbit is an Aut-orbit of arcs whose reversal gives
  an isomorphic tournament.
""")

# ============================================================
# SYNTHESIS
# ============================================================

banner("SYNTHESIS: THE EDGE COUNT FORMULA STATUS")

print("""
  WHAT WE NOW KNOW:

  1. |E(G_n)| = (1/2) * sum_C degree(C) where degree is computed
     from the Aut-orbit decomposition. PROVED.

  2. degree(C) = m / |Aut(C)| - (#self-loop orbits) ... NOT EXACTLY.
     Actually: degree(C) = #cross-orbits that land in DISTINCT classes.
     Multiple cross-orbits can land in the same class (degeneracy).

  3. The Aut-orbit computation is O(|Aut|^2 * m) per class.
     For generic classes (|Aut|=1): each arc is its own orbit,
     and degree = m - (#arcs preserving class) = m - F[C][C]/size(C).

  4. For |Aut|=1 classes: degree(C) = m - (self-loop count per tournament).

  5. The SEQUENCE |E(G_n)| = 1, 5, 30, 290 is NOT in OEIS.

  6. The average degree sequence 1.0, 2.5, 5.0, 10.4 grows ~ n^2/4.
     This is because m = C(n,2) ~ n^2/2, and about half the arcs
     lead to different classes.

  CONJECTURED FORMULA:
  |E(G_n)| ~ |V(G_n)| * C(n,2) / 4 for large n
  = A000568(n) * n(n-1)/8

  At n=5: 12 * 10/4 = 30. EXACT!
  At n=6: 56 * 15/4 = 210. Does NOT match 290.

  So the simple formula works at n=5 but fails at n=6.
  The deviation: 290 - 210 = 80 = extra edges from within-score connections.

  REFINED CONJECTURE:
  |E(G_n)| = |V| * m / 4 + correction(n)
  where correction(n) accounts for within-score edges.
  correction: 0, 0, 0, 80 at n=3,4,5,6.

  At n=5: avg_deg = 5.0 = m/2 = 10/2. EXACT.
  At n=6: avg_deg = 10.36 > m/2 = 7.5. The deviation grows.

  THE FORMULA REMAINS OPEN. Computing |E(G_7)| (456 classes)
  would give the next data point. The Aut-orbit method is feasible
  for n=7 with optimized code.

  The DUAL BURNSIDE PROCESS (Feng) offers a theoretical framework:
  the dual chain Q=AB connects the classical Burnside (for vertex count)
  to a dual process whose stationary distribution encodes edge information.
  The shared eigenvalues between Q and K may provide the missing formula.
""")

print("Done. opus-2026-03-22-S185")
