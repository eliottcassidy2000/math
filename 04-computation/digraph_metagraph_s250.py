#!/usr/bin/env python3
"""
digraph_metagraph_s250.py -- opus-2026-03-23-S250

EXTENDING THE MERGED META-GRAPH TO GENERAL DIRECTED GRAPHS

The tournament meta-graph G_n = Q_m / S_n lives on the m-cube because
every pair of vertices has exactly ONE arc (tournament = complete digraph).
What happens when we relax this?

THREE GENERALIZATIONS:

1. ORIENTED GRAPHS: each pair has at most one arc (0 or 1 direction).
   The "flip" operation: add, remove, or reverse an arc.
   State space: {0, 1, 2}^{C(n,2)} = 3^m (three states per pair:
   no arc, forward, backward). The meta-graph lives on the 3-ary cube!

2. DIGRAPHS (allowing both directions): each pair has 0, 1, or 2 arcs.
   State space: {00, 01, 10, 11}^{C(n,2)} = 4^m = 2^{2m} = 2^{n(n-1)}.
   This is the full directed graph space.

3. DIGRAPHS WITH SELF-LOOPS: each vertex can also have a self-loop.
   State space: 2^{n^2} (one bit per (i,j) entry, including i=i).

The TOURNAMENT case is the SIMPLEST: 2^m = 2^{C(n,2)}.
It's the Boolean hypercube because each pair has exactly 2 choices.

THE KEY QUESTION: Which of our results survive generalization?

Author: opus-2026-03-23-S250
"""

import sys
from math import comb, factorial
from itertools import permutations, product
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

# ============================================================
banner("GENERALIZATION 1: ORIENTED GRAPHS")
# ============================================================

print("""
  An ORIENTED GRAPH on n vertices: each pair {i,j} has:
  - No arc (state 0)
  - Arc i→j (state 1)
  - Arc j→i (state 2)

  Total oriented graphs: 3^{C(n,2)} = 3^m.

  The FLIP GRAPH: two oriented graphs are adjacent if they differ
  at exactly one pair. But now each pair has 3 states, so each
  vertex has 2m neighbors (2 choices for each of m pairs).

  This is the HAMMING GRAPH H(m, 3) — the m-dimensional 3-ary cube.

  The isomorphism class graph: G_n^{orient} = H(m, 3) / S_n.

  WHAT SURVIVES from tournaments:
  - Burnside formula for V = #{iso classes}: same structure, different Fix
  - The S_n action on pairs is the SAME
  - The "staircase" becomes a staircase with 3-valued cells (ternary tiling)

  WHAT'S NEW:
  - The state "no arc" is a NEUTRAL state. It's the 0 in the 3-valued encoding.
  - The complement map: swap state 1↔2, leave 0 fixed. Still order 2.
  - But now there's also a "deletion" operation: map to state 0.
  - The 3 states form the GROUP Z_3? No — {0, 1, 2} under flip is Z_3 only
    if we define flip cyclically. Under arc reversal: {1, 2} swap, 0 is fixed.

  THE 2 AND 3 CONNECTION:
  - In tournaments: m cells, each binary (2 states) → 2^m configurations
  - In oriented graphs: m cells, each TERNARY (3 states) → 3^m configurations
  - The tournament theory uses the ORDER-2 generator (complement)
  - The oriented graph theory ADDS the ORDER-3 generator (ternary cycle on states)!

  THIS IS THE LITERAL "2 AND 3" OF PSL(2,Z):
  - Tournament = the 2-structure (binary)
  - Oriented graph = the 3-structure (ternary)
  - The modular group PSL(2,Z) = <S, ST> with orders 2 and 3
    governs the interplay between these two.
""")

# Compute oriented graph iso classes at small n
for n in [3, 4]:
    m = comb(n, 2)
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

    # Each pair has 3 states: 0 (no arc), 1 (i→j), 2 (j→i)
    # Total: 3^m oriented graphs

    def canon_orient(state):
        """Canonical form of an oriented graph state tuple under S_n."""
        best = None
        for perm in permutations(range(n)):
            new_state = []
            for (i,j) in pairs:
                pi, pj = perm[i], perm[j]
                if pi < pj:
                    idx = pairs.index((pi, pj))
                    new_state.append(state[idx])
                else:
                    idx = pairs.index((pj, pi))
                    s = state[idx]
                    new_state.append(2 if s == 1 else (1 if s == 2 else 0))
            new_state = tuple(new_state)
            if best is None or new_state < best:
                best = new_state
        return best

    classes = set()
    for state in product(range(3), repeat=m):
        c = canon_orient(state)
        classes.add(c)

    nc = len(classes)
    print("  n=%d: 3^%d = %d oriented graphs, %d iso classes" % (n, m, 3**m, nc))

    # Compare with tournaments
    n_tourn = 2**m
    tc = set()
    for bits in range(n_tourn):
        state = tuple(1 if (bits >> k) & 1 else 2 for k in range(m))
        tc.add(canon_orient(state))
    print("       %d tournaments, %d tournament iso classes" % (n_tourn, len(tc)))
    print("       Ratio: %d / %d = %.2f (oriented/tournament classes)" %
          (nc, len(tc), nc / len(tc)))

# ============================================================
banner("GENERALIZATION 2: GENERAL DIGRAPHS")
# ============================================================

print("""
  A GENERAL DIGRAPH on n vertices: each ordered pair (i,j) with i≠j has
  independently an arc or not. Total: 2^{n(n-1)}.

  This is MUCH larger than tournaments (2^{C(n,2)}).
  At n=4: 2^12 = 4096 digraphs vs 2^6 = 64 tournaments.

  The flip operation: toggle one arc (add or remove).
  The meta-graph: H(n(n-1), 2) = the n(n-1)-dimensional hypercube,
  quotiented by S_n.

  WHAT SURVIVES:
  - Burnside counting for iso classes (A000273 for digraphs)
  - The S_n action on the n(n-1) arc positions
  - The complement map (toggle ALL arcs)

  WHAT'S NEW:
  - S_n acts on ORDERED pairs (i,j), not unordered {i,j}.
  - The "complement" of arc (i→j) is "no arc from i to j",
    not "arc (j→i)". So complement in digraphs ≠ reversal.
  - We can define MULTIPLE merge operations:
    a) Merge by COMPLEMENT (toggle all)
    b) Merge by REVERSAL (reverse all arcs) — gives the "transpose" digraph
    c) Merge by BOTH (complement + reversal)

  THE DIRECTED-GRAPH STAIRCASE:
  For tournaments: the staircase δ_{n-2} has C(n,2) cells.
  For digraphs: the matrix has n(n-1) cells (the full off-diagonal).
  The staircase becomes the FULL SQUARE minus the diagonal.

  The tournament staircase is the UPPER TRIANGLE of this square.
  The tournament constraint "exactly one arc per pair" means:
  cell (i,j) and cell (j,i) are COMPLEMENTARY (sum to 1).
  In digraphs, they're INDEPENDENT.

  So: DIGRAPH THEORY = tournament theory WITHOUT the complementarity constraint.
  The tournament is a SECTION of the digraph space:
  digraphs ∩ {A[i][j] + A[j][i] = 1 for all i<j} = tournaments.
""")

# Count digraph iso classes at n=3
n = 3
total_digraphs = 2**(n*(n-1))
print("  Digraph iso classes at n=%d:" % n)
print("  Total digraphs: 2^%d = %d" % (n*(n-1), total_digraphs))

# Enumerate and canonicalize
ordered_pairs = [(i,j) for i in range(n) for j in range(n) if i != j]
m_dig = len(ordered_pairs)

def canon_digraph(bits, n, ordered_pairs):
    best = None
    for perm in permutations(range(n)):
        new_bits = 0
        for k, (i,j) in enumerate(ordered_pairs):
            pi, pj = perm[i], perm[j]
            new_k = ordered_pairs.index((pi, pj))
            if bits & (1 << k):
                new_bits |= (1 << new_k)
        if best is None or new_bits < best:
            best = new_bits
    return best

dig_classes = set()
for bits in range(total_digraphs):
    dig_classes.add(canon_digraph(bits, n, ordered_pairs))

print("  Digraph iso classes: %d" % len(dig_classes))
# OEIS A000273: 1, 1, 3, 16, 218, 9608, ...
# At n=3: should be 16
print("  Expected (A000273): 16")

# Now build the digraph flip graph and count edges
edges = set()
for bits in range(total_digraphs):
    c1 = canon_digraph(bits, n, ordered_pairs)
    for k in range(m_dig):
        flipped = bits ^ (1 << k)
        c2 = canon_digraph(flipped, n, ordered_pairs)
        if c1 != c2:
            edges.add((min(c1, c2), max(c1, c2)))

n_dig_classes = len(dig_classes)
n_dig_edges = len(edges)
print("  Digraph meta-graph: V=%d, E=%d" % (n_dig_classes, n_dig_edges))
print("  Average degree: %.2f (m_dig = %d)" % (2*n_dig_edges/n_dig_classes, m_dig))

# Compare with tournament meta-graph
print("\n  Tournament meta-graph at n=3: V=2, E=1")
print("  Ratio: digraph V/tournament V = %d/2 = %.1f" % (n_dig_classes, n_dig_classes/2))
print("  Ratio: digraph E/tournament E = %d/1 = %.1f" % (n_dig_edges, n_dig_edges))

# ============================================================
banner("THE THREE LEVELS OF GENERALIZATION")
# ============================================================

print("""
  LEVEL 0: TOURNAMENTS (complete, antisymmetric)
  ================================================
  State space: {0,1}^{C(n,2)} = the m-cube Q_m
  Meta-graph: G_n = Q_m / S_n
  Merge: Z_2 (complement = reverse all arcs)
  Merged: G_n / Z_2

  Properties:
  - V = A000568(n) (tournament iso classes)
  - edge_orbits = T_n/2 + (n-2)!
  - H is always ODD
  - β_2 = 0 for all tournaments

  LEVEL 1: ORIENTED GRAPHS (at most one arc per pair)
  ===================================================
  State space: {0,1,2}^{C(n,2)} = H(m, 3) (3-ary cube)
  Meta-graph: G_n^{orient} = H(m, 3) / S_n
  Merge: Z_2 (reverse all arcs, 0 stays fixed)
  Merged: G_n^{orient} / Z_2

  Properties:
  - V = A001174(n) (oriented graph iso classes: 1, 2, 7, 42, 582, ...)
  - Tournaments embed as the "complete" layer
  - The "empty graph" (all 0s) is a single fixed point
  - The staircase cells become TERNARY (0, 1, 2)

  WHAT CHANGES:
  - The Hamiltonian path count H generalizes to "path count in digraph"
  - Arcs can be ABSENT, not just reversed
  - The "apex" arc can be absent, forward, or backward (3 states)
  - The recursive decomposition has 3 apex states instead of 2

  LEVEL 2: GENERAL DIGRAPHS (any arcs, including mutual)
  ======================================================
  State space: {0,1}^{n(n-1)} = Q_{n(n-1)} (the big cube)
  Meta-graph: G_n^{dig} = Q_{n(n-1)} / S_n
  Merge: multiple options (complement, transpose, both)
  Merged: G_n^{dig} / (S_n × merge_group)

  Properties:
  - V = A000273(n) (digraph iso classes: 1, 1, 3, 16, 218, 9608, ...)
  - Tournaments embed as a SECTION (constraint: A[i][j]+A[j][i]=1)
  - Mutual arcs allowed (A[i][j]=A[j][i]=1)
  - The staircase becomes the FULL SQUARE

  WHAT CHANGES:
  - The staircase is no longer a TRIANGLE but a SQUARE
  - The complement map changes meaning (delete arcs vs reverse)
  - TWO merge operations: complement AND transpose
  - The Schur-Weyl decomposition involves more irreps
  - The recursive structure: n(n-1) arcs, not just C(n,2)

  LEVEL 3: DIGRAPHS WITH SELF-LOOPS (the most general)
  =====================================================
  State space: {0,1}^{n^2} = Q_{n^2}
  Meta-graph: G_n^{full} = Q_{n^2} / S_n
  V = A000595(n) (binary relations: 1, 2, 10, 104, 3044, ...)

  This is the meta-graph of ALL BINARY RELATIONS on [n].
""")

# ============================================================
banner("WHAT SURVIVES AT EACH LEVEL")
# ============================================================

print("""
  FORMULA / PROPERTY          | TOURN | ORIENT | DIGRAPH | FULL
  ----------------------------|-------|--------|---------|-----
  V = Burnside formula        |  YES  |  YES   |  YES    | YES
  T_n = arc-orbit count       |  YES  |  YES   |  YES    | YES
  edge_orbits = T_n/2 + corr  |  ???  |  ???   |  ???    | ???
  Case A = T_n/2              |  YES  |  YES?  |  YES?   | YES?
  Case B = (n-2)!             |  ???  |  ???   |  ???    | ???
  H always odd                |  YES  |  NO    |  NO     | NO
  Clique complex beta_1       |  2,15 |  ???   |  ???    | ???
  Staircase = triangle        |  YES  |  YES   |  NO     | NO
  3-color decomposition       |  YES  |  ???   |  ???    | ???
  f(n) = arc neutrality       |  YES  |  ???   |  ???    | ???

  THE HIERARCHY:
  FULL ⊃ DIGRAPH ⊃ ORIENTED ⊃ TOURNAMENT

  Each level adds constraints:
  FULL → DIGRAPH: no self-loops (A[i][i] = 0)
  DIGRAPH → ORIENTED: at most one arc per pair (A[i][j]*A[j][i] = 0)
  ORIENTED → TOURNAMENT: exactly one arc per pair (A[i][j]+A[j][i] = 1)
""")

# ============================================================
banner("THE 2 AND 3 AT EACH LEVEL")
# ============================================================

print("""
  THE MODULAR GROUP STRUCTURE:

  At LEVEL 0 (tournaments):
  - Each cell has 2 states → the 2-cube per cell
  - The complement swaps states (order 2)
  - The near-automorphism involves 3-cycles (order 3)
  - PSL(2,Z) = <S(2), ST(3)> governs the structure

  At LEVEL 1 (oriented graphs):
  - Each cell has 3 states → the 3-ary cube per cell
  - The reversal swaps states 1↔2, fixes 0 (order 2)
  - The cyclic rotation 0→1→2→0 has order 3
  - PSL(2,Z) appears LITERALLY: the 2 and 3 ARE the cell states!
  - S = reversal (order 2), ST = cyclic (order 3)
  - The modular group IS the symmetry of the 3-state cell!

  At LEVEL 2 (digraphs):
  - Each ordered pair has 2 states → back to binary
  - But pairs (i,j) and (j,i) are INDEPENDENT
  - The transpose swaps (i,j) ↔ (j,i): order 2
  - The complement toggles each arc: order 2
  - Together: Klein 4-group Z_2 × Z_2 (not PSL(2,Z))

  THE DEEP INSIGHT:
  Oriented graphs are the NATURAL HOME of PSL(2,Z) because
  the 3 states per cell match the 3 orbifold points:
  - State 0 (no arc) = the cusp (transitive limit)
  - State 1 (forward) = the i-point (complement)
  - State 2 (backward) = the ω-point (cycle)

  Tournaments = the oriented graphs with NO state-0 cells.
  They live in the "non-cuspidal" part of the modular curve.
  The cusp is excluded, forcing every cell to be either 1 or 2.

  This is why tournament theory is governed by orders 2 and 3 ONLY:
  the cusp (order ∞) is absent. Removing the cusp collapses PSL(2,Z)
  to the finite group <S, ST | S^2 = (ST)^3 = 1> = the modular group
  modulo the cusp stabilizer.
""")

# ============================================================
banner("CONCRETE COMPUTATION: ORIENTED GRAPH META-GRAPH")
# ============================================================

# At n=3: compute the oriented graph meta-graph
n = 3; m = comb(n, 2)
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

print("  ORIENTED GRAPH META-GRAPH at n=3:\n")

# Enumerate all 3^3 = 27 oriented graphs
orient_classes = {}
for state in product(range(3), repeat=m):
    c = canon_orient(state)
    if c not in orient_classes:
        orient_classes[c] = len(orient_classes)

n_oc = len(orient_classes)
print("  Iso classes: %d" % n_oc)

# Build adjacency (flip one pair's state to an adjacent state)
oc_map = {}
for state in product(range(3), repeat=m):
    oc_map[state] = orient_classes[canon_orient(state)]

orient_edges = set()
for state in product(range(3), repeat=m):
    c1 = oc_map[state]
    for k in range(m):
        for new_val in range(3):
            if new_val == state[k]: continue
            new_state = list(state)
            new_state[k] = new_val
            new_state = tuple(new_state)
            c2 = oc_map[new_state]
            if c1 != c2:
                orient_edges.add((min(c1,c2), max(c1,c2)))

n_oe = len(orient_edges)
print("  Edges: %d" % n_oe)
print("  Average degree: %.2f" % (2*n_oe/n_oc))

# List the classes with their properties
print("\n  Classes:")
for c_key, cid in sorted(orient_classes.items(), key=lambda x: x[1]):
    state = c_key
    n_arcs = sum(1 for s in state if s > 0)
    is_tournament = all(s > 0 for s in state)
    # Count labeled copies
    copies = sum(1 for s in product(range(3), repeat=m) if canon_orient(s) == c_key)
    print("    %d: state=%s, arcs=%d, is_tournament=%s, copies=%d" %
          (cid, state, n_arcs, is_tournament, copies))

# How many are tournaments?
tourn_classes = [c for c, cid in orient_classes.items() if all(s > 0 for s in c)]
print("\n  Tournament classes among oriented: %d out of %d" %
      (len(tourn_classes), n_oc))

# ============================================================
banner("THE EDGE ORBIT FORMULA: DOES IT GENERALIZE?")
# ============================================================

print("""
  For TOURNAMENTS: edge_orbits = T_n/2 + (n-2)!

  For ORIENTED GRAPHS: each pair has 3 states, each "flip" changes
  one pair's state. An edge orbit of the Hamming graph H(m,3) under S_n:

  edge_orbits^{orient} = (1/n!) * sum_σ Fix_orient(σ)

  where Fix_orient(σ) counts edges of H(m,3) fixed by σ.

  An edge (state, state') differs in exactly one pair's value.
  σ fixes this edge iff σ maps (state, state') to itself (both fixed)
  or swaps them.

  The SAME Case A / Case B decomposition applies:
  Case A: σ fixes both endpoints
  Case B: σ swaps them

  For Case A: Fix_orient_A(σ) = Fix_OG(σ) × C(fix(σ), 2) × 2
  (the ×2 because each fixed pair has 2 possible flips: s→s', s→s'')

  For Case B: analogous to tournaments but with 3-state structure.

  CONJECTURE: edge_orbits^{orient} = T_n^{orient}/2 + correction(n)
  where T_n^{orient} counts arc-orbits in oriented graph space
  and the correction is related to (n-2)! × (some factor involving 3).

  THE 3 ENTERS: in tournaments, (n-2)! came from pair-swapping.
  In oriented graphs, the pair-swap has RICHER structure because
  the "empty" state introduces a new symmetry class.
  The correction might be (n-2)! × 2 or (n-2)! × 3.
""")

# Actually compute for n=3
print("  Direct computation at n=3:\n")

# T_n^{orient} = arc-orbits in oriented graph space
# = #{S_n-orbits on (oriented_graph, pair, new_state) triples}
# Wait: an "arc" in the Hamming graph is a pair (state, k, new_val).
# The S_n orbit of such a triple gives an arc-orbit.

# Let me count edge_orbits directly
# From above: n_oe = edges in the meta-graph (distinct class pairs)
# But we also need self-loops and multi-edges

total_orient_edges_labeled = 0
for state in product(range(3), repeat=m):
    for k in range(m):
        for new_val in range(3):
            if new_val == state[k]: continue
            new_state = list(state); new_state[k] = new_val; new_state = tuple(new_state)
            if tuple(sorted([state, new_state])) not in set():  # count all
                total_orient_edges_labeled += 1

# Each labeled edge counted twice (once per direction)
total_orient_edges_labeled //= 2

# Orient T_n = arc-orbits = (1/n!) * sum_σ Fix_{(OG, arc)}(σ)
# Simple: T_n = total_labeled_arcs / n! ... no, Burnside needed.

# For n=3: S_3 has 6 elements
fix_arc_sum = 0
for perm in permutations(range(n)):
    for state in product(range(3), repeat=m):
        # Check if perm fixes this oriented graph
        new_state = []
        for (i,j) in pairs:
            pi, pj = perm[i], perm[j]
            if pi < pj:
                idx = pairs.index((pi, pj))
                new_state.append(state[idx])
            else:
                idx = pairs.index((pj, pi))
                s = state[idx]
                new_state.append(2 if s == 1 else (1 if s == 2 else 0))
        if tuple(new_state) == state:
            # perm fixes this OG. Count how many arcs are fixed.
            for k in range(m):
                i, j = pairs[k]
                pi, pj = perm[i], perm[j]
                if pi < pj and pairs.index((pi,pj)) == k:
                    # This pair is mapped to itself
                    fix_arc_sum += 2  # 2 flips per fixed pair
                # Also check if pair is mapped to itself with reversal
                # More careful: a "flip" at pair k changes state[k].
                # For σ to fix this flip: σ must fix the pair {i,j}
                # AND the flip type (which new_val) must be preserved.

# This is getting complex. Let me just count edge_orbits directly via Burnside.
fix_edge_sum = 0
for perm in permutations(range(n)):
    def apply_orient(state):
        new_state = []
        for (i,j) in pairs:
            pi, pj = perm[i], perm[j]
            if pi < pj:
                idx = pairs.index((pi, pj))
                new_state.append(state[idx])
            else:
                idx = pairs.index((pj, pi))
                s = state[idx]
                new_state.append(2 if s == 1 else (1 if s == 2 else 0))
        return tuple(new_state)

    for state in product(range(3), repeat=m):
        for k in range(m):
            for new_val in range(3):
                if new_val == state[k]: continue
                new_state = list(state); new_state[k] = new_val
                new_state = tuple(new_state)
                if new_state > state:  # unordered edge, count once
                    ps = apply_orient(state)
                    pn = apply_orient(new_state)
                    if (ps == state and pn == new_state) or (ps == new_state and pn == state):
                        fix_edge_sum += 1

orient_edge_orbits = fix_edge_sum // factorial(n)
print("  n=3: orient edge_orbits = %d" % orient_edge_orbits)
print("  n=3: orient E (meta-graph edges) = %d" % n_oe)
print("  n=3: orient gap_orbits = %d" % (orient_edge_orbits - n_oe))

# Compare with tournament:
print("\n  Tournament at n=3: edge_orbits = 3, E = 1, gap = 2")
print("  Oriented at n=3: edge_orbits = %d, E = %d, gap = %d" %
      (orient_edge_orbits, n_oe, orient_edge_orbits - n_oe))

print("\nDone. opus-2026-03-23-S250")
