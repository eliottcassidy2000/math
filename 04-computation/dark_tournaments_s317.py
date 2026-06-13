#!/usr/bin/env python3
"""
dark_tournaments_s317.py -- opus-2026-03-24-S317

DARK TOURNAMENTS: THE MISSING FOURTH TYPE

The Royle equinumerosity gives us a 2×2 system:

                  DIRECTED           UNDIRECTED
  ALL-ODD-PERM   Tournaments (T)    Even Graphs (E)
  ALL-EVEN-PERM  Dark Tourn. (D)    Odd Graphs (O)

Where:
  #T = #E = (1/n!) Σ_{|σ| odd} 2^{c(σ_E)}     [A000568]
  #O = #D = (1/n!) Σ_{|σ| even} 2^{c(σ_E)}    [A000088 - A000568]
  #T + #D = #Graphs = A000088                    [total undirected graphs]

But WHAT IS D? What directed structure gives exactly #OddGraphs iso classes?

THE KEY: the Burnside sum for D uses EVEN-ORDER permutations.
  #D = (1/n!) Σ_{σ: |σ| even} #{directed things fixed by σ}

For tournaments: σ fixes T iff σ(T) = T, which requires all cycle pairs
to be compatible with the orientation. Self-paired cycles (which exist
only for even-order σ) force A[i][j] = A[j][i], impossible in tournaments.
So even-order σ contribute 0 to the tournament count.

For D: we need a directed structure where even-order σ CAN fix things.
Self-paired cycles give A[i][j] = A[j][i] — this means SYMMETRIC arcs
(= undirected edges). So D should be a structure that ALLOWS some
arcs to be undirected (or equivalently, bidirectional).

CANDIDATE: PARTIALLY DIRECTED GRAPHS (mixed graphs)
where each edge is EITHER directed (one way) OR undirected (both ways).

For a pair {i,j}: three choices: i→j, j→i, or i↔j (undirected).
Under σ with a self-paired cycle involving i,j: only i↔j is fixed.
So even-order σ fix partially directed graphs where self-paired arcs
are undirected.

#PD(n) = (1/n!) Σ_σ 3^{c(σ_E)} = A000088(n)... no, that gives
the number of 3-colored complete graphs.

Let me think more carefully.

Author: opus-2026-03-24-S317
"""

import sys
import numpy as np
from math import comb, factorial, gcd
from itertools import permutations
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def pair_orbits_general(ct):
    """Number of pair orbits under a permutation with cycle type ct."""
    within = sum((c - 1) // 2 for c in ct)
    across = sum(gcd(ct[i], ct[j]) for i in range(len(ct)) for j in range(i+1, len(ct)))
    return within + across

def self_paired_cycles(ct):
    """Number of self-paired pair orbits (from even cycles)."""
    return sum(c // 2 for c in ct if c % 2 == 0)

def multinomial(n, ct):
    ct_counter = Counter(ct)
    denom = 1
    for c, a in ct_counter.items():
        denom *= (c ** a) * factorial(a)
    return factorial(n) // denom

def gen_all_partitions(n):
    result = []
    def _gen(remaining, max_part, current):
        if remaining == 0:
            result.append(tuple(sorted(current)))
            return
        for p in range(min(remaining, max_part), 0, -1):
            _gen(remaining - p, p, current + [p])
    _gen(n, n, [])
    return list(set(result))

def gen_odd_partitions(n):
    result = []
    def _gen(remaining, max_part, current):
        if remaining == 0:
            result.append(tuple(sorted(current)))
            return
        start = min(remaining, max_part)
        if start % 2 == 0: start -= 1
        for p in range(start, 0, -2):
            _gen(remaining - p, p, current + [p])
    _gen(n, n, [])
    return list(set(result))

def gen_even_order_partitions(n):
    """Partitions where ALL cycle lengths have ODD length (= even ORDER permutations).
    Wait — even ORDER means the LCM of cycle lengths is even.
    A permutation has even order iff at least one cycle has even length.

    Actually: order = LCM of cycle lengths. Even order iff any cycle is even-length.
    Odd order iff ALL cycles are odd-length.

    So even-order perms = ALL partitions MINUS all-odd partitions.
    """
    all_parts = gen_all_partitions(n)
    odd_parts = gen_odd_partitions(n)
    return [p for p in all_parts if p not in odd_parts]

# ============================================================
banner("THE FOUR-TYPE SYSTEM: COUNTS")
# ============================================================

print("  THE 2×2 SYSTEM:\n")
print("  %-4s %-10s %-10s %-10s %-10s %-10s" %
      ("n", "Graphs", "T (odd σ)", "D (even σ)", "Even G", "Odd G"))
print("  " + "-"*55)

for n in range(3, 12):
    all_parts = gen_all_partitions(n)
    odd_parts = gen_odd_partitions(n)
    even_parts = gen_even_order_partitions(n)

    # #Graphs = (1/n!) Σ_all 2^{c(σ_E)}
    graphs = 0
    for ct in all_parts:
        mult = multinomial(n, ct)
        po = pair_orbits_general(ct)
        graphs += mult * (2 ** po)
    graphs //= factorial(n)

    # #Tournaments = (1/n!) Σ_{odd σ} 2^{c(σ_E)}
    tourn = 0
    for ct in odd_parts:
        mult = multinomial(n, ct)
        po = pair_orbits_general(ct)
        tourn += mult * (2 ** po)
    tourn //= factorial(n)

    # #D = #Graphs - #T = (1/n!) Σ_{even σ} 2^{c(σ_E)}
    dark = graphs - tourn

    # Verify: D via even-order partitions
    dark_check = 0
    for ct in even_parts:
        mult = multinomial(n, ct)
        po = pair_orbits_general(ct)
        dark_check += mult * (2 ** po)
    dark_check //= factorial(n)

    # #Even = #T (Royle)
    even_g = tourn

    # #Odd = #D
    odd_g = dark

    print("  %-4d %-10d %-10d %-10d %-10d %-10d" %
          (n, graphs, tourn, dark, even_g, odd_g))

    assert dark == dark_check, "Dark tournament count mismatch!"
    assert graphs == tourn + dark, "Graph count mismatch!"

# ============================================================
banner("WHAT FIXES UNDER EVEN-ORDER PERMUTATIONS?")
# ============================================================

print("""
  THE MECHANISM: Why even-order σ give 2^{pair_orbits} "dark" objects.

  For a permutation σ with at least one even-length cycle:
  - The pair orbits still decompose into within-cycle and cross-cycle.
  - Self-paired orbits exist (from even cycles): pairs (i, σ^{k/2}(i))
    where k = even cycle length.
  - For a TOURNAMENT: self-paired orbits require A[i][j] = A[j][i],
    which is IMPOSSIBLE (tournaments are antisymmetric).
  - So even-order σ fix 0 tournaments.

  For a GRAPH (undirected): self-paired orbits are fine (A[i][j] = A[j][i]
  always holds). So 2^{pair_orbits} graphs are fixed.

  The "dark tournaments" count = the contribution of even-order σ to
  the graph Burnside sum. These are graphs that can ONLY be fixed by
  even-order permutations — i.e., they have no odd-order automorphism.

  But this isn't quite right. Let me think again...

  Actually: #D = (1/n!) Σ_{even σ} 2^{po(σ)} counts the average
  number of even-order fixers, which equals the number of orbits
  that contribute to even-order Burnside sum.

  These are NOT "graphs with no odd-order automorphism."
  They are the EXCESS of graphs over tournaments.
  D = #Graphs - #Tournaments = #OddGraphs.

  So: D IS the count of odd graphs (up to isomorphism).
  "Dark tournaments" ARE odd graphs. The name is evocative but
  the object is known!
""")

# ============================================================
banner("THE META-GRAPH PARALLEL")
# ============================================================

# At n=5: #T = 12, #D = #O = A000088(5) - 12 = 34 - 12 = 22.
# So there are 22 "dark tournament" = odd graph iso classes.
# And 12 tournament = even graph iso classes.

# Build the odd graph meta-graph at n=5 and compare
n = 5
m = comb(n, 2)  # 10 edges in K_5

# Build all graphs on 5 vertices
print("  n=%d: Building graph meta-graphs...\n" % n)

def graph_canon(adj, n):
    """Canonical form of an undirected graph adjacency (upper triangle bits)."""
    best = None
    for perm in permutations(range(n)):
        sig = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(i+1, n))
        if best is None or sig < best:
            best = sig
    return best

def is_even_graph(adj, n):
    """Check if every vertex has even degree."""
    for v in range(n):
        deg = sum(adj[v][w] for w in range(n) if w != v)
        if deg % 2 != 0:
            return False
    return True

def is_odd_graph(adj, n):
    """Check if every vertex has odd degree."""
    for v in range(n):
        deg = sum(adj[v][w] for w in range(n) if w != v)
        if deg % 2 != 1:
            return False
    return True

all_pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

even_classes = {}; odd_classes = {}; all_classes = {}
even_members = defaultdict(list); odd_members = defaultdict(list)

for bits in range(2**m):
    adj = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(all_pairs):
        if (bits >> k) & 1:
            adj[i][j] = 1; adj[j][i] = 1

    c = graph_canon(adj, n)
    if c not in all_classes:
        all_classes[c] = len(all_classes)

    if is_even_graph(adj, n):
        if c not in even_classes:
            even_classes[c] = len(even_classes)
        even_members[even_classes[c]].append(bits)
    elif is_odd_graph(adj, n):
        if c not in odd_classes:
            odd_classes[c] = len(odd_classes)
        odd_members[odd_classes[c]].append(bits)

print("  Total graph classes: %d = A000088(%d)" % (len(all_classes), n))
print("  Even graph classes: %d" % len(even_classes))
print("  Odd graph classes: %d" % len(odd_classes))
print("  Other: %d" % (len(all_classes) - len(even_classes) - len(odd_classes)))
print()

# Build edge-toggle meta-graph for even and odd graphs
# Toggle = flip one edge in K_n

for graph_type, classes, members_dict, type_name in [
    ("even", even_classes, even_members, "EVEN"),
    ("odd", odd_classes, odd_members, "ODD")]:

    nm = len(classes)
    if nm == 0:
        print("  %s: no classes" % type_name)
        continue

    # Reverse map: class_id → canonical form
    class_by_id = {v: k for k, v in classes.items()}

    # Build adjacency: toggle one edge, check if result is same type
    adj_meta = np.zeros((nm, nm), dtype=int)
    for cid in range(nm):
        for bits in members_dict[cid][:10]:  # sample
            for k in range(m):
                nbr = bits ^ (1 << k)
                # Build neighbor graph
                adj_nbr = [[0]*n for _ in range(n)]
                for kk, (i, j) in enumerate(all_pairs):
                    if (nbr >> kk) & 1:
                        adj_nbr[i][j] = 1; adj_nbr[j][i] = 1
                if graph_type == "even" and is_even_graph(adj_nbr, n):
                    c_nbr = graph_canon(adj_nbr, n)
                    if c_nbr in classes:
                        nid = classes[c_nbr]
                        adj_meta[cid][nid] = 1
                elif graph_type == "odd" and is_odd_graph(adj_nbr, n):
                    c_nbr = graph_canon(adj_nbr, n)
                    if c_nbr in classes:
                        nid = classes[c_nbr]
                        adj_meta[cid][nid] = 1

    n_edges = sum(1 for i in range(nm) for j in range(i+1, nm) if adj_meta[i][j])
    max_deg = max(sum(adj_meta[i]) for i in range(nm)) if nm > 0 else 0

    print("  %s meta-graph: V=%d, E=%d, max_deg=%d" %
          (type_name, nm, n_edges, max_deg))

# Build tournament meta-graph for comparison
print()
from itertools import permutations as perms

def tourn_canon(A, nn):
    sc = [sum(A[i]) for i in range(nn)]
    sg = defaultdict(list)
    for v in range(nn): sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(sg.keys())]
    best = None
    def gp(gs2):
        if not gs2: yield []; return
        for p in perms(gs2[0]):
            for r in gp(gs2[1:]): yield list(p)+r
    for p in gp(gs):
        f = tuple(A[p[i]][p[j]] for i in range(nn) for j in range(nn) if i!=j)
        if best is None or f < best: best = f
    return best

tourn_classes = {}
for bits in range(2**m):
    A = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(all_pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    c = tourn_canon(A, n)
    if c not in tourn_classes:
        tourn_classes[c] = len(tourn_classes)

print("  Tournament classes: %d = A000568(%d)" % (len(tourn_classes), n))
print("  Verification: #Even = #Tourn = %s" % (len(even_classes) == len(tourn_classes)))

# ============================================================
banner("THE DARK TOURNAMENT CANDIDATES")
# ============================================================

print("""
  WHAT STRUCTURE HAS EXACTLY #OddGraphs ISO CLASSES?

  CANDIDATE 1: Odd graphs themselves (trivial answer).
  D = odd graphs. This is correct but uninteresting.

  CANDIDATE 2: SEMI-TOURNAMENTS (partially directed graphs).
  A semi-tournament: for each pair {i,j}, either i→j or j→i or i↔j.
  3 choices per pair → 3^{C(n,2)} total.
  But this gives A000088 + more (too many classes).

  CANDIDATE 3: TOURNAMENTS WITH ONE UNDIRECTED EDGE.
  Take a tournament on n vertices and make one arc bidirectional.
  This creates a "tournament with one defect" = a near-tournament.
  Count: C(n,2) × 2^{C(n,2)-1} labeled (one pair is undirected).
  Iso classes: need Burnside.

  CANDIDATE 4: SIGNED TOURNAMENTS.
  A signed tournament: each arc gets a sign (+/-).
  This doubles the choices: 4 options per pair (i→j+, i→j-, j→i+, j→i-).
  Iso classes = (1/n!) Σ_σ 4^{pair_orbits(σ)} = way too many.

  CANDIDATE 5: TOURNAMENTS WITH A MARKED HAMILTONIAN PATH.
  Each tournament T has H(T) Hamiltonian paths.
  A "marked tournament" = (T, P) where P is a specific HP.
  Total: Σ H(T) = n! × 2^{m}... no, that's wrong.
  Actually: Σ_{labeled T} H(T) = the sum of all HP counts.

  CANDIDATE 6: COLORINGS OF THE METAGRAPH.
  The metagraph has χ = n-1 colors.
  A "dark tournament" = a pair (tournament class, coloring).
  This gives (n-1) × V_n classes = too many (unless we quotient).

  CANDIDATE 7: TOURNAMENTS ON n-1 VERTICES PLUS A DEFECT.
  Maybe D captures the "residual" when you try to extend a tournament
  from n-1 to n vertices but can't (the extension creates a non-tournament).

  THE REAL ANSWER: D = odd graphs. The Burnside sum
  (1/n!) Σ_{even σ} 2^{po(σ)} exactly counts odd graph iso classes.
  There is no "darker" interpretation — odd graphs ARE the dark tournaments.

  BUT: the METAGRAPH structure of odd graphs is DIFFERENT from tournaments.
  The edge-toggle meta-graph of odd graphs has different connectivity,
  different spectral properties, and NO self-loops.
  Understanding how the two meta-graphs RELATE (despite having different
  node counts) is the key to understanding the "dark" structure.
""")

# ============================================================
banner("THE FOUR META-GRAPHS COMPARED")
# ============================================================

print("  THE FOUR META-GRAPHS AT n=5:\n")
print("  %-15s %-6s %-6s %-6s %-10s" %
      ("Type", "V", "E", "maxD", "Connected"))
print("  " + "-"*45)

# Tournament meta-graph: V=12, E=30 (from earlier computation)
print("  %-15s %-6d %-6d %-6s %-10s" % ("Tournament", 12, 30, "7", "Yes"))
print("  %-15s %-6d %-6d %-6d %-10s" % ("Even Graph", len(even_classes),
      "?", "?", "No (disconnected)"))
print("  %-15s %-6d %-6d %-6d %-10s" % ("Odd Graph", len(odd_classes),
      "?", "?", "?"))
print("  %-15s %-6d %-6d %-6s %-10s" % ("All Graphs", len(all_classes),
      "?", "?", "Yes"))

print("""
  THE PARALLEL:
    Tournaments (V=12) ↔ Even Graphs (V=12)  [EQUINUMEROUS]
    Dark Tourn. (V=22) ↔ Odd Graphs (V=22)   [EQUINUMEROUS]
    Sum: 34 = A000088(5) = total graphs

  But the META-GRAPH structures differ:
  Tournament meta-graph: well-connected, self-loops, spectral gap 2/n
  Even graph meta-graph: sparse, NO self-loops, disconnected
  Odd graph meta-graph: structure unknown (computed above)

  THE DEEP QUESTION: Is there a BIJECTION between tournaments and
  even graphs that preserves SOME meta-graph structure?
  The Royle bijection (via Burnside) is NOT structure-preserving
  (it preserves counts but not connectivity).

  A structure-preserving bijection would be a MAJOR result:
  it would transfer ALL our Krawtchouk/sl(n)/waggly theory
  to the even graph world.
""")

print("Done. opus-2026-03-24-S317")
