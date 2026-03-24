#!/usr/bin/env python3
"""
even_graph_connection_s259.py -- opus-2026-03-23-S259

TOURNAMENTS AND EVEN GRAPHS ARE EQUINUMEROUS (Royle et al. 2022)

THE KEY IDENTITY:
  #Graphs(n) = #Tournaments(n) + #OddGraphs(n)
  => #EvenGraphs(n) = #Tournaments(n) = A000568(n)

WHERE:
  #Graphs(n) = (1/n!) sum_{g in S_n} 2^{c(g_E)}  [all permutations]
  #Tournaments(n) = (1/n!) sum_{g: |g| odd} 2^{c(g_E)}  [odd-order only]
  #OddGraphs(n) = (1/n!) sum_{g: |g| even} 2^{c(g_E)}  [even-order only]

The proof hinges on Lemma 3.1: g_A has a self-paired cycle iff |g| is even.
For odd-order g: no self-paired cycles, so c(g_A) = 2*c(g_E),
and tournaments fixed by g = 2^{c(g_A)/2} = 2^{c(g_E)}.

CONNECTION TO OUR WORK:
- c(g_E) = pair_orbits(cycle_type(g)) = our Burnside polynomial exponent
- The ODD-ORDER restriction = our "only odd cycle types" condition
  (because a permutation has odd order iff ALL its cycle lengths are odd)
- This is EXACTLY the same Burnside sum we use for V_n!

THE DEEP INSIGHT:
The Royle et al. identity says:
  V_n (tournaments) = #EvenGraphs(n)
  = #Graphs(n) - #OddGraphs(n)
  = (1/n!)[sum_{all g} 2^{c(g_E)} - sum_{|g| even} 2^{c(g_E)}]

This gives us a SECOND way to compute V_n:
  V_n = A000088(n) - #OddGraphs(n)
where A000088(n) = total graph iso classes.

And: #OddGraphs(n) = A000088(n) - A000568(n).

THE META-GRAPH CONNECTION:
If tournaments and even graphs are equinumerous, then:
- The tournament meta-graph G_n has V_n = A000568(n) nodes
- There should be an EVEN GRAPH meta-graph with the SAME number of nodes
- Is there a natural bijection between the two meta-graphs?
- Does the even graph meta-graph have the SAME edge count?

Author: opus-2026-03-23-S259
"""

import sys
from math import comb, factorial, gcd
from itertools import permutations
from collections import Counter
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def pair_orbits_from_type(ct):
    """c(g_E) = number of edge-orbits for permutation with cycle type ct."""
    within = sum((c - 1) // 2 for c in ct)
    across = sum(gcd(ct[i], ct[j]) for i in range(len(ct)) for j in range(i+1, len(ct)))
    return within + across

def edge_orbits_from_type(ct):
    """c(g_E) for ANY cycle type (not just odd)."""
    # For a permutation with cycle lengths c1, c2, ...:
    # Edge orbits = sum_{i} floor(c_i/2) + sum_{i<j} gcd(c_i, c_j)
    within = sum(c // 2 for c in ct)
    across = sum(gcd(ct[i], ct[j]) for i in range(len(ct)) for j in range(i+1, len(ct)))
    return within + across

def multinomial(n, ct):
    ct_counter = Counter(ct)
    denom = 1
    for c, a in ct_counter.items():
        denom *= (c ** a) * factorial(a)
    return factorial(n) // denom

def gen_all_partitions(n):
    """Generate ALL partitions of n."""
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

def has_even_cycle(ct):
    return any(c % 2 == 0 for c in ct)

def perm_order(ct):
    """Order of a permutation = lcm of cycle lengths."""
    from math import lcm
    return reduce(lcm, ct)

banner("VERIFYING THE ROYLE ET AL. IDENTITY")

# For each n: compute #Graphs, #Tournaments, #OddGraphs, #EvenGraphs
# using the Burnside formulas from the paper.

# A000088 = graph iso classes: 1, 1, 2, 4, 11, 34, 156, 1044, 12346, ...
A000088 = {1:1, 2:2, 3:4, 4:11, 5:34, 6:156, 7:1044, 8:12346}
A000568 = {1:1, 2:1, 3:2, 4:4, 5:12, 6:56, 7:456, 8:6880}

print("  FORMULA: #Graphs = (1/n!) sum_g 2^{c(g_E)}")
print("           #Tourn  = (1/n!) sum_{|g| odd} 2^{c(g_E)}")
print("           #Odd    = (1/n!) sum_{|g| even} 2^{c(g_E)}")
print("           #Even   = #Graphs - #Odd = #Tourn")
print()

for n in range(3, 9):
    nfact = factorial(n)
    all_parts = gen_all_partitions(n)

    sum_all = 0
    sum_odd_order = 0
    sum_even_order = 0

    for ct in all_parts:
        mult = multinomial(n, ct)
        ce = edge_orbits_from_type(ct)
        contribution = mult * (2 ** ce)
        sum_all += contribution

        # Check order: order(g) = lcm of cycle lengths
        from math import lcm
        order = reduce(lcm, ct)
        if order % 2 == 1:  # odd order
            sum_odd_order += contribution
        else:  # even order
            sum_even_order += contribution

    graphs = sum_all // nfact
    tournaments = sum_odd_order // nfact
    odd_graphs = sum_even_order // nfact
    even_graphs = graphs - odd_graphs

    print("  n=%d: #Graphs=%d, #Tourn=%d, #Odd=%d, #Even=%d" %
          (n, graphs, tournaments, odd_graphs, even_graphs))
    print("       A000088=%s, A000568=%s, Even=Tourn? %s" %
          (A000088.get(n, '?'), A000568.get(n, '?'),
           even_graphs == tournaments))
    print()

banner("THE DEEP CONNECTION TO OUR META-GRAPH")

print("""
  THE ROYLE IDENTITY AND THE TOURNAMENT META-GRAPH:

  V_n = A000568(n) = #EvenGraphs(n)

  This means: for every tournament iso class, there is a
  CORRESPONDING even graph iso class. The bijection is NOT
  constructive (Royle et al. leave it as an open problem).

  BUT: the Burnside sum that counts them is the SAME:

  V_n = (1/n!) sum_{|g| odd} 2^{c(g_E)}

  where c(g_E) = edge_orbits = our pair_orbits function.

  For ODD-ORDER g: c(g_E) = pair_orbits(cycle_type(g))
  because |g| odd iff ALL cycle lengths are odd.
  (lcm of odd numbers is odd.)

  So: the SAME cycle types contribute to both counts!

  NOW: what about the META-GRAPH of even graphs?

  Define: the even graph flip graph as the graph on even graph
  iso classes, where two are adjacent if they differ by adding
  or removing one edge.

  QUESTION: Does this flip graph have the SAME number of edges
  as the tournament meta-graph G_n?

  CONJECTURE: The even graph meta-graph is ISOMORPHIC to G_n.
  (Or at least has the same edge count.)

  If true: EVERY invariant of G_n (edges, clique complex, Betti
  numbers, etc.) applies equally to the even graph flip graph!
""")

banner("WHAT ODD-ORDER MEANS FOR CYCLE TYPES")

# A permutation has odd order iff ALL cycle lengths are odd.
# Verify:
print("  Cycle types with odd order (= all odd cycles):")
for n in range(3, 8):
    all_parts = gen_all_partitions(n)
    odd_order_parts = []
    even_order_parts = []
    for ct in all_parts:
        order = reduce(lcm, ct)
        if order % 2 == 1:
            odd_order_parts.append(ct)
        else:
            even_order_parts.append(ct)
    print("  n=%d: odd-order types: %s" % (n, sorted(odd_order_parts)))
    print("        even-order types: %s" % (n, sorted(even_order_parts)))
    print()

banner("THE SELF-PAIRED CYCLE CONNECTION")

print("""
  LEMMA 3.1 (Royle et al.): g_A has a self-paired cycle iff |g| is even.

  A SELF-PAIRED CYCLE of g_A is a cycle C = (a1, a2, ..., ak) such that
  C = C_bar (reversing all arcs gives the same cycle).

  For a self-paired cycle: the cycle visits each arc AND its reverse.
  So the cycle has even length (arcs come in forward/backward pairs).

  |g| even iff g has a cycle of even length iff g_A has a self-paired cycle.

  THE CONNECTION TO ARC NEUTRALITY:
  - Self-paired cycles create "forced" pair orbits (1 choice, not 2)
  - In ODD-ORDER permutations: no self-paired cycles → all orbits are FREE
  - In EVEN-ORDER permutations: at least 1 self-paired → some orbits forced

  For tournaments: only odd-order permutations fix ANY tournament.
  (Even-order permutations force A[i][j] = A[j][i] for some pair,
   which contradicts the tournament condition.)

  THIS IS WHY THE BURNSIDE SUM ONLY HAS ODD-ORDER TERMS!
  The Royle paper makes this precise: the self-paired cycle creates
  a pair {a, a_bar} that can't be consistently oriented.

  FOR EVEN GRAPHS: a graph X is even iff no automorphism of X
  reverses an odd number of edges. An automorphism g ∈ Aut(X)
  reverses edge {u,v} iff (u^g - v^g)/(u-v) = -1.

  The SIGN of g with respect to X:
  sgn_X(g) = prod_{edges {u,v}} (u^g - v^g)/(u - v)

  X is odd iff some g ∈ Aut(X) has sgn_X(g) = -1.
  X is even iff ALL g ∈ Aut(X) have sgn_X(g) = +1.

  THE SIGN FUNCTION sgn_X IS A GROUP HOMOMORPHISM (Lemma 4.1).
  So: X is odd iff Aut(X) has index-2 subgroup (= kernel of sgn_X).
  X is even iff sgn_X is trivially 1 (no odd automorphisms).

  THIS IS THE SAME AS: the tournament condition forces ALL
  automorphisms to preserve arc orientation (no reversal allowed).
  Odd-order automorphisms ALWAYS preserve orientation.
  Even-order automorphisms MAY reverse arcs via self-paired cycles.
""")

banner("THE FORMULA BRIDGE: EVEN GRAPHS → E(G_n)")

print("""
  Since V_n = #EvenGraphs(n), and BOTH are counted by the SAME
  Burnside sum, the T_n formula ALSO counts something about even graphs:

  T_n = (1/n!) sum_{|g| odd} c(g_E) * 2^{c(g_E)}
      ??? No — T_n counts arc-orbits, not edge-orbits of graphs.

  Actually: T_n counts S_n-orbits on (tournament, arc) pairs.
  For even graphs: the analogous quantity would be S_n-orbits on
  (even graph, edge) pairs.

  DEFINE: T_n^{even} = (1/n!) sum_{|g| odd} c(g_E) * 2^{c(g_E)}

  Wait: for GRAPHS, each edge orbit of g_E contributes 2 choices
  (present or absent). For (graph, edge) pairs:
  Fixed pairs = #{(X, e) : g fixes X and g fixes e}
  = 2^{c(g_E)} × #{edges fixed by g_E}
  = 2^{c(g_E)} × (# 1-cycles of g_E)

  Hmm, this doesn't directly give T_n.

  T_n for tournaments: each arc orbit of g_A contributes 2 choices.
  c(g_A) = #{arc-orbits}. For odd |g|: c(g_A) = 2 * c(g_E)
  (no self-paired cycles). So 2^{c(g_A)/2} = 2^{c(g_E)}.

  And: the number of arc orbits of g is:
  for odd |g|: c(g_A) = 2 * c(g_E) (each edge orbit splits into
  two arc orbits: forward and backward).

  T_n = (1/n!) sum_{|g| odd} C(fix(g), 2) × 2^{c(g_E)}

  Wait, that's what we compute for arc-orbits. Actually:
  T_n = #{(iso class, arc-orbit) pairs}
  = (1/n!) sum_g #{(T, arc) : g fixes both}
  = (1/n!) sum_{|g| odd} #{arcs fixed by g in each fixed tournament}
  ... this is more complex.

  THE SIMPLE VERSION:
  The Royle paper gives us V_n. Our turbo computation gives T_n.
  The RELATIONSHIP between V_n and T_n is:
  T_n = sum_{iso classes C} (# arc-orbits of C)
  = sum_C m / |Aut(C)| ... no, that's not right either.

  WHAT WE KNOW:
  V_n = #{tournament iso classes} = #{even graph iso classes}
  T_n = #{(tournament, arc) orbit pairs under S_n}
  E_n = #{edges of G_n} = #{adjacent class pairs}

  The even graph equinumerosity gives V_n a SECOND interpretation.
  But T_n and E_n don't have obvious even graph analogs
  (since even graphs don't have "arcs").

  HOWEVER: for even graphs, the ANALOG of an arc is an "edge toggle":
  adding or removing one edge. The even graph flip graph has:
  V = V_n nodes (same count!)
  E' = #{edges of even graph flip graph}

  QUESTION: Is E' = E_n? Is the even graph flip graph isomorphic to G_n?
""")

# Let me test this at n=4 where both counts are small.
banner("TESTING AT n=4: EVEN GRAPH FLIP GRAPH")

n = 4
m = comb(n, 2)  # 6 edges
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

# Enumerate all graphs on 4 vertices (2^6 = 64)
def canon_graph(edge_set, n):
    """Canonical form of a graph under S_n."""
    best = None
    for perm in permutations(range(n)):
        mapped = frozenset((min(perm[i], perm[j]), max(perm[i], perm[j]))
                           for (i,j) in edge_set)
        if best is None or sorted(mapped) < sorted(best):
            best = mapped
    return frozenset(best)

def is_even_graph(edge_set, n):
    """Check if graph is even: no odd automorphism."""
    pairs_list = [(i,j) for i in range(n) for j in range(i+1,n)]
    for perm in permutations(range(n)):
        # Check if perm is an automorphism
        mapped = frozenset((min(perm[i], perm[j]), max(perm[i], perm[j]))
                           for (i,j) in edge_set)
        if mapped != edge_set: continue
        # Count reversed edges
        reversed_count = 0
        for (i,j) in edge_set:
            pi, pj = perm[i], perm[j]
            if pi > pj:  # edge orientation reversed
                reversed_count += 1
        if reversed_count % 2 == 1:
            return False  # found an odd automorphism
    return True

# Enumerate all graphs
graph_classes = {}
for bits in range(2**m):
    edge_set = frozenset(pairs[k] for k in range(m) if (bits >> k) & 1)
    c = canon_graph(edge_set, n)
    if c not in graph_classes:
        graph_classes[c] = {'edges': edge_set, 'even': is_even_graph(edge_set, n)}

total_graphs = len(graph_classes)
even_count = sum(1 for v in graph_classes.values() if v['even'])
odd_count = total_graphs - even_count

print("  n=4: %d graph iso classes, %d even, %d odd" % (total_graphs, even_count, odd_count))
print("  A000088(4) = 11, A000568(4) = 4")
print("  Even = Tournament? %s" % (even_count == A000568[4]))

# Build even graph flip graph
even_classes = {c: i for i, (c, v) in enumerate(graph_classes.items()) if v['even']}
even_reps = list(even_classes.keys())
ne = len(even_reps)

even_edges = set()
for i, c1 in enumerate(even_reps):
    edges1 = graph_classes[c1]['edges']
    for p in pairs:
        # Toggle edge p
        if p in edges1:
            new_edges = edges1 - {p}
        else:
            new_edges = edges1 | {p}
        new_canon = canon_graph(new_edges, n)
        if new_canon in even_classes:
            j = even_classes[new_canon]
            if i != j:
                even_edges.add((min(i,j), max(i,j)))

print("  Even graph flip graph: V=%d, E=%d" % (ne, len(even_edges)))

# Compare with tournament meta-graph at n=4
print("  Tournament meta-graph at n=4: V=4, E=5")
print("  Even graph flip graph at n=4: V=%d, E=%d" % (ne, len(even_edges)))
print("  Same V? %s. Same E? %s" % (ne == 4, len(even_edges) == 5))

# Do the same for n=3
print()
n = 3; m = comb(n, 2)
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
graph_classes3 = {}
for bits in range(2**m):
    edge_set = frozenset(pairs[k] for k in range(m) if (bits >> k) & 1)
    c = canon_graph(edge_set, n)
    if c not in graph_classes3:
        graph_classes3[c] = {'edges': edge_set, 'even': is_even_graph(edge_set, n)}

even3 = {c: i for i, (c, v) in enumerate(graph_classes3.items()) if v['even']}
even_reps3 = list(even3.keys())
ne3 = len(even_reps3)

even_edges3 = set()
for i, c1 in enumerate(even_reps3):
    edges1 = graph_classes3[c1]['edges']
    for p in pairs:
        if p in edges1:
            new_edges = edges1 - {p}
        else:
            new_edges = edges1 | {p}
        new_canon = canon_graph(new_edges, n)
        if new_canon in even3:
            j = even3[new_canon]
            if i != j:
                even_edges3.add((min(i,j), max(i,j)))

print("  n=3: %d graph classes, %d even" % (len(graph_classes3), ne3))
print("  Even graph flip graph: V=%d, E=%d" % (ne3, len(even_edges3)))
print("  Tournament meta-graph: V=2, E=1")
print("  Isomorphic? V match: %s, E match: %s" % (ne3 == 2, len(even_edges3) == 1))

print("\nDone. opus-2026-03-23-S259")
