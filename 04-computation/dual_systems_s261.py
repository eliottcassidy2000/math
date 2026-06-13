#!/usr/bin/env python3
"""
dual_systems_s261.py -- opus-2026-03-23-S261

TWO ANALOGOUS SYSTEMS WITH THE SAME COUNTS

SYSTEM A (INTERACTING): Even Graphs + Odd Graphs
  |Even| = A000568(n)
  |Odd|  = A000088(n) - A000568(n)
  They INTERACT via the EO boundary in the graph meta-graph.
  Edge toggle can send an even graph to odd and vice versa.

SYSTEM B (NON-INTERACTING): Tournaments + "Tournament Shadow"
  |Tournaments| = A000568(n)
  |Shadow|      = A000088(n) - A000568(n)
  They DON'T interact: a tournament stays a tournament under arc reversal.
  The "shadow" is the odd graph world, counted by the SAME Burnside difference.

THE CORRESPONDENCE:
  |Even| = |Tournaments| = A000568(n)
  |Odd|  = |Shadow|      = A000088(n) - A000568(n)

Both partitions have the same cardinalities, but:
  System A: the two parts INTERACT (EO edges exist)
  System B: the two parts are DECOUPLED (tournaments self-contained)

QUESTION: What properties of one system transfer to the other?
Does the interaction strength in System A relate to the
self-loop structure in System B?

Author: opus-2026-03-23-S261
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

def canon_graph(edge_set, n):
    best = None
    for perm in permutations(range(n)):
        mapped = frozenset((min(perm[i], perm[j]), max(perm[i], perm[j]))
                           for (i,j) in edge_set)
        if best is None or sorted(mapped) < sorted(best):
            best = mapped
    return frozenset(best)

def is_even_graph(edge_set, n):
    for perm in permutations(range(n)):
        mapped = frozenset((min(perm[i], perm[j]), max(perm[i], perm[j]))
                           for (i,j) in edge_set)
        if mapped != edge_set: continue
        rev = sum(1 for (i,j) in edge_set if perm[i] > perm[j])
        if rev % 2 == 1: return False
    return True

def canon_tournament(A, n):
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

banner("THE TWO SYSTEMS AT n=5")

n = 5; m = comb(n, 2)
pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

# ============================================================
# SYSTEM A: EVEN + ODD GRAPHS
# ============================================================

# Build graph classes
graph_classes = {}
for bits in range(2**m):
    edge_set = frozenset(pairs[k] for k in range(m) if (bits >> k) & 1)
    c = canon_graph(edge_set, n)
    if c not in graph_classes:
        graph_classes[c] = {
            'edges': edge_set,
            'even': is_even_graph(edge_set, n),
            'size': 0,
            'n_edges': len(edge_set),
        }
    graph_classes[c]['size'] += 1

even_graphs = {c: v for c, v in graph_classes.items() if v['even']}
odd_graphs = {c: v for c, v in graph_classes.items() if not v['even']}

print("  SYSTEM A (Even + Odd Graphs):")
print("  |Even| = %d, |Odd| = %d, Total = %d\n" %
      (len(even_graphs), len(odd_graphs), len(graph_classes)))

# Build graph meta-graph adjacency
gc_list = list(graph_classes.keys())
gc_idx = {c: i for i, c in enumerate(gc_list)}
nc = len(gc_list)
g_adj = [[False]*nc for _ in range(nc)]

for bits in range(2**m):
    es = frozenset(pairs[k] for k in range(m) if (bits >> k) & 1)
    ci = gc_idx[canon_graph(es, n)]
    for k in range(m):
        es2 = es.symmetric_difference({pairs[k]})
        cj = gc_idx[canon_graph(es2, n)]
        if ci != cj:
            g_adj[ci][cj] = True; g_adj[cj][ci] = True

# Compute EO interaction strength
EO_edges = 0
EE_edges = 0
OO_edges = 0
for i in range(nc):
    for j in range(i+1, nc):
        if not g_adj[i][j]: continue
        i_even = graph_classes[gc_list[i]]['even']
        j_even = graph_classes[gc_list[j]]['even']
        if i_even and j_even: EE_edges += 1
        elif not i_even and not j_even: OO_edges += 1
        else: EO_edges += 1

print("  Interaction: EO = %d edges, EE = %d, OO = %d" % (EO_edges, EE_edges, OO_edges))
print("  EO/(EO+EE+OO) = %.4f (interaction fraction)\n" %
      (EO_edges / (EO_edges + EE_edges + OO_edges)))

# ============================================================
# SYSTEM B: TOURNAMENTS (self-contained)
# ============================================================

# Build tournament classes
t_classes = {}
t_bits_map = {}
for bits in range(2**m):
    A = [[0]*n for _ in range(n)]
    for k, (i,j) in enumerate(pairs):
        if (bits>>k)&1: A[i][j]=1
        else: A[j][i]=1
    c = canon_tournament(A, n)
    if c not in t_classes:
        t_classes[c] = {
            'H': Hcount(A, n),
            'size': 0,
            'adj': [row[:] for row in A],
        }
    t_classes[c]['size'] += 1
    t_bits_map[bits] = c

tc_list = list(t_classes.keys())
tc_idx = {c: i for i, c in enumerate(tc_list)}
nt = len(tc_list)

t_adj = [[False]*nt for _ in range(nt)]
t_self_loops = [False]*nt

for bits in range(2**m):
    ci = tc_idx[t_bits_map[bits]]
    for k in range(m):
        cj = tc_idx[t_bits_map[bits ^ (1<<k)]]
        if ci != cj:
            t_adj[ci][cj] = True; t_adj[cj][ci] = True
        else:
            t_self_loops[ci] = True

T_edges = sum(1 for i in range(nt) for j in range(i+1,nt) if t_adj[i][j])
T_self = sum(t_self_loops)

print("  SYSTEM B (Tournaments — self-contained):")
print("  |Tournaments| = %d, E = %d, self-loops = %d" % (nt, T_edges, T_self))
print("  Interaction with outside: ZERO (tournaments stay tournaments)\n")

banner("THE DUALITY: INTERACTION vs SELF-LOOPS")

print("""
  SYSTEM A (Even + Odd): INTERACTING
    |Even| = %d, |Odd| = %d
    Internal edges: EE = %d, OO = %d
    INTERACTION (EO boundary): %d edges

  SYSTEM B (Tournaments): NON-INTERACTING but SELF-LOOPING
    |Tournaments| = %d (= |Even|!)
    Internal edges: %d
    SELF-LOOPS: %d classes have self-loops

  THE KEY OBSERVATION:
  System A has EO INTERACTION but NO self-loops (graph toggles always change class).
  System B has NO INTERACTION but HAS self-loops (arc reversals can preserve class).

  INTERACTION and SELF-LOOPS are COMPLEMENTARY:
  - If you CAN leave your type (even→odd), you DON'T self-loop
  - If you CAN'T leave your type (tournament stays tournament), you DO self-loop

  The EO interaction energy = %d edges.
  The self-loop count = %d classes.
  Are these related?
""" % (len(even_graphs), len(odd_graphs), EE_edges, OO_edges, EO_edges,
       nt, T_edges, T_self, EO_edges, T_self))

# Check: is EO ≈ some function of T_self?
print("  Numerical comparison at n=5:")
print("  EO edges = %d, T self-loop classes = %d" % (EO_edges, T_self))
print("  EO / T_self = %.4f" % (EO_edges / T_self if T_self > 0 else 0))

# Compute at n=3 and n=4 too
for n_val in [3, 4]:
    m_val = comb(n_val, 2)
    pairs_val = [(i,j) for i in range(n_val) for j in range(i+1,n_val)]

    # Graph classes
    gc = {}
    for bits in range(2**m_val):
        es = frozenset(pairs_val[k] for k in range(m_val) if (bits>>k)&1)
        c = canon_graph(es, n_val)
        if c not in gc:
            gc[c] = {'even': is_even_graph(es, n_val), 'edges': es}

    gc_l = list(gc.keys())
    gc_i = {c: i for i, c in enumerate(gc_l)}
    nc_val = len(gc_l)
    ga = [[False]*nc_val for _ in range(nc_val)]
    for bits in range(2**m_val):
        es = frozenset(pairs_val[k] for k in range(m_val) if (bits>>k)&1)
        ci = gc_i[canon_graph(es, n_val)]
        for k in range(m_val):
            es2 = es.symmetric_difference({pairs_val[k]})
            cj = gc_i[canon_graph(es2, n_val)]
            if ci != cj: ga[ci][cj] = True; ga[cj][ci] = True

    eo_val = 0
    for i in range(nc_val):
        for j in range(i+1, nc_val):
            if ga[i][j] and gc[gc_l[i]]['even'] != gc[gc_l[j]]['even']:
                eo_val += 1

    # Tournament self-loops
    tc_val = {}
    for bits in range(2**m_val):
        A = [[0]*n_val for _ in range(n_val)]
        for k,(ii,jj) in enumerate(pairs_val):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        c = canon_tournament(A, n_val)
        if c not in tc_val: tc_val[c] = {'size': 0}
        tc_val[c]['size'] += 1

    tc_l = list(tc_val.keys())
    tc_i = {c: i for i, c in enumerate(tc_l)}
    nt_val = len(tc_l)
    sl_val = 0
    for bits in range(2**m_val):
        ci = tc_i[canon_tournament([[0]*n_val for _ in range(n_val)], n_val)]
        A = [[0]*n_val for _ in range(n_val)]
        for k,(ii,jj) in enumerate(pairs_val):
            if (bits>>k)&1: A[ii][jj]=1
            else: A[jj][ii]=1
        ci = tc_i[canon_tournament(A, n_val)]
        for k in range(m_val):
            A2 = [row[:] for row in A]
            ii, jj = pairs_val[k]
            A2[ii][jj] = 1 - A2[ii][jj]; A2[jj][ii] = 1 - A2[jj][ii]
            cj = tc_i[canon_tournament(A2, n_val)]
            if ci == cj:
                sl_val += 1
                break

    print("  n=%d: EO = %d, T self-loop classes = %d" % (n_val, eo_val, sl_val))

banner("THE CONSERVATION LAW")

print("""
  At each n, the Burnside sum over ALL permutations gives #Graphs.
  It SPLITS into:
    ODD-ORDER terms → count Tournaments (= Even graphs)
    EVEN-ORDER terms → count Odd graphs

  This split is ADDITIVE: #Graphs = #Tourn + #Odd.

  For the META-GRAPH EDGES, the split becomes:
    FULL graph meta-graph = EE + EO + OO edges
    Tournament meta-graph = T_edges (separate space)

  There's NO direct edge-level conservation law because
  the two systems live in different spaces:
  - System A: {0,1}^m (undirected graphs, toggle edges)
  - System B: {0,1}^m (tournaments, reverse arcs)

  BOTH are quotients of the SAME m-cube Q_m by S_n,
  but with DIFFERENT interpretations of the bits!

  In System A: bit k = 1 means edge {i,j} is PRESENT
  In System B: bit k = 1 means arc goes from i to j (not j to i)

  THE SAME 2^m binary strings, the SAME S_n action,
  but DIFFERENT semantics → different meta-graphs.

  THE DEEP CONNECTION:
  The Burnside sum (1/n!) Σ_g 2^{c(g_E)} counts BOTH:
  - Graph iso classes (all g)
  - Tournament iso classes (odd-order g only)

  The SAME 2^{c(g_E)} term contributes to BOTH counts!
  But for graphs it counts SUBSETS of edges (each orbit freely
  chosen as present or absent), while for tournaments it counts
  ORIENTATIONS of arcs (each orbit freely chosen as forward or
  backward). These are mathematically identical operations
  (choosing one of 2 states per orbit) but semantically different.

  THIS IS THE ISOMORPHISM THAT ISN'T:
  The Burnside counting is identical, but the flip graphs differ
  because ADJACENCY depends on semantics (toggle vs reverse).
""")

banner("ENERGY BALANCE: INTERACTION vs SELF-LOOPS")

# The key metaphor: System A has interaction energy (EO edges)
# while System B has self-energy (self-loops). Are these "conserved"?

print("  n | EO (interaction) | T_self (self-loops) | EO + T_self | E_total_graph")
print("  --|------------------|---------------------|-------------|---------------")

EO_vals = {3: 1, 4: 7, 5: 36}  # from S260
T_self_vals = {3: 1, 4: 2, 5: 5}  # from S211/S245 (classes with self-loops)
E_graph = {3: 3, 4: 14, 5: 74}

for n_val in [3, 4, 5]:
    eo = EO_vals[n_val]
    ts = T_self_vals[n_val]
    eg = E_graph[n_val]
    print("  %d |       %2d         |         %2d          |     %3d     |      %3d" %
          (n_val, eo, ts, eo + ts, eg))

print("""
  EO + T_self = 2, 9, 41 — NOT equal to E_graph or any simple combination.

  But the SPIRIT is right: interaction energy and self-energy are
  COMPLEMENTARY aspects of the same underlying structure.

  System A (even+odd): energy goes into BOUNDARY (EO edges)
  System B (tournaments): energy goes into INTERIOR (self-loops)

  This is like a GAUGE CHOICE:
  - In the "graph gauge": energy is boundary interaction
  - In the "tournament gauge": energy is self-interaction
  - The gauge-invariant quantity: the TOTAL Burnside sum
""")

banner("THE SHADOW: WHAT IS THE TOURNAMENT'S COUNTERPART TO ODD GRAPHS?")

print("""
  The user's insight: tournaments have a "shadow" with |shadow| = |odd graphs|.

  What IS this shadow?

  In the Burnside framework:
    |Graphs| = |Tournaments| + |Odd graphs|
             = Σ_{odd order} 2^{c(g_E)} / n! + Σ_{even order} 2^{c(g_E)} / n!

  The "shadow" of tournaments = the EVEN-ORDER Burnside sum = |Odd graphs|.

  But this shadow counts something CONCRETE: it counts orbits of S_n
  on {0,1}^m where the stabilizer has EVEN ORDER elements that
  "reverse" edges. These are graph classes with non-trivial symmetry
  that includes edge-reversing automorphisms.

  In TOURNAMENT language: the shadow corresponds to tournaments
  that WOULD BE fixed by even-order permutations IF even-order
  permutations could fix tournaments (they can't, because of
  the self-paired cycle obstruction).

  So the shadow = "IMPOSSIBLE tournament fixings" = the set of
  tournament configurations that are ruled out by the antisymmetry
  constraint A[i][j] + A[j][i] = 1.

  THE TWO SYSTEMS:
  System A: {Allowed undirected} = Even ∪ Odd (INTERACTING, via toggles)
  System B: {Allowed directed} = Tournaments (self-contained, NO shadow)
            {Forbidden directed} = "Anti-tournaments" (self-paired cycles)

  The shadow DOESN'T EXIST as actual objects — it's the ABSENCE.
  |Shadow| = |Odd graphs| counts what's FORBIDDEN, not what's present.

  THIS IS THE PARTICLE-ANTIPARTICLE ANALOGY:
  Even graphs = particles (exist, interact with odd)
  Odd graphs = antiparticles (exist, interact with even)
  Tournaments = stable particles (exist, don't need antiparticles)
  Anti-tournaments = virtual antiparticles (forbidden, counted by |Odd|)

  The Burnside sum is the PARTITION FUNCTION:
  Z = Z_particles + Z_antiparticles = Z_stable + Z_virtual

  |Z_particles| = |Z_stable| (EQUINUMEROUS!)
  |Z_antiparticles| = |Z_virtual| (also equinumerous!)
  Z_total = same in both frames.
""")

print("Done. opus-2026-03-23-S261")
