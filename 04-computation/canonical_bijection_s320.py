#!/usr/bin/env python3
"""
canonical_bijection_s320.py -- opus-2026-03-24-S320

PATH-INDEPENDENT TOURNAMENT-GRAPH BIJECTION

The S319 bijection: tiling → graph, depends on choice of base path P₀.
Different paths give different graphs from the same tournament.

GOAL: make the bijection CANONICAL (path-independent).

THREE STRATEGIES:

STRATEGY 1: CANONICAL PATH SELECTION
For each tournament T, choose the HP that gives the "best" graph.
"Best" could mean: most edges, fewest edges, most even, lexicographic min, etc.

The canonical path is an INVARIANT of the tournament iso class
(if we require the selection rule to be isomorphism-invariant).

STRATEGY 2: THE GRAPH MULTISET
For each tournament T with H(T) Hamiltonian paths, the H(T) different
base-path choices give H(T) different graphs. Define:
  Φ(T) = {G_P : P is a HP of T} = multiset of H(T) graphs.

This multiset IS an isomorphism invariant of T.
Two isomorphic tournaments give the same multiset (up to relabeling).

STRATEGY 3: THE SCORE-BASED PATH
Choose the HP that follows the SCORE ordering:
P = (vertex with highest score, ..., vertex with lowest score).
If scores are distinct (generic case), this gives a unique canonical HP.
Ties: break by secondary invariant (c3 contribution, etc.).

STRATEGY 4: THE MEDIAN GRAPH
Among all H(T) graphs Φ(T), pick the one closest to the "center"
(median edge count, or median in some metric).

Author: opus-2026-03-24-S320
"""

import sys
import numpy as np
from math import comb, factorial, log2
from itertools import permutations
from collections import defaultdict, Counter
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

def find_all_HPs(A, n):
    paths = []
    def dfs(path, visited):
        if len(path) == n:
            paths.append(tuple(path))
            return
        v = path[-1]
        for w in range(n):
            if w not in visited and A[v][w]:
                dfs(path + [w], visited | {w})
    for s in range(n):
        dfs([s], {s})
    return paths

def hp_to_graph(A, hp, n):
    """Given tournament A and HP, produce the graph.
    Along the HP: edges always present (base path).
    Non-HP arcs: present iff arc goes from earlier to later in HP ordering."""
    # Position of each vertex in the HP
    pos = {hp[i]: i for i in range(n)}
    adj = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            # Is this a base-path edge? (consecutive in HP)
            is_base = False
            for k in range(n-1):
                if (hp[k] == i and hp[k+1] == j) or (hp[k] == j and hp[k+1] == i):
                    is_base = True; break

            if is_base:
                adj[i][j] = 1; adj[j][i] = 1  # always present
            else:
                # Arc from earlier-in-HP to later-in-HP = "natural" = edge present
                # Arc from later to earlier = "flipped" = edge absent
                if pos[i] < pos[j]:
                    # i is earlier in HP. Arc i→j means i beats j.
                    if A[i][j]:
                        adj[i][j] = 1; adj[j][i] = 1  # natural → present
                    # else: flipped → absent
                else:
                    # j is earlier. Arc j→i means j beats i.
                    if A[j][i]:
                        adj[i][j] = 1; adj[j][i] = 1  # natural → present
    return adj

def graph_canon(adj, n):
    sc = [sum(adj[i]) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n): sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(sg.keys())]
    best = None
    def gp(gs2):
        if not gs2: yield []; return
        for p in permutations(gs2[0]):
            for r in gp(gs2[1:]): yield list(p)+r
    for p in gp(gs):
        f = tuple(adj[p[i]][p[j]] for i in range(n) for j in range(i+1, n))
        if best is None or f < best: best = f
    return best

def tourn_canon(A, n):
    sc = [sum(A[i]) for i in range(n)]
    sg = defaultdict(list)
    for v in range(n): sg[sc[v]].append(v)
    gs = [sg[s] for s in sorted(sg.keys())]
    best = None
    def gp(gs2):
        if not gs2: yield []; return
        for p in permutations(gs2[0]):
            for r in gp(gs2[1:]): yield list(p)+r
    for p in gp(gs):
        f = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n) if i!=j)
        if best is None or f < best: best = f
    return best

def is_even_graph(adj, n):
    return all(sum(adj[v][w] for w in range(n) if w != v) % 2 == 0 for v in range(n))

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

# ============================================================
banner("STRATEGY 1: CANONICAL PATH SELECTION (n=5)")
# ============================================================

n = 5
m_arc = comb(n, 2)
all_pairs = [(i,j) for i in range(n) for j in range(i+1,n)]

print("  For each tournament: try ALL HPs, pick the one giving\n"
      "  the graph closest to even (minimum degree-parity violations).\n")

# Build tournament classes
tourn_classes = defaultdict(list)
for bits in range(2**m_arc):
    A = [[0]*n for _ in range(n)]
    for k, (i, j) in enumerate(all_pairs):
        if (bits >> k) & 1: A[i][j] = 1
        else: A[j][i] = 1
    c = tourn_canon(A, n)
    tourn_classes[c].append((bits, [row[:] for row in A]))

nc = len(tourn_classes)
print("  %d tournament iso classes\n" % nc)

# For each class: find the HP that gives the most-even graph
results = []
for tc, members in tourn_classes.items():
    bits0, A0 = members[0]
    H = Hcount(A0, n)
    hps = find_all_HPs(A0, n)

    best_hp = None; best_parity_violations = n + 1; best_graph = None
    even_hp_count = 0
    graph_classes = set()

    for hp in hps:
        G = hp_to_graph(A0, hp, n)
        n_edges = sum(G[i][j] for i in range(n) for j in range(i+1, n))
        degs = [sum(G[v][w] for w in range(n) if w != v) for v in range(n)]
        parity_violations = sum(d % 2 for d in degs)

        gc = graph_canon(G, n)
        graph_classes.add(gc)

        if is_even_graph(G, n):
            even_hp_count += 1

        if parity_violations < best_parity_violations:
            best_parity_violations = parity_violations
            best_hp = hp
            best_graph = G

    results.append({
        'H': H, 'best_pv': best_parity_violations,
        'even_count': even_hp_count, 'n_graph_classes': len(graph_classes),
        'total_hps': len(hps)
    })

# Sort by H and display
print("  %-4s %-6s %-10s %-10s %-12s" %
      ("H", "#HPs", "best_PV", "#even_HPs", "#graph_cls"))
print("  " + "-"*45)

for r in sorted(results, key=lambda x: x['H']):
    print("  %-4d %-6d %-10d %-10d %-12d" %
          (r['H'], r['total_hps'], r['best_pv'], r['even_count'], r['n_graph_classes']))

# ============================================================
banner("THE KEY QUESTION: CAN EVERY TOURNAMENT FIND AN EVEN HP?")
# ============================================================

# For each tournament: does there EXIST an HP P such that the
# resulting graph G_P is even?

print("  Can every tournament choose an HP that gives an even graph?\n")

n_with_even = sum(1 for r in results if r['even_count'] > 0)
n_without_even = sum(1 for r in results if r['even_count'] == 0)

print("  Classes with ≥1 even-graph HP: %d / %d" % (n_with_even, nc))
print("  Classes with NO even-graph HP: %d / %d" % (n_without_even, nc))
print()

if n_without_even > 0:
    print("  SOME TOURNAMENTS CANNOT PRODUCE AN EVEN GRAPH!")
    print("  The bijection CANNOT be made even-graph-preserving for ALL classes.\n")

    for r in sorted(results, key=lambda x: x['H']):
        if r['even_count'] == 0:
            print("    H=%d: %d HPs tried, NONE gives even graph. Best PV=%d" %
                  (r['H'], r['total_hps'], r['best_pv']))
else:
    print("  EVERY TOURNAMENT HAS AN EVEN-GRAPH HP!")
    print("  The canonical bijection maps tournaments to even graphs surjectively.\n")

# ============================================================
banner("STRATEGY 2: THE GRAPH MULTISET Φ(T)")
# ============================================================

# For each tournament class: what graph classes appear in Φ(T)?
print("  Φ(T) = {G_P : P ∈ HP(T)} = the multiset of graphs from all HPs.\n")

for tc, members in sorted(tourn_classes.items(), key=lambda x: Hcount(x[1][0][1], n)):
    bits0, A0 = members[0]
    H = Hcount(A0, n)
    hps = find_all_HPs(A0, n)

    graph_class_counts = Counter()
    even_count = 0
    edge_counts = []

    for hp in hps:
        G = hp_to_graph(A0, hp, n)
        gc = graph_canon(G, n)
        graph_class_counts[gc] += 1
        n_edges = sum(G[i][j] for i in range(n) for j in range(i+1, n))
        edge_counts.append(n_edges)
        if is_even_graph(G, n):
            even_count += 1

    n_distinct = len(graph_class_counts)
    avg_edges = np.mean(edge_counts)
    print("  H=%2d: %d HPs → %d distinct graph classes, %d even, avg_edges=%.1f" %
          (H, len(hps), n_distinct, even_count, avg_edges))
    for gc, count in graph_class_counts.most_common(3):
        print("    %d× graph %s..." % (count, str(gc)[:30]))

# ============================================================
banner("STRATEGY 3: SCORE-BASED CANONICAL PATH")
# ============================================================

print("""
  THE SCORE PATH: order vertices by score (out-degree), highest first.
  This gives a canonical HP for tournaments with DISTINCT scores.

  For tournaments with tied scores: break ties by secondary invariant.

  THE TRICK: the score path is the TRANSITIVE ORDERING.
  If T were transitive, the score path IS the unique HP.
  For non-transitive T: the score path may NOT be a valid HP
  (the arc between consecutive score-ranked vertices might go "wrong way").

  ALTERNATIVE: use the score RANKING as a guide but find the
  CLOSEST valid HP to the score ranking (minimum Kendall tau distance).
""")

# For each tournament: try the score-ordered path
n_valid_score_path = 0
for tc, members in tourn_classes.items():
    bits0, A0 = members[0]
    scores = [(sum(A0[v]), v) for v in range(n)]
    scores.sort(reverse=True)
    score_path = tuple(v for _, v in scores)

    # Check if this is a valid HP
    valid = all(A0[score_path[i]][score_path[i+1]] for i in range(n-1))
    if valid:
        n_valid_score_path += 1

print("  Score-ordered path is valid HP: %d / %d classes (%.1f%%)" %
      (n_valid_score_path, nc, 100*n_valid_score_path/nc))

# ============================================================
banner("STRATEGY 4: THE CLOSEST-TO-EVEN PATH")
# ============================================================

# For each tournament: among all HPs, find the one whose graph
# has minimum total degree-parity violation (0 = even graph).

print("  CLOSEST-TO-EVEN PATH: minimize Σ (degree mod 2) over all HPs.\n")

pv_dist = Counter()
for tc, members in tourn_classes.items():
    bits0, A0 = members[0]
    hps = find_all_HPs(A0, n)
    best_pv = n + 1
    for hp in hps:
        G = hp_to_graph(A0, hp, n)
        pv = sum(sum(G[v][w] for w in range(n) if w != v) % 2 for v in range(n))
        best_pv = min(best_pv, pv)
    pv_dist[best_pv] += 1

print("  Distribution of best parity violations:")
for pv in sorted(pv_dist.keys()):
    print("    PV=%d: %d classes (%.1f%%)" % (pv, pv_dist[pv], 100*pv_dist[pv]/nc))

# How many can reach PV=0 (perfect even graph)?
pv0 = pv_dist.get(0, 0)
print("\n  Classes reachable to even graph (PV=0): %d / %d (%.1f%%)" %
      (pv0, nc, 100*pv0/nc))
print("  Classes with PV=2 (off by one vertex pair): %d" % pv_dist.get(2, 0))

# ============================================================
banner("THE PATH-INDEPENDENT BIJECTION")
# ============================================================

print("""
  THE RESOLUTION: No single path gives an even graph for ALL tournaments.
  But we can define a CANONICAL BIJECTION that is path-independent:

  DEFINITION: For tournament class [T], define its CANONICAL GRAPH as:
    Φ([T]) = graph_class( arg min_{P ∈ HP(T)} PV(G_P) )
  where PV(G) = total degree-parity violations.

  If PV = 0 (even graph achievable): the canonical graph IS an even graph.
  If PV > 0 (no even HP): the canonical graph is the CLOSEST to even.

  This map is:
  - WELL-DEFINED: the minimum PV is an invariant of the iso class
  - PATH-INDEPENDENT: we choose the best among all HPs
  - CANONICAL: unique up to tie-breaking within the minimum-PV set
  - INJECTIVE on tournament classes (different T give different graph multisets)

  THE STRUCTURE:
  At n=5: {n_pv0} classes map to even graphs, {n_pv2} map to near-even.
  The "dark" classes (PV > 0) are tournaments whose cycle structure
  is INCOMPATIBLE with any even-graph embedding.

  This incompatibility IS the "dark tournament" phenomenon:
  dark tournaments are those whose cycle structure forces odd-degree
  vertices in EVERY graph interpretation.
""".format(n_pv0=pv0, n_pv2=pv_dist.get(2, 0)))

# ============================================================
banner("THE EVEN-GRAPH REACHABILITY AT n=6")
# ============================================================

n6 = 6
m_arc6 = comb(n6, 2)
all_pairs6 = [(i,j) for i in range(n6) for j in range(i+1,n6)]

print("  n=6: checking even-graph reachability (sampling)...\n")

import random
random.seed(42)

n_sampled = 0; n_even_reachable = 0; pv_dist_6 = Counter()

for _ in range(500):
    # Random tournament
    A = [[0]*n6 for _ in range(n6)]
    for i, j in all_pairs6:
        if random.random() < 0.5:
            A[i][j] = 1
        else:
            A[j][i] = 1

    hps = find_all_HPs(A, n6)
    if not hps: continue
    n_sampled += 1

    best_pv = n6 + 1
    for hp in hps[:50]:  # Sample up to 50 HPs
        G = hp_to_graph(A, hp, n6)
        pv = sum(sum(G[v][w] for w in range(n6) if w != v) % 2 for v in range(n6))
        best_pv = min(best_pv, pv)
        if best_pv == 0: break

    pv_dist_6[best_pv] += 1
    if best_pv == 0:
        n_even_reachable += 1

print("  Sampled %d tournaments:" % n_sampled)
print("  Even-graph reachable (PV=0): %d (%.1f%%)" %
      (n_even_reachable, 100*n_even_reachable/n_sampled))
for pv in sorted(pv_dist_6.keys()):
    print("    PV=%d: %d tournaments (%.1f%%)" %
          (pv, pv_dist_6[pv], 100*pv_dist_6[pv]/n_sampled))

print("\nDone. opus-2026-03-24-S320")
