#!/usr/bin/env python3
"""
even_graph_bijection_s20dt.py — Find bijection between tournaments and even graphs
kind-pasteur-2026-03-23-S20dt

A graph G is ODD if there exists an orientation O and automorphism sigma
such that sigma reverses an ODD number of oriented edges.
A graph G is EVEN if it is NOT odd.

Royle et al. (2023): #even graphs on n vertices = #tournaments on n = A000568(n).
NO bijection is known. Let us find one.

Strategy: enumerate both sets at n=3,4,5, compute invariants of each,
and find a natural correspondence.
"""

import sys
from math import comb, factorial
from itertools import permutations, combinations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  FINDING BIJECTION: EVEN GRAPHS <-> TOURNAMENTS")
print("  kind-pasteur-2026-03-23-S20dt")
print("=" * 80)

def analyze(n):
    t0 = time.time()
    ALL_EDGES = [(i,j) for i in range(n) for j in range(i+1,n)]
    m = len(ALL_EDGES)
    perms = list(permutations(range(n)))

    # ============================================================
    # TOURNAMENTS
    # ============================================================
    def make_tourn(bits):
        A = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(ALL_EDGES):
            if bits & (1 << k): A[i][j] = 1
            else: A[j][i] = 1
        return A

    def canon_tourn(A):
        best = None
        for p in perms:
            s = tuple(A[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or s < best: best = s
        return best

    def Hdp(A):
        dp = {}
        for v in range(n): dp[(1<<v, v)] = 1
        for S in range(1, 1<<n):
            for v in range(n):
                if not (S & (1<<v)): continue
                val = dp.get((S,v), 0)
                if val == 0: continue
                for u in range(n):
                    if S & (1<<u): continue
                    if A[v][u]:
                        dp[(S|(1<<u), u)] = dp.get((S|(1<<u), u), 0) + val
        return sum(dp.get(((1<<n)-1, v), 0) for v in range(n))

    tourn_classes = {}
    for bits in range(1 << m):
        A = make_tourn(bits)
        cn = canon_tourn(A)
        if cn not in tourn_classes:
            H = Hdp(A)
            sc = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
            c3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
                     if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]))
            aut = sum(1 for p in perms if all(A[p[i]][p[j]] == A[i][j]
                      for i in range(n) for j in range(n)))
            tourn_classes[cn] = {'H': H, 'scores': sc, 'c3': c3, 'aut': aut,
                                 'orbit': factorial(n)//aut}

    # ============================================================
    # EVEN GRAPHS
    # ============================================================
    def graph_from_bits(bits):
        """Undirected graph as edge set."""
        edges = set()
        for k, (i,j) in enumerate(ALL_EDGES):
            if bits & (1 << k):
                edges.add((i,j))
        return frozenset(edges)

    def graph_adj(bits):
        A = [[0]*n for _ in range(n)]
        for k, (i,j) in enumerate(ALL_EDGES):
            if bits & (1 << k):
                A[i][j] = 1; A[j][i] = 1
        return A

    def canon_graph(bits):
        best = None
        for p in perms:
            nb = 0
            for k, (i,j) in enumerate(ALL_EDGES):
                pi, pj = min(p[i],p[j]), max(p[i],p[j])
                nk = ALL_EDGES.index((pi, pj))
                if (bits >> k) & 1: nb |= (1 << nk)
            if best is None or nb < best: best = nb
        return best

    def graph_auts(bits):
        """Find all automorphisms of graph defined by bits."""
        A = graph_adj(bits)
        return [p for p in perms if all(A[p[i]][p[j]] == A[i][j]
                for i in range(n) for j in range(n))]

    def is_odd_graph(bits):
        """Check if graph is odd: exists orientation + aut reversing odd # edges."""
        edge_list = [(i,j) for k, (i,j) in enumerate(ALL_EDGES) if bits & (1 << k)]
        if not edge_list:
            return False  # empty graph is even

        auts = graph_auts(bits)

        # Try all orientations (2^|edges| of them)
        num_edges = len(edge_list)
        if num_edges > 15:  # skip if too many
            return None

        for orient_bits in range(1 << num_edges):
            # Build oriented edge list
            oriented = []
            for k, (i,j) in enumerate(edge_list):
                if orient_bits & (1 << k):
                    oriented.append((j, i))  # reverse
                else:
                    oriented.append((i, j))  # forward

            for sigma in auts:
                if sigma == list(range(n)):
                    continue  # skip identity
                # Count reversed edges
                reversed_count = 0
                for u, v in oriented:
                    su, sv = sigma[u], sigma[v]
                    # The oriented edge u->v maps to su->sv.
                    # Find which original edge this corresponds to.
                    orig_edge = (min(su,sv), max(su,sv))
                    if orig_edge not in [(min(a,b),max(a,b)) for a,b in edge_list]:
                        continue  # shouldn't happen if sigma is automorphism
                    # What's the orientation of this edge in the original?
                    idx = edge_list.index(orig_edge)
                    orig_u, orig_v = oriented[idx]
                    # sigma maps u->v to su->sv. The edge is orig_edge.
                    # Direction preserved if su->sv matches orig orientation.
                    if su == orig_u and sv == orig_v:
                        pass  # preserved
                    elif su == orig_v and sv == orig_u:
                        reversed_count += 1  # reversed
                    else:
                        # sigma maps to the same edge but different endpoints
                        # This means (su,sv) = (orig_v, orig_u), so reversed
                        reversed_count += 1

                if reversed_count % 2 == 1:
                    return True  # odd!

        return False  # even

    # Classify all graphs
    graph_classes = {}
    for bits in range(1 << m):
        cn = canon_graph(bits)
        if cn in graph_classes:
            continue
        A = graph_adj(bits)
        deg_seq = tuple(sorted(sum(A[i][j] for j in range(n)) for i in range(n)))
        num_edges = bin(bits).count('1')
        aut = len(graph_auts(bits))
        odd = is_odd_graph(bits)
        graph_classes[cn] = {'bits': bits, 'deg': deg_seq, 'edges': num_edges,
                            'aut': aut, 'odd': odd, 'even': not odd}

    even_graphs = {cn: info for cn, info in graph_classes.items() if info['even']}
    odd_graphs = {cn: info for cn, info in graph_classes.items() if info['odd']}

    print(f"\nn={n}: {len(graph_classes)} graphs, {len(even_graphs)} even, {len(odd_graphs)} odd")
    print(f"       {len(tourn_classes)} tournaments")
    print(f"       Match? {len(even_graphs) == len(tourn_classes)}")
    print(f"       Time: {time.time()-t0:.1f}s")

    # Print both sets
    print(f"\n  EVEN GRAPHS:")
    print(f"  {'edges':>5} {'deg_seq':>20} {'|Aut|':>5} {'orbit':>5}")
    for cn in sorted(even_graphs.keys(), key=lambda c: even_graphs[c]['edges']):
        info = even_graphs[cn]
        print(f"  {info['edges']:5d} {str(info['deg']):>20} {info['aut']:5d} {factorial(n)//info['aut']:5d}")

    print(f"\n  TOURNAMENTS:")
    print(f"  {'H':>5} {'scores':>20} {'|Aut|':>5} {'orbit':>5} {'c3':>3}")
    for cn in sorted(tourn_classes.keys(), key=lambda c: tourn_classes[c]['H']):
        info = tourn_classes[cn]
        print(f"  {info['H']:5d} {str(info['scores']):>20} {info['aut']:5d} {info['orbit']:5d} {info['c3']:3d}")

    # Look for matching invariants
    print(f"\n  POTENTIAL BIJECTION INVARIANTS:")
    even_auts = sorted([even_graphs[cn]['aut'] for cn in even_graphs])
    tourn_auts = sorted([tourn_classes[cn]['aut'] for cn in tourn_classes])
    print(f"  Even graph |Aut| values: {even_auts}")
    print(f"  Tournament |Aut| values: {tourn_auts}")
    print(f"  |Aut| multisets match? {even_auts == tourn_auts}")

    even_orbits = sorted([factorial(n)//even_graphs[cn]['aut'] for cn in even_graphs])
    tourn_orbits = sorted([tourn_classes[cn]['orbit'] for cn in tourn_classes])
    print(f"  Even orbit sizes: {even_orbits}")
    print(f"  Tournament orbit sizes: {tourn_orbits}")
    print(f"  Orbit multisets match? {even_orbits == tourn_orbits}")

    return even_graphs, tourn_classes


# ============================================================================
# MAIN
# ============================================================================

for n in [3, 4, 5]:
    analyze(n)

print(f"\n  DONE.")
print("=" * 80)
