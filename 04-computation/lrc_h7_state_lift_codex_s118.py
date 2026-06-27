#!/usr/bin/env python3
"""S118: the H=7 / binary-arc transfer target for LRC14.

This is a small exact guardrail script for the owner's digraph prompt.

The useful facts are:

1. For a connected simple graph G, I(G,2)=7 forces G=K3.
2. For tournaments, H(T)=I(Omega(T),2), and the K3 conflict atom is forbidden.
3. For arbitrary directed graphs, H=7 is not forbidden; the binary state must
   be a complete orientation state, not merely edge present/absent.

The LRC14 proof route therefore needs a state-lift theorem from an apex-7
over-cover atom to a tournament-realizable conflict graph.  If the lifted graph
has connected partition function 7, it is K3, hence forbidden.
"""

from collections import Counter, defaultdict
from itertools import combinations


def independence_at_2(n, edges):
    edge_set = {tuple(sorted(e)) for e in edges}
    total = 0
    for mask in range(1 << n):
        ok = True
        size = mask.bit_count()
        bits = [i for i in range(n) if (mask >> i) & 1]
        for a, b in combinations(bits, 2):
            if (a, b) in edge_set:
                ok = False
                break
        if ok:
            total += 2**size
    return total


def connected(n, edges):
    if n == 0:
        return False
    adj = [[] for _ in range(n)]
    for a, b in edges:
        adj[a].append(b)
        adj[b].append(a)
    seen = {0}
    stack = [0]
    while stack:
        v = stack.pop()
        for w in adj[v]:
            if w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) == n


def graph_edges_from_mask(n, mask):
    pairs = list(combinations(range(n), 2))
    return [pairs[i] for i in range(len(pairs)) if (mask >> i) & 1]


def audit_independence_preimage(max_n=6):
    by_n = {}
    witnesses = []
    for n in range(1, max_n + 1):
        conn_count = 0
        val_counts = Counter()
        for mask in range(1 << (n * (n - 1) // 2)):
            edges = graph_edges_from_mask(n, mask)
            if not connected(n, edges):
                continue
            conn_count += 1
            val = independence_at_2(n, edges)
            val_counts[val] += 1
            if val == 7:
                witnesses.append((n, tuple(edges)))
        by_n[n] = (conn_count, val_counts.get(7, 0))
    return by_n, witnesses


def hamiltonian_path_count(adj):
    n = len(adj)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            cur = dp[mask][last]
            if cur == 0:
                continue
            for nxt in range(n):
                if (mask >> nxt) & 1:
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += cur
    return sum(dp[full])


def tournament_adj_from_mask(n, mask):
    adj = [[False] * n for _ in range(n)]
    bit = 0
    for i, j in combinations(range(n), 2):
        if (mask >> bit) & 1:
            adj[i][j] = True
        else:
            adj[j][i] = True
        bit += 1
    return adj


def audit_tournament_H(max_n=6):
    spectra = {}
    counts = {}
    for n in range(1, max_n + 1):
        c = Counter()
        for mask in range(1 << (n * (n - 1) // 2)):
            h = hamiltonian_path_count(tournament_adj_from_mask(n, mask))
            c[h] += 1
        spectra[n] = sorted(c)
        counts[n] = c
    return spectra, counts


def digraph_adj_from_mask(n, mask):
    arcs = [(i, j) for i in range(n) for j in range(n) if i != j]
    adj = [[False] * n for _ in range(n)]
    for bit, (i, j) in enumerate(arcs):
        if (mask >> bit) & 1:
            adj[i][j] = True
    return adj, arcs


def find_arbitrary_digraph_with_H(target=7, n=4):
    for mask in range(1 << (n * (n - 1))):
        adj, arcs = digraph_adj_from_mask(n, mask)
        if hamiltonian_path_count(adj) == target:
            edges = [(i, j) for bit, (i, j) in enumerate(arcs) if (mask >> bit) & 1]
            return edges
    return None


def main():
    print("S118 H=7 / BINARY-ARC STATE-LIFT AUDIT")
    print("=" * 72)
    print()

    print("1. Connected graph preimage of I(G,2)=7")
    by_n, witnesses = audit_independence_preimage()
    for n, (conn_count, seven_count) in by_n.items():
        print(f"   n={n}: connected labelled graphs={conn_count}, I(G,2)=7 count={seven_count}")
    print(f"   witnesses: {witnesses}")
    print("   proof shortcut: if n>=4, I(G,2)>=1+2n>=9;")
    print("   for n=3, the only connected graph with no independent pair is K3.")
    print()

    print("2. Tournament Hamiltonian-path spectra for complete binary orientations")
    spectra, counts = audit_tournament_H()
    for n, vals in spectra.items():
        seven = counts[n].get(7, 0)
        print(f"   n={n}: H-values={vals}; H=7 labelled count={seven}")
    print("   This finite audit is only a checksum; the project theorem is the all-n")
    print("   obstruction H(T)=I(Omega(T),2), with the K3 conflict atom forbidden.")
    print()

    print("3. Guardrail: arbitrary digraphs do realize H=7")
    edges = find_arbitrary_digraph_with_H()
    print(f"   first 4-vertex present/absent digraph with H=7 has arcs={edges}")
    print("   Therefore the LRC transfer cannot use 'digraph' in the loose sense.")
    print("   It needs a complete binary orientation / tournament-realizable packet.")
    print()

    print("4. LRC14 state-lift theorem shape")
    print("   If a primitive top-balanced LRC14 counterexample yields a connected")
    print("   binary packet graph G_LRC with partition function I(G_LRC,2)=7, then")
    print("   G_LRC=K3.  If the packet graph is a tournament odd-cycle conflict")
    print("   graph, THM-201/THM-343 forbid it as a whole atom or component.")
    print("   So that lift would prove the counterexample impossible.  The missing")
    print("   work is the lift from apex-7 over-cover cells to the")
    print("   tournament-realizable conflict packet, not the H=7 obstruction.")
    print("   Caveat: K3 is not forbidden as an arbitrary subgraph; n=8 realizes")
    print("   H=63 with Omega=K31.")


if __name__ == "__main__":
    main()
