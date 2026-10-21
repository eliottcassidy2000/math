"""tribonacci_exception_check_20261020.py -- opus-2026-10-20-S47
Identify the classes at the ends of closed blue lines without a three-move black detour
(blue_line_black_detours_20261020.cpp exceptions), by VF2 isomorphism against the tribonacci node
T*_n (transitive tournament with its Hamiltonian path reversed: arcs i->i+1 and j->i for j >= i+2)
and against the transitive tournament.  Usage: python3 tribonacci_exception_check_20261020.py n tiling [n tiling ...]
"""
import sys
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher

def cells(n):
    return [(i, j) for i in range(n) for j in range(i + 2, n)]

def tournament(n, t):
    G = nx.DiGraph(); G.add_nodes_from(range(n))
    for i in range(n - 1): G.add_edge(i, i + 1)
    for k, (i, j) in enumerate(cells(n)):
        if t >> k & 1: G.add_edge(i, j)
        else: G.add_edge(j, i)
    return G

def kappa(n, t):
    return t ^ ((1 << len(cells(n))) - 1)

def iso(G, H):
    return DiGraphMatcher(G, H).is_isomorphic()

def main():
    args = list(map(int, sys.argv[1:]))
    for n, t in zip(args[::2], args[1::2]):
        T = tournament(n, t); K = tournament(n, kappa(n, t))
        star = tournament(n, 0); trans = tournament(n, (1 << len(cells(n))) - 1)
        H = lambda G: sum(1 for _ in ()) if False else None
        sc = lambda G: sorted(d for _, d in G.out_degree())
        print(f"n={n} tiling={t}: closed member scores {sc(T)} iso T*_n: {iso(T, star)}, iso T_n: {iso(T, trans)}; "
              f"open partner scores {sc(K)} iso T*_n: {iso(K, star)}, iso T_n: {iso(K, trans)}; "
              f"partner is converse of closed member: {iso(K, T.reverse())}; closed member self-converse: {iso(T, T.reverse())}")

if __name__ == "__main__":
    main()
