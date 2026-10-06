"""Trees with an oriented diametral path as the double cover carrying the tournament converse,
and the Berggren leg swap on the three branches.

Session opus-2026-10-05-S11. Companion of trees5_tournaments4_converse_20261005.py.

For a tree Y let OD(Y) be the set of oriented diametral paths of Y modulo Aut(Y), and rev the
reversal. Prints, for N = 4..8, the number of oriented classes, how many are reversal-fixed and
how many swapped pairs, next to the tournament numbers (classes, self-converse, pairs) at n = N-1.
Expected: Z/2-set isomorphism at N = 4 (2 fixed) and N = 5 (2 fixed + 1 pair), failure from N = 6.

Also verifies with sympy that the leg swap S = (a b) conjugates the Berggren matrices by the
transposition (A C) and fixes B.

Run:  python -X utf8 trees5_tournaments4_oriented_diametral_20261005.py [--save]
"""
from __future__ import annotations

import os
import sys
import warnings

warnings.filterwarnings("ignore")
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
import sympy as sp

OUT = []


def say(s=""):
    OUT.append(s)
    print(s)


def oriented_diametral_classes(T):
    D = nx.diameter(T)
    dist = dict(nx.all_pairs_shortest_path_length(T))
    paths = [tuple(nx.shortest_path(T, u, v)) for u in T.nodes for v in T.nodes if dist[u][v] == D]
    auts = list(GraphMatcher(T, T).isomorphisms_iter())
    orbits, seen = [], set()
    for p in paths:
        if p in seen:
            continue
        orb = {tuple(a[x] for x in p) for a in auts}
        seen |= orb
        orbits.append(frozenset(orb))
    idx = {q: i for i, orb in enumerate(orbits) for q in orb}
    pairing = {i: idx[tuple(reversed(next(iter(orb))))] for i, orb in enumerate(orbits)}
    return orbits, pairing


def tree_label(T):
    N = T.number_of_nodes()
    deg = tuple(sorted((d for _, d in T.degree()), reverse=True))
    if deg[0] == N - 1:
        return "star"
    if deg[0] == 2:
        return "path"
    if N == 5:
        return "fork"
    return "deg=%s" % (deg,)


TOURN = {3: (2, 2, 0), 4: (4, 2, 1), 5: (12, 8, 2), 6: (56, 12, 22), 7: (456, 88, 184)}


def main():
    say("# Oriented diametral paths: the double cover of trees that carries the converse")
    say()
    say("| N | tree classes | oriented classes | reversal-fixed | swapped pairs | n=N-1 tournament classes | self-converse | converse pairs | Z/2-sets isomorphic |")
    say("|---|---|---|---|---|---|---|---|---|")
    detail = []
    results = {}
    for N in range(4, 9):
        total = fixed = pairs = 0
        trees = list(nx.nonisomorphic_trees(N))
        for T in trees:
            orbits, pairing = oriented_diametral_classes(T)
            k = len(orbits)
            f = sum(1 for i in pairing if pairing[i] == i)
            total += k
            fixed += f
            pairs += (k - f) // 2
            if N <= 6:
                detail.append("  N=%d %-16s diam=%d  oriented=%d fixed=%d pairs=%d" % (N, tree_label(T), nx.diameter(T), k, f, (k - f) // 2))
        tc, sc, pr = TOURN[N - 1]
        iso = (total, fixed, pairs) == (tc, sc, pr)
        results[N] = (total, fixed, pairs, iso)
        say("| %d | %d | %d | %d | %d | %d | %d | %d | %s |" % (N, len(trees), total, fixed, pairs, tc, sc, pr, "YES" if iso else "no"))
    say()
    for line in detail:
        say(line)
    assert results[4][:3] == (2, 2, 0) and results[4][3]
    assert results[5][:3] == (4, 2, 1) and results[5][3]
    assert not results[6][3] and not results[7][3]
    say()
    say("At N = 5: star and path are reversal-fixed, the fork has two orientations exchanged by reversal:")
    say("exactly the pattern (TT4 fixed, STRONG fixed, vortex pair swapped) of the four 4-tournaments under converse.")
    say("At N = 6, 7, 8 the Z/2-sets differ from the 5-, 6-, 7-tournaments (counts above).")
    say()
    # Berggren
    A = sp.Matrix([[1, -2, 2], [2, -1, 2], [2, -2, 3]])
    B = sp.Matrix([[1, 2, 2], [2, 1, 2], [2, 2, 3]])
    C = sp.Matrix([[-1, 2, 2], [-2, 1, 2], [-2, 2, 3]])
    S = sp.Matrix([[0, 1, 0], [1, 0, 0], [0, 0, 1]])
    ok = (S * A * S == C) and (S * B * S == B) and (S * C * S == A)
    assert ok
    say("Berggren matrices A, B, C and the leg swap S = (a b): S A S = C, S B S = B, S C S = A  -> %s" % ok)
    say("The leg swap acts on the three Berggren branches as the transposition (A C) fixing B: one fixed branch and one swapped pair,")
    say("the opposite pattern to the tournament converse on the merged grades (trivial) and on the four classes (two fixed, one pair).")
    v = sp.Matrix([3, 4, 5])
    say("(3,4,5) children: A -> %s  B -> %s  C -> %s" % ((A * v).T.tolist()[0], (B * v).T.tolist()[0], (C * v).T.tolist()[0]))
    say()
    say("ALL ASSERTIONS PASSED")
    if "--save" in sys.argv:
        here = os.path.dirname(os.path.abspath(__file__))
        path = os.path.normpath(os.path.join(here, "..", "..", "05-knowledge", "results", "trees5_tournaments4_oriented_diametral_20261005.out"))
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(OUT) + "\n")
        print("saved", path)


if __name__ == "__main__":
    main()
