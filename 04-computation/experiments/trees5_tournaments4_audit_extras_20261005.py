"""Session-lead re-verification of the workflow auditors' new findings (opus-2026-10-05-S11).

Independent of the lens/auditor scratch code; reuses only the class builder of
trees5_tournaments4_converse_20261005.py (same folder).

Checks
  A. centre-marked trees with the centre-swap involution: (marked, fixed, orbits) for N = 4..8
     versus (A000568, A002785, merged) at n = N-1  (expect equality iff N <= 5)
  B. arc-flip label law  delta c3 = s_u - s_v - 1  on every (class, arc), n <= 7
  C. hostiles: no-source-no-sink non-strong classes, every-vertex-on-a-3-cycle non-strong classes,
     every-arc-on-a-3-cycle => strong, n <= 7 (and n = 8 with --n8)
  D. palindromic slack path (self-complementary score sequence) versus self-converse: first class-level
     failure at n = 6 (12 of 24), none at n <= 5
  E. credit grammar S = e | H S G S G S: no letter-level symmetry preserves the language (m <= 5);
     the planar mirror is an involution of the language with slot permutation (1 3)(2)
  F. converse-invariant orientations of the three 5-vertex trees: star 1, fork 0, path 2
  G. unit-step restriction is one-sided: merged flip metagraph at n = 4 is K3, every edge rotation of
     a 5-tree changes the diameter by exactly 1

Run:  python -X utf8 trees5_tournaments4_audit_extras_20261005.py [--n8] [--save]
"""
from __future__ import annotations

import itertools
import os
import sys
import warnings

warnings.filterwarnings("ignore")
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher, GraphMatcher

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from trees5_tournaments4_converse_20261005 import (  # noqa: E402
    tournament_classes, scores, three_cycles, largest_strong, is_self_converse, tree_classes, ahu)

OUT = []


def say(s=""):
    OUT.append(s)
    print(s)


A000568 = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}
A002785 = {1: 1, 2: 1, 3: 2, 4: 2, 5: 8, 6: 12, 7: 88, 8: 176}


def main(do_n8):
    nmax = 8 if do_n8 else 7
    TC = tournament_classes(nmax)

    # ------------------------------------------------------------------ A
    say("## A. Centre-marked trees with the centre-swap involution")
    say()
    say("| N | trees | centre-marked | swap-fixed | orbits | n=N-1: classes | self-converse | merged | equal |")
    say("|---|---|---|---|---|---|---|---|---|")
    resA = {}
    for N in range(4, 9):
        trees = tree_classes(N)
        marked = fixed = 0
        for T in trees:
            centres = nx.center(T)
            if len(centres) == 1:
                marked += 1
                fixed += 1
            else:
                c1, c2 = centres
                # marks at c1 and c2 are equivalent iff an automorphism exchanges them
                same = ahu(T, c1) == ahu(T, c2)
                if same:
                    marked += 1
                    fixed += 1
                else:
                    marked += 2
        orbits = (marked + fixed) // 2
        n = N - 1
        eq = (marked, fixed, orbits) == (A000568[n], A002785[n], (A000568[n] + A002785[n]) // 2)
        resA[N] = eq
        say("| %d | %d | %d | %d | %d | %d | %d | %d | %s |" % (N, len(trees), marked, fixed, orbits, A000568[n], A002785[n], (A000568[n] + A002785[n]) // 2, "YES" if eq else "no"))
        assert orbits == len(trees)
    assert resA[4] and resA[5] and not resA[6] and not resA[7]
    say()
    say("Centre-marked trees reproduce (classes, self-converse, merged) exactly for N <= 5 (n <= 4) and fail from N = 6.")
    say("At N = 5: star and path have a unique centre (fixed); the fork is bicentral with inequivalent centres (a swapped pair).")

    # ------------------------------------------------------------------ B
    say()
    say("## B. Arc-flip label law delta c3 = s_u - s_v - 1")
    say()
    checked = 0
    for n in range(3, min(nmax, 7) + 1):
        for G in TC[n].reps:
            c = three_cycles(G)
            for u, v in list(G.edges):
                H = G.copy()
                H.remove_edge(u, v)
                H.add_edge(v, u)
                assert three_cycles(H) - c == G.out_degree(u) - G.out_degree(v) - 1
                checked += 1
    say("verified on %d (class, arc) pairs, n <= %d. Proof: c3 = C(n,3) - sum C(s_w,2); the flip changes s_u -> s_u - 1, s_v -> s_v + 1, so" % (checked, min(nmax, 7)))
    say("delta c3 = C(s_u,2) + C(s_v,2) - C(s_u-1,2) - C(s_v+1,2) = (s_u - 1) - s_v.")

    # ------------------------------------------------------------------ C
    say()
    say("## C. Hostiles: positive local readings without strong connectivity")
    say()
    say("| n | no source, no sink, not strong | ... mod converse | every vertex on a 3-cycle, not strong | every arc on a 3-cycle, not strong |")
    say("|---|---|---|---|---|")
    rowsC = {}
    for n in range(3, nmax + 1):
        C = TC[n]
        conv = {}
        nss = []
        v3 = 0
        a3 = 0
        for i, G in enumerate(C.reps):
            if largest_strong(G) == n:
                continue
            od = dict(G.out_degree())
            if min(od.values()) >= 1 and max(od.values()) <= n - 2:
                nss.append(i)
            # vertex on a 3-cycle
            on3 = {w: False for w in G.nodes}
            arc_on3 = {e: False for e in G.edges}
            for a, b, c in itertools.permutations(G.nodes, 3):
                if a < b and a < c and G.has_edge(a, b) and G.has_edge(b, c) and G.has_edge(c, a):
                    on3[a] = on3[b] = on3[c] = True
                    arc_on3[(a, b)] = arc_on3[(b, c)] = arc_on3[(c, a)] = True
            if all(on3.values()):
                v3 += 1
            if all(arc_on3.values()):
                a3 += 1
        # mod converse
        seen = set()
        nss_mc = 0
        for i in nss:
            if i in seen:
                continue
            nss_mc += 1
            j = C.find(C.reps[i].reverse(copy=True))
            seen.add(i)
            seen.add(j)
        rowsC[n] = (len(nss), nss_mc, v3, a3)
        say("| %d | %d | %d | %d | %d |" % (n, len(nss), nss_mc, v3, a3))
        assert a3 == 0
    assert [rowsC[n][0] for n in range(3, 8)] == [0, 0, 0, 1, 3]
    assert [rowsC[n][2] for n in range(3, 8)] == [0, 0, 0, 1, 2]
    if do_n8:
        say("n = 8 census: no-source-no-sink non-strong = %d (Theorem A of the synthesis predicted 16), every-vertex-on-3-cycle non-strong = %d (predicted 13)" % (rowsC[8][0], rowsC[8][2]))
    say()
    say("Minimal order 6, unique witness C3 => C3 (scores 1,1,1,4,4,4). 'Every arc on a 3-cycle' forces strong for all n: an arc between")
    say("distinct strong components lies on no cycle, so all arcs are internal and there is one component. (PROVED; census agrees.)")

    # ------------------------------------------------------------------ D
    say()
    say("## D. Palindromic slack (self-complementary score sequence) versus self-converse, class level")
    say()
    say("| n | classes with self-complementary scores | of which self-converse | failures |")
    say("|---|---|---|---|")
    resD = {}
    for n in range(3, min(nmax, 7) + 1):
        pal = sc = 0
        first = None
        for G in TC[n].reps:
            s = scores(G)
            if tuple(sorted(n - 1 - x for x in s)) == s:
                pal += 1
                if is_self_converse(G):
                    sc += 1
                elif first is None:
                    first = (s, sorted(G.edges))
        resD[n] = (pal, sc, first)
        say("| %d | %d | %d | %d |" % (n, pal, sc, pal - sc))
    assert all(resD[n][0] == resD[n][1] for n in range(3, 6))
    assert resD[6][0] - resD[6][1] == 12 and resD[6][0] == 24
    say()
    say("First class-level failure at n = 6: scores %s, arcs %s" % (resD[6][2][0], resD[6][2][1]))
    say("(self-converse => palindromic slack always; the converse holds at the class level iff n <= 5; at the score level for all n by Eplett 1979, CITED)")

    # ------------------------------------------------------------------ E
    say()
    say("## E. Credit grammar S = e | H S G S G S: letter-level symmetries and the planar mirror")
    say()

    def balanced_words(m):
        # all words with m H (+2) and 2m G (-1), nonnegative prefix sums, total 0
        res = []
        for pos in itertools.combinations(range(3 * m), m):
            w = ["G"] * (3 * m)
            for p in pos:
                w[p] = "H"
            bal = 0
            ok = True
            for ch in w:
                bal += 2 if ch == "H" else -1
                if bal < 0:
                    ok = False
                    break
            if ok and bal == 0:
                res.append("".join(w))
        return res

    def parse(w):
        # returns nested tuple (S1, S2, S3) or () for empty, using first returns
        if not w:
            return ()
        assert w[0] == "H"
        bal = 2
        i = 1
        start = 1
        parts = []
        while len(parts) < 2:
            if w[i] == "G" and bal == (2 if len(parts) == 0 else 1):
                parts.append(w[start:i])
                start = i + 1
                bal -= 1
            else:
                bal += 2 if w[i] == "H" else -1
            i += 1
        parts.append(w[start:])
        return tuple(parse(p) for p in parts)

    def unparse(t):
        if t == ():
            return ""
        return "H" + unparse(t[0]) + "G" + unparse(t[1]) + "G" + unparse(t[2])

    def mirror(t):
        if t == ():
            return ()
        return (mirror(t[2]), mirror(t[1]), mirror(t[0]))

    from math import comb
    for m in range(1, 6):
        W = balanced_words(m)
        assert len(W) == comb(3 * m, m) // (2 * m + 1)
        Wset = set(W)
        rev_ok = sum(1 for w in W if w[::-1] in Wset)
        swap_ok = sum(1 for w in W if w.translate(str.maketrans("HG", "GH")) in Wset)
        both_ok = sum(1 for w in W if w[::-1].translate(str.maketrans("HG", "GH")) in Wset)
        mir = [unparse(mirror(parse(w))) for w in W]
        assert all(x in Wset for x in mir)
        assert all(unparse(mirror(mirror(parse(w)))) == w for w in W)
        fixed = sum(1 for w, x in zip(W, mir) if w == x)
        say("m=%d: %d words; reversal preserves %d, H<->G swap preserves %d, both %d; planar mirror is an involution of the language with %d fixed trees" % (
            m, len(W), rev_ok, swap_ok, both_ok, fixed))
        assert rev_ok == 0 and swap_ok == 0 and both_ok == 0
    say("Mirror-fixed counts 1,1,2,3,7 for m=1..5 (A047749 'symmetric ternary trees' starts 1,1,1,2,3,7 with m=0). The mirror permutes the")
    say("three slots by (1 3)(2) and reverses the entry-balance chain +2,+1,0: orbit type (1,2), like the Berggren leg swap.")

    # ------------------------------------------------------------------ F
    say()
    say("## F. Converse-invariant orientations of the three 5-vertex trees")
    say()
    for T in tree_classes(5):
        E = list(T.edges)
        classes = []
        for bits in range(1 << len(E)):
            D = nx.DiGraph()
            D.add_nodes_from(T.nodes)
            for k, (a, b) in enumerate(E):
                D.add_edge(a, b) if (bits >> k) & 1 else D.add_edge(b, a)
            if any(DiGraphMatcher(D, X).is_isomorphic() for X in classes):
                continue
            classes.append(D)
        inv = sum(1 for D in classes if DiGraphMatcher(D, D.reverse(copy=True)).is_isomorphic())
        deg = tuple(sorted((d for _, d in T.degree()), reverse=True))
        name = {4: "star", 3: "fork", 2: "path"}[deg[0]]
        say("%s: %d orientation classes, %d converse-invariant, %d modulo converse" % (name, len(classes), inv, (len(classes) + inv) // 2))
        if name == "fork":
            assert inv == 0
        if name == "star":
            assert inv == 1
        if name == "path":
            assert inv == 2
    say("The fork, the tree at the merged (converse-pair) grade, is the only 5-vertex tree with no converse-invariant orientation:")
    say("its unique degree-3 vertex has out-degree in {0,1,2,3}, and reversal sends out-degree k to 3-k, never equal (PROVED).")

    # ------------------------------------------------------------------ G
    say()
    say("## G. The unit-step restriction is one-sided")
    say()
    C4 = TC[4]
    adj = set()
    for i, G in enumerate(C4.reps):
        for u, v in list(G.edges):
            H = G.copy()
            H.remove_edge(u, v)
            H.add_edge(v, u)
            j = C4.find(H)
            if j != i:
                adj.add((min(i, j), max(i, j)))
    conv = [C4.find(G.reverse(copy=True)) for G in C4.reps]
    merged_edges = {tuple(sorted((min(i, conv[i]), min(j, conv[j])))) for i, j in adj if min(i, conv[i]) != min(j, conv[j])}
    say("merged flip metagraph at n = 4: %d nodes, %d edges (K3)" % (len({min(i, conv[i]) for i in range(4)}), len(merged_edges)))
    assert len(merged_edges) == 3
    trees5 = tree_classes(5)
    deltas = set()
    for T in trees5:
        for e in list(T.edges):
            U = T.copy()
            U.remove_edge(*e)
            A, B = list(nx.connected_components(U))
            for a in A:
                for b in B:
                    if {a, b} == set(e):
                        continue
                    V = U.copy()
                    V.add_edge(a, b)
                    if any(GraphMatcher(V, T).is_isomorphic() for _ in [0]):
                        continue
                    deltas.add(abs(nx.diameter(V) - nx.diameter(T)))
    say("5-vertex trees: every class-changing edge rotation changes the diameter by exactly %s" % sorted(deltas))
    assert deltas == {1}
    say("So 'P3 = P3' needs the unit-c3 restriction on the tournament side only; the owner's merged metagraph is the triangle K3.")

    say()
    say("ALL ASSERTIONS PASSED (n_max = %d)" % nmax)
    if "--save" in sys.argv:
        here = os.path.dirname(os.path.abspath(__file__))
        path = os.path.normpath(os.path.join(here, "..", "..", "05-knowledge", "results", "trees5_tournaments4_audit_extras_20261005.out"))
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(OUT) + "\n")
        print("saved", path)


if __name__ == "__main__":
    main("--n8" in sys.argv)
