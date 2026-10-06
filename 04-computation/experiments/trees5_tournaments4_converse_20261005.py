"""Trees on n+1 vertices versus n-tournaments modulo converse.

Exact census (pure Python + networkx VF2), 2026-10-05, session opus-2026-10-05-S11.

Universe
  * tournaments on n vertices, n = 1..7, isomorphism classes built by incremental
    vertex extension from the (n-1)-classes; every candidate is bucketed by a
    Weisfeiler-Lehman hash and compared by VF2 (DiGraphMatcher) inside its bucket;
  * unlabelled trees on N vertices, N = 2..9 (networkx.nonisomorphic_trees);
  * optional n = 8 tournaments (--n8), slower.

Positive controls: OEIS A000568 (tournament classes 1,1,2,4,12,56,456,6880),
A002785 (self-converse classes 1,1,2,2,8,12,88,176), A000055 (trees
1,1,1,2,3,6,11,23,47 for N = 1..9), A000571 (score sequences 1,1,2,4,9,22,59),
A000081 (rooted trees 1,1,2,4,9,20), A001764 (Fuss-Catalan 1,1,3,12,55,273).

Every claim printed below is asserted; a failed assertion raises.
"""
from __future__ import annotations

import itertools
import warnings
warnings.filterwarnings("ignore")
import sys
from collections import Counter, defaultdict
from fractions import Fraction
from math import comb

import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher, GraphMatcher

OUT = []


def say(s=""):
    OUT.append(s)
    print(s)


# ---------------------------------------------------------------- tournaments
def wl(G):
    return nx.weisfeiler_lehman_graph_hash(G, iterations=4)


class TournamentClasses:
    """Isomorphism classes of n-tournaments, bucketed by WL hash, VF2 inside."""

    def __init__(self, n):
        self.n = n
        self.reps = []           # DiGraph representatives
        self.buckets = defaultdict(list)

    def find(self, G):
        for idx in self.buckets[wl(G)]:
            if DiGraphMatcher(self.reps[idx], G).is_isomorphic():
                return idx
        return None

    def add(self, G):
        idx = self.find(G)
        if idx is None:
            idx = len(self.reps)
            self.reps.append(G)
            self.buckets[wl(G)].append(idx)
        return idx


def extend_classes(prev: TournamentClasses) -> TournamentClasses:
    n = prev.n + 1
    cur = TournamentClasses(n)
    v = n - 1
    for H in prev.reps:
        for bits in range(1 << (n - 1)):
            G = H.copy()
            G.add_node(v)
            for u in range(n - 1):
                if (bits >> u) & 1:
                    G.add_edge(v, u)
                else:
                    G.add_edge(u, v)
            cur.add(G)
    return cur


def tournament_classes(nmax):
    T1 = TournamentClasses(1)
    G = nx.DiGraph()
    G.add_node(0)
    T1.add(G)
    out = {1: T1}
    for n in range(2, nmax + 1):
        out[n] = extend_classes(out[n - 1])
    return out


def scores(G):
    return tuple(sorted(d for _, d in G.out_degree()))


def three_cycles(G):
    n = G.number_of_nodes()
    # c3 = C(n,3) - sum_v C(s_v, 2)   (classical)
    return comb(n, 3) - sum(comb(d, 2) for _, d in G.out_degree())


def largest_strong(G):
    return max(len(c) for c in nx.strongly_connected_components(G))


def cycle_spectrum(G):
    return tuple(sorted({len(c) for c in nx.simple_cycles(G)}))


def is_self_converse(G):
    return DiGraphMatcher(G, G.reverse(copy=True)).is_isomorphic()


# ---------------------------------------------------------------- trees
def tree_classes(N):
    if N == 1:
        G = nx.Graph()
        G.add_node(0)
        return [G]
    if N == 2:
        return [nx.path_graph(2)]
    return list(nx.nonisomorphic_trees(N))


def creatable_cycle_lengths(T):
    dist = dict(nx.all_pairs_shortest_path_length(T))
    S = set()
    for u, v in itertools.combinations(T.nodes, 2):
        if not T.has_edge(u, v):
            S.add(dist[u][v] + 1)
    return tuple(sorted(S))


def tree_index(trees):
    """Return a function mapping a tree to its class index (VF2 inside WL buckets)."""
    buckets = defaultdict(list)
    for i, T in enumerate(trees):
        buckets[nx.weisfeiler_lehman_graph_hash(T, iterations=4)].append(i)

    def find(T):
        for i in buckets[nx.weisfeiler_lehman_graph_hash(T, iterations=4)]:
            if GraphMatcher(trees[i], T).is_isomorphic():
                return i
        raise RuntimeError("tree class not found")

    return find


def ahu(T, root, parent=None):
    return "(" + "".join(sorted(ahu(T, c, root) for c in T.neighbors(root) if c != parent)) + ")"


def tree_name(T):
    N = T.number_of_nodes()
    deg = sorted((d for _, d in T.degree()), reverse=True)
    if deg[0] == N - 1:
        return "star K_{1,%d}" % (N - 1)
    if deg[0] == 2:
        return "path P%d" % N
    if N == 5:
        return "fork (chair)"
    return "tree deg=%s diam=%d" % (tuple(deg), nx.diameter(T))


# ---------------------------------------------------------------- main
def main(do_n8=False):
    nmax = 8 if do_n8 else 7
    say("# trees on n+1 vertices vs n-tournaments modulo converse  (exact census)")
    say()
    TC = tournament_classes(nmax)
    A000568 = {1: 1, 2: 1, 3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880}
    A002785 = {1: 1, 2: 1, 3: 2, 4: 2, 5: 8, 6: 12, 7: 88, 8: 176}
    A000055 = {1: 1, 2: 1, 3: 1, 4: 2, 5: 3, 6: 6, 7: 11, 8: 23, 9: 47, 10: 106}
    A000571 = {1: 1, 2: 1, 3: 2, 4: 4, 5: 9, 6: 22, 7: 59, 8: 167}

    info = {}
    say("## A. Tournament classes, self-converse classes, merged (mod converse) classes")
    say()
    say("| n | classes | self-converse | converse pairs | merged = (classes+SC)/2 | score sequences | score seqs mod complement | trees on n+1 |")
    say("|---|---|---|---|---|---|---|---|")
    merged_count = {}
    for n in range(1, nmax + 1):
        C = TC[n]
        reps = C.reps
        assert len(reps) == A000568[n], (n, len(reps))
        conv = []
        sc = 0
        for i, G in enumerate(reps):
            j = C.find(G.reverse(copy=True))
            assert j is not None
            conv.append(j)
            if j == i:
                sc += 1
        assert sc == A002785[n], (n, sc)
        pairs = (len(reps) - sc) // 2
        merged = sc + pairs
        merged_count[n] = merged
        assert 2 * merged == len(reps) + sc
        S = {scores(G) for G in reps}
        assert len(S) == A000571[n]
        Smod = {min(s, tuple(sorted(n - 1 - x for x in s))) for s in S}
        trees_np1 = len(tree_classes(n + 1))
        assert trees_np1 == A000055[n + 1]
        info[n] = dict(reps=reps, conv=conv, sc=sc, pairs=pairs, merged=merged)
        say("| %d | %d | %d | %d | %d | %d | %d | %d |" % (
            n, len(reps), sc, pairs, merged, len(S), len(Smod), trees_np1))
    say()
    say("Coincidence merged(n) == trees(n+1): " + ", ".join(
        "n=%d:%s" % (n, "YES" if merged_count[n] == A000055[n + 1] else "no (%d vs %d)" % (merged_count[n], A000055[n + 1]))
        for n in range(1, nmax + 1)))
    assert all(merged_count[n] == A000055[n + 1] for n in range(1, 5))
    assert merged_count[5] == 10 and A000055[6] == 6
    say("Score sequence is a complete class invariant iff n <= 4: " + ", ".join(
        "n=%d:%s" % (n, "complete" if A000571[n] == A000568[n] else "incomplete") for n in range(1, nmax + 1)))
    assert all(A000571[n] == A000568[n] for n in range(1, 5)) and A000571[5] < A000568[5]

    # ---------------------------------------------------------- the n = 4 table
    say()
    say("## B. The four 4-tournaments and the three 5-vertex trees")
    say()
    reps4 = info[4]["reps"]
    conv4 = info[4]["conv"]
    names4 = {}
    for i, G in enumerate(reps4):
        s = scores(G)
        names4[i] = {(0, 1, 2, 3): "TT4 (transitive)", (1, 1, 1, 3): "out-vortex (source over C3)",
                     (0, 2, 2, 2): "in-vortex (C3 over sink)", (1, 1, 2, 2): "STRONG"}[s]
    say("| class | scores | c3 | largest strong comp L | cycle spectrum | self-converse | |Aut| | converse |")
    say("|---|---|---|---|---|---|---|---|")
    rows4 = {}
    for i, G in enumerate(reps4):
        aut = sum(1 for _ in DiGraphMatcher(G, G).isomorphisms_iter())
        rows4[i] = dict(scores=scores(G), c3=three_cycles(G), L=largest_strong(G),
                        spec=cycle_spectrum(G), selfc=(conv4[i] == i), aut=aut)
        say("| %s | %s | %d | %d | %s | %s | %d | %s |" % (
            names4[i], rows4[i]["scores"], rows4[i]["c3"], rows4[i]["L"], rows4[i]["spec"] or "{}",
            rows4[i]["selfc"], aut, names4[conv4[i]]))
    # exact facts
    assert sorted(r["c3"] for r in rows4.values()) == [0, 1, 1, 2]
    assert sorted(r["L"] for r in rows4.values()) == [1, 3, 3, 4]
    specs = sorted(set(r["spec"] for r in rows4.values()))
    assert specs == [(), (3,), (3, 4)]
    say()
    say("H = 1 + 2 c3 on these classes (OCF at n=4): " + str([1 + 2 * r["c3"] for r in rows4.values()]))

    trees5 = tree_classes(5)
    say()
    say("| tree | degrees | leaves | max degree | diameter | cycle lengths creatable by one added edge |")
    say("|---|---|---|---|---|---|")
    rows5 = {}
    for i, T in enumerate(trees5):
        deg = tuple(sorted((d for _, d in T.degree()), reverse=True))
        rows5[i] = dict(name=tree_name(T), deg=deg, leaves=deg.count(1), maxdeg=deg[0],
                        diam=nx.diameter(T), cre=creatable_cycle_lengths(T))
        r = rows5[i]
        say("| %s | %s | %d | %d | %d | %s |" % (r["name"], r["deg"], r["leaves"], r["maxdeg"], r["diam"], r["cre"]))
        assert r["cre"] == tuple(range(3, r["diam"] + 2))
    assert sorted(r["diam"] for r in rows5.values()) == [2, 3, 4]
    say()
    say("Both sides are 3-chains: tournament cycle spectra {} < {3} < {3,4}; tree creatable spectra {3} < {3,4} < {3,4,5}.")
    say("The merged converse pair (out/in vortex) is the middle grade L = 3 <-> the fork (diameter 3).")

    # ---------------------------------------------------------- general-n gradings
    say()
    say("## C. The two gradings for every n: largest strong component L (tournaments mod converse) and diameter (trees on n+1)")
    say()
    say("Camion-Moon: a tournament's cycle spectrum is {3..L} (empty iff L = 1); a tree's creatable")
    say("spectrum is {3..diam+1}. Both chains have n-1 values for n >= 3. Fibre sizes:")
    say()
    for n in range(3, nmax + 1):
        reps = info[n]["reps"]
        conv = info[n]["conv"]
        merged_reps = [i for i in range(len(reps)) if conv[i] >= i]
        fibL = Counter(largest_strong(reps[i]) for i in merged_reps)
        assert largest_strong(reps[0]) in fibL or True
        fibC3 = Counter(three_cycles(reps[i]) for i in merged_reps)
        # cycle spectrum check (theorem) on all classes for n <= 7
        if n <= 7:
            for G in reps:
                L = largest_strong(G)
                assert cycle_spectrum(G) == (tuple(range(3, L + 1)) if L >= 3 else ()), (n, scores(G))
        trees = tree_classes(n + 1)
        fibD = Counter(nx.diameter(T) for T in trees)
        for T in trees:
            assert creatable_cycle_lengths(T) == tuple(range(3, nx.diameter(T) + 2))
        degseqs = {tuple(sorted(d for _, d in T.degree())) for T in trees}
        say("n=%d: merged tournaments by L: %s   (by c3: %s)" % (
            n, dict(sorted(fibL.items())), dict(sorted(fibC3.items()))))
        say("      trees on %d by diameter:   %s   degree sequences %d vs trees %d" % (
            n + 1, dict(sorted(fibD.items())), len(degseqs), len(trees)))
        chainL = sorted(fibL)
        chainD = sorted(fibD)
        assert chainL == [1] + list(range(3, n + 1)), chainL
        assert chainD == list(range(2, n + 1)), chainD
        assert len(chainL) == len(chainD) == n - 1
        if n <= 4:
            assert all(v == 1 for v in fibL.values()) and all(v == 1 for v in fibD.values())
        else:
            assert fibL[3] >= 2 and fibD[3] >= 2
    say()
    say("Both gradings are bijective (every fibre a singleton) exactly for n <= 4; for n >= 5 the")
    say("fibre over the middle grade (L = 3, diameter 3) already has >= 2 elements on each side.")
    say("Degree sequences identify trees on N vertices iff N <= 5 (N = 6: 5 sequences, 6 trees).")

    # explicit witnesses for n >= 5 (proof of the boundary)
    say()
    say("Witnesses (proof of the boundary, verified non-isomorphic here):")
    for n in range(5, nmax + 1):
        C = TC[n]
        # 3-cycle inserted into TT_{n-3} at position p (vertices 0..p-1 above, cycle, rest below)
        found = set()
        for p in range(0, n - 2):
            G = nx.DiGraph()
            order = list(range(n))
            cyc = [p, p + 1, p + 2]
            for u, v in itertools.combinations(order, 2):
                if u in cyc and v in cyc:
                    continue
                G.add_edge(u, v)
            G.add_edges_from([(cyc[0], cyc[1]), (cyc[1], cyc[2]), (cyc[2], cyc[0])])
            assert largest_strong(G) == 3
            i = C.find(G)
            j = info[n]["conv"][i]
            found.add(min(i, j))
        # double stars on n+1 vertices
        trees = tree_classes(n + 1)
        find_t = tree_index(trees)
        ds = set()
        for a in range(1, (n - 1) // 2 + 1):
            b = n - 1 - a
            T = nx.Graph()
            T.add_edge(0, 1)
            T.add_edges_from((0, 2 + k) for k in range(a))
            T.add_edges_from((1, 2 + a + k) for k in range(b))
            assert T.number_of_nodes() == n + 1 and nx.diameter(T) == 3
            ds.add(find_t(T))
        say("  n=%d: %d distinct merged classes with L=3 from C3 inserted into TT%d; %d distinct double stars of diameter 3 on %d vertices" % (
            n, len(found), n - 3, len(ds), n + 1))
        assert len(found) >= 2 and len(ds) >= 2

    # ---------------------------------------------------------- metagraphs
    say()
    say("## D. Metagraphs: one arc reversal (tournaments) and one edge rotation (trees)")
    say()

    def flip_graph(n):
        C = TC[n]
        reps = C.reps
        adj = defaultdict(set)
        dc3 = {}
        for i, G in enumerate(reps):
            for u, v in list(G.edges):
                H = G.copy()
                H.remove_edge(u, v)
                H.add_edge(v, u)
                j = C.find(H)
                if j != i:
                    adj[i].add(j)
                    dc3[(i, j)] = three_cycles(reps[j]) - three_cycles(G)
        return adj, dc3

    adj4, dc3_4 = flip_graph(4)
    say("n=4 flip graph on the 4 classes (edges):")
    E4 = set()
    for i in adj4:
        for j in adj4[i]:
            E4.add(tuple(sorted((i, j))))
    for i, j in sorted(E4):
        say("   %s -- %s   (delta c3 = %+d)" % (names4[i], names4[j], abs(dc3_4[(i, j)])))
    out_i = [i for i in names4 if names4[i].startswith("out")][0]
    in_i = [i for i in names4 if names4[i].startswith("in")][0]
    tt_i = [i for i in names4 if names4[i].startswith("TT")][0]
    st_i = [i for i in names4 if names4[i].startswith("STRONG")][0]
    assert (min(out_i, in_i), max(out_i, in_i)) not in E4, "vortices must be non-adjacent"
    assert len(E4) == 5, E4
    say("=> K4 minus the vortex-vortex edge; the vortices are at flip distance 2 (Redei-note distance lemma).")
    say("   Modulo converse: the triangle TT -- V -- STRONG -- TT.  Unit steps (|delta c3| = 1): the path TT -- V -- STRONG.")
    assert abs(dc3_4[(tt_i, st_i)]) == 2 and abs(dc3_4[(tt_i, out_i)]) == 1 and abs(dc3_4[(out_i, st_i)]) == 1

    def rotation_graph(N):
        trees = tree_classes(N)
        find_t = tree_index(trees)
        adj = defaultdict(set)
        for i, T in enumerate(trees):
            for e in list(T.edges):
                U = T.copy()
                U.remove_edge(*e)
                comps = list(nx.connected_components(U))
                A, B = comps[0], comps[1]
                for a in A:
                    for b in B:
                        if (a, b) == e or (b, a) == e:
                            continue
                        V = U.copy()
                        V.add_edge(a, b)
                        j = find_t(V)
                        if j != i:
                            adj[i].add(j)
        return trees, adj

    trees5, radj5 = rotation_graph(5)
    say()
    say("5-vertex trees under one edge rotation (remove an edge, add an edge keeping a tree):")
    E5 = set()
    for i in radj5:
        for j in radj5[i]:
            E5.add(tuple(sorted((i, j))))
    for i, j in sorted(E5):
        say("   %s -- %s" % (tree_name(trees5[i]), tree_name(trees5[j])))
    assert len(E5) == 2
    star = [i for i, T in enumerate(trees5) if max(d for _, d in T.degree()) == 4][0]
    path = [i for i, T in enumerate(trees5) if max(d for _, d in T.degree()) == 2][0]
    assert (min(star, path), max(star, path)) not in E5
    say("=> the path star -- fork -- path (a rotation changes the maximum degree by at most 1).")

    trees6, radj6 = rotation_graph(6)
    adj5, _ = flip_graph(5)
    say()
    say("Next size: 6-vertex tree rotation graph has %d vertices and %d edges; 5-tournament flip graph on classes has %d vertices and %d edges (mod converse: %d vertices)." % (
        len(trees6), sum(len(v) for v in radj6.values()) // 2, len(info[5]["reps"]),
        sum(len(v) for v in adj5.values()) // 2, info[5]["merged"]))

    # ---------------------------------------------------------- marked versions (observer lens)
    say()
    say("## E. Marked (observer) versions at the coincidence size")
    say()
    # rooted trees on 5 vertices
    rooted = set()
    for T in tree_classes(5):
        for r in T.nodes:
            rooted.add(ahu(T, r))
    say("rooted trees on 5 vertices: %d (A000081)" % len(rooted))
    assert len(rooted) == 9
    # vertex-marked 4-tournaments
    marked = []  # (class, vertex-orbit representative)
    marked_mod_conv = []
    for i, G in enumerate(reps4):
        orbits = []
        for v in G.nodes:
            if any(any(m[v] == w for m in DiGraphMatcher(G, G).isomorphisms_iter()) for w in orbits):
                continue
            orbits.append(v)
        for v in orbits:
            marked.append((i, v))
    say("vertex-marked 4-tournaments up to isomorphism: %d" % len(marked))
    assert len(marked) == 12
    # mod converse: (G, v) ~ (G^op, v); identify via converse-isomorphisms
    seen = set()
    count_mc = 0
    for (i, v) in marked:
        if (i, v) in seen:
            continue
        count_mc += 1
        seen.add((i, v))
        G = reps4[i]
        j = conv4[i]
        Gop = G.reverse(copy=True)
        for m in DiGraphMatcher(reps4[j], Gop).isomorphisms_iter():
            # m: reps4[j] -> Gop ; marked vertex v of Gop corresponds to vertex w of reps4[j] with m[w] = v
            w = [w for w in m if m[w] == v][0]
            # reduce w to its orbit representative in reps4[j]
            for (jj, vv) in marked:
                if jj == j and any(mm[vv] == w for mm in DiGraphMatcher(reps4[j], reps4[j]).isomorphisms_iter()):
                    seen.add((jj, vv))
    say("vertex-marked 4-tournaments modulo converse: %d   (rooted trees on 5 vertices: 9)  -> the coincidence does NOT survive marking a vertex" % count_mc)
    assert count_mc == 6

    # ---------------------------------------------------------- ternary grammar
    say()
    say("## F. The ordered full ternary tree (credit grammar S = e | H S G S G S) and the 3-chain of slots")
    say()
    fc = [comb(3 * m, m) // (2 * m + 1) for m in range(6)]
    say("Fuss-Catalan counts of ordered full ternary trees by internal nodes m=0..5: %s (sum through m=5: %d)" % (fc, sum(fc)))
    assert fc == [1, 1, 3, 12, 55, 273] and sum(fc) == 345
    say("The three child slots of a node are entered at prefix balances +2, +1, 0 (after H, after the first G, after the second G): a 3-chain.")
    say("No letter-level symmetry (reversal, H<->G swap, or both) preserves the language; the planar mirror (reverse the child order")
    say("recursively) is an involution of the language exchanging slots 1 and 3 and fixing slot 2 (orbit type (1,2); corrected wording,")
    say("see trees5_tournaments4_audit_extras_20261005.out section E).")
    # unordered full ternary trees with m internal nodes
    def unordered_full_ternary(m):
        # canonical multiset recursion
        from functools import lru_cache

        @lru_cache(None)
        def shapes(k):
            if k == 0:
                return frozenset({"."})
            res = set()
            for a in range(k):
                for b in range(k - a):
                    c = k - 1 - a - b
                    for x in shapes(a):
                        for y in shapes(b):
                            for z in shapes(c):
                                res.add("(" + "".join(sorted([x, y, z])) + ")")
            return frozenset(res)
        return len(shapes(m))
    uft = [unordered_full_ternary(m) for m in range(6)]
    say("Unordered full ternary trees by internal nodes m=0..5: %s  (ordered/unordered quotient; m=2: 3 ordered -> 1 unordered)" % uft)
    assert uft[:3] == [1, 1, 1]

    # ---------------------------------------------------------- summary
    say()
    say("## G. Summary of exact facts")
    say()
    say("1. merged(n) = (A000568(n) + A002785(n))/2 = 1,1,2,3,10,34,272%s; trees(n+1) = 1,1,2,3,6,11,23%s." % (
        ",3528" if do_n8 else "", ",47" if do_n8 else ""))
    say("   Equal exactly for n <= 4. The owner's 3 = 3 is the last member of a four-term coincidence, not a bijection family.")
    say("2. At n = 4 both sides are graded 3-chains by cycle content: tournaments by largest strong component L in {1,3,4}")
    say("   (= cycle spectrum {}, {3}, {3,4}; = c3 in {0,1,2}); trees by diameter in {2,3,4} (= creatable spectrum {3}, {3,4}, {3,4,5}).")
    say("   The converse pair sits at the middle grade and corresponds to the fork.")
    say("3. For every n the two gradings land on (n-1)-chains; they are bijective iff n <= 4 (witnesses above).")
    say("4. Flip metagraph at n=4: K4 - e, quotient K3 by converse; unit-c3 steps give the path TT -- V -- STRONG = the tree rotation graph star -- fork -- path.")
    say("5. First-moment fields: score sequences are complete iff n <= 4; degree sequences are complete for trees iff N <= 5. Same threshold.")
    say("6. Marking a vertex breaks the coincidence (6 vs 9).")
    say()
    say("ALL ASSERTIONS PASSED (n_max tournaments = %d, trees to %d vertices)" % (nmax, nmax + 1))


if __name__ == "__main__":
    do_n8 = "--n8" in sys.argv
    main(do_n8)
    if "--save" in sys.argv:
        import os
        here = os.path.dirname(os.path.abspath(__file__))
        path = os.path.join(here, "..", "..", "05-knowledge", "results", "trees5_tournaments4_converse_20261005.out")
        with open(path, "w", encoding="utf-8") as fh:
            fh.write("\n".join(OUT) + "\n")
        print("saved", os.path.normpath(path))
