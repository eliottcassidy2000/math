"""Referee: THM-4835 part 3 (end-edge map for trees, root-swap map for functional orientations)
and the corollary F_n identity.  Run: python maps.py"""
import itertools
import random
import time
from math import comb, factorial

import networkx as nx
import numpy as np

from core import automorphisms


# ------------------------------------------------------------ end-edge map (trees)
def end_edge_map(n, edges):
    """Return dict: ordered non-adjacent pair (u,w) -> ordered pair of edges (e_u, e_w)."""
    G = nx.Graph()
    G.add_nodes_from(range(n))
    G.add_edges_from(edges)
    E = {frozenset(e) for e in edges}
    mp = {}
    for u in range(n):
        for w in range(n):
            if u == w or frozenset((u, w)) in E:
                continue
            path = nx.shortest_path(G, u, w)
            assert len(path) >= 3
            mp[(u, w)] = (frozenset(path[:2]), frozenset(path[-2:]))
    return mp


def check_end_edge(n, edges, auts):
    mp = end_edge_map(n, edges)
    E = [frozenset(e) for e in edges]
    # bijection onto ordered pairs of distinct edges, swap compatible
    imgs = set(mp.values())
    target = {(a, b) for a in E for b in E if a != b}
    if imgs != target or len(imgs) != len(mp):
        return False
    for (u, w), (a, b) in mp.items():
        if mp[(w, u)] != (b, a):
            return False
    for g in auts:
        ge = lambda e: frozenset(g[x] for x in e)
        for (u, w), (a, b) in mp.items():
            if mp[(g[u], g[w])] != (ge(a), ge(b)):
                return False
    return True


def tree_checks(nmax=9):
    for n in range(3, nmax + 1):
        t0 = time.time()
        ok = 0
        tot = 0
        maxaut = 0
        for T in nx.nonisomorphic_trees(n):
            edges = list(T.edges())
            arcs = edges + [(b, a) for a, b in edges]  # undirected graph automorphisms
            auts = automorphisms(n, arcs)
            maxaut = max(maxaut, len(auts))
            tot += 1
            ok += check_end_edge(n, edges, auts)
        print(f"trees n={n}: {ok}/{tot} end-edge maps are Aut(G)-equivariant signed bijections "
              f"(max |Aut(G)| = {maxaut}) [{time.time() - t0:.1f}s]", flush=True)


# ---------------------------------------------------- root-swap map (functional)
def root_swap_map(n, r, f):
    """f: dict v -> f(v) for v != r.  Map ordered pairs of V\\{r} to ordered non-edges."""
    E = {frozenset((v, f[v])) for v in f}
    mp = {}
    for a in range(n):
        for b in range(n):
            if a == b or r in (a, b):
                continue
            if frozenset((a, b)) in E:
                if f.get(a) == b and b != r:          # a -> f(a)  <->  r -> a
                    mp[(a, b)] = (r, a)
                elif f.get(b) == a and a != r:        # f(b) <- b, i.e. (f(b), b) <-> (b, r)
                    mp[(a, b)] = (b, r)
                else:
                    raise RuntimeError
            else:
                mp[(a, b)] = (a, b)
    return mp, E


def check_root_swap(n, r, f, auts):
    mp, E = root_swap_map(n, r, f)
    target = {(a, b) for a in range(n) for b in range(n) if a != b and frozenset((a, b)) not in E}
    if set(mp.values()) != target or len(set(mp.values())) != len(mp):
        return False
    for (a, b), (c, d) in mp.items():
        if mp[(b, a)] != (d, c):
            return False
    for g in auts:
        for (a, b), (c, d) in mp.items():
            if mp[(g[a], g[b])] != (g[c], g[d]):
                return False
    return True


def functional_digraphs(n):
    """All labelled (r, f): f: V\\{r} -> V, no loops, no 2-cycles."""
    for r in range(n):
        others = [v for v in range(n) if v != r]
        for img in itertools.product(range(n), repeat=n - 1):
            f = dict(zip(others, img))
            if any(f[v] == v for v in others):
                continue
            if any(f[v] != r and f.get(f[v]) == v for v in others):
                continue
            yield r, f


def F_formula(n):
    m = n - 1
    return n * sum((-1) ** k * factorial(m) * m ** (m - 2 * k) //
                   (factorial(k) * 2 ** k * factorial(m - 2 * k)) for k in range(m // 2 + 1))


def functional_checks():
    for n in range(3, 7):
        t0 = time.time()
        cnt = ok = 0
        for r, f in functional_digraphs(n):
            cnt += 1
            arcs = [(v, f[v]) for v in f]
            auts = automorphisms(n, arcs)
            ok += check_root_swap(n, r, f, auts)
        print(f"functional n={n}: {cnt} labelled (r,f) [formula F_n={F_formula(n)}], "
              f"root-swap equivariant signed bijection in {ok}/{cnt} [{time.time() - t0:.1f}s]",
              flush=True)
    n = 7
    cnt = sum(1 for _ in functional_digraphs(n))
    print(f"functional n=7: {cnt} labelled (r,f), formula F_7 = {F_formula(7)}")


def random_checks(seed=20261027, reps=100):
    rng = random.Random(seed)
    for n in (8, 9, 10):
        okt = okf = 0
        for _ in range(reps):
            # random labelled tree (Pruefer) with random orientation; Aut(D) <= Aut(G)
            T = nx.from_prufer_sequence([rng.randrange(n) for _ in range(n - 2)])
            edges = [(a, b) if rng.random() < 0.5 else (b, a) for a, b in T.edges()]
            auts = automorphisms(n, edges)
            okt += check_end_edge(n, [tuple(e) for e in edges], auts)
            # random functional orientation: root r, f with no loops / 2-cycles
            while True:
                r = rng.randrange(n)
                f = {v: rng.choice([w for w in range(n) if w != v]) for v in range(n) if v != r}
                if not any(f[v] != r and f.get(f[v]) == v for v in f):
                    break
            arcs = [(v, f[v]) for v in f]
            okf += check_root_swap(n, r, f, automorphisms(n, arcs))
        print(f"random n={n}: end-edge ok {okt}/{reps}, root-swap ok {okf}/{reps}", flush=True)


def score_identity():
    for n in range(3, 8):
        t0 = time.time()
        pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
        P = len(pairs)
        codes = np.arange(2 ** P, dtype=np.int64)
        scores = np.zeros((2 ** P, n), dtype=np.int64)
        for k, (i, j) in enumerate(pairs):
            bit = (codes >> k) & 1      # 1: i -> j
            scores[:, i] += bit
            scores[:, j] += 1 - bit
        tot = 0
        for r in range(n):
            prod = np.ones(2 ** P, dtype=object) if n > 7 else np.ones(2 ** P, dtype=np.int64)
            for v in range(n):
                if v != r:
                    prod = prod * scores[:, v]
            tot += int(prod.sum())
        rhs = F_formula(n) * 2 ** comb(n - 1, 2)
        print(f"score identity n={n}: LHS={tot} F_n*2^C(n-1,2)={rhs} equal={tot == rhs} "
              f"F_n={F_formula(n)} [{time.time() - t0:.1f}s]", flush=True)


if __name__ == "__main__":
    tree_checks(9)
    functional_checks()
    random_checks()
    score_identity()
