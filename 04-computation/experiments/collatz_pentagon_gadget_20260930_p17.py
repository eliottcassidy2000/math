#!/usr/bin/env python3
"""collatz_pentagon_gadget_20260930_p17.py -- the hatted-icosahedron test (Theorem A of the thirteenth note) for P_17 and the
Clebsch graph without building C_5^2 in full: for an induced icosahedron I of C_5(G), its twelve link pentagons are vertices of
C_5^2(G); a hat is an induced pentagon q of C_5(G) sharing an edge with exactly four link pentagons that form a tadpole in I.
Candidate hats are found as the induced pentagons through the edges of the twelve links.  Validated on P_13 (108 hats).
(session collatz-posets-zeta5-20260927, opus, 2026-09-30, direction D39.)
"""
import time, sys
from collatz_pentagon_gadget_20260930 import graph_from_edges, circulant, pentagon_graph, find_icosahedra, is_induced_icosahedron

T0 = time.time()


def paley(q):
    sq = {(x * x) % q for x in range(1, q)}
    return graph_from_edges(q, [(i, j) for i in range(q) for j in range(i + 1, q) if (j - i) % q in sq])


def clebsch():
    return graph_from_edges(16, [(i, j) for i in range(16) for j in range(i + 1, 16) if bin(i ^ j).count("1") in (1, 4)])


def pentagons_through_edge(adj, u, v, limit=200000):
    """induced pentagons containing the edge uv: paths u-v-c-d-e-u with no chords."""
    out = set()
    for c in adj[v]:
        if c == u or c in adj[u]:
            continue
        for d in adj[c]:
            if d in (u, v) or d in adj[u] or d in adj[v]:
                continue
            for e in adj[d]:
                if e in (u, v, c) or e not in adj[u] or e in adj[v] or e in adj[c]:
                    continue
                cyc = (u, v, c, d, e)
                m = min(cyc); i = cyc.index(m); r = cyc[i:] + cyc[:i]
                out.add(min(r, (r[0],) + r[:0:-1]))
                if len(out) > limit:
                    return out
    return out


def cyc_edges(p):
    return set((min(p[i], p[(i + 1) % 5]), max(p[i], p[(i + 1) % 5])) for i in range(5))


def link_pentagon(adj, v, vs):
    N = adj[v] & vs
    start = next(iter(N)); cyc = [start]; prev = None
    while len(cyc) < 5:
        nxt = [u for u in adj[cyc[-1]] & N if u != prev and u not in cyc]
        prev = cyc[-1]; cyc.append(nxt[0])
    return tuple(cyc)


def hats_of_lifted(C1, I):
    vs = set(I)
    links = {v: link_pentagon(C1, v, vs) for v in I}
    link_edges = {v: cyc_edges(links[v]) for v in I}
    link_keys = set()
    for p in links.values():
        m = min(p); i = p.index(m); r = p[i:] + p[:i]
        link_keys.add(min(r, (r[0],) + r[:0:-1]))
    cands = set()
    for v in I:
        for (a, b) in link_edges[v]:
            cands |= pentagons_through_edge(C1, a, b)
    cands -= link_keys
    hats = []
    # adjacency in the lifted icosahedron C_5(I): two link pentagons are adjacent iff they share an edge of C_5(G)
    lifted_adj = {v: {w for w in I if w != v and link_edges[v] & link_edges[w]} for v in I}
    for q in cands:
        qe = cyc_edges(q)
        touched = [v for v in I if qe & link_edges[v]]
        if len(touched) != 4:
            continue
        T = set(touched)
        inner = sum(1 for x in T for y in T if x < y and y in lifted_adj[x])
        degs = sorted(len(lifted_adj[x] & T) for x in T)
        if inner == 4 and degs == [1, 2, 2, 3]:
            hats.append(q)
    return hats, len(cands)


def find_induced_copy(pattern, host, budget=300, avoid=()):
    """an induced copy of `pattern` in `host` as a vertex list (or None), by backtracking in BFS order of the pattern."""
    n = len(pattern); order = []; seen = set()
    for s0 in range(n):
        if s0 in seen:
            continue
        seen.add(s0); q = [s0]
        while q:
            x = q.pop(0); order.append(x)
            for y in sorted(pattern[x]):
                if y not in seen:
                    seen.add(y); q.append(y)
    mapping = {}; used = set(avoid); H = len(host); t1 = time.time()
    def bt(i):
        if time.time() - t1 > budget:
            return None
        if i == n:
            return True
        x = order[i]
        prev = [y for y in pattern[x] if y in mapping]
        cand = set.intersection(*[host[mapping[y]] for y in prev]) if prev else set(range(H))
        for w in sorted(cand):
            if w in used or len(host[w]) < len(pattern[x]):
                continue
            if all(((y in pattern[x]) == (mapping[y] in host[w])) for y in mapping):
                mapping[x] = w; used.add(w)
                r = bt(i + 1)
                if r:
                    return True
                if r is None:
                    return None
                del mapping[x]; used.discard(w)
        return False
    if bt(0):
        return [mapping[x] for x in range(n)]
    return None


def run(name, G, max_icos=6, budget=420, pole_search=True):
    print("== %s ==" % name, flush=True)
    C1, pents = pentagon_graph(G)
    print(" %s: %d induced pentagons; C_5: %d vertices, degrees %d..%d (%.0fs)" % (name, len(pents), len(C1), min(len(a) for a in C1), max(len(a) for a in C1), time.time() - T0), flush=True)
    icos = find_icosahedra(C1, range(len(C1)), max_found=max_icos, time_budget=T0 and (time.time() - T0 + budget)) if pole_search else []
    if not icos:
        from collatz_pentagon_gadget_20260930 import ICO
        copy = find_induced_copy(ICO, C1, budget=budget)
        icos = [copy] if copy else []
        print(" (pole search unavailable; generic induced-subgraph search used)", flush=True)
    print(" induced icosahedra found in C_5(%s): %d (%.0fs)" % (name, len(icos), time.time() - T0), flush=True)
    for I in icos:
        h, nc = hats_of_lifted(C1, I)
        print("  lifted icosahedron of %s: %d candidate pentagons through its link edges, %d tadpole hats%s (%.0fs)" % (
            I[:4], nc, len(h), "  -> EXPANSION CERTIFICATE" if h else "", time.time() - T0), flush=True)
        if h:
            return True
    return False


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    if which in ("all", "p13"):
        run("P_13 (validation; expect 108 hats)", circulant(13, [1, 3, 4, 9, 10, 12]), max_icos=2)
    if which in ("all", "p17"):
        run("P_17", paley(17), max_icos=4)
    if which in ("all", "clebsch"):
        run("Clebsch", clebsch(), max_icos=4, pole_search=False)
    print("total %.0fs" % (time.time() - T0))
