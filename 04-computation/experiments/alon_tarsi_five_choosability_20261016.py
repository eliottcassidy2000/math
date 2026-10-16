#!/usr/bin/env python3
"""Alon-Tarsi certificates of five-choosability.

For a graph G and an outdegree target d (d_v <= 4, sum d_v = |E|), the Combinatorial Nullstellensatz gives:
if the coefficient of prod_v x_v^(d_v) in the graph polynomial P_G = prod_(uv in E, u<v) (x_u - x_v) is nonzero,
then G is L-colourable for every list assignment with |L(v)| >= d_v + 1, in particular five-choosable.  The
coefficient is computed exactly (Python integers) by a frontier dynamic programme over the edges: each edge is
oriented towards the endpoint whose exponent it raises (sign -1 when it raises x_v for u < v); a vertex leaves the
frontier after its last edge, and only states where it has used exactly d_v survive.

A failed six-colour exterior is not five-choosable (THM-4661 2(c)), so a nonzero coefficient certifies that G is not
a failed exterior.  Session opus-2026-10-16-S38.
Usage: python alon_tarsi_five_choosability_20261016.py graphs.json [tries]   (graphs.json: list of {"edges": ...})
"""
import sys, json, random, time
import networkx as nx


def frontier_order(G):
    """Vertex order with a small frontier (greedy: next vertex with most neighbours already placed)."""
    order = []
    placed = set()
    start = min(G.nodes(), key=lambda v: G.degree(v))
    order.append(start); placed.add(start)
    while len(order) < G.number_of_nodes():
        cand = max((v for v in G.nodes() if v not in placed),
                   key=lambda v: (sum(1 for u in G[v] if u in placed), -G.degree(v)))
        order.append(cand); placed.add(cand)
    return order


def coefficient(G, d, order):
    pos = {v: i for i, v in enumerate(order)}
    edges = sorted(((min(u, v, key=pos.get), max(u, v, key=pos.get)) for u, v in G.edges()),
                   key=lambda e: (pos[e[1]], pos[e[0]]))
    last = {}
    for i, (u, v) in enumerate(edges):
        last[u] = i; last[v] = i
    # state: tuple of (vertex, used) for active vertices, sorted by vertex -> coefficient
    states = {(): 1}
    active = []
    maxfront = 0
    for i, (u, v) in enumerate(edges):
        for w in (u, v):
            if w not in active:
                active.append(w)
        idx = {w: k for k, w in enumerate(active)}
        new = {}
        sgn_u = 1 if u < v else -1        # term x_u has sign +1 iff u is the smaller label in (x_min - x_max)
        for st, c in states.items():
            used = list(st) + [0] * (len(active) - len(st))
            # orient raising x_u
            if used[idx[u]] < d[u]:
                nu = used.copy(); nu[idx[u]] += 1
                key = tuple(nu)
                new[key] = new.get(key, 0) + c * sgn_u
            if used[idx[v]] < d[v]:
                nv = used.copy(); nv[idx[v]] += 1
                key = tuple(nv)
                new[key] = new.get(key, 0) - c * sgn_u
        # retire vertices whose last edge is i
        for w in (u, v):
            if last[w] == i:
                k = active.index(w)
                filt = {}
                for st, c in new.items():
                    if st[k] == d[w] and c != 0:
                        key = st[:k] + st[k + 1:]
                        filt[key] = filt.get(key, 0) + c
                new = filt
                active.pop(k)
        states = {k: c for k, c in new.items() if c != 0}
        maxfront = max(maxfront, len(states))
        if not states:
            return 0, maxfront
    return states.get((), 0), maxfront


def random_target(G, cap=4, rng=random):
    """Random outdegree sequence realisable by an orientation, all outdegrees <= cap (orientation by local search)."""
    E = list(G.edges())
    out = {v: 0 for v in G.nodes()}
    ori = []
    for u, v in E:
        if rng.random() < 0.5:
            ori.append((u, v)); out[u] += 1
        else:
            ori.append((v, u)); out[v] += 1
    for _ in range(100000):
        bad = [i for i, (a, b) in enumerate(ori) if out[a] > cap]
        if not bad:
            return out
        i = rng.choice(bad)
        a, b = ori[i]
        ori[i] = (b, a); out[a] -= 1; out[b] += 1
    return None


if __name__ == '__main__':
    graphs = json.load(open(sys.argv[1]))
    tries = int(sys.argv[2]) if len(sys.argv) > 2 else 20
    rng = random.Random(20261016)
    results = []
    for gi, g in enumerate(graphs):
        G = nx.Graph(); G.add_edges_from(map(tuple, g['edges']))
        order = frontier_order(G)
        found = None
        t0 = time.time()
        for k in range(tries):
            d = random_target(G, 4, rng)
            if d is None:
                continue
            c, mf = coefficient(G, d, order)
            if c != 0:
                found = dict(outdegrees=[d[v] for v in sorted(G.nodes())], coefficient=str(c), max_states=mf)
                break
        results.append(dict(graph=gi, n=G.number_of_nodes(), m=G.number_of_edges(), certificate=found,
                            tries=k + 1, seconds=round(time.time() - t0, 1)))
        print(f'graph {gi}: n={G.number_of_nodes()} m={G.number_of_edges()}: '
              f'{"5-choosable (AT coefficient " + found["coefficient"] + ")" if found else "no certificate"} '
              f'after {k+1} targets, {time.time()-t0:.1f} s', flush=True)
    json.dump(results, open(sys.argv[1].replace('.json', '_alon_tarsi.json'), 'w'), indent=1)
