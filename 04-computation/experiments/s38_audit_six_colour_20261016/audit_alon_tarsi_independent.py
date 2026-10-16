#!/usr/bin/env python3
"""S38 audit (independent referee): re-check of the n = 22 models of THM-4661 statement 5.

1. Each stored model (05-knowledge/results/six_colour_exterior_models_n22_20261016.json) is checked against the
   necessary conditions (N1)-(N5) with its stored t and high flags (my own checks, networkx + SAT for (N5)).
2. The Alon-Tarsi coefficient of prod x_v^(d_v) in P_G = prod_(uv in E, u<v by label) (x_u - x_v) is recomputed for
   the stored outdegree targets with an independent implementation: vertices are processed in a breadth-first order
   (not the S38 greedy order), the state is the vector of partial exponents of the vertices that still have
   unprocessed edges, and an edge contributes +1 (exponent to the smaller label) or -1 (to the larger label).  The
   value must coincide with the stored coefficient (the coefficient of a fixed monomial is an invariant of P_G).
3. Sanity of the method on small graphs against brute-force expansion with sympy (K4 minus an edge, C5, the Petersen
   graph's coefficient for a 3-regular orientation target, a random 9-vertex graph).
"""
import json, os, random, itertools, time
import networkx as nx
from pysat.solvers import Solver

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
OUTJ = os.path.join(HERE, 'audit_alon_tarsi_independent.json')


def log(*a):
    print(*a, flush=True)


def at_coefficient(G, d):
    order = list(nx.bfs_tree(G, max(G.nodes(), key=G.degree)).nodes())
    order += [v for v in G.nodes() if v not in order]
    pos = {v: i for i, v in enumerate(order)}
    # process edges grouped by the later endpoint (in BFS order)
    edges = sorted(G.edges(), key=lambda e: (max(pos[e[0]], pos[e[1]]), min(pos[e[0]], pos[e[1]])))
    remaining = {v: G.degree(v) for v in G.nodes()}
    active = []                       # list of vertices with partial exponents in the state
    states = {(): 1}
    for (u, v) in edges:
        for w in (u, v):
            if w not in active:
                active.append(w)
                states = {k + (0,): c for k, c in states.items()}
        iu, iv = active.index(u), active.index(v)
        lo, hi = (u, v) if u < v else (v, u)
        ilo, ihi = active.index(lo), active.index(hi)
        new = {}
        for k, c in states.items():
            if k[ilo] < d[lo]:
                kk = list(k); kk[ilo] += 1; kk = tuple(kk)
                new[kk] = new.get(kk, 0) + c          # factor (x_lo - x_hi): x_lo with +1
            if k[ihi] < d[hi]:
                kk = list(k); kk[ihi] += 1; kk = tuple(kk)
                new[kk] = new.get(kk, 0) - c          # x_hi with -1
        states = {k: c for k, c in new.items() if c}
        for w in (u, v):
            remaining[w] -= 1
            if remaining[w] == 0:
                i = active.index(w)
                states2 = {}
                for k, c in states.items():
                    if k[i] == d[w]:
                        kk = k[:i] + k[i + 1:]
                        states2[kk] = states2.get(kk, 0) + c
                states = {k: c for k, c in states2.items() if c}
                active.pop(i)
        if not states:
            return 0
    return states.get((), 0)


def brute_coefficient(G, d):
    import sympy as sp
    xs = {v: sp.Symbol(f'x{v}') for v in G.nodes()}
    P = 1
    for u, v in G.edges():
        lo, hi = (u, v) if u < v else (v, u)
        P = sp.expand(P * (xs[lo] - xs[hi]))
    mono = 1
    for v in G.nodes():
        mono *= xs[v] ** d[v]
    return sp.Poly(P, *xs.values()).coeff_monomial(mono)


def orientation_target(G, rng):
    out = {v: 0 for v in G.nodes()}
    for u, v in G.edges():
        if rng.random() < 0.5:
            out[u] += 1
        else:
            out[v] += 1
    return out


# ---- 3. sanity on small graphs
rng = random.Random(1)
small = [nx.complete_graph(4), nx.cycle_graph(5), nx.petersen_graph(), nx.gnp_random_graph(9, 0.45, seed=3)]
small[0].remove_edge(0, 1)
san = []
for G in small:
    for _ in range(3):
        d = orientation_target(G, rng)
        a = at_coefficient(G, d); b = brute_coefficient(G, d)
        san.append((a, int(b)))
        assert a == b, (a, b)
log(f'sanity against sympy expansion on small graphs: {len(san)} cases agree, e.g. {san[:6]}')

# ---- 1. and 2.
models = json.load(open(os.path.join(ROOT, '05-knowledge/results/six_colour_exterior_models_n22_20261016.json')))
cert = json.load(open(os.path.join(ROOT, '05-knowledge/results/six_colour_exterior_models_n22_20261016_alon_tarsi.json')))
assert len(models) == len(cert) == 36
res = []
for gi, (mdl, ct) in enumerate(zip(models, cert)):
    G = nx.Graph(); G.add_nodes_from(range(22)); G.add_edges_from(map(tuple, mdl['edges']))
    t, hi = mdl['t'], mdl['high']
    deg = dict(G.degree())
    n1 = sum(nx.triangles(G).values()) == 0 and all(len(set(G[a]) & set(G[b])) <= 2 for a, b in itertools.combinations(G, 2))
    n2 = all(deg[v] >= 6 - t[v] and (hi[v] == (deg[v] >= 7 - t[v])) for v in G) and all(deg[v] <= 7 for v in G)
    low = [v for v in G if not hi[v]]
    L = G.subgraph(low)
    c4 = any(True for cyc in nx.simple_cycles(L.to_directed()) if len(cyc) == 4) if L.number_of_edges() else False
    n3 = (not c4) and (4 * L.number_of_edges() <= 5 * (len(low) - 1) if low else True)
    n4 = any(hi)
    with Solver(name='minisat22') as S:
        vv = lambda v, c: 3 * v + c + 1
        for v in G:
            S.add_clause([vv(v, c) for c in range(3)])
        for a, b in G.edges():
            for c in range(3):
                S.add_clause([-vv(a, c), -vv(b, c)])
        n5 = not S.solve()
    d = {v: ct['certificate']['outdegrees'][v] for v in range(22)}
    assert sum(d.values()) == G.number_of_edges() and max(d.values()) <= 4
    t0 = time.time()
    c = at_coefficient(G, d)
    rec = dict(graph=gi, N1=n1, N2=n2, N3=n3, N4=n4, N5_non3col=n5, edges=G.number_of_edges(),
               coefficient=c, stored=int(ct['certificate']['coefficient']), agree=(c == int(ct['certificate']['coefficient'])),
               seconds=round(time.time() - t0, 1))
    res.append(rec)
    log(json.dumps(rec))
log(f'all models satisfy N1-N5: {all(r["N1"] and r["N2"] and r["N3"] and r["N4"] and r["N5_non3col"] for r in res)}; '
    f'all coefficients reproduced and nonzero: {all(r["agree"] and r["coefficient"] != 0 for r in res)}')
json.dump(dict(sanity=san, models=res), open(OUTJ, 'w'), indent=1)
log('wrote', OUTJ)
