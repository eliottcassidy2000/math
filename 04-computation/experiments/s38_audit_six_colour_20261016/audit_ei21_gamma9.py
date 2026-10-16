#!/usr/bin/env python3
"""S38 audit (independent referee): THM-4661 statement 4(c) and THM-4662 statement 5.

(a) From the House of Graphs adjacency lists stored in ei21_placement_q_zeta55_20261012.py (EI21, HoG 51378) and
    ei19_fish_placement_q_zeta55_20261012.py (EI17, EI19, fish): check that EI21 contains the 4-cycles
    (10,13,15,14), (13,15,19,17), (14,15,20,16), (15,19,18,20) and the edges 10-18, 16-17 (the Gamma_9 at hub 15).
(b) Count all injective copies (monomorphisms) of Gamma_9 in EI17, EI19, EI21 and the fish, by hub.
(c) Symbolic check of the algebra: with a, b, c, d unit complex numbers, the four rhombus completions and the two
    closing edges give |P|^2 + |Q|^2 = 1, Re(P conj Q) = 0 for P = a - d, Q = b - c, and the three square relations of
    C3 with X = |P|^2 (via (x conj y - conj x y)^2 = 4 Re(x conj y)^2 - 4 |x|^2 |y|^2).
(d) A real realization of Gamma_9 (one-parameter family): numeric construction and check of all 14 unit edges and
    injectivity.
"""
import itertools, json, math, os, re, cmath
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
OUTJ = os.path.join(HERE, 'audit_ei21_gamma9.json')
out = {}


def log(*a):
    print(*a, flush=True)


def load_list(path, name):
    src = open(path).read()
    i = src.index(name + ' = [')
    j = i + len(name + ' = ')
    depth = 0
    for k in range(j, len(src)):
        if src[k] == '[':
            depth += 1
        elif src[k] == ']':
            depth -= 1
            if depth == 0:
                return eval(src[j:k + 1])
    raise ValueError


EI21 = load_list(os.path.join(ROOT, '04-computation/experiments/ei21_placement_q_zeta55_20261012.py'), 'EI21')
src19 = open(os.path.join(ROOT, '04-computation/experiments/ei19_fish_placement_q_zeta55_20261012.py')).read()
blk = src19[src19.index('GRAPHS = {'):]
GR = eval(blk[len('GRAPHS = '):blk.index('}') + 1])


def graph(adj):
    G = nx.Graph()
    for i, nb in enumerate(adj):
        for j in nb:
            assert i in adj[j]
            G.add_edge(i, j)
    return G


G21 = graph(EI21)
log(f'EI21: {G21.number_of_nodes()} vertices, {G21.number_of_edges()} edges, triangles {sum(nx.triangles(G21).values())//3}')
cyc = [(10, 13, 15, 14), (13, 15, 19, 17), (14, 15, 20, 16), (15, 19, 18, 20)]
ok_c = all(G21.has_edge(c[i], c[(i + 1) % 4]) for c in cyc for i in range(4))
ok_e = G21.has_edge(10, 18) and G21.has_edge(16, 17)
ok_h = all(G21.has_edge(15, v) for v in (13, 14, 19, 20))
log(f'EI21 hub-15 structure: four 4-cycles {ok_c}, closing edges 10-18, 16-17 {ok_e}, hub edges {ok_h}')
out['EI21_hub15'] = dict(cycles=ok_c, closing=ok_e, hub=ok_h)

# Gamma_9
g9 = nx.Graph()
g9.add_edges_from([('h', 'a'), ('h', 'b'), ('h', 'c'), ('h', 'd'), ('a', 'wab'), ('b', 'wab'), ('a', 'wac'), ('c', 'wac'),
                   ('b', 'wbd'), ('d', 'wbd'), ('c', 'wcd'), ('d', 'wcd'), ('wab', 'wcd'), ('wac', 'wbd')])
assert g9.number_of_nodes() == 9 and g9.number_of_edges() == 14
census = {}
for name, G in (('EI17', graph(GR['EI17'])), ('EI19', graph(GR['EI19'])), ('EI21', G21), ('fish', graph(GR['FISH']))):
    GM = GraphMatcher(G, g9)
    hubs = {}
    cnt = 0
    for mp in GM.subgraph_monomorphisms_iter():
        inv = {v: k for k, v in mp.items()}
        hubs[inv['h']] = hubs.get(inv['h'], 0) + 1
        cnt += 1
    census[name] = dict(monomorphisms=cnt, hubs=hubs)
    log(f'{name}: {cnt} injective copies of Gamma_9 (labelled), hubs {hubs}')
out['gamma9_census'] = {k: dict(monomorphisms=v['monomorphisms'], hubs={str(a): b for a, b in v['hubs'].items()})
                        for k, v in census.items()}

# (c) symbolic algebra
import sympy as sp
th = sp.symbols('t1:5', real=True)
a, b, c, d = [sp.exp(sp.I * t) for t in th]
P, Q = a - d, b - c
wab, wac, wbd, wcd = a + b, a + c, b + d, c + d
e1 = sp.simplify(sp.expand((wab - wcd) * sp.conjugate(wab - wcd)) - sp.expand((P + Q) * sp.conjugate(P + Q)))
e2 = sp.simplify(sp.expand((wac - wbd) * sp.conjugate(wac - wbd)) - sp.expand((P - Q) * sp.conjugate(P - Q)))
log(f'closing edges are |P+Q| and |P-Q|: residuals {e1}, {e2}')
# from |P+Q|^2 = |P-Q|^2 = 1: |P|^2 + |Q|^2 = 1 and Re(P conj Q) = 0 (linear algebra)
X, Y, R = sp.symbols('X Y R', real=True)       # X = |P|^2, Y = |Q|^2, R = Re(P conj Q)
sol = sp.solve([X + Y + 2 * R - 1, X + Y - 2 * R - 1], [Y, R], dict=True)[0]
log(f'solution: {sol}')
# square relations: (x conj y - conj x y)^2 = 4 Re^2 - 4|x|^2|y|^2 ; for unit a, d with |a-d|^2 = X: Re(a conj d) = 1 - X/2
rel1 = sp.expand(4 * (1 - X / 2) ** 2 - 4)           # (a conj d - conj a d)^2
rel2 = sp.expand(4 * (1 - (1 - X) / 2) ** 2 - 4)     # (b conj c - conj b c)^2, |b-c|^2 = 1 - X
rel0 = sp.expand(4 * 0 - 4 * X * (1 - X))            # (P conj Q - conj P Q)^2
out['C3_relations'] = dict(z0=str(sp.factor(-rel0 / 4 * 11 / 11)), z1=str(sp.factor(-rel1)), z2=str(sp.factor(-rel2)))
log(f'-(PQbar - conj)^2/4 = {sp.factor(-rel0/4)} (11 z0^2 = X(1-X)); -(a dbar - conj)^2 = {sp.factor(-rel1)} '
    f'(11 z1^2 = X(4-X)); -(b cbar - conj)^2 = {sp.factor(-rel2)} (11 z2^2 = (1-X)(X+3))')
assert sp.expand(-rel0 / 4 - X * (1 - X)) == 0 and sp.expand(-rel1 - X * (4 - X)) == 0 and sp.expand(-rel2 - (1 - X) * (X + 3)) == 0

# (d) a real realization: choose X, angle of P; a, d unit with a - d = P; b, c unit with b - c = Q, Q = i s P/|P| ...
def realize(X, alpha, beta):
    # P = sqrt(X) e^{i alpha}; Q = sqrt(1-X) i e^{i alpha}; a, d unit with a - d = P: a = e^{i alpha}(sqrt(X)/2 + i r),
    # d = a - P, r = sqrt(1 - X/4); b, c unit with b - c = Q: b = i e^{i alpha}(sqrt(1-X)/2) + e^{i alpha} s... general:
    Pv = math.sqrt(X) * cmath.exp(1j * alpha)
    Qv = math.sqrt(1 - X) * 1j * cmath.exp(1j * alpha) * (1 if beta >= 0 else -1)
    def unit_pair(D):
        m = abs(D)
        r = math.sqrt(1 - m * m / 4)
        u = D / m
        p1 = u * (m / 2 + 1j * r)
        return p1, p1 - D
    aa, dd = unit_pair(Pv)
    bb, cc = unit_pair(Qv)
    pos = {'h': 0, 'a': aa, 'b': bb, 'c': cc, 'd': dd, 'wab': aa + bb, 'wac': aa + cc, 'wbd': bb + dd, 'wcd': cc + dd}
    errs = max(abs(abs(pos[u] - pos[v]) - 1) for u, v in g9.edges())
    mind = min(abs(pos[u] - pos[v]) for u, v in itertools.combinations(pos, 2))
    return errs, mind


samples = [realize(Xv, al, 1) for Xv in (0.2, 0.37, 0.5, 0.81) for al in (0.0, 0.7)]
log(f'real realizations of Gamma_9: max edge error {max(s[0] for s in samples):.2e}, min vertex distance '
    f'{min(s[1] for s in samples):.3f}')
out['real_realizations'] = dict(max_edge_error=max(s[0] for s in samples), min_vertex_distance=min(s[1] for s in samples))
json.dump(out, open(OUTJ, 'w'), indent=1, default=str)
log('wrote', OUTJ)
