"""Independent audit (S70 audit) of THM-4852 (2), (4) and the note's sections 4-6: fraction anatomy over 2..29, the
register-flow graph against the n = 5 merged tiling line metagraph (parsed directly from THM-4729's n = 5 table),
alternative graphs, and baselines for the numerology.  Written from scratch.
"""
import glob
import os
import re
from fractions import Fraction
from itertools import combinations, permutations
from math import comb, lcm, prod

import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..', '..'))
P = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
OWNER = [(17, 91), (78, 85), (19, 51), (23, 38), (29, 33), (77, 29), (95, 23), (77, 19), (1, 17), (11, 13), (13, 11),
         (15, 2), (1, 7), (55, 1)]            # Conway 1987 Theorem 1 = Wikipedia/OEIS = the owner's list
GUY = OWNER[:11] + [(15, 14), (15, 2), (55, 1)]  # Guy 1983


def say(*a):
    print(*a, flush=True)


def fac(n):
    e = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            e[d] = e.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        e[n] = e.get(n, 0) + 1
    return e


def defect(N):
    """repository convention (aliquot_66_five_fifteen_20261028/PROOFS.md): F = tau-2, U = omega-[prime],
    S = 2^omega - 1 - [squarefree]; returns F - U - S."""
    e = fac(N)
    tau = prod(v + 1 for v in e.values())
    om = len(e)
    prime = om == 1 and sum(e.values()) == 1
    sqf = all(v == 1 for v in e.values())
    return (tau - 2) - (om - prime) - (2 ** om - 1 - sqf)


def shape(N):
    e = sorted(fac(N).values(), reverse=True)
    return '1' if N == 1 else '*'.join(('p', 'q', 'r', 's')[i] + (f'^{v}' if v > 1 else '') for i, v in enumerate(e))


# ------------------------------------------------------------------ anatomy
say('=' * 100)
say('ANATOMY over 2..29 (owner = Conway 1987 list; Guy 1983 list differs in rules L, M)')
say('=' * 100)
# balanced shapes by brute force
bal = sorted({shape(N) for N in range(2, 200000) if defect(N) == 0})
say(f'  balanced shapes (F = U + S) among 2 <= N < 200000: {bal}')
assert bal == ['p', 'p*q^2*r', 'p^2*q*r', 'p^3'] or set(bal) == {'p', 'p^3', 'p^2*q*r'}
for name, L in (('owner', OWNER), ('Guy', GUY)):
    say(f'  -- {name} list')
    net = {p: 0 for p in P}
    for k, (a, b) in enumerate(L):
        ea, eb = fac(a), fac(b)
        ch = []
        for p in P:
            d = ea.get(p, 0) - eb.get(p, 0)
            net[p] += d
            if d:
                ch.append(f'{p}{"+" if d > 0 else "-"}')
        da = '-' if a == 1 else defect(a)
        db = '-' if b == 1 else defect(b)
        say(f'    {chr(65 + k)} {a:>3}/{b:<3} shapes {shape(a):>7} / {shape(b):<7} change {", ".join(ch):<22} defect {da} / {db}')
    pr = prod(Fraction(a, b) for a, b in L)
    say(f'    net exponent change {dict((p, v) for p, v in net.items() if v)}; product of the 14 fractions = {pr} = {fac(pr.numerator)} / {fac(pr.denominator)}')
    say(f'    lcm(denominators) = {lcm(*[b for a, b in L])}; pure fractions {[(a, b) for a, b in L if a == 1 or b == 1]}')
    parts = sorted(({a for a, b in L} | {b for a, b in L}) - {1})
    say(f'    distinct parts > 1: {parts}; balanced ones: {[x for x in parts if defect(x) == 0]}')
po = prod(Fraction(a, b) for a, b in OWNER)
pg = prod(Fraction(a, b) for a, b in GUY)
assert po == Fraction(55 ** 2, 2 * 17 ** 2) and pg / po == Fraction(15, 2)
assert lcm(*[b for a, b in OWNER]) == lcm(*[b for a, b in GUY]) == prod(P) == 6469693230
say(f'  owner product 55^2/(2*17^2) = {po}; Guy/owner = {pg / po}; lcm = 29# = {prod(P)} (both)')
say(f'  owner set {{7, 17, 55}} vs balanced shapes {{p, p^3, p^2qr}}: defects 7 -> {defect(7)}, 17 -> {defect(17)}, 55 -> {defect(55)} (55 = pq is NOT balanced)')

# ------------------------------------------------------------------ flow graphs
say('')
say('=' * 100)
say('REGISTER-FLOW GRAPH (arc p -> q: p | den, q | num of one fraction) and alternatives')
say('=' * 100)


def flow(L):
    arcs = set()
    for a, b in L:
        for p in P:
            if b % p == 0:
                for q in P:
                    if a % q == 0:
                        arcs.add((p, q))
    return arcs


def cooc(L):
    G = nx.Graph()
    G.add_nodes_from(P)
    for a, b in L:
        S = [p for p in P if (a * b) % p == 0]
        G.add_edges_from(combinations(S, 2))
    return G


def same_side(L):
    G = nx.Graph()
    G.add_nodes_from(P)
    for a, b in L:
        for side in (a, b):
            S = [p for p in P if side % p == 0]
            G.add_edges_from(combinations(S, 2))
    return G


graphs = {}
for name, L in (('owner', OWNER), ('Guy', GUY)):
    arcs = flow(L)
    D = nx.DiGraph()
    D.add_nodes_from(P)
    D.add_edges_from(arcs)
    G = nx.Graph(D)
    loops = [e for e in arcs if e[0] == e[1]]
    say(f'  flow[{name}]: {len(arcs)} arcs, {G.number_of_edges()} simple pairs, loops {loops}; '
        f'(out,in) = {[(p, D.out_degree(p), D.in_degree(p)) for p in P]}')
    say(f'      simple degree sequence {sorted(d for _, d in G.degree())}')
    graphs[f'flow[{name}]'] = G
    for gname, H in ((f'cooccurrence[{name}]', cooc(L)), (f'same-side[{name}]', same_side(L))):
        say(f'  {gname}: {H.number_of_edges()} edges, degree sequence {sorted(d for _, d in H.degree())}')
        graphs[gname] = H

# ------------------------------------------------------------------ metagraph from THM-4729
say('')
say('=' * 100)
say('n = 5 MERGED TILING LINE METAGRAPH, parsed from THM-4729 (n = 5 table)')
say('=' * 100)
f4729 = glob.glob(os.path.join(ROOT, '01-canon', 'theorems', 'THM-4729-*.md'))[0]
txt = open(f4729, encoding='utf-8').read()
mb = re.search(r'^blue : (.*)$', txt, re.M).group(1)
mk = re.search(r'^black: (.*)$', txt, re.M).group(1)
blue = [tuple(int(v) for v in e.split('-')) for e in mb.split(', ')]
black = []
for e in mk.split(', '):
    ed, mult = e.split(' x')
    black += [tuple(int(v) for v in ed.split('-'))] * int(mult)
tab = re.findall(r'^\s+(\d)\s+\|\s+(\w+)\s+\|\s+(\w+)\s+\|\s+(\d+)\s+\|\s+(\d+)\s+\|\s+(\d+)\s+\|', txt, re.M)
nodes = {int(r[0]): dict(type=r[1], strong=r[2], H=int(r[3]), Aut=int(r[4]), t=int(r[5])) for r in tab[:10]}
say(f'  parsed {len(blue)} blue and {len(black)} black lines; node (H, Aut, t): {[(k, v["H"], v["Aut"], v["t"]) for k, v in sorted(nodes.items())]}')
tl = sorted(v['t'] for v in nodes.values())
say(f'  tilings per node sorted {tl}, sum {sum(tl)} (= 2^6 = 64 tilings, 32 lines = {len(blue) + len(black)})')
M = nx.MultiGraph()
M.add_nodes_from(range(10))
M.add_edges_from(blue + black)
S = nx.Graph()
S.add_nodes_from(range(10))
S.add_edges_from(e for e in blue + black if e[0] != e[1])
loops = sorted({e for e in black + blue if e[0] == e[1]})
pure_blue = [v for v in range(10) if not any(v in e for e in black)]
say(f'  simple non-loop pairs {S.number_of_edges()}, loops {loops}, pure-blue nodes {pure_blue}')
say(f'  simple degree sequence {sorted(d for _, d in S.degree())}; multigraph degrees (loop = 2) {sorted(d for _, d in M.degree())}')
assert len(blue) == 8 and len(black) == 24 and S.number_of_edges() == 18 and pure_blue == [0, 3, 9]
# blue/black overlap on the same pair?
ov = {frozenset(e) for e in blue} & {frozenset(e) for e in black}
say(f'  pairs carrying both a blue and a black line: {sorted(tuple(sorted(e)) for e in ov)}')

say('')
say('  bijection tests against the metagraph simple graph S (18 edges):')
for gname, G in graphs.items():
    Gr = nx.relabel_nodes(G, {p: i for i, p in enumerate(P)})
    iso = nx.is_isomorphic(Gr, S)
    m1 = GraphMatcher(S, Gr).subgraph_is_monomorphic() if Gr.number_of_edges() <= S.number_of_edges() else False
    m2 = GraphMatcher(Gr, S).subgraph_is_monomorphic() if S.number_of_edges() <= Gr.number_of_edges() else False
    say(f'    {gname:<24} edges {G.number_of_edges():>2}: isomorphic {iso}; {gname} embeds in S {m1}; S embeds in {gname} {m2}')

# weak (homomorphic) bijections into the co-occurrence graph: how cheap are they?
C = nx.relabel_nodes(graphs['cooccurrence[owner]'], {p: i for i, p in enumerate(P)})
n_emb, n_owner = 0, 0
for m in GraphMatcher(C, S).subgraph_monomorphisms_iter():   # m: C-node -> S-node, S edges map onto C edges
    inv = {v: P[k] for k, v in m.items()}                     # metagraph node -> prime
    n_emb += 1
    if {inv[0], inv[3], inv[9]} in ({7, 17, 5}, {7, 17, 11}):
        n_owner += 1
say(f'  adjacency-preserving bijections metagraph -> co-occurrence graph: {n_emb}; of these {n_owner} send the '
    f'pure-blue nodes {{TT5, 1->C3->1, T5}} onto {{7, 17, 5 or 11}} (the owner\'s 7, 17, 55)')

# a graph on the FRACTIONS (the owner's {7, 17, 55} are fraction parts): rule-succession graph of the orbit
# (Conway's flowchart), 14 nodes; merging Conway's loop pairs AB, DG, EF and the inverse pair JK gives 10 macro-rules
def succession(L, nphase):
    num = [fac(a) for a, b in L]
    den = [fac(b) for a, b in L]
    x = {2: 1}
    prev, arcs, entries = None, set(), 0
    while entries <= nphase:
        if not any(x.get(s, 0) for s in (11, 13, 17, 19, 23, 29)) and prev is not None and prev != 11 and prev != 12:
            entries += 1
        for k in range(14):
            if all(x.get(p, 0) >= e for p, e in den[k].items()):
                break
        for p, e in den[k].items():
            x[p] -= e
        for p, e in num[k].items():
            x[p] = x.get(p, 0) + e
        if prev is not None:
            arcs.add((prev, k))
        prev = k
    return arcs


succ = succession(OWNER, 120)
assert succ == succession(OWNER, 40)          # the arc set is already complete after 40 phases
Gs = nx.Graph()
Gs.add_nodes_from(range(14))
Gs.add_edges_from((a, b) for a, b in succ if a != b)
say(f'  rule-succession graph (owner; same arc set after 40 and 120 phases): {len(succ)} arcs on 14 rules, {Gs.number_of_edges()} simple pairs, '
    f'self-loops at {sorted(chr(65 + a) for a, b in succ if a == b)}')
macro = {0: 'AB', 1: 'AB', 2: 'C', 3: 'DG', 6: 'DG', 4: 'EF', 5: 'EF', 7: 'H', 8: 'I', 9: 'JK', 10: 'JK', 11: 'L', 12: 'M', 13: 'N'}
Gm = nx.Graph()
Gm.add_nodes_from(sorted(set(macro.values())))
Gm.add_edges_from((macro[a], macro[b]) for a, b in succ if macro[a] != macro[b])
say(f'  merged into 10 macro-rules {sorted(Gm.nodes)}: {Gm.number_of_edges()} pairs, degrees {sorted(d for _, d in Gm.degree())}; '
    f'pure macro-rules I, M, N; isomorphic to S: {nx.is_isomorphic(Gm, S)}; '
    f'embeds in S: {GraphMatcher(S, Gm).subgraph_is_monomorphic()}; S embeds in it: {GraphMatcher(Gm, S).subgraph_is_monomorphic()}')

# how special is "no adjacency bijection"? random graphs with the flow graph's edge count
import random
rng = random.Random(5)
hits = 0
for _ in range(2000):
    R = nx.gnm_random_graph(10, 18, seed=rng.randrange(10 ** 9))
    hits += nx.is_isomorphic(R, S)
say(f'  baseline: random G(10, 18) graphs isomorphic to S: {hits} / 2000 (a negative bijection test carries no information)')

# T*_5 = node 7?  R_5 = base path 5->4->3->2->1 plus every tile oriented upward (THM-4729 L4 / THM-787)
R5 = nx.DiGraph()
R5.add_nodes_from(range(1, 6))
for a in range(1, 6):
    for b in range(1, a):
        if a == b + 1:
            R5.add_edge(a, b)      # path arc a -> a-1
        else:
            R5.add_edge(b, a)      # tile oriented upward
Ts = nx.DiGraph()
Ts.add_nodes_from(range(5))
for i in range(5):
    for j in range(i + 1, 5):
        Ts.add_edge(i, j)
Ts.remove_edge(0, 2); Ts.add_edge(2, 0)
Ts.remove_edge(2, 4); Ts.add_edge(4, 2)


def ham_paths(D):
    return sum(1 for perm in permutations(D.nodes) if all(D.has_edge(perm[i], perm[i + 1]) for i in range(len(perm) - 1)))


def n_aut(D):
    return sum(1 for _ in nx.algorithms.isomorphism.DiGraphMatcher(D, D).isomorphisms_iter())


say(f'  node 7 check: R_5 (end of TT5\'s blue line 0-7) has (H, |Aut|) = ({ham_paths(R5)}, {n_aut(R5)}); '
    f'T*_5 (THM-4732) has ({ham_paths(Ts)}, {n_aut(Ts)}); isomorphic {nx.is_isomorphic(R5, Ts)}')

# ------------------------------------------------------------------ numerology baselines
say('')
say('=' * 100)
say('NUMEROLOGY BASELINES')
say('=' * 100)
cyc = [17, 25, 37, 55, 41, 61, 91]
parts_all = {a for a, b in OWNER} | {b for a, b in OWNER}
odd_range = list(range(17, 92, 2))
odd_parts = [x for x in odd_range if x in parts_all]
hit = [x for x in cyc if x in parts_all]
Nn, K, n = len(odd_range), len(odd_parts), len(cyc)
pge = sum(comb(K, j) * comb(Nn - K, n - j) for j in range(len(hit), n + 1)) / comb(Nn, n)
say(f'  3x-1 cycle {sorted(cyc)} hits PRIMEGAME parts {hit}; odd parts in [17,91]: {odd_parts} ({K} of {Nn} odd numbers)')
say(f'  P(a random 7-set of odd numbers in [17,91] contains >= {len(hit)} parts) = {pge:.3f}  (expected hits {n * K / Nn:.2f})')
OGG = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]
assert prod(OGG) == prod(p for p in OGG if p <= 31) * 41 * 47 * 59 * 71
say(f'  Ogg primes {OGG}: registers = first ten; unused {OGG[10:]}; product of the 15 = {prod(OGG)}')
say(f'  any finite prime set S works as a base for THM-4852 (3) (B = prod S); "P#" and "Ogg" are just two choices')
say('DONE')
