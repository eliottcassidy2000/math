#!/usr/bin/env python3
"""How small can a vertex-minimal failed six-colour exterior (THM-4661) be?  A SAT search over the combinatorial
necessary conditions, rooted at a high vertex of maximum degree.

Necessary conditions on such an exterior H with n vertices (THM-4661 statement 2), all used here:
  (N1) triangle-free; any two vertices have at most two common neighbours;
  (N2) every vertex b has list size 6 - t(b), t(b) in {0,1}, and deg(b) >= 6 - t(b); b is low if deg(b) = 6 - t(b),
       high if deg(b) >= 7 - t(b); the counting lemma (2) gives deg(b) <= D = floor((n-1)/3);
  (N3) every component of the low subgraph is a Gallai tree whose blocks are edges or odd cycles of length >= 5;
       used through two consequences: no 4-cycle inside the low subgraph, and 4|E(low)| <= 5(#low - 1) when some
       vertex is low;
  (N4) some vertex is high;
  (N5) H is not three-colourable (THM-4661 2(c)); imposed lazily: each three-colouring of a candidate is cut off by
       the clause "some pair coloured alike by it is an edge" (counterexample-guided refinement).
Root: z = a high vertex of maximum degree d among high vertices (so every high vertex has degree <= d; d in [6, D];
low vertices have degree |L| <= 6, so every vertex has degree <= d).
N(z) = d independent vertices; every neighbour u of z has degree in [5, D]; the second neighbourhood W (m vertices)
consists of vertices adjacent to one or two vertices of N(z), no two of them adjacent to the same pair; the
remaining n - 1 - d - m vertices are far (distance >= 3).  The attachment pattern of W on N(z) (a simple "pair
graph" on N(z) plus single attachments) is enumerated up to isomorphism (pair graphs from the networkx graph atlas,
single attachments up to the pair graph's automorphisms); the edges inside W and among far vertices and between them
are SAT variables.

UNSAT for every pattern means that no n-vertex graph satisfies (N1)-(N5): every failed exterior has more than n
vertices.  Verdicts are CaDiCaL verdicts (FINITE-EXACT as solver results; not independently proof-checked).
A model satisfying (N1)-(N5) would not by itself be a failed exterior; it is reported for a separate list test.
Session opus-2026-10-16-S38.  Usage: python six_colour_exterior_twenty_sat_20261016.py N [N2 ...]
"""
import sys, itertools, json, time
import networkx as nx
from networkx.algorithms.isomorphism import GraphMatcher
from pysat.formula import IDPool
from pysat.card import CardEnc, EncType
from pysat.solvers import Solver

ATLAS = nx.graph_atlas_g()


def patterns(n, d, D):
    """Attachment patterns (m, pairs, singles) for a root of degree d, up to isomorphism."""
    out = []
    gs = [g for g in ATLAS if g.number_of_nodes() == d]
    for m in range(0, n - d):
        if 1 + d + m > n:
            break
        for g in gs:
            p = g.number_of_edges()
            s = m - p
            if s < 0:
                continue
            pdeg = dict(g.degree())
            auts = [mp for mp in GraphMatcher(g, g).isomorphisms_iter()]
            seen = set()
            for S in itertools.combinations_with_replacement(range(d), s):
                att = {u: pdeg[u] for u in range(d)}
                for x in S:
                    att[x] += 1
                if any(not (4 <= att[u] <= d - 1) for u in range(d)):
                    continue
                key = min(tuple(sorted(mp[x] for x in S)) for mp in auts)
                if key in seen:
                    continue
                seen.add(key)
                out.append((m, [list(e) for e in g.edges()], list(key)))
    return out


def three_colour(G):
    with Solver(name='cadical195') as S:
        vv = lambda v, c: 3 * v + c + 1
        for v in G.nodes():
            S.add_clause([vv(v, c) for c in range(3)])
        for x, y in G.edges():
            for c in range(3):
                S.add_clause([-vv(x, c), -vv(y, c)])
        if not S.solve():
            return None
        m = S.get_model()
        return {v: next(c for c in range(3) if m[vv(v, c) - 1] > 0) for v in G.nodes()}


def solve_pattern(n, d, D, m, pairs, singles, budget_rounds=20000):
    Z = 0
    NZ = list(range(1, d + 1))
    W = list(range(d + 1, d + 1 + m))
    F = list(range(d + 1 + m, n))
    pool = IDPool()
    cls = []
    fixed = {}
    for v in range(1, n):
        fixed[frozenset((Z, v))] = v in NZ
    for u, w in itertools.combinations(NZ, 2):
        fixed[frozenset((u, w))] = False
    att = {w: [] for w in W}
    for idx, (x, y) in enumerate(pairs):
        att[W[idx]] = [NZ[x], NZ[y]]
    for k, s in enumerate(singles):
        att[W[len(pairs) + k]] = [NZ[s]]
    for u in NZ:
        for w in W:
            fixed[frozenset((u, w))] = u in att[w]
        for f in F:
            fixed[frozenset((u, f))] = False
    var = {}
    rest = W + F
    for x, y in itertools.combinations(rest, 2):
        var[frozenset((x, y))] = pool.id(('e', x, y))
    for f in F:                                   # far vertices are not adjacent to... (no constraint beyond N(z))
        pass

    def lit(x, y):
        k = frozenset((x, y))
        if k in fixed:
            return fixed[k]
        return var[k]

    V = list(range(n))
    for x, y, zz in itertools.combinations(V, 3):
        ls = [lit(x, y), lit(y, zz), lit(x, zz)]
        if any(l is False for l in ls):
            continue
        c = [-l for l in ls if l is not True]
        if not c:
            return 'UNSAT-fixed', 0
        cls.append(c)
    for x, y in itertools.combinations(V, 2):
        terms = []
        nfix = 0
        for k in V:
            if k in (x, y):
                continue
            a, b = lit(x, k), lit(y, k)
            if a is False or b is False:
                continue
            if a is True and b is True:
                nfix += 1
            elif a is True:
                terms.append(b)
            elif b is True:
                terms.append(a)
            else:
                cv = pool.id(('c', x, y, k))
                cls += [[-cv, a], [-cv, b], [cv, -a, -b]]
                terms.append(cv)
        if nfix > 2:
            return 'UNSAT-fixed', 0
        if len(terms) > 2 - nfix:
            cls += CardEnc.atmost(terms, bound=2 - nfix, vpool=pool, encoding=EncType.seqcounter).clauses
    t = {v: pool.id(('t', v)) for v in V}
    h = {v: pool.id(('h', v)) for v in V}
    T = pool.id(('true',))
    cls.append([T])
    for v in V:
        inc = [lit(v, u) for u in V if u != v]
        fc = sum(1 for l in inc if l is True)
        vs = [l for l in inc if l is not True and l is not False]

        def atleast(k, sel):
            need = k - fc
            if need <= 0:
                return
            if need > len(vs):
                cls.append([-sel]); return
            for c in CardEnc.atleast(vs, bound=need, vpool=pool, encoding=EncType.seqcounter).clauses:
                cls.append(c + [-sel])

        def atmost(k, sel):
            cap = k - fc
            if cap < 0:
                cls.append([-sel]); return
            if cap >= len(vs):
                return
            for c in CardEnc.atmost(vs, bound=cap, vpool=pool, encoding=EncType.seqcounter).clauses:
                cls.append(c + [-sel])

        atleast(5, T)
        atmost(d, T)                    # low vertices have degree <= 6 <= d, high ones <= d by the choice of z
        s6 = pool.id(('s6', v)); atleast(6, s6); cls.append([t[v], s6])
        s7 = pool.id(('s7', v)); atleast(7, s7)
        cls.append([-h[v], t[v], s7]); cls.append([-h[v], -t[v], s6])
        a5 = pool.id(('a5', v)); atmost(5, a5)
        a6 = pool.id(('a6', v)); atmost(6, a6)
        cls.append([h[v], -t[v], a5]); cls.append([h[v], t[v], a6])
        ad = pool.id(('ad', v)); atmost(d, ad); cls.append([-h[v], ad])        # high vertices have degree <= d
    cls.append([h[Z]])
    for a_, b_, c_, d_ in itertools.permutations(V, 4):
        if not (a_ < b_ and a_ < c_ and a_ < d_ and b_ < d_):
            continue
        ls = [lit(a_, b_), lit(b_, c_), lit(c_, d_), lit(d_, a_)]
        if any(l is False for l in ls):
            continue
        cls.append([-l for l in ls if l is not True] + [h[a_], h[b_], h[c_], h[d_]])
    lits = []
    for x, y in itertools.combinations(V, 2):
        l = lit(x, y)
        if l is False:
            continue
        le = pool.id(('le', x, y))
        if l is True:
            cls += [[-le, -h[x]], [-le, -h[y]], [le, h[x], h[y]]]
        else:
            cls += [[-le, l], [-le, -h[x]], [-le, -h[y]], [le, -l, h[x], h[y]]]
        lits += [le] * 4
    lits += [h[v] for v in V for _ in range(5)]
    # the Gallai count needs at least one low vertex: 4|E(low)| + 5 #high + 5 [some vertex low] <= 5n
    some_low = pool.id(('somelow',))
    cls.append([-some_low] + [-h[v] for v in V])
    for v in V:
        cls.append([some_low, h[v]])
    lits += [some_low] * 5
    cls += CardEnc.atmost(lits, bound=5 * n, vpool=pool, encoding=EncType.totalizer).clauses
    rounds = 0
    with Solver(name='cadical195', bootstrap_with=cls) as S:
        while True:
            if not S.solve():
                return 'UNSAT', rounds
            M = set(l for l in S.get_model() if l > 0)
            G = nx.Graph()
            G.add_nodes_from(V)
            for x, y in itertools.combinations(V, 2):
                l = lit(x, y)
                if l is True or (l is not False and l in M):
                    G.add_edge(x, y)
            col = three_colour(G)
            if col is None:
                return 'SAT-non3col', dict(edges=sorted(map(sorted, G.edges())), t=[int(t[v] in M) for v in V],
                                          high=[int(h[v] in M) for v in V])
            S.add_clause([var[frozenset((x, y))] for x, y in itertools.combinations(rest, 2) if col[x] == col[y]])
            rounds += 1
            if rounds > budget_rounds:
                return 'UNDECIDED', rounds


allres = {}
for n in [int(a) for a in sys.argv[1:]]:
    D = (n - 1) // 3
    t0 = time.time()
    res = []
    for d in range(6, D + 1):
        if d > 7:
            raise SystemExit('pair graphs from the atlas need d <= 7')
        pats = patterns(n, d, D)
        print(f'n = {n}, root degree {d} (max degree {D}): {len(pats)} attachment patterns', flush=True)
        for i, (m, P, S) in enumerate(pats):
            t1 = time.time()
            verdict, info = solve_pattern(n, d, D, m, P, S)
            res.append(dict(root_degree=d, m=m, pairs=P, singles=S, verdict=verdict,
                            info=info, seconds=round(time.time() - t1, 2)))
            if verdict != 'UNSAT':
                print(f'   pattern {i}: m={m}, pairs={P}, singles={S}: {verdict} {info}', flush=True)
    verdicts = sorted(set(r['verdict'] for r in res))
    rounds = sum(r['info'] for r in res if isinstance(r['info'], int))
    print(f'n = {n}: {len(res)} patterns, verdicts {verdicts}, non-3-colourability cuts used {rounds}, '
          f'{time.time()-t0:.1f} s', flush=True)
    allres[n] = res
json.dump(allres, open('05-knowledge/results/six_colour_exterior_twenty_sat_20261016.json', 'w'), indent=1)
print('wrote 05-knowledge/results/six_colour_exterior_twenty_sat_20261016.json')
