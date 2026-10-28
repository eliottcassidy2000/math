#!/usr/bin/env python3
"""Census of oriented (n-1)-edge carriers D on n vertices and their trijectivity.

For an oriented carrier D (underlying graph G with n vertices, n-1 edges), the tournaments containing D are the
orientations of the complement pair set Pi (|Pi| = C(n-1,2)).  A = Aut(D) acts on Pi with a sign (whether g swaps the
two ends of a pair).  D is LINEARLY TRIJECTIVE if the signed A-set Pi is isomorphic to the signed pair set of some
A-set X with |X| = n-1; then {T contains D} = Tour(X) A-equivariantly (an affine F2 isomorphism).
D is WEAKLY TRIJECTIVE if the A-set {T contains D} is isomorphic to Tour(X) for some X (Burnside marks agree on all
subgroups) without a pair-level identification.
Usage: python carrier_census.py NMAX
"""
import sys, itertools, json
from collections import Counter
import networkx as nx


def carriers(n):
    out = []
    for g in nx.graph_atlas_g():
        if g.number_of_nodes() == n and g.number_of_edges() == n - 1:
            out.append(sorted(tuple(sorted(e)) for e in g.edges()))
    return out


def aut_graph(n, edges):
    E = set(edges)
    deg = [0] * n
    for a, b in E:
        deg[a] += 1; deg[b] += 1
    res = []
    for p in itertools.permutations(range(n)):
        if any(deg[p[i]] != deg[i] for i in range(n)):
            continue
        ok = True
        for a, b in E:
            x, y = p[a], p[b]
            if (min(x, y), max(x, y)) not in E:
                ok = False; break
        if ok:
            res.append(p)
    return res


def arcs_of(edges, mask):
    return frozenset((a, b) if not (mask >> i) & 1 else (b, a) for i, (a, b) in enumerate(edges))


def mask_of(edges, arcs):
    m = 0
    for i, (a, b) in enumerate(edges):
        if (b, a) in arcs:
            m |= 1 << i
    return m


def orientation_reps(n, edges, autG):
    seen = set(); reps = []
    for mask in range(1 << len(edges)):
        if mask in seen:
            continue
        arcs = arcs_of(edges, mask)
        orb = set()
        for p in autG:
            orb.add(mask_of(edges, frozenset((p[a], p[b]) for a, b in arcs)))
        seen |= orb
        reps.append(min(orb))
    return reps


class Group:
    def __init__(self, elems):
        self.el = [tuple(e) for e in elems]
        self.idx = {e: i for i, e in enumerate(self.el)}
        m = len(self.el[0])
        self.mul = [[self.idx[tuple(a[b[k]] for k in range(m))] for b in self.el] for a in self.el]
        self.inv = [self.idx[tuple(sorted(range(m), key=lambda k: a[k]))] for a in self.el]
        self.e = self.idx[tuple(range(m))]

    def conj_set(self, c, S):
        ci = self.inv[c]
        return frozenset(self.mul[self.mul[c][s]][ci] for s in S)

    def closure(self, gens):
        S = {self.e} | set(gens); frontier = list(S)
        while frontier:
            new = []
            for a in frontier:
                for b in list(S):
                    for x in (self.mul[a][b], self.mul[b][a]):
                        if x not in S:
                            S.add(x); new.append(x)
            frontier = new
        return frozenset(S)

    def subgroups(self):
        subs = {frozenset([self.e])}
        cyc = {self.closure([g]) for g in range(len(self.el))}
        subs |= cyc
        changed = True
        while changed:
            changed = False
            cur = list(subs)
            for H in cur:
                for C in cyc:
                    if C <= H:
                        continue
                    J = self.closure(list(H | C))
                    if J not in subs:
                        subs.add(J); changed = True
        return subs

    def canon_sub(self, H):
        return min(tuple(sorted(self.conj_set(c, H))) for c in range(len(self.el)))


def signed_types(G, points, act, pairs):
    """Multiset of A-orbit types (stabiliser, reversing subset) up to conjugacy, for an A-set of unordered pairs."""
    seen = set(); types = Counter()
    for (a, b) in pairs:
        key = frozenset((a, b))
        if key in seen:
            continue
        orb = set(); S = []; R = []
        for i in range(len(G.el)):
            x, y = act(i, a), act(i, b)
            orb.add(frozenset((x, y)))
            if (x, y) == (a, b):
                S.append(i)
            elif (x, y) == (b, a):
                S.append(i); R.append(i)
        seen |= orb
        S = frozenset(S); R = frozenset(R)
        canon = min((tuple(sorted(G.conj_set(c, S))), tuple(sorted(G.conj_set(c, R)))) for c in range(len(G.el)))
        types[canon] += 1
    return types


def coset_action(G, H):
    """Points = left cosets gH; returns (points list, act(i, coset_index) -> coset_index)."""
    cos = []; where = {}
    for g in range(len(G.el)):
        c = frozenset(G.mul[g][h] for h in H)
        if c not in where:
            where[c] = len(cos); cos.append(c)
    table = [[where[frozenset(G.mul[i][x] for x in c)] for c in cos] for i in range(len(G.el))]
    return len(cos), table


def partitions_into(target, parts):
    """multisets (as sorted tuples of part ids) of the given parts (id, size) summing to target."""
    parts = sorted(parts, key=lambda t: t[1])
    res = []
    def rec(start, rem, cur):
        if rem == 0:
            res.append(tuple(cur)); return
        for j in range(start, len(parts)):
            pid, sz = parts[j]
            if sz <= rem:
                rec(j, rem - sz, cur + [pid])
    rec(0, target, [])
    return res


def marks(G, act, pairs, subgroups):
    """For each subgroup K: number of K-fixed orientations of the pair set (2^#orbits, or 0 if some K-elt reverses)."""
    out = {}
    for K in subgroups:
        rev = False; orbits = 0; seen = set()
        for (a, b) in pairs:
            if frozenset((a, b)) in seen:
                continue
            orbits += 1
            orb = set()
            for k in K:
                x, y = act(k, a), act(k, b)
                orb.add(frozenset((x, y)))
                if (x, y) == (b, a):
                    rev = True
            seen |= orb
        out[K] = 0 if rev else orbits
    return out


def classify(n, edges, mask, autG, max_order=200):
    arcs = arcs_of(edges, mask)
    A = [p for p in autG if frozenset((p[a], p[b]) for a, b in arcs) == arcs]
    G = Group(A)
    E = set(edges)
    Pi = [(a, b) for a, b in itertools.combinations(range(n), 2) if (a, b) not in E]
    outdeg = Counter(a for a, b in arcs); indeg = Counter(b for a, b in arcs)
    is_tree = nx.is_tree(nx.Graph(edges)) if len(set(v for e in edges for v in e)) == n else False
    Gx = nx.Graph(); Gx.add_nodes_from(range(n)); Gx.add_edges_from(edges)
    is_tree = nx.is_connected(Gx)
    functional = all(outdeg[v] <= 1 for v in range(n)) and sum(1 for v in range(n) if outdeg[v] == 0) == 1
    cofunctional = all(indeg[v] <= 1 for v in range(n)) and sum(1 for v in range(n) if indeg[v] == 0) == 1
    act_V = lambda i, x: G.el[i][x]
    piT = signed_types(G, range(n), act_V, Pi)
    # Burnside marked count
    cnt = 0
    for i in range(len(G.el)):
        K = G.closure([i])
        rev = False; orbits = 0; seen = set()
        for (a, b) in Pi:
            if frozenset((a, b)) in seen:
                continue
            orbits += 1; orb = set()
            for k in K:
                x, y = G.el[k][a], G.el[k][b]
                orb.add(frozenset((x, y)))
                if (x, y) == (b, a):
                    rev = True
            seen |= orb
        cnt += 0 if rev else 2 ** orbits
    marked = cnt // len(G.el)
    # fixed-subgraph criterion: for every subgroup K (cyclic first, then all), e(G[Fix K]) = |Fix K| - 1
    def fix_ok(K):
        F = [v for v in range(n) if all(G.el[k][v] == v for k in K)]
        if not F:
            return None
        eF = sum(1 for a, b in edges if a in F and b in F)
        return eF == len(F) - 1
    crit_cyc = all(fix_ok(G.closure([i])) is not False for i in range(len(G.el)))
    fixes_vertex = any(all(G.el[k][v] == v for k in range(len(G.el))) for v in range(n))
    res = {'order': len(G.el), 'tree': is_tree, 'functional': functional, 'cofunctional': cofunctional,
           'marked_types': marked, 'linear': None, 'weak': None, 'crit_cyclic': crit_cyc, 'A_fixes_vertex': fixes_vertex}
    if len(G.el) > max_order:
        return res
    subs = list(G.subgroups())
    res['crit_all'] = all(fix_ok(K) is not False for K in subs)
    reps = {}
    for H in subs:
        reps.setdefault(G.canon_sub(H), H)
    parts = [(cid, len(G.el) // len(H)) for cid, H in reps.items()]
    lin = False; weak = False; witness = None
    piMarks = marks(G, act_V, Pi, subs)
    for combo in partitions_into(n - 1, parts):
        # build X as disjoint union of coset spaces
        tables = []; offset = 0; size = 0
        for cid in combo:
            m, tab = coset_action(G, reps[cid]); tables.append((offset, m, tab)); offset += m
        def actX(i, x, tables=tables):
            for off, m, tab in tables:
                if off <= x < off + m:
                    return off + tab[i][x - off]
            raise ValueError
        pairsX = list(itertools.combinations(range(offset), 2))
        if not lin and signed_types(G, range(offset), actX, pairsX) == piT:
            lin = True; weak = True; witness = [len(G.el) // len(reps[c]) for c in combo]
            break
        if not weak:
            mX = marks(G, actX, pairsX, subs)
            if all((piMarks[K] == 0) == (mX[K] == 0) and (piMarks[K] == 0 or piMarks[K] == mX[K]) for K in subs):
                weak = True; witness = witness or ('weak', [len(G.el) // len(reps[c]) for c in combo])
    res.update({'linear': lin, 'weak': weak, 'witness': witness})
    return res


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    summary = {}
    for n in range(3, nmax + 1):
        rows = []
        for edges in carriers(n):
            autG = aut_graph(n, edges)
            for mask in orientation_reps(n, edges, autG):
                r = classify(n, edges, mask, autG)
                r['edges'] = edges; r['arcs'] = sorted(arcs_of(edges, mask))
                rows.append(r)
        c = Counter((r['tree'], r['functional'] or r['cofunctional'], r['linear'], r['weak']) for r in rows)
        agree = Counter((r['linear'], r.get('crit_all'), r['crit_cyclic'], r['A_fixes_vertex']) for r in rows)
        print(f"n={n}: (linear, crit_all_subgroups, crit_cyclic, A fixes a vertex) -> {dict(agree)}", flush=True)
        print(f"n={n}: {len(rows)} oriented carriers; (tree, func/cofunc, linear, weak) -> count:", flush=True)
        for k, v in sorted(c.items(), key=lambda kv: str(kv[0])):
            print("    ", k, v)
        bad = [r for r in rows if r['linear'] is False]
        for r in bad:
            print(f"   NON-LINEAR: |A|={r['order']} tree={r['tree']} weak={r['weak']} marked={r['marked_types']} arcs={r['arcs']}")
        summary[n] = rows
    json.dump({str(k): v for k, v in summary.items()}, open(f'carrier_census_n{nmax}.json', 'w'), default=str)


if __name__ == '__main__':
    main()
