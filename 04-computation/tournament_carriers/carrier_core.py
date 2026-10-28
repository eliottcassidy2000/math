#!/usr/bin/env python3
"""Shared machinery for oriented carriers (kind-pasteur-2026-10-27-S8).

An oriented carrier D on vertex set range(n): a set of arcs (a, b) with no 2-cycles.  Its complement pair set Pi is
the set of unordered pairs that carry no arc.  A = Aut(D) acts on Pi with a sign (whether an element swaps the two ends
of a pair).  For j >= 1 and an A-set X of size n - j, the carrier is LINEARLY j-JECTIVE if Pi is isomorphic to the
signed pair set of X (j = 1: "trijective" in the sense of the S3 trijection); WEAKLY j-JECTIVE if only the A-sets of
orientations (tournaments containing D, resp. Tour(X)) are isomorphic (equal Burnside marks).
"""
import itertools
from collections import Counter
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher


def automorphisms(n, arcs):
    D = nx.DiGraph(); D.add_nodes_from(range(n)); D.add_edges_from(arcs)
    out = []
    for m in DiGraphMatcher(D, D).isomorphisms_iter():
        out.append(tuple(m[i] for i in range(n)))
    return out


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
            for H in list(subs):
                for C in cyc:
                    if not C <= H:
                        J = self.closure(list(H | C))
                        if J not in subs:
                            subs.add(J); changed = True
        return subs

    def canon_sub(self, H):
        return min(tuple(sorted(self.conj_set(c, H))) for c in range(len(self.el)))


def signed_types(G, act, pairs):
    seen = set(); types = Counter()
    for (a, b) in pairs:
        if frozenset((a, b)) in seen:
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
        types[min((tuple(sorted(G.conj_set(c, S))), tuple(sorted(G.conj_set(c, R)))) for c in range(len(G.el)))] += 1
    return types


def coset_action(G, H):
    cos = []; where = {}
    for g in range(len(G.el)):
        c = frozenset(G.mul[g][h] for h in H)
        if c not in where:
            where[c] = len(cos); cos.append(c)
    return len(cos), [[where[frozenset(G.mul[i][x] for x in c)] for c in cos] for i in range(len(G.el))]


def partitions_into(target, parts):
    parts = sorted(parts, key=lambda t: t[1]); res = []
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
    out = {}
    for K in subgroups:
        rev = False; orbits = 0; seen = set()
        for (a, b) in pairs:
            if frozenset((a, b)) in seen:
                continue
            orbits += 1; orb = set()
            for k in K:
                x, y = act(k, a), act(k, b)
                orb.add(frozenset((x, y)))
                if (x, y) == (b, a):
                    rev = True
            seen |= orb
        out[K] = None if rev else orbits
    return out


def analyse(n, arcs, j=1, max_order=400, want_weak=True):
    """Return dict with A order, Maxwell criterion, linear / weak j-jectivity, Burnside marked count."""
    arcs = [tuple(a) for a in arcs]
    edges = {frozenset(a) for a in arcs}
    A = automorphisms(n, arcs)
    G = Group(A)
    Pi = [(a, b) for a, b in itertools.combinations(range(n), 2) if frozenset((a, b)) not in edges]
    actV = lambda i, x: G.el[i][x]
    # Burnside count of marked types
    tot = 0
    for i in range(len(G.el)):
        K = G.closure([i])
        m = marks(G, actV, Pi, [K])[K]
        tot += 0 if m is None else 2 ** m
    res = {'order': len(G.el), 'marked_types': tot // len(G.el), 'fixes_vertex': any(
        all(G.el[k][v] == v for k in range(len(G.el))) for v in range(n))}
    if len(G.el) > max_order:
        res.update({'linear': None, 'weak': None, 'maxwell': None})
        return res
    subs = list(G.subgroups())
    def maxwell_ok(K):
        F = {v for v in range(n) if all(G.el[k][v] == v for k in K)}
        eF = sum(1 for e in edges if e <= F)
        return eF == j * len(F) - j * (j + 1) // 2
    res['maxwell'] = all(maxwell_ok(K) for K in subs)
    reps = {}
    for H in subs:
        reps.setdefault(G.canon_sub(H), H)
    parts = [(cid, len(G.el) // len(H)) for cid, H in reps.items()]
    piT = signed_types(G, actV, Pi)
    piM = marks(G, actV, Pi, subs) if want_weak else None
    lin = False; weak = False
    for combo in partitions_into(n - j, parts):
        tables = []; off = 0
        for cid in combo:
            m, tab = coset_action(G, reps[cid]); tables.append((off, m, tab)); off += m
        def actX(i, x, tables=tables):
            for o, m, tab in tables:
                if o <= x < o + m:
                    return o + tab[i][x - o]
        pairsX = list(itertools.combinations(range(off), 2))
        if not lin and signed_types(G, actX, pairsX) == piT:
            lin = True; weak = True
            break
        if want_weak and not weak:
            mX = marks(G, actX, pairsX, subs)
            if all(piM[K] == mX[K] for K in subs):
                weak = True
    res.update({'linear': lin, 'weak': weak if want_weak else None})
    return res
