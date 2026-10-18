#!/usr/bin/env python3
"""Additive pre-screen for A + uC + vC + (uv)C with the reflected contact family (S41).

golden_product_direction_census_20261017.py found that the coordinate uv adds exactly one family of cross edges for the
golden source's coefficient sets, identified by golden_product_relations_20261017.py as q (0, 1, 1, conj d) in the
coordinates (1, u, v, uv) (the reflected contact: |u + v + conj(d) uv| = |conj(d + u + v)| = 1). Together with the
v-contact q (d, 1, 1, 0), lemma G1 and the census give the complete edge set for these coefficient sets. This script runs
the additive pre-screen c_A(a) + c_1(b) + c_2(c) + c_3(l) mod 6 for all 570 port patterns.
Usage: python golden_forcing_search_uv_20261017.py
Session opus-2026-10-17-S41.
"""
import os, sys, json
from pysat.solvers import Solver
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from golden_boundary_20261017c import PORT_EDGES, patterns
from golden_forcing_search_20261017 import Geo, log


def main():
    geo = Geo(); f = geo.f; P = list(geo.ports); M = geo.M; d = geo.d
    dbar = geo.bar(d)
    def uniq(lst):
        out = []; seen = set()
        for x in lst:
            if x not in seen:
                seen.add(x); out.append(x)
        return out
    A = uniq(P + [f.add(p, m) for p in P for m in M]); C = uniq(M + P[:5]); L = C
    EA = geo.unit_edges(A); EC = geo.unit_edges(C)
    IA = {x: i for i, x in enumerate(A)}; IC = {x: i for i, x in enumerate(C)}
    units = {}
    for b, b2 in EC:
        for s, t in ((b, b2), (b2, b)):
            units.setdefault(f.sub(C[t], C[s]), []).append((s, t))
    mixV, mixR = set(), set()
    for q, cp in units.items():
        dq = f.mul(d, q); dbq = f.mul(dbar, q)
        ap = [(a, IA[f.add(A[a], dq)]) for a in range(len(A)) if f.add(A[a], dq) in IA]
        lp = [(l, IC[f.add(L[l], dbq)]) for l in range(len(L)) if f.add(L[l], dbq) in IC]
        for (s, t) in cp:
            for (c, c2) in cp:
                for (a, a2) in ap:
                    mixV.add(((a, s, c), (a2, t, c2)))        # v-contact, uv-layer fixed
                for (l, l2) in lp:
                    mixR.add(((s, c, l), (t, c2, l2)))        # reflected contact, A fixed
    log('|A| %d |C| %d; v-contact patterns %d, reflected-contact patterns %d' % (len(A), len(C), len(mixV), len(mixR)))
    pats = patterns(8, PORT_EDGES)
    nA, nC = len(A), len(C)
    top = [0]
    def new():
        top[0] += 1; return top[0]
    X = [[new() for _ in range(6)] for _ in range(nA)]
    Y = [[[new() for _ in range(6)] for _ in range(nC)] for _ in range(3)]
    cl = []
    for V in [X] + Y:
        for row in V:
            cl.append(row[:])
            for i in range(6):
                for k in range(i):
                    cl.append([-row[i], -row[k]])
    for a, b in EA:
        for k in range(6):
            cl.append([-X[a][k], -X[b][k]])
    for V in Y:
        for a, b in EC:
            for k in range(6):
                cl.append([-V[a][k], -V[b][k]])
        cl.append([V[0][0]])
    cache = {}
    def diff(V, tag, i, i2):
        key = (tag, i, i2)
        if key not in cache:
            Dv = [new() for _ in range(6)]
            for x in range(6):
                for y in range(6):
                    cl.append([-V[i2][x], -V[i][y], Dv[(x - y) % 6]])
            cache[key] = Dv
        return cache[key]
    def tri(D1, D2, D3):
        for x in range(6):
            for y in range(6):
                cl.append([-D1[x], -D2[y], -D3[(-x - y) % 6]])
    for (a, s, c), (a2, t, c2) in mixV:
        tri(diff(X, 'A', a, a2), diff(Y[0], 'u', s, t), diff(Y[1], 'v', c, c2))
    for (s, c, l), (t, c2, l2) in mixR:
        tri(diff(Y[0], 'u', s, t), diff(Y[1], 'v', c, c2), diff(Y[2], 'w', l, l2))
    blocked = []
    with Solver(name='g3', bootstrap_with=cl) as sol:
        for pat in pats:
            if not sol.solve(assumptions=[X[i][c] for i, c in enumerate(pat)]):
                blocked.append(pat)
    log('additive pre-screen with the reflected contact: blocked %d of %d' % (len(blocked), len(pats)))
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(dict(mixV=len(mixV), mixR=len(mixR), blocked=len(blocked)), open(os.path.join(
        here, '..', '..', '05-knowledge', 'results', 'golden_forcing_search_uv_20261017.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
