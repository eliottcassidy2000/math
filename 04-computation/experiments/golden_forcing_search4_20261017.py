#!/usr/bin/env python3
"""Two coupled contact families: additive pre-screen and SAT for A + uC + vC + Y C, Y = e^2 = -Z/conj(Z) (S41).

golden_product_direction_census_20261017.py found that Y = e^2 (an element of the base field F) creates cross edges with
the golden source. golden_product_relations_20261017.py identifies them as one contact family q (d, 1, 0, conj d) in the
coordinates (1, u, v, Y): indeed d + u + conj(d) Y = Z conj(u)/conj(Z) is a unit. Together with clock_duality's v-contact
q (d, 1, 1, 0), both families share the same (A, u)-step dq, q. For the coefficient sets of the golden source the census
is complete (every unit difference with nonzero Y-part is a pure Y-axis unit or in this family), and G1
(marked_contacts_root_current_20261017e.md section 5) classifies the
rest. This script builds the 4-coordinate source combinatorially, runs the additive pre-screen
c_A(a) + c_1(b) + c_2(c) + c_3(l) mod 6, and runs the full SAT problem (relative six-core first) on the patterns the
pre-screen blocks, or on all patterns with --full.
Usage: python golden_forcing_search4_20261017.py [--full] [--ylayer C|pent|M]
Session opus-2026-10-17-S41.
"""
import os, sys, json, time
from fractions import Fraction as Q
from pysat.solvers import Solver
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from escaped_ports_20261017d import need
from golden_boundary_20261017c import PORT_EDGES, patterns
from golden_forcing_search_20261017 import Geo, log

T0 = time.time()


def main():
    full = '--full' in sys.argv
    ylayer = sys.argv[sys.argv.index('--ylayer') + 1] if '--ylayer' in sys.argv else 'C'
    geo = Geo(); f = geo.f; P = list(geo.ports); M = geo.M; d = geo.d
    dbar = geo.bar(d)
    def uniq(lst):
        out = []; seen = set()
        for x in lst:
            if x not in seen:
                seen.add(x); out.append(x)
        return out
    A = uniq(P + [f.add(p, m) for p in P for m in M])
    C = uniq(M + P[:5])
    Ly = {'C': C, 'pent': uniq([f.zero] + P[:5]), 'M': uniq(M)}[ylayer]
    EA = geo.unit_edges(A); EC = geo.unit_edges(C); EL = geo.unit_edges(Ly)
    IA = {x: i for i, x in enumerate(A)}; IC = {x: i for i, x in enumerate(C)}; IL = {x: i for i, x in enumerate(Ly)}
    units = {}
    for b, b2 in EC:
        for s, t in ((b, b2), (b2, b)):
            units.setdefault(f.sub(C[t], C[s]), []).append((s, t))
    mixV, mixY = set(), set()
    for q, cpairs in units.items():
        dq = f.mul(d, q); dbq = f.mul(dbar, q)
        apairs = [(a, IA[f.add(A[a], dq)]) for a in range(len(A)) if f.add(A[a], dq) in IA]
        lpairs = [(l, IL[f.add(Ly[l], dbq)]) for l in range(len(Ly)) if f.add(Ly[l], dbq) in IL]
        for (a, a2) in apairs:
            for (s, t) in cpairs:
                for (c, c2) in cpairs:
                    mixV.add(((a, s, c), (a2, t, c2)))          # v-contact, any Y-layer index
                for (l, l2) in lpairs:
                    mixY.add(((a, s, l), (a2, t, l2)))          # Y-contact, any v-layer index
    log('|A| %d |C| %d |L_Y| %d; edges A %d C %d L %d; v-contact patterns %d, Y-contact patterns %d'
        % (len(A), len(C), len(Ly), len(EA), len(EC), len(EL), len(mixV), len(mixY)))
    pats = patterns(8, PORT_EDGES)
    nA, nC, nL = len(A), len(C), len(Ly)

    def additive(shared):
        top = [0]
        def new():
            top[0] += 1; return top[0]
        X = [[new() for _ in range(6)] for _ in range(nA)]
        Y1 = [[new() for _ in range(6)] for _ in range(nC)]
        Y2 = Y1 if shared else [[new() for _ in range(6)] for _ in range(nC)]
        Y3 = [[new() for _ in range(6)] for _ in range(nL)]
        cl = []
        for V in ([X, Y1, Y3] if shared else [X, Y1, Y2, Y3]):
            for row in V:
                cl.append(row[:])
                for i in range(6):
                    for k in range(i):
                        cl.append([-row[i], -row[k]])
        for (V, E) in ((X, EA), (Y1, EC), (Y2, EC), (Y3, EL)):
            for a, b in E:
                for k in range(6):
                    cl.append([-V[a][k], -V[b][k]])
        cl += [[Y1[0][0]], [Y2[0][0]], [Y3[0][0]]]
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
            tri(diff(X, 'A', a, a2), diff(Y1, '1', s, t), diff(Y2, '2' if not shared else '1', c, c2))
        for (a, s, l), (a2, t, l2) in mixY:
            tri(diff(X, 'A', a, a2), diff(Y1, '1', s, t), diff(Y3, '3', l, l2))
        blocked = []
        with Solver(name='g3', bootstrap_with=cl) as sol:
            for pat in pats:
                if not sol.solve(assumptions=[X[i][c] for i, c in enumerate(pat)]):
                    blocked.append(pat)
        return blocked
    bs = additive(True); bf = additive(False)
    log('additive pre-screen: blocked with c_1 = c_2 (shared): %d; with independent coefficient colourings: %d (of %d)'
        % (len(bs), len(bf), len(pats)))
    out = dict(ylayer=ylayer, A=nA, C=nC, L=nL, mixV=len(mixV), mixY=len(mixY), blocked_shared=len(bs), blocked=len(bf), blocked_list=bf)
    if bf or full:
        test = bf if bf else pats
        idx = lambda a, b, c, l: ((a * nC + b) * nC + c) * nL + l
        n = nA * nC * nC * nL
        edges = set()
        for a, a2 in EA:
            for b in range(nC):
                for c in range(nC):
                    for l in range(nL):
                        edges.add((idx(a, b, c, l), idx(a2, b, c, l)))
        for b, b2 in EC:
            for a in range(nA):
                for c in range(nC):
                    for l in range(nL):
                        edges.add((idx(a, b, c, l), idx(a, b2, c, l)))
                        edges.add((idx(a, c, b, l), idx(a, c, b2, l)))
        for l, l2 in EL:
            for a in range(nA):
                for b in range(nC):
                    for c in range(nC):
                        edges.add((idx(a, b, c, l), idx(a, b, c, l2)))
        for (a, s, c), (a2, t, c2) in mixV:
            for l in range(nL):
                edges.add((idx(a, s, c, l), idx(a2, t, c2, l)))
        for (a, s, l), (a2, t, l2) in mixY:
            for c in range(nC):
                edges.add((idx(a, s, c, l), idx(a2, t, c, l2)))
        nb = [set() for _ in range(n)]
        for x, y in edges:
            nb[x].add(y); nb[y].add(x)
        ports = [idx(i, 0, 0, 0) for i in range(8)]
        live = set(range(n))
        while True:
            batch = [v for v in live if v not in ports and len(nb[v] & live) <= 5]
            if not batch:
                break
            for v in batch:
                live.remove(v)
        log('full graph: %d vertices, %d edges, relative six-core %d' % (n, len(edges), len(live)))
        var = {v: i for i, v in enumerate(sorted(live))}
        blocked = []
        with Solver(name='g3') as sol:
            for v, i in var.items():
                sol.add_clause([6 * i + c + 1 for c in range(6)])
            for x, y in edges:
                if x in var and y in var:
                    for c in range(6):
                        sol.add_clause([-(6 * var[x] + c + 1), -(6 * var[y] + c + 1)])
            for pat in test:
                if not sol.solve(assumptions=[6 * var[ports[i]] + c + 1 for i, c in enumerate(pat)]):
                    blocked.append(pat)
        log('full SAT: blocked %d of %d tested' % (len(blocked), len(test)))
        out.update(full_vertices=n, full_edges=len(edges), six_core=len(live), full_blocked=blocked)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'golden_forcing_search4_%s_20261017.json' % ylayer), 'w'), indent=1)


if __name__ == '__main__':
    main()
