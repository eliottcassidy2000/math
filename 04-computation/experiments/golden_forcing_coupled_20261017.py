#!/usr/bin/env python3
"""All three known golden contact families in one source: additive pre-screen with CaDiCaL (S42).

S41 (golden_additive_boundary_20261017.md sections 3c-4) found three families of mixed unit edges for the golden source
in the coordinates (1, u, v, uv, Y), Y = e^2 = -Z/conj(Z):
    v-contact        q (d, 1, 1, 0, 0)        (clock_duality's G1 contact)
    reflected contact q (0, 1, 1, conj d, 0)  (uv family: |u + v + conj(d) uv| = 1)
    e^2 contact      q (d, 1, 0, 0, conj d)   (|d + u + conj(d) e^2| = 1)
Each was tested separately (all additively free). This script puts all three into one five-coordinate source
A' + uC' + vC' + (uv)W + e^2 L (W = L = C'), i.e. the v- and e^2-contacts share their (A, u)-step and the reflected and
v-contacts share their (u, v)-step, and runs the additive pre-screen

    c(a, b, c, w, l) = c_A(a) + c_1(b) + c_2(c) + c_3(w) + c_4(l)  (mod 6)

for all 570 canonical port patterns (shared c_1 = c_2, and independent). No census of differences with both uv- and
e^2-parts has been made, so the graph built here is a SUBGRAPH of the unit-distance graph on these addresses: a blocked
pattern would still have to be certified on the full graph, while freedom here does not prove freedom of the full graph.
Usage: python golden_forcing_coupled_20261017.py AREC CREC     (AREC in base, dM, contact1; CREC in base, MM)
Session opus-2026-10-17-S42.
"""
import os, sys, json, time
from pysat.solvers import Solver
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from golden_boundary_20261017c import PORT_EDGES, patterns
from golden_forcing_search_20261017 import Geo, recipe, log


def main():
    arec, crec = sys.argv[1], sys.argv[2]
    geo = Geo(); f = geo.f; d = geo.d; dbar = geo.bar(d)
    A = recipe(geo, arec)[0]
    C = recipe(geo, 'MM')[1] if crec == 'MM' else recipe(geo, 'base')[1]
    EA = geo.unit_edges(A); EC = geo.unit_edges(C)
    IA = {x: i for i, x in enumerate(A)}; IC = {x: i for i, x in enumerate(C)}
    units = {}
    for b, b2 in EC:
        for s, t in ((b, b2), (b2, b)):
            units.setdefault(f.sub(C[t], C[s]), []).append((s, t))
    V, R, E2 = set(), set(), set()
    for q, cp in units.items():
        dq = f.mul(d, q); dbq = f.mul(dbar, q)
        ap = [(a, IA[f.add(A[a], dq)]) for a in range(len(A)) if f.add(A[a], dq) in IA]
        lp = [(l, IC[f.add(C[l], dbq)]) for l in range(len(C)) if f.add(C[l], dbq) in IC]
        for (s, t) in cp:
            for (a, a2) in ap:
                for (c, c2) in cp:
                    V.add(((a, s, c), (a2, t, c2)))
                for (l, l2) in lp:
                    E2.add(((a, s, l), (a2, t, l2)))
            for (c, c2) in cp:
                for (l, l2) in lp:
                    R.add(((s, c, l), (t, c2, l2)))
    log('A %s |A| = %d (%d edges); C %s |C| = %d (%d edges); unit directions %d; patterns: v-contact %d, reflected %d, e^2 %d'
        % (arec, len(A), len(EA), crec, len(C), len(EC), len(units), len(V), len(R), len(E2)))
    pats = patterns(8, PORT_EDGES)
    nA, nC = len(A), len(C)
    res = {}
    for shared in (True, False):
        top = [0]
        def new():
            top[0] += 1; return top[0]
        X = [[new() for _ in range(6)] for _ in range(nA)]
        L1 = [[new() for _ in range(6)] for _ in range(nC)]
        L2 = L1 if shared else [[new() for _ in range(6)] for _ in range(nC)]
        L3 = [[new() for _ in range(6)] for _ in range(nC)]
        L4 = [[new() for _ in range(6)] for _ in range(nC)]
        layers = [X, L1, L3, L4] if shared else [X, L1, L2, L3, L4]
        cl = []
        for Vv in layers:
            for row in Vv:
                cl.append(row[:])
                for i in range(6):
                    for k in range(i):
                        cl.append([-row[i], -row[k]])
        for a, b in EA:
            for k in range(6):
                cl.append([-X[a][k], -X[b][k]])
        for Vv in layers[1:]:
            for a, b in EC:
                for k in range(6):
                    cl.append([-Vv[a][k], -Vv[b][k]])
            cl.append([Vv[0][0]])
        cache = {}
        def diff(Vv, tag, i, i2):
            key = (tag, i, i2)
            if key not in cache:
                Dv = [new() for _ in range(6)]
                for x in range(6):
                    for y in range(6):
                        cl.append([-Vv[i2][x], -Vv[i][y], Dv[(x - y) % 6]])
                cache[key] = Dv
            return cache[key]
        def tri(D1, D2, D3):
            for x in range(6):
                for y in range(6):
                    cl.append([-D1[x], -D2[y], -D3[(-x - y) % 6]])
        t2 = '1' if shared else '2'
        for (a, s, c), (a2, t, c2) in V:
            tri(diff(X, 'A', a, a2), diff(L1, '1', s, t), diff(L2, t2, c, c2))
        for (s, c, l), (t, c2, l2) in R:
            tri(diff(L1, '1', s, t), diff(L2, t2, c, c2), diff(L3, '3', l, l2))
        for (a, s, l), (a2, t, l2) in E2:
            tri(diff(X, 'A', a, a2), diff(L1, '1', s, t), diff(L4, '4', l, l2))
        t0 = time.time()
        blocked = []
        with Solver(name='cadical195', bootstrap_with=cl) as sol:
            base_ok = sol.solve()
            for pat in pats:
                if not sol.solve(assumptions=[X[i][c] for i, c in enumerate(pat)]):
                    blocked.append(pat)
        key = 'shared' if shared else 'independent'
        res[key] = dict(clauses=len(cl), variables=top[0], satisfiable_without_ports=bool(base_ok), blocked=len(blocked),
                        blocked_list=blocked)
        log('additive pre-screen (%s coefficient colourings): %d clauses; blocked %d of %d patterns [%.0f s]'
            % (key, len(cl), len(blocked), len(pats), time.time() - t0))
    out = dict(A_recipe=arec, C_recipe=crec, A=nA, C=nC, edges_A=len(EA), edges_C=len(EC), directions=len(units),
               v_contact=len(V), reflected=len(R), e2_contact=len(E2), results=res)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results',
                                     'golden_forcing_coupled_%s_%s_20261017.json' % (arec, crec)), 'w'), indent=1)


if __name__ == '__main__':
    main()
