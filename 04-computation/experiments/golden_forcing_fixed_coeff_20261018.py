#!/usr/bin/env python3
"""Sufficient additive test with FIXED coefficient colourings (fast), for large golden sources (S42).

The full additive pre-screen (golden_forcing_coupled_20261017.py) solves for c_A, c_1, c_2, c_3, c_4 jointly; for the
contact2 recipe (|A'| = 10,597) that CSP has about 1.2e7 clauses. Here the coefficient colourings are fixed first:
  1. solve the coefficient-only CSP (layers u, v, w, l over C', proper on the unit graph of C', value 0 at 0, the
     reflected-contact constraints), shared c_1 = c_2, enumerating several solutions and keeping the one minimising the
     total size of the contact shift sets;
  2. with these fixed, every mixed edge becomes a forbidden-difference constraint on c_A alone:
       v-contact  (a -> a + dq):  c_A(a+dq) - c_A(a) not in -S^V_q,  S^V_q = {D_q c_1(b) + D_q c_2(c)},
       e^2 contact (a -> a + dq): c_A(a+dq) - c_A(a) not in -S^E_q,  S^E_q = {D_q c_1(b) + D_(conj(d) q) c_4(l)},
     plus properness on the unit graph of A'; and the 570 port patterns are tested.
A pattern solved here extends additively (hence extends); a pattern not solved is inconclusive (other coefficient
colourings may work). Same subgraph caveat as the coupled script.
Usage: python golden_forcing_fixed_coeff_20261018.py AREC [NSOL]
Session opus-2026-10-17-S42.
"""
import os, sys, json, time, random
from pysat.solvers import Solver
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from golden_boundary_20261017c import PORT_EDGES, patterns
from golden_forcing_search_20261017 import Geo, recipe, log


def coverage(A, EA, qs, apairs, cand_list, argv):
    nA = len(A)
    X = [[6 * a + c + 1 for c in range(6)] for a in range(nA)]
    base = []
    for a in range(nA):
        base.append(X[a][:])
        for i in range(6):
            for k in range(i):
                base.append([-X[a][i], -X[a][k]])
    for a, b in EA:
        for k in range(6):
            base.append([-X[a][k], -X[b][k]])
    pats = patterns(8, PORT_EDGES)
    union = set(); bestcov = 0; covs = []
    for ci, (cost, cu, cll, SV, SE) in enumerate(cand_list):
        cl = list(base)
        for qi in range(len(qs)):
            forb = {(-s) % 6 for s in SV[qi]} | {(-s) % 6 for s in SE[qi]}
            for (a, a2) in apairs[qi]:
                for x in range(6):
                    for delta in forb:
                        cl.append([-X[a][x], -X[a2][(x + delta) % 6]])
        ok = set()
        with Solver(name='cadical195', bootstrap_with=cl) as sol:
            for pi, pat in enumerate(pats):
                if sol.solve(assumptions=[X[i][c] for i, c in enumerate(pat)]):
                    ok.add(pi)
        covs.append(len(ok)); union |= ok; bestcov = max(bestcov, len(ok))
    log('coverage over %d coefficient colourings: best single %d of %d patterns; union %d; mean %.1f'
        % (len(cand_list), bestcov, len(pats), len(union), sum(covs) / len(covs)))


def main():
    arec = sys.argv[1]; nsol = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    geo = Geo(); f = geo.f; d = geo.d; dbar = geo.bar(d)
    A = recipe(geo, arec)[0]; C = recipe(geo, 'base')[1]
    EA = geo.unit_edges(A); EC = geo.unit_edges(C)
    IA = {x: i for i, x in enumerate(A)}; IC = {x: i for i, x in enumerate(C)}
    units = {}
    for b, b2 in EC:
        for s, t in ((b, b2), (b2, b)):
            units.setdefault(f.sub(C[t], C[s]), []).append((s, t))
    qs = list(units)
    apairs = {}; lpairs = {}
    for qi, q in enumerate(qs):
        dq = f.mul(d, q); dbq = f.mul(dbar, q)
        apairs[qi] = [(a, IA[f.add(A[a], dq)]) for a in range(len(A)) if f.add(A[a], dq) in IA]
        lpairs[qi] = [(l, IC[f.add(C[l], dbq)]) for l in range(len(C)) if f.add(C[l], dbq) in IC]
    log('A %s |A| = %d (%d edges); |C| = %d (%d edges); directions %d; d-step pairs in A %d'
        % (arec, len(A), len(EA), len(C), len(EC), len(qs), sum(len(v) for v in apairs.values())))
    nC = len(C)
    # 1. coefficient-only CSP: layers u (= v, shared), w, l
    top = [0]
    def new():
        top[0] += 1; return top[0]
    Lu = [[new() for _ in range(6)] for _ in range(nC)]
    Lw = [[new() for _ in range(6)] for _ in range(nC)]
    Ll = [[new() for _ in range(6)] for _ in range(nC)]
    cl = []
    for V in (Lu, Lw, Ll):
        for row in V:
            cl.append(row[:])
            for i in range(6):
                for k in range(i):
                    cl.append([-row[i], -row[k]])
        for a, b in EC:
            for k in range(6):
                cl.append([-V[a][k], -V[b][k]])
        cl.append([V[0][0]])
    # reflected contact: Delta_q c_u(b) + Delta_q c_u(c) + Delta_{dbar q} c_w(l) != 0 for all q-edges b, c and dbar q-pairs l
    for qi, q in enumerate(qs):
        for (s, t) in units[q]:
            for (c, c2) in units[q]:
                for (l, l2) in lpairs[qi]:
                    for x1 in range(6):
                        for y1 in range(6):
                            for x2 in range(6):
                                for y2 in range(6):
                                    for x3 in range(6):
                                        y3 = (x3 - ((y1 - x1) + (y2 - x2))) % 6   # makes the total 0
                                        cl.append([-Lu[s][x1], -Lu[t][y1], -Lu[c][x2], -Lu[c2][y2], -Lw[l][x3], -Lw[l2][y3]])
    log('coefficient CSP: %d clauses' % len(cl))
    best = None
    cand_list = []
    random.seed(1)
    with Solver(name='cadical195', bootstrap_with=cl) as sol:
        for it in range(nsol):
            if not sol.solve():
                break
            model = set(v for v in sol.get_model() if v > 0)
            cu = [next(k for k in range(6) if Lu[b][k] in model) for b in range(nC)]
            cll = [next(k for k in range(6) if Ll[b][k] in model) for b in range(nC)]
            SV = {}; SE = {}
            for qi, q in enumerate(qs):
                Dq = {(cu[t] - cu[s]) % 6 for (s, t) in units[q]}
                SV[qi] = {(x + y) % 6 for x in Dq for y in Dq}
                Dl = {(cll[l2] - cll[l]) % 6 for (l, l2) in lpairs[qi]}
                SE[qi] = {(x + y) % 6 for x in Dq for y in Dl}
            cost = sum(len(v) for v in SV.values()) + sum(len(v) for v in SE.values())
            if best is None or cost < best[0]:
                best = (cost, cu, cll, SV, SE)
            cand_list.append((cost, cu, cll, SV, SE))
            # block this solution's u- and l-layer colourings to enumerate another
            sol.add_clause([-Lu[b][cu[b]] for b in range(nC)] + [-Ll[b][cll[b]] for b in range(nC)])
    if '--coverage' in sys.argv:
        coverage(A, EA, qs, apairs, cand_list, sys.argv)
        return
    cost, cu, cll, SV, SE = best
    log('best coefficient colouring after %d solutions: total shift-set size %d; c_u = %s, c_l = %s' % (nsol, cost, cu, cll))
    # 2. first-coordinate CSP with fixed shifts
    nA = len(A)
    X = [[6 * a + c + 1 for c in range(6)] for a in range(nA)]
    clA = []
    for a in range(nA):
        clA.append(X[a][:])
        for i in range(6):
            for k in range(i):
                clA.append([-X[a][i], -X[a][k]])
    for a, b in EA:
        for k in range(6):
            clA.append([-X[a][k], -X[b][k]])
    for qi in range(len(qs)):
        forb = {(-s) % 6 for s in SV[qi]} | {(-s) % 6 for s in SE[qi]}
        for (a, a2) in apairs[qi]:
            for x in range(6):
                for delta in forb:
                    clA.append([-X[a][x], -X[a2][(x + delta) % 6]])
    log('first-coordinate CSP: %d clauses' % len(clA))
    pats = patterns(8, PORT_EDGES)
    unsolved = []
    t0 = time.time()
    with Solver(name='cadical195', bootstrap_with=clA) as sol:
        for pat in pats:
            if not sol.solve(assumptions=[X[i][c] for i, c in enumerate(pat)]):
                unsolved.append(pat)
    log('fixed-coefficient additive test: %d of %d patterns extend; %d inconclusive [%.0f s]'
        % (len(pats) - len(unsolved), len(pats), len(unsolved), time.time() - t0))
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(dict(A_recipe=arec, A=nA, edges_A=len(EA), shift_cost=cost, c_u=cu, c_l=cll,
                   extend=len(pats) - len(unsolved), inconclusive=len(unsolved), inconclusive_list=unsolved),
              open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'golden_forcing_fixed_coeff_%s_20261018.json' % arec), 'w'), indent=1)


if __name__ == '__main__':
    main()
