#!/usr/bin/env python3
"""Can the reflected (uv) contact force odd reflection? Additive pre-screen with an enlarged uv-layer (S42).

Mechanism (preprint_connections_solfeggio_digits_20261017.md section 7, L6). In the source A' + uC' + vC' + (uv)W the
reflected contact joins (b, c, w) to (b + q, c + q, w + conj(d) q), so additive colourings must satisfy
    Delta_q c_1(b) + Delta_q c_2(c) + Delta_(conj(d) q) c_3(w) != 0  (mod 6).
If every proper colouring c_3 of W showed all five nonzero differences on the conj(d) q-steps of W, this would force
Delta_q c_1(b) + Delta_q c_2(c) = 0 for all q-edges, i.e. homomorphism-like coefficient colourings with c_2 = -c_1 (odd
reflection, zero contact shifts), and the v-contact would reduce to plain 6-colouring of the golden two-distance graph
on A', which blocks 465 of the 570 patterns. This script enlarges W to W = C' + conj(d) ({0} + U(C')), U = unit
differences, computes its unit graph exactly, and runs the additive pre-screen with the v-contact and reflected contact
(coordinates a, b, c, w), recording also the number of distinct conj(d) q-differences in a solver colouring of W.
Usage: python golden_reflection_forcing_20261017.py [WREC]   (WREC in C, Cd, Cdd)
Session opus-2026-10-17-S42.
"""
import os, sys, json, time
from pysat.solvers import Solver
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from golden_boundary_20261017c import PORT_EDGES, patterns
from golden_forcing_search_20261017 import Geo, recipe, log


def main():
    wrec = sys.argv[1] if len(sys.argv) > 1 else 'Cd'
    geo = Geo(); f = geo.f; d = geo.d; dbar = geo.bar(d)
    A, C = recipe(geo, 'base')
    def uniq(lst):
        out = []; seen = set()
        for x in lst:
            if x not in seen:
                seen.add(x); out.append(x)
        return out
    U = sorted({f.sub(b, a) for a in C for b in C if a != b and geo.is_unit(f.sub(b, a))})
    if wrec == 'C':
        W = list(C)
    elif wrec == 'Cd':
        W = uniq(list(C) + [f.add(c, f.mul(dbar, q)) for c in C for q in U])
    elif wrec == 'Cdd':
        W1 = uniq(list(C) + [f.add(c, f.mul(dbar, q)) for c in C for q in U])
        W = uniq(W1 + [f.add(w, f.mul(dbar, q)) for w in W1 for q in U])
    else:
        raise ValueError(wrec)
    EA = geo.unit_edges(A); EC = geo.unit_edges(C); EW = geo.unit_edges(W)
    IA = {x: i for i, x in enumerate(A)}; IW = {x: i for i, x in enumerate(W)}
    units = {}
    for b, b2 in EC:
        for s, t in ((b, b2), (b2, b)):
            units.setdefault(f.sub(C[t], C[s]), []).append((s, t))
    V, R = set(), set()
    wpairs_by_q = {}
    for q, cp in units.items():
        dq = f.mul(d, q); dbq = f.mul(dbar, q)
        ap = [(a, IA[f.add(A[a], dq)]) for a in range(len(A)) if f.add(A[a], dq) in IA]
        wp = [(w, IW[f.add(W[w], dbq)]) for w in range(len(W)) if f.add(W[w], dbq) in IW]
        wpairs_by_q[q] = wp
        for (s, t) in cp:
            for (c, c2) in cp:
                for (a, a2) in ap:
                    V.add(((a, s, c), (a2, t, c2)))
                for (w, w2) in wp:
                    R.add(((s, c, w), (t, c2, w2)))
    log('W %s: |W| = %d with %d unit edges; |A| %d (%d edges), |C| %d (%d edges); v-contact %d, reflected %d; '
        'conj(d)q-pairs in W per direction: min %d max %d'
        % (wrec, len(W), len(EW), len(A), len(EA), len(C), len(EC), len(V), len(R),
           min(len(v) for v in wpairs_by_q.values()), max(len(v) for v in wpairs_by_q.values())))
    pats = patterns(8, PORT_EDGES)
    nA, nC, nW = len(A), len(C), len(W)
    res = {}
    for shared in (True, False):
        top = [0]
        def new():
            top[0] += 1; return top[0]
        X = [[new() for _ in range(6)] for _ in range(nA)]
        L1 = [[new() for _ in range(6)] for _ in range(nC)]
        L2 = L1 if shared else [[new() for _ in range(6)] for _ in range(nC)]
        L3 = [[new() for _ in range(6)] for _ in range(nW)]
        cl = []
        for Vv in ([X, L1, L3] if shared else [X, L1, L2, L3]):
            for row in Vv:
                cl.append(row[:])
                for i in range(6):
                    for k in range(i):
                        cl.append([-row[i], -row[k]])
        for (Vv, E) in ((X, EA), (L1, EC), (L2, EC), (L3, EW)):
            for a, b in E:
                for k in range(6):
                    cl.append([-Vv[a][k], -Vv[b][k]])
        cl += [[L1[0][0]], [L2[0][0]], [L3[0][0]]]
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
        for (s, c, w), (t, c2, w2) in R:
            tri(diff(L1, '1', s, t), diff(L2, t2, c, c2), diff(L3, '3', w, w2))
        blocked = []
        ndiff = None
        with Solver(name='cadical195', bootstrap_with=cl) as sol:
            if sol.solve():
                model = set(v for v in sol.get_model() if v > 0)
                col3 = [next(k for k in range(6) if L3[w][k] in model) for w in range(nW)]
                ndiff = {}
                for q, wp in wpairs_by_q.items():
                    ndiff[str(len(ndiff))] = sorted({(col3[w2] - col3[w]) % 6 for (w, w2) in wp})
            for pat in pats:
                if not sol.solve(assumptions=[X[i][c] for i, c in enumerate(pat)]):
                    blocked.append(pat)
        key = 'shared' if shared else 'independent'
        sizes = sorted(len(v) for v in ndiff.values()) if ndiff else []
        res[key] = dict(blocked=len(blocked), blocked_list=blocked, w_difference_sets=ndiff)
        log('additive pre-screen (%s): blocked %d of %d; sizes of conj(d)q-difference sets in one W-colouring: %s'
            % (key, len(blocked), len(pats), sizes))
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(dict(W_recipe=wrec, W=nW, edges_W=len(EW), v_contact=len(V), reflected=len(R), results=res),
              open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'golden_reflection_forcing_%s_20261017.json' % wrec), 'w'),
              indent=1)


if __name__ == '__main__':
    main()
