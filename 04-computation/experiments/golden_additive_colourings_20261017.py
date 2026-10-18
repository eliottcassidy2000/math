#!/usr/bin/env python3
"""Why every port pattern extends: additive colourings of clock_duality's shared golden contact source (S41).

The source G = A + uC + vC of shared_golden_contact_20261017g.md (5400 vertices, 30856 edges) consists of three axis
families (Cartesian-sum edges) and 56 mixed edges in the contact direction q(d, 1, 1). For a Cartesian sum, every
colouring of the first factor extends: colour (a, b, c) by c_A(a) + c_1(b) + c_2(c) mod 6. This script asks, for each of
the 570 canonical six-colour port patterns, whether such an ADDITIVE colouring exists that also respects the mixed edges,
i.e. c_A(a') - c_A(a) + c_1(b') - c_1(b) + c_2(c') - c_2(c) != 0 mod 6 on every mixed edge (a,b,c) ~ (a',b',c').

If additive colourings exist for every pattern, the boundary relation is unrestricted for a structural reason (the
mixed edges are cohomologically trivial for Z/6), and a forcing attachment must create a cocycle obstruction. The script
also reports the mixed-edge difference types, verifies every additive colouring against the full edge list, and counts
how many patterns survive when the group Z/6 is replaced by the 'rigid' ansatz c_1 = c_2 (shared golden coefficients).
Input: the receipt 04-computation/experiments/root_positive_transport_20261017g/golden_shared.json (clock_duality_oct15).
Usage: python golden_additive_colourings_20261017.py [receipt]
Session opus-2026-10-17-S41.
"""
import sys, os, json
from collections import Counter
from pysat.solvers import Solver
from pysat.card import CardEnc


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(here, 'root_positive_transport_20261017g', 'golden_shared.json')
    R = json.load(open(path))
    addr = [tuple(a) for a in R['source_addresses']]
    edges = [tuple(e) for e in R['edges']]
    mixed = set(tuple(e) for e in R['mixed_edges'])
    nA, nC = R['first_axis_vertices'], R['coefficient_vertices']
    pats = R['patterns']
    print('vertices %d, edges %d, mixed %d, |A| = %d, |C| = %d, patterns %d' % (len(addr), len(edges), len(mixed), nA, nC, len(pats)))
    EA, EU, EV, MX = set(), set(), set(), []
    for e in edges:
        (a, b, c), (a2, b2, c2) = addr[e[0]], addr[e[1]]
        if e in mixed:
            MX.append(((a, b, c), (a2, b2, c2)))
        elif b == b2 and c == c2:
            EA.add(tuple(sorted((a, a2))))
        elif a == a2 and c == c2:
            EU.add(tuple(sorted((b, b2))))
        elif a == a2 and b == b2:
            EV.add(tuple(sorted((c, c2))))
        else:
            raise RuntimeError('unclassified non-mixed edge %s' % (e,))
    print('axis graphs: A %d edges, u-axis C %d edges, v-axis C %d edges; mixed %d' % (len(EA), len(EU), len(EV), len(MX)))
    types = Counter((a2 != a, b2 != b, c2 != c) for (a, b, c), (a2, b2, c2) in MX)
    print('mixed-edge coordinate-change types:', dict(types))

    def solve(pat, shared=False):
        top = [0]
        def new():
            top[0] += 1; return top[0]
        X = [[new() for _ in range(6)] for _ in range(nA)]
        Y = [[new() for _ in range(6)] for _ in range(nC)]
        Z = Y if shared else [[new() for _ in range(6)] for _ in range(nC)]
        cl = []
        for V in (X, Y) if shared else (X, Y, Z):
            for row in V:
                cl.append(row[:])
                for i in range(6):
                    for j in range(i):
                        cl.append([-row[i], -row[j]])
        for a, a2 in EA:
            for k in range(6):
                cl.append([-X[a][k], -X[a2][k]])
        for b, b2 in EU:
            for k in range(6):
                cl.append([-Y[b][k], -Y[b2][k]])
        for c, c2 in EV:
            for k in range(6):
                cl.append([-Z[c][k], -Z[c2][k]])
        cl.append([Y[0][0]])
        if not shared:
            cl.append([Z[0][0]])
        for i, col in enumerate(pat):
            cl.append([X[i][col]])               # ports are (i, 0, 0)

        def diff(P, i1, i2):
            if i1 == i2:
                return None
            D = [new() for _ in range(6)]
            for i in range(6):
                for j in range(6):
                    cl.append([-P[i2][i], -P[i1][j], D[(i - j) % 6]])
            return D
        for (a, b, c), (a2, b2, c2) in MX:
            ds = [d for d in (diff(X, a, a2), diff(Y, b, b2), diff(Z, c, c2)) if d is not None]
            if len(ds) == 1:
                cl.append([-ds[0][0]])
            elif len(ds) == 2:
                for i in range(6):
                    cl.append([-ds[0][i], -ds[1][(-i) % 6]])
            else:
                for i in range(6):
                    for j in range(6):
                        cl.append([-ds[0][i], -ds[1][j], -ds[2][(-i - j) % 6]])
        with Solver(name='g3', bootstrap_with=cl) as s:
            if not s.solve():
                return None
            m = set(v for v in s.get_model() if v > 0)
            cA = [next(k for k in range(6) if X[a][k] in m) for a in range(nA)]
            c1 = [next(k for k in range(6) if Y[b][k] in m) for b in range(nC)]
            c2 = [next(k for k in range(6) if Z[c][k] in m) for c in range(nC)]
            return cA, c1, c2

    ok = 0; ok_shared = 0; failures = []
    for pat in pats:
        sol = solve(pat)
        if sol is None:
            failures.append(pat); continue
        cA, c1, c2 = sol
        col = [(cA[a] + c1[b] + c2[c]) % 6 for (a, b, c) in addr]
        assert col[:8] == list(pat)
        assert all(col[x] != col[y] for x, y in edges), 'additive colouring violates an edge'
        ok += 1
        if solve(pat, shared=True) is not None:
            ok_shared += 1
    print('patterns with an additive Z/6 colouring (verified on all %d edges): %d of %d' % (len(edges), ok, len(pats)))
    print('patterns with an additive colouring with c_1 = c_2 (shared coefficient colouring): %d of %d' % (ok_shared, len(pats)))
    json.dump(dict(vertices=len(addr), edges=len(edges), mixed=len(MX), axis=[len(EA), len(EU), len(EV)],
                   mixed_types={str(k): v for k, v in types.items()}, additive=ok, additive_shared=ok_shared,
                   patterns=len(pats), failures=failures),
              open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'golden_additive_colourings_20261017.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
