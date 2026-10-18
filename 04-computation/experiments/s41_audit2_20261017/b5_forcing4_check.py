#!/usr/bin/env python3
"""Referee check B5 (second S41 audit): golden_forcing_search4_20261017.py, source A + uC + vC + e^2 C.

1. Classify EVERY unit difference of the source from first principles: float census of the whole box
   (dA, d1, d2, d3) in (D(A) u 0) x (D(C) u 0)^3 minus 0 at the intended embedding, every candidate identified exactly in
   K as an axis unit, the v-contact q(d,1,1,0) or the e^2-contact q(d,1,0,conj d); anything else is reported.
2. Build the full 54,000-vertex graph from those unit vectors by exact address lookup, count edges by type.
3. Rebuild the session's mixV / mixY sets with the session's own construction (copied from search4 lines 37-56, run on
   its Geo object; its main() is NOT called because it writes into 05-knowledge/results) and compare them with the
   referee's contact edges after mapping the session's elements to the referee's representation.
4. Port graph and canonical patterns recomputed; for each of the 570 patterns solve the additive problem with a DIRECT
   encoding (no auxiliary difference variables; one clause per forbidden 6-tuple), shared (c1 = c2) and independent,
   and verify every resulting colouring explicitly on all edges of the referee's graph and on the literal ports.
Usage: python b5_forcing4_check.py
"""
import os, sys, time, itertools, json
import numpy as np
from numba import njit
from pysat.solvers import Solver
import audit2_field as F

T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
EXP = os.path.dirname(HERE)


def log(*a):
    print('[%5.0f s]' % (time.time() - T0), *a, flush=True)


@njit(cache=True)
def scan(XA, X12, X3, tol):
    out = []
    mind = 1e9
    for i3 in range(X3.shape[0]):
        for i12 in range(X12.shape[0]):
            y = X12[i12] + X3[i3]
            for iA in range(XA.shape[0]):
                z = XA[iA] + y
                dev = abs(z.real * z.real + z.imag * z.imag - 1.0)
                if dev < tol:
                    out.append((iA, i12, i3))
                elif dev < mind:
                    mind = dev
    return out, mind


def main():
    S = F.golden_sets()
    A, C, P, d = S['A'], S['C'], S['P'], S['d']
    dbar = F.conj(d)
    nA, nC = len(A), len(C)
    DA = [F.K(0)] + F.diffs(A); DC = [F.K(0)] + F.diffs(C)
    G = F.golden_numbers()
    u, v, e = G['u'], G['v'], G['e']
    Y = e * e
    XA = np.array([F.to_complex(x) for x in DA])
    XC = [F.to_complex(x) for x in DC]
    X12 = np.array([u * a + v * b for a in XC for b in XC])
    X3 = np.array([Y * x for x in XC])
    cand, mind = scan(XA, X12, X3, 1e-8)
    log('box census for Y = e^2 with d3 in D(C) u {0}: %d x %d x %d x %d = %d differences; %d candidates (incl. 0), next deviation %.3e'
        % (len(DA), len(DC), len(DC), len(DC), len(DA) * len(DC) ** 3, len(cand), mind))
    types = {}
    units = []
    for iA, i12, i3 in cand:
        i1, i2 = divmod(i12, len(DC))
        dA, d1, d2, d3 = DA[iA], DC[i1], DC[i2], DC[i3]
        z = [F.is_zero(t) for t in (dA, d1, d2, d3)]
        if all(z):
            continue
        if z.count(False) == 1:
            q = [t for t, zz in zip((dA, d1, d2, d3), z) if not zz][0]
            ty = 'axis' + str(z.index(False)) if F.is_unit(q) else 'other'
        elif F.is_unit(d1) and dA == F.mul(d, d1) and d2 == d1 and z[3]:
            ty = 'v-contact'
        elif F.is_unit(d1) and dA == F.mul(d, d1) and z[2] and d3 == F.mul(dbar, d1):
            ty = 'e2-contact'
        else:
            ty = 'other'
        types[ty] = types.get(ty, 0) + 1
        units.append((ty, dA, d1, d2, d3))
    log('unit difference vectors by type:', types)
    assert 'other' not in types
    # 2. full graph
    IA = {x: i for i, x in enumerate(A)}; IC = {x: i for i, x in enumerate(C)}
    def shift(X, I, dlt):
        return [I.get(F.add(x, dlt)) for x in X]
    cacheA, cacheC = {}, {}
    vid = lambda a, b, c, l: ((a * nC + b) * nC + c) * nC + l
    edges = set(); bytype = {}
    contact_pat = {'v-contact': set(), 'e2-contact': set()}
    for ty, dA, d1, d2, d3 in units:
        sA = cacheA.setdefault(dA, shift(A, IA, dA))
        s1 = cacheC.setdefault(d1, shift(C, IC, d1)); s2 = cacheC.setdefault(d2, shift(C, IC, d2)); s3 = cacheC.setdefault(d3, shift(C, IC, d3))
        pa = [(a, t) for a, t in enumerate(sA) if t is not None]
        p1 = [(a, t) for a, t in enumerate(s1) if t is not None]
        p2 = [(a, t) for a, t in enumerate(s2) if t is not None]
        p3 = [(a, t) for a, t in enumerate(s3) if t is not None]
        cnt = 0
        for (a, a2), (b, b2), (c, c2), (l, l2) in itertools.product(pa, p1, p2, p3):
            x, y = vid(a, b, c, l), vid(a2, b2, c2, l2)
            ed = (min(x, y), max(x, y))
            if ed not in edges:
                edges.add(ed); cnt += 1
                if ty == 'v-contact':
                    contact_pat[ty].add(tuple(sorted(((a, b, c), (a2, b2, c2)))))
                elif ty == 'e2-contact':
                    contact_pat[ty].add(tuple(sorted(((a, b, l), (a2, b2, l2)))))
        bytype[ty] = bytype.get(ty, 0) + cnt
    log('referee graph: %d vertices, %d edges; by type %s; undirected contact patterns: v %d, e^2 %d'
        % (nA * nC ** 3, len(edges), bytype, len(contact_pat['v-contact']), len(contact_pat['e2-contact'])))
    # 3. session construction (search4 lines 37-56) on the session's Geo object
    sys.path.insert(0, EXP)
    from golden_forcing_search_20261017 import Geo
    geo = Geo(); f = geo.f
    def to_mine(x):
        idx = {bb: i for i, bb in enumerate(f.basis)}
        assert all(x[idx[(a, b, c)]] == 0 for (a, b, c) in f.basis if c != 0)
        return (tuple(x[idx[(a, 0, 0)]] for a in range(8)), tuple(x[idx[(a, 1, 0)]] for a in range(8)))
    Ps = list(geo.ports); Ms = geo.M; ds = geo.d; dbs = geo.bar(ds)
    def uniq(lst):
        out = []; seen = set()
        for x in lst:
            if x not in seen:
                seen.add(x); out.append(x)
        return out
    As = uniq(Ps + [f.add(p, m) for p in Ps for m in Ms]); Cs = uniq(Ms + Ps[:5]); Ls = Cs
    same_order = [to_mine(x) for x in As] == A and [to_mine(x) for x in Cs] == C
    log('session A, C equal the referee A, C element by element (same order):', same_order)
    EA = geo.unit_edges(As); EC = geo.unit_edges(Cs); EL = geo.unit_edges(Ls)
    IAs = {x: i for i, x in enumerate(As)}; ILs = {x: i for i, x in enumerate(Ls)}
    unitsS = {}
    for b, b2 in EC:
        for s_, t_ in ((b, b2), (b2, b)):
            unitsS.setdefault(f.sub(Cs[t_], Cs[s_]), []).append((s_, t_))
    mixV, mixY = set(), set()
    for q, cpairs in unitsS.items():
        dq = f.mul(ds, q); dbq = f.mul(dbs, q)
        apairs = [(a, IAs[f.add(As[a], dq)]) for a in range(len(As)) if f.add(As[a], dq) in IAs]
        lpairs = [(l, ILs[f.add(Ls[l], dbq)]) for l in range(len(Ls)) if f.add(Ls[l], dbq) in ILs]
        for (a, a2) in apairs:
            for (s_, t_) in cpairs:
                for (c, c2) in cpairs:
                    mixV.add(((a, s_, c), (a2, t_, c2)))
                for (l, l2) in lpairs:
                    mixY.add(((a, s_, l), (a2, t_, l2)))
    mixVu = {tuple(sorted(x)) for x in mixV}; mixYu = {tuple(sorted(x)) for x in mixY}
    log('session: |EA| %d |EC| %d |EL| %d; directed mixV %d mixY %d; undirected %d, %d'
        % (len(EA), len(EC), len(EL), len(mixV), len(mixY), len(mixVu), len(mixYu)))
    log('session mixV == referee v-contact patterns:', mixVu == contact_pat['v-contact'],
        '; session mixY == referee e^2-contact patterns:', mixYu == contact_pat['e2-contact'])
    # referee axis edges inside A and C, compared with the session's
    myEA = sorted({tuple(sorted((a, IA[F.add(A[a], dA)]))) for ty, dA, d1, d2, d3 in units if ty == 'axis0' for a in range(nA) if F.add(A[a], dA) in IA})
    myEC = sorted({tuple(sorted((b, IC[F.add(C[b], d1)]))) for ty, dA, d1, d2, d3 in units if ty == 'axis1' for b in range(nC) if F.add(C[b], d1) in IC})
    log('axis edges: referee A %d C %d; equal to session EA, EC:' % (len(myEA), len(myEC)),
        myEA == sorted(tuple(sorted(x)) for x in EA), myEC == sorted(tuple(sorted(x)) for x in EC))
    sess_total = len(EA) * nC ** 3 + 3 * len(EC) * nA * nC ** 2 + len(mixVu) * nC + len(mixYu) * nC
    log('session edge count implied by its construction: %d; referee graph: %d' % (sess_total, len(edges)))
    # 4. ports and patterns
    port_edges = sorted({(i, k) for i in range(8) for k in range(i + 1, 8) if F.is_unit(F.sub(P[i], P[k]))})
    pats = []
    def walk(acc):
        if len(acc) == 8:
            pats.append(tuple(acc)); return
        for col in range(min(5, max(acc, default=-1) + 1) + 1):
            b = acc + [col]
            if all(b[i] != b[k] for i, k in port_edges if max(i, k) < len(b)):
                walk(b)
    walk([])
    from golden_boundary_20261017c import PORT_EDGES, patterns
    log('port graph edges (referee):', port_edges, '; equal to PORT_EDGES:', port_edges == sorted(tuple(sorted(x)) for x in PORT_EDGES),
        '; canonical patterns:', len(pats), '; equal to session list:', pats == patterns(8, PORT_EDGES))
    assert all(A[i] == P[i] for i in range(8)) and F.is_zero(C[0])
    E = np.array(sorted(edges), dtype=np.int64)
    aIdx = np.repeat(np.arange(nA), nC ** 3); bIdx = np.tile(np.repeat(np.arange(nC), nC ** 2), nA)
    cIdx = np.tile(np.repeat(np.arange(nC), nC), nA * nC); lIdx = np.tile(np.arange(nC), nA * nC ** 2)
    vmix = sorted(contact_pat['v-contact']); ymix = sorted(contact_pat['e2-contact'])
    for shared in (True, False):
        top = [0]
        def new():
            top[0] += 1; return top[0]
        X = [[new() for _ in range(6)] for _ in range(nA)]
        Y1 = [[new() for _ in range(6)] for _ in range(nC)]
        Y2 = Y1 if shared else [[new() for _ in range(6)] for _ in range(nC)]
        Y3 = [[new() for _ in range(6)] for _ in range(nC)]
        cl = []
        for V in ([X, Y1, Y3] if shared else [X, Y1, Y2, Y3]):
            for row in V:
                cl.append(row[:])
                for i in range(6):
                    for k in range(i):
                        cl.append([-row[i], -row[k]])
        for a, a2 in myEA:
            for k in range(6):
                cl.append([-X[a][k], -X[a2][k]])
        for V in (Y1, Y2, Y3):
            for b, b2 in myEC:
                for k in range(6):
                    cl.append([-V[b][k], -V[b2][k]])
        cl += [[Y1[0][0]], [Y2[0][0]], [Y3[0][0]]]
        tuples6 = [t for t in itertools.product(range(6), repeat=6) if ((t[1] - t[0]) + (t[3] - t[2]) + (t[5] - t[4])) % 6 == 0]
        for (a, b, c), (a2, b2, c2) in vmix:
            for t in tuples6:
                cl.append([-X[a][t[0]], -X[a2][t[1]], -Y1[b][t[2]], -Y1[b2][t[3]], -Y2[c][t[4]], -Y2[c2][t[5]]])
        for (a, b, l), (a2, b2, l2) in ymix:
            for t in tuples6:
                cl.append([-X[a][t[0]], -X[a2][t[1]], -Y1[b][t[2]], -Y1[b2][t[3]], -Y3[l][t[4]], -Y3[l2][t[5]]])
        blocked = 0; verified = 0
        with Solver(name='cd15', bootstrap_with=cl) as sol:
            for pat in pats:
                if not sol.solve(assumptions=[X[i][c] for i, c in enumerate(pat)]):
                    blocked += 1; continue
                m = set(x for x in sol.get_model() if x > 0)
                cA = np.array([next(k for k in range(6) if X[a][k] in m) for a in range(nA)])
                c1 = np.array([next(k for k in range(6) if Y1[b][k] in m) for b in range(nC)])
                c2 = np.array([next(k for k in range(6) if Y2[b][k] in m) for b in range(nC)])
                c3 = np.array([next(k for k in range(6) if Y3[b][k] in m) for b in range(nC)])
                col = (cA[aIdx] + c1[bIdx] + c2[cIdx] + c3[lIdx]) % 6
                ok = bool(np.all(col[E[:, 0]] != col[E[:, 1]])) and tuple(int(col[vid(i, 0, 0, 0)]) for i in range(8)) == pat
                verified += ok
        log('additive (%s): clauses %d; blocked %d of %d; colourings verified on all %d referee edges and literal ports: %d'
            % ('c1 = c2' if shared else 'independent', len(cl), blocked, len(pats), len(edges), verified))


if __name__ == '__main__':
    main()
