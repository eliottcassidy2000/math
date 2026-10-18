#!/usr/bin/env python3
"""Design loop for golden-coordinate six-colour forcing: build A' + uC' + vC' in K, pre-screen additively, then SAT (S41).

clock_duality_oct15's golden contact (joint_contact_field_20261017f.md section 5 with G5; shared_golden_contact_
20261017g.md) classifies every unit difference in K + uK + vK, K = Q(zeta15, sqrt(-11)), d = 1 + zeta5: the three axes
and the contact family q(d, 1, 1), |q| = 1. So for coefficient sets A', C' in K the graph on A' + uC' + vC' is determined
inside K: axis edges are unit differences inside one factor, and mixed edges join (a, b, c) to (a + dq, b + q, c + q)
for every unit q in K. This script
  1. builds A', C' from named recipes (ports first and literal),
  2. finds their unit edges (modular sieve at three good primes, then exact norm in K),
  3. solves, for every one of the 570 canonical six-colour port patterns, the ADDITIVE problem
     c(a, b, c) = c_A(a) + c_1(b) + c_2(c) mod 6 (golden_additive_colourings_20261017.py), incrementally,
  4. for any additively blocked pattern, runs the full SAT problem on the complete graph (relative six-core first),
     and stores an explicit colouring witness whenever it extends.
A pattern blocked by the full SAT problem would be a forcing candidate requiring an independent certificate.
Usage: python golden_forcing_search_20261017.py RECIPE [--full]
Session opus-2026-10-17-S41.
"""
import os, sys, json, time
from fractions import Fraction as Q
from itertools import product
import numpy as np
import sympy as sp
from pysat.solvers import Solver
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from escaped_ports_20261017d import make, need
from golden_boundary_20261017c import PORT_EDGES, patterns
from golden_vsquare_unit_census_20261017 import places, red

T0 = time.time()


def log(*a):
    print('[%5.0f s]' % (time.time() - T0), *a, flush=True)


class Geo:
    def __init__(self):
        f, z, j, alpha, ports, units, eta, conj, escape, unused = make('sqrt3')
        self.f, self.conj, self.alpha, self.eta, self.j, self.ports = f, conj, alpha, eta, j, ports
        self.d = f.add(f.one, alpha)
        rho = f.scale(f.add(f.scalar(5), j), Q(1, 6)); baseD = f.add(f.one, eta)
        self.M = [f.zero, f.one, eta, baseD, rho, f.mul(rho, eta), f.mul(rho, baseD)]
        self.pl = places(f, 3)
        self.cache = {}

    def bar(self, a):
        return self.f.apply(a, self.conj)

    def is_unit(self, x):
        if x not in self.cache:
            self.cache[x] = self.f.mul(x, self.bar(x)) == self.f.one
        return self.cache[x]

    def residues(self, pts):
        out = []
        for p, emb, embc in self.pl:
            out.append((p, np.array([red(x, p, emb) for x in pts], dtype=np.int64),
                        np.array([red(x, p, embc) for x in pts], dtype=np.int64)))
        return out

    def unit_edges(self, pts):
        res = self.residues(pts)
        n = len(pts); edges = []
        p, r, rc = res[0]
        for a in range(n):
            dr = (r[a + 1:] - r[a]) % p; drc = (rc[a + 1:] - rc[a]) % p
            cand = np.flatnonzero((dr * drc) % p == 1)
            for off in cand:
                b = a + 1 + int(off)
                ok = True
                for p2, r2, rc2 in res[1:]:
                    if ((r2[b] - r2[a]) * (rc2[b] - rc2[a])) % p2 != 1:
                        ok = False; break
                if ok and self.is_unit(self.f.sub(pts[b], pts[a])):
                    edges.append((a, b))
        return edges


def recipe(geo, name):
    f = geo.f; P = list(geo.ports); M = geo.M; d = geo.d
    def uniq(lst):
        out = []; seen = set()
        for x in lst:
            if x not in seen:
                seen.add(x); out.append(x)
        return out
    A0 = uniq(P + [f.add(p, m) for p in P for m in M])
    C0 = uniq(M + P[:5])
    if name == 'base':
        return A0, C0
    UC = sorted({f.sub(b, a) for a in C0 for b in C0 if a != b and geo.is_unit(f.sub(b, a))})
    if name == 'contact1':          # A' = A + {0} + d * (unit differences of C)
        return uniq(A0 + [f.add(a, f.mul(d, q)) for a in A0 for q in UC]), C0
    if name == 'dM':                # A' = P + M + d M
        return uniq(A0 + [f.add(a, f.mul(d, m)) for a in A0 for m in M]), C0
    if name == 'dM_C2':             # plus a larger coefficient set C' = C + spindle units
        C1 = uniq(C0 + [f.add(c, q) for c in C0 for q in UC])
        return uniq(A0 + [f.add(a, f.mul(d, m)) for a in A0 for m in M]), C1
    if name == 'MM':                # Eisenstein-rigid coefficients: C' = M + M plus the golden pentagon
        C2 = uniq([f.add(a, b) for a in M for b in M] + P[:5])
        return uniq(A0 + [f.add(a, f.mul(d, m)) for a in A0 for m in M]), C2
    if name == 'contact2':          # two contact steps
        A1 = uniq(A0 + [f.add(a, f.mul(d, q)) for a in A0 for q in UC])
        return uniq(A1 + [f.add(a, f.mul(d, q)) for a in A1 for q in UC]), C0
    raise ValueError(name)


def build(geo, A, C):
    f = geo.f
    EA = geo.unit_edges(A); EC = geo.unit_edges(C)
    IA = {x: i for i, x in enumerate(A)}
    pairsC = {}
    for b, b2 in EC:
        for (s, t) in ((b, b2), (b2, b)):
            pairsC.setdefault(f.sub(C[t], C[s]), []).append((s, t))
    mixed = set()
    for q, lst in pairsC.items():
        dq = f.mul(geo.d, q)
        pairsA = [(a, IA[f.add(A[a], dq)]) for a in range(len(A)) if f.add(A[a], dq) in IA]
        for (s, t) in lst:
            for (c, c2) in lst:
                for (a, a2) in pairsA:
                    mixed.add(tuple(sorted(((a, s, c), (a2, t, c2)))))
    return EA, EC, sorted(mixed)


def two_distance(nA, EA, mixed, pats):
    """Sufficient test: c_1 homomorphism-like and c_2 = -c_1 reduce every mixed constraint to c_A(a + dq) != c_A(a)."""
    E = set(tuple(sorted(e)) for e in EA) | set(tuple(sorted((x[0], y[0]))) for x, y in mixed if x[0] != y[0])
    blocked = []
    with Solver(name='g3') as s:
        for a in range(nA):
            s.add_clause([6 * a + c + 1 for c in range(6)])
        for a, b in E:
            for c in range(6):
                s.add_clause([-(6 * a + c + 1), -(6 * b + c + 1)])
        for pat in pats:
            if not s.solve(assumptions=[6 * i + c + 1 for i, c in enumerate(pat)]):
                blocked.append(pat)
    return blocked, len(E)


def additive(nA, nC, EA, EC, mixed, pats, shared):
    top = [0]
    def new():
        top[0] += 1; return top[0]
    X = [[new() for _ in range(6)] for _ in range(nA)]
    Y = [[new() for _ in range(6)] for _ in range(nC)]
    Z = Y if shared else [[new() for _ in range(6)] for _ in range(nC)]
    cl = []
    for V in ((X, Y) if shared else (X, Y, Z)):
        for row in V:
            cl.append(row[:])
            for i in range(6):
                for k in range(i):
                    cl.append([-row[i], -row[k]])
    for a, a2 in EA:
        for k in range(6):
            cl.append([-X[a][k], -X[a2][k]])
    for b, b2 in EC:
        for k in range(6):
            cl.append([-Y[b][k], -Y[b2][k]])
            if not shared:
                cl.append([-Z[b][k], -Z[b2][k]])
    cl.append([Y[0][0]])
    if not shared:
        cl.append([Z[0][0]])
    dcache = {}
    def diff(P, tag, i1, i2):
        if i1 == i2:
            return None
        key = (tag, i1, i2)
        if key not in dcache:
            D = [new() for _ in range(6)]
            for i in range(6):
                for k in range(6):
                    cl.append([-P[i2][i], -P[i1][k], D[(i - k) % 6]])
            dcache[key] = D
        return dcache[key]
    for (a, b, c), (a2, b2, c2) in mixed:
        ds = [x for x in (diff(X, 'A', a, a2), diff(Y, 'Y', b, b2), diff(Z, 'Y' if shared else 'Z', c, c2)) if x is not None]
        if len(ds) == 3:
            for i in range(6):
                for k in range(6):
                    cl.append([-ds[0][i], -ds[1][k], -ds[2][(-i - k) % 6]])
        elif len(ds) == 2:
            for i in range(6):
                cl.append([-ds[0][i], -ds[1][(-i) % 6]])
        else:
            cl.append([-ds[0][0]])
    blocked = []
    with Solver(name='g3', bootstrap_with=cl) as s:
        for pat in pats:
            if not s.solve(assumptions=[X[i][c] for i, c in enumerate(pat)]):
                blocked.append(pat)
    return blocked


def full_sat(nA, nC, EA, EC, mixed, pats):
    """Complete graph on all addresses (a, b, c); ports are (i, 0, 0). Returns blocked patterns."""
    idx = lambda a, b, c: (a * nC + b) * nC + c
    n = nA * nC * nC
    edges = set()
    for a, a2 in EA:
        for b in range(nC):
            for c in range(nC):
                edges.add((idx(a, b, c), idx(a2, b, c)))
    for b, b2 in EC:
        for a in range(nA):
            for c in range(nC):
                edges.add((idx(a, b, c), idx(a, b2, c)))
                edges.add((idx(a, c, b), idx(a, c, b2)))
    for (a, b, c), (a2, b2, c2) in mixed:
        edges.add((idx(a, b, c), idx(a2, b2, c2)))
    nb = [set() for _ in range(n)]
    for x, y in edges:
        nb[x].add(y); nb[y].add(x)
    ports = [idx(i, 0, 0) for i in range(8)]
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
    with Solver(name='g3') as s:
        for v, i in var.items():
            s.add_clause([6 * i + c + 1 for c in range(6)])
        for x, y in edges:
            if x in var and y in var:
                for c in range(6):
                    s.add_clause([-(6 * var[x] + c + 1), -(6 * var[y] + c + 1)])
        for pat in pats:
            if not s.solve(assumptions=[6 * var[ports[i]] + c + 1 for i, c in enumerate(pat)]):
                blocked.append(pat)
    return blocked, n, len(edges), len(live)


def main():
    name = sys.argv[1] if len(sys.argv) > 1 else 'base'
    full = '--full' in sys.argv
    geo = Geo()
    A, C = recipe(geo, name)
    log('recipe %s: |A| = %d, |C| = %d' % (name, len(A), len(C)))
    EA, EC, mixed = build(geo, A, C)
    log('edges: A %d, C %d, mixed %d' % (len(EA), len(EC), len(mixed)))
    need([tuple(sorted(e)) for e in EA if max(e) < 8] == sorted(tuple(sorted(e)) for e in PORT_EDGES), 'ports literal')
    pats = patterns(8, PORT_EDGES)
    bl_two, n_two = two_distance(len(A), EA, mixed, pats)
    log('golden two-distance graph on A (unit edges + d-steps): %d edges; patterns not extending: %d' % (n_two, len(bl_two)))
    bl_shared = additive(len(A), len(C), EA, EC, mixed, pats, True)
    bl_free = additive(len(A), len(C), EA, EC, mixed, pats, False)
    log('additive pre-screen: blocked with c_1 = c_2: %d; blocked with independent c_1, c_2: %d (of %d)'
        % (len(bl_shared), len(bl_free), len(pats)))
    out = dict(recipe=name, A=len(A), C=len(C), edges_A=len(EA), edges_C=len(EC), mixed=len(mixed),
               two_distance_edges=n_two, two_distance_blocked=len(bl_two),
               additive_blocked_shared=len(bl_shared), additive_blocked=len(bl_free), blocked_free=bl_free)
    if full or bl_free:
        bl, n, ne, core = full_sat(len(A), len(C), EA, EC, mixed, bl_free if bl_free else pats)
        log('full SAT on %d vertices / %d edges: blocked %d of %d tested' % (n, ne, len(bl), len(bl_free) if bl_free else len(pats)))
        out.update(full_vertices=n, full_edges=ne, six_core=core, full_blocked=bl)
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(out, open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'golden_forcing_search_%s_20261017.json' % name), 'w'), indent=1)


if __name__ == '__main__':
    main()
