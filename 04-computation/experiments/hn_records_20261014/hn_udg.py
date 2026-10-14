#!/usr/bin/env python3
"""Exact unit-distance graphs over a real multiquadratic field, SAT colourability and core-based minimisation.

Field.  F = Q(sqrt p : p in PRIMES), a real multiquadratic field of degree 2^len(PRIMES).  An element is a tuple of
2^k Fractions; coordinate b (a bitmask) multiplies sqrt(prod of PRIMES[i] with bit i set in b).  A point of the plane is
a pair (x, y) of elements, i.e. a complex number x + iy with x, y in F.

Exactness.  Edges are found by a floating-point candidate search (|d - 1| < 1e-9) and every candidate is then decided by
exact arithmetic, |p - q|^2 == 1 in F.  Every exact unit pair passes the float filter (float error ~ 1e-13), so the edge
set is complete for the given vertex set, and every edge is an exact unit distance.

Colouring.  k-colourability is a CNF with variables x[v][c] and one selector s[v] per vertex: (-s_v | x_v1 | ... | x_vk)
and (-x_uc | -x_vc) for every edge and colour.  A vertex whose selector is off is free (all its colour variables can be
false).  Unsatisfiable cores over the selectors give non-k-colourable induced subgraphs; repeated deletion gives a
vertex-critical one.
"""
from fractions import Fraction as Fr
from itertools import combinations
import math
import random

import numpy as np
from scipy.spatial import cKDTree
from pysat.solvers import Solver

import os
PRIMES = tuple(int(x) for x in os.environ.get('HN_PRIMES', '3,5,11').split(','))
K = len(PRIMES)
DIM = 1 << K
BASIS_FLOAT = []
for b in range(DIM):
    v = 1
    for i in range(K):
        if b >> i & 1:
            v *= PRIMES[i]
    BASIS_FLOAT.append(math.sqrt(v))
BASIS_FLOAT = np.array(BASIS_FLOAT)
# product table: sqrt(S) * sqrt(T) = (prod of primes in S & T) * sqrt(S ^ T)
MULC = [[1] * DIM for _ in range(DIM)]
for a in range(DIM):
    for b in range(DIM):
        c = 1
        for i in range(K):
            if (a & b) >> i & 1:
                c *= PRIMES[i]
        MULC[a][b] = c

ZERO = tuple(Fr(0) for _ in range(DIM))
ONE = tuple(Fr(1) if b == 0 else Fr(0) for b in range(DIM))


def el(**kw):
    """Field element from keyword radicands, e.g. el(r1=Fr(1,2), r3=Fr(1,2)) = 1/2 + sqrt3/2."""
    out = [Fr(0)] * DIM
    for key, val in kw.items():
        r = int(key[1:])
        b = 0
        rem = r
        for i, p in enumerate(PRIMES):
            if rem % p == 0:
                rem //= p
                b |= 1 << i
        assert rem == 1, f"radicand {r} not a product of distinct {PRIMES}"
        out[b] += Fr(val)
    return tuple(out)


def fadd(a, b):
    return tuple(x + y for x, y in zip(a, b))


def fsub(a, b):
    return tuple(x - y for x, y in zip(a, b))


def fneg(a):
    return tuple(-x for x in a)


def fmul(a, b):
    out = [Fr(0)] * DIM
    for i in range(DIM):
        ai = a[i]
        if not ai:
            continue
        for j in range(DIM):
            bj = b[j]
            if bj:
                out[i ^ j] += ai * bj * MULC[i][j]
    return tuple(out)


def ffloat(a):
    return float(sum(float(a[b]) * BASIS_FLOAT[b] for b in range(DIM)))


# complex numbers z = (x, y)
def cadd(z, w):
    return (fadd(z[0], w[0]), fadd(z[1], w[1]))


def csub(z, w):
    return (fsub(z[0], w[0]), fsub(z[1], w[1]))


def cmul(z, w):
    return (fsub(fmul(z[0], w[0]), fmul(z[1], w[1])), fadd(fmul(z[0], w[1]), fmul(z[1], w[0])))


def cconj(z):
    return (z[0], fneg(z[1]))


def cnorm(z):
    return fadd(fmul(z[0], z[0]), fmul(z[1], z[1]))


def cfloat(z):
    return (ffloat(z[0]), ffloat(z[1]))


C0 = (ZERO, ZERO)
C1 = (ONE, ZERO)


def is_unit(z):
    return cnorm(z) == ONE


def cpow(z, n):
    if n < 0:
        z, n = cconj(z), -n   # unit z: z^-1 = conj(z)
        assert is_unit(z)
    out = C1
    for _ in range(n):
        out = cmul(out, z)
    return out


# ---------------------------------------------------------------- point sets
def pset(points):
    return frozenset(points)


def rotate(S, z):
    return pset(cmul(p, z) for p in S)


def translate(S, t):
    return pset(cadd(p, t) for p in S)


def minkowski(A, B):
    return pset(cadd(a, b) for a in A for b in B)


def union(*sets):
    out = set()
    for S in sets:
        out |= S
    return pset(out)


# ---------------------------------------------------------------- graphs
def unit_graph(points, tol=1e-9):
    """Return (pts list, float array, edges) of the exact unit-distance graph on the given points."""
    pts = sorted(points, key=lambda p: cfloat(p))
    xy = np.array([cfloat(p) for p in pts])
    tree = cKDTree(xy)
    cand = tree.query_pairs(1.0 + 1e-6, output_type='ndarray')
    if len(cand) == 0:
        return pts, xy, []
    d = np.hypot(*(xy[cand[:, 0]] - xy[cand[:, 1]]).T)
    cand = cand[np.abs(d - 1.0) < tol]
    edges = []
    for i, j in cand:
        if is_unit(csub(pts[i], pts[j])):
            edges.append((int(i), int(j)))
    return pts, xy, sorted(edges)


def check_float_gap(xy, edges, tol=1e-9):
    """Sanity: no non-edge pair is within tol of unit distance (would indicate a float/exact disagreement)."""
    tree = cKDTree(xy)
    cand = tree.query_pairs(1.0 + 1e-6, output_type='ndarray')
    d = np.hypot(*(xy[cand[:, 0]] - xy[cand[:, 1]]).T)
    near = cand[np.abs(d - 1.0) < tol]
    E = set(edges)
    return sum(1 for i, j in near if (int(i), int(j)) not in E)


def triangles_through(n, edges, limit=None):
    adj = [set() for _ in range(n)]
    for i, j in edges:
        adj[i].add(j)
        adj[j].add(i)
    out = []
    for i, j in edges:
        for k in adj[i] & adj[j]:
            if k > j:
                out.append((i, j, k))
                if limit and len(out) >= limit:
                    return out
    return out


class Colourer:
    """Incremental k-colourability oracle with vertex selectors."""

    def __init__(self, n, edges, k=4, solver='cadical153', fixed_triangle=None):
        self.n, self.k = n, k
        self.edges = list(edges)
        self.solver = Solver(name=solver)
        x = lambda v, c: v * k + c + 1
        self.x = x
        self.sel = [n * k + v + 1 for v in range(n)]
        for v in range(n):
            self.solver.add_clause([-self.sel[v]] + [x(v, c) for c in range(k)])
        for (u, v) in edges:
            for c in range(k):
                self.solver.add_clause([-x(u, c), -x(v, c)])
        self.fixed = ()
        if fixed_triangle is not None:
            # symmetry breaking: valid only while these vertices stay active
            self.fixed = tuple(fixed_triangle)
            for c, v in enumerate(fixed_triangle):
                self.solver.add_clause([x(v, c)])

    def colourable(self, active):
        act = set(active) | set(self.fixed)
        ok = self.solver.solve(assumptions=[self.sel[v] for v in sorted(act)])
        if ok:
            return True, None
        core = self.solver.get_core() or []
        sel_index = {s: v for v, s in enumerate(self.sel)}
        return False, sorted(set(sel_index[s] for s in core if s in sel_index) | set(self.fixed))

    def colouring(self, active):
        act = set(active) | set(self.fixed)
        if not self.solver.solve(assumptions=[self.sel[v] for v in sorted(act)]):
            return None
        model = set(l for l in self.solver.get_model() if l > 0)
        col = {}
        for v in act:
            cs = [c for c in range(self.k) if self.x(v, c) in model]
            col[v] = cs[0]
        return col

    def close(self):
        self.solver.delete()


def shrink_to_core(col, active, rounds=50):
    """Iterate core extraction until stable. Returns (colourable?, active)."""
    active = sorted(active)
    for _ in range(rounds):
        ok, core = col.colourable(active)
        if ok:
            return True, active
        if len(core) >= len(active):
            return False, active
        active = core
    return False, active


def vertex_critical(col, active, order=None, rng=None):
    """Greedy deletion with core shrinking; the result is non-colourable and every single deletion makes it colourable."""
    ok, active = shrink_to_core(col, active)
    assert not ok
    active = set(active)
    cand = list(order) if order is not None else sorted(active)
    if rng is not None:
        rng.shuffle(cand)
    for v in cand:
        if v not in active or v in col.fixed:
            continue
        trial = sorted(active - {v})
        ok, core = col.colourable(trial)
        if not ok:
            active = set(core)
    # confirm criticality
    return sorted(active)


def is_vertex_critical(col, active):
    active = sorted(active)
    ok, _ = col.colourable(active)
    if ok:
        return False
    for v in active:
        ok, _ = col.colourable([u for u in active if u != v])
        if not ok:
            return False
    return True


def induced_edges(edges, keep):
    keep = set(keep)
    return [(i, j) for (i, j) in edges if i in keep and j in keep]


# ---------------------------------------------------------------- standard pieces
SQ3_2 = el(r3=Fr(1, 2))
HALF = el(r1=Fr(1, 2))
HEX = [(ONE, ZERO), (HALF, SQ3_2), (fneg(HALF), SQ3_2), (fneg(ONE), ZERO), (fneg(HALF), fneg(SQ3_2)), (HALF, fneg(SQ3_2))]
H = pset([C0] + HEX)                                              # the hexagonal wheel W7
ETA = (el(r33=Fr(1, 6)), el(r3=Fr(1, 6))) if 11 in PRIMES and 3 in PRIMES else None  # exp((i/2) arccos(5/6))
RHO = (el(r1=Fr(7, 8)), el(r15=Fr(1, 8))) if 5 in PRIMES else None  # exp(i arccos(7/8)); needs sqrt5
I_OVER_SQ3 = (ZERO, el(r3=Fr(1, 3)))                              # i / sqrt3 (not unit)
I_SQ3 = (ZERO, el(r3=Fr(1)))                                      # i sqrt3
assert (ETA is None or is_unit(ETA)) and (RHO is None or is_unit(RHO))
if ETA is not None:
    assert cmul(ETA, ETA) == (el(r1=Fr(5, 6)), el(r11=Fr(1, 6)))  # eta^2 = (5 + i sqrt11)/6


def H_rot(alphas):
    return union(*[rotate(H, cpow(ETA, a)) for a in alphas])
