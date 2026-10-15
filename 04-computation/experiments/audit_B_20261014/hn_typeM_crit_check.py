#!/usr/bin/env python3
"""Audit B: independent rebuild of Parts' type-M graph L727 u rho S361 and check of a vertex-critical subgraph
(crit_*.json of the session).

Representation (different from the session's real coordinates in Q(sqrt3, sqrt5, sqrt11)): a point is a complex
number in K = Q(sqrt-3, sqrt-11, sqrt5), stored as 8 rationals on the basis e_b (b = bitmask over the generators
g0 = sqrt-3, g1 = sqrt-11, g2 = sqrt5), e_a e_b = (prod of g_i^2 over common bits) e_(a xor b), with g0^2 = -3,
g1^2 = -11, g2^2 = 5.  Complex conjugation: g0 -> -g0, g1 -> -g1, g2 -> g2.  |z|^2 = z conj(z).

Pieces (as described in the session's note, after Parts arXiv:2010.12665):
  H = {0} u {6th roots of unity},  eta = (sqrt33 + i sqrt3)/6,  H^m = U_{|a|<=m} eta^a H,
  V25 = H^1 u (i/sqrt3) H,  V31 = H^1 u i sqrt3 (eta^-1 H u eta H),  V37 = H^1 u (i/sqrt3) H^1,
  V49 = H^2 u (i/sqrt3) H^1,  L727 = V49 + V37,  S361 = V31 + V25,  G = L727 u rho S361,  rho = (7 + i sqrt15)/8.
Edges: ALL pairs decided exactly (no floating-point filter at all).

Then for the crit JSON: map its exact coordinates into K, check membership in G, recompute the induced unit graph,
compare with the listed edges, test non-4-colourability with two solvers other than CaDiCaL, and vertex-criticality
by exhibiting and checking a proper 4-colouring of G - v for every v.

Usage: python3 -u hn_typeM_crit_check.py crit.json [all|unsat|vc|g]  (g = glucose4 UNSAT check + vertex-criticality)
"""
import sys
import json
import time
from fractions import Fraction as Fr
from math import lcm
import cmath

import numpy as np
from pysat.solvers import Solver

T0 = time.time()


def log(*a):
    print(*a, f"[{time.time() - T0:.1f}s]", flush=True)


SQ = (-3, -11, 5)
VAL = []
for b in range(8):
    v = 1 + 0j
    for i in range(3):
        if b >> i & 1:
            v *= cmath.sqrt(SQ[i])
    VAL.append(v)
MC = [[1] * 8 for _ in range(8)]
for a in range(8):
    for b in range(8):
        c = 1
        for i in range(3):
            if (a & b) >> i & 1:
                c *= SQ[i]
        MC[a][b] = c
CONJS = [(-1) ** ((b & 1) + (b >> 1 & 1)) for b in range(8)]


def mul(x, y):
    out = [Fr(0)] * 8
    for a in range(8):
        if x[a]:
            for b in range(8):
                if y[b]:
                    out[a ^ b] += x[a] * y[b] * MC[a][b]
    return tuple(out)


def add(x, y):
    return tuple(u + v for u, v in zip(x, y))


def conj(x):
    return tuple(c * s for c, s in zip(x, CONJS))


def val(x):
    return sum(complex(float(c)) * VAL[b] for b, c in enumerate(x))


def E(**kw):
    out = [Fr(0)] * 8
    for k, v in kw.items():
        out[int(k[1:], 2)] += Fr(v)
    return tuple(out)


ONE = E(e000=1)
ZERO = tuple([Fr(0)] * 8)
w6 = E(e000=Fr(1, 2), e001=Fr(1, 2))                       # (1 + sqrt-3)/2 = exp(i pi/3)
eta = E(e011=Fr(-1, 6), e001=Fr(1, 6))                     # (sqrt33 + i sqrt3)/6, sqrt33 = -sqrt-3 sqrt-11
rho = E(e000=Fr(7, 8), e101=Fr(1, 8))                      # (7 + i sqrt15)/8, i sqrt15 = sqrt-3 sqrt5
i_over_sqrt3 = E(e001=Fr(1, 3))
i_sqrt3 = E(e001=1)
assert abs(val(w6) - cmath.exp(1j * cmath.pi / 3)) < 1e-12
assert abs(val(eta) - (33 ** 0.5 + 1j * 3 ** 0.5) / 6) < 1e-12
assert abs(val(rho) - (7 + 1j * 15 ** 0.5) / 8) < 1e-12
for u in (w6, eta, rho):
    assert mul(u, conj(u)) == ONE
assert mul(eta, eta) == E(e000=Fr(5, 6), e010=Fr(1, 6))     # (5 + i sqrt11)/6


def pw(u, n):
    if n < 0:
        u, n = conj(u), -n
    out = ONE
    for _ in range(n):
        out = mul(out, u)
    return out


Hs = frozenset([ZERO] + [pw(w6, k) for k in range(6)])


def rot(S, u):
    return frozenset(mul(p, u) for p in S)


def H_(alphas):
    out = set()
    for a in alphas:
        out |= rot(Hs, pw(eta, a))
    return frozenset(out)


def mink(A, B):
    return frozenset(add(a, b) for a in A for b in B)


H1, H2 = H_([-1, 0, 1]), H_([-2, -1, 0, 1, 2])
V25 = H1 | rot(Hs, i_over_sqrt3)
V31 = H1 | rot(H_([-1, 1]), i_sqrt3)
V37 = H1 | rot(H1, i_over_sqrt3)
V49 = H2 | rot(H1, i_over_sqrt3)
L = mink(V49, V37)
S = mink(V31, V25)
G = L | rot(S, rho)
log(f"|H1| = {len(H1)}, |H2| = {len(H2)}, |V25| = {len(V25)}, |V31| = {len(V31)}, |V37| = {len(V37)}, "
    f"|V49| = {len(V49)}, |L727| = {len(L)}, |S361| = {len(S)}, |G| = {len(G)}")

pts = sorted(G, key=lambda p: (val(p).real, val(p).imag))
D = 1
for p in pts:
    for c in p:
        D = lcm(D, c.denominator)
Xi = np.array([[int(c * D) for c in p] for p in pts], dtype=np.int64)
MCa = np.array(MC, dtype=np.int64)
CS = np.array(CONJS, dtype=np.int64)


def norm_rows(Dd):
    """exact z conj(z) for integer rows (over D): returns integer rows over D^2"""
    C = Dd * CS
    out = np.zeros_like(Dd)
    for a in range(8):
        for b in range(8):
            out[:, a ^ b] += MCa[a, b] * Dd[:, a] * C[:, b]
    return out


def unit_pairs(Xs):
    n = len(Xs)
    target = np.zeros(8, dtype=np.int64)
    target[0] = D * D
    I, J = [], []
    for i in range(n - 1):
        dd = Xs[i + 1:] - Xs[i]
        ok = np.all(norm_rows(dd) == target, axis=1)
        jj = np.nonzero(ok)[0] + i + 1
        I.extend([i] * len(jj))
        J.extend(jj.tolist())
    return I, J


I, J = unit_pairs(Xi)
log(f"L727 u rho S361: |V| = {len(pts)}, exact unit edges (all pairs, no float filter) = {len(I)}")

if len(sys.argv) > 2 and sys.argv[2] == 'full':
    # the full type-M graph L727 u rho S361, two solvers other than CaDiCaL, triangle through the origin fixed
    n0 = len(pts)
    orig = [i for i, p in enumerate(pts) if p == ZERO][0]
    adj = [set() for _ in range(n0)]
    for a, b in zip(I, J):
        adj[a].add(b); adj[b].add(a)
    u = min(adj[orig]); w = min(adj[orig] & adj[u])
    for name in ('glucose4', 'kissat404'):
        s_ = Solver(name=name)
        x = lambda v, c: v * 4 + c + 1
        for v in range(n0):
            s_.add_clause([x(v, c) for c in range(4)])
        for a, b in zip(I, J):
            for c in range(4):
                s_.add_clause([-x(a, c), -x(b, c)])
        for c, v in enumerate((orig, u, w)):
            s_.add_clause([x(v, c)])
        t = time.time()
        log(f"L727 u rho S361 4-colourable ({name}, triangle {orig, u, w} fixed): {s_.solve()}  ({time.time() - t:.1f}s)")
        s_.delete()
    sys.exit(0)

# ------------------------------------------------------------- crit JSON
crit = json.load(open(sys.argv[1]))
primes = crit['spec']['primes']
assert primes == [3, 5, 11], primes
# session basis index b (bits over 3, 5, 11): real radical sqrt(prod).  Map x (real) and i*y into K.
#  real: 1 -> e000, sqrt5 -> e100, sqrt33 -> -e011, sqrt165 -> -e111   (sqrt3, sqrt11, sqrt15, sqrt55 must vanish)
#  i*(real): i sqrt3 -> e001, i sqrt11 -> e010, i sqrt15 -> e101, i sqrt55 -> e110 (others must vanish in y)
REAL = {0: (0b000, 1), 2: (0b100, 1), 5: (0b011, -1), 7: (0b111, -1)}      # session bits: 1=3, 2=5, 4=11
IMAG = {1: (0b001, 1), 4: (0b010, 1), 3: (0b101, 1), 6: (0b110, 1)}


def to_K(xy):
    x = [Fr(s) for s in xy[0]]
    y = [Fr(s) for s in xy[1]]
    out = [Fr(0)] * 8
    for b, c in enumerate(x):
        if c:
            assert b in REAL, ("x has a radical outside K", b)
            t, s = REAL[b]
            out[t] += s * c
    for b, c in enumerate(y):
        if c:
            assert b in IMAG, ("y has a radical outside K", b)
            t, s = IMAG[b]
            out[t] += s * c
    return tuple(out)


cv = [to_K(v) for v in crit['vertices']]
Gset = set(G)
inG = sum(1 for p in cv if p in Gset)
fromL = sum(1 for p in cv if p in L)
log(f"crit: n = {crit['n']}, listed edges {crit['edges']}; vertices in G: {inG}/{len(cv)} "
    f"(in L727: {fromL}, others from rho S361: {len(cv) - fromL}); distinct: {len(set(cv))}")
Xc = np.array([[int(c * D) for c in p] for p in cv], dtype=np.int64)
Ic, Jc = unit_pairs(Xc)
mine = set(zip(Ic, Jc))
listed = set((min(a, b), max(a, b)) for a, b in crit['edge_list'])
log(f"crit induced unit edges recomputed: {len(mine)}; equal to the listed edge list: {mine == listed}")
n = len(cv)
Ec = sorted(mine)


def solve(active, name, want_model=False):
    k = 4
    s = Solver(name=name)
    x = lambda v, c: v * k + c + 1
    act = set(active)
    for v in act:
        s.add_clause([x(v, c) for c in range(k)])
    for a, b in Ec:
        if a in act and b in act:
            for c in range(k):
                s.add_clause([-x(a, c), -x(b, c)])
    r = s.solve()
    col = None
    if r and want_model:
        m = set(l for l in s.get_model() if l > 0)
        col = {v: [c for c in range(k) if x(v, c) in m][0] for v in act}
    s.delete()
    return r, col


MODE = sys.argv[2] if len(sys.argv) > 2 else 'all'
if MODE in ('all', 'unsat', 'g'):
    for name in (('glucose4',) if MODE == 'g' else ('glucose4', 'minisat22', 'lingeling')):
        t = time.time()
        r, _ = solve(range(n), name)
        log(f"crit graph 4-colourable ({name}): {r}  ({time.time() - t:.1f}s)")
if MODE == 'unsat':
    sys.exit(0)

bad = []
for v in range(n):
    act = [u for u in range(n) if u != v]
    r, col = solve(act, 'glucose4', want_model=True)
    if not r:
        bad.append(v)
        continue
    mono = sum(1 for a, b in Ec if a != v and b != v and col[a] == col[b])
    assert mono == 0
    if v % 50 == 0:
        log(f"   ... G - v checked for v < {v + 1}")
log(f"vertex-critical: every G - v has a CHECKED proper 4-colouring: {not bad}  (failures: {bad[:10]})")
