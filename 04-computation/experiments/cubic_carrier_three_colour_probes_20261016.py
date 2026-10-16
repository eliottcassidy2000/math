#!/usr/bin/env python3
"""Three-colourability probes of finite unit-distance graphs in the cubic carrier L' = Q(zeta55, theta),
2 theta^3 + 2 theta^2 - 2 theta - 1 = 0 (clock_duality_oct15's G7), with the width-three unit
v = (theta - theta xi + xi^2)/(theta - theta xi^-1 + xi^-2), xi = zeta11, and its conjugates sigma_m(v),
m in (Z/11)^x (sigma_m: zeta11 -> zeta11^m, theta fixed; sigma_(-m)(v) = conj(v) = 1/v).

Exact arithmetic: elements of L' are coefficient arrays over the basis zeta55^i theta^k (0 <= i < 40, 0 <= k < 3);
points are integer combinations of directions, stored as integer vectors scaled by the common denominator, so point
identity is exact.  All unit pairs of a probe (not only generator steps) are found by a three-embedding filter
(every embedding of the CM field L' preserves norm one) and then each is confirmed by the exact norm
(x - y)(bar x - bar y) = 1 in L'.  Three-colourability by CaDiCaL; a colouring is re-verified on every exact edge;
a non-three-colourable probe would be reported (none occurred).  Also: the number of unit edges that are not
generator steps and the Q-ranks of the direction sets (the ranks explain why the probes are Cartesian: their direction
sets are linearly independent), and short integer relations coupling rotated conjugates zeta^k v_r to the roots of
unity, found by kernel-LLL on the saturated relation lattice ([10^12 A | I]); an earlier version reduced a rational
null-space basis, whose integer span has index about 10^107 to 10^192 in the relation lattice, and reported lengths
about 1000 times too long (MISTAKE-622).
FINITE-EXACT for the probes as run.  Triangle-free: checked on every probe.

Session opus-2026-10-16-S38.  Usage: python cubic_carrier_three_colour_probes_20261016.py [--big]
"""
import sys, json, time, math, itertools
import numpy as np
import sympy as sp
from flint import fmpq_poly, fmpq
from scipy.spatial import cKDTree
from pysat.solvers import Solver

BIG = '--big' in sys.argv
z, th, j = sp.symbols('z th j')
PHI55 = sp.Poly(sp.cyclotomic_poly(55, z), z)
NZ = 40

# ---------- exact field ----------
PHI = fmpq_poly([int(c) for c in reversed(PHI55.all_coeffs())])


def zpoly_from_sympy(expr):
    P = sp.Poly(sp.expand(expr), z)
    R = P.rem(PHI55) if P.degree() >= NZ else P
    cs = R.all_coeffs()[::-1]
    return [sp.Rational(c) for c in cs] + [sp.Rational(0)] * (NZ - len(cs))


def to_fmpq(coeffs):
    return fmpq_poly([fmpq(int(sp.Rational(c).p), int(sp.Rational(c).q)) for c in coeffs])


def zconj(p):
    """zeta -> zeta^-1 on a polynomial of degree < 40, reduced mod Phi55."""
    cs = [p[i] for i in range(NZ)]
    out = [fmpq(0)] * 55
    for i, c in enumerate(cs):
        out[(-i) % 55] += c
    return fmpq_poly(out) % PHI


def mul(a, b):
    c = [fmpq_poly(0)] * 5
    for i in range(3):
        for k in range(3):
            c[i + k] = c[i + k] + a[i] * b[k]
    c = [x % PHI for x in c]
    half = fmpq(1, 2)
    return [c[0] + c[3] * half - c[4] * half, c[1] + c[3] - c[4] * half, c[2] - c[3] + c[4] * 2]


def conj(a):
    return [zconj(x) for x in a]


def is_norm_one(a):
    n = mul(a, conj(a))
    return n[0] == 1 and n[1] == 0 and n[2] == 0


# ---------- the unit v (exact, from the G7 receipt) ----------
coef_str = json.load(open('05-knowledge/results/elliptic_cubic_direction_20261016c.json'))['unit_coefficients']
coef = [sp.sympify(s.replace('theta', 'th')) for s in coef_str]
gauss = sum(sp.legendre_symbol(i, 11) * z**(5 * i) for i in range(1, 11))     # sqrt(-11) as a Gauss sum in zeta11
v_expr = sum(c.subs(j, gauss) * z**(5 * h) for h, c in enumerate(coef))
v_expr = sp.expand(v_expr)
# split by theta-degree
vt = [sp.expand(v_expr).coeff(th, k) for k in range(3)]
assert sp.expand(v_expr - sum(vt[k] * th**k for k in range(3))) == 0
v_el = [to_fmpq(zpoly_from_sympy(vt[k])) for k in range(3)]
assert is_norm_one(v_el), 'v must have norm one'


def sigma(a, m):
    """zeta55 -> zeta55^m."""
    out = []
    for p in a:
        cs = [fmpq(0)] * 55
        for i in range(NZ):
            cs[(i * m) % 55] += p[i]
        out.append(fmpq_poly(cs) % PHI)
    return out


def root_of_unity(k):
    cs = [fmpq(0)] * 55
    cs[k % 55] = fmpq(1)
    return [fmpq_poly(cs) % PHI, fmpq_poly(0), fmpq_poly(0)]


def neg(a):
    return [-x for x in a]


# ---------- integer vectors ----------
DEN = 1541


def to_vec(a):
    out = np.zeros(3 * NZ, dtype=np.int64)
    for k in range(3):
        for i in range(NZ):
            c = a[k][i] * DEN
            assert c.q == 1, 'denominator must divide 1541'
            out[k * NZ + i] = int(c.p)
    return out


def from_vec(vec):
    return [fmpq_poly([fmpq(int(vec[k * NZ + i]), DEN) for i in range(NZ)]) for k in range(3)]


# embeddings: zeta55 -> exp(2 pi i a/55), theta -> one of the three real roots
troots = sorted(float(r) for r in sp.real_roots(2 * th**3 + 2 * th**2 - 2 * th - 1))
EMB = [(1, troots[2]), (7, troots[0]), (13, troots[1])]
BASIS = np.array([[np.exp(2j * np.pi * a * i / 55) * t**k for k in range(3) for i in range(NZ)] for (a, t) in EMB]).T


def emb(vecs):
    return (vecs.astype(np.float64) @ BASIS) / DEN


# ---------- directions ----------
ms = [m for m in range(1, 55) if math.gcd(m, 55) == 1 and m % 5 == 1]        # sigma acts on zeta11 only
vconj = {}
for m in ms:
    vconj[m % 11] = sigma(v_el, m)
assert len(vconj) == 10


def dirset(spec):
    dirs = []
    if 'mu10' in spec:
        dirs += [root_of_unity(11 * k) for k in range(5)]
    if 'mu22' in spec:
        dirs += [root_of_unity(5 * k) for k in range(11)]
    if 'mu110' in spec:
        dirs += [root_of_unity(k) for k in range(55)]
    if 'v' in spec:
        dirs += [vconj[r] for r in sorted(vconj)]
    if 'v5' in spec:          # v-conjugates times fifth roots of unity
        dirs += [mul(vconj[r], root_of_unity(11 * k)) for r in sorted(vconj) for k in range(1, 5)]
    dirs = dirs + [neg(d) for d in dirs]
    vecs = []
    seen = set()
    for d in dirs:
        vv = to_vec(d)
        t = tuple(vv)
        if t not in seen:
            seen.add(t)
            vecs.append(vv)
    return np.array(vecs)


def ball(D, radius, cap):
    pts = {tuple(np.zeros(3 * NZ, dtype=np.int64)): 0}
    order = [np.zeros(3 * NZ, dtype=np.int64)]
    frontier = [order[0]]
    for r in range(radius):
        new = []
        for p in frontier:
            for d in D:
                q = p + d
                t = tuple(q)
                if t not in pts:
                    pts[t] = len(order)
                    order.append(q)
                    new.append(q)
                    if len(order) >= cap:
                        return np.array(order)
        frontier = new
    return np.array(order)


def unit_edges(P):
    X = emb(P)
    R6 = np.concatenate([X.real, X.imag], axis=1)
    tree = cKDTree(R6)
    pairs = tree.query_pairs(math.sqrt(3) + 1e-7, output_type='ndarray')
    if len(pairs) == 0:
        return np.zeros((0, 2), dtype=np.int64), 0
    W = X[pairs[:, 0]] - X[pairs[:, 1]]
    ok = np.all(np.abs(np.abs(W) - 1) < 1e-8, axis=1)
    cand = pairs[ok]
    exact = []
    for (a, b) in cand:
        if is_norm_one(from_vec(P[a] - P[b])):
            exact.append((a, b))
    return np.array(exact, dtype=np.int64).reshape(-1, 2), len(cand)


def three_colour(n, E):
    with Solver(name='cadical195') as S:
        var = lambda v, c: 3 * v + c + 1
        for v in range(n):
            S.add_clause([var(v, c) for c in range(3)])
        for a, b in E.tolist():
            for c in range(3):
                S.add_clause([-var(a, c), -var(b, c)])
        S.add_clause([var(0, 0)])
        ok = S.solve()
        if not ok:
            return None
        m = S.get_model()
        col = [next(c for c in range(3) if m[var(v, c) - 1] > 0) for v in range(n)]
    return col


def triangles(n, E):
    adj = [set() for _ in range(n)]
    for a, b in E:
        adj[a].add(b); adj[b].add(a)
    return sum(len(adj[a] & adj[b]) for a, b in E) // 3


PROBES = [('mu10+v', ['mu10', 'v'], 3, 60000), ('mu22+v', ['mu22', 'v'], 3, 60000),
          ('mu110+v', ['mu110', 'v'], 2, 60000), ('mu10+v+v5', ['mu10', 'v', 'v5'], 2, 60000)]
if BIG:
    PROBES += [('mu10+v r4', ['mu10', 'v'], 4, 120000)]
results = []
for name, spec, rad, cap in PROBES:
    t0 = time.time()
    D = dirset(spec)
    P = ball(D, rad, cap)
    E, ncand = unit_edges(P)
    tri = triangles(len(P), E)
    col = three_colour(len(P), E)
    ok = col is not None and all(col[a] != col[b] for a, b in E.tolist())
    Dset = set(map(tuple, D))
    nongen = sum(1 for a, b in E.tolist() if tuple(P[a] - P[b]) not in Dset)
    res = dict(probe=name, directions=int(len(D)), radius=rad, vertices=int(len(P)), exact_unit_edges=int(len(E)),
               numeric_candidates=int(ncand), non_generator_unit_edges=int(nongen), triangles=int(tri),
               three_colourable=col is not None, colouring_verified=bool(ok), seconds=round(time.time() - t0, 1))
    results.append(res)
    print(json.dumps(res), flush=True)
    assert tri == 0

# ---------- why the probes are Cartesian (ranks), and how short the coupling relations really are ----------
from flint import fmpz_mat
rk = lambda rows: fmpz_mat([[int(x) for x in r] for r in rows]).rank()
Zrows = [to_vec(root_of_unity(k)) for k in range(55)]
Vrows = [to_vec(vconj[r]) for r in sorted(vconj)]
VZrows = [to_vec(mul(vconj[r], root_of_unity(k))) for r in sorted(vconj) for k in range(55)]
M5rows = [to_vec(root_of_unity(11 * k)) for k in range(5)]
VM5rows = [to_vec(mul(vconj[r], root_of_unity(11 * k))) for r in sorted(vconj) for k in range(5)]
ranks = dict(roots_of_unity=rk(Zrows), roots_plus_v_conjugates=rk(Zrows + Vrows),
             roots_plus_all_zeta_k_v_m=rk(Zrows + VZrows), mu5_plus_v_m_mu5=rk(M5rows + VM5rows))
print('Q-ranks:', ranks, flush=True)
coupling = []
one = root_of_unity(0)
for trip in [(1, 3, 4), (2, 9, 3), (1, 10, 2), (1, 2, 3)]:
    R = []
    for w in [one] + [vconj[r] for r in trip]:
        R += [to_vec(mul(w, root_of_unity(k))) for k in range(40)]
    nr, nc = len(R), len(R[0])
    SC = 10 ** 12
    Mk = fmpz_mat([[SC * int(x) for x in R[i]] + [1 if j == i else 0 for j in range(nr)] for i in range(nr)])
    B = Mk.lll()
    rels = []
    for r in range(B.nrows()):
        if all(int(B[r, c]) == 0 for c in range(nc)):
            vec = [int(B[r, nc + j]) for j in range(nr)]
            if any(vec[40:]):
                rels.append(vec)
    # verify exactly: sum_i vec_i * R_i = 0
    best = min(rels, key=lambda v: sum(abs(x) for x in v))
    assert all(sum(best[i] * int(R[i][c]) for i in range(nr)) == 0 for c in range(nc))
    coupling.append(dict(v_conjugates=list(trip), kernel_vectors=len(rels), shortest_L1=sum(abs(x) for x in best),
                         max_coefficient=max(abs(x) for x in best),
                         relation={'one': best[:40], 'v_%d' % trip[0]: best[40:80], 'v_%d' % trip[1]: best[80:120],
                                   'v_%d' % trip[2]: best[120:160]}))
    print('coupling relation', {k: v for k, v in coupling[-1].items() if k != 'relation'}, flush=True)
json.dump(dict(probes=results, ranks=ranks, coupling_relations_kernel_lll=coupling),
          open('05-knowledge/results/cubic_carrier_three_colour_probes_20261016.json', 'w'), indent=1)
print('wrote 05-knowledge/results/cubic_carrier_three_colour_probes_20261016.json')
