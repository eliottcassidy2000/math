"""Swap reciprocity of the coset theta series of M = {rho, z}^perp (session opus-2026-10-22-S53; HYP-9363 structure).

Notation.
- rho is the Leech cusp of II_{25,1}. z is the cusp of a Niemeier lattice N, with z . rho = -h.
- M = {rho, z}^perp is positive definite of rank 24 with discriminant group D = (Z/h)^2. D is hyperbolic: it is the
  direct sum of two isotropic lines, H_L = Lambda_M/M and H_N = N_M/M.
- Write gamma = alpha g_L + beta g_N with alpha, beta in Z/h. The swap sigma(alpha, beta) = (beta, alpha) is an
  isometry of D, but it is not induced by any isometry of M (it would force theta_Lambda = theta_N).
- E_gamma = theta_(M + gamma) - theta_(M + sigma gamma).

Root-count form. For a root r of II_{25,1} with (r . rho, r . z) = (-a, -b), the projection of r to M has norm
2 + 2ab/h and lies in the class (a, b) mod h, up to the fixed identifications. Swap reciprocity is the statement that
the number of roots with (-r.rho, -r.z) = (a, b) equals the number with (b, a) whenever h divides neither a nor b.

Hypothesis tested here (it implies HYP-9363):
- E vanishes off the axes {alpha = 0} u {beta = 0};
- on an axis, E_(alpha, 0) is a multiple of Delta that depends only on the order of alpha.

Method: enumerate every vector of M^* with |x|^2 <= B, classify it by its class in D, and compare the coefficient
tables of theta_gamma and theta_(sigma gamma).
Needs cypari and numba. Usage: python swap_reciprocity.py B type [type ...]"""
import sys, time, math
import numpy as np
from fractions import Fraction
from cypari import pari
from numba import njit

B = float(sys.argv[1]) if len(sys.argv) > 1 else 4.0
ONLY = sys.argv[2:] if len(sys.argv) > 2 else ['a24', 'd24', 'a2^12', 'a8^3', 'd4^6', 'e8^3']
pari.allocatemem(2 * 10 ** 9)
T0 = time.time()
FAIL = []


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


# the cusp table is read from all_holes.py (same directory)
src = open(__file__.replace('swap_reciprocity.py', 'all_holes.py'), encoding='utf-8').read()
TABLE = src.split('TABLE = """\\\n')[1].split('"""')[0]
n = 25


def lor2(x, y):
    return int(sum(int(a) * int(b) for a, b in zip(x[:-1], y[:-1])) - int(x[-1]) * int(y[-1]))


rows = {}
for line in TABLE.strip('\n').split('\n'):
    t = line.split()
    co = [int(x) for x in t[6:]]
    rows[t[4]] = (int(t[0]), int(t[1]), [co[1 + 24 - i] for i in range(n)] + [co[0]])

gens = []
for i in range(n - 1):
    v = [0] * (n + 1); v[i] = 2; v[i + 1] = -2; gens.append(v)
v = [0] * (n + 1); v[n - 2] = 2; v[n - 1] = 2; gens.append(v)
gens.append([1] * (n + 1))
v = [0] * (n + 1); v[n - 1] = 2; v[n] = 2; gens.append(v)
v = [0] * (n + 1); v[0] = 4; gens.append(v)
A = pari.matrix(n + 1, len(gens), [gens[j][i] for i in range(n + 1) for j in range(len(gens))])
Hb = pari.mathnf(A)
basis = [[int(Hb[i, j]) for i in range(n + 1)] for j in range(n + 1)]
rho2 = [2 * i for i in range(n)] + [140]


def perp_basis(u2, w2):
    F = pari.matrix(2, n + 1, [lor2(b, u2) // 4 for b in basis] + [lor2(b, w2) // 4 for b in basis])
    K = pari.matkerint(F)
    vecs = [[sum(int(K[i, j]) * basis[i][t] for i in range(n + 1)) for t in range(n + 1)] for j in range(24)]
    G = pari.matrix(24, 24, [lor2(a, b) // 4 for a in vecs for b in vecs])
    U = pari.qflllgram(G)
    vecs = [[sum(int(U[i, j]) * vecs[i][t] for i in range(24)) for t in range(n + 1)] for j in range(24)]
    return vecs, [[lor2(a, b) // 4 for b in vecs] for a in vecs]


def kernel_perp1(u2):
    F = pari.matrix(1, n + 1, [lor2(b, u2) // 4 for b in basis])
    K = pari.matkerint(F)
    return [[sum(int(K[i, j]) * basis[i][t] for i in range(n + 1)) for t in range(n + 1)] for j in range(25)]


@njit(cache=True)
def enum_dual(Q, Bd, Bnd):
    """All integer vectors a with a^T Q a <= Bnd (Q = Gram of M^* scaled by Bd^2 to be integral up to Bd^2), half-space."""
    d = Q.shape[0]
    q = Q.astype(np.float64).copy()
    for i in range(d):
        for j in range(i + 1, d):
            q[j, i] = q[i, j]
            q[i, j] = q[i, j] / q[i, i]
        for k in range(i + 1, d):
            for l in range(k, d):
                q[k, l] -= q[k, i] * q[i, l]
    out = []
    a = np.zeros(d, dtype=np.int64)
    T = np.zeros(d + 1)
    cen = np.zeros(d)
    hi = np.zeros(d, dtype=np.int64)
    i = d - 1
    rad = np.sqrt(Bnd / q[i, i])
    a[i] = np.int64(np.ceil(-rad))
    hi[i] = np.int64(np.floor(rad))
    while True:
        if a[i] > hi[i]:
            i += 1
            if i == d:
                break
            a[i] += 1
            continue
        t = a[i] - cen[i]
        T[i] = T[i + 1] + q[i, i] * t * t
        if T[i] > Bnd + 1e-7:
            a[i] += 1
            continue
        if i == 0:
            out.append(a.copy())
            a[i] += 1
            continue
        j = i - 1
        s = 0.0
        for l in range(j + 1, d):
            s += q[j, l] * a[l]
        cen[j] = -s
        rad = np.sqrt(max(Bnd + 1e-7 - T[i], 0.0) / q[j, j])
        a[j] = np.int64(np.ceil(cen[j] - rad))
        hi[j] = np.int64(np.floor(cen[j] + rad))
        T[j + 1] = T[i]
        i = j
    return out


for name in ONLY:
    t1 = time.time()
    h, n0, z2 = rows[name]
    Mv, GM = perp_basis(rho2, z2)
    GMp = pari.matrix(24, 24, [v for row in GM for v in row])
    detM = int(pari.matdet(GMp))
    Ginv = pari.matsolve(GMp, pari.matid(24))
    Gi = [[Fraction(str(Ginv[i, j])) for j in range(24)] for i in range(24)]
    # dual basis vectors (coordinates in Mv): columns of G^-1; Gram of dual = G^-1
    den = 1
    for row in Gi:
        for x in row:
            den = den * x.denominator // math.gcd(den, x.denominator)
    QD = np.array([[int(x * den) for x in row] for row in Gi], dtype=np.int64)   # den * Gram(M^*)
    # lifts: Lambda_M = proj of rho^perp, N_M = proj of z^perp (projection along span(rho, z))
    def proj_M(x2):         # x2 doubled; returns coordinates (Fractions) of the M-component in basis Mv
        g = [Fraction(lor2(a, x2), 4) for a in Mv]
        return [sum(Gi[i][j] * g[j] for j in range(24)) for i in range(24)]
    rp = kernel_perp1(rho2)
    zp = kernel_perp1(z2)
    def order(c):
        o = 1
        for t in c:
            o = o * t.denominator // math.gcd(o, t.denominator)
        return o
    gL = max((proj_M(x) for x in rp), key=order)
    gN = max((proj_M(x) for x in zp), key=order)
    check(order(gL) == h and order(gN) == h, f'{name}: generators of H_L, H_N have order h')
    # pairing of a dual vector (coords y in dual basis -> M-basis coords Gi y) with gL, gN:  (Gi y) . G . gX = y . gX
    def cls(y):
        al = sum(Fraction(int(y[i])) * gN[i] for i in range(24))      # y in dual-coordinates: x = sum y_i e_i^*, x . v = sum y_i v_i
        be = sum(Fraction(int(y[i])) * gL[i] for i in range(24))
        return (int((al * h) % h), int((be * h) % h))
    Bnd = B * den
    vecs = enum_dual(QD, den, Bnd)
    table = {}
    for y in vecs:
        nrm = Fraction(int(y @ QD @ y), int(den))          # |x|^2 with x = sum y_i e_i^*
        c = cls(y)
        table.setdefault(c, {}).setdefault(nrm, 0)
        table[c][nrm] += 1
    # reciprocity test
    ok_off = True
    bad = []
    axis = {}
    for al in range(h):
        for be in range(h):
            if al == 0 or be == 0:
                continue
            t1_ = table.get((al, be), {})
            t2_ = table.get((be, al), {})
            if t1_ != t2_:
                ok_off = False
                bad.append(((al, be), sorted(t1_.items())[:3], sorted(t2_.items())[:3]))
    for al in range(1, h):
        ta, tb = table.get((al, 0), {}), table.get((0, al), {})
        keys = sorted(set(ta) | set(tb))
        axis[al] = [(k, ta.get(k, 0) - tb.get(k, 0)) for k in keys if ta.get(k, 0) != tb.get(k, 0)]
    print(f'== {name} (h = {h}, 24+k = {n0}): det M = {detM}, |M^* vectors with |x|^2 <= {B}| = {len(vecs)}, classes seen = {len(table)} [{time.time() - t1:.0f}s]')
    check(detM == h * h, f'{name}: det M = h^2')
    check(ok_off, f'{name}: theta_(alpha,beta) = theta_(beta,alpha) for all {(h - 1) ** 2} off-axis classes, |x|^2 <= {B}' + ('' if ok_off else f'; first failures {bad[:2]}'))
    from math import gcd
    byorder = {}
    for al, diffs in axis.items():
        d = h // gcd(al, h)
        byorder.setdefault(d, set()).add(tuple(diffs))
    print('     axis differences theta_(a,0) - theta_(0,a) by order of a (norm: diff):', {d: sorted(v)[:1] for d, v in sorted(byorder.items())})
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f' ({time.time() - T0:.0f}s)')
