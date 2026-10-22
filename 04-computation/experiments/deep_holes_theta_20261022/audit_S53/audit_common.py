"""Shared helpers for the S53 audit (independent re-implementation where possible; the II_{25,1} set-up and
Borcherds' cusp table are read from ../all_holes.py, which is not modified).

Conventions: vectors of II_{25,1} are stored in doubled coordinates (26 integers, the last one timelike), so that
lor2(x, y) = 4 x.y. The Leech lattice is {rho, rho + r0}^perp with r0 = e_0 - e_1 a Leech root."""
import os, math
import numpy as np
from cypari import pari
from numba import njit

HERE = os.path.dirname(os.path.abspath(__file__))
PARENT = os.path.dirname(HERE)
n = 25
pari.allocatemem(2 * 10 ** 9)


def lor2(x, y):
    return int(sum(int(a) * int(b) for a, b in zip(x[:-1], y[:-1])) - int(x[-1]) * int(y[-1]))


src = open(os.path.join(PARENT, 'all_holes.py'), encoding='utf-8').read()
TABLE = src.split('TABLE = """\\\n')[1].split('"""')[0]
ROWS = []
for line in TABLE.strip('\n').split('\n'):
    t = line.split()
    co = [int(x) for x in t[6:]]
    ROWS.append((t[4], int(t[0]), int(t[1]), [co[1 + 24 - i] for i in range(n)] + [co[0]]))

gens = []
for i in range(n - 1):
    v = [0] * (n + 1); v[i] = 2; v[i + 1] = -2; gens.append(v)
v = [0] * (n + 1); v[n - 2] = 2; v[n - 1] = 2; gens.append(v)
gens.append([1] * (n + 1))
v = [0] * (n + 1); v[n - 1] = 2; v[n] = 2; gens.append(v)
v = [0] * (n + 1); v[0] = 4; gens.append(v)
_Am = pari.matrix(n + 1, len(gens), [gens[j][i] for i in range(n + 1) for j in range(len(gens))])
_Hb = pari.mathnf(_Am)
BASIS = [[int(_Hb[i, j]) for i in range(n + 1)] for j in range(n + 1)]
RHO2 = [2 * i for i in range(n)] + [140]
RHOP2 = list(RHO2); RHOP2[0] += 2; RHOP2[1] -= 2

assert lor2(RHO2, RHO2) == 0
# check unimodularity of the II basis
_G = pari.matrix(n + 1, n + 1, [lor2(a, b) // 4 for a in BASIS for b in BASIS])
assert abs(int(pari.matdet(_G))) == 1


def perp(forms, dim):
    F = pari.matrix(len(forms), n + 1, [lor2(b, f) // 4 for f in forms for b in BASIS])
    K = pari.matkerint(F)
    return [[sum(int(K[i, j]) * BASIS[i][t] for i in range(n + 1)) for t in range(n + 1)] for j in range(dim)]


def lll(vecs):
    d = len(vecs)
    G = pari.matrix(d, d, [lor2(a, b) // 4 for a in vecs for b in vecs])
    U = pari.qflllgram(G)
    return [[sum(int(U[i, j]) * vecs[i][t] for i in range(d)) for t in range(n + 1)] for j in range(d)]


LAM = lll(perp([RHO2, RHOP2], 24))
GLAM = [[lor2(a, b) // 4 for b in LAM] for a in LAM]
GLAM_I = np.array(GLAM, dtype=np.int64)
GLAM_F = GLAM_I.astype(np.float64)
GLAM_P = pari.matrix(24, 24, [x for row in GLAM for x in row])
assert int(pari.matdet(GLAM_P)) == 1
assert int(pari.qfminim(GLAM_P, 2, 0, 0)[0]) == 0        # rootless
assert int(pari.qfminim(GLAM_P, 4, 0, 0)[0]) == 196560    # Leech


def leech_coords(x2):
    """coordinates (rational, as a pari column) of the Leech projection of x2 in the LAM basis."""
    sol = pari.matsolve(GLAM_P, pari.matrix(24, 1, [lor2(a, x2) // 4 for a in LAM]))
    return sol


@njit(cache=True)
def shell_vecs(Gf, Gi, Yd, den, R2):
    """all integer a with |a + Yd/den|^2 <= R2 (Gram Gf); returns exact den^2 |a + Yd/den|^2 and the vectors a.
    Floating-point Fincke-Pohst with a 1e-9 safety margin; exact integer norms are recomputed."""
    d = Gf.shape[0]
    q = Gf.copy()
    for i in range(d):
        for j in range(i + 1, d):
            q[j, i] = q[i, j]
            q[i, j] = q[i, j] / q[i, i]
        for k in range(i + 1, d):
            for l in range(k, d):
                q[k, l] -= q[k, i] * q[i, l]
    y = -Yd.astype(np.float64) / den
    vals = []
    vecs = []
    a = np.zeros(d, dtype=np.int64)
    T = np.zeros(d + 1)
    cen = np.zeros(d)
    hi = np.zeros(d, dtype=np.int64)
    i = d - 1
    cen[i] = y[i]
    rad = np.sqrt(R2 + 1e-9) / np.sqrt(q[i, i])
    a[i] = np.int64(np.ceil(cen[i] - rad))
    hi[i] = np.int64(np.floor(cen[i] + rad))
    while True:
        if a[i] > hi[i]:
            i += 1
            if i == d:
                break
            a[i] += 1
            continue
        t = a[i] - cen[i]
        T[i] = T[i + 1] + q[i, i] * t * t
        if T[i] > R2 + 1e-9:
            a[i] += 1
            continue
        if i == 0:
            s = 0
            for u in range(d):
                zu = den * a[u] + Yd[u]
                if zu == 0:
                    continue
                acc = 0
                for w in range(d):
                    acc += Gi[u, w] * (den * a[w] + Yd[w])
                s += zu * acc
            vals.append(s)
            vecs.append(a.copy())
            a[i] += 1
            continue
        j = i - 1
        sacc = 0.0
        for l in range(j + 1, d):
            sacc += q[j, l] * (a[l] - y[l])
        cen[j] = y[j] - sacc
        rad = np.sqrt(max(R2 + 1e-9 - T[i], 0.0) / q[j, j])
        a[j] = np.int64(np.ceil(cen[j] - rad))
        hi[j] = np.int64(np.floor(cen[j] + rad))
        T[j + 1] = T[i]
        i = j
    return vals, vecs


def bez(vals):
    g, coef = 0, [0] * len(vals)
    for j, val in enumerate(vals):
        if val == 0:
            continue
        if g == 0:
            g, coef = val, [0] * len(vals); coef[j] = 1; continue
        a0, b0, a1, b1, r0_, r1_ = 1, 0, 0, 1, g, val
        while r1_:
            qq = r0_ // r1_
            r0_, r1_ = r1_, r0_ - qq * r1_
            a0, a1 = a1, a0 - qq * a1
            b0, b1 = b1, b0 - qq * b1
        coef = [a0 * c for c in coef]; coef[j] += b0; g = r0_
    return [-c * g for c in coef]


def h_plus_one_point(z2, h):
    """The (h+1)-point c'' of a Niemeier cusp z (Borcherds' table): returns (Ym, m) with m c'' = sum Ym_i LAM_i.
    Construction as in fricke_h_plus_one.py: z' with z.z' = -1, the N-basis of {z, z'}^perp, the projected Weyl
    vector, the closest point of N + rho_N/h, and z'' = z' + x + (x^2/2) z."""
    zv = [lor2(b, z2) // 4 for b in BASIS]
    cf = bez(zv)
    y2 = [sum(cf[j] * BASIS[j][t] for j in range(n + 1)) for t in range(n + 1)]
    yy = lor2(y2, y2) // 4
    zp2 = [a + (yy // 2) * b for a, b in zip(y2, z2)]
    assert lor2(zp2, zp2) == 0 and lor2(zp2, z2) // 4 == -1
    Nv = lll(perp([z2, zp2], 24))
    GN = [[lor2(a, b) // 4 for b in Nv] for a in Nv]
    zrho = lor2(z2, RHO2) // 4
    rzp = lor2(RHO2, zp2) // 4
    rhoN2 = [r + rzp * a + zrho * b for r, a, b in zip(RHO2, z2, zp2)]
    GNp = pari.matrix(24, 24, [v for row in GN for v in row])
    sol = pari.matsolve(GNp, pari.matrix(24, 1, [lor2(a, rhoN2) // 4 for a in Nv]))
    Yr = np.array([int(sol[i, 0]) for i in range(24)], dtype=np.int64)
    Gi = np.array(GN, dtype=np.int64)
    sv, sa = shell_vecs(Gi.astype(np.float64), Gi, -Yr, h, 2.0 * (h + 1) / h + 1e-6)
    best = min(int(s) for s in sv)
    a0 = [a for s, a in zip(sv, sa) if int(s) == best][0]
    x2 = [sum(int(a0[i]) * Nv[i][t] for i in range(24)) for t in range(n + 1)]
    xx = lor2(x2, x2) // 4
    zpp2 = [a + b + (xx // 2) * c for a, b, c in zip(zp2, x2, z2)]
    assert lor2(zpp2, zpp2) == 0 and lor2(zpp2, z2) // 4 == -1
    m = -lor2(zpp2, RHO2) // 4
    zr = lor2(zpp2, RHO2) // 4
    zrp_ = lor2(zpp2, RHOP2) // 4
    lz2 = [a + zrp_ * b + zr * c for a, b, c in zip(zpp2, RHO2, RHOP2)]
    sol = leech_coords(lz2)
    Ym = np.array([int(sol[i, 0]) for i in range(24)], dtype=np.int64)
    return Ym, m


def deep_hole_point(z2, h):
    """The deep hole c of type z: h c = proj_Lambda(z); returns (Yh, h') with h' = -z.rho."""
    m = -lor2(z2, RHO2) // 4
    zr = lor2(z2, RHO2) // 4
    zrp_ = lor2(z2, RHOP2) // 4
    lz2 = [a + zrp_ * b + zr * c for a, b, c in zip(z2, RHO2, RHOP2)]
    sol = leech_coords(lz2)
    return np.array([int(sol[i, 0]) for i in range(24)], dtype=np.int64), m


def coset_profile(Y, den, R2=2.0):
    """norm profile {norm*den^2: count} of the coset Lambda + (Y/den) (sign irrelevant) up to norm R2."""
    vals, _ = shell_vecs(GLAM_F, GLAM_I, np.ascontiguousarray(Y, dtype=np.int64), den, R2 + 1e-6)
    prof = {}
    for s in vals:
        prof[int(s)] = prof.get(int(s), 0) + 1
    return prof
