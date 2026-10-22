"""Fricke twining at the (h+1)-points of the Leech lattice (session opus-2026-10-22-S53; THM-4759 (10)).

For every Niemeier lattice N let z'' be the cusp of type N with z.z'' = -1 and z''.rho = -(h+1), and let
c'' = proj_Lambda(z'')/(h+1), a point of order e = h + 1. The script checks:
  - Lambda + c'' has a unique vector of norm 2 - 2/e and exactly 24 vectors of norm 2;
  - Fricke invariance of Z = Theta_Lambda^[c'']/eta^24, i.e. Theta^[c'']/Delta = theta_(Lambda+c'')(e tau)/Delta(e tau),
    through q^(NORM/2 - 1);
  - the proper multiples j c'' (j | e, 1 < j < e) are deep holes of norm 2, so they give no further poles.
With the principal-part argument of THM-4759 (7), Z - 24 is then the Gamma_0(e)+e Hauptmodul.
Needs cypari and numba. Usage: python fricke_h_plus_one.py [NORM=6] [types...]"""
import sys, time, math, re
import numpy as np
from cypari import pari
from numba import njit

NORM = int(sys.argv[1]) if len(sys.argv) > 1 else 6
ONLY = sys.argv[2:]
pari.allocatemem(2 * 10 ** 9)
FAIL = []


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


src = open(__file__.replace('fricke_h_plus_one.py', 'all_holes.py'), encoding='utf-8').read()
TABLE = src.split('TABLE = """\\\n')[1].split('"""')[0]
n = 25


def lor2(x, y):
    return int(sum(int(a) * int(b) for a, b in zip(x[:-1], y[:-1])) - int(x[-1]) * int(y[-1]))


rows = []
for line in TABLE.strip('\n').split('\n'):
    t = line.split()
    co = [int(x) for x in t[6:]]
    rows.append((t[4], int(t[0]), int(t[1]), [co[1 + 24 - i] for i in range(n)] + [co[0]]))
gens = []
for i in range(n - 1):
    v = [0] * (n + 1); v[i] = 2; v[i + 1] = -2; gens.append(v)
v = [0] * (n + 1); v[n - 2] = 2; v[n - 1] = 2; gens.append(v)
gens.append([1] * (n + 1))
v = [0] * (n + 1); v[n - 1] = 2; v[n] = 2; gens.append(v)
v = [0] * (n + 1); v[0] = 4; gens.append(v)
Am = pari.matrix(n + 1, len(gens), [gens[j][i] for i in range(n + 1) for j in range(len(gens))])
Hb = pari.mathnf(Am)
basis = [[int(Hb[i, j]) for i in range(n + 1)] for j in range(n + 1)]
rho2 = [2 * i for i in range(n)] + [140]


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


def perp_basis(u2, w2):
    F = pari.matrix(2, n + 1, [lor2(b, u2) // 4 for b in basis] + [lor2(b, w2) // 4 for b in basis])
    K = pari.matkerint(F)
    vecs = [[sum(int(K[i, j]) * basis[i][t] for i in range(n + 1)) for t in range(n + 1)] for j in range(24)]
    G = pari.matrix(24, 24, [lor2(a, b) // 4 for a in vecs for b in vecs])
    U = pari.qflllgram(G)
    vecs = [[sum(int(U[i, j]) * vecs[i][t] for i in range(24)) for t in range(n + 1)] for j in range(24)]
    return vecs, [[lor2(a, b) // 4 for b in vecs] for a in vecs]


@njit(cache=True)
def shell_vecs(Gf, Gi, Yd, den, R2):
    """all a with |a + Yd/den|^2 <= R2: returns den^2 |a + Yd/den|^2 and a."""
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



rhop2 = list(rho2); rhop2[0] += 2; rhop2[1] -= 2
F = pari.matrix(2, n + 1, [lor2(b, rho2) // 4 for b in basis] + [lor2(b, rhop2) // 4 for b in basis])
K = pari.matkerint(F)
lam = [[sum(int(K[i, j]) * basis[i][t] for i in range(n + 1)) for t in range(n + 1)] for j in range(24)]
GL = pari.matrix(24, 24, [lor2(a, b) // 4 for a in lam for b in lam])
U = pari.qflllgram(GL)
lam = [[sum(int(U[i, j]) * lam[i][t] for i in range(24)) for t in range(n + 1)] for j in range(24)]
GmL = [[lor2(a, b) // 4 for b in lam] for a in lam]
GpL = pari.matrix(24, 24, [x for row in GmL for x in row])
GL_i = np.array(GmL, dtype=np.int64)
GL_f = GL_i.astype(np.float64)
import cmath


@njit(cache=True)
def twisted_counts(Gf, Gi, GY, h, Bnd):
    d = Gf.shape[0]
    q = Gf.copy()
    for i in range(d):
        for j in range(i + 1, d):
            q[j, i] = q[i, j]
            q[i, j] = q[i, j] / q[i, i]
        for k in range(i + 1, d):
            for l in range(k, d):
                q[k, l] -= q[k, i] * q[i, l]
    M = Bnd // 2
    out = np.zeros((M + 1, h), dtype=np.int64)
    a = np.zeros(d, dtype=np.int64)
    T = np.zeros(d + 1)
    cen = np.zeros(d)
    hi = np.zeros(d, dtype=np.int64)
    i = d - 1
    rad = np.sqrt(Bnd / q[i, i])
    a[i] = np.int64(np.ceil(-rad - 1e-9))
    hi[i] = np.int64(np.floor(rad + 1e-9))
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
            s = 0
            for u in range(d):
                if a[u] == 0:
                    continue
                acc = 0
                for w in range(d):
                    acc += Gi[u, w] * a[w]
                s += a[u] * acc
            r = 0
            for u in range(d):
                r += a[u] * GY[u]
            out[s // 2, r % h] += 1
            a[i] += 1
            continue
        j = i - 1
        sacc = 0.0
        for l in range(j + 1, d):
            sacc += q[j, l] * a[l]
        cen[j] = -sacc
        rad = np.sqrt(max(Bnd + 1e-7 - T[i], 0.0) / q[j, j])
        a[j] = np.int64(np.ceil(cen[j] - rad))
        hi[j] = np.int64(np.floor(cen[j] + rad))
        T[j + 1] = T[i]
        i = j
    return out


def inv_delta(L, step=1):
    s = [0] * L; s[0] = 1
    for mm in range(1, L // step + 2):
        for _ in range(24):
            for t in range(step * mm, L):
                s[t] += s[t - step * mm]
    return s

for name, h, n0, z2 in rows:
    if ONLY and name not in ONLY:
        continue
    t1 = time.time()
    zv = [lor2(b, z2) // 4 for b in basis]
    cf = bez(zv)
    y2 = [sum(cf[j] * basis[j][t] for j in range(n + 1)) for t in range(n + 1)]
    yy = lor2(y2, y2) // 4
    zp2 = [a + (yy // 2) * b for a, b in zip(y2, z2)]
    Nv, GN = perp_basis(z2, zp2)
    zrho = lor2(z2, rho2) // 4
    rzp = lor2(rho2, zp2) // 4
    rhoN2 = [r + rzp * a + zrho * b for r, a, b in zip(rho2, z2, zp2)]
    GNp = pari.matrix(24, 24, [v for row in GN for v in row])
    sol = pari.matsolve(GNp, pari.matrix(24, 1, [lor2(a, rhoN2) // 4 for a in Nv]))
    Yr = np.array([int(sol[i, 0]) for i in range(24)], dtype=np.int64)
    Gi = np.array(GN, dtype=np.int64)
    Gf = Gi.astype(np.float64)
    sv, sa = shell_vecs(Gf, Gi, -Yr, h, 2.0 * (h + 1) / h + 1e-6)
    best = min(int(s) for s in sv)
    a0 = [a for s, a in zip(sv, sa) if int(s) == best][0]
    # z'' = z' + x + (x^2/2) z with x = a0 (N-coordinates) makes the projected Weyl vector the true one
    x2 = [sum(int(a0[i]) * Nv[i][t] for i in range(24)) for t in range(n + 1)]
    xx = lor2(x2, x2) // 4
    zpp2 = [a + b + (xx // 2) * c for a, b, c in zip(zp2, x2, z2)]
    ok_null = lor2(zpp2, zpp2) == 0 and lor2(zpp2, z2) // 4 == -1
    m = -lor2(zpp2, rho2) // 4
    # Leech-side point c'' = proj_Lambda(z'')/m
    zr = lor2(zpp2, rho2) // 4
    zrp_ = lor2(zpp2, rhop2) // 4
    lz2 = [a + zrp_ * b + zr * c for a, b, c in zip(zpp2, rho2, rhop2)]          # m c'' (doubled)
    sol = pari.matsolve(GpL, pari.matrix(24, 1, [lor2(a, lz2) // 4 for a in lam]))
    Ym = np.array([int(sol[i, 0]) for i in range(24)], dtype=np.int64)
    svL, saL = shell_vecs(GL_f, GL_i, Ym, m, 2.0 + 1e-6)
    vals = np.array([int(s) for s in svL], dtype=np.int64)
    mn = int(vals.min()) if len(vals) else None
    nmin = int((vals == mn).sum()) if len(vals) else 0
    n2 = int((vals == 2 * m * m).sum()) if len(vals) else 0
    fr = mn == 2 * m * m - 2 * m
    e = m
    cnt = twisted_counts(GL_f, GL_i, (GL_i @ Ym).astype(np.int64), e, NORM)
    M = NORM // 2
    th = [sum(int(cnt[k][r]) * cmath.exp(2j * math.pi * r / e) for r in range(e)) for k in range(M + 1)]
    ok_real = all(abs(zz.imag) < 1e-6 and abs(zz.real - round(zz.real)) < 1e-6 for zz in th)
    th = [round(zz.real) for zz in th]
    idl = inv_delta(M + 1)
    left = [sum(th[i] * idl[k - i] for i in range(0, k + 1)) for k in range(M + 1)]
    svR, saR = shell_vecs(GL_f, GL_i, Ym, e, 2.0 * (M + e) / e + 1e-6)
    expo = {}
    okint = True
    for s in svR:
        num = e * int(s)
        if num % (2 * e * e):
            okint = False
            continue
        kk = num // (2 * e * e)
        expo[kk] = expo.get(kk, 0) + 1
    idle = inv_delta(e * (M + 2), e)
    right = []
    for k in range(M + 1):
        tot = 0
        for ex, cn in expo.items():
            jj = (-1 + k) - (ex - e)
            if 0 <= jj < len(idle):
                tot += cn * idle[jj]
        right.append(tot)
    mult_ok = True
    for jm in [jm for jm in range(2, e) if e % jm == 0]:
        svj, saj = shell_vecs(GL_f, GL_i, jm * Ym, e, 2.0 + 1e-6)
        vj = [int(s) for s in svj]
        if not vj or min(vj) != 2 * e * e:
            mult_ok = False
    ok = ok_null and m == h + 1 and fr and nmin == 1 and n2 == 24 and ok_real and okint and left == right and mult_ok
    series = ' + '.join(f'{c}q^{k - 1}' for k, c in enumerate(left))
    check(ok, f'{name:9s} e = h+1 = {e}: unique norm 2-2/e vector, 24 of norm 2; proper multiples deep: {mult_ok}; '
              f'Fricke W_{e}: {left == right}; Z = {series} [{time.time() - t1:.0f}s]')
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
