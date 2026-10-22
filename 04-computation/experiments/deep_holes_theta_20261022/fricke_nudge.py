"""The Fricke point is the deep hole nudged towards a second-shell Leech point (session opus-2026-10-22-S53; THM-4759 (10)).

For every Niemeier type, with c'' the (h+1)-point of leech_h_plus_one.py, the script checks:
  - z'' = z + r for a Leech root r (r^2 = 2, r.rho = -1) with r.z = -1, i.e. z'' is the reflection of the deep-hole
    cusp z in a wall of the Leech chamber;
  - c'' = c + (lambda - c)/(h+1), where lambda is the unique nearest Leech point to c'' and c = ((h+1)c'' - lambda)/h
    is a deep hole of the same type: 24 + k vertices at norm 2, and |lambda - c|^2 = 2 + 2/h (second shell);
  - exactly 24 of the 24 + k hole vertices stay at distance sqrt 2 from c''; the other k move away (the k affine
    nodes).
Usage: python fricke_nudge.py"""
import sys, time, math, re
import numpy as np
from cypari import pari
from numba import njit

ONLY = sys.argv[1:]
pari.allocatemem(2 * 10 ** 9)
FAIL = []


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


src = open(__file__.replace('fricke_nudge.py', 'all_holes.py'), encoding='utf-8').read()
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
    # wall: r = z'' - z
    r2 = [a - b for a, b in zip(zpp2, z2)]
    ok_wall = lor2(r2, r2) // 4 == 2 and lor2(r2, rho2) // 4 == -1 and lor2(r2, z2) // 4 == -1
    # unique nearest Leech point to c'' : vectors a with |a + Ym/m|^2 minimal -> lambda = -a
    svL, saL = shell_vecs(GL_f, GL_i, Ym, m, 2.0 + 1e-6)
    vals = [int(s) for s in svL]
    mn = min(vals)
    amin = [a for s, a in zip(svL, saL) if int(s) == mn]
    lam_c = -amin[0]                                    # lambda in Leech coordinates; c'' - lambda = (Ym + m a)/m
    # c = ((h+1) c'' - lambda)/h  -> h-scaled coordinates: Yc = Ym - lam_c  (since (h+1)c'' = Ym in units 1/1? careful)
    # (h+1) c'' has coordinates Ym (integers); so h c = Ym - lambda, i.e. coordinates of h c are Yc = Ym - lam_c
    Yc = Ym - lam_c
    svc, sac = shell_vecs(GL_f, GL_i, Yc, h, 2.0 + 2.0 / h + 1e-6)
    vc = [int(s) for s in svc]
    nv = sum(1 for s in vc if s == 2 * h * h)
    n2nd = sum(1 for s in vc if s == 2 * h * h + 2 * h)
    deep = min(vc) == 2 * h * h and nv == n0
    lam_dist = int((h * lam_c - Yc) @ GL_i @ (h * lam_c - Yc))          # h^2 |lambda - c|^2
    second = lam_dist == 2 * h * h + 2 * h
    # vertices of c (vectors of Lambda + c of norm 2): which stay at distance sqrt 2 from c''?
    stay = 0
    for s, a in zip(svc, sac):
        if int(s) != 2 * h * h:
            continue
        X = h * a + Yc                                   # h (vertex - ... ) : vector of Lambda + c scaled by h
        # point of Lambda + c'': same lattice translate a, i.e. a + Ym/m ; norm m^2 |a + Ym/m|^2
        Xpp = m * a + Ym
        if int(Xpp @ GL_i @ Xpp) == 2 * m * m:
            stay += 1
    k = n0 - 24
    check(ok_wall and len(amin) == 1 and deep and second and stay == 24,
          f'{name:9s} h={h:2d}: z\'\' = z + r (Leech wall, r.z = -1): {ok_wall}; c\'\' = c + (lambda - c)/(h+1) with lambda the unique nearest point, '
          f'c a deep hole ({nv} vertices): {deep}, |lambda - c|^2 = 2 + 2/h: {second}; vertices kept at sqrt 2: {stay} (of {n0}; k = {k} move away)')
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
