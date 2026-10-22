"""The point rho_N/(h+1) of a Niemeier lattice and the inverse of Fricke descent (session opus-2026-10-22-S53; HYP-9368).

Venkov's |rho_N|^2 = 2h(h+1) gives |rho_N/(h+1)|^2 = 2 - 2/(h+1), which is exactly the Fricke norm 2 - 2/p for
p = h + 1. For every Niemeier lattice N the script:
  1. normalises rho_N: the closest-vector problem for rho/h selects a true Weyl vector of norm 2h(h+1);
  2. finds the minimal norm of the coset N + rho_N/(h+1) and how often it is attained;
  3. builds the (h+1)-neighbour N'' = N_0 + Z w (N_0 = {x : x.rho/(h+1) in Z}, w in rho/(h+1) + N of even norm) and
     counts its roots.
When N'' is rootless it is the Leech lattice. Then rho_N/(h+1) is the image, under the inverse orbifold, of a Fricke
point of Lambda of order h + 1.
Result: this happens for ALL 23 Niemeier lattices (a second holy construction, THM-4759 (10)). In every case the
minimum is attained once, at rho_N/(h+1) itself. An earlier version of this script predicted that only the ten
Fricke descents (A1^24, A2^12, A3^8, A4^6, D4^6, D6^4, A7^2D5^2, A8^3, D8^3, D12^2) would qualify; that prediction was
refuted.
Needs cypari and numba. Usage: python weyl_h_plus_one.py [types...]"""
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


src = open(__file__.replace('weyl_h_plus_one.py', 'all_holes.py'), encoding='utf-8').read()
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
    Yr = np.array([int(sol[i, 0]) for i in range(24)], dtype=np.int64)      # rho_N + h x for some x in N
    Gi = np.array(GN, dtype=np.int64)
    Gf = Gi.astype(np.float64)
    # 1. true Weyl vector: rho = Yr - h y with y a closest point of N to Yr/h (norm 2h(h+1))
    sv, sa = shell_vecs(Gf, Gi, -Yr, h, 2.0 * (h + 1) / h + 1e-6)    # vectors a with |a - Yr/h|^2 small
    best = min(int(s) for s in sv)
    a0 = [a for s, a in zip(sv, sa) if int(s) == best][0]
    rho = Yr - h * a0                    # h (a0 - Yr/h) = h a0 - Yr; rho = Yr - h a0
    ok_rho = int(rho @ Gi @ rho) == 2 * h * (h + 1)
    # 2. the coset N + rho/(h+1)
    p = h + 1
    sv, sa = shell_vecs(Gf, Gi, rho, p, 2.0 + 1e-6)
    vals = np.array([int(s) for s in sv], dtype=np.int64)
    mn = int(vals.min())
    nmin = int((vals == mn).sum())
    fricke_norm = mn == 2 * p * p - 2 * p                  # p^2 (2 - 2/p)
    # 3. the p-neighbour: w = rho/p + lam with |w|^2 even; roots = norm-2 x in N + j rho/p with (x - j w).(rho/p) in Z
    W = None
    for lam_c in [np.zeros(24, dtype=np.int64)] + [np.eye(24, dtype=np.int64)[i] * s for i in range(24) for s in (1, -1)]:
        Wv = rho + p * lam_c
        nw = int(Wv @ Gi @ Wv)
        if nw % (p * p) == 0 and (nw // (p * p)) % 2 == 0:
            W = Wv; break
    if W is None:
        for i in range(24):
            for k in range(24):
                lam_c = np.eye(24, dtype=np.int64)[i] + np.eye(24, dtype=np.int64)[k]
                Wv = rho + p * lam_c
                nw = int(Wv @ Gi @ Wv)
                if nw % (p * p) == 0 and (nw // (p * p)) % 2 == 0:
                    W = Wv; break
            if W is not None:
                break
    nroots = 0
    for j in range(1, p):
        sv, sa = shell_vecs(Gf, Gi, j * rho, p, 2.0 + 1e-6)
        for s, a in zip(sv, sa):
            if int(s) != 2 * p * p:
                continue
            X = p * a + j * rho                             # p x
            if int((X - j * W) @ Gi @ rho) % (p * p) == 0:
                nroots += 1
    leech = nroots == 0
    check(ok_rho and leech and fricke_norm and nmin == 1,
          f'{name:9s} h={h:2d}: N + rho/{p}: min {mn}/{p * p} (x{nmin}){" = 2 - 2/" + str(p) if fricke_norm else ""}; '
          f'{p}-neighbour has {nroots} roots{" -> LEECH (inverse Fricke descent)" if leech else ""} '
          f'[{time.time() - t1:.0f}s]')
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
