"""Multiples of Leech deep holes, enumerated on the Leech side (session opus-2026-10-22-S53; THM-4759 (5)).

For each deep hole c of type X (Coxeter number h, so c has order h modulo Lambda) and each divisor d of h with d < h,
enumerate Lambda + d c up to norm 2 exactly (integer arithmetic, float pruning with slack). Prediction:
  - if d is a bad prime of every component of X and gcd(d, h/d) = 1: the minimal norm is 2 - 2d/h, attained once
    (d c is an interior point of a Voronoi cell, not a hole);
  - otherwise the minimal norm is 2, so d c is again a deep hole; its order is e = h/d, so its type has Coxeter number
    e, and its number of vertices 24 + k' identifies the type.
Independent of the Niemeier-side enumeration in multiple_cosets.py (the two agree below norm 2 by THM-4758).
Needs cypari and numba. Usage: python deep_hole_multiples.py [types...]"""
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


src = open(__file__.replace('deep_hole_multiples.py', 'all_holes.py'), encoding='utf-8').read()
TABLE = src.split('TABLE = """\\\n')[1].split('"""')[0]
n = 25


def lor2(x, y):
    return int(sum(int(a) * int(b) for a, b in zip(x[:-1], y[:-1])) - int(x[-1]) * int(y[-1]))


rows = []
for line in TABLE.strip('\n').split('\n'):
    t = line.split()
    co = [int(x) for x in t[6:]]
    rows.append((t[4], int(t[0]), int(t[1]), [co[1 + 24 - i] for i in range(n)] + [co[0]]))
TYPES = {}
for name, h, n0, z2 in rows:
    TYPES.setdefault((h, n0), []).append(name)          # (Coxeter number, number of hole vertices) -> type
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
rhop2 = list(rho2); rhop2[0] += 2; rhop2[1] -= 2
F = pari.matrix(2, n + 1, [lor2(b, rho2) // 4 for b in basis] + [lor2(b, rhop2) // 4 for b in basis])
K = pari.matkerint(F)
lam = [[sum(int(K[i, j]) * basis[i][t] for i in range(n + 1)) for t in range(n + 1)] for j in range(24)]
GL = pari.matrix(24, 24, [lor2(a, b) // 4 for a in lam for b in lam])
U = pari.qflllgram(GL)
lam = [[sum(int(U[i, j]) * lam[i][t] for i in range(24)) for t in range(n + 1)] for j in range(24)]
Gm = [[lor2(a, b) // 4 for b in lam] for a in lam]
Gp = pari.matrix(24, 24, [x for row in Gm for x in row])
Gi = np.array(Gm, dtype=np.int64)
Gf = np.array(Gm, dtype=np.float64)
assert int(pari.matdet(Gp)) == 1 and int(pari.qfminim(Gp, 2, 0, 0)[0]) == 0      # unimodular, rootless: Leech


@njit(cache=True)
def shell(Gf, Gi, Yd, h, R2):
    """values h^2 |a + Yd/h|^2 for lattice vectors a with |a + Yd/h|^2 <= R2."""
    d = Gf.shape[0]
    q = Gf.copy()
    for i in range(d):
        for j in range(i + 1, d):
            q[j, i] = q[i, j]
            q[i, j] = q[i, j] / q[i, i]
        for k in range(i + 1, d):
            for l in range(k, d):
                q[k, l] -= q[k, i] * q[i, l]
    y = -Yd.astype(np.float64) / h
    vals = []
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
                zu = h * a[u] + Yd[u]
                if zu == 0:
                    continue
                acc = 0
                for w in range(d):
                    acc += Gi[u, w] * (h * a[w] + Yd[w])
                s += zu * acc
            vals.append(s)
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
    return vals


BADP = {'a': lambda r: [], 'd': lambda r: [2], 'e': lambda r: [2, 3] if r in (6, 7) else [2, 3, 5]}
for name, h, n0, z2 in rows:
    if ONLY and name not in ONLY:
        continue
    t1 = time.time()
    comps = []
    for t, rr, ex in re.findall(r'([ade])(\d+)(?:\^(\d+))?', name):
        comps += [(t, int(rr))] * (int(ex) if ex else 1)
    zrho = lor2(z2, rho2) // 4
    zrp = lor2(z2, rhop2) // 4
    lz2 = [a + zrp * b + zrho * c for a, b, c in zip(z2, rho2, rhop2)]      # h c (doubled coordinates)
    sol = pari.matsolve(Gp, pari.matrix(24, 1, [lor2(a, lz2) // 4 for a in lam]))
    Yh = np.array([int(sol[i, 0]) for i in range(24)], dtype=np.int64)        # coordinates of h c
    allok = True
    out = []
    for d in [d for d in range(1, h) if h % d == 0]:
        vals = np.array(shell(Gf, Gi, d * Yh, h, 2.0 + 1e-6), dtype=np.int64)
        H2 = h * h
        mn = int(vals.min())
        nmin = int((vals == mn).sum())
        shallow = all(d in BADP[t](rr) for t, rr in comps) and math.gcd(d, h // d) == 1
        if d == 1:
            ok = mn == 2 * H2 and nmin == n0
            desc = f'c: deep hole, {nmin} vertices'
        elif shallow:
            ok = mn == 2 * H2 - 2 * d * h and nmin == 1
            desc = f'{d}c: interior point, min {mn}/{H2} = 2 - 2*{d}/{h} (x{nmin})'
        else:
            e = h // d
            typ = TYPES.get((e, nmin), ['?'])
            ok = mn == 2 * H2 and typ != ['?']
            desc = f'{d}c: deep hole of order {e}, {nmin} vertices -> {"/".join(typ)}'
        allok &= ok
        out.append(desc)
    check(allok, f'{name:9s} h={h:2d}: ' + '; '.join(out) + f' [{time.time() - t1:.0f}s]')
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
