"""The Kneser neighbour of the Leech lattice at a Fricke point (session opus-2026-10-22-S53; THM-4759 (7)).

For a shallow multiple v = d c (order e) the inner automorphism g_v = e^{2 pi i v(0)} has type 0. Its cyclic
orbifold is the lattice VOA of the e-neighbour

    Lambda' = Lambda_0 + Z w,   Lambda_0 = {x in Lambda : x.v in Z},

where w is in v + Lambda and |w|^2 is even. The roots of Lambda' are the norm-2 vectors x of Lambda + j v, j = 1..e-1,
with (x - j w).v in Z. This script collects them exactly and reads off the root system of Lambda'. That is the
Niemeier lattice obtained by orbifolding V_Lambda along the Fricke point.
Needs cypari and numba. Usage: python fricke_neighbors.py [types...]"""
import sys, time, math, re
from fractions import Fraction as Fr
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


src = open(__file__.replace('fricke_neighbors.py', 'all_holes.py'), encoding='utf-8').read()
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
rhop2 = list(rho2); rhop2[0] += 2; rhop2[1] -= 2
F = pari.matrix(2, n + 1, [lor2(b, rho2) // 4 for b in basis] + [lor2(b, rhop2) // 4 for b in basis])
K = pari.matkerint(F)
lam = [[sum(int(K[i, j]) * basis[i][t] for i in range(n + 1)) for t in range(n + 1)] for j in range(24)]
GL = pari.matrix(24, 24, [lor2(a, b) // 4 for a in lam for b in lam])
U = pari.qflllgram(GL)
lam = [[sum(int(U[i, j]) * lam[i][t] for i in range(24)) for t in range(n + 1)] for j in range(24)]
Gm = [[lor2(a, b) // 4 for b in lam] for a in lam]
Gp = pari.matrix(24, 24, [x for row in Gm for x in row])
G = np.array(Gm, dtype=np.int64)
Gf = G.astype(np.float64)


@njit(cache=True)
def shell_vecs(Gf, Gi, Yd, den, R2):
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


def root_system(R, h_scale):
    """R: list of integer vectors (scaled by h_scale) of norm 2; returns sorted list of component sizes."""
    m = len(R)
    Rm = np.array(R, dtype=np.int64)
    IP = Rm @ G @ Rm.T                                  # h_scale^2 * inner products
    adj = IP != 0
    seen = [False] * m
    comps = []
    for s0 in range(m):
        if seen[s0]:
            continue
        stack = [s0]; seen[s0] = True; c = 0
        while stack:
            u = stack.pop(); c += 1
            for w_ in np.nonzero(adj[u])[0]:
                if not seen[w_]:
                    seen[w_] = True; stack.append(int(w_))
        comps.append(c)
    return sorted(comps, reverse=True)


def name_of(sizes, hh):
    out = []
    for s in sizes:
        if s == hh * (hh - 1):
            out.append(f'A{hh - 1}')
        elif hh % 2 == 0 and s == 2 * (hh // 2 + 1) * (hh // 2):
            out.append(f'D{hh // 2 + 1}')
        elif (s, hh) == (72, 12):
            out.append('E6')
        elif (s, hh) == (126, 18):
            out.append('E7')
        elif (s, hh) == (240, 30):
            out.append('E8')
        else:
            out.append(f'?({s})')
    from collections import Counter
    c = Counter(out)
    return ''.join(f'{k}^{v}' if v > 1 else k for k, v in sorted(c.items()))


BADP = {'a': lambda r: [], 'd': lambda r: [2], 'e': lambda r: [2, 3] if r in (6, 7) else [2, 3, 5]}
for name, h, n0, z2 in rows:
    comps = []
    for t, rr, ex in re.findall(r'([ade])(\d+)(?:\^(\d+))?', name):
        comps += [(t, int(rr))] * (int(ex) if ex else 1)
    shallow = [d for d in range(2, h) if h % d == 0 and math.gcd(d, h // d) == 1
               and all(d in BADP[t](rr) for t, rr in comps)]
    if not shallow or (ONLY and name not in ONLY):
        continue
    t1 = time.time()
    zrho = lor2(z2, rho2) // 4
    zrp = lor2(z2, rhop2) // 4
    lz2 = [a + zrp * b + zrho * c for a, b, c in zip(z2, rho2, rhop2)]
    sol = pari.matsolve(Gp, pari.matrix(24, 1, [lor2(a, lz2) // 4 for a in lam]))
    Yh = np.array([int(sol[i, 0]) for i in range(24)], dtype=np.int64)        # coordinates of h c
    for d in shallow:
        e = h // d
        Yv = d * Yh                                     # coordinates of h v  (v = d c), so v = Yv / h
        # w = v + lam0 with |w|^2 even: search lam0 among small Leech vectors
        vv = Fr(int(Yv @ G @ Yv), h * h)
        best_w = None
        for lam_c in [np.zeros(24, dtype=np.int64)] + [np.eye(24, dtype=np.int64)[i] * s for i in range(24) for s in (1, -1)]:
            Wv = Yv + h * lam_c
            nw = Fr(int(Wv @ G @ Wv), h * h)
            if nw.denominator == 1 and nw.numerator % 2 == 0:
                best_w = Wv; break
        if best_w is None:                              # try sums of two basis vectors
            for i in range(24):
                for k in range(24):
                    lam_c = np.eye(24, dtype=np.int64)[i] + np.eye(24, dtype=np.int64)[k]
                    Wv = Yv + h * lam_c
                    nw = Fr(int(Wv @ G @ Wv), h * h)
                    if nw.denominator == 1 and nw.numerator % 2 == 0:
                        best_w = Wv; break
                if best_w is not None:
                    break
        roots = []
        for j in range(1, e):
            sv, sa = shell_vecs(Gf, G, (j * Yv) % (h * 10 ** 6) if False else j * Yv, h, 2.0 + 1e-6)
            for s, a in zip(sv, sa):
                if int(s) != 2 * h * h:
                    continue
                X = h * a + j * Yv                       # h x
                # (x - j w) . v in Z  <=>  (X - j W) . Yv / h^2 in Z
                ipv = int((X - j * best_w) @ G @ Yv)
                if ipv % (h * h) == 0:
                    roots.append(X)
        nr = len(roots)
        hh = nr // 24 if nr % 24 == 0 else None
        sizes = root_system(roots, h) if roots else []
        typ = name_of(sizes, hh) if hh else '?'
        check(hh is not None and '?' not in typ,
              f'{name:9s} v = {d}c (order {e}): neighbour Lambda\' has {nr} roots, Coxeter number {hh}: {typ} [{time.time() - t1:.0f}s]')
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
