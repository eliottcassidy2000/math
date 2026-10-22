"""Order-6 Fricke points of the Leech lattice and their orbifold neighbours (session opus-2026-10-22-S53; HYP-9368).

An order-6 point v is v = a - b, with a = 3v of order 2 and b = -2v of order 3 (CRT). If v is a Fricke point, its
proper multiples are deep holes: a is an A1^24 hole and b an A2^12 hole. This script fixes the A1^24 hole a and samples
A2^12 holes b = y/3, where y is a norm-18 vector whose class mod 3 Lambda has minimal norm 18. It tests v = a - b for
the Fricke property (unique minimal vector, of norm 5/3) and builds the 6-neighbour of Lambda at each Fricke sample.
HYP-9368 predicts that every such neighbour is A4^6, with 120 roots; that is dim V'_1 = 24*6, the composite-order step
missing from the prime-order argument.
Usage: python order6_fricke_census.py [TRIALS=20000]"""
import sys, time, math, re
from fractions import Fraction as Fr
import numpy as np
from cypari import pari
from numba import njit

TRIALS = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
pari.allocatemem(2 * 10 ** 9)
FAIL = []


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


src = open(__file__.replace('order6_fricke_census.py', 'all_holes.py'), encoding='utf-8').read()
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



rowd = {name: (h, n0, z2) for name, h, n0, z2 in rows}


def hole_coords(name):
    h, n0, z2 = rowd[name]
    zrho = lor2(z2, rho2) // 4
    zrp = lor2(z2, rhop2) // 4
    lz2 = [a + zrp * b + zrho * c for a, b, c in zip(z2, rho2, rhop2)]
    sol = pari.matsolve(Gp, pari.matrix(24, 1, [lor2(a, lz2) // 4 for a in lam]))
    return np.array([int(sol[i, 0]) for i in range(24)], dtype=np.int64)       # coordinates of h c


Ya = hole_coords('a1^24')            # 2a  (a = Ya/2)
mv = pari.qfminim(Gp, 4, 200000, 0)
L4 = np.array([[int(mv[2][i, j]) for i in range(24)] for j in range(int(mv[0]) // 2)], dtype=np.int64)
L4 = np.vstack([L4, -L4])
GL4 = L4 @ G
print(f'       |Lambda_4| = {len(L4)}', flush=True)
rng = np.random.default_rng(20261024)


def min_shell(Y, den, bound):
    sv, sa = shell_vecs(Gf, G, Y, den, bound)
    vals = [int(s) for s in sv]
    return vals, sa


def neighbour_roots(Yv, e):
    """roots of the e-neighbour Lambda_0 + Z w at v = Yv/e (Yv = coordinates of e v)."""
    W = None
    cands = [np.zeros(24, dtype=np.int64)] + [np.eye(24, dtype=np.int64)[i] * s for i in range(24) for s in (1, -1)]
    for i in range(24):
        for k in range(i, 24):
            cands.append(np.eye(24, dtype=np.int64)[i] + np.eye(24, dtype=np.int64)[k])
    for lam_c in cands:
        Wv = Yv + e * lam_c
        nw = int(Wv @ G @ Wv)
        if nw % (e * e) == 0 and (nw // (e * e)) % 2 == 0:
            W = Wv; break
    if W is None:
        return None
    nr = 0
    for j in range(1, e):
        sv, sa = shell_vecs(Gf, G, j * Yv, e, 2.0 + 1e-6)
        for s, a in zip(sv, sa):
            if int(s) != 2 * e * e:
                continue
            X = e * a + j * Yv
            if int((X - j * W) @ G @ Yv) % (e * e) == 0:
                nr += 1
    return nr


found = {}
n_holes = 0
t0 = time.time()
for trial in range(TRIALS):
    # random norm-18 vector y = x1 + x2 + x3 with x_i in Lambda_4 and pairwise inner products 1
    x1 = L4[rng.integers(0, len(L4))]
    ip1 = GL4 @ x1
    cand2 = np.nonzero(ip1 == 1)[0]
    x2 = L4[cand2[rng.integers(0, len(cand2))]]
    ip2 = GL4 @ x2
    cand3 = np.nonzero((ip1 == 1) & (ip2 == 1))[0]
    if len(cand3) == 0:
        continue
    x3 = L4[cand3[rng.integers(0, len(cand3))]]
    y = x1 + x2 + x3                                   # norm 12 + 2*3 = 18
    vals, _ = min_shell(y, 3, 2.0 + 1e-6)              # Lambda + y/3 up to norm 2
    if not vals or min(vals) != 18 or sum(1 for s in vals if s == 18) != 36:
        continue                                       # b = y/3 must be an A2^12 deep hole (36 vertices)
    n_holes += 1
    # v = a - b, coordinates of 6 v: 6a - 6b = 3 Ya - 2 y
    Yv = 3 * Ya - 2 * y
    vals6, _ = min_shell(Yv, 6, 5.0 / 3 + 1e-6)
    if not vals6:
        continue
    mn = min(vals6)
    if mn != 60 or sum(1 for s in vals6 if s == 60) != 1:
        continue                                       # not a Fricke point (need unique norm 5/3 = 60/36)
    nr = neighbour_roots(Yv, 6)
    found[nr] = found.get(nr, 0) + 1
print(f'       {TRIALS} trials, {n_holes} A2^12 holes b, Fricke points v = a - b by neighbour root count: {found} '
      f'[{time.time() - t0:.0f}s]', flush=True)
check(len(found) >= 1 and set(found) == {120},
      'every sampled order-6 Fricke point has the 6-neighbour with 120 roots (Coxeter number 5, i.e. A4^6)')
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
