"""Fricke descent = roots of height divisible by p (session opus-2026-10-22-S53; THM-4759 (8)).

For each Niemeier lattice N and each common bad prime p || h, the roots of N whose height alpha.rho_N is divisible
by p form a root system X_p, the fixed points of the principal element of order p. The script checks that:
  - |X_p| = 24 (h/p - 1), i.e. sum over components of 2 sum_i floor(m_i/p) = l (h/p - 1), which is equivalent to
    sum_i (m_i mod p) = l p/2;
  - every component of X_p has Coxeter number h/p - 1, and the type agrees with the Fricke descent found on the Leech
    side (fricke_neighbors.out): D_n -> D_{n/2}^2, E7 -> A7, E8 -> D8 / A8 / A4^2, E6 -> A2^3.
This identifies the Fricke descent with the classical fixed-point subsystem of exp(2 pi i rho/p), i.e. with
Borel-de Siebenthal.
Needs cypari and numba. Usage: python height_descent.py"""
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


src = open(__file__.replace('height_descent.py', 'all_holes.py'), encoding='utf-8').read()
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



from collections import Counter
EXPECT = {('d4^6', 2): 'A1^24', ('d6^4', 2): 'A3^8', ('e6^4', 3): 'A2^12', ('d8^3', 2): 'D4^6', ('d10e7^2', 2): 'A7^2D5^2',
          ('d12^2', 2): 'D6^4', ('d16e8', 2): 'D8^3', ('e8^3', 2): 'D8^3', ('e8^3', 3): 'A8^3', ('e8^3', 5): 'A4^6',
          ('d24', 2): 'D12^2'}


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
    rho = Yr - h * a0
    # roots of N
    mv = pari.qfminim(GNp, 2, 10 ** 5, 0)
    R = np.array([[int(mv[2][i, j]) for i in range(24)] for j in range(int(mv[0]) // 2)], dtype=np.int64)
    R = np.vstack([R, -R])
    ht = R @ Gi @ rho
    ok_ht = len(R) == 24 * h and int(np.abs(ht).max()) == h - 1 and int((ht > 0).sum()) == 12 * h
    for p in shallow:
        Xp = R[ht % p == 0]
        IP = Xp @ Gi @ Xp.T
        adj = IP != 0
        seen = [False] * len(Xp)
        sizes = []
        for s0 in range(len(Xp)):
            if seen[s0]:
                continue
            stack = [s0]; seen[s0] = True; c = 0
            while stack:
                u = stack.pop(); c += 1
                for w_ in np.nonzero(adj[u])[0]:
                    if not seen[w_]:
                        seen[w_] = True; stack.append(int(w_))
            sizes.append(c)
        hp = h // p - 1
        typ = name_of(sorted(sizes, reverse=True), hp)
        exp_typ = EXPECT[(name, p)]
        same = sorted(typ.replace('^', ' ').split()) == sorted(exp_typ.replace('^', ' ').split()) or typ == exp_typ
        check(ok_ht and len(Xp) == 24 * hp and typ == exp_typ,
              f'{name:9s} p={p}: roots of height divisible by {p}: {len(Xp)} = 24(h/p - 1) = {24 * hp}; type {typ} '
              f'(Fricke descent on the Leech side: {exp_typ})')
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
