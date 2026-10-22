"""Shallow multiples of the Weyl vector: bad primes and Atkin-Lehner involutions (session opus-2026-10-22-S53).

For a Niemeier lattice N (root system X, Coxeter number h, Weyl vector rho_N, k components, glue |N/Q|) and a divisor
d of h, enumerate the vectors of the coset N + d rho_N/h of norm <= 2 + 2/h exactly (integer arithmetic; float
pruning with slack). Checked for all 23 types (THM-4759 (5)):
  d = 1:          min norm 2 + 2/h, attained |N/Q| times (Borcherds, thesis Lemma 2.2.7);
  1 < d < h:      if d is a bad prime of every component of X and gcd(d, h/d) = 1, then the minimum is 2 - 2d/h,
                  attained by a single vector, and there are phi(d) T(d) vectors of norm 2, where T(d) is the
                  number of components with bad prime d; otherwise the minimum is >= 2;
  and the first case happens exactly when the Atkin-Lehner involution W_{h/d} fixes T^X = 1/eta_{pi_X}
  (exact test on the eta exponents: r_{delta * e} = r_delta with delta * e = delta e / gcd(delta, e)^2).
Bad primes: D_n: 2; E6, E7: 2, 3; E8: 2, 3, 5; A_n: none.
Needs cypari and numba. Usage: python multiple_cosets.py [types...]"""
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


src = open(__file__.replace('multiple_cosets.py', 'all_holes.py'), encoding='utf-8').read()
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
def shell(Gf, Gi, Yd, h, R2):
    """vectors a (N-coords) with |a + Yd/h|^2 <= R2 (Yd = coordinates of d rho_N); returns h^2 |a + Yd/h|^2 values and a."""
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
                zu = h * a[u] + Yd[u]
                if zu == 0:
                    continue
                acc = 0
                for w in range(d):
                    acc += Gi[u, w] * (h * a[w] + Yd[w])
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


def phi(m):
    return sum(1 for j in range(1, m + 1) if math.gcd(j, m) == 1)


BADP = {'a': lambda r: [], 'd': lambda r: [2], 'e': lambda r: [2, 3] if r in (6, 7) else [2, 3, 5]}


def coxeter_frame(t, r):
    if t == 'a':
        return {r + 1: 1, 1: -1}
    if t == 'd':
        out = {}
        for a_, m_ in ((2 * r - 2, 1), (2, 1), (r - 1, -1), (1, -1)):
            out[a_] = out.get(a_, 0) + m_
        return out
    return {6: {12: 1, 2: 1, 3: 1, 6: -1, 4: -1, 1: -1},
            7: {18: 1, 3: 1, 2: 1, 9: -1, 6: -1, 1: -1},
            8: {30: 1, 5: 1, 3: 1, 2: 1, 15: -1, 10: -1, 6: -1, 1: -1}}[r]


for name, h, n0, z2 in rows:
    if ONLY and name not in ONLY:
        continue
    t1 = time.time()
    comps = []
    for t, rr, ex in re.findall(r'([ade])(\d+)(?:\^(\d+))?', name):
        comps += [(t, int(rr))] * (int(ex) if ex else 1)
    T = {}
    pi = {}
    for t, rr in comps:
        for p in BADP[t](rr):
            T[p] = T.get(p, 0) + 1
        for a_, m_ in coxeter_frame(t, rr).items():
            pi[a_] = pi.get(a_, 0) + m_
    pi = {a_: m_ for a_, m_ in pi.items() if m_}
    glue2 = 1
    for a_, m_ in pi.items():
        glue2 *= a_ ** m_ if m_ > 0 else 1
    den2 = 1
    for a_, m_ in pi.items():
        den2 *= a_ ** (-m_) if m_ < 0 else 1
    glue = math.isqrt(glue2 // den2)                     # |N/Q| = prod a^(m_a/2)
    assert glue * glue * den2 == glue2
    # Atkin-Lehner involutions fixing 1/eta_pi (exact exponent test)
    S = [e for e in range(1, h + 1) if h % e == 0 and math.gcd(e, h // e) == 1
         and all(pi.get(dl * e // math.gcd(dl, e) ** 2, 0) == m_ for dl, m_ in pi.items())]
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
    Gf = np.array(GN, dtype=np.float64)
    allok = True    # Yr = rho_N + h x for some x in N (choice of z'); only the coset N + rho_N/h matters
    out = []
    for d in [d for d in range(1, h) if h % d == 0]:
        vals, vecs = shell(Gf, Gi, d * Yr, h, 2.0 + 2.0 / h + 1e-6)
        vals = np.array(vals, dtype=np.int64)
        H2 = h * h
        mn = int(vals.min()) if len(vals) else None
        nmin = int((vals == mn).sum()) if len(vals) else 0
        n2 = int((vals == 2 * H2).sum()) if len(vals) else 0
        if d == 1:
            ok = mn == 2 * H2 + 2 * h and nmin == glue
            kind = 'Weyl coset'
        else:
            al = (h // d) in S
            badall = all(d in BADP[t](rr) for t, rr in comps) and math.gcd(d, h // d) == 1
            if al != badall:
                ok = False
            elif al:
                ok = mn == 2 * H2 - 2 * d * h and nmin == 1 and n2 == phi(d) * T.get(d, 0)
            else:
                ok = mn is None or mn >= 2 * H2
            kind = 'SHALLOW (W_%d)' % (h // d) if al else 'deep'
        allok &= ok
        out.append((d, kind, mn, nmin, n2))
    desc = '; '.join(f'd={d}: {kind}, min {"> 2 + 2/h" if mn is None else str(mn) + "/" + str(h * h)} (x{nm}), norm-2 {n2}'
                     for d, kind, mn, nm, n2 in out)
    check(allok, f'{name:9s} h={h:2d} S_X={S} T={T}: {desc} [{time.time() - t1:.0f}s]')
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
