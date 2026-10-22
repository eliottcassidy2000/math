"""Twined theta of a Niemeier lattice by its Weyl vector vs Coxeter-element Frame shapes (session opus-2026-10-22-S53).

Conjecture (tested here). For every Niemeier lattice N with roots, with Coxeter number h, Weyl vector rho_N and root
system X,

    Theta_N^[rho_N/h](tau) := sum_{x in N} e^{2 pi i x.rho_N/h} q^{|x|^2/2} = eta(tau)^24 / eta_{pi_X}(tau),

where pi_X is the product over the components of the Frame shape of the Coxeter element, i.e. of
prod_{exponents m} (x - e^{2 pi i m/h}) written as prod_a (x^a - 1)^{m_a}, and eta_pi = prod_a eta(a tau)^{m_a}.
For pure A-types this is the S-transform of Garvan-Kim-Stanton (THM-4758 (3')). Together with THM-4758 it gives

    Tr_{V_Lambda}(e^{2 pi i c(0)} q^{L_0 - 1}) = 24 + k + 1/eta_{pi_X}(tau).

The functions 1/eta_{pi_X} are compared with the Hauptmoduln of the umbral lambencies (Cheng-Duncan-Harvey).
Needs cypari and numba. Usage: python coxeter_frame_twined.py [NORM=6] [types...]"""
import sys, time, math, cmath
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


src = open(__file__.replace('coxeter_frame_twined.py', 'all_holes.py'), encoding='utf-8').read()
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


def coxeter_frame(comp):
    """Frame shape (dict a -> m_a) of the Coxeter element of an irreducible simply-laced root system."""
    t, r = comp[0], int(comp[1:])
    if t == 'a':
        return {r + 1: 1, 1: -1}
    if t == 'd':
        hh = 2 * r - 2
        fr = {hh: 1, 2: 1, r - 1: -1, 1: -1}
        out = {}
        for a_, m_ in fr.items():
            out[a_] = out.get(a_, 0) + m_
        return out
    if t == 'e':
        return {6: {12: 1, 2: 1, 3: 1, 6: -1, 4: -1, 1: -1},
                7: {18: 1, 3: 1, 2: 1, 9: -1, 6: -1, 1: -1},
                8: {30: 1, 5: 1, 3: 1, 2: 1, 15: -1, 10: -1, 6: -1, 1: -1}}[r]
    raise ValueError(comp)


def parse_components(name):
    import re
    out = []
    for t, r, e in re.findall(r'([ade])(\d+)(?:\^(\d+))?', name):
        out += [t + r] * (int(e) if e else 1)
    return out


def eta_product_series(frame, L):
    """q-expansion (coefficients of q^0..L-1, after removing q^{sum a m_a /24}) of prod_a eta(a tau)^{m_a}."""
    s = [0] * L; s[0] = 1
    for a_, m_ in frame.items():
        for mm in range(1, L // a_ + 2):
            for _ in range(abs(m_)):
                if m_ > 0:
                    s = [s[t] - (s[t - a_ * mm] if t >= a_ * mm else 0) for t in range(L)]
                else:
                    for t in range(a_ * mm, L):
                        s[t] += s[t - a_ * mm]
    shift = sum(a_ * m_ for a_, m_ in frame.items())   # q^(shift/24)
    return s, shift


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
    cnt = twisted_counts(np.array(GN, dtype=np.float64), Gi, (Gi @ Yr).astype(np.int64), h, NORM)
    th = [sum(int(cnt[m][r]) * cmath.exp(2j * math.pi * r / h) for r in range(h)) for m in range(NORM // 2 + 1)]
    # prediction eta^24 / eta_pi: eta^24 = q * prod(1-q^n)^24 ; 1/eta_pi = q^{-shift/24} prod(...)^{-m}
    frame = {}
    for comp in parse_components(name):
        for a_, m_ in coxeter_frame(comp).items():
            frame[a_] = frame.get(a_, 0) + m_
    frame = {a_: m_ for a_, m_ in frame.items() if m_}
    L = NORM // 2 + 3
    inv = {a_: -m_ for a_, m_ in frame.items()}
    inv[1] = inv.get(1, 0) + 24
    ser, shift = eta_product_series(inv, L)       # eta^24 / eta_pi = q^{shift/24} * ser
    ok_shift = shift == 0                        # expect weight-12 series starting at q^0
    pred = ser[:NORM // 2 + 1]
    ok = ok_shift and all(abs(th[m].imag) < 1e-6 and abs(th[m].real - pred[m]) < 1e-6 for m in range(NORM // 2 + 1))
    frame_str = ' '.join(f'{a_}^{m_}' for a_, m_ in sorted(frame.items(), key=lambda x: -x[0]))
    check(ok, f'{name:9s} h={h:2d}: Theta_N^[rho/h] = eta^24/eta_pi, pi = {frame_str}; enumeration {[round(t.real) for t in th]} vs {pred} [{time.time() - t1:.0f}s]')
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
