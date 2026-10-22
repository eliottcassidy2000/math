"""Shallow multiples of Leech deep holes twine the Leech VOA by Fricke Hauptmoduln (session opus-2026-10-22-S53).

Conjecture tested (HYP-9368). Let c be a deep hole of type X with Coxeter number h. Let d be a common bad prime of X
with d || h (THM-4759 (5a)), and put v = d c, an element of order e = h/d. Then

    Z(tau) := Tr_{V_Lambda}(e^{2 pi i v(0)} q^{L_0 - 1}) = Theta_Lambda^[v](tau) / eta(tau)^24

is invariant under the Fricke involution tau -> -1/(e tau). Its only poles on X_0(e) are simple poles at infinity
and at 0; the one at 0 is the unique shallow vector of Lambda + v. Hence Z - 24 is the normalised Hauptmodul of
Gamma_0(e)+, a monstrous Fricke group.

Fricke invariance is equivalent to the exact q-series identity

    Theta_Lambda^[v](tau) / Delta(tau) = theta_{Lambda + v}(e tau) / Delta(e tau),

since Z(-1/(e tau)) = theta_{Lambda+v}(e tau)/eta(e tau)^24 by Poisson summation. The left side is computed from
Leech vectors of norm <= NORM weighted by e^{2 pi i x.v}, the right side from the vectors of Lambda + v. The script
compares the two Laurent series through q^(NORM/2 - 1).
Needs cypari and numba. Usage: python fricke_multiples.py [NORM=8] [types...]"""
import sys, time, math, cmath, re
import numpy as np
from cypari import pari
from numba import njit

NORM = int(sys.argv[1]) if len(sys.argv) > 1 else 8
ONLY = sys.argv[2:]
pari.allocatemem(2 * 10 ** 9)
FAIL = []


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


src = open(__file__.replace('fricke_multiples.py', 'all_holes.py'), encoding='utf-8').read()
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
Gi = np.array(Gm, dtype=np.int64)
Gf = np.array(Gm, dtype=np.float64)


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


@njit(cache=True)
def shell(Gf, Gi, Yd, h, R2):
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


def inv_delta(L, step=1):
    """coefficients of q^{-step} prod (1 - q^{step n})^{-24} = 1/Delta(q^step) as Laurent series: list c[k] for
    exponents -step + k, k = 0..L-1."""
    s = [0] * L; s[0] = 1
    for m in range(1, L // step + 2):
        for _ in range(24):
            for t in range(step * m, L):
                s[t] += s[t - step * m]
    return s


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
    cnt = twisted_counts(Gf, Gi, (Gi @ Yh).astype(np.int64), h, NORM)
    M = NORM // 2                     # Theta^[v] known through q^M
    for d in shallow:
        e = h // d
        # left: Theta^[d c] / Delta, Laurent coefficients for exponents -1 .. M-1
        th = [sum(int(cnt[m][r]) * cmath.exp(2j * math.pi * d * r / h) for r in range(h)) for m in range(M + 1)]
        ok_real = all(abs(z.imag) < 1e-6 and abs(z.real - round(z.real)) < 1e-6 for z in th)
        th = [round(z.real) for z in th]
        idl = inv_delta(M + 1)                                     # 1/Delta = q^-1 (1 + 24 q + ...)
        left = [sum(th[i] * idl[k - i] for i in range(0, k + 1)) for k in range(M + 1)]   # exponents -1 + k
        # right: theta_{Lambda + d c}(e tau) / Delta(e tau)
        R2 = 2.0 * (M + e) / e + 1e-6
        vals = np.array(shell(Gf, Gi, d * Yh, h, R2), dtype=np.int64)
        H2 = h * h
        # exponent of x in theta(e tau): e |x|^2 / 2 = e * vals / (2 h^2)
        expo = {}
        okint = True
        for s in vals:
            num = e * int(s)
            if num % (2 * H2):
                okint = False
                continue
            k = num // (2 * H2)
            expo[k] = expo.get(k, 0) + 1
        idle = inv_delta(e * (M + 2), e)                         # 1/Delta(q^e) = q^-e (1 + 24 q^e + ...)
        right = []
        for k in range(M + 1):                                   # exponent -1 + k
            tot = 0
            for ex, cn in expo.items():
                j = (-1 + k) - (ex - e)                          # need q^{ex - e} * q^{j}, j >= 0 index into idle
                if 0 <= j < len(idle):
                    tot += cn * idle[j]
            right.append(tot)
        same = ok_real and okint and left == right
        series = ' + '.join(f'{c}q^{k - 1}' for k, c in enumerate(left))
        check(same, f'{name:9s} v = {d}c (order {e}): Theta^[v]/Delta = theta_(Lambda+v)(e tau)/Delta(e tau) through '
                    f'q^{M - 1} (Fricke W_{e} invariance); Z = {series} [{time.time() - t1:.0f}s]')
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
