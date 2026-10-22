"""The deep-hole automorphism's trace on the Leech VOA is 24 + a Hauptmodul (session opus-2026-10-22-S53; THM-4758).

Applying tau -> -1/tau to THM-4758 (2) (Poisson summation; Lambda and N are unimodular) gives

    Theta_Lambda^[c](tau) - Theta_N^[rho_N/h](tau) = (24 + k) Delta,
    where  Theta_L^[x](tau) = sum_{y in L} e^{2 pi i y.x} q^{|y|^2/2}.

For pure A-types A_(h-1)^k, the S-transform of Garvan-Kim-Stanton gives Theta_N^[rho/h] = eta(tau)^(hk)/eta(h tau)^k.
Hence

    Theta_Lambda^[c] / eta^24 = 24 + T_h(tau),   T_h = (eta(tau)/eta(h tau))^(24/(h-1)) + 24/(h-1),

where T_h is the normalised Hauptmodul of Gamma_0(h). For h = 2, 3, 4, 5, 7, 9, 13 it is the Monster McKay-Thompson
series T_2B, T_3B, T_4C, T_5B, T_7B, T_9B, T_13B. In VOA language: Tr_{V_Lambda}(e^{2 pi i c(0)} q^{L_0 - 1})
= 24 + T_h(tau). This trace equals the trace of a fixed-point-free Conway element of Frame shape h^k / 1^k plus 24 + k;
the difference is supported on weight 1 only.

Checked here by direct enumeration of Leech vectors (characters e(lambda . c)) for the A1^24, A2^12, A4^6, A6^4 and
A24 holes. Needs cypari and numba."""
import sys, time, math, cmath
import numpy as np
from cypari import pari
from numba import njit

NORM = int(sys.argv[1]) if len(sys.argv) > 1 else 8
TYPES = sys.argv[2:] or ['a1^24', 'a2^12', 'a4^6', 'a6^4', 'a24']
pari.allocatemem(2 * 10 ** 9)
FAIL = []


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


src = open(__file__.replace('twined_hauptmodul.py', 'all_holes.py'), encoding='utf-8').read()
TABLE = src.split('TABLE = """\\\n')[1].split('"""')[0]
n = 25


def lor2(x, y):
    return int(sum(int(a) * int(b) for a, b in zip(x[:-1], y[:-1])) - int(x[-1]) * int(y[-1]))


rows = {}
for line in TABLE.strip('\n').split('\n'):
    t = line.split()
    co = [int(x) for x in t[6:]]
    rows[t[4]] = (int(t[0]), int(t[1]), [co[1 + 24 - i] for i in range(n)] + [co[0]])
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


@njit(cache=True)
def twisted_counts(Gf, Gi, GY, h, Bnd):
    """sum over Leech vectors y (coordinates a) with |y|^2 <= Bnd of e(y . c), c = Yh/h in the dual-coordinate sense:
    y . c = (a^T G Yh)/h. Returns counts[norm/2][residue mod h]."""
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
    pass
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


def series_prod(L, terms):     # prod over (step, exponent) of prod_n (1 - q^(step n))^e, coefficients 0..L-1
    s = [0] * L; s[0] = 1
    for step, e in terms:
        for mm in range(1, L // step + 2):
            for _ in range(abs(e)):
                if e > 0:
                    s = [s[t] - (s[t - step * mm] if t >= step * mm else 0) for t in range(L)]
                else:
                    for t in range(step * mm, L):
                        s[t] += s[t - step * mm]
    return s


Gf = np.array(Gm, dtype=np.float64)
for name in TYPES:
    t1 = time.time()
    h, n0, z2 = rows[name]
    kc = 24 // (h - 1)
    zrho = lor2(z2, rho2) // 4
    zrp = lor2(z2, rhop2) // 4
    lz2 = [a + zrp * b + zrho * c for a, b, c in zip(z2, rho2, rhop2)]      # h c (doubled coordinates)
    g = pari.matrix(24, 1, [lor2(a, lz2) // 4 for a in lam])
    sol = pari.matsolve(Gp, g)
    Yh = np.array([int(sol[i, 0]) for i in range(24)], dtype=np.int64)        # coordinates of h c
    cnt = twisted_counts(Gf, Gi, (Gi @ Yh).astype(np.int64), h, NORM)
    # Theta^[c] coefficients: sum_r cnt[m][r] e(r/h)
    th = []
    for m in range(NORM // 2 + 1):
        z = sum(int(cnt[m][r]) * cmath.exp(2j * math.pi * r / h) for r in range(h))
        th.append(z)
    # prediction: (24 + k) Delta + eta^(24 + k)/eta(h tau)^k, the eta-quotient starting at q^0
    L = NORM // 2 + 2
    D = series_prod(L, [(1, 24)])
    E = series_prod(L, [(1, 24 + kc), (h, -kc)])
    pred = [E[0]] + [(24 + kc) * D[m - 1] + E[m] for m in range(1, L - 1)]
    ok = all(abs(th[m].imag) < 1e-6 and abs(th[m].real - pred[m]) < 1e-6 for m in range(NORM // 2 + 1))
    check(ok, f'{name}: Theta_Lambda^[c] = (24+{kc}) Delta + eta^{24 + kc}/eta({h} tau)^{kc} through q^{NORM // 2}: '
              f'{[round(t.real) for t in th]} vs {pred[:NORM // 2 + 1]}  [{time.time() - t1:.0f}s]')
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
