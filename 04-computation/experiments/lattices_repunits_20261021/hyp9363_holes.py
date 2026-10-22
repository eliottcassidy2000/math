"""HYP-9363: theta series of the Leech lattice around its pure-A deep holes (session opus-2026-10-21-S52).

Conjecture. Let c be a deep hole of type A_(h-1)^k, where k(h-1) = 24, so h is one of 2, 3, 4, 5, 7, 9, 13, 25. Then

    theta_(Lambda + c)(tau) = h^(k/2) (eta(tau)^h / eta(tau/h))^k + (24 + k) Delta(tau)
                            = theta_(N + rho_N/h)(tau) + (24 + k) Delta(tau).

The second form uses two facts:
- h^(k/2) = |glue code| of N = N(A_(h-1)^k);
- theta_(A_(h-1) + rho/h) = eta(tau)^h / eta(tau/h), the h-core generating function (Garvan-Kim-Stanton 1990).

The coefficient 24 + k is the number of hole vertices, and it equals 24h/(h-1) = (roots of N)/(h-1).

Modes:
- `A1 KMAX`: c = v/2 for a Leech vector v of norm 8; the A_1^24 holes (SPLAG ch. 23).
- `A2 KMAX`: c = (a + b)/3 for Leech vectors with |a|^2 = |b|^2 = |a - b|^2 = 6, chosen so that c has 36 nearest points
  at distance sqrt 2 and none closer: the centroid of an A~_2 triangle of a type-A_2^12 hole.
- `A24 KMAX`: c = the A_24 hole of THM-4752 (nu = (1/2^25 | 5/2)); this repeats hyp9363_deep_hole_theta.py.

Leech points near c are enumerated exactly (Fincke-Pohst, exact integer norms), and points closer than sqrt 2 are
counted separately (there must be none). Needs cypari and numba.
Prints ALL CHECKS PASSED or the first disagreement."""
import sys, time
import numpy as np
from cypari import pari
from numba import njit

MODE = sys.argv[1] if len(sys.argv) > 1 else 'A2'
KMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 4
pari.allocatemem(2 * 10 ** 9)
T0 = time.time()
FAIL = []


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


n = 25


def lor2(x, y):
    return int(sum(int(a) * int(b) for a, b in zip(x[:-1], y[:-1])) - int(x[-1]) * int(y[-1]))


gens = []
for i in range(n - 1):
    v = [0] * (n + 1); v[i] = 2; v[i + 1] = -2; gens.append(v)
v = [0] * (n + 1); v[n - 2] = 2; v[n - 1] = 2; gens.append(v)
gens.append([1] * (n + 1))
v = [0] * (n + 1); v[n - 1] = 2; v[n] = 2; gens.append(v)
v = [0] * (n + 1); v[0] = 4; gens.append(v)
A = pari.matrix(n + 1, len(gens), [gens[j][i] for i in range(n + 1) for j in range(len(gens))])
Hb = pari.mathnf(A)
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
check(int(pari.matdet(Gp)) == 1 and int(pari.qfminim(Gp, 2, 0)[0]) == 0 and all(Gm[i][i] % 2 == 0 for i in range(24)),
      'Leech lattice Lambda = {rho, rho\'}^perp rebuilt (even, unimodular, rootless)')
Gf = np.array(Gm, dtype=np.float64)
Gi = np.array(Gm, dtype=np.int64)


def ip(x, y):
    return sum(x[a] * Gm[a][b] * y[b] for a in range(24) for b in range(24))


@njit(cache=True)
def enumerate_hole(Gf, Gi, Yv, h, R2, KMAX):
    d = Gf.shape[0]
    q = Gf.copy()
    for i in range(d):
        for j in range(i + 1, d):
            q[j, i] = q[i, j]
            q[i, j] = q[i, j] / q[i, i]
        for k in range(i + 1, d):
            for l in range(k, d):
                q[k, l] -= q[k, i] * q[i, l]
    y = Yv.astype(np.float64) / h
    counts = np.zeros(KMAX + 1, dtype=np.int64)
    closer = 0
    a = np.zeros(d, dtype=np.int64)
    T = np.zeros(d + 1)
    cen = np.zeros(d)
    hi = np.zeros(d, dtype=np.int64)
    i = d - 1
    cen[i] = y[i]
    rad = np.sqrt(max(R2 + 1e-9 - T[i + 1], 0.0) / q[i, i])
    a[i] = np.int64(np.ceil(cen[i] - rad))
    hi[i] = np.int64(np.floor(cen[i] + rad))
    total = 0
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
                zu = h * a[u] - Yv[u]
                if zu == 0:
                    continue
                acc = 0
                for w in range(d):
                    acc += Gi[u, w] * (h * a[w] - Yv[w])
                s += zu * acc
            if s < 2 * h * h:
                closer += 1
            elif (s - 2 * h * h) % (2 * h) == 0:
                kk = (s - 2 * h * h) // (2 * h)
                if kk <= KMAX:
                    counts[kk] += 1
            total += 1
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
    return counts, closer, total


if MODE == 'A1':
    h, kc = 2, 24
    sv = pari.qfminim(Gp, 8, 400, 2)[2]
    v8 = next([int(sv[i, j]) for i in range(24)] for j in range(int(pari.matsize(sv)[1]))
              if ip([int(sv[i, j]) for i in range(24)], [int(sv[i, j]) for i in range(24)]) == 8)
    Y = v8
    print(f'     c = v/2, v = {v8} (norm 8)')
elif MODE == 'A2':
    h, kc = 3, 12
    sv = pari.qfminim(Gp, 6, 60000, 2)[2]
    six = [[int(sv[i, j]) for i in range(24)] for j in range(int(pari.matsize(sv)[1]))]
    six = [x for x in six if ip(x, x) == 6]
    Y = None
    tried = 0
    for a_ in six[:40]:
        for b_ in six:
            for sgn in (1, -1):
                bb = [sgn * t for t in b_]
                if ip(a_, bb) != 3:
                    continue
                cand = [p + q_ for p, q_ in zip(a_, bb)]          # 3 c = a + b
                cnt, cl, tot = enumerate_hole(Gf, Gi, np.array(cand, dtype=np.int64), 3, 2.0, 0)
                tried += 1
                if cl == 0 and cnt[0] == 36:
                    Y = cand
                    break
            if Y is not None:
                break
        if Y is not None:
            break
    check(Y is not None, f'found an A_2^12 deep hole c = (a + b)/3 (36 nearest points, none closer) after {tried} candidate triangles')
elif MODE == 'A24':
    h, kc = 25, 1
    nu2 = [1] * n + [5]
    lnu2 = [a - 25 * b - 25 * c_ for a, b, c_ in zip(nu2, rho2, rhop2)]
    g = pari.matrix(24, 1, [lor2(a, lnu2) // 4 for a in lam])
    Ycol = pari.matsolve(Gp, g)
    Y = [int(Ycol[i, 0]) for i in range(24)]
else:
    raise SystemExit('mode must be A1, A2 or A24')

Yv = np.array(Y, dtype=np.int64)
R2 = 2.0 + 2.0 * KMAX / h
counts, closer, total = enumerate_hole(Gf, Gi, Yv, h, R2, KMAX)
print(f'     {MODE}: enumerated {total} Leech points with |lambda - c|^2 <= {R2:.2f}; closer than sqrt 2: {closer} [{time.time() - T0:.0f}s]')

L = h + KMAX + 2


def mul(a, b):
    return [sum(a[i] * b[t - i] for i in range(t + 1)) for t in range(L)]


def prod_pow(step, e):           # prod_n (1 - Q^(step n))^e, e may be negative
    s = [0] * L; s[0] = 1
    for m in range(1, L // step + 2):
        for _ in range(abs(e)):
            if e > 0:
                s = [s[t] - (s[t - step * m] if t >= step * m else 0) for t in range(L)]
            else:
                for t in range(step * m, L):
                    s[t] += s[t - step * m]
    return s


# theta = sum_k N_k Q^(h + k), Q = q^(1/h)
D = prod_pow(h, 24)                                     # (24 + kc) Q^h prod(1 - Q^(hn))^24
E = mul(prod_pow(h, h * kc), prod_pow(1, -kc))           # h^(kc/2) Q^(h+1) prod(1-Q^(hn))^(h kc) / prod(1-Q^n)^kc
glue = int(round(h ** (kc / 2)))
pred = [(24 + kc) * D[k] + (glue * E[k - 1] if k >= 1 else 0) for k in range(KMAX + 1)]
print(f'     prediction: {glue} (eta^{h}/eta(tau/{h}))^{kc} + {24 + kc} Delta')
print('      k   N_k (exact)    prediction')
for k in range(KMAX + 1):
    print(f'     {k:3d}  {int(counts[k]):13d}  {pred[k]:13d}' + ('' if counts[k] == pred[k] else '   <-- differs'))
check(closer == 0 and int(counts[0]) == 24 + kc, f'no Leech point closer than sqrt 2; {24 + kc} at distance sqrt 2 (vertices of A~_{h - 1}^{kc})')
check(all(int(counts[k]) == pred[k] for k in range(KMAX + 1)),
      f'theta_(Lambda + c) = {glue} (eta^{h}/eta(tau/{h}))^{kc} + {24 + kc} Delta through q^(1 + {KMAX}/{h})')
print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f' ({time.time() - T0:.0f}s)')
