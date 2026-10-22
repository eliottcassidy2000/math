"""Order-4 Fricke points of the Leech lattice (session opus-2026-10-22-S53; HYP-9368 at composite order 4).

A point v of order 4 is a Fricke point when three things hold: min |Lambda + v|^2 = 2 - 2/4 = 3/2, attained once; 2v
is a deep hole (min norm 2); and 3v = -v. With y = 4x for the minimal vector x, these say that y has norm 24, that y
is the unique minimum of y + 4 Lambda, and that y/2 is a deep hole, i.e. y + 2 Lambda has minimal norm 8. For norm-24
y in a type-4 class mod 2 Lambda, uniqueness reduces to max_{mu in Lambda_4} y.mu <= 7: norm-6 vectors cannot reach
y.lambda = 12 unless y = 2 lambda.

The samples are y = u + 2w, with u = mu1 + mu2 (mu_i in Lambda_4, mu1.mu2 = 0, so |u|^2 = 8) and w = mu3 + mu4 of
norm 6 or 8 with u.w = 4 - |w|^2. Then |y|^2 = 24 and y = u mod 2 Lambda. (Taking w in Lambda_4 forces y.w = 8, so
such samples are never Fricke.) The orbit invariant is the distribution of y.mu over
Lambda_4. HYP-9368 predicts a single invariant among the Fricke samples: the (h+1)-point of A2^12, which is the
shallow multiple 3c of E6^4. Usage: python order4_fricke_census.py [SAMPLES=3000]"""
import sys, time, math
import numpy as np
from cypari import pari
from numba import njit

NS = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
pari.allocatemem(2 * 10 ** 9)
FAIL = []


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


src = open(__file__.replace('order4_fricke_census.py', 'all_holes.py'), encoding='utf-8').read()
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
G = np.array(Gm, dtype=np.int64)

# Lambda_4 (one of each +- pair suffices for the distribution up to sign; we keep all for clarity)
mv = pari.qfminim(Gp, 4, 200000, 0)
assert int(mv[0]) == 196560
L4 = np.array([[int(mv[2][i, j]) for i in range(24)] for j in range(196560 // 2)], dtype=np.int64)
L4 = np.vstack([L4, -L4])
GL4 = L4 @ G                                   # rows: lambda^T G, so y . lambda = (lambda^T G) y



def invariant(y):
    ip = GL4 @ y
    vals_, counts = np.unique(ip, return_counts=True)
    return tuple(zip(vals_.tolist(), counts.tolist()))


# the E6^4 shallow multiple 3c (order 4) as reference: its unique minimal vector x0, y0 = 4 x0
h, n0, z2 = rows['e6^4']
zrho = lor2(z2, rho2) // 4
zrp = lor2(z2, rhop2) // 4
lz2 = [a + zrp * b + zrho * c for a, b, c in zip(z2, rho2, rhop2)]
sol = pari.matsolve(Gp, pari.matrix(24, 1, [lor2(a, lz2) // 4 for a in lam]))
Yh = np.array([int(sol[i, 0]) for i in range(24)], dtype=np.int64)        # coordinates of 12 c


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


Y3c = 3 * Yh                                   # 12 * (3c) = coordinates; Lambda + 3c <-> a + Y3c/12
sv, sa = shell_vecs(G.astype(np.float64), G, Y3c, 12, 2.0)
vals = sorted(zip([int(x) for x in sv], sa), key=lambda t: t[0])
best = vals[0][0]
mins = [12 * a + Y3c for s, a in vals if s == best]
check(best == 216 and len(mins) == 1, f'E6^4: Lambda + 3c has a unique vector of norm {best}/144 = 3/2')
y0 = mins[0] // 3                              # 4 x0 = (12 a + 36 c)/3
assert int(y0 @ G @ y0) == 24
inv0 = invariant(y0)
check(max(t[0] for t in inv0) <= 7, f'E6^4 reference y0 = 4 x0: max y0.mu over Lambda_4 = {max(t[0] for t in inv0)}; distribution {inv0}')

rng = np.random.default_rng(20261023)
seen = {}
fricke = {}
tried = 0
t0 = time.time()
while sum(seen.values()) < NS and tried < 2000 * NS:
    tried += 1
    i1, i2, i3, i4 = rng.integers(0, len(L4), 4)
    m1, m2 = L4[i1], L4[i2]
    if int(m1 @ G @ m2) != 0:
        continue
    u = m1 + m2                                  # norm 8
    m3, m4 = L4[i3], L4[i4]
    p34 = int(m3 @ G @ m4)
    if p34 not in (-1, 0):
        continue
    w = m3 + m4                                  # norm 6 (p34 = -1) or 8 (p34 = 0)
    nw = int(w @ G @ w)
    if int(u @ G @ w) != 4 - nw:                 # |u + 2w|^2 = 8 + 4 u.w + 4|w|^2 = 24
        continue
    y = u + 2 * w
    assert int(y @ G @ y) == 24
    inv = invariant(y)
    seen[inv] = seen.get(inv, 0) + 1
    if max(t[0] for t in inv) <= 7:
        fricke[inv] = fricke.get(inv, 0) + 1
print(f'       {sum(seen.values())} samples of norm-24 vectors in type-4 classes mod 2 Lambda ({time.time() - t0:.0f}s); '
      f'{len(seen)} invariants; Fricke (max y.mu <= 7): {sum(fricke.values())} samples, {len(fricke)} invariant(s)', flush=True)
for inv, c in sorted(fricke.items(), key=lambda t: -t[1]):
    print(f'       Fricke invariant (x{c}): {inv}', flush=True)
check(len(fricke) == 1 and list(fricke)[0] == inv0,
      'every sampled order-4 Fricke point has the invariant of the E6^4 shallow multiple 3c (= the A2^12 (h+1)-point)')
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
