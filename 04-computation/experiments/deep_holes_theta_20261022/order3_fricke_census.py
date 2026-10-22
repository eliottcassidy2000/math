"""Which order-3 points of the Leech lattice twine V_Lambda by T_3A + 24? (session opus-2026-10-22-S53; HYP-9368)

For v of order 3 modulo Lambda, the trace F_v = Theta_Lambda^[v]/eta^24 is a modular function for Gamma_0(3). Its
poles are at infinity and at 0; the pole at 0 comes from the vectors of Lambda + v of norm < 2. By the argument in
THM-4759 (7), F_v - 24 is the Gamma_0(3)+ Hauptmodul T_3A exactly when Lambda + v has a unique vector of norm 4/3.
Writing y = 3x for that vector, y is a norm-12 vector of Lambda that is the unique minimum of its class
y + 3 Lambda. A norm-12 vector y is such a minimum iff no lambda in Lambda_4 has y.lambda = 6.

HYP-9368 predicts that all such y form a single Co_0-orbit, the class of the shallow multiple 2c of the D4^6 deep
hole. This script samples norm-12 vectors of Lambda (Fincke-Pohst from several randomised bases). For each sample it
records k(y) = #{lambda in Lambda_4 : y.lambda = 6} and the full distribution of y.lambda over Lambda_4, an orbit
invariant. It reports the distinct invariants among the k = 0 vectors, and checks that the D4^6 vector has k = 0.
Needs cypari and numba. Usage: python order3_fricke_census.py [SAMPLES_PER_BASIS=400] [BASES=8]
Result: two orbit invariants, k = 2 (sums of two norm-4 vectors with product 2, stabiliser U6(2):2) and k = 0, the
invariant of the M24-fixed vector (2^24)/sqrt 8. So the D4^6 shallow multiple 2c is the M24-invariant point
(2^24)/(3 sqrt 8) modulo Lambda, and it is the only order-3 point twining V_Lambda by T_3A + 24 (given the classical
two-orbit decomposition of the type-6 vectors)."""
import sys, time, math
import numpy as np
from cypari import pari
from numba import njit

SPB = int(sys.argv[1]) if len(sys.argv) > 1 else 400
NB = int(sys.argv[2]) if len(sys.argv) > 2 else 8
pari.allocatemem(2 * 10 ** 9)
FAIL = []


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


src = open(__file__.replace('order3_fricke_census.py', 'all_holes.py'), encoding='utf-8').read()
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

# the D4^6 shallow vector: y = 3 x0 with x0 the unique norm-4/3 vector of Lambda + 2c
h, n0, z2 = rows['d4^6']
zrho = lor2(z2, rho2) // 4
zrp = lor2(z2, rhop2) // 4
lz2 = [a + zrp * b + zrho * c for a, b, c in zip(z2, rho2, rhop2)]
sol = pari.matsolve(Gp, pari.matrix(24, 1, [lor2(a, lz2) // 4 for a in lam]))
Yh = np.array([int(sol[i, 0]) for i in range(24)], dtype=np.int64)        # coordinates of 6 c


@njit(cache=True)
def shell_vecs(Gf, Gi, Yd, den, R2):
    """all a in Z^24 with |a + Yd/den|^2 <= R2; returns den^2 |a + Yd/den|^2 and a."""
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


Y2c = 2 * Yh                                   # coordinates of 6 * (2c) = 12 c; x in Lambda + 2c <-> a + Y2c/6
sv, sa = shell_vecs(G.astype(np.float64), G, Y2c, 6, 2.0)
vals = sorted(zip([int(x) for x in sv], sa), key=lambda t: t[0])
best = vals[0][0]
ymin = [6 * a + Y2c for s, a in vals if s == best]
check(best * 9 == 12 * 36 and len(ymin) == 1, f'D4^6: Lambda + 2c has a unique vector of norm {best}/36 = 4/3')
y0 = ymin[0] // 2                              # 3 x0 in Lambda coordinates (6 a + 12 c)/2 = 3 (a + 2c)
assert int(y0 @ G @ y0) == 12


def invariant(y):
    ip = GL4 @ y
    vals_, counts = np.unique(ip, return_counts=True)
    return tuple(zip(vals_.tolist(), counts.tolist()))


inv0 = invariant(y0)
k0 = dict(inv0).get(6, 0)
check(k0 == 0, f'D4^6 shallow vector y0 = 3 x0: k(y0) = #{{lambda in Lambda_4 : y0.lambda = 6}} = {k0}; '
               f'inner-product distribution {inv0}')


@njit(cache=True)
def sample_norm(Q, target, maxn):
    """first maxn vectors a with a^T Q a == target (Fincke-Pohst, Q integral PSD)."""
    d = Q.shape[0]
    q = Q.astype(np.float64).copy()
    for i in range(d):
        for j in range(i + 1, d):
            q[j, i] = q[i, j]
            q[i, j] = q[i, j] / q[i, i]
        for k in range(i + 1, d):
            for l in range(k, d):
                q[k, l] -= q[k, i] * q[i, l]
    out = np.zeros((maxn, d), dtype=np.int64)
    cnt = 0
    a = np.zeros(d, dtype=np.int64)
    T = np.zeros(d + 1)
    cen = np.zeros(d)
    hi = np.zeros(d, dtype=np.int64)
    Bnd = target + 1e-6
    i = d - 1
    rad = np.sqrt(Bnd / q[i, i])
    a[i] = np.int64(np.ceil(-rad))
    hi[i] = np.int64(np.floor(rad))
    while cnt < maxn:
        if a[i] > hi[i]:
            i += 1
            if i == d:
                break
            a[i] += 1
            continue
        t = a[i] - cen[i]
        T[i] = T[i + 1] + q[i, i] * t * t
        if T[i] > Bnd:
            a[i] += 1
            continue
        if i == 0:
            s = 0
            for u in range(d):
                acc = 0
                for w in range(d):
                    acc += Q[u, w] * a[w]
                s += a[u] * acc
            if s == target:
                out[cnt] = a
                cnt += 1
            a[i] += 1
            continue
        j = i - 1
        sacc = 0.0
        for l in range(j + 1, d):
            sacc += q[j, l] * a[l]
        cen[j] = -sacc
        rad = np.sqrt(max(Bnd - T[i], 0.0) / q[j, j])
        a[j] = np.int64(np.ceil(cen[j] - rad))
        hi[j] = np.int64(np.floor(cen[j] + rad))
        T[j + 1] = T[i]
        i = j
    return out[:cnt]


rng = np.random.default_rng(20261022)
seen = {}
total = 0
t0 = time.time()
# (i) Fincke-Pohst samples from randomised bases (biased towards short sums of basis vectors)
for b in range(NB):
    Mx = np.eye(24, dtype=np.int64)
    for _ in range(60):
        i, j = rng.choice(24, 2, replace=False)
        Mx[:, i] += int(rng.integers(-2, 3)) * Mx[:, j]
    Qb = Mx.T @ G @ Mx
    Ul = pari.qflllgram(pari.matrix(24, 24, [int(x) for x in Qb.flatten()]))
    Uln = np.array([[int(Ul[i, j]) for j in range(24)] for i in range(24)], dtype=np.int64)
    Mb = Mx @ Uln
    Qb = Mb.T @ G @ Mb
    S = sample_norm(Qb, 12, SPB)
    for a in S:
        inv = invariant(Mb @ a)
        seen[inv] = seen.get(inv, 0) + 1
        total += 1
# (ii) random sums lambda_1 + lambda_2 + lambda_3 of norm-4 vectors with norm 12 (pairwise products summing to 0)
trials = 0
while trials < 40 * SPB:
    trials += 1
    i1, i2, i3 = rng.integers(0, len(L4), 3)
    y = L4[i1] + L4[i2] + L4[i3]
    if int(y @ G @ y) != 12:
        continue
    inv = invariant(y)
    seen[inv] = seen.get(inv, 0) + 1
    total += 1
print(f'       sampled {total} norm-12 vectors from {NB} randomised bases in {time.time() - t0:.0f}s; '
      f'{len(seen)} distinct inner-product distributions over Lambda_4 (orbit invariants)', flush=True)
for inv, c in sorted(seen.items(), key=lambda t: -t[1]):
    k = dict(inv).get(6, 0)
    print(f'       k = {k:5d}  (x{c:5d})  distribution {inv}', flush=True)
unique_classes = [inv for inv in seen if dict(inv).get(6, 0) == 0]
check(len(unique_classes) == 1 and unique_classes[0] == inv0 and len(seen) == 2,
      f'two invariants among sampled norm-12 vectors (k = 0 and k = 2), and every k = 0 sample has the D4^6 invariant')
# orbit sizes: Co_0 has two orbits on norm-12 vectors with stabilisers M24 and U6(2):2; their sizes add up exactly
co0, m24, u62 = 8315553613086720000, 244823040, 9196830720
check(co0 % m24 == 0 and co0 % (2 * u62) == 0 and co0 // m24 + co0 // (2 * u62) == 34417656000,
      f'|Co_0|/|M24| + |Co_0|/|U6(2):2| = {co0 // m24} + {co0 // (2 * u62)} = 34417656000 = number of norm-12 vectors')
# the M24-fixed vector (2^24)/sqrt 8 has the same Lambda_4 inner-product distribution as y0 (hand count in HYP-9368)
hand = ((-5, 24), (-4, 759), (-3, 6072), (-2, 21528), (-1, 43056), (0, 53682), (1, 43056), (2, 21528), (3, 6072), (4, 759), (5, 24))
check(inv0 == hand, 'the D4^6 shallow vector has the inner-product distribution of the M24-fixed vector (2^24)/sqrt 8')
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}')
