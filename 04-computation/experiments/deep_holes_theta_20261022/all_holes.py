"""HYP-9363 for all 23 deep-hole types of the Leech lattice (session opus-2026-10-22-S53).

Conjecture (HYP-9363). For a deep hole c of type R, write N for the Niemeier lattice, h for its Coxeter number,
rho_N for its Weyl vector and k for the number of components of R. Then

    theta_(Lambda + c) - theta_(N + rho_N/h) = (24 + k) Delta.

Coefficientwise this says, for every j >= 0,

    N_j - N'_j = (24 + k) tau(1 + j/h)  if h | j,  and  0  otherwise,

where
- N_j  = #{lambda in Lambda : |lambda - c|^2 = 2 + 2j/h}      (= #{Leech roots r : r . z = -j}),
- N'_j = #{x in N : |x - rho_N/h|^2 = 2 + 2j/h}               (= #{roots r : r . z = -1, r . rho = -j}).

Data. The cusp vectors z (primitive null vectors of II_{25,1} with z^perp/z = N and z . rho = -h) are taken from
Borcherds' table of norm-0 vectors (https://math.berkeley.edu/~reb/lattices/norm0). His coordinates are doubled, with
time first and the 25 space coordinates ordered so that rho = (70 | 24, ..., 0). They are converted here to the
tournament coordinates of THM-4752, where rho = (0, ..., 24 | 70). The table's columns 2-4 give the numbers of Leech
roots at levels 0, 1, 2, i.e. N_0 = 24 + k, N_1, N_2; they are re-derived here, not assumed.

Method.
- Leech side: Lambda = {rho, rho'}^perp with rho' = rho + e_0 - e_1, and c = (z + (z.rho') rho + (z.rho) rho')/h.
- N side: N = {z, z'}^perp with z' = y + (y^2/2) z for any y in II with y . z = -1, and
  rho_N = rho + (rho.z') z + (rho.z) z'.
- Both sides are enumerated exactly (Fincke-Pohst, integer norms) up to |.|^2 <= 2 + 2 J/h with J = JFAC * h.
Needs cypari and numba. Usage: python all_holes.py [JFAC=2.0] [types...]. Prints ALL CHECKS PASSED or the
disagreements."""
import sys, time
import numpy as np
from cypari import pari
from numba import njit

JFAC = float(sys.argv[1]) if len(sys.argv) > 1 else 2.0
ONLY = sys.argv[2:]
pari.allocatemem(2 * 10 ** 9)
T0 = time.time()
FAIL = []


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


TABLE = """\
  2 48 4096 97152  a1^24 244823040 68 24 22 22 20 20 18 18 16 16 14 14 12 12 10 10  8  8  6  6  4  4  2  2  0  0
  3 36 729 8748  a2^12 190080 48 16 16 16 14 14 14 12 12 10 10 10  8  8  8  6  6  6  4  4  4  2  2  2  0  0
  4 32 256 2048  a3^8 2688 36 12 12 12 10 10 10 10  8  8  8  8  6  6  6  6  4  4  4  2  2  2  2  0  0  0
  5 30 125 750  a4^6 240 30 10 10 10  8  8  8  8  8  6  6  6  6  6  4  4  4  4  4  2  2  2  2  2  0  0
  6 29 72 360  a5^4d4 48 24  8  8  8  8  6  6  6  6  6  6  4  4  4  4  4  2  2  2  2  2  2  0  0  0  0
  6 30 64 384  d4^6 2160 32 12 10 10 10 10  8  8  8  8  6  6  6  6  4  4  4  4  2  2  2  2  0  0  0  0
  7 28 49 196  a6^4 24 21  7  7  7  7  5  5  5  5  5  5  5  3  3  3  3  3  3  3  1  1  1  1  1  1 -1
  8 28 32 128  a7^2d5^2 8 18  6  6  6  6  6  4  4  4  4  4  4  4  2  2  2  2  2  2  2  2  0  0  0  0  0
  9 27 27 81  a8^3 12 17  5  5  5  5  5  5  5  5  3  3  3  3  3  3  3  3  3  1  1  1  1  1  1  1 -1
 10 27 20 60  a9^2d6 4 14  4  4  4  4  4  4  4  4  4  4  2  2  2  2  2  2  2  2  2  0  0  0  0  0  0
 10 28 16 64  d6^4 24 20  8  6  6  6  6  6  6  4  4  4  4  4  4  2  2  2  2  2  2  0  0  0  0  0  0
 12 27 12 36  a11d7e6 2 12  4  4  4  4  4  4  2  2  2  2  2  2  2  2  2  2  2  2  0  0  0  0  0  0  0
 12 28  9 36  e6^4 48 24  8  8  8  8  8  8  6  6  4  4  4  4  4  4  2  2  2  2  2  2  0  0  0  0  0
 13 26 13 26  a12^2 4 11  3  3  3  3  3  3  3  3  3  3  3  3  1  1  1  1  1  1  1  1  1  1  1  1 -1
 14 27  8 24  d8^3 6 14  6  4  4  4  4  4  4  4  4  2  2  2  2  2  2  2  2  0  0  0  0  0  0  0  0
 16 26  8 16  a15d9 2  8  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  2  0  0  0  0  0  0  0  0  0
 18 26  6 12  a17e7 2  9  3  3  3  3  3  3  3  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1 -1
 18 27  4 12  d10e7^2 2 12  4  4  4  4  4  4  4  2  2  2  2  2  2  2  2  0  0  0  0  0  0  0  0  0  0
 22 26  4  8  d12^2 2  8  4  2  2  2  2  2  2  2  2  2  2  2  2  0  0  0  0  0  0  0  0  0  0  0  0
 25 25  5  5  a24 2  5  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1  1
 30 26  2  4  d16e8 1  6  2  2  2  2  2  2  2  2  2  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 30 27  1  3  e8^3 6  18  6  6  6  6  6  6  6  6  2  2  2  2  2  2  2  2  2  0  0  0  0  0  0  0  0
 46 25  2  2  d24 1  2  2  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
"""
n = 25


def lor2(x, y):          # doubled coordinates -> 4 * (true inner product)
    return int(sum(int(a) * int(b) for a, b in zip(x[:-1], y[:-1])) - int(x[-1]) * int(y[-1]))


rows = []
for line in TABLE.strip('\n').split('\n'):
    t = line.split()
    h, n0, n1, n2, name, aut = int(t[0]), int(t[1]), int(t[2]), int(t[3]), t[4], int(t[5])
    co = [int(x) for x in t[6:]]
    assert len(co) == 26
    z2 = [co[1 + 24 - i] for i in range(n)] + [co[0]]          # doubled, my coordinates
    rows.append((name, h, n0, n1, n2, z2))

# II_{25,1} basis (doubled coordinates)
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
inII = lambda x2: (all(t % 2 == 0 for t in x2) or all(t % 2 == 1 for t in x2)) and ((sum(x2[:-1]) - x2[-1]) // 2) % 2 == 0
rho2 = [2 * i for i in range(n)] + [140]
rhop2 = list(rho2); rhop2[0] += 2; rhop2[1] -= 2


def sublattice_perp(u2, w2):
    F = pari.matrix(2, n + 1, [lor2(b, u2) // 4 for b in basis] + [lor2(b, w2) // 4 for b in basis])
    K = pari.matkerint(F)
    vecs = [[sum(int(K[i, j]) * basis[i][t] for i in range(n + 1)) for t in range(n + 1)] for j in range(24)]
    G = pari.matrix(24, 24, [lor2(a, b) // 4 for a in vecs for b in vecs])
    U = pari.qflllgram(G)
    vecs = [[sum(int(U[i, j]) * vecs[i][t] for i in range(24)) for t in range(n + 1)] for j in range(24)]
    Gm = [[lor2(a, b) // 4 for b in vecs] for a in vecs]
    return vecs, Gm


def coords_of(vecs, Gm, x2):       # integer coordinates of a lattice vector x (doubled coords) in the basis vecs
    Gp = pari.matrix(24, 24, [v for row in Gm for v in row])
    g = pari.matrix(24, 1, [lor2(a, x2) // 4 for a in vecs])
    sol = pari.matsolve(Gp, g)
    out = [sol[i, 0] for i in range(24)]
    assert all(pari.denominator(t) == 1 for t in out)
    return [int(t) for t in out]


def bezout_minus_one(vals):
    """Integer coefficients a with sum a_j vals_j = -1 (gcd of vals is 1 since z is primitive in a unimodular lattice)."""
    g, coef = 0, [0] * len(vals)
    for j, v in enumerate(vals):
        if v == 0:
            continue
        if g == 0:
            g, coef = v, [0] * len(vals)
            coef[j] = 1
            continue
        # extended gcd of g and v
        a0, b0, a1, b1, r0, r1 = 1, 0, 0, 1, g, v
        while r1:
            qq = r0 // r1
            r0, r1 = r1, r0 - qq * r1
            a0, a1 = a1, a0 - qq * a1
            b0, b1 = b1, b0 - qq * b1
        coef = [a0 * c for c in coef]
        coef[j] += b0
        g = r0
    assert abs(g) == 1
    return [-c * g for c in coef]


@njit(cache=True)
def enumerate_center(Gf, Gi, Yv, h, R2, JMAX):
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
    counts = np.zeros(JMAX + 1, dtype=np.int64)
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
                if kk <= JMAX:
                    counts[kk] += 1
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
    return counts, closer


def tau_list(m):
    pr = [1] + [0] * m
    for k in range(1, m + 1):
        for _ in range(24):
            pr = [pr[i] - (pr[i - k] if i >= k else 0) for i in range(m + 1)]
    return [0] + pr[:m]          # tau(1..m) at index 1..m


TAU = tau_list(12)
lam, GL = sublattice_perp(rho2, rhop2)
GLp = pari.matrix(24, 24, [v for row in GL for v in row])
check(int(pari.matdet(GLp)) == 1 and int(pari.qfminim(GLp, 2, 0)[0]) == 0, 'Leech lattice {rho, rho\'}^perp (det 1, rootless)')
GLf, GLi = np.array(GL, dtype=np.float64), np.array(GL, dtype=np.int64)

summary = []
for name, h, n0, n1, n2, z2 in rows:
    if ONLY and name not in ONLY:
        continue
    t1 = time.time()
    zrho = lor2(z2, rho2) // 4
    prim = not inII([t // 2 for t in z2]) if all(t % 2 == 0 for t in z2) else True
    ok_z = lor2(z2, z2) == 0 and inII(z2) and zrho == -h and prim
    # Leech side: c = (z + (z.rho') rho + (z.rho) rho') / h
    zrp = lor2(z2, rhop2) // 4
    lz2 = [a + zrp * b + zrho * c for a, b, c in zip(z2, rho2, rhop2)]          # = h c (doubled)
    Yc = coords_of(lam, GL, lz2)
    JMAX = int(round(JFAC * h))
    R2 = 2.0 + 2.0 * JMAX / h
    NL, closeL = enumerate_center(GLf, GLi, np.array(Yc, dtype=np.int64), h, R2, JMAX)
    # N side: y with y.z = -1 via extended gcd on the basis pairings, z' = y + (y^2/2) z
    pz = [lor2(b, z2) // 4 for b in basis]
    acoef = bezout_minus_one(pz)                                   # sum a_j pz_j = -1
    y2 = [sum(acoef[j] * basis[j][t] for j in range(n + 1)) for t in range(n + 1)]
    assert lor2(y2, z2) // 4 == -1
    yy = lor2(y2, y2) // 4
    zp2 = [a + (yy // 2) * b for a, b in zip(y2, z2)]
    assert lor2(zp2, zp2) == 0 and lor2(zp2, z2) == -4
    Nv, GN = sublattice_perp(z2, zp2)
    GNp = pari.matrix(24, 24, [v for row in GN for v in row])
    roots_N = int(pari.qfminim(GNp, 2, 0)[0])            # pairs +-
    rzp = lor2(rho2, zp2) // 4
    rhoN2 = [r + rzp * a + zrho * b for r, a, b in zip(rho2, z2, zp2)]            # rho_N (doubled)
    Yr = coords_of(Nv, GN, rhoN2)
    NN, closeN = enumerate_center(np.array(GN, dtype=np.float64), np.array(GN, dtype=np.int64), np.array(Yr, dtype=np.int64), h, R2, JMAX)
    # compare
    diffs = []
    good = True
    for j in range(JMAX + 1):
        pred = (n0 * TAU[1 + j // h]) if j % h == 0 else 0
        dj = int(NL[j]) - int(NN[j])
        diffs.append(dj)
        good &= dj == pred
    ok_tab = int(NL[0]) == n0 and (JMAX < 1 or int(NL[1]) == n1) and (JMAX < 2 or int(NL[2]) == n2)
    ok_struct = closeL == 0 and closeN == 0 and int(NN[0]) == 0 and roots_N == 24 * h and int(pari.matdet(GNp)) == 1
    check(ok_z and ok_tab and ok_struct and good,
          f'{name:9s} h={h:2d}: N_0..N_2 = {list(NL[:3])} (table {n0}, {n1}, {n2}); |R(N)| = {roots_N} = 24h; '
          f'N_j - N\'_j = {n0} tau(1 + j/h) [h | j], 0 otherwise, for j <= {JMAX} [{time.time() - t1:.0f}s]')
    if not good:
        print('       j, N_j, N\'_j, diff:', [(j, int(NL[j]), int(NN[j]), diffs[j]) for j in range(JMAX + 1) if diffs[j] != ((n0 * TAU[1 + j // h]) if j % h == 0 else 0)][:12])
    summary.append((name, h, n0, good))

print()
print('summary:', [(nm, h, ok) for nm, h, n0, ok in summary])
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f' ({time.time() - T0:.0f}s)')
