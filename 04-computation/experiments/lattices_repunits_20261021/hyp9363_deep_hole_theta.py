"""HYP-9363 test: the theta series of the Leech lattice around its A_24 deep hole (session opus-2026-10-21-S52).

In the tournament coordinates of THM-4752 (II_{25,1}, rho = (0,...,24 | 70), nu = (1/2^25 | 5/2)) the Leech roots
(r^2 = 2, r . rho = -1) are r_lambda = lambda + rho' + (lambda^2/2 - 1) rho, lambda in Lambda = {rho, rho'}^perp, with
rho' = rho + e_0 - e_1 (null, rho . rho' = -1).  Then r_lambda . nu = 25 (1 - |lambda - c|^2 / 2) with
c = nu/25 - rho - rho', so N_k = #{Leech roots with r . nu = -k} = #{lambda : |lambda - c|^2 = 2 + 2k/25}: the theta
series of the coset Lambda - c around the A_24 deep hole c is  theta(tau) = sum_k N_k q^(1 + k/25).
The independent audit (audit_S52, A5) found N_0 = 25 and N_k = 5 p(k-1) for 1 <= k <= 22 (p = partitions).
Pre-registered candidate (written before any count beyond k = 22 was computed):
    theta(tau) = 25 eta(tau)^24 + 5 eta(tau)^25 / eta(tau/25)
(both terms of weight 12; in Q = q^(1/25): 25 Q^25 prod(1 - Q^25n)^24 + 5 Q^26 prod(1 - Q^25n)^25 / prod(1 - Q^n)),
which agrees with 5 p(k-1) for k <= 24 and predicts the first departures N_25 = 5 p(24) - 600 = 7275,
N_26 = 5 (p(25) - 25) = 12055.
Here Lambda is built from scratch (PARI HNF / integer kernel / LLL; checked even, unimodular, rootless = Leech) and the
points of Lambda near c are enumerated exactly (Fincke-Pohst with exact integer norms) up to k = KMAX.
Needs cypari and numba.  Prints ALL CHECKS PASSED or the first disagreement."""
import sys, time
from fractions import Fraction
import numpy as np
from cypari import pari
from numba import njit

KMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 40
pari.allocatemem(2 * 10 ** 9)
T0 = time.time()
FAIL = []


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


n = 25
J = np.diag([1] * n + [-1])


def lor2(x, y):          # x, y in doubled coordinates (integers); returns the true inner product * 4
    return int(sum(int(a) * int(b) for a, b in zip(x[:-1], y[:-1])) - int(x[-1]) * int(y[-1]))


# ---------------- II_{25,1} in doubled coordinates
gens = []
for i in range(n - 1):
    v = [0] * (n + 1); v[i] = 2; v[i + 1] = -2; gens.append(v)
v = [0] * (n + 1); v[n - 2] = 2; v[n - 1] = 2; gens.append(v)            # e_23 + e_24
gens.append([1] * (n + 1))                                                 # (1/2)^26
v = [0] * (n + 1); v[n - 1] = 2; v[n] = 2; gens.append(v)                  # e_24 + e_25
v = [0] * (n + 1); v[0] = 4; gens.append(v)                                # 2 e_0
A = pari.matrix(n + 1, len(gens), [gens[j][i] for i in range(n + 1) for j in range(len(gens))])
Hb = pari.mathnf(A)                                                        # columns: a basis (doubled coords)
Hm = [[int(Hb[i, j]) for j in range(n + 1)] for i in range(n + 1)]
basis = [[Hm[i][j] for i in range(n + 1)] for j in range(n + 1)]
G = [[Fraction(lor2(a, b), 4) for b in basis] for a in basis]
detG = pari.matdet(pari.matrix(n + 1, n + 1, [int(G[i][j]) for i in range(n + 1) for j in range(n + 1)]))
inII = lambda x2: (all(t % 2 == 0 for t in x2) or all(t % 2 == 1 for t in x2)) and ((sum(x2[:-1]) - x2[-1]) // 2) % 2 == 0
check(all(g.denominator == 1 for row in G for g in row) and all(G[i][i] % 2 == 0 for i in range(n + 1))
      and int(detG) == -1 and all(inII(b) for b in basis),
      'II_{25,1}: an HNF basis of 26 lattice vectors with even integral Gram matrix of determinant -1')

rho2 = [2 * i for i in range(n)] + [140]
rhop2 = list(rho2); rhop2[0] += 2; rhop2[1] -= 2                         # rho' = rho + e_0 - e_1
nu2 = [1] * n + [5]
check(lor2(rho2, rho2) == 0 and lor2(rhop2, rhop2) == 0 and lor2(rho2, rhop2) == -4 and inII(rhop2),
      "rho' = rho + e_0 - e_1 is null with rho . rho' = -1 (II = Lambda + U)")

# ---------------- Lambda = {rho, rho'}^perp: integer kernel, Gram, LLL
F = pari.matrix(2, n + 1, [lor2(b, rho2) // 4 for b in basis] + [lor2(b, rhop2) // 4 for b in basis])
K = pari.matkerint(F)
Kc = [[int(K[i, j]) for j in range(24)] for i in range(n + 1)]
lam = [[sum(Kc[i][j] * basis[i][t] for i in range(n + 1)) for t in range(n + 1)] for j in range(24)]
GL = pari.matrix(24, 24, [lor2(a, b) // 4 for a in lam for b in lam])
U = pari.qflllgram(GL)
Uc = [[int(U[i, j]) for j in range(24)] for i in range(24)]
lam = [[sum(Uc[i][j] * lam[i][t] for i in range(24)) for t in range(n + 1)] for j in range(24)]
Gm = [[lor2(a, b) // 4 for b in lam] for a in lam]
Gp = pari.matrix(24, 24, [x for row in Gm for x in row])
short = pari.qfminim(Gp, 2, 0)
check(all(lor2(a, b) % 4 == 0 for a in lam for b in lam) and all(Gm[i][i] % 2 == 0 for i in range(24))
      and int(pari.matdet(Gp)) == 1 and int(short[0]) == 0 and all(lor2(a, rho2) == 0 and lor2(a, rhop2) == 0 for a in lam),
      'Lambda = {rho, rho\'}^perp: even, unimodular, rank 24, no vectors of norm 2: the Leech lattice')

# ---------------- the deep hole c = nu/25 - rho - rho' in Lambda (x) Q
# lambda_nu = 25 c = nu - 25 rho - 25 rho'  (doubled coordinates)
lnu2 = [a - 25 * b - 25 * c_ for a, b, c_ in zip(nu2, rho2, rhop2)]
check(lor2(lnu2, rho2) == 0 and lor2(lnu2, rhop2) == 0 and lor2(lnu2, lnu2) == 4 * 1250,
      'lambda_nu = 25 c = nu - 25 rho - 25 rho\' lies in Lambda, |25 c|^2 = 1250, so |c|^2 = 2 (lambda = 0 is a vertex of the hole)')
g = pari.matrix(24, 1, [lor2(a, lnu2) // 4 for a in lam])
Ycol = pari.matsolve(Gp, g)                        # coordinates of 25 c in the LLL basis
Y = [int(Ycol[i, 0]) for i in range(24)]
check(all(Ycol[i, 0] == Y[i] for i in range(24)) and
      [sum(Y[j] * lam[j][t] for j in range(24)) for t in range(n + 1)] == lnu2,
      '25 c has integral coordinates in the Leech basis')

# ---------------- exact enumeration of lambda with |lambda - c|^2 <= 2 + 2 KMAX/25
Gf = np.array(Gm, dtype=np.float64)
Gi = np.array(Gm, dtype=np.int64)
Yv = np.array(Y, dtype=np.int64)
R2 = 2.0 + 2.0 * KMAX / 25.0


@njit(cache=True)
def enumerate_hole(Gf, Gi, Yv, R2, KMAX):
    d = Gf.shape[0]
    # Fincke-Pohst form Q(x) = sum_i q[i,i] (x_i + sum_{j>i} q[i,j] x_j)^2
    q = Gf.copy()
    for i in range(d):
        for j in range(i + 1, d):
            q[j, i] = q[i, j]
            q[i, j] = q[i, j] / q[i, i]
        for k in range(i + 1, d):
            for l in range(k, d):
                q[k, l] -= q[k, i] * q[i, l]
    y = Yv.astype(np.float64) / 25.0
    counts = np.zeros(KMAX + 1, dtype=np.int64)
    a = np.zeros(d, dtype=np.int64)
    T = np.zeros(d + 1)
    cen = np.zeros(d)
    hi = np.zeros(d, dtype=np.int64)
    i = d - 1
    # initialise level d-1
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
            # exact norm: 625 |lambda - c|^2 = (25 a - Y)^T G (25 a - Y)
            s = 0
            for u in range(d):
                zu = 25 * a[u] - Yv[u]
                if zu == 0:
                    continue
                acc = 0
                for w in range(d):
                    acc += Gi[u, w] * (25 * a[w] - Yv[w])
                s += zu * acc
            if (s - 1250) % 50 == 0:
                kk = (s - 1250) // 50
                if 0 <= kk <= KMAX:
                    counts[kk] += 1
            total += 1
            a[i] += 1
            continue
        # descend
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
    return counts, total


counts, total = enumerate_hole(Gf, Gi, Yv, R2, KMAX)
print(f'     enumerated {total} Leech points with |lambda - c|^2 <= {R2:.2f} [{time.time() - T0:.0f}s]')

# partitions
P = [1] + [0] * (KMAX + 60)
for part in range(1, len(P)):
    for m in range(part, len(P)):
        P[m] += P[m - part]
# candidate: 25 prod(1 - Q^25n)^24  +  5 Q prod(1 - Q^25n)^25 / prod(1 - Q^n)
L = KMAX + 1
A24 = [0] * L; A24[0] = 1
for nn in range(1, L // 25 + 2):
    for _ in range(24):
        A24 = [A24[m] - (A24[m - 25 * nn] if m >= 25 * nn else 0) for m in range(L)]
A25 = [0] * L; A25[0] = 1
for nn in range(1, L // 25 + 2):
    for _ in range(25):
        A25 = [A25[m] - (A25[m - 25 * nn] if m >= 25 * nn else 0) for m in range(L)]
B = [sum(A25[i] * P[m - i] for i in range(m + 1)) for m in range(L)]          # prod(1-Q^25n)^25 / prod(1-Q^n)
cand = [25 * A24[k] + (5 * B[k - 1] if k >= 1 else 0) for k in range(L)]
naive = [25] + [5 * P[k - 1] for k in range(1, L)]
print('      k   N_k (exact)   5 p(k-1)   candidate')
for k in range(L):
    print(f'     {k:3d}  {int(counts[k]):11d}  {naive[k]:9d}  {cand[k]:10d}' + ('' if counts[k] == cand[k] else '   <-- differs'))
check(list(counts[:23]) == naive[:23], 'N_k = 25, 5 p(k-1) for k <= 22 (reproduces the audit enumeration of Leech roots)')
first_naive = next((k for k in range(L) if counts[k] != naive[k]), None)
print(f'     first k with N_k != 5 p(k-1): {first_naive}')
check(all(int(counts[k]) == cand[k] for k in range(L)),
      f'theta_(Lambda - c) = 25 eta(tau)^24 + 5 eta(tau)^25/eta(tau/25) through k = {KMAX} (q^(1 + {KMAX}/25))')
# sanity: the Leech-root correspondence on a few explicit vectors
print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f' ({time.time() - T0:.0f}s)')
