#!/usr/bin/env python3
"""S37 audit 2, script b2: THM-4658 statement 6 (character inheritance) and the remainder identity of HYP-9285.

Exact census t(n) for n <= 2^22 + 2 (simple sequential recursion: t(2m) = -t(m); odd n iterated until below n).
 (1) Statement 6 for every Dirichlet character psi mod 9 and mod 27 (lambda = psi(2)), the twisted splits
     phi_j(n) = [3^j || n] chi(n/3^j) (j = 0, 1, 2, chi the character mod 3), and the untwisted multiples of 3^j:
     C_phi(N) + lambda C_phi(floor(N/2)) - O_phi(N) = 0 and O_phi(N) = -sum_(m <= 3 ceil(N/2), m = 2 (3)) phi((2m-1)/3) t(m)
     for every N <= 2^20.  Negative control: phi = 1_(n = 1 mod 9) (not a halving eigenfunction).
 (2) Octave-averaged ratio <D_j/sqrt N>/<O_j/sqrt N> on the grid of 64 points per octave, 2^10..2^20, versus 2 + sqrt 2.
 (3) H(N) = -[S(N) + 3S(floor(2N/3)) + 3S(floor(N/3))] for every N <= 2^22; the same identity for a halving-null
     sequence (random odd signs, t(2m) = -t(m)) as a control; the algebra 1 + 3(2/3)^r + 3(1/3)^r = 3(2/3)^r D(r).
"""
import math, cmath
import numpy as np
import numba as nb
import mpmath as mp

NMAX = (1 << 22) + 8


@nb.njit
def census(nmax):
    t = np.zeros(nmax + 1, dtype=np.int64)
    t[1] = 1
    for n in range(2, nmax + 1):
        if n % 2 == 0:
            t[n] = -t[n // 2]
        else:
            x = n
            s = 0
            while x >= n:
                x = x // 2 if x % 2 == 0 else (3 * x + 1) // 2
                s += 1
            t[n] = t[x] if s % 2 == 0 else -t[x]
    return t


t = census(NMAX)
# spot check against full trajectories
def tau_par(n):
    s = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2
        s += 1
    return s & 1
rng = np.random.default_rng(7)
assert all(t[n] == 1 - 2 * tau_par(int(n)) for n in rng.integers(1, NMAX, 3000))
print(f'census t(n), n <= {NMAX}: S(2^20) = {int(t[1:2**20 + 1].sum())}')

NN = 1 << 20
n_all = np.arange(NMAX + 1)


def dlog_table(q):
    m = 3 ** q
    dl = -np.ones(m, dtype=np.int64)
    x = 1
    for j in range(2 * 3 ** (q - 1)):
        dl[x] = j
        x = (2 * x) % m
    return dl


def check(phi, lam, name):
    # phi: complex array over n = 0..NMAX
    f = phi * t
    C = np.cumsum(f)                       # C[N] = sum_(n<=N) phi(n) t(n)
    odd = np.where(n_all % 2 == 1, f, 0)
    O = np.cumsum(odd)
    N = np.arange(1, NN + 1)
    e1 = np.abs(C[N] + lam * C[N // 2] - O[N]).max()
    # second expression: m = 2 mod 3, m <= 3 ceil(N/2)
    g = np.zeros(NMAX + 1, dtype=complex)
    m = np.arange(2, NMAX + 1, 3)
    g[m] = phi[(2 * m - 1) // 3] * t[m]
    Gc = np.cumsum(g)
    up = 3 * ((N + 1) // 2)
    e2 = np.abs(O[N] + Gc[up]).max()
    return e1, e2, C, O


print('(1) statement 6: max over N <= 2^20 of |C(N) + lambda C(N/2) - O(N)| and |O(N) - second expression|')
for q in (2, 3):
    m = 3 ** q
    dl = dlog_table(q)
    order = 2 * 3 ** (q - 1)
    for k in range(order):
        phi = np.zeros(NMAX + 1, dtype=complex)
        r = n_all % m
        unit = dl[r] >= 0
        phi[unit] = np.exp(2j * np.pi * k * dl[r[unit]] / order)
        lam = cmath.exp(2j * math.pi * k / order)        # psi(2)
        e1, e2, _, _ = check(phi, lam, f'psi_{k} mod {m}')
        if k in (0, 1, order // 2) or e1 > 1e-6 or e2 > 1e-6:
            print(f'  psi_{k} mod {m}: psi(2) = {lam:.4f}; errors {e1:.1e}, {e2:.1e}')
    print(f'  all {order} characters mod {m} checked')
chi = np.array([0, 1, -1])
res2 = {}
for j in range(3):
    phi = np.zeros(NMAX + 1, dtype=complex)
    v = n_all.copy()
    sel = np.zeros(NMAX + 1, dtype=bool)
    sel[1:] = (n_all[1:] % 3 ** j == 0) & ((n_all[1:] // 3 ** j) % 3 != 0)
    phi[sel] = chi[(n_all[sel] // 3 ** j) % 3]
    e1, e2, C, O = check(phi, -1, f'phi_{j}')
    print(f'  twisted split phi_{j} (lambda = -1): errors {e1:.1e}, {e2:.1e}')
    res2[j] = (C, O)
    phi_u = np.zeros(NMAX + 1, dtype=complex)
    phi_u[sel] = 1.0
    e1, e2, _, _ = check(phi_u, 1, f'untwisted_{j}')
    print(f'  untwisted [3^{j} || n] (lambda = +1): errors {e1:.1e}, {e2:.1e}')
phi_bad = np.where(n_all % 9 == 1, 1.0 + 0j, 0)
e1, e2, _, _ = check(phi_bad, 1, 'bad')
e1m, _, _, _ = check(phi_bad, -1, 'bad')
print(f'  negative control 1_(n = 1 mod 9): |C + C(N/2) - O| max {e1:.0f}, |C - C(N/2) - O| max {e1m:.0f} (second expression {e2:.1e})')

print('(2) octave-averaged ratio <D_j>/<O_j> (64 log-uniform points per octave, 2^10..2^20)')
grid = np.array(sorted(set(int(round(2 ** (10 + i / 64))) for i in range(0, 10 * 64))))
for j in range(3):
    C, O = res2[j]
    d = (C[grid].real / np.sqrt(grid)).mean()
    o = (O[grid].real / np.sqrt(grid)).mean()
    print(f'  j = {j}: <D_j> = {d:+.4f}, <O_j> = {o:+.4f}, ratio {d / o:.3f} (2 + sqrt 2 = {2 + math.sqrt(2):.3f})')

print('(3) remainder identity H(N) = -[S(N) + 3S(floor(2N/3)) + 3S(floor(N/3))]')
w = np.where(n_all % 3 == 2, 2, -1)
w[0] = 0
S = np.cumsum(t)
H = np.cumsum(w * t)
N = np.arange(1, 1 << 22)
err = np.abs(H[N] + S[N] + 3 * S[(2 * N) // 3] + 3 * S[N // 3]).max()
print(f'  Collatz: max |H(N) + S(N) + 3S(2N/3) + 3S(N/3)| over N < 2^22 = {err}')
# halving null control: t0(n) = (-1)^v2(n) s(odd part), s random
rs = np.random.default_rng(11).integers(0, 2, NMAX + 1) * 2 - 1
tn = np.zeros(NMAX + 1, dtype=np.int64)
for n in range(1, NMAX + 1):
    v, mm = 0, n
    while mm % 2 == 0:
        mm //= 2; v += 1
    tn[n] = rs[mm] * (-1) ** v
Sn = np.cumsum(tn); Hn = np.cumsum(w * tn)
errn = np.abs(Hn[N] + Sn[N] + 3 * Sn[(2 * N) // 3] + 3 * Sn[N // 3])
print(f'  halving null: identity fails at {int((errn > 0).sum())} of {len(N)} N (max deviation {int(errn.max())}); '
      f'it needs the odd-step relation t(2k-1) = -t(3k-1)')
mp.mp.dps = 40
worst = 0
for r in (mp.mpc('0.99820922', '240.2478'), mp.mpc('0.3', '-7'), mp.mpc('2.5', '0.1')):
    Dv = 1 + mp.power(2, -r) + mp.power(mp.mpf(3) / 2, r) / 3
    lhs = 1 + 3 * mp.power(mp.mpf(2) / 3, r) + 3 * mp.power(mp.mpf(1) / 3, r)
    worst = max(worst, abs(lhs - 3 * mp.power(mp.mpf(2) / 3, r) * Dv))
print(f'  algebra 1 + 3(2/3)^r + 3(1/3)^r - 3(2/3)^r D(r): max {mp.nstr(worst, 3)} at three test points')
