#!/usr/bin/env python3
"""S37 audit 2, script b6c: direct check of Lemma B2 of collatz_amplitude_boundary_20261013c.md (the finite identity on
which every residue certificate rests), for t(n) = (-1)^tau(n) (a = b = -1, kappa = 0):
    -(3/2)^s C_N(s) - D(s) F_N(s) = -(F_N - F_L) - 2^-s (F_N - F_J),   C_N = E_N + R_N,
with E_N = sum_(n<=N) (1[n = 2 mod 3] - 1/3) t(n) n^-s, R_N = sum_(n<=N, n = 2 mod 3) t(n) [(n - 1/2)^-s - n^-s],
L = floor((2N-1)/3), J = floor(L/2), at random complex s and many N (mpmath, 30 digits).  Also the consequence used by
the certificates: at a zero rho of D, C_N(rho) = (2/3)^rho [(F_N - F_L) + 2^-rho (F_N - F_J)].
"""
import random
import mpmath as mp

mp.mp.dps = 30


def tau_par(n):
    s = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2
        s += 1
    return s & 1


NMAX = 3000
t = [0] + [1 - 2 * tau_par(n) for n in range(1, NMAX + 1)]
rng = random.Random(5)
worst = mp.mpf(0)
for trial in range(12):
    s = mp.mpc(rng.uniform(0.3, 2.0), rng.uniform(-60, 60))
    pw = [mp.mpf(0)] + [mp.power(n, -s) for n in range(1, NMAX + 1)]
    F = [mp.mpf(0)] * (NMAX + 1); E = [mp.mpf(0)] * (NMAX + 1); R = [mp.mpf(0)] * (NMAX + 1)
    for n in range(1, NMAX + 1):
        F[n] = F[n - 1] + t[n] * pw[n]
        E[n] = E[n - 1] + ((1 if n % 3 == 2 else 0) - mp.mpf(1) / 3) * t[n] * pw[n]
        R[n] = R[n - 1] + (t[n] * (mp.power(n - mp.mpf(1) / 2, -s) - pw[n]) if n % 3 == 2 else 0)
    Dv = 1 + mp.power(2, -s) + mp.power(mp.mpf(3) / 2, s) / 3
    for N in rng.sample(range(2, NMAX + 1), 40):
        L = (2 * N - 1) // 3; J = L // 2
        lhs = -mp.power(mp.mpf(3) / 2, s) * (E[N] + R[N]) - Dv * F[N]
        rhs = -(F[N] - F[L]) - mp.power(2, -s) * (F[N] - F[J])
        worst = max(worst, abs(lhs - rhs) / max(1, abs(rhs)))
print(f'Lemma B2: max relative |lhs - rhs| over 12 random s and 480 random N <= {NMAX}: {mp.nstr(worst, 3)}')
rho = mp.findroot(lambda z: 1 + mp.power(2, -z) + mp.power(mp.mpf(3) / 2, z) / 3, mp.mpc('0.6939', '22.2532'))
pw = [mp.mpf(0)] + [mp.power(n, -rho) for n in range(1, NMAX + 1)]
F = [mp.mpf(0)] * (NMAX + 1); C = [mp.mpf(0)] * (NMAX + 1)
for n in range(1, NMAX + 1):
    F[n] = F[n - 1] + t[n] * pw[n]
    C[n] = C[n - 1] + ((1 if n % 3 == 2 else 0) - mp.mpf(1) / 3) * t[n] * pw[n] + \
        (t[n] * (mp.power(n - mp.mpf(1) / 2, -rho) - pw[n]) if n % 3 == 2 else 0)
w = 0
for N in (100, 1000, 2999, 3000):
    L = (2 * N - 1) // 3; J = L // 2
    w = max(w, abs(C[N] - mp.power(mp.mpf(2) / 3, rho) * ((F[N] - F[L]) + mp.power(2, -rho) * (F[N] - F[J]))))
print(f'at rho = {mp.nstr(rho, 10)}: max |C_N - (2/3)^rho[(F_N - F_L) + 2^-rho (F_N - F_J)]| = {mp.nstr(w, 3)}')
