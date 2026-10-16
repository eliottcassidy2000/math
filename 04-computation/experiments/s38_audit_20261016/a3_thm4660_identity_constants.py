#!/usr/bin/env python3
"""S38 audit, item 2: THM-4660 at the level of finite sums and its explicit constants.

(1) The exact cut-off identity behind F D = G (Lemma B2 of collatz_amplitude_boundary_20261013c.md with a = b = -1,
    kappa = 0): for every s and N,
        -(3/2)^s [F_H,N/3 + E_N] - D(s) F_N = -(F_N - F_L) - 2^(-s) (F_N - F_J),
    F_N = sum_(n<=N) t(n) n^(-s), F_H,N = sum_(n<=N) w_n t(n) n^(-s), E_N = sum_(m<=N, m=2 (3)) t(m)[(m - 1/2)^(-s) - m^(-s)],
    L = floor((2N-1)/3), J = floor(L/2); exact t(n) for n <= 2^22 (own descent code), checked at eight s including two
    zeros of D.  Also the even/odd split itself (t(2m) = -t(m), t(2m-1) = -t(3m-1) for all m <= 2^21).
(2) Lemma C constants: 1 - f(3/2) on Re s = 3/2, max |D| on the Jensen disk, the Jensen ratio, the rectangle corners;
    the enlarged rectangle needed for the good-height step.
(3) The Riesz kernel: (1/2 pi i) int_(2) y^s 2 ds/(s(s+1)(s+2)) = (1 - 1/y)^2 [y >= 1] at a few y (mpmath quadrature).
Usage: python a3_thm4660_identity_constants.py
"""
import os, json, math
import numpy as np
import numba as nb
import mpmath as mp

HERE = os.path.dirname(os.path.abspath(__file__))
out = {}
NN = 1 << 22


@nb.njit
def parities(nmax):
    p = np.zeros(nmax + 1, np.int8)
    for n in range(2, nmax + 1):
        if n % 2 == 0:
            p[n] = 1 - p[n // 2]
        else:
            x = n
            s = 0
            while x >= n:
                x = x // 2 if x % 2 == 0 else (3 * x + 1) // 2
                s += 1
            p[n] = (s + p[x]) & 1
    return p


par = parities(NN)
t = (1 - 2 * par.astype(np.int64))
t[0] = 0
n = np.arange(NN + 1, dtype=np.float64)
m = np.arange(1, NN // 2 + 1)
ok_even = bool(np.all(t[2 * m] == -t[m]))
mo = np.arange(1, (NN + 1) // 3)
ok_odd = bool(np.all(t[2 * mo - 1] == -t[3 * mo - 1]))
out['t(2m) = -t(m) for m <= 2^21'] = ok_even
out['t(2m-1) = -t(3m-1) for 3m-1 <= 2^22'] = ok_odd
w = np.where(np.arange(NN + 1) % 3 == 2, 2, -1)
L2, L15 = math.log(2), math.log(1.5)


def D(s):
    return 1 + 2.0 ** (-s) + 1.5 ** s / 3


def Fsum(coef, s, N):
    k = np.arange(1, N + 1)
    return np.sum(coef[1:N + 1] * np.exp(-s * np.log(k.astype(np.float64))))


def Esum(s, N):
    k = np.arange(2, N + 1, 3)
    kf = k.astype(np.float64)
    return np.sum(t[k] * (np.exp(-s * np.log(kf - 0.5)) - np.exp(-s * np.log(kf))))


rows = []
zs = [complex(0.99820922076525681, 240.24783014069124), complex(0.69391829382601, 22.253194573399)]
for s in [3 + 5j, 2 + 30j, 1.2 - 7j, 0.7 + 100j, 0.55 + 1000j, 0.3 + 3j] + zs:
    for N in (1000, 12345, NN):
        Lc = (2 * N - 1) // 3
        Jc = Lc // 2
        FN, FL, FJ = Fsum(t, s, N), Fsum(t, s, Lc), Fsum(t, s, Jc)
        FH = Fsum(t * w, s, N)
        EN = Esum(s, N)
        lhs = -(1.5 ** s) * (FH / 3 + EN) - D(s) * FN
        rhs = -(FN - FL) - 2.0 ** (-s) * (FN - FJ)
        rows.append([str(s), N, abs(lhs - rhs), abs(rhs)])
out['B2_identity_rows'] = rows
out['B2_identity_max_abs_error'] = max(r[2] for r in rows)
print('(1) even/odd relations:', ok_even, ok_odd, '; cut-off identity max |lhs - rhs| over 24 (s, N):',
      out['B2_identity_max_abs_error'], flush=True)
for r in rows:
    print('    s =', r[0], 'N =', r[1], ' |lhs - rhs| =', f'{r[2]:.2e}', ' |rhs| =', f'{r[3]:.3e}')
# (2) Lemma C constants
f15 = 2 ** -1.5 + 1.5 ** 1.5 / 3
out['lemmaC_lower_bound_on_Re_3/2'] = 1 - f15
out['lemmaC_max_D_on_disk_bound'] = 1 + 2.0 ** 1 + 1.5 ** 4 / 3
out['lemmaC_jensen_ratio_r1.6'] = math.log(5 / (1 - f15)) / math.log(2.5 / 1.6)
out['lemmaC_corner_distance_rectangle'] = math.hypot(1.5, 0.5)
eta = 1 / 48
out['lemmaC_corner_distance_enlarged_rectangle'] = math.hypot(1.5 + eta, 0.5 + eta)
out['lemmaC_jensen_ratio_r1.61'] = math.log(5 / (1 - f15)) / math.log(2.5 / 1.61)
print('(2) Lemma C:', {k: v for k, v in out.items() if k.startswith('lemmaC')}, flush=True)
# (3) Riesz kernel
mp.mp.dps = 30
kr = []
for y in (0.5, 1.0, 1.7, 10.0):
    f = lambda tt: (mp.mpf(y) ** (2 + 1j * tt) * 2 / ((2 + 1j * tt) * (3 + 1j * tt) * (4 + 1j * tt)))
    val = mp.quad(f, [-mp.inf, -50, 0, 50, mp.inf]) / (2 * mp.pi)
    kr.append([y, float(mp.re(val)), (1 - 1 / y) ** 2 if y >= 1 else 0.0])
out['riesz_kernel_check'] = kr
print('(3) Riesz kernel (y, integral, (1 - 1/y)^2 [y >= 1]):', kr, flush=True)
json.dump(out, open(os.path.join(HERE, 'a3_thm4660_identity_constants.json'), 'w'), indent=1, default=str)
