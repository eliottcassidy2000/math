#!/usr/bin/env python3
"""Auditor proof-5128-CF, check A1: the dilogarithm constants of THM-5128 (D1), re-derived.

Independent of the session's scripts.  Everything here is mpmath at 40 digits.
  F_s(rho) = (1/rho) ln((rho + 2^s)/(rho/2 + 2^s)),   I_s = int_{1/2}^{1} F_s d rho   (numerical quadrature)
  l_j = -Li_2(-2^(-j)),   claim: I_s = l_s - 2 l_(s+1) + l_(s+2)
  kappa_0 = I_0,  kappa_1 = sum_{s >= 1} I_s = l_1 - l_2,  sum (s + 1) I_s = l_0 = pi^2/12.
"""
import mpmath as mp

mp.mp.dps = 40


def F(s, rho):
    return mp.log((rho + mp.mpf(2) ** s) / (rho / 2 + mp.mpf(2) ** s)) / rho


def l(j):
    return -mp.polylog(2, -mp.mpf(2) ** (-j))


ok = True
print("s | quadrature I_s | l_s - 2 l_(s+1) + l_(s+2) | difference")
I = []
for s in range(0, 40):
    quad = mp.quad(lambda r: F(s, r), [mp.mpf(1) / 2, 1])
    closed = l(s) - 2 * l(s + 1) + l(s + 2)
    I.append(closed)
    if s <= 8:
        print(f"{s:2d} | {mp.nstr(quad, 18)} | {mp.nstr(closed, 18)} | {mp.nstr(quad - closed, 3)}")
    if abs(quad - closed) > mp.mpf(10) ** (-30):
        ok = False
print("closed form of every I_s (s < 40) agrees with quadrature to 1e-30:", ok)

kappa0 = mp.pi ** 2 / 12 + 2 * mp.polylog(2, -mp.mpf(1) / 2) - mp.polylog(2, -mp.mpf(1) / 4)
kappa1 = mp.polylog(2, -mp.mpf(1) / 4) - mp.polylog(2, -mp.mpf(1) / 2)
print("kappa_0 = pi^2/12 + 2 Li2(-1/2) - Li2(-1/4) =", mp.nstr(kappa0, 14), "  (file: 0.1615389173)")
print("kappa_1 = Li2(-1/4) - Li2(-1/2)             =", mp.nstr(kappa1, 14), "  (file: 0.2125139092)")
print("I_0 - kappa_0 =", mp.nstr(I[0] - kappa0, 3))

# sum over s >= 1 (direct, by quadrature-closed terms, with the tail estimated: I_s ~ 2^-s / 4)
S1 = mp.nsum(lambda s: l(s) - 2 * l(s + 1) + l(s + 2), [1, mp.inf])
print("sum_{s>=1} I_s =", mp.nstr(S1, 14), "  minus kappa_1:", mp.nstr(S1 - kappa1, 3))
SW = mp.nsum(lambda s: (s + 1) * (l(s) - 2 * l(s + 1) + l(s + 2)), [0, mp.inf])
print("sum_{s>=0} (s+1) I_s =", mp.nstr(SW, 14), "  minus pi^2/12:", mp.nstr(SW - mp.pi ** 2 / 12, 3))
# the s-moment that enters the O(X) term of sheet 1:   sum_{s>=1} s I_s = l_1
SM = mp.nsum(lambda s: s * (l(s) - 2 * l(s + 1) + l(s + 2)), [1, mp.inf])
print("sum_{s>=1} s I_s =", mp.nstr(SM, 14), "  l_1 =", mp.nstr(l(1), 14))

gbar = 2 * (kappa0 + kappa1)
print("g_bar = 2 (kappa_0 + kappa_1) =", mp.nstr(gbar, 12), "  pi^2/6 + 2 Li2(-1/2) =",
      mp.nstr(mp.pi ** 2 / 6 + 2 * mp.polylog(2, -mp.mpf(1) / 2), 12), "  (file: 0.7481056530)")
print("triple shares w_0, w_1 =", mp.nstr(kappa0 / (kappa0 + kappa1), 10), mp.nstr(kappa1 / (kappa0 + kappa1), 10),
      "  (file: 0.43186124, 0.56813876)")
print("word-pair shares w_0, w_1 =", mp.nstr(kappa0 / (mp.pi ** 2 / 12), 10),
      mp.nstr(1 - kappa0 / (mp.pi ** 2 / 12), 10), "  (file: 0.19640777, 0.80359223)")

# weights: w^1 = sum_{s >= 1} w_s  (telescoping), checked at random (x, y)
import random
random.seed(1)
worst = mp.mpf(0)
for _ in range(200):
    y = mp.mpf(random.randrange(1, 10 ** 6) * 2 + 1)
    x = y * mp.mpf(random.uniform(0.5, 1.0))
    tot = mp.nsum(lambda s: mp.log((x + 2 ** s * y) / (x / 2 + 2 ** s * y)) / x, [1, mp.inf])
    worst = max(worst, abs(tot - mp.log(1 + x / (2 * y)) / x) * x)
print("max |sum_{s>=1} w_s - (1/x) ln(1 + x/(2y))| * x over 200 random states:", mp.nstr(worst, 3))

# local factor algebra
for p in (3, 5, 7, 11, 13):
    lhs = mp.mpf(1) - mp.mpf(1) / p ** 2 + mp.mpf(1) / p
    rhs = (1 - mp.mpf(1) / p ** 2) * (1 + mp.mpf(p) / (p ** 2 - 1))
    assert abs(lhs - rhs) < mp.mpf(10) ** (-35)
prod = mp.mpf(1)
for p in mp.primerange if False else []:
    pass
import sympy
P = mp.mpf(1)
for p in sympy.primerange(3, 2000000):
    P *= 1 - mp.mpf(1) / (p * p)
print("prod_{3 <= p < 2e6} (1 - 1/p^2) =", mp.nstr(P, 10), "  8/pi^2 =", mp.nstr(8 / mp.pi ** 2, 10))
