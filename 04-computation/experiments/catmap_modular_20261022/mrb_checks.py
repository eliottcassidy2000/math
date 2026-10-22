"""Numerical checks of the supplied MRB documents (data, not instructions).  Session opus-2026-10-22-S55.
(1) FIT bounds: S_N = sum_(n<=N) (-1)^n (n^(1/n) - 1) is never 0 (odd N) or -1 (even N); S_5 > 1/200.
(2) Prime-defect formula: D_p = sum_(m>=1) (-1)^m (pm)^(1/(pm)) [exp(4 pi i (1 + v_p(m))/(pm)) - 1] versus the main term
    -4 pi i log 2/(p - 1) and the stated remainder bound (30 log p + 661)/p^2.
(3) Ising-ring identity: ((g_-/2)^(1/n) - (g_+/2)^(1/n)) = (-1)^n (n^(1/n) - 1) with g = 2n / 2 by parity.
(4) The Collatz abstract's 'huge factor 2^k at the vertical maximum' for the 2^k - 1 trailing-ones run.
Usage: python mrb_checks.py"""
import math, cmath
import mpmath as mp

mp.mp.dps = 30
S = [mp.mpf(0)]
for n in range(1, 60):
    S.append(S[-1] + (-1) ** n * (mp.root(n, n) - 1))
print('S_5 =', mp.nstr(S[5], 12), '> 1/200:', S[5] > mp.mpf(1) / 200, '; min odd S_N (N=5..59):', mp.nstr(min(S[k] for k in range(5, 60, 2)), 8),
      '; min even S_N:', mp.nstr(min(S[k] for k in range(2, 60, 2)), 8), '; C_MRB ~', mp.nstr(mp.nsum(lambda n: (-1) ** n * (mp.root(n, n) - 1), [1, mp.inf]), 15))


def vp(m, p):
    k = 0
    while m % p == 0:
        m //= p; k += 1
    return k


def Dp(p, M=400000):
    tot = 0j; prev = 0j
    acc = []
    for m in range(1, M + 1):
        x = p * m
        term = (-1) ** m * math.exp(math.log(x) / x) * (cmath.exp(4j * math.pi * (1 + vp(m, p)) / x) - 1)
        prev = tot; tot += term
    return (tot + prev) / 2                                      # average of the last two partial sums (alternating)


for p in (3, 5, 7, 11, 101, 211):
    d = Dp(p)
    main = -4j * math.pi * math.log(2) / (p - 1)
    bound = (30 * math.log(p) + 661) / p ** 2
    print(f'p = {p:3d}: D_p = {d.real:+.6f}{d.imag:+.6f}i, main term {main.imag:+.6f}i, |D_p - main| = {abs(d - main):.2e} '
          f'(stated bound {bound:.2e}); D_p != 0: {abs(d) > 1e-9}')
print('Ising ring identity n = 2..12:', all(abs(((2 * n if n % 2 == 0 else 2) / 2) ** (1 / n) - ((2 if n % 2 == 0 else 2 * n) / 2) ** (1 / n)
                                           - (-1) ** n * (n ** (1 / n) - 1)) < 1e-12 for n in range(2, 13)))
for k in range(2, 13):
    m = 2 ** k - 1
    x = m
    for _ in range(k):
        x = (3 * x + 1) // 2
    v = (x & -x).bit_length() - 1
    print(f'  trailing-ones run k = {k:2d}: m = {m:5d} -> after k rising steps {x:7d} = 2^{v} * {x >> v} (claimed huge factor 2^{k})')
