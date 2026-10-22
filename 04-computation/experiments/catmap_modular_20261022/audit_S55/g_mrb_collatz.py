"""Audit G (cont.): synthesis-note sections 5-6: Ising ring at exp(-2K) = sqrt5, Collatz LTE claim, MRB partial sums
and the trace argument, independent recomputation of the prime defect D_p.  Auditor for opus-2026-10-22-S55.
Usage: python g_mrb_collatz.py"""
import math, cmath
import sympy as sp
import mpmath as mp
import numpy as np

# ---------------- Ising ring
x = sp.sqrt(5)                                      # x = exp(-2K)
K = -sp.log(5) / 4
Tm = sp.Matrix([[sp.exp(K), sp.exp(-K)], [sp.exp(-K), sp.exp(K)]])
Tflip = sp.Matrix([[sp.exp(-K), sp.exp(K)], [sp.exp(K), sp.exp(-K)]])          # one antiperiodic (sign-flipped) bond
# exact algebra in Q(sqrt5) with the normalised matrices e^(-K) T = [[1, x], [x, 1]], e^(-K) T' = [[x, 1], [1, x]]
Tn = sp.Matrix([[1, x], [x, 1]]); Tf = sp.Matrix([[x, 1], [1, x]])
ok_P = all(sp.expand((Tn ** n).trace() - 2 ** n * sp.lucas(n)) == 0 for n in range(1, 21))
ok_A = all(sp.expand((Tn ** (n - 1) * Tf).trace() - 2 ** n * sp.sqrt(5) * sp.fibonacci(n)) == 0 for n in range(1, 21))
ev = sorted((Tm * sp.exp(-K)).eigenvals().keys(), key=lambda z: float(z))
print('Ising ring, exp(-2K) = sqrt5 (K = -ln5/4 < 0, antiferromagnetic):  e^(-nK) Z_periodic = 2^n L_n (n<=20):', ok_P,
      '; e^(-nK) Z_antiperiodic = 2^n sqrt5 F_n:', ok_A, '; eigenvalue ratio', sp.nsimplify(sp.simplify(ev[0] / ev[1])),
      '= -1/phi^2:', sp.simplify(ev[0] / ev[1] + 1 / ((1 + sp.sqrt(5)) / 2) ** 2) == 0)
print('   unnormalised: Z_periodic(n) = 5^(-n/4) 2^n L_n, e.g. n = 2:', sp.nsimplify(sp.simplify((Tm ** 2).trace())))

# ---------------- Collatz LTE
def v2(n):
    return (n & -n).bit_length() - 1
bad = [k for k in range(1, 2001) if v2(3 ** k - 1) != (1 if k % 2 else 2 + v2(k))]
print('\nv_2(3^k - 1) = 1 (k odd), 2 + v_2(k) (k even), k = 1..2000: failures', bad)
traj_ok = True
for k in range(1, 80):
    y = 2 ** k - 1
    for j in range(k):
        traj_ok &= (y % 2 == 1)
        y = (3 * y + 1) // 2
    traj_ok &= (y == 3 ** k - 1)
print('2^k - 1 -> (k steps of (3n+1)/2, all on odd values) -> 3^k - 1, k = 1..79:', traj_ok)
print('k with v_2(3^k - 1) = k (k <= 2000):', [k for k in range(1, 2001) if v2(3 ** k - 1) == k], '(k = 1 is 1 -> 2; k = 4 is 15 -> 80)')

# ---------------- MRB partial sums S_N = sum_{n<=N} (-1)^n (n^(1/n) - 1)
mp.mp.dps = 30
S = [mp.mpf(0)]
for n in range(1, 60):
    S.append(S[-1] + (-1) ** n * (mp.root(n, n) - 1))
odd = {N: S[N] for N in range(1, 60, 2)}; even = {N: S[N] for N in range(2, 60, 2)}
print('\nS_1, S_3, S_5, S_7 =', [mp.nstr(S[N], 8) for N in (1, 3, 5, 7)])
print('min odd S_N, 3 <= N <= 59:', mp.nstr(min(odd[N] for N in odd if N >= 3), 8), 'at N =', min((odd[N], N) for N in odd if N >= 3)[1],
      '; min odd S_N, 5 <= N <= 59:', mp.nstr(min(odd[N] for N in odd if N >= 5), 8),
      '; min |S_N| odd 3 <= N <= 59:', mp.nstr(min(abs(odd[N]) for N in odd if N >= 3), 8), 'at N =', min((abs(odd[N]), N) for N in odd if N >= 3)[1])
print('min even S_N, N <= 59:', mp.nstr(min(even.values()), 8), 'at N =', min((v, N) for N, v in even.items())[1],
      '; S_5 > 1/200:', S[5] > mp.mpf(1) / 200)
# trace argument: minimal polynomial of n^(1/n) has zero x^(d-1) coefficient
xs = sp.Symbol('x')
tr0 = []
for n in range(2, 31):
    mpoly = sp.Poly(sp.minimal_polynomial(sp.Integer(n) ** sp.Rational(1, n), xs), xs)
    cs = mpoly.all_coeffs()
    tr0.append(len(cs) > 2 and cs[1] == 0 or (len(cs) == 2 and cs[1] != 0 and False))
print('n^(1/n), n = 2..30: minimal polynomial has zero trace coefficient (and degree >= 2):', all(tr0))

# ---------------- prime defect D_p (formula as encoded in mrb_checks.py; source PDF not in the repo)
def vp_arr(m, p):
    v = np.zeros_like(m)
    mm = m.copy()
    while True:
        mask = (mm % p == 0)
        if not mask.any():
            break
        v[mask] += 1
        mm[mask] //= p
    return v


def Dp(p, Mmax=4_000_000):
    m = np.arange(1, Mmax + 1, dtype=np.int64)
    x = (p * m).astype(float)
    v = vp_arr(m, p).astype(float)
    terms = (-1.0) ** m * np.exp(np.log(x) / x) * (np.exp(4j * np.pi * (1 + v) / x) - 1)
    # sum over complete blocks of 2p consecutive m (one full period of sign and p-divisibility), then average the
    # partial sums at the last two block ends (alternating-block tail)
    cs = np.cumsum(terms)
    B = 2 * p
    ends = cs[B - 1::B]
    return (ends[-1] + ends[-2]) / 2, ends[-1] - ends[-2]


for p in (3, 5, 7, 11, 101, 211):
    d, last = Dp(p)
    main = -4j * math.pi * math.log(2) / (p - 1)
    bound = (30 * math.log(p) + 661) / p ** 2
    print(f'p = {p:3d}: D_p = {d.real:+.6f}{d.imag:+.6f}i (last block {abs(last):.1e}), |D_p - main| = {abs(d - main):.3e}, '
          f'bound {bound:.3e}, within: {abs(d - main) < bound}')
