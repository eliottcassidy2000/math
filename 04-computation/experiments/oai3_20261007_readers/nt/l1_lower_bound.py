#!/usr/bin/env python3
"""Explicit lower bound for L(1,chi), chi real primitive mod q, assuming only
   L(sigma,chi) != 0 for 7/8 < sigma < 1   (openai/math #003 gives this for every chi).

Goldfeld-type argument (Goldfeld 1975; Iwaniec-Kowalski Thm 5.28 proof) with explicit constants:
  F = zeta*L = sum a(n) n^-s, a = 1*chi >= 0, a(1) = 1.
  S(x) = sum a(n) n^-beta e^{-n/x} >= e^{-1/x}
       = L(1,chi) Gamma(1-beta) x^{1-beta} + F(beta) - [beta+delta>1] F(beta-1)/x + I,
  F(beta) = zeta(beta) L(beta,chi) < 0  (L(beta,chi) > 0 since no zero in [beta,1), beta > 7/8),
  |F(beta-1)| <= |zeta(beta-1)| (q/pi)^{3/2-beta} R_a(1-beta,0) zeta(2-beta),
  |I| <= x^{-beta-delta} q^{1/2+delta} J_a(beta,delta),
  J_a = (1/2pi) pi^{-1-2delta} zeta(1+delta)^2 int R_0(delta,t) R_a(delta,t) |Gamma(-beta-delta+it)| dt,
  R_a(d,t) = |Gamma((1+d+a+it)/2)| / |Gamma((-d+a+it)/2)|   (functional equation + |L(1+d+it)| <= zeta(1+d)).
Hence L(1,chi) >= (e^{-1/x} - E1 - E2) / (Gamma(1-beta) x^{1-beta}).
NUMERICAL (mpmath quadrature at 30 digits); the reported constant is then shaved by 1%.
"""
import mpmath as mp
mp.mp.dps = 30

def R(a, d, t):
    return abs(mp.gamma((1 + d + a + 1j * t) / 2)) / abs(mp.gamma((-d + a + 1j * t) / 2))

_Jcache = {}
def J(a, beta, d):
    key = (a, float(beta), float(d))
    if key in _Jcache: return _Jcache[key]
    f = lambda t: R(0, d, t) * R(a, d, t) * abs(mp.gamma(-beta - d + 1j * t))
    integral = 2 * mp.quad(f, [0, 0.5, 2, 8, 20, 40, 80])
    val = integral * mp.zeta(1 + d) ** 2 * mp.pi ** (-1 - 2 * d) / (2 * mp.pi)
    _Jcache[key] = val
    return val

def lower(q, a, beta, d, logx):
    x = mp.e ** logx
    E2 = mp.e ** (-(beta + d) * logx) * mp.mpf(q) ** (mp.mpf(1) / 2 + d) * J(a, beta, d)
    E1 = 0
    if beta + d > 1:
        Fm = abs(mp.zeta(beta - 1)) * (mp.mpf(q) / mp.pi) ** (mp.mpf(3) / 2 - beta) * R(a, 1 - beta, 0) * mp.zeta(2 - beta)
        E1 = Fm / x
    num = mp.e ** (-1 / x) - E1 - E2
    if num <= 0: return mp.mpf(0)
    return num / (mp.gamma(1 - beta) * mp.e ** ((1 - beta) * logx))

def best(q, a):
    bestv = (mp.mpf(0), None)
    for d in [0.08, 0.12, 0.16, 0.2, 0.25, 0.3, 0.35, 0.4]:
        for beta in [0.876, 0.89, 0.9, 0.92, 0.94, 0.96, 0.97, 0.98, 0.985, 0.99]:
            if beta + d == 1: continue
            # scan log x
            lq = mp.log(q)
            for k in range(0, 61):
                logx = lq * (0.3 + 0.02 * k)
                v = lower(q, a, mp.mpf(beta), mp.mpf(d), logx)
                if v > bestv[0]: bestv = (v, (d, beta, float(logx / lq)))
    return bestv

if __name__ == "__main__":
    import sys
    for a in (1, 0):
        print(f"parity a={a} ({'odd: imaginary quadratic' if a else 'even: real quadratic'})")
        for k in [4, 5, 6, 8, 10, 11, 12, 13, 14, 16, 20, 30, 50, 100]:
            q = mp.mpf(10) ** k
            v, par = best(q, a)
            c = v * mp.log(q)
            print(f"  q=1e{k:<3d} L(1,chi) >= {mp.nstr(v, 6):>10s}   (= {mp.nstr(c, 5)}/log q)  params d,beta,logx/logq = {par}")
            sys.stdout.flush()
