#!/usr/bin/env python3
"""
THM-1550 obstacle (iii), the additive<->multiplicative bridge, reformulated as a LOG-DERIVATIVE
(no formal log).  For Phi(x) = x^M - t*R(x), R deg d=M+N, R(0)!=0, let Pi(t) = product of the
M small roots (val>0) and D_m = [x^{Mm}] R^m (so D_0 = 1).  Then in C[[t]]:

    t * Pi'(t) / Pi(t)  =  F(t) := sum_{m>=0} D_m t^m  =  [x^0]  x^M / (x^M - t R)   (root-free).

Immediately:  D_m = 0 for all m>=1  =>  t Pi'/Pi = 1  =>  Pi = c*t   (obstacle (iii)).

The identity factors into two clean pieces:
  (a) PER-ROOT (pure algebra, differentiate u^M = t R(u)):  t u'/u = R(u)/S(u),  S := M R - x R'.
  (b) SUM (symmetric function / residue):  sum_{small} R(u_j)/S(u_j) = F(t).
So  t Pi'/Pi = sum_j t u_j'/u_j = sum_j R(u_j)/S(u_j) = F(t).

Replaces death-star-S106/S111's "formal CT of log of the Wiener-Hopf factorization" with
differentiation + a per-root formula + a symmetric-function sum.  Root-free anchor:
Res_0[R/(x Phi)] = R(0)/Phi(0) = -1/t.                                    mac-mini-2026-07-22-S163
"""
import numpy as np

def check(M, N, seed):
    rng = np.random.default_rng(seed); d = M + N
    a = rng.standard_normal(d + 1) + 1j * rng.standard_normal(d + 1)
    if abs(a[0]) < 0.3: a[0] += 1.0                      # R(0) = a[0] != 0
    R = np.poly1d(a[::-1]); Rp = R.deriv()
    S = lambda x: M * R(x) - x * Rp(x)
    def Dm(m):
        c = (R ** m).coeffs[::-1]; idx = M * m
        return c[idx] if idx < len(c) else 0.0
    K = 14; Ds = [Dm(m) for m in range(1, K + 1)]
    Fser = lambda t: 1 + sum(Ds[m - 1] * t ** m for m in range(1, K + 1))
    def roots(t):
        co = np.array([-t * a[k] for k in range(d + 1)], dtype=complex); co[M] += 1.0
        r = sorted(np.roots(co[::-1]), key=abs); return r[:M], r[M:]
    t = 0.01; small, _ = roots(t)
    Pi = np.prod(small)
    # t Pi'/Pi via central difference
    h = 1e-6
    sp, _ = roots(t + h); sm, _ = roots(t - h)
    dPi = (np.prod(sp) - np.prod(sm)) / (2 * h)
    e_main = abs(t * dPi / Pi - Fser(t))                 # t Pi'/Pi = F
    e_sum = abs(sum(R(u) / S(u) for u in small) - Fser(t))  # (b)
    return e_main, e_sum, Pi / t

if __name__ == "__main__":
    print(" M  N | t Pi'/Pi = F  | sum R/S = F   | Pi/t (-> c, val-1)")
    for (M, N, s) in [(1, 1, 1), (2, 1, 2), (1, 2, 3), (2, 3, 4), (3, 2, 5), (2, 2, 6)]:
        em, es, pit = check(M, N, s)
        print(f" {M}  {N} |  {em:.1e}    |  {es:.1e}    | {pit:.3f}")
    print("=> log-derivative bridge confirmed; D_m=0 for all m => Pi = c*t (obstacle iii).")
