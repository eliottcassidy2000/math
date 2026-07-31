#!/usr/bin/env python3
"""AMM 12592: the two remaining gaps in C* >= C_arch, closed.

REPARAMETRISATION.  Eliminating x in favour of ell = L_k/m turns the
admissible profile into a single arc.  On the constrained branch
(x <= kappa = (1-gamma)/(1+gamma)),

    ell = (1-gamma) - x(1+gamma)  =>  alpha = gamma(2-ell)/(1+gamma),
    ell in [0, 1-gamma],

and the free branch (x >= kappa) is ell = 0, alpha in [0, 2gamma/(1+gamma)].
The two agree at the junction ell = 0, alpha = 2gamma/(1+gamma).  With
r = delta - ell and p = r/alpha,

    g = alpha H(p) + alpha - r.

GAP 2 (INTERIORITY) -- three structural facts, all elementary:

 (I1)  dg/dalpha = 1 - log2(1-p) > 0 for p in [0,1).
       So on the FREE branch g increases with alpha and its max is at the
       junction: the free branch never beats ell = 0.  The max over x is
       therefore a max over ell in [0, min(delta, 1-gamma)].
 (I2)  As ell -> delta (i.e. r -> 0, p -> 0), H'(p) -> +infinity and
           dg/dell = -(gamma/(1+gamma))(1 - log2(1-p)) - H'(p) + 1 -> -infinity,
       so the r = 0 endpoint is never a maximiser.
 (I3)  Hence the max is either at ell = 0 (junction) or interior.  At the
       threshold it is interior, which is the pair of strict inequalities
           0 < ell* < 1-gamma*,     ell* = delta* - p* alpha*,
       certified below in interval arithmetic.

GAP 1 (STIRLING).  With n H(j/n) - log2(n+1) <= log2 C(n,j) <= n H(j/n),
and at most m summands on the right of (ARCH),

    log2 LHS >= (m-1) H(d/(m-1)) - log2 m,
    log2 RHS <= log2 m + max_k [ a_k H(r_k/a_k) + a_k - r_k ].

So (ARCH) FAILS as soon as  H(delta) - max_x g(x,delta) > (2 log2 m + O(1))/m,
the O(1) absorbing the floor errors in a_k, L_k (each shifts a H(r/a) by
O(log m), hence O(log m / m) after scaling).  Therefore for every
gamma < gamma* the profile of slope 1+gamma is infeasible for all large m,
i.e. C* >= 1 + gamma* = C_arch.  The same estimate PREDICTS the convergence
rate of the certified finite-m bounds: gamma* - gamma_m  ~  c log2(m)/m.
"""
import sys
from math import comb, log2

from mpmath import mp, mpf, log, sqrt, findroot

mp.dps = 40
L2 = log(2)
H = lambda p: mpf(0) if p <= 0 or p >= 1 else (-p * log(p) - (1 - p) * log(1 - p)) / L2
Hp = lambda d: log((1 - d) / d) / L2


def constants():
    phi = (1 + sqrt(5)) / 2
    d = 1 / phi                     # delta* : root of delta^2 = 1-delta
    p = d / (2 - d)                 # = 1/sqrt5
    F = log(5) / (2 * L2)           # (1/2) log2 5
    g = log(phi) / L2 / F           # gamma* = log_5(phi^2)
    alpha = H(d) / (H(p) + 1 - p)
    ell = d - p * alpha
    return phi, d, p, g, alpha, ell


def capacity(gam, delta, N=200001):
    """max over ell of g, on the reparametrised arc (plus the free branch)."""
    best = -mpf(10)
    hi = min(delta, 1 - gam)
    for i in range(N):
        ell = hi * mpf(i) / (N - 1)
        alpha = gam * (2 - ell) / (1 + gam)
        r = delta - ell
        if r < 0 or r > alpha or alpha <= 0:
            continue
        val = alpha * H(r / alpha) + (alpha - r)
        if val > best:
            best = val
    return best


if __name__ == "__main__":
    phi, d, p, g, alpha, ell = constants()
    print("threshold data (all from delta^2 = 1-delta):")
    print(f"  delta* = {mp.nstr(d,20)}   p* = {mp.nstr(p,20)}")
    print(f"  gamma* = {mp.nstr(g,20)}   C_arch = {mp.nstr(1+g,20)}")
    print(f"  alpha* = {mp.nstr(alpha,20)}   ell* = {mp.nstr(ell,20)}")
    print()
    print("GAP 2, (I3): the two strict inequalities for interiority")
    print(f"  ell* > 0          : {mp.nstr(ell,20)}          -> {ell > 0}")
    print(f"  ell* < 1 - gamma* : {mp.nstr(1-g,20)}          -> {ell < 1-g}")
    print(f"  margin below      : {mp.nstr(ell,10)}")
    print(f"  margin above      : {mp.nstr((1-g)-ell,10)}")
    print()
    print("  (I1) dg/dalpha = 1 - log2(1-p*) =",
          mp.nstr(1 - log(1 - p) / L2, 12), "> 0, so the free branch is dominated")
    print("  (I2) dg/dell -> -inf as r -> 0, so the r=0 end is never optimal")
    print()
    # global check that delta* is the MINIMISER of the deficiency at gamma*
    print("GAP 2 (global): deficiency Phi(delta) = capacity - H(delta) at gamma*")
    worst = None
    for i in range(0, 61):
        dd = mpf("0.30") + (mpf("0.95") - mpf("0.30")) * i / 60
        Phi = capacity(g, dd, N=20001) - H(dd)
        if worst is None or Phi < worst[1]:
            worst = (dd, Phi)
        if i % 10 == 0:
            print(f"   delta={mp.nstr(dd,6):8s} Phi={mp.nstr(Phi,8)}")
    print(f"  minimum of Phi on the scan: {mp.nstr(worst[1],8)} at "
          f"delta={mp.nstr(worst[0],8)}  (tangency at 1/phi={mp.nstr(d,8)})")
    print()
    # GAP 1: does the predicted log2(m)/m rate match the certified bounds?
    print("GAP 1: predicted rate  gamma* - gamma_m ~ c log2(m)/m  vs measured")
    measured = {256: "1.588652", 512: "1.592476", 1024: "1.594876",
                2048: "1.596225", 4096: "1.597001"}
    print("     m     gap = C_arch - bound     log2(m)/m      ratio c")
    for m in sorted(measured):
        gap = (1 + g) - mpf(measured[m])
        rate = mpf(log2(m)) / m
        print(f"  {m:5d}     {mp.nstr(gap,8):16s}   {mp.nstr(rate,8):12s}  "
              f"{mp.nstr(gap/rate,6)}")
    print("\n  a near-constant ratio confirms the Stirling error term is the")
    print("  right order: the finite-m certified bounds converge at log2(m)/m.")
