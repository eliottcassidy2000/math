#!/usr/bin/env python3
"""Asymptotic threshold of the archimedean necessary condition (AMM 12592).

(ARCH) says, for a depth profile a_k and every d,

    C(m-1,d)  <=  sum_{k : 0 <= d-L_k <= a_k} C(a_k, d-L_k) 2^{a_k-d+L_k}.

Scale k = xm, d = delta m, a_k = alpha(x) m, L_k = ell(x) m.  With
C = 1+gamma the target ratio, the profile a_k = min(m-1-k, gamma(m+k)) gives

    kappa = (1-gamma)/(1+gamma),
    x <= kappa :  alpha = gamma(1+x),   ell = (1-gamma) - x(1+gamma),
    x >= kappa :  alpha = 1-x,          ell = 0.

Taking log_2 and dividing by m, C(m-1,d) -> H(delta) and each summand ->
alpha H(r/alpha) + (alpha - r) with r = delta - ell(x).  So (ARCH) survives
in the limit iff

    for all delta in [0,1]:   H(delta)  <=  max_x [ alpha H(r/alpha) + alpha - r ].

The least gamma for which this holds is the asymptotic archimedean floor
C_arch = 1 + gamma*: no balanced block scheme of slope below C_arch can exist.
"""
from math import log2

GRID_X = 4001
GRID_D = 2001


def H(p):
    if p <= 0 or p >= 1:
        return 0.0
    return -p * log2(p) - (1 - p) * log2(1 - p)


def profile_scaled(gamma, x):
    kappa = (1 - gamma) / (1 + gamma)
    if x <= kappa:
        return gamma * (1 + x), (1 - gamma) - x * (1 + gamma)
    return 1 - x, 0.0


def capacity(gamma, delta):
    """max_x [ alpha H(r/alpha) + alpha - r ]  over admissible x."""
    best = -1e9
    for i in range(GRID_X):
        x = i / (GRID_X - 1)
        alpha, ell = profile_scaled(gamma, x)
        if alpha <= 0:
            continue
        r = delta - ell
        if r < -1e-12 or r > alpha + 1e-12:
            continue
        r = min(max(r, 0.0), alpha)
        val = alpha * H(r / alpha) + (alpha - r)
        if val > best:
            best = val
    return best


def feasible(gamma):
    """Does (ARCH) survive for every delta?"""
    worst = None
    for j in range(GRID_D):
        delta = j / (GRID_D - 1)
        need = H(delta)
        have = capacity(gamma, delta)
        if have < need - 1e-9:
            if worst is None or need - have > worst[1]:
                worst = (delta, need - have)
    return worst is None, worst


def threshold():
    lo, hi = 0.0, 1.0
    for _ in range(40):
        mid = (lo + hi) / 2
        ok, _ = feasible(mid)
        if ok:
            hi = mid
        else:
            lo = mid
    return hi


if __name__ == "__main__":
    g = threshold()
    print(f"gamma* = {g:.6f}")
    print(f"C_arch = 1 + gamma* = {1+g:.6f}")
    ok, worst = feasible(g + 1e-4)
    print(f"just above: feasible={ok}")
    ok2, worst2 = feasible(g - 1e-3)
    print(f"just below: feasible={ok2}  worst delta={worst2[0]:.4f} "
          f"deficit={worst2[1]:.5f}" if worst2 else "")
    print()
    print(" gamma   C=1+gamma   ARCH survives?   binding delta")
    for gg in (0.50, 0.55, 0.58, 0.60, 0.62, 0.63, 0.64, 0.65, 0.70, 0.80):
        ok, w = feasible(gg)
        b = "-" if w is None else f"{w[0]:.3f}"
        print(f" {gg:5.3f}   {1+gg:8.4f}   {str(ok):14s}   {b}")
