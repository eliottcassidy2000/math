#!/usr/bin/env python3
"""epsilon-fair extractors: does the archimedean floor survive approximate fairness?

RELAXATION.  epsilon-fair means |P_p(H) - 1/2| <= epsilon for all p, i.e.
|sum_m D_m(p)| <= 2 epsilon.  Since D_m = sum_j c_j p^{2m-j}q^j and the
Bernstein functions binom(2m,j)p^{2m-j}q^j are a partition of unity,

    |D_m(p)| <= max_j |c_j| / binom(2m,j),

so it SUFFICES to allow |c_j| <= 2 epsilon binom(m,j) at the relative layers.
The within-block identity then becomes

    sum_k E_k(u)(1+u)^{L_k} = +- u^{m-1} + Gamma(u),
    Gamma(u) = sum_j c'_j u^{j-1},   |c'_j| <= 2 epsilon binom(m,j).

Expanding at u = -1 and using sum_j binom(m,j) binom(j,d) = binom(m,d) 2^{m-d},

    |[v^d] Gamma(-1+v)| <= 2 epsilon binom(m,d) 2^{m-d},

so (ARCH) weakens to

    binom(m-1,d) - 2 epsilon binom(m,d) 2^{m-d}
        <= sum_{k: 0<=d-L_k<=a_k} binom(a_k,d-L_k) 2^{a_k-d+L_k}.   (ARCH-eps)

THE SCALING IS DECISIVE.  With d = delta m, log_2 of the main term is
m H(delta) but log_2 of the correction is m[H(delta) + 1 - delta], which is
EXPONENTIALLY LARGER for every delta < 1.  So for any FIXED epsilon > 0 the
left side goes negative and (ARCH-eps) is vacuous once m is large: the
archimedean floor is a phenomenon of EXACT fairness only.  It survives
exactly when 2 epsilon 2^{m(1-delta*)} = O(1), i.e.

    epsilon  <~  2^{-m(1-delta*)} = 2^{-m/phi^2},    1 - 1/phi = 1/phi^2
                                                    = 0.3819660...

PREDICTION: writing epsilon = 2^{-beta m}, the certified floor survives for
beta > 1/phi^2 and collapses for beta < 1/phi^2.  That critical exponent is
tested below.
"""
import sys
from fractions import Fraction
from math import comb

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(1000000)


def profile(m, C):
    a = []
    for k in range(m):
        n = m + k
        cap = Fraction(C) * n - n - 1
        if cap < 0:
            return None
        a.append(min(m - 1 - k, int(cap)))
    return a


def arch_eps_ok(m, a, num, den):
    """(ARCH-eps) with epsilon = num/den (exact rationals)."""
    L = [m - 1 - k - a[k] for k in range(m)]
    for d in range(m):
        lhs = Fraction(comb(m - 1, d)) - 2 * Fraction(num, den) * comb(m, d) * 2 ** (m - d)
        if lhs <= 0:
            continue                       # slack already swamps this level
        rhs = 0
        for k in range(m):
            r = d - L[k]
            if 0 <= r <= a[k]:
                rhs += comb(a[k], r) << (a[k] - r)
            if rhs >= lhs:
                break
        if lhs > rhs:
            return False
    return True


def floor_at(m, num, den):
    """Largest candidate ratio refuted by (ARCH-eps)."""
    cands = sorted({Fraction(m + k + 1 + aa, m + k)
                    for k in range(m) for aa in range(m - k)})
    lo, hi, best = 0, len(cands) - 1, None
    while lo <= hi:
        mid = (lo + hi) // 2
        a = profile(m, cands[mid])
        ok = False if a is None else arch_eps_ok(m, a, num, den)
        if ok:
            hi = mid - 1
        else:
            best = cands[mid]
            lo = mid + 1
    return best


if __name__ == "__main__":
    PHI2INV = 0.3819660112501051          # 1/phi^2 = 1 - 1/phi
    print("epsilon = 2^{-beta m};  predicted critical beta* = 1/phi^2 = "
          f"{PHI2INV:.7f}")
    print()
    print("   m   beta      epsilon=2^-(beta m)   certified floor rho >")
    for m in (64, 128, 256):
        exact = floor_at(m, 0, 1)
        print(f"  {m:4d}   exact     0                    {float(exact):.6f}")
        for beta in (0.25, 0.34, 0.38, 0.42, 0.50, 0.70):
            e = int(beta * m)
            f = floor_at(m, 1, 2 ** e)
            tag = "" if f is None else f"{float(f):.6f}"
            print(f"  {m:4d}   {beta:.2f}      2^-{e:<4d}              "
                  f"{tag if f else 'COLLAPSED (nothing refuted)'}")
        print()
    print("""READING.  The floor is EXACTLY INTACT for beta above 1/phi^2 (at m=256,
beta = 0.42 and 0.50 both reproduce the exact-fairness value 1.588652) and
DEGRADES CONTINUOUSLY below it (0.38 -> 1.586441, 0.34 -> 1.579710,
0.25 -> 1.556391).  The onset sits right at the predicted 1/phi^2 = 0.38197.
It is a smooth degradation, not a cliff: the floor does not vanish, it erodes.
So the golden ratio appears twice -- delta* = 1/phi is the binding degree
fraction and 1 - delta* = 1/phi^2 is the critical slack exponent, the two
tied by the same quadratic delta^2 = 1 - delta.  Caveat: this is a statement
about the lower-bound METHOD, not a construction; it shows the bound is not
robust to slack, not that C*(epsilon) itself drops.""")
