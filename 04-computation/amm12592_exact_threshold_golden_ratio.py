#!/usr/bin/env python3
"""Solve the archimedean threshold exactly (stationarity), test delta* = 1/phi.

At the threshold three conditions hold simultaneously:
  (1) d g/d x = 0        (interior max over the stratum coordinate x)
  (2) g = H(delta)       (capacity touches the requirement)
  (3) d g/d delta = H'(delta)   (tangency in delta)

With p = r/alpha, (3) reduces to (1-p)/p = 2(1-delta)/delta, i.e.

        p = delta / (2 - delta),

and (1) reduces, using log2((1-p)/p) = 1 + H'(delta), to

        gamma = -H'(delta) / [ H(p) + (1+H'(delta))(1-p) ].

Then (2) gives alpha = H(delta)/[H(p)+1-p], while the definition of r gives
alpha = (2-delta)/[ (1+gamma)/gamma - p ].  Equating the two determines
delta.
"""
from mpmath import mp, mpf, log, sqrt, findroot

mp.dps = 40
LOG2 = log(2)


def H(p):
    if p <= 0 or p >= 1:
        return mpf(0)
    return (-p * log(p) - (1 - p) * log(1 - p)) / LOG2


def Hp(d):                      # derivative of H
    return log((1 - d) / d) / LOG2


def residual(delta):
    delta = mpf(delta)
    p = delta / (2 - delta)
    hp = Hp(delta)
    gamma = -hp / (H(p) + (1 + hp) * (1 - p))
    alpha1 = H(delta) / (H(p) + 1 - p)
    alpha2 = (2 - delta) / ((1 + gamma) / gamma - p)
    return alpha1 - alpha2


if __name__ == "__main__":
    d = findroot(residual, mpf("0.618"))
    p = d / (2 - d)
    hp = Hp(d)
    gamma = -hp / (H(p) + (1 + hp) * (1 - p))
    alpha = H(d) / (H(p) + 1 - p)
    x = alpha / gamma - 1
    phi = (1 + sqrt(5)) / 2

    print(f"delta*  = {mp.nstr(d, 25)}")
    print(f"1/phi   = {mp.nstr(1/phi, 25)}")
    print(f"delta* - 1/phi = {mp.nstr(d - 1/phi, 10)}")
    print()
    print(f"p*      = {mp.nstr(p, 25)}")
    print(f"1/sqrt5 = {mp.nstr(1/sqrt(5), 25)}")
    print(f"p* - 1/sqrt5   = {mp.nstr(p - 1/sqrt(5), 10)}")
    print()
    print(f"gamma*  = {mp.nstr(gamma, 25)}")
    print(f"C_arch  = {mp.nstr(1 + gamma, 25)}")
    print(f"alpha*  = {mp.nstr(alpha, 25)}")
    print(f"x*      = {mp.nstr(x, 25)}")
    print()
    # if delta* = 1/phi were exact, p* would be 1/sqrt5 exactly
    print("closed-form candidate assuming delta* = 1/phi:")
    d2 = 1 / phi
    p2 = d2 / (2 - d2)
    hp2 = Hp(d2)
    g2 = -hp2 / (H(p2) + (1 + hp2) * (1 - p2))
    print(f"  p(1/phi) = {mp.nstr(p2, 25)}   1/sqrt5 = {mp.nstr(1/sqrt(5), 25)}")
    print(f"  gamma    = {mp.nstr(g2, 25)}   C = {mp.nstr(1+g2, 25)}")
    print(f"  residual at 1/phi = {mp.nstr(residual(d2), 10)}  (0 iff exact)")
