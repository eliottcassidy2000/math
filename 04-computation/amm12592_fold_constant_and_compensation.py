#!/usr/bin/env python3
"""Working the three objects of THM-3009 sec 11.

OBJECT 1 -- THE DEFICIENCY FOLD gives the convergence constant.
Near the threshold, Phi_gamma(delta) = max_x g - H(delta) is a fold:
Phi ~ -A(gamma*-gamma) + (B/2)(delta-delta*)^2.  (ARCH) fails iff the depth
A(gamma*-gamma) exceeds the Stirling error E, so the certified bound obeys

    gamma* - gamma_m  =  E / A.

By the envelope theorem A = dg/dgamma at the maximiser, and only alpha
depends on gamma at fixed ell:  alpha = gamma(2-ell)/(1+gamma), so

    A = (2-ell*)/(1+gamma*)^2 * (1 - log2(1-p*)),

using dg/dalpha = 1 - log2(1-p).  With the SHARP Stirling constant -- the
central term is 2^{nH}/sqrt(2 pi n p(1-p)) and the O(m) summands of (ARCH)
form a Gaussian of width ~sqrt(m), so the sum costs +1/2 log2 m against the
LHS's -1/2 log2 m, leaving E = (1/2) log2(m)/m -- the prediction is

    c = (gamma* - gamma_m) / (log2(m)/m) = 1/(2A).

The crude bound E = 2 log2(m)/m of THM-3009 sec 10.3 would give c = 2/A,
about five times too large; the measured c ~ 0.30-0.34 discriminates between
them.
"""
from mpmath import mp, mpf, log, sqrt

mp.dps = 40
L2 = log(2)
H = lambda p: mpf(0) if p <= 0 or p >= 1 else (-p * log(p) - (1 - p) * log(1 - p)) / L2


def data():
    phi = (1 + sqrt(5)) / 2
    d = 1 / phi
    p = d / (2 - d)
    F = log(5) / (2 * L2)
    g = log(phi) / L2 / F
    alpha = H(d) / (H(p) + 1 - p)
    ell = d - p * alpha
    return phi, d, p, g, alpha, ell


if __name__ == "__main__":
    phi, d, p, g, alpha, ell = data()
    dgda = 1 - log(1 - p) / L2                      # dg/dalpha
    A = (2 - ell) / (1 + g) ** 2 * dgda             # dPhi/dgamma at the max
    print("OBJECT 1 -- the fold constant")
    print(f"  ell*            = {mp.nstr(ell,18)}")
    print(f"  dg/dalpha       = 1 - log2(1-p*) = {mp.nstr(dgda,18)}")
    print(f"  A = dPhi/dgamma = {mp.nstr(A,18)}")
    print()
    print(f"  crude  E = 2 log2 m / m   ->  c = 2/A   = {mp.nstr(2/A,12)}")
    print(f"  sharp  E = (1/2)log2 m/m  ->  c = 1/(2A) = {mp.nstr(1/(2*A),12)}")
    print()
    meas = {256: "1.588652", 512: "1.592476", 1024: "1.594876",
            2048: "1.596225", 4096: "1.597001"}
    print("     m      measured c        (c_sharp - c_m)")
    cs = 1 / (2 * A)
    from math import log2 as l2
    for m in sorted(meas):
        gap = (1 + g) - mpf(meas[m])
        cm = gap / (mpf(l2(m)) / m)
        print(f"  {m:5d}      {mp.nstr(cm,10):14s}    {mp.nstr(cs-cm,8)}")
    print(f"\n  predicted c = 1/(2A) = {mp.nstr(cs,10)}")
    print("  measured c rises 0.2987 -> 0.3367 monotonically toward it;")
    print("  the residual is the O(1/log m) correction from the fold's width.")

    print("\nOBJECT 2 -- compensation depth (structural)")
    print("""  Let C*_B be the least slope over schemes balanced on ratio-B blocks
  (THM-3007: balanced blocks are exactly [2^a, 2^g), so B = 2^j).  Coarser
  balance is FEWER constraints, hence
        C*_1 >= C*_2 >= C*_4 >= ... >= C*  (global balance only).
  C_arch = log_5(5 phi^2) is a lower bound for C*_1 ONLY.  The whole gap
  between 1.598 and the true C* is the value of cross-shell compensation.
  First decisive test: is C*_2 < C*_1?  For a ratio-B block the middle
  composition range N <= k <= l carries FREE symmetric deficits
  delta_j = delta_{N+l-j} (the two branches couple), which the ratio-2 block
  does not have -- so the extra freedom is real and localised.""")

    print("\nOBJECT 3 -- syllable depth (verdict)")
    print("""  Refining the balance requirement from n_1 to (n_1,...,n_j) partitions
  the shells FINER, i.e. imposes strictly MORE constraints.  So depth-j
  balance can only raise the achievable slope, never lower it, and it is not
  a weaker hypothesis either -- it is a strengthening of exactly the
  assumption one wants to remove.  Dead end in both directions; its only use
  is as bookkeeping.  Recorded as a negative.""")
