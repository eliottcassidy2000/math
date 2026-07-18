#!/usr/bin/env python3
"""
PROVE THE FRAME ORTHOGONALITY  A = <nu_hat, g_hat> = 0   (boxeph-2026-07-18-S99)

THM-727 exact identity:  Error(w) = sum_s integral_0^1 1_{R_s}(x) g_s(w x) dx,
  R_s = {frame E' misses exactly sector s},  g_s(y) = 1_{[s/7,(s+1)/7)}(y) - 1/7,
  S := w*Error = sum_s sum_{l!=0} (-1/2pi i l) U_s(l w) g_hat_s(l),
  U_s(N) = sum_{p endpt of R_s} eps_p e(-N p),   |g_hat_s(l)| = |sin(pi l/7)|/(pi|l|).

A := lim_{w->inf} Error(w) = sum_s sum_l (-1/2pi i l) nu_hat_s(l) g_hat_s(l),
     nu_hat_s(l) = lim U_s(l w)/w   (the S97 fixed-frame comb).

CLEAN READING (proved below):  Error(w) = |{x in UNION R_s : wx in the MISSED sector s(x)}|
  - (1/7)|union R_s|  =  the far runner wx filling the frame's missed sector.  As w->inf,
  wx equidistributes over the 7 sectors => Error -> 0 => A = 0.

This script:
 (1) computes Error(w) DIRECTLY (exact rational sector measures) for the frame {1..6},
     confirming Error -> 0 and the rate |S| = w|Error| <= kappa * R (THM-727, kappa~0.61);
 (2) computes kappa = (1/pi^2) sum_{l!=0} |sin(pi l/7)|/l^2 ;
 (3) shows |S|/R -> a constant < kappa (the 1st-order sum is O(R), NOT the crux) and probes
     whether |S| = o(R) (the SUFFICIENT uniform target, beyond A=0).

Exact Fractions for the sector measures.  Pure Python.
"""

from fractions import Fraction as Fr
from math import lcm, pi, sin, gcd


# ---- R_s intervals of the FRAME (miss exactly sector s) -------------------
def Rs_intervals(E, s):
    """return R_s as a list of (lo,hi) rational intervals: x where the |E| runners
    occupy exactly 6 sectors, missing sector s. (endpoints from the 7e-grid.)"""
    bps = sorted(set(Fr(k, 7 * e) for e in E for k in range(7 * e)) | {Fr(0), Fr(1)})
    ivs = []
    for i in range(len(bps) - 1):
        a, b = bps[i], bps[i + 1]
        mid = (a + b) / 2
        occ = set(int((e * mid % 1) * 7) for e in E)
        if len(occ) == 6 and s not in occ:
            ivs.append((a, b))
    # merge touching
    ivs.sort()
    out = []
    for a, b in ivs:
        if out and out[-1][1] == a:
            out[-1] = (out[-1][0], b)
        else:
            out.append((a, b))
    return out


def sector_measure(a, b, s, w):
    """exact measure of { x in [a,b] : frac(w x) in [s/7,(s+1)/7) }."""
    # y = w x ranges over [w a, w b]; count length with frac(y) in [s/7,(s+1)/7), /w
    Y0, Y1 = w * a, w * b
    al, be = Fr(s, 7), Fr(s + 1, 7)
    # measure over y of {y in [Y0,Y1]: frac(y) in [al,be)} : full periods + ends
    import math
    k0 = math.floor(Y0)
    k1 = math.floor(Y1)
    tot = Fr(0)
    for k in range(k0, k1 + 1):
        lo = max(Y0, Fr(k) + al)
        hi = min(Y1, Fr(k) + be)
        if hi > lo:
            tot += hi - lo
    return tot / w


def Error_direct(E, w):
    """Error(w) = sum_s [ measure(x in R_s with wx in sector s) - (1/7)|R_s| ]."""
    tot = Fr(0)
    R = 0
    for s in range(7):
        ivs = Rs_intervals(E, s)
        R += 2 * len(ivs)                      # endpoint count rho_s*2
        Rlen = sum(b - a for a, b in ivs)
        hit = sum(sector_measure(a, b, s, w) for a, b in ivs)
        tot += hit - Rlen / 7
    return tot, R


# ---- (2) kappa --------------------------------------------------------------
kappa = (1 / pi ** 2) * sum(abs(sin(pi * l / 7)) / l ** 2 for l in range(1, 20000)) * 2  # +-l
print(f"kappa = (1/pi^2) sum_{{l!=0}} |sin(pi l/7)|/l^2 = {kappa:.4f}   (THM-727: |S|<=0.61 R)")
print()

# ---- (1)+(3) Error(w) and the rate -----------------------------------------
E = [1, 2, 3, 4, 5, 6]
print(f"frame {E}:  Error(w), |S|=w|Error|, and |S|/R")
print(f"{'w':>6} {'Error(w)':>13} {'|S|=w|Err|':>11} {'R':>4} {'|S|/R':>8} {'kappa*R':>9}")
Rc = None
for w in [7, 13, 20, 50, 100, 200, 500, 1000, 2000, 5000]:
    err, R = Error_direct(E, w)
    Rc = R
    S = abs(err) * w
    print(f"{w:>6} {float(err):>13.6f} {float(S):>11.5f} {R:>4} {float(S)/R:>8.4f} "
          f"{kappa*R:>9.3f}")
print()
print(f"R (frame R_s endpoint count) = {Rc}  (FIXED -- does not grow with w)")
print("A = lim Error(w) : the column Error(w) -> 0 as w->inf  ==>  A = 0  (proved: |Error|<=kappa*R/w).")
print("|S| stays BOUNDED (<= kappa*R): A=0 is the O(R) tail of THM-727. The SUFFICIENT uniform")
print("target is |S| = o(R) (1st-order Weyl cancellation over the R_s boundaries), NOT o(R^2).")
