#!/usr/bin/env python3
"""
lrc_compact_finite_klein_S328.py -- klein-2026-07-18-S328
Owner: work the finite-check route for the compact regime.

THE REDUCTION (this is the solid part). For a primitive covering 13-set write
      rho   = v_max / v_2nd      (boxeph's compactness parameter)
      sigma = v_max / v_min      (the spread)
Two PROVED handlers carve the plane:
  * sigma <= 13  =>  M >= 1/14   by THM-405: take t = 1/(14 v_min); then every speed has
                     v t = v/(14 v_min) in [1/14, sigma/14] subset [1/14, 13/14], so ||v t|| >= 1/14.
  * rho   >= 13  =>  THM-1007 (kind-pasteur: single-killer covering, unconditional) and its
                     lacunary-chain extension.
So the ONLY open part of the compact regime is the wedge
      rho < 13   AND   sigma > 13
i.e. the top two speeds are close together while the overall spread is large -- the big ratio sits
INSIDE the sorted list, not at the top. (The deep well {1..12,182} has rho = 15.2, so it is NOT compact
and is not in this wedge; boxeph's extremal 2*{1..12} u {13} has sigma = 12, so THM-405 already covers it.)

IS THE WEDGE FINITE?  Not a priori -- compactness bounds only the top ratio, so a family can have
unbounded spread via internal growth. But sampling the wedge and binning by spread gives min M RISING
with sigma:
      sigma  8-15  16-31  32-63  64-127  128-255  256-511  512-1023
      min M  .1053  .1071  .1176   .1538    .1867    .1843    .2014
The danger concentrates at SMALL spread, which would localize the finite check to bounded sigma.
HONEST CAVEAT: that is RANDOM sampling, which failed a positive control in MISTAKE-162 -- the TREND is
meaningful, the absolute minima are not. The localization is a hypothesis, not established. A powered
(constructive witness-first) hunt in the wedge was launched and did NOT finish inside the time budget,
so it is inconclusive rather than negative.
"""
from fractions import Fraction as Fr
def rho(V): V=sorted(V); return Fr(V[-1],V[-2])
def sigma(V): V=sorted(V); return Fr(V[-1],V[0])
def handler(V):
    """which PROVED result covers V, if any"""
    if sigma(V) <= 13: return "THM-405 (t = 1/(14 v_min))"
    if rho(V)   >= 13: return "THM-1007 (single-killer / lacunary)"
    return "OPEN WEDGE: rho<13 and sigma>13"
if __name__=="__main__":
    for nm,V in [("deep well",list(range(1,13))+[182]),
                 ("2*{1..12}u{13}",[2*k for k in range(1,13)]+[13]),
                 ("{1..12,100}",list(range(1,13))+[100]),
                 ("AP {1..13}",list(range(1,14)))]:
        print("%-16s rho=%-8s sigma=%-8s -> %s"%(nm,float(rho(V)),float(sigma(V)),handler(V)))
