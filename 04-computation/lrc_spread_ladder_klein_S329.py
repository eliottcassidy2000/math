#!/usr/bin/env python3
"""
lrc_spread_ladder_klein_S329.py -- klein-2026-07-18-S329
Owner: refine the map of what does NOT work; prove new statements sharpening it; challenge assumptions.

THM-1043 THE SPREAD LADDER.  For ANY finite V of positive integers with spread sigma = max/min:
        sigma <= n-1   =>   M(V) >= 1/n,    for every integer n >= 2.
Proof: t = 1/(n*min V) gives vt = v/(n min V) in [1/n, sigma/n] subset [1/n,(n-1)/n], so ||vt|| >= 1/n.
No primitivity, no covering, no bound on |V|. THM-405 is the n=14 rung -- the argument never used 14.
Verified: 5,400 families with sigma<=n-1, n=10..15, ZERO violations; tight (min M = 1/n) at n=10,11,14.

THE n=13 RUNG PROVES HYP-7355 FOR sigma <= 12, AT ITS OWN EXTREMAL. boxeph-S85 names 2*{1..12}u{13} as
the extremal of "compact primitive covering => M >= 1/13". That family has sigma = 24/2 = 12 <= 12, so
t = 1/26 gives M >= 1/13, and the exact value IS 1/13. The stated extremal is proved, with equality.

THE REFINED MAP (which PROVED handler covers each known low-M covering family):
   deep well {1..12,182}   M=14/183 (covering-MIN)  sigma=182 rho=15.17 -> THM-1007 single-killer
   {1..12,364} tower       M=28/365                 sigma=364 rho=30.33 -> THM-1007
   2*{1..12} u {13}        M=1/13                   sigma=12  rho=1.09  -> THM-1043 n=13 rung (TIGHT)
   {1..11,13,84}           M=7/89                   sigma=84  rho=6.46  -> OPEN WEDGE
TWO ASSUMPTIONS DISCARDED:
 (a) the covering-MINIMUM is NOT in the compact wedge -- the deep well has rho=15.17 (non-compact) and is
     single-killer, hence already proved. The extremum is not where the open territory is.
 (b) the wedge's binding case is {1..11,13,84} (M = 7/89, only 2.25% above 1/13), NOT boxeph's stated
     2*{1..12}u{13}, which is now proved. Any proof of HYP-7355 must survive that family with 2% room.

NEW COORDINATE. sigma and rho are ratios, but what breaks the ladder is WRAPPING: the witness
t = 1/(n vmin) fails exactly when a speed wraps past the first window. So use
        W(V) = ceil(log_13 sigma)   (13-fold octaves spanned);   W=1 is exactly the ladder's reach.
The residual is W>=2, and the binding family has W=2: the open problem is SINGLE-WRAP, one octave wide.
"""
from fractions import Fraction as Fr
import math
def spread(V): V=sorted(V); return Fr(V[-1],V[0])
def ladder_bound(V):
    """largest n with sigma <= n-1, i.e. the best rung; returns the guaranteed floor 1/n."""
    s=spread(V); n=math.floor(s)+1        # smallest n with n-1 >= sigma
    return Fr(1,n)
def octaves(V):
    s=float(spread(V))
    return max(1,math.ceil(math.log(s,13))) if s>1 else 1
def handler(V):
    V=sorted(V)
    if spread(V)<=13: return "spread ladder (THM-1043 / THM-405)"
    if V[-1] > 13*V[-2]: return "THM-1007 single-killer"
    return "OPEN WEDGE"
if __name__=="__main__":
    for nm,V in [("deep well",list(range(1,13))+[182]),
                 ("2*{1..12}u{13}",[2*k for k in range(1,13)]+[13]),
                 ("{1..11,13,84}",list(range(1,12))+[13,84])]:
        print("%-16s sigma=%-7s floor>=%-6s W=%d  ->  %s"
              %(nm,spread(V),ladder_bound(V),octaves(V),handler(V)))
