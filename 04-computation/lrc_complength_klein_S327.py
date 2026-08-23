#!/usr/bin/env python3
"""
lrc_complength_klein_S327.py -- klein-2026-07-18-S327
Owner: work the component-length distribution of the good set (the missing input named in S326).

THM-1042. Positive-length components of G_B(1/14) have endpoints
(14k +- 1)/(14v), v in B, so every positive component length is a rational
with denominator dividing 14*v*v' (e.g. 1/588, 588 = 14*6*7).  The count
below is the number of positive-length arcs; closed safe sets can additionally
have isolated equality-wall components.  Exact table for B = {1..k}:

   B        mu        #comps  L_max    1/L_max   next speed k+1
   {1..3}   0.69048     4     5/21       4.2        4
   {1..4}   0.61905     6     9/56       6.2        5
   {1..5}   0.50476    10     4/35       8.8        6
   {1..6}   0.45714    12     1/12      12.0        7
   {1..7}   0.33469    18     3/49      16.3        8
   {1..8}   0.26582    20     5/112     22.4        9
   {1..9}   0.18107    20     2/63      31.5       10
   {1..10}  0.13798    20     3/140     46.7       11
   {1..11}  0.05633    14     1/77      77.0       12

An additive step charges a PROPORTIONAL loss, valid only when components exceed the incoming arc period
1/w; otherwise one period spans a whole component. So base B admits w only if w > 1/L_max(B). But
1/L_max(k) > k+1 at EVERY row and the gap widens (ratios 1.05,1.24,1.47,1.71,2.04,2.49,3.15,4.25,6.42):
1/L_max grows superlinearly, the next speed linearly.

=> AN ADDITIVE CERTIFICATE CAN NEVER ABSORB A CONSECUTIVE SPEED, so it fails on every family with a
small-integer core (AP, deep well, GW, every covering family with a consecutive run). Verified on the deep
well: every initial split is blocked, and the blockers are exactly the consecutive speeds.

EXPLAINS UNIFORMLY: (1) THM-1015's large-killer thresholds 65.7..347.5 ARE the 1/L_max of its bases;
(2) the S326 measure-recursion died at w=8 because 8 < 16.3 = 1/L_max({1..7}); (3) the S314 radius-3
fragmentation wall is the same short components seen from the Hamming side. Changing the state variable
(interval -> measure+components) removes the r<7 cap but not this: both price a speed against component
scale. NEGATIVE about the METHOD, not about LRC(14).
"""
from fractions import Fraction as Fr
d=Fr(1,14)
def components(B):
    """exact components of G_B(1/14) as lengths (Fractions)"""
    cuts={Fr(0),Fr(1)}
    for v in B:
        for k in range(0,v+1):
            for e in (Fr(k,v)-d/v, Fr(k,v)+d/v):
                if 0<=e<=1: cuts.add(e)
    cs=sorted(cuts); segs=[]
    for i in range(len(cs)-1):
        a,b=cs[i],cs[i+1]
        if b<=a: continue
        mid=(a+b)/2
        if all(min((v*mid)%1,1-((v*mid)%1))>=d for v in B): segs.append([a,b])
    m=[]
    for a,b in segs:
        if m and m[-1][1]==a: m[-1][1]=b
        else: m.append([a,b])
    return [b-a for a,b in m]
def admits(B,w):
    """the additive step is valid only if w exceeds 1/L_max(B)"""
    L=components(B)
    return bool(L) and w > 1/max(L)
if __name__=="__main__":
    print(" k | mu       | N  | L_max   | 1/L_max | k+1 | blocked?")
    for k in range(3,12):
        B=list(range(1,k+1)); L=components(B)
        Lm=max(L)
        print(" %2d | %.5f | %2d | %-7s | %7.1f | %3d | %s"
              %(k,float(sum(L)),len(L),Lm,float(1/Lm),k+1,"YES" if (k+1)<=1/Lm else "no"))
