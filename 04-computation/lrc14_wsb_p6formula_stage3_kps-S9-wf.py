#!/usr/bin/env python3
"""
lrc14_wsb_p6formula_stage3_kps-S9-wf.py  (kind-pasteur-2026-06-19-S9)

Stage 3.  Stage 2 found the clean exact formula  p_6(AP_k) = 1/(7(k-1))  (verified k=7..14):
the probability that ALL six sectors {1..6} are empty for consec_k.  And L_y(k=8) =
P(N=0) + (1/10)P(N=3) + P(N=6), with g(0)=g(6)=1 the two unit-weight terms.

GOALS
 (G1) PROVE/EXPLAIN p_6(AP_k)=1/(7(k-1)).  N=6 means sectors {1..6} all empty, i.e. ALL of
      frac(e x), e=1..k-1, land in sector 0 = [0,1/7).  Equivalently in theta=7x coords:
      floor(e theta) mod 7 == 0 for all e=1..k-1, with theta in [0,7).  Characterize the theta-set.
 (G2) Is p_6 MAXIMIZED by consec over same-k sets?  And p_0?  And p_0+p_6?  (the two g=1 terms)
      Exhaustive k=8 bank.  This is the "two-extreme-events" sub-question.
 (G3) N=6 combinatorial meaning for GENERAL E: all e in E\{0} have frac(e x) in [0,1/7).
      = {x : ||e x|| has the same integer part region} -- this is exactly the LRC 'almost lonely'
      cluster event!  meas{x: all frac(e_i x) in [0,1/7)} -- relate to the lonely runner measure
      directly.  Compute and compare to 1/(7(k-1)) bound;  is consec the max?
 (G4) The map: P(N=6)(E) = meas{x: e x mod 1 in [0,1/7) for all e in E}.  By scaling y=7x... no,
      sectors are fixed.  Actually P(N=6) = meas{x: max_e frac(e x) < 1/7 AND min... } -- it's the
      measure of the set where the WHOLE cluster orbit is confined to the first 1/7 arc.  This is a
      classic 'all points in a short arc' = a SMALL-DENOMINATOR / continued-fraction event.
"""
import sys, itertools, functools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)

def law_N(E):
    E = sorted(set(E)); bps = set([F(0), F(1)])
    for e in E:
        if e == 0: continue
        for j in range(7):
            c = F(j, 7)
            for m in range(e): bps.add((c + m) / e)
    bps = sorted(z for z in bps if F(0) <= z < F(1))
    p = [F(0)] * 7
    for i in range(len(bps)):
        x0 = bps[i]; x1 = bps[i + 1] if i + 1 < len(bps) else F(1)
        if x1 <= x0: continue
        xm = (x0 + x1) / 2
        hit = set()
        for e in E:
            fr = (e * xm) % 1
            hit.add((fr.numerator * 7) // fr.denominator)
        nempty = len([j for j in range(1, 7) if j not in hit])
        p[nempty] += x1 - x0
    return p

def measAllInArc(E, a, b):
    """meas{ x in [0,1): frac(e x) in [a,b) for all e in E\{0} }."""
    E = sorted(set(e for e in E if e != 0))
    bps = set([F(0), F(1)])
    for e in E:
        for t in (a, b):
            for m in range(e): bps.add((t + m) / e)
    bps = sorted(z for z in bps if F(0) <= z < F(1))
    tot = F(0)
    for i in range(len(bps)):
        x0 = bps[i]; x1 = bps[i+1] if i+1<len(bps) else F(1)
        if x1<=x0: continue
        xm=(x0+x1)/2
        if all(a <= (e*xm)%1 < b for e in E): tot += x1-x0
    return tot

print("="*78)
print("(G1) characterize the theta-set for N=6 (all in sector 0) for consec_k; verify p6=1/(7(k-1))")
print("="*78)
# N=6 for consec_k: floor(e theta) mod 7 == 0 for all e=1..k-1, theta in [0,7).
# In x coords: frac(e x) in [0,1/7) for all e=1..k-1.
for k in (7,8,9,10):
    m = measAllInArc(list(range(k)), F(0), F(1,7))
    p6 = law_N(list(range(k)))[6]
    print(f"  k={k}: meas{{all frac(e x)<1/7, e=1..{k-1}}}={m}  p6(law)={p6}  equal={m==p6}  "
          f"=1/(7(k-1))? {m==F(1,7*(k-1))}")
# WHY 1/(7(k-1)): frac((k-1)x)<1/7 is the binding constraint; for x in [0, 1/(7(k-1))) all
# frac(e x)=e x < (k-1)x < 1/7.  But there are also other intervals near rationals.  Let's list them.
print("\n  the N=6 set for consec_8 (all frac(e x)<1/7, e=1..7) as a union of intervals:")
def arc_intervals(E,a,b):
    E=sorted(set(e for e in E if e!=0)); bps=set([F(0),F(1)])
    for e in E:
        for t in (a,b):
            for mm in range(e): bps.add((t+mm)/e)
    bps=sorted(z for z in bps if F(0)<=z<F(1)); out=[]
    for i in range(len(bps)):
        x0=bps[i]; x1=bps[i+1] if i+1<len(bps) else F(1)
        if x1<=x0: continue
        xm=(x0+x1)/2
        if all(a<=(e*xm)%1<b for e in E): out.append((x0,x1))
    # merge
    out.sort(); mg=[]
    for a2,b2 in out:
        if mg and a2<=mg[-1][1]: mg[-1]=(mg[-1][0],max(mg[-1][1],b2))
        else: mg.append((a2,b2))
    return mg
ivs = arc_intervals(list(range(8)),F(0),F(1,7))
for a2,b2 in ivs: print(f"    [{a2}, {b2})  len={b2-a2}")
print(f"    total={sum(b2-a2 for a2,b2 in ivs)} = 1/(7*7)={F(1,49)}")

print()
print("="*78)
print("(G2) does consec MAXIMIZE p_0, p_6, p_0+p_6 over same-k bank? (the g=1 terms)")
print("="*78)
k=8; pc=law_N(list(range(k)))
bank=[[0]+list(r) for r in itertools.combinations(range(1,13),k-1)]
b0=b6=b06=0
maxp0=maxp6=maxp06=(F(0),None)
for E in bank:
    p=law_N(E)
    if p[0]>pc[0]+F(1,10**12): b0+=1
    if p[6]>pc[6]+F(1,10**12): b6+=1
    if p[0]+p[6]>pc[0]+pc[6]+F(1,10**12): b06+=1
    if p[0]>maxp0[0]: maxp0=(p[0],E)
    if p[6]>maxp6[0]: maxp6=(p[6],E)
    if p[0]+p[6]>maxp06[0]: maxp06=(p[0]+p[6],E)
print(f"  consec p0={float(pc[0]):.5f} p6={float(pc[6]):.5f} p0+p6={float(pc[0]+pc[6]):.5f}")
print(f"  shapes beating consec:  p0:{b0}   p6:{b6}   p0+p6:{b06}")
print(f"  global max p0={float(maxp0[0]):.5f} at {maxp0[1]}")
print(f"  global max p6={float(maxp6[0]):.5f} at {maxp6[1]}")
print(f"  global max p0+p6={float(maxp06[0]):.5f} at {maxp06[1]}")

print()
print("="*78)
print("(G3) P(N=6)(E) = meas{all frac(e x) in [0,1/7)} -- is consec the max over same-k?")
print("="*78)
# already covered by p6 above; confirm via direct arc measure and find the maximizer over a WIDER box
k=8
best=(F(0),None); over_consec=0
consec_p6 = measAllInArc(list(range(k)),F(0),F(1,7))
for r in itertools.combinations(range(1,15),k-1):
    E=[0]+list(r)
    m=measAllInArc(E,F(0),F(1,7))
    if m>best[0]: best=(m,E)
    if m>consec_p6+F(1,10**12): over_consec+=1
print(f"  consec_8 P(N=6)={consec_p6}={float(consec_p6):.6f}")
print(f"  over a wider box maxE<=14: max P(N=6)={best[0]}={float(best[0]):.6f} at {best[1]}")
print(f"  shapes with P(N=6)>consec: {over_consec}  (0 => consec maximizes the all-in-1/7-arc measure)")

print("\nDONE stage 3.")
