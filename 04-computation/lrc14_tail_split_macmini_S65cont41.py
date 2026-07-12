#!/usr/bin/env python3
"""cont.41 extension of klein THM-717: does POS = 6T1+4T2+2T3 SPLIT? Test whether
T1, T2, T3 are EACH individually minimized at consec. If so, POS >= 6T1min+4T2min+2T3min
and if all mins at consec => POS >= POS(consec) = 4717/882, PROVING klein's bound from
three separate (single-tail) covering statements. T1 = P(N>=1) = 1-p0 = meas(S7)."""
from fractions import Fraction as F
import random
def Ts(E):
    pts=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(e):
            for s in range(8): pts.add(F(m,e)+F(s,7*e))
    pts=sorted(x for x in pts if 0<=x<=1); p=[F(0)]*8
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2; hit=set()
        for e in E:
            if e==0: hit.add(0); continue
            fr=e*mid-(e*mid).__floor__(); hit.add(int(fr*7))
        p[7-len(hit)]+=b-a
    return [sum(p[n] for n in range(j,8)) for j in range(4)]  # T0,T1,T2,T3
cons=list(range(1,10)); Tc=Ts(cons)
print(f"consec {{1..9}} tails: T1={float(Tc[1]):.5f} T2={float(Tc[2]):.5f} T3={float(Tc[3]):.5f}")
print(f"  POS(consec) = 6T1+4T2+2T3 = {float(6*Tc[1]+4*Tc[2]+2*Tc[3]):.5f} = 4717/882={float(F(4717,882)):.5f}")
print("\nQ: is EACH T_j individually minimized at consec? (adversarial MIN per tail)")
random.seed(411)
def climb(jt):
    E=sorted(random.sample(range(0,30),9)); best=Ts(E)[jt]
    for _ in range(500):
        i=random.randrange(9); nd=random.choice([-2,-1,1,2])
        E2=E[:]; E2[i]+=nd
        if len(set(E2))!=9 or min(E2)<0 or max(E2)>55: continue
        v=Ts(sorted(E2))[jt]
        if v<best: best=v; E=sorted(E2)
    return E,best
for jt in [1,2,3]:
    gmin=(F(9),None)
    for t in range(6):
        E,v=climb(jt)
        if v<gmin[0]: gmin=(v,E)
    below = gmin[0]<Tc[jt]
    print(f"  T{jt}: consec={float(Tc[jt]):.5f}, adversarial min={float(gmin[0]):.5f} at {gmin[1]}")
    print(f"       => consec {'IS' if not below else 'is NOT'} the T{jt}-minimizer  {'(splits!)' if not below else '(*** T%d min elsewhere: %s)'%(jt,gmin[1])}")
# structural: is T1 = meas(S7) minimized at consec? (the known best-coverer question)
print("\nT1 = meas(S7) = 1 - p0: consec maximizes p0 (best coverer)?")
p0c = 1-float(Tc[1]); print(f"  p0(consec) = {p0c:.5f}")
