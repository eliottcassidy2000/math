#!/usr/bin/env python3
"""cont.39: the CLOSED-FORM lower bound J >= 6(1-p0). Since 0 in E hits sector 0 always,
N (empty among the 6 nonzero sectors) <= 6, and N(7-N) >= 6 for N>=1, = 0 for N=0.
So J = E[N(7-N)] >= 6*P(N>=1) = 6(1-p0), p0 = P(all 6 nonzero sectors hit).
CLAIM to test: p0 <= 19/91 uniformly => J >= 6*72/91 = 432/91 (the threshold) CLOSES k=9.
CRUX: low-mu (good-coverer) families have HIGH p0 -- do any exceed 19/91?"""
from fractions import Fraction as F
import random
def p0_and_J(E):
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
    p0=p[0]  # N=0 = all 7 sectors hit (incl 0). Since 0 always hits sector 0, N=0 <=> 6 nonzero all hit
    J=sum(p[n]*n*(7-n) for n in range(8))
    return p0, J
thr=F(432,91); p0max=F(19,91)
print(f"threshold J>=432/91={float(thr):.4f}; simple bound needs p0<=19/91={float(p0max):.4f}")
print(f"{'family':30s} {'p0':>8s} {'6(1-p0)':>9s} {'J':>8s}  bound-closes?")
def show(nm,E):
    p0,J=p0_and_J(E); lb=6*(1-p0)
    ok = lb>=thr
    print(f"  {nm:28s} {float(p0):8.4f} {float(lb):9.4f} {float(J):8.4f}  {'YES' if ok else '*** NO (p0=%0.4f>0.2088)'%float(p0)}")
    return p0,ok
show("consec {1..9}",list(range(1,10)))
show("consec {0..8}",list(range(9)))
show("mod7",[1,8,15,22,29,36,43,50,57])
# the crux: LOW-mu families = good coverers = high p0
random.seed(39)
worst_p0=(F(0),None); nbreak=0; ntot=0
def climb_lowmu():
    E=sorted(random.sample(range(0,30),9)); best=p0_and_J(E)[0]
    for _ in range(500):
        i=random.randrange(9); nd=random.choice([-2,-1,1,2])
        E2=E[:]; E2[i]+=nd
        if len(set(E2))!=9 or min(E2)<0 or max(E2)>60: continue
        v=p0_and_J(sorted(E2))[0]
        if v>best: best=v; E=sorted(E2)
    return E,best
print("\nMAX-p0 hill-climbs (the families that BREAK the simple bound if p0>19/91):")
for t in range(8):
    E,p0=climb_lowmu(); ntot+=1
    _,J=p0_and_J(E); closes = J>=thr
    br = p0>p0max
    if br: nbreak+=1
    if p0>worst_p0[0]: worst_p0=(p0,E)
    print(f"  run {t}: p0={float(p0):.4f} J={float(J):.4f}  {'BOUND BREAKS (p0>0.2088) but real J %s' % ('OK' if closes else 'FAILS!') if br else 'bound holds'}")
print(f"\nMAX p0 found: {worst_p0[0]}={float(worst_p0[0]):.4f} at {worst_p0[1]}")
print(f"simple bound p0<=19/91 holds uniformly: {worst_p0[0]<=p0max}  ({nbreak}/{ntot} runs exceed)")
