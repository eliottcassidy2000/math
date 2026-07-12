#!/usr/bin/env python3
"""cont.41: extend klein THM-717 -- the (POS) bound 6T1+4T2+2T3 >= 4717/882.
POS = E[g(N)], g=(0,6,10,12,12,12,12,12) NON-DECREASING. Non-decreasing => POS is
MINIMIZED by stochastically-smallest N = best coverage = low mu. So the POS-minimizer
should be a LOW-mu (well-covering) family, NOT the saddle consec. TEST: does any low-mu
family have POS < 4717/882 (= correcting klein's 'extremal at consec')?"""
from fractions import Fraction as F
import random
def Tvec(E):
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
    T=[sum(p[n] for n in range(j,8)) for j in range(8)]  # T[j]=P(N>=j)
    mu=sum(j*p[j] for j in range(8))
    POS=6*T[1]+4*T[2]+2*T[3]
    J=sum(p[n]*n*(7-n) for n in range(8))
    return float(mu),POS,J,p
target=F(4717,882)
print(f"klein POS bound: 6T1+4T2+2T3 >= 4717/882 = {float(target):.4f} (extremal at consec, per THM-717)")
print(f"{'family':22s} {'mu':>6s} {'POS':>9s} {'J':>8s}  vs 4717/882")
def show(nm,E):
    mu,POS,J,p=Tvec(E)
    flag = "*** BELOW klein bound!" if POS<target else ""
    print(f"  {nm:20s} {mu:6.3f} {float(POS):9.4f} {float(J):8.4f}  {flag}")
    return POS
show("consec {1..9}",list(range(1,10)))
show("consec {0..8}",list(range(9)))
# low-mu families from cont.38 (good coverers) -- do they undercut POS?
show("lowmu A",[6,9,12,15,17,18,21,22,24])
show("lowmu B",[1,7,10,13,14,16,18,19,28])
show("lowmu C",[2,3,6,10,14,18,22,26,31])
# adversarial: MINIMIZE POS directly
random.seed(41)
def climb_pos():
    E=sorted(random.sample(range(0,32),9)); best=Tvec(E)[1]
    for _ in range(700):
        i=random.randrange(9); nd=random.choice([-2,-1,1,2])
        E2=E[:]; E2[i]+=nd
        if len(set(E2))!=9 or min(E2)<0 or max(E2)>60: continue
        v=Tvec(sorted(E2))[1]
        if v<best: best=v; E=sorted(E2)
    return E,best
print("\nadversarial MIN-POS hill-climbs (undercut 4717/882 => correct klein):")
gmin=(F(99),None)
for t in range(10):
    E,POS=climb_pos()
    if POS<gmin[0]: gmin=(POS,E)
    mu=Tvec(E)[0]
    print(f"  run {t}: POS={float(POS):.4f} mu={mu:.3f}  {'*** BELOW' if POS<target else 'ok'} at {E}")
print(f"\nGLOBAL min POS found: {float(gmin[0]):.4f} at {gmin[1]}  (klein bound 4717/882={float(target):.4f})")
print(f"klein POS bound HOLDS (min >= 4717/882): {gmin[0]>=target}  {'-- consec NOT the minimizer!' if gmin[1] not in ([1,2,3,4,5,6,7,8,9],list(range(9))) and gmin[0]<target else ''}")
