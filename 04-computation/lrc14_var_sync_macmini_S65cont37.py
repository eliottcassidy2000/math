#!/usr/bin/env python3
"""THM-713: the variance-synchronization decomposition. Var(N) = DIAG + COV where
DIAG = Sum_s A_s(1-A_s) (A_s = P(s empty)), COV = Sum_{s!=s'} [A_{ss'} - A_s A_s'].
Test: is consec's Var-max lead carried by COV (three-gap synchronization)?"""
from fractions import Fraction as F
import random
def stats(E):
    pts=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(e):
            for s in range(8): pts.add(F(m,e)+F(s,7*e))
    pts=sorted(x for x in pts if 0<=x<=1)
    A=[F(0)]*7; AA=[[F(0)]*7 for _ in range(7)]
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2; hit=set()
        for e in E:
            if e==0: hit.add(0); continue
            fr=e*mid-(e*mid).__floor__(); hit.add(int(fr*7))
        emp=[s for s in range(7) if s not in hit]
        for s in emp:
            A[s]+=b-a
            for s2 in emp:
                if s2!=s: AA[s][s2]+=b-a
    DIAG=sum(A[s]*(1-A[s]) for s in range(7))
    COV=sum(AA[s][s2]-A[s]*A[s2] for s in range(7) for s2 in range(7) if s2!=s)
    return DIAG,COV,DIAG+COV
k=9
fams=[("consec {1..9}",list(range(1,10))),("consec {0..8}",list(range(9))),
      ("best-mean coverer-ish",[2,4,5,6,7,8,10,12,14][:9]),
      ("near-AP",[1,2,3,4,5,6,7,8,26]),("mod7",[1,8,15,22,29,36,43,50,57]),
      ("spread",[1,3,7,12,19,25,31,38,44])]
random.seed(6)
for _ in range(12): fams.append(("rnd",sorted(random.sample(range(0,40),9))))
print(f"{'family':24s} {'DIAG':>8s} {'COV':>8s} {'Var(N)':>8s}")
rows=[]
for nm,E in fams:
    D,C,V=stats(E); rows.append((V,nm,D,C))
    print(f"{nm:24s} {float(D):8.4f} {float(C):8.4f} {float(V):8.4f}")
rows.sort(reverse=True)
print()
print(f"Var-max: {rows[0][1]} (V={float(rows[0][0]):.4f}, DIAG={float(rows[0][2]):.4f}, COV={float(rows[0][3]):.4f})")
dspread=max(float(r[2]) for r in rows)-min(float(r[2]) for r in rows)
cspread=max(float(r[3]) for r in rows)-min(float(r[3]) for r in rows)
print(f"DIAG spread across families: {dspread:.4f};  COV spread: {cspread:.4f}")
print(f"SYNCHRONIZATION VERDICT: consec's lead is {'COV-carried' if cspread>2*dspread else 'MIXED'}")
