#!/usr/bin/env python3
"""THM-703: the empty-moment majorant reduction. q(N) = 1 - 2N/3 + N(N-1)/6 majorizes
g(N) = [N=0] + [N=1]/3 on {0..7} (check: q(0)=1, q(1)=1/3, q(2)=0, q(3)=0, q(4)=1/3...).
Hence Phi = E[g(N_empty)] <= 1 - (2/3)m1 + (1/6)m2 pointwise-integrated: PROVED.
TEST: does the moment bound close (<= cap_{k+1}) at TRUE moment values across the landscape?
Where it closes, the global lemma reduces to the two-moment (pair-correlation) inequalities."""
from fractions import Fraction as F
import random, itertools
CAP={8:F(2243,5880),9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7),13:F(1)}
def pvec(E):
    pts=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(e):
            for s in range(8): pts.add(F(m,e)+F(s,7*e))
    pts=sorted(p for p in pts if 0<=p<=1); p=[F(0)]*8
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2; hit=set()
        for e in E:
            if e==0: hit.add(0); continue
            fr=e*mid-(e*mid).__floor__(); hit.add(int(fr*7))
        p[7-len(hit)]+=b-a
    return p
# majorant check on {0..7}
q=lambda N: 1-F(2,3)*N+F(1,6)*N*(N-1)
g=lambda N: 1 if N==0 else (F(1,3) if N==1 else 0)
assert all(q(N)>=g(N) for N in range(8)), [ (N,q(N)) for N in range(8)]
print("majorant q >= g on {0..7}: PROVED (min slack at N=2,3: q=0)")
print()
print("landscape: Phi_true vs moment-bound 1-(2/3)m1+(1/6)m2 vs cap_(k+1)")
random.seed(3)
rows=[]
def test(E,label):
    p=pvec(E); k=len(E)
    Phi=p[0]+p[1]/3; m1=sum(j*p[j] for j in range(8)); m2=sum(j*(j-1)*p[j] for j in range(8))
    mb=1-F(2,3)*m1+F(1,6)*m2; cap=CAP.get(k+1)
    closes = mb<=cap if cap else None
    print(f"  {label:28s} k={k}: Phi={float(Phi):.4f} mbound={float(mb):.4f} cap={float(cap):.4f} {'CLOSES' if closes else 'GAP +%.4f'%float(mb-cap)}")
    return closes
test(list(range(8)),"consec {0..7}")
test(list(range(1,9)),"consec {1..8}")
test([1,2,3,5,7,8,9,11],"argmin-mP-style k=8")
test(list(range(9)),"consec k=9")
test(list(range(10)),"consec k=10")
nclose=0; ntot=0
for _ in range(25):
    k=random.choice([8,9]); E=sorted(random.sample(range(0,30),k)); ntot+=1
    p=pvec(E); Phi=p[0]+p[1]/3
    m1=sum(j*p[j] for j in range(8)); m2=sum(j*(j-1)*p[j] for j in range(8))
    mb=1-F(2,3)*m1+F(1,6)*m2
    if mb<=CAP[k+1]: nclose+=1
print(f"  random spread<=30 k=8,9: {nclose}/{ntot} close at true moments")
