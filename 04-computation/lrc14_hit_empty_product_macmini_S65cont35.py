#!/usr/bin/env python3
"""THM-711: the hit-empty product form. J(E) := 6 m1 - m2 = E[N(7-N)] (exact identity:
E[N(7-N)] = 7m1 - E[N^2] = 7m1 - (m2 + m1) = 6m1 - m2). The k=9 deg-2 base holds iff
J >= 12(1 - cap_10) = 12*1611/4004... compute exact. HUNT: adversarial 9-sets, all classes,
for the global inf -- if inf > threshold, the base needs NO spread bound."""
from fractions import Fraction as F
import random, itertools
cap10=F(55,91)
THR=12*(1-cap10)
print(f"threshold 12(1-cap_10) = {THR} = {float(THR):.4f}")
def J(E):
    pts=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(e):
            for s in range(8): pts.add(F(m,e)+F(s,7*e))
    pts=sorted(x for x in pts if 0<=x<=1); tot=F(0)
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2; hit=set()
        for e in E:
            if e==0: hit.add(0); continue
            fr=e*mid-(e*mid).__floor__(); hit.add(int(fr*7))
        N=7-len(hit); tot+=(b-a)*N*(7-N)
    return tot
worst=(None,None)
def rec(E,label):
    global worst
    j=J(E)
    if worst[0] is None or j<worst[0]: worst=(j,list(E)+[label])
    return j
print(f"  consec {{0..8}}: J = {float(rec(list(range(9)),'consec0')):.4f}")
print(f"  consec {{1..9}}: J = {float(rec(list(range(1,10)),'consec1')):.4f}")
print(f"  doubling 12->24 style {{1..8,18}}: J = {float(rec([1,2,3,4,5,6,7,8,18],'dbl')):.4f}")
print(f"  near-AP {{1..8,26}}: J = {float(rec([1,2,3,4,5,6,7,8,26],'nearAP')):.4f}")
print(f"  mod7-aligned {{1,8,15,22,29,36,43,50,57}}: J = {float(rec([1,8,15,22,29,36,43,50,57],'mod7')):.4f}")
random.seed(9)
for _ in range(40):
    rec(sorted(random.sample(range(0,45),9)),'rnd')
for _ in range(12):
    base=sorted(random.sample(range(0,12),8)); far=random.randint(60,400)
    rec(base+[far],'farmix')
print(f"  +52 adversarial (random, far-mix): done")
print(f"  GLOBAL WORST: J = {worst[0]} = {float(worst[0]):.4f} at {worst[1]}")
print(f"  UNCONDITIONAL? worst >= {float(THR):.4f}: {'YES -- margin '+format(float(worst[0]-THR),'+.4f') if worst[0]>=THR else '*** NO -- VIOLATION ***'}")
