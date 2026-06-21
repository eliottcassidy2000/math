#!/usr/bin/env python3
"""Does far-smoothing (the 'death kernel' on the miss-count N) decrease L_yK8=10q0+q3+10q6?
death(q)_t = q_t*(1-t/7) + q_{t+1}*(t+1)/7  (a far runner deletes a missed sector w.p. N/7).
Compute L_yK8(death(q)) - L_yK8(q) on actual miss-distributions; if always <=0, far smoothing
lowers gK8 (the decorrelated concentration extremality). Also iterate (r far runners)."""
import sys, itertools
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def sector_of(p): return int((p%1)*7)
def missdist(E):
    E=sorted(set(E)); b=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(0,7*e+1): b.add(F(m,7*e))
    b=sorted(b); q=[F(0)]*7
    for i in range(len(b)-1):
        x0,x1=b[i],b[i+1]
        if x1<=x0: continue
        t=7-len(set(sector_of(e*((x0+x1)/2)) for e in E))
        if t<=6: q[t]+=x1-x0
    return q
def death(q):
    q2=[F(0)]*7
    for t in range(7):
        q2[t]+=q[t]*(F(7-t,7))           # no deletion
        if t+1<=6: q2[t]+=q[t+1]*(F(t+1,7)) # deletion from t+1
    return q2
def Lyk8(q): return 10*q[0]+q[3]+10*q[6]
print("L_yK8(death(q)) - L_yK8(q) for various bases (death = one far-runner smoothing):")
bad=0
for E in [(0,1,2,3,4,5,6,7),(0,1,2,3,4,5,6,7,8),(0,2,4,6,8,10,12,14),(0,1,2,4,8),(0,1,2,3,4,5,6,7,8,9)]:
    q=missdist(E); dq=death(q); diff=Lyk8(dq)-Lyk8(q)
    print(f"  E={str(E):<28} q={[float(x) for x in q]}")
    print(f"     L_yK8={float(Lyk8(q)):.4f} -> death {float(Lyk8(dq)):.4f}  diff={float(diff):+.5f}  {'(decreases OK)' if diff<=0 else '(INCREASES!)'}")
    if diff>0: bad+=1
print(f"\n  death increases L_yK8 in {bad} cases. The exact condition: 10q1+4q4 <= 3q3+60q6 (from diff<=0).")
# check the algebraic condition
print("  diff = (10q1 - 3q3 + 4q4 - 60q6)/7 ; far smoothing lowers gK8 IFF 10q1+4q4<=3q3+60q6.")
