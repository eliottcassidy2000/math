#!/usr/bin/env python3
"""ADVERSARIAL verification of the gK8 wide-closure route (claude-opus HYP-2809/2812):
(1) gK8 dual validity: 10*p0 <= L_yK8=10q0+q3+10q6 (trivial, q3,q6>=0) -- confirm.
(2) finite check: L_yK8(consec_k) vs 10*cap_k -- the margin.
(3) is consec the L_yK8-MAX? broad adversarial search (wide, dilated, doublet, single-far)."""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(7)
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
def Lyk8(q): return 10*q[0]+q[3]+10*q[6]
caps={8:F(2243,5880),9:F(1979,4004),10:F(4,7),11:F(5,7),12:F(6,7)}
print("(2) L_yK8(consec_k) vs 10*cap_k:")
for k in (8,9,10,11,12):
    cons=tuple(range(k)); Lc=Lyk8(missdist(cons)); cap10=10*caps[k]
    print(f"  k={k}: L_yK8(consec)={float(Lc):.4f}  10*cap={float(cap10):.4f}  margin={float(cap10-Lc):.4f}  {'OK' if Lc<cap10 else 'FAIL'}")
print("\n(3) adversarial: any config with L_yK8 > L_yK8(consec_k)?")
for k in (8,9,10):
    cons=tuple(range(k)); Lc=Lyk8(missdist(cons)); viol=0; worst=F(0); N=0
    fams=[]
    # dilated consec, doublets, single-far, random wide, spread
    for d in (2,3,5): fams.append(tuple(d*i for i in range(k)))
    for M in (15,18,21,30): fams.append(tuple(range(k-2))+(M,M+1))
    for f in (15,21,50,100): fams.append(tuple(range(k-1))+(f,))
    for _ in range(60): fams.append(tuple(sorted([0]+random.sample(range(1,80),k-1))))
    for E in fams:
        L=Lyk8(missdist(E)); N+=1
        if L>Lc+F(1,10**6): viol+=1
        if L>worst: worst=L
    print(f"  k={k}: consec L_yK8={float(Lc):.4f}; {viol}/{N} configs exceed it; max found={float(worst):.4f} {'(consec is MAX, gK8 route sound)' if viol==0 else '(EXCEEDED!)'}")
