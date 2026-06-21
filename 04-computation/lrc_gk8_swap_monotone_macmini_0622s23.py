#!/usr/bin/env python3
"""The CORRECT concentration-extremality local move: SAME-SIZE swap (replace a base runner by a far
one). Does L_yK8=10q0+q3+10q6 DECREASE under base->far swap? If so + connectivity => bounded max.
Compare consec_k vs consec_{k-1} u {far}, and general swaps."""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(2)
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
def q0q6(q): return q[0]+q[6]
print("SAME-SIZE swap: consec_k vs consec_{k-1} u {far}; does L_yK8 (and q0+q6) decrease?")
for k in (8,9,10):
    base=tuple(range(k)); Lb=Lyk8(missdist(base)); qb=missdist(base)
    print(f"\n k={k}: consec_k={base}  L_yK8={float(Lb):.4f} q0+q6={float(q0q6(qb)):.4f}")
    inc_L=0; inc_q=0; tot=0
    for far in [15,16,18,21,30,50,100]:
        swap=tuple(range(k-1))+(far,)   # drop k-1, add far
        qs=missdist(swap); Ls=Lyk8(qs); tot+=1
        dL=Ls-Lb; dq=q0q6(qs)-q0q6(qb)
        if dL>0: inc_L+=1
        if dq>0: inc_q+=1
        print(f"    swap k-1->far={far:>4}: L_yK8={float(Ls):.4f} (d={float(dL):+.4f}) q0+q6={float(q0q6(qs)):.4f} (d={float(dq):+.4f}) {'L UP!' if dL>0 else ''}")
    print(f"    => swap raised L_yK8 in {inc_L}/{tot}, raised q0+q6 in {inc_q}/{tot} (want 0/0: bounded is max)")
# random general swaps (replace any base elt by far)
print("\nrandom same-size base->far swaps (does L_yK8 ever increase above consec?):")
viol=0;N=0
for k in (8,9):
    cons=tuple(range(k)); Lc=Lyk8(missdist(cons))
    for _ in range(40):
        # random size-k config with a far element
        base=sorted(random.sample(range(1,14), k-2)); far=random.choice([15,16,18,20,21,25,30])
        E=tuple([0]+base+[far]); L=Lyk8(missdist(E)); N+=1
        if L>Lc+F(1,10**6): viol+=1
print(f"  L_yK8(wide) > L_yK8(consec_k) in {viol}/{N} random wide configs (0 => consec is the L_yK8-max)")
