#!/usr/bin/env python3
"""CONCENTRATION via VARIANCE: q0+q6 = P(N in {0,6}) = mass at the EXTREMES of miss-count N.
N = sum_{j=1..6} 1[sector j empty]. consec clusters runners => sector-empty indicators POSITIVELY
correlated => high Var(N) => extreme N. Test: does consec maximize Var(N)? is q0+q6 monotone in
moments of N? Compute E[N], Var(N), q0+q6, L_yK8 for consec vs wide; look for the clean inequality."""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(3)
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
def moments(q):
    EN=sum(t*q[t] for t in range(7)); EN2=sum(t*t*q[t] for t in range(7))
    return EN, EN2-EN*EN  # mean, variance
def Lyk8(q): return 10*q[0]+q[3]+10*q[6]
print(f"{'config':<34}{'E[N]':>8}{'Var(N)':>9}{'q0+q6':>8}{'q3':>7}{'L_yK8':>8}")
def show(name,E):
    q=missdist(E); m,v=moments(q); print(f"{name:<34}{float(m):>8.4f}{float(v):>9.4f}{float(q[0]+q[6]):>8.4f}{float(q[3]):>7.4f}{float(Lyk8(q)):>8.4f}")
for k in (8,9,10):
    print(f"--- k={k} ---")
    show(f"consec_{k}", tuple(range(k)))
    show(f"dilated 2*consec", tuple(2*i for i in range(k)))
    show(f"doublet base+{15,16}", tuple(range(k-2))+(15,16))
    show(f"single-far +21", tuple(range(k-1))+(21,))
    show(f"wide AP step large", tuple(13*i for i in range(k)))
# Does consec MAXIMIZE Var(N)? E[N] is ~constant (=6*(6/7)^? decorr)? check
print("\nIs E[N] ~ invariant and Var(N) maximized by consec? (=> consec extremal by 2nd moment)")
for k in (9,):
    cons=tuple(range(k)); _,vc=moments(missdist(cons)); EN_c,_=moments(missdist(cons))
    maxv=vc; argmax='consec'; N=0; viol_var=0
    EN_range=[EN_c,EN_c]
    for _ in range(80):
        E=tuple(sorted([0]+random.sample(range(1,60),k-1))); m,v=moments(missdist(E)); N+=1
        EN_range[0]=min(EN_range[0],m); EN_range[1]=max(EN_range[1],m)
        if v>vc: viol_var+=1
        if v>maxv: maxv=v; argmax=str(E)
    print(f"  k={k}: Var(consec)={float(vc):.4f}; {viol_var}/{N} wide exceed Var(consec); max Var={float(maxv):.4f} @ {argmax[:40]}")
    print(f"        E[N] range over wide = [{float(EN_range[0]):.3f}, {float(EN_range[1]):.3f}] (consec E[N]={float(EN_c):.3f})")
