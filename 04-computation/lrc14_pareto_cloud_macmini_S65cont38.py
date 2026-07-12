#!/usr/bin/env python3
"""Slice 1: the (mu, Var) Pareto cloud at k=9. J = E[N(7-N)] = mu(7-mu) - Var
(mu=E[N]). Lower-J frontier = MAX Var per mu. Question: does a one-parameter family
(consec + near-AP deformations) span the frontier => finite-dimensional extremal."""
from fractions import Fraction as F
import random
def muVarJ(E):
    pts=set([F(0),F(1)])
    for e in E:
        if e==0: continue
        for m in range(e):
            for s in range(8): pts.add(F(m,e)+F(s,7*e))
    pts=sorted(x for x in pts if 0<=x<=1)
    p=[F(0)]*8
    for a,b in zip(pts,pts[1:]):
        mid=(a+b)/2; hit=set()
        for e in E:
            if e==0: hit.add(0); continue
            fr=e*mid-(e*mid).__floor__(); hit.add(int(fr*7))
        p[7-len(hit)]+=b-a
    mu=sum(j*p[j] for j in range(8))
    m2=sum(j*(j-1)*p[j] for j in range(8))
    EN2=sum(j*j*p[j] for j in range(8))
    Var=EN2-mu*mu
    J=6*mu-m2   # = mu(7-mu)-Var
    return float(mu),float(Var),float(J),J
random.seed(38)
cloud=[]
# structured: consec-shifts, near-APs (one detune), mod-7, dilated APs
structured=[]
for sh in range(0,20): structured.append(("cshift%d"%sh,list(range(sh,sh+9))))
for det in range(1,15):
    base=list(range(1,9)); structured.append(("nearAP+%d"%det,base+[8+det]))
for step in [7,14]:
    structured.append(("mod%d"%step,[1+step*i for i in range(9)]))
for d in [2,3]:
    structured.append(("dilAP%d"%d,[d*i for i in range(1,10)]))
for nm,E in structured:
    mu,Va,J,Jx=muVarJ(E); cloud.append((mu,Va,J,nm,E))
# random spread
for _ in range(400):
    E=sorted(random.sample(range(0,50),9)); mu,Va,J,Jx=muVarJ(E); cloud.append((mu,Va,J,"rnd",E))
# global J-min
Jmin=min(cloud,key=lambda r:r[2])
print(f"GLOBAL J-min: {Jmin[3]} {Jmin[4]}  J={Jmin[2]:.4f} (mu={Jmin[0]:.3f}, Var={Jmin[1]:.3f})")
# frontier: max Var per mu-bin
print("\nLOWER-J FRONTIER (max Var per mu-bin => min J):")
bins={}
for mu,Va,J,nm,E in cloud:
    b=round(mu,1)
    if b not in bins or Va>bins[b][1]: bins[b]=(mu,Va,J,nm,E)
for b in sorted(bins):
    mu,Va,J,nm,E=bins[b]
    print(f"  mu~{b}: Var_max={Va:.4f} J={J:.4f}  {nm} {E}")
# which family types populate the frontier?
print("\nfrontier family types:", sorted(set(bins[b][3] for b in bins)))
# consec vs the frontier at each mu it hits
print("\nis consec ON the frontier? (consec (mu,Var) vs frontier max at that mu-bin)")
for sh in [0,1,2]:
    mu,Va,J,Jx=muVarJ(list(range(sh,sh+9)))
    b=round(mu,1); fmax=bins.get(b,(0,0))[1]
    print(f"  cshift{sh}: mu={mu:.3f} Var={Va:.4f} vs frontier {fmax:.4f}  {'ON' if abs(Va-fmax)<1e-9 else 'below by %.4f'%(fmax-Va)}")
