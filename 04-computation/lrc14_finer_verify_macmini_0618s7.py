#!/usr/bin/env python3
"""
lrc14_finer_verify_macmini_0618s7.py  (mac-mini-2026-06-18-S7)

Verify the FINER-COVER certificate route:
 (1) closure for the dangerous rows k=10,11 (L=14): M_L + C_L*W(consec) <= cap_k?
 (2) the universal W bound: W(E) <= W(consec_k) for all integer E (the gap-monotonicity).
     RAW version W_raw=Sum 1/|(e_j-e_i)(e_l-e_j)(e_l-e_i)| is term-by-term <= consec (PROVED,
     since e_j-e_i>=j-i). PRIMITIVE version W=Sum 1/h_primitive: check consec maximizes it.
 (3) W0 = (cap_k - M_L)/C_L vs W(consec): the closure margin (want W(consec) << W0).
"""
import sys, itertools, random
from math import gcd
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(618700)
SEV=F(1,7)
def measSL(E, L):
    E=sorted(set(E)); bps=set([F(0),F(1)]); ends=set()
    for j in range(L):
        ends.add(F(j,L)%1); ends.add((F(j,L)+SEV)%1)
    for e in E:
        if e==0: continue
        for t in ends:
            for m in range(e): bps.add((t+m)/e)
    bps=sorted(z for z in bps if 0<=z<=1); total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2; cov=[False]*L
        for e in E:
            p=(e*xm)%1; jhi=int(p*L); jlo=int((p-SEV)*L)
            for jj in range(jlo+1,jhi+1): cov[jj%L]=True
        if all(cov): total+=x1-x0
    return total
def Wprim(E):
    E=sorted(set(E)); w=0.0
    for i,j,l in itertools.combinations(range(len(E)),3):
        a,b,c=E[j]-E[l],E[l]-E[i],E[i]-E[j]
        g=gcd(gcd(abs(a),abs(b)),abs(c)) or 1
        h=abs((a//g)*(b//g)*(c//g))
        if h>0: w+=1.0/h
    return w
def Wraw(E):
    E=sorted(set(E)); w=0.0
    for i,j,l in itertools.combinations(range(len(E)),3):
        h=abs((E[j]-E[i])*(E[l]-E[j])*(E[l]-E[i]))
        if h>0: w+=1.0/h
    return w
def M_iid(k,L,trials=120000):
    import random as R; cnt=0; sev=1.0/7
    for _ in range(trials):
        pts=[0.0]+[R.random() for _ in range(k-1)]; cov=[False]*L
        for p in pts:
            jhi=int(p*L); jlo=int((p-sev)*L)
            for jj in range(jlo+1,jhi+1): cov[jj%L]=True
        if all(cov): cnt+=1
    return cnt/trials
H1=F(1,14)
def danger(u):
    iv=[]
    for j in range(u):
        c=F(j,u); a=(c-H1/u)%1; b=(c+H1/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def mgm(iv):
    iv=sorted(iv); o=[]
    for a,b in iv:
        if o and a<=o[-1][1]: o[-1]=(o[-1][0],max(o[-1][1],b))
        else: o.append((a,b))
    return o
def measGP(P):
    dz=mgm([iv for u in P for iv in danger(u)]); s=F(0); prev=F(0)
    for a,b in dz:
        if a>prev: s+=a-prev
        prev=max(prev,b)
    if prev<1: s+=1-prev
    return s
import functools
@functools.lru_cache(None)
def cap(k): return min(measGP(list(P)) for P in itertools.combinations(range(1,14),13-k))

print("(2) UNIVERSAL W bound: W(E) <= W(consec_k)?  (raw=PROVED term-by-term; primitive=check)")
for k in (8,9,10,11):
    Wrc=Wraw(list(range(k))); Wpc=Wprim(list(range(k)))
    vr=vp=0; mxp=Wpc
    for _ in range(3000):
        sp=random.randint(k-1,4*k); E=sorted(set([0,sp]+random.sample(range(1,sp),k-2)))
        if len(E)!=k: continue
        if Wraw(E)>Wrc+1e-9: vr+=1
        wp=Wprim(E)
        if wp>Wpc+1e-9: vp+=1; mxp=max(mxp,wp)
    print(f"  k={k}: W_raw(consec)={Wrc:.4f} (raw-violations {vr}); "
          f"W_prim(consec)={Wpc:.4f} (prim-violations {vp}, max seen {mxp:.4f})")

print("\n(1,3) finer-cover closure k=10,11 (L=14): bound = M_L + C_L*W(consec) vs cap_k")
for k in (10,11):
    capk=float(cap(k)); Wc=Wprim(list(range(k))); ML=M_iid(k,14)
    bank=[list(range(k))]
    for extra in itertools.combinations(range(1,k+4),k-1):
        E=[0]+list(extra)
        if max(E)<=k+3: bank.append(E)
    worst=0.0
    for E in bank:
        if len(set(E))!=k: continue
        corr=float(measSL(E,14))-ML; w=Wprim(E)
        if w>1e-9: worst=max(worst,corr/w)
    bound=ML+worst*Wc; W0=(capk-ML)/worst if worst>0 else 999
    print(f"  k={k}: M_L={ML:.4f} C_L*={worst:.5f} bound={bound:.4f} cap={capk:.4f} "
          f"{'CLOSES' if bound<=capk else 'no'}; W0={(capk-ML)/worst:.1f} vs W(consec)={Wc:.2f}")
print("\nDONE.")
