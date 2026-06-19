#!/usr/bin/env python3
"""
lrc14_finer_certificate_v2_macmini_0618s7.py  (mac-mini-2026-06-18-S7)

CLEAN test: does the FINER-COVER (L arcs of length 1/7) close the crude certificate
meas(S_L(E)) <= cap_k WITHOUT the extremal/finite-residual that S7 needed?
Fixes v1: TRUE iid main term M_L via Monte-Carlo (k-1 uniform points, all L arcs hit);
THOROUGH bank including the perforated shapes that gave S7's worst ratio.
Report: M_L, worst C_L*=max corr_L/W over bank, bound=M_L+C_L*W(consec), cap_k, CLOSES?
AND the DIRECT check max_E meas(S_L(E)) <= cap_k + whether consec maximizes meas(S_L)
over the thorough bank (codex found AP-beaters for S7 at k=12,13 — check finer L).
"""
import sys, itertools, random
from math import gcd
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(61870)
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

def M_iid(k, L, trials=200000):
    """true iid main: prob that k-1 uniform points (plus 0) cover all L arcs of length 1/7."""
    import random as R; cnt=0; sev=1.0/7
    for _ in range(trials):
        pts=[0.0]+[R.random() for _ in range(k-1)]
        cov=[False]*L
        for p in pts:
            jhi=int(p*L); jlo=int((p-sev)*L)
            for jj in range(jlo+1,jhi+1): cov[jj%L]=True
        if all(cov): cnt+=1
    return cnt/trials

def W(E):
    E=sorted(set(E)); w=0.0
    for i,j,l in itertools.combinations(range(len(E)),3):
        a,b,c=E[j]-E[l],E[l]-E[i],E[i]-E[j]
        g=gcd(gcd(abs(a),abs(b)),abs(c)) or 1
        h=abs((a//g)*(b//g)*(c//g))
        if h>0: w+=1.0/h
    return w
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

def thorough_bank(k):
    bank=[list(range(k))]
    # all single + double perforations within {0..k+3}
    for extra in itertools.combinations(range(1,k+4), k-1):
        E=[0]+list(extra)
        if max(E)-min(E)<=k+3: bank.append(E)
    return [E for E in bank if len(set(E))==k]

for k in (8,9):
    capk=float(cap(k)); Wc=W(list(range(k)))
    bank=thorough_bank(k)
    print(f"\n===== k={k}: cap_k={capk:.4f}, W(consec)={Wc:.3f}, |bank|={len(bank)} =====")
    for L in (7,14,21):
        ML=M_iid(k,L)
        worst=0.0; wargE=None; mx=(measSL(list(range(k)),L),'consec'); consmax=True
        for E in bank:
            sL=float(measSL(E,L)); w=W(E)
            if sL>float(mx[0])+1e-12: mx=(sL,E); consmax=False
            corr=sL-ML
            if w>1e-9 and corr/w>worst: worst=corr/w; wargE=E
        bound=ML+worst*Wc
        print(f"  L={L:2d}: M_L(iid)~{ML:.4f}  worstC_L={worst:.5f}  bound={bound:.4f}  cap={capk:.4f}  "
              f"{'CLOSES (cert-only!)' if bound<=capk else 'no (needs extremal)'} | "
              f"max meas(S_L)={float(mx[0]):.4f} {'(consec)' if consmax else 'BEATEN by '+str(mx[1])}")
print("\nDONE.")
