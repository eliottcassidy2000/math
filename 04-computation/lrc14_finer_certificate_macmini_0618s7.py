#!/usr/bin/env python3
"""
lrc14_finer_certificate_macmini_0618s7.py  (mac-mini-2026-06-18-S7)

Does the FINER-COVER certificate close where S7 did not?
Certificate: meas(S_L(E)) = M_L(k) + corr_L(E), corr_L(E) <= C_L * W(E) (relation-height split,
same W=Sum_triples 1/height since same offsets). Then meas(S_L) <= M_L + C_L*W(consec_k) for
ALL E (consec maximizes W). CLOSES iff M_L + C_L*W(consec_k) <= cap_k.
S7 FAILED: M_7 + C_7**W(consec_8) = 0.0245 + 0.359 = 0.384 > cap_8=0.381.
Test L=14,21,28 for k=8,9,10: M_L (dissoc proxy), C_L*=max corr_L/W over a bank, the bound,
and whether consec maximizes meas(S_L).
"""
import sys, itertools, random
from math import gcd
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
random.seed(6187)
SEV=F(1,7); H=F(1,14)

def measSL(E, L):
    """fast: point p covers arcs j with L(p-1/7) < j <= Lp (a contiguous run mod L);
       all L arcs covered iff union of runs = Z/L."""
    E=sorted(set(E))
    bps=set([F(0),F(1)]); ends=set()
    for j in range(L):
        ends.add(F(j,L)%1); ends.add((F(j,L)+SEV)%1)
    for e in E:
        if e==0: continue
        for t in ends:
            for m in range(e): bps.add((t+m)/e)
    bps=sorted(z for z in bps if 0<=z<=1)
    total=F(0)
    for i in range(len(bps)-1):
        x0,x1=bps[i],bps[i+1]
        if x1<=x0: continue
        xm=(x0+x1)/2
        covered=[False]*L
        for e in E:
            p=(e*xm)%1
            # j in (L*(p-1/7), L*p]  -> integer j, mod L
            jhi=int(p*L)               # floor(L p) (p<1 so <L); arc j=jhi has j/L<=p
            jlo=int((p-SEV)*L)         # floor(L(p-1/7))
            # j ranges jlo+1 .. jhi (mod L)
            for jj in range(jlo+1, jhi+1):
                covered[jj % L]=True
            if all(covered): break
        if all(covered): total+=x1-x0
    return total

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

def dissoc(k): return [0]+[2**i for i in range(k-1)]   # {0,1,2,4,8,...} high-height proxy

for k in (8,9,10):
    capk=cap(k); Wc=W(list(range(k)))
    print(f"\n===== k={k}: cap_k={float(capk):.4f}, W(consec)={Wc:.3f} =====")
    # bank: consec + perforated (drop one) + random bounded-spread
    bank=[list(range(k))]
    for drop in range(1,k):
        bank.append([0]+[i for i in range(1,k+1) if i!=drop][:k-1])
    for _ in range(40):
        sp=random.randint(k-1,k+4); E=sorted(set([0,sp]+random.sample(range(1,sp),k-2)))
        if len(E)==k: bank.append(E)
    for L in (7,14,21):
        ML=float(measSL(dissoc(k),L))               # high-height main proxy
        maxratio=0.0; consmax=True; sc=float(measSL(list(range(k)),L))
        for E in bank:
            sL=float(measSL(E,L)); w=W(E)
            if sL>sc+1e-9: consmax=False
            corr=sL-ML
            if w>1e-9 and corr/w>maxratio: maxratio=corr/w
        bound=ML+maxratio*Wc
        print(f"  L={L:2d}: M_L={ML:.4f}  C_L*={maxratio:.5f}  meas(S_L,consec)={sc:.4f}  "
              f"bound=M_L+C_L*W(consec)={bound:.4f}  cap={float(capk):.4f}  "
              f"{'CLOSES' if bound<=float(capk) else 'no'}  consec-max-S_L:{consmax}")
print("\nDONE.")
