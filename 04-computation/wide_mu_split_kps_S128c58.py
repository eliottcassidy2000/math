#!/usr/bin/env python3
"""wide_mu_split_kps_S128c58.py -- kind-pasteur S128 cont.58.
Split the COVERING WIDE CLUSTER census by mu = min(core).  Print data only."""
import sys
from fractions import Fraction as F
from math import ceil, gcd
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def la(v,q):
    r=v%q; return min(r,q-r)
def band_hits(P,K,qlo=15,qhi=600):
    out=[]
    for q in range(qlo,qhi+1):
        E=[]; ok=True
        for v in list(P)+list(K):
            e=la(v,q)
            if e==0: ok=False; break
            E.append(e)
        if not ok: continue
        emin,emax=min(E),max(E)
        if emax>13*emin: continue
        lo=F(q,14*emin); hi=F(13*q,14*emax); a=ceil(lo)
        if a>hi or a<=0: continue
        out.append((q,a,emin,emax))
    return out
def first_witness(P,K,qhi=600):
    V=list(P)+list(K)
    for q in range(2,qhi+1):
        for a in range(1,q//2+1):
            if gcd(a,q)!=1: continue
            t=F(a,q)
            if min(nd(v*t) for v in V)>=F(1,14): return (q,a,min(nd(v*t) for v in V))
    return None
from collections import defaultdict
stat=defaultdict(lambda:[0,0]); uncert=[]
for drop in range(1,13):
    P=[x for x in range(1,13) if x!=drop]
    mu,M=min(P),max(P); span=M-mu
    for a13 in range(1,40):
        k1=13*a13
        if k1<=13*M: continue
        for b14 in range(1,40):
            k2=14*b14
            if k2<=13*M or k2==k1: continue
            if abs(k2-k1)<=span: continue
            K=sorted([k1,k2])
            if len(set(P+K))!=13 or K[0]<=13*M: continue
            h=band_hits(P,K)
            stat[mu][0]+=1
            if h: stat[mu][1]+=1
            else: uncert.append((drop,P,K,abs(k2-k1)))
print("### covering WIDE clusters, band-certified split by mu = min(core) ###")
print("  mu   tested   band-certified   rate")
for mu in sorted(stat):
    t,c=stat[mu]; print("  %-3d  %-8d %-16d %.4f"%(mu,t,c,c/t))
print()
print("### the uncertified ones, in full ###")
print("  #uncertified: %d"%len(uncert))
seen=set()
for drop,P,K,D in uncert:
    key=(drop,tuple(K))
    if key in seen: continue
    seen.add(key)
    if len(seen)>12: break
    w=first_witness(P,K)
    print("  drop=%-3d core=%s"%(drop,P))
    print("      K=%s spread=%d ; 13|%d ; 14|%d ; first true witness: %s"%(
        K,D,K[0] if K[0]%13==0 else K[1], K[0] if K[0]%14==0 else K[1], w))
    if w:
        q,a,mn=w
        E=sorted(set(la(v,q) for v in list(P)+K))
        print("      at that q: residues=%s  emin=%d emax=%d  ratio=%.2f  (band needs emax<=13*emin AND an integer in [q/14emin,13q/14emax])"%(
            E,min(E),max(E),max(E)/min(E)))
        lo=F(q,14*min(E)); hi=F(13*q,14*max(E))
        print("      band interval at that q = [%s, %s] = [%.3f, %.3f] ; a used = %d ; a in band: %s"%(
            lo,hi,float(lo),float(hi),a,lo<=a<=hi))
print("DONE")
