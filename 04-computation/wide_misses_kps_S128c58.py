#!/usr/bin/env python3
"""wide_misses_kps_S128c58.py -- kind-pasteur S128 cont.58.
(A) Reproduce and CHARACTERIZE the 3 wide-cluster misses of cont.57.
(B) Note: the cont.57 sweep did NOT impose COVERING.  Covering is the actual hypothesis
    of 'covering => M > 1/14', and with a core inside {1..12} it forces the KILLERS to
    carry the 13- and 14-divisibility.  Build the covering wide clusters and test them.
PRINT DATA ONLY -- no interpretation baked in."""
import sys, random
from fractions import Fraction as F
from math import ceil, gcd
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def la(v,q):
    r=v%q; return min(r,q-r)
def band_scan(P,K,qlo=15,qhi=400):
    hits=[]
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
        t=F(a,q); mn=min(nd(v*t) for v in list(P)+list(K))
        if mn>=F(1,14): hits.append((q,a,mn))
    return hits
def any_witness(P,K,qhi=2000):
    """ANY lonely witness a/q with q<=qhi (not just band-shaped)"""
    V=list(P)+list(K); best=(F(0),None)
    for q in range(2,qhi+1):
        for a in range(1,q//2+1):
            if gcd(a,q)!=1: continue
            t=F(a,q); mn=min(nd(v*t) for v in V)
            if mn>best[0]: best=(mn,(q,a))
            if mn>=F(1,14): return mn,(q,a)
    return best
print("### (A) reproduce the cont.57 wide-cluster misses ###")
random.seed(575)
misses=[]
for r in [2,3,4]:
    for lo_mult,hi_mult,lab in [(1,3,'a'),(3,10,'b'),(10,60,'c')]:
        for _ in range(500):
            P=sorted(random.sample(range(1,13),13-r))
            mu,M=min(P),max(P); span=M-mu
            if span<2: continue
            D=random.randint(span*lo_mult+1, span*hi_mult)
            v1=random.randint(13*M+1,13*M+3000)
            offs=sorted(random.sample(range(1,D+1),r-2)+[D]) if r>=2 else []
            K=[v1]+[v1+o for o in offs]
            if len(set(P+K))!=13: continue
            if not band_scan(P,K): misses.append((r,lab,P,K))
print("  misses found: %d"%len(misses))
for r,lab,P,K in misses:
    mu,M=min(P),max(P)
    print("   r=%d band=%s core=%s (mu=%d M=%d)  K=%s  spread=%d"%(r,lab,P,mu,M,K,K[-1]-K[0]))
    print("     core contains 1: %s ; 13|any: %s ; 14|any: %s"%(
        1 in P, any(v%13==0 for v in list(P)+K), any(v%14==0 for v in list(P)+K)))
    mn,qa=any_witness(P,K)
    print("     best witness with q<=2000: min||vt||=%s at %s  (>=1/14: %s)"%(mn,qa,mn>=F(1,14)))
print()
print("### (B) COVERING wide clusters: core in {1..12}, killers carry 13| and 14| ###")
print("  core-drop  mu M  killers        spread  D>M-mu  #band-q  first(q,a)   min||vt||")
rows=0; certified=0; nocert=[]
for drop in range(1,13):
    P=[x for x in range(1,13) if x!=drop]
    mu,M=min(P),max(P); span=M-mu
    for a13 in range(1,40):
        k1=13*a13
        if k1<=13*M: continue
        for b14 in range(1,40):
            k2=14*b14
            if k2<=13*M or k2==k1: continue
            D=abs(k2-k1)
            if D<=span: continue          # WIDE only
            K=sorted([k1,k2])
            if len(set(P+K))!=13: continue
            if K[0]<=13*M: continue
            rows+=1
            hits=band_scan(P,K)
            if hits:
                certified+=1
                if rows<=14:
                    print("  drop%-3d   %-2d %-2d %-14s %-7d %-7s %-8d %-12s %s"%(
                        drop,mu,M,str(K),D,D>span,len(hits),str(hits[0][:2]),hits[0][2]))
            else:
                nocert.append((P,K,D))
print("  ...")
print("  TOTAL covering wide clusters: %d ; band-certified: %d ; NOT: %d"%(rows,certified,len(nocert)))
if nocert:
    print("  uncertified examples (first 8):")
    for P,K,D in nocert[:8]:
        mn,qa=any_witness(P,K)
        print("    core=%s K=%s spread=%d  best q<=2000: %s at %s (>=1/14: %s)"%(
            [min(P),max(P)],K,D,mn,qa,mn>=F(1,14)))
print("DONE")
