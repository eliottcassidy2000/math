#!/usr/bin/env python3
"""wide_cluster_residual_kps_S128c57.py -- kind-pasteur S128 cont.57.
SIZING THE ONE REMAINING GAP.  THM-1032 closes clustered killers with spread D <= M-mu.
What does the residual (D > M-mu, non-lacunary) actually look like?  Print data only."""
import sys, random
from fractions import Fraction as F
from math import ceil
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def la(v,q):
    r=v%q; return min(r,q-r)
def band_scan(P,K,qhi=400):
    """does SOME modulus certify?  return (q,a,minnorm,count_valid)"""
    first=None; cnt=0
    for q in range(15,qhi+1):
        E=[]; ok=True
        for v in list(P)+K:
            e=la(v,q)
            if e==0: ok=False; break
            E.append(e)
        if not ok: continue
        emin,emax=min(E),max(E)
        if emax>13*emin: continue
        lo=F(q,14*emin); hi=F(13*q,14*emax)
        a=ceil(lo)
        if a>hi or a<=0: continue
        t=F(a,q); mn=min(nd(v*t) for v in list(P)+K)
        if mn>=F(1,14):
            cnt+=1
            if first is None: first=(q,a,mn)
    return first,cnt
random.seed(575)
print("### wide clusters (spread D > M-mu): does the BAND still find a modulus? ###")
print("  r  D-range        tested  band-certified  median#valid-q  first-q median")
for r in [2,3,4]:
    for lo_mult,hi_mult,lab in [(1,3,'(M-mu, 3(M-mu)]'),(3,10,'(3(M-mu),10(M-mu)]'),(10,60,'(10x,60x]')]:
        tot=0; ok=0; cnts=[]; qs=[]
        for _ in range(500):
            P=sorted(random.sample(range(1,13),13-r))
            mu,M=min(P),max(P); span=M-mu
            if span<2: continue
            D=random.randint(span*lo_mult+1, span*hi_mult)
            v1=random.randint(13*M+1,13*M+3000)
            offs=sorted(random.sample(range(1,D+1),r-2)+[D]) if r>=2 else []
            K=[v1]+[v1+o for o in offs]
            if len(set(P+K))!=13: continue
            tot+=1
            first,cnt=band_scan(P,K)
            if first: ok+=1; cnts.append(cnt); qs.append(first[0])
        if tot:
            cnts.sort(); qs.sort()
            print("  %d  %-18s %5d   %5d           %-6s        %s"%(
                r,lab,tot,ok,cnts[len(cnts)//2] if cnts else '-',qs[len(qs)//2] if qs else '-'))
print()
print("### control: spread <= M-mu (THM-1032 territory), same scan ###")
tot=0; ok=0
for _ in range(500):
    P=sorted(random.sample(range(1,13),11)); mu,M=min(P),max(P); span=M-mu
    if span<1: continue
    D=random.randint(1,span); v1=random.randint(13*M+1,13*M+3000)
    K=[v1,v1+D]
    if len(set(P+K))!=13: continue
    tot+=1
    first,cnt=band_scan(P,K)
    if first: ok+=1
print("  band-certified %d / %d  (THM-1032 also gives these an EXPLICIT q with no scan)"%(ok,tot))
print("DONE")
