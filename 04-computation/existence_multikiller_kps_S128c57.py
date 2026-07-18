#!/usr/bin/env python3
"""existence_multikiller_kps_S128c57.py -- kind-pasteur S128 cont.57.
HOW FAR DOES q = v1 + M REACH?  The d<=5 two-killer hypothesis is NOT needed.

CLAIM (spread form).  V = P u {v_1<...<v_r}, |P| = 13-r core speeds, mu=min P, M=max P,
v_1 > 13M, and SPREAD D := v_r - v_1 <= M - mu.  Put q = v_1 + M.  Then each killer
v_i = v_1 + delta_i with delta_i <= D < M, so v_i < q and q - v_i = M - delta_i in [mu, M].
Hence E = P u {M - delta_i} lies in [mu,M]: e_min = mu, e_max = M, ratio inherited.
Band nonempty+integral since q >= 14M+1 and M <= 12 mu.  ==> certificate for ANY r, ANY
spread up to M - mu.  For a core inside {1..12}: M - mu is 10 or 11, so this covers every
clustered killer cluster of width <= 10 -- strictly more than the d <= 5 stratum."""
import sys, random
from fractions import Fraction as F
from math import ceil
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def la(v,q):
    r=v%q; return min(r,q-r)
def cert(P,K):
    mu=min(P); M=max(P); q=K[0]+M; a=ceil(F(q,14*mu))
    E=sorted(set(list(P)+[la(v,q) for v in K])); emin,emax=min(E),max(E)
    if emin==0: return None,'zero-residue'
    if a>F(13*q,14*emax): return None,'no-integer'
    V=sorted(list(P)+K); t=F(a,q)
    return (q,a,min(nd(v*t) for v in V),emin,emax),None
print("### r killers, spread D <= M-mu : predicted ALWAYS certified ###")
print("  r  spread-cap  tested  certified  min-margin(min||vt||)")
random.seed(571)
for r in [2,3,4,5]:
    for capfrac in ['<=M-mu','=M-mu+1(over)','=M-mu+3(over)']:
        tot=0; ok=0; worst=None
        for _ in range(1500):
            ncore=13-r
            P=sorted(random.sample(range(1,13),ncore))
            mu,M=min(P),max(P); span=M-mu
            if span<r: continue
            cap = span if capfrac=='<=M-mu' else span+(1 if 'M-mu+1' in capfrac else 3)
            v1=random.randint(13*M+1,13*M+3000)
            offs=sorted(random.sample(range(1,cap+1),r-1)) if cap>=r-1 else None
            if offs is None: continue
            if capfrac!='<=M-mu' and offs[-1]<=span: offs[-1]=cap   # force over-spread
            K=[v1]+[v1+o for o in offs]
            if len(set(P+K))!=13: continue
            tot+=1
            res,err=cert(P,K)
            if res and res[2]>=F(1,14):
                ok+=1
                if worst is None or res[2]<worst: worst=res[2]
        if tot: print("  %d  %-14s %5d   %5d      %s"%(r,capfrac,tot,ok,worst))
print()
print("### the d<=5 two-killer stratum is a SUBSET: exhaustive re-check with spread up to 11 ###")
tot=0; ok=0
for k in range(1,13):
    P=[x for x in range(1,13) if x!=k]; mu,M=min(P),max(P)
    for v1 in range(13*M+1,13*M+800):
        for d in range(1,M-mu+1):
            K=[v1,v1+d]
            if len(set(P+K))!=13: continue
            tot+=1
            res,err=cert(P,K)
            if res and res[2]>=F(1,14): ok+=1
print("  spread d = 1..M-mu (not just 1..5): certified %d / %d"%(ok,tot))
print()
print("### sharpness: what happens exactly at spread = M-mu+1 ? ###")
P=list(range(2,13)); mu,M=2,12
for d in range(9,15):
    v1=13*M+7; K=[v1,v1+d]
    res,err=cert(P,K)
    tag = ("min||vt||=%s ok=%s"%(res[2],res[2]>=F(1,14))) if res else err
    print("  spread %2d (M-mu=%d): %s"%(d,M-mu,tag))
print("DONE")
