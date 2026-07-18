#!/usr/bin/env python3
"""c_generalized_band_kps_S128c58.py -- kind-pasteur S128 cont.58.
THE c-GENERALIZED CLUSTER BAND.  THM-1032 took q = v1 + M.  General: q = v1 + c gives
killer residues  e_i = |delta_i - c|  (delta_i = v_i - v1), because v_i < q => q - v_i =
c - delta_i, and v_i > q => v_i - q = delta_i - c.  THM-1032 is c = M (residues M-delta_i,
needs delta_i <= M-mu).  Letting c FLOAT: need dist(c, delta_i) in [mu, M] for all i.
For r=2, Delta = {0,D}: c in [mu,M] cap ([D-M,D-mu] u [D+mu,D+M]), nonempty iff
   2*mu <= D <= 2*M   (the reflected branch)   OR   D <= M-mu  (THM-1032's branch).
PREDICTION: explicit reach extends from D <= M-mu to D <= 2M.  Test it.  PRINT DATA ONLY."""
import sys, random
from fractions import Fraction as F
from math import ceil
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def la(v,q):
    r=v%q; return min(r,q-r)
def c_cert(P,K,Ucap=None):
    """try every c; return first (c,q,a,minnorm) whose band works"""
    mu=min(P); M=max(P); v1=K[0]; deltas=[v-v1 for v in K]
    U = M if Ucap is None else Ucap
    for c in range(1, U*3+2):
        E=list(P)
        good=True
        for d in deltas:
            e=abs(d-c)
            if e<mu or e>U: good=False; break
            E.append(e)
        if not good: continue
        q=v1+c
        emin,emax=min(E),max(E)
        if 13*emin<=emax: continue
        lo=F(q,14*emin); hi=F(13*q,14*emax); a=ceil(lo)
        if a>hi or a<=0: continue
        # verify residues really are |d-c|
        if any(la(v,q)!=abs((v-v1)-c) for v in K): continue
        t=F(a,q); mn=min(nd(v*t) for v in list(P)+list(K))
        if mn>=F(1,14): return (c,q,a,mn)
    return None
print("### reach test: r=2 killers, core = 11-subset of {1..12}, spread D swept ###")
print("  mu M   D-range        tested  THM-1032(c=M)  c-general(U=M)  c-general(U=22)")
random.seed(58)
for drop in [1,6,12]:
    P=[x for x in range(1,13) if x!=drop]; mu,M=min(P),max(P)
    for lo,hi in [(1,M-mu),(M-mu+1,2*M),(2*M+1,3*M),(3*M+1,6*M)]:
        tot=0;c1=0;c2=0;c3=0
        for _ in range(260):
            D=random.randint(lo,hi)
            v1=random.randint(13*M+1,13*M+2000); K=[v1,v1+D]
            if len(set(P+K))!=13: continue
            tot+=1
            # THM-1032 branch: c = M exactly
            if D<=M-mu: c1+=1
            if c_cert(P,K): c2+=1
            if c_cert(P,K,Ucap=22): c3+=1
        if tot: print("  %-2d %-2d  [%3d,%4d]      %4d    %5d          %5d           %5d"%(
            mu,M,lo,hi,tot,c1,c2,c3))
print()
print("### exhaustive: does c-general cover ALL D <= 2M ? (core-drop x D x v1) ###")
tot=0; ok=0; fails=[]
for drop in range(1,13):
    P=[x for x in range(1,13) if x!=drop]; mu,M=min(P),max(P)
    for D in range(1,2*M+1):
        for v1 in range(13*M+1,13*M+60):
            K=[v1,v1+D]
            if len(set(P+K))!=13: continue
            tot+=1
            if c_cert(P,K): ok+=1
            else: fails.append((drop,mu,M,D,v1))
print("  tested %d ; c-general certified %d ; failures %d"%(tot,ok,len(fails)))
if fails:
    from collections import Counter
    print("  failure D-values:",sorted(Counter(f[3] for f in fails).items())[:12])
    print("  failure mu-values:",sorted(Counter(f[1] for f in fails).items()))
    for f in fails[:4]: print("    drop=%d mu=%d M=%d D=%d v1=%d"%f)
print()
print("### sharpness: first D where c-general (U=M) fails, per core ###")
for drop in range(1,13):
    P=[x for x in range(1,13) if x!=drop]; mu,M=min(P),max(P)
    firstfail=None
    for D in range(1,3*M+6):
        v1=13*M+7; K=[v1,v1+D]
        if len(set(P+K))!=13: continue
        if not c_cert(P,K): firstfail=D; break
    print("  drop=%-3d mu=%-2d M=%-2d  2M=%-3d  first failing D = %s"%(drop,mu,M,2*M,firstfail))
print("DONE")
