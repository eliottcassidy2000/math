#!/usr/bin/env python3
"""dmin_stratum_kps_S128c55.py -- kind-pasteur S128 cont.55.
MEASURE the d_min <= 5 stratum (killer blocks with two speeds within 5).
Question (data only): is this stratum DANGEROUS (M near 1/14) or merely uncertifiable
by (BG-K)?  Exact M on covering families with near-equal killers."""
import sys, random
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)
random.seed(55)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def M_exact(V):
    cand=set()
    for v in V:
        for a in range(1,2*v): cand.add(F(a,2*v))
    for i in range(len(V)):
        for j in range(i+1,len(V)):
            for s in (V[i]+V[j],abs(V[i]-V[j])):
                if s:
                    for a in range(1,s): cand.add(F(a,s))
    b=F(0); bt=None
    for t in cand:
        if 0<t<1:
            m=min(nd(v*t) for v in V)
            if m>b: b=m; bt=t
    return b,bt
def is_covering(V): return all(any(v%q==0 for v in V) for q in range(2,15))
def dmin_block(K):
    return min(K[j]-K[i] for i in range(len(K)) for j in range(i+1,len(K)))
print("exact M for families with NEAR-EQUAL killers (d_min <= 5), core max m, killers > 13m")
print("  %-30s %-6s %-4s %-12s %-9s %-6s"%("V","d_min","cov","M (exact)","M float","M*14"))
rows=[]
tries=0
while len(rows)<12 and tries<80000:
    tries+=1
    ncore=random.choice([9,10,11])
    core=sorted(random.sample(range(1,13), ncore))
    m=max(core); j=13-ncore
    if j<2: continue
    base=random.randint(13*m+1, 13*m+120)
    d=random.randint(1,5)
    K=[base, base+d]
    while len(K)<j: K.append(base+random.randint(6,80))
    K=sorted(set(K))
    if len(K)!=j: continue
    V=sorted(set(core+K))
    if len(V)!=13: continue
    if reduce(gcd,V)!=1: continue
    dm=dmin_block(K)
    if dm>5: continue
    Mx,t=M_exact(V)
    cov=is_covering(V)
    rows.append((float(Mx),dm,cov,tuple(V),Mx))
    print("  %-30s %-6d %-4s %-12s %-9.6f %-6.3f"%(str(list(V))[:30],dm,"Y" if cov else "n",Mx,float(Mx),float(Mx)*14))
print()
if rows:
    rows.sort()
    print("  min M over stratum: %.6f  (= %.3f x threshold 1/14)"%(rows[0][0], rows[0][0]*14))
    print("  max M over stratum: %.6f  (= %.3f x)"%(rows[-1][0], rows[-1][0]*14))
    covrows=[r for r in rows if r[2]]
    print("  covering members: %d ; their min M: %s"%(len(covrows), ("%.6f (%.3f x)"%(covrows[0][0],covrows[0][0]*14)) if covrows else "none"))
    print("  any at or below 1/14: %d"%len([r for r in rows if r[0]<=1/14+1e-12]))
print("DONE")
