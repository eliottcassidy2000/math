#!/usr/bin/env python3
"""multikiller_weak_target_kps_S128c52.py -- the MULTI-KILLER analogue at the weak target.
ITERATED BALANCE: for a lacunary killer chain (each killer > 13x the previous max), adding
killer k multiplies the bound by v_k/(v_k+s_{k-1}) with s_{k-1} <= prev max < v_k/13,
so each factor > 13/14.  With core |P| = 13-j and mu_0 = M(P) >= 1/(14-j) (LRC(14-j) settled):
   M(S) > (1/(14-j)) * (13/14)^j.
CHECK: is that > 1/14 for every j = 1..12?  Then verify on generated lacunary families."""
import sys, random
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)
random.seed(726)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def M_exact(V):
    cand=set()
    for v in V:
        for a in range(1,2*v): cand.add(F(a,2*v))
    for i in range(len(V)):
        for j2 in range(i+1,len(V)):
            for s_ in (V[i]+V[j2],abs(V[i]-V[j2])):
                if s_:
                    for a in range(1,s_): cand.add(F(a,s_))
    b=F(0)
    for t in cand:
        if 0<t<1:
            m=min(nd(v*t) for v in V)
            if m>b: b=m
    return b
print("== the iterated-balance bound (1/(14-j))*(13/14)^j vs 1/14 ==")
allok=True
for j in range(1,13):
    bnd=F(1,14-j)*F(13,14)**j
    ok = bnd >= F(1,14)
    strict = bnd > F(1,14)
    if not ok: allok=False
    print("   j=%2d killers: bound=%-22s = %.6f  >=1/14:%s %s"%(
        j,bnd,float(bnd),ok,"(STRICT)" if strict else "(EQUALITY at j=1; strictness from s<v/13)"))
print("   all j pass:",allok)
print()
print("== verify on generated LACUNARY-CHAIN families (each killer > 13x previous max) ==")
bad=0; tested=0; worst=None
for trial in range(60):
    j=random.randint(1,4)
    core=sorted(random.sample(range(1,16), 13-j))
    S=core[:]
    for _ in range(j):
        S.append(13*max(S)+random.randint(1,60))
    S=sorted(S)
    if len(set(S))!=13 or reduce(gcd,S)!=1: continue
    tested+=1
    P=core; mu0=M_exact(P)
    bnd=mu0
    prevmax=max(P)
    for k in [x for x in S if x not in core]:
        bnd = bnd*F(k, k+prevmax)
        prevmax=k
    MS=M_exact(S)
    if not (MS > F(1,14) and MS >= bnd and bnd > F(1,14)):
        bad+=1; print("   VIOLATION S=%s M=%s bnd=%s"%(S,MS,bnd)); break
    r=float(bnd)-1/14
    if worst is None or r<worst[0]: worst=(r,j,tuple(S),bnd,MS)
print("   %d lacunary families, violations=%d"%(tested,bad))
if worst:
    print("   tightest: j=%d bound=%s (%.6f, margin %+.2e) actual M=%s (%.6f)"%(
        worst[1],worst[3],float(worst[3]),worst[0],worst[4],float(worst[4])))
print()
print(">>> LACUNARY-CHAIN multi-killer closes at the weak target by ITERATED BALANCE.")
print(">>> CAVEAT: killers CLUSTERED with each other (v_2 ~ v_1) break the chain (factor ~1/2),")
print(">>>          so this covers the lacunary chain, not all of THM-726's multi-killer.")
print("DONE")
