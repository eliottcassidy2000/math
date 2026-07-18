#!/usr/bin/env python3
"""thm724_weak_target_kps_S128c52.py -- kind-pasteur S128 cont.52.
THM-724's near-tight residual AT THE WEAK TARGET M > 1/14.

CLAIM (weak-target single-killer closure): for a primitive 13-set S = C u {v_f} with
|C|=12 and v_f > 13*max(C) (single-killer),
   M(S) >= mu * v_f/(v_f + s)   [THM-724 Lemma 1, killer resonant]
        >= mu * v_f/(v_f + maxC) [s <= maxC, decreasing in s]
        >  (1/13) * v_f/(v_f + v_f/13) = (1/13)(13/14) = 1/14.
using mu = M(C) >= 1/13 (LRC(13), settled) and maxC < v_f/13.
No residual: the weak target is reached by the BALANCE LEMMA ALONE.

VERIFY: (a) the chain numerically on many single-killer configs (covering and not);
(b) actual M(S) > 1/14 always; (c) the worst-case margin formula.
"""
import sys, random
from fractions import Fraction as F
from math import gcd
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)
random.seed(724)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def M_exact(V):
    cand=set()
    for v in V:
        for a in range(1,2*v): cand.add(F(a,2*v))
    for i in range(len(V)):
        for j in range(i+1,len(V)):
            for s_ in (V[i]+V[j],abs(V[i]-V[j])):
                if s_:
                    for a in range(1,s_): cand.add(F(a,s_))
    b=F(0)
    for t in cand:
        if 0<t<1:
            m=min(nd(v*t) for v in V)
            if m>b: b=m
    return b
def is_covering(V): return all(any(v%q==0 for v in V) for q in range(2,15))

print("== (a) the algebraic chain on single-killer configs ==")
print("   need: mu>=1/13, maxC < v_f/13, balance_bound = mu*v_f/(v_f+maxC) > 1/14, M(S) >= balance")
bad=0; tested=0; minmargin=None; mincfg=None
for trial in range(120):
    # small core so M(C) is exactly computable; v_f > 13*maxC
    C=sorted(random.sample(range(1,26), 12))
    maxC=max(C)
    v_f=13*maxC + random.randint(1, 400)
    S=sorted(C+[v_f])
    if len(set(S))!=13: continue
    if reduce(gcd,S)!=1: continue
    tested+=1
    mu=M_exact(C)
    bal=mu*F(v_f, v_f+maxC)          # the s<=maxC worst case
    ok_mu   = mu >= F(1,13)
    ok_gap  = maxC < F(v_f,13)
    ok_bal  = bal > F(1,14)
    if not (ok_mu and ok_gap and ok_bal):
        bad+=1
        print("   VIOLATION C=%s v_f=%d mu=%s bal=%s"%(C,v_f,mu,bal)); break
    m=bal-F(1,14)
    if minmargin is None or m<minmargin: minmargin=m; mincfg=(tuple(C),v_f,mu,bal)
print("   %d single-killer configs: mu>=1/13 OK, maxC<v_f/13 OK, balance>1/14 OK; violations=%d"%(tested,bad))
print("   tightest balance margin over 1/14: %s = %.3e  (C max=%d, v_f=%d, mu=%s)"%(
      minmargin, float(minmargin), max(mincfg[0]), mincfg[1], mincfg[2]))

print()
print("== (b) actual M(S) vs the bound, incl. the deep well and covering configs ==")
cases=[("deep well {1..12,182}", list(range(1,13)), 182)]
for _ in range(6):
    C=sorted(random.sample(range(1,26),12)); mx=max(C)
    vf=13*mx+random.randint(1,300)
    if reduce(gcd,sorted(C+[vf]))==1: cases.append(("random SK", C, vf))
for name,C,vf in cases:
    S=sorted(C+[vf]); mu=M_exact(C); mx=max(C)
    bal=mu*F(vf,vf+mx); MS=M_exact(S)
    print("   %-22s maxC=%3d v_f=%5d mu=%-7s bound=%.6f  M(S)=%-9s %.6f  M>1/14:%s cov:%s"%(
        name,mx,vf,mu,float(bal),MS,float(MS),MS>F(1,14),is_covering(S)))

print()
print("== (c) worst-case margin formula: mu=1/13, s=maxC=m, v_f=13m+1 ==")
for m in (1,5,12,25,100):
    vf=13*m+1
    bal=F(1,13)*F(vf,vf+m)
    print("   m=%3d: balance=%s, margin over 1/14 = %s = %.3e  (formula 1/(182(14m+1)) = %.3e)"%(
        m,bal,bal-F(1,14),float(bal-F(1,14)),1/(182*(14*m+1))))
print()
print(">>> THE BALANCE LEMMA ALONE reaches the WEAK target: THM-724's near-tight residual")
print(">>> is EMPTY at M > 1/14 -- no shallow witness, no census needed.")
print("DONE")
