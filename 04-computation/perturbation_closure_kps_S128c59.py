#!/usr/bin/env python3
"""perturbation_closure_kps_S128c59.py -- kind-pasteur S128 cont.59.
THE DICHOTOMY, MADE COMPLETE.
LARGE-KILLER HORN (perturbation).  For any core P inside {1..12}, t0 = 1/13 has
||p/13|| >= 1/13 for every p <= 12 (13 is prime > 12).  For |t - 1/13| <= rho,
   ||p t|| >= 1/13 - 12*rho ,  which is >= 1/14 exactly when rho <= 1/2184.
So the core-safe set CONTAINS the interval I = [1/13 - 1/2184, 1/13 + 1/2184], length
   ell = 1/1092.
A killer k is unsafe on a union of intervals of length 1/(7k) spaced 1/k apart, so
   |Bad_k n I| <= ell/7 + 2/(7k).
For r killers the union misses I whenever  r*ell/7 + 2r/(7 k_min) < ell , i.e.
   k_min > 2r / (ell * (7 - r)).      r = 2  =>  k_min > 4*1092/5 = 873.6.
So TWO killers, both >= 874, are certified by a t inside I -- no modulus search at all.
SMALL-KILLER HORN: killers in (13M, 874) is a FINITE set -- check it exhaustively with the
small-modulus criterion via bitmasks.  PRINT DATA ONLY."""
import sys
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def la(r,q):
    r%=q; return min(r,q-r)
print("### (1) the core-safe interval claim ###")
rho=F(1,2184)
lo,hi=F(1,13)-rho,F(1,13)+rho
worst=None
for drop in range(1,13):
    P=[x for x in range(1,13) if x!=drop]
    # exact: min over p of ||p t|| at the two endpoints and at 1/13
    for t in [lo,hi,F(1,13)]:
        m=min(nd(p*t) for p in P)
        if worst is None or m<worst[0]: worst=(m,drop,t)
print("  I = [1/13 - 1/2184, 1/13 + 1/2184] , length ell = 1/1092")
print("  min over 12 cores and over the endpoints/centre of min_p ||p t|| = %s (%.6f)"%(worst[0],float(worst[0])))
print("  threshold 1/14 = %.6f  -> interval valid: %s"%(1/14,worst[0]>=F(1,14)))
print("  (attained at core-drop %d, t = %s)"%(worst[1],worst[2]))
print()
print("### (2) large-killer horn: does a good t in I really exist for k >= 874? ###")
def find_t_in_I(P,K,steps=200000):
    """scan I on a fine grid for a t with every speed >= 1/14"""
    import fractions
    N=steps
    for j in range(N+1):
        t=lo+(hi-lo)*F(j,N)
        if min(nd(v*t) for v in list(P)+list(K))>=F(1,14): return t
    return None
import random
random.seed(59)
tested=0; found=0; fails=[]
for _ in range(220):
    drop=random.randint(1,12); P=[x for x in range(1,13) if x!=drop]
    k1=random.randint(874,60000); k2=random.randint(874,60000)
    if k1==k2 or len(set(P+[k1,k2]))!=13: continue
    tested+=1
    t=find_t_in_I(P,[k1,k2],4000)
    if t: found+=1
    else: fails.append((drop,k1,k2))
print("  random pairs with both killers >= 874 : %d tested, %d certified by a t inside I"%(tested,found))
if fails: print("  grid missed (grid is coarse, not a refutation):",fails[:4])
print()
print("  adversarial: the lcm-built killers that defeat EVERY modulus")
from math import gcd
def lcmf(a,b): return a*b//gcd(a,b)
L=1
for q in range(15,31): L=lcmf(L,q)
for drop in [1,6]:
    P=[x for x in range(1,13) if x!=drop]
    t=find_t_in_I(P,[L,2*L],20000)
    print("   drop=%d  K=(lcm(15..30), 2*lcm)  ~10^%d : t in I found: %s"%(
        drop,len(str(L))-1,"YES" if t else "no (grid too coarse)"))
    if t: print("        t=%s  min||vt||=%s"%(t,min(nd(v*t) for v in P+[L,2*L])))
print()
print("### (3) small-killer horn: EXHAUSTIVE finite check, killers in (13M, 874) ###")
QS=[(q,a) for q in range(15,41) for a in range(1,q)]
idx={qa:i for i,qa in enumerate(QS)}
print("  (q,a) pairs in play: %d"%len(QS))
def kmask(k):
    m=0
    for i,(q,a) in enumerate(QS):
        if la(k*a,q)>=-(-q//14): m|=(1<<i)
    return m
KMAX=874
masks={k:kmask(k) for k in range(1,KMAX)}
tot=0; cert=0; bad=[]
for drop in range(1,13):
    P=[x for x in range(1,13) if x!=drop]; M=max(P)
    cm=0
    for i,(q,a) in enumerate(QS):
        if all(la(p*a,q)>=-(-q//14) for p in P): cm|=(1<<i)
    for k1 in range(13*M+1,KMAX):
        m1=cm & masks[k1]
        if m1==0:
            # k1 alone defeats every (q,a); record and continue
            pass
        for k2 in range(k1+1,KMAX):
            V=P+[k1,k2]
            if len(set(V))!=13: continue
            if not all(any(v%q==0 for v in V) for q in range(2,15)): continue
            tot+=1
            if m1 & masks[k2]: cert+=1
            else: bad.append((drop,k1,k2))
print("  covering families with both killers < %d : %d"%(KMAX,tot))
print("  certified by the small-modulus criterion (q<=40) : %d"%cert)
print("  UNCERTIFIED : %d"%len(bad))
if bad:
    print("  examples:",bad[:8])
print("DONE")
