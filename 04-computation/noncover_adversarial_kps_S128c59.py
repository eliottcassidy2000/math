#!/usr/bin/env python3
"""noncover_adversarial_kps_S128c59.py -- kind-pasteur S128 cont.59.
BEFORE trying to PROVE the non-covering lemma, ATTACK it.
A modulus q is TRIVIALLY killing if q | k1 or q | k2 (that killer is bad for every a).
So if a killer is divisible by every q in the search range, EVERY modulus is killing and
the small-modulus criterion has no witness at all.  Take k1 = lcm(15..Q).  Is the resulting
family a legal trapped-core-shape family?  PRINT DATA ONLY."""
import sys
from math import gcd
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def la(r,q):
    r%=q; return min(r,q-r)
def lcm(a,b): return a*b//gcd(a,b)
def killing(P,k1,k2,q):
    thr=-(-q//14)
    for a in range(1,q):
        if all(la(p*a,q)>=thr for p in P) and la(k1*a,q)>=thr and la(k2*a,q)>=thr:
            return False
    return True
def legal(P,k1,k2):
    V=sorted(P+[k1,k2])
    if len(set(V))!=13: return False,"not 13 distinct"
    if not all(any(v%q==0 for v in V) for q in range(2,15)): return False,"not covering"
    # compressed: every speed within 13x of SOME other
    for i,v in enumerate(V):
        if not any(j!=i and v<=13*V[j] for j in range(len(V))): return False,"not compressed"
    # gap family: some ratio > 13
    if not any(a>13*b for a in V for b in V): return False,"no gap"
    return True,"legal trapped-core shape"
print("### the adversarial construction: k1 = lcm(15..Q), k2 = 2*k1 ###")
P=list(range(2,13))   # core {2..12}, mu=2 M=12
for Q in [20,25,30,40,60]:
    L=1
    for q in range(15,Q+1): L=lcm(L,q)
    k1,k2=L,2*L
    ok,why=legal(P,k1,k2)
    nk=sum(1 for q in range(15,Q+1) if killing(P,k1,k2,q))
    tot=Q-15+1
    print("  Q=%-3d  lcm(15..Q) has %d digits ; legal: %-5s (%s)"%(Q,len(str(L)),ok,why))
    print("        killing at %d of %d moduli in [15,%d]  -> criterion has a witness: %s"%(
        nk,tot,Q,nk<tot))
print()
print("### is the criterion's FAILURE fatal?  scale separation of these families ###")
for Q in [20,30]:
    L=1
    for q in range(15,Q+1): L=lcm(L,q)
    print("  Q=%-3d killers L, 2L with L ~ 10^%d ; core max = 12"%(Q,len(str(L))-1))
    print("        killer/core scale ratio = L/12 ~ 10^%d ; killer-killer ratio = 2 (comparable)"%(len(str(L//12))-1))
print()
print("### how SMALL can a killer pair be and still kill every modulus in [15,Q]? ###")
print("  Q    min max(k1,k2) killing all of [15,Q]      example")
for Q in [18,20,22,24,26,28]:
    best=None
    for k1 in range(157,4001):
        if best and k1>best[0]: break
        for k2 in range(157,4001):
            if k1==k2: continue
            if best and max(k1,k2)>=best[0]: continue
            ok,_=legal(P,k1,k2)
            if not ok: continue
            if all(killing(P,k1,k2,q) for q in range(15,Q+1)):
                if best is None or max(k1,k2)<best[0]: best=(max(k1,k2),k1,k2)
    print("  %-4d %-38s %s"%(Q,best[0] if best else "> 4000 (none found)",
                             ("(%d,%d)"%(best[1],best[2])) if best else "-"))
print("DONE")
