#!/usr/bin/env python3
"""exact_support_decomposition_kps_S128c39.py -- kind-pasteur S128 cont.39.
THE EXACT SIGNED SUPPORT DECOMPOSITION: the M_s relation masses of THM-935 computed as
EXACT RATIONALS by inclusion-exclusion of subset overlap measures:
  Mhat(A) := mu(cap_{v in A} B_v)  (exact sweep per subset)
  M(A)    := Mhat(A) - sum_{S subsetneq A, |S|>=2} M(S) * (2lam)^{|A|-|S|} - (2lam)^{|A|}
  M_s     := sum_{|A|=s} M(A)
Referee: B5(direct sweep) == B5eq + sum_s (-1)^s E_s^{(5)} M_s  EXACTLY (rationals).
Outputs: per-packet decomposition, signed-vs-absolute slack, top blocker subsets."""
import sys
from fractions import Fraction as F
from itertools import combinations
from math import comb, gcd
sys.stdout.reconfigure(line_buffering=True)
LAM=F(1,14); TWO=2*LAM

def sweep_mu_all(speeds):
    """exact mu(cap of ALL bad sets of these speeds) via integer sweep on 14*lcm grid."""
    L=1
    for v in speeds: L=L*v//gcd(L,v)
    D=14*L
    ev=[]
    for v in speeds:
        w=L//v
        for j in range(v):
            lo=(14*j-1)*w; hi=(14*j+1)*w
            if lo<0: ev.append((0,1)); ev.append((hi,-1)); ev.append((D+lo,1)); ev.append((D,-1))
            else: ev.append((lo,1)); ev.append((hi,-1))
    ev.sort()
    k=len(speeds); depth=0; last=0; tot=0
    for x,d in ev:
        if depth==k: tot+=x-last
        depth+=d
        if depth==k: last=x
    return F(tot,D)

def depth_B5(speeds):
    """direct exact B5 via the full depth-spectrum sweep."""
    L=1
    for v in speeds: L=L*v//gcd(L,v)
    # full 13-speed lcm too big -> use Fraction sweep instead
    ev=[]
    for v in speeds:
        for j in range(v):
            lo=F(14*j-1,14*v); hi=F(14*j+1,14*v)
            if lo<0: ev.append((F(0),1)); ev.append((hi,-1)); ev.append((lo+1,1)); ev.append((F(1),-1))
            else: ev.append((lo,1)); ev.append((hi,-1))
    ev.sort()
    n=len(speeds); mu=[F(0)]*(n+1); depth=0; last=F(0)
    for x,d in ev:
        if x>last: mu[depth]+=x-last; last=x
        depth+=d
    if F(1)>last: mu[depth]+=F(1)-last
    S=[sum(comb(d,k)*mu[d] for d in range(n+1)) for k in range(6)]
    return sum((-1)**k*S[k] for k in range(6))

def decompose(name,pk):
    n=len(pk)
    M={}   # frozenset -> exact mass
    for s in range(2,6):
        for A in combinations(pk,s):
            mhat=sweep_mu_all(list(A))
            corr=TWO**s
            for r in range(2,s):
                for Ssub in combinations(A,r):
                    corr+=M[frozenset(Ssub)]*TWO**(s-r)
            M[frozenset(A)]=mhat-corr
        print("  %s: level %d done"%(name,s),flush=True)
    Ms={s: sum(M[frozenset(A)] for A in combinations(pk,s)) for s in range(2,6)}
    E={2:F(24,343),3:F(24,49),4:F(-2,7),5:F(1)}
    B5eq=sum((-1)**k*comb(13,k)*LAM.__mul__(2)**k for k in range(6))
    pred=B5eq+sum((-1)**s*E[s]*Ms[s] for s in range(2,6))
    direct=depth_B5(pk)
    ok=(pred==direct)
    print("  %s: EXACT IDENTITY %s"%(name,"HOLDS" if ok else "FAILS"))
    print("    M_2..M_5 = %s"%["%+.6f"%float(Ms[s]) for s in range(2,6)])
    absdebt=sum(abs(E[s])*abs(Ms[s]) for s in range(2,6))
    signed=sum((-1)**s*E[s]*Ms[s] for s in range(2,6))
    print("    B5 = %.6f = B5eq %.6f + signed %.6f ; ABSOLUTE debt bound %.6f (slack factor %.2f)"%(
        float(direct),float(B5eq),float(signed),float(absdebt),float(absdebt/abs(signed)) if signed!=0 else float('inf')))
    tops=sorted(M.items(),key=lambda kv:-abs(kv[1]))[:5]
    print("    top blocker subsets: %s"%[ (tuple(sorted(k)),"%+.5f"%float(v)) for k,v in tops])
    return Ms,direct,pred

cert=[307,425,541,671,800,944,1087,1413,1943,2147,2570,3056,3310]
opus30Z=[420,450,510,570,690,870,1230,1770,2370,3210,4170,5190,7230]
import random
random.seed(7)
rnd=sorted(random.sample(range(300,4001),13))
for name,pk in [("certified",cert),("opus30Z",opus30Z),("random",rnd)]:
    decompose(name,pk)
print("DONE")
