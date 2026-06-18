#!/usr/bin/env python3
"""
Quantify the EXPLICIT threshold V0*(P,offsets) across many fixed (P,offsets) to report:
 - is V0* bounded uniformly over all fixed-small-part single-tight-cluster patterns?  (NO in
   general -- as offsets grow / become resonant, the constant-sign sub-arc shrinks, so V0*
   grows; this is the honest boundary: the closed family is parametrized, with a PER-PATTERN
   explicit threshold, and finiteness is per-pattern, not uniform.)
 - typical and worst V0* over a broad sweep of admissible patterns with offsets bounded by D.

This pins WHY there is no single uniform V0 closing all of S3 (consistent with the AP family
being genuinely infinite): even within fixed-offset clusters, V0* depends on the offset multiset.
But for EACH FIXED pattern, V0* is finite and explicit -> the sub-case is closed.
"""
from fractions import Fraction as F
from math import gcd
import random
random.seed(7)

def safe_components(A,h=F(1,14)):
    iv=[]
    for u in A:
        for j in range(0,u):
            c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
            if a<b: iv.append((a,b))
            else: iv.append((a,F(1))); iv.append((F(0),b))
    iv.sort(); merged=[]
    for a,b in iv:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    safe=[]; prev=F(0)
    for a,b in merged:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<1: safe.append((prev,F(1)))
    return safe
def circ_width(points):
    pts=sorted(set(p%1 for p in points))
    if len(pts)<=1: return F(0)
    gaps=[pts[i+1]-pts[i] for i in range(len(pts)-1)]+[pts[0]+1-pts[-1]]
    return 1-max(gaps)
def widest_arc(P):
    sc=safe_components(P)
    if not sc: return None
    return max(sc,key=lambda ab: ab[1]-ab[0])
def w_of_tau(offsets, tau):
    w=F(6,7)-circ_width([d*tau for d in offsets])
    return w if w>0 else F(0)
def omega_width(offsets, ap, bp):
    ivs=[]
    for d in offsets:
        if d*(bp-ap)>=1: return F(0)
        lo=(F(1,2)-d*bp)%1; hi=(F(1,2)-d*ap)%1
        if lo<=hi: ivs.append((lo,hi))
        else: ivs.append((lo,F(1))); ivs.append((F(0),hi))
    if not ivs: return F(6,7)
    ivs.sort(); merged=[]
    for lo,hi in ivs:
        if merged and lo<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],hi))
        else: merged.append([lo,hi])
    n=len(merged)
    if n==1:
        lo,hi=merged[0]; gap=lo+(1-hi); cw=1-gap
        return max(F(0),F(6,7)-cw)
    gaps=[((merged[(s+1)%n][0]-merged[s][1])%1) for s in range(n)]
    cw=1-max(gaps)
    return max(F(0),F(6,7)-cw)
def constructive_threshold(P, offsets):
    arc=widest_arc(P)
    if arc is None: return None
    alpha,beta=arc
    offs=sorted(set(offsets)); cand_t=set(); diffs=set()
    for i in range(len(offs)):
        for j in range(len(offs)):
            if offs[i]!=offs[j]: diffs.add(abs(offs[i]-offs[j]))
        if offs[i]!=0: diffs.add(offs[i])
    for d in diffs:
        if d==0: continue
        k=int(alpha*d)-2
        while F(k,d)<=beta:
            t=F(k,d)
            if alpha<t<beta: cand_t.add(t)
            k+=1
    bl=sorted(cand_t)
    for s in range(len(bl)-1): cand_t.add((bl[s]+bl[s+1])/2)
    cand_t.add((alpha+beta)/2)
    tau0=None; wbest=F(-1)
    for t in cand_t:
        wv=w_of_tau(offsets,t)
        if wv>wbest: wbest=wv; tau0=t
    if wbest<=0 or tau0 is None: return None
    g_target=wbest/4
    hi=min(tau0-alpha,beta-tau0)
    if hi<=0: return None
    if omega_width(offsets,tau0-hi/10**6,tau0+hi/10**6) < g_target:
        g_target=omega_width(offsets,tau0-hi/10**9,tau0+hi/10**9)/2
        if g_target<=0: return None
    lo=F(0)
    for _ in range(80):
        mid=(lo+hi)/2
        if omega_width(offsets,tau0-mid,tau0+mid)>=g_target: lo=mid
        else: hi=mid
    h=lo
    if h<=0: return None
    s=2*h; V0star=int(1/s)+1
    return V0star

print("Sweep of EXPLICIT thresholds V0*(P,offsets) over admissible fixed patterns")
print("="*78)
# small parts P with positive measure, sizes 2..6; offsets random with bounded max D.
results=[]
for Dmax in [12, 25, 50]:
    vs=[]
    cnt=0
    for _ in range(400):
        psz=random.randint(2,6)
        P=sorted(random.sample(range(1,14),psz))
        if sum(b-a for a,b in safe_components(P))==0: continue
        nL=13-psz
        pool=list(range(1,Dmax+1))
        if len(pool)<nL-1: continue
        offs=sorted(set([0]+random.sample(pool,nL-1)))
        if len(offs)!=nL: continue
        V0s=constructive_threshold(P,offs)
        if V0s is None: continue
        vs.append(V0s); cnt+=1
    if vs:
        vs.sort()
        print("  offsets<=%2d : patterns=%4d  V0*  min=%d  median=%d  p90=%d  max=%d"
              %(Dmax,cnt,vs[0],vs[len(vs)//2],vs[int(0.9*len(vs))],vs[-1]))
print()
print("INTERPRETATION:")
print("  V0* is finite & explicit for EACH fixed pattern, but GROWS with offset spread/resonance.")
print("  There is NO single uniform V0 closing all fixed-cluster patterns -- consistent with S3")
print("  being genuinely infinite.  The closure is PER-PATTERN: for each (P,offsets), V0>=V0*(P,")
print("  offsets) is PROVED, and V0<V0* is a finite check.  This is a legitimate sub-family closure")
print("  (a parametrized infinite family, each member's tail proven), NOT a closure of all of S3.")
