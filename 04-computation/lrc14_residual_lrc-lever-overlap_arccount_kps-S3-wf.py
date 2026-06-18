#!/usr/bin/env python3
"""
lrc14_residual_lrc-lever-overlap_arccount_kps-S3-wf

THE ARC-COUNTING / TRANSVERSALITY core.

Established (decouple script): for every S3 set, the cluster L'=L\\{Vmax} alone has
W(G_{L'}) > 1/(7 Vmax), ratio 5.5..8.2.  But meas(G_P)+meas(G_{L'})<1 so a pure measure
pigeonhole fails.  The win must come from: G_{L'} has SEVERAL wide arcs spread across
[0,1), and at least one lies in a P-safe region.

Strategy to PROVE:  Count the L'-safe arcs of width > 1/(7 Vmax).  Call this N_wide.
The P-danger set D_P has measure m_DP and consists of teeth.  If N_wide arcs are spread
so that they cannot ALL be hit by P's teeth, one survives.

Specifically: a wide L'-arc of width w is KILLED by P only if NO sub-interval of width
1/(7Vmax) inside it is P-safe, i.e. P-danger covers all of it minus margins.  A single
P-tooth of a runner u<=13 has width 1/(7u) >= 1/91.  Compare to 1/(7Vmax) (tiny).
So one P-tooth can "block" a wide L'-arc only over its own width.

KEY measurable claims:
  (K1) Number N_wide of L'-arcs with width > 1/(7Vmax).
  (K2) Number of "maximal P-safe sub-windows" that each have width >= 1/(7Vmax) AND lie
       inside some L'-arc.  If >=1 always, C(S) holds.
  (K3) The DECISIVE count:  among the L'-safe arcs (each itself a union the cluster carves),
       how many are ENTIRELY inside G_P?  And what is the total L'-safe measure that lies
       in G_P (= meas(G_A))?  Compare to 1/(7Vmax).

Also test the cleaner sufficient condition:
  (SUF) meas(G_{L'}) > meas(D_P).   Then by measure, G_{L'} cannot be covered by D_P, so
        G_{L'} cap G_P has positive measure -> SOME common safe point -> but we need a wide
        ARC, not just a point.  Strengthen: if meas(G_{L'}) > meas(D_P) AND every L'-arc has
        width > 1/(7Vmax)... no.  Need arc-level.

  (SUF2) The cluster L' produces arcs that come in a near-periodic comb of period ~1/V0 with
        the wide arcs of width ~6/(7V0) appearing in a bounded pattern.  Count how many of
        the wide L'-arcs survive intersecting D_P.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random, sys

def flush(*a):
    print(*a); sys.stdout.flush()
C=F(1,14)

def nrm(x):
    r=x-int(x); r=r+1 if r<0 else r
    return r if r<=F(1,2) else 1-r

def safe_components(A,h=C):
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

def measure_safe(A): return sum(b-a for a,b in safe_components(A))
def widest_arc(arclist):
    if not arclist: return F(0)
    ws=[b-a for a,b in arclist]
    if arclist[0][0]==0 and arclist[-1][1]==1 and len(arclist)>1:
        ws.append(arclist[0][1]+(1-arclist[-1][0]))
    return max(ws)
def Wwidth(A): return widest_arc(safe_components(A))
def intersect_arclists(L1,L2):
    out=[]; i=j=0
    while i<len(L1) and j<len(L2):
        a,b=L1[i]; c,d=L2[j]
        lo=max(a,c); hi=min(b,d)
        if lo<hi: out.append((lo,hi))
        if b<d: i+=1
        else: j+=1
    return out
def is_covering(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def gcd_all(S): return reduce(gcd,S)

def make_S3_set(rng, max_attempts=20000):
    for _ in range(max_attempts):
        small=set(rng.sample(range(1,14), rng.randint(4,9)))
        V0=rng.randint(20,400); s=rng.randint(8,45)
        clsize=rng.randint(2, 13-len(small))
        cluster=set(); tries=0
        while len(cluster)<clsize and tries<200:
            cluster.add(V0+rng.randint(0,s)); tries+=1
        S=small|cluster
        while len(S)<13:
            S.add(rng.randint(V0,V0+s) if rng.random()<0.5 else rng.randint(1,13))
        if len(S)!=13: continue
        S=set(list(S))
        if len(S)!=13 or gcd_all(S)!=1 or not is_covering(S): continue
        Vmin=min(S); Vmax=max(S); k=sum(1 for v in S if v>13)
        if k>=2 and Vmax>=13*Vmin: return S
    return None

def run(n_sets=300, seed=23):
    rng=random.Random(seed)
    flush("="*70)
    flush("ARC-COUNT / TRANSVERSALITY: how many wide L'-arcs, how many survive D_P")
    flush("="*70)
    got=0; attempts=0
    SUF_holds=0      # meas(G_L') > meas(D_P)
    Nwide_min=10**9
    surviving_min=10**9
    surviving_data=[]
    suf_ratio_min=None
    while got<n_sets and attempts<n_sets*40:
        attempts+=1
        S=make_S3_set(rng)
        if S is None: continue
        got+=1
        Vmax=max(S)
        P=set(v for v in S if v<=13); L=set(v for v in S if v>13)
        Lp=L-{Vmax}
        thresh=F(1,7*Vmax)
        GLp = safe_components(Lp) if Lp else [(F(0),F(1))]
        GP  = safe_components(P)  if P  else [(F(0),F(1))]
        mLp=measure_safe(Lp) if Lp else F(1)
        mDP=1-(measure_safe(P) if P else F(1))
        # SUF
        ratio = mLp/mDP if mDP>0 else F(10**9)
        if suf_ratio_min is None or ratio<suf_ratio_min: suf_ratio_min=ratio
        if mLp>mDP: SUF_holds+=1
        # count wide L'-arcs
        wide=[(a,b) for a,b in GLp if b-a>thresh]
        # wrap arc
        if GLp and GLp[0][0]==0 and GLp[-1][1]==1 and len(GLp)>1:
            if GLp[0][1]+(1-GLp[-1][0])>thresh and (GLp[-1][0],GLp[0][1]) not in [(a,b) for a,b in [GLp[-1],GLp[0]]]:
                pass
        Nwide=len(wide)
        Nwide_min=min(Nwide_min,Nwide)
        # how many wide L'-arcs contain a P-safe sub-window of width>=thresh
        inter=intersect_arclists(GLp,GP)
        # surviving = number of intersection arcs of width>thresh
        survivors=[(a,b) for a,b in inter if b-a>thresh]
        # wrap survivor
        if inter and inter[0][0]==0 and inter[-1][1]==1 and len(inter)>1:
            if inter[0][1]+(1-inter[-1][0])>thresh: survivors.append(('wrap',))
        nsurv=len(survivors)
        surviving_min=min(surviving_min,nsurv)
        surviving_data.append((Nwide,nsurv,float(ratio)))
    flush(f"\n  S3 sets: {got}")
    flush(f"  SUF (meas(G_L') > meas(D_P)) holds:  {SUF_holds}/{got}")
    flush(f"    min ratio meas(G_L')/meas(D_P) = {float(suf_ratio_min):.4f}  exact={suf_ratio_min}")
    flush(f"  min # wide L'-arcs (width>1/(7Vmax)): {Nwide_min}")
    flush(f"  min # SURVIVING common arcs (width>1/(7Vmax)): {surviving_min}")
    avg_wide=sum(d[0] for d in surviving_data)/len(surviving_data)
    avg_surv=sum(d[1] for d in surviving_data)/len(surviving_data)
    flush(f"  avg # wide L'-arcs: {avg_wide:.1f}   avg # survivors: {avg_surv:.1f}")

if __name__=="__main__":
    run(300)
