#!/usr/bin/env python3
"""
lrc14_residual_lrc-lever-overlap_struct_kps-S3-wf

STRUCTURE PROBE for the overlap angle. We locate WHERE the witness arc of
G_A = G_P cap G_{L'} sits and what controls its width.

Two competing geometric facts to reconcile (from the head-start):
  (a) The cluster L' = {V0+d : d in offsets} has danger teeth of width 1/(7u) each,
      spaced 1/u apart.  Near a given tau, ALL cluster members u land near the same
      fractional value (they "collapse"), so a CLUSTER-SAFE window is where they ALL
      land in a common gap.  Those windows have tau-width ~ (1 - small)/(7 V0) but are
      spaced ~1/V0 apart (fine comb).
  (b) The small part P has WIDE safe arcs (width O(1/7)).

The widest arc of the intersection is a CLUSTER-safe window that ALSO lies in a
P-safe arc.  Question: is the widest such window provably > 1/(7 Vmax)?

This probe extracts, per S3 set:
  - the witness tau* (where widest arc of G_A sits), as a Fraction
  - the LEFT and RIGHT boundary of that arc: which runner's tooth bounds each side
  - whether both bounds come from CLUSTER members (head-start claim) or from P
  - the cluster-collapse picture: spread of {frac(u tau*): u in L'} -- should be tiny
  - the "single-gap" index the cluster sits in at tau*
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

def widest_arc_info(A):
    """Return (width, lo, hi) of the widest safe arc (with circular wrap)."""
    sc=safe_components(A)
    if not sc: return F(0),None,None
    best=(F(0),None,None)
    for a,b in sc:
        if b-a>best[0]: best=(b-a,a,b)
    # wrap
    if sc[0][0]==0 and sc[-1][1]==1 and len(sc)>1:
        w=sc[0][1]+(1-sc[-1][0])
        if w>best[0]: best=(w,sc[-1][0],sc[0][1]+1)  # wrapping arc (hi>1)
    return best

def bounding_runner(A, tau_edge, side):
    """Which runner u in A has a tooth boundary exactly at tau_edge?
       side='right' means tau_edge is the RIGHT end of a safe arc => a tooth starts here
       (||u tau||=1/14 with u tau just entering danger). Return list of such u."""
    res=[]
    for u in A:
        val=nrm(u*tau_edge)
        if val==C:  # exactly on the tooth boundary
            res.append(u)
    return res

def probe(n_sets=120, seed=7):
    rng=random.Random(seed)
    flush("="*70)
    flush("STRUCTURE PROBE: witness arc of G_{S\\Vmax}: location, bounds, cluster collapse")
    flush("="*70)
    got=0; attempts=0
    cluster_bounds_both=0
    cluster_bounds_one=0
    cluster_bounds_none=0
    collapse_spreads=[]
    while got<n_sets and attempts<n_sets*40:
        attempts+=1
        S=make_S3_set(rng)
        if S is None: continue
        got+=1
        Vmax=max(S); A=set(S)-{Vmax}
        w,lo,hi=widest_arc_info(A)
        mid=(lo+hi)/2
        if mid>=1: mid-=1
        # identify cluster = large runners (>13) within 50 of max of A
        Amax=max(A); larges=[u for u in A if u>13]
        Lp=set(u for u in larges if Amax-u<=50) if larges else set()
        P=A-Lp
        # bounding runners at the two edges
        loe = lo if lo<1 else lo-1
        hie = hi if hi<1 else hi-1
        bl=bounding_runner(A,loe,'left')
        br=bounding_runner(A,hie,'right')
        # are bounds in cluster?
        bl_cl = any(u in Lp for u in bl)
        br_cl = any(u in Lp for u in br)
        if bl_cl and br_cl: cluster_bounds_both+=1
        elif bl_cl or br_cl: cluster_bounds_one+=1
        else: cluster_bounds_none+=1
        # cluster collapse at mid
        if Lp:
            fr=[nrm(u*mid) for u in Lp]
            collapse_spreads.append(max(fr)-min(fr))
    flush(f"\n  S3 sets probed: {got}")
    flush(f"  widest arc bounded by CLUSTER on BOTH sides: {cluster_bounds_both}")
    flush(f"  bounded by cluster on ONE side:              {cluster_bounds_one}")
    flush(f"  bounded by cluster on NEITHER side:          {cluster_bounds_none}")
    if collapse_spreads:
        flush(f"  cluster collapse spread max_u frac(u*mid)-min_u frac(u*mid): "
              f"min={float(min(collapse_spreads)):.5f} max={float(max(collapse_spreads)):.5f}")

if __name__=="__main__":
    probe(120)
