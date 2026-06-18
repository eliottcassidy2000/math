#!/usr/bin/env python3
"""
lrc14_residual_lrc-lever-overlap_finegrid_kps-S3-wf

THE FINE-GRID PIGEONHOLE — the correct rigorous mechanism.

The survivor exists because G_{L'} is a FINE comb: its safe arcs sample [0,1) densely enough
that EVERY sufficiently-wide P-safe arc must contain at least one L'-safe sub-arc of width>thresh.

RIGOROUS PIGEONHOLE (candidate THEOREM).
  Let G = max gap of the cluster comb, defined as the largest tau-interval CONTAINING NO point
  whose entire thresh-neighborhood is L'-safe.  More usable: let
     Gmax = the largest gap between consecutive L'-safe-arcs' "thresh-cores"
            (a thresh-core of an L'-arc [a,b] is [a, b-thresh], the set of left-endpoints of
             length-thresh windows inside it).
  If a P-safe arc [p,q] has width q-p > Gmax + thresh, then it must contain a full length-thresh
  window that is L'-safe -> survivor.  PROVABLE pigeonhole IF we can bound Gmax.

  Equivalent cleaner statement.  Let D_{L'} = cluster danger (complement of G_{L'}).
  Cluster comb fineness:  the longest run of cluster-danger is <= (one cluster tooth-block).
  By THM-526 arc-width lemma applied to the cluster, the longest cluster-DANGER run is exactly
  the longest danger run of the SLOWEST cluster member... no -- it's the union.

KEY QUESTION to settle exactly:
  (Q) What is  Gmax = longest gap in [0,1) where NO length-thresh L'-safe window starts?
      Compare to  WPmax = widest P-safe arc, and to  WP2 = 2nd widest, etc.
  If  Gmax + thresh < WPmax  for ALL S3  -> clean pigeonhole proof of survivor in WIDEST P-arc.
  We saw E3 (a cruder version) failed 60/400.  Test the SHARP version with thresh-cores.

  ALSO test the "covering by several P-arcs" version: the total P-safe measure is spread over
  multiple arcs; pigeonhole over the SUM.  Specifically:
     If  sum over P-safe arcs of  max(0, width - (Gmax+thresh))  > 0,  survivor exists.
  Test whether this SUM is always > 0 (it should match survivor existence).

Report the exact Gmax distribution and which pigeonhole form is universal.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random, sys

def flush(*a):
    print(*a); sys.stdout.flush()
C=F(1,14)

def teeth_intervals(u,h=C):
    iv=[]
    for j in range(0,u):
        c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv
def merged_intervals(iv):
    iv=sorted(iv); merged=[]
    for a,b in iv:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    return merged
def danger_merged(A,h=C):
    iv=[]
    for u in A: iv+=teeth_intervals(u,h)
    return merged_intervals(iv)
def safe_from_danger(merged):
    safe=[]; prev=F(0)
    for a,b in merged:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<1: safe.append((prev,F(1)))
    return safe
def safe_components(A,h=C): return safe_from_danger(danger_merged(A,h))
def measure_safe(A): return sum(b-a for a,b in safe_components(A))
def is_covering(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def gcd_all(S): return reduce(gcd,S)

def make_S3_set(rng, max_attempts=20000, spread_max=80, vlo=14, vhi=400):
    for _ in range(max_attempts):
        small=set(rng.sample(range(1,14), rng.randint(4,9)))
        V0=rng.randint(vlo,vhi); s=rng.randint(6,spread_max)
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

def thresh_cores(GLp, thresh):
    """For each L'-safe arc [a,b] with b-a>thresh, the core [a, b-thresh] of valid left-endpoints.
       Return merged list of cores (where a length-thresh L'-safe window can START)."""
    cores=[]
    for a,b in GLp:
        if b-a>thresh: cores.append((a,b-thresh))
    # handle wrap arc
    if GLp and GLp[0][0]==0 and GLp[-1][1]==1 and len(GLp)>1:
        w=GLp[0][1]+(1-GLp[-1][0])
        if w>thresh:
            cores.append((GLp[-1][0], GLp[0][1]+1-thresh))  # may exceed 1
    return merged_intervals(cores)

def max_gap_complement(cores):
    """Largest gap in [0,1) NOT covered by cores (circularly). cores may have entries with hi>1."""
    # normalize to [0,1) by splitting wrap
    norm=[]
    for a,b in cores:
        if b<=1: norm.append((a,b))
        else:
            norm.append((a,F(1)));
            bb=b-1
            if bb>0: norm.append((F(0),bb))
    norm=merged_intervals(norm)
    if not norm: return F(1)
    # gaps between consecutive, plus wrap
    gaps=[]
    for i in range(len(norm)-1):
        gaps.append(norm[i+1][0]-norm[i][1])
    gaps.append(norm[0][0]+1-norm[-1][1])  # wrap
    return max(gaps) if gaps else F(0)

def run(n_sets=600, seed=101):
    rng=random.Random(seed)
    flush("="*70)
    flush("FINE-GRID PIGEONHOLE: Gmax (cluster comb max gap of thresh-cores) vs P-safe arcs")
    flush("="*70)
    got=0; attempts=0
    pig_widest=0        # Gmax+thresh < widest P-arc
    pig_sum=0           # sum over P-arcs of max(0, w-(Gmax+thresh)) > 0
    survivor=0
    Gmax_data=[]
    worst=None
    while got<n_sets and attempts<n_sets*50:
        attempts+=1
        S=make_S3_set(rng)
        if S is None: continue
        got+=1
        Vmax=max(S)
        P=frozenset(v for v in S if v<=13); L=sorted(v for v in S if v>13)
        Lp=[v for v in L if v!=Vmax]
        thresh=F(1,7*Vmax)
        GLp=safe_components(Lp) if Lp else [(F(0),F(1))]
        GP=safe_components(P) if P else [(F(0),F(1))]
        cores=thresh_cores(GLp,thresh)
        Gmax=max_gap_complement(cores)
        Gmax_data.append(Gmax)
        # widest P-arc (with wrap)
        wp=[b-a for a,b in GP]
        if GP and GP[0][0]==0 and GP[-1][1]==1 and len(GP)>1:
            wp.append(GP[0][1]+(1-GP[-1][0]))
        WPmax=max(wp) if wp else F(0)
        if Gmax+thresh < WPmax: pig_widest+=1
        # sum form
        ssum=sum(max(F(0), w-(Gmax+thresh)) for w in wp)
        if ssum>0: pig_sum+=1
        # actual survivor
        # widest of intersection
        inter=[]; i=j=0
        L1=[(a,b) for a,b in GLp]; L2=GP
        while i<len(L1) and j<len(L2):
            a,b=L1[i]; c,d=L2[j]
            lo=max(a,c); hi=min(b,d)
            if lo<hi: inter.append((lo,hi))
            if b<d: i+=1
            else: j+=1
        wsurv=max((b-a for a,b in inter),default=F(0))
        if inter and inter[0][0]==0 and inter[-1][1]==1 and len(inter)>1:
            wsurv=max(wsurv, inter[0][1]+(1-inter[-1][0]))
        if wsurv>thresh: survivor+=1
        if worst is None or (Gmax+thresh-WPmax) > worst[0]:
            worst=(Gmax+thresh-WPmax,float(Gmax),float(WPmax),float(thresh),sorted(S))
    flush(f"\n  S3 sets: {got}")
    flush(f"  survivor EXISTS:                            {survivor}/{got}")
    flush(f"  PIGEONHOLE widest (Gmax+thresh<WPmax):      {pig_widest}/{got}")
    flush(f"  PIGEONHOLE sum   (sum max(0,w-Gmax-thresh)>0): {pig_sum}/{got}")
    flush(f"  Gmax: min={float(min(Gmax_data)):.5f} max={float(max(Gmax_data)):.5f} "
          f"avg={float(sum(Gmax_data)/len(Gmax_data)):.5f}")
    if worst:
        flush(f"\n  worst (Gmax+thresh-WPmax={worst[0]}):  Gmax={worst[1]:.5f} WPmax={worst[2]:.5f}")
        flush(f"    S={worst[4]}")

if __name__=="__main__":
    run(600)
