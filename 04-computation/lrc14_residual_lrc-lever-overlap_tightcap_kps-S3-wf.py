#!/usr/bin/env python3
"""
lrc14_residual_lrc-lever-overlap_tightcap_kps-S3-wf

TIGHTEN the destruction-capacity bound.

A wide L'-arc [a,b] (width w > thresh = 1/(7Vmax)) yields a common arc with G_P UNLESS
D_P (small-part danger) covers all but margins, i.e. there is no length-thresh P-safe window
inside [a,b].  We over-counted with "2 arcs per D_P component".

Precise destruction model:
  A length-thresh P-safe window inside [a,b] fails to exist  <=>  D_P teeth tile [a,b]
  leaving every gap < thresh.  Since thresh is TINY (1/(7Vmax)) and a single P-tooth of
  runner u<=13 has width 1/(7u) >= 1/91 >> thresh, ONE P-tooth covers a length 1/(7u)
  chunk.  An L'-arc of width w survives a single P-tooth as long as w - (tooth overlap) > thresh
  on at least one side.

KEY exact quantities to settle the right capacity:
  (Q1) For each S3 set: ACTUAL number K of wide L'-arcs DESTROYED (no thresh-window survives).
       Compare K to #comp(D_P) and to measure-based predictions.
  (Q2) Is it ALWAYS true that K <= #comp(D_P)?  (each connected P-danger component destroys
       at most ONE wide L'-arc).  If so, then  N_wide > #comp(D_P)  =>  survivor exists.  CLEAN.
  (Q3) Does the WEAKER inequality N_wide > #comp(D_P) hold for ALL S3 (we saw 399/400)?
       Examine the one failure; is it a genuine S3 set, and does C(S) still hold via another v?

  (Q4) Better capacity: total P-danger MEASURE m_DP.  An L'-arc of width w needs P-danger to
       cover length > w - thresh to be destroyed (roughly).  So #destroyed <= m_DP / (min arc
       width).  Since min wide-arc width ~ thresh... not clean.  Try: #destroyed-arcs * thresh
       <= m_DP + (#comp boundary effects).

This script:
  - computes K (true destroyed count) exactly
  - checks K <= #comp(D_P) on all S3
  - examines the failure set of N_wide > #comp(D_P)
  - reports the ROBUST clean inequality that holds universally on the sample.
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

def wide_arcs(arclist,thresh):
    """List wide arcs (a,b) width>thresh, treating wrap arc as separate entry (a,b) with b may>1."""
    res=[(a,b) for a,b in arclist if b-a>thresh]
    if arclist and arclist[0][0]==0 and arclist[-1][1]==1 and len(arclist)>1:
        if arclist[0][1]+(1-arclist[-1][0])>thresh:
            res.append((arclist[-1][0], arclist[0][1]+1))  # wrapping
    return res

def arc_survives(arc, DP_merged, thresh):
    """Does wide L'-arc [a,b] contain a P-safe subwindow of width>thresh?
       i.e. does (arc minus DP) have a piece of width>thresh?"""
    a,b=arc
    # build DP within [a,b] (handle wrap by shifting)
    # represent arc on real line; map DP teeth into [a,b] including +1 copies for wrap
    pieces=[]
    prev=a
    # gather DP intervals intersecting [a,b], also +1 shifted
    segs=[]
    for c,d in DP_merged:
        for shift in (F(0),F(1)):
            cc,dd=c+shift,d+shift
            lo=max(a,cc); hi=min(b,dd)
            if lo<hi: segs.append((lo,hi))
    segs=merged_intervals(segs)
    for c,d in segs:
        if c>prev and c-prev>thresh: return True
        prev=max(prev,d)
    if b>prev and b-prev>thresh: return True
    return False

def run(n_sets=400, seed=41):
    rng=random.Random(seed)
    flush("="*70)
    flush("TIGHT CAPACITY: true destroyed count K, test K <= #comp(D_P)")
    flush("="*70)
    got=0; attempts=0
    K_le_comp=0
    Nwide_gt_comp=0
    survivor_exists=0
    worst_margin=None
    failure_examples=[]
    Kvals=[]
    while got<n_sets and attempts<n_sets*40:
        attempts+=1
        S=make_S3_set(rng)
        if S is None: continue
        got+=1
        Vmax=max(S)
        P=set(v for v in S if v<=13); L=set(v for v in S if v>13)
        Lp=L-{Vmax}
        thresh=F(1,7*Vmax)
        GLp=safe_components(Lp) if Lp else [(F(0),F(1))]
        DP=danger_merged(P) if P else []
        ncomp=len(DP)
        warcs=wide_arcs(GLp,thresh)
        Nwide=len(warcs)
        survivors=[arc for arc in warcs if arc_survives(arc,DP,thresh)]
        K=Nwide-len(survivors)   # destroyed
        Kvals.append(K)
        if K<=ncomp: K_le_comp+=1
        if Nwide>ncomp: Nwide_gt_comp+=1
        if len(survivors)>=1: survivor_exists+=1
        else:
            failure_examples.append(sorted(S))
        m=Nwide-ncomp
        if worst_margin is None or m<worst_margin[0]:
            worst_margin=(m,Nwide,ncomp,K,len(survivors),sorted(S))
    flush(f"\n  S3 sets: {got}")
    flush(f"  survivor (common arc>thresh) EXISTS:  {survivor_exists}/{got}")
    flush(f"  K (destroyed) <= #comp(D_P):          {K_le_comp}/{got}")
    flush(f"  N_wide > #comp(D_P):                  {Nwide_gt_comp}/{got}")
    flush(f"  max K (destroyed arcs): {max(Kvals)}   avg K: {sum(Kvals)/len(Kvals):.2f}")
    if worst_margin:
        flush(f"\n  tightest N_wide-#comp: {worst_margin[0]}  "
              f"N_wide={worst_margin[1]} #comp={worst_margin[2]} K={worst_margin[3]} surv={worst_margin[4]}")
        flush(f"    S={worst_margin[5]}")
    if failure_examples:
        flush(f"\n  *** NO-SURVIVOR cases: {len(failure_examples)} ***")
        for fx in failure_examples[:5]: flush(f"    {fx}")

if __name__=="__main__":
    run(400)
