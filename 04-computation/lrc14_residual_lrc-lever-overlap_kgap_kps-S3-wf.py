#!/usr/bin/env python3
"""
lrc14_residual_lrc-lever-overlap_kgap_kps-S3-wf

EXTEND the Collapse Lemma: from single-gap (one cluster band) to the FULL cluster safe set,
using the cluster's OWN exact safe arcs (not the over-rigid single-gap approximation), and
search for an ELEMENTARY witness valid on the WHOLE of S3.

The single-gap construction failed because the true witness band straddles several gap indices.
The CORRECT object is G_{L'} itself.  The Collapse Lemma's value was: G_{L'}-arcs are wide and
many.  The remaining task: prove ONE of them lands in G_P.

ELEMENTARY UNIVERSAL ATTEMPT (the real lever).  Use the FIRST cluster member w0 only as a
*ruler*.  Its safe arcs I_0,...,I_{w0-1} are EXACTLY periodic (period 1/w0, each width 6/(7w0)).
Within each I_j, ALL of L' is safe on the sub-window
   K_j = { tau in I_j : frac(u tau) in [1/14,13/14] for all u in L' }.
K_j is NONEMPTY and we can compute its width exactly. Over the w0 periods, the K_j are NEARLY
congruent (three-distance): the cluster offset pattern d_u*tau drifts by d_u/w0 per period, so
the family {K_j} takes only a BOUNDED number (<= 3, by three-distance) of distinct widths/shapes.

CLAIM TO TEST:  At least a constant fraction of the K_j have width > thresh, and they are
equidistributed mod the P-danger, so >=1 lands in G_P.  We test the EXACT counts:
  - n_K = #{j : width(K_j) > thresh}
  - n_K_Psafe = #{j : K_j has a length-thresh sub-window in G_P}
  - whether n_K_Psafe >= 1 on ALL S3 (this is the TRUE single-runner-ruler coverage, should
    beat the rigid single-gap version since K_j uses the cluster's real safe window not the
    monotone bound).

ALSO: the THREE-DISTANCE count of distinct K_j shapes (validates the bounded-pattern claim).
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

def Kj_in_period(j, w0, Lp, thresh, GP):
    """The L'-safe sub-window(s) inside the w0-safe arc I_j = ((j+1/14)/w0,(j+13/14)/w0).
       Compute exactly by intersecting I_j with G_{Lp}. Return list of arcs and whether any
       length-thresh piece is also P-safe."""
    lo=F(14*j+1,14*w0); hi=F(14*j+13,14*w0)
    if hi>=1: hi=F(1)  # clamp (period j=w0-1 may approach 1)
    if lo>=hi: return [], False
    # G_{Lp} restricted to [lo,hi]
    GLp=safe_components(Lp) if Lp else [(F(0),F(1))]
    # intersect [lo,hi] with GLp and GP
    pieces=[]
    for a,b in GLp:
        x=max(lo,a); y=min(hi,b)
        if x<y: pieces.append((x,y))
    # P-safe check: any piece has a length-thresh sub-window in GP
    psafe=False
    for x,y in pieces:
        if y-x<=thresh: continue
        # intersect [x,y] with GP, check a >thresh piece
        for a,b in GP:
            xx=max(x,a); yy=min(y,b)
            if yy-xx>thresh: psafe=True; break
        if psafe: break
    return pieces, psafe

def run(n_sets=500, seed=111):
    rng=random.Random(seed)
    flush("="*70)
    flush("K_j RULER (single cluster member w0 as ruler, exact cluster window)")
    flush("="*70)
    got=0; attempts=0
    ruler_certified=0   # some K_j P-safe
    residual=[]
    nK_data=[]
    while got<n_sets and attempts<n_sets*50:
        attempts+=1
        S=make_S3_set(rng)
        if S is None: continue
        got+=1
        Vmax=max(S)
        P=frozenset(v for v in S if v<=13); L=sorted(v for v in S if v>13)
        Lp=[v for v in L if v!=Vmax]
        w0=min(Lp); thresh=F(1,7*Vmax)
        GP=safe_components(P) if P else [(F(0),F(1))]
        found=False; nK=0
        for j in range(0,w0):
            pieces,psafe=Kj_in_period(j,w0,Lp,thresh,GP)
            if any(y-x>thresh for x,y in pieces): nK+=1
            if psafe: found=True; break
        nK_data.append(nK)
        if found: ruler_certified+=1
        else: residual.append(sorted(S))
    flush(f"\n  S3 sets: {got}")
    flush(f"  RULER certified (some K_j P-safe):  {ruler_certified}/{got}  "
          f"({100*ruler_certified/got:.1f}%)")
    flush(f"  RESIDUAL: {len(residual)}")
    for r in residual[:8]: flush(f"    {r}")

if __name__=="__main__":
    run(500, seed=111)
