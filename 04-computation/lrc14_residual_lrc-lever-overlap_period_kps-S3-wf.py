#!/usr/bin/env python3
"""
lrc14_residual_lrc-lever-overlap_period_kps-S3-wf

THE ALMOST-PERIODIC / BOUNDED-PATTERN core — the heart of the lrc-lever angle.

We test the central structural claim that would give a RIGOROUS proof:

CLAIM (cluster comb).  Let L' = {V0 + d : d in Delta} (Delta = bounded offset set,
0 = min, max offset = s).  Consider the cluster-safe set
   G_{L'} = {tau : ||u tau|| >= 1/14  for all u in L'}.
Then in EACH "super-cell" of length 1/g (g = gcd-related period) the safe set has a
BOUNDED, V0-independent NUMBER of wide arcs, and the WIDEST is ~ (1 - eps)/(7 V0) wide,
with eps depending only on Delta.  As V0 grows, G_{L'} -> a fine comb that EQUIDISTRIBUTES.

To PROVE C(S) we want:  among the >=24 wide arcs of G_{L'}, at least one avoids D_P.
The rigorous lever: G_{L'} is "spread" — its wide arcs are NOT clustered; they sample
[0,1) almost uniformly.  Quantify by:

  (E1) The MAXIMAL GAP between consecutive wide-arc CENTERS of G_{L'} (a "dispersion").
       If max-gap < (min P-safe arc length), pigeonhole forces a wide L'-arc into each
       P-safe arc -> survivor.
  (E2) meas(G_P) (total P-safe measure) and the # of P-safe arcs.  A wide L'-arc center
       in a P-safe arc, far enough from its ends, survives.

  (E3) THE CLEAN SUFFICIENT THEOREM candidate:
       If (max gap between consecutive G_{L'} wide-arc midpoints) < (widest P-safe arc) - 2*thresh,
       then some wide L'-arc midpoint lands in the widest P-safe arc with thresh-margin -> survivor.
       Test this inequality exactly on S3.

  (E4) Even cleaner — relate to meas(D_P):  the wide L'-arc midpoints form a set X of >=24
       points.  D_P has measure m_DP < 1.  If X equidistributes, expect (1-m_DP)*|X| of them
       in G_P.  Test the ACTUAL count of wide-L'-midpoints landing in G_P, and whether it is
       always >= 1 (it must be, since survivor exists).  Then find the SHARPEST provable lower
       bound: #(midpoints in G_P) >= |X| * (1 - m_DP) - (boundary slack)?

This pins down WHICH dispersion statistic gives a clean, provable inequality.
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
def measure_safe(A): return sum(b-a for a,b in safe_components(A))
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

def wide_arc_mids(arclist,thresh):
    mids=[]
    for a,b in arclist:
        if b-a>thresh: mids.append(((a+b)/2, b-a))
    if arclist and arclist[0][0]==0 and arclist[-1][1]==1 and len(arclist)>1:
        w=arclist[0][1]+(1-arclist[-1][0])
        if w>thresh:
            m=(arclist[-1][0]+arclist[0][1]+1)/2
            if m>=1: m-=1
            mids.append((m,w))
    return sorted(mids)

def max_gap(points):
    """max circular gap between consecutive points in [0,1)."""
    if not points: return F(1)
    pts=sorted(points)
    g=F(0)
    for i in range(len(pts)-1):
        g=max(g,pts[i+1]-pts[i])
    g=max(g, pts[0]+1-pts[-1])  # wrap
    return g

def in_safe(t, safe):
    for a,b in safe:
        if a<=t<b: return True
    return False

def run(n_sets=400, seed=53):
    rng=random.Random(seed)
    flush("="*70)
    flush("DISPERSION: max-gap of wide L'-arc midpoints vs widest P-safe arc")
    flush("="*70)
    got=0; attempts=0
    E3_holds=0
    mids_in_GP_min=10**9
    worst_E3=None
    maxgap_data=[]; widestP_data=[]
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
        GP=safe_components(P) if P else [(F(0),F(1))]
        mids=wide_arc_mids(GLp,thresh)
        midpts=[m for m,w in mids]
        mg=max_gap(midpts)
        widestP=max((b-a for a,b in GP), default=F(0))
        # wrap for P
        if GP and GP[0][0]==0 and GP[-1][1]==1 and len(GP)>1:
            widestP=max(widestP, GP[0][1]+(1-GP[-1][0]))
        maxgap_data.append(mg); widestP_data.append(widestP)
        # E3 test
        if mg < widestP - 2*thresh: E3_holds+=1
        else:
            if worst_E3 is None or (widestP-2*thresh-mg) < worst_E3[0]:
                worst_E3=(widestP-2*thresh-mg,float(mg),float(widestP),sorted(S))
        # how many wide L'-midpoints land in G_P
        cnt=sum(1 for m in midpts if in_safe(m,GP))
        mids_in_GP_min=min(mids_in_GP_min,cnt)
    flush(f"\n  S3 sets: {got}")
    flush(f"  E3  (maxgap < widestP - 2 thresh):  {E3_holds}/{got}")
    flush(f"  min #(wide L'-midpoints inside G_P): {mids_in_GP_min}")
    flush(f"  maxgap of midpoints: min={float(min(maxgap_data)):.5f} max={float(max(maxgap_data)):.5f}")
    flush(f"  widest P-safe arc:   min={float(min(widestP_data)):.5f} max={float(max(widestP_data)):.5f}")
    if worst_E3:
        flush(f"\n  worst E3 (slack {worst_E3[0]:.5f}): maxgap={worst_E3[1]:.5f} widestP={worst_E3[2]:.5f}")
        flush(f"    S={worst_E3[3]}")

if __name__=="__main__":
    run(400)
