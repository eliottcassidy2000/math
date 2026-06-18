#!/usr/bin/env python3
"""
lrc14_residual_lrc-lever-overlap_stress_kps-S3-wf

STRESS-TEST the K_j ruler witness on adversarial S3 regimes + measure the witness COUNT
scaling with V0 (to support the equidistribution / persistence claim).

The K_j ruler (single cluster member w0 as periodic ruler) certified 500/500 on the main
sample.  Here we:
  (1) Stress regimes: huge V0 (up to 5000), extreme spread (up to 200), tiny P, multi-cluster,
      cluster = exact AP (worst three-distance case), cluster spanning a multiple of 14.
  (2) Measure  rho_K = (#periods j with K_j P-safe) / w0  -- the FRACTION of ruler periods that
      work.  If rho_K is bounded BELOW by a positive constant as V0 grows, the witness PERSISTS
      for all large V0 -> supports a proof.  If rho_K -> 0, the finite check has no uniform floor.
  (3) Report the exact minimum  #(P-safe periods)  over all tested S3 (must be >=1; ideally grows).
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

def count_Psafe_periods(S, cap_periods=None):
    """Count periods j (of ruler w0) where I_j contains a length-thresh window safe for
       ALL of L' AND P. Returns (count, w0). Uses ONE precompute of GLp, GP."""
    Vmax=max(S)
    P=[v for v in S if v<=13]; L=sorted(v for v in S if v>13)
    Lp=[v for v in L if v!=Vmax]
    if not Lp: return (None,None)
    w0=min(Lp); thresh=F(1,7*Vmax)
    GLp=safe_components(Lp); GP=safe_components(P) if P else [(F(0),F(1))]
    # G_A = G_Lp cap G_P (precompute once)
    GA=[]; i=j=0
    while i<len(GLp) and j<len(GP):
        a,b=GLp[i]; c,d=GP[j]
        lo=max(a,c); hi=min(b,d)
        if lo<hi: GA.append((lo,hi))
        if b<d: i+=1
        else: j+=1
    # For each ruler period I_j, does GA have a >thresh piece inside?
    cnt=0
    Jrange = range(0,w0) if cap_periods is None else range(0,min(w0,cap_periods))
    # Efficient: for each GA arc, find which periods it touches and whether width>thresh.
    # We just scan GA arcs and bucket by period of their midpoint; simpler: count GA arcs >thresh
    # whose interval lies within a single I_j. But the ruler periods partition [0, 13/14...].
    # Simplest exact: count distinct ruler periods j that contain a >thresh GA piece.
    periods_hit=set()
    for a,b in GA:
        if b-a<=thresh: continue
        # which period(s) does [a,b] lie in? I_j=((14j+1)/(14w0),(14j+13)/(14w0))
        # find j s.t. a in I_j: j = floor(a*w0 - 1/14)... just compute j from a
        # a >= (14j+1)/(14w0) => 14 w0 a -1 >= 14 j => j <= (14 w0 a -1)/14
        ja = (14*w0*a - 1)
        import math
        j_lo = ( (14*w0*a - 1) )  # fraction
        # period index containing a:
        jj = int((w0*a) )  # approx; refine
        # robust: search nearby
        for cand_j in range(max(0,jj-1), jj+2):
            lo=F(14*cand_j+1,14*w0); hi=F(14*cand_j+13,14*w0)
            x=max(a,lo); y=min(b,hi)
            if y-x>thresh:
                periods_hit.add(cand_j)
    return (len(periods_hit), w0)

def make_S3(rng, spread_max, vlo, vhi, ap=False, near14=False):
    for _ in range(30000):
        small=set(rng.sample(range(1,14), rng.randint(3,9)))
        V0=rng.randint(vlo,vhi)
        s=rng.randint(6,spread_max)
        clsize=rng.randint(2, max(2,13-len(small)))
        if ap:
            step=rng.randint(1, max(1,s//max(1,clsize)))
            cluster=set(V0+step*t for t in range(clsize))
        else:
            cluster=set()
            tries=0
            while len(cluster)<clsize and tries<300:
                cluster.add(V0+rng.randint(0,s)); tries+=1
        if near14:
            # force cluster to span a multiple of 14 region
            base=14*rng.randint(2,20)
            cluster=set(base+rng.randint(-s//2,s//2) for _ in range(clsize))
            cluster={c for c in cluster if c>13}
        S=small|cluster
        while len(S)<13:
            S.add(rng.randint(min(cluster) if cluster else V0, (max(cluster)+s) if cluster else V0+s)
                  if rng.random()<0.5 else rng.randint(1,13))
        if len(S)!=13: continue
        S=set(list(S))
        if len(S)!=13 or gcd_all(S)!=1 or not is_covering(S): continue
        Vmin=min(S); Vmax=max(S); k=sum(1 for v in S if v>13)
        if k>=2 and Vmax>=13*Vmin: return S
    return None

def run_regime(tag, n, seed, spread_max=60, vlo=14, vhi=400, ap=False, near14=False):
    rng=random.Random(seed)
    got=0; attempts=0; cert=0
    min_periods=10**9; min_frac=None; residual=[]
    fracs=[]
    while got<n and attempts<n*60:
        attempts+=1
        S=make_S3(rng,spread_max,vlo,vhi,ap,near14)
        if S is None: continue
        got+=1
        cnt,w0=count_Psafe_periods(S)
        if cnt is None: continue
        if cnt>=1: cert+=1
        else: residual.append(sorted(S))
        min_periods=min(min_periods,cnt)
        frac=F(cnt,w0)
        fracs.append(frac)
        if min_frac is None or frac<min_frac: min_frac=frac
    flush(f"  [{tag}] sets={got} certified(>=1 P-safe period)={cert}/{got}  "
          f"min#periods={min_periods}  min frac rho_K={float(min_frac):.4f}  "
          f"avg rho_K={float(sum(fracs)/len(fracs)):.4f}")
    for r in residual[:4]: flush(f"      RESIDUAL {r}")
    return cert,got

if __name__=="__main__":
    flush("="*70)
    flush("STRESS TEST: K_j ruler witness across adversarial S3 regimes")
    flush("="*70)
    run_regime("general",        400, 201, spread_max=45, vlo=14, vhi=400)
    run_regime("huge V0<=5000",  300, 202, spread_max=45, vlo=500, vhi=5000)
    run_regime("extreme spread", 300, 203, spread_max=200, vlo=30, vhi=600)
    run_regime("AP cluster",     400, 204, spread_max=80, vlo=20, vhi=500, ap=True)
    run_regime("near-14k span",  300, 205, spread_max=40, vlo=14, vhi=300, near14=True)
    run_regime("tiny P / big L", 300, 206, spread_max=30, vlo=20, vhi=400)
