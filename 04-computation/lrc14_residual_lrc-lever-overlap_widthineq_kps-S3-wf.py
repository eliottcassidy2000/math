#!/usr/bin/env python3
"""
lrc14_residual_lrc-lever-overlap_widthineq_kps-S3-wf

VERIFY the WIDTH-INEQ as a clean PROVED fact and find its exact hypothesis.

CLAIM (WIDTH-INEQ).  For any S3 covering 13-set, with v=Vmax, A=S\{v}, P=small part,
L'=cluster minus top, w0=min L', wT=max L', the FIRST-GAP cluster window
   J_0 = [ 1/(14 w0) , 13/(14 wT) ]
has width  (13 w0 - wT)/(14 w0 wT)  and this EXCEEDS thresh = 1/(7 Vmax) = 2/(14 Vmax)  iff
   (13 w0 - wT) * Vmax  >  2 * w0 * wT.                                       (*)

When does (*) hold UNCONDITIONALLY on S3?  We have Vmax >= wT (Vmax is the max of S, wT=max L'
which is <= Vmax; actually Vmax could equal wT+something or be a separate large runner).  Two cases:
  - If Vmax is IN the cluster top region: wT < Vmax (since we removed Vmax from L).  Then
    (13 w0 - wT) Vmax >= (13 w0 - wT) wT.  (*) reduces to needing 13 w0 - wT > 2 w0, i.e.
    11 w0 > wT, i.e. cluster spread ratio wT/w0 < 11.
  - The cluster spread: wT = w0 + s.  wT/w0 < 11  <=>  s < 10 w0.  For ANY cluster (s <= ~tens,
    w0 >= 14) this is MASSIVELY satisfied.  So (*) holds whenever 11 w0 > wT AND Vmax >= wT,
    i.e. essentially always on S3.

We VERIFY (*) holds on a large S3 sample, find its exact min margin, and confirm the clean
hypothesis  "11 w0 > wT"  (cluster internal spread ratio < 11) is SUFFICIENT and HOLDS on S3.
We ALSO check: is wT/w0 < 11 ALWAYS true for an S3 cluster?  The cluster is the set of large
runners; if there are large runners spread over a factor >11 they wouldn't be 'a cluster'.
Define L' more carefully as the GROUP near Vmax; test whether some S3 sets have wT/w0 >= 11
(a 'wide cluster' -- the genuinely hard remainder).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import random, sys

def flush(*a):
    print(*a); sys.stdout.flush()
C=F(1,14)

def is_covering(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def gcd_all(S): return reduce(gcd,S)

def make_S3(rng, spread_max, vlo, vhi):
    for _ in range(30000):
        small=set(rng.sample(range(1,14), rng.randint(3,9)))
        V0=rng.randint(vlo,vhi); s=rng.randint(1,spread_max)
        clsize=rng.randint(2, max(2,13-len(small)))
        cluster=set(); tries=0
        while len(cluster)<clsize and tries<300:
            cluster.add(V0+rng.randint(0,s)); tries+=1
        S=small|cluster
        while len(S)<13:
            S.add(rng.randint(V0,V0+s) if rng.random()<0.5 else rng.randint(1,13))
        if len(S)!=13: continue
        S=set(list(S))
        if len(S)!=13 or gcd_all(S)!=1 or not is_covering(S): continue
        Vmin=min(S); Vmax=max(S); k=sum(1 for v in S if v>13)
        if k>=2 and Vmax>=13*Vmin: return sorted(S)
    return None

def run(n=2000, seed=501, spread_max=200, vlo=14, vhi=2000):
    rng=random.Random(seed)
    got=0; attempts=0
    star_holds=0
    spread11=0  # wT/w0 < 11 (cluster = ALL large runners)
    min_margin=None
    wide_cluster=[]
    while got<n and attempts<n*60:
        attempts+=1
        S=make_S3(rng,spread_max,vlo,vhi)
        if S is None: continue
        got+=1
        Vmax=max(S)
        L=sorted(v for v in S if v>13); Lp=[v for v in L if v!=Vmax]
        if not Lp: continue
        w0=min(Lp); wT=max(Lp)
        # (*) : (13 w0 - wT) Vmax > 2 w0 wT
        lhs=(13*w0-wT)*Vmax; rhs=2*w0*wT
        if lhs>rhs: star_holds+=1
        margin=F(lhs-rhs, rhs) if rhs>0 else F(0)
        if min_margin is None or margin<min_margin: min_margin=margin
        if 11*w0>wT: spread11+=1
        else: wide_cluster.append((S,w0,wT,F(wT,w0)))
    flush(f"\n  [spread<={spread_max}, V0 in [{vlo},{vhi}]] S3 sets={got}")
    flush(f"    WIDTH-INEQ (*) holds:        {star_holds}/{got}")
    flush(f"    cluster spread wT/w0 < 11:   {spread11}/{got}")
    flush(f"    min margin of (*) = {float(min_margin):.4f}  (exact {min_margin})")
    flush(f"    'wide cluster' (wT/w0>=11): {len(wide_cluster)}")
    for S,w0,wT,r in wide_cluster[:5]:
        flush(f"      {S}  w0={w0} wT={wT} ratio={float(r):.2f}")

if __name__=="__main__":
    flush("="*70)
    flush("WIDTH-INEQ verification across wide regimes")
    flush("="*70)
    run(2000, 501, 200, 14, 2000)
    run(1500, 502, 500, 14, 3000)   # very loose clusters
    run(1500, 503, 30, 14, 60)      # small V0 boundary
