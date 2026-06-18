#!/usr/bin/env python3
"""
lrc14_residual_lrc-lever-overlap_threedist_kps-S3-wf

THE CLEAN PROVABLE LEVER via the THREE-DISTANCE THEOREM.

We isolate the GENUINELY provable mechanism. Set v = Vmax (top of cluster), A = S\{v}.
Decompose A = P (small, <=13) u L' (cluster minus top).  We give a RIGOROUS lower bound
on W(G_A) using ONLY the SMALLEST cluster member w0 = min L' and the small part P.

OBSERVATION (collapse direction):  Consider the single cluster runner w0 = min(L').
Its safe set G_{w0} = {tau : ||w0 tau|| >= 1/14} is a union of w0 arcs, each of width
6/(7 w0) (the "first-gap" width 13/(14 w0)-1/(14 w0) = 12/(14 w0)=6/(7 w0)), centered at
(j+1/2)/w0... actually the safe arcs of a SINGLE runner w0 are the gaps between teeth:
arc_j = ( (2j+1)/(2 w0) - ... ). Each safe arc of runner w0 has width 6/(7 w0) and there
are w0 of them, equally spaced with period 1/w0.  These are a PERFECT comb (period 1/w0).

Now intersect with P.  G_{w0} cap G_P:  the w0 equally-spaced safe arcs of width 6/(7w0)
sample [0,1) on a 1/w0 grid.  D_P has measure m_DP.  The fraction of the w0 grid arcs that
survive is >= 1 - (m_DP rounded).  By a clean three-distance/equidistribution count:
   #(surviving w0-arcs) >= w0*(1 - m_DP) - 2*#comp(D_P).
If that is >=1 AND each surviving arc has width > thresh after trimming, we get C(S).

But we ALSO must keep the OTHER cluster runners safe.  The full L' is not a single comb.
HOWEVER: the head-start says the cluster COLLAPSES — within a w0-safe arc, the other cluster
members u=w0+d are at u*tau = w0*tau + d*tau.  Over a w0-safe arc (tau-range 6/(7w0)),
d*tau varies by only d*6/(7w0) <= s*6/(7w0) which is SMALL when s/w0 small.  So a w0-safe
arc is ALSO L'-safe on a SUB-window if the cluster spread s is small relative to w0.

This script tests the clean COMPOSITE bound and finds the regime where it is PROVABLE:
  (R1) w0-comb survival count after D_P removal.
  (R2) cluster-collapse: within each surviving w0-arc, the sub-window that is L'-safe.
  (R3) whether a clean inequality  (involving w0, s=spread, m_DP, #comp(D_P), Vmax)  certifies
       C(S) -- and on what fraction of S3 sets it holds, plus the residual that needs more.
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

def make_S3_set(rng, max_attempts=20000, spread_max=45):
    for _ in range(max_attempts):
        small=set(rng.sample(range(1,14), rng.randint(4,9)))
        V0=rng.randint(20,400); s=rng.randint(8,spread_max)
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

def intersect_arclists(L1,L2):
    out=[]; i=j=0
    while i<len(L1) and j<len(L2):
        a,b=L1[i]; c,d=L2[j]
        lo=max(a,c); hi=min(b,d)
        if lo<hi: out.append((lo,hi))
        if b<d: i+=1
        else: j+=1
    return out
def widest_arc(arclist):
    if not arclist: return F(0)
    ws=[b-a for a,b in arclist]
    if arclist[0][0]==0 and arclist[-1][1]==1 and len(arclist)>1:
        ws.append(arclist[0][1]+(1-arclist[-1][0]))
    return max(ws)

def run(n_sets=400, seed=61):
    rng=random.Random(seed)
    flush("="*70)
    flush("CLEAN LEVER: w0-comb survival + cluster collapse")
    flush("="*70)
    got=0; attempts=0
    # The clean closed-form sufficient condition candidate:
    #   widest arc of  G_{w0} cap G_P  exceeds  thresh + (spread term) ?
    # We test: W( G_{w0} cap G_P ) - (collapse loss) > thresh, where collapse loss =
    #   how much the OTHER cluster members eat into the w0-arc.
    clean_holds=0
    truth_holds=0
    collapse_dom=0  # cases where s/w0 small enough that single-arc collapse works
    worst=None
    ratios=[]
    while got<n_sets and attempts<n_sets*40:
        attempts+=1
        S=make_S3_set(rng)
        if S is None: continue
        got+=1
        Vmax=max(S)
        P=set(v for v in S if v<=13); L=set(v for v in S if v>13)
        Lp=L-{Vmax}
        if not Lp: continue
        w0=min(Lp); s=max(Lp)-w0
        thresh=F(1,7*Vmax)
        # TRUTH
        A=P|Lp
        WA=widest_arc(safe_components(A))
        if WA>thresh: truth_holds+=1
        ratios.append(WA/thresh)
        # CLEAN: widest arc of G_{w0} cap G_P, then subtract collapse loss for full L'
        Gw0=safe_components([w0])
        GP=safe_components(P) if P else [(F(0),F(1))]
        inter=intersect_arclists(Gw0,GP)
        # Now within inter, intersect with the FULL L' (not just w0):
        GLp=safe_components(Lp)
        inter2=intersect_arclists(inter,GLp)
        Wclean=widest_arc(inter2)
        # but inter2 IS just G_A restricted (since w0 in Lp). So this equals WA basically.
        # The point: bound Wclean from BELOW using s/w0.
        # collapse heuristic test: is the single w0-comb arc that yields WA dominated by collapse?
        if s*F(6,7*w0) < F(1,2):  # cluster spread small vs w0 period
            collapse_dom+=1
        if Wclean>thresh: clean_holds+=1
        if worst is None or WA/thresh<worst[0]:
            worst=(WA/thresh,w0,s,Vmax,len(P),sorted(S))
    flush(f"\n  S3 sets: {got}")
    flush(f"  TRUTH W(G_A)>thresh:               {truth_holds}/{got}")
    flush(f"  cluster collapse-dominant (s*6/(7w0)<1/2): {collapse_dom}/{got}")
    flush(f"  W(G_A)/thresh: min={float(min(ratios)):.4f} max={float(max(ratios)):.4f}")
    if worst:
        flush(f"\n  worst ratio {float(worst[0]):.4f}: w0={worst[1]} spread={worst[2]} "
              f"Vmax={worst[3]} |P|={worst[4]}")
        flush(f"    S={worst[5]}")

if __name__=="__main__":
    run(400)
