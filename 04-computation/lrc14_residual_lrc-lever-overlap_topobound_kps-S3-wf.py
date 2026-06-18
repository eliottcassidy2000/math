#!/usr/bin/env python3
"""
lrc14_residual_lrc-lever-overlap_topobound_kps-S3-wf

THE TOPOLOGICAL COUNTING BOUND — search for the universal rigorous inequality.

We have on all S3 (tested): a common arc (G_P cap G_{L'}) of width > thresh EXISTS.
We want a PROVABLE inequality forcing this.

Setup: G_{L'} = union of "L'-arcs" (safe arcs of the cluster minus top).  D_P = small-part
danger = union of "P-teeth components".  A common arc of width>thresh fails to exist iff
EVERY L'-arc, after removing D_P, has all pieces of width <= thresh.

RIGOROUS DESTRUCTION COUNT.  Consider the L'-arcs of width > thresh (there are N_wide of them).
A wide L'-arc [a,b] is "killed" (yields no >thresh common subarc) only if D_P intersects it
such that no gap >thresh remains.  Count how D_P can do this:
  - The boundaries of D_P components inside (a,b) chop it.  If a D_P component lies ENTIRELY
    inside (a,b) it splits one arc into two pieces (both must be <=thresh to help kill).
  - To kill an arc of width w, need to cover it down to <=thresh-pieces: requires the
    D_P-measure inside the arc >= w - (#pieces)*thresh.

CLEANEST TRUE INVARIANT to test (Euler-characteristic style):
  Let B = total number of D_P-component ENDPOINTS that fall strictly inside L'-arcs.
  Each such endpoint can create at most one new "cut".  A wide L'-arc survives unless it is
  cut into all-small pieces.  The number of L'-arcs that can be killed is <= B (each kill
  consumes >=1 interior endpoint... or rather a component spanning an arc).

  Test the inequality:   N_wide  >  (number of D_P components that intersect SOME L'-arc).
  And the SHARPER:        meas(G_{L'} cap G_P) > 0  AND  it contains an arc>thresh.

ALSO probe the SHARPEST clean closed form using ONLY scalar data:
  Define  rho = meas(D_P)  (small-part danger measure, <1).
  Define  the L'-arcs have a MINIMUM width  wmin_L' and there are  n_arcs  of them with total
  measure meas(G_{L'}).
  A clean theorem:  if  meas(G_{L'}) - rho  >  thresh * (n_arcs + #comp(D_P)) ... test variants.

GOAL: find ONE scalar inequality, provable from (Vmax, w0, wT, P), that implies survivor.
We grid-search candidate inequalities and report which hold 100% on a large S3 sample.
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

def wide_arcs(arclist,thresh):
    res=[(a,b) for a,b in arclist if b-a>thresh]
    if arclist and arclist[0][0]==0 and arclist[-1][1]==1 and len(arclist)>1:
        if arclist[0][1]+(1-arclist[-1][0])>thresh:
            res.append((arclist[-1][0],arclist[0][1]+1))
    return res

def run(n_sets=600, seed=91):
    rng=random.Random(seed)
    flush("="*70)
    flush("TOPO-BOUND search for a universal scalar inequality => survivor")
    flush("="*70)
    got=0; attempts=0
    # candidate inequalities, count how many hold AND are tight to survivor existence
    cand_names=[
        "N_wide_Lp > comp_DP",                 # # wide L'-arcs > # P-danger components
        "N_wide_Lp > 2*comp_DP",
        "meas_GLp > meas_DP",                  # measure pigeonhole (known to fail)
        "N_wide_Lp - comp_DP >= 1",
        "N_wide_Lp > comp_DP + arcs_in_DP",    # subtract L'-arcs fully inside a D_P comp
    ]
    holds={k:0 for k in cand_names}
    survivor_exists=0
    min_Nwide_minus_comp=10**9
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
        DP=danger_merged(P) if P else []
        comp_DP=len(DP)
        warcs=wide_arcs(GLp,thresh)
        N_wide=len(warcs)
        mGLp=measure_safe(Lp) if Lp else F(1)
        mDP=1-(measure_safe(P) if P else F(1))
        # L'-arcs fully inside a D_P component (these are auto-killed but small)
        arcs_in_DP=0
        for a,b in warcs:
            aa=a if a<1 else a-1; bb=b if b<=1 else b-1
            for c,d in DP:
                if c<=aa and bb<=d: arcs_in_DP+=1; break
        # survivor: intersect
        GP=safe_from_danger(DP) if DP else [(F(0),F(1))]
        # widest of intersection
        inter=[]; i=j=0
        L1=[(a,b) for a,b in GLp]; L2=GP
        while i<len(L1) and j<len(L2):
            a,b=L1[i]; c,d=L2[j]
            lo=max(a,c); hi=min(b,d)
            if lo<hi: inter.append((lo,hi))
            if b<d: i+=1
            else: j+=1
        wsurv=max((b-a for a,b in inter), default=F(0))
        if inter and inter[0][0]==0 and inter[-1][1]==1 and len(inter)>1:
            wsurv=max(wsurv, inter[0][1]+(1-inter[-1][0]))
        surv = wsurv>thresh
        if surv: survivor_exists+=1
        vals={
            "N_wide_Lp > comp_DP": N_wide>comp_DP,
            "N_wide_Lp > 2*comp_DP": N_wide>2*comp_DP,
            "meas_GLp > meas_DP": mGLp>mDP,
            "N_wide_Lp - comp_DP >= 1": N_wide-comp_DP>=1,
            "N_wide_Lp > comp_DP + arcs_in_DP": N_wide>comp_DP+arcs_in_DP,
        }
        for k,v in vals.items():
            if v: holds[k]+=1
        m=N_wide-comp_DP
        if m<min_Nwide_minus_comp:
            min_Nwide_minus_comp=m; worst=(m,N_wide,comp_DP,sorted(S),float(thresh))
    flush(f"\n  S3 sets: {got}")
    flush(f"  survivor EXISTS: {survivor_exists}/{got}")
    for k in cand_names:
        flush(f"  [{holds[k]:>4}/{got}]  {k}")
    flush(f"\n  min (N_wide - comp_DP) = {min_Nwide_minus_comp}")
    if worst:
        flush(f"   at N_wide={worst[1]} comp_DP={worst[2]}")
        flush(f"   S={worst[3]}")

if __name__=="__main__":
    run(600)
