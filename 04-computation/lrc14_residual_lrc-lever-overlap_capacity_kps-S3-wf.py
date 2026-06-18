#!/usr/bin/env python3
"""
lrc14_residual_lrc-lever-overlap_capacity_kps-S3-wf

THE COUNTING ARGUMENT, made rigorous-candidate.

Picture established:
  - G_{L'} (cluster minus top) is a COMB of many wide safe arcs (>=24, avg 267), each of
    width > 1/(7 Vmax), spread across [0,1).
  - D_P (small-part danger) has teeth: for each u in P (u<=13), 2u teeth each of width 1/(7u),
    so 2u tooth-intervals; total over P at most sum 2u*(1/(7u)) = 2|P|/7 measure, but the
    NUMBER of teeth is sum_{u in P} 2u <= 2*(sum of P) which is small (speeds <=13).

DESTRUCTION CAPACITY:  a wide L'-arc [a,b] (width w>thresh) is "killed" only if D_P leaves no
length-thresh sub-window safe inside it.  A single P-tooth (one interval of D_P) of width
1/(7u) can, by sitting in the middle, split a wide arc but to FULLY kill an arc of width w you
need D_P to cover all but <thresh on each P-safe piece.  Cleanest rigorous bound:

  An L'-arc SURVIVES (yields a common arc > thresh) UNLESS every length-thresh window inside it
  meets D_P.  Since D_P is a union of teeth, the number of L'-arcs that D_P can possibly
  intersect at all is at most (number of D_P teeth) + (number of L'-arcs whose interior a tooth
  lands in).  A single tooth (interval) intersects at most ... it can lie inside one L'-arc, or
  straddle a boundary touching 2.  So #(L'-arcs touched by D_P) <= 2 * #teeth(D_P).

  #teeth(D_P) = sum_{u in P} 2u  (u teeth in each half? -- count exactly here).

So if N_wide(L') > 2 * #teeth(D_P), then at least one wide L'-arc is UNTOUCHED by D_P, hence
lies entirely in G_P, hence is a common arc of width > thresh => C(S).  PROVED conditionally.

THIS SCRIPT measures, exactly:
  - N_wide = number of L'-safe arcs of width > 1/(7 Vmax)
  - Tcap   = 2 * (number of D_P teeth) = 2 * sum_{u in P} (#teeth of u)
  - whether N_wide > Tcap   (the clean sufficient inequality)
  - AND the refined count: actual #(wide L'-arcs untouched by D_P).
Report worst cases and whether the clean inequality always holds on S3.
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
    """Danger teeth of runner u: {tau: ||u tau|| < h}. Return list of (lo,hi) on [0,1),
       splitting wrap. Each tooth centered at j/u, radius h/u."""
    iv=[]
    for j in range(0,u):
        c=F(j,u); a=(c-h/u)%1; b=(c+h/u)%1
        if a<b: iv.append((a,b))
        else: iv.append((a,F(1))); iv.append((F(0),b))
    return iv

def danger_merged(A,h=C):
    iv=[]
    for u in A: iv+=teeth_intervals(u,h)
    iv.sort(); merged=[]
    for a,b in iv:
        if merged and a<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],b))
        else: merged.append((a,b))
    return merged

def safe_from_danger(merged):
    safe=[]; prev=F(0)
    for a,b in merged:
        if a>prev: safe.append((prev,a))
        prev=max(prev,b)
    if prev<1: safe.append((prev,F(1)))
    return safe

def safe_components(A,h=C):
    return safe_from_danger(danger_merged(A,h))

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

def num_teeth_DP(P,h=C):
    """Number of CONNECTED danger components of D_P (after merging overlaps).
       This is the right 'capacity' since merged teeth count once."""
    return len(danger_merged(P,h))

def count_wide_arcs(arclist,thresh):
    cnt=0
    for a,b in arclist:
        if b-a>thresh: cnt+=1
    # wrap
    if arclist and arclist[0][0]==0 and arclist[-1][1]==1 and len(arclist)>1:
        if arclist[0][1]+(1-arclist[-1][0])>thresh: cnt+=1
    return cnt

def run(n_sets=400, seed=31):
    rng=random.Random(seed)
    flush("="*70)
    flush("CAPACITY COUNTING: N_wide(L') vs destruction capacity of D_P")
    flush("="*70)
    got=0; attempts=0
    clean_ineq_holds=0   # N_wide > 2 * #merged-teeth(D_P)
    clean_ineq_holds_strict=0  # N_wide > #merged-teeth(D_P) (weaker capacity model)
    untouched_min=10**9
    worst=None
    Nwide_min=10**9; cap_max=0
    for_record=[]
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
        Nwide=count_wide_arcs(GLp,thresh)
        ncomp_DP=num_teeth_DP(P) if P else 0   # # merged danger components of P
        cap2=2*ncomp_DP
        Nwide_min=min(Nwide_min,Nwide); cap_max=max(cap_max,cap2)
        if Nwide>cap2: clean_ineq_holds+=1
        if Nwide>ncomp_DP: clean_ineq_holds_strict+=1
        # actual untouched count: wide L'-arcs that contain a length-thresh P-safe window
        DP=danger_merged(P) if P else []
        # for each wide L'-arc, is there a P-safe subwindow of width>thresh inside?
        GP=safe_from_danger(DP) if DP else [(F(0),F(1))]
        # intersect
        untouched=0
        i=j=0
        L1=[(a,b) for a,b in GLp]; L2=GP
        # count wide arcs of intersection
        inter=[]
        ii=jj=0
        while ii<len(L1) and jj<len(L2):
            a,b=L1[ii]; c,d=L2[jj]
            lo=max(a,c); hi=min(b,d)
            if lo<hi: inter.append((lo,hi))
            if b<d: ii+=1
            else: jj+=1
        untouched=count_wide_arcs(inter,thresh)
        untouched_min=min(untouched_min,untouched)
        if worst is None or Nwide-cap2 < worst[0]:
            worst=(Nwide-cap2,Nwide,cap2,ncomp_DP,len(P),sorted(S),float(thresh))
        for_record.append((Nwide,cap2,ncomp_DP,untouched))
    flush(f"\n  S3 sets: {got}")
    flush(f"  CLEAN INEQ  N_wide(L') > 2*#comp(D_P):   {clean_ineq_holds}/{got}")
    flush(f"  WEAKER INEQ N_wide(L') > #comp(D_P):     {clean_ineq_holds_strict}/{got}")
    flush(f"  min N_wide={Nwide_min}, max 2*#comp(D_P)={cap_max}")
    flush(f"  min #(untouched wide L'-arcs surviving): {untouched_min}")
    if worst:
        flush(f"\n  tightest margin (N_wide - 2*#comp): {worst[0]}")
        flush(f"    N_wide={worst[1]} cap2={worst[2]} #comp(D_P)={worst[3]} |P|={worst[4]}")
        flush(f"    S={worst[5]}")
    # distribution
    margins=[r[0]-r[1] for r in for_record]
    flush(f"\n  N_wide - 2*#comp(D_P): min={min(margins)} max={max(margins)} "
          f"avg={sum(margins)/len(margins):.1f}")

if __name__=="__main__":
    run(400)
