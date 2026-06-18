#!/usr/bin/env python3
"""
lrc14_residual_lrc-lever-overlap_decouple_kps-S3-wf

THE DECOUPLING TEST.

We want C(S): W(S\{Vmax}) > 1/(7 Vmax).   A = S\{Vmax} = P u L'.
G_A = G_P cap G_{L'}.

Three candidate LOWER bounds for W(G_A), in increasing strength of decoupling:

  (T0) W(G_A)  itself (the truth).
  (T1) Does the cluster alone already over-deliver:  W(G_{L'}) >> 1/(7 Vmax)?
       If G_{L'} has MANY wide arcs spread around [0,1), at least one must avoid the
       (few, narrow) dangerous regions of P.  -> a covering/pigeonhole on arc count.
  (T2) MEASURE pigeonhole:  meas(G_{L'} cap G_P) >= meas(G_{L'}) + meas(G_P) - 1.
       If that's positive AND G_{L'}'s arcs are individually wide, get a wide common arc.

CRITICAL quantities to measure exactly for each S3 set:
  - meas(G_P), widest arc of G_P, NUMBER of P-danger teeth in [0,1) and their total width
  - meas(G_{L'}), widest arc of G_{L'}, number of L'-safe arcs, their widths
  - the actual W(G_A) and whether T1/T2 would have certified it

Also: the small part P is a covering-ish system with all speeds <= 13 (mostly).
Its danger set D_P = {tau: some u in P has ||u tau|| < 1/14} has measure
   meas(D_P) <= sum_{u in P} 2u * (1/(7u)) = sum_{u in P} 2/7 = 2|P|/7.
That's a HARD a-priori bound independent of the speeds!  With |P|<=11, meas(D_P)<=22/7>1
(useless), but P here is SMALL speeds, |P| typically 4-8, and the teeth OVERLAP a lot.
Measure meas(D_P) exactly and compare to the naive 2|P|/7.
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

def measure_safe(A): return sum(b-a for a,b in safe_components(A))
def widest_arc(arclist):
    if not arclist: return F(0)
    ws=[b-a for a,b in arclist]
    if arclist[0][0]==0 and arclist[-1][1]==1 and len(arclist)>1:
        ws.append(arclist[0][1]+(1-arclist[-1][0]))
    return max(ws)
def Wwidth(A): return widest_arc(safe_components(A))
def intersect_arclists(L1,L2):
    out=[]; i=j=0
    while i<len(L1) and j<len(L2):
        a,b=L1[i]; c,d=L2[j]
        lo=max(a,c); hi=min(b,d)
        if lo<hi: out.append((lo,hi))
        if b<d: i+=1
        else: j+=1
    return out
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

def decompose_PL(S):
    """P = speeds <=13 ; L = speeds >13 (the large part). Then L' = L\\{Vmax}."""
    S=set(S); P=set(v for v in S if v<=13); L=set(v for v in S if v>13)
    return P,L

def run(n_sets=300, seed=11):
    rng=random.Random(seed)
    flush("="*70)
    flush("DECOUPLING TEST: P (<=13) vs L (>13). Can cluster alone certify?")
    flush("="*70)
    got=0; attempts=0
    T1_certifies=0  # W(G_{L'}) > 1/(7Vmax) alone
    T2_certifies=0  # meas(G_{L'})+meas(G_P)-1 > 0 AND a common arc > thresh
    truth_holds=0
    measDP_vs_naive=[]
    WLp_ratios=[]; meas_excess=[]
    worst_truth=None
    while got<n_sets and attempts<n_sets*40:
        attempts+=1
        S=make_S3_set(rng)
        if S is None: continue
        got+=1
        Vmax=max(S); P,L=decompose_PL(S)
        Lp=L-{Vmax}
        thresh=F(1,7*Vmax)
        # truth
        A=P|Lp
        WA=Wwidth(A)
        if WA>thresh: truth_holds+=1
        else:
            if worst_truth is None or WA/thresh < worst_truth[0]:
                worst_truth=(WA/thresh,sorted(S))
        # T1: cluster alone
        WLp=Wwidth(Lp) if Lp else F(1)
        WLp_ratios.append(WLp/thresh)
        if WLp>thresh: T1_certifies+=1
        # measures
        mP=measure_safe(P) if P else F(1)
        mLp=measure_safe(Lp) if Lp else F(1)
        excess=mP+mLp-1
        meas_excess.append(excess)
        if excess>0:
            inter=intersect_arclists(safe_components(P) if P else [(F(0),F(1))],
                                     safe_components(Lp) if Lp else [(F(0),F(1))])
            if widest_arc(inter)>thresh: T2_certifies+=1
        # meas(D_P) vs naive 2|P|/7
        mDP = 1-mP
        naive = F(2*len(P),7)
        measDP_vs_naive.append((mDP,naive,len(P)))
    flush(f"\n  S3 sets: {got}")
    flush(f"  TRUTH  W(G_A)>1/(7Vmax):           {truth_holds}/{got}")
    flush(f"  T1 (cluster L' alone certifies):   {T1_certifies}/{got}")
    flush(f"  T2 (measure-excess + common arc):  {T2_certifies}/{got}")
    if WLp_ratios:
        flush(f"\n  W(G_L')/thresh:  min={float(min(WLp_ratios)):.4f}  max={float(max(WLp_ratios)):.4f}")
    if meas_excess:
        pos=sum(1 for e in meas_excess if e>0)
        flush(f"  meas(G_P)+meas(G_L')-1 > 0 :  {pos}/{got}  "
              f"(min={float(min(meas_excess)):.4f} max={float(max(meas_excess)):.4f})")
    if measDP_vs_naive:
        ratios=[float(a/b) for a,b,_ in measDP_vs_naive if b>0]
        flush(f"  meas(D_P) / (2|P|/7):  min={min(ratios):.4f} max={max(ratios):.4f} "
              f"(how slack the union bound is)")
        mx=max(measDP_vs_naive, key=lambda t:t[0])
        flush(f"  max meas(D_P)={float(mx[0]):.4f} at |P|={mx[2]}  (meas of P-danger set)")
    if worst_truth:
        flush(f"\n  worst truth ratio: {float(worst_truth[0]):.4f}  S={worst_truth[1]}")

if __name__=="__main__":
    run(300)
