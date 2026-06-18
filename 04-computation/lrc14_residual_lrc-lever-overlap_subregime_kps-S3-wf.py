#!/usr/bin/env python3
"""
lrc14_residual_lrc-lever-overlap_subregime_kps-S3-wf

CHARACTERIZE the PROVED sub-regime of the Collapse Lemma with a CLEAN inequality, and
quantify the residual honestly.

The single-gap witness J_j = [(j+1/14)/w0, (j+13/14)/wT] (wT=max L') is nonempty with
width > 0 iff  (j+13/14)/wT > (j+1/14)/w0  iff  (j+13/14) w0 > (j+1/14) wT  iff
   13/14 w0 - 1/14 wT  >  j (wT - w0) = j*s.
So small j works; the LARGEST usable j is j_max = floor( (13 w0 - wT)/(14 s) ) (for s>0).
For the window to even exist need 13 w0 > wT, i.e. wT/w0 < 13  (cluster spread ratio < 13).
Since wT = w0 + s, this is s < 12 w0, ALWAYS true for clusters (s<=~80, w0>=14).  Good.

Width of J_0 = 13/(14 wT) - 1/(14 w0) = (13 w0 - wT)/(14 w0 wT).
For this to beat thresh = 1/(7 Vmax) = 2/(14 Vmax):
   (13 w0 - wT)/(14 w0 wT) > 2/(14 Vmax)
   (13 w0 - wT) Vmax > 2 w0 wT.                              (WIDTH-INEQ)

So the FIRST-GAP window J_0 alone beats threshold whenever WIDTH-INEQ holds.  But J_0 sits
near tau ~ 1/(14 w0) (close to 0) where P (small speeds) is DANGEROUS.  We need a J_j inside
a P-safe arc.  There are (j_max+1) candidate windows J_0..J_{j_max} on a ~1/w0 grid.

THE CLEAN THEOREM candidate (single-gap lever):
  (T) If  (#candidate windows J_j of width>thresh)  >  (capacity of D_P to block them),
      then one J_j is P-safe -> C(S).
  Concretely the windows J_j for j=0..j_max are located at left ends (j+1/14)/w0, i.e. an
  arithmetic progression with step ~1/w0 (plus drift).  They sample [0, (j_max)/w0] ~ [0, 13/14].

This script:
  1. Computes j_max and the number of width>thresh windows N_J.
  2. Computes how many J_j are P-safe (the certified count) exactly.
  3. Tests a clean SUFFICIENT inequality:  N_J  >  2*(sum over u in P of u)/1 ... no.
     Better: the windows J_j sit on a grid of step delta_j = (J_{j+1}.left - J_j.left).
     The widest P-safe arc has width WP.  If a run of consecutive J_j of length > 1/WP-ish
     all lie in one P-danger comb... we just COUNT.
  4. Reports: fraction of S3 where >=1 J_j is P-safe (the TRUE single-gap coverage), and the
     residual = sets needing the MULTI-gap (straddling) witness.

ALSO: the honest residual handler. For residual sets, we VERIFY (computationally, exact) that
C(S) still holds via the FULL safe-set computation W(S\{Vmax})>thresh, and report that the
gap is ONLY in the closed-form proof, not in the truth.
"""
from fractions import Fraction as F
from math import gcd, floor
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
def widest_arc(arclist):
    if not arclist: return F(0)
    ws=[b-a for a,b in arclist]
    if arclist[0][0]==0 and arclist[-1][1]==1 and len(arclist)>1:
        ws.append(arclist[0][1]+(1-arclist[-1][0]))
    return max(ws)
def is_covering(S): return all(any(v%q==0 for v in S) for q in range(2,15))
def gcd_all(S): return reduce(gcd,S)

def P_safe_interval(p,q,GP):
    for a,b in GP:
        if a<=p and q<=b: return True
    return False

def make_S3_set(rng, max_attempts=20000, spread_max=45, vlo=20, vhi=400):
    for _ in range(max_attempts):
        small=set(rng.sample(range(1,14), rng.randint(4,9)))
        V0=rng.randint(vlo,vhi); s=rng.randint(8,spread_max)
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

def analyze(S):
    Vmax=max(S)
    P=frozenset(v for v in S if v<=13); L=sorted(v for v in S if v>13)
    Lp=[v for v in L if v!=Vmax]
    w0=min(Lp); wT=max(Lp)
    s=wT-w0
    thresh=F(1,7*Vmax)
    GP=safe_components(P) if P else [(F(0),F(1))]
    # j_max for nonempty J_j width
    if s>0:
        jmax = (13*w0 - wT)//(14*s)  # floor; window nonempty while (j+13/14)w0>(j+1/14)wT
        # recompute exactly: need (14j+13)w0 > (14j+1)wT  <=> 13w0 - wT > 14 j s
        jmax = (13*w0 - wT - 1)//(14*s) if (13*w0-wT)>0 else -1
    else:
        jmax = w0-1
    N_J=0; N_J_wide=0; N_J_psafe=0
    psafe_widths=[]
    for j in range(0, max(0,jmax+1)):
        p=F(14*j+1,14*w0); q=F(14*j+13,14*wT)
        if q<=p or q>=1: continue
        N_J+=1
        if q-p>thresh:
            N_J_wide+=1
            if P_safe_interval(p,q,GP):
                N_J_psafe+=1; psafe_widths.append(q-p)
    WIDTH_INEQ = (13*w0-wT)*Vmax > 2*w0*wT  # J_0 width>thresh
    truth = widest_arc(safe_components(list(P)+Lp))>thresh
    return dict(Vmax=Vmax,w0=w0,wT=wT,s=s,thresh=thresh,jmax=jmax,
                N_J=N_J,N_J_wide=N_J_wide,N_J_psafe=N_J_psafe,
                WIDTH_INEQ=WIDTH_INEQ,truth=truth,
                mDP=1-measure_safe(P) if P else F(0))

def run(n_sets=500, seed=81, spread_max=45, vlo=20, vhi=400, tag=""):
    rng=random.Random(seed)
    flush("="*70)
    flush(f"SINGLE-GAP SUBREGIME  [{tag}]  spread<= {spread_max}, V0 in [{vlo},{vhi}]")
    flush("="*70)
    got=0; attempts=0
    width_ineq_holds=0
    single_gap_certified=0  # >=1 J_j P-safe
    truth_holds=0
    residual=[]
    NJ_when_cert=[]
    while got<n_sets and attempts<n_sets*50:
        attempts+=1
        S=make_S3_set(rng,spread_max=spread_max,vlo=vlo,vhi=vhi)
        if S is None: continue
        got+=1
        d=analyze(S)
        if d['truth']: truth_holds+=1
        if d['WIDTH_INEQ']: width_ineq_holds+=1
        if d['N_J_psafe']>=1:
            single_gap_certified+=1
            NJ_when_cert.append(d['N_J_wide'])
        else:
            residual.append((sorted(S),d))
    flush(f"\n  S3 sets: {got}")
    flush(f"  TRUTH W(G_A)>thresh:                       {truth_holds}/{got}")
    flush(f"  WIDTH-INEQ (J_0 width>thresh):             {width_ineq_holds}/{got}")
    flush(f"  SINGLE-GAP certified (some J_j P-safe):    {single_gap_certified}/{got}  "
          f"({100*single_gap_certified/got:.1f}%)")
    flush(f"  RESIDUAL (needs straddling/multi-gap):     {len(residual)}")
    # verify ALL residual still TRUE
    res_true=sum(1 for _,d in residual if d['truth'])
    flush(f"  of residual, C(S) STILL TRUE (exact full check): {res_true}/{len(residual)}")
    if NJ_when_cert:
        flush(f"  when certified: #wide J_j windows available min={min(NJ_when_cert)} max={max(NJ_when_cert)}")
    return single_gap_certified, got, residual

if __name__=="__main__":
    run(500, seed=81, spread_max=45, tag="GENERAL")
    run(400, seed=82, spread_max=14, tag="TIGHT<=14")
    run(400, seed=83, spread_max=80, tag="LOOSE<=80")
