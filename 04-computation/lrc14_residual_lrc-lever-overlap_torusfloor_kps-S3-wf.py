#!/usr/bin/env python3
"""
lrc14_residual_lrc-lever-overlap_torusfloor_kps-S3-wf

THE LIMITING TORUS FLOOR — confirm rho*(Delta,P) > 0 provably, and find a UNIFORM positive
lower bound over ALL admissible (Delta, P).

In the limit w0 -> infty, on ruler period I_j the cluster band is governed by tau ~ x/w0 with
x = (j + phi), phi in (1/14,13/14).  As j sweeps, the relevant phases are:
   - cluster: the window where ALL {(w0+d) tau} in [1/14,13/14].  With tau -> 0, (w0+d)tau
     spans an interval of length (s)*tau ~ s*x/w0.  In the SCALED variable, set tau itself the
     parameter; the cluster-safe condition on the *continuum* is just G_{L'} which is exactly
     periodic-ish.  The clean limiting object: the LEFT-ENDPOINT phase of the ruler window
     theta_j = frac(w0 tau_j) ... we instead compute rho* by a DIRECT high-resolution exact
     fraction at very large w0, AND by the analytic "independent-phases" lower bound below.

ANALYTIC LOWER BOUND (the provable floor).  Consider one ruler period I_j (width 6/(7 w0)).
For tau in I_j, set y = w0 tau - j in (1/14, 13/14).  Then:
  - cluster member w0+d:  (w0+d) tau = j + y + d*tau.  Safe iff y + d*tau in [1/14,13/14] mod 1.
    Since d*tau in [0, s*tau], for tau in I_j this is a small shift delta_d = d*tau in
    [0, s*(j+13/14)/w0^2 *w0] ... bounded by s/(w0) * (j/w0) <= s/w0 (since j<w0).  Actually
    delta_d = d*tau <= s * (j+1)/w0 ... grows with j.  Handle via the EXACT computation.
  - small runner u in P:  u*tau safe iff frac(u tau) in [1/14,13/14].  u tau = u(j+y)/w0.

The cluster constraint costs at most  (sum_d shift) but for the WINDOW to be nonempty we need
the common intersection over d of [1/14 - d tau, 13/14 - d tau] (mod the gap) to have length
> thresh.  Worst case shrink = s*tau_max.  Over period I_j, tau ~ (j+1/2)/w0, so the cluster
eats  ~ s*(j+1/2)/w0 * (factor)  -- but the gap has length 6/(7w0)*w0 = 6/7 in y-units...

Rather than a messy hand-bound, we COMPUTE the limiting rho* by exact evaluation at a LADDER
of primes w0 and confirm the floor, scanning the SPACE of (Delta, P) for the MINIMUM rho*.

GOAL: report  inf over tested (Delta,P)  of  lim_{w0} rho_K  -- if it is a positive number
bounded away from 0, the equidistribution proof has a uniform floor.  Also report the SMALLEST
#good ever seen at the largest w0 (must stay >= 1, ideally grows).
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
from itertools import combinations
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

def good_periods(P, Lp, Vmax):
    w0=min(Lp); thresh=F(1,7*Vmax)
    GLp=safe_components(Lp); GP=safe_components(P) if P else [(F(0),F(1))]
    GA=[]; i=j=0
    while i<len(GLp) and j<len(GP):
        a,b=GLp[i]; c,d=GP[j]
        lo=max(a,c); hi=min(b,d)
        if lo<hi: GA.append((lo,hi))
        if b<d: i+=1
        else: j+=1
    periods=set()
    for a,b in GA:
        if b-a<=thresh: continue
        jj=int(w0*a)
        for cj in range(max(0,jj-1),jj+2):
            lo=F(14*cj+1,14*w0); hi=F(14*cj+13,14*w0)
            x=max(a,lo); y=min(b,hi)
            if y-x>thresh: periods.add(cj)
    return len(periods), w0

def rho_limit(P, Delta, w0=3203):
    s=max(Delta)
    Lp=[w0+d for d in Delta]
    Vmax=max(Lp)+14
    cnt,_=good_periods(P,Lp,Vmax)
    return F(cnt,w0), cnt

def run():
    flush("="*70)
    flush("TORUS FLOOR: minimum limiting rho* over a broad scan of (Delta, P)")
    flush("="*70)
    rng=random.Random(303)
    min_rho=None; min_cnt=10**9; worst=None
    worst_cnt=None
    n=0
    for trial in range(250):
        # random Delta (offset pattern), small spread, and random P from 1..13
        ssize=rng.randint(2,8)
        s=rng.randint(1,60)
        Delta=sorted(set([0]+[rng.randint(0,s) for _ in range(ssize-1)]))
        if max(Delta)==0: continue
        psize=rng.randint(3,11)
        P=sorted(rng.sample(range(1,14),psize))
        rho,cnt=rho_limit(P,Delta)
        n+=1
        if min_rho is None or rho<min_rho:
            min_rho=rho; worst=(P,Delta,float(rho),cnt)
        if cnt<min_cnt:
            min_cnt=cnt; worst_cnt=(P,Delta,cnt)
    flush(f"\n  scanned {n} (Delta,P) pairs at w0=3203")
    flush(f"  MIN rho* = {float(min_rho):.5f}   (= #good/w0)")
    flush(f"    at P={worst[0]} Delta={worst[1]} cnt={worst[3]}")
    flush(f"  MIN #good periods = {min_cnt}")
    flush(f"    at P={worst_cnt[0]} Delta={worst_cnt[1]}")
    # The single most adversarial known pattern: P = ALL of 1..13, Delta = {0,1}
    flush("\n  Adversarial extreme: P=range(1,14), Delta={0,1}:")
    for w0 in [101,401,1601,6401,12809]:
        Lp=[w0,w0+1]; Vmax=w0+1+14
        cnt,_=good_periods(list(range(1,14)),Lp,Vmax)
        flush(f"    w0={w0:>6}  #good={cnt:>5}  rho_K={float(F(cnt,w0)):.5f}")

if __name__=="__main__":
    run()
