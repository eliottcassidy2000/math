#!/usr/bin/env python3
"""
ANGLE H part 4 — LINK 2 end-to-end: the slow-fast / global-witness reduction.

The chain claims:   rho*_{1/7}(P,E) = meas(G_P cap {maxgap{frac(e_i x)}>1/7}) > 0
                   ==>  M(S) >= 1/14  for S = P ∪ {Vmax - e_i}  (Vmax large).

We test this DIRECTLY and ADVERSARIALLY at FINITE Vmax:
  1. Pick (P,E) admissible (|P|+|E|=13, 0 in E).
  2. For a sweep of Vmax, build S = P ∪ {Vmax - e : e in E}.
  3. Compute exact M(S).  Check M(S) >= 1/14.
  4. Independently compute rho*_{1/7}(P,E) and check rho*>0.
  5. The KEY adversarial question: find (P,E,Vmax) with rho*_{1/7}>0 BUT M(S)<1/14
     (would break the SUFFICIENCY direction = the reduction is WRONG), OR with
     rho*=0 but it's actually a covering set that IS lonely (reduction not
     necessary — benign) vs NOT lonely (genuine LRC counterexample - catastrophe).

Also DIRECTLY re-derive the threshold:  a single point tau is safe for runner u
iff ||u tau|| >= 1/14.  In slow-fast: u = Vmax - e, tau = (j + phi/(... )) ...
We avoid the asymptotic and just compute M(S) exactly for the reconstructed S,
then ask whether the reduction's PREDICTION (rho*>0 => M>=1/14) is correct.
"""
from fractions import Fraction as F
from itertools import combinations
import sys
sys.path.insert(0, "04-computation")
from lrc14_adversarial_chain_macmini_0618sH import (
    M_of_S, measure_safe, mu_theta, frac_dist, frac_dist_signed, circular_maxgap)

def rho_star_17(P, E):
    """meas(G_P cap {x: circular maxgap{frac(e x):e in E} > 1/7}), exact."""
    # combine breakpoints of G_P and of the maxgap indicator
    P = sorted(set(P)); E = sorted(set(E))
    half = F(1,14); theta = F(1,7)
    bps = set([F(0), F(1)])
    for p in P:
        for k in range(0, p+1):
            for s in (half, -half):
                t = (F(k)+s)/p
                if F(0)<=t<=F(1): bps.add(t)
    diffs = set()
    for a,b in combinations(E,2):
        d=abs(a-b)
        if d: diffs.add(d)
    for d in diffs:
        for n in range(-1, d+2):
            bps.add(F(n,d))
            for s in (theta,-theta,1-theta,theta-1):
                t=(F(n)+s)/d
                if F(0)<t<F(1): bps.add(t)
    bps = sorted(b for b in bps if F(0)<=b<=F(1))
    total=F(0)
    for i in range(len(bps)-1):
        a,b=bps[i],bps[i+1]
        if b==a: continue
        mid=(a+b)/2
        if all(frac_dist(p*mid)>=half for p in P):
            pts=[frac_dist_signed(e*mid) for e in E]
            if circular_maxgap(pts)>theta:
                total+=(b-a)
    return total

print("="*70)
print("LINK 2: reconstruct S=P∪{Vmax-e}, exact M(S), test rho*>0 ==> M>=1/14")
print("="*70)

# admissible test families (P,E), |P|+|E|=13
tests = [
  # (P, E) with |P|+|E|=13
  ([1,2,3,4,5],            [0,1,2,3,4,5,6,7]),         # k=8 consec cluster
  ([1,5,7,8,9],            [0,1,2,3,4,5,6,7]),         # k=8, the binding P for cap_8
  ([1,2,3,12,13],          [0,2,3,4,5,6,7,8]),         # k=8 perforated (THM-530 row)
  ([1,2,3,13],             [0,1,2,3,4,5,6,7,8]),       # k=9 consec
  ([1,2,3],                [0,1,2,3,4,5,6,7,8,9]),     # k=10 consec
  ([1,2,3,6,12,13],        [0,2,3,4,5,6,8]),           # k=7 the via-max-zero witness
]
thr = lambda k: F(1) - min(measure_safe(list(P)) for P in combinations(range(1,14),13-k))

for P, E in tests:
    k = len(E)
    rs = rho_star_17(P, E)
    mu = mu_theta(E, F(1,7))
    gP = measure_safe(P)
    print(f"\nP={P}  E={E}  (k={k}, |P|={len(P)})")
    print(f"  meas(G_P)={float(gP):.4f}  mu_1/7(E)={float(mu):.4f}  rho*_1/7={float(rs):.4f}")
    # reconstruct S for a sweep of Vmax (must exceed max e so all speeds>13 ideally)
    emax = max(E)
    fails = []
    Mvals = []
    for Vmax in range(max(20, emax+14), max(20,emax+14)+24):
        S = sorted(set(P) | {Vmax - e for e in E})
        if len(S) != 13:      # collisions: skip (not primitive 13-set)
            continue
        if any(v<=0 for v in S): continue
        import math
        g = 0
        for v in S: g = math.gcd(g, v)
        # gcd!=1 -> scale; M is scale-invariant so fine, but skip non-primitive label
        M, arg = M_of_S(S)
        Mvals.append((Vmax, M))
        if M < F(1,14):
            fails.append((Vmax, S, M))
    if Mvals:
        mn = min(m for _,m in Mvals)
        print(f"  Vmax sweep: min M(S) over {len(Mvals)} sets = {mn} = {float(mn):.5f}"
              f"  (1/14={float(F(1,14)):.5f})")
        if fails:
            print(f"  *** {len(fails)} sets with M<1/14! e.g. {fails[0][1]} M={fails[0][2]}")
            # is it covering? does rho*>0 still? -> would be a genuine break
        else:
            print(f"  all M(S) >= 1/14 (reduction prediction holds for these Vmax)")
