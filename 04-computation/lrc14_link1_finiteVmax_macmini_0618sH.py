#!/usr/bin/env python3
"""
ANGLE H part 5 — LINK 1: the finite-Vmax control + the (un-enumerated) small-Vmax check.

THM-527-A: rho_K := #{good Vmax-periods j}/Vmax  ->  rho*(P,E) as Vmin->inf,
with rho_K = rho* + O(#arcs/Vmax).  For a witness we need rho_K > 0, i.e.
at least ONE good period j.  Claim: Vmax > C/rho* gives rho_K>0, and Vmax<=V0
is "a finite check".  THE GAP: that finite check is asserted, never enumerated
for k>=3.  Also: is the O(1/Vmax) error rigorously a LOWER bound (so #good>0
is guaranteed once Vmax large), or could ALL good periods vanish at some finite
Vmax even though rho*>0?

We TEST the actual good-period count directly:
  period j is GOOD iff its safe gap I_j contains a point safe for ALL of S.
  Equivalently (THM-527-A): x_j := frac of the slow coordinate at period j,
  and we need x_j in G_P and maxgap{frac(e_i x_j)}>1/7.  BUT this is the
  ASYMPTOTIC surrogate.  The HONEST test: directly count, for the reconstructed
  S at finite Vmax, how many Vmax-safe periods I_j actually contain an S-safe
  point.  If that count is ALWAYS >=1 (so M(S)>=1/14 via criterion C at v=Vmax),
  the reduction delivers.  We look for the SMALLEST good-count and whether it
  ever hits 0 at finite Vmax with rho*>0.
"""
from fractions import Fraction as F
from itertools import combinations
import sys, math
sys.path.insert(0, "04-computation")
from lrc14_adversarial_chain_macmini_0618sH import (
    M_of_S, measure_safe, frac_dist, frac_dist_signed, circular_maxgap)

def safe_for_all(tau, S, half=F(1,14)):
    return all(frac_dist(v*tau) >= half for v in S)

def count_good_periods(S, Vmax):
    """For each Vmax-safe period I_j = ((14j+1)/(14Vmax),(14j+13)/(14Vmax)),
    j=0..Vmax-1, decide if it contains a point safe for ALL of S.  Exact:
    inside I_j the only runners that vary are the others; we sweep breakpoints
    of ||u tau|| for u in S within I_j and test midpoints.  Return #good."""
    good = 0
    others = [u for u in S if u != Vmax]
    half = F(1,14)
    for j in range(Vmax):
        a = F(14*j+1, 14*Vmax); b = F(14*j+13, 14*Vmax)
        # breakpoints of ||u tau||=1/14 for u in others, within (a,b)
        bps = {a, b}
        for u in others:
            # u tau = k +- 1/14 -> tau=(k +-1/14)/u
            klo = int(a*u) - 1; khi = int(b*u) + 1
            for k in range(klo, khi+1):
                for s in (half, -half):
                    t = (F(k)+s)/u
                    if a < t < b: bps.add(t)
        bl = sorted(bps)
        ok = False
        for i in range(len(bl)-1):
            mid = (bl[i]+bl[i+1])/2
            if safe_for_all(mid, S):
                ok = True; break
        if ok: good += 1
    return good

# rho* via the slow-fast surrogate (reuse link2's exact computation)
sys.path.insert(0, "04-computation")
from lrc14_link2_globalwitness_macmini_0618sH import rho_star_17

print("="*70)
print("LINK 1: actual good-period COUNT vs Vmax (does it ever hit 0?)")
print("="*70)

tests = [
  ([1,5,7,8,9],            [0,1,2,3,4,5,6,7]),   # k=8, binding P (smallest meas G_P)
  ([1,2,3],                [0,1,2,3,4,5,6,7,8,9]), # k=10 consec
  ([1,2,3,12,13],          [0,2,3,4,5,6,7,8]),   # k=8 perforated (mu_1/7=1)
]
for P, E in tests:
    k = len(E); emax = max(E)
    rs = rho_star_17(P, E)
    print(f"\nP={P} E={E} (k={k}): rho*_1/7 = {float(rs):.4f}")
    print(f"  {'Vmax':>6} {'#good':>6} {'rho_K':>9} {'M(S)':>10}  primitive?")
    mincount = None
    for Vmax in [40, 60, 80, 120, 200, 300, 500]:
        S = sorted(set(P) | {Vmax - e for e in E})
        if len(S) != 13:
            print(f"  {Vmax:>6}  (collision, skip)"); continue
        g = 0
        for v in S: g = math.gcd(g, v)
        gc = count_good_periods(S, Vmax)
        rho_K = F(gc, Vmax)
        M, _ = M_of_S(S)
        prim = "yes" if g==1 else f"gcd={g}"
        flag = "  <<< ZERO GOOD!" if gc==0 else ""
        print(f"  {Vmax:>6} {gc:>6} {float(rho_K):>9.4f} {float(M):>10.5f}  {prim}{flag}")
        mincount = gc if mincount is None else min(mincount, gc)
    print(f"  min #good over sweep = {mincount}  (need >=1 for criterion-C witness)")
