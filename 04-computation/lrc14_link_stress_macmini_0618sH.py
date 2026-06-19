#!/usr/bin/env python3
"""
ANGLE H part 2 — STRESS THE INDIVIDUAL LINKS.

LINK 2: the global-witness reduction.  Claim (THM-530/531):
   In the slow-fast picture (large Vmax), there is a point safe for ALL of S
   iff x in G_P AND the cluster co-offset phases {frac(e_i x)} leave a circular
   gap > 1/7.  We re-derive the 1/7 and check the half-width arithmetic.

   The teeth: each cluster runner u_i = Vmax - e_i.  Near a Vmax-safe gap, in
   the fast phase phi=frac(Vmax tau), runner u_i is safe iff phi avoids a tooth
   of HALF-WIDTH 1/14 centered at frac(e_i x) (x = slow coordinate).  A single
   safe phi exists iff the union of teeth (each full width 2*(1/14)=1/7) does
   NOT cover the circle, i.e. the centers {frac(e_i x)} leave a circular gap
   > 1/7 (so a phi sits >1/14 from both neighboring centers).  => threshold 1/7.

   ADVERSARIAL TEST: build ACTUAL covering 13-sets S = P ∪ {Vmax - e_i}, compute
   the TRUE M(S) exactly, and check the implication:
     (x in G_P at the relevant slow value AND maxgap>1/7)  ==>  M(S) >= 1/14.
   More importantly check the CONVERSE risk: is it possible that mu_{1/7}(E)
   large but M(S) < 1/14 for a reconstructed S?  (would break the reduction).

LINK 3: N(E) subset S7(E) exact?  maxgap <= 1/7  ==> every sector hit?
   A 1/7-net (maxgap<=1/7) — does it necessarily hit EVERY closed-open sector
   [j/7,(j+1)/7)?  Boundary subtlety: a point exactly at a 1/7 boundary.

LINK "covering used?": is the bound claimed for ALL integer E?  Find an
   integer E with mu_{1/7}(E) < thr_k, or meas(S7(E)) > cap_k.  If none exists
   the claim is genuinely universal (covering not needed for the crux); if one
   exists, covering/primitivity MUST be invoked or the reduction has a gap.
"""
from fractions import Fraction as F
from itertools import combinations, product
import random, sys
sys.path.insert(0, "04-computation")
from lrc14_adversarial_chain_macmini_0618sH import (
    M_of_S, measure_safe, mu_theta, meas_S7, M7, frac_dist, frac_dist_signed,
    circular_maxgap)

print("="*70)
print("LINK 3: N(E) subset S7(E).  Does maxgap<=1/7 force every sector hit?")
print("="*70)
# A 1/7-net: points with circular maxgap <= 1/7.  Claim: hits every sector.
# Proof check: sector S_j = [j/7,(j+1)/7) has length 1/7. If a point set has
# maxgap <= 1/7 then between any two cyclically adjacent points the gap is
# <= 1/7.  Could a sector of length 1/7 contain NO point yet maxgap<=1/7?
# If sector (length exactly 1/7) is empty, the two points flanking it are
# the nearest on each side; the gap spanning the empty sector is >= 1/7, and
# = 1/7 only if both flanking points sit exactly ON the sector boundaries
# j/7 and (j+1)/7.  Then point at (j+1)/7 is in sector j+1 (half-open), point
# at j/7 is in sector j... wait j/7 is the LEFT boundary of S_j => in S_j.
# So if a point sits exactly at j/7 it's IN S_j.  The only way S_j empty with
# maxgap=1/7 is points exactly at j/7 boundaries but then j/7 in S_j.
# CONCLUSION: maxgap<=1/7 (closed) ==> every sector hit, EXCEPT possibly a
# measure-zero set of x.  Verify empirically with exact x near boundaries.
def check_link3(E, samples=2000):
    # sample many x; whenever maxgap<=1/7 verify S7 holds (all sectors hit)
    bad = 0
    for _ in range(samples):
        x = F(random.randint(1, 10**6), 10**6)
        pts = sorted(frac_dist_signed(e*x) for e in E)
        mg = circular_maxgap(pts)
        if mg <= F(1,7):
            hit = set()
            for e in E:
                fx = frac_dist_signed(e*x)
                sec = int(fx*7)
                if sec == 7: sec = 6
                hit.add(sec)
            if len(hit) != 7:
                bad += 1
    return bad

for k in [8, 10, 13]:
    E = list(range(k))
    b = check_link3(E)
    print(f"  consec k={k}: N subset S7 violations (random exact x) = {b}")
# also exact-boundary x where a point sits at j/7
print("  exact boundary x test: x=1/7 with E={0..7}:")
x = F(1,7); E = list(range(8))
pts = sorted(frac_dist_signed(e*x) for e in E)
print(f"    points = {[str(p) for p in pts]}")
print(f"    maxgap = {circular_maxgap(pts)} (vs 1/7={F(1,7)})")
hit = set()
for e in E:
    fx = frac_dist_signed(e*x); sec=int(fx*7); sec=6 if sec==7 else sec; hit.add(sec)
print(f"    sectors hit = {sorted(hit)}  -> {'ALL 7' if len(hit)==7 else 'MISSING'}")

print()
print("="*70)
print("LINK 'covering used?': is meas(S7(E))<=cap_k claimed for ALL integer E?")
print("Hunt for an integer E (0 in E, |E|=k) with meas(S7(E)) > cap_k.")
print("="*70)
cap = {8:F(2243,5880),9:F(1979,4004),10:F(55,91),11:F(66,91),12:F(6,7),13:F(1)}
# NOTE cap_k = min_{|P|=13-k} meas(G_P): for k, |P|=13-k.
# k=8 -> |P|=5 -> 2243/5880; k=9->|P|=4->1979/4004; k=10->|P|=3->55/91;
# k=11->|P|=2->66/91; k=12->|P|=1->6/7; k=13->|P|=0->1.
worst = {}
for k in [8,9,10,11,12,13]:
    capk = cap[k]
    consec = meas_S7(list(range(k)))
    mx = consec; argmx = tuple(range(k))
    # exhaustive over k-subsets of {0..B} containing 0, small box
    B = {8:14,9:13,10:12,11:13,12:13,13:14}[k]
    cnt=0; exceed=0
    for rest in combinations(range(1, B+1), k-1):
        E = (0,)+rest
        v = meas_S7(list(E))
        cnt+=1
        if v > mx:
            mx = v; argmx = E
        if v > capk:
            exceed += 1
    worst[k]=(consec, mx, argmx, capk, exceed)
    print(f"  k={k}: cap={float(capk):.4f}  meas_S7(consec)={float(consec):.4f}"
          f"  MAX_S7={float(mx):.4f} at {argmx}  (#>cap={exceed}/{cnt}, box<= {B})")
