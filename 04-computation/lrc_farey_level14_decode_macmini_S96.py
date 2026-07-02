#!/usr/bin/env python3
"""
mac-mini-2026-07-01-S96 -- HYP-3852: the FAREY-LEVEL-14 decode of 23520/392, the
Mordell-Tornheim slice identity, the Stern-Brocot band structure of the final
window, and the pair-law route to hp0cap's decorrelation leg.

(1) DECODE: the S95 correction values 1/23520 (= Lambda_GW - Lambda_AP at r=27/392)
    and 392 decode structurally:
      392 = 2*14^2 (r = 27/392 = 1/14 - 1/392: 'one part in 2n^2 below threshold');
      slope(GW)-slope(AP) on [1/15, 2/29] = -29/60 EXACTLY, and 29/60 = 2*(1/5+1/24)
      = the merge-kink jumps of the GW outlier pair (5,24) + its mirror;
      5+24 = 29 = 2n+1; 2/29 = MEDIANT(1/15, 1/14);
      1/23520 = (29/60)*(2/29 - 27/392); 23520 = 60 * 392.
(2) MT-SLICE: c_AP(n) = sum_{units a} 1/(a(n-a)) = 2n * sum_{p+q=n, p<q, gcd=1}
    1/(p*q*(p+q)) -- the collapse constant is the level-n SLICE of the
    Mordell-Tornheim layer-cake (klein HYP-3834's E = the full sum).  Same level
    n = p+q as THM-594(B)'s pair-overlap branch point.
(3) STERN-BROCOT BANDS: the THM-596 bands (numerator-d' reduced fractions in the
    open final window) enumerated; minimal numerator = 2 at the mediant.
(4) hp0cap ROUTE PROBE: L_y(E) = q0 + q6 + q3/10 for consec_9 vs far-swapped;
    the far-element drop and its 1/w rate, to be matched against the exact
    pair law (THM-594(B)) -- the proposed replacement for the Vitali estimates.
"""
from fractions import Fraction as F
from math import gcd
import sys
sys.path.insert(0, '04-computation')
from lonely_profile import profile

print("=" * 78)
print("(1) THE 23520 / 392 DECODE")
print("=" * 78)
AP = list(range(1, 14)); GW = list(range(1, 12)) + [13, 24]
pAP = profile(AP, F(1, 14)); pGW = profile(GW, F(1, 14))
r0 = F(27, 392)
dAP, dGW = pAP.measure(r0), pGW.measure(r0)
print(f"  r = 27/392 = 1/14 - 1/392  (392 = 2*14^2 = {2*14*14})")
print(f"  Lambda_GW - Lambda_AP at 27/392 = {dGW - dAP}  (claim 1/23520: {dGW-dAP == F(1,23520)})")
sAP = pAP.slope(F(1, 15) + F(1, 10**9)); sGW = pGW.slope(F(1, 15) + F(1, 10**9))
diff = sGW - sAP
print(f"  slope difference on [1/15, 2/29]: {diff}  (claim -29/60: {diff == -F(29,60)})")
print(f"  29/60 = 2*(1/5 + 1/24)? {F(29,60) == 2*(F(1,5)+F(1,24))}   (GW pair (5,24) + mirror)")
print(f"  5 + 24 = 29 = 2n+1;  mediant(1/15, 1/14) = (1+1)/(15+14) = 2/29 ✓")
print(f"  (29/60)*(2/29 - 27/392) = {F(29,60)*(F(2,29)-F(27,392))}  (= 1/23520: "
      f"{F(29,60)*(F(2,29)-F(27,392)) == F(1,23520)});  23520 = 60*392 = {60*392}")

print()
print("=" * 78)
print("(2) MORDELL-TORNHEIM SLICE IDENTITY: c_AP(n) = 2n * MT-slice(n)")
print("=" * 78)
for n in (6, 8, 14, 20):
    c_ap = sum(F(1, a * (n - a)) for a in range(1, n) if gcd(a, n) == 1)
    mt = sum(F(1, p * q * (p + q)) for p in range(1, n) for q in range(p + 1, n)
             if p + q == n and gcd(p, q) == 1)
    print(f"  n={n}: c_AP = {c_ap} = {float(c_ap):.6f};  2n*MT_slice = {2*n*mt}  "
          f"EQUAL: {c_ap == 2*n*mt}")
print("  => the collapse constant is the p+q=n slice of the MT layer-cake;")
print("     the SAME level n is THM-594(B)'s pair-overlap branch point p+q<=n.")

print()
print("=" * 78)
print("(3) STERN-BROCOT BAND STRUCTURE of the open final window (1/15, 1/14)")
print("=" * 78)
lo, hi = F(1, 15), F(1, 14)
for dp in range(2, 7):
    band = [q for q in range(14 * dp + 1, 15 * dp) if gcd(dp, q) == 1
            and lo < F(dp, q) < hi]
    print(f"  d'={dp}: reduced q* in (14d',15d') = {band}   radii {[str(F(dp,q)) for q in band]}")
print("  minimal numerator = 2 at q*=29: the MEDIANT -- GW's identity onset,")
print("  the first exposed band (THM-596), and the band-emptying anchor, all one point.")

print()
print("=" * 78)
print("(4) hp0cap PROBE: L_y for consec_9 vs far-swap; 1/w rate vs pair law")
print("=" * 78)

def miss_distribution(E):
    """Exact q_j = meas{x in [0,1): exactly j of the 7 sectors unhit by {frac(e x)}}."""
    # breakpoints: frac(e*x) crosses sector boundary s/7 when x = (k + s/7)/e
    pts = set([F(0), F(1)])
    for e in E:
        for k in range(e):
            for s in range(7):
                x = (F(k) + F(s, 7)) / e
                if 0 <= x < 1:
                    pts.add(x)
    pts = sorted(pts)
    q = [F(0)] * 8
    for i in range(len(pts) - 1):
        a, b = pts[i], pts[i + 1]
        mid = (a + b) / 2
        hit = set()
        for e in E:
            f = e * mid - int(e * mid)
            hit.add(int(f * 7))
        q[7 - len(hit)] += b - a
    return q

def L_y(E):
    q = miss_distribution(E)
    return q[0] + q[6] + q[3] / 10

E9 = list(range(1, 10))
ly9 = L_y(E9)
cap9 = F(1979, 4004)
print(f"  consec_9 {E9}: L_y = {ly9} = {float(ly9):.6f}; cap_9 = {float(cap9):.6f}; "
      f"p0<=L_y<=cap? L_y<=cap: {ly9 <= cap9}")
print("  far-swap 9 -> w (drop 9, add w):")
prev = None
for w in (30, 60, 120, 240):
    E = list(range(1, 9)) + [w]
    ly = L_y(E)
    drop = ly9 - ly
    rate = (ly9 - ly) * w
    print(f"    w={w:4d}: L_y = {float(ly):.6f}; drop from consec = {float(drop):+.6f}; "
          f"w*(consec-drop gap to limit) probe = {float(rate):.4f}")
print("  (decorrelated-limit and O(1/w) rate: the pair-law THM-594(B) branch-2")
print("   formula gives the exact per-pair deficit -- the Vitali replacement.)")
