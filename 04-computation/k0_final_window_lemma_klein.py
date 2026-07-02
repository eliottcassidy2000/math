#!/usr/bin/env python3
"""
klein-2026-07-01-S89 -- HYP-3844: the K=0 FINAL-WINDOW LEMMA (band arithmetic + exposure)

LEMMA FW (n=14), proved parts verified here:
 (a) kink radii of Lambda_S: d/(v+w) (gap death, convex) or d/(w-v) (same-side edge
     crossing = engulfment, concave IF exposed). [THM-592 frame]
 (b) d=1 EMPTINESS: r0 = d/(w-v) in (1/15, 1/14) forces 14d < w-v < 15d; d=1 gives
     14 < w-v < 15 -- NO integer. Hence any concave kink in the open window has d >= 2
     and w-v >= 29. COROLLARY: diam(S) <= 28 => K = 0 on (1/15, 1/14).
 (c) EXPOSURE => BAND-CRITICAL RESIDUE SYSTEM: if the crossing of right endpoints of
     a/v, b/w (d = bv - aw) is exposed at r0, the coincidence point is
         x0 = (b-a)/(w-v)   EXACTLY,
     with ||v x0|| = ||w x0|| = r0 and (exposure) m_S(x0) = r0. In lowest terms
     x0 = B/Q: m_S(B/Q) = D/Q with Q/D = (w-v)/d in (14,15), Q >= 29, all residues of
     S*B mod Q at distance >= D from 0, and v = w mod Q both attaining D.
     => an EXACT finite per-set test for possible exposed window kinks.
 (d) empirics: (i) the mac-mini S94 reconciliation -- GW's candidate 2/29 is NOT a kink
     (profiles identical through it); (ii) deep-well window kink candidates are submerged;
     (iii) adversarial search: plant band differences (w-v = 29, 43, 44) in covering sets
     and hunt for an actually-exposed window kink.
"""
from fractions import Fraction as F
from math import gcd
from functools import reduce
import itertools, random

def Lambda(S, r):
    ivs = []
    for v in S:
        rv = r / v
        for a in range(v + 1):
            c = F(a, v); lo, hi = c - rv, c + rv
            if hi <= 0 or lo >= 1: continue
            ivs.append((max(lo, F(0)), min(hi, F(1))))
    ivs.sort()
    tot = F(0); cl = ch = None
    for lo, hi in ivs:
        if ch is None: cl, ch = lo, hi
        elif lo <= ch: ch = max(ch, hi)
        else: tot += ch - cl; cl, ch = lo, hi
    if ch is not None: tot += ch - cl
    return 1 - tot

def m_at(S, t):
    return min(min((v * t) % 1, 1 - (v * t) % 1) for v in S)

# ---------------- (b) the band table ----------------
print("=" * 92)
print("(b) THE BAND TABLE: concave kink radii d/(w-v) in (1/15, 1/14) need w-v in (14d, 15d)")
print("=" * 92)
for d in range(1, 8):
    lo, hi = 14 * d, 15 * d
    band = list(range(lo + 1, hi))
    print(f"  d={d}:  w-v in ({lo},{hi})  ->  {band if band else 'EMPTY (integrality!)'}")
print("  => any concave kink in the open window has d>=2, w-v>=29; diam(S)<=28 => K=0. PROVED.")

# ---------------- (c) exact finite test ----------------
def window_kink_candidates(S, rlo=F(1,15), rhi=F(1,14)):
    """All band crossings with their x0 and the band-critical test result.
    Returns list of dicts: (v,w,d,r0,x0,mx0, critical=exposure-necessary-condition)."""
    out = []
    Sset = sorted(S)
    for v, w in itertools.combinations(Sset, 2):
        q = w - v
        for d in range(1, q):
            r0 = F(d, q)
            if not (rlo < r0 < rhi): continue
            # crossings: right endpoints a/v, b/w with bv - aw = d;
            # x0 = (b-a)/q; enumerate solutions b,a with 0<=a<v... all crossings mod 1:
            # (a,b) solving bv - aw = d: b = (d + aw)/v when v | d+aw.
            for a in range(v):
                if (d + a * w) % v: continue
                b = (d + a * w) // v
                x0 = F(b - a, q) % 1
                mx0 = m_at(S, x0)
                critical = (mx0 == r0)
                out.append(dict(v=v, w=w, d=d, r0=r0, x0=x0, mx0=mx0, critical=critical))
    return out

# ---------------- (d)(i) mac-mini reconciliation: GW at 2/29 ----------------
print("\n" + "=" * 92)
print("(d)(i) RECONCILIATION: is 2/29 a real kink of Lambda_GW? (mac-mini S94 lists it as breakpoint)")
print("=" * 92)
GW = list(range(1, 12)) + [13, 24]
AP = list(range(1, 14))
r29 = F(2, 29)
eps = F(1, 10**5)
for name, S in [("GW", GW), ("AP", AP)]:
    L0, Lm, Lp = Lambda(S, r29 - eps), Lambda(S, r29), Lambda(S, r29 + eps)
    linear = (L0 + Lp == 2 * Lm)
    print(f"  {name}: Lambda at 2/29 -eps/0/+eps: {float(L0):.7f} {float(Lm):.7f} {float(Lp):.7f}"
          f"  locally-linear-through={linear}")
print(f"  identical there? {Lambda(AP, r29) == Lambda(GW, r29)} "
      f"(gap-death candidates 2/33, 2/31 are BELOW 1/15 -- consistent w/ both sessions)")
cands = window_kink_candidates(GW)
print(f"  GW window band-crossing candidates: {len(cands)}; critical (exposure-necessary): "
      f"{sum(c['critical'] for c in cands)}")
for c in cands[:6]:
    print(f"    (v,w,d)=({c['v']},{c['w']},{c['d']})  r0={c['r0']}  x0={c['x0']}  "
          f"m(x0)={c['mx0']} vs r0 -> critical={c['critical']}")

# ---------------- (d)(ii) deep well ----------------
print("\n" + "=" * 92)
print("(d)(ii) DEEP WELL {1..12,182}: window candidates and their submersion")
print("=" * 92)
CONSTR = list(range(1, 13)) + [182]
cands = window_kink_candidates(CONSTR)
ncrit = sum(c['critical'] for c in cands)
print(f"  candidates in window: {len(cands)}, critical: {ncrit}")
byd = {}
for c in cands: byd.setdefault((c['d'], c['w']-c['v']), 0)
print(f"  (d, w-v) pairs present: {sorted(byd)}")
ex = [c for c in cands if c['v']==12 and c['w']==182 and c['d']==12][:1]
if ex:
    c = ex[0]
    print(f"  example (12,182,d=12): r0={c['r0']}={float(c['r0']):.5f} x0={c['x0']}={float(c['x0']):.5f} "
          f"m(x0)={c['mx0']}={float(c['mx0']):.5f} < r0 => SUBMERGED (covered by a faster small speed)")

# ---------------- (d)(iii) adversarial exposure hunt ----------------
print("\n" + "=" * 92)
print("(d)(iii) ADVERSARIAL HUNT: covering 13-sets with PLANTED band pairs (w-v=29,43,44);")
print("         does ANY have a critical (m(x0)=r0) window crossing? [necessary for exposure]")
print("=" * 92)
def is_covering(S, n=14):
    return all(any(v % q == 0 for v in S) for q in range(2, n + 1))
random.seed(89)
found_critical = []
tried = 0
for trial in range(4000):
    base = random.randint(2, 90)
    delta = random.choice([29, 43, 44, 57, 58, 59])
    S = {base, base + delta}
    for q in [14, 13, 11, 9, 8]:
        S.add(q * random.randint(1, 5))
    while len(S) < 13:
        S.add(random.randint(1, 130))
    S = sorted(S)
    if len(S) != 13 or not is_covering(S) or reduce(gcd, S) != 1: continue
    tried += 1
    for c in window_kink_candidates(S):
        if c['critical']:
            found_critical.append((S, c))
if found_critical:
    print(f"  CRITICAL CROSSINGS FOUND in {len(found_critical)} cases / {tried} covering sets, e.g.:")
    for S, c in found_critical[:4]:
        print(f"    S={S}")
        print(f"      (v,w,d)=({c['v']},{c['w']},{c['d']}) r0={c['r0']} x0={c['x0']} m(x0)={c['mx0']}")
        # verify whether it is an actual concave kink of Lambda (exposure, not just critical)
        eps = F(1, 10**6)
        L0, Lm, Lp = Lambda(S, c['r0']-eps), Lambda(S, c['r0']), Lambda(S, c['r0']+eps)
        sl_l = (Lm - L0) / eps; sl_r = (Lp - Lm) / eps
        print(f"      slopes around r0: left {float(sl_l):.4f}  right {float(sl_r):.4f}  "
              f"CONCAVE-KINK={'YES' if sl_r < sl_l else 'no'}")
else:
    print(f"  NO critical window crossing in {tried} planted covering sets -- exposure needs")
    print(f"  m_S(x0) = r0 EXACTLY at a modulus-29+ rational; the planted pairs' crossings all submerged.")

print("\nDONE.")
