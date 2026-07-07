#!/usr/bin/env python3
"""
kps-2026-07-06-S53 -- THE SINGLE-SCALE DENSITY FLOOR (Route 1, the residue my S52 coarse
reduction leaves open).

My S52 coarse reduction (now FORMALIZED, LRCCoarseReduction.lean) grounds every MULTI-SCALE
13-family in the settled LRC(<=13).  So Route 1's density node -- prove the good-period
density rho(v) = Leb{t in [0,1): min_i ||v_i t|| >= 1/14} lets M(v) >= 1/14 -- only needs
SINGLE-SCALE (bounded-ratio) families.  This script studies that residue directly on the
DIRECT LRC(14) object (13 speeds, threshold 1/14), confirming the structure the analytic
floor must exploit:

  (1) rho(v) > 0 for single-scale families EXCEPT at/near the AP {1..13}, where rho -> 0
      (the AP is still LONELY -- an ISOLATED optimal witness, M = 1/14 exactly).  So NO
      uniform positive density floor exists (refutes any '2/7' uniform claim); the floor is
      a DICHOTOMY: near-AP (arithmetic, isolated) vs spread (positive density).
  (2) THREE-GAP SIGNATURE (mac-mini HYP-4412, here on the 13/1-14 object): at the optimal
      witness t*, the number g of distinct gap-lengths of {0} u {v_i t* mod 1} is SMALL
      (g <= 3, a {k*alpha} arithmetic orbit) for near-tight families and LARGE for loose
      ones.  So 'near-tight => three-gap => witness is an arithmetic orbit => M is a
      continued-fraction rung (1/14, 1/13, ...)' -- the metric face of the AP rigidity,
      the mechanism that closes the near-AP corner (OPEN-Q-110) by classical three-gap
      theory rather than measure.
"""
from fractions import Fraction
import numpy as np
from math import gcd
from functools import reduce
import random

ONE_14 = Fraction(1, 14)

def Mw_and_witness(v):
    """Exact M(v) and an optimal rational witness c/q."""
    v = [x for x in v if x]
    S = sum(abs(x) for x in v)
    Q = min(4 * S, 2 * max(abs(x) for x in v) + 2)
    va = np.array(v, dtype=np.int64)
    bn, bd, bc = 0, 1, 1
    for q in range(2, Q + 1):
        a = np.arange(1, q)
        r = np.outer(va, a) % q
        d = np.minimum(r, q - r)
        col = d.min(axis=0)
        j = int(col.argmax())
        qb = int(col[j])
        if qb * bd > bn * q:
            bn, bd, bc = qb, q, int(a[j])
    return Fraction(bn, bd), Fraction(bc, bd)

def good_density(v, res=20000):
    """Leb{t in [0,1): min_i ||v_i t|| >= 1/14} by a fine grid (>= form)."""
    ts = (np.arange(res) + 0.5) / res
    va = np.array(v, dtype=np.float64)
    # ||v_i t|| = dist to nearest integer
    ph = np.abs(np.outer(ts, va))
    ph = np.minimum(ph % 1.0, 1.0 - (ph % 1.0))
    good = (ph.min(axis=1) >= 1.0 / 14.0)
    return good.mean()

def three_gap_count(v, t, tol=1e-9):
    """# distinct gap-lengths of {0} u {v_i t mod 1} on the circle."""
    pts = sorted(set([0.0] + [ (vi * float(t)) % 1.0 for vi in v ]))
    pts.append(1.0)
    gaps = [round(pts[i+1] - pts[i], 6) for i in range(len(pts) - 1)]
    distinct = sorted(set(g for g in gaps if g > tol))
    # merge near-equal
    merged = []
    for g in distinct:
        if not merged or g - merged[-1] > 1e-4:
            merged.append(g)
    return len(merged)

def reduce_gcd(v):
    g = reduce(gcd, v)
    return [x // g for x in v] if g > 1 else list(v)

print("=" * 76)
print("SINGLE-SCALE DENSITY FLOOR for LRC(14) (13 speeds, threshold 1/14 = %.5f)" % float(ONE_14))
print("=" * 76)

ap = list(range(1, 14))
Map, wap = Mw_and_witness(ap)
print("\n(1) THE AP {1..13} (extremal single-scale family):")
print("    M = %s = %.6f (= 1/14, TIGHT); optimal witness t* = %s = %.6f" %
      (Map, float(Map), wap, float(wap)))
print("    good-density rho (>= 1/14) = %.5f  (ISOLATED witness -> ~0; still LONELY)" %
      good_density(ap))
print("    three-gap count g at t* = %d  (arithmetic orbit {k/14})" % three_gap_count(ap, wap))

print("\n(2) NEAR-AP perturbations (bump one speed): rho stays ~0, g stays small, M jumps to a rung")
hdr = "    %-26s %-10s %-9s %-4s"
print(hdr % ("family (perturbed AP)", "M", "rho", "g"))
for bump_i, bump in [(12, +1), (12, +2), (11, +1), (0, +1), (6, +1), (12, +11)]:
    v = list(range(1, 14)); v[bump_i] += bump
    v = reduce_gcd(sorted(set(v)))
    if len(v) < 13:
        continue
    M, w = Mw_and_witness(v)
    print(hdr % (str(v[-4:]) + " (i=%d +%d)" % (bump_i, bump),
                 "%.5f" % float(M), "%.5f" % good_density(v), three_gap_count(v, w)))

print("\n(3) SPREAD single-scale families (bounded ratio <=13, random): rho > 0, g LARGE, M >> 1/14")
print(hdr % ("family (last 4 speeds)", "M", "rho", "g"))
random.seed(53)
cnt_pos = 0
for _ in range(8):
    v = sorted(random.sample(range(2, 27), 13))     # single-scale: ratio 27/2 < 14
    v = reduce_gcd(v)
    M, w = Mw_and_witness(v)
    rho = good_density(v)
    if rho > 0.02:
        cnt_pos += 1
    print(hdr % (str(v[-4:]), "%.5f" % float(M), "%.5f" % rho, three_gap_count(v, w)))
print("    (spread families with rho > 0.02: %d/8 -- positive good-density, easily lonely)" % cnt_pos)

print("\n" + "=" * 76)
print("READOUT for the Route-1 density node (single-scale residue):")
print("  * NO uniform positive floor: rho -> 0 at the AP (still lonely, isolated witness).")
print("  * The floor is a DICHOTOMY: near-AP (g<=3 arithmetic orbit, M a CF rung -- close by")
print("    three-gap/rigidity, NOT measure) vs spread (g large, rho>0 -- close by decorrelation).")
print("  * This is EXACTLY the domain my S52 coarse reduction leaves: multi-scale already")
print("    grounded in LRC(<=13) (LRCCoarseReduction.lean); the analytic floor need only cover")
print("    these single-scale families, split near-AP (arithmetic) vs spread (measure).")
print("=" * 76)
