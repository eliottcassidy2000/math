#!/usr/bin/env python3
"""
kps-2026-07-06-S54 -- CENSUS of the bottom of the single-scale 13-family M-spectrum, and
verification of opus-S130's mu_{1/7} density floor.

Route 1 (bound M(v) >= 1/14 directly).  My S52/S53 coarse reduction grounds MULTI-SCALE
families in LRC(<=13); the residue is SINGLE-SCALE (bounded-ratio) families, where the
density floor is a dichotomy: near-AP (rigidity) + spread (decorrelation).  This script
does the owner's CENSUS on the near-AP corner:

  PART 1.  Enumerate all primitive 13-subsets of {1..max} (single-scale, ratio <= max),
  compute exact M, and bucket the families with M <= threshold.  Question: what is the
  BOTTOM of the spectrum -- is there a gap above 1/14 (an LRC14-analogue of (G))?  What
  parametric families populate the low buckets?

  PART 2.  Verify opus-S130's floor mu_{1/7}(E) = meas{x: circular max-gap {frac(e x)} > 1/7},
  reproduce mu_{1/7}(AP_13) = 477/1078 ~ 0.4425, and adversarially test AP-minimality.
"""
from fractions import Fraction
import numpy as np
from math import gcd
from functools import reduce
from itertools import combinations
import random

def Mw(v):
    v = [x for x in v if x]
    S = sum(abs(x) for x in v)
    Q = min(4 * S, 2 * max(abs(x) for x in v) + 2)
    va = np.array(v, dtype=np.int64)
    bn, bd = 0, 1
    for q in range(2, Q + 1):
        a = np.arange(1, q)
        r = np.outer(va, a) % q
        d = np.minimum(r, q - r)
        col = d.min(axis=0)
        qb = int(col.max())
        if qb * bd > bn * q:
            bn, bd = qb, q
    return Fraction(bn, bd)

def primitive(v):
    g = reduce(gcd, v)
    return tuple(x // g for x in v) if g > 1 else tuple(v)

print("=" * 78)
print("PART 1 -- CENSUS: bottom of the single-scale 13-family M-spectrum")
print("  (primitive 13-subsets of {1..max}; threshold M <= 0.085; 1/14=%.5f 2/27=%.5f 1/13=%.5f)"
      % (1/14, 2/27, 1/13))
print("=" * 78)

THRESH = Fraction(85, 1000)
from collections import defaultdict
buckets = defaultdict(list)     # M -> list of families
total = 0
for mx in range(13, 21):
    # 13-subsets of {1..mx} with max element = mx (so we don't recount smaller-max ones)
    for combo in combinations(range(1, mx), 12):
        v = combo + (mx,)
        if reduce(gcd, v) != 1:
            continue
        total += 1
        M = Mw(v)
        if M <= THRESH:
            buckets[M].append(v)
print("enumerated %d primitive 13-families (max<=20); %d with M<=%.3f" %
      (total, sum(len(x) for x in buckets.values()), float(THRESH)))
print()
print("  M (exact)      ~value   count   smallest example                 structure note")
for M in sorted(buckets)[:16]:
    fams = buckets[M]
    ex = min(fams)
    # structure heuristics
    note = ""
    if ex == tuple(range(1, 14)):
        note = "= AP {1..13} (TIGHT, M=1/14)"
    elif list(ex[:-1]) == list(range(1, 13)):
        note = "{1..12, %d}" % ex[-1]
    elif list(ex[:11]) == list(range(1, 12)) and ex[11] == 13:
        note = "{1..11,13, %d}" % ex[-1]
    print("  %-13s %.5f  %5d   %-32s %s" % (str(M), float(M), len(fams), str(ex), note))

print()
print("  -> the GAP structure just above 1/14:")
vals = sorted(buckets)
for i, M in enumerate(vals[:6]):
    lo = float(M); nxt = float(vals[i+1]) if i+1 < len(vals) else None
    print("     M=%s=%.6f  (%d families)%s" % (M, float(M), len(buckets[M]),
          "" if nxt is None else "   gap to next: %.6f" % (nxt-lo)))

print()
print("=" * 78)
print("PART 2 -- opus-S130 density floor mu_{1/7}(E) = meas{x: maxgap{frac(e x)} > 1/7}")
print("=" * 78)

def mu_maxgap(E, thr=1.0/7, res=200000):
    """meas{x in [0,1): circular max-gap of {frac(e_i x)} > thr} by fine sampling."""
    xs = (np.arange(res) + 0.5) / res
    Ea = np.array(E, dtype=np.float64)
    ph = np.mod(np.outer(xs, Ea), 1.0)          # res x k
    ph.sort(axis=1)
    gaps = np.diff(ph, axis=1)
    wrap = (ph[:, 0] + 1.0) - ph[:, -1]
    maxgap = np.maximum(gaps.max(axis=1), wrap)
    return (maxgap > thr).mean()

print("\n(reproduce opus's exact AP table via sampling, res=200k):")
print("   k   mu_1/7(AP={1..k})   opus-exact")
opus = {8:691/735, 9:247/294, 10:38/49, 11:1381/2205, 12:13823/24255, 13:477/1078}
for k in range(8, 14):
    ap = list(range(1, k+1))
    print("   %-3d %.5f            %.5f" % (k, mu_maxgap(ap), opus[k]))

print("\nAP-MINIMALITY test (k=13): is mu_1/7(AP) the MINIMUM over 13-sets?")
ap13 = list(range(1, 14)); base = mu_maxgap(ap13)
print("   mu_1/7(AP_13) = %.5f (= 477/1078 = %.5f)" % (base, 477/1078))
random.seed(54); beat = 0; worst_below = base
for _ in range(200):
    E = random.sample(range(1, 200), 13)
    m = mu_maxgap(E, res=40000)
    if m < base - 0.003:
        beat += 1; worst_below = min(worst_below, m)
# structured adversaries
for name, E in [("geometric", [2**i for i in range(13)]),
                ("doubled AP", [2*i for i in range(1,14)]),
                ("dilated AP*7", [7*i for i in range(1,14)]),
                ("{1..12,26}", list(range(1,13))+[26]),
                ("primes", [2,3,5,7,11,13,17,19,23,29,31,37,41])]:
    m = mu_maxgap(E, res=40000)
    flag = " <-- BELOW AP!" if m < base - 0.005 else ""
    print("   %-14s mu_1/7 = %.5f%s" % (name, m, flag))
print("   random 13-sets beating AP (by >0.003): %d/200%s" %
      (beat, "" if beat == 0 else "  (min found %.5f)" % worst_below))

print("\n" + "=" * 78)
print("READOUT: the census reveals the near-tight parametric families (Part 1) and")
print("confirms/tests opus's AP-minimal density floor (Part 2). Structural targets below.")
print("=" * 78)
