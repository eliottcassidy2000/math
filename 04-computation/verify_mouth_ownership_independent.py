#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
INDEPENDENT adversarial verification of the drop-6 MOUTH-OWNERSHIP GEOMETRY claim.
Written from scratch, exact Fraction arithmetic, no reuse of the engines under test.

We re-derive everything independently:
  - lonely measure meas(G_C) = meas{t in [0,1): ||d t|| > 1/14 for all d in C}
  - the safe components (mouths) and their exact endpoints
  - which speed(s) wall each component
  - the determinant numerator of each R->L gap
  - the effect of deleting tower bits {1,2,4,8}
  - the two-tail min row geometry
"""
from fractions import Fraction
from functools import reduce
from math import gcd

THETA = Fraction(1, 14)
THR1  = Fraction(7, 858)      # drop-6 mouth, global min
THR2  = Fraction(426, 35035)  # AP one-hole second value

# ---------------------------------------------------------------------------
# Independent danger-arc / safe-component engine.
# Danger arcs of speed d: { t : ||d t|| <= 1/14 } = union_m [ (m/d) - 1/(14d), (m/d) + 1/(14d) ].
# We tag each safe endpoint with which (speed, m, side) walls it.
# ---------------------------------------------------------------------------

def danger_arcs_tagged(d):
    """Return list of (lo, hi, left_wall_info, right_wall_info) wrapped into [0,1).
       Each arc on the real line is [ m/d - 1/(14d), m/d + 1/(14d) ].
       left wall point = m/d - 1/(14d) = (14m-1)/(14d)
       right wall point = m/d + 1/(14d) = (14m+1)/(14d)
       We carry the wall label (d, m, 'L'/'R')."""
    out = []
    den = 14 * d
    for m in range(0, d + 1):
        lo = Fraction(14*m - 1, den)
        hi = Fraction(14*m + 1, den)
        # wall labels: lo is an 'L' wall (left edge of danger -> right edge of a safe gap on its left),
        # hi is an 'R' wall.  But for the SAFE gap, the gap's LEFT endpoint = a danger 'hi' (R wall),
        # and the gap's RIGHT endpoint = a danger 'lo' (L wall).  We store both.
        out.append((lo, hi, (d, m, 'L'), (d, m, 'R')))
    return out

def safe_components(core):
    """Exact safe components with endpoint wall ownership.
       Returns list of (lo, hi, owners_at_lo, owners_at_hi) where owners are sets of (speed,m,side)."""
    # Collect all danger intervals clipped to [0,1), and a map point -> set of wall tags.
    intervals = []
    wall_at = {}  # point -> set of (d, m, side)
    for d in core:
        if d == 0:
            continue
        for lo, hi, lt, rt in danger_arcs_tagged(d):
            # clip to [0,1)
            a = max(lo, Fraction(0))
            b = min(hi, Fraction(1))
            if a < b:
                intervals.append((a, b))
            # record wall points (only those landing in [0,1])
            if 0 <= lo <= 1:
                wall_at.setdefault(lo, set()).add(lt)
            if 0 <= hi <= 1:
                wall_at.setdefault(hi, set()).add(rt)
    # merge danger intervals
    intervals.sort()
    merged = []
    for a, b in intervals:
        if merged and a <= merged[-1][1]:
            if b > merged[-1][1]:
                merged[-1] = (merged[-1][0], b)
        else:
            merged.append((a, b))
    # complement = safe gaps
    comps = []
    cursor = Fraction(0)
    for a, b in merged:
        if a > cursor:
            comps.append((cursor, a))
        if b > cursor:
            cursor = b
    if cursor < 1:
        comps.append((cursor, Fraction(1)))
    # attach owners
    out = []
    for lo, hi in comps:
        out.append((lo, hi, wall_at.get(lo, set()), wall_at.get(hi, set())))
    return out

def meas(core):
    return sum((hi - lo for lo, hi, _, _ in safe_components(core)), Fraction(0))

def speeds_walling(core):
    """For each safe component, which speeds OWN its two endpoints (the binding walls)."""
    res = []
    for lo, hi, ol, oh in safe_components(core):
        lo_speeds = sorted(set(s for (s, m, side) in ol))
        hi_speeds = sorted(set(s for (s, m, side) in oh))
        res.append((lo, hi, hi - lo, lo_speeds, hi_speeds, ol, oh))
    return res

def det_numerator(lo, hi, ol, oh):
    """lo is an R-wall of (v, a): lo = (14a+1)/(14v); hi is an L-wall of (w,b): hi=(14b-1)/(14w).
       hi-lo numerator over 14vw = v(14b-1) - w(14a+1)."""
    # find the R wall owning lo and the L wall owning hi
    rwalls = [(v, a) for (v, a, side) in ol if side == 'R']
    lwalls = [(w, b) for (w, b, side) in oh if side == 'L']
    # canonical: smallest (speed, tooth)
    v, a = min(rwalls)
    w, b = min(lwalls)
    return v * (14*b - 1) - w * (14*a + 1), (v, a), (w, b)

def core_from(holes, tails):
    base = [d for d in range(1, 14) if d not in set(holes)]
    return tuple(sorted(base + list(tails)))

# ---------------------------------------------------------------------------
print("="*78)
print("INDEPENDENT VERIFICATION — drop-6 mouth ownership geometry")
print("="*78)

drop6 = core_from({6}, [])
print(f"\ndrop-6 core = {drop6}")
print(f"meas(drop-6) = {meas(drop6)} = {float(meas(drop6)):.9f}  (THR1=7/858={float(THR1):.9f})")
assert meas(drop6) == THR1, "drop-6 measure mismatch!"
print("  CHECK: meas == 7/858  OK")

print("\n--- (1) The four mouths: endpoints, length, walling speeds, determinant ---")
comps = speeds_walling(drop6)
print(f"number of safe components = {len(comps)}")
all_wall_speeds = set()
for lo, hi, length, los, his, ol, oh in comps:
    detn, rw, lw = det_numerator(lo, hi, ol, oh)
    all_wall_speeds |= set(los) | set(his)
    print(f"  [{lo}, {hi}] len={length}={float(length):.6f}")
    print(f"     left endpoint walled by speeds {los}, right by {his}")
    print(f"     R{rw}->L{lw}  det_num={detn}")
print(f"\nUnion of all mouth-bounding speeds = {sorted(all_wall_speeds)}")
tower = {1, 2, 4, 8}
print(f"Tower {{1,2,4,8}} intersect mouth-walls = {sorted(all_wall_speeds & tower)}")

print("\n--- (2) Deleting each tower bit: measure + mouth survival ---")
for s in sorted(tower):
    c11 = tuple(sorted(set(drop6) - {s}))
    m = meas(c11)
    cmps = speeds_walling(c11)
    # check the 4 original mouths survive: do the same 4 intervals appear?
    orig_intervals = {(lo, hi) for lo, hi, *_ in comps}
    new_intervals  = {(lo, hi) for lo, hi, *_ in cmps}
    survived = orig_intervals & new_intervals
    print(f"  del {s}: core(11)={c11}")
    print(f"     meas={m} = {float(m):.9f}  >=THR2? {m>=THR2}  (THR2={float(THR2):.9f})")
    print(f"     original 4 mouths surviving exactly: {len(survived)}/4")

print("\n--- (3) Two-tail min row {1:-4, 5:+2, 23:+2} ---")
# core (1,2,3,5,7,8,9,11,12,13,20,46): holes {4,6,10}, tails {20,46}
two_tail = core_from({4, 6, 10}, [20, 46])
print(f"  core = {two_tail}")
mt = meas(two_tail)
print(f"  meas = {mt} = {float(mt):.9f}  >=THR2? {mt>=THR2}")
# which original mouths survive?
cmps2 = {(lo, hi) for lo, hi, *_ in speeds_walling(two_tail)}
orig = [(lo, hi, length) for lo, hi, length, *_ in comps]
print("  original mouth survival in two-tail core:")
det3_surv = Fraction(0)
for lo, hi, length, los, his, ol, oh in comps:
    detn, rw, lw = det_numerator(lo, hi, ol, oh)
    surv = (lo, hi) in cmps2
    print(f"     mouth [{lo},{hi}] det={detn} len={length}: {'SURVIVES' if surv else 'destroyed'}")

print("\n--- (4) Exhaustive 1-tail check: only below-THR2 1-tail row ---")
below = []
for holes in [(h1, h2) for h1 in range(1,14) for h2 in range(h1+1,14)]:
    for tail in range(14, 60):
        c = core_from(set(holes), [tail])
        if len(c) != 12:
            continue
        if reduce(gcd, c) != 1:
            continue
        m = meas(c)
        if m < THR2:
            below.append((m, holes, tail))
below.sort()
print(f"  below-THR2 one-tail rows (holes pair + 1 tail, tail<=59): {len(below)}")
for m, h, t in below:
    print(f"     meas={m}={float(m):.9f} holes={h} tail={t}")
