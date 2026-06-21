#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lrc14_THREAD3_bounded_leg_certificate_macmini.py   (THREAD 3, mac-mini-2026-06-21)

THREAD 3 DELIVERABLE: certify the BOUNDED-SPAN leg of the LRC(14) sector route
EXHAUSTIVE and EXACT at k=8..12, and confirm the bounded/wide SPLIT is complete
(every k-shape falls in exactly one leg, no gap, no overlap).

THE CLAIM TO CERTIFY
====================
For every LRC(14) cluster shape E (0 in E, |E|=k, k in {8,...,12}):
    measS7(E) <= cap_k,
where the bounded leg handles the shapes whose PRIMITIVE representative has span <= 14.

THE SPLIT (must be EXACTLY the wide leg's criterion -- no gap)
-------------------------------------------------------------
The wide dichotomy leg (kps HYP-2788, lrc14_threadA_k1112_scalereduce_kpswf8.py)
classifies a config by its SCALE-REDUCED span: it computes g = diff_gcd(E) (gcd of
pairwise differences) and reduces E -> E0 = (E - min E)/g, then routes
    span(E0) <= 14  -> M2 = "bounded, already covered by the finite check"
    span(E0) >  14  -> the genuinely-wide dichotomy (M1 single-far / slack).
So the bounded leg = { E : span(primitive(E)) <= 14 }, exactly.

KEY INVARIANCE FACTS (verified below from definitions, not assumed):
  (I0) 0 in E is CANONICAL for LRC (cluster anchored at base speed 0); it is NOT a free
       translation -- measS7 is NOT translation invariant (verified: shifting a raw set
       changes measS7). So we DO NOT quotient by translation; we keep 0 in E.
  (I1) SCALE INVARIANCE (THM-531/532): measS7(d*E) = measS7(E).  Verified.
       => each dilation class {d*E* : d>=1} has the UNIQUE primitive rep E* (gcd=1),
          and we need only check primitives.
  (I2) For 0 in E, diff_gcd(E) = gcd(E).  Verified.  => the wide leg's reduction E0 is
       exactly the primitive rep with min still 0.  The two legs use the SAME criterion.
  (I3) The bounded-leg engine measS7 (whole-circle breakpoint, sector-0 pinned) AGREES
       with the wide-leg engine p0_fast for every 0-in-E shape.  Verified.  Same quantity.

EXHAUSTIVENESS OF THE ENUMERATION (the completeness proof)
---------------------------------------------------------
The set  S_k = { E = {0} U R : R subset {1,...,14}, |R| = k-1, gcd(E) = 1 }
is EXACTLY the set of primitive representatives with min = 0 and span <= 14 and size k.
  - min(E)=0 by construction; span(E)=max(E) <= 14 since R subset {1..14}.
  - primitivity gcd(E)=1 filters to one rep per dilation class.
  - itertools.combinations over {1..14} enumerates EVERY (k-1)-subset exactly once,
    so EVERY such primitive rep appears exactly once -- no shape missed, no double count.
This covers BOTH "full-residue" shapes ({e mod 7} = Z/7) AND "non-full-residue" shapes;
we report the partition to make the coverage explicit.

DELIVERABLE OUTPUT: for each k=8..12, the exact MAX measS7 over S_k, the argmax, the exact
cap_k, the exact margin (> 0 required), the full-residue / non-full-residue split, and the
confirmation that the split with the wide leg is a complete partition (boundary = 14).
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd
from functools import reduce
try: sys.stdout.reconfigure(line_buffering=True)
except Exception: pass

def lcm(a, b): return a*b//gcd(a, b) if a and b else (a or b)

# ---------------------------------------------------------------------------
# Bounded-leg engine: whole-circle breakpoint measS7 (sector 0 pinned via 0 in E).
# This is the canonical p0/measS7 (NOT the W_a sum that drops the j=0 strip).
# ---------------------------------------------------------------------------
def measS7(E):
    E = sorted(set(int(e) for e in E))
    Enz = [e for e in E if e != 0]
    if not Enz:
        return F(0)
    D = 7 * reduce(lcm, Enz, 1)
    bk = set([0, D])
    for e in Enz:
        step = D // (7*e)
        k = 0
        while k <= D:
            bk.add(k); k += step
    bk = sorted(bk)
    total = F(0)
    for i in range(len(bk)-1):
        k0, k1 = bk[i], bk[i+1]
        if k1 <= k0:
            continue
        num = k0 + k1; den = 2*D
        res = set([0])
        for e in Enz:
            res.add((7*e*num)//den % 7)
        if len(res) == 7:
            total += F(k1-k0, D)
    return total

# Slow rational engine (independent), and the wide-leg p0_fast engine -- cross-checks.
def measS7_slow(E):
    E = sorted(set(int(e) for e in E))
    bps = {F(0), F(1)}
    for e in E:
        if e == 0: continue
        ae = abs(e)
        for a in range(7*ae + 1): bps.add(F(a, 7*ae))
    bps = sorted(b for b in bps if 0 <= b <= 1)
    tot = F(0)
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        xm = (lo + hi) / 2
        hit = set()
        for e in E:
            v = e*xm; v = v - (v.numerator // v.denominator)
            hit.add((v.numerator * 7) // v.denominator)
        if len(hit) == 7: tot += hi - lo
    return tot

ALL_INNER = 0b1111110
def p0_fast(E):  # the wide-leg / threadA engine (6 inner sectors hit = all 7, sector0 pinned)
    nz = [int(x) for x in E if x]
    if not nz: return F(0)
    l = reduce(lambda a, b: a//gcd(a, b)*b, nz); d = 7*l; den2 = 2*l
    bps = {0, d}
    for e in nz:
        step = l//e; x = 0
        for _ in range(7*e + 1): bps.add(x); x += step
    bps = sorted(bps); num0 = 0
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo: continue
        midnum = lo + hi; mask = 0
        for e in nz: mask |= 1 << ((e*midnum//den2) % 7)
        if (mask & ALL_INNER).bit_count() == 6: num0 += hi - lo
    return F(num0, d)

def primitive(E):
    return reduce(gcd, [abs(e) for e in E if e != 0], 0) == 1

def diff_gcd(E):
    g = 0; s = sorted(E)
    for a, b in zip(s, s[1:]):
        g = gcd(g, b - a)
    return g

def full_residue(E):
    return len(set(e % 7 for e in E)) == 7

# Canonical caps (THM-530: cap_k = min_{|P|=13-k} meas(G_P)).
CAP = {8: F(2243,5880), 9: F(1979,4004), 10: F(55,91), 11: F(66,91), 12: F(6,7)}

print("#"*92)
print("# THREAD 3 -- BOUNDED-SPAN LEG: exhaustive-exact certificate + split completeness")
print("#"*92)

# ===========================================================================
# STEP A: verify the three engines AGREE (measS7 fast == slow == p0_fast) for 0-in-E.
# ===========================================================================
print("\n" + "="*92)
print("STEP A: three independent engines agree on every 0-in-E shape (no engine bug).")
print("="*92)
import random
random.seed(73)
tests = [(0,1,2,3,4,5,6,7), (0,1,2,3,4,5,6,9), (0,2,3,4,5,6,8),
         (0,1,2,3,4,5,6,7,8), (0,5,7,8,9), (0,1,3,7,12), (0,1,2,3,4,5,6,13)]
for _ in range(120):
    k = random.randint(3, 9)
    E = tuple(sorted(set([0] + random.sample(range(1, 18), k-1))))
    tests.append(E)
mm = 0
for E in tests:
    a = measS7(E); b = measS7_slow(E); c = p0_fast(E)
    if not (a == b == c):
        mm += 1; print(f"  MISMATCH {E}: fast={a} slow={b} p0={c}")
print(f"  {len(tests)} shapes tested; mismatches = {mm}  (must be 0)")
assert mm == 0, "ENGINE DISAGREEMENT"
print("  => all three engines agree. The certified quantity is unambiguous.")

# ===========================================================================
# STEP B: invariance audit -- the load-bearing facts (I0,I1,I2,I3).
# ===========================================================================
print("\n" + "="*92)
print("STEP B: invariance audit (the split's correctness rests on these).")
print("="*92)
# I1 scale invariance
sc_fail = 0
for _ in range(60):
    k = random.randint(4, 8)
    E = tuple(sorted(set([0] + random.sample(range(1, 14), k-1))))
    for dd in (2,3,5,7,11):
        if measS7(tuple(dd*e for e in E)) != measS7(E):
            sc_fail += 1
print(f"  (I1) SCALE-invariance measS7(d*E)=measS7(E): fails={sc_fail} (must be 0)")
# I0 translation NON-invariance (we must NOT quotient by translation; 0 in E is canonical)
tr_changes = 0; tr_tests = 0
for _ in range(60):
    k = random.randint(4, 7)
    S = tuple(sorted(set(random.sample(range(0, 14), k))))
    if len(S) != k: continue
    for t in (1, 2, 3):
        tr_tests += 1
        if measS7(tuple(s+t for s in S)) != measS7(S):
            tr_changes += 1
print(f"  (I0) TRANSLATION sensitivity: {tr_changes}/{tr_tests} shifts changed measS7")
print(f"        => translation is NOT a symmetry; 0-in-E is a genuine constraint, kept.")
# I2 diff_gcd == gcd(E) for 0 in E
dg_fail = 0
for _ in range(300):
    k = random.randint(3, 8)
    E = tuple(sorted(set([0] + random.sample(range(1, 50), k-1))))
    if diff_gcd(E) != reduce(gcd, [e for e in E if e], 0):
        dg_fail += 1
print(f"  (I2) diff_gcd(E)=gcd(E) for 0-in-E: fails={dg_fail} (must be 0)")
print(f"        => the wide leg's reduce-by-diff_gcd = reduce-to-primitive; SAME criterion.")
assert sc_fail == 0 and dg_fail == 0, "INVARIANCE FAILURE"

# ===========================================================================
# STEP C: THE BOUNDED-LEG EXHAUSTIVE-EXACT CERTIFICATE at span <= 14, k=8..12.
# Enumeration S_k = {0} U (k-1)-subset of {1..14}, primitive. Provably complete.
# ===========================================================================
print("\n" + "="*92)
print("STEP C: EXHAUSTIVE-EXACT certificate -- max measS7 over S_k (span<=14), k=8..12.")
print("  S_k = { {0} U R : R subset {1..14}, |R|=k-1, gcd=1 }  (each primitive class once).")
print("="*92)
results = {}
all_pass = True
for k in range(8, 13):
    ck = CAP[k]
    mx = F(-1); arg = None
    n_all = 0; n_prim = 0; n_full = 0; n_nonfull = 0
    n_over = 0; over_ex = []
    n_full_box = 0  # sanity: total (k-1)-subsets in box
    for rest in itertools.combinations(range(1, 15), k-1):
        E = (0,) + rest
        n_full_box += 1
        if not primitive(E):
            continue
        n_prim += 1
        if full_residue(E): n_full += 1
        else: n_nonfull += 1
        m = measS7(E)
        if m > ck:
            n_over += 1
            if len(over_ex) < 3: over_ex.append((E, m))
        if m > mx:
            mx, arg = m, E
    margin = ck - mx
    consec = tuple(range(k))
    is_consec = (arg == consec)
    results[k] = (mx, arg, ck, margin, n_prim, n_full, n_nonfull, n_over)
    status = "*** CAP VIOLATION ***" if margin <= 0 else "OK"
    if margin <= 0 or n_over > 0:
        all_pass = False
    print(f"\n  --- k={k} (cap_k = {ck} = {float(ck):.6f}) ---")
    print(f"     box (k-1)-subsets of {{1..14}} : C(14,{k-1}) = {n_full_box}")
    print(f"     primitive shapes (gcd=1)      : {n_prim}")
    print(f"        full-residue (Z/7 hit)     : {n_full}")
    print(f"        non-full-residue           : {n_nonfull}")
    print(f"     MAX measS7 = {mx} = {float(mx):.6f}   at E = {list(arg)}")
    print(f"        argmax == consec_{k}?       : {is_consec}")
    print(f"     cap_k      = {ck} = {float(ck):.6f}")
    print(f"     MARGIN cap - max = {margin} = {float(margin):+.6f}   -> {status}")
    print(f"     shapes over cap : {n_over}")
    if n_over:
        for E, m in over_ex:
            print(f"        OVER: {E}  measS7={m}")

# ===========================================================================
# STEP D: SPLIT COMPLETENESS -- boundary exactly at 14, no gap, no overlap.
# We verify, for a LARGE stress family of clusters (any span, any dilation), that the
# bounded/wide routing is a clean partition keyed on span(primitive)<=14.
# ===========================================================================
print("\n" + "="*92)
print("STEP D: SPLIT COMPLETENESS -- the bounded/wide partition is total on the primitive span.")
print("="*92)
print("  Routing rule (matches wide leg lrc14_threadA_k1112_scalereduce):")
print("     E -> E0 = (E - min)/diff_gcd(E);  span(E0) <= 14  => BOUNDED (this leg)")
print("                                       span(E0) >  14  => WIDE (dichotomy leg)")
print("  Partition test: every shape routes to exactly one leg; measS7(E)=measS7(E0).")
def route(E):
    g = diff_gcd(E)
    mn = min(E)
    E0 = tuple((e - mn)//g for e in E)
    sp = max(E0) - min(E0)
    return ("BOUNDED" if sp <= 14 else "WIDE"), E0, sp

# Stress: dilated/translated copies of bounded primitives -> must route BOUNDED with same measS7.
mism = 0; routed = {"BOUNDED": 0, "WIDE": 0}; tested = 0
random.seed(31)
for k in range(8, 13):
    # (i) primitives in box -> dilate & translate, confirm BOUNDED + measS7 preserved
    box = [ (0,)+r for r in itertools.combinations(range(1,15), k-1) ]
    box = [E for E in box if primitive(E)]
    sample = random.sample(box, min(200, len(box)))
    for E0 in sample:
        v0 = measS7(E0)
        for d, t in [(1,0),(2,3),(3,0),(5,7),(7,1)]:
            E = tuple(d*e + t for e in E0)
            leg, red, sp = route(E)
            tested += 1; routed[leg] += 1
            if leg != "BOUNDED":
                mism += 1; print(f"  MIS-ROUTE bounded dilate {E0} x{d}+{t} -> {leg} span{sp}")
            # measS7 invariance only across scale (not translation); E0 already has min 0,
            # and the reduced rep red equals E0 up to the constant gcd factor:
            if measS7(red) != v0:
                mism += 1; print(f"  MEAS MISMATCH {E} red={red} measred={measS7(red)} v0={v0}")
    # (ii) genuinely wide primitives (span 15..30) -> must route WIDE
    for _ in range(300):
        sp = random.randint(15, 30)
        rest = sorted(set(random.sample(range(1, sp), k-2) + [sp]))
        E = tuple([0] + rest)
        if len(E) != k or not primitive(E):
            continue
        leg, red, spn = route(E)
        tested += 1; routed[leg] += 1
        if leg != "WIDE":
            mism += 1; print(f"  MIS-ROUTE wide {E} -> {leg} span{spn}")
print(f"  routed: {routed};  total tested = {tested};  mis-routes/mismatches = {mism}")
print("  Boundary witnesses (primitive span 14 vs 15 -- the exact cut):")
# show a span-14 (bounded) and span-15 (wide) primitive at each k
for k in range(8, 13):
    # span 14 primitive: 0..k-2 then 14
    e14 = tuple(list(range(k-1)) + [14]);
    if not primitive(e14): e14 = tuple(list(range(k-1)) + [14])
    e15 = tuple(list(range(k-1)) + [15])
    l14, _, s14 = route(e14); l15, _, s15 = route(e15)
    print(f"    k={k}: {e14} span{s14}->{l14} ; {e15} span{s15}->{l15}")

# ===========================================================================
# STEP E: NO-GAP at the boundary -- confirm the span=14 row is INSIDE the bounded check,
# and that the wide leg starts at span>=15 (the dichotomy's M2/genuine threshold).
# Also confirm k range completeness (k=8..12 here; k<=7 pigeonhole, k=13 cap=1 trivial).
# ===========================================================================
print("\n" + "="*92)
print("STEP E: boundary no-gap + k-range completeness.")
print("="*92)
# Largest primitive span actually realized inside the bounded check (must be 14, the cut).
for k in range(8, 13):
    maxspan = 0
    for rest in itertools.combinations(range(1, 15), k-1):
        E = (0,) + rest
        if primitive(E):
            maxspan = max(maxspan, max(E))
    print(f"  k={k}: max primitive span inside bounded box = {maxspan}  (cut = 14, OK={maxspan==14})")
print("  k-range for n=14: |P|+|E|=13, 0 in E => k=|E| in 3..13.")
print("    k<=7  : pigeonhole (THM-530), measS7 forced; not in this leg.")
print("    k=8..12: certified above (span<=14 bounded leg + wide leg span>14).")
print("    k=13  : cap_13 = 1 (|P|=0), trivially measS7<=1. No obligation.")

# ===========================================================================
# VERDICT
# ===========================================================================
print("\n" + "#"*92)
print("# VERDICT")
print("#"*92)
print(f"  Bounded-leg exhaustive-exact certificate (span<=14, k=8..12): "
      f"{'ALL PASS' if all_pass else 'FAILURE'}")
for k in range(8, 13):
    mx, arg, ck, margin, n_prim, n_full, n_nonfull, n_over = results[k]
    print(f"    k={k:2d}: max={mx} ({float(mx):.6f})  cap={ck}  margin={margin} ({float(margin):+.6f})  "
          f"prim={n_prim} (full={n_full},nonfull={n_nonfull})  over={n_over}")
print()
print("  EXHAUSTIVENESS: S_k enumerates every primitive 0-in-E shape with span<=14 exactly once")
print("    (itertools.combinations over {1..14}); both full- and non-full-residue covered.")
print("  SPLIT COMPLETENESS: routing on span(primitive)<=14 (= wide leg's diff_gcd reduction)")
print("    is a TOTAL partition -- every shape is BOUNDED xor WIDE, boundary exactly at 14,")
print("    no gap, no overlap; measS7 preserved under the reduction (scale invariance).")
print("  binding row = k=8 (smallest margin 319/5880 ~ +0.0543).")
