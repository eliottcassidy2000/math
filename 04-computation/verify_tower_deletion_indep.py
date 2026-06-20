#!/usr/bin/env python3
"""
INDEPENDENT adversarial verification of the LRC(14) tower-deletion lemma (THM-545).
Claim: deleting any dyadic-1 tower bit 2^a (a in {0,1,2,3}, i.e. delete 1,2,4,8)
from a primitive AP-tail 12-core C forces meas(G_C) >= thr2 = 426/35035.

Written from scratch with fractions.Fraction. NOT importing kps engine.
"""
from fractions import Fraction
import itertools
from math import gcd
from functools import reduce

THETA = Fraction(1, 14)
THR1  = Fraction(7, 858)
THR2  = Fraction(426, 35035)

def lonely_measure(C, theta=THETA, want_comps=False):
    """
    meas{ t in [0,1): ||d*t|| > theta for all d in C }, exact.
    ||x|| = distance to nearest integer. ||d t|| <= theta  <=>  t within theta/d of some m/d.
    Danger set = union over d, m of [m/d - theta/d, m/d + theta/d] (mod 1).
    Safe = complement. Engine #2: build danger as half-open events on a sorted breakpoint sweep.
    """
    # collect raw intervals, reduce mod 1 by adding integer translates that meet [0,1)
    events = []  # (point, +1 start / -1 end)
    raw = []
    for d in C:
        if d == 0:
            continue
        w = theta / d
        for m in range(0, d + 1):
            lo = Fraction(m, d) - w
            hi = Fraction(m, d) + w
            raw.append((lo, hi))
    # wrap into [0,1)
    pieces = []
    for lo, hi in raw:
        for s in (-1, 0, 1):
            a = lo + s
            b = hi + s
            a = max(a, Fraction(0))
            b = min(b, Fraction(1))
            if a < b:
                pieces.append((a, b))
    # sweep-line union via events
    for a, b in pieces:
        events.append((a, 1))
        events.append((b, -1))
    events.sort()
    covered = Fraction(0)
    depth = 0
    last = None
    comps = 0
    in_run = False
    for pt, kind in events:
        if depth > 0 and last is not None:
            covered += pt - last
        prev_depth = depth
        depth += kind
        if prev_depth == 0 and depth > 0:
            comps += 1  # new covered run starts
        last = pt
    safe = Fraction(1) - covered
    if want_comps:
        # number of SAFE components = number of gaps between covered runs (on circle->line [0,1))
        # Recompute cleanly: build merged covered intervals.
        merged = []
        pcs = sorted(pieces)
        cur = Fraction(-1)
        for a, b in pcs:
            if a > cur:
                merged.append([a, b]); cur = b
            elif b > cur:
                merged[-1][1] = b; cur = b
        # safe components = maximal gaps inside [0,1)
        nc = 0
        prev = Fraction(0)
        for a, b in merged:
            if a > prev:
                nc += 1
            prev = max(prev, b)
        if prev < Fraction(1):
            nc += 1
        return safe, nc
    return safe

def carry_profile(C):
    d = {}
    for e in C:
        if e == 0:
            continue
        m, a = e, 0
        while m % 2 == 0:
            m //= 2; a += 1
        d[m] = d.get(m, 0) + 2 ** a
    return dict(sorted(d.items()))

def ap_tail_core(holes, tails):
    return tuple(sorted([d for d in range(1, 14) if d not in holes] + list(tails)))

def ceil_frac(x):
    return x.numerator // x.denominator + (1 if x.numerator % x.denominator else 0)

# ---------------------------------------------------------------------------
# Cross-check engine against a SECOND, totally different measure method:
# directly enumerate the safe set as complement using the kps merge method.
def lonely_measure_v2(C, theta=THETA):
    arcs = []
    for d in C:
        w = theta / d
        for m in range(0, d + 1):
            arcs.append((Fraction(m, d) - w, Fraction(m, d) + w))
    segs = []
    for lo, hi in arcs:
        for s in (-1, 0, 1):
            a2 = max(lo + s, Fraction(0)); b2 = min(hi + s, Fraction(1))
            if a2 < b2:
                segs.append((a2, b2))
    segs.sort()
    cur = Fraction(-1); u = []
    for a, b in segs:
        if a > cur:
            u.append([a, b]); cur = b
        elif b > cur:
            u[-1][1] = b; cur = b
    return Fraction(1) - sum(b - a for a, b in u)

if __name__ == "__main__":
    print("=== STEP 0: engine self-consistency (two independent methods agree) ===")
    tests = [ap_tail_core({6}, []), ap_tail_core({6,10},[20]),
             ap_tail_core({4,6,10},[20,46]), ap_tail_core({12}, [])]
    for C in tests:
        v1 = lonely_measure(C)
        v2 = lonely_measure_v2(C)
        assert v1 == v2, (C, v1, v2)
    print("  two-method agreement: OK on all canonical cores")

    print("\n=== STEP 1: canonical values ===")
    d6 = lonely_measure(ap_tail_core({6}, []))
    exc = lonely_measure(ap_tail_core({6,10}, [20]))
    d12 = lonely_measure(ap_tail_core({12}, []))
    print(f"  drop-6 meas        = {d6}   (canon 7/858={THR1})   match={d6==THR1}")
    print(f"  {{5:+2}} exc meas    = {exc}  (canon 3859/420420)   match={exc==Fraction(3859,420420)}")
    print(f"  exc == 7/858+1/980 = {exc==THR1+Fraction(1,980)}")
    print(f"  drop-12 meas       = {d12}  (canon 426/35035={THR2})  match={d12==THR2}")
    print(f"  thr2 = 426/35035 = {float(THR2):.8f}")

    # =====================================================================
    print("\n=== STEP 2: VALIDATE the comb inequality empirically ===")
    # Claim: meas(G_B \ D_r) >= (6/7)*M_B - 2*c_B/(7*r) for new speed r.
    # D_r = union of arcs [m/r-1/(14r), m/r+1/(14r)]; each arc has length 1/(7r).
    # Test on many (base B, new speed r) pairs that the inequality HOLDS exactly.
    import random
    viol = 0; tested = 0
    rng = random.Random(12345)
    for a in [0,1,2,3]:
        bit = 2**a
        for h in [d for d in range(1,14) if d != bit]:
            B = tuple(d for d in range(1,14) if d not in (bit,h))
            M, c = lonely_measure(B, want_comps=True)
            for r in range(14, 90):
                if r in B: continue
                C = tuple(sorted(B + (r,)))
                if len(C) != 12: continue
                LHS = lonely_measure(C)            # meas(G_B \ D_r) since adding speed r
                RHS = Fraction(6,7)*M - Fraction(2*c, 7*r)
                tested += 1
                if LHS < RHS:
                    viol += 1
                    if viol <= 5:
                        print(f"  VIOLATION a={a} h={h} r={r}: LHS={float(LHS):.6f} < RHS={float(RHS):.6f}")
    print(f"  comb inequality tested on {tested} (B,r) pairs; violations={viol}")
    print(f"  => comb inequality {'HOLDS (no violation)' if viol==0 else 'FAILS'}")

    # =====================================================================
    print("\n=== STEP 3: independent re-run of k=1 tower-deletion proof ===")
    # For each tower bit, comb-certify the tail r>=R, exact-check residue r in [14,R).
    all_ok = True
    for a in [0,1,2,3]:
        bit = 2**a
        worstR = 0; ncores = 0; below = []
        for h in [d for d in range(1,14) if d != bit]:
            B = tuple(d for d in range(1,14) if d not in (bit,h))
            M, c = lonely_measure(B, want_comps=True)
            denom = 6*M - 7*THR2
            assert denom > 0, (a,h,M)
            # comb tail certified for all r where (6/7)M - 2c/(7r) >= thr2
            #   <=> r >= 2c / (6M - 7*thr2)
            R = ceil_frac(Fraction(2*c) / denom)
            worstR = max(worstR, R)
            for r in range(14, R):
                if r in B: continue
                C = tuple(sorted(B + (r,)))
                if len(C) != 12: continue
                if reduce(gcd, C) != 1: continue
                ncores += 1
                L = lonely_measure(C)
                if L < THR2:
                    below.append((L, h, r))
        status = "PROVED" if not below else f"{len(below)} BELOW"
        if below: all_ok = False
        print(f"  a={a} (delete {bit}): comb cutoff R={worstR}, residue cores={ncores}, "
              f"below_thr2={len(below)}  => {status}")
        for L,h,r in sorted(below)[:5]:
            print(f"      BELOW meas={L} hole={h} tail={r}")
    print(f"  k=1 layer all four bits: {'ALL PROVED' if all_ok else 'COUNTEREXAMPLE FOUND'}")
