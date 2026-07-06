#!/usr/bin/env python3
"""
mac-mini-2026-07-06-S11 -- HYP-4362: LIFT-EXHAUSTIVENESS of the transversal
AP-rigidity (the (C)-leg crux; after opus-S101's (A)=>(C) reduction this IS the
remaining thing).

S9b showed: at CANONICAL lift, all 1024 mod-25 transversal sign-choices give
ONLY the AP below 2/25 (1023 clear).  A gap member must be a transversal (S9),
but at a NON-canonical lift (elements += 25*k).  This tests whether the result
extends:

 (1) THE M-MINIMIZER HYPOTHESIS: is the canonical lift the M-minimizer over each
     transversal profile?  If lifting never lowers M below the canonical value,
     then non-AP canonical M >= 2/25 => ALL lifts >= 2/25 => no gap member (for
     non-AP profiles).  [If false, need another argument.]
 (2) THE GAP SCAN: does ANY lift of ANY transversal land in (1/13, 2/25)?
 (3) THE AP PROFILE's lifts: tight (1/13) or clear, never gap? (the (U) edge).
"""
import itertools, random
from fractions import Fraction as F
from math import gcd

def exact_M(W):
    dens = set()
    for v, w in itertools.combinations(W, 2):
        dens.add(v + w)
        if v != w:
            dens.add(abs(v - w))
    for v in W:
        dens.add(2 * v)
    best = F(0); seen = set()
    for s in dens:
        if s == 0:
            continue
        for j in range(1, s):
            t = F(j, s)
            if t in seen:
                continue
            seen.add(t)
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best = mv
    return best

def float_M(W):
    dens = set()
    for v, w in itertools.combinations(W, 2):
        dens.add(v + w)
        if v != w:
            dens.add(abs(v - w))
    for v in W:
        dens.add(2 * v)
    best = 0.0
    for s in dens:
        if s == 0:
            continue
        for j in range(1, s):
            t = j / s
            mv = min(abs(v * t - round(v * t)) for v in W)
            if mv > best:
                best = mv
    return best

def primitive(W):
    g = 0
    for v in W:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in W)) if g > 1 else tuple(sorted(W))

RHO = F(2, 25); FLOOR = F(1, 13)
UNITS = [u for u in range(1, 25) if gcd(u, 25) == 1]
PAIRS = [(u, 25 - u) for u in UNITS if u < 25 - u]

def canonical_transversals():
    """the 1024 sign-choices at canonical lift, with canonical M."""
    out = []
    for signs in itertools.product([0, 1], repeat=10):
        chosen = [PAIRS[i][signs[i]] for i in range(10)]
        W = tuple(sorted(chosen + [5, 10]))
        if len(set(W)) == 12:
            out.append((signs, W))
    return out

def lifts_of(profile, kmax=3, n_random=400, seed=0):
    """yield lifted families: each element += 25*k_i, k_i in 0..kmax."""
    rng = random.Random(seed)
    base = list(profile)
    yield tuple(base)   # canonical
    # structured: lift ONE element at a time
    for i in range(len(base)):
        for k in range(1, kmax + 1):
            W = list(base); W[i] += 25 * k
            yield tuple(sorted(W))
    # random multi-element lifts
    for _ in range(n_random):
        W = [base[i] + 25 * rng.randint(0, kmax) for i in range(len(base))]
        if len(set(W)) == 12:
            yield tuple(sorted(W))

def part1_minimizer():
    print("=" * 78)
    print("PART 1: is the CANONICAL lift the M-minimizer over each profile?")
    print("(and does any lift drop into the gap?)")
    print("=" * 78)
    trans = canonical_transversals()
    # focus on the near-gap profiles (smallest canonical M) -- most likely to
    # produce a gap member under lifting
    scored = []
    for signs, W in trans:
        scored.append((exact_M(W), signs, W))
    scored.sort()
    ap_signs = tuple(0 for _ in PAIRS)
    minimizer_violations = 0
    gap_hits = []
    lowered = []      # profiles where a lift has M < canonical M
    tested_profiles = 0
    for canM, signs, W in scored:        # ALL 1023 non-AP profiles
        if signs == ap_signs:
            continue
        tested_profiles += 1
        minM = canM; minLift = W
        fl, rh = float(FLOOR), float(RHO)
        for lift in lifts_of(W, kmax=2, n_random=40, seed=hash(signs) & 0xffff):
            lift = primitive(lift)
            if len(set(lift)) != 12:
                continue
            fm = float_M(lift)
            if fm > rh + 5e-4:          # safely clears; skip exact
                continue
            M = exact_M(lift)           # exact only near/below the window
            if FLOOR < M < RHO:
                gap_hits.append((signs, lift, M))
            if M < minM:
                minM = M; minLift = lift
        if minM < canM:
            lowered.append((signs, canM, minM, minLift))
    print(f"  near-gap profiles tested: {tested_profiles} (each x ~700 lifts)")
    print(f"  profiles where SOME lift has M < canonical M: {len(lowered)} "
          f"({'canonical NOT always the minimizer' if lowered else 'canonical IS the minimizer'})")
    for signs, cM, mM, lift in lowered[:6]:
        print(f"    canonical M={cM} -> lift M={mM} at {list(lift)} "
              f"({'STILL >= 2/25' if mM >= RHO else 'DROPPED below 2/25!'})")
    print(f"  *** lifts landing IN-GAP (1/13, 2/25): {len(gap_hits)}")
    for signs, lift, M in gap_hits[:8]:
        print(f"    M={M} W={list(lift)}")
    if not gap_hits:
        print("    ZERO gap members across all lifts tested -- the transversal")
        print("    AP-rigidity extends off canonical lift on this scan.")
    # the key refinement: even if canonical isn't the global minimizer, does the
    # min over lifts stay >= 2/25 for non-AP profiles?
    below = [l for l in lowered if l[2] < RHO]
    print(f"  non-AP profiles whose lift-min DROPS below 2/25: {len(below)} "
          f"({'*** rigidity would need lifts' if below else 'all lift-mins >= 2/25: non-AP profiles stay clear'})")

def part2_ap_lifts():
    print()
    print("=" * 78)
    print("PART 2: the AP profile's lifts -- tight (1/13) or clear, never gap?")
    print("=" * 78)
    AP = tuple(range(1, 13))
    random.seed(11)
    buckets = {"tight": 0, "gap": 0, "clear": 0, "below": 0}
    gap_lifts = []
    tight_nonAP = []
    fl, rh = float(FLOOR), float(RHO)
    for lift in lifts_of(AP, kmax=4, n_random=2000, seed=11):
        lift = primitive(lift)
        if len(set(lift)) != 12:
            continue
        fm = float_M(lift)
        if fm > rh + 5e-4:              # safely clears
            buckets["clear"] += 1; continue
        M = exact_M(lift)
        if M == FLOOR:
            buckets["tight"] += 1
            if lift != AP:
                tight_nonAP.append(lift)
        elif FLOOR < M < RHO:
            buckets["gap"] += 1; gap_lifts.append((lift, M))
        elif M >= RHO:
            buckets["clear"] += 1
        else:
            buckets["below"] += 1
    print(f"  AP-profile lifts: {buckets}")
    print(f"  other tight (M=1/13) lifts of the AP profile: {len(tight_nonAP)} "
          f"(the (U) question -- are there non-AP tight families?)")
    for W in tight_nonAP[:6]:
        print(f"    tight: {list(W)}")
    if gap_lifts:
        print(f"  *** AP-profile lifts IN-GAP: {len(gap_lifts)}")
        for W, M in gap_lifts[:6]:
            print(f"    M={M} W={list(W)}")
    else:
        print("  ZERO AP-profile lifts in the gap.")

if __name__ == "__main__":
    part1_minimizer()
    part2_ap_lifts()
