#!/usr/bin/env python3
"""
LRC(14) covering-bound investigation (kind-pasteur-2026-06-27).

Integrates mac-mini-S59's REDIRECT (proof = covering bound, census is the easy
case's extremal) + codex-S122 THM-571 (lift sieve closes |M14|>=7) + kps gamma-trick.

Goal: gain a *comprehensive, checkable* understanding of the |M14|<=6 residual.

Three probes:
  Part 1. Confirm the reduction: every 14-free set has M >= 1/14 at t=1/14 (trivial half).
  Part 2. The covering MARGIN: minimize M over covering sets (one multiple of 14 + 12 others).
          If covering sets are bounded away from 1/14, the bound is *strictly* easier than
          the tight-locus census (which lives at exactly 1/14). This is the quantitative
          test of mac-mini's "the proof is strictly weaker than the census".
  Part 3. The optimized-gamma 14-lift sieve on |M14|<=6: for a Q-safe fine phase gamma,
          is there a lift m in {0..13} with all of R safe? Tabulate forbidden-lift counts
          by speed type (coprime-to-14 / even / 7-type) -- locate exactly where the sieve
          can fail.

All M computed EXACTLY via the crossing/peak candidate method (Fractions).
"""
from fractions import Fraction as F
from itertools import combinations
import random

# ---------- exact distance-to-nearest-integer ----------
def nrm(x):
    """||x|| = distance from x to nearest integer, as a Fraction in [0,1/2]."""
    x = x - int(x)            # frac in (-1,1)
    if x < 0:
        x += 1                # frac in [0,1)
    return min(x, 1 - x)

def mindist(speeds, t):
    """min_i ||v_i t|| at a rational t."""
    return min(nrm(v * t) for v in speeds)

# ---------- exact M(S) = max_t min_i ||v_i t|| ----------
def M_exact(speeds):
    """
    Exact loneliness M(S). The lower envelope f(t)=min_i||v_i t|| is piecewise
    linear; its local maxima occur where two sawtooths cross (t=k/(v_i+v_j) or
    k/(v_i-v_j)) or at a single sawtooth peak (t=(2k+1)/(2 v_i)). Evaluate f on
    all such candidates in (0,1/2]; return the max. Exact and complete.
    """
    speeds = sorted(set(speeds))
    cand = set()
    n = len(speeds)
    # peaks
    for v in speeds:
        for k in range(0, 2 * v):
            num = 2 * k + 1
            t = F(num, 2 * v)
            if 0 < t <= F(1, 2):
                cand.add(t)
    # crossings of pairs (sum and difference denominators)
    for i in range(n):
        for j in range(i + 1, n):
            a, b = speeds[i], speeds[j]
            for D in (a + b, b - a):
                if D <= 0:
                    continue
                for k in range(1, D):
                    t = F(k, D)
                    if 0 < t <= F(1, 2):
                        cand.add(t)
    best = F(0)
    for t in cand:
        d = mindist(speeds, t)
        if d > best:
            best = d
    return best

# ---------- helpers ----------
def is_covering(S):
    return any(v % 14 == 0 for v in S)

def gcd14(b):
    import math
    return math.gcd(b, 14)   # 1, 2, 7, or 14

ONE14 = F(1, 14)

# ============================================================
print("=" * 70)
print("LRC(14) COVERING-BOUND INVESTIGATION  (kind-pasteur S31 cont.)")
print("=" * 70)

# ---------- Part 1: the trivial half (reduction check) ----------
print("\n[Part 1] Reduction: every 14-FREE set has M >= 1/14 at t = 1/14.")
random.seed(11)
bad = 0
tested = 0
for _ in range(3000):
    S = random.sample(range(1, 60), 13)
    if is_covering(S):
        continue
    tested += 1
    # t = 1/14 : ||v/14|| >= 1/14 because v not = 0 mod 14
    d = mindist(S, F(1, 14))
    if d < ONE14:
        bad += 1
        print("   COUNTEREXAMPLE (should be impossible):", S, d)
print(f"   tested {tested} random 14-free sets; t=1/14 fails on {bad}.")
print("   => 14-free half of LRC(14) is TRIVIAL (THM-523). Confirmed." if bad == 0
      else "   => UNEXPECTED FAILURE")

# ---------- Part 2: the covering margin ----------
print("\n[Part 2] Covering MARGIN: min M over covering sets (a multiple of 14 present).")
print("         If min M is bounded away from 1/14, the bound is *strictly* weaker")
print("         than the tight-locus census {AP,GW} (which sits exactly at 1/14).")

def min_M_over_family(gen, label):
    best = None
    arg = None
    cnt = 0
    for S in gen:
        if len(set(S)) != 13:
            continue
        if not is_covering(S):
            continue
        cnt += 1
        m = M_exact(S)
        if best is None or m < best:
            best = m
            arg = tuple(sorted(S))
    print(f"   {label}: {cnt} sets, min M = {best} = {float(best):.5f}"
          + (f"  at {arg}" if arg else ""))
    return best, arg

# Family A: 13 CONSECUTIVE integers a..a+12 that contain a multiple of 14
def consec_cov():
    for a in range(1, 200):
        S = list(range(a, a + 13))
        if is_covering(S):
            yield S
mA = min_M_over_family(consec_cov(), "consecutive blocks (covering)")

# Family B: AP {1..13} with ONE entry replaced by a multiple of 14 (14,28,42,...)
def ap_swap_cov():
    base = list(range(1, 14))
    for drop in range(13):
        for mult in (14, 28, 42, 56, 70, 84, 98, 112, 126, 140):
            S = base[:drop] + base[drop + 1:] + [mult]
            yield S
mB = min_M_over_family(ap_swap_cov(), "AP with one speed -> multiple of 14")

# Family C: small bounded covering sets, exhaustive-ish over {1..20} containing 14
def small_cov():
    pool = [v for v in range(1, 21) if v != 14]
    for combo in combinations(pool, 12):
        yield list(combo) + [14]
# this is C(19,12) ~ 50k; sample to keep it quick but representative
def small_cov_sample(k=40000):
    pool = [v for v in range(1, 21) if v != 14]
    rng = random.Random(7)
    seen = 0
    while seen < k:
        combo = rng.sample(pool, 12)
        seen += 1
        yield combo + [14]
mC = min_M_over_family(small_cov_sample(), "random 12-subsets of {1..20}\\{14} + {14}")

# Family D: targeted near-AP covering (consecutive with a gap, includes 14)
def near_ap_cov():
    rng = random.Random(3)
    for _ in range(20000):
        # take a window and force a multiple of 14 in it
        a = rng.randint(2, 8)
        S = list(range(a, a + 12))      # 12 consecutive
        # append a multiple of 14 near the window
        mult = 14 * rng.randint(1, 3)
        if mult in S:
            continue
        yield S + [mult]
mD = min_M_over_family(near_ap_cov(), "12-consecutive + a small multiple of 14")

covmins = [m for m in (mA[0], mB[0], mC[0], mD[0]) if m is not None]
overall = min(covmins) if covmins else None
print(f"\n   OVERALL min M over all covering families tried: {overall} = {float(overall):.5f}")
print(f"   1/14 = {float(ONE14):.5f} ; 1/13 = {float(F(1,13)):.5f} ; 1/8 = {float(F(1,8)):.5f}")
if overall is not None:
    if overall > ONE14:
        gap = overall - ONE14
        print(f"   => covering sets observed STRICTLY above 1/14 (margin >= {gap} = {float(gap):.5f}).")
        print(f"      Margin vs 1/13? min >= 1/13 : {overall >= F(1,13)}")
    else:
        print("   => some covering set hits 1/14 (no margin) -- bound is genuinely tight on covering sets.")

# ---------- Part 3: optimized-gamma 14-lift sieve on |M14| <= 6 ----------
print("\n[Part 3] Optimized-gamma 14-lift sieve on |M14| <= 6 (the residual branch).")
print("         S = 14Q | R.  Fix Q-safe fine phase gamma (||q*gamma||>=1/14).")
print("         lifts t_m=(gamma+m)/14, m=0..13: 14Q safe at ALL; count R-forbidden lifts.")
print("         A surviving lift => M(S) >= 1/14 for THIS gamma.  Scan gamma over a grid.")

def lift_sieve_ok(S, gamma_den=4200):
    """
    Try to certify M(S)>=1/14 by the 14-lift sieve with an optimized gamma.
    Q = multiples of 14 (scaled: q=v/14). R = rest. Need gamma with all q safe,
    then a lift m with all R safe. Scan gamma on a fine rational grid.
    Returns (certified, best_gamma, surviving_m, type_counts).
    """
    Q = [v // 14 for v in S if v % 14 == 0]
    R = [v for v in S if v % 14 != 0]
    # type census of R
    tc = {1: 0, 2: 0, 7: 0}
    for b in R:
        tc[gcd14(b)] += 1
    for g in range(1, gamma_den):
        gamma = F(g, gamma_den)
        # Q-safe at this fine phase?  ||q*gamma|| >= 1/14
        if any(nrm(q * gamma) < ONE14 for q in Q):
            continue
        # try the 14 lifts
        for m in range(14):
            t = (gamma + m) / 14
            if all(nrm(b * t) >= ONE14 for b in R):
                return True, gamma, m, tc
    return False, None, None, tc

# test on the covering families restricted to |M14|<=6 (here typically |M14|=1)
print("   Testing the sieve on the covering minimizers + structured |M14|<=6 sets:")
test_sets = []
for m, arg in (mA, mB, mC, mD):
    if arg:
        test_sets.append(arg)
# add some explicit small covering sets with 1..3 multiples of 14
test_sets += [
    tuple(sorted([14] + list(range(1, 13)))),          # {1..12,14}, |M14|=1
    tuple(sorted([14, 28] + list(range(1, 12)))),      # |M14|=2
    tuple(sorted([14, 28, 42] + list(range(1, 11)))),  # |M14|=3
    tuple(sorted(list(range(2, 15)))),                 # {2..14}, |M14|=1 (mac-mini's M=1/8)
]
fails = 0
for S in test_sets:
    nM14 = sum(1 for v in S if v % 14 == 0)
    ok, g, m, tc = lift_sieve_ok(list(S))
    Mtrue = M_exact(list(S))
    status = "OK" if ok else "SIEVE-FAILS"
    if not ok:
        fails += 1
    print(f"   |M14|={nM14} R-types(1/2/7)={tc[1]}/{tc[2]}/{tc[7]}  M={float(Mtrue):.4f}  "
          f"sieve={status}" + (f" (gamma={g}, m={m})" if ok else "") + f"  S={S}")
print(f"\n   lift-sieve failures: {fails}/{len(test_sets)}")
print("   (A sieve-fail does NOT mean M<1/14 -- only that THIS certificate didn't fire.)")

print("\n" + "=" * 70)
print("DONE.")
print("=" * 70)
