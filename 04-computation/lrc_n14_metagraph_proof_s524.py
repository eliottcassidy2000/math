#!/usr/bin/env python3
"""LRC@14 proof attempt via metagraph walk and CRT quotient.

opus-2026-06-01-S524

Strategy: LRC@14 = the observer outdegree walk d(t) on {0,...,13} reaches 13.

THE CREATIVE TWEAK: instead of the full 14-vertex pointed tournament,
collapse to the observer's OUTDEGREE d(t) = #{runners at distance >= 1/14}.
This is a walk on {0,...,13} with steps of ±1.

Key structural tools:
1. THM-387 gap monotonicity: between wrap-arounds, d is non-decreasing
   from the RIGHT and non-increasing from the LEFT
2. CRT 14=2*7: partition runners into 7 mod-7 classes of size 2 (plus {7})
3. THM-369 sieve: primitive speed sets must cover all moduli
4. The "race" at each wrap-around: d needs to reach 13 before the next
   incoming blocker

New idea: the QUOTIENT TOURNAMENT on 7 CRT classes.
Each class {a, a+7} forms a 2-runner sub-problem.
The observer beats class C iff BOTH runners in C are safe.
This gives a tournament on 7+1=8 vertices (observer + 7 classes).
LRC = observer is source in this 8-vertex quotient.
A000568(8) = 6880, but the CRT quotient constrains which classes appear.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import gcd
from functools import reduce
from collections import Counter, defaultdict


ONE = Fraction(1)
ZERO = Fraction(0)


def frac(x: Fraction) -> Fraction:
    return x - Fraction(x.numerator // x.denominator)


def dist0(x: Fraction) -> Fraction:
    f = frac(x)
    return min(f, ONE - f)


# ═══════════════════════════════════════════════════════════════
# PART 1: Observer outdegree walk d(t) for n=14
# ═══════════════════════════════════════════════════════════════

def outdegree_walk_n14(speeds, verbose=False):
    """Trace the observer outdegree d(t) for a speed set at n=14.

    d(t) = #{i : ||v_i t|| >= 1/14}
    LRC: does d(t) = 13 for some t?

    Returns the walk trajectory and key statistics.
    """
    n = 14
    thr = Fraction(1, 14)

    # Compute threshold walls: times when ||v_i t|| = 1/n
    walls = set()
    walls.add(ZERO)
    for v in speeds:
        for a in [1, n - 1]:
            for k in range(v):
                t = Fraction(k * n + a, v * n)
                if ZERO <= t < ONE:
                    walls.add(t)
        for k in range(v):
            walls.add(Fraction(k, v))

    walls = sorted(walls)
    walls_ext = walls + [ONE]

    # Trace d(t)
    max_d = 0
    max_d_t = ZERO
    d_values = []
    d_hist = Counter()

    for idx in range(len(walls)):
        t_start = walls_ext[idx]
        t_end = walls_ext[idx + 1]
        width = t_end - t_start
        if width <= 0:
            continue

        t_mid = (t_start + t_end) / 2
        d = sum(1 for v in speeds if dist0(Fraction(v) * t_mid) >= thr)
        d_values.append((float(t_mid), d))
        d_hist[d] += 1

        if d > max_d:
            max_d = d
            max_d_t = t_mid

    # Check walls too
    wall_max_d = 0
    for t in walls:
        d = sum(1 for v in speeds if dist0(Fraction(v) * t) >= thr)
        if d > wall_max_d:
            wall_max_d = d

    return {
        "max_d_open": max_d,
        "max_d_wall": wall_max_d,
        "max_d_t": float(max_d_t),
        "lonely": max_d == 13 or wall_max_d == 13,
        "d_hist": dict(d_hist),
        "d_values": d_values,
        "num_cells": len(d_values),
    }


def part1_outdegree_walks():
    """Trace d(t) for various n=14 speed sets."""
    print("=" * 70)
    print("PART 1: Observer outdegree walk d(t) at n=14")
    print("=" * 70)
    print()

    n = 14
    test_sets = [
        ("initial", tuple(range(1, 14))),
        ("primes", (1, 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37)),
        ("evens+1", (1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24)),
        ("coprime_small", (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15)),
        ("hard_sporadic", (1, 3, 4, 5, 9, 11, 13, 15, 17, 19, 23, 25, 29)),
    ]

    for name, speeds in test_sets:
        g = reduce(gcd, speeds)
        if g > 1:
            speeds = tuple(s // g for s in speeds)

        result = outdegree_walk_n14(speeds)
        print(f"  {name:20s}: max_d(open)={result['max_d_open']:2d}  "
              f"max_d(wall)={result['max_d_wall']:2d}  "
              f"lonely={result['lonely']}  cells={result['num_cells']}")
        print(f"    d histogram: {result['d_hist']}")

        # Show the high-d region
        high_d = [(t, d) for t, d in result["d_values"] if d >= result["max_d_open"] - 1]
        if high_d:
            print(f"    near-max region ({len(high_d)} cells):")
            for t, d in high_d[:5]:
                print(f"      t≈{t:.6f}  d={d}")
        print()


# ═══════════════════════════════════════════════════════════════
# PART 2: CRT quotient tournament (14 = 2 × 7)
# ═══════════════════════════════════════════════════════════════

def crt_quotient_n14(speeds, t):
    """Build the CRT quotient tournament for n=14.

    Partition 13 runners into mod-7 classes:
    Class 0: speed 7 (singleton)
    Class k (k=1..6): speeds {k, k+7} (pairs)

    The observer beats class C iff ALL runners in C are safe.
    Class C beats class D: majority vote of pairwise runner comparisons.

    Returns: observer outdeg in quotient, quotient tournament.
    """
    n = 14
    thr = Fraction(1, 14)

    # Group speeds by mod 7
    classes = defaultdict(list)
    for v in speeds:
        classes[v % 7].append(v)

    # Observer vs each class: observer beats class iff all runners safe
    obs_beats = {}
    for cls, runners in sorted(classes.items()):
        all_safe = all(dist0(Fraction(v) * t) >= thr for v in runners)
        obs_beats[cls] = all_safe

    obs_quotient_outdeg = sum(1 for v in obs_beats.values() if v)

    return {
        "obs_quotient_outdeg": obs_quotient_outdeg,
        "obs_beats": obs_beats,
        "classes": dict(classes),
        "num_classes": len(classes),
    }


def part2_crt_quotient():
    """CRT quotient analysis for n=14."""
    print("=" * 70)
    print("PART 2: CRT quotient tournament (14 = 2×7)")
    print("=" * 70)
    print()

    speeds = tuple(range(1, 14))
    n = 14
    thr = Fraction(1, 14)

    # Show class structure
    classes = defaultdict(list)
    for v in speeds:
        classes[v % 7].append(v)
    print("Mod-7 classes:")
    for cls, runners in sorted(classes.items()):
        print(f"  class {cls}: speeds {runners}")
    print()

    # Trace the quotient outdegree over time
    num_samples = 2000
    max_qd = 0
    max_qd_t = 0.0
    qd_hist = Counter()

    for s in range(1, num_samples + 1):
        t = Fraction(s, num_samples)
        result = crt_quotient_n14(speeds, t)
        qd = result["obs_quotient_outdeg"]
        qd_hist[qd] += 1
        if qd > max_qd:
            max_qd = qd
            max_qd_t = float(t)

    # Also check wall times
    walls = set()
    for v in speeds:
        for a in [1, n - 1]:
            for k in range(v):
                tt = Fraction(k * n + a, v * n)
                if ZERO <= tt < ONE:
                    walls.add(tt)
        for k in range(v):
            walls.add(Fraction(k, v))

    wall_max_qd = 0
    wall_max_t = ZERO
    for t in sorted(walls):
        result = crt_quotient_n14(speeds, t)
        qd = result["obs_quotient_outdeg"]
        if qd > wall_max_qd:
            wall_max_qd = qd
            wall_max_t = t

    num_classes = len(classes)
    print(f"Quotient tournament: observer + {num_classes} CRT classes = {num_classes+1} vertices")
    print(f"Observer needs to beat all {num_classes} classes for loneliness.")
    print(f"  max quotient outdeg (samples): {max_qd} at t≈{max_qd_t:.6f}")
    print(f"  max quotient outdeg (walls): {wall_max_qd} at t={float(wall_max_t):.6f}")
    print(f"  quotient lonely: {max_qd == num_classes or wall_max_qd == num_classes}")
    print(f"  quotient outdeg histogram: {dict(qd_hist)}")
    print()

    # Show what happens at the lonely time t=1/14
    t_lonely = Fraction(1, 14)
    result = crt_quotient_n14(speeds, t_lonely)
    print(f"At t=1/14 (regular polygon):")
    print(f"  observer quotient outdeg: {result['obs_quotient_outdeg']}")
    for cls, beats in sorted(result["obs_beats"].items()):
        runners = result["classes"][cls]
        dists = [float(dist0(Fraction(v) * t_lonely)) for v in runners]
        print(f"  class {cls} ({runners}): all_safe={beats}  distances={[f'{d:.4f}' for d in dists]}")
    print()

    # THE KEY: the quotient has only 7 classes.
    # LRC says all 7 must be "beaten" (all runners safe).
    # This is LRC on a 7-vertex problem (sort of).
    # But each "runner" in the quotient is actually 2 real runners,
    # so "beating" a class is harder than beating a single runner.
    print("THE QUOTIENT REDUCTION:")
    print(f"  LRC@14 reduces to: observer beats all {num_classes} CRT classes")
    print(f"  Each class is a 2-runner sub-problem (except class 0: 1 runner)")
    print(f"  Class beaten iff BOTH runners safe: ||at||>=1/14 AND ||(a+7)t||>=1/14")
    print(f"  For class k: this is ||kt||>=1/14 AND ||(k+7)t||>=1/14")
    print()

    # For each class, what's the measure of "class safe" time?
    print("Class-safe measures (initial segment):")
    for cls, runners in sorted(classes.items()):
        safe_count = 0
        for s in range(1, 10001):
            t = Fraction(s, 10000)
            if all(dist0(Fraction(v) * t) >= thr for v in runners):
                safe_count += 1
        print(f"  class {cls} ({runners}): safe {safe_count}/10000 = {safe_count/100:.1f}%")
    print()


# ═══════════════════════════════════════════════════════════════
# PART 3: The "pair constraint" — when are both runners safe?
# ═══════════════════════════════════════════════════════════════

def pair_constraint_analysis():
    """For each CRT pair {a, a+7}, analyze when BOTH are safe.

    ||at|| >= 1/14 AND ||(a+7)t|| >= 1/14

    Write θ = {at} and φ = {7t}. Then {(a+7)t} = {θ + φ} (mod 1).
    Both safe iff: θ ∈ [1/14, 13/14] AND {θ+φ} ∈ [1/14, 13/14]

    The forbidden region for (θ,φ) is:
    F = {(θ,φ) : θ ∈ [0,1/14)∪(13/14,1] OR {θ+φ} ∈ [0,1/14)∪(13/14,1]}

    The safe region has measure (12/14)^2 - [overlap correction] ≈ (6/7)^2 ≈ 73.5%
    """
    print("=" * 70)
    print("PART 3: Pair constraint — both runners in CRT class safe")
    print("=" * 70)
    print()

    n = 14
    thr = Fraction(1, 14)

    # For each pair (a, a+7), compute the safe measure
    for a in range(1, 7):
        v1, v2 = a, a + 7

        safe_count = 0
        num_samples = 100000
        for s in range(1, num_samples + 1):
            t = Fraction(s, num_samples)
            d1 = dist0(Fraction(v1) * t)
            d2 = dist0(Fraction(v2) * t)
            if d1 >= thr and d2 >= thr:
                safe_count += 1

        pair_safe = safe_count / num_samples

        # Compare with (6/7)^2 (independent model)
        indep = (12 / 14) ** 2
        print(f"  pair ({v1},{v2}): safe={pair_safe:.4f}  "
              f"independent={(indep):.4f}  ratio={pair_safe/indep:.4f}")

    # The singleton class 0 (speed 7)
    safe_count = 0
    for s in range(1, num_samples + 1):
        t = Fraction(s, num_samples)
        if dist0(Fraction(7) * t) >= thr:
            safe_count += 1
    print(f"  singleton (7): safe={safe_count/num_samples:.4f}  (theory: 12/14={12/14:.4f})")
    print()

    # CRITICAL: the ALL-classes-safe measure
    speeds = tuple(range(1, 14))
    all_safe = 0
    for s in range(1, num_samples + 1):
        t = Fraction(s, num_samples)
        if all(dist0(Fraction(v) * t) >= thr for v in speeds):
            all_safe += 1
    print(f"  ALL classes safe (initial segment): {all_safe}/{num_samples}")
    print(f"    = {all_safe/num_samples:.6f} (wall-only for initial segment)")
    print()

    # Product of pair-safe measures (independence assumption)
    pair_safes = []
    for a in range(1, 7):
        v1, v2 = a, a + 7
        sc = sum(1 for s in range(1, 10001)
                 if dist0(Fraction(v1) * Fraction(s, 10000)) >= thr and
                    dist0(Fraction(v2) * Fraction(s, 10000)) >= thr)
        pair_safes.append(sc / 10000)
    singleton_safe = sum(1 for s in range(1, 10001)
                         if dist0(Fraction(7) * Fraction(s, 10000)) >= thr) / 10000
    pair_safes.append(singleton_safe)

    product = 1.0
    for ps in pair_safes:
        product *= ps
    print(f"  product of class-safe measures: {product:.6f}")
    print(f"  (if classes were independent, this would be the lonely measure)")
    print()


# ═══════════════════════════════════════════════════════════════
# PART 4: The walk-on-quotient proof approach
# ═══════════════════════════════════════════════════════════════

def quotient_walk_proof():
    """Attempt to prove LRC@14 via the CRT quotient walk.

    The quotient tournament has 8 vertices: observer + 7 CRT classes.
    Observer beats class C iff both runners in C are safe.

    The walk changes this quotient at walls where:
    1. A runner crosses the 1/14 threshold (changes observer-class edge)
    2. Runners within a class swap order (doesn't change quotient directly)

    The quotient outdegree d_Q(t) counts how many classes are fully safe.
    LRC iff d_Q reaches 7.

    KEY STRUCTURAL FACT: for the quotient to reach 7, we need
    ALL 7 classes to be SIMULTANEOUSLY safe. This is a conjunction of
    7 conditions, each involving 2 runners (except class 0: 1 runner).

    PROOF IDEA: Use the regular polygon fact. At t=k/14 (gcd(k,14)=1),
    all runners form a regular 14-gon, so all are at distance exactly 1/14.
    This is the BOUNDARY of the safe region.

    For non-initial speed sets: the trajectory might not pass through
    the regular polygon, but it must pass NEAR it (THM-387 gap race).
    The quotient structure should help bound how close.
    """
    print("=" * 70)
    print("PART 4: Quotient walk proof approach for n=14")
    print("=" * 70)
    print()

    n = 14
    thr = Fraction(1, 14)

    # For the initial segment, the proof is trivial:
    # at t=1/14, all positions are k/14, all distances are k/14 or (14-k)/14
    # and min distance is 1/14 (boundary).
    print("INITIAL SEGMENT PROOF:")
    print(f"  At t=1/14: positions = {{k/14 : k=0,...,13}}")
    print(f"  Distances: min_k ||k/14|| = 1/14 = threshold")
    print(f"  This is a BOUNDARY lonely time (wall). QED for initial segment.")
    print()

    # For GENERAL speed sets: need to show the walk reaches d=13 or wall d=13.
    # Key: any PRIMITIVE speed set with 13 speeds has a common structure.

    # SIEVE ARGUMENT (THM-369):
    # A primitive set must contain speeds covering all moduli 2,...,13.
    # In particular, gcd(v_1,...,v_13) = 1 ensures the lattice generated
    # by the speeds is Z (or a sublattice of Z with index 1).

    # DENSITY ARGUMENT:
    # The total "unsafe measure" per runner is 2/14 = 1/7.
    # For 13 runners, the total unsafe measure is 13/7 ≈ 1.86.
    # But by the THREE-DISTANCE THEOREM, the gaps between consecutive
    # positions of a single runner on the circle take at most 3 values.
    # For v runners of speed v: the v positions are equally spaced at 1/v.
    # The close zone has width 2/14 around 0, containing at most
    # ceil(2v/14) = ceil(v/7) of the v positions.

    # PAIR INDEPENDENCE:
    # For coprime a, a+7: the pair (||at||, ||(a+7)t||) is equidistributed
    # on [0,1/2]^2 as t varies over [0,1). The safe region for the pair
    # is [1/14, 1/2]^2 (both distances >= 1/14). Its measure is
    # (1/2 - 1/14)^2 / (1/2)^2 = (5/14)^2 / (1/4) = 25/49 ≈ 51%.
    # Wait, that's not right. Let me recalculate.

    # For a single runner: P(||vt|| >= 1/14) = 1 - 2/14 = 12/14 = 6/7
    # For a pair (a, a+7) with gcd(a, a+7) = gcd(a, 7):
    #   If gcd(a,7)=1 (a=1,2,3,4,5,6): the pair is "coprime mod 7"
    #   and P(both safe) ≈ (6/7)^2 ≈ 73.5% (approximately independent)

    # But we need ALL 7 classes safe SIMULTANEOUSLY.
    # If they were independent: P(all safe) ≈ (6/7)^6 * (6/7) = (6/7)^7 ≈ 37.6%
    # So if independence held, the lonely measure would be ~37.6%, which is HUGE.

    print("INDEPENDENCE MODEL:")
    indep_all = (6 / 7) ** 7
    print(f"  If all 7 CRT classes were independent:")
    print(f"    P(all safe) ≈ (6/7)^7 = {indep_all:.4f} = {indep_all*100:.1f}%")
    print(f"  This would give a LARGE lonely measure.")
    print(f"  The actual lonely measure is MUCH smaller (often wall-only)")
    print(f"  because the classes are NOT independent — they're correlated")
    print(f"  through the shared time parameter t.")
    print()

    # THE REAL PROOF TARGET:
    # Show that the correlation structure cannot make P(all safe) = 0.
    # Even though classes are correlated, the correlation cannot be
    # SO negative as to make the intersection empty.

    print("PROOF TARGET:")
    print(f"  Show: for any primitive {{v_1,...,v_13}}, the set")
    print(f"    L = {{t : all ||v_i t|| >= 1/14}} is nonempty.")
    print()
    print(f"  EQUIVALENT (via CRT 14=2×7):")
    print(f"    L = intersection of 7 sets, one per CRT class.")
    print(f"    Each class-safe set has measure >= (6/7)^2 ≈ 73.5% (for coprime pairs)")
    print(f"    or 6/7 ≈ 85.7% (for the singleton class).")
    print()
    print(f"  The intersection is empty only if the 7 class-safe sets are")
    print(f"  'perfectly anti-correlated' — covering each other's holes exactly.")
    print(f"  The metagraph walk structure should prevent this.")
    print()


# ═══════════════════════════════════════════════════════════════
# PART 5: Exhaustive check of many n=14 speed sets
# ═══════════════════════════════════════════════════════════════

def exhaustive_n14_check():
    """Check LRC@14 for many primitive speed sets.

    Can't be truly exhaustive (infinite), but check structured families:
    1. Initial segment + perturbations
    2. Random-looking coprime sets
    3. Adversarial attempts (clustered speeds, large speeds)
    """
    print("=" * 70)
    print("PART 5: Checking LRC@14 for many speed sets")
    print("=" * 70)
    print()

    n = 14
    checked = 0
    lonely_found = 0
    max_open_d = 0
    hardest_set = None

    test_sets = []

    # Family 1: Initial segment perturbations
    base = list(range(1, 14))
    for i in range(13):
        for delta in [-1, 1, 2, -2]:
            perturbed = base.copy()
            perturbed[i] = base[i] + delta
            if perturbed[i] > 0 and len(set(perturbed)) == 13:
                s = tuple(sorted(perturbed))
                if reduce(gcd, s) == 1 and s not in test_sets:
                    test_sets.append(s)

    # Family 2: Coprime sets with small speeds
    for combo in combinations(range(1, 20), 13):
        if reduce(gcd, combo) == 1:
            test_sets.append(combo)
        if len(test_sets) > 300:
            break

    # Family 3: Sets with large max speed
    for combo in combinations(range(1, 25), 13):
        if reduce(gcd, combo) == 1 and max(combo) >= 20:
            test_sets.append(combo)
        if len(test_sets) > 500:
            break

    # Deduplicate
    test_sets = list(set(test_sets))[:500]

    print(f"Testing {len(test_sets)} speed sets...")

    for speeds in test_sets:
        result = outdegree_walk_n14(speeds)
        checked += 1

        if result["lonely"]:
            lonely_found += 1
        else:
            print(f"  *** NOT LONELY: {speeds} ***")

        if result["max_d_open"] > max_open_d or (
            result["max_d_open"] == max_open_d and
            result["max_d_wall"] < (13 if hardest_set else 14)
        ):
            if not result["lonely"] or result["max_d_open"] < 13:
                max_open_d = result["max_d_open"]
                hardest_set = speeds

    print(f"\nResults: {lonely_found}/{checked} lonely ({100*lonely_found/checked:.1f}%)")
    print(f"Hardest set (lowest max open d): {hardest_set}")
    print(f"  max open d = {max_open_d}")
    if hardest_set:
        result = outdegree_walk_n14(hardest_set)
        print(f"  d histogram: {result['d_hist']}")
        print(f"  wall max d: {result['max_d_wall']}")
    print()


# ═══════════════════════════════════════════════════════════════
# PART 6: Proof sketch
# ═══════════════════════════════════════════════════════════════

def proof_sketch():
    """Summarize the proof strategy and what's still needed."""
    print("=" * 70)
    print("PART 6: Proof sketch for LRC@14")
    print("=" * 70)
    print()

    print("PROVED INGREDIENTS:")
    print("  THM-381: Lonely iff observer is source in marked tournament")
    print("  THM-384: Lonely iff both adjacent gaps >= 1/14")
    print("  THM-387: Gap flow LS→LL→SL (directed, no SL→LL)")
    print("  THM-369: Sieve — primitive set covers all moduli")
    print("  HYP-1996: Only 2*Fib(11)=178 circular iso-classes for 13 runners")
    print()

    print("THE CRT QUOTIENT STRATEGY:")
    print("  14 = 2 × 7. Group 13 runners into 7 mod-7 classes.")
    print("  LRC = all 7 classes simultaneously safe.")
    print("  Each class is a 2-runner sub-problem (or 1-runner singleton).")
    print()

    print("WHAT'S NEEDED:")
    print("  1. Show that the CRT class-safe sets can't perfectly tile.")
    print("     (Their union measure is ~7*(6/7)^2 ≈ 5.3 >> 1, so massive overlap)")
    print("  2. Use the gap-race structure (THM-387) on the quotient.")
    print("  3. The quotient walk on {0,...,7} (quotient outdeg) must reach 7.")
    print("  4. Key: each class-safe set has measure ≈ (6/7)^2 ≈ 73.5%.")
    print("     The intersection of 7 such sets (measure ≈ 37.6% if independent)")
    print("     cannot be empty if the correlation is bounded.")
    print()

    print("POTENTIAL PROOF ROUTE:")
    print("  A. CORRELATION BOUND: for any two CRT classes C,D with coprime")
    print("     speed representatives, the correlation between safe_C and safe_D")
    print("     is bounded: |corr| < 1 - epsilon(n).")
    print("  B. INCLUSION-EXCLUSION: with bounded correlations and each class")
    print("     safe ≈ 73.5% of the time, the intersection is positive.")
    print("  C. CONTINUITY: positive measure implies the set is nonempty.")
    print()

    print("THE CREATIVE INSIGHT:")
    print("  The CRT quotient reduces LRC@14 to a 7-DIMENSIONAL coverage problem.")
    print("  Each dimension is a CRT class. The 'holes' in each class-safe set")
    print("  are thin strips of width 2/14 on the torus. The question is whether")
    print("  7 thin strips can cover the entire torus. Since each strip has")
    print("  width 2/14 ≈ 14.3%, and there are only 7 of them (total width < 100%),")
    print("  there must be GAPS in the coverage — and those gaps are the lonely times.")
    print()
    print("  Wait — the 7 class-UNSAFE sets each have measure ≈ (1-(6/7)^2) ≈ 26.5%.")
    print("  Their union has measure ≤ 7*26.5% = 185% >> 100%.")
    print("  So the union CAN cover [0,1). The question is whether it ALWAYS does.")
    print()
    print("  But with the CRT structure, the unsafe sets have ARITHMETIC spacing.")
    print("  They can't tile perfectly because the gcd constraints from primitivity")
    print("  force misalignment between different classes' unsafe intervals.")
    print()


# ═══════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════

def main():
    print("LRC@14 Proof via Metagraph Walk + CRT Quotient — opus-S524")
    print()

    part1_outdegree_walks()
    part2_crt_quotient()
    pair_constraint_analysis()
    quotient_walk_proof()
    exhaustive_n14_check()
    proof_sketch()


if __name__ == "__main__":
    main()
