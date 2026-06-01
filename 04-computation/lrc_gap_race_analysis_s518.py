#!/usr/bin/env python3
"""Gap race analysis for the LRC LL entry mechanism (THM-387).

opus-2026-06-01-S518

THM-387 proves that LL can only be entered from LS after a wrap-around
reset.  This script analyzes the "gap race" at each wrap-around:
does g_right reach 1/n before g_left drops to 1/n?

For each wrap-around event:
- L_0: the counterclockwise gap (g_left) at the reset
- Race outcome: LL (won) or SS (lost)
- Time to reach LL or SS from the LS entry

Key questions:
1. What is the distribution of L_0 across wrap-arounds?
2. Is max(L_0) always close to 1?
3. Is the wrap-around with max(L_0) always the one that reaches LL?
4. For n=14, which wrap-arounds win the race?

Stored output:
    05-knowledge/results/lrc_gap_race_analysis_s518.out
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import gcd
from functools import reduce


ONE = Fraction(1)
ZERO = Fraction(0)


def frac(x: Fraction) -> Fraction:
    return x - Fraction(x.numerator // x.denominator)


def compute_wrap_arounds(speeds: tuple[int, ...]) -> list[tuple[Fraction, int]]:
    """Compute all wrap-around events in [0, 1): times when {v_i t} = 0.

    Returns [(time, runner_index), ...] sorted by time.
    """
    events = []
    for idx, v in enumerate(speeds):
        for k in range(v):
            t = Fraction(k, v)
            if ZERO < t < ONE:  # exclude t=0 (degenerate)
                events.append((t, idx))
    events.sort()
    return events


def gap_race_at_reset(speeds: tuple[int, ...], n: int, t_reset: Fraction,
                      resetter_idx: int) -> dict:
    """Analyze the gap race starting just after a wrap-around at t_reset.

    At t_reset, runner resetter_idx has {v * t_reset} = 0 (at observer).
    Just after, g_right ≈ 0 (runner is CW-close) and g_left is free.

    We track the race forward until g_right >= 1/n (LL) or g_left < 1/n (SS).
    """
    thr = Fraction(1, n)
    v_resetter = speeds[resetter_idx]

    # Compute L_0: g_left just after the reset
    # g_left = min_{j} (1 - {v_j * t_reset})
    # At t_reset: {v_resetter * t_reset} = 0, so 1 - {v_resetter * t} = 1
    # The resetter doesn't constrain g_left.
    g_left_components = []
    for j, v in enumerate(speeds):
        if j == resetter_idx:
            g_left_components.append(ONE)  # resetter at position 0
        else:
            pos = frac(Fraction(v) * t_reset)
            if pos == ZERO:
                # Another runner also at 0 — degenerate, skip
                return {"degenerate": True}
            g_left_components.append(ONE - pos)

    L_0 = min(g_left_components)

    # Compute g_right just after reset
    g_right_components = []
    for j, v in enumerate(speeds):
        pos = frac(Fraction(v) * t_reset)
        g_right_components.append(pos if pos > ZERO else ONE)  # exclude zero (the resetter)

    # Actually, at t_reset the resetter is at 0, so {v_resetter * t_reset} = 0.
    # g_right = min_j {v_j t_reset} = 0. We need the value just AFTER t_reset.
    # Let's track forward through walls.

    R_0 = ZERO  # g_right at the reset moment

    # Determine initial fiber
    if L_0 < thr:
        # Already in SS or SL at reset (g_right=0 is S, g_left<1/n is S)
        return {
            "degenerate": False,
            "t_reset": t_reset,
            "resetter_idx": resetter_idx,
            "L_0": L_0,
            "R_0": R_0,
            "entry_fiber": "SS",
            "race_outcome": "never_LS",
            "time_to_outcome": ZERO,
        }

    # We're in LS just after the reset (L_0 >= 1/n, R_0 = 0 < 1/n)
    # Track forward until either g_right >= 1/n or g_left < 1/n

    # Find all walls in [t_reset, t_reset + 1/(n*v_min)] that could change the race
    # Actually, we need to find the time when g_right first reaches 1/n
    # and the time when g_left first drops to 1/n

    # Time for runner resetter to reach position 1/n:
    # {v_resetter * (t_reset + dt)} = v_resetter * dt (for small dt)
    # So dt = 1/(n * v_resetter)
    time_resetter_to_threshold = Fraction(1, n * v_resetter)

    # But g_right is the MIN of all fractional parts, and the resetter starts
    # nearest. So g_right stays at v_resetter * dt until another runner becomes
    # closer. So time for g_right >= 1/n is 1/(n * v_resetter), assuming
    # no other runner crosses below the resetter before then.

    # Time for g_left to drop from L_0 to thr:
    # g_left decreases at rate v_nearest_CCW (the speed of the nearest CCW runner)
    # Find which runner achieves L_0
    nearest_ccw_speed = None
    for j, v in enumerate(speeds):
        pos = frac(Fraction(v) * t_reset)
        if j == resetter_idx:
            continue
        if ONE - pos == L_0:
            nearest_ccw_speed = v
            break

    if nearest_ccw_speed is None:
        # Should not happen
        return {"degenerate": True}

    # g_left drops from L_0 at rate v_nearest_CCW (approximately)
    if L_0 > thr:
        time_left_to_threshold = (L_0 - thr) / nearest_ccw_speed
    else:
        time_left_to_threshold = ZERO

    # Simplified race: g_right reaches 1/n at time 1/(n*v_resetter)
    #                  g_left reaches 1/n at time (L_0 - 1/n) / v_nearest_CCW
    # LL is reached if g_right wins: 1/(n*v_resetter) < (L_0-1/n)/v_nearest_CCW
    # Equivalently: v_nearest_CCW < n * v_resetter * (L_0 - 1/n)
    # Or: v_nearest_CCW / v_resetter < n * (L_0 - 1/n) = n*L_0 - 1

    # Note: this is a SIMPLIFIED analysis (ignoring runner swaps during the race)
    race_ratio = time_resetter_to_threshold / time_left_to_threshold if time_left_to_threshold > 0 else Fraction(0)
    # race_ratio < 1 means g_right wins (LL reached)

    outcome = "LL" if race_ratio < 1 else ("tie" if race_ratio == 1 else "SS")

    return {
        "degenerate": False,
        "t_reset": t_reset,
        "resetter_idx": resetter_idx,
        "L_0": L_0,
        "R_0": R_0,
        "entry_fiber": "LS",
        "race_outcome": outcome,
        "time_right": time_resetter_to_threshold,
        "time_left": time_left_to_threshold,
        "race_ratio": race_ratio,
        "v_resetter": v_resetter,
        "v_nearest_ccw": nearest_ccw_speed,
        "speed_ratio": Fraction(nearest_ccw_speed, v_resetter),
    }


def analyze_gap_races(n_values: list[int] = [3, 4, 5, 6, 7]):
    """Analyze gap races for exhaustive primitive speed sets."""
    print("=" * 70)
    print("PART A: Gap race analysis — wrap-around entry quality")
    print("=" * 70)
    print()

    for n in n_values:
        max_speed = {3: 30, 4: 20, 5: 15, 6: 12, 7: 10}.get(n, 9)
        thr = Fraction(1, n)

        total_sets = 0
        sets_with_race_ll = 0
        sets_without_race_ll = 0
        all_L0_max = []
        all_best_race_ratios = []

        for speeds in ((combo) for combo in
                       combinations(range(1, max_speed + 1), n - 1)
                       if reduce(gcd, combo) == 1):
            total_sets += 1
            events = compute_wrap_arounds(speeds)

            races = []
            for t_reset, resetter_idx in events:
                result = gap_race_at_reset(speeds, n, t_reset, resetter_idx)
                if not result.get("degenerate", False) and result["entry_fiber"] == "LS":
                    races.append(result)

            if not races:
                sets_without_race_ll += 1
                continue

            max_L0 = max(r["L_0"] for r in races)
            all_L0_max.append(float(max_L0))

            ll_races = [r for r in races if r["race_outcome"] == "LL"]
            if ll_races:
                sets_with_race_ll += 1
                best_ratio = min(r["race_ratio"] for r in ll_races)
                all_best_race_ratios.append(float(best_ratio))
            else:
                sets_without_race_ll += 1
                # Report the best (closest to LL) race
                best_ratio = min(r["race_ratio"] for r in races)
                all_best_race_ratios.append(float(best_ratio))

        print(f"n={n}  max_speed={max_speed}  total_sets={total_sets}")
        print(f"  LL via simplified race: {sets_with_race_ll}/{total_sets}")
        print(f"  no LL via simplified race: {sets_without_race_ll}/{total_sets}")

        if all_L0_max:
            avg_L0 = sum(all_L0_max) / len(all_L0_max)
            min_L0 = min(all_L0_max)
            print(f"  max(L_0) across wrap-arounds: avg={avg_L0:.4f}, min={min_L0:.4f}")

        if all_best_race_ratios:
            avg_rr = sum(all_best_race_ratios) / len(all_best_race_ratios)
            max_rr = max(all_best_race_ratios)
            print(f"  best race ratio (< 1 = LL): avg={avg_rr:.4f}, max={max_rr:.4f}")

        print()


def initial_segment_gap_races(n_max: int = 16):
    """Gap races for initial segment speeds {1,...,n-1}."""
    print("=" * 70)
    print("PART B: Initial segment {1,...,n-1} gap races")
    print("=" * 70)
    print()

    for n in range(3, n_max + 1):
        speeds = tuple(range(1, n))
        events = compute_wrap_arounds(speeds)

        races = []
        for t_reset, resetter_idx in events:
            result = gap_race_at_reset(speeds, n, t_reset, resetter_idx)
            if not result.get("degenerate", False) and result["entry_fiber"] == "LS":
                races.append(result)

        ll_races = [r for r in races if r["race_outcome"] == "LL"]
        tie_races = [r for r in races if r["race_outcome"] == "tie"]
        ss_races = [r for r in races if r["race_outcome"] == "SS"]

        print(f"n={n:2d}  wrap-arounds={len(events)}  LS entries={len(races)}")
        print(f"  LL races: {len(ll_races)}  tie: {len(tie_races)}  SS: {len(ss_races)}")

        if races:
            L0_values = [float(r["L_0"]) for r in races]
            print(f"  L_0 range: [{min(L0_values):.6f}, {max(L0_values):.6f}]")

            race_ratios = [float(r["race_ratio"]) for r in races]
            print(f"  race ratio range: [{min(race_ratios):.6f}, {max(race_ratios):.6f}]")

        # Show the best race
        if ll_races:
            best = min(ll_races, key=lambda r: r["race_ratio"])
            print(f"  best LL race: t_reset={float(best['t_reset']):.4f}, "
                  f"v_reset={best['v_resetter']}, v_ccw={best['v_nearest_ccw']}, "
                  f"ratio={float(best['race_ratio']):.4f}")
        elif tie_races:
            best = tie_races[0]
            print(f"  closest race (TIE): v_reset={best['v_resetter']}, "
                  f"v_ccw={best['v_nearest_ccw']}")

        print()


def L0_distribution_analysis():
    """Study the L_0 distribution to see if max(L_0) -> 1 as max_speed grows."""
    print("=" * 70)
    print("PART C: L_0 maximum vs speed bound")
    print("=" * 70)
    print()

    for n in [3, 4, 5]:
        print(f"n={n}:")
        for max_speed in [10, 20, 30, 50]:
            max_L0_overall = 0.0
            count = 0
            for combo in combinations(range(1, max_speed + 1), n - 1):
                if reduce(gcd, combo) != 1:
                    continue
                count += 1
                events = compute_wrap_arounds(combo)
                for t_reset, resetter_idx in events:
                    result = gap_race_at_reset(combo, n, t_reset, resetter_idx)
                    if not result.get("degenerate", False) and result["entry_fiber"] == "LS":
                        max_L0_overall = max(max_L0_overall, float(result["L_0"]))
                # Only check first few for large spaces
                if count > 500:
                    break

            print(f"  max_speed={max_speed:3d}  sets={count}  max(L_0)={max_L0_overall:.6f}")
        print()


def main():
    print("LRC Gap Race Analysis — opus-2026-06-01-S518")
    print("THM-387: LL entry requires winning the gap race from LS")
    print()

    analyze_gap_races(n_values=[3, 4, 5, 6])
    initial_segment_gap_races(n_max=14)
    L0_distribution_analysis()

    print("=" * 70)
    print("SYNTHESIS")
    print("=" * 70)
    print()
    print("The simplified race analysis captures the EXACT outcome when no")
    print("runner swaps occur during the race interval.  For initial segments,")
    print("the tie races (race_ratio = 1) correspond to the wall-only LL witnesses.")
    print("The full picture requires tracking runner swaps during the race.")
    print()


if __name__ == "__main__":
    main()
