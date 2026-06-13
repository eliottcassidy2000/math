#!/usr/bin/env python3
"""n=6 inside debt as function of the 3-pair joint state.

opus-2026-06-01-S533

For n=6: the tiling model has 10 tiles with 3 independent tiles:
  (4,2), (5,3), (6,1)  — disjoint vertex sets {4,2}, {5,3}, {6,1}

These cover ALL 6 vertices. Each can be aligned (0) or anti-aligned (1).
The joint state (a,b,c) ∈ {0,1}³ gives 8 configurations.

Tile (6,1) is the APEX — the observer-to-sink arc.

The goal: compute the inside debt for each of the 8 joint states, and
find the mod-6 three-channel analogue of "a+b+c odd" from n=4.

For each primitive speed set, at each time t:
1. Record the joint state (a,b,c) of the 3 independent tiles
2. Record the observer outdegree (lonely status)
3. Aggregate: which (a,b,c) states have lonely times?

The "parity law" = a condition on (a,b,c) that determines loneliness.
"""

from __future__ import annotations

from fractions import Fraction
from math import gcd
from functools import reduce
from itertools import combinations
from collections import Counter, defaultdict


ONE = Fraction(1)
ZERO = Fraction(0)


def frac(x: Fraction) -> Fraction:
    return x - Fraction(x.numerator // x.denominator)


def dist0(x: Fraction) -> Fraction:
    f = frac(x)
    return min(f, ONE - f)


# ═══════════════════════════════════════════════════════════════
# The 3 independent tiles for n=6
# ═══════════════════════════════════════════════════════════════

# Base path: 6→5→4→3→2→1 (1-indexed, vertex 6 = observer for LRC)
# Tiles (non-adjacent base-path arcs):
# (6,4), (6,3), (6,2), (6,1), (5,3), (5,2), (5,1), (4,2), (4,1), (3,1)
#
# Independent tiles: (4,2), (5,3), (6,1)
# These use vertices {4,2}, {5,3}, {6,1} — all 6 covered.
#
# In LRC with observer=vertex 6, runners at vertices 1..5 with speeds v_1..v_5:
#   Tile (4,2): arc between runner_4 and runner_2 (half-turn)
#   Tile (5,3): arc between runner_5 and runner_3 (half-turn)
#   Tile (6,1): arc between observer and runner_1 (threshold 1/6)
#
# Tile aligned = "higher vertex beats lower" in base-path order.
# (4,2) aligned: runner_4 → runner_2 in half-turn
# (5,3) aligned: runner_5 → runner_3
# (6,1) aligned: observer → runner_1 (runner_1 is safe = lonely-compatible)

def compute_tile_state(speeds, n, t):
    """Compute the state of the 3 independent tiles at time t.

    Returns (a, b, c) where:
      a = 1 if tile (4,2) is aligned, else 0
          i.e., frac((v_4 - v_2) * t) ∈ (0, 1/2) [half-turn: runner 4 beats runner 2]
      b = 1 if tile (5,3) is aligned
          i.e., frac((v_5 - v_3) * t) ∈ (0, 1/2)
      c = 1 if tile (6,1) is aligned (observer beats runner 1)
          i.e., ||v_1 * t|| ≥ 1/n

    Note: using 0-indexed speeds. speeds[0]=v_1, ..., speeds[4]=v_5.
    In the base path: vertex 1 has speed v_1 (slowest), vertex 5 has v_5 (fastest).
    """
    thr = Fraction(1, n)

    # Tile (4,2): runner at speed v_4 vs runner at speed v_2
    # In 0-indexed speeds: v_4 = speeds[3], v_2 = speeds[1]
    diff_42 = frac(Fraction(speeds[3] - speeds[1]) * t)
    a = 1 if ZERO < diff_42 < Fraction(1, 2) else 0

    # Tile (5,3): runner at speed v_5 vs runner at speed v_3
    diff_53 = frac(Fraction(speeds[4] - speeds[2]) * t)
    b = 1 if ZERO < diff_53 < Fraction(1, 2) else 0

    # Tile (6,1): observer vs runner_1 = apex tile
    # Aligned = observer beats runner_1 = runner_1 is safe
    c = 1 if dist0(Fraction(speeds[0]) * t) >= thr else 0

    return (a, b, c)


def compute_lonely(speeds, n, t):
    """Check if time t is lonely."""
    thr = Fraction(1, n)
    return all(dist0(Fraction(v) * t) >= thr for v in speeds)


def compute_obs_outdeg(speeds, n, t):
    """Observer outdegree = number of safe runners."""
    thr = Fraction(1, n)
    return sum(1 for v in speeds if dist0(Fraction(v) * t) >= thr)


# ═══════════════════════════════════════════════════════════════
# PART 1: Joint state distribution for initial segment
# ═══════════════════════════════════════════════════════════════

def initial_segment_analysis():
    """For n=6, speeds=(1,2,3,4,5): compute joint state vs loneliness."""
    print("=" * 70)
    print("PART 1: n=6 initial segment — 3-channel joint state")
    print("=" * 70)
    print()

    n = 6
    speeds = (1, 2, 3, 4, 5)

    # Compute all walls
    thr = Fraction(1, n)
    walls = set()
    walls.add(ZERO)
    for v in speeds:
        for a in [1, n - 1]:
            for k in range(v):
                walls.add(Fraction(k * n + a, v * n))
        for k in range(v):
            walls.add(Fraction(k, v))
    for i, vi in enumerate(speeds):
        for j, vj in enumerate(speeds):
            if i >= j:
                continue
            diff = abs(vi - vj)
            for k in range(diff):
                walls.add(Fraction(k, diff))
                t = Fraction(2 * k + 1, 2 * diff)
                if ZERO <= t < ONE:
                    walls.add(t)

    walls = sorted(walls)
    walls_ext = walls + [ONE]

    # For each cell: record joint state and lonely status
    state_data = defaultdict(lambda: {"measure": Fraction(0), "lonely_measure": Fraction(0),
                                       "obs_outdeg": Counter(), "count": 0})

    for idx in range(len(walls)):
        t_start = walls_ext[idx]
        t_end = walls_ext[idx + 1]
        width = t_end - t_start
        if width <= 0:
            continue

        t_mid = (t_start + t_end) / 2
        state = compute_tile_state(speeds, n, t_mid)
        lonely = compute_lonely(speeds, n, t_mid)
        obs_od = compute_obs_outdeg(speeds, n, t_mid)

        state_data[state]["measure"] += width
        if lonely:
            state_data[state]["lonely_measure"] += width
        state_data[state]["obs_outdeg"][obs_od] += 1
        state_data[state]["count"] += 1

    # Also check wall points
    wall_lonely_states = Counter()
    for t in walls:
        if compute_lonely(speeds, n, t):
            state = compute_tile_state(speeds, n, t)
            wall_lonely_states[state] += 1

    # Report
    print(f"speeds = {speeds}, n = {n}")
    print(f"cells: {sum(d['count'] for d in state_data.values())}")
    print()

    print(f"{'state':>10} {'measure':>10} {'lonely_m':>10} {'cells':>6} {'outdeg_dist':>30} {'parity':>7}")
    print("-" * 85)

    for state in sorted(state_data.keys()):
        d = state_data[state]
        a, b, c = state
        parity = (a + b + c) % 2
        outdeg_str = str(dict(d["obs_outdeg"]))
        print(f"  ({a},{b},{c})  {float(d['measure']):10.6f} {float(d['lonely_measure']):10.6f} "
              f"{d['count']:6d} {outdeg_str:>30} {'odd' if parity else 'even':>7}")

    print()
    print("Wall lonely states:", dict(wall_lonely_states))
    print()

    # Check: which states have lonely measure > 0 or wall lonely?
    print("LONELY STATES:")
    for state in sorted(state_data.keys()):
        d = state_data[state]
        a, b, c = state
        has_lonely = d["lonely_measure"] > 0 or wall_lonely_states.get(state, 0) > 0
        parity = (a + b + c) % 2
        print(f"  ({a},{b},{c}) parity={'odd' if parity else 'even'}: "
              f"{'LONELY' if has_lonely else 'not lonely'}  "
              f"(open={float(d['lonely_measure']):.6f}, wall={wall_lonely_states.get(state, 0)})")

    print()


# ═══════════════════════════════════════════════════════════════
# PART 2: Exhaustive speed sets at n=6
# ═══════════════════════════════════════════════════════════════

def exhaustive_n6():
    """For all primitive 5-speed sets at n=6: which joint states have lonely times?"""
    print("=" * 70)
    print("PART 2: Exhaustive n=6 — joint state vs loneliness")
    print("=" * 70)
    print()

    n = 6
    max_speed = 12
    thr = Fraction(1, n)

    # For each speed set, find which joint states are visited and which are lonely
    state_has_lonely = defaultdict(int)  # state → # speed sets where it's lonely
    state_visited = defaultdict(int)  # state → # speed sets that visit it
    total_sets = 0

    for combo in combinations(range(1, max_speed + 1), 5):
        if reduce(gcd, combo) != 1:
            continue
        total_sets += 1
        speeds = combo

        # Quick wall computation
        walls = set()
        walls.add(ZERO)
        for v in speeds:
            for a in [1, n - 1]:
                for k in range(v):
                    walls.add(Fraction(k * n + a, v * n))
            for k in range(v):
                walls.add(Fraction(k, v))
        for i, vi in enumerate(speeds):
            for j, vj in enumerate(speeds):
                if i >= j:
                    continue
                diff = abs(vi - vj)
                for k in range(diff):
                    walls.add(Fraction(k, diff))
                    t = Fraction(2 * k + 1, 2 * diff)
                    if ZERO <= t < ONE:
                        walls.add(t)
        walls = sorted(walls)
        walls_ext = walls + [ONE]

        visited = set()
        lonely_states = set()

        # Check cells
        for idx in range(len(walls)):
            if walls_ext[idx + 1] <= walls_ext[idx]:
                continue
            t_mid = (walls_ext[idx] + walls_ext[idx + 1]) / 2
            state = compute_tile_state(speeds, n, t_mid)
            visited.add(state)
            if compute_lonely(speeds, n, t_mid):
                lonely_states.add(state)

        # Check walls
        for t in walls:
            state = compute_tile_state(speeds, n, t)
            visited.add(state)
            if compute_lonely(speeds, n, t):
                lonely_states.add(state)

        for s in visited:
            state_visited[s] += 1
        for s in lonely_states:
            state_has_lonely[s] += 1

    print(f"total primitive sets: {total_sets}")
    print()

    print(f"{'state':>10} {'visited':>8} {'lonely':>8} {'lonely%':>8} {'parity':>7} {'a+b+c':>6}")
    print("-" * 55)

    for state in sorted(set(state_visited.keys()) | set(state_has_lonely.keys())):
        a, b, c = state
        vis = state_visited[state]
        lon = state_has_lonely.get(state, 0)
        pct = 100 * lon / vis if vis > 0 else 0
        parity = (a + b + c) % 2
        print(f"  ({a},{b},{c})  {vis:8d} {lon:8d} {pct:7.1f}% {'odd' if parity else 'even':>7} {a+b+c:>6}")

    print()

    # Check if parity determines loneliness
    odd_lonely = sum(1 for s in state_has_lonely if (s[0]+s[1]+s[2]) % 2 == 1)
    even_lonely = sum(1 for s in state_has_lonely if (s[0]+s[1]+s[2]) % 2 == 0)
    odd_visited = sum(1 for s in state_visited if (s[0]+s[1]+s[2]) % 2 == 1)
    even_visited = sum(1 for s in state_visited if (s[0]+s[1]+s[2]) % 2 == 0)

    print(f"PARITY CHECK:")
    print(f"  odd (a+b+c≡1 mod 2): {odd_lonely}/{odd_visited} states have lonely times")
    print(f"  even (a+b+c≡0 mod 2): {even_lonely}/{even_visited} states have lonely times")
    print()

    # Check mod-3 and mod-6 patterns
    for modulus in [2, 3, 6]:
        print(f"  mod-{modulus} pattern of (a+b+c):")
        for r in range(modulus):
            states_in_r = [s for s in state_visited if (s[0]+s[1]+s[2]) % modulus == r]
            lonely_in_r = [s for s in state_has_lonely if (s[0]+s[1]+s[2]) % modulus == r]
            print(f"    (a+b+c)≡{r} mod {modulus}: {len(lonely_in_r)}/{len(states_in_r)} states lonely")
    print()

    # Look for more complex patterns
    print("CHANNEL-SPECIFIC PATTERNS:")
    for channel in ["a", "b", "c"]:
        idx = {"a": 0, "b": 1, "c": 2}[channel]
        for val in [0, 1]:
            states_with = [s for s in state_visited if s[idx] == val]
            lonely_with = [s for s in state_has_lonely if s[idx] == val]
            print(f"  {channel}={val}: {len(lonely_with)}/{len(states_with)} states lonely")
    print()

    # The APEX channel (c): c=1 means runner_1 is safe (apex aligned)
    # Does c=1 make loneliness more likely?
    print("APEX CHANNEL (c = observer beats runner_1):")
    c1_lonely = sum(state_has_lonely.get(s, 0) for s in state_visited if s[2] == 1)
    c1_total = sum(state_visited.get(s, 0) for s in state_visited if s[2] == 1)
    c0_lonely = sum(state_has_lonely.get(s, 0) for s in state_visited if s[2] == 0)
    c0_total = sum(state_visited.get(s, 0) for s in state_visited if s[2] == 0)
    print(f"  c=1 (apex aligned): {c1_lonely} lonely visits out of {c1_total}")
    print(f"  c=0 (apex anti): {c0_lonely} lonely visits out of {c0_total}")
    print()


# ═══════════════════════════════════════════════════════════════
# PART 3: The debt as a function of joint state
# ═══════════════════════════════════════════════════════════════

def debt_by_joint_state():
    """Compute the resonance debt for each (a,b,c) configuration.

    For each joint state: the 7 dependent tiles can vary.
    Compute the AVERAGE lonely measure across all speed sets
    that visit each joint state.
    """
    print("=" * 70)
    print("PART 3: Inside debt by joint state (averaged over speed sets)")
    print("=" * 70)
    print()

    n = 6
    max_speed = 12
    thr = Fraction(1, n)
    f0 = (n - 2) / n
    m = n - 1
    outside_credit = f0 ** m

    # For each state: accumulate lonely measures
    state_lonely_sum = defaultdict(float)
    state_measure_sum = defaultdict(float)
    state_count = defaultdict(int)

    for combo in combinations(range(1, max_speed + 1), 5):
        if reduce(gcd, combo) != 1:
            continue
        speeds = combo

        walls = set()
        walls.add(ZERO)
        for v in speeds:
            for a in [1, n - 1]:
                for k in range(v):
                    walls.add(Fraction(k * n + a, v * n))
            for k in range(v):
                walls.add(Fraction(k, v))
        for i, vi in enumerate(speeds):
            for j, vj in enumerate(speeds):
                if i >= j:
                    continue
                diff = abs(vi - vj)
                for k in range(diff):
                    walls.add(Fraction(k, diff))
                    t = Fraction(2 * k + 1, 2 * diff)
                    if ZERO <= t < ONE:
                        walls.add(t)
        walls = sorted(walls)
        walls_ext = walls + [ONE]

        # Total lonely measure for this speed set
        total_lonely = Fraction(0)
        state_lonely = defaultdict(lambda: Fraction(0))
        state_meas = defaultdict(lambda: Fraction(0))

        for idx in range(len(walls)):
            if walls_ext[idx + 1] <= walls_ext[idx]:
                continue
            t_mid = (walls_ext[idx] + walls_ext[idx + 1]) / 2
            width = walls_ext[idx + 1] - walls_ext[idx]

            state = compute_tile_state(speeds, n, t_mid)
            lonely = compute_lonely(speeds, n, t_mid)

            state_meas[state] += width
            if lonely:
                state_lonely[state] += width
                total_lonely += width

        for state in state_meas:
            state_measure_sum[state] += float(state_meas[state])
            state_lonely_sum[state] += float(state_lonely[state])
            state_count[state] += 1

    print(f"{'state':>10} {'avg_measure':>12} {'avg_lonely':>12} {'lonely_rate':>12} {'sets':>6} {'a+b+c%2':>8}")
    print("-" * 70)

    for state in sorted(state_count.keys()):
        a, b, c = state
        avg_m = state_measure_sum[state] / state_count[state]
        avg_l = state_lonely_sum[state] / state_count[state]
        rate = avg_l / avg_m if avg_m > 0 else 0
        parity = (a + b + c) % 2

        print(f"  ({a},{b},{c})  {avg_m:12.6f} {avg_l:12.6f} {rate:12.6f} "
              f"{state_count[state]:6d} {'odd' if parity else 'even':>8}")

    print()

    # DEBT ANALYSIS: lonely rate vs outside credit
    print("DEBT ANALYSIS:")
    print(f"  outside credit = ((n-2)/n)^{m} = {outside_credit:.6f}")
    print()
    for state in sorted(state_count.keys()):
        a, b, c = state
        avg_l = state_lonely_sum[state] / state_count[state]
        avg_m = state_measure_sum[state] / state_count[state]
        rate = avg_l / avg_m if avg_m > 0 else 0
        debt = outside_credit - rate
        print(f"  ({a},{b},{c}): lonely_rate={rate:.6f}, debt={debt:.6f}, "
              f"debt/credit={debt/outside_credit:.4f}")
    print()


# ═══════════════════════════════════════════════════════════════
# PART 4: Search for the parity law
# ═══════════════════════════════════════════════════════════════

def search_parity_law():
    """Try to find the mod-6 rule that determines loneliness from (a,b,c)."""
    print("=" * 70)
    print("PART 4: Search for the three-channel parity law")
    print("=" * 70)
    print()

    n = 6
    max_speed = 12
    thr = Fraction(1, n)

    # Collect: for each speed set, which states are lonely?
    lonely_profiles = Counter()  # frozenset of lonely states → count

    total_sets = 0
    for combo in combinations(range(1, max_speed + 1), 5):
        if reduce(gcd, combo) != 1:
            continue
        total_sets += 1
        speeds = combo

        walls = set()
        walls.add(ZERO)
        for v in speeds:
            for a in [1, n - 1]:
                for k in range(v):
                    walls.add(Fraction(k * n + a, v * n))
            for k in range(v):
                walls.add(Fraction(k, v))
        for i, vi in enumerate(speeds):
            for j, vj in enumerate(speeds):
                if i >= j:
                    continue
                diff = abs(vi - vj)
                for k in range(diff):
                    walls.add(Fraction(k, diff))
                    t = Fraction(2 * k + 1, 2 * diff)
                    if ZERO <= t < ONE:
                        walls.add(t)
        walls = sorted(walls)
        walls_ext = walls + [ONE]

        lonely_states = set()
        for idx in range(len(walls)):
            if walls_ext[idx + 1] <= walls_ext[idx]:
                continue
            t_mid = (walls_ext[idx] + walls_ext[idx + 1]) / 2
            if compute_lonely(speeds, n, t_mid):
                lonely_states.add(compute_tile_state(speeds, n, t_mid))
        for t in walls:
            if compute_lonely(speeds, n, t):
                lonely_states.add(compute_tile_state(speeds, n, t))

        lonely_profiles[frozenset(lonely_states)] += 1

    print(f"total speed sets: {total_sets}")
    print(f"distinct lonely profiles: {len(lonely_profiles)}")
    print()

    print("Most common lonely profiles:")
    for profile, count in lonely_profiles.most_common(10):
        states = sorted(profile)
        parities = [(s[0]+s[1]+s[2]) % 2 for s in states]
        print(f"  {count:4d} sets: lonely at states {states}")
        print(f"         parities: {parities}")
    print()

    # The UNIVERSAL lonely state: which (a,b,c) is lonely for ALL speed sets?
    all_lonely = None
    for profile, count in lonely_profiles.items():
        if all_lonely is None:
            all_lonely = set(profile)
        else:
            all_lonely &= set(profile)

    print(f"States lonely for ALL speed sets: {sorted(all_lonely) if all_lonely else 'NONE'}")
    print()

    # For each state: what fraction of speed sets is it lonely?
    state_lonely_frac = {}
    for state in [(a,b,c) for a in [0,1] for b in [0,1] for c in [0,1]]:
        count = sum(cnt for profile, cnt in lonely_profiles.items() if state in profile)
        state_lonely_frac[state] = count / total_sets
        a, b, c = state
        print(f"  ({a},{b},{c}) a+b+c={a+b+c}: lonely in {count}/{total_sets} = {100*count/total_sets:.1f}% of sets")
    print()

    # THE PARITY LAW SEARCH: find a boolean function f(a,b,c) such that
    # f=1 => lonely for all speed sets (or a high fraction)
    print("BOOLEAN FUNCTION SEARCH:")
    print("  Looking for f(a,b,c) → {lonely always, lonely sometimes, never lonely}")
    print()

    for desc, func in [
        ("a+b+c odd", lambda a,b,c: (a+b+c) % 2 == 1),
        ("a+b+c even", lambda a,b,c: (a+b+c) % 2 == 0),
        ("c=1 (apex)", lambda a,b,c: c == 1),
        ("a=b", lambda a,b,c: a == b),
        ("a≠b (channels differ)", lambda a,b,c: a != b),
        ("a+b+c ≥ 2", lambda a,b,c: a+b+c >= 2),
        ("a+b+c ≤ 1", lambda a,b,c: a+b+c <= 1),
        ("a·b·c = 1 (all aligned)", lambda a,b,c: a*b*c == 1),
        ("c=1 and a=b", lambda a,b,c: c==1 and a==b),
        ("c=1 and a≠b", lambda a,b,c: c==1 and a!=b),
        ("c=1 and a+b odd", lambda a,b,c: c==1 and (a+b)%2==1),
    ]:
        matching = [s for s in state_lonely_frac if func(*s)]
        if not matching:
            continue
        min_frac = min(state_lonely_frac[s] for s in matching)
        max_frac = max(state_lonely_frac[s] for s in matching)
        avg_frac = sum(state_lonely_frac[s] for s in matching) / len(matching)
        print(f"  {desc:25s}: {len(matching)} states, "
              f"lonely range [{min_frac:.1%}, {max_frac:.1%}], avg {avg_frac:.1%}")

    print()


def main():
    print("n=6 Three-Channel Parity Law — opus-S533")
    print()

    initial_segment_analysis()
    exhaustive_n6()
    debt_by_joint_state()
    search_parity_law()

    print("=" * 70)
    print("SYNTHESIS")
    print("=" * 70)
    print()


if __name__ == "__main__":
    main()
