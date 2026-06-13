#!/usr/bin/env python3
"""Two-gap dynamics for the LRC source-gap criterion (THM-384).

opus-2026-06-01-S518

THM-384 reduces LRC to: both observer-adjacent gaps >= 1/n simultaneously.
This script studies the DYNAMICS of the two-gap pair (g_left, g_right) to
understand WHY the (long,long) state must always be visited.

Key questions:
1. What is the time-average of g_left + g_right?  Is it always >= 2/n?
2. How does the two-gap trajectory traverse the four fibers
   (SS, SL, LS, LL) where S = short (<1/n), L = long (>=1/n)?
3. What is the "exchange mechanism" — when one gap shrinks below 1/n,
   does the other grow?
4. Can we prove a pigeonhole/measure bound forcing a visit to LL?

Approach:
- For each primitive speed set, compute ALL wall times (rational points
  where the circular order changes or a runner crosses the 1/n threshold)
- Track (g_left, g_right) through each cell
- Compute the measure of time spent in each fiber
- Look for structural invariants

Tournament Analysis declaration:
    vertices: representative primitive speed sets at fixed n
    pairwise observable: open LL measure, then exchange ratio, then gap-sum integral
    switch/gauge: larger open LL measure wins; ties use exchange/gap-sum
    tie Hamiltonian path: lexicographic order on speed sets

Stored output:
    05-knowledge/results/lrc_two_gap_dynamics_s518.out
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations
from math import gcd
from collections import Counter, defaultdict
from functools import reduce
from typing import NamedTuple


# ─── Utilities ───────────────────────────────────────────────────────────

ONE = Fraction(1)
ZERO = Fraction(0)


def frac(x: Fraction) -> Fraction:
    """Fractional part in [0, 1)."""
    return x - Fraction(x.numerator // x.denominator)


def dist0(x: Fraction) -> Fraction:
    """Circular distance from 0: min({x}, 1-{x})."""
    f = frac(x)
    return min(f, ONE - f)


def is_primitive(speeds: tuple[int, ...]) -> bool:
    """Speed set is primitive iff gcd of all speeds is 1."""
    return reduce(gcd, speeds) == 1


def primitive_speed_sets(n: int, max_speed: int | None = None):
    """Generate primitive speed sets with n-1 distinct positive speeds.

    For small n, exhaustive up to max_speed.
    """
    k = n - 1  # number of runners
    if max_speed is None:
        max_speed = 3 * n  # reasonable bound for exploration
    for combo in combinations(range(1, max_speed + 1), k):
        if reduce(gcd, combo) == 1:
            yield combo


# ─── Wall computation ────────────────────────────────────────────────────

class Cell(NamedTuple):
    """A cell in the runner clock with its gap data."""
    t_start: Fraction
    t_end: Fraction
    g_left: Fraction   # counterclockwise gap at midpoint
    g_right: Fraction  # clockwise gap at midpoint
    fiber: str         # "LL", "LS", "SL", "SS"


def compute_walls(speeds: tuple[int, ...], n: int) -> list[Fraction]:
    """Compute all wall times in [0, 1) where circular structure changes.

    Walls come from:
    1. LRC threshold crossings: ||v_i t|| = 1/n  =>  {v_i t} = 1/n or (n-1)/n
    2. Runner-observer ties: {v_i t} = 0  =>  t = k/v_i
    3. Runner-runner order crossings (needed for gap tracking)

    For the two-gap problem, we need walls where:
    - A runner enters/exits the clockwise forbidden zone [0, 1/n)
    - A runner enters/exits the counterclockwise forbidden zone (1-1/n, 1)
    - The identity of the nearest clockwise/counterclockwise runner changes
    """
    thr = Fraction(1, n)
    walls = set()

    for v in speeds:
        # LRC threshold walls: {v*t} = 1/n or (n-1)/n
        # {v*t} = a/n  =>  v*t ≡ a/n (mod 1)  =>  t = (k*n + a)/(v*n)
        for a in [1, n - 1]:
            for k in range(v):
                t = Fraction(k * n + a, v * n)
                if ZERO <= t < ONE:
                    walls.add(t)

        # Observer ties: {v*t} = 0  =>  t = k/v
        for k in range(v):
            t = Fraction(k, v)
            if ZERO <= t < ONE:
                walls.add(t)

    # Runner-runner order crossings (when two runners swap circular position)
    for i, vi in enumerate(speeds):
        for j, vj in enumerate(speeds):
            if i >= j:
                continue
            diff = vi - vj
            if diff == 0:
                continue
            # {vi*t} = {vj*t}  =>  (vi-vj)*t ≡ 0 (mod 1)  =>  t = k/|diff|
            d = abs(diff)
            for k in range(d):
                t = Fraction(k, d)
                if ZERO <= t < ONE:
                    walls.add(t)

    walls.add(ZERO)
    return sorted(walls)


def compute_gaps_at(speeds: tuple[int, ...], n: int, t: Fraction) -> tuple[Fraction, Fraction]:
    """Compute (g_left, g_right) at time t.

    g_right = min fractional part (clockwise distance from observer)
    g_left = min (1 - fractional part) (counterclockwise distance from observer)

    If any runner is tied with observer, that gap is 0.
    """
    positions = [frac(Fraction(v) * t) for v in speeds]

    # Clockwise gap: smallest positive fractional part
    # If any position is 0, clockwise gap is 0
    cw_distances = []
    ccw_distances = []
    for p in positions:
        if p == ZERO:
            return (ZERO, ZERO)  # tied with observer
        cw_distances.append(p)        # clockwise distance
        ccw_distances.append(ONE - p) # counterclockwise distance

    g_right = min(cw_distances)
    g_left = min(ccw_distances)
    return (g_left, g_right)


def fiber_label(g_left: Fraction, g_right: Fraction, thr: Fraction) -> str:
    """Label the two-gap fiber."""
    left_long = g_left >= thr
    right_long = g_right >= thr
    return ("L" if left_long else "S") + ("L" if right_long else "S")


def analyze_speed_set(speeds: tuple[int, ...], n: int) -> dict:
    """Full two-gap analysis for a speed set."""
    thr = Fraction(1, n)
    walls = compute_walls(speeds, n)

    # Add a sentinel for periodicity
    walls_extended = walls + [ONE]

    cells = []
    fiber_measures = Counter()
    total_g_sum = ZERO  # integral of g_left + g_right
    max_g_sum = ZERO
    min_g_sum = ONE  # will be updated
    lonely_measure = ZERO

    for idx in range(len(walls)):
        t_start = walls_extended[idx]
        t_end = walls_extended[idx + 1]
        width = t_end - t_start

        if width <= 0:
            continue

        # Evaluate at midpoint
        t_mid = (t_start + t_end) / 2
        g_left, g_right = compute_gaps_at(speeds, n, t_mid)
        fib = fiber_label(g_left, g_right, thr)

        cells.append(Cell(t_start, t_end, g_left, g_right, fib))
        fiber_measures[fib] += width

        g_sum = g_left + g_right
        total_g_sum += g_sum * width
        max_g_sum = max(max_g_sum, g_sum)
        min_g_sum = min(min_g_sum, g_sum)

        if fib == "LL":
            lonely_measure += width

    # Check wall values too for boundary witnesses
    wall_ll_count = 0
    for t in walls:
        g_left, g_right = compute_gaps_at(speeds, n, t)
        if fiber_label(g_left, g_right, thr) == "LL":
            wall_ll_count += 1

    return {
        "speeds": speeds,
        "n": n,
        "num_cells": len(cells),
        "num_walls": len(walls),
        "fiber_measures": dict(fiber_measures),
        "lonely_measure": lonely_measure,
        "avg_g_sum": total_g_sum,  # this is the integral, not avg
        "max_g_sum": max_g_sum,
        "min_g_sum": min_g_sum,
        "wall_ll_count": wall_ll_count,
        "cells": cells,
        "has_open_ll": lonely_measure > 0,
        "has_wall_ll": wall_ll_count > 0,
        "has_any_ll": lonely_measure > 0 or wall_ll_count > 0,
    }


def gap_sum_integral(speeds: tuple[int, ...], n: int) -> Fraction:
    """Integral of g_left + g_right without storing the full cell movie."""
    walls = compute_walls(speeds, n)
    walls_extended = walls + [ONE]
    total_g_sum = ZERO

    for idx in range(len(walls)):
        t_start = walls_extended[idx]
        t_end = walls_extended[idx + 1]
        width = t_end - t_start
        if width <= 0:
            continue

        t_mid = (t_start + t_end) / 2
        g_left, g_right = compute_gaps_at(speeds, n, t_mid)
        total_g_sum += (g_left + g_right) * width

    return total_g_sum


# ─── Fiber transition analysis ──────────────────────────────────────────

def fiber_transitions(cells: list[Cell]) -> Counter:
    """Count transitions between fibers in the cell sequence."""
    transitions = Counter()
    for i in range(len(cells) - 1):
        transitions[(cells[i].fiber, cells[i + 1].fiber)] += 1
    # Wraparound
    if len(cells) >= 2:
        transitions[(cells[-1].fiber, cells[0].fiber)] += 1
    return transitions


def gap_exchange_analysis(cells: list[Cell], thr: Fraction) -> dict:
    """Analyze the exchange mechanism between g_left and g_right.

    Key question: when g_left drops below 1/n, does g_right tend to increase?
    """
    exchange_events = 0
    anti_exchange_events = 0

    for i in range(len(cells) - 1):
        c1, c2 = cells[i], cells[i + 1]
        # Left gap change
        dl = c2.g_left - c1.g_left
        # Right gap change
        dr = c2.g_right - c1.g_right

        if dl != 0 and dr != 0:
            if (dl > 0 and dr < 0) or (dl < 0 and dr > 0):
                exchange_events += 1  # anti-correlated
            else:
                anti_exchange_events += 1  # co-moving

    return {
        "exchange_events": exchange_events,
        "co_moving_events": anti_exchange_events,
        "exchange_ratio": (
            exchange_events / (exchange_events + anti_exchange_events)
            if (exchange_events + anti_exchange_events) > 0
            else None
        ),
    }


# ─── Initial segment analysis ───────────────────────────────────────────

def initial_segment_analysis(n_max: int = 20):
    """Analyze {1,...,n-1} speeds for n = 3 to n_max."""
    print("=" * 70)
    print("PART A: Initial segment speeds {1,...,n-1}")
    print("=" * 70)
    print()

    for n in range(3, n_max + 1):
        speeds = tuple(range(1, n))
        result = analyze_speed_set(speeds, n)

        thr = Fraction(1, n)
        # Find the exact lonely time and gap values
        lonely_times = []
        for cell in result["cells"]:
            if cell.fiber == "LL":
                lonely_times.append((cell.t_start, cell.t_end, cell.g_left, cell.g_right))

        # Also check wall times
        walls = compute_walls(speeds, n)
        wall_lonely = []
        for t in walls:
            gl, gr = compute_gaps_at(speeds, n, t)
            if fiber_label(gl, gr, thr) == "LL":
                wall_lonely.append((t, gl, gr))

        print(f"n={n:2d}  speeds={speeds}")
        print(f"  cells={result['num_cells']:5d}  walls={result['num_walls']:5d}")
        print(f"  fiber measures: {', '.join(f'{k}={float(v):.6f}' for k,v in sorted(result['fiber_measures'].items()))}")
        print(f"  lonely measure (LL): {float(result['lonely_measure']):.8f}")
        print(f"  avg(g_left+g_right) = {float(result['avg_g_sum']):.6f}  (fair share = {float(Fraction(2,n)):.6f})")
        print(f"  max(g_left+g_right) = {float(result['max_g_sum']):.6f}")
        print(f"  min(g_left+g_right) = {float(result['min_g_sum']):.6f}")

        if lonely_times:
            print(f"  open LL cells: {len(lonely_times)}")
            for ts, te, gl, gr in lonely_times[:3]:
                print(f"    [{float(ts):.6f}, {float(te):.6f}]  g_left={float(gl):.6f}  g_right={float(gr):.6f}")

        if wall_lonely:
            print(f"  wall LL points: {len(wall_lonely)}")
            for t, gl, gr in wall_lonely[:3]:
                print(f"    t={float(t):.6f}  g_left={float(gl):.6f}  g_right={float(gr):.6f}")

        if not lonely_times and not wall_lonely:
            print(f"  *** NO LL STATE FOUND — THIS WOULD BE A COUNTEREXAMPLE ***")

        print()

        # Stop if computation gets too heavy
        if n >= 12 and result["num_walls"] > 50000:
            print(f"  (stopping initial segment at n={n} due to wall count)")
            break


# ─── Exhaustive small-n analysis ─────────────────────────────────────────

def exhaustive_small_n(n_values: list[int] = [3, 4, 5, 6]):
    """Exhaustive analysis of all primitive speed sets at small n."""
    print("=" * 70)
    print("PART B: Exhaustive primitive speed sets")
    print("=" * 70)
    print()

    for n in n_values:
        k = n - 1
        max_speed = {3: 20, 4: 14, 5: 11, 6: 9, 7: 8, 8: 8}.get(n, 8)

        total = 0
        ll_count = 0
        no_ll_count = 0
        min_lonely_meas = ONE
        max_lonely_meas = ZERO
        fiber_stats = Counter()
        avg_g_sum_total = ZERO
        exchange_ratios = []

        worst_set = None

        for speeds in primitive_speed_sets(n, max_speed):
            result = analyze_speed_set(speeds, n)
            total += 1

            if result["has_any_ll"]:
                ll_count += 1
            else:
                no_ll_count += 1
                print(f"  *** NO LL: n={n} speeds={speeds} ***")

            if result["lonely_measure"] < min_lonely_meas and result["has_any_ll"]:
                min_lonely_meas = result["lonely_measure"]
                worst_set = speeds
            max_lonely_meas = max(max_lonely_meas, result["lonely_measure"])

            for fib, meas in result["fiber_measures"].items():
                fiber_stats[fib] += meas

            avg_g_sum_total += result["avg_g_sum"]

            # Exchange analysis
            thr = Fraction(1, n)
            exch = gap_exchange_analysis(result["cells"], thr)
            if exch["exchange_ratio"] is not None:
                exchange_ratios.append(exch["exchange_ratio"])

        print(f"n={n}  max_speed={max_speed}  total_sets={total}")
        print(f"  LL found: {ll_count}/{total}  ({100*ll_count/total:.1f}%)")
        if no_ll_count > 0:
            print(f"  *** {no_ll_count} sets with NO LL state — possible counterexamples! ***")
        else:
            print(f"  ALL sets visit LL — LRC holds for all tested.")

        if total > 0:
            avg_lonely = float(fiber_stats.get("LL", 0)) / total
            avg_g_sum = float(avg_g_sum_total) / total
            fair_share = float(Fraction(2, n))
            print(f"  avg LL measure: {avg_lonely:.6f}")
            print(f"  avg g_sum integral: {avg_g_sum:.6f}  (fair share: {fair_share:.6f})")
            print(f"  min lonely measure (open): {float(min_lonely_meas):.8f}  set={worst_set}")
            print(f"  max lonely measure: {float(max_lonely_meas):.8f}")

            if exchange_ratios:
                avg_exch = sum(exchange_ratios) / len(exchange_ratios)
                print(f"  avg exchange ratio: {avg_exch:.4f}  (1.0 = perfect anti-correlation)")

        print()


# ─── Gap-sum integral theory ─────────────────────────────────────────────

def gap_sum_integral_theory(n_max: int = 12):
    """Test the conjecture: integral of (g_left + g_right) >= 2/n for all
    primitive speed sets.

    If true, combined with continuity, this would force a visit to the
    (long,long) state.
    """
    print("=" * 70)
    print("PART C: Gap-sum integral >= 2/n conjecture")
    print("=" * 70)
    print()

    for n in range(3, n_max + 1):
        k = n - 1
        max_speed = {3: 28, 4: 18, 5: 13, 6: 11, 7: 9, 8: 8, 9: 8, 10: 7, 11: 7, 12: 7}.get(n, 7)
        fair_share = Fraction(2, n)

        total = 0
        violations = 0
        min_ratio = Fraction(10)  # will be updated
        min_ratio_set = None
        max_ratio = ZERO

        for speeds in primitive_speed_sets(n, max_speed):
            total += 1

            integral = gap_sum_integral(speeds, n)
            ratio = integral / fair_share if fair_share > 0 else Fraction(10)

            if ratio < min_ratio:
                min_ratio = ratio
                min_ratio_set = speeds
            max_ratio = max(max_ratio, ratio)

            if integral < fair_share:
                violations += 1

        print(f"n={n:2d}  sets={total}  violations={violations}")
        print(f"  min ratio (integral / fair_share): {float(min_ratio):.6f}  set={min_ratio_set}")
        print(f"  max ratio: {float(max_ratio):.6f}")

        if violations > 0:
            print(f"  *** CONJECTURE VIOLATED: {violations} sets have integral < 2/n ***")
        else:
            print(f"  integral >= 2/n holds for all tested sets.")
        print()

    print("Targeted exact gap-sum counterexamples outside the bounded scan:")
    for n, speeds in [
        (5, (5, 11, 12, 17)),
        (6, (5, 7, 8, 9, 13)),
    ]:
        result = analyze_speed_set(speeds, n)
        fair_share = Fraction(2, n)
        ratio = result["avg_g_sum"] / fair_share
        print(
            f"  n={n} speeds={speeds}: integral={float(result['avg_g_sum']):.6f} "
            f"fair={float(fair_share):.6f} ratio={float(ratio):.6f} "
            f"LL={float(result['lonely_measure']):.8f} wallLL={result['wall_ll_count']}"
        )
    print()


# ─── n=14 sanity check ──────────────────────────────────────────────────

def n14_periodic_sanity_check():
    """Two-gap dynamics for the n=14 initial row and its doubled traversal.

    This is not the full S514 hard-row ladder.  It only checks that doubling
    all speeds preserves the initial source-gap dynamics up to a twofold
    traversal of the circle.
    """
    print("=" * 70)
    print("PART D: n=14 initial/doubled periodic sanity check")
    print("=" * 70)
    print()

    n = 14
    thr = Fraction(1, n)

    # Initial segment
    rows = [
        ("initial", tuple(range(1, 14))),
        ("doubled initial (2x traversal)", tuple(2 * i for i in range(1, 14))),
    ]

    for name, speeds in rows:
        # Check primitivity
        g = reduce(gcd, speeds)
        if g > 1:
            print(f"  {name}: gcd={g}, reducing to primitive form")
            speeds = tuple(s // g for s in speeds)
            n_eff = n  # threshold stays 1/n for the original problem
        else:
            n_eff = n

        print(f"\n--- {name}: speeds={speeds[:6]}... (n={n_eff}) ---")

        # Limit wall computation for large speed sets
        max_wall_speed = max(speeds)
        est_walls = sum(2 * n_eff + v for v in speeds) + sum(abs(vi - vj) for vi in speeds for vj in speeds if vi > vj)
        print(f"  estimated wall count: ~{est_walls}")

        if est_walls > 100000:
            print(f"  (too many walls, using sampling instead)")
            # Sample approach
            num_samples = 10000
            ll_hits = 0
            g_sum_total = 0.0
            max_g_sum = 0.0
            min_g_sum = 1.0

            for s in range(num_samples):
                t = Fraction(s, num_samples)
                positions = [frac(Fraction(v) * t) for v in speeds]
                if any(p == ZERO for p in positions):
                    continue

                cw = min(positions)
                ccw = min(ONE - p for p in positions)
                g_sum = float(cw + ccw)
                g_sum_total += g_sum
                max_g_sum = max(max_g_sum, g_sum)
                min_g_sum = min(min_g_sum, g_sum)

                if cw >= thr and ccw >= thr:
                    ll_hits += 1

            print(f"  sampled LL hits: {ll_hits}/{num_samples} ({100*ll_hits/num_samples:.2f}%)")
            print(f"  avg g_sum: {g_sum_total/num_samples:.6f}  (fair: {float(Fraction(2,n_eff)):.6f})")
            print(f"  max g_sum: {max_g_sum:.6f}  min g_sum: {min_g_sum:.6f}")
        else:
            result = analyze_speed_set(speeds, n_eff)
            print(f"  cells={result['num_cells']}  walls={result['num_walls']}")
            print(f"  fibers: {', '.join(f'{k}={float(v):.6f}' for k,v in sorted(result['fiber_measures'].items()))}")
            print(f"  lonely measure: {float(result['lonely_measure']):.8f}")
            print(f"  avg g_sum: {float(result['avg_g_sum']):.6f}  (fair: {float(Fraction(2,n_eff)):.6f})")


# ─── Tournament Analysis summary ────────────────────────────────────────

def hamiltonian_path_count(adj: tuple[tuple[int, ...], ...]) -> int:
    n = len(adj)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for start in range(n):
        dp[1 << start][start] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            for nxt in range(n):
                if not ((mask >> nxt) & 1) and adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += count
    return sum(dp[full])


def directed_triangle_count(adj: tuple[tuple[int, ...], ...]) -> int:
    total = 0
    n = len(adj)
    for a, b, c in combinations(range(n), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            total += 1
    return total


def scc_sizes(adj: tuple[tuple[int, ...], ...]) -> tuple[int, ...]:
    n = len(adj)

    def reach(start: int, reverse: bool = False) -> set[int]:
        seen = {start}
        stack = [start]
        while stack:
            cur = stack.pop()
            for nxt in range(n):
                has_edge = adj[nxt][cur] if reverse else adj[cur][nxt]
                if has_edge and nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        return seen

    remaining = set(range(n))
    sizes = []
    while remaining:
        start = min(remaining)
        comp = reach(start) & reach(start, reverse=True)
        sizes.append(len(comp))
        remaining -= comp
    return tuple(sorted(sizes, reverse=True))


def score_histogram(adj: tuple[tuple[int, ...], ...]) -> str:
    counts = Counter(sum(row) for row in adj)
    return " ".join(f"{score}:{counts[score]}" for score in sorted(counts))


def exchange_ratio(result: dict) -> float:
    data = gap_exchange_analysis(result["cells"], Fraction(1, result["n"]))
    return float(data["exchange_ratio"] or 0.0)


def choose_representatives(n: int, max_speed: int) -> list[tuple[int, ...]]:
    """Choose speed sets that expose different two-gap dynamics."""
    records = []
    fair_share = Fraction(2, n)
    for speeds in primitive_speed_sets(n, max_speed):
        result = analyze_speed_set(speeds, n)
        records.append(
            {
                "speeds": speeds,
                "result": result,
                "ll": result["lonely_measure"],
                "gap_ratio": result["avg_g_sum"] / fair_share,
                "exchange": exchange_ratio(result),
            }
        )

    selectors = [
        min(records, key=lambda r: (r["ll"], r["speeds"])),
        max(records, key=lambda r: (r["ll"], tuple(-x for x in r["speeds"]))),
        min(records, key=lambda r: (r["gap_ratio"], r["speeds"])),
        max(records, key=lambda r: (r["gap_ratio"], tuple(-x for x in r["speeds"]))),
        min(records, key=lambda r: (r["exchange"], r["speeds"])),
        max(records, key=lambda r: (r["exchange"], tuple(-x for x in r["speeds"]))),
    ]
    violators = [r for r in records if r["gap_ratio"] < 1]
    if violators:
        selectors.append(min(violators, key=lambda r: (r["gap_ratio"], r["speeds"])))

    out: list[tuple[int, ...]] = []
    for record in selectors:
        speeds = record["speeds"]
        if speeds not in out:
            out.append(speeds)
    return out


def two_gap_tournament(n: int, max_speed: int) -> None:
    """Tournament Analysis on representative speed sets.

    Pairwise observable: `(LL measure, exchange ratio, gap-sum integral)`.
    Switch: larger open LL measure wins; then larger exchange ratio; then
    larger gap-sum integral.  Remaining ties follow lexicographic order.
    """
    vertices = choose_representatives(n, max_speed)
    metrics = {speeds: analyze_speed_set(speeds, n) for speeds in vertices}
    exch = {speeds: exchange_ratio(metrics[speeds]) for speeds in vertices}

    adj = [[0] * len(vertices) for _ in vertices]
    for i, j in combinations(range(len(vertices)), 2):
        a, b = vertices[i], vertices[j]
        key_a = (metrics[a]["lonely_measure"], Fraction.from_float(exch[a]), metrics[a]["avg_g_sum"])
        key_b = (metrics[b]["lonely_measure"], Fraction.from_float(exch[b]), metrics[b]["avg_g_sum"])
        if key_a > key_b:
            winner = i
        elif key_b > key_a:
            winner = j
        else:
            winner = i if a < b else j
        loser = j if winner == i else i
        adj[winner][loser] = 1

    adj_t = tuple(tuple(row) for row in adj)
    print("=" * 70)
    print(f"PART G: Tournament Analysis of representative n={n} speed sets")
    print("=" * 70)
    print()
    print("Pairwise observable: open LL measure, exchange ratio, gap-sum integral.")
    print("Switch: larger open LL measure wins; ties use exchange/gap-sum, then lexicographic path.")
    print()
    print(f"vertices={len(vertices)}  H={hamiltonian_path_count(adj_t)}  "
          f"c3={directed_triangle_count(adj_t)}  SCCs={scc_sizes(adj_t)}  "
          f"score=[{score_histogram(adj_t)}]")
    print("ranked vertices by score:")
    scores = [sum(row) for row in adj_t]
    for idx in sorted(range(len(vertices)), key=lambda k: (-scores[k], vertices[k])):
        speeds = vertices[idx]
        result = metrics[speeds]
        print(
            f"  score={scores[idx]:>2}  speeds={speeds}  "
            f"LL={float(result['lonely_measure']):.8f}  "
            f"wallLL={result['wall_ll_count']:>2}  "
            f"gapInt={float(result['avg_g_sum']):.6f}  "
            f"exchange={exch[speeds]:.4f}"
        )
    print()


# ─── Fiber sequence structure ────────────────────────────────────────────

def fiber_sequence_structure(n_values: list[int] = [3, 4, 5, 6]):
    """Analyze the structure of the fiber transition sequence.

    The THM-384 walk visits fibers SS, SL, LS, LL in sequence as t goes 0->1.
    Prove: the sequence must contain LL.

    Key structural question: can we have a cycle SS->SL->SS->SL->... that
    never hits LL? What constraints prevent this?
    """
    print("=" * 70)
    print("PART E: Fiber transition structure")
    print("=" * 70)
    print()

    for n in n_values:
        max_speed = {3: 20, 4: 14, 5: 11, 6: 9}.get(n, 9)

        print(f"n={n}:")
        transition_totals = Counter()
        fiber_seq_lengths = Counter()
        total = 0

        for speeds in primitive_speed_sets(n, max_speed):
            result = analyze_speed_set(speeds, n)
            total += 1

            # Get fiber sequence
            fiber_seq = [c.fiber for c in result["cells"]]

            # Count unique fiber types visited
            unique_fibers = set(fiber_seq)
            fiber_seq_lengths[len(unique_fibers)] += 1

            # Count transitions
            trans = fiber_transitions(result["cells"])
            for key, count in trans.items():
                transition_totals[key] += count

        # Report transitions
        print(f"  total sets: {total}")
        print(f"  unique fiber counts: {dict(fiber_seq_lengths)}")
        print(f"  transitions (total across all sets):")
        for (f1, f2), count in sorted(transition_totals.items()):
            print(f"    {f1} -> {f2}: {count}")

        # Key question: are there forbidden transitions?
        all_possible = [(a, b) for a in ["SS", "SL", "LS", "LL"] for b in ["SS", "SL", "LS", "LL"]]
        missing = [(a, b) for a, b in all_possible if transition_totals.get((a, b), 0) == 0]
        if missing:
            print(f"  MISSING transitions: {missing}")
        print()


# ─── Pigeonhole bound ────────────────────────────────────────────────────

def pigeonhole_bound(n_values: list[int] = [3, 4, 5, 6, 7, 8]):
    """Test a pigeonhole argument for why LL must be visited.

    Argument sketch:
    - The total measure of "clockwise close" times ∪_i {t : {v_i t} < 1/n}
      is at most (n-1)/n < 1 (each runner contributes measure 1/n).
    - So the complement (all runners clockwise-far) has measure >= 1/n.
    - Similarly for counterclockwise.
    - The question is whether these two complements OVERLAP.

    We compute:
    - mu(B+) = measure of ∪_i {t : {v_i t} < 1/n}
    - mu(B-) = measure of ∪_i {t : 1-{v_i t} < 1/n}
    - mu(B+ ∩ B-) = measure of times where both a CW-close and CCW-close runner exist
    - mu(B+^c ∩ B-^c) = mu(LL) = lonely measure

    Key identity: mu(LL) = 1 - mu(B+) - mu(B-) + mu(B+ ∩ B-)

    So mu(LL) >= 0 iff mu(B+ ∩ B-) >= mu(B+) + mu(B-) - 1.
    We need mu(LL) > 0 (or = 0 with wall LL).
    """
    print("=" * 70)
    print("PART F: Pigeonhole measure decomposition")
    print("=" * 70)
    print()

    for n in n_values:
        max_speed = {3: 20, 4: 14, 5: 11, 6: 9, 7: 8, 8: 8}.get(n, 8)
        thr = Fraction(1, n)

        total = 0
        min_overlap = ONE
        min_overlap_set = None

        print(f"n={n}:")

        for speeds in primitive_speed_sets(n, max_speed):
            walls = compute_walls(speeds, n)
            walls_ext = walls + [ONE]

            mu_bp = ZERO  # measure of B+ (some runner CW-close)
            mu_bm = ZERO  # measure of B- (some runner CCW-close)
            mu_both = ZERO  # measure of B+ ∩ B- (both present)
            mu_ll = ZERO  # measure of LL (neither present)

            for idx in range(len(walls)):
                t_start = walls_ext[idx]
                t_end = walls_ext[idx + 1]
                width = t_end - t_start
                if width <= 0:
                    continue

                t_mid = (t_start + t_end) / 2
                positions = [frac(Fraction(v) * t_mid) for v in speeds]

                has_cw_close = any(p < thr and p > ZERO for p in positions) or any(p == ZERO for p in positions)
                has_ccw_close = any(ONE - p < thr and p > ZERO for p in positions) or any(p == ZERO for p in positions)

                if has_cw_close:
                    mu_bp += width
                if has_ccw_close:
                    mu_bm += width
                if has_cw_close and has_ccw_close:
                    mu_both += width
                if not has_cw_close and not has_ccw_close:
                    mu_ll += width

            total += 1

            # Verify identity: mu_ll = 1 - mu_bp - mu_bm + mu_both
            check = ONE - mu_bp - mu_bm + mu_both
            assert abs(check - mu_ll) < Fraction(1, 10**10), f"Identity failed: {check} != {mu_ll}"

            # The "overlap" is the key: mu_both
            # For LL > 0, need mu_both > mu_bp + mu_bm - 1
            overlap_excess = mu_both - (mu_bp + mu_bm - ONE)

            if overlap_excess < min_overlap:
                min_overlap = overlap_excess
                min_overlap_set = speeds

        print(f"  total sets: {total}")
        print(f"  min overlap excess (mu_both - (mu_bp + mu_bm - 1)): {float(min_overlap):.8f}")
        print(f"    at speeds: {min_overlap_set}")
        print(f"  NOTE: overlap excess = mu(LL), so this should always be >= 0")
        print()


# ─── Main ────────────────────────────────────────────────────────────────

def main():
    print("LRC Two-Gap Dynamics — opus-2026-06-01-S518")
    print("THM-384 reduction: lonely iff g_left >= 1/n AND g_right >= 1/n")
    print()

    # Part A: Initial segment
    initial_segment_analysis(n_max=14)

    # Part B: Exhaustive small n
    exhaustive_small_n(n_values=[3, 4, 5, 6])

    # Part C: Gap-sum integral conjecture
    gap_sum_integral_theory(n_max=8)

    # Part D: n=14 initial/doubled sanity check
    n14_periodic_sanity_check()

    # Part E: Fiber transition structure
    fiber_sequence_structure(n_values=[3, 4, 5])

    # Part F: Pigeonhole decomposition
    pigeonhole_bound(n_values=[3, 4, 5, 6])

    # Part G: Tournament Analysis of representative rows
    two_gap_tournament(n=5, max_speed=11)

    print("=" * 70)
    print("SYNTHESIS")
    print("=" * 70)
    print()
    print("1. HYP-1981/THM-381 turns LRC into source reachability; THM-384")
    print("   reduces that source condition to the two observer-adjacent gaps.")
    print("2. In bounded exhaustive windows for n=3..6, every primitive speed set")
    print("   visits the closed LL source-gap fiber.  Initial segments are the")
    print("   extremal wall-only cases: open LL measure is 0 but wall LL points exist.")
    print("3. The naive integral proof route fails.  The conjecture")
    print("   integral(g_left+g_right) >= 2/n has explicit violations at n=5 and n=6.")
    print("4. The dynamic signal is exchange, not average size: gap changes are")
    print("   strongly anti-correlated in the tested windows, and LL appears as a")
    print("   bridge from one-sided-long LS into one-sided-long SL.")
    print("5. The inclusion-exclusion/pigeonhole identity is useful bookkeeping but")
    print("   tautological by itself: the overlap excess is exactly the LL measure.")
    print("6. Tournament Analysis on representative speed sets is transitive in this")
    print("   audit, ranking source-corridor abundance rather than exposing a cyclic")
    print("   obstruction.  That pushes the proof burden back to labelled endpoint")
    print("   pressure and compactified walk constraints.")
    print()


if __name__ == "__main__":
    main()
