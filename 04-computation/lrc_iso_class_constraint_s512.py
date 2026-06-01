#!/usr/bin/env python3
"""
lrc_iso_class_constraint_s512.py

codex-2026-06-01 S512

Explore the user's hypothesis that the Lonely Runner Conjecture is ultimately a
restriction on tournament isomorphism classes, in the A000568 sense.

The point of this script is not to pretend that raw unmarked tournament classes
already prove LRC.  It tests the exact opposite question carefully:

    Which lift of the runner problem to a tournament-class object is strong
    enough that "visited class lies in a forced subset" becomes a witness?

We compare several fibers over ordinary A000568 classes.

* phase_half:
    vertices are runners; arcs use the half-turn circular comparator.

* phase_marked_observer:
    same tournament, but the stationary observer is a colored vertex.

* phase_safe_colored:
    same tournament, with colors recording observer / safe moving / unsafe
    moving.  This is tautological, but it gives an upper bound on what a
    colored class proof can do.

* gap_rank_marked_adjacent:
    vertices are the circular gaps between consecutive runners; arcs rank gaps
    by length; the two gaps adjacent to the observer are colored.  This is the
    most geometric non-tautological lift in this audit.

* gap_threshold_fiber:
    the same gap-rank tournament, but gap colors also remember whether a gap is
    at least 1/n.  This is the first honest threshold-decorated fiber.

* pair_deficit_origin_marked:
    vertices are unordered runner-pairs; arcs compare the LRC danger deficit
    max(0, 1/n - dist(i,j)); pair-cells incident to the observer are colored.

* pair_deficit_threshold_fiber:
    the same pair-cell tournament, with colors recording whether each pair-cell
    has positive danger deficit.

The sampled clock states include both open cells and exact wall times for the
half-turn clock and the LRC endpoint clock.  Thus equality witnesses such as
the initial-segment time t=1/n are not missed.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
from math import comb, factorial, gcd, prod
from typing import Callable


Adj = tuple[tuple[int, ...], ...]
ClassKey = tuple[tuple[int, ...], tuple[int, ...]]


@dataclass(frozen=True)
class State:
    speeds: tuple[int, ...]
    time: Fraction
    good: bool
    classes: dict[str, ClassKey]


def frac(value: Fraction) -> Fraction:
    return value - (value.numerator // value.denominator)


def circular_distance(a: Fraction, b: Fraction) -> Fraction:
    diff = frac(a - b)
    return min(diff, 1 - diff)


def positions(speeds: tuple[int, ...], t: Fraction) -> tuple[Fraction, ...]:
    return tuple(frac(Fraction(speed) * t) for speed in speeds)


def primitive_speed_sets(total_n: int, max_speed: int) -> list[tuple[int, ...]]:
    out: list[tuple[int, ...]] = []
    for moving in combinations(range(1, max_speed + 1), total_n - 1):
        g = 0
        for speed in moving:
            g = gcd(g, speed)
        if g == 1:
            out.append((0,) + moving)
    return out


def lrc_good(speeds: tuple[int, ...], t: Fraction) -> bool:
    threshold = Fraction(1, len(speeds))
    pos = positions(speeds, t)
    return all(circular_distance(pos[0], pos[i]) >= threshold for i in range(1, len(speeds)))


def event_times(speeds: tuple[int, ...]) -> list[Fraction]:
    total_n = len(speeds)
    out: set[Fraction] = {Fraction(0), Fraction(1)}

    # Half-turn phase tournament walls.
    for i, j in combinations(range(total_n), 2):
        d = abs(speeds[i] - speeds[j])
        if d == 0:
            continue
        for m in range(2 * d + 1):
            out.add(Fraction(m, 2 * d))

    # LRC endpoint walls: ||v t|| = 1/n.
    for speed in speeds[1:]:
        for m in range(speed):
            out.add(Fraction(m * total_n + 1, total_n * speed))
            out.add(Fraction(m * total_n + total_n - 1, total_n * speed))

    return sorted(time for time in out if 0 <= time <= 1)


def sampled_times(speeds: tuple[int, ...]) -> list[Fraction]:
    walls = event_times(speeds)
    out: set[Fraction] = set(walls)
    for left, right in zip(walls, walls[1:]):
        if left < right:
            out.add((left + right) / 2)
    return sorted(out)


def adjacency_from_winners(n: int, winner: Callable[[int, int], int]) -> Adj:
    adj = [[0] * n for _ in range(n)]
    for i, j in combinations(range(n), 2):
        w = winner(i, j)
        loser = j if w == i else i
        adj[w][loser] = 1
    return tuple(tuple(row) for row in adj)


@lru_cache(maxsize=None)
def canonical_colored(adj: Adj, colors: tuple[int, ...] | None = None) -> ClassKey:
    n = len(adj)
    if colors is None:
        colors = (0,) * n
    best: ClassKey | None = None
    for perm in permutations(range(n)):
        perm_colors = tuple(colors[perm[i]] for i in range(n))
        bits = tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(i + 1, n))
        key = (perm_colors, bits)
        if best is None or key < best:
            best = key
    assert best is not None
    return best


def hamiltonian_path_count(adj: Adj) -> int:
    n = len(adj)
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(1 << n)]
    for start in range(n):
        dp[1 << start][start] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp[mask][last]
            if not val:
                continue
            for nxt in range(n):
                if not ((mask >> nxt) & 1) and adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += val
    return sum(dp[full])


def score_sequence(adj: Adj) -> tuple[int, ...]:
    return tuple(sorted(sum(row) for row in adj))


def triangle_count(adj: Adj) -> int:
    total = 0
    for a, b, c in combinations(range(len(adj)), 3):
        if (adj[a][b] and adj[b][c] and adj[c][a]) or (
            adj[a][c] and adj[c][b] and adj[b][a]
        ):
            total += 1
    return total


def phase_half_adj(speeds: tuple[int, ...], t: Fraction) -> Adj:
    def winner(i: int, j: int) -> int:
        gap = frac(Fraction(speeds[i] - speeds[j]) * t)
        if gap == 0 or gap == Fraction(1, 2):
            return i
        if 0 < gap < Fraction(1, 2):
            return i
        return j

    return adjacency_from_winners(len(speeds), winner)


def safe_colors(speeds: tuple[int, ...], t: Fraction) -> tuple[int, ...]:
    threshold = Fraction(1, len(speeds))
    pos = positions(speeds, t)
    colors = [2]
    for idx in range(1, len(speeds)):
        colors.append(1 if circular_distance(pos[0], pos[idx]) >= threshold else 0)
    return tuple(colors)


def gap_rank_adj_and_colors(
    speeds: tuple[int, ...], t: Fraction
) -> tuple[Adj, tuple[int, ...], tuple[Fraction, ...]]:
    pos = positions(speeds, t)
    order = sorted(range(len(speeds)), key=lambda idx: (pos[idx], idx))
    gaps: list[Fraction] = []
    colors: list[int] = []
    n = len(order)
    for k, idx in enumerate(order):
        nxt = order[(k + 1) % n]
        length = frac(pos[nxt] - pos[idx])
        gaps.append(length)
        colors.append(1 if idx == 0 or nxt == 0 else 0)

    def winner(i: int, j: int) -> int:
        if gaps[i] > gaps[j]:
            return i
        if gaps[j] > gaps[i]:
            return j
        return i

    return adjacency_from_winners(n, winner), tuple(colors), tuple(gaps)


def pair_deficit_adj_and_colors(
    speeds: tuple[int, ...], t: Fraction
) -> tuple[Adj, tuple[int, ...], tuple[Fraction, ...]]:
    labels = tuple(combinations(range(len(speeds)), 2))
    threshold = Fraction(1, len(speeds))
    pos = positions(speeds, t)
    deficits: list[Fraction] = []
    colors: list[int] = []
    for i, j in labels:
        deficits.append(max(Fraction(0), threshold - circular_distance(pos[i], pos[j])))
        colors.append(1 if i == 0 else 0)

    def winner(a: int, b: int) -> int:
        if deficits[a] > deficits[b]:
            return a
        if deficits[b] > deficits[a]:
            return b
        return a

    return adjacency_from_winners(len(labels), winner), tuple(colors), tuple(deficits)


def state_classes(speeds: tuple[int, ...], t: Fraction) -> dict[str, ClassKey]:
    phase = phase_half_adj(speeds, t)
    gap_adj, gap_colors, gaps = gap_rank_adj_and_colors(speeds, t)
    pair_adj, pair_colors, deficits = pair_deficit_adj_and_colors(speeds, t)
    threshold = Fraction(1, len(speeds))
    gap_threshold_colors = tuple(
        (2 if gap_colors[idx] else 0) + (1 if gaps[idx] >= threshold else 0)
        for idx in range(len(gaps))
    )
    pair_threshold_colors = tuple(
        (2 if pair_colors[idx] else 0) + (1 if deficits[idx] == 0 else 0)
        for idx in range(len(deficits))
    )
    return {
        "phase_half": canonical_colored(phase),
        "phase_marked_observer": canonical_colored(phase, (1,) + (0,) * (len(speeds) - 1)),
        "phase_safe_colored": canonical_colored(phase, safe_colors(speeds, t)),
        "gap_rank_marked_adjacent": canonical_colored(gap_adj, gap_colors),
        "gap_threshold_fiber": canonical_colored(gap_adj, gap_threshold_colors),
        "pair_deficit_origin_marked": canonical_colored(pair_adj, pair_colors),
        "pair_deficit_threshold_fiber": canonical_colored(pair_adj, pair_threshold_colors),
    }


def odd_partitions(n: int) -> list[tuple[int, ...]]:
    out: list[tuple[int, ...]] = []

    def rec(remaining: int, max_part: int, parts: list[int]) -> None:
        if remaining == 0:
            out.append(tuple(parts))
            return
        start = min(max_part, remaining)
        if start % 2 == 0:
            start -= 1
        for part in range(start, 0, -2):
            rec(remaining - part, part, parts + [part])

    rec(n, n, [])
    return out


def a000568_term(parts: tuple[int, ...]) -> Fraction:
    counts = Counter(parts)
    exponent = (
        sum(counts[r] * counts[s] * gcd(r, s) for r in counts for s in counts)
        - sum(counts.values())
    ) // 2
    denominator = prod((part ** mult) * factorial(mult) for part, mult in counts.items())
    return Fraction(2**exponent, denominator)


def a000568_burnside(n: int) -> int:
    return int(sum(a000568_term(parts) for parts in odd_partitions(n)))


def all_tournament_class_count(n: int) -> int:
    classes: set[ClassKey] = set()
    pairs = tuple(combinations(range(n), 2))
    for mask in range(1 << len(pairs)):
        adj = [[0] * n for _ in range(n)]
        for bit, (i, j) in enumerate(pairs):
            if (mask >> bit) & 1:
                adj[i][j] = 1
            else:
                adj[j][i] = 1
        classes.add(canonical_colored(tuple(tuple(row) for row in adj)))
    return len(classes)


def analyze_total_n(total_n: int, max_speed: int) -> None:
    speed_sets = primitive_speed_sets(total_n, max_speed)
    lift_stats: dict[str, dict[ClassKey, Counter[str]]] = defaultdict(lambda: defaultdict(Counter))
    speed_visits: dict[str, dict[tuple[int, ...], set[ClassKey]]] = defaultdict(lambda: defaultdict(set))
    witness_times: dict[tuple[int, ...], list[Fraction]] = {}
    state_count = 0

    for speeds in speed_sets:
        good_times: list[Fraction] = []
        for t in sampled_times(speeds):
            good = lrc_good(speeds, t)
            if good:
                good_times.append(t)
            classes = state_classes(speeds, t)
            state_count += 1
            bucket = "good" if good else "bad"
            for lift, key in classes.items():
                lift_stats[lift][key][bucket] += 1
                speed_visits[lift][speeds].add(key)
        if good_times:
            witness_times[speeds] = good_times

    print(f"TOTAL n={total_n} (observer + {total_n - 1} moving), max_speed={max_speed}")
    print("-" * 100)
    print(f"primitive speed sets: {len(speed_sets)}")
    print(f"sampled exact states: {state_count}")
    print(f"speed sets with sampled LRC witness: {len(witness_times)}/{len(speed_sets)}")
    print()
    print("Lift separation by class/fiber:")
    print("  lift                         classes  good-only  bad-only  mixed  certifies speeds")
    for lift in sorted(lift_stats):
        stats = lift_stats[lift]
        good_only = {key for key, c in stats.items() if c["good"] and not c["bad"]}
        bad_only = {key for key, c in stats.items() if c["bad"] and not c["good"]}
        mixed = {key for key, c in stats.items() if c["good"] and c["bad"]}
        certified = sum(
            1 for speeds in speed_sets if speed_visits[lift][speeds] & good_only
        )
        print(
            f"  {lift:<28} {len(stats):>7} {len(good_only):>10} "
            f"{len(bad_only):>9} {len(mixed):>6} {certified:>8}/{len(speed_sets)}"
        )
    print()

    print("Representative tight/equality witnesses:")
    for speeds in list(speed_sets)[:4]:
        goods = witness_times.get(speeds, [])
        text = ", ".join(str(t) for t in goods[:5]) if goods else "-"
        print(f"  speeds={speeds}: {text}")
    initial = tuple(range(total_n))
    if initial in witness_times:
        print(f"  initial segment {initial}: first witnesses {witness_times[initial][:6]}")
    print()


def circular_menu(total_n: int, max_speed: int) -> None:
    speed_sets = primitive_speed_sets(total_n, max_speed)
    open_classes: dict[ClassKey, tuple[int, tuple[int, ...], int, tuple[int, ...]]] = {}
    compact_classes: dict[ClassKey, tuple[int, tuple[int, ...], int, tuple[int, ...]]] = {}
    for speeds in speed_sets:
        walls = event_times(speeds)
        open_times = [(left + right) / 2 for left, right in zip(walls, walls[1:]) if left < right]
        for t in open_times:
            adj = phase_half_adj(speeds, t)
            key = canonical_colored(adj)
            if key not in open_classes:
                open_classes[key] = (
                    hamiltonian_path_count(adj),
                    score_sequence(adj),
                    triangle_count(adj),
                    speeds,
                )
        for t in set(open_times) | set(walls):
            adj = phase_half_adj(speeds, t)
            key = canonical_colored(adj)
            if key not in compact_classes:
                compact_classes[key] = (
                    hamiltonian_path_count(adj),
                    score_sequence(adj),
                    triangle_count(adj),
                    speeds,
                )

    print(f"Circular half-turn menu probe for total n={total_n}")
    print("-" * 100)
    print(
        f"open-cell phase classes: {len(open_classes)}; "
        f"wall-compactified classes: {len(compact_classes)}; "
        f"A000568({total_n})={a000568_burnside(total_n)}"
    )
    print("  open-cell class details:")
    for h, score, tri, speeds in sorted(open_classes.values()):
        print(f"  H={h:<5} score={score} c3={tri:<2} first_speed={speeds}")
    if len(compact_classes) != len(open_classes):
        extra = sorted(set(compact_classes) - set(open_classes))
        print(f"  tie-wall adds {len(extra)} extra phase classes under the fixed tie path")
    print()


def print_a000568() -> None:
    print("A000568 base layer")
    print("=" * 100)
    print("n  Burnside  enumerated<=5  odd partitions")
    for n in range(1, 8):
        enum = all_tournament_class_count(n) if n <= 5 else None
        parts = ", ".join("+".join(map(str, p)) for p in odd_partitions(n)[:5])
        print(f"{n:>2} {a000568_burnside(n):>9} {str(enum or '-'):>13}  {parts}")
    print()


def main() -> None:
    print("LRC AS A CONSTRAINED TOURNAMENT-CLASS PROBLEM (S512)")
    print("=" * 100)
    print(
        "Raw A000568 classes are the base.  LRC needs a runner-clock walk in a\n"
        "decorated fiber over that base: observer mark, safe mask, gap marks, or\n"
        "pair-cell deficit marks.  The question is which decoration is minimal."
    )
    print()

    print_a000568()
    for total_n, max_speed in ((3, 16), (4, 10), (5, 8)):
        circular_menu(total_n, max_speed)
    for total_n, max_speed in ((3, 16), (4, 10)):
        analyze_total_n(total_n, max_speed)

    print("SYNTHESIS")
    print("=" * 100)
    print("1. The open half-turn runner clock lands in a circular-menu subset of")
    print("   A000568 classes.  Equality walls form a tie-path compactification")
    print("   that can add many boundary classes, so tight LRC states need the")
    print("   boundary fiber, not only the open circular menu.")
    print("2. For LRC itself, unmarked phase class is too coarse: good and bad")
    print("   states live in the same A000568 class.  This is the first warning.")
    print("3. Pure gap-rank and pair-deficit rank are still too coarse.  The")
    print("   first separating objects are threshold-decorated fibers over the")
    print("   A000568 base: gap-rank plus 1/n gap colors, or pair-deficit plus")
    print("   zero/nonzero danger colors.")
    print("4. A future n=14 proof by this route would likely show that the hard")
    print("   quotient-ladder walk is trapped in the circular/pair-cell fiber and")
    print("   must hit a good-only decorated class before closing.")


if __name__ == "__main__":
    main()
