#!/usr/bin/env python3
"""Probe the LRC <-> A000568 tournament-isoclass analogy.

codex-2026-06-01-S510

The user hypothesis is that LRC is ultimately analogous to a problem on
tournament isomorphism classes and A000568.  This script makes that precise in
three finite ways:

1. Use the odd-cycle Burnside formula to print A000568(n).
2. Turn small runner systems into exact closed walks through tournament
   isomorphism classes, with an observer-marked fiber over each unmarked class.
3. For n=14 and n=18 LRC rows, where exact isomorphism reduction is too large,
   record coarse phase fingerprints.  These are lower-resolution shadows of
   the same A000568 chamber walk.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, permutations
from math import factorial, gcd


ONE = Fraction(1, 1)
HALF = Fraction(1, 2)


@lru_cache(maxsize=None)
def pair_list(n: int) -> tuple[tuple[int, int], ...]:
    return tuple((i, j) for i in range(n) for j in range(i + 1, n))


@lru_cache(maxsize=None)
def pair_index(n: int) -> dict[tuple[int, int], int]:
    return {pair: idx for idx, pair in enumerate(pair_list(n))}


@lru_cache(maxsize=None)
def all_perms(n: int) -> tuple[tuple[int, ...], ...]:
    return tuple(permutations(range(n)))


@lru_cache(maxsize=None)
def observer_perms(n: int) -> tuple[tuple[int, ...], ...]:
    return tuple((0,) + rest for rest in permutations(range(1, n)))


def circle(value: Fraction) -> Fraction:
    return value % ONE


def dist0(value: Fraction) -> Fraction:
    value = circle(value)
    return min(value, ONE - value)


def odd_partitions(n: int) -> list[tuple[tuple[int, int], ...]]:
    """Partitions of n into odd parts, returned as ((part, multiplicity), ...)."""
    out: list[tuple[tuple[int, int], ...]] = []
    odd_parts = list(range(n if n % 2 else n - 1, 0, -2))

    def rec(idx: int, remaining: int, parts: list[tuple[int, int]]) -> None:
        if remaining == 0:
            out.append(tuple(parts))
            return
        if idx >= len(odd_parts):
            return
        part = odd_parts[idx]
        rec(idx + 1, remaining, parts)
        for mult in range(1, remaining // part + 1):
            rec(idx + 1, remaining - mult * part, parts + [(part, mult)])

    rec(0, n, [])
    return out


def davis_exponent(parts: tuple[tuple[int, int], ...]) -> int:
    """Davis exponent for tournament Burnside fixed sets."""
    total = 0
    for i, (ki, mi) in enumerate(parts):
        total += mi * ((ki - 1) // 2)
        total += mi * (mi - 1) // 2 * ki
        for kj, mj in parts[i + 1 :]:
            total += mi * mj * gcd(ki, kj)
    return total


def z_lambda(parts: tuple[tuple[int, int], ...]) -> int:
    total = 1
    for part, mult in parts:
        total *= part**mult * factorial(mult)
    return total


def a000568(n: int) -> int:
    """Number of unlabeled tournaments on n vertices."""
    if n <= 1:
        return 1
    nfact = factorial(n)
    total = 0
    for parts in odd_partitions(n):
        total += (nfact // z_lambda(parts)) * (2 ** davis_exponent(parts))
    return total // nfact


def phase_bits(speeds: tuple[int, ...], t: Fraction) -> int:
    """Half-turn phase tournament; bit 1 for i->j when i<j."""
    bits = 0
    for idx, (i, j) in enumerate(pair_list(len(speeds))):
        f = circle(Fraction(speeds[i] - speeds[j]) * t)
        if f == 0 or f == HALF:
            winner_is_i = True
        else:
            winner_is_i = f < HALF
        if winner_is_i:
            bits |= 1 << idx
    return bits


def edge(bits: int, n: int, i: int, j: int) -> bool:
    if i == j:
        raise ValueError("no loop edge")
    if i < j:
        return bool((bits >> pair_index(n)[(i, j)]) & 1)
    return not bool((bits >> pair_index(n)[(j, i)]) & 1)


@lru_cache(maxsize=None)
def row_masks(bits: int, n: int) -> tuple[int, ...]:
    masks = [0] * n
    for i, j in pair_list(n):
        if edge(bits, n, i, j):
            masks[i] |= 1 << j
        else:
            masks[j] |= 1 << i
    return tuple(masks)


def scores(bits: int, n: int) -> tuple[int, ...]:
    return tuple(mask.bit_count() for mask in row_masks(bits, n))


def score_sequence(bits: int, n: int) -> tuple[int, ...]:
    return tuple(sorted(scores(bits, n)))


def score_hist(bits: int, n: int) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(scores(bits, n)).items()))


def directed_triangles(bits: int, n: int) -> int:
    total = 0
    for a, b, c in combinations(range(n), 3):
        if (edge(bits, n, a, b) and edge(bits, n, b, c) and edge(bits, n, c, a)) or (
            edge(bits, n, a, c) and edge(bits, n, c, b) and edge(bits, n, b, a)
        ):
            total += 1
    return total


def largest_scc(bits: int, n: int) -> int:
    masks = row_masks(bits, n)
    graph = [[j for j in range(n) if (masks[i] >> j) & 1] for i in range(n)]
    reverse = [[] for _ in range(n)]
    for i, row in enumerate(graph):
        for j in row:
            reverse[j].append(i)

    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for nxt in graph[v]:
            if nxt not in seen:
                dfs(nxt)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)

    best = 0
    seen.clear()
    for start in reversed(order):
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        size = 0
        while stack:
            v = stack.pop()
            size += 1
            for nxt in reverse[v]:
                if nxt not in seen:
                    seen.add(nxt)
                    stack.append(nxt)
        best = max(best, size)
    return best


@lru_cache(maxsize=None)
def hamiltonian_paths(bits: int, n: int) -> int:
    masks = row_masks(bits, n)
    full = (1 << n) - 1

    @lru_cache(maxsize=None)
    def dp(mask: int, last: int) -> int:
        if mask == full:
            return 1
        total = 0
        available = full ^ mask
        for nxt in range(n):
            if (available >> nxt) & 1 and (masks[last] >> nxt) & 1:
                total += dp(mask | (1 << nxt), nxt)
        return total

    return sum(dp(1 << start, start) for start in range(n))


def encode_under_perm(bits: int, n: int, perm: tuple[int, ...]) -> int:
    out = 0
    for idx, (i, j) in enumerate(pair_list(n)):
        if edge(bits, n, perm[i], perm[j]):
            out |= 1 << idx
    return out


@lru_cache(maxsize=None)
def canonical(bits: int, n: int) -> int:
    return min(encode_under_perm(bits, n, perm) for perm in all_perms(n))


@lru_cache(maxsize=None)
def marked_canonical(bits: int, n: int) -> int:
    return min(encode_under_perm(bits, n, perm) for perm in observer_perms(n))


def wall_times(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    walls: set[Fraction] = set()
    for i, j in combinations(range(len(speeds)), 2):
        d = abs(speeds[i] - speeds[j])
        if d == 0:
            continue
        for m in range(2 * d):
            walls.add(Fraction(m, 2 * d))
    return tuple(sorted(walls))


def clock_midpoints(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    walls = list(wall_times(speeds))
    if not walls:
        return (Fraction(1, 2),)
    walls.append(ONE)
    return tuple((a + b) / 2 for a, b in zip(walls, walls[1:]))


def observer_margin(speeds: tuple[int, ...], t: Fraction) -> Fraction:
    n = len(speeds)
    threshold = Fraction(1, n)
    return min(dist0(Fraction(speed) * t) for speed in speeds[1:]) / threshold


def safe_for_observer(speeds: tuple[int, ...], t: Fraction) -> bool:
    return observer_margin(speeds, t) >= 1


def split_interval(lo: Fraction, hi: Fraction) -> list[tuple[Fraction, Fraction]]:
    """Split a circular interval into [0,1] pieces."""
    width = hi - lo
    if width >= ONE:
        return [(Fraction(0), ONE)]
    lo = circle(lo)
    hi = lo + width
    if hi <= ONE:
        return [(lo, hi)]
    return [(lo, ONE), (Fraction(0), hi - ONE)]


def forbidden_intervals(speeds: tuple[int, ...]) -> tuple[tuple[Fraction, Fraction], ...]:
    """Open forbidden arcs for the marked observer, represented by closures."""
    n = len(speeds)
    intervals: list[tuple[Fraction, Fraction]] = []
    for speed in speeds[1:]:
        radius = Fraction(1, n * speed)
        for k in range(speed):
            center = Fraction(k, speed)
            intervals.extend(split_interval(center - radius, center + radius))
    return tuple(intervals)


def merge_intervals(
    intervals: tuple[tuple[Fraction, Fraction], ...],
) -> tuple[tuple[Fraction, Fraction], ...]:
    if not intervals:
        return ()
    merged: list[list[Fraction]] = []
    for lo, hi in sorted(intervals):
        if not merged or lo > merged[-1][1]:
            merged.append([lo, hi])
        elif hi > merged[-1][1]:
            merged[-1][1] = hi
    return tuple((lo, hi) for lo, hi in merged)


def safe_gaps(speeds: tuple[int, ...]) -> tuple[tuple[Fraction, Fraction], ...]:
    """Positive-length safe intervals in the complement of forbidden arcs."""
    merged = merge_intervals(forbidden_intervals(speeds))
    if not merged:
        return ((Fraction(0), ONE),)
    if len(merged) == 1 and merged[0] == (Fraction(0), ONE):
        return ()

    gaps: list[tuple[Fraction, Fraction]] = []
    for (_, hi), (next_lo, _) in zip(merged, merged[1:]):
        if hi < next_lo:
            gaps.append((hi, next_lo))
    last_hi = merged[-1][1]
    first_lo = merged[0][0]
    if last_hi < ONE or first_lo > 0:
        gaps.append((last_hi, first_lo + ONE))
    return tuple(gaps)


def midpoint_on_circle(gap: tuple[Fraction, Fraction]) -> Fraction:
    return circle((gap[0] + gap[1]) / 2)


def safe_witnesses(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    """Exact witness samples: safe endpoints plus one midpoint per positive gap."""
    candidates = {Fraction(0)}
    for lo, hi in forbidden_intervals(speeds):
        candidates.add(circle(lo))
        candidates.add(circle(hi))
    for gap in safe_gaps(speeds):
        candidates.add(midpoint_on_circle(gap))
    return tuple(sorted(t for t in candidates if safe_for_observer(speeds, t)))


def sum_abs_diffs(speeds: tuple[int, ...]) -> int:
    return sum(abs(a - b) for a, b in combinations(speeds, 2))


def compressed_walk(values: list[int]) -> list[tuple[int, int]]:
    if not values:
        return []
    out = []
    last = values[0]
    count = 1
    for value in values[1:]:
        if value == last:
            count += 1
        else:
            out.append((last, count))
            last = value
            count = 1
    out.append((last, count))
    return out


@dataclass(frozen=True)
class SmallSummary:
    label: str
    n: int
    cells: int
    a_count: int
    unmarked: int
    marked: int
    safe_samples: int
    safe_unmarked: int
    safe_marked: int
    mixed_safe_unmarked: int
    mixed_observer_score_unmarked: int
    h_values: tuple[int, ...]
    class_runs: tuple[tuple[int, int], ...]


def exact_small_summary(label: str, speeds: tuple[int, ...]) -> SmallSummary:
    n = len(speeds)
    mids = clock_midpoints(speeds)
    witnesses = safe_witnesses(speeds)
    a_count = a000568(n)
    class_ids: dict[int, int] = {}
    class_walk: list[int] = []
    unmarked_safe: dict[int, set[bool]] = defaultdict(set)
    unmarked_observer_scores: dict[int, set[int]] = defaultdict(set)
    safe_unmarked: set[int] = set()
    safe_marked: set[int] = set()
    marked: set[int] = set()
    h_values: set[int] = set()

    for t in mids:
        bits = phase_bits(speeds, t)
        can = canonical(bits, n)
        mcan = marked_canonical(bits, n)
        class_ids.setdefault(can, len(class_ids))
        class_walk.append(class_ids[can])
        marked.add(mcan)
        unmarked_safe[can].add(False)
        unmarked_observer_scores[can].add(scores(bits, n)[0])
        h_values.add(hamiltonian_paths(bits, n))

    for t in witnesses:
        bits = phase_bits(speeds, t)
        can = canonical(bits, n)
        mcan = marked_canonical(bits, n)
        safe_unmarked.add(can)
        safe_marked.add(mcan)
        unmarked_safe[can].add(True)
        unmarked_observer_scores[can].add(scores(bits, n)[0])

    return SmallSummary(
        label=label,
        n=n,
        cells=len(mids),
        a_count=a_count,
        unmarked=len(class_ids),
        marked=len(marked),
        safe_samples=len(witnesses),
        safe_unmarked=len(safe_unmarked),
        safe_marked=len(safe_marked),
        mixed_safe_unmarked=sum(1 for values in unmarked_safe.values() if len(values) > 1),
        mixed_observer_score_unmarked=sum(
            1 for values in unmarked_observer_scores.values() if len(values) > 1
        ),
        h_values=tuple(sorted(h_values)),
        class_runs=tuple(compressed_walk(class_walk)[:18]),
    )


def class_fiber_report(label: str, speeds: tuple[int, ...]) -> None:
    n = len(speeds)
    by_class: dict[int, dict[str, object]] = {}
    local_ids: dict[int, int] = {}
    for t in clock_midpoints(speeds):
        bits = phase_bits(speeds, t)
        can = canonical(bits, n)
        local_ids.setdefault(can, len(local_ids))
        bucket = by_class.setdefault(
            can,
            {
                "cells": 0,
                "safe": 0,
                "scores0": set(),
                "margins": [],
                "H": hamiltonian_paths(bits, n),
                "score": score_sequence(bits, n),
            },
        )
        bucket["cells"] = int(bucket["cells"]) + 1
        if safe_for_observer(speeds, t):
            bucket["safe"] = int(bucket["safe"]) + 1
        bucket["scores0"].add(scores(bits, n)[0])  # type: ignore[union-attr]
        bucket["margins"].append(observer_margin(speeds, t))  # type: ignore[union-attr]

    safe_by_class: Counter[int] = Counter()
    for t in safe_witnesses(speeds):
        safe_by_class[canonical(phase_bits(speeds, t), n)] += 1

    print(f"\nMarked fiber sample: {label}")
    print("  id  H     score-seq             cells safeW observer-score  margin-range")
    for can, bucket in sorted(by_class.items(), key=lambda item: (item[1]["H"], item[0])):
        margins = bucket["margins"]  # type: ignore[assignment]
        margin_range = f"{float(min(margins)):.3f}..{float(max(margins)):.3f}"
        score_values = ",".join(str(v) for v in sorted(bucket["scores0"]))  # type: ignore[arg-type]
        print(
            f"  {local_ids[can]:2d}  {bucket['H']:<5} "
            f"{str(bucket['score']):<21} {bucket['cells']:>5} "
            f"{safe_by_class[can]:>5}  {score_values:<14} {margin_range}"
        )


def ladder(n: int, scale: int, skip: int) -> tuple[int, ...]:
    speeds = tuple(sorted({1} | {scale * q for q in range(1, n) if q != skip}))
    if len(speeds) != n - 1:
        raise ValueError((n, scale, skip, speeds))
    return speeds


def phase_fingerprint(bits: int, n: int) -> tuple[object, ...]:
    out = scores(bits, n)
    return (
        score_hist(bits, n),
        directed_triangles(bits, n),
        out[0],
        largest_scc(bits, n),
        sum(1 for value in out if value == n - 1),
        sum(1 for value in out if value == 0),
    )


@dataclass(frozen=True)
class LargeSummary:
    label: str
    n: int
    cells: int
    a_count: int
    fingerprints: int
    safe_samples: int
    safe_fingerprints: int
    observer_scores: tuple[int, ...]
    c3_range: tuple[int, int]
    largest_scc_range: tuple[int, int]


def coarse_large_summary(label: str, total_n: int, moving_speeds: tuple[int, ...]) -> LargeSummary:
    speeds = (0,) + moving_speeds
    fingerprints: set[tuple[object, ...]] = set()
    safe_fingerprints: set[tuple[object, ...]] = set()
    observer_scores: set[int] = set()
    c3_values: list[int] = []
    scc_values: list[int] = []
    mids = clock_midpoints(speeds)

    for t in mids:
        bits = phase_bits(speeds, t)
        fp = phase_fingerprint(bits, total_n)
        fingerprints.add(fp)
        observer_scores.add(scores(bits, total_n)[0])
        c3_values.append(fp[1])  # type: ignore[arg-type]
        scc_values.append(fp[3])  # type: ignore[arg-type]

    witnesses = safe_witnesses(speeds)
    for t in witnesses:
        bits = phase_bits(speeds, t)
        safe_fingerprints.add(phase_fingerprint(bits, total_n))

    return LargeSummary(
        label=label,
        n=total_n,
        cells=len(mids),
        a_count=a000568(total_n),
        fingerprints=len(fingerprints),
        safe_samples=len(witnesses),
        safe_fingerprints=len(safe_fingerprints),
        observer_scores=tuple(sorted(observer_scores)),
        c3_range=(min(c3_values), max(c3_values)),
        largest_scc_range=(min(scc_values), max(scc_values)),
    )


def print_a000568_table() -> None:
    print("=" * 78)
    print("A000568 as the chamber count for unmarked tournament shapes")
    print("=" * 78)
    print(" n | A000568(n) | odd partitions used")
    print("---+------------+--------------------")
    for n in range(1, 19):
        print(f"{n:2d} | {a000568(n):>10d} | {len(odd_partitions(n)):>18d}")


def print_small_clock_section() -> None:
    print("\n" + "=" * 78)
    print("Exact small runner clocks as walks through A000568 chambers")
    print("=" * 78)
    families = [
        ("N5 initial 0..4", (0, 1, 2, 3, 4)),
        ("N5 prime-like", (0, 2, 3, 5, 7)),
        ("N5 lopsided", (0, 1, 4, 9, 11)),
        ("N6 initial 0..5", (0, 1, 2, 3, 4, 5)),
        ("N6 coprime spread", (0, 1, 5, 7, 11, 13)),
        ("N7 initial 0..6", (0, 1, 2, 3, 4, 5, 6)),
        ("N7 prime-like", (0, 2, 3, 5, 7, 11, 13)),
    ]

    print(
        "case                  N cells diffsum unmarked/A  marked safe-samp  "
        "safe-classes mixed-safe mixed-score H-values"
    )
    print("-" * 118)
    for label, speeds in families:
        s = exact_small_summary(label, speeds)
        print(
            f"{label:<22} {s.n:>1} {s.cells:>5} {sum_abs_diffs(speeds):>7} "
            f"{s.unmarked:>4}/{s.a_count:<4} {s.marked:>6} "
            f"{s.safe_samples:>5} {s.safe_unmarked:>3}/{s.safe_marked:<3} "
            f"{s.mixed_safe_unmarked:>10} {s.mixed_observer_score_unmarked:>11} "
            f"{s.h_values}"
        )
        run_text = " ".join(f"{cid}^{count}" for cid, count in s.class_runs)
        print(f"  local class run prefix: {run_text}")

    class_fiber_report("N5 initial 0..4", (0, 1, 2, 3, 4))
    class_fiber_report("N6 initial 0..5", (0, 1, 2, 3, 4, 5))


def print_large_shadow_section() -> None:
    print("\n" + "=" * 78)
    print("n=14 and n=18 LRC rows: coarse shadows of huge A000568 walks")
    print("=" * 78)
    rows = [
        ("n14 initial", 14, tuple(range(1, 14))),
        ("n14 row-parent", 14, ladder(14, 7, 6)),
        ("n14 gate", 14, ladder(14, 14, 6)),
        ("n18 initial", 18, tuple(range(1, 18))),
        ("n18 row-parent", 18, ladder(18, 9, 8)),
        ("n18 gate", 18, ladder(18, 18, 8)),
    ]

    print(
        "case             N cells     A000568(N)         fp safe-samp  "
        "safe-fp obs-score-range c3-range         SCC-range"
    )
    print("-" * 122)
    for label, n, speeds in rows:
        s = coarse_large_summary(label, n, speeds)
        obs_range = f"{min(s.observer_scores)}..{max(s.observer_scores)}"
        c3_range = f"{s.c3_range[0]}..{s.c3_range[1]}"
        scc_range = f"{s.largest_scc_range[0]}..{s.largest_scc_range[1]}"
        print(
            f"{label:<16} {s.n:>2} {s.cells:>5} {s.a_count:>22d} "
            f"{s.fingerprints:>8} {s.safe_samples:>10} {s.safe_fingerprints:>7} "
            f"{obs_range:<15} {c3_range:<16} {scc_range}"
        )

    print(
        "\nInterpretation: fingerprints are coarse, so they are a lower-resolution "
        "shadow of visited isomorphism classes.  But cells are an upper bound on "
        "visited classes, and cells << A000568(N) already says LRC rows trace a "
        "microscopic walk in the full chamber space."
    )


def print_convoluted_paths() -> None:
    print("\n" + "=" * 78)
    print("Creative paths suggested by the data")
    print("=" * 78)
    paths = [
        (
            "unmarked chamber problem",
            "A000568(N) is the number of unlabelled tournament chambers; an LRC "
            "speed row is a closed walk through a tiny subset of these chambers.",
        ),
        (
            "marked fiber problem",
            "The observer is not quotientable.  The LRC object is a marked fiber "
            "over A000568: one unmarked chamber can carry safe and unsafe "
            "observer placements.",
        ),
        (
            "target-avoidance problem",
            "A counterexample becomes a closed arithmetic walk avoiding all "
            "observer-safe marked chambers while respecting endpoint-wall "
            "protection.",
        ),
        (
            "Burnside resonance problem",
            "A000568 skips even permutation cycles; runner clocks collapse wall "
            "events when speed differences share denominators.  Both are "
            "quotient counts controlled by cycle parity and gcd data.",
        ),
        (
            "q-deformed pressure problem",
            "The repo's A(n,q) viewpoint suggests weighting chambers by endpoint "
            "pressure, with q=2 the unmarked tournament count and higher q "
            "remembering labelled protection debt.",
        ),
        (
            "Royle-even shadow",
            "The A000568/even-graph equinumerosity may make endpoint-safe masks "
            "an even-graph shadow of the marked tournament clock.",
        ),
        (
            "H as a chamber coordinate",
            "H is not the whole metric, but within the circular clock image it is "
            "a useful coordinate on chambers; endpoint-aware gauges add the "
            "missing marked coordinates.",
        ),
        (
            "proof by sparse image",
            "For large N, LRC rows occupy cells many orders below A000568(N).  A "
            "proof may only need the circular arithmetic image, not all "
            "tournament classes.",
        ),
    ]
    for idx, (name, text) in enumerate(paths, 1):
        print(f"{idx}. {name}: {text}")


def main() -> None:
    print("LRC / A000568 isomorphism-class analogy (codex-2026-06-01-S510)")
    print_a000568_table()
    print_small_clock_section()
    print_large_shadow_section()
    print_convoluted_paths()


if __name__ == "__main__":
    main()
