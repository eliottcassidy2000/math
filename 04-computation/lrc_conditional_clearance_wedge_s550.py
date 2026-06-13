#!/usr/bin/env python3
"""Exact checks for THM-395 and the LRC conditional-clearance cascade.

The tournament half verifies that backward-wedge mass is exactly three times
the directed-triangle count, and that zero wedge debt is equivalent to
transitivity.  The LRC half prints exact safe-measure factors for sample
runner rows, interpreting a cascade as a product of conditional clearances.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations


def fmt_fraction(value: Fraction | None) -> str:
    if value is None:
        return "undefined"
    if value.denominator == 1:
        return str(value.numerator)
    return f"{value.numerator}/{value.denominator}"


def tournament_edge(n: int, mask: int):
    """Return an edge predicate for the tournament encoded by mask."""
    pair_index: dict[tuple[int, int], int] = {}
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            pair_index[(i, j)] = bit
            bit += 1

    def edge(i: int, j: int) -> bool:
        if i == j:
            raise ValueError("loop requested")
        if i < j:
            return bool((mask >> pair_index[(i, j)]) & 1)
        return not bool((mask >> pair_index[(j, i)]) & 1)

    return edge


def directed_triangle_count(n: int, edge) -> int:
    count = 0
    for a, b, c in combinations(range(n), 3):
        if edge(a, b) and edge(b, c) and edge(c, a):
            count += 1
        elif edge(a, c) and edge(c, b) and edge(b, a):
            count += 1
    return count


def backward_wedge_sum(n: int, edge) -> int:
    total = 0
    for x in range(n):
        for y in range(n):
            if x == y or not edge(x, y):
                continue
            for z in range(n):
                if z == x or z == y:
                    continue
                if edge(z, x) and edge(y, z):
                    total += 1
    return total


def is_transitive(n: int, edge) -> bool:
    for x in range(n):
        for y in range(n):
            if x == y or not edge(x, y):
                continue
            for z in range(n):
                if z == x or z == y:
                    continue
                if edge(y, z) and not edge(x, z):
                    return False
    return True


def score_histogram(n: int, edge) -> tuple[int, ...]:
    scores = []
    for i in range(n):
        scores.append(sum(1 for j in range(n) if i != j and edge(i, j)))
    return tuple(sorted(scores))


def scc_sizes(n: int, edge) -> tuple[int, ...]:
    graph = [[j for j in range(n) if i != j and edge(i, j)] for i in range(n)]
    reverse = [[j for j in range(n) if i != j and edge(j, i)] for i in range(n)]
    seen: set[int] = set()
    order: list[int] = []

    def dfs(v: int) -> None:
        seen.add(v)
        for w in graph[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in range(n):
        if v not in seen:
            dfs(v)

    seen.clear()
    sizes = []

    def rdfs(v: int) -> int:
        seen.add(v)
        size = 1
        for w in reverse[v]:
            if w not in seen:
                size += rdfs(w)
        return size

    for v in reversed(order):
        if v not in seen:
            sizes.append(rdfs(v))
    return tuple(sorted(sizes, reverse=True))


def hamiltonian_path_count(n: int, edge) -> int:
    total = 0
    for path in permutations(range(n)):
        if all(edge(path[i], path[i + 1]) for i in range(n - 1)):
            total += 1
    return total


def exhaustive_tournament_checks(max_n: int = 6) -> list[str]:
    lines = []
    for n in range(3, max_n + 1):
        total = 1 << (n * (n - 1) // 2)
        c3_hist: Counter[int] = Counter()
        transitive_count = 0
        zero_wedge_count = 0
        max_wedge = 0
        score_samples: Counter[tuple[int, ...]] = Counter()
        for mask in range(total):
            edge = tournament_edge(n, mask)
            c3 = directed_triangle_count(n, edge)
            wedge = backward_wedge_sum(n, edge)
            transitive = is_transitive(n, edge)
            zero_wedge = wedge == 0
            if wedge != 3 * c3:
                raise AssertionError((n, mask, wedge, c3))
            if transitive != zero_wedge:
                raise AssertionError((n, mask, transitive, zero_wedge))
            c3_hist[c3] += 1
            score_samples[score_histogram(n, edge)] += 1
            transitive_count += int(transitive)
            zero_wedge_count += int(zero_wedge)
            max_wedge = max(max_wedge, wedge)
        common_scores = ", ".join(
            f"{hist}:{count}" for hist, count in score_samples.most_common(4)
        )
        lines.append(
            "n={n}: tournaments={total}, transitive={transitive_count}, "
            "zero_wedge={zero_wedge_count}, max_wedge={max_wedge}, "
            "c3_hist={c3_hist}, common_score_hists={common_scores}".format(
                n=n,
                total=total,
                transitive_count=transitive_count,
                zero_wedge_count=zero_wedge_count,
                max_wedge=max_wedge,
                c3_hist=dict(sorted(c3_hist.items())),
                common_scores=common_scores,
            )
        )
    return lines


def cyclic_regular_edge(n: int):
    winners = set(range(1, (n + 1) // 2))

    def edge(i: int, j: int) -> bool:
        if i == j:
            raise ValueError("loop requested")
        return ((j - i) % n) in winners

    return edge


def transitive_edge(n: int):
    def edge(i: int, j: int) -> bool:
        if i == j:
            raise ValueError("loop requested")
        return i < j

    return edge


def sample_fingerprints() -> list[str]:
    samples = [
        ("transitive_5", 5, transitive_edge(5)),
        ("directed_cycle_3", 3, cyclic_regular_edge(3)),
        ("regular_cycle_5", 5, cyclic_regular_edge(5)),
    ]
    lines = []
    for name, n, edge in samples:
        c3 = directed_triangle_count(n, edge)
        wedge = backward_wedge_sum(n, edge)
        lines.append(
            "{name}: scores={scores}, c3={c3}, wedge={wedge}, "
            "scc_sizes={scc}, hamiltonian_paths={hp}".format(
                name=name,
                scores=score_histogram(n, edge),
                c3=c3,
                wedge=wedge,
                scc=scc_sizes(n, edge),
                hp=hamiltonian_path_count(n, edge),
            )
        )
    return lines


def split_circle_interval(center: Fraction, radius: Fraction) -> list[tuple[Fraction, Fraction]]:
    one = Fraction(1, 1)
    start = center - radius
    end = center + radius
    if start < 0:
        return [(start + one, one), (Fraction(0), end)]
    if end > one:
        return [(start, one), (Fraction(0), end - one)]
    return [(start, end)]


def merge_intervals(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    intervals = sorted((a, b) for a, b in intervals if b > a)
    merged: list[list[Fraction]] = []
    for start, end in intervals:
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        elif end > merged[-1][1]:
            merged[-1][1] = end
    return [(a, b) for a, b in merged]


def safe_measure(n: int, speeds: list[int]) -> Fraction:
    intervals: list[tuple[Fraction, Fraction]] = []
    for speed in speeds:
        radius = Fraction(1, n * speed)
        for residue in range(speed):
            center = Fraction(residue, speed)
            intervals.extend(split_circle_interval(center, radius))
    covered = sum((b - a for a, b in merge_intervals(intervals)), Fraction(0))
    return Fraction(1) - covered


def clearance_rows(n: int, speeds: list[int]) -> list[tuple[int, int, Fraction, Fraction | None]]:
    rows = []
    prefix: list[int] = []
    previous = Fraction(1)
    for k, speed in enumerate(speeds, start=1):
        prefix.append(speed)
        measure = safe_measure(n, prefix)
        factor = None if previous == 0 else measure / previous
        rows.append((k, speed, measure, factor))
        previous = measure
    return rows


def clearance_report() -> list[str]:
    cases = [
        ("AP_n5", 5, [1, 2, 3, 4]),
        ("generic_n5", 5, [1, 2, 4, 7]),
        ("AP_n7", 7, [1, 2, 3, 4, 5, 6]),
        ("lacunary_n7", 7, [1, 2, 4, 8, 16, 32]),
    ]
    lines = []
    for name, n, speeds in cases:
        lines.append(f"{name}: n={n}, speeds={speeds}")
        for k, speed, measure, factor in clearance_rows(n, speeds):
            lines.append(
                "  k={k}, speed={speed}, mu(F_k)={measure}, P_k={factor}".format(
                    k=k,
                    speed=speed,
                    measure=fmt_fraction(measure),
                    factor=fmt_fraction(factor),
                )
            )
    return lines


def main() -> None:
    print("THM-395 backward-wedge transitivity verifier")
    print("pairwise observable: oriented tournament edge X->Y")
    print("switch/gauge: the tournament orientation; no ties in these exact checks")
    print("tie Hamiltonian path: not used here because all pair data are strict")
    print()
    print("exhaustive tournament checks")
    for line in exhaustive_tournament_checks():
        print(line)
    print()
    print("sample tournament fingerprints")
    for line in sample_fingerprints():
        print(line)
    print()
    print("LRC conditional-clearance products")
    for line in clearance_report():
        print(line)


if __name__ == "__main__":
    main()
