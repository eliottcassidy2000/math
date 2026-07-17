#!/usr/bin/env python3
"""Exact ratio-layer bound for the continuous LRC(14) pair covariance.

For coprime positive magnitudes ``a,b`` define

    P(a,b) = (B2((a-b)/14) - B2((a+b)/14)) / (a*b).

The script bounds the total negative weight carried by any graph on thirteen
positive magnitudes.  The layer cake telescopes through seven clique/Turan
caps and two sharp top-layer structural caps; the much larger degree ledger is
not needed.  All arithmetic affecting the result uses ``Fraction``.

Tournament-analysis note: the pair observable is symmetric.  Thresholding it
produces an undirected multiplicative compatibility graph on oriented primitive
ratios; orienting by numerical order would be a gauge with no proof content.
The useful fingerprints are degree, clique number, edge cap, and threshold
flips, not a tie Hamiltonian path.  The challenged vertex assumption is that
vertices must be runners: quotienting by an anchor makes primitive ratios the
vertices.  This preserves all threshold-clique constraints but destroys the
absolute scale, which the covariance itself does not use.

This is an exact computational certificate, not yet a Lean theorem.  The
remaining formal bridges are a replayable certificate for the seven finite
clique upper bounds and the continuum-to-clean-grid discrepancy estimate.
"""

from fractions import Fraction as F
from math import gcd, isqrt


def bernoulli(residue: int) -> F:
    point = F(residue % 14, 14)
    return point * point - point + F(1, 6)


def pair_covariance(first: int, second: int) -> F:
    return (bernoulli(first - second) - bernoulli(first + second)) / (
        first * second
    )


def turan_edges(vertex_count: int, clique_cap: int) -> int:
    quotient, remainder = divmod(vertex_count, clique_cap)
    part_sizes = [quotient + 1] * remainder + [quotient] * (
        clique_cap - remainder
    )
    return (vertex_count * vertex_count - sum(size * size for size in part_sizes)) // 2


TIERS = {
    2: F(6, 637),
    3: F(2, 441),
    4: F(5, 2646),
    5: F(1, 1764),
    6: F(1, 3136),
    7: F(5, 37632),
    8: F(1, 9996),
    9: F(1, 14112),
}
CUT = TIERS[9]


def primitive_negative_ratios(threshold: F) -> dict[tuple[int, int], F]:
    product_limit = int(F(12, 49) / threshold) + 1
    ratios: dict[tuple[int, int], F] = {}
    for first in range(1, isqrt(product_limit) + 1):
        for second in range(first + 1, product_limit // first + 1):
            if gcd(first, second) != 1:
                continue
            covariance = pair_covariance(first, second)
            assert abs(covariance) <= F(12, 49 * first * second)
            if covariance < 0 and -covariance >= threshold:
                ratios[(first, second)] = -covariance
                ratios[(second, first)] = -covariance
    return ratios


def ratio_graph(threshold: F) -> tuple[list[tuple[int, int]], list[int]]:
    allowed = primitive_negative_ratios(threshold)
    vertices = sorted(allowed)
    adjacency = [0] * len(vertices)
    for first_index, (a, b) in enumerate(vertices):
        for second_index in range(first_index):
            c, d = vertices[second_index]
            numerator = a * d
            denominator = b * c
            common = gcd(numerator, denominator)
            if (numerator // common, denominator // common) in allowed:
                adjacency[first_index] |= 1 << second_index
                adjacency[second_index] |= 1 << first_index
    return vertices, adjacency


def maximum_anchored_clique(
    vertices: list[tuple[int, int]], adjacency: list[int]
) -> tuple[int, list[tuple[int, int]]]:
    """Exact bitset branch-and-bound; the omitted anchor contributes one."""

    best: list[int] = []

    def expand(current: list[int], candidates: int) -> None:
        nonlocal best
        if len(current) + candidates.bit_count() <= len(best):
            return
        ordered: list[int] = []
        color_bounds: list[int] = []
        uncolored = candidates
        color = 0
        while uncolored:
            color += 1
            available = uncolored
            while available:
                bit = available & -available
                vertex = bit.bit_length() - 1
                ordered.append(vertex)
                color_bounds.append(color)
                uncolored &= ~bit
                available &= ~bit
                available &= ~adjacency[vertex]
        for position in range(len(ordered) - 1, -1, -1):
            if len(current) + color_bounds[position] <= len(best):
                return
            vertex = ordered[position]
            bit = 1 << vertex
            if candidates & bit:
                expand(current + [vertex], candidates & adjacency[vertex])
                candidates &= ~bit
        if len(current) > len(best):
            best = current[:]

    expand([], (1 << len(vertices)) - 1)
    witness = [vertices[index] for index in best]
    for position, first in enumerate(best):
        for second in best[:position]:
            assert adjacency[first] & (1 << second)
    return 1 + len(best), witness


def main() -> None:
    oriented = primitive_negative_ratios(CUT)
    weights = sorted(value for (a, b), value in oriented.items() if a < b)
    levels = sorted(set(weights))
    assert len(weights) == 1489
    assert len(levels) == 482
    assert levels[-1] == TIERS[2]

    print("METHOD exact Fraction ratio layers; deterministic bitset clique audit")
    print("VERTICES oriented primitive ratios; symmetric threshold graph, no tournament gauge")
    print(f"RATIOS {len(weights)} LEVELS {len(levels)} CUT {CUT}")
    print("CLIQUE TIERS")
    for clique_size in range(3, 10):
        threshold = TIERS[clique_size]
        next_level = next(level for level in levels if level > threshold)
        tier_vertices, tier_adjacency = ratio_graph(threshold)
        tier_clique, tier_witness = maximum_anchored_clique(
            tier_vertices, tier_adjacency
        )
        next_vertices, next_adjacency = ratio_graph(next_level)
        next_clique, _ = maximum_anchored_clique(next_vertices, next_adjacency)
        assert tier_clique == clique_size
        assert next_clique == clique_size - 1
        edge_count = sum(bits.bit_count() for bits in tier_adjacency) // 2
        print(
            f"  size={clique_size} threshold={threshold} next={next_level} "
            f"oriented={len(tier_vertices)} edges={edge_count} "
            f"next_cap={next_clique} witness={tier_witness}"
        )

    total = 78 * TIERS[9]
    total += sum(
        turan_edges(13, clique_cap) * (TIERS[clique_cap] - TIERS[clique_cap + 1])
        for clique_cap in range(3, 9)
    )
    total += 42 * (F(4, 539) - TIERS[3])
    total += 22 * (F(5, 588) - F(4, 539))
    total += 12 * (TIERS[2] - F(5, 588))

    target = F(13, 30)
    expected = F(175847381, 411675264)
    assert total == expected
    assert total < target
    print("LAYER_CAPS tail=78; Turan k=8..3; top=42,22,12")
    print(f"TOTAL {total} = {float(total):.15f}")
    print(f"TARGET {target} = {float(target):.15f}")
    print(f"MARGIN {target - total} = {float(target - total):.15f}")

    path_only = total + 2 * (F(5, 588) - F(4, 539))
    path_expected = F(176738453, 411675264)
    assert path_only == path_expected
    assert path_only < target
    print("PATH_ONLY_CAPS tail=78; Turan k=8..3; top=42,24,12")
    print(f"PATH_ONLY_TOTAL {path_only} = {float(path_only):.15f}")
    print(
        f"PATH_ONLY_MARGIN {target - path_only} = "
        f"{float(target - path_only):.15f}"
    )
    print("STATUS PASS; finite clique replay and continuum-to-grid Lean bridges remain")


if __name__ == "__main__":
    main()
