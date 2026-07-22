#!/usr/bin/env python3
"""Exact referee for THM-2117's scalar/full-Toeplitz separation.

On the guard half-fiber, a terminal ``a*p+n*q`` becomes the shifted comb

    || |n|*beta + (a mod 2)/2 || < 1/14.

The fixed eight-comb family below passes every small clock from THM-2105,
has maximum Hunter-tree intersection mass below 1/7, and satisfies every
scalar Fourier necessary inequality from THM-2115.  Nevertheless, exact
interval arithmetic exhibits an open phase gap.  The paper proof then uses
Fejer convergence to turn that gap into a negative finite Toeplitz quadratic.

Tournament-analysis declaration
--------------------------------
Runners are not the faithful vertices here.  The relevant finite objects are
LCM-closed divisor packets (for scalar harmonics), pair-overlap edges (for
Hunter), and Fourier modes (for Toeplitz positivity).  Orienting two packets
by their scalar size produces a transitive tie-broken tournament, but destroys
simultaneous phase and therefore cannot preserve positive semidefiniteness.
The exact divisibility hypergraph and Toeplitz moment sequence are retained.
"""

from __future__ import annotations

from fractions import Fraction as F
from itertools import combinations
from math import gcd, lcm, pi, sin


FAMILY = (
    (840, 0),
    (4228, 1),
    (4984, 1),
    (5691, 1),
    (7070, 1),
    (4256, 1),
    (6195, 1),
    (4445, 1),
)
RADIUS = F(1, 14)
WIDTH = F(1, 7)
SAFE_BETA = F(223163, 6285230)
EXPECTED_GAP = (F(1103, 31115), F(176, 4949))


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def frac(x: F) -> F:
    return x - (x.numerator // x.denominator)


def circle_distance(x: F) -> F:
    y = frac(x)
    return min(y, 1 - y)


def bernoulli_two(x: F) -> F:
    y = frac(x)
    return y * y - y + F(1, 6)


def pair_mass(left: tuple[int, int], right: tuple[int, int]) -> F:
    """THM-641 for two width-1/7 windows centered at 0 or 1/2."""
    v, e = left
    w, f = right
    common = gcd(v, w)
    q, r = v // common, w // common
    a = F(e, 2) - RADIUS
    b = F(f, 2) - RADIUS
    phi = q * b - r * a
    correction = (
        -bernoulli_two(phi + q * WIDTH)
        + bernoulli_two(phi)
        + bernoulli_two(phi + (q - r) * WIDTH)
        - bernoulli_two(phi - r * WIDTH)
    )
    return WIDTH * WIDTH + correction / (2 * q * r)


def maximum_spanning_tree(
    weights: dict[tuple[int, int], F], size: int
) -> tuple[F, tuple[tuple[int, int, F], ...]]:
    parent = list(range(size))

    def root(vertex: int) -> int:
        while parent[vertex] != vertex:
            parent[vertex] = parent[parent[vertex]]
            vertex = parent[vertex]
        return vertex

    total = F(0)
    chosen: list[tuple[int, int, F]] = []
    for (i, j), weight in sorted(
        weights.items(), key=lambda item: (item[1], item[0]), reverse=True
    ):
        first, second = root(i), root(j)
        if first == second:
            continue
        parent[first] = second
        total += weight
        chosen.append((i, j, weight))
        if len(chosen) == size - 1:
            break
    require(len(chosen) == size - 1, "Kruskal tree did not span")
    return total, tuple(chosen)


def danger_intervals(v: int, e: int) -> list[tuple[F, F]]:
    """Closed interval representatives; endpoints do not affect Haar mass."""
    intervals: list[tuple[F, F]] = []
    local_radius = F(1, 14 * v)
    for integer in range(v):
        center = frac(F(2 * integer - e, 2 * v))
        left, right = center - local_radius, center + local_radius
        if left < 0:
            intervals.append((0, right))
            intervals.append((left + 1, F(1)))
        elif right > 1:
            intervals.append((left, F(1)))
            intervals.append((0, right - 1))
        else:
            intervals.append((left, right))
    return intervals


def merged_danger_union() -> tuple[tuple[tuple[F, F], ...], F]:
    intervals = sorted(
        interval
        for speed, parity in FAMILY
        for interval in danger_intervals(speed, parity)
    )
    merged: list[list[F]] = []
    for left, right in intervals:
        if not merged or left > merged[-1][1]:
            merged.append([left, right])
        elif right > merged[-1][1]:
            merged[-1][1] = right
    mass = sum((right - left for left, right in merged), F(0))
    return tuple((left, right) for left, right in merged), mass


def safe_gaps(merged: tuple[tuple[F, F], ...]) -> tuple[tuple[F, F], ...]:
    require(merged and merged[0][0] == 0 and merged[-1][1] == 1,
            "danger union must meet the circle cut")
    return tuple(
        (merged[index][1], merged[index + 1][0])
        for index in range(len(merged) - 1)
        if merged[index][1] < merged[index + 1][0]
    )


def scalar_packet_census() -> tuple[F, int, tuple[int, ...], int, int]:
    """Maximize sum(v/L) over every LCM-closed packet with >=2 speeds."""
    speeds = tuple(speed for speed, _ in FAMILY)
    records: set[tuple[int, tuple[int, ...]]] = set()
    for size in range(2, len(speeds) + 1):
        for seed in combinations(speeds, size):
            modulus = lcm(*seed)
            closure = tuple(speed for speed in speeds if modulus % speed == 0)
            records.add((modulus, closure))
    scored = [
        (sum((F(speed, modulus) for speed in closure), F(0)), modulus, closure)
        for modulus, closure in records
        if len(closure) >= 2
    ]
    require(scored, "missing multi-divisor packets")
    score, modulus, closure = max(scored)
    maximal_count = sum(value == score for value, _, _ in scored)
    return score, modulus, closure, len(records), maximal_count


def main() -> None:
    require(len(FAMILY) == 8, "rank-eight terminal family")
    require(len({speed for speed, _ in FAMILY}) == 8, "distinct speeds")
    require(all(parity in (0, 1) for _, parity in FAMILY), "binary shifts")

    # The first comb alone blinds every THM-2105 clock m=2,...,7.
    clock_rows: list[tuple[int, tuple[int, ...]]] = []
    for modulus_half in range(2, 8):
        modulus = 2 * modulus_half
        covered = {
            residue
            for residue in range(modulus)
            if any(
                (speed * residue + modulus_half * parity) % modulus == 0
                for speed, parity in FAMILY
            )
        }
        require(covered == set(range(modulus)), f"clock {modulus_half} failed")
        clock_rows.append((modulus_half, tuple(sorted(covered))))

    # Exact THM-641 pair matrix and strongest Hunter spanning tree.
    weights: dict[tuple[int, int], F] = {}
    for i, j in combinations(range(len(FAMILY)), 2):
        forward = pair_mass(FAMILY[i], FAMILY[j])
        backward = pair_mass(FAMILY[j], FAMILY[i])
        require(forward == backward, f"pair formula asymmetry at {(i, j)}")
        require(0 <= forward <= WIDTH, f"invalid pair mass at {(i, j)}")
        weights[(i, j)] = forward
    tree_mass, tree = maximum_spanning_tree(weights, len(FAMILY))
    require(tree_mass < F(1, 7), "Hunter tree unexpectedly excludes the family")

    # All frequencies with at least two divisors reduce to a finite LCM census.
    (
        packet_score,
        packet_modulus,
        packet,
        packet_count,
        maximal_packet_count,
    ) = scalar_packet_census()
    require(packet_score == F(67, 472), "unexpected maximal divisor score")
    require(packet_modulus == 49560, "unexpected maximal packet modulus")
    require(packet == (840, 6195), "unexpected maximal divisor packet")
    require(maximal_packet_count == 1, "maximal divisor packet is not unique")
    require(packet_score < F(1, 7), "multi-divisor scalar envelope is too large")
    # For a singleton packet with r=k*v, |sin(pi*k/7)|/(pi*k) is at most
    # sin(pi/7)/pi < 1/7.  For a multi-packet, |hat H(r)| is at most
    # (1/pi)*sum(v/r) <= (1/pi)*67/472 < 1/7.
    singleton_numeric = sin(pi / 7) / pi
    require(singleton_numeric < 1 / 7, "singleton analytic envelope failed")

    # Literal exact safe phase and the complete circle interval union.
    phase_rows: list[tuple[int, int, F, F]] = []
    for speed, parity in FAMILY:
        distance = circle_distance(speed * SAFE_BETA + F(parity, 2))
        margin = distance - RADIUS
        require(margin > 0, f"safe phase fails at speed {speed}")
        phase_rows.append((speed, parity, distance, margin))

    merged, danger_mass = merged_danger_union()
    gaps = safe_gaps(merged)
    largest_length = max(right - left for left, right in gaps)
    require(EXPECTED_GAP in gaps, "distinguished safe gap is missing")
    require(EXPECTED_GAP[1] - EXPECTED_GAP[0] == largest_length,
            "distinguished safe gap is not maximal")
    largest_gap = EXPECTED_GAP
    require(SAFE_BETA == sum(EXPECTED_GAP, F(0)) / 2, "safe beta is not midpoint")
    require(danger_mass < 1, "exact danger union unexpectedly covers")

    print("LRC14 SIGNED HALF-FIBER SCALAR/TOEPLITZ SEPARATION -- exact audit")
    print(f"family (speed, parity): {FAMILY}")
    print("small clocks m=2..7: PASS (speed 840 is universal)")
    for modulus_half, residues in clock_rows:
        print(f"  m={modulus_half}: {residues}")
    print(f"maximum Hunter tree mass: {tree_mass}")
    print(f"Hunter slack 1/7-tree: {F(1, 7) - tree_mass}")
    print("Hunter tree edges (indices, exact mass):")
    for i, j, weight in tree:
        print(f"  ({i},{j}) {weight}")
    print(f"LCM closures audited: {packet_count}")
    print(
        "largest multi-divisor score: "
        f"{packet_score} at r={packet_modulus}, packet={packet}"
    )
    print(f"multi-divisor gap to 1/7: {F(1, 7) - packet_score}")
    print(f"singleton analytic envelope sin(pi/7)/pi: {singleton_numeric:.15f} < 1/7")
    print(f"strict-safe beta: {SAFE_BETA}")
    for speed, parity, distance, margin in phase_rows:
        print(f"  speed={speed} parity={parity}: distance={distance}, margin={margin}")
    print(f"merged danger components: {len(merged)}")
    print(f"safe components: {len(gaps)}")
    print(f"danger union mass: {danger_mass}")
    print(f"safe mass: {1 - danger_mass}")
    print(
        "largest safe gap: "
        f"({largest_gap[0]}, {largest_gap[1]}), length={largest_gap[1]-largest_gap[0]}"
    )
    print("full Toeplitz separation: Fejer convergence on the open safe gap")
    print("PASS")


if __name__ == "__main__":
    main()
