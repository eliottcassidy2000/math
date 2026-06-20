#!/usr/bin/env python3
"""Euler parity carriers for the LRC14 collar route.

The prompt asks to place two equinumerosities beside the LRC14 proof:

  * odd partitions = distinct partitions:
        prod_n (1+x^n) = prod_{m odd} 1/(1-x^m)
    This is a 2-adic carry theorem: a distinct part 2^a*m contributes
    2^a copies of the odd part m.

  * tournaments and even graphs are equinumerous.  At the labelled suspension
    level, an arbitrary graph on n-1 vertices is completed to an even graph on
    n vertices by adding exactly the apex edges needed to make all degrees even.

Both say the same structural thing: free binary choices plus parity completion
can be moved into a constrained carrier.  This script audits that carrier on
the exact LRC14 AP-collar rows already isolated by THM-541/543/544.
"""

from __future__ import annotations

from collections import Counter
from fractions import Fraction
from functools import lru_cache
from itertools import combinations
from math import comb


N = 80
TARGET_DENOM = 14
AP_SECOND = Fraction(426, 35035)
DROP6_COLLAR = Fraction(7, 858)
DROP6_MOUTHS = (
    (Fraction(29, 182), Fraction(9, 56)),
    (Fraction(29, 168), Fraction(27, 154)),
    (Fraction(127, 154), Fraction(139, 168)),
    (Fraction(47, 56), Fraction(153, 182)),
)


def multiply_distinct(coeffs: list[int], part: int) -> list[int]:
    out = coeffs[:]
    for n in range(part, len(coeffs)):
        out[n] += coeffs[n - part]
    return out


def distinct_coeffs(n: int) -> list[int]:
    coeffs = [1] + [0] * n
    for part in range(1, n + 1):
        coeffs = multiply_distinct(coeffs, part)
    return coeffs


def odd_coeffs(n: int) -> list[int]:
    coeffs = [1] + [0] * n
    for part in range(1, n + 1, 2):
        for total in range(part, n + 1):
            coeffs[total] += coeffs[total - part]
    return coeffs


def distinct_partitions(n: int, max_part: int | None = None) -> list[tuple[int, ...]]:
    if n == 0:
        return [tuple()]
    if max_part is None:
        max_part = n
    out: list[tuple[int, ...]] = []
    for part in range(min(n, max_part), 0, -1):
        for tail in distinct_partitions(n - part, min(part - 1, n - part)):
            out.append((part,) + tail)
    return out


def odd_partitions(n: int, max_part: int | None = None) -> list[tuple[int, ...]]:
    if n == 0:
        return [tuple()]
    if max_part is None:
        max_part = n if n % 2 else n - 1
    out: list[tuple[int, ...]] = []
    start = min(n, max_part)
    if start % 2 == 0:
        start -= 1
    for part in range(start, 0, -2):
        for tail in odd_partitions(n - part, min(part, n - part)):
            out.append((part,) + tail)
    return out


def odd_part_and_power(n: int) -> tuple[int, int]:
    power = 0
    while n % 2 == 0:
        n //= 2
        power += 1
    return n, power


def distinct_to_odd(partition: tuple[int, ...]) -> tuple[int, ...]:
    odd_parts: list[int] = []
    for part in partition:
        odd, power = odd_part_and_power(part)
        odd_parts.extend([odd] * (1 << power))
    return tuple(sorted(odd_parts, reverse=True))


def odd_to_distinct(partition: tuple[int, ...]) -> tuple[int, ...]:
    counts = Counter(partition)
    parts: list[int] = []
    for odd, multiplicity in counts.items():
        bit = 0
        while multiplicity:
            if multiplicity & 1:
                parts.append(odd * (1 << bit))
            multiplicity >>= 1
            bit += 1
    return tuple(sorted(parts, reverse=True))


def verify_euler_bijection(n: int) -> tuple[int, int, bool, bool]:
    distinct = distinct_partitions(n)
    odd = odd_partitions(n)
    image = {distinct_to_odd(part) for part in distinct}
    inverse_image = {odd_to_distinct(part) for part in odd}
    return len(distinct), len(odd), image == set(odd), inverse_image == set(distinct)


def suspension_even_graph(base_bits: int, n: int) -> tuple[tuple[int, int], ...]:
    """Complete a graph on vertices 0..n-2 to an even graph on 0..n-1."""

    base_edges = [(i, j) for i in range(n - 1) for j in range(i + 1, n - 1)]
    edges: set[tuple[int, int]] = set()
    degree = [0] * n
    for bit, (i, j) in enumerate(base_edges):
        if (base_bits >> bit) & 1:
            edges.add((i, j))
            degree[i] += 1
            degree[j] += 1
    apex = n - 1
    for v in range(n - 1):
        if degree[v] % 2:
            edges.add((v, apex))
            degree[v] += 1
            degree[apex] += 1
    assert all(d % 2 == 0 for d in degree)
    return tuple(sorted(edges))


def verify_even_graph_suspension(n: int) -> tuple[int, int, bool]:
    seen = set()
    base_count = 1 << comb(n - 1, 2)
    for bits in range(base_count):
        seen.add(suspension_even_graph(bits, n))
    return base_count, len(seen), len(seen) == base_count


def merge(intervals: list[tuple[Fraction, Fraction]]) -> list[tuple[Fraction, Fraction]]:
    out: list[tuple[Fraction, Fraction]] = []
    for lo, hi in sorted(intervals):
        if lo >= hi:
            continue
        if out and lo <= out[-1][1]:
            if hi > out[-1][1]:
                out[-1] = (out[-1][0], hi)
        else:
            out.append((lo, hi))
    return out


@lru_cache(maxsize=None)
def danger_arcs(speed: int) -> tuple[tuple[Fraction, Fraction], ...]:
    arcs: list[tuple[Fraction, Fraction]] = []
    denom = TARGET_DENOM * speed
    for tooth in range(speed):
        left = TARGET_DENOM * tooth - 1
        right = TARGET_DENOM * tooth + 1
        if left < 0:
            arcs.append((Fraction(0), Fraction(right, denom)))
            arcs.append((Fraction(denom + left, denom), Fraction(1)))
        else:
            arcs.append((Fraction(left, denom), Fraction(right, denom)))
    return tuple(arcs)


def safe_components(core: tuple[int, ...]) -> tuple[tuple[Fraction, Fraction], ...]:
    danger: list[tuple[Fraction, Fraction]] = []
    for speed in core:
        danger.extend(danger_arcs(speed))
    safe: list[tuple[Fraction, Fraction]] = []
    cursor = Fraction(0)
    for lo, hi in merge(danger):
        if lo > cursor:
            safe.append((cursor, lo))
        if hi > cursor:
            cursor = hi
    if cursor < 1:
        safe.append((cursor, Fraction(1)))
    return tuple(safe)


def measure(intervals: tuple[tuple[Fraction, Fraction], ...]) -> Fraction:
    return sum((hi - lo for lo, hi in intervals), Fraction(0))


def intersect_measure(
    xs: tuple[tuple[Fraction, Fraction], ...],
    ys: tuple[tuple[Fraction, Fraction], ...],
) -> Fraction:
    total = Fraction(0)
    i = j = 0
    while i < len(xs) and j < len(ys):
        lo = max(xs[i][0], ys[j][0])
        hi = min(xs[i][1], ys[j][1])
        if lo < hi:
            total += hi - lo
        if xs[i][1] < ys[j][1]:
            i += 1
        else:
            j += 1
    return total


def odd_carry_profile(core: tuple[int, ...]) -> dict[int, int]:
    profile: Counter[int] = Counter()
    for speed in core:
        odd, power = odd_part_and_power(speed)
        profile[odd] += 1 << power
    return dict(sorted(profile.items()))


def profile_delta(core: tuple[int, ...], base: tuple[int, ...]) -> dict[int, int]:
    current = Counter(odd_carry_profile(core))
    reference = Counter(odd_carry_profile(base))
    keys = sorted(set(current) | set(reference))
    return {key: current[key] - reference[key] for key in keys if current[key] != reference[key]}


def fmt_profile(profile: dict[int, int]) -> str:
    return "{" + ", ".join(f"{k}:{v}" for k, v in profile.items()) + "}"


def main() -> None:
    print("LRC14 Euler parity carrier scout")
    print("exact arithmetic only")
    print()

    distinct = distinct_coeffs(N)
    odd = odd_coeffs(N)
    assert distinct == odd
    print("Euler product identity")
    print(f"  verified coefficients through x^{N}: prod(1+x^n) == prod_odd 1/(1-x^m)")
    print("  first coefficients:")
    print("   " + ", ".join(f"p_odd_distinct({i})={distinct[i]}" for i in range(15)))
    print()

    print("binary carry bijection samples")
    for n in [14, 20, 30]:
        d_count, o_count, forward_ok, backward_ok = verify_euler_bijection(n)
        print(
            f"  n={n}: distinct={d_count}, odd={o_count}, "
            f"forward={forward_ok}, backward={backward_ok}"
        )
    sample = (20, 10, 8, 3, 1)
    print(f"  sample distinct {sample} -> odd {distinct_to_odd(sample)}")
    print(f"  inverse -> {odd_to_distinct(distinct_to_odd(sample))}")
    print()

    print("labelled tournament/even-graph suspension")
    print("  arbitrary graph on n-1 vertices -> even graph on n vertices by apex parity completion")
    for n in range(3, 8):
        free, image, injective = verify_even_graph_suspension(n)
        print(
            f"  n={n}: free_bits=C({n-1},2), count={free}, "
            f"even_graphs_with_chosen_apex={image}, injective={injective}"
        )
    print("  formula for all n: count = 2^C(n-1,2), with inverse deleting the chosen apex")
    print()

    drop6 = tuple(v for v in range(1, 14) if v != 6)
    rows = [
        ("THM-541 drop-6 collar", drop6),
        ("AP second drop-12", tuple(v for v in range(1, 14) if v != 12)),
        ("THM-543 one-tail exception", (1, 2, 3, 4, 5, 7, 8, 9, 11, 12, 13, 20)),
        ("THM-544 corrected two-tail minimum", (1, 2, 3, 5, 7, 8, 9, 11, 12, 13, 20, 46)),
        ("three-scale HYP-2655 witness shadow", (1, 2, 30, 31, 32, 60, 61, 62, 63, 64, 65, 66)),
    ]

    print("LRC14 AP-collar rows as odd-part carry profiles")
    print("  profile records odd_part -> binary multiplicity sum from speed=2^a*odd_part")
    for name, core in rows:
        comps = safe_components(tuple(sorted(core)))
        safe = measure(comps)
        old = intersect_measure(comps, DROP6_MOUTHS)
        delta = profile_delta(tuple(sorted(core)), drop6)
        new_odd = sorted(set(odd_carry_profile(core)) - set(odd_carry_profile(drop6)))
        print(f"  {name}")
        print(f"    core={tuple(sorted(core))}")
        print(f"    safe={safe}, old_mouth_survivor={old}, above_second={safe >= AP_SECOND}")
        print(f"    odd_carry_profile={fmt_profile(odd_carry_profile(tuple(sorted(core))))}")
        print(f"    delta_vs_drop6={fmt_profile(delta)}, new_odd_parts={new_odd}")
    print()

    print("carrier conclusion")
    print("  Euler odd/distinct: binary free choices become odd-shell multiplicities.")
    print("  even-graph suspension: binary free choices become parity-completed cycle carriers.")
    print("  HYP-2656 CRT/halving: the odd base is rigid; dyadic levels refine it.")
    print("  LRC14 reading: bounded AP-tail rows should be routed by odd-shell carry plus")
    print("  exact mouth survival before using scalar measure. The unique one-tail")
    print("  below-second exception is a pure carry inside an existing odd shell.")
    print("  The corrected two-tail minimum introduces a new odd shell and damages")
    print("  the old mouth, but already pays the second-value threshold exactly.")
    print()

    print("Tournament Analysis")
    print("  vertices: parity carrier lenses")
    print("  pairwise observable: which lens preserves exact safe-measure routing")
    print("  switch/gauge: carry profile before scalarizing to safe measure")
    print("  Hamiltonian path:")
    print("    odd_shell_carry > mouth_survivor > state_word_damage")
    print("    > apex_parity_completion > raw_speed_set")
    print("  directed 3-cycles: 0 (declared proof-lens order)")
    print()
    print("PASS: Euler parity carrier scout complete.")


if __name__ == "__main__":
    main()
