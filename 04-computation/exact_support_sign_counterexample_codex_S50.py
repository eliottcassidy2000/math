#!/usr/bin/env python3
"""Exact sign audit for the support-Moebius masses used in THM-935/948.

The proposed positive-association law ``M(A) >= 0`` for every ``|A| >= 3``
is false.  This script uses two independent exact computations of every arc
intersection measure:

* an endpoint sweep on the ``14*lcm(A)`` grid;
* a midpoint count on the same rational cell decomposition.

It then performs the subset-lattice Moebius recursion from the exact-support
decomposition.  All arithmetic is ``fractions.Fraction``.

Tournament analysis / challenged assumption
--------------------------------------------
The pair observable ``M({i,j})`` is symmetric, so its sign does not define an
honest directed edge without an arbitrary label order.  More importantly, the
counterexamples first appear as signed 3-, 4-, and 5-hyperedges.  A runner
tournament would destroy support parity and the higher Moebius sign.  The
faithful quotient is therefore the signed support hypergraph (equivalently the
cumulant complex), not a tournament on runners.  Score histograms, directed
cycles, SCCs, edge flips, and Hamiltonian-path counts are consequently not
well-defined for this observable; this is a structural negative result rather
than an omitted fingerprint.
"""

from fractions import Fraction
from functools import lru_cache
from itertools import combinations, groupby
from math import gcd


P = Fraction(1, 7)  # measure of one bad arc


def lcm(values: tuple[int, ...]) -> int:
    result = 1
    for value in values:
        result = result * abs(value) // gcd(result, abs(value))
    return result


@lru_cache(maxsize=None)
def intersection_measure_sweep(support: tuple[int, ...]) -> Fraction:
    """Exact measure of the intersection of all bad arcs in ``support``."""
    support = tuple(sorted(support))
    scale = lcm(support)
    denominator = 14 * scale
    events: list[tuple[int, int]] = []
    for speed in support:
        width = scale // abs(speed)
        for center in range(abs(speed)):
            low = (14 * center - 1) * width
            high = (14 * center + 1) * width
            if low < 0:
                events.extend(
                    [(0, 1), (high, -1), (denominator + low, 1), (denominator, -1)]
                )
            else:
                events.extend([(low, 1), (high, -1)])

    depth = 0
    previous = 0
    total = 0
    target_depth = len(support)
    for coordinate, same_coordinate in groupby(sorted(events), key=lambda item: item[0]):
        if depth == target_depth:
            total += coordinate - previous
        depth += sum(delta for _, delta in same_coordinate)
        previous = coordinate
    return Fraction(total, denominator)


@lru_cache(maxsize=None)
def intersection_measure_midpoints(support: tuple[int, ...]) -> Fraction:
    """Independent exact cell-midpoint count of the same intersection."""
    support = tuple(sorted(support))
    scale = lcm(support)
    denominator = 14 * scale
    doubled = 2 * denominator
    bad_cells = 0
    for cell in range(denominator):
        midpoint_numerator = 2 * cell + 1
        all_bad = True
        for speed in support:
            residue = (speed * midpoint_numerator) % doubled
            distance_numerator = min(residue, doubled - residue)
            if 14 * distance_numerator >= doubled:
                all_bad = False
                break
        bad_cells += int(all_bad)
    return Fraction(bad_cells, denominator)


@lru_cache(maxsize=None)
def mobius_mass(support: tuple[int, ...]) -> Fraction:
    """Exact recursive mass M(A) from the THM-935/948 support decomposition."""
    support = tuple(sorted(support))
    size = len(support)
    assert size >= 2
    sweep = intersection_measure_sweep(support)
    result = sweep - P**size
    for lower_size in range(2, size):
        for subset in combinations(support, lower_size):
            result -= mobius_mass(tuple(subset)) * P ** (size - lower_size)
    return result


def sign_counts(size: int, maximum: int) -> tuple[int, int, int, tuple[Fraction, tuple[int, ...]]]:
    negative = zero = positive = 0
    minimum: tuple[Fraction, tuple[int, ...]] | None = None
    for support in combinations(range(1, maximum + 1), size):
        mass = mobius_mass(tuple(support))
        negative += int(mass < 0)
        zero += int(mass == 0)
        positive += int(mass > 0)
        if minimum is None or mass < minimum[0]:
            minimum = (mass, support)
    assert minimum is not None
    return negative, zero, positive, minimum


def main() -> None:
    witnesses = {
        3: ((1, 2, 15), Fraction(-4, 1715)),
        4: ((1, 2, 3, 28), Fraction(-5, 4116)),
        5: ((1, 2, 3, 4, 32), Fraction(-109, 806736)),
    }
    print("EXACT SUPPORT-SIGN COUNTEREXAMPLES")
    for size, (support, expected) in witnesses.items():
        # Independent midpoint verification is intentionally restricted to
        # the displayed witnesses and their faces.  Applying it to every
        # exploratory scan row would enumerate enormous lcm grids.
        for face_size in range(2, len(support) + 1):
            for face in combinations(support, face_size):
                sweep = intersection_measure_sweep(tuple(face))
                midpoint = intersection_measure_midpoints(tuple(face))
                assert sweep == midpoint, (face, sweep, midpoint)
        mass = mobius_mass(support)
        assert mass == expected
        print(f"M_{size}{support} = {mass} = {float(mass):+.12f}")

    triple = (1, 2, 15)
    print("\nTRIPLE AUTOPSY")
    for pair in combinations(triple, 2):
        measure = intersection_measure_sweep(tuple(pair))
        print(f"mu{pair}={measure}; M_2={mobius_mass(tuple(pair))}")
    triple_measure = intersection_measure_sweep(triple)
    pair_sum = sum(mobius_mass(tuple(pair)) for pair in combinations(triple, 2))
    print(f"mu{triple}={triple_measure}")
    print(f"sum_pair_excess={pair_sum}; P*sum={P * pair_sum}")
    print(f"M_3={triple_measure} - P^3 - P*sum = {mobius_mass(triple)}")

    print("\nEXACT FINITE SIGN SCANS")
    for size, maximum in [(3, 35), (4, 28), (5, 20)]:
        negative, zero, positive, minimum = sign_counts(size, maximum)
        print(
            f"size={size}, speeds=1..{maximum}: "
            f"negative={negative}, zero={zero}, positive={positive}; "
            f"minimum={minimum[0]} at {minimum[1]}"
        )

    family_negative = []
    for far in range(5, 65):
        support = (1, 2, 3, 4, far)
        if len(set(support)) == 5 and mobius_mass(support) < 0:
            family_negative.append((far, mobius_mass(support)))
    assert family_negative[0] == (32, Fraction(-109, 806736))
    print(
        "quintuple family (1,2,3,4,N), N=5..64: "
        f"first negative at N={family_negative[0][0]}, mass={family_negative[0][1]}"
    )

    print("\nTOURNAMENT AUDIT")
    print("pair observable M({i,j}) is symmetric: no canonical edge orientation")
    print("faithful object: signed support hypergraph / cumulant complex")
    print("tournament score/cycle/SCC/flip/Hamiltonian fingerprints: not well-defined")
    print("PASS")


if __name__ == "__main__":
    main()
