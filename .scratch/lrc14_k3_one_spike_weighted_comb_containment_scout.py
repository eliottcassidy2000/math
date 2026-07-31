#!/usr/bin/env python3
"""Exact hostile scout for the positioned c=3 optional-spike event.

For an aligned triple A let R_A be its compact safe set.  For a spike
denominator d=7a+r and reduced numerator c, the high-count bit is one on

    B(c, parity(a), r)
      = {u : ||c*u + a/2|| < r/14}.

A non-guaranteed one-spike cover requires R_A subset B.  This scout checks
that compact/open containment exactly from rational boundary cells.  It is
a finite hostile probe, not a uniform theorem in unbounded A or c.

There is also a rigorous analytic tail sidecar.  The complement of B has c
components and total mass delta=(7-r)/7.  Applying the sharp one-comb
interval discrepancy to every component shows that containment forces

    c * sum_{a in A} 1/a >= 2(7-r)/3.

This bounds genuinely high-scale triples but does not make the remaining
projective ratio bank finite by itself.
"""

from argparse import ArgumentParser
from bisect import bisect_right
from fractions import Fraction as Q
from itertools import combinations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def norm(value):
    residue = value % 1
    return min(residue, 1 - residue)


def danger(multiplier, u):
    return norm(multiplier * u) < Q(1, 14)


def high(numerator, parity, remainder, u):
    center = Q(parity, 2)
    return norm(numerator * u + center) < Q(remainder, 14)


def equality_boundaries(multiplier, radius, center=Q(0)):
    """Boundaries of ||multiplier*u+center||<radius in [0,1]."""
    points = {Q(0), Q(1)}
    for integer in range(-1, multiplier + 2):
        for sign in (-1, 1):
            u = (Q(integer) - center + sign * radius) / multiplier
            if 0 <= u <= 1:
                points.add(u)
    return tuple(sorted(points))


def safe_atoms(A):
    boundaries = tuple(
        sorted(
            {Q(0), Q(1)}
            | {
                point
                for multiplier in A
                for point in equality_boundaries(multiplier, Q(1, 14))
            }
        )
    )
    safe_points = tuple(
        point
        for point in boundaries
        if point < 1 and not any(danger(a, point) for a in A)
    )
    safe_cells = tuple(
        (left, right)
        for left, right in zip(boundaries, boundaries[1:])
        if left < right
        and not any(danger(a, (left + right) / 2) for a in A)
    )
    require(safe_points or safe_cells, ("aligned triple has empty safe set", A))
    return safe_points, safe_cells


def containment_witness(atoms, numerator, parity, remainder):
    """Return None for containment, otherwise an exact safe point outside B."""
    safe_points, safe_cells = atoms
    for point in safe_points:
        if not high(numerator, parity, remainder, point):
            return point
    boundaries = equality_boundaries(
        numerator,
        Q(remainder, 14),
        Q(parity, 2),
    )
    for left, right in safe_cells:
        # A strict B-boundary in the interior is itself safe and outside B.
        index = bisect_right(boundaries, left)
        if index < len(boundaries) and boundaries[index] < right:
            return boundaries[index]
        midpoint = (left + right) / 2
        if not high(numerator, parity, remainder, midpoint):
            return midpoint
    return None


def direct_partition_control(A, numerator, parity, remainder):
    """Compare the atom test with a literal common-boundary sweep."""
    boundaries = tuple(
        sorted(
            {Q(0), Q(1)}
            | {
                point
                for multiplier in A
                for point in equality_boundaries(multiplier, Q(1, 14))
            }
            | set(
                equality_boundaries(
                    numerator,
                    Q(remainder, 14),
                    Q(parity, 2),
                )
            )
        )
    )
    samples = set(boundaries[:-1])
    samples.update(
        (left + right) / 2
        for left, right in zip(boundaries, boundaries[1:])
        if left < right
    )
    literal = all(
        any(danger(a, point) for a in A)
        or high(numerator, parity, remainder, point)
        for point in samples
    )
    atom = containment_witness(safe_atoms(A), numerator, parity, remainder) is None
    require(literal == atom, ("atom/partition disagreement", A, numerator, parity, remainder))
    return 1


def main(A_max, c_max):
    cases = ((0, 5), (0, 6), (1, 5), (1, 6))
    containments = {case: [] for case in cases}
    triples = 0
    pair_tests = 0
    partition_controls = 0
    minimum_witness = {}
    for A in combinations(range(1, A_max + 1), 3):
        atoms = safe_atoms(A)
        triples += 1
        for numerator in range(1, c_max + 1):
            for case in cases:
                parity, remainder = case
                witness = containment_witness(
                    atoms, numerator, parity, remainder
                )
                pair_tests += 1
                if witness is None:
                    containments[case].append((A, numerator))
                elif case not in minimum_witness:
                    minimum_witness[case] = (A, numerator, witness)
                if A[-1] <= 6 and numerator <= 6:
                    partition_controls += direct_partition_control(
                        A, numerator, parity, remainder
                    )

    # Boundary sensitivity: R_{1,2,3} is contained in the *closed* wide
    # d=13,c=1 band but not in its strict-open version because u=1/14 stays.
    require(
        containment_witness(safe_atoms((1, 2, 3)), 1, 1, 6) == Q(1, 14),
        "strict shifted-boundary control changed",
    )
    # A second comb can cover those boundary points, so shifted containments
    # are genuinely possible and must not be silently generalized away.
    require(
        containment_witness(safe_atoms((1, 2, 14)), 1, 1, 6) is None,
        "positive shifted containment control changed",
    )

    minimum_tariff = {}
    for case in cases:
        _parity, remainder = case
        tariff_floor = Q(2 * (7 - remainder), 3)
        for A, numerator in containments[case]:
            tariff = numerator * sum(Q(1, a) for a in A)
            require(
                tariff >= tariff_floor,
                ("analytic containment tariff failed", case, A, numerator, tariff),
            )
            if case not in minimum_tariff or tariff < minimum_tariff[case]:
                minimum_tariff[case] = tariff

    print("LRC14 k3 one-spike positioned weighted-comb containment scout")
    print(f"A_max={A_max}")
    print(f"c_max={c_max}")
    print(f"aligned_triples={triples}")
    print(f"case_tests={pair_tests}")
    print(f"independent_partition_controls={partition_controls}")
    for case in cases:
        print(f"case={case},containment_count={len(containments[case])}")
        print(f"case={case},first_containments={containments[case][:20]}")
        print(f"case={case},minimum_noncontainment_witness={minimum_witness.get(case)}")
        print(
            f"case={case},analytic_tariff_floor={Q(2 * (7 - case[1]), 3)},"
            f"minimum_containment_tariff={minimum_tariff.get(case)}"
        )
    print("scope=finite hostile scout only; A and numerator are unbounded")
    print("all_exact_controls=PASS")


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--A-max", type=int, default=20)
    parser.add_argument("--c-max", type=int, default=20)
    args = parser.parse_args()
    main(args.A_max, args.c_max)
