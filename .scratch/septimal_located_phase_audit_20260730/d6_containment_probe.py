#!/usr/bin/env python3
"""Exact hostile probe for the D=168,c=4 lone d=6 located-phase gate.

For a reduced numerator c in {1,5}, let V_c be the phases at which the
strict d=6 mask selects one residue class (the high spike bit).  A proposed
located-phase closure would need to rule out R_{a,b} subset V_c uniformly.
Equivalently, the compact low set C_c=T\V_c must not be contained in the
open union of the two aligned danger combs D_a union D_b.  We test this
containment exactly on the complete rational endpoint partition.
"""

from fractions import Fraction as Q
from math import gcd


def norm(value):
    value %= 1
    return min(value, 1 - value)


def in_danger(a, u):
    return norm(a * u) < Q(1, 14)


def selected_count(d, numerator, u):
    return sum(
        norm(Q(numerator) * (residue + u) / d) < Q(1, 14)
        for residue in range(d)
    )


def boundaries_for_danger(a):
    points = {Q(0), Q(1)}
    for integer in range(a + 1):
        for sign in (-1, 1):
            point = (Q(integer) + Q(sign, 14)) / a
            if 0 <= point <= 1:
                points.add(point)
    return points


def boundaries_for_mask(d, numerator):
    points = {Q(0), Q(1)}
    for residue in range(d):
        for integer in range(-1, numerator + 2):
            for sign in (-1, 1):
                point = Q(d) * (Q(integer) + Q(sign, 14)) / numerator - residue
                if 0 <= point <= 1:
                    points.add(point)
    return points


def low_covered(a, b, numerator):
    points = sorted(
        boundaries_for_danger(a)
        | boundaries_for_danger(b)
        | boundaries_for_mask(6, numerator)
    )
    samples = set(points[:-1])
    samples.update((x + y) / 2 for x, y in zip(points, points[1:]))
    hostile = []
    for u in sorted(samples):
        low = selected_count(6, numerator, u) == 0
        covered = in_danger(a, u) or in_danger(b, u)
        if low and not covered:
            hostile.append(u)
    return not hostile, tuple(hostile[:4])


def low_intersection_measure(a, numerator):
    """Exact Haar measure of C_numerator intersect D_a."""
    answer = Q(0)
    radius = Q(1, 14 * a)
    for j in range(numerator):
        left = (Q(j) + Q(3, 7)) / numerator
        right = (Q(j) + Q(4, 7)) / numerator
        # Include the two wrap representatives k=0,a.
        for k in range(a + 1):
            center = Q(k, a)
            overlap = min(right, center + radius) - max(left, center - radius)
            if overlap > 0:
                answer += overlap
    return answer


def main():
    first = None
    witnesses = []
    max_aligned = 30
    for numerator in (1, 5):
        for b in range(2, max_aligned + 1):
            for a in range(1, b):
                covered, hostile = low_covered(a, b, numerator)
                if covered:
                    witness = (numerator, a, b)
                    witnesses.append(witness)
                    if first is None:
                        first = witness
    print("d6 located-phase exact containment probe")
    print(f"universe=numerator in {{1,5}}, 1<=a<b<={max_aligned}")
    print(f"containment_witness_count={len(witnesses)}")
    print(f"first_witness={first}")
    print(f"first_20={witnesses[:20]}")
    for numerator in (1, 5):
        print(
            f"numerator_{numerator}_first=",
            next((w for w in witnesses if w[0] == numerator), None),
        )
    for numerator in (1, 5):
        ranked = sorted(
            ((low_intersection_measure(a, numerator), a) for a in range(1, 1001)),
            reverse=True,
        )
        print(f"numerator_{numerator}_top_intersections={ranked[:20]}")
        distinct_pair_bound = ranked[0][0] + next(
            mass for mass, a in ranked if a != ranked[0][1]
        )
        print(
            f"numerator_{numerator}_top_distinct_pair_sum={distinct_pair_bound} "
            f"target={Q(1,7)}"
        )


if __name__ == "__main__":
    main()
