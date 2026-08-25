#!/usr/bin/env python3
"""Standalone exact audit of THM-4030's d=4 affine defect-lattice boundary.

All arithmetic is exact.
The two independent decision paths are:

1. direct wall-cell subdivision in the pack phase y;
2. finite enumeration of the affine centre lattice (A,B,C).
"""

from __future__ import annotations

from fractions import Fraction as F
from math import gcd
import sys


sys.stdout.reconfigure(newline="\n")

CHECKS = 0


def require(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(label)


def distance(speed: int, phase: F) -> F:
    residue = (speed * phase) % 1
    return min(residue, 1 - residue)


def bad_label_sets(y: F, even: int, first: int, second: int) -> tuple[tuple[int, ...], ...]:
    return tuple(
        tuple(
            label
            for label in range(4)
            if distance(speed, (y + label) / 4) < F(1, 14)
        )
        for speed in (even, first, second)
    )


def labels_are_covered(y: F, even: int, first: int, second: int) -> bool:
    bad = bad_label_sets(y, even, first, second)
    return set().union(*map(set, bad)) == set(range(4))


def direct_covering_phase(even: int, first: int, second: int) -> F | None:
    """Decide phase-cover existence by exact wall-cell subdivision."""
    walls = {F(0), F(1)}
    for speed in (even, first, second):
        for label in range(4):
            # On 0<=y<=1 and 0<=label<=3, speed*(y+label)/4 lies
            # between 0 and speed.  The padded range is therefore complete.
            for integer in range(-1, speed + 2):
                for sign in (-1, 1):
                    wall = F(4, speed) * (F(integer) + sign * F(1, 14)) - label
                    if 0 <= wall <= 1:
                        walls.add(wall)
    ordered = sorted(walls)
    probes = set(ordered)
    probes.update((left + right) / 2 for left, right in zip(ordered, ordered[1:]))
    for y in sorted(probes):
        if labels_are_covered(y, even, first, second):
            return y
    return None


def defects(r: int, first: int, second: int, A: int, B: int, C: int) -> tuple[int, int, int]:
    n12 = 2 * second * A - first * (2 * B - second)
    n1e = 4 * r * A - first * (2 * C - r)
    n2e = 4 * r * B - 2 * second * C - second * r
    return n12, n1e, n2e


def certificate_point(
    r: int,
    first: int,
    second: int,
    A: int,
    B: int,
    C: int,
) -> F | None:
    """Return a point in the three open real intervals, if one exists."""
    centres = (F(A, first), F(B, second) - F(1, 2), F(C, 2 * r) - F(1, 4))
    radii = (F(1, 14 * first), F(1, 14 * second), F(1, 28 * r))
    left = max(centre - radius for centre, radius in zip(centres, radii))
    right = min(centre + radius for centre, radius in zip(centres, radii))
    if left < right:
        return (left + right) / 2
    return None


def oriented_lattice_certificate(r: int, first: int, second: int):
    """Finite exact decision for first@z, second@(z+1/2), 2r@(z+1/4).

    Simultaneous integral translation normalizes 0<=A<first.  If the
    intervals meet, pairwise radius bounds force 0<=B<=2*second and
    0<=C<=3*r, so the displayed box is complete, not a heuristic cutoff.
    """
    require(r > 0 and first > 0 and second > 0, "positive speeds")
    require(r % 2 == first % 2 == second % 2 == 1, "odd d=4 parameters")
    for A in range(first):
        for B in range(2 * second + 1):
            for C in range(3 * r + 1):
                n12, n1e, n2e = defects(r, first, second, A, B, C)
                if not 7 * abs(n12) < first + second:
                    continue
                if not 7 * abs(n1e) < first + 2 * r:
                    continue
                if not 7 * abs(n2e) < second + 2 * r:
                    continue
                point = certificate_point(r, first, second, A, B, C)
                require(point is not None, "open-interval Helly reconstruction")
                require(2 * r * n12 - second * n1e + first * n2e == 0, "defect identity")
                require(n12 % 2 and n1e % 2 and n2e % 2, "defects are odd")
                require(n12 % gcd(first, second) == 0, "gcd(first,second) divides defect")
                require(n1e % gcd(first, r) == 0, "gcd(first,r) divides defect")
                require(n2e % gcd(second, r) == 0, "gcd(second,r) divides defect")
                y = (4 * point) % 1
                require(labels_are_covered(y, 2 * r, first, second), "certificate reconstructs labels")
                return (A, B, C, n12, n1e, n2e, point, y)
    return None


def lattice_certificate(even: int, first: int, second: int):
    require(even % 4 == 2, "unique even exception has order two")
    r = even // 2
    forward = oriented_lattice_certificate(r, first, second)
    if forward is not None:
        return ("first_then_second", forward)
    reverse = oriented_lattice_certificate(r, second, first)
    if reverse is not None:
        return ("second_then_first", reverse)
    return None


def oriented_defect_certificate(r: int, first: int, second: int):
    """Eliminate A,B,C: bounded odd defects plus one relation are exact."""
    bound12 = (first + second - 1) // 7
    bound1e = (first + 2 * r - 1) // 7
    bound2e = (second + 2 * r - 1) // 7
    for n12 in range(-bound12, bound12 + 1):
        if n12 == 0 or n12 % 2 == 0 or n12 % gcd(first, second):
            continue
        for n1e in range(-bound1e, bound1e + 1):
            if n1e == 0 or n1e % 2 == 0 or n1e % gcd(first, r):
                continue
            numerator = second * n1e - 2 * r * n12
            if numerator % first:
                continue
            n2e = numerator // first
            if not -bound2e <= n2e <= bound2e:
                continue
            if n2e == 0 or n2e % 2 == 0 or n2e % gcd(second, r):
                continue

            # Construct A,B,C from the defect triple.  The local-CRT lemma in
            # the proof note guarantees an A; the finite residue loop is an
            # independent exact realization rather than an assumption.
            k12 = (n12 - first * second) // 2
            k1e = (n1e - first * r) // 2
            k2e = (n2e + second * r) // 2
            for A in range(first):
                if (second * A - k12) % first:
                    continue
                if (2 * r * A - k1e) % first:
                    continue
                B = (second * A - k12) // first
                C = (2 * r * A - k1e) // first
                require(2 * r * B - second * C == k2e, "defect CRT third equation")
                require(defects(r, first, second, A, B, C) == (n12, n1e, n2e), "defect realization")
                return (n12, n1e, n2e, A, B, C)
            raise AssertionError(("unrealized defect triple", r, first, second, n12, n1e, n2e))
    return None


def defect_certificate(even: int, first: int, second: int):
    r = even // 2
    forward = oriented_defect_certificate(r, first, second)
    if forward is not None:
        return ("first_then_second", forward)
    reverse = oriented_defect_certificate(r, second, first)
    if reverse is not None:
        return ("second_then_first", reverse)
    return None


def direct_clearance(speeds: tuple[int, ...], phase: F) -> F:
    return min(distance(speed, phase) for speed in speeds)


def main() -> None:
    # Independent direct-phase versus affine-lattice equivalence.
    equivalence_profiles = 0
    for r in range(1, 12, 2):
        odds = tuple(range(1, 20, 2))
        for index, first in enumerate(odds):
            for second in odds[index + 1 :]:
                direct = direct_covering_phase(2 * r, first, second)
                lattice = lattice_certificate(2 * r, first, second)
                defect = defect_certificate(2 * r, first, second)
                require((direct is not None) == (lattice is not None), f"equivalence {(r, first, second)}")
                require((direct is not None) == (defect is not None), f"defect equivalence {(r, first, second)}")
                equivalence_profiles += 1

    # The three cheap necessary gates are exact consequences of odd,
    # gcd-divisible nonzero defects.  They are deliberately not treated as
    # sufficient: (2,7,11) is the first small compatibility hostile.
    pairwise_gate_checks = 0
    first_false_converse = None
    for r in range(1, 12, 2):
        odds = tuple(range(1, 20, 2))
        for index, first in enumerate(odds):
            for second in odds[index + 1 :]:
                gate = (
                    first + second > 7 * gcd(first, second)
                    and first + 2 * r > 7 * gcd(first, r)
                    and second + 2 * r > 7 * gcd(second, r)
                )
                lattice = lattice_certificate(2 * r, first, second)
                if lattice is not None:
                    require(gate, f"lattice implies gcd gates {(r, first, second)}")
                elif gate and first_false_converse is None:
                    first_false_converse = (2 * r, first, second)
                pairwise_gate_checks += 1
    require(first_false_converse == (2, 7, 11), "first pairwise false converse")

    # Minimal positive signal beyond the old d=2 odd-pair gate.
    # (2,3,5) has 3+5>7 but cannot cover; (2,5,7) isolates the strict-open
    # endpoint 5+2=7 while its other two gcd gates are strict.
    require(3 + 5 > 7 * gcd(3, 5), "old odd-pair gate passes")
    require(lattice_certificate(2, 3, 5) is None, "small new positive control")
    require(direct_covering_phase(2, 3, 5) is None, "small positive direct control")
    require(5 + 7 > 7 * gcd(5, 7), "endpoint control odd pair strict")
    require(5 + 2 == 7 * gcd(5, 1), "quarter-shift endpoint equality")
    require(7 + 2 > 7 * gcd(7, 1), "other quarter gate strict")
    require(lattice_certificate(2, 5, 7) is None, "open endpoint is safe")

    # Global minimum hostile by maximum exception speed.  The exact scan is
    # backed by the proof note: below 9 the only pairwise survivor (6,5,7)
    # would require three +/-1 defects, incompatible with 6N12-7N1e+5N2e=0.
    minimal = None
    for maximum in range(1, 10):
        for even in range(2, maximum + 1, 4):
            odds = tuple(x for x in range(1, maximum + 1, 2))
            for index, first in enumerate(odds):
                for second in odds[index + 1 :]:
                    phase = direct_covering_phase(even, first, second)
                    if phase is not None:
                        minimal = (maximum, even, first, second, phase)
                        break
                if minimal is not None:
                    break
            if minimal is not None:
                break
        if minimal is not None:
            break
    require(minimal == (9, 2, 7, 9, F(6, 49)), "minimal geometric hostile")
    require(bad_label_sets(F(6, 49), 2, 7, 9) == ((0, 2), (1,), (3,)), "minimal label partition")
    minimal_certificate = oriented_lattice_certificate(1, 7, 9)
    require(minimal_certificate[:6] == (2, 7, 1, 1, 1, 1), "minimal unit-defect certificate")
    require(oriented_defect_certificate(1, 7, 9)[:3] == (-1, -1, -1), "minimal defect-only certificate")
    require(2 * 1 - 9 * 1 + 7 * 1 == 0, "minimal additive relation 2+7=9")

    # Canonical H={1,...,10} selector hostile.  Exact search confirms this is
    # the smallest-max exception triple covering one of the ten H-optimal
    # phases y=k/11.  The full row is nevertheless 1/11-safe elsewhere.
    h10_minimal = None
    for maximum in range(1, 12):
        for even in range(2, maximum + 1, 4):
            odds = tuple(x for x in range(1, maximum + 1, 2))
            for index, first in enumerate(odds):
                for second in odds[index + 1 :]:
                    for k in range(1, 11):
                        y = F(k, 11)
                        if labels_are_covered(y, even, first, second):
                            h10_minimal = (maximum, even, first, second, y)
                            break
                    if h10_minimal is not None:
                        break
                if h10_minimal is not None:
                    break
            if h10_minimal is not None:
                break
        if h10_minimal is not None:
            break
    require(h10_minimal == (11, 2, 9, 11, F(1, 11)), "minimal H10 selector hostile")
    require(bad_label_sets(F(1, 11), 2, 9, 11) == ((0, 2), (3,), (1,)), "H10 label partition")

    divisible_body = (8, 12, 20, 24, 28, 32, 36, 40)
    pair = (4, 16)
    hostile_body = divisible_body + (2, 9, 11)
    hostile_row = hostile_body + pair
    require(gcd(*hostile_body) == 1, "typed hostile body primitive")
    require(len(set(hostile_row)) == 13, "typed hostile row distinct")
    require(direct_clearance(hostile_row, F(21, 22)) == F(1, 11), "hostile row itself lonely")

    positive_body = divisible_body + (2, 5, 7)
    positive_row = positive_body + pair
    require(gcd(*positive_body) == 1, "typed endpoint body primitive")
    require(len(set(positive_row)) == 13, "typed endpoint row distinct")
    require(direct_clearance(positive_row, F(3, 11)) == F(1, 11), "typed endpoint explicit lift")

    print("LRC14_D4_AFFINE_DEFECT_LATTICE_BOUNDARY_THM4030_EXACT_AUDIT")
    print("scope=THM4004_c4_equals_8_selector_boundary;LRC14=OPEN")
    print("iff=direct_phase_cover=A_B_C_affine_certificate=bounded_odd_defect_relation")
    print("relation=(2r)N_ab-bN_ar+aN_br=0;all_N_odd_nonzero")
    print("gcd_gates=a+b>7gcd(a,b);a+2r>7gcd(a,r);b+2r>7gcd(b,r)")
    print(f"equivalence_profiles={equivalence_profiles};pairwise_gate_checks={pairwise_gate_checks}")
    print(f"pairwise_not_sufficient={first_false_converse}")
    print(f"minimal_geometric_hostile={minimal};certificate={minimal_certificate[:6]}")
    print(f"minimal_H10_hostile={h10_minimal};safe_phase=21/22;clearance=1/11")
    print("strict_endpoint_positive=(2,5,7);new_small_positive=(2,3,5)")
    print(f"checks={CHECKS}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
