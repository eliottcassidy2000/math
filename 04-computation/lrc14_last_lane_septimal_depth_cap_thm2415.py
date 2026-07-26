#!/usr/bin/env python3
"""Exact companion for THM-2415.

The proof itself is symbolic.  This companion checks its finite root-word
identities, exact Haar-measure constants, a bank of hostile transition masks,
the orbit-hitting cutoff, and the sharp M=2 stopping control using only
standard-library rational arithmetic.
"""

from fractions import Fraction
from math import floor


def frac(x: Fraction) -> Fraction:
    return x - floor(x)


def circle_norm(x: Fraction) -> Fraction:
    r = frac(x)
    return min(r, 1 - r)


def danger(v: int, x: Fraction, width: int = 1) -> bool:
    """Open D_(v,width)={x: ||v*x|| < width/14}."""

    return circle_norm(v * x) < Fraction(width, 14)


def partition_measure(coefficients, predicate) -> Fraction:
    """Exact measure of a Boolean combination of open rational combs.

    Endpoints have measure zero.  On every complementary cell the predicate
    is constant, so testing its midpoint gives the exact Haar measure.
    """

    cuts = {Fraction(0), Fraction(1)}
    for v in coefficients:
        radius = Fraction(1, 14 * v)
        for k in range(v):
            center = Fraction(k, v)
            cuts.add(frac(center - radius))
            cuts.add(frac(center + radius))
    points = sorted(cuts)
    total = Fraction(0)
    for left, right in zip(points, points[1:]):
        mid = (left + right) / 2
        if predicate(mid):
            total += right - left
    return total


def transition_measure(u: int, w: int, d: int) -> Fraction:
    high = (7**d) * w
    coefficients = (u, 13 * u, high)
    return partition_measure(
        coefficients,
        lambda x: danger(13 * u, x)
        and not danger(u, x)
        and not danger(high, x),
    )


def root_word(v: int, y: Fraction):
    return {
        j
        for j in range(7)
        if danger(v, (y + j) / 7, width=2)
    }


def guard_word_in_h_gauge(r: Fraction):
    """Root word of E_H after relabelling j by H*j modulo seven."""

    return {
        k
        for k in range(7)
        if circle_norm((r + k) / 7) < Fraction(1, 7)
    }


def drifted_guard_word_in_h_gauge(r: Fraction):
    return {
        k
        for k in range(7)
        if circle_norm(13 * (r + k) / 7) < Fraction(1, 7)
    }


def grid_hits_interval(n: int, shift_index: int) -> int:
    """Count an n-grid shifted by shift_index/(13*n) in (6/13,7/13)."""

    shift = Fraction(shift_index, 13 * n)
    return sum(
        Fraction(6, 13) < frac(shift + Fraction(j, n)) < Fraction(7, 13)
        for j in range(n)
    )


def main() -> None:
    # The two unit combs touch only along their central 1/91 cell.
    inter = partition_measure(
        (1, 13), lambda x: danger(1, x) and danger(13, x)
    )
    diff = partition_measure(
        (1, 13), lambda x: danger(13, x) and not danger(1, x)
    )
    assert inter == Fraction(1, 91)
    assert diff == Fraction(12, 91)

    one_root_base_mass = 7 * diff
    high_safe_base_mass = Fraction(6, 7)
    base_intersection_floor = (
        one_root_base_mass + high_safe_base_mass - 1
    )
    transition_floor = base_intersection_floor / 7
    assert one_root_base_mass == Fraction(12, 13)
    assert base_intersection_floor == Fraction(71, 91)
    assert transition_floor == Fraction(71, 637)

    # Hostile bank: unrelated u,w can make the lower bound non-sharp, but
    # never violate it.
    tested = []
    for u in range(1, 14):
        if u % 7 == 0:
            continue
        for w in range(1, 14):
            if w % 7 == 0:
                continue
            for d in (1, 2):
                mass = transition_measure(u, w, d)
                assert mass >= transition_floor
                tested.append((mass, u, w, d))
    min_mass, min_u, min_w, min_d = min(tested)

    # Exact affine guard transport.  Every open thirteenth digit chamber has
    # W_13H={b,b+1}; chamber b=6 is exactly the collision with W_H.
    base_guard = {0, 6}
    for b in range(13):
        r = Fraction(2 * b + 1, 26)
        assert guard_word_in_h_gauge(r) == base_guard
        expected = {b % 7, (b + 1) % 7}
        assert drifted_guard_word_in_h_gauge(r) == expected
    r6 = Fraction(13, 26)
    assert Fraction(6, 13) < r6 < Fraction(7, 13)
    assert drifted_guard_word_in_h_gauge(r6) == base_guard

    # Every shifted 7^a grid hits the phase chamber once a>=2.  Enumerating
    # shifts at all combinatorial endpoint chambers is an exact finite check.
    phase_minima = {}
    for n in (7, 49, 343):
        counts = [grid_hits_interval(n, s) for s in range(13 * n)]
        phase_minima[n] = min(counts)
        assert min(counts) >= n // 13
    assert phase_minima[7] == 0
    assert phase_minima[49] >= 1

    # Sharp M=2 stopping control for the phase-forcing mechanism.
    q = 7
    d = 49
    h = 1
    y = Fraction(1, 91)
    assert danger(13 * q, y)
    assert not danger(q, y)
    assert not danger(d, y)
    orbit_phases = [frac(h * (y + Fraction(j, 7))) for j in range(7)]
    assert not any(
        Fraction(6, 13) < r < Fraction(7, 13) for r in orbit_phases
    )

    print("theorem=THM-2415")
    print("arithmetic=Fraction-only")
    print(f"mu(D_13u_inter_D_u)={inter}")
    print(f"mu(D_13u_minus_D_u)={diff}")
    print(f"one-root-base-mass={one_root_base_mass}")
    print(f"high-safe-base-mass={high_safe_base_mass}")
    print(f"base-intersection-floor={base_intersection_floor}")
    print(f"transition-floor={transition_floor}")
    print(f"hostile-bank-size={len(tested)}")
    print(
        "hostile-bank-min="
        f"{min_mass} at (u,w,d)=({min_u},{min_w},{min_d})"
    )
    print(
        "guard-transport="
        "W_H={0,6}; W_13H(b)={b,b+1}; b=6 gives equality"
    )
    print(
        "phase-grid-minima="
        + ",".join(f"N{n}:{phase_minima[n]}" for n in sorted(phase_minima))
    )
    print(
        "M2-hostile="
        "Q=7,D=49,H=1,y=1/91 lies in transition set but misses phase chamber"
    )
    print("status=PASS")


if __name__ == "__main__":
    main()
