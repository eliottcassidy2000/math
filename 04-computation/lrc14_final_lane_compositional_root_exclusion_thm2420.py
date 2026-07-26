#!/usr/bin/env python3
"""Exact companion for THM-2420.

Checks the two composed 7/13-root disintegrations which exclude depths
M=2 and M=1 in the primitive final septimal lane.  It also checks the
uniform narrow-comb/phase half-overlap lemma and a hostile showing that
the inherited 13-unit hypothesis on H is load-bearing.

Only standard-library integer and Fraction arithmetic is used.
"""

from fractions import Fraction
from math import ceil, gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def frac(x):
    return x % 1


def circle_norm(x):
    r = frac(x)
    return min(r, 1 - r)


def danger(v, x):
    return circle_norm(v * x) < Fraction(1, 14)


def phase_chamber(h, x):
    r = frac(h * x)
    return Fraction(6, 13) < r < Fraction(7, 13)


def phase_far(h, x):
    return circle_norm(h * x) > Fraction(3, 13)


def norm_boundaries(v, radius):
    """Boundary points of ||v*x|| < radius on the circle."""

    v = abs(v)
    out = set()
    for j in range(v):
        out.add(frac((Fraction(j) + radius) / v))
        out.add(frac((Fraction(j) - radius) / v))
    return out


def phase_boundaries(v, left, right):
    """Boundary points of left < {v*x} < right."""

    v = abs(v)
    out = set()
    for j in range(v):
        out.add((Fraction(j) + left) / v)
        out.add((Fraction(j) + right) / v)
    return out


def exact_measure(predicate, boundary_sets):
    points = {Fraction(0), Fraction(1)}
    for boundary_set in boundary_sets:
        points.update(boundary_set)
    ordered = sorted(points)
    total = Fraction(0)
    for left, right in zip(ordered, ordered[1:]):
        if left == right:
            continue
        midpoint = (left + right) / 2
        if predicate(midpoint):
            total += right - left
    return total


def m1_measures(h, u, v, g):
    """Exact M=1 phase measures; g=L-1>=1."""

    high = 13 * 7**g * v
    j_measure = exact_measure(
        lambda x: danger(13 * u, x)
        and not danger(high, x)
        and phase_chamber(h, x),
        (
            norm_boundaries(13 * u, Fraction(1, 14)),
            norm_boundaries(high, Fraction(1, 14)),
            phase_boundaries(h, Fraction(6, 13), Fraction(7, 13)),
        ),
    )
    narrow_overlap = exact_measure(
        lambda x: danger(u, x)
        and danger(13 * u, x)
        and phase_chamber(h, x),
        (
            norm_boundaries(u, Fraction(1, 14)),
            norm_boundaries(13 * u, Fraction(1, 14)),
            phase_boundaries(h, Fraction(6, 13), Fraction(7, 13)),
        ),
    )
    transition_phase = exact_measure(
        lambda x: danger(13 * u, x)
        and not danger(u, x)
        and not danger(high, x)
        and phase_chamber(h, x),
        (
            norm_boundaries(13 * u, Fraction(1, 14)),
            norm_boundaries(u, Fraction(1, 14)),
            norm_boundaries(high, Fraction(1, 14)),
            phase_boundaries(h, Fraction(6, 13), Fraction(7, 13)),
        ),
    )
    return j_measure, narrow_overlap, transition_phase


def m2_measures(h, u, v, g):
    """Exact M=2 composed-root measures; g=L-2>=1."""

    high_terminal = 13 * 7**g * v
    transformed = exact_measure(
        lambda t: danger(13 * u, t)
        and not danger(u, t)
        and not danger(high_terminal, t)
        and phase_far(h, t),
        (
            norm_boundaries(13 * u, Fraction(1, 14)),
            norm_boundaries(u, Fraction(1, 14)),
            norm_boundaries(high_terminal, Fraction(1, 14)),
            norm_boundaries(h, Fraction(3, 13)),
        ),
    )

    q_divided = 7 * u
    high_original = 13 * 7 ** (g + 1) * v
    landed = exact_measure(
        lambda y: danger(13 * q_divided, y)
        and not danger(q_divided, y)
        and not danger(high_original, y)
        and phase_chamber(h, y),
        (
            norm_boundaries(13 * q_divided, Fraction(1, 14)),
            norm_boundaries(q_divided, Fraction(1, 14)),
            norm_boundaries(high_original, Fraction(1, 14)),
            phase_boundaries(h, Fraction(6, 13), Fraction(7, 13)),
        ),
    )
    return transformed, landed


def useful_m2_roots(s, h, u, v, g):
    high = 13 * 7**g * v
    roots = [(s + j) / 13 for j in range(13)]
    return sum(
        danger(13 * u, t)
        and not danger(u, t)
        and not danger(high, t)
        and phase_far(h, t)
        for t in roots
    )


def main():
    e_mass = Fraction(6, 49)
    m2_transformed_floor = Fraction(36, 637)
    m2_landed_floor = Fraction(36, 4459)
    m1_j_mass = Fraction(6, 637)
    overlap_ceiling = Fraction(1, 182)
    m1_landed_floor = Fraction(5, 1274)

    # The elementary numerical inequalities in the half-overlap lemma.
    half_ratio_checks = 0
    for b in range(2, 10001):
        if b % 13 == 0:
            continue
        require(
            Fraction(ceil(Fraction(b, 13)), b) <= Fraction(1, 2),
            f"grid half-ratio failed at b={b}",
        )
        half_ratio_checks += 1
    for a in range(17, 10001):
        length = Fraction(a, 91)
        require(
            (length + 1) / 13 <= length / 2,
            f"one-comb periodic bound failed at a={a}",
        )

    # Independent exact interval replay of the narrow overlap.
    overlap_checks = 0
    max_overlap = Fraction(0)
    for a in range(1, 81):
        if a % 13 == 0:
            continue
        for b in range(1, 81):
            if b % 13 == 0 or gcd(a, b) != 1:
                continue
            overlap = exact_measure(
                lambda x: circle_norm(b * x) < Fraction(1, 182)
                and phase_chamber(a, x),
                (
                    norm_boundaries(b, Fraction(1, 182)),
                    phase_boundaries(
                        a, Fraction(6, 13), Fraction(7, 13)
                    ),
                ),
            )
            require(
                overlap <= overlap_ceiling,
                f"narrow half-overlap failed at (a,b)=({a},{b})",
            )
            max_overlap = max(max_overlap, overlap)
            overlap_checks += 1
    require(max_overlap == overlap_ceiling, "half-overlap ceiling not sharp")

    m1_cases = (
        (1, 1, 1, 1),
        (5, 2, 1, 1),
        (8, 3, 2, 1),
        (11, 5, 3, 2),
        (2, 8, 5, 1),
        (17, 11, 2, 1),
    )
    m1_min = None
    for h, u, v, g in m1_cases:
        require(h % 7 and h % 13 and u % 7 and u % 13 and v % 7, "bad M1 case")
        j_measure, overlap, transition_phase = m1_measures(h, u, v, g)
        require(
            j_measure == m1_j_mass,
            f"M1 root-section mass failed at {(h,u,v,g)}: {j_measure}",
        )
        require(
            overlap <= overlap_ceiling,
            f"M1 overlap failed at {(h,u,v,g)}: {overlap}",
        )
        require(
            transition_phase >= m1_landed_floor,
            f"M1 landed floor failed at {(h,u,v,g)}: {transition_phase}",
        )
        m1_min = (
            transition_phase
            if m1_min is None
            else min(m1_min, transition_phase)
        )

    m2_cases = (
        (1, 1, 1, 1),
        (5, 2, 3, 1),
        (8, 3, 2, 1),
        (11, 5, 4, 1),
        (2, 8, 5, 2),
    )
    m2_min = None
    for h, u, v, g in m2_cases:
        require(h % 7 and h % 13 and u % 7 and u % 13 and v % 7, "bad M2 case")
        transformed, landed = m2_measures(h, u, v, g)
        require(
            transformed >= m2_transformed_floor,
            f"M2 transformed floor failed at {(h,u,v,g)}: {transformed}",
        )
        require(
            landed == transformed / 7,
            f"M2 7-root landing failed at {(h,u,v,g)}",
        )
        require(
            landed >= m2_landed_floor,
            f"M2 landed floor failed at {(h,u,v,g)}: {landed}",
        )
        m2_min = landed if m2_min is None else min(m2_min, landed)

    # Direct 13-root controls on a denominator coprime to every speed.
    root_bases_checked = 0
    for h, u, v, g in m2_cases:
        for numerator in range(1, 1009):
            s = Fraction(numerator, 1009)
            if not danger(u, s) or danger(7**g * v, s):
                continue
            roots = [(s + j) / 13 for j in range(13)]
            require(all(danger(13 * u, t) for t in roots), "13u constancy failed")
            require(
                all(not danger(13 * 7**g * v, t) for t in roots),
                "high 13-root constancy failed",
            )
            require(sum(danger(u, t) for t in roots) == 1, "one-root law failed")
            require(sum(phase_far(h, t) for t in roots) == 7, "phase 7-root law failed")
            require(
                useful_m2_roots(s, h, u, v, g) >= 6,
                "six-useful-root floor failed",
            )
            root_bases_checked += 1

    # The 13-unit condition on H is load-bearing.
    hostile_s = Fraction(11, 1009)
    require(danger(1, hostile_s) and not danger(7, hostile_s), "bad hostile base")
    hostile_useful = useful_m2_roots(hostile_s, 13, 1, 1, 1)
    require(hostile_useful == 0, "13-divisible-H hostile did not collapse")

    print("theorem=THM-2420")
    print("arithmetic=Fraction-only")
    print(f"different-septimal-depth-safe-mass={e_mass}")
    print("M2-useful-root-floor=6-of-13")
    print(f"M2-transformed-band-floor={m2_transformed_floor}")
    print(f"M2-landed-phase-floor={m2_landed_floor}")
    print(f"M2-exact-cases={len(m2_cases)},minimum={m2_min}")
    print(f"M1-root-section-mass={m1_j_mass}")
    print(f"narrow-phase-overlap-ceiling={overlap_ceiling}")
    print(f"M1-landed-phase-floor={m1_landed_floor}")
    print(f"M1-exact-cases={len(m1_cases)},minimum={m1_min}")
    print(f"half-ratio-checks={half_ratio_checks}")
    print(f"independent-overlap-checks={overlap_checks},maximum={max_overlap}")
    print(f"direct-13-root-bases-checked={root_bases_checked}")
    print(f"hostile-13-divisible-H-useful-roots={hostile_useful}")
    print("conclusion=primitive-final-lane-empty")
    print("scope=no scalar-row decrement and no LRC(14) conclusion")
    print("status=PASS")


if __name__ == "__main__":
    main()
