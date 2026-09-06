#!/usr/bin/env python3
"""Independent integer-grid audit of the q=2 sharp physical-union cap.

This implementation does not import the Fraction interval engine.  Every
periodic tooth endpoint is lifted to one exact common integer denominator,
and all intersections are accumulated as integer lengths before one final
Fraction conversion.
"""

from __future__ import annotations

from fractions import Fraction as Q
from math import gcd, lcm


IntInterval = tuple[int, int]
TERNARY_TARGET = Q(72, 539)
FULL_TARGET = Q(214, 1449)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def centered7(n: int) -> int:
    r = n % 7
    return r if r <= 3 else r - 7


def pair_mass(p: int, q: int) -> Q:
    defect = centered7((q - p) // 2) ** 2 - centered7((q + p) // 2) ** 2
    return Q(2, 49) * (1 + Q(defect, p * q))


def pair_defect(p: int, q: int) -> Q:
    defect = centered7((q - p) // 2) ** 2 - centered7((q + p) // 2) ** 2
    return Q(defect, p * q)


def exceptional_rays(
    neutral: Q, product_limit: int, ternary_units: bool
) -> list[tuple[int, int]]:
    answer: list[tuple[int, int]] = []
    for p in range(1, product_limit, 2):
        if p * p >= product_limit:
            break
        if ternary_units and p % 3 == 0:
            continue
        for q in range(p + 2, ((product_limit - 1) // p) + 1, 2):
            if (ternary_units and q % 3 == 0) or gcd(p, q) != 1:
                continue
            if pair_defect(p, q) > neutral:
                answer.append((p, q))
    return answer


def merge(arcs: list[IntInterval]) -> list[IntInterval]:
    answer: list[list[int]] = []
    for left, right in sorted(arcs):
        if left >= right:
            continue
        if not answer or left > answer[-1][1]:
            answer.append([left, right])
        else:
            answer[-1][1] = max(answer[-1][1], right)
    return [(left, right) for left, right in answer]


def periodic_danger(speed: int, denominator: int) -> list[IntInterval]:
    require(denominator % (14 * speed) == 0, (speed, denominator))
    step = denominator // speed
    radius = denominator // (14 * speed)
    arcs: list[IntInterval] = []
    for k in range(speed):
        left, right = k * step - radius, k * step + radius
        if left < 0:
            arcs.extend(((0, right), (denominator + left, denominator)))
        elif right > denominator:
            arcs.extend(((left, denominator), (0, right - denominator)))
        else:
            arcs.append((left, right))
    return merge(arcs)


def half_shift(arcs: list[IntInterval], denominator: int) -> list[IntInterval]:
    shift = denominator // 2
    answer: list[IntInterval] = []
    for left, right in arcs:
        left += shift
        right += shift
        if left >= denominator:
            answer.append((left - denominator, right - denominator))
        elif right > denominator:
            answer.extend(((left, denominator), (0, right - denominator)))
        else:
            answer.append((left, right))
    return merge(answer)


def subtract(first: list[IntInterval], second: list[IntInterval]) -> list[IntInterval]:
    answer: list[IntInterval] = []
    for left, right in first:
        pieces = [(left, right)]
        for cut_left, cut_right in second:
            following: list[IntInterval] = []
            for piece_left, piece_right in pieces:
                if cut_right <= piece_left or cut_left >= piece_right:
                    following.append((piece_left, piece_right))
                else:
                    if piece_left < cut_left:
                        following.append((piece_left, cut_left))
                    if cut_right < piece_right:
                        following.append((cut_right, piece_right))
            pieces = following
        answer.extend(pieces)
    return merge(answer)


def intersection_length(first: list[IntInterval], second: list[IntInterval]) -> int:
    i = j = 0
    total = 0
    while i < len(first) and j < len(second):
        total += max(0, min(first[i][1], second[j][1]) - max(first[i][0], second[j][0]))
        if first[i][1] < second[j][1]:
            i += 1
        elif second[j][1] < first[i][1]:
            j += 1
        else:
            i += 1
            j += 1
    return total


def intersection(first: list[IntInterval], second: list[IntInterval]) -> list[IntInterval]:
    i = j = 0
    answer: list[IntInterval] = []
    while i < len(first) and j < len(second):
        left = max(first[i][0], second[j][0])
        right = min(first[i][1], second[j][1])
        if left < right:
            answer.append((left, right))
        if first[i][1] < second[j][1]:
            i += 1
        elif second[j][1] < first[i][1]:
            j += 1
        else:
            i += 1
            j += 1
    return merge(answer)


def exposure_grid(p: int, q: int) -> tuple[int, list[IntInterval], Q, Q]:
    denominator = 14 * p * q
    a = merge(periodic_danger(p, denominator) + periodic_danger(q, denominator))
    b = half_shift(a, denominator)
    r = subtract(b, a)
    sigma = Q(intersection_length(a, b), denominator)
    mu_r = Q(sum(right - left for left, right in r), denominator)
    require(sigma == pair_mass(p, q), (p, q, sigma, pair_mass(p, q)))
    return denominator, r, mu_r, sigma


def lifted_exposure(
    base: list[IntInterval], base_denominator: int, t: int, denominator: int
) -> list[IntInterval]:
    source_denominator = base_denominator * t
    require(denominator % source_denominator == 0, (source_denominator, denominator))
    scale = denominator // source_denominator
    # k-major order is sorted and the components are disjoint a.e.
    return [
        ((left + k * base_denominator) * scale, (right + k * base_denominator) * scale)
        for k in range(t)
        for left, right in base
    ]


def physical_mass(
    p: int,
    q: int,
    t: int,
    c: int,
    base_denominator: int,
    r: list[IntInterval],
    sigma: Q,
) -> Q:
    denominator = lcm(base_denominator * t, 14 * c)
    rt = lifted_exposure(r, base_denominator, t, denominator)
    dc = periodic_danger(c, denominator)
    return sigma + 2 * Q(intersection_length(rt, dc), denominator)


def audit_ray(
    p: int, q: int, target: Q, ternary_units: bool
) -> dict[str, object]:
    base_denominator, r, mu_r, sigma = exposure_grid(p, q)
    eta = (target - sigma) / 2 - mu_r / 7
    require(eta > 0, ((p, q), "nonpositive asymptotic margin", eta))
    cutoff = Q(len(r), 3) / eta
    require(cutoff.denominator > 1, ((p, q), "integral cutoff", cutoff))
    checked = 0
    maximum = Q(-1)
    maximizers: set[tuple[int, int, int]] = set()
    for t in range(1, cutoff.numerator // cutoff.denominator + 1, 2):
        if ternary_units and t % 3 == 0:
            continue
        upper_c = (cutoff.numerator - 1) // (cutoff.denominator * t)
        for c in range(1, upper_c + 1, 2):
            if (
                (ternary_units and c % 3 == 0)
                or gcd(t, c) != 1
                or c in (t * p, t * q)
            ):
                continue
            mass = physical_mass(p, q, t, c, base_denominator, r, sigma)
            require(mass <= target, ((p, q), t, c, mass))
            checked += 1
            shape = tuple(sorted((t * p, t * q, c)))
            if mass > maximum:
                maximum, maximizers = mass, {shape}
            elif mass == maximum:
                maximizers.add(shape)
    return {
        "ray": (p, q),
        "base_denominator": base_denominator,
        "components": len(r),
        "mu_r": mu_r,
        "sigma": sigma,
        "eta": eta,
        "cutoff": cutoff,
        "checked": checked,
        "maximum": maximum,
        "maximizers": sorted(maximizers),
    }


def main() -> None:
    print("LRC14_Q2_THREE_ODD_TAIL_TWO_DOMAIN_SHARP_CAP_INDEPENDENT")
    print("status=INDEPENDENT_INTEGER_GRID_ALL_HEIGHT_FINITE_REDUCTION_AUDIT")
    general_pair_box = [
        (pair_mass(p, q), p, q)
        for p in range(1, 17, 2)
        for q in range(p + 2, 18, 2)
        if gcd(p, q) == 1 and p * q < 17
    ]
    ternary_pair_box = [
        (pair_mass(p, q), p, q)
        for p in range(1, 33, 2)
        for q in range(p + 2, 34, 2)
        if p % 3 and q % 3 and gcd(p, q) == 1 and p * q < 33
    ]
    require(max(general_pair_box) == (Q(4, 63), 1, 9), general_pair_box)
    require(max(ternary_pair_box) == (Q(4, 77), 1, 11), ternary_pair_box)
    ternary_product_33 = [
        (p, q)
        for p in range(1, 34, 2)
        for q in range(p + 2, 34, 2)
        if p % 3 and q % 3 and gcd(p, q) == 1 and p * q == 33
    ]
    require(ternary_product_33 == [], ternary_product_33)
    ternary_neutral_edges = [
        (p, q)
        for p in range(1, 100, 2)
        for q in range(p + 2, 100, 2)
        if p % 3 and q % 3 and gcd(p, q) == 1 and p * q <= 99
        and pair_defect(p, q) == Q(1, 11)
    ]
    require(ternary_neutral_edges == [], ternary_neutral_edges)
    print("PAIR_CAP full_odd=4/63 equality=(1,9) odd_3_unit=4/77 equality=(1,11) ternary_pq33_boundary=excluded ternary_neutral_equal_edges=none")

    domains = (
        ("full_odd", FULL_TARGET, Q(128, 621), 44, False,
         [(1, 9), (1, 11), (1, 23), (3, 11)], (1, 9, 23)),
        ("odd_3_unit", TERNARY_TARGET, Q(1, 11), 99, True,
         [(1, 11), (1, 23), (1, 25), (1, 37), (5, 11)], (1, 7, 11)),
    )
    for name, target, neutral, product_limit, unit_filter, expected_rays, expected_equality in domains:
        rays = exceptional_rays(neutral, product_limit, unit_filter)
        require(rays == expected_rays, (name, rays))
        print(f"DOMAIN name={name} target={target} neutral={neutral} exceptional_rays={rays}")
        rows = [audit_ray(*ray, target, unit_filter) for ray in rays]
        for row in rows:
            print(name.upper() + "_FAMILY " + " ".join(f"{key}={value}" for key, value in row.items()))
        require(max(row["maximum"] for row in rows) == target, (name, rows))
        equality = sorted(
            {
                shape
                for row in rows
                if row["maximum"] == target
                for shape in row["maximizers"]
            }
        )
        require(equality == [expected_equality], (name, equality))
        print(f"{name.upper()}_primitive_equality_shapes={equality}")
        print(f"{name.upper()}_finite_presentations={sum(row['checked'] for row in rows)}")

    for name, tails, physical_count, physical_width, quotient_count, quotient_width in (
        ("full_odd", (1, 9, 23), 16, Q(10, 483), 8, Q(20, 483)),
        ("odd_3_unit", (1, 7, 11), 12, Q(1, 77), 6, Q(2, 77)),
    ):
        denominator = lcm(*(14 * speed for speed in tails))
        danger_union = merge(
            [arc for speed in tails for arc in periodic_danger(speed, denominator)]
        )
        physical = intersection(danger_union, half_shift(danger_union, denominator))
        require(len(physical) == physical_count, (name, physical))
        require(Q(max(right - left for left, right in physical), denominator) == physical_width, (name, physical))
        quotient_arcs: list[IntInterval] = []
        for left, right in physical:
            left, right = 2 * left, 2 * right
            if left >= denominator:
                quotient_arcs.append((left - denominator, right - denominator))
            elif right > denominator:
                quotient_arcs.extend(((left, denominator), (0, right - denominator)))
            else:
                quotient_arcs.append((left, right))
        quotient = merge(quotient_arcs)
        require(len(quotient) == quotient_count, (name, quotient))
        require(Q(max(right - left for left, right in quotient), denominator) == quotient_width, (name, quotient))
        print(
            f"{name.upper()}_ENDPOINT physical_components={physical_count} "
            f"physical_max_width={physical_width} quotient_components={quotient_count} "
            f"quotient_max_width={quotient_width}"
        )
    print(f"HAAR_GATES full_odd={FULL_TARGET} odd_3_unit={TERNARY_TARGET}")
    print("PASS")


if __name__ == "__main__":
    main()
