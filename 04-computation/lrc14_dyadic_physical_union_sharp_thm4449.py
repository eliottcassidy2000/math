#!/usr/bin/env python3
"""All-height finite reduction for the q=2 three-odd-tail physical union.

The analytic input is documented in REPORT.md.  This exact certificate:

* derives both domains' exceptional pair rays from seventh-rounding;
* constructs each one-sided exposure set R=(A+1/2) minus A;
* derives the BV/Fourier product cutoff exactly; and
* checks every coprime primitive (scale, third-tail) presentation below it.

All arithmetic used for theorem-facing comparisons is fractions.Fraction.
"""

from __future__ import annotations

from fractions import Fraction as Q
from itertools import combinations
from math import gcd


Interval = tuple[Q, Q]
TERNARY_TARGET = Q(72, 539)
FULL_TARGET = Q(214, 1449)
BASELINE = Q(2, 49)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def centered7(n: int) -> int:
    residue = n % 7
    return residue if residue <= 3 else residue - 7


def pair_defect(p: int, q: int) -> Q:
    """Seventh-rounding correction delta with coprime odd p<q."""
    require(0 < p < q and p % 2 and q % 2 and gcd(p, q) == 1, (p, q))
    d = centered7((p + q) // 2)
    e = centered7((q - p) // 2)
    return Q(e * e - d * d, p * q)


def pair_mass(p: int, q: int) -> Q:
    return BASELINE * (1 + pair_defect(p, q))


def reduce_pair(a: int, b: int) -> tuple[int, int]:
    g = gcd(a, b)
    p, q = a // g, b // g
    return (p, q) if p < q else (q, p)


def merge(intervals: list[Interval]) -> list[Interval]:
    answer: list[list[Q]] = []
    for left, right in sorted(intervals):
        if left >= right:
            continue
        if not answer or left > answer[-1][1]:
            answer.append([left, right])
        else:
            answer[-1][1] = max(answer[-1][1], right)
    return [(left, right) for left, right in answer]


def danger(speed: int) -> list[Interval]:
    arcs: list[Interval] = []
    radius = Q(1, 14 * speed)
    for k in range(speed):
        left, right = Q(k, speed) - radius, Q(k, speed) + radius
        if left < 0:
            arcs.extend(((Q(0), right), (1 + left, Q(1))))
        elif right > 1:
            arcs.extend(((left, Q(1)), (Q(0), right - 1)))
        else:
            arcs.append((left, right))
    return merge(arcs)


def shift_half(intervals: list[Interval]) -> list[Interval]:
    arcs: list[Interval] = []
    for left, right in intervals:
        left, right = left + Q(1, 2), right + Q(1, 2)
        if left >= 1:
            arcs.append((left - 1, right - 1))
        elif right > 1:
            arcs.extend(((left, Q(1)), (Q(0), right - 1)))
        else:
            arcs.append((left, right))
    return merge(arcs)


def intersect(first: list[Interval], second: list[Interval]) -> list[Interval]:
    i = j = 0
    answer: list[Interval] = []
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


def subtract(first: list[Interval], second: list[Interval]) -> list[Interval]:
    answer: list[Interval] = []
    for left, right in first:
        pieces = [(left, right)]
        for cut_left, cut_right in second:
            next_pieces: list[Interval] = []
            for piece_left, piece_right in pieces:
                if cut_right <= piece_left or cut_left >= piece_right:
                    next_pieces.append((piece_left, piece_right))
                else:
                    if piece_left < cut_left:
                        next_pieces.append((piece_left, cut_left))
                    if cut_right < piece_right:
                        next_pieces.append((cut_right, piece_right))
            pieces = next_pieces
        answer.extend(pieces)
    return merge(answer)


def measure(intervals: list[Interval]) -> Q:
    return sum((right - left for left, right in intervals), Q(0))


def literal_physical(tails: tuple[int, int, int]) -> Q:
    union = merge([arc for speed in tails for arc in danger(speed)])
    return measure(intersect(union, shift_half(union)))


def preimage(intervals: list[Interval], multiplier: int) -> list[Interval]:
    # k-major order is already circularly sorted; avoiding a Fraction sort is
    # important in the densest (1,37) finite box.
    return [
        ((left + k) / multiplier, (right + k) / multiplier)
        for k in range(multiplier)
        for left, right in intervals
    ]


def exceptional_rays(
    neutral: Q, product_limit: int, ternary_units: bool
) -> list[tuple[int, int]]:
    answer: list[tuple[int, int]] = []
    # delta>neutral and delta<=9/(pq) force pq<product_limit.
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


def exposure(p: int, q: int) -> tuple[list[Interval], Q, Q]:
    a = merge(danger(p) + danger(q))
    base_bad = intersect(a, shift_half(a))
    r = subtract(shift_half(a), a)
    sigma = measure(base_bad)
    require(sigma == pair_mass(p, q), ((p, q), sigma, pair_mass(p, q)))
    require(measure(r) == measure(a) - sigma, (p, q, "exposure mass"))
    return r, measure(r), sigma


def finite_family(
    p: int, q: int, target: Q, ternary_units: bool
) -> dict[str, object]:
    r, mu_r, sigma = exposure(p, q)
    components = len(r)
    eta = (target - sigma) / 2 - mu_r / 7
    require(eta > 0, ((p, q), "nonpositive asymptotic margin", eta))
    cutoff = Q(components, 3) / eta
    require(cutoff.denominator > 1, ((p, q), "integral cutoff", cutoff))

    checked = 0
    maximum = Q(-1)
    maximizers: set[tuple[int, int, int]] = set()
    minimum_gap = Q(1)
    maximum_product = 0
    r_preimages: dict[int, list[Interval]] = {}
    d_cache: dict[int, list[Interval]] = {}

    # Products at least cutoff satisfy the analytic BV/Fourier estimate.
    # Enumerate exactly every remaining primitive presentation gcd(t,c)=1.
    for t in range(1, cutoff.numerator // cutoff.denominator + 1, 2):
        if ternary_units and t % 3 == 0:
            continue
        upper_c = (cutoff.numerator - 1) // (cutoff.denominator * t)
        if upper_c < 1:
            continue
        rt = r_preimages.setdefault(t, preimage(r, t))
        for c in range(1, upper_c + 1, 2):
            if (
                (ternary_units and c % 3 == 0)
                or gcd(t, c) != 1
                or c in (t * p, t * q)
            ):
                continue
            dc = d_cache.setdefault(c, danger(c))
            intersection_mass = measure(intersect(rt, dc))
            physical = sigma + 2 * intersection_mass
            checked += 1
            maximum_product = max(maximum_product, t * c)
            gap = target - physical
            require(gap >= 0, ((p, q), t, c, physical, "counterexample"))
            if physical > maximum:
                maximum = physical
                maximizers = {tuple(sorted((t * p, t * q, c)))}
            elif physical == maximum:
                maximizers.add(tuple(sorted((t * p, t * q, c))))
            minimum_gap = min(minimum_gap, gap)

            # Exact algebraic identity against a literal full danger union.
            # Run this independent route on every small presentation and on
            # all equality candidates; it is deliberately not needed for the
            # all-height reduction itself.
            if (t <= 5 and c <= 25) or gap == 0:
                full_a = merge(danger(t * p) + danger(t * q) + danger(c))
                literal = measure(intersect(full_a, shift_half(full_a)))
                require(literal == physical, ((p, q), t, c, literal, physical))

    return {
        "ray": (p, q),
        "sigma": sigma,
        "r_components": components,
        "mu_r": mu_r,
        "eta": eta,
        "cutoff": cutoff,
        "finite_presentations": checked,
        "maximum_product": maximum_product,
        "maximum": maximum,
        "maximizers": sorted(maximizers),
        "minimum_gap": minimum_gap,
    }


def main() -> None:
    print("LRC14_Q2_THREE_ODD_TAIL_TWO_DOMAIN_SHARP_CAP_PRIMARY")
    print("status=PROVED_RELATIVE_TO_BV_FOURIER_LEMMA_PLUS_FINITE_EXACT")
    # Exact pair-cap controls, including equality rays.
    general_pair_box = [
        (pair_mass(p, q), p, q)
        for p in range(1, 17, 2)
        for q in range(p + 2, 18, 2)
        if gcd(p, q) == 1 and p * q < 17
    ]
    require(max(general_pair_box) == (Q(4, 63), 1, 9), general_pair_box)
    ternary_pair_box = [
        (pair_mass(p, q), p, q)
        for p in range(1, 33, 2)
        for q in range(p + 2, 34, 2)
        if p % 3 and q % 3 and gcd(p, q) == 1 and p * q < 33
    ]
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

    domains = [
        (
            "full_odd",
            FULL_TARGET,
            Q(128, 621),
            44,
            False,
            [(1, 9), (1, 11), (1, 23), (3, 11)],
            (1, 9, 23),
        ),
        (
            "odd_3_unit",
            TERNARY_TARGET,
            Q(1, 11),
            99,
            True,
            [(1, 11), (1, 23), (1, 25), (1, 37), (5, 11)],
            (1, 7, 11),
        ),
    ]
    for name, target, neutral, product_limit, unit_filter, expected_rays, equality in domains:
        rays = exceptional_rays(neutral, product_limit, unit_filter)
        require(rays == expected_rays, (name, rays))
        print(f"DOMAIN name={name} target={target} neutral={neutral} exceptional_rays={rays}")
        results = [finite_family(*ray, target, unit_filter) for ray in rays]
        for row in results:
            print(name.upper() + "_FAMILY " + " ".join(f"{key}={value}" for key, value in row.items()))
        equality_shapes = sorted(
            {
                shape
                for row in results
                for shape in row["maximizers"]
                if row["maximum"] == target
            }
        )
        require(equality_shapes == [equality], (name, equality_shapes))
        require(max(row["maximum"] for row in results) == target, (name, results))
        print(f"{name.upper()}_primitive_equality_shapes={equality_shapes}")

    # Outside-domain hostile to the final constant and to silently dropping
    # the ternary-unit hypothesis.
    hostile = (1, 7, 9)
    hostile_mass = literal_physical(hostile)
    require(hostile_mass == Q(20, 147), hostile_mass)
    require(hostile_mass > TERNARY_TARGET, (hostile_mass, TERNARY_TARGET))
    first_hostiles = [
        tails
        for tails in combinations(range(1, 10, 2), 3)
        if literal_physical(tails) > TERNARY_TARGET
    ]
    require(first_hostiles == [hostile], first_hostiles)
    print(f"non_3_unit_hostile={hostile} mass={hostile_mass} excess={hostile_mass-TERNARY_TARGET}")

    # The stronger experimental envelope mu(F)<=4/49+max_pair Sigma also
    # genuinely needs the ternary-unit domain.  This is its first failure by
    # maximum tail height among positive odd triples.
    stronger_hostile = (5, 9, 11)
    stronger_mass = literal_physical(stronger_hostile)
    stronger_rhs = Q(4, 49) + max(
        pair_mass(*reduce_pair(a, b))
        for a, b in ((5, 9), (5, 11), (9, 11))
    )
    require(stronger_mass == Q(446, 3465), stronger_mass)
    require(stronger_rhs == Q(346, 2695), stronger_rhs)
    require(stronger_mass > stronger_rhs, (stronger_mass, stronger_rhs))
    stronger_first = []
    for tails in combinations(range(1, 12, 2), 3):
        pair_ceiling = Q(4, 49) + max(
            pair_mass(*reduce_pair(a, b)) for a, b in combinations(tails, 2)
        )
        if literal_physical(tails) > pair_ceiling:
            stronger_first.append(tails)
    require(stronger_first == [stronger_hostile], stronger_first)
    print(
        f"stronger_envelope_non_3_unit_hostile={stronger_hostile} "
        f"mass={stronger_mass} rhs={stronger_rhs} excess={stronger_mass-stronger_rhs}"
    )

    # Quotient/physical address distinction at equality.  Doubling maps the
    # physical pullback two-to-one to the quotient owner-cut set.
    for name, tails, physical_count, physical_width, quotient_count, quotient_width in (
        ("full_odd", (1, 9, 23), 16, Q(10, 483), 8, Q(20, 483)),
        ("odd_3_unit", (1, 7, 11), 12, Q(1, 77), 6, Q(2, 77)),
    ):
        equal_a = merge([arc for speed in tails for arc in danger(speed)])
        physical_cells = intersect(equal_a, shift_half(equal_a))
        require(len(physical_cells) == physical_count, (name, physical_cells))
        require(max(r - l for l, r in physical_cells) == physical_width, (name, physical_cells))
        doubled: list[Interval] = []
        for left, right in physical_cells:
            left, right = 2 * left, 2 * right
            if left >= 1:
                doubled.append((left - 1, right - 1))
            elif right > 1:
                doubled.extend(((left, Q(1)), (Q(0), right - 1)))
            else:
                doubled.append((left, right))
        quotient_cells = merge(doubled)
        require(len(quotient_cells) == quotient_count, (name, quotient_cells))
        require(max(r - l for l, r in quotient_cells) == quotient_width, (name, quotient_cells))
        print(
            f"{name.upper()}_ENDPOINT physical_components={physical_count} "
            f"physical_max_width={physical_width} quotient_components={quotient_count} "
            f"quotient_max_width={quotient_width}"
        )
    print(f"HAAR_GATES full_odd={FULL_TARGET} odd_3_unit={TERNARY_TARGET}")
    print("PASS")


if __name__ == "__main__":
    main()
