#!/usr/bin/env python3
"""Exact audits for the LRC(14) signed-wall septimal quadrature lane.

The physical signed ``u``-walls are labelled by

    (n, sigma),  n in Z/uZ, sigma in {-1,+1},
    x(n,sigma) = (14*n+sigma)/(14*u) mod 1.

Danger is always strict: ``||s*x|| < 1/14``.  The script independently
checks seven pieces of a proposed strengthening of THM-4348:

1. an exact generalized-CRT formula for arbitrary intersections and unions
   of speed masks on the two sign sheets;
2. a probability measure supported on a 7-adic lift of two base walls, under
   which every odd speed has danger mass at most 1/7;
3. the exact three-case mass of the even anchor, according to the relative
   7-adic valuations of ``h`` and ``u``;
4. the sharp seven-tail partition, re-anchoring controls, and explicit
   hostile rows delimiting what wall masks alone can prove.
5. exact ``7^{-r}`` intersections for masks at pairwise distinct lower
   7-valuations;
6. the common-height cyclic-AP owner formula, its complete first-shell
   classification, and finite-exact deeper-shell cover data; and
7. the first-shell coarse-defect escape lemma, including the exact ``8/171``
   permutation locus of the smallest wall tiler and its physical safe lift.

All arithmetic is integer or ``Fraction`` exact.  This is a lane-specific
research audit, not a proof of LRC(14).
"""

from __future__ import annotations

from fractions import Fraction
from functools import cache, reduce
from itertools import combinations, product
from math import gcd, lcm
from random import Random


SIGNS = (-1, 1)


def v7(number: int) -> int:
    """The exponent of 7 in a positive integer."""
    assert number > 0
    exponent = 0
    while number % 7 == 0:
        number //= 7
        exponent += 1
    return exponent


def centered_residue(number: int, modulus: int) -> int:
    """Return e with -modulus/2 < e <= modulus/2 and e = number mod modulus."""
    assert modulus > 0 and modulus % 2 == 0
    residue = number % modulus
    if residue > modulus // 2:
        residue -= modulus
    return residue


@cache
def wall_points(u: int) -> tuple[tuple[int, int], ...]:
    assert u > 0 and u % 2 == 1
    return tuple((n, sigma) for sigma in SIGNS for n in range(u))


def is_dangerous(u: int, speed: int, point: tuple[int, int]) -> bool:
    """Strict danger on an addressed wall, without floating point."""
    n, sigma = point
    numerator = speed * (14 * n + sigma)
    return abs(centered_residue(numerator, 14 * u)) < u


@cache
def direct_mask(u: int, speed: int) -> frozenset[tuple[int, int]]:
    return frozenset(point for point in wall_points(u) if is_dangerous(u, speed, point))


def is_dangerous_fraction(u: int, speed: int, point: tuple[int, int]) -> bool:
    """Independent rational-circle implementation of strict wall danger."""
    n, sigma = point
    phase = (Fraction(speed * (14 * n + sigma), 14 * u)) % 1
    return min(phase, 1 - phase) < Fraction(1, 14)


@cache
def quotient_classes(u: int, speed: int, sigma: int) -> tuple[int, tuple[int, ...]]:
    """Return q and the distinct n mod q classes in the strict mask.

    Put d=gcd(u,s), q=u/d, S=s/d.  Every strict centered residue y obeys

        -q < y < q,  y = S*sigma (mod 14),
        n = S^{-1} ((y-S*sigma)/14) (mod q).
    """
    assert u > 0 and u % 2 == 1 and speed > 0 and sigma in SIGNS
    divisor = gcd(u, speed)
    q = u // divisor
    reduced_speed = speed // divisor
    inverse = pow(reduced_speed, -1, q) if q > 1 else 0
    classes: list[int] = []
    for y in range(-q + 1, q):
        if (y - reduced_speed * sigma) % 14:
            continue
        z = (y - reduced_speed * sigma) // 14
        residue = (inverse * z) % q if q > 1 else 0
        classes.append(residue)
    assert len(classes) == len(set(classes))
    return q, tuple(sorted(classes))


@cache
def formula_mask(u: int, speed: int) -> frozenset[tuple[int, int]]:
    answer: set[tuple[int, int]] = set()
    for sigma in SIGNS:
        q, classes = quotient_classes(u, speed, sigma)
        class_set = set(classes)
        answer.update((n, sigma) for n in range(u) if n % q in class_set)
    return frozenset(answer)


def crt_compatible(moduli: tuple[int, ...], residues: tuple[int, ...]) -> bool:
    return all(
        (a - b) % gcd(q, r) == 0
        for (q, a), (r, b) in combinations(zip(moduli, residues), 2)
    )


@cache
def crt_intersection_count(u: int, speeds: tuple[int, ...]) -> int:
    """Exact arbitrary-fold intersection count on X_u via generalized CRT."""
    assert speeds
    total = 0
    for sigma in SIGNS:
        data = tuple(quotient_classes(u, speed, sigma) for speed in speeds)
        moduli = tuple(q for q, _ in data)
        lift_count = u // reduce(lcm, moduli, 1)
        for residues in product(*(classes for _, classes in data)):
            if crt_compatible(moduli, residues):
                total += lift_count
    return total


@cache
def direct_intersection_count(u: int, speeds: tuple[int, ...]) -> int:
    return len(set.intersection(*(set(direct_mask(u, speed)) for speed in speeds)))


def crt_union_count(u: int, speeds: tuple[int, ...]) -> int:
    """Inclusion-exclusion using the exact CRT intersections."""
    total = 0
    for size in range(1, len(speeds) + 1):
        sign = 1 if size % 2 else -1
        total += sign * sum(
            crt_intersection_count(u, subset)
            for subset in combinations(speeds, size)
        )
    return total


def direct_union_count(u: int, speeds: tuple[int, ...]) -> int:
    return len(set().union(*(direct_mask(u, speed) for speed in speeds)))


def septimal_support(u: int) -> tuple[tuple[int, int], ...]:
    """The 2*7^a support points, where u=7^a*v and 7 does not divide v."""
    assert u > 0 and u % 2 == 1
    a = v7(u)
    v = u // 7**a
    inverse_14 = pow(14, -1, v) if v > 1 else 0
    base = {sigma: (-sigma * inverse_14) % v if v > 1 else 0 for sigma in SIGNS}
    support = tuple(
        (base[sigma] + j * v, sigma)
        for sigma in SIGNS
        for j in range(7**a)
    )
    assert len(support) == 2 * 7**a == len(set(support))
    return support


def support_mass(u: int, speed: int) -> Fraction:
    support = septimal_support(u)
    hits = sum(is_dangerous(u, speed, point) for point in support)
    return Fraction(hits, len(support))


def anchor_expected_mass(u: int, h: int) -> Fraction:
    if v7(h) > v7(u):
        return Fraction(1)
    if v7(h) == v7(u):
        return Fraction()
    return Fraction(1, 7)


def physical_clearance(speeds: tuple[int, ...], time: Fraction) -> Fraction:
    distances = []
    for speed in speeds:
        residue = (speed * time) % 1
        distances.append(min(residue, 1 - residue))
    return min(distances)


def support_intersection_counts(
    u: int,
    speeds: tuple[int, ...],
    *,
    fraction_path: bool = False,
) -> tuple[int, int]:
    """Return exact hit counts on the negative and positive support sheets."""
    predicate = is_dangerous_fraction if fraction_path else is_dangerous
    counts = []
    support = septimal_support(u)
    for sigma in SIGNS:
        counts.append(
            sum(
                all(predicate(u, speed, point) for speed in speeds)
                for point in support
                if point[1] == sigma
            )
        )
    return tuple(counts)  # type: ignore[return-value]


UNIT_REPRESENTATIVES_14 = (-5, -3, -1, 1, 3, 5)


def unit_representative_14(speed: int) -> int:
    """The unique representative of a 7-unit odd residue in (-7,7)."""
    matches = tuple(
        delta for delta in UNIT_REPRESENTATIVES_14 if (delta - speed) % 14 == 0
    )
    assert len(matches) == 1
    return matches[0]


@cache
def owner_block_formula(k: int, reduced_speed: int, sigma: int = 1) -> frozenset[int]:
    """The common-height owner block in Z/7^k Z from the cyclic-AP formula."""
    assert k >= 1 and reduced_speed > 0 and reduced_speed % 2 and reduced_speed % 7
    assert sigma in SIGNS
    modulus = 7**k
    delta = unit_representative_14(reduced_speed)
    shift = (delta - reduced_speed) // 14
    inverse = pow(reduced_speed, -1, modulus)
    first = (-modulus - delta) // 14 + 1
    last = (modulus - delta - 1) // 14
    positive = frozenset(
        (inverse * (ell + shift)) % modulus for ell in range(first, last + 1)
    )
    assert len(positive) == modulus // 7
    if sigma == 1:
        return positive
    return frozenset((-address) % modulus for address in positive)


@cache
def owner_block_direct(k: int, reduced_speed: int, sigma: int = 1) -> frozenset[int]:
    """Independent centered-residue owner block on one sign sheet."""
    modulus = 7**k
    return frozenset(
        n for n in range(modulus) if is_dangerous(modulus, reduced_speed, (n, sigma))
    )


def bit_mask(addresses: frozenset[int]) -> int:
    answer = 0
    for address in addresses:
        answer |= 1 << address
    return answer


def unique_owner_masks(k: int) -> tuple[int, tuple[tuple[int, int], ...]]:
    """Return (U, ((bit mask, least positive speed representative),...))."""
    modulus = 7**k
    representatives: dict[int, int] = {}
    for speed in range(1, 14 * modulus, 2):
        if speed % 7:
            block = bit_mask(owner_block_formula(k, speed))
            representatives.setdefault(block, speed)
    return modulus, tuple(representatives.items())


def exact_owner_partitions(
    modulus: int,
    masks: tuple[tuple[int, int], ...],
    required_speed: int | None = None,
) -> tuple[tuple[int, ...], ...]:
    """Enumerate exact seven-block partitions by first-uncovered-point search."""
    full = (1 << modulus) - 1
    through_point: list[list[int]] = [[] for _ in range(modulus)]
    for index, (mask, _) in enumerate(masks):
        remaining = mask
        while remaining:
            low_bit = remaining & -remaining
            through_point[low_bit.bit_length() - 1].append(index)
            remaining -= low_bit

    answers: set[tuple[int, ...]] = set()

    def search(covered: int, chosen: tuple[int, ...]) -> None:
        if covered == full:
            if len(chosen) == 7:
                answers.add(tuple(sorted(masks[index][1] for index in chosen)))
            return
        if len(chosen) >= 7:
            return
        remaining = full ^ covered
        first_open = (remaining & -remaining).bit_length() - 1
        for index in through_point[first_open]:
            mask = masks[index][0]
            if not (mask & covered):
                search(covered | mask, chosen + (index,))

    if required_speed is None:
        search(0, ())
    else:
        required_index = next(
            index for index, (_, speed) in enumerate(masks) if speed == required_speed
        )
        search(masks[required_index][0], (required_index,))
    return tuple(sorted(answers))


def unsigned_residue(number: int, modulus: int) -> int:
    residue = number % modulus
    return min(residue, (-residue) % modulus)


def is_parallel_lift(partition: tuple[int, ...], modulus: int) -> bool:
    """Check the observed +/-R mod 2U parallel-lift normal form."""
    two_u = 2 * modulus
    common = unsigned_residue(partition[0], two_u)
    if any(unsigned_residue(speed, two_u) != common for speed in partition):
        return False
    oriented = tuple(
        speed % (14 * modulus)
        if speed % two_u == common
        else (-speed) % (14 * modulus)
        for speed in partition
    )
    return set(oriented) == {
        (common + 2 * modulus * j) % (14 * modulus) for j in range(7)
    }


def lift_owner_set_fraction(speed: int, y: Fraction) -> tuple[int, ...]:
    """Strictly dangerous lifts (y+j)/7, using rational-circle arithmetic."""
    answer = []
    for j in range(7):
        phase = (Fraction(speed, 7) * (y + j)) % 1
        if min(phase, 1 - phase) < Fraction(1, 14):
            answer.append(j)
    return tuple(answer)


def lift_owner_set_integer(speed: int, y: Fraction) -> tuple[int, ...]:
    """Independent modular-integer implementation of the same lift owners."""
    numerator, denominator = y.numerator, y.denominator
    phase_denominator = 7 * denominator
    answer = []
    for j in range(7):
        residue = (speed * (numerator + j * denominator)) % phase_denominator
        distance_numerator = min(residue, phase_denominator - residue)
        if 14 * distance_numerator < phase_denominator:
            answer.append(j)
    return tuple(answer)


def owner_breakpoints(speeds: tuple[int, ...]) -> tuple[Fraction, ...]:
    """All lift-owner transition points on [0,1], plus the two endpoints."""
    points = {Fraction(), Fraction(1)}
    for speed in speeds:
        assert speed % 2 and speed % 7
        points.update(Fraction(2 * m + 1, 2 * speed) for m in range(speed))
    return tuple(sorted(points))


def permutation_locus(
    speeds: tuple[int, ...],
) -> tuple[Fraction, tuple[tuple[Fraction, Fraction, tuple[int, ...]], ...], int]:
    """Exact open-cell sweep of the seven-lift owner-permutation locus."""
    boundaries = owner_breakpoints(speeds)
    intervals = []
    checked_points = 0
    for left, right in zip(boundaries, boundaries[1:]):
        midpoint = (left + right) / 2
        owners_fraction = tuple(lift_owner_set_fraction(speed, midpoint) for speed in speeds)
        owners_integer = tuple(lift_owner_set_integer(speed, midpoint) for speed in speeds)
        assert owners_fraction == owners_integer
        checked_points += 1
        if all(len(owner) == 1 for owner in owners_fraction):
            labels = tuple(owner[0] for owner in owners_fraction)
            if len(set(labels)) == 7:
                intervals.append((left, right, labels))
    for boundary in boundaries[:-1]:
        assert tuple(lift_owner_set_fraction(speed, boundary) for speed in speeds) == tuple(
            lift_owner_set_integer(speed, boundary) for speed in speeds
        )
        checked_points += 1
    measure = sum((right - left for left, right, _ in intervals), Fraction())
    return measure, tuple(intervals), checked_points


def audit_mask_formula() -> dict[str, int]:
    cases = 0
    endpoint_hits = 0
    for u in range(1, 100, 2):
        for speed in range(1, 14 * u + 1):
            direct = direct_mask(u, speed)
            formula = formula_mask(u, speed)
            assert direct == formula
            cases += 1
            endpoint_hits += sum(
                abs(centered_residue(speed * (14 * n + sigma), 14 * u)) == u
                for n, sigma in wall_points(u)
            )
    return {"cases": cases, "strict_endpoint_occurrences": endpoint_hits}


def audit_crt_calculus() -> dict[str, int]:
    intersection_cases = 0
    union_cases = 0
    maximum_fold = 0
    # Exhaust every 1--4 subset of a compact residue bank for each odd u.
    # The generalized-CRT formula itself is unbounded; this audit universe is
    # deliberately small enough for three independent interpreter replays.
    for u in range(1, 30, 2):
        bank = tuple(range(1, min(14 * u, 9) + 1))
        for size in range(1, min(4, len(bank)) + 1):
            for speeds in combinations(bank, size):
                assert crt_intersection_count(u, speeds) == direct_intersection_count(u, speeds)
                intersection_cases += 1
                maximum_fold = max(maximum_fold, size)
        # Unions are more expensive; use every triple plus deterministic 4--6 folds.
        for speeds in combinations(bank, min(3, len(bank))):
            assert crt_union_count(u, speeds) == direct_union_count(u, speeds)
            union_cases += 1

    rng = Random(4348)
    for size in range(4, 7):
        for _ in range(80):
            u = rng.randrange(1, 40, 2)
            speeds = tuple(rng.sample(range(1, 14 * u + 1), size))
            assert crt_intersection_count(u, speeds) == direct_intersection_count(u, speeds)
            assert crt_union_count(u, speeds) == direct_union_count(u, speeds)
            intersection_cases += 1
            union_cases += 1
            maximum_fold = max(maximum_fold, size)
    return {
        "intersection_cases": intersection_cases,
        "union_cases": union_cases,
        "maximum_fold": maximum_fold,
    }


def audit_support_and_anchor() -> dict[str, int]:
    odd_speed_cases = 0
    anchor_cases = 0
    lift_fibres = 0
    exact_one_seventh_odd = 0
    valuation_safe_odd = 0
    independent_anchor_cases = 0
    for u in range(1, 400, 2):
        support = septimal_support(u)
        a = v7(u)
        v = u // 7**a

        # The support is precisely the uniform iterated preimage of its two
        # base addresses under (n,sigma) -> (n mod w,sigma).
        for exponent in range(a):
            w = v * 7**exponent
            lower = set(septimal_support(w))
            upper = set(septimal_support(7 * w))
            projected = {(n % w, sigma) for n, sigma in upper}
            assert projected == lower
            for point in lower:
                children = [(n, sigma) for n, sigma in upper if (n % w, sigma) == point]
                assert len(children) == 7
                lift_fibres += 1

        for speed in range(1, min(14 * u, 1200) + 1, 2):
            mass = support_mass(u, speed)
            expected = Fraction(1, 7) if v7(speed) < a else Fraction()
            assert mass == expected
            if expected:
                exact_one_seventh_odd += 1
            else:
                valuation_safe_odd += 1
            odd_speed_cases += 1

        for h in range(1, min(5 * u, 1600) + 1):
            observed = support_mass(u, 2 * h)
            expected = anchor_expected_mass(u, h)
            assert observed == expected
            anchor_cases += 1
            if u < 120 and h <= min(3 * u, 400):
                fraction_hits = sum(
                    is_dangerous_fraction(u, 2 * h, point) for point in support
                )
                assert Fraction(fraction_hits, len(support)) == expected
                independent_anchor_cases += 1
    return {
        "odd_speed_cases": odd_speed_cases,
        "anchor_cases": anchor_cases,
        "independent_fraction_anchor_cases": independent_anchor_cases,
        "lift_fibres": lift_fibres,
        "odd_exact_one_seventh": exact_one_seventh_odd,
        "odd_valuation_safe": valuation_safe_odd,
    }


def audit_sharp_seven_partition() -> dict[str, int]:
    cases = 0
    for u in range(7, 800, 14):
        speeds = tuple(1 + 2 * j * u for j in range(7))
        masks = tuple(direct_mask(u, speed) for speed in speeds)
        support = set(septimal_support(u))
        restricted = tuple(mask & support for mask in masks)
        assert all(len(mask) == len(support) // 7 for mask in restricted)
        assert sum(map(len, restricted)) == len(support)
        assert set().union(*restricted) == support
        assert all(not (left & right) for left, right in combinations(restricted, 2))
        # Independent Fraction path: each support point has one and only one
        # strict hit, with no centered-residue helper involved.
        assert all(
            sum(is_dangerous_fraction(u, speed, point) for speed in speeds) == 1
            for point in support
        )
        cases += 1
    return {"u_cases": cases}


def audit_distinct_valuation_intersections() -> dict[str, int]:
    """Audit the exact 7^{-r} law on a declared composite-u universe."""
    direct_cases = 0
    independent_fraction_cases = 0
    unit_multipliers = (1, 3, 5, 9)
    coprime_parts = (1, 3, 5, 11)
    for a in range(1, 5):
        for v in coprime_parts:
            u = 7**a * v
            for fold in range(1, a + 1):
                for valuations in combinations(range(a), fold):
                    for units in product(unit_multipliers, repeat=fold):
                        speeds = tuple(
                            unit * 7**valuation
                            for unit, valuation in zip(units, valuations)
                        )
                        expected_sheet_count = 7 ** (a - fold)
                        assert support_intersection_counts(u, speeds) == (
                            expected_sheet_count,
                            expected_sheet_count,
                        )
                        direct_cases += 1
                        # A second Fraction-circle path covers every fold
                        # through height three at one nontrivial base quotient.
                        if v == 11 and a <= 3:
                            assert support_intersection_counts(
                                u, speeds, fraction_path=True
                            ) == (expected_sheet_count, expected_sheet_count)
                            independent_fraction_cases += 1
    return {
        "direct_cases": direct_cases,
        "independent_fraction_cases": independent_fraction_cases,
        "maximum_fold": 4,
    }


def audit_owner_blocks_and_finite_rigidity() -> dict[str, object]:
    """Audit the proved owner formula and explicitly scoped finite data."""
    expected_classes = {
        0: {1, 97, 3, 95, 5, 93},
        1: {13, 85, 33, 65, 39, 59},
        2: {17, 81, 27, 71, 37, 61},
        3: {9, 89, 25, 73, 41, 57},
        4: {19, 79, 31, 67, 43, 55},
        5: {11, 87, 29, 69, 47, 51},
        6: {15, 83, 23, 75, 45, 53},
    }
    observed_classes: dict[int, set[int]] = {label: set() for label in range(7)}
    formula_direct_checks = 0
    for speed in range(1, 98, 2):
        if speed % 7 == 0:
            continue
        for sigma in SIGNS:
            assert owner_block_formula(1, speed, sigma) == owner_block_direct(1, speed, sigma)
            formula_direct_checks += 1
        delta = unit_representative_14(speed)
        kappa = (pow(speed, -1, 7) * ((delta - speed) // 14)) % 7
        assert owner_block_formula(1, speed, 1) == frozenset({kappa})
        observed_classes[kappa].add(speed)
    assert observed_classes == expected_classes

    canonical_tiler = (1, 9, 11, 13, 15, 17, 19)
    canonical_owners = tuple(next(iter(owner_block_formula(1, speed))) for speed in canonical_tiler)
    assert canonical_owners == (0, 3, 5, 1, 6, 2, 4)
    for sigma in SIGNS:
        blocks = tuple(owner_block_formula(1, speed, sigma) for speed in canonical_tiler)
        assert set().union(*map(set, blocks)) == set(range(7))
        assert all(not (left & right) for left, right in combinations(blocks, 2))

    finite_data: dict[int, tuple[int, int]] = {}
    for k, expected_mask_count, expected_partition_count in (
        (2, 147, 21),
        (3, 1029, 147),
    ):
        modulus, masks = unique_owner_masks(k)
        assert len(masks) == expected_mask_count
        # Independent centered-residue comparison for every distinct block,
        # on both sign sheets.
        for _, speed in masks:
            for sigma in SIGNS:
                assert owner_block_formula(k, speed, sigma) == owner_block_direct(
                    k, speed, sigma
                )
                formula_direct_checks += 1
        partitions = exact_owner_partitions(modulus, masks)
        assert len(partitions) == expected_partition_count
        assert all(is_parallel_lift(partition, modulus) for partition in partitions)
        assert len(
            {unsigned_residue(partition[0], 2 * modulus) for partition in partitions}
        ) == expected_partition_count
        finite_data[k] = (len(masks), len(partitions))

    # At k=4 this is deliberately a targeted completion through the block
    # containing speed 1, not an exhaustive classification of all covers.
    modulus, masks = unique_owner_masks(4)
    assert len(masks) == 7203
    speed_one_mask = next(mask for mask, speed in masks if speed == 1)
    disjoint_partners = sum(not (mask & speed_one_mask) for mask, _ in masks)
    assert disjoint_partners == 46
    completions = exact_owner_partitions(modulus, masks, required_speed=1)
    expected_completion = (1, 4801, 4803, 9603, 9605, 14405, 14407)
    assert completions == (expected_completion,)
    assert is_parallel_lift(expected_completion, modulus)
    for speed in expected_completion:
        for sigma in SIGNS:
            assert owner_block_formula(4, speed, sigma) == owner_block_direct(
                4, speed, sigma
            )
            formula_direct_checks += 1
    assert all(
        sum(is_dangerous_fraction(modulus, speed, point) for speed in expected_completion) == 1
        for point in wall_points(modulus)
    )

    return {
        "k1_unit_residues": sum(map(len, observed_classes.values())),
        "k1_class_sizes": tuple(len(observed_classes[label]) for label in range(7)),
        "k1_canonical_tiler": canonical_tiler,
        "k1_positive_owners": canonical_owners,
        "formula_direct_checks": formula_direct_checks,
        "k2_unique_masks_and_covers": finite_data[2],
        "k3_unique_masks_and_covers": finite_data[3],
        "k4_target_masks": len(masks),
        "k4_speed1_disjoint_partners": disjoint_partners,
        "k4_speed1_completions": len(completions),
        "k4_completion": expected_completion,
    }


def coarse_danger_boundaries(speed: int) -> set[Fraction]:
    """Endpoints of {y: ||speed*y||<1/14} on [0,1)."""
    return {
        Fraction((14 * m + sign) % (14 * speed), 14 * speed)
        for m in range(speed)
        for sign in SIGNS
    }


def audit_coarse_defect_and_retained_hostile() -> dict[str, object]:
    """Audit the exact first-shell defect closure and its wall-only hostile."""
    lower = (1, 9, 11, 13, 15, 17, 19)
    upper = (7, 21, 35, 63, 77)
    anchor = 14
    row = (anchor,) + lower + upper
    support = wall_points(7)

    # On the critical walls, the seven lower tails partition and every
    # retained upper-cone speed is safe: n_a=5 and there is no deeper tail.
    assert all(
        sum(is_dangerous(7, speed, point) for speed in lower) == 1
        for point in support
    )
    assert all(
        sum(is_dangerous_fraction(7, speed, point) for speed in lower) == 1
        for point in support
    )
    assert all(not is_dangerous(7, speed, point) for speed in (anchor,) + upper for point in support)

    witness = Fraction(15, 182)
    clearance = physical_clearance(row, witness)
    binders = tuple(
        speed for speed in row if physical_clearance((speed,), witness) == clearance
    )
    assert clearance == Fraction(1, 14)
    assert binders == (13,)

    permutation_measure, intervals, owner_points_checked = permutation_locus(lower)
    expected_intervals = (
        (Fraction(1, 18), Fraction(3, 38), (0, 3, 5, 1, 6, 2, 4)),
        (Fraction(35, 38), Fraction(17, 18), (6, 3, 1, 5, 0, 4, 2)),
    )
    assert intervals == expected_intervals
    assert permutation_measure == Fraction(8, 171)
    defect_measure = 1 - permutation_measure
    assert defect_measure == Fraction(163, 171) > Fraction(6, 7)

    # Same X_7 wall classification, different physical owner evolution: this
    # parallel tiler lies beyond the scalar defect-measure certificate.
    parallel_tiler = (1, 15, 29, 43, 57, 71, 85)
    for sigma in SIGNS:
        blocks = tuple(owner_block_formula(1, speed, sigma) for speed in parallel_tiler)
        assert set().union(*map(set, blocks)) == set(range(7))
        assert all(not (left & right) for left, right in combinations(blocks, 2))
    parallel_measure, _, parallel_points_checked = permutation_locus(parallel_tiler)
    assert parallel_measure == Fraction(7122, 46835) > Fraction(1, 7)

    # Independent pointwise audits of the coarse-defect lift implication.
    rng = Random(4370)
    profiles = []
    odd_bank = tuple(range(1, 40, 2))
    for _ in range(40):
        profiles.append((rng.randrange(1, 31), tuple(rng.sample(odd_bank, 5))))
    lift_points_checked = 0
    for h_reduced, upper_reduced in profiles:
        coarse_speeds = (2 * h_reduced,) + upper_reduced
        physical_row = (14 * h_reduced,) + lower + tuple(
            7 * speed for speed in upper_reduced
        )
        points = set(owner_breakpoints(lower))
        for speed in coarse_speeds:
            points.update(coarse_danger_boundaries(speed))
        ordered = sorted(points)
        test_points = set(ordered[:-1])
        test_points.update((left + right) / 2 for left, right in zip(ordered, ordered[1:]))
        for y in test_points:
            killed = set().union(
                *(set(lift_owner_set_fraction(speed, y)) for speed in lower)
            )
            coarse_bad = any(
                physical_clearance((speed,), y) < Fraction(1, 14)
                for speed in coarse_speeds
            )
            fibre_bad = all(
                physical_clearance(physical_row, (y + j) / 7) < Fraction(1, 14)
                for j in range(7)
            )
            assert fibre_bad == (len(killed) == 7 or coarse_bad)
            if len(killed) < 7 and not coarse_bad:
                free_lift = min(set(range(7)) - killed)
                time = (y + free_lift) / 7
                assert physical_clearance(physical_row, time) >= Fraction(1, 14)
            lift_points_checked += 1

    return {
        "critical_wall_lower_partition_points": len(support),
        "critical_exact_height_tails": len(upper),
        "critical_deeper_tails": 0,
        "actual_safe_time": str(witness),
        "actual_clearance": str(clearance),
        "actual_binders": binders,
        "permutation_intervals": tuple(
            (str(left), str(right), labels) for left, right, labels in intervals
        ),
        "permutation_measure": str(permutation_measure),
        "defect_measure": str(defect_measure),
        "owner_points_checked": owner_points_checked,
        "parallel_tiler": parallel_tiler,
        "parallel_permutation_measure": str(parallel_measure),
        "parallel_owner_points_checked": parallel_points_checked,
        "coarse_profiles": len(profiles),
        "coarse_lift_points_checked": lift_points_checked,
    }


def audit_reanchor_and_hostiles() -> dict[str, object]:
    # Headline valuation control with no tail at the anchor's exact height.
    # Here nu_7(h)=1, six tails lie below it, and six lie strictly above it.
    # The auxiliary wall scale u=7 need not itself occur as a row speed.
    valuation_tails = (1, 3, 5, 9, 11, 13, 49, 147, 245, 343, 441, 539)
    valuation_u = 7
    assert all(v7(speed) != v7(420) for speed in valuation_tails)
    assert sum(v7(speed) < v7(420) for speed in valuation_tails) == 6
    valuation_safe = [
        point
        for point in wall_points(valuation_u)
        if not any(is_dangerous(valuation_u, speed, point) for speed in (840,) + valuation_tails)
    ]
    assert valuation_safe
    valuation_point = valuation_safe[0]
    valuation_time = Fraction(14 * valuation_point[0] + valuation_point[1], 14 * valuation_u) % 1
    valuation_clearance = physical_clearance((840,) + valuation_tails, valuation_time)
    assert valuation_clearance == Fraction(13, 98)

    # A physical 2+12 row in the h=420 resonance branch at the original u=1.
    # Re-anchor at v=7: ten of the eleven other odd tails are odd multiples
    # of 7 and hence wall-safe; only speed 1 is a nonmultiple.
    h = 420
    tails = (1, 7, 21, 35, 63, 77, 91, 105, 119, 133, 147, 161)
    reanchor = 7
    others = tuple(speed for speed in tails if speed != reanchor)
    multiples = tuple(speed for speed in others if speed % reanchor == 0)
    nonmultiples = tuple(speed for speed in others if speed % reanchor)
    assert len(tails) == 12 and len(set(tails)) == 12
    assert v7(reanchor) == v7(h)
    assert len(multiples) == 10 and len(nonmultiples) == 1
    assert all((speed // reanchor) % 2 == 1 for speed in multiples)
    safe_points = [
        point
        for point in wall_points(reanchor)
        if not any(is_dangerous(reanchor, speed, point) for speed in (2 * h,) + tails)
    ]
    assert safe_points
    point = safe_points[0]
    n, sigma = point
    time = Fraction(14 * n + sigma, 14 * reanchor) % 1
    clearance = physical_clearance((2 * h,) + tails, time)
    assert clearance >= Fraction(1, 14)
    binders = tuple(
        speed
        for speed in (2 * h,) + tails
        if physical_clearance((speed,), time) == clearance
    )
    assert binders == (7, 91, 105)

    # A strict wall-cover hostile: an anchor plus four odd residual masks can
    # already cover X_5. Adding 14u to a speed preserves this wall mask.
    hostile_u = 5
    hostile_h = 7
    hostile_residuals = (71, 81, 83, 87)
    assert all(speed > 2 * hostile_h for speed in hostile_residuals)
    assert direct_union_count(hostile_u, (2 * hostile_h,) + hostile_residuals) == 2 * hostile_u
    for speed in hostile_residuals:
        assert direct_mask(hostile_u, speed) == direct_mask(hostile_u, speed % (14 * hostile_u))

    # A primitive anchor-plus-five cover at the canonical physical h=420
    # scale. It blocks any universal improvement from five odd residuals to
    # four without extra structure. Each tail can be lifted independently by
    # a multiple of 14u while preserving its wall mask.
    five_u = 11
    five_representatives = (5, 23, 61, 49, 65)
    five_lifted = tuple(speed + 14 * five_u * (j + 1) for j, speed in enumerate(five_representatives))
    assert gcd(2 * h, *five_lifted, five_u) == 1
    assert direct_union_count(five_u, (2 * h,) + five_lifted) == 2 * five_u
    assert all(
        direct_mask(five_u, speed) == direct_mask(five_u, representative)
        for speed, representative in zip(five_lifted, five_representatives)
    )

    return {
        "no_matching_height_lower_tails": 6,
        "no_matching_height_safe_time": str(valuation_time),
        "no_matching_height_clearance": str(valuation_clearance),
        "reanchor_h": h,
        "reanchor_tail_count": len(tails),
        "reanchor_multiples": len(multiples),
        "reanchor_nonmultiples": len(nonmultiples),
        "safe_wall": point,
        "safe_time": str(time),
        "clearance": str(clearance),
        "binding_speeds": binders,
        "anchor_plus_four_cover": (hostile_u, hostile_h, hostile_residuals),
        "anchor_plus_five_cover": (five_u, h, five_lifted),
    }


def main() -> None:
    print("LRC14 SIGNED-WALL SEPTIMAL QUADRATURE / REANCHOR AUDIT")
    print("SCOPE=FINITE-EXACT CONTROLS; STRICT danger ||s*x||<1/14; LRC14_OPEN")
    print("mask_formula", audit_mask_formula())
    print("crt_calculus", audit_crt_calculus())
    print("septimal_measure", audit_support_and_anchor())
    print("sharp_seven_partition", audit_sharp_seven_partition())
    print("distinct_valuation_intersections", audit_distinct_valuation_intersections())
    print("owner_blocks_and_finite_rigidity", audit_owner_blocks_and_finite_rigidity())
    print("coarse_defect_and_retained_hostile", audit_coarse_defect_and_retained_hostile())
    print("reanchor_and_hostiles", audit_reanchor_and_hostiles())
    print("PASS")


if __name__ == "__main__":
    main()
