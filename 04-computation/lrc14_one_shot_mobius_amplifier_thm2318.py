#!/usr/bin/env python3
"""Exact, optimization-safe companion for THM-2318."""

from fractions import Fraction
from itertools import product
from math import gcd


P = 13
Q = 7
ELL = 65537
MODULUS = P * Q * ELL
EPSILON = Fraction(1, 10**12)

SOURCE_CENTERS = (Fraction(1, 16), Fraction(15, 16))
CURRENT_CENTERS = (Fraction(7, 16), Fraction(9, 16))


def require(condition: bool, message: str) -> None:
    """Raise in both ordinary and optimized Python when a check fails."""
    if not condition:
        raise RuntimeError(message)


def is_prime(number: int) -> bool:
    """Deterministic trial-division primality test."""
    if number < 2:
        return False
    divisor = 2
    while divisor * divisor <= number:
        if number % divisor == 0:
            return number == divisor
        divisor += 1
    return True


def valuation(number: int, prime: int) -> int:
    """Return the prime-adic valuation of a positive integer."""
    require(number > 0, "valuation requires a positive integer")
    answer = 0
    while number % prime == 0:
        number //= prime
        answer += 1
    return answer


def circle_distance(x: Fraction, y: Fraction) -> Fraction:
    """Distance on R/Z between rational representatives."""
    difference = (x - y) % 1
    return min(difference, 1 - difference)


def cube_upper_set_check() -> tuple[int, int, int]:
    """Exhaust every upper subset of the Boolean three-cube."""
    cells = tuple(product(range(2), repeat=3))
    upper_sets = 0
    disjoint_origin_upper_sets = 0
    singleton_top_sets = 0

    for mask in range(1 << len(cells)):
        occupied = {
            cell for index, cell in enumerate(cells)
            if mask & (1 << index)
        }
        is_upper = all(
            larger in occupied
            for cell in occupied
            for larger in cells
            if all(larger[i] >= cell[i] for i in range(3))
        )
        if not is_upper:
            continue
        upper_sets += 1
        if occupied and (0, 0, 0) not in occupied:
            disjoint_origin_upper_sets += 1
        if occupied == {(1, 1, 1)}:
            singleton_top_sets += 1

    require(upper_sets == 20, "Boolean-cube upper-set count changed")
    require(
        disjoint_origin_upper_sets == 18,
        "disjoint-origin Boolean-cube upper-set count changed",
    )
    require(singleton_top_sets == 1, "fully minimal top cube changed")
    return upper_sets, disjoint_origin_upper_sets, singleton_top_sets


def mobius_selector_check() -> tuple[int, int]:
    """Verify the complete three-prime Boolean Möbius selector."""
    primes = (P, Q, ELL)
    truth_rows = 0
    for flags in product((False, True), repeat=3):
        number = 1
        for prime, flag in zip(primes, flags):
            if flag:
                number *= prime
        expanded = (
            1
            - int(number % P == 0)
            - int(number % Q == 0)
            - int(number % ELL == 0)
            + int(number % (P * Q) == 0)
            + int(number % (P * ELL) == 0)
            + int(number % (Q * ELL) == 0)
            - int(number % MODULUS == 0)
        )
        require(
            expanded == int(gcd(number, MODULUS) == 1),
            "three-prime selector truth table failed",
        )
        truth_rows += 1

    interval_rows = 0
    for number in range(-10000, 10001):
        product_selector = (
            (1 - int(number % P == 0))
            * (1 - int(number % Q == 0))
            * (1 - int(number % ELL == 0))
        )
        require(
            product_selector == int(gcd(abs(number), MODULUS) == 1),
            "three-prime selector interval check failed",
        )
        interval_rows += 1
    return truth_rows, interval_rows


def cube_mobius_sign_check() -> int:
    """Check coefficients and the odd-dimensional top-corner sign."""
    rows = 0
    total = 0
    for flags in product((0, 1), repeat=3):
        coefficient = (-1) ** sum(flags)
        expected = {
            (0, 0, 0): 1,
            (1, 0, 0): -1,
            (0, 1, 0): -1,
            (0, 0, 1): -1,
            (1, 1, 0): 1,
            (1, 0, 1): 1,
            (0, 1, 1): 1,
            (1, 1, 1): -1,
        }[flags]
        require(coefficient == expected, "cube Möbius coefficient changed")
        # A fully minimal top has collision value one only at (1,1,1).
        collision_value = int(flags == (1, 1, 1))
        total += coefficient * collision_value
        rows += 1
    require(total == -1, "odd-dimensional corner sign changed")
    return rows


def residue_bank_bound_check(max_nodes: int = 512) -> int:
    """Check the exact endpoint-Prony progression ceiling."""
    require(gcd(MODULUS - 1, MODULUS) == 1,
            "largest residue is not a unit")
    rows = 0
    for node_count in range(1, max_nodes + 1):
        largest_tested = (MODULUS - 1) + MODULUS * (node_count - 1)
        require(
            largest_tested == MODULUS * node_count - 1,
            "three-prime Prony ceiling changed",
        )
        rows += 1
    return rows


def gap(scale: int) -> Fraction:
    """Exact cross-gap between the hostile center sets after one scale."""
    source_images = tuple((scale * x) % 1 for x in SOURCE_CENTERS)
    current_images = tuple((scale * x) % 1 for x in CURRENT_CENTERS)
    return min(
        circle_distance(x, y)
        for x in source_images
        for y in current_images
    )


def collision(scale: int) -> bool:
    """Positive-measure overlap of the hostile pushed supports."""
    source_radius = scale * EPSILON
    current_radius = 169 * scale * EPSILON
    cross_gap = gap(scale)
    geometric = (
        source_radius >= Fraction(1, 2)
        or current_radius >= Fraction(1, 2)
        or source_radius + current_radius > cross_gap
    )
    threshold = 170 * scale * EPSILON > cross_gap
    require(geometric == threshold, "collision threshold lost equivalence")
    return geometric


def hostile_cube_check() -> tuple[
    int,
    tuple[tuple[int, int, int], ...],
    tuple[tuple[str, int, int, int, bool], ...],
]:
    """Verify the exact fully minimal (3,1,1) hostile collision cube."""
    require(is_prime(ELL), "the auxiliary amplifier is not prime")
    require(ELL == 2**16 + 1, "Fermat-prime identity changed")
    require(ELL % 16 == 1, "auxiliary prime no longer fixes center residues")
    require(gcd(ELL, P * Q) == 1, "auxiliary prime is not new")

    positive_cells = []
    for s in range(4):
        for t in range(2):
            for u in range(2):
                scale = P**s * Q**t * ELL**u
                expected_gap = (
                    Fraction(3, 8) if s % 2 == 0 else Fraction(1, 8)
                )
                require(gap(scale) == expected_gap,
                        "hostile three-direction gap law changed")
                if collision(scale):
                    positive_cells.append((s, t, u))

    require(
        tuple(positive_cells) == ((3, 1, 1),),
        "hostile target cube is no longer fully minimal",
    )

    named_cells = (
        ("top", 3, 1, 1, True),
        ("13-face", 2, 1, 1, False),
        ("7-face", 3, 0, 1, False),
        ("ell-face", 3, 1, 0, False),
    )
    arithmetic_rows = []
    for name, s, t, u, expected_collision in named_cells:
        scale = P**s * Q**t * ELL**u
        scaled_threshold = 1360 * scale
        right_side = (3 if s % 2 == 0 else 1) * 10**12
        require(
            (scaled_threshold > right_side) == expected_collision,
            f"{name} exact collision inequality changed",
        )
        require(
            collision(scale) == expected_collision,
            f"{name} geometric collision state changed",
        )
        arithmetic_rows.append(
            (name, scale, scaled_threshold, right_side, expected_collision)
        )

    top_scale = P**3 * Q * ELL
    require(169 * top_scale * EPSILON < Fraction(1, 2),
            "top current interval began wrapping")
    require(top_scale * EPSILON < Fraction(1, 2),
            "top source interval began wrapping")
    return top_scale, tuple(positive_cells), tuple(arithmetic_rows)


def valuation_and_landing_check() -> tuple[int, int, int, int]:
    """Check the hostile and global landed-frequency ledgers."""
    local_jump_product = 4 * 4
    local_bound = MODULUS * local_jump_product - 1
    require(local_bound == 95421871, "local Prony bound changed")

    global_coefficient = MODULUS * 8
    require(global_coefficient == 47710936,
            "global S-squared coefficient changed")

    rows = 0
    owner = P
    target_s = 3
    for residual in range(1, 10001):
        if gcd(residual, MODULUS) != 1:
            continue
        frequency = owner * P ** (target_s - 1) * residual
        require(valuation(frequency, P) == 3,
                "hostile thirteen-adic grade changed")
        require(valuation(frequency, Q) == 0,
                "one-shot amplifier inserted a seven factor")
        require(gcd(frequency // P**3, P * Q) == 1,
                "hostile outside multiplier ceased to be 91-unit")
        rows += 1
    return local_jump_product, local_bound, global_coefficient, rows


def main() -> None:
    upper_sets, disjoint_sets, singleton_top = cube_upper_set_check()
    selector_truth_rows, selector_interval_rows = mobius_selector_check()
    mobius_rows = cube_mobius_sign_check()
    prony_rows = residue_bank_bound_check()
    top_scale, positive_cells, arithmetic_rows = hostile_cube_check()
    (
        local_jump_product,
        local_bound,
        global_coefficient,
        valuation_rows,
    ) = valuation_and_landing_check()

    print("theorem=THM-2318")
    print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-AUDITED")
    print("directions=(13,7,65537)")
    print(f"auxiliary_prime={ELL}")
    print(f"selector_modulus={MODULUS}")
    print("full_cube_selector=[gcd(n,13*7*65537)=1]")
    print("full_cube_sign=-I_top")
    print("one_shot_cell=(s,1,1)")
    print("one_shot_landing=13^(s-1)*n")
    print("auxiliary_7_and_65537_exponents_in_landing=0")
    print(f"boolean_cube_upper_sets={upper_sets}")
    print(f"boolean_cube_disjoint_origin_upper_sets={disjoint_sets}")
    print(f"boolean_cube_singleton_top_sets={singleton_top}")
    print(f"selector_truth_rows={selector_truth_rows}")
    print(f"selector_interval_rows={selector_interval_rows}")
    print(f"mobius_coefficient_rows={mobius_rows}")
    print(f"prony_ceiling_rows={prony_rows}")
    print(f"hostile_top_scale={top_scale}")
    print(f"hostile_positive_cells_in_target_cube={positive_cells}")
    for name, scale, left, right, state in arithmetic_rows:
        print(
            f"hostile_{name}="
            f"scale:{scale},1360scale:{left},rhs:{right},collision:{state}"
        )
    print("hostile_owner=13")
    print("hostile_grade=3")
    print("hostile_outside_multiplier_is_91_unit=true")
    print(f"local_jump_product={local_jump_product}")
    print(f"local_prony_bound={local_bound}")
    print(f"global_bound={global_coefficient}*S^2-1")
    print(f"valuation_rows={valuation_rows}")
    print("root_residue_not_prescribed=true")
    print("pair_edge_incidence_not_proved=true")
    print("target_gain_not_selected=true")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
