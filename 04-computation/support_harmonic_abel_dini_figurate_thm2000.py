#!/usr/bin/env python3
"""Exact/numerically certified referee for THM-2000.

Paper topology and classical convergence theorems stay in the theorem text.
This program freezes the finite Abel identities, support/multiplicity collision
ledger, dyadic block bounds, polygonal and master-figurate algebra, exact
simplex tails, Abel--Dini integer lifts, Egyptian refinements, Sylvester
remainders, and the tournament maximum-cyclic-triangle reciprocal constant.
"""

from __future__ import annotations

import ast
import hashlib
from fractions import Fraction
from itertools import combinations
from math import ceil, comb, factorial, floor, isqrt
from pathlib import Path

import mpmath as mp


F = Fraction
mp.mp.dps = 60


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def harmonic(number: int) -> F:
    return sum((F(1, index) for index in range(1, number + 1)), F(0))


def reciprocal_mass(values: list[int], *, deduplicate: bool) -> F:
    positive = [value for value in values if value > 0]
    if deduplicate:
        positive = sorted(set(positive))
    return sum((F(1, value) for value in positive), F(0))


def collision_tax(values: list[int]) -> F:
    return reciprocal_mass(values, deduplicate=False) - reciprocal_mass(
        values, deduplicate=True
    )


def support_collision_audit() -> tuple[int, tuple[tuple[str, F], ...]]:
    prefixes = {
        "factorial_0_to_10": [factorial(index) for index in range(11)],
        "fibonacci_1_to_18": [
            1,
            1,
            2,
            3,
            5,
            8,
            13,
            21,
            34,
            55,
            89,
            144,
            233,
            377,
            610,
            987,
            1597,
            2584,
        ],
        "catalan_0_to_10": [comb(2 * index, index) // (index + 1) for index in range(11)],
        "legacy_A000568_prefix": [1, 1, 1, 2, 4, 12, 56, 456, 6880],
        "legacy_A002854_prefix": [1, 1, 2, 3, 7, 16, 54, 243],
        "legacy_A000571_prefix": [1, 1, 1, 2, 4, 9, 22, 59],
    }
    expected = {
        "factorial_0_to_10": F(1),
        "fibonacci_1_to_18": F(1),
        "catalan_0_to_10": F(1),
        "legacy_A000568_prefix": F(2),
        "legacy_A002854_prefix": F(1),
        "legacy_A000571_prefix": F(2),
    }
    rows: list[tuple[str, F]] = []
    for name, values in prefixes.items():
        tax = collision_tax(values)
        require(tax == expected[name], f"collision tax changed for {name}")
        require(
            reciprocal_mass(values, deduplicate=False)
            == reciprocal_mass(values, deduplicate=True) + tax,
            f"collision identity failed for {name}",
        )
        rows.append((name, tax))

    labeled = [2 ** comb(index, 2) for index in range(1, 12)]
    switching = [2 ** comb(index - 1, 2) for index in range(1, 13)]
    require(set(labeled) == set(switching), "tournament/switching supports differ")
    require(
        reciprocal_mass(switching, deduplicate=False)
        - reciprocal_mass(labeled, deduplicate=False)
        == 1,
        "indexed tournament/switching collision tax is not one",
    )
    require(
        reciprocal_mass(switching, deduplicate=True)
        == reciprocal_mass(labeled, deduplicate=True),
        "support normalization did not remove the duplicate one",
    )
    return len(rows) + 2, tuple(rows)


def finite_abel_identity(support: list[int], cutoff: int) -> tuple[F, F]:
    values = sorted({value for value in support if 1 <= value <= cutoff})
    lhs = sum((F(1, value) for value in values), F(0))
    boundary = F(len(values), cutoff)
    integral = F(0)
    for index, value in enumerate(values, start=1):
        next_value = values[index] if index < len(values) else cutoff
        if next_value > value:
            integral += index * (F(1, value) - F(1, next_value))
    return lhs, boundary + integral


def abel_and_block_audit() -> tuple[int, int]:
    supports = [
        list(range(1, 31)),
        [index * index for index in range(1, 31)],
        [comb(index, 2) for index in range(2, 32)],
        [2**index for index in range(12)],
        [1, 2, 3, 6],
        [3, 6, 10, 15, 21, 28, 36],
    ]
    abel_rows = 0
    block_rows = 0
    for support in supports:
        cutoff = max(support) + 17
        lhs, rhs = finite_abel_identity(support, cutoff)
        require(lhs == rhs, "finite Abel--Stieltjes identity failed")
        abel_rows += 1

        maximum = max(support)
        exponent = 0
        while 2**exponent <= maximum:
            block = [
                value
                for value in set(support)
                if 2**exponent <= value < 2 ** (exponent + 1)
            ]
            mass = sum((F(1, value) for value in block), F(0))
            count = len(block)
            require(
                F(count, 2 ** (exponent + 1)) <= mass <= F(count, 2**exponent),
                "dyadic occupancy sandwich failed",
            )
            block_rows += 1
            exponent += 1

    require(harmonic(3) < 2 < harmonic(4), "harmonic crossing of two changed")
    require(
        reciprocal_mass([1, 2, 3, 6], deduplicate=True) == 2,
        "finite mass-two support changed",
    )
    return abel_rows, block_rows


def dirichlet_valuation_audit() -> int:
    rows = 0
    sets = [
        {1, 2, 5, 13},
        {2, 3, 5, 8, 13},
        {1, 4, 9, 16},
        {3, 6, 10, 15},
    ]
    for left, right in combinations(sets, 2):
        for exponent in (1, 2, 3, 4):
            def profile(support: set[int]) -> F:
                return sum((F(1, value**exponent) for value in support), F(0))

            require(
                profile(left | right) + profile(left & right)
                == profile(left) + profile(right),
                "Dirichlet support valuation failed",
            )
            rows += 1
    return rows


def polygonal(side_count: int, index: int) -> int:
    return ((side_count - 2) * index * index - (side_count - 4) * index) // 2


def polygonal_count(side_count: int, height: int) -> int:
    a = side_count - 2
    b = side_count - 4
    discriminant = b * b + 8 * a * height
    estimate = (b + isqrt(discriminant)) // (2 * a)
    while polygonal(side_count, estimate + 1) <= height:
        estimate += 1
    while estimate > 0 and polygonal(side_count, estimate) > height:
        estimate -= 1
    return estimate


def polygonal_arithmetic_audit() -> tuple[int, int, tuple[str, ...]]:
    partial_fraction_rows = 0
    count_rows = 0
    for side_count in range(3, 41):
        for index in range(1, 161):
            value = polygonal(side_count, index)
            if side_count != 4:
                shifted = F(index - 1) + F(2, side_count - 2)
                rhs = F(2, side_count - 4) * (F(1, shifted) - F(1, index))
                require(F(1, value) == rhs, "polygonal harmonic-clock split failed")
                partial_fraction_rows += 1
        for height in range(1, 2000, 17):
            count = polygonal_count(side_count, height)
            require(
                polygonal(side_count, count) <= height
                and polygonal(side_count, count + 1) > height,
                "polygonal counting formula failed",
            )
            count_rows += 1

    for side_count in range(3, 21):
        for cutoff in (1, 2, 7, 31):
            direct = sum(
                (F(1, polygonal(side_count, index)) for index in range(1, cutoff + 1)),
                F(0),
            )
            if side_count == 4:
                require(
                    direct == sum((F(1, index * index) for index in range(1, cutoff + 1)), F(0)),
                    "square partial sum changed",
                )
            else:
                clock = F(2, side_count - 4) * (
                    sum(
                        (
                            F(1, F(index - 1) + F(2, side_count - 2))
                            for index in range(1, cutoff + 1)
                        ),
                        F(0),
                    )
                    - harmonic(cutoff)
                )
                require(direct == clock, "finite polygonal clock identity failed")

    def mass(side_count: int) -> mp.mpf:
        if side_count == 4:
            return mp.pi**2 / 6
        return mp.mpf(2) / (side_count - 4) * (
            -mp.euler - mp.digamma(mp.mpf(2) / (side_count - 2))
        )

    closed_rows = {
        3: mp.mpf(2),
        4: mp.pi**2 / 6,
        5: 3 * mp.log(3) - mp.pi / mp.sqrt(3),
        6: 2 * mp.log(2),
        8: mp.mpf(3) * mp.log(3) / 4 + mp.sqrt(3) * mp.pi / 12,
        10: mp.log(2) + mp.pi / 6,
        14: 2 * mp.log(2) / 5 + 3 * mp.log(3) / 10 + mp.sqrt(3) * mp.pi / 10,
    }
    summaries: list[str] = []
    for side_count, closed in closed_rows.items():
        require(abs(mass(side_count) - closed) < mp.mpf("1e-52"), "paper closed form failed")
        summaries.append(f"S_{side_count}={mp.nstr(closed, 24)}")

    left = mass(4 - 10 ** -8)
    right = mass(4 + 10 ** -8)
    require(
        abs(left - mp.pi**2 / 6) < mp.mpf("1e-7")
        and abs(right - mp.pi**2 / 6) < mp.mpf("1e-7"),
        "square confluent limit failed",
    )
    return partial_fraction_rows, count_rows, tuple(summaries)


def master_figurate(side_count: int, dimension: int, index: int) -> int:
    return (side_count - 2) * comb(index + dimension - 2, dimension) + comb(
        index + dimension - 2, dimension - 1
    )


def master_figurate_audit() -> tuple[int, int]:
    factor_rows = 0
    simplex_rows = 0
    for side_count in range(3, 13):
        t = side_count - 2
        for dimension in range(2, 9):
            for index in range(1, 101):
                value = master_figurate(side_count, dimension, index)
                beta = F(
                    factorial(index - 1) * factorial(dimension - 2),
                    factorial(index + dimension - 2),
                )
                factor = F(
                    dimension * (dimension - 1),
                    dimension + t * (index - 1),
                ) * beta
                require(F(1, value) == factor, "master beta factorization failed")
                require(
                    value
                    >= master_figurate(3, dimension, index)
                    == comb(index + dimension - 1, dimension),
                    "master simplex domination failed",
                )
                factor_rows += 1

            cutoff = 80
            simplex_sum = sum(
                (F(1, master_figurate(3, dimension, index)) for index in range(1, cutoff + 1)),
                F(0),
            )
            exact = F(dimension, dimension - 1) * (
                1 - F(1, comb(cutoff + dimension - 1, dimension - 1))
            )
            require(simplex_sum == exact, "simplex exact tail identity failed")
            simplex_rows += 1
    return factor_rows, simplex_rows


def ceil_fraction(value: F) -> int:
    return -(-value.numerator // value.denominator)


def abel_dini_integer_lift_audit() -> tuple[int, tuple[str, ...]]:
    running = F(0)
    previous = {1: 0, 2: 0}
    masses = {1: F(0), 2: F(0)}
    rows = 0
    for index in range(1, 121):
        running += F(1, index)
        for power in (1, 2):
            lifted = ceil_fraction(F(index) * running**power)
            require(lifted > previous[power], "integer Abel--Dini lift is not increasing")
            previous[power] = lifted
            masses[power] += F(1, lifted)
            rows += 1

    block_summaries: list[str] = []
    for power in (1, 2):
        energies = []
        for cutoff in (16, 64, 256):
            energy = sum(
                (F(floor(F(2**level, level**power)), 2**level) for level in range(2, cutoff + 1)),
                F(0),
            )
            energies.append(float(energy))
        require(energies[0] < energies[1] < energies[2], "block energy lost monotonicity")
        block_summaries.append(
            f"p={power}:lift_mass_120={float(masses[power]):.12f};"
            f"block_energy={','.join(f'{value:.9f}' for value in energies)}"
        )
    return rows, tuple(block_summaries)


def fibonacci_table(maximum_index: int) -> list[int]:
    values = [0, 1]
    while len(values) <= maximum_index:
        values.append(values[-1] + values[-2])
    return values


def is_fibbinary(value: int) -> bool:
    return value >= 0 and value & (value >> 1) == 0


def is_moser_de_bruijn(value: int) -> bool:
    if value < 0:
        return False
    while value:
        value, digit = divmod(value, 4)
        if digit > 1:
            return False
    return True


def moser_value(mask: int) -> int:
    value = 0
    place = 1
    while mask:
        if mask & 1:
            value += place
        mask >>= 1
        place *= 4
    return value


def automatic_support_audit() -> tuple[int, int, tuple[str, ...]]:
    """Check the fibbinary/Moser Mahler decompositions and exact block laws."""

    word_cutoff = 1 << 18
    for value in range(word_cutoff):
        fibbinary_split = (
            (value % 2 == 0 and is_fibbinary(value // 2))
            or (value % 4 == 1 and is_fibbinary(value // 4))
        )
        moser_split = value % 4 in (0, 1) and is_moser_de_bruijn(value // 4)
        require(is_fibbinary(value) == fibbinary_split, "fibbinary Mahler split failed")
        require(is_moser_de_bruijn(value) == moser_split, "Moser Mahler split failed")
        require(
            not is_moser_de_bruijn(value) or is_fibbinary(value),
            "Moser support escaped fibbinary support",
        )

    fibonacci = fibonacci_table(26)
    block_rows = 0
    fibbinary_values = [value for value in range(1 << 19) if is_fibbinary(value)]
    for exponent in range(19):
        cumulative = sum(value < 2**exponent for value in fibbinary_values)
        block = sum(2**exponent <= value < 2 ** (exponent + 1) for value in fibbinary_values)
        require(cumulative == fibonacci[exponent + 2], "fibbinary cumulative count failed")
        require(block == fibonacci[exponent + 1], "fibbinary dyadic block count failed")
        block_rows += 1

    for exponent in range(16):
        cumulative_values = [moser_value(mask) for mask in range(1 << exponent)]
        block_values = [
            moser_value(mask) for mask in range(1 << exponent, 1 << (exponent + 1))
        ]
        require(len(set(cumulative_values)) == 2**exponent, "Moser cumulative count failed")
        require(
            all(0 <= value < 4**exponent for value in cumulative_values),
            "Moser cumulative block escaped",
        )
        require(len(set(block_values)) == 2**exponent, "Moser block count failed")
        require(
            min(block_values) == 4**exponent
            and max(block_values) == (4 ** (exponent + 1) - 1) // 3,
            "Moser block endpoints failed",
        )
        block_rows += 1

    tail_summaries: list[str] = []
    for exponent in (8, 12, 18):
        fib_lower = F(fibonacci[exponent + 3], 2**exponent)
        fib_upper = 2 * fib_lower
        finite_upper = sum(
            (F(fibonacci[level + 1], 2**level) for level in range(exponent, 23)),
            F(0),
        )
        remaining_upper = F(fibonacci[26], 2**22)
        require(
            finite_upper + remaining_upper == fib_upper,
            "fibbinary tail generating-function identity failed",
        )
        tail_summaries.append(
            f"fibbinary_K={exponent}:tail_between={fib_lower}..{fib_upper}"
        )
    for exponent in (8, 12, 16):
        moser_lower = F(3, 2 ** (exponent + 1))
        moser_upper = F(1, 2 ** (exponent - 1))
        finite_upper = sum(
            (F(1, 2**level) for level in range(exponent, 23)),
            F(0),
        )
        remaining_upper = F(1, 2**22)
        require(
            finite_upper + remaining_upper == moser_upper,
            "Moser tail geometric identity failed",
        )
        tail_summaries.append(
            f"moser_R={exponent}:tail_between={moser_lower}..{moser_upper}"
        )
    return word_cutoff, block_rows, tuple(tail_summaries)


def egyptian_and_sylvester_audit() -> tuple[int, tuple[int, ...], F]:
    refinement_rows = 0
    for denominator in range(1, 500):
        require(
            F(1, denominator)
            == F(1, denominator + 1) + F(1, denominator * (denominator + 1)),
            "Egyptian refinement failed",
        )
        refinement_rows += 1

    base = {1} | {2**index for index in range(1, 10)}
    for mask in range(1 << 8):
        support = set(base)
        expected_size = len(base)
        for bit in range(8):
            if mask & (1 << bit):
                denominator = 2 ** (bit + 1)
                support.remove(denominator)
                support.add(denominator + 1)
                support.add(denominator * (denominator + 1))
                expected_size += 1
        require(len(support) == expected_size, "Egyptian refinement collided")
        require(
            reciprocal_mass(sorted(support), deduplicate=True)
            == reciprocal_mass(sorted(base), deduplicate=True),
            "Egyptian refinement changed finite mass",
        )

    terms = [2]
    for _ in range(5):
        terms.append(terms[-1] * terms[-1] - terms[-1] + 1)
    partial = F(0)
    for index, term in enumerate(terms[:-1]):
        partial += F(1, term)
        require(
            1 - partial == F(1, terms[index + 1] - 1),
            "Sylvester exact remainder failed",
        )
    return refinement_rows + (1 << 8), tuple(terms), partial


def tournament_sequence_audit() -> tuple[str, str, str]:
    odd_partial = mp.nsum(
        lambda index: mp.mpf(6) / (index * (index + 1) * (2 * index + 1)),
        [1, mp.inf],
    )
    even_partial = mp.nsum(
        lambda index: mp.mpf(3) / (index * (index - 1) * (index + 1)),
        [2, mp.inf],
    )
    odd_closed = 18 - 24 * mp.log(2)
    even_closed = mp.mpf(3) / 4
    require(abs(odd_partial - odd_closed) < mp.mpf("1e-52"), "odd max-c3 mass failed")
    require(abs(even_partial - even_closed) < mp.mpf("1e-52"), "even max-c3 mass failed")

    q = mp.mpf(1) / 2
    theta_sum = mp.nsum(lambda index: q ** (index * (index + 1) / 2), [0, mp.inf])
    theta_product = mp.nprod(
        lambda index: (1 - q ** (2 * index)) / (1 - q ** (2 * index - 1)),
        [1, mp.inf],
    )
    require(abs(theta_sum - theta_product) < mp.mpf("1e-52"), "Ramanujan psi product failed")
    return (
        f"max_c3_odd={mp.nstr(odd_closed, 25)}",
        f"max_c3_all={mp.nstr(odd_closed + even_closed, 25)}",
        f"tournament_theta={mp.nstr(theta_sum, 25)}",
    )


def source_has_no_assert_nodes() -> int:
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    count = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    require(count == 0, "optimization-sensitive assert node found")
    return count


def main() -> None:
    no_asserts = source_has_no_assert_nodes()
    collision_rows, collision_summary = support_collision_audit()
    abel_rows, block_rows = abel_and_block_audit()
    valuation_rows = dirichlet_valuation_audit()
    polygon_rows, polygon_count_rows, polygon_summary = polygonal_arithmetic_audit()
    master_rows, simplex_rows = master_figurate_audit()
    dini_rows, dini_summary = abel_dini_integer_lift_audit()
    automatic_rows, automatic_block_rows, automatic_summary = automatic_support_audit()
    egyptian_rows, sylvester_terms, sylvester_partial = egyptian_and_sylvester_audit()
    tournament_summary = tournament_sequence_audit()

    print("THM-2000 SUPPORT-HARMONIC / ABEL-DINI / FIGURATE EXACT AUDIT")
    print(f"Python assert nodes = {no_asserts}")
    print(f"support/collision rows = {collision_rows}: {collision_summary}")
    print(f"finite Abel rows / dyadic block rows = {abel_rows}/{block_rows}")
    print(f"Dirichlet valuation rows = {valuation_rows}")
    print(f"polygonal split/count rows = {polygon_rows}/{polygon_count_rows}")
    print(f"polygonal closed rows = {polygon_summary}")
    print(f"master beta rows / simplex-tail rows = {master_rows}/{simplex_rows}")
    print(f"integer Abel-Dini lift rows = {dini_rows}")
    print(f"Abel-Dini block summaries = {dini_summary}")
    print(f"automatic language rows / block rows = {automatic_rows}/{automatic_block_rows}")
    print(f"automatic exact tail summaries = {automatic_summary}")
    print(f"Egyptian refinement rows = {egyptian_rows}")
    print(f"Sylvester terms = {sylvester_terms}")
    print(f"Sylvester checked partial/remainder = {sylvester_partial}/{1-sylvester_partial}")
    print(f"tournament exact rows = {tournament_summary}")
    print("support_profile=D_A(s)=sum_(m in A)m^(-s); mass=D_A(1)")
    print("Abel=sum_(m<=X,m in A)1/m=A(X)/X+integral_1^X A(t)/t^2 dt")
    print("dyadic=converges iff sum_k #A[2^k,2^(k+1))/2^k converges")
    print("Bertrand=n*log(n)*...*log_r(n)^p converges iff p>1")
    print("master_axes=s=3 gives simplices; d=2 gives polygons")
    print("automatic_dimensions=fibbinary log_2(phi); Moser-de-Bruijn 1/2")
    print("scalar_guardrail=triangular_support and {1,2,3,6} both have mass 2")
    print("STATUS=PASS")
    print(f"source_sha256={hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}")


if __name__ == "__main__":
    main()
