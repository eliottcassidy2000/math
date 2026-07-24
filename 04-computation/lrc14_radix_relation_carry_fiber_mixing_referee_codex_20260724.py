#!/usr/bin/env python3
"""Exact referee for THM-2163 and the strengthened THM-2161 hostile pair."""

from collections import defaultdict
from itertools import product
from math import gcd, lcm


HEIGHT = 29
L = 360_360

A = (
    1,
    2,
    3,
    4,
    5,
    6,
    1_441_447,
    44_324_288,
    1_331_890_569,
    39_958_158_250,
    1_198_739_341_811,
    35_962_180_254_012,
    1_078_865_413_025_413,
)

B = (
    1,
    2,
    3,
    4,
    5,
    6,
    360_367,
    13_693_688,
    412_972_569,
    12_384_492_130,
    371_535_484_331,
    11_146_064_529_612,
    1_078_865_413_025_413,
)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def dot(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return sum(x * y for x, y in zip(left, right))


def centered_abs(residue: int, modulus: int) -> int:
    residue %= modulus
    return min(residue, modulus - residue)


def modulus_margin(row: tuple[int, ...], modulus: int) -> int:
    return max(
        min(centered_abs(multiplier * value, modulus) for value in row)
        for multiplier in range(1, modulus)
    )


def relation_box_core_count() -> int:
    counts: dict[int, int] = {0: 1}
    for speed in range(1, 7):
        new_counts: defaultdict[int, int] = defaultdict(int)
        for subtotal, multiplicity in counts.items():
            for coefficient in range(-HEIGHT, HEIGHT + 1):
                new_counts[subtotal + coefficient * speed] += multiplicity
        counts = dict(new_counts)
    return counts[0]


def radix_digits(values: tuple[int, ...], base: int, length: int) -> list[tuple[int, ...]]:
    return [
        tuple((value // (base**j)) % base for value in values)
        for j in range(length)
    ]


def transition_path(
    coefficients: tuple[int, ...],
    digits: list[tuple[int, ...]],
    base: int,
) -> tuple[bool, list[int]]:
    carry = 0
    path = [carry]
    for digit in digits:
        numerator = carry + dot(coefficients, digit)
        if numerator % base:
            return False, path
        carry = numerator // base
        path.append(carry)
    return True, path


def reconstruct(digits: list[tuple[int, ...]], base: int) -> tuple[int, ...]:
    dimension = len(digits[0])
    return tuple(
        sum((base**j) * digits[j][i] for j in range(len(digits)))
        for i in range(dimension)
    )


def exhaustive_converse_referee() -> int:
    accepted = 0
    for base in (2, 3, 4):
        dimension, length = 2, 3
        coefficient_words = [
            word
            for word in product(range(-2, 3), repeat=dimension)
            if any(word)
        ]
        for flat_digits in product(range(base), repeat=dimension * length):
            digits = [
                tuple(flat_digits[j * dimension : (j + 1) * dimension])
                for j in range(length)
            ]
            values = reconstruct(digits, base)
            for coefficients in coefficient_words:
                valid, path = transition_path(coefficients, digits, base)
                if valid and path[-1] == 0:
                    require(dot(coefficients, values) == 0, "converse reconstruction")
                    accepted += 1
                if dot(coefficients, values) == 0:
                    require(valid and path[-1] == 0, "forward digit path")
                    l1 = sum(abs(value) for value in coefficients)
                    require(all(abs(carry) < l1 for carry in path), "strict carry bound")
    return accepted


def verify_actual_decomposition() -> None:
    examples = (
        ((7, 19, 31), (2, -3, 1)),   # 14-57+31=-12, not a relation
        ((5, 8, 11), (3, -6, 3)),    # 15-48+33=0
        ((1, 2, 3, 4), (2, -1, 0, 0)),
    )
    for values, coefficients in examples:
        if dot(coefficients, values) != 0:
            continue
        for base in (2, 3, 5):
            power = 1
            length = 0
            while power <= max(values):
                power *= base
                length += 1
            digits = radix_digits(values, base, length)
            valid, path = transition_path(coefficients, digits, base)
            require(valid and path[0] == path[-1] == 0, "actual carry endpoints")
            l1 = sum(abs(value) for value in coefficients)
            require(all(abs(carry) < l1 for carry in path), "actual strict carry bound")
            for j, carry in enumerate(path):
                scale = base**j
                quotient = tuple(value // scale for value in values)
                remainder = tuple(value % scale for value in values)
                require(carry == dot(coefficients, remainder) // scale, "remainder carry")
                require(carry == -dot(coefficients, quotient), "quotient carry")


def main() -> None:
    require(lcm(*range(2, 14)) == L, "small-modulus lcm")
    for name, row in (("A", A), ("B", B)):
        require(len(row) == 13 and len(set(row)) == 13, f"{name} distinctness")
        require(gcd(*row) == 1, f"{name} primitivity")
        require(row[:6] == tuple(range(1, 7)), f"{name} retained core")
        require(
            all(row[index] != index + 1 for index in range(6, 13)),
            f"{name} defect exactly seven",
        )
        require(
            all(value % L == index + 1 for index, value in enumerate(row)),
            f"{name} labelled CRT residues",
        )
        for index in range(6, 13):
            require(
                row[index] > HEIGHT * sum(row[:index]),
                f"{name} coefficient-superincreasing index {index + 1}",
            )

    require(max(A) == max(B), "same scalar maximum")
    require(A[-1] == B[-1] == 1_078_865_413_025_413, "common final owner")

    margins = []
    ap = tuple(range(1, 14))
    for modulus in range(2, 14):
        margin_a = modulus_margin(A, modulus)
        margin_b = modulus_margin(B, modulus)
        margin_ap = modulus_margin(ap, modulus)
        require(margin_a == margin_b == margin_ap == 0, f"modulus {modulus} bank")
        margins.append(margin_a)

    margin17_a = modulus_margin(A, 17)
    margin17_b = modulus_margin(B, 17)
    require(margin17_a == 0, "A mod-17 margin")
    require(margin17_b == 2, "B mod-17 margin")
    require(
        min(centered_abs(2 * value, 17) for value in B) == 2,
        "B mod-17 witness multiplier",
    )

    core_count = relation_box_core_count()
    require(core_count > 1, "nontrivial retained relation box")

    verify_actual_decomposition()
    accepted_toy_paths = exhaustive_converse_referee()
    require(accepted_toy_paths > 0, "toy carry census")

    print("THM-2163 exact radix relation-carry / CRT-fiber referee")
    print(f"L=lcm(2,...,13)={L}; height={HEIGHT}")
    print(f"A and B: primitive, distinct, defect=7, same max={max(A)}")
    print("labelled residues agree with (1,...,13) for every modulus 2..13")
    print(f"coefficient-superincreasing tail verified at all 7 tail indices")
    print(f"common full height-29 relation box has {core_count} vectors, all tail coefficients zero")
    print(f"small-bank margins q=2..13: {tuple(margins)}")
    print(f"adaptive split: m_17(A)={margin17_a}, m_17(B)={margin17_b}; B gives M>=2/17")
    print(f"exhaustive radix converse/forward toy paths accepted={accepted_toy_paths}")
    print("carry endpoints, divisibility, quotient identity, and strict l1 bound all verified")
    print("ALL EXACT ASSERTIONS PASSED")


if __name__ == "__main__":
    main()
