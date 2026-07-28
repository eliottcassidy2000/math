#!/usr/bin/env python3
"""Exact companion for THM-2828.

Audits the adjacent-difference cubic tensor and the lower-prefix cone
extension of THM-2824 using exact integers and fractions only.
"""

from fractions import Fraction
from functools import lru_cache
from math import comb, factorial


MAX_TENSOR_INDEX = 24
MAX_B = 16
MAX_C = 24


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


@lru_cache(maxsize=None)
def moment(indices):
    return factorial(sum(indices)) // product(factorial(i) for i in indices)


def product(values):
    out = 1
    for value in values:
        out *= value
    return out


def add_sparse(*polynomials):
    out = {}
    for polynomial in polynomials:
        for index, value in polynomial.items():
            out[index] = out.get(index, Fraction(0)) + value
            if out[index] == 0:
                del out[index]
    return out


def scale_sparse(polynomial, scalar):
    return {
        index: scalar * value
        for index, value in polynomial.items()
        if scalar * value
    }


def d(index):
    return {index + 1: Fraction(1), index: Fraction(-1)}


def bilinear(left, right):
    return sum(
        a * b * moment(tuple(sorted((i, j))))
        for i, a in left.items()
        for j, b in right.items()
    )


def trilinear(left, middle, right):
    return sum(
        a * b * c * moment(tuple(sorted((i, j, k))))
        for i, a in left.items()
        for j, b in middle.items()
        for k, c in right.items()
    )


def cubic_closed(i, j, k):
    total = i + j + k
    e2 = i * j + i * k + j * k
    e3 = i * j * k
    multinomial = Fraction(
        factorial(total),
        factorial(i) * factorial(j) * factorial(k),
    )
    numerator = 2 * (total + 1) ** 2 + total * e2 - e3
    denominator = (i + 1) * (j + 1) * (k + 1)
    return multinomial * Fraction(numerator, denominator)


def atomic_d(i, b, c):
    v = {c: Fraction(1), b: Fraction(-1)}
    return (
        2 * trilinear(v, v, v) * bilinear(d(i), v)
        - 3 * trilinear(d(i), v, v) * bilinear(v, v)
    )


def cone_fixtures(b):
    fixtures = [
        tuple(Fraction(1 if i == 0 else 0) for i in range(b)),
        tuple(Fraction(1 if i == b - 1 else 0) for i in range(b)),
        tuple(Fraction(1) for _ in range(b)),
        tuple(Fraction(i + 1) for i in range(b)),
        tuple(Fraction((3 * i + b) % 5) for i in range(b)),
        tuple(Fraction((i + 1) % 3, i + 2) for i in range(b)),
    ]
    return tuple(values for values in fixtures if any(values))


def verify_prefix_coefficients(lambdas, u):
    b = len(lambdas)
    coefficients = tuple(u.get(i, Fraction(0)) for i in range(b + 1))
    prefix = Fraction(0)
    for i in range(b):
        prefix += coefficients[i]
        require(prefix == -lambdas[i], f"prefix coordinate drift at {(b, i)}")
        require(prefix <= 0, f"prefix cone sign drift at {(b, i)}")
    prefix += coefficients[b]
    require(prefix == 0, f"total coefficient drift at b={b}")


def main():
    tensor_checks = 0
    quadratic_checks = 0
    cone_checks = 0
    atom_checks = 0
    gaussian_eligible = 0
    smallest_cubic = None
    smallest_cubic_index = None

    for i in range(MAX_TENSOR_INDEX + 1):
        for j in range(MAX_TENSOR_INDEX + 1):
            for k in range(MAX_TENSOR_INDEX + 1):
                direct = trilinear(d(i), d(j), d(k))
                closed = cubic_closed(i, j, k)
                require(direct == closed, f"cubic tensor formula drift at {(i, j, k)}")
                require(closed > 0, f"cubic tensor positivity drift at {(i, j, k)}")
                if smallest_cubic is None or closed < smallest_cubic:
                    smallest_cubic = closed
                    smallest_cubic_index = (i, j, k)
                tensor_checks += 1

    for i in range(MAX_TENSOR_INDEX + 1):
        for j in range(MAX_TENSOR_INDEX + 1):
            value = bilinear(d(i), d(j))
            require(value == comb(i + j, i), f"Pascal fixed-point drift at {(i, j)}")
            quadratic_checks += 1

    for b in range(1, MAX_B + 1):
        for c in range(b + 1, MAX_C + 1):
            v = {c: Fraction(1), b: Fraction(-1)}
            atoms = tuple(atomic_d(i, b, c) for i in range(b))
            for i, value in enumerate(atoms):
                require(value >= 0, f"inherited atomic sign drift at {(i, b, c)}")
                require(
                    (value == 0) == (i == b - 1 and c == b + 1),
                    f"inherited atomic equality drift at {(i, b, c)}",
                )
                atom_checks += 1

            for lambdas in cone_fixtures(b):
                u = {}
                for i, value in enumerate(lambdas):
                    u = add_sparse(u, scale_sparse(d(i), value))
                verify_prefix_coefficients(lambdas, u)
                require(sum(u.values()) == 0, f"mean-zero coefficient drift at {(b, c)}")

                t111 = trilinear(u, u, u)
                tensor_sum = sum(
                    lambdas[i] * lambdas[j] * lambdas[k] * cubic_closed(i, j, k)
                    for i in range(b)
                    for j in range(b)
                    for k in range(b)
                )
                require(t111 == tensor_sum and t111 > 0, f"cone cubic drift at {(b, c)}")

                g11 = bilinear(u, u)
                g12 = bilinear(u, v)
                g22 = bilinear(v, v)
                t122 = trilinear(u, v, v)
                t222 = trilinear(v, v, v)
                d_total = 2 * t222 * g12 - 3 * t122 * g22
                atomic_sum = sum(lambdas[i] * atoms[i] for i in range(b))
                require(d_total == atomic_sum, f"cone atomic linearity drift at {(b, c)}")
                require(d_total >= 0, f"cone orientation sign drift at {(b, c)}")
                require(
                    g11 > 0 and g22 > 0 and g11 * g22 - g12 * g12 > 0,
                    f"cone Gram positivity drift at {(b, c)}",
                )

                i2 = (
                    3 * t122 * g11 * g22
                    - 2 * t222 * g12 * g11
                    - t111 * g22 * g22
                )
                require(
                    i2 == -g11 * d_total - t111 * g22 * g22,
                    f"cone I2 identity drift at {(b, c)}",
                )
                require(i2 < 0, f"cone divisibility obstruction drift at {(b, c)}")
                if lambdas[0] == 0:
                    gaussian_eligible += 1
                cone_checks += 1

    print("THM-2828 LOWER-PREFIX FACTORIAL CONE")
    print(
        f"tensor_max_index={MAX_TENSOR_INDEX}; cubic_tensor_checks={tensor_checks}; "
        f"pascal_fixed_point_checks={quadratic_checks}"
    )
    print(
        f"smallest_cubic_tensor={smallest_cubic}; "
        f"index={smallest_cubic_index}; all_strict_positive=True"
    )
    print(
        f"cone_bank_b<={MAX_B}; c<={MAX_C}; cone_directions={cone_checks}; "
        f"inherited_atom_checks={atom_checks}"
    )
    print(
        f"prefix_equivalence=all; cubic_orientation=all_strict; "
        f"mixed_orientation=all_nonnegative; I2=all_strict_negative"
    )
    print(
        f"gaussian_eligible_lambda0_zero={gaussian_eligible}; "
        "scope=many lower slots in one prefix cone, not arbitrary radial coefficients"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
