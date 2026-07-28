#!/usr/bin/env python3
"""Exact companion for THM-2824.

The script verifies the binary quadratic/cubic reduction for normalized
factorial monomials f_n=s^n/n!, the strict cubic orientation t_111>0, the
Wronskian formula, and the remaining atomic inequality through c<=50.

The finite atomic scan is evidence, not a proof of the open universal
inequality.  The script uses only integers and fractions, has explicit
exception gates, and has no scratch dependency.
"""

from fractions import Fraction
from functools import lru_cache
from math import factorial


MAX_C = 50


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


@lru_cache(maxsize=None)
def moment(indices):
    total = sum(indices)
    denominator = 1
    for index in indices:
        denominator *= factorial(index)
    return factorial(total) // denominator


def bilinear(left, right):
    return sum(
        a_value * b_value * moment(tuple(sorted((a_index, b_index))))
        for a_index, a_value in left.items()
        for b_index, b_value in right.items()
    )


def trilinear(left, middle, right):
    return sum(
        a_value
        * b_value
        * c_value
        * moment(tuple(sorted((a_index, b_index, c_index))))
        for a_index, a_value in left.items()
        for b_index, b_value in middle.items()
        for c_index, c_value in right.items()
    )


def add_sparse(*polynomials):
    out = {}
    for polynomial in polynomials:
        for index, value in polynomial.items():
            out[index] = out.get(index, 0) + value
            if out[index] == 0:
                del out[index]
    return out


def scale_sparse(polynomial, scalar):
    return {
        index: scalar * value
        for index, value in polynomial.items()
        if scalar * value
    }


def multiply_sparse(left, right):
    out = {}
    for i, a in left.items():
        for j, b in right.items():
            out[i + j] = out.get(i + j, Fraction(0)) + a * b
    return {index: value for index, value in out.items() if value}


def derivative_sparse(polynomial):
    return {
        index - 1: index * value
        for index, value in polynomial.items()
        if index
    }


def normalized_monomial(index):
    return {index: Fraction(1, factorial(index))}


def factorial_ratio(high, low):
    require(high >= low, "factorial ratio order drift")
    return factorial(high) // factorial(low)


def invariants(a, b, c):
    f_a = {a: 1}
    f_b = {b: 1}
    f_c = {c: 1}
    u = add_sparse(f_b, scale_sparse(f_a, -1))
    v = add_sparse(f_c, scale_sparse(f_b, -1))

    g11 = bilinear(u, u)
    g12 = bilinear(u, v)
    g22 = bilinear(v, v)
    t111 = trilinear(u, u, u)
    t112 = trilinear(u, u, v)
    t122 = trilinear(u, v, v)
    t222 = trilinear(v, v, v)

    i1 = (
        3 * t112 * g11 * g22
        - t222 * g11 * g11
        - 2 * t111 * g12 * g22
    )
    i2 = (
        3 * t122 * g11 * g22
        - 2 * t222 * g12 * g11
        - t111 * g22 * g22
    )
    d_total = 2 * t222 * g12 - 3 * t122 * g22
    return (
        u,
        v,
        g11,
        g12,
        g22,
        t111,
        t112,
        t122,
        t222,
        i1,
        i2,
        d_total,
    )


def atomic_d(i, b, c):
    d_i = {i + 1: 1, i: -1}
    v = {c: 1, b: -1}
    g22 = bilinear(v, v)
    t222 = trilinear(v, v, v)
    return (
        2 * t222 * bilinear(d_i, v)
        - 3 * trilinear(d_i, v, v) * g22
    )


def verify_wronskian(a, b, c):
    p = b - a
    r = c - b
    a_scale = factorial_ratio(b, a)
    b_scale = factorial_ratio(c, b)

    f_a = normalized_monomial(a)
    x = {p: Fraction(1, a_scale)}
    y = {r: Fraction(1, b_scale)}
    one = {0: Fraction(1)}
    u = multiply_sparse(f_a, add_sparse(x, scale_sparse(one, -1)))
    v = multiply_sparse(
        multiply_sparse(f_a, x),
        add_sparse(y, scale_sparse(one, -1)),
    )
    actual = add_sparse(
        multiply_sparse(u, derivative_sparse(v)),
        scale_sparse(multiply_sparse(derivative_sparse(u), v), -1),
    )

    prefactor = {
        2 * a + p - 1: Fraction(1, factorial(a) ** 2 * a_scale)
    }
    xy = multiply_sparse(x, y)
    bracket = add_sparse(
        {0: Fraction(p)},
        scale_sparse(xy, r),
        scale_sparse(y, -(p + r)),
    )
    expected = multiply_sparse(prefactor, bracket)
    require(actual == expected, f"Wronskian identity drift at {(a, b, c)}")

    # B^(1/r)>A^(1/p), written without radicals.
    require(
        b_scale**p > a_scale**r,
        f"root-mean orientation drift at {(a, b, c)}",
    )


def verify_t111_ratio(a, b, t111):
    r = b - a
    f_values = []
    for k in range(4):
        numerator = factorial(3 * a + k * r)
        denominator = factorial(a) ** 3 * factorial_ratio(b, a) ** k
        f_values.append(Fraction(numerator, denominator))
    reconstructed = f_values[3] - 3 * f_values[2] + 3 * f_values[1] - f_values[0]
    require(reconstructed == t111, f"t111 reconstruction drift at {(a, b)}")
    ratios = tuple(f_values[k + 1] / f_values[k] for k in range(3))
    products = tuple(
        Fraction(
            factorial(3 * a + (k + 1) * r),
            factorial(3 * a + k * r) * factorial_ratio(b, a),
        )
        for k in range(3)
    )
    require(ratios == products, f"ratio product drift at {(a, b)}")
    if a == 0:
        require(ratios[0] == 1, f"a=0 rho0 drift at {(a, b)}")
        require(ratios[1] >= 2 and ratios[2] >= 3, "a=0 ratio floor drift")
    else:
        require(
            ratios[0] > 1 and ratios[1] > 2**r and ratios[2] >= 3**r,
            f"a>0 ratio floor drift at {(a, b)}",
        )
    require(t111 > 0, f"strict cubic orientation failed at {(a, b)}")


def main():
    support_count = 0
    atomic_count = 0
    atomic_positive = 0
    atomic_zero = []
    d_total_zero = []
    i1_negative = 0
    i2_negative = 0
    divisibility_identity_checks = 0
    telescoping_checks = 0
    wronskian_checks = 0
    ratio_checks = 0

    smallest_positive_atomic = None
    smallest_positive_atomic_index = None

    for c in range(2, MAX_C + 1):
        for b in range(1, c):
            for i in range(b):
                value = atomic_d(i, b, c)
                require(value >= 0, f"atomic orientation failed at {(i, b, c)}")
                atomic_count += 1
                if value == 0:
                    atomic_zero.append((i, b, c))
                else:
                    atomic_positive += 1
                    if smallest_positive_atomic is None or value < smallest_positive_atomic:
                        smallest_positive_atomic = value
                        smallest_positive_atomic_index = (i, b, c)

            for a in range(b):
                (
                    _u,
                    _v,
                    g11,
                    g12,
                    g22,
                    t111,
                    t112,
                    t122,
                    t222,
                    i1,
                    i2,
                    d_total,
                ) = invariants(a, b, c)
                require(
                    g11 > 0 and g22 > 0 and g11 * g22 - g12 * g12 > 0,
                    f"Gram positivity drift at {(a, b, c)}",
                )
                verify_t111_ratio(a, b, t111)
                ratio_checks += 1
                verify_wronskian(a, b, c)
                wronskian_checks += 1

                require(
                    i2 == -g11 * d_total - t111 * g22 * g22,
                    f"I2/D identity drift at {(a, b, c)}",
                )
                divisibility_identity_checks += 1

                atomic_sum = sum(atomic_d(i, b, c) for i in range(a, b))
                require(
                    d_total == atomic_sum,
                    f"atomic telescoping drift at {(a, b, c)}",
                )
                telescoping_checks += 1
                require(d_total >= 0, f"finite D total failed at {(a, b, c)}")
                require(i2 < 0, f"finite common-null obstruction failed at {(a, b, c)}")
                if i1 < 0:
                    i1_negative += 1
                if i2 < 0:
                    i2_negative += 1
                if d_total == 0:
                    d_total_zero.append((a, b, c))
                support_count += 1

    expected_count = (MAX_C + 1) * MAX_C * (MAX_C - 1) // 6
    expected_zero = tuple((i, i + 1, i + 2) for i in range(MAX_C - 1))
    require(support_count == expected_count, "support triple count drift")
    require(atomic_count == expected_count, "atomic triple count drift")
    require(tuple(atomic_zero) == expected_zero, "atomic equality boundary drift")
    require(tuple(d_total_zero) == expected_zero, "total equality boundary drift")
    require(i1_negative == support_count, "finite I1 sign drift")
    require(i2_negative == support_count, "finite I2 sign drift")

    print("THM-2824 ARBITRARY THREE-SLOT FACTORIAL ATOMIC ORIENTATION")
    print(f"scan_max_c={MAX_C}; support_triples={support_count}; atomic_triples={atomic_count}")
    print(
        "gram_positive=all; t111_positive=all; "
        f"ratio_proof_checks={ratio_checks}; wronskian_checks={wronskian_checks}"
    )
    print(
        f"divisibility_identity_checks={divisibility_identity_checks}; "
        f"atomic_telescoping_checks={telescoping_checks}"
    )
    print(
        f"atomic_nonnegative={atomic_count}/{atomic_count}; "
        f"atomic_positive={atomic_positive}; atomic_zero={len(atomic_zero)}"
    )
    print(
        "atomic_zero_boundary=(i,b,c)=(j,j+1,j+2); "
        f"j=0..{MAX_C - 2}"
    )
    print(
        f"smallest_positive_atomic={smallest_positive_atomic}; "
        f"index={smallest_positive_atomic_index}"
    )
    print(
        f"I1_negative={i1_negative}/{support_count}; "
        f"I2_negative={i2_negative}/{support_count}; "
        f"D_total_zero={len(d_total_zero)}"
    )
    print(
        "scope=finite exact evidence only for universal atomic D_i>=0; "
        "the quadratic/cubic divisibility reduction, t111 orientation, "
        "Wronskian formula, and conditional implication are exact"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
