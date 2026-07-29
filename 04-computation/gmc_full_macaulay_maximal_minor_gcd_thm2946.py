#!/usr/bin/env python3
"""Exact controls for THM-2946.

The theorem itself is the universal Fitting/resultant argument.  This
companion checks its Hilbert/Koszul invoice and an independent six-point
specialization with one transverse common root.
"""

from hashlib import sha256
from itertools import product
from math import comb

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def monomials_of_degree(variables, degree):
    x, y, z = variables
    return [
        x**i * y**j * z ** (degree - i - j)
        for i in range(degree, -1, -1)
        for j in range(degree - i, -1, -1)
    ]


def hilbert_complete_intersection(maximum_degree):
    degrees = (2, 3, 4)
    values = []
    for degree in range(maximum_degree + 1):
        value = 0
        for mask in range(1 << len(degrees)):
            shift = sum(
                degrees[index]
                for index in range(len(degrees))
                if mask & (1 << index)
            )
            sign = -1 if mask.bit_count() % 2 else 1
            value += sign * (
                comb(degree - shift + 2, 2) if degree >= shift else 0
            )
        values.append(value)
    return values


def coefficient_vector(polynomial, target_basis, variables):
    expanded = sp.Poly(sp.expand(polynomial), *variables)
    return [
        expanded.coeff_monomial(monomial)
        for monomial in target_basis
    ]


def macaulay_matrix(Q, C, F, variables):
    target = monomials_of_degree(variables, 7)
    columns = []
    for form, multiplier_degree in ((Q, 5), (C, 4), (F, 3)):
        for multiplier in monomials_of_degree(variables, multiplier_degree):
            columns.append(
                coefficient_vector(form * multiplier, target, variables)
            )
    return sp.Matrix.hstack(*(sp.Matrix(column) for column in columns))


def lf_digest(text):
    return sha256(text.replace("\r\n", "\n").encode()).hexdigest()


def main():
    print("THM-2946 full Macaulay maximal-minor gcd controls")

    dimensions = tuple(comb(degree + 2, 2) for degree in (5, 4, 3, 7))
    require(dimensions == (21, 15, 10, 36), "graded dimensions changed")
    hilbert = hilbert_complete_intersection(9)
    require(
        hilbert == [1, 3, 5, 6, 5, 3, 1, 0, 0, 0],
        "complete-intersection Hilbert series changed",
    )
    koszul_kernel = comb(4, 2) + comb(3, 2) + comb(2, 2)
    require(koszul_kernel == 10, "degree-seven Koszul kernel changed")
    require(sum(dimensions[:3]) - dimensions[3] == 10, "rank invoice changed")

    x, y, z, t = sp.symbols("x y z t")
    variables = (x, y, z)
    Q = x**2 - z**2
    C = y**3 - z**3
    F = z**3 * (x + 2 * y - 3 * z) + t * z**4

    matrix_zero = macaulay_matrix(Q, C, F.subs(t, 0), variables)
    matrix_one = macaulay_matrix(Q, C, F.subs(t, 1), variables)
    require(matrix_zero.shape == (36, 46), "Macaulay shape changed")
    require(matrix_zero.rank() == 35, "simple-root rank is not 35")
    require(matrix_one.rank() == 36, "off-resultant rank is not 36")

    omega = sp.symbols("omega")
    first_norm = sp.resultant(
        omega**3 - 1,
        t - 2 + 2 * omega,
        omega,
    )
    second_norm = sp.resultant(
        omega**3 - 1,
        t - 4 + 2 * omega,
        omega,
    )
    require(
        sp.expand(first_norm) == sp.expand((t - 2) ** 3 + 8),
        "first three-point norm changed",
    )
    require(
        sp.expand(second_norm) == sp.expand((t - 4) ** 3 + 8),
        "second three-point norm changed",
    )
    six_point_norm = sp.factor(first_norm * second_norm)
    expected_norm = (
        t
        * (t**2 - 6 * t + 12)
        * (t**3 - 12 * t**2 + 48 * t - 56)
    )
    require(
        sp.expand(six_point_norm - expected_norm) == 0,
        "six-point norm factorization changed",
    )
    require(
        expected_norm.subs(t, 0) == 0
        and sp.diff(expected_norm, t).subs(t, 0) == -672,
        "simple resultant valuation changed",
    )

    # At t=0, the chosen point contributes zero and the product of the
    # other five values is the derivative of the six-point product.
    require(
        sp.diff(expected_norm, t).subs(t, 0) == -672,
        "five-point nonzero product changed",
    )

    u, v = sp.symbols("u v")
    require(sp.gcd(sp.Poly(u, u, v), sp.Poly(v, u, v)).total_degree() == 0,
            "universal specialization hostile changed")
    require(sp.gcd(sp.Poly(u, u), sp.Poly(u, u)).as_expr() == u,
            "specialized gcd hostile changed")

    record = (
        f"hilbert={','.join(str(value) for value in hilbert)};"
        f"ranks={matrix_zero.rank()},{matrix_one.rank()};"
        f"norm={sp.sstr(expected_norm)}"
    )
    print("macaulay_shape=36x46;domain=21+15+10")
    print("hilbert_series=1,3,5,6,5,3,1,0,0,0")
    print("socle_degree=6;degree7_kernel=10;degree7_rank=36")
    print("specialization_ranks=t0:35;t1:36")
    print(
        "six_point_norm="
        "t*(t^2-6*t+12)*(t^3-12*t^2+48*t-56)"
    )
    print("simple_root=t0;derivative=-672")
    print("specialization_can_gain_gcd=true")
    print("control_digest=" + lf_digest(record))
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
