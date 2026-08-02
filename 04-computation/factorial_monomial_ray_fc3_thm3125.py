#!/usr/bin/env python3
"""Exact transfer controls for THM-3125's three-variable monomial rays."""

from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import factorial, prod


V_MAX = 3
D_MAX = 4
N_MAX = 7


def rising(shape: Fraction, n: int) -> Fraction:
    return prod(shape + j for j in range(n))


def factorial_ray_weight(v, d, n):
    return prod(factorial(d * vi * n) for vi in v)


def gauss_ray_weight(v, d, n):
    positive_q = [d * vi for vi in v if vi > 0]
    scale = prod(q**q for q in positive_q)
    answer = Fraction(scale**n, 1)
    for q in positive_q:
        for k in range(1, q + 1):
            answer *= rising(Fraction(k, q), n)
    return answer


def normalized_product_moment(weights, *vectors):
    """Moment of a product of sparse vectors in f_j=t^j/W_j coordinates."""
    states = {0: Fraction(1)}
    for vector in vectors:
        nxt = {}
        for old_degree, old_coef in states.items():
            for degree, coef in vector.items():
                new_degree = old_degree + degree
                nxt[new_degree] = nxt.get(new_degree, Fraction(0)) + old_coef * coef
        states = nxt
    return sum(coef * weights[degree] for degree, coef in states.items())


def orientation_invariants(v, d):
    weights = [Fraction(factorial_ray_weight(v, d, n), 1) for n in range(7)]
    u = {1: Fraction(1, weights[1]), 0: Fraction(-1)}
    z = {2: Fraction(1, weights[2]), 1: Fraction(-1, weights[1])}
    g11 = normalized_product_moment(weights, u, u)
    g12 = normalized_product_moment(weights, u, z)
    g22 = normalized_product_moment(weights, z, z)
    t111 = normalized_product_moment(weights, u, u, u)
    t112 = normalized_product_moment(weights, u, u, z)
    t122 = normalized_product_moment(weights, u, z, z)
    t222 = normalized_product_moment(weights, z, z, z)
    i1 = 3 * t112 * g11 * g22 - t222 * g11**2 - 2 * t111 * g12 * g22
    i2 = 3 * t122 * g11 * g22 - 2 * t222 * g12 * g11 - t111 * g22**2
    return i1, i2


def main():
    vectors = [v for v in product(range(V_MAX + 1), repeat=3) if any(v)]

    gauss_cells = 0
    gauss_ok = True
    gauss_hash = sha256()
    for v in vectors:
        for d in range(1, D_MAX + 1):
            for n in range(N_MAX + 1):
                lhs = factorial_ray_weight(v, d, n)
                rhs = gauss_ray_weight(v, d, n)
                gauss_ok &= lhs == rhs
                gauss_cells += 1
                gauss_hash.update(f"{v}|{d}|{n}|{lhs}|{rhs}\n".encode())

    orientation_cells = 0
    i1_negative = 0
    i2_negative = 0
    orientation_hash = sha256()
    for v in vectors:
        for d in range(1, D_MAX + 1):
            i1, i2 = orientation_invariants(v, d)
            orientation_cells += 1
            i1_negative += i1 < 0
            i2_negative += i2 < 0
            orientation_hash.update(f"{v}|{d}|{i1}|{i2}\n".encode())

    hostile_v = (1, 2, 3)
    hostile_d = 2
    hostile_n = 3
    correct = gauss_ray_weight(hostile_v, hostile_d, hostile_n)
    positive_q = [hostile_d * vi for vi in hostile_v]
    scale = prod(q**q for q in positive_q)
    wrong_scale = correct / scale

    print(
        f"gauss_cells={gauss_cells} v_box=0..{V_MAX}^3_nonzero d=1..{D_MAX} "
        f"n=0..{N_MAX} all_exact={gauss_ok} digest={gauss_hash.hexdigest()}"
    )
    print(
        f"orientation_cells={orientation_cells} strict_negative_I1={i1_negative} "
        f"strict_negative_I2={i2_negative} digest={orientation_hash.hexdigest()}"
    )
    print(
        "wrong_geometric_scale_hostile="
        f"v={hostile_v},d={hostile_d},n={hostile_n},detected={wrong_scale != correct}"
    )
    print("zero_exponent_vector_excluded=", (0, 0, 0) not in vectors, sep="")


if __name__ == "__main__":
    main()
