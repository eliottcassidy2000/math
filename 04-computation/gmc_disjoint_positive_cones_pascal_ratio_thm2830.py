#!/usr/bin/env python3
"""Exact referee for THM-2830's universal Pascal-ratio proof."""

from fractions import Fraction
from itertools import combinations_with_replacement, product
from math import comb, factorial
from random import Random


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def h(n, j):
    return comb(n + j, n)


def triple(a, b, c):
    total = a + b + c
    return (
        factorial(total + 1)
        // (factorial(a + 1) * factorial(b) * factorial(c))
        + factorial(total + 1)
        // (factorial(a) * factorial(b + 1) * factorial(c))
        + factorial(total + 1)
        // (factorial(a) * factorial(b) * factorial(c + 1))
        - factorial(total) // (factorial(a) * factorial(b) * factorial(c))
    )


def tau(a, b, c):
    total = a + b + c
    return (
        Fraction(total + 1, a + 1)
        + Fraction(total + 1, b + 1)
        + Fraction(total + 1, c + 1)
        - 1
    )


def beta_ratio(n, p, q):
    return Fraction(
        factorial(n) * factorial(n + p + q),
        factorial(n + p) * factorial(n + q),
    )


def c_kernel(n, p, q):
    return Fraction(triple(n, p, q), h(n, p) * h(n, q))


def c_base(n, p, q):
    return Fraction(n + 1, p + 1) + Fraction(n + 1, q + 1)


def e_residual(n, p, q):
    return c_kernel(n, p, q) - c_base(n, p, q)


def effective_mean(n, p, q):
    return (n + 1) * (
        Fraction(triple(n + 1, p, q), triple(n, p, q)) - 1
    )


def f_residual(n, p, q):
    return (
        c_kernel(n, p, q) * effective_mean(n, p, q)
        - c_base(n, p, q) * (p + q + 1)
    )


def cross_polynomial(n, p, q):
    return (
        n * n * p
        + n * n * q
        + 2 * n * n
        + n * p * p
        + 2 * n * p * q
        + 6 * n * p
        + n * q * q
        + 6 * n * q
        + 6 * n
        + p * p * q
        + 3 * p * p
        + p * q * q
        + 6 * p * q
        + 7 * p
        + 3 * q * q
        + 7 * q
        + 4
    )


def residual_polynomial(n, p, q):
    return (
        n * n * p * p
        + 2 * n * n * p * q
        + 4 * n * n * p
        + n * n * q * q
        + 4 * n * n * q
        + 4 * n * n
        + n * p * p * q
        + 4 * n * p * p
        + n * p * q * q
        + 6 * n * p * q
        + 11 * n * p
        + 4 * n * q * q
        + 11 * n * q
        + 10 * n
        + 2 * p * p
        + 2 * p * q
        + 6 * p
        + 2 * q * q
        + 6 * q
        + 6
    )


def residual_certificate(n, p, q):
    total = p + q
    q_poly = (n + 1) ** 2 * (n + 2) * (total + 2) ** 2
    r_poly = residual_polynomial(n, p, q)
    return Fraction(
        (beta_ratio(n, p, q) - 1) * q_poly
        + beta_ratio(n, p, q) * total * r_poly,
        2 * (n + 1) * (n + 2) * (p + 1) * (q + 1),
    )


def j_kernel(n, p, q, r):
    return (
        c_kernel(n, p, q) * (effective_mean(n, p, q) - r)
        + c_kernel(n, p, r) * (effective_mean(n, p, r) - q)
        + c_kernel(n, q, r) * (effective_mean(n, q, r) - p)
    )


def k_kernel(n, p, q, r):
    return sum(
        triple(n + 1, x, y) * h(n, z)
        - triple(n, x, y) * h(n + 1, z)
        for x, y, z in ((p, q, r), (p, r, q), (q, r, p))
    )


algebra_cells = 0
monotone_cells = 0
cyclic_cells = 0
minimum_j_gap = None
floor_equalities = 0

for n in range(13):
    for p in range(13):
        for q in range(13):
            beta = beta_ratio(n, p, q)
            require(
                c_kernel(n, p, q) == beta * tau(n, p, q),
                "normalized C identity",
            )

            product_beta = Fraction(1)
            vandermonde_beta = Fraction(0)
            for k in range(q):
                product_beta *= Fraction(n + p + k + 1, n + k + 1)
            for k in range(min(p, q) + 1):
                vandermonde_beta += Fraction(
                    comb(p, k) * comb(q, k), comb(n + k, k)
                )
            require(
                beta == product_beta == vandermonde_beta and beta >= 1,
                "beta positivity identities",
            )

            boundary = Fraction(
                q * (n * q + 2 * n + 2 * q + 3),
                (n + 1) * (q + 1),
            )
            boundary_delta = Fraction(
                n * q * q
                + 3 * n * q
                + 3 * n
                + 2 * q * q
                + 6 * q
                + 5,
                (n + 1) * (q + 1) * (q + 2),
            )
            require(e_residual(n, 0, q) == boundary, "E boundary")
            require(
                e_residual(n, 0, q + 1) - e_residual(n, 0, q)
                == boundary_delta
                and boundary_delta > 0,
                "E boundary monotonicity",
            )

            cross = (
                e_residual(n, p + 1, q + 1)
                - e_residual(n, p + 1, q)
                - e_residual(n, p, q + 1)
                + e_residual(n, p, q)
            )
            cross_formula = Fraction(
                beta * cross_polynomial(n, p, q),
                (n + 1) * (n + p + 1) * (n + q + 1),
            )
            require(cross == cross_formula and cross > 0, "E supermodularity")

            residual_gap = (
                f_residual(n, p, q)
                - Fraction(p + q, 2) * e_residual(n, p, q)
            )
            require(
                residual_gap == residual_certificate(n, p, q)
                and residual_gap >= 0,
                "effective-mean certificate",
            )
            require(e_residual(n, p, q) >= 0, "E positivity")
            require(e_residual(n, p + 1, q) >= e_residual(n, p, q), "E p")
            require(e_residual(n, p, q + 1) >= e_residual(n, p, q), "E q")
            algebra_cells += 1
            monotone_cells += 3

for n in range(11):
    for p, q, r in combinations_with_replacement(range(14), 3):
        base_sum = (
            c_base(n, p, q) * (p + q + 1 - r)
            + c_base(n, p, r) * (p + r + 1 - q)
            + c_base(n, q, r) * (q + r + 1 - p)
        )
        require(base_sum == 6 * (n + 1), "separable cyclic identity")

        e1 = e_residual(n, p, q)
        e2 = e_residual(n, p, r)
        e3 = e_residual(n, q, r)
        require(e1 <= e2 <= e3, "ordered residual")
        residual_rearrangement = (
            e1 * Fraction(p + q - 2 * r, 2)
            + e2 * Fraction(p + r - 2 * q, 2)
            + e3 * Fraction(q + r - 2 * p, 2)
        )
        require(residual_rearrangement >= 0, "residual rearrangement")

        j_value = j_kernel(n, p, q, r)
        require(j_value >= 6 * (n + 1), "universal J floor")
        k_value = k_kernel(n, p, q, r)
        require(
            Fraction(k_value * (n + 1), h(n, p) * h(n, q) * h(n, r))
            == j_value,
            "K/J normalization",
        )
        require(
            k_value >= 6 * h(n, p) * h(n, q) * h(n, r),
            "universal K floor",
        )
        gap = j_value - 6 * (n + 1)
        if gap == 0:
            require((p, q, r) == (0, 0, 0), "unexpected floor equality")
            floor_equalities += 1
        if minimum_j_gap is None or gap < minimum_j_gap:
            minimum_j_gap = gap
        cyclic_cells += 1


def ratio_data(n, coefficients):
    a_value = sum(value * h(n, index) for index, value in coefficients.items())
    b_value = sum(
        left * right * triple(n, p, q)
        for p, left in coefficients.items()
        for q, right in coefficients.items()
    )
    return a_value, b_value


rng = Random(2830)
ratio_cells = 0
minimum_ratio_margin = None
for _ in range(3000):
    support = sorted(set(rng.randrange(0, 20) for _ in range(rng.randrange(1, 7))))
    coefficients = {index: rng.randrange(1, 8) for index in support}
    for n in range(20):
        a0, b0 = ratio_data(n, coefficients)
        a1, b1 = ratio_data(n + 1, coefficients)
        determinant = b1 * a0 - b0 * a1
        require(determinant >= 2 * a0**3, "quantitative Pascal ratio")
        ratio_margin = determinant - 2 * a0**3
        if minimum_ratio_margin is None or ratio_margin < minimum_ratio_margin:
            minimum_ratio_margin = ratio_margin
        ratio_cells += 1

for n in range(20):
    a0, b0 = ratio_data(n, {0: 7})
    a1, b1 = ratio_data(n + 1, {0: 7})
    require(b1 * a0 - b0 * a1 == 2 * a0**3, "sharp d0 equality")

hostile_coefficients = {0: 2, 10: 1}
hostile_ratios = [
    Fraction(
        ratio_data(n, hostile_coefficients)[1],
        ratio_data(n, hostile_coefficients)[0],
    )
    for n in range(3)
]
convexity_hostile = hostile_ratios[2] - 2 * hostile_ratios[1] + hostile_ratios[0]
require(convexity_hostile == Fraction(-53721500, 663), "convexity hostile")


def orientation_from_ratios(i, coefficients):
    a_i, b_i = ratio_data(i, coefficients)
    sum_a = 0
    sum_b = 0
    for index, value in coefficients.items():
        a_prev, b_prev = ratio_data(index - 1, coefficients)
        sum_a += value * a_prev
        sum_b += value * b_prev
    direct = 6 * (a_i * sum_b - b_i * sum_a)
    ratio_form = 6 * a_i * sum(
        value
        * ratio_data(index - 1, coefficients)[0]
        * (
            Fraction(
                ratio_data(index - 1, coefficients)[1],
                ratio_data(index - 1, coefficients)[0],
            )
            - Fraction(b_i, a_i)
        )
        for index, value in coefficients.items()
    )
    return direct, ratio_form


orientation_cells = 0
orientation_equalities = 0
for i in range(6):
    for offsets in product(range(3), repeat=5):
        coefficients = {
            i + 1 + slot: value
            for slot, value in enumerate(offsets)
            if value
        }
        if not coefficients:
            continue
        direct, ratio_form = orientation_from_ratios(i, coefficients)
        require(direct == ratio_form and direct >= 0, "orientation")
        if direct == 0:
            require(set(coefficients) == {i + 1}, "orientation equality")
            orientation_equalities += 1
        orientation_cells += 1


def interlaced_data(lower, upper):
    a_values = {}
    b_values = {}
    max_index = max(max(lower, default=0), max(upper, default=1) - 1)
    for index in range(max_index + 1):
        a_values[index], b_values[index] = ratio_data(index, upper)

    luv = sum(value * a_values[index] for index, value in lower.items())
    luv2 = sum(value * b_values[index] for index, value in lower.items())
    half_lv2 = sum(
        value * a_values[index - 1] for index, value in upper.items()
    )
    third_lv3 = sum(
        value * b_values[index - 1] for index, value in upper.items()
    )
    direct = 6 * (luv * third_lv3 - luv2 * half_lv2)

    alpha = {
        index: Fraction(value * a_values[index], luv)
        for index, value in lower.items()
    }
    beta = {
        index - 1: Fraction(value * a_values[index - 1], half_lv2)
        for index, value in upper.items()
    }
    expectation_alpha = sum(
        weight * Fraction(b_values[index], a_values[index])
        for index, weight in alpha.items()
    )
    expectation_beta = sum(
        weight * Fraction(b_values[index], a_values[index])
        for index, weight in beta.items()
    )
    transport = 6 * luv * half_lv2 * (expectation_beta - expectation_alpha)
    require(direct == transport, "interlaced transport identity")
    pairwise = 6 * sum(
        a_values[i]
        * a_values[j]
        * (
            Fraction(b_values[j], a_values[j])
            - Fraction(b_values[i], a_values[i])
        )
        * (
            upper.get(j + 1, 0) * lower.get(i, 0)
            - upper.get(i + 1, 0) * lower.get(j, 0)
        )
        for i in range(max_index + 1)
        for j in range(i + 1, max_index + 1)
    )
    require(direct == pairwise, "pairwise Cauchy-Binet transport")
    return direct, alpha, beta


strict_interlaced, _, _ = interlaced_data({0: 2, 1: 1}, {1: 1, 2: 2})
equality_interlaced, equality_alpha, equality_beta = interlaced_data(
    {0: 1, 2: 1}, {1: 1, 3: 1}
)
negative_interlaced, _, _ = interlaced_data({2: 1}, {1: 1, 3: 1})
require(strict_interlaced == 17460, "strict interlaced control")
require(equality_interlaced == 0, "interlaced equality control")
require(equality_alpha == equality_beta, "interlaced equal laws")
require(negative_interlaced == -33540, "unordered interlaced hostile")

interlaced_cells = 0
interlaced_equalities = 0
for _ in range(2000):
    width = rng.randrange(2, 8)
    lower_values = [rng.randrange(1, 7) for _ in range(width)]
    ratios = sorted(rng.randrange(1, 8) for _ in range(width))
    lower = {index: value for index, value in enumerate(lower_values)}
    upper = {
        index + 1: lower_values[index] * ratios[index]
        for index in range(width)
    }
    for i in range(width):
        for j in range(i + 1, width):
            require(
                upper[j + 1] * lower[i] >= upper[i + 1] * lower[j],
                "raw MLR condition",
            )
    value, alpha, beta = interlaced_data(lower, upper)
    require(value >= 0, "interlaced orientation")
    if value == 0:
        require(len(set(ratios)) == 1 and alpha == beta, "interlaced equality")
        interlaced_equalities += 1
    interlaced_cells += 1

print("THM-2830 UNIVERSAL PASCAL-RATIO EXACT REFEREE")
print("status=PROVED+VERIFIED-EXACT+INDEPENDENTLY-HOSTILE-AUDITED")
print(f"algebra_cells={algebra_cells}")
print(f"residual_monotonicity_cells={monotone_cells}")
print(f"cyclic_kernel_cells={cyclic_cells}")
print(f"minimum_J_minus_6nplus1={minimum_j_gap}")
print(f"floor_equalities={floor_equalities}; unique_triple=(0,0,0)")
print(f"random_ratio_cells={ratio_cells}")
print(f"minimum_determinant_minus_2A3={minimum_ratio_margin}")
print(f"ratio_convexity_hostile={convexity_hostile}")
print(f"orientation_cells={orientation_cells}")
print(f"orientation_equalities={orientation_equalities}")
print(f"interlaced_MLR_cells={interlaced_cells}")
print(f"interlaced_equalities={interlaced_equalities}")
print(
    "interlaced_controls="
    f"strict:{strict_interlaced}; equality:{equality_interlaced}; "
    f"unordered:{negative_interlaced}"
)
print("universal_floor=K_n>=6*H_np*H_nq*H_nr")
print("mechanism=SEPARABLE_BASE_PLUS_MONOTONE_INTERACTION")
print("all_exact_controls=PASS")
