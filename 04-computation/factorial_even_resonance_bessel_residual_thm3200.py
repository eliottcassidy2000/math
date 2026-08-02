#!/usr/bin/env python3
"""Exact controls for THM-3200's even-resonance Bessel residual proof."""

from fractions import Fraction
from math import comb, isqrt


EVEN_D_MAX = 160
DIRECT_D_MAX = 10
REINDEX_D_MAX = 20
HOSTILE_ODD = (15, 33, 35)
FACTORIAL = [1]
for integer in range(1, 2 * EVEN_D_MAX + 1):
    FACTORIAL.append(FACTORIAL[-1] * integer)


def vp(integer: int, prime: int):
    if integer == 0:
        return None
    answer = 0
    while integer % prime == 0:
        answer += 1
        integer //= prime
    return answer


def unit_residue(integer: int, prime: int) -> int:
    while integer % prime == 0:
        integer //= prime
    return integer % prime


def odd_double_factorial(index: int) -> int:
    answer = 1
    for integer in range(1, 2 * index, 2):
        answer *= integer
    return answer


def moment_coefficient(n: int, d: int, j: int) -> int:
    """[v^j] L((d-t+v*t^2)^n), from the direct multinomial sum."""
    return comb(n, j) * sum(
        comb(n - j, ell)
        * d ** (n - j - ell)
        * (-1) ** ell
        * FACTORIAL[2 * j + ell]
        for ell in range(n - j + 1)
    )


def moment_coefficients(n: int, d: int):
    return [moment_coefficient(n, d, j) for j in range(n + 1)]


def direct_bivariate_moment(n: int, d: int):
    """Independent repeated-product expansion for the small controls."""
    polynomial = {(0, 0): 1}
    base = {(0, 0): d, (1, 0): -1, (2, 1): 1}
    for _ in range(n):
        product = {}
        for (t0, v0), coefficient0 in polynomial.items():
            for (t1, v1), coefficient1 in base.items():
                key = (t0 + t1, v0 + v1)
                product[key] = product.get(key, 0) + coefficient0 * coefficient1
        polynomial = product
    answer = [0] * (n + 1)
    for (t_degree, v_degree), coefficient in polynomial.items():
        answer[v_degree] += coefficient * FACTORIAL[t_degree]
    return answer


def reindexed_coefficient(n: int, d: int, j: int) -> int:
    """The k=#constant-factors form in THM-3200 equation (10)."""
    normalized = sum(
        Fraction((-d) ** k, FACTORIAL[k])
        * comb(n + j - k, 2 * j)
        * odd_double_factorial(j)
        for k in range(n - j + 1)
    )
    answer = (-1) ** (n - j) * FACTORIAL[n] * 2**j * normalized
    assert answer.denominator == 1
    return answer.numerator


def bessel_coefficient(n: int, j: int) -> int:
    return comb(n + j, 2 * j) * odd_double_factorial(j)


def bessel_polynomial(n: int):
    return [bessel_coefficient(n, j) for j in range(n + 1)]


def lower_hull(coefficients, prime: int):
    points = [
        (index, vp(coefficient, prime))
        for index, coefficient in enumerate(coefficients)
        if coefficient != 0
    ]
    hull = []
    for point in points:
        while len(hull) >= 2:
            old_slope = Fraction(
                hull[-1][1] - hull[-2][1], hull[-1][0] - hull[-2][0]
            )
            new_slope = Fraction(
                point[1] - hull[-1][1], point[0] - hull[-1][0]
            )
            if old_slope >= new_slope:
                hull.pop()
            else:
                break
        hull.append(point)
    return hull


def trim(polynomial, prime: int):
    answer = [coefficient % prime for coefficient in polynomial]
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return answer


def polynomial_remainder(dividend, divisor, prime: int):
    answer = trim(dividend, prime)
    divisor = trim(divisor, prime)
    while len(answer) >= len(divisor) and answer != [0]:
        shift = len(answer) - len(divisor)
        multiplier = answer[-1] * pow(divisor[-1], -1, prime) % prime
        for index, coefficient in enumerate(divisor):
            answer[index + shift] -= multiplier * coefficient
        answer = trim(answer, prime)
    return answer


def polynomial_gcd(left, right, prime: int):
    left = trim(left, prime)
    right = trim(right, prime)
    while right != [0]:
        left, right = right, polynomial_remainder(left, right, prime)
    inverse = pow(left[-1], -1, prime)
    return [(coefficient * inverse) % prime for coefficient in left]


def residual_edges(coefficients, prime: int):
    hull = lower_hull(coefficients, prime)
    edges = {}
    for left, right in zip(hull, hull[1:]):
        i0, y0 = left
        i1, y1 = right
        slope = Fraction(y1 - y0, i1 - i0)
        step = slope.denominator
        residual = []
        for index in range(i0, i1 + 1, step):
            valuation = vp(coefficients[index], prime)
            on_edge = valuation is not None and valuation == y0 + slope * (index - i0)
            residual.append(
                unit_residue(coefficients[index], prime) if on_edge else 0
            )
        edges[slope] = trim(residual, prime)
    return hull, edges


def prime_divisors(integer: int):
    answer = []
    candidate = 2
    while candidate <= isqrt(integer):
        if integer % candidate == 0:
            answer.append(candidate)
            while integer % candidate == 0:
                integer //= candidate
        candidate += 1
    if integer > 1:
        answer.append(integer)
    return answer


def even_controls():
    direct = True
    reindexed = True
    normalized_reduction = True
    exact_hulls = True
    residual_gcds = True
    pair_count = 0
    coefficient_count = 0
    even_values = list(range(4, EVEN_D_MAX + 1, 2))

    max_n = EVEN_D_MAX - 1
    bessel = [bessel_polynomial(n) for n in range(max_n + 1)]
    exact_recurrence = True
    for n in range(1, max_n):
        expected = [0] * (n + 2)
        for j, coefficient in enumerate(bessel[n]):
            expected[j + 1] += (2 * n + 1) * coefficient
        for j, coefficient in enumerate(bessel[n - 1]):
            expected[j] += coefficient
        exact_recurrence &= expected == bessel[n + 1]

    for d in even_values:
        for n in (d - 2, d - 1):
            pair_count += 1
            coefficients = moment_coefficients(n, d)
            coefficient_count += len(coefficients)
            if d <= DIRECT_D_MAX:
                direct &= coefficients == direct_bivariate_moment(n, d)
            if d <= REINDEX_D_MAX:
                reindexed &= coefficients == [
                    reindexed_coefficient(n, d, j) for j in range(n + 1)
                ]
            h_n = vp(FACTORIAL[n], 2)
            residual = []
            for j, coefficient in enumerate(coefficients):
                scale = 2 ** (h_n + j)
                normalized_reduction &= coefficient % scale == 0
                residual.append((coefficient // scale) % 2)
            expected_residual = [coefficient % 2 for coefficient in bessel[n]]
            normalized_reduction &= residual == expected_residual
            exact_hulls &= lower_hull(coefficients, 2) == [
                (0, h_n),
                (n, h_n + n),
            ]
        left = [coefficient % 2 for coefficient in bessel[d - 2]]
        right = [coefficient % 2 for coefficient in bessel[d - 1]]
        residual_gcds &= polynomial_gcd(left, right, 2) == [1]

    return {
        "even_values": even_values,
        "pair_count": pair_count,
        "coefficient_count": coefficient_count,
        "direct": direct,
        "reindexed": reindexed,
        "normalized_reduction": normalized_reduction,
        "exact_hulls": exact_hulls,
        "exact_recurrence": exact_recurrence,
        "residual_gcds": residual_gcds,
    }


def hostile_controls():
    rows = []
    for d in HOSTILE_ODD:
        left_coefficients = moment_coefficients(d - 2, d)
        right_coefficients = moment_coefficients(d - 1, d)
        for prime in prime_divisors(d):
            left_hull, left_edges = residual_edges(left_coefficients, prime)
            right_hull, right_edges = residual_edges(right_coefficients, prime)
            common = sorted(set(left_edges).intersection(right_edges))
            data = []
            for slope in common:
                common_gcd = polynomial_gcd(
                    left_edges[slope], right_edges[slope], prime
                )
                data.append(
                    (
                        str(slope),
                        slope.denominator % prime != 0,
                        common_gcd,
                    )
                )
            rows.append((d, prime, left_hull, right_hull, data))
    expected_degrees = {
        (15, 3): [1, 1],
        (15, 5): [2],
        (33, 3): [1, 1],
        (33, 11): [2],
        (35, 5): [1, 1],
        (35, 7): [2],
    }
    assert {
        (d, prime): [len(common_gcd) - 1 for _, _, common_gcd in data]
        for d, prime, _, _, data in rows
    } == expected_degrees
    return rows


def main():
    controls = even_controls()
    even_values = controls["even_values"]
    print(
        f"even_universe={even_values[0]}..{even_values[-1]} "
        f"count={len(even_values)} polynomial_instances={controls['pair_count']} "
        f"coefficients={controls['coefficient_count']}"
    )
    print(f"independent_bivariate_controls_d_le_{DIRECT_D_MAX}={controls['direct']}")
    print(f"exact_reindexed_controls_d_le_{REINDEX_D_MAX}={controls['reindexed']}")
    print(f"normalized_bessel_reduction_mod_2={controls['normalized_reduction']}")
    print(f"single_slope_one_newton_hulls={controls['exact_hulls']}")
    print(f"integer_bessel_three_term_recurrence={controls['exact_recurrence']}")
    print(f"consecutive_residual_gcds_one={controls['residual_gcds']}")
    print("hostile_odd_controls=begin")
    for d, prime, left_hull, right_hull, data in hostile_controls():
        faces = ",".join(
            f"slope={slope}:tame={tame}:gcd={common_gcd}"
            for slope, tame, common_gcd in data
        )
        print(
            f"d={d} p={prime} A_hull={left_hull} B_hull={right_hull} "
            f"shared=[{faces}]"
        )
    print("hostile_odd_controls=end")


if __name__ == "__main__":
    main()
