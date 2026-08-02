#!/usr/bin/env python3
"""Exact controls for THM-3204's odd-prime-power residual theorem."""

from fractions import Fraction
from math import comb, isqrt


PRIME_POWER_MAX = 400
LAGUERRE_MAX = 49
RATIO_MAX = 125
FACTORIAL = [1]
for integer in range(1, 2 * PRIME_POWER_MAX + 1):
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


def integer_binomial(top: int, bottom: int) -> int:
    """Generalized binomial for integral top and nonnegative bottom."""
    if bottom < 0:
        return 0
    numerator = 1
    for offset in range(bottom):
        numerator *= top - offset
    return numerator // FACTORIAL[bottom]


def generalized_laguerre_at_minus_d(m: int, alpha: int, d: int):
    """L_m^(alpha)(-d) in the convention used by THM-3204."""
    return sum(
        Fraction(integer_binomial(m + alpha, m - k) * d**k, FACTORIAL[k])
        for k in range(m + 1)
    )


def moment_coefficient(n: int, d: int, j: int) -> int:
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


def bessel_coefficient(n: int, j: int) -> int:
    return comb(n + j, 2 * j) * odd_double_factorial(j)


def terminal_coefficient(n: int, j: int) -> int:
    return (-1) ** (n - j) * FACTORIAL[n] * 2**j * bessel_coefficient(n, j)


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


def odd_primes(limit: int):
    answer = []
    for candidate in range(3, limit + 1, 2):
        if all(candidate % divisor for divisor in range(3, isqrt(candidate) + 1, 2)):
            answer.append(candidate)
    return answer


def prime_power_rows():
    rows = []
    for prime in odd_primes(isqrt(PRIME_POWER_MAX)):
        exponent = 2
        value = prime**exponent
        while value <= PRIME_POWER_MAX:
            rows.append((value, prime, exponent))
            exponent += 1
            value *= prime
    return sorted(rows)


def predicted_hull(d: int, prime: int, n: int, h: int):
    if n == d - 1:
        return [(0, h), (d - 1, h + (d - 1) // (prime - 1))]
    answer = [
        (0, h),
        (1, h),
        (d - prime + 1, h + (d - prime) // (prime - 1)),
    ]
    if prime > 3:
        answer.append((d - 2, h + (d - 1) // (prime - 1)))
    return answer


def ratio_valuation(d: int, prime: int, exponent: int, n: int, j: int, k: int):
    m = n - j
    big_n = n + j
    return (
        exponent * k
        - vp(FACTORIAL[k], prime)
        + vp(comb(m, k), prime)
        - vp(comb(big_n, k), prime)
    )


def main_controls():
    rows = prime_power_rows()
    direct = True
    laguerre = True
    ratio_nonnegative = True
    ratio_equalities = True
    exact_hulls = True
    scalar_prefix = True
    residual_gcds = True
    digit_face_sets = True
    coefficient_count = 0
    hostile_lifts = {3: [], 5: []}

    for d, prime, exponent in rows:
        coefficients = {}
        for n in (d - 2, d - 1):
            coefficients[n] = moment_coefficients(n, d)
            coefficient_count += len(coefficients[n])
            if d == 9:
                direct &= coefficients[n] == direct_bivariate_moment(n, d)
            if d <= LAGUERRE_MAX:
                for j, coefficient in enumerate(coefficients[n]):
                    laguerre_value = generalized_laguerre_at_minus_d(
                        n - j, -n - j - 1, d
                    )
                    reconstructed = (
                        Fraction(FACTORIAL[n] * FACTORIAL[2 * j], FACTORIAL[j])
                        * laguerre_value
                    )
                    laguerre &= reconstructed.denominator == 1
                    laguerre &= reconstructed.numerator == coefficient

            h = vp(FACTORIAL[n], prime)
            exact_hulls &= lower_hull(coefficients[n], prime) == predicted_hull(
                d, prime, n, h
            )

            bessel_face = []
            for j in range(n + 1):
                terminal_valuation = h + vp(bessel_coefficient(n, j), prime)
                actual_valuation = vp(coefficients[n][j], prime)
                if actual_valuation == terminal_valuation:
                    bessel_face.append(j)
            if n == d - 1:
                expected_face = [0] + [
                    d - prime**rho for rho in range(exponent - 1, -1, -1)
                ]
            else:
                expected_face = [1] + [
                    d - prime**rho + 1 for rho in range(exponent - 1, 0, -1)
                ]
                expected_face = [0] + expected_face + [d - 2]
            # Equality with terminal valuation can occur away from the lower
            # face; compare only indices on the predicted lower-face lines.
            hull, edges = residual_edges(coefficients[n], prime)
            on_lower_faces = set()
            for left, right in zip(hull, hull[1:]):
                slope = Fraction(right[1] - left[1], right[0] - left[0])
                for j in range(left[0], right[0] + 1):
                    if vp(coefficients[n][j], prime) == left[1] + slope * (j - left[0]):
                        on_lower_faces.add(j)
            digit_face_sets &= sorted(on_lower_faces) == sorted(set(expected_face))

            if d <= RATIO_MAX:
                delta = d - n
                equalities = []
                for j in range(n + 1):
                    for k in range(1, n - j + 1):
                        valuation = ratio_valuation(d, prime, exponent, n, j, k)
                        ratio_nonnegative &= valuation >= 0
                        if valuation == 0:
                            equalities.append((j, k))
                ratio_equalities &= equalities == [(delta, 1)]

        slope = Fraction(1, prime - 1)
        _, left_edges = residual_edges(coefficients[d - 2], prime)
        _, right_edges = residual_edges(coefficients[d - 1], prime)
        left_residual = left_edges[slope]
        right_residual = right_edges[slope]
        scalar = (-2) % prime
        scalar_prefix &= left_residual == [
            scalar * coefficient % prime for coefficient in right_residual[:-1]
        ]
        scalar_prefix &= right_residual[-1] != 0
        residual_gcds &= polynomial_gcd(left_residual, right_residual, prime) == [1]

        if prime == 3:
            n, j = d - 1, 1
            hostile_lifts[3].append(
                (d, vp(coefficients[n][j], prime) - vp(terminal_coefficient(n, j), prime))
            )
        if prime == 5:
            n, j = d - 2, 2
            hostile_lifts[5].append(
                (d, vp(coefficients[n][j], prime) - vp(terminal_coefficient(n, j), prime))
            )

    return {
        "rows": rows,
        "coefficient_count": coefficient_count,
        "direct": direct,
        "laguerre": laguerre,
        "ratio_nonnegative": ratio_nonnegative,
        "ratio_equalities": ratio_equalities,
        "exact_hulls": exact_hulls,
        "digit_face_sets": digit_face_sets,
        "scalar_prefix": scalar_prefix,
        "residual_gcds": residual_gcds,
        "hostile_lifts": hostile_lifts,
    }


def d15_hostile():
    d = 15
    answer = []
    for prime in (3, 5):
        left_hull, left_edges = residual_edges(moment_coefficients(d - 2, d), prime)
        right_hull, right_edges = residual_edges(moment_coefficients(d - 1, d), prime)
        common = sorted(set(left_edges).intersection(right_edges))
        data = [
            (
                str(slope),
                slope.denominator % prime != 0,
                polynomial_gcd(left_edges[slope], right_edges[slope], prime),
            )
            for slope in common
        ]
        answer.append((prime, left_hull, right_hull, data))
    return answer


def main():
    controls = main_controls()
    rows = controls["rows"]
    print(
        f"prime_power_universe=odd_powers_le_{PRIME_POWER_MAX} count={len(rows)} "
        f"values={','.join(str(value) for value, _, _ in rows)} "
        f"coefficients={controls['coefficient_count']}"
    )
    print(f"independent_bivariate_d9={controls['direct']}")
    print(f"negative_parameter_laguerre_identity_d_le_{LAGUERRE_MAX}={controls['laguerre']}")
    print(f"ratio_valuations_nonnegative_d_le_{RATIO_MAX}={controls['ratio_nonnegative']}")
    print(f"ratio_equality_only_j_delta_k1={controls['ratio_equalities']}")
    print(f"predicted_newton_hulls={controls['exact_hulls']}")
    print(f"digit_sum_face_index_sets={controls['digit_face_sets']}")
    print(f"scalar_prefix_residual_identity={controls['scalar_prefix']}")
    print(f"common_residual_gcds_one={controls['residual_gcds']}")
    print(f"hostile_p3_nDminus1_j1_lifts={controls['hostile_lifts'][3]}")
    print(f"hostile_p5_nDminus2_j2_lifts={controls['hostile_lifts'][5]}")
    for prime, left_hull, right_hull, data in d15_hostile():
        faces = ",".join(
            f"slope={slope}:tame={tame}:gcd={common_gcd}"
            for slope, tame, common_gcd in data
        )
        print(
            f"hostile_d15 p={prime} A_hull={left_hull} B_hull={right_hull} "
            f"shared=[{faces}]"
        )


if __name__ == "__main__":
    main()
