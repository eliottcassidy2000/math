from __future__ import annotations

from collections import Counter
from math import comb


ROLES = (2, 4, 6, 8)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise AssertionError(message)


def generalized_binom(n: int, k: int) -> int:
    numerator = 1
    for j in range(k):
        numerator *= n - j
    return numerator // factorial(k)


def factorial(k: int) -> int:
    answer = 1
    for j in range(2, k + 1):
        answer *= j
    return answer


def factorization(n: int) -> list[tuple[int, int]]:
    factors: list[tuple[int, int]] = []
    p = 2
    while p * p <= n:
        if n % p == 0:
            exponent = 0
            while n % p == 0:
                n //= p
                exponent += 1
            factors.append((p, exponent))
        p += 1
    if n > 1:
        factors.append((n, 1))
    return factors


def floor_log_p(k: int, p: int) -> int:
    exponent = 0
    power = p
    while power <= k:
        exponent += 1
        power *= p
    return exponent


def exact_period(m: int, k: int) -> int:
    period = 1
    for p, a in factorization(m):
        period *= p ** (a + floor_log_p(k, p))
    return period


def circular_convolution(left: list[int], right: list[int]) -> list[int]:
    m = len(left)
    answer = [0] * m
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            answer[(i + j) % m] += a * b
    return answer


def center_value_mod_odd_m(k: int, m: int) -> int:
    r = k // 2
    return ((-1) ** r * comb(2 * r, r) * pow(pow(16, r, m), -1, m)) % m


def check_polynomial_reflection_and_difference() -> int:
    checks = 0
    for k in range(1, 13):
        for n in range(-30, 31):
            value = generalized_binom(n, k)
            reflected = generalized_binom(k - 1 - n, k)
            require(reflected == (-1) ** k * value, ("reflection", k, n))
            require(generalized_binom(n + 1, k) - value == generalized_binom(n, k - 1),
                    ("Pascal difference", k, n))
            odd_increment = generalized_binom(n, k - 1)
            reflected_increment = generalized_binom(k - 2 - n, k - 1)
            require(reflected_increment == (-1) ** (k - 1) * odd_increment,
                    ("increment reflection", k, n))
            checks += 3
    return checks


def check_modular_fibres(limit: int = 100) -> tuple[int, dict[int, tuple[int, list[int]]]]:
    checks = 0
    controls: dict[int, tuple[int, list[int]]] = {}
    for m in range(2, limit + 1):
        histograms: list[list[int]] = []
        for k in ROLES:
            period = exact_period(m, k)
            require((period % 2) == (m % 2), ("period parity", m, k))
            histogram = [0] * m
            for n in range(period):
                value = comb(n, k) % m
                reflected_n = (k - 1 - n) % period
                require(comb(reflected_n, k) % m == value,
                        ("modular reflection", m, k, n))
                histogram[value] += 1
                checks += 1
            histograms.append(histogram)

            if m % 2 == 0:
                require(all(count % 2 == 0 for count in histogram),
                        ("even one-role histogram", m, k))
            else:
                center = ((k - 1) * pow(2, -1, period)) % period
                distinguished = comb(center, k) % m
                closed_form = center_value_mod_odd_m(k, m)
                require(distinguished == closed_form, ("fixed value", m, k))
                require([i for i, count in enumerate(histogram) if count % 2]
                        == [distinguished], ("odd one-role histogram", m, k))

        total = [1] + [0] * (m - 1)
        for histogram in histograms:
            total = circular_convolution(total, histogram)

        if m % 2 == 0:
            require(all(count % 16 == 0 for count in total),
                    ("complete even histogram", m))
        else:
            distinguished = sum(center_value_mod_odd_m(k, m) for k in ROLES) % m
            closed_form = (-3453 * pow(32768, -1, m)) % m
            require(distinguished == closed_form, ("Sun fixed target", m))
            require([i for i, count in enumerate(total) if count % 2]
                    == [distinguished], ("complete odd histogram", m))

        if m in (3, 11, 33):
            controls[m] = (distinguished, total)
    return checks, controls


def main() -> None:
    polynomial_checks = check_polynomial_reflection_and_difference()
    modular_checks, controls = check_modular_fibres()
    print(f"polynomial_checks={polynomial_checks}")
    print(f"modular_lane_checks={modular_checks}")
    for m in (3, 11, 33):
        distinguished, total = controls[m]
        odd_targets = [i for i, count in enumerate(total) if count % 2]
        print(
            f"m={m} distinguished={distinguished} odd_targets={odd_targets} "
            f"minimum_v2={min((count & -count).bit_length() - 1 for count in total)}"
        )
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
