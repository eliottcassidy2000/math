#!/usr/bin/env python3
"""Exact companion for THM-2802 normalized unit-bispectrum reconstruction."""

from collections import Counter
from itertools import product


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def prime_divisors(n):
    divisors = []
    candidate = 2
    while candidate * candidate <= n:
        if n % candidate == 0:
            divisors.append(candidate)
            while n % candidate == 0:
                n //= candidate
        candidate += 1
    if n > 1:
        divisors.append(n)
    return tuple(divisors)


def root_of_exact_order(order, prime):
    require((prime - 1) % order == 0, "field does not split the cycle")
    for value in range(2, prime):
        if pow(value, order, prime) != 1:
            continue
        if all(
            pow(value, order // divisor, prime) != 1
            for divisor in prime_divisors(order)
        ):
            return value
    raise RuntimeError("no root of requested order")


def dft(coefficients, prime, root):
    n = len(coefficients)
    return tuple(
        sum(
            coefficient * pow(root, -frequency * position, prime)
            for position, coefficient in enumerate(coefficients)
        )
        % prime
        for frequency in range(n)
    )


def normalized_bispectrum(fourier, prime):
    n = len(fourier)
    require(all(fourier), "bispectrum requires a group-algebra unit")
    inverses = tuple(pow(value, -1, prime) for value in fourier)
    return tuple(
        fourier[k]
        * fourier[l]
        * inverses[(k + l) % n]
        * inverses[0]
        % prime
        for k in range(n)
        for l in range(n)
    )


def scalar_translate_fourier(fourier, scalar, shift, prime, root):
    return tuple(
        scalar * pow(root, -frequency * shift, prime) * value % prime
        for frequency, value in enumerate(fourier)
    )


def exhaustive_unit_census(order, prime):
    root = root_of_exact_order(order, prime)
    counts = Counter()
    nonzero = range(1, prime)
    for fourier in product(nonzero, repeat=order):
        counts[normalized_bispectrum(fourier, prime)] += 1

    expected_orbit = order * (prime - 1)
    expected_classes = (prime - 1) ** order // expected_orbit
    require(
        len(counts) == expected_classes,
        f"C{order}/F{prime} bispectrum class census",
    )
    require(
        set(counts.values()) == {expected_orbit},
        f"C{order}/F{prime} bispectrum fibre is not one scalar-translation orbit",
    )
    return len(counts), expected_orbit


def c13_controls():
    order = 13
    prime = 53
    root = root_of_exact_order(order, prime)

    semantic = tuple(1 if index in (6, 7) else 0 for index in range(order))
    semantic_fourier = dft(semantic, prime, root)
    require(all(semantic_fourier), "two-point semantic chain is not a unit")
    semantic_bispectrum = normalized_bispectrum(semantic_fourier, prime)

    invariant_checks = 0
    for scalar in range(1, prime):
        for shift in range(order):
            transformed = scalar_translate_fourier(
                semantic_fourier, scalar, shift, prime, root
            )
            require(
                normalized_bispectrum(transformed, prime) == semantic_bispectrum,
                "normalized bispectrum lost scalar/origin invariance",
            )
            invariant_checks += 1

    support_a = {0, 1, 3, 9}
    support_b = {1, 2, 5, 7}
    vector_a = tuple(int(index in support_a) for index in range(order))
    vector_b = tuple(int(index in support_b) for index in range(order))
    fourier_a = dft(vector_a, prime, root)
    fourier_b = dft(vector_b, prime, root)
    require(all(fourier_a) and all(fourier_b), "homometric control is not a unit")
    power_a = tuple(
        fourier_a[k] * fourier_a[-k % order] % prime for k in range(order)
    )
    power_b = tuple(
        fourier_b[k] * fourier_b[-k % order] % prime for k in range(order)
    )
    require(power_a == power_b, "homometric power spectra separated")
    require(
        normalized_bispectrum(fourier_a, prime)
        != normalized_bispectrum(fourier_b, prime),
        "normalized bispectrum failed to separate the homometric pair",
    )

    projective_translates = {
        scalar_translate_fourier(fourier_a, scalar, shift, prime, root)
        for scalar in range(1, prime)
        for shift in range(order)
    }
    require(
        fourier_b not in projective_translates,
        "homometric control unexpectedly became a scalar translate",
    )
    return invariant_checks, power_a[0], power_a[1]


def main():
    c3_classes, c3_orbit = exhaustive_unit_census(3, 7)
    c5_classes, c5_orbit = exhaustive_unit_census(5, 11)
    invariant_checks, homometric_zero_power, homometric_nonzero_power = (
        c13_controls()
    )

    print("THM-2802 NORMALIZED UNIT BISPECTRUM")
    print("status=VERIFIED-EXACT abstract splitting-field companion")
    print(
        "identity=B_f(k,l)=fhat(k)fhat(l)/"
        "(fhat(k+l)fhat(0)); equality iff scalar-translation orbit"
    )
    print(
        f"exhaustive_C3_F7_classes={c3_classes} "
        f"orbit_size={c3_orbit} units={6**3}"
    )
    print(
        f"exhaustive_C5_F11_classes={c5_classes} "
        f"orbit_size={c5_orbit} units={10**5}"
    )
    print(
        f"C13_F53_scalar_origin_invariance={invariant_checks}/"
        f"{52 * 13}"
    )
    print(
        "homometric_C13_pair=0139_vs_1257 "
        f"power_zero={homometric_zero_power} "
        f"power_nonzero={homometric_nonzero_power} "
        "normalized_bispectrum=distinct"
    )
    print(
        "scope=projective coefficient/origin reconstruction for units; "
        "no positive carrier, frame naturality, semantic attachment, or LRC14"
    )
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
