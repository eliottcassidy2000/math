#!/usr/bin/env python3
"""Independent numerical Fourier hostile checks for THM-4061."""

import cmath
import math


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def phi(q):
    return sum(math.gcd(a, q) == 1 for a in range(1, q))


def stern_sign(a, q):
    inverse = pow(a, -1, q)
    return 1 if (a + inverse) % 2 == 0 else -1


def packet_sum(q):
    return sum(
        stern_sign(a, q) for a in range(1, q) if math.gcd(a, q) == 1
    )


def hyperbola_count(q):
    m = (q - 1) // 2
    return sum(
        pow(4 * r, -1, q) <= m
        for r in range(1, m + 1)
        if math.gcd(r, q) == 1
    )


def harmonic(n):
    return sum(1.0 / j for j in range(1, n + 1))


def fourier_reconstruction(p):
    zeta = cmath.exp(2j * math.pi / p)
    coefficients = []
    for h in range(p):
        coefficients.append(
            sum(((-1) ** x) * zeta ** (-h * x) for x in range(p)) / p
        )
    value = 0j
    for h in range(p):
        for k in range(p):
            kloosterman = sum(
                zeta ** (h * a + k * pow(a, -1, p)) for a in range(1, p)
            )
            value += coefficients[h] * coefficients[k] * kloosterman
    return value, coefficients


def main():
    checked = 0
    for q in range(3, 2002, 2):
        packet = packet_sum(q)
        n_box = hyperbola_count(q)
        require(packet == 4 * n_box - phi(q), f"hyperbola identity q={q}")
        checked += 1
    print("odd_hyperbola_rows", checked)

    primes = []
    for p in range(3, 80, 2):
        if all(p % d for d in range(2, math.isqrt(p) + 1)):
            primes.append(p)
            reconstructed, coefficients = fourier_reconstruction(p)
            packet = packet_sum(p)
            require(abs(reconstructed.real - packet) < 1e-9, f"Fourier real p={p}")
            require(abs(reconstructed.imag) < 1e-9, f"Fourier imag p={p}")
            for h, coefficient in enumerate(coefficients):
                expected = 1.0 / (p * abs(math.cos(math.pi * h / p)))
                require(
                    abs(abs(coefficient) - expected) < 1e-12,
                    f"coefficient p={p},h={h}",
                )
            hstar = 2 * harmonic(p - 1) - harmonic((p - 1) // 2)
            coefficient_mass = sum(abs(coefficients[h]) for h in range(1, p))
            require(
                coefficient_mass <= hstar + 1e-12,
                f"harmonic bound p={p}",
            )
            rhs = (
                2 * math.sqrt(p) * hstar * hstar
                + 2 * hstar / p
                + (p - 1) / (p * p)
            )
            require(abs(packet) <= rhs, f"final bound p={p}")
    print("prime_fourier_rows", len(primes), "last_prime", primes[-1])
    print("sample_S", tuple((p, packet_sum(p)) for p in primes))


if __name__ == "__main__":
    main()
