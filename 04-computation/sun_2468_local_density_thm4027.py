#!/usr/bin/env python3
"""Exact local densities and archimedean scale for Sun's 2-4-6-8 sum."""

from __future__ import annotations

from fractions import Fraction
from math import comb, factorial, gamma, gcd


N = 896_315_812_331_399
DEGREES = (2, 4, 6, 8)


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def factor_prime_powers(q: int) -> list[tuple[int, int]]:
    out = []
    p = 2
    while p * p <= q:
        if q % p == 0:
            a = 0
            while q % p == 0:
                q //= p
                a += 1
            out.append((p, a))
        p += 1
    if q > 1:
        out.append((q, 1))
    return out


def lcm(a: int, b: int) -> int:
    return a // gcd(a, b) * b


def binomial_period(r: int, q: int) -> int:
    """A certified period of k -> C(k,r) modulo q.

    For p^a, Vandermonde and
      v_p(C(p^b,j)) = b-v_p(j), 0<j<p^b,
    show that p^(a+floor(log_p r)) is a period.  CRT gives the lcm.
    """
    ans = 1
    for p, a in factor_prime_powers(q):
        e = 0
        power = 1
        while power * p <= r:
            power *= p
            e += 1
        ans = lcm(ans, p ** (a + e))
    return ans


def distribution(r: int, q: int) -> tuple[list[int], int]:
    period = binomial_period(r, q)
    counts = [0] * q
    vals = [comb(k, r) % q for k in range(2 * period + r + 1)]
    require(
        vals[: period + r + 1] == vals[period : 2 * period + r + 1],
        "binomial period failure",
    )
    for k in range(period):
        counts[comb(k, r) % q] += 1
    require(sum(counts) == period, "histogram mass drift")
    return counts, period


def cyclic_convolve(a: list[int], b: list[int]) -> list[int]:
    q = len(a)
    require(len(b) == q, "convolution modulus mismatch")
    out = [0] * q
    for i, ai in enumerate(a):
        if not ai:
            continue
        for j, bj in enumerate(b):
            if bj:
                out[(i + j) % q] += ai * bj
    return out


def local_profile(q: int) -> tuple[list[Fraction], list[int]]:
    total = [1] + [0] * (q - 1)
    denominator = 1
    periods = []
    for r in DEGREES:
        dist, period = distribution(r, q)
        total = cyclic_convolve(total, dist)
        denominator *= period
        periods.append(period)
    require(sum(total) == denominator, "convolution mass drift")
    factors = [Fraction(q * count, denominator) for count in total]
    return factors, periods


def show_profile(q: int) -> None:
    factors, periods = local_profile(q)
    residue = N % q
    minimum = min(factors)
    minimizers = [i for i, f in enumerate(factors) if f == minimum]
    minimizer_text = str(minimizers)
    if len(minimizers) > 8:
        minimizer_text = f"{minimizers[:8]}... count={len(minimizers)}"
    maximum = max(factors)
    print(
        f"q={q:3d} Nmod={residue:3d} periods={periods} "
        f"factor(N)={factors[residue]}={float(factors[residue]):.12f} "
        f"min={minimum}={float(minimum):.12f}@{minimizer_text} "
        f"max={float(maximum):.12f} all_positive={all(f > 0 for f in factors)}"
    )


def main() -> None:
    reciprocal_sum = sum(Fraction(1, r) for r in DEGREES)
    volume_constant = 1.0
    for r in DEGREES:
        volume_constant *= gamma(1.0 + 1.0 / r) * factorial(r) ** (1.0 / r)
    volume_constant /= gamma(1.0 + float(reciprocal_sum))
    shell_constant = float(reciprocal_sum) * volume_constant
    shell_constant_direct = 1.0
    for r in DEGREES:
        shell_constant_direct *= gamma(1.0 + 1.0 / r) * factorial(r) ** (1.0 / r)
    shell_constant_direct /= gamma(float(reciprocal_sum))
    require(
        abs(shell_constant - shell_constant_direct) < 1e-12,
        "Gamma-constant identity failure",
    )
    print(f"N={N}")
    print(f"reciprocal_sum={reciprocal_sum} excess={reciprocal_sum - 1}")
    print(f"cumulative_volume_constant={volume_constant:.15f}")
    print(f"formal_shell_constant={shell_constant:.15f}")
    print(f"formal_shell_mean_at_N={shell_constant * N ** (1 / 24):.15f}")
    print()
    for q in (
        2, 4, 8, 16, 32,
        3, 9, 27,
        5, 25, 125,
        7, 49,
        11, 121,
        13, 169,
        17, 289,
        19, 361,
        33,
        99,
    ):
        show_profile(q)

    f3, _ = local_profile(3)
    f9, _ = local_profile(9)
    f11, _ = local_profile(11)
    f33, _ = local_profile(33)
    f99, _ = local_profile(99)
    require(
        f33[N % 33] == f3[N % 3] * f11[N % 11],
        "CRT factorization failure at 33",
    )
    require(
        f99[N % 99] == f9[N % 9] * f11[N % 11],
        "CRT factorization failure at 99",
    )
    print()
    print(f"CRT-check factor_33(N)=factor_3(N)*factor_11(N)={f33[N % 33]}")
    print(f"CRT-check factor_99(N)=factor_9(N)*factor_11(N)={f99[N % 99]}")
    print(
        "archimedean_times_9_11="
        f"{shell_constant * N ** (1 / 24) * float(f99[N % 99]):.15f}"
    )


if __name__ == "__main__":
    main()
