#!/usr/bin/env python3
"""Exact finite controls for the fixed-modulus average-order theorem.

For q=99 this counts all admissible tuples with total at most X, separated by
the total residue modulo q.  It compares those exact counts with the rigorous
leading constants predicted by anisotropic Riemann sums.  This is an
average-order statement; it makes no pointwise circle-method assertion.
"""

from __future__ import annotations

from fractions import Fraction
from hashlib import sha256
from math import comb, factorial, gamma, gcd, isqrt


Q = 99
DEGREES = (2, 4, 6, 8)
STARTS = (2, 3, 5, 7)
TARGET_RESIDUE = 53
WORST_RESIDUE = 86


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def factor_prime_powers(q: int) -> tuple[tuple[int, int], ...]:
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
    return tuple(out)


def floor_log_p(k: int, p: int) -> int:
    e = 0
    while k >= p:
        k //= p
        e += 1
    return e


def period(q: int, k: int) -> int:
    answer = 1
    for p, a in factor_prime_powers(q):
        component = p ** (a + floor_log_p(k, p))
        answer = answer // gcd(answer, component) * component
    return answer


def cyclic_convolution(left: tuple[int, ...], right: tuple[int, ...]) -> tuple[int, ...]:
    answer = [0] * len(left)
    for i, a in enumerate(left):
        if a:
            for j, b in enumerate(right):
                if b:
                    answer[(i + j) % len(left)] += a * b
    return tuple(answer)


def local_factors(q: int) -> tuple[Fraction, ...]:
    total = (1,) + (0,) * (q - 1)
    universe = 1
    for k in DEGREES:
        orbit = period(q, k)
        hist = [0] * q
        for x in range(orbit):
            hist[comb(x, k) % q] += 1
        total = cyclic_convolution(total, tuple(hist))
        universe *= orbit
    return tuple(Fraction(q * count, universe) for count in total)


def binomial_values(k: int, start: int, limit: int) -> tuple[int, ...]:
    values = []
    x = start
    while comb(x, k) <= limit:
        values.append(comb(x, k))
        x += 1
    return tuple(values)


def triangular_top(limit: int) -> int:
    if limit < 1:
        return 1
    top = (1 + isqrt(1 + 8 * limit)) // 2
    while comb(top, 2) > limit:
        top -= 1
    while comb(top + 1, 2) <= limit:
        top += 1
    return top


def triangular_prefix_data(q: int) -> tuple[tuple[int, ...], tuple[tuple[int, ...], ...]]:
    orbit = period(q, 2)
    sequence = tuple(comb(w, 2) % q for w in range(2, 2 + orbit))
    prefixes = [[0] * q]
    for residue in sequence:
        row = prefixes[-1].copy()
        row[residue] += 1
        prefixes.append(row)
    return tuple(prefixes[-1]), tuple(tuple(row) for row in prefixes)


def cumulative_residue_counts(limit: int) -> tuple[int, ...]:
    c4 = binomial_values(4, 3, limit - 1)
    c6 = binomial_values(6, 5, limit - 1)
    c8 = binomial_values(8, 7, limit - 1)
    orbit = period(Q, 2)
    cycle_hist, prefixes = triangular_prefix_data(Q)
    answer = [0] * Q
    for t8 in c8:
        for t6 in c6:
            partial = t8 + t6
            if partial >= limit:
                break
            for t4 in c4:
                base = partial + t4
                if base >= limit:
                    break
                top = triangular_top(limit - base)
                number_of_w = top - 1  # w=2,...,top
                cycles, remainder = divmod(number_of_w, orbit)
                prefix = prefixes[remainder]
                shift = base % Q
                for triangular_residue in range(Q):
                    multiplicity = cycles * cycle_hist[triangular_residue] + prefix[triangular_residue]
                    answer[(shift + triangular_residue) % Q] += multiplicity
    return tuple(answer)


def constants() -> tuple[Fraction, float]:
    reciprocal_sum = sum((Fraction(1, k) for k in DEGREES), Fraction(0))
    numerator = 1.0
    for k in DEGREES:
        numerator *= gamma(1 + 1 / k) * factorial(k) ** (1 / k)
    return reciprocal_sum, numerator / gamma(1 + float(reciprocal_sum))


def main() -> None:
    sigma = local_factors(Q)
    require(sigma[TARGET_RESIDUE] == Fraction(544, 1089), "target sigma")
    require(sigma[WORST_RESIDUE] == Fraction(496, 1089), "worst sigma")
    require(min(sigma) == sigma[WORST_RESIDUE], "worst residue")
    reciprocal_sum, volume = constants()
    require(reciprocal_sum == Fraction(25, 24), "critical exponent")

    rows = []
    for limit in (10_000, 100_000, 1_000_000, 10_000_000):
        counts = cumulative_residue_counts(limit)
        total = sum(counts)
        predicted_total = volume * limit ** float(reciprocal_sum)
        predicted_target = predicted_total * float(sigma[TARGET_RESIDUE]) / Q
        predicted_worst = predicted_total * float(sigma[WORST_RESIDUE]) / Q
        rows.append(
            (
                limit,
                total,
                counts[TARGET_RESIDUE],
                counts[WORST_RESIDUE],
                total / predicted_total,
                counts[TARGET_RESIDUE] / predicted_target,
                counts[WORST_RESIDUE] / predicted_worst,
            )
        )

    payload = "\n".join(
        f"{limit}:{total}:{target}:{worst}"
        for limit, total, target, worst, _, _, _ in rows
    )
    print("SUN_2468_FIXED_MODULUS_AVERAGE_ORDER_CONTROLS")
    print(f"q={Q} sigma_target={sigma[TARGET_RESIDUE]} sigma_worst={sigma[WORST_RESIDUE]}")
    print(f"reciprocal_sum={reciprocal_sum} volume_constant={volume:.15f}")
    for limit, total, target, worst, ratio_total, ratio_target, ratio_worst in rows:
        print(
            f"X={limit} total={total} residue_{TARGET_RESIDUE}={target} "
            f"residue_{WORST_RESIDUE}={worst} ratio_total={ratio_total:.12f} "
            f"ratio_target={ratio_target:.12f} ratio_worst={ratio_worst:.12f}"
        )
    print(f"semantic_sha256={sha256(payload.encode('ascii')).hexdigest()}")
    print("theorem_scope=FIXED_Q_CUMULATIVE_RESIDUE_CLASS_ASYMPTOTIC")
    print("pointwise_inference=NONE")
    print("PASS")


if __name__ == "__main__":
    main()
