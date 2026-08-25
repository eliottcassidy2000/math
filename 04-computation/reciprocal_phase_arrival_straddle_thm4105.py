#!/usr/bin/env python3
"""Exact hostile audit for reciprocal exponentiation and LRC arrival.

This is a research scout.  It uses only integer and Fraction arithmetic and
keeps every truth gate active under ``python -O``.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations, product
from math import gcd


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def egcd(a: int, b: int) -> tuple[int, int, int]:
    if b == 0:
        return abs(a), 1 if a >= 0 else -1, 0
    g, x1, y1 = egcd(b, a % b)
    return g, y1, x1 - (a // b) * y1


def bezout_many(values: tuple[int, ...]) -> tuple[int, tuple[int, ...]]:
    g = 0
    coeffs: list[int] = []
    for value in values:
        new_g, scale_old, scale_new = egcd(g, value)
        coeffs = [scale_old * c for c in coeffs] + [scale_new]
        g = new_g
    require(sum(c * v for c, v in zip(coeffs, values)) == g, "bad Bezout")
    return g, tuple(coeffs)


def circle_norm(x: Fraction) -> Fraction:
    residue = x % 1
    return min(residue, 1 - residue)


def pair_defect_mod(speeds: tuple[int, ...], phases: tuple[int, ...], q: int) -> bool:
    return all(
        (speeds[j] * phases[i] - speeds[i] * phases[j]) % q == 0
        for i in range(len(speeds))
        for j in range(i + 1, len(speeds))
    )


def modular_descent_audit() -> tuple[int, int, int]:
    states = 0
    compatible = 0
    reconstructed = 0
    for size in range(2, 5):
        for speeds in combinations(range(1, 8), size):
            if gcd(*speeds) != 1:
                continue
            g, coeffs = bezout_many(speeds)
            require(g == 1, "primitive tuple did not have Bezout gcd one")
            for q in range(2, 9):
                for phases in product(range(q), repeat=size):
                    states += 1
                    if not pair_defect_mod(speeds, phases, q):
                        continue
                    compatible += 1
                    time = sum(c * h for c, h in zip(coeffs, phases)) % q
                    require(
                        all((speed * time - phase) % q == 0 for speed, phase in zip(speeds, phases)),
                        f"descent failure: {speeds=} {q=} {phases=}",
                    )
                    reconstructed += 1

    # Primitivity is load-bearing.  These reciprocal equations hold, but the
    # phase vector is not on the (2,4)-orbit modulo two.
    speeds = (2, 4)
    phases = (0, 1)
    q = 2
    require(pair_defect_mod(speeds, phases, q), "nonprimitive hostile lost compatibility")
    require(
        not any(
            all((speed * time - phase) % q == 0 for speed, phase in zip(speeds, phases))
            for time in range(q)
        ),
        "nonprimitive hostile unexpectedly arrived",
    )
    require(
        not pair_defect_mod((1, 2), phases, q),
        "primitive normalization should expose the hostile defect",
    )
    return states, compatible, reconstructed


def quantitative_arrival_audit() -> tuple[int, Fraction]:
    gates = 0
    largest_tariff = Fraction(0)
    for size in range(2, 5):
        for speeds in combinations(range(1, 7), size):
            if gcd(*speeds) != 1:
                continue
            _, coeffs = bezout_many(speeds)
            tariff = sum(abs(c) for c in coeffs)
            largest_tariff = max(largest_tariff, Fraction(tariff))
            vmax = max(speeds)
            for q in (5, 7):
                for phases_mod_q in product(range(q), repeat=size):
                    phases = tuple(Fraction(h, q) for h in phases_mod_q)
                    defects = {
                        (i, j): circle_norm(speeds[j] * phases[i] - speeds[i] * phases[j])
                        for i in range(size)
                        for j in range(size)
                    }
                    defect = max(defects.values(), default=Fraction(0))
                    time = sum(c * phase for c, phase in zip(coeffs, phases)) % 1
                    arrival_at_reconstruction = max(
                        circle_norm(phase - speed * time)
                        for speed, phase in zip(speeds, phases)
                    )
                    require(
                        arrival_at_reconstruction <= tariff * defect,
                        f"upper arrival tariff failed: {speeds=} {phases_mod_q=}",
                    )
                    for time_num in range(q):
                        test_time = Fraction(time_num, q)
                        error = max(
                            circle_norm(phase - speed * test_time)
                            for speed, phase in zip(speeds, phases)
                        )
                        require(
                            defect <= 2 * vmax * error,
                            f"lower arrival tariff failed: {speeds=} {phases_mod_q=}",
                        )
                        gates += 1
    return gates, largest_tariff


def reciprocal_commutator(a: int, b: int) -> Fraction:
    return Fraction(a**b, b**a)


def natural_power_order_audit() -> tuple[tuple[tuple[int, int], ...], int, int]:
    ties: list[tuple[int, int]] = []
    for a in range(1, 41):
        for b in range(a + 1, 41):
            if a**b == b**a:
                ties.append((a, b))
    require(ties == [(2, 4)], f"unexpected natural reciprocal-power ties: {ties}")

    strict_triples = 0
    directed_cycles = 0
    for a, b, c in combinations(range(1, 21), 3):
        ab = (a**b > b**a) - (a**b < b**a)
        bc = (b**c > c**b) - (b**c < c**b)
        ca = (c**a > a**c) - (c**a < a**c)
        if 0 in (ab, bc, ca):
            continue
        strict_triples += 1
        if ab == bc == ca:
            directed_cycles += 1
    require(directed_cycles == 0, "reciprocal-power orientation formed a directed cycle")

    triangle_gates = 0
    for a, b, c in combinations(range(2, 10), 3):
        weighted_curl = (
            reciprocal_commutator(a, b) ** c
            * reciprocal_commutator(b, c) ** a
            * reciprocal_commutator(c, a) ** b
        )
        require(weighted_curl == 1, f"weighted triangle law failed: {(a, b, c)}")
        triangle_gates += 1
    return tuple(ties), strict_triples, triangle_gates


def straddle_factorization_audit() -> tuple[int, tuple[int, int, int, int, Fraction, Fraction]]:
    gates = 0
    for u in range(1, 17):
        for w in range(1, 17):
            if u == w:
                continue
            q = u + w
            for alpha in range(0, 19):
                for beta in range(0, 19):
                    p = alpha + beta
                    determinant = u * beta - w * alpha
                    if determinant <= 0 or 2 * determinant > q:
                        continue
                    time = Fraction(p, q)
                    margin = Fraction(determinant, q)
                    require(u * time - alpha == margin, "positive straddle lift failed")
                    require(w * time - beta == -margin, "negative straddle lift failed")

                    # (u^beta / w^alpha)^q
                    #   = (u^w / w^u)^p * (uw)^determinant.
                    lhs = Fraction(u**beta, w**alpha) ** q
                    rhs = reciprocal_commutator(u, w) ** p * Fraction(u * w) ** determinant
                    require(lhs == rhs, "reciprocal exponentiation factorization failed")
                    require(beta * q == w * p + determinant, "u-exponent ledger failed")
                    require(alpha * q == u * p - determinant, "w-exponent ledger failed")
                    gates += 1

    ap13 = (1, 13, 0, 1, Fraction(1, 14), Fraction(1, 14))
    u, w, alpha, beta, time, margin = ap13
    require(u * time - alpha == margin, "AP13 lower owner failed")
    require(w * time - beta == -margin, "AP13 upper owner failed")
    return gates, ap13


def abc_reciprocal_power_carrier_audit() -> tuple[int, int, int]:
    coprime_pairs = 0
    mixed_gcd_one = 0
    odd_gcd_two = 0
    for a in range(2, 19):
        for b in range(a + 1, 19):
            if gcd(a, b) != 1:
                continue
            x = a**b
            y = b**a
            require(x != y, "coprime off-diagonal reciprocal powers tied")
            height = max(x, y)
            gap = abs(x - y)
            total = x + y
            require(gcd(x, y) == 1, "primitive powered terms not coprime")
            require(height >= 2 ** max(a, b), "height/log carrier bound failed")
            require(a * b <= max(a, b) ** 2, "base product bound failed")
            common = gcd(total, gap)
            if (a * b) % 2 == 0:
                require(common == 1, "mixed-parity sum/gap gcd should be one")
                mixed_gcd_one += 1
            else:
                require(common == 2, "odd sum/gap gcd should be two")
                odd_gcd_two += 1
            coprime_pairs += 1

    # Every oriented LRC determinant packet becomes a primitive ABC triple
    # only after its actual product gcd is removed.
    straddle_packets = 0
    for u in range(1, 13):
        for w in range(1, 13):
            for alpha in range(1, 13):
                for beta in range(1, 13):
                    x = u * beta
                    y = w * alpha
                    if x == y:
                        continue
                    d = abs(x - y)
                    scale = gcd(x, y)
                    xx, yy, dd = x // scale, y // scale, d // scale
                    require(gcd(xx, yy) == gcd(xx, dd) == gcd(yy, dd) == 1, "bad ABC normalization")
                    require(max(xx, yy) == min(xx, yy) + dd, "bad normalized additive packet")
                    straddle_packets += 1
    return coprime_pairs, mixed_gcd_one + odd_gcd_two, straddle_packets


def main() -> None:
    states, compatible, reconstructed = modular_descent_audit()
    quantitative_gates, largest_tariff = quantitative_arrival_audit()
    ties, strict_triples, triangle_gates = natural_power_order_audit()
    straddle_gates, ap13 = straddle_factorization_audit()
    power_pairs, parity_gates, abc_packets = abc_reciprocal_power_carrier_audit()

    print("RECIPROCAL EXPONENTIATION / LRC ARRIVAL EXACT AUDIT")
    print(f"modular states={states} compatible={compatible} reconstructed={reconstructed}")
    print("nonprimitive hostile: v=(2,4), q=2, phase=(0,1) compatible but not physical")
    print(f"quantitative pointwise gates={quantitative_gates} largest_sample_Bezout_l1={largest_tariff}")
    print(f"natural off-diagonal ties={ties}")
    print(f"strict natural triples={strict_triples} directed_cycles=0")
    print(f"weighted reciprocal triangle gates={triangle_gates}")
    print(f"LRC straddle factorization gates={straddle_gates}")
    print(f"AP13 straddle={ap13}")
    print(f"coprime reciprocal-power pairs={power_pairs} parity_gcd_gates={parity_gates}")
    print(f"normalized LRC determinant ABC packets={abc_packets}")
    print("PASS")


if __name__ == "__main__":
    main()


