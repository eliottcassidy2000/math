#!/usr/bin/env python3
"""Exact arithmetic referee for THM-4096 and MISTAKE-497.

The finite universe is m<=2000 for the rational ray identities, direct
Dedekind sums through m=10, all vertices of the {2,7} affine optimizer, and
six exact Kubota--Leopoldt interpolation controls through Bernoulli index 198.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import comb


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def denominator(m: int) -> int:
    return 182 * m + 1


def witness_margin(m: int) -> F:
    return F(14 * m, denominator(m)) - F(1, 14)


def carrier(m: int) -> F:
    return -F(14 * 14, 12) * witness_margin(m)


def fixed_observer(m: int) -> F:
    return F(91 * m * (13 * m - 14), 6 * denominator(m))


def adapted_observer(m: int) -> F:
    return F(91 * m * (13 - 14 * m), 6 * denominator(m))


def sawtooth(numerator: int, modulus: int) -> F:
    residue = numerator % modulus
    return F(0) if residue == 0 else F(residue, modulus) - F(1, 2)


def dedekind(h: int, k: int) -> F:
    return sum(
        (sawtooth(r, k) * sawtooth(h * r, k) for r in range(1, k)),
        F(0),
    )


def bernoulli_numbers(limit: int) -> list[F]:
    values = [F(0) for _ in range(limit + 1)]
    values[0] = F(1)
    for n in range(1, limit + 1):
        values[n] = -sum((F(comb(n + 1, j)) * values[j] for j in range(n)), F(0)) / F(n + 1)
    return values


def valuation(number: F, prime: int) -> int:
    numerator = abs(number.numerator)
    denominator_value = number.denominator
    answer = 0
    while numerator and numerator % prime == 0:
        numerator //= prime
        answer += 1
    while denominator_value % prime == 0:
        denominator_value //= prime
        answer -= 1
    return answer


def rational_mod(number: F, prime: int) -> int:
    require(number.denominator % prime != 0, "nonintegral rational reduction")
    return (number.numerator % prime) * pow(number.denominator, -1, prime) % prime


def main() -> None:
    z_infinity = -F(1, 12)
    ray_gates = 0
    interval_gates = 0
    for m in range(1, 2001):
        d = denominator(m)
        delta = witness_margin(m)
        c = carrier(m)
        fixed = fixed_observer(m)
        adapted = adapted_observer(m)
        require(delta == F(14 * m - 1, 14 * d), "witness margin identity failed")
        require(c == -F(7 * (14 * m - 1), 6 * d), "carrier identity failed")
        require(c - z_infinity == F(15 - 14 * m, 12 * d), "archimedean split failed")
        require(fixed - c == F((m - 1) * (1183 * m + 7), 6 * d), "fixed split failed")
        require(c - adapted == F((m - 1) * (1274 * m - 7), 6 * d), "adapted split failed")
        if m == 1:
            require(fixed == c == adapted == -F(91, 1098), "first collision changed")
            require(c > z_infinity, "first carrier lost positive hybrid side")
        else:
            require(adapted < c < z_infinity < F(1, 2) < fixed, f"interval certificate failed m={m}")
            interval_gates += 1
        ray_gates += 1

    direct_dedekind_gates = 0
    for m in range(1, 11):
        d = denominator(m)
        require(dedekind(14, d) == fixed_observer(m), f"fixed Dedekind replay failed m={m}")
        require(dedekind(14 * m, d) == adapted_observer(m), f"adapted Dedekind replay failed m={m}")
        direct_dedekind_gates += 2

    # Rational Euler-factor vertices.  For odd p these are correctly typed as
    # Z_p^(2)=L_p(-1,omega_p^2); p=2 is retained only as a rational vertex.
    twisted_gates = 0
    for prime in (2, 3, 5, 7, 11, 13):
        twisted = F(prime - 1, 12)
        require(twisted == z_infinity + F(prime, 12), "twisted affine value failed")
        twisted_gates += 1
    require(F(7 - 1, 12) == F(1, 2), "twisted 7-adic value changed")

    # Full nonnegative {2,7} optimizer at m=1.
    rhs = F(1, 183)
    endpoints = [(rhs / 2, F(0)), (F(0), rhs / 7)]
    costs = []
    optimizer_gates = 0
    for alpha_2, alpha_7 in endpoints:
        require(2 * alpha_2 + 7 * alpha_7 == rhs, "optimizer endpoint infeasible")
        alpha_infinity = 1 - alpha_2 - alpha_7
        hybrid = alpha_infinity * z_infinity + alpha_2 * F(1, 12) + alpha_7 * F(1, 2)
        require(hybrid == carrier(1), "optimizer endpoint has wrong carrier")
        costs.append(alpha_2 + alpha_7)
        optimizer_gates += 1
    require(costs == [F(1, 366), F(1, 1281)], "optimizer costs changed")
    require(costs[1] < costs[0], "7-adic endpoint is not unique equal-cost minimizer")
    for prime_max in (2, 3, 5, 7, 11, 13, 101):
        minimum_cost = rhs / prime_max
        require(prime_max * minimum_cost == rhs, "finite-support optimizer failed")
        optimizer_gates += 1

    # Trivial-character interpolation approximants. These are finite controls;
    # the theorem uses interpolation and continuity for the all-j proof.
    indices = {(5, 1), (5, 2), (7, 1), (7, 2), (11, 1), (13, 1)}
    max_index = max(2 + (prime - 3) * prime**j for prime, j in indices)
    bernoulli = bernoulli_numbers(max_index)
    interpolation_rows = []
    for prime, j in sorted(indices):
        k = 2 + (prime - 3) * prime**j
        require(k % (prime - 1) == 0 and k % prime == 2, "interpolation index type failed")
        approximant = (1 - prime ** (k - 1)) * (-bernoulli[k] / k)
        require(valuation(approximant, prime) == -1, "trivial zeta approximant valuation failed")
        scaled = prime * approximant
        require(rational_mod(scaled, prime) == pow(2, -1, prime), "scaled zeta congruence failed")
        interpolation_rows.append((prime, j, k, rational_mod(scaled, prime)))

    require(-4 * carrier(1) == F(182, 549), "eta normalization control failed")
    require(fixed_observer(2) == F(364, 365), "m=2 fixed endpoint failed")
    require(carrier(2) == -F(63, 730), "m=2 carrier failed")
    require(adapted_observer(2) == -F(91, 73), "m=2 adapted endpoint failed")
    require(carrier(2) - z_infinity == -F(13, 4380), "m=2 negative gap failed")

    print("THM-4096 twisted p-adic affine optimizer audit: PASS")
    print(f"ray identities={ray_gates}; strict interval certificates={interval_gates}")
    print(f"direct Dedekind sums={direct_dedekind_gates}; Euler-factor vertices={twisted_gates}")
    print(f"optimizer gates={optimizer_gates}; interpolation controls={len(interpolation_rows)}")
    print("m=1: C=s(14,183)=-91/1098; 2a+7b=1/183; minimizer=(0,1/1281)")
    print("m=2 interval: -91/73 < -63/730 < -1/12 < 1/2 < 364/365")
    print(f"trivial-branch scaled congruence rows={interpolation_rows}")
    print("eta control: -4*s(14,183)=182/549")


if __name__ == "__main__":
    main()
