#!/usr/bin/env python3
"""Independent exact-arithmetic hostile audit for THM-4096.

This audit deliberately avoids the primary referee's sawtooth summation and
Bernoulli recurrence.  It reconstructs the LRC witness from folded residues,
computes Dedekind sums by the Euclidean reciprocity algorithm, computes
Bernoulli numbers by the Akiyama--Tanigawa transform, and certifies the affine
optimizer by exact exchange identities.

Finite audit universe: all 1 <= m <= 4096 for the phase and affine sign,
1 <= m <= 512 for both reciprocity computations, the complete {2,7}
optimizer segment through exact endpoint/exchange formulas, seven finite
trivial-character interpolation controls, and finite-prime optimizer controls
through 97.  The all-m and all-finite-prime conclusions use the displayed
cleared-denominator/exchange identities, not extrapolation from this universe.
"""

from __future__ import annotations

from fractions import Fraction as Q
from math import gcd


def require(condition: bool, message: str) -> None:
    """Optimization-safe assertion."""
    if not condition:
        raise RuntimeError(message)


def folded_distance(numerator: int, denominator: int) -> Q:
    residue = numerator % denominator
    return Q(min(residue, denominator - residue), denominator)


def deep_well_phase(m: int) -> tuple[Q, Q, tuple[int, ...]]:
    """Reconstruct the THM-2057 phase and its active speeds from residues."""
    denominator = 182 * m + 1
    time = Q(14 * m, denominator)
    speeds = tuple(range(1, 13)) + (182 * m,)
    distances = tuple(folded_distance(speed * time.numerator, denominator) for speed in speeds)
    minimum = min(distances)
    owners = tuple(speed for speed, distance in zip(speeds, distances) if distance == minimum)
    return time, minimum, owners


def dedekind_by_reciprocity(h: int, k: int) -> Q:
    """Exact Dedekind sum via Euclid and reciprocity, not sawtooth enumeration."""
    require(k > 0, "Dedekind modulus must be positive")
    h %= k
    require(gcd(h, k) == 1, "Dedekind arguments must be coprime")
    if h == 0:
        require(k == 1, "zero numerator is coprime only at modulus one")
        return Q(0)
    if h == 1:
        return Q((k - 1) * (k - 2), 12 * k)
    reciprocal = -Q(1, 4) + (Q(h, k) + Q(k, h) + Q(1, h * k)) / 12
    return reciprocal - dedekind_by_reciprocity(k % h, h)


def bernoulli_table(limit: int) -> list[Q]:
    """Akiyama--Tanigawa transform (the B_1=+1/2 convention)."""
    work: list[Q] = []
    answer: list[Q] = []
    for m in range(limit + 1):
        work.append(Q(1, m + 1))
        for j in range(m, 0, -1):
            work[j - 1] = j * (work[j - 1] - work[j])
        answer.append(work[0])
    return answer


def valuation(number: Q, prime: int) -> int:
    require(number != 0, "valuation of zero was not requested")
    numerator = abs(number.numerator)
    denominator = number.denominator
    result = 0
    while numerator % prime == 0:
        numerator //= prime
        result += 1
    while denominator % prime == 0:
        denominator //= prime
        result -= 1
    return result


def residue_of_integral_rational(number: Q, prime: int) -> int:
    require(number.denominator % prime != 0, "rational is not p-integral")
    return number.numerator * pow(number.denominator, -1, prime) % prime


def affine_hybrid(z_infinity: Q, coefficients: dict[int, Q]) -> Q:
    """Rational affine hybrid; p=2 denotes the Euler-factor vertex (7)."""
    alpha_infinity = 1 - sum(coefficients.values(), Q(0))
    require(alpha_infinity >= 0, "negative archimedean coefficient")
    return alpha_infinity * z_infinity + sum(
        (alpha * Q(prime - 1, 12) for prime, alpha in coefficients.items()),
        Q(0),
    )


def main() -> None:
    z_infinity = -Q(1, 12)
    phase_rows = 0
    negative_side_rows = 0

    for m in range(1, 4097):
        denominator = 182 * m + 1
        time, lonely_value, owners = deep_well_phase(m)
        require(time == Q(14 * m, denominator), f"phase reconstruction failed at m={m}")
        require(lonely_value == time, f"phase value failed at m={m}")
        require(owners == (1, 182 * m), f"active owner pair changed at m={m}")

        margin = lonely_value - Q(1, 14)
        carrier = -Q(49, 3) * margin
        # These two cleared-denominator identities are the all-m algebraic
        # certificate for the margin and its archimedean displacement.
        require(14 * denominator * margin == 14 * m - 1, f"margin polynomial failed at m={m}")
        affine_rhs = 12 * (carrier - z_infinity)
        require(denominator * affine_rhs == 15 - 14 * m, f"affine polynomial failed at m={m}")
        if m == 1:
            require(affine_rhs == Q(1, 183), "first affine right side changed")
        else:
            require(affine_rhs < 0, f"nonnegative hybrid obstruction failed at m={m}")
            negative_side_rows += 1
        phase_rows += 1

    reciprocity_rows = 0
    first_fixed = Q(0)
    first_adapted = Q(0)
    for m in range(1, 513):
        denominator = 182 * m + 1
        fixed = dedekind_by_reciprocity(14, denominator)
        adapted = dedekind_by_reciprocity(14 * m, denominator)
        expected_fixed = Q(91 * m * (13 * m - 14), 6 * denominator)
        expected_adapted = Q(91 * m * (13 - 14 * m), 6 * denominator)
        require(fixed == expected_fixed, f"fixed reciprocity formula failed at m={m}")
        require(adapted == expected_adapted, f"adapted reciprocity formula failed at m={m}")

        _, lonely_value, _ = deep_well_phase(m)
        carrier = -Q(49, 3) * (lonely_value - Q(1, 14))
        require(
            fixed - carrier == Q((m - 1) * (1183 * m + 7), 6 * denominator),
            f"fixed/carrier split failed at m={m}",
        )
        require(
            carrier - adapted == Q((m - 1) * (1274 * m - 7), 6 * denominator),
            f"carrier/adapted split failed at m={m}",
        )
        if m == 1:
            first_fixed, first_adapted = fixed, adapted
            require(fixed == carrier == adapted, "triple collision failed")
        else:
            require(adapted < carrier < z_infinity < Q(1, 2) < fixed, f"observer order failed at m={m}")
        reciprocity_rows += 2

    # For odd p, equation (1.1) with chi=omega_p^2 and n=2 has
    # chi*omega_p^{-2}=1, hence Z_p^(2)=(1-p) zeta(-1).
    bernoulli = bernoulli_table(198)
    require(bernoulli[2] == Q(1, 6), "independent B_2 computation failed")
    complex_zeta_minus_one = -bernoulli[2] / 2
    require(complex_zeta_minus_one == z_infinity, "zeta(-1) normalization failed")
    twisted_primes = (3, 5, 7, 11, 13, 17)
    twisted_values: dict[int, Q] = {}
    for prime in twisted_primes:
        value = (1 - prime) * complex_zeta_minus_one
        require(value == Q(prime - 1, 12), f"twisted interpolation algebra failed at p={prime}")
        require(value - z_infinity == Q(prime, 12), f"twisted affine displacement failed at p={prime}")
        twisted_values[prime] = value
    require(twisted_values[7] == Q(1, 2), "7-adic twisted value changed")
    two_euler_vertex = Q(2 - 1, 12)
    require(two_euler_vertex == Q(1, 12), "p=2 rational Euler-factor vertex changed")
    require(two_euler_vertex - z_infinity == Q(2, 12), "p=2 affine displacement changed")

    # Complete {2,7} feasible segment at m=1.  The exchange identity
    # cost-rhs/7=(5/7)alpha_2 proves uniqueness of the equal-cost minimizer.
    rhs = Q(1, 183)
    optimizer_samples = 0
    costs: list[Q] = []
    for numerator in range(0, 101):
        alpha_7 = Q(numerator, 100) * rhs / 7
        alpha_2 = (rhs - 7 * alpha_7) / 2
        require(alpha_2 >= 0 and alpha_7 >= 0, "optimizer sample left feasible segment")
        require(2 * alpha_2 + 7 * alpha_7 == rhs, "optimizer equation failed")
        require(
            affine_hybrid(z_infinity, {2: alpha_2, 7: alpha_7})
            == -Q(91, 1098),
            "hybrid failed to recover first carrier",
        )
        cost = alpha_2 + alpha_7
        require(cost - rhs / 7 == Q(5, 7) * alpha_2, "optimizer exchange identity failed")
        costs.append(cost)
        optimizer_samples += 1
    require(costs[0] == Q(1, 366), "2-adic endpoint cost changed")
    require(costs[-1] == Q(1, 1281), "7-adic endpoint cost changed")
    require(all(left > right for left, right in zip(costs, costs[1:])), "segment cost is not strictly decreasing")

    strict_costs: list[Q] = []
    for scale in (2, 3, 5, 11, 101, 1009):
        alpha_2 = rhs / (2 * scale)
        alpha_7 = (rhs - 2 * alpha_2) / 7
        require(alpha_2 > 0 and alpha_7 > 0, "strict optimizer control hit boundary")
        cost = alpha_2 + alpha_7
        require(cost > rhs / 7, "strict-positive optimizer attained forbidden endpoint")
        strict_costs.append(cost)
    require(all(left > right for left, right in zip(strict_costs, strict_costs[1:])), "strict costs do not approach infimum")

    # General finite-prime equal-cost exchange: cost-rhs/q_max is a sum of
    # nonnegative terms alpha_p(1-p/q_max), with equality only at q_max.
    finite_primes = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97)
    finite_optimizer_rows = 0
    previous_single_cost: Q | None = None
    for end in range(1, len(finite_primes) + 1):
        bank = finite_primes[:end]
        largest = bank[-1]
        single_cost = rhs / largest
        if previous_single_cost is not None:
            require(single_cost < previous_single_cost, "unrestricted prime controls did not improve")
        previous_single_cost = single_cost
        weights = {prime: Q(index + 1) for index, prime in enumerate(bank)}
        normalizer = sum((prime * weight for prime, weight in weights.items()), Q(0))
        coefficients = {prime: rhs * weight / normalizer for prime, weight in weights.items()}
        require(sum((prime * alpha for prime, alpha in coefficients.items()), Q(0)) == rhs, "finite-bank feasibility failed")
        cost = sum(coefficients.values(), Q(0))
        exchange = sum((alpha * (1 - Q(prime, largest)) for prime, alpha in coefficients.items()), Q(0))
        require(cost - single_cost == exchange, "finite-bank exchange identity failed")
        require(cost >= single_cost, "finite-bank maximum-prime optimizer failed")
        require(affine_hybrid(z_infinity, coefficients) == -Q(91, 1098), "finite-bank hybrid value failed")
        finite_optimizer_rows += 1

    # Odd-prime trivial-character branch: k_j is divisible by p-1 and tends
    # p-adically to 2.  (At p=2 the valuation is -2, so it is intentionally
    # outside this v_p=-1 correction.)
    interpolation_indices = ((3, 1), (5, 1), (5, 2), (7, 1), (7, 2), (11, 1), (13, 1))
    interpolation_rows: list[tuple[int, int, int, int]] = []
    for prime, depth in interpolation_indices:
        k = 2 + (prime - 3) * prime**depth
        require(k % (prime - 1) == 0, "trivial interpolation parity failed")
        require(k % prime == 2 % prime, "interpolation did not approach 2")
        approximant = (1 - prime ** (k - 1)) * (-bernoulli[k] / k)
        require(valuation(approximant, prime) == -1, "odd-prime trivial approximant valuation failed")
        scaled = prime * approximant
        residue = residue_of_integral_rational(scaled, prime)
        require(residue == pow(2, -1, prime), "odd-prime scaled residue failed")
        interpolation_rows.append((prime, depth, k, residue))

    seven_approximation = (1 - 7**29) * (-bernoulli[30] / 30)
    require(valuation(seven_approximation - Q(1, 2), 7) == -1, "old 7-adic mistyping was not rejected")
    eta_signature = -4 * first_fixed
    require(first_fixed == first_adapted == -Q(91, 1098), "first Dedekind value changed")
    require(eta_signature == Q(182, 549), "eta normalization failed")
    require(eta_signature != first_fixed, "hostile missing-factor control did not separate")

    _, second_lonely_value, _ = deep_well_phase(2)
    second_carrier = -Q(49, 3) * (second_lonely_value - Q(1, 14))
    second_fixed = dedekind_by_reciprocity(14, 365)
    second_adapted = dedekind_by_reciprocity(28, 365)
    require(second_adapted == -Q(91, 73), "m=2 adapted endpoint changed")
    require(second_carrier == -Q(63, 730), "m=2 carrier changed")
    require(second_carrier - z_infinity == -Q(13, 4380), "m=2 obstruction gap changed")
    require(second_fixed == Q(364, 365), "m=2 fixed endpoint changed")

    print("THM-4096 independent hostile audit: PASS")
    print(f"finite universe: phase/affine m=1..{phase_rows}; negative-side rows={negative_side_rows}")
    print(f"Euclidean-reciprocity Dedekind evaluations={reciprocity_rows}")
    print("active phase owners for every audited m: {1,182m}; value=14m/(182m+1)")
    print(f"rational vertices: p=2 Euler-factor -> {two_euler_vertex}; p=7 twisted -> {twisted_values[7]}")
    print(f"{{2,7}} optimizer samples={optimizer_samples}; endpoint costs=1/366,1/1281; unique minimizer=(0,1/1281)")
    print(f"finite-prime exchange rows={finite_optimizer_rows}; last single-prime cost={previous_single_cost}")
    print(f"odd-prime trivial-branch controls={interpolation_rows}")
    print("hostiles rejected: trivial zeta_7(-1)=1/2; eta_sig=s; nonnegative hybrid at m=2")
    print("m=1: s(14,183)=s(14,183 adapted)=C=-91/1098; eta_sig=182/549")
    print("m=2: -91/73 < -63/730 < -1/12 < 1/2 < 364/365; C-z_infinity=-13/4380")


if __name__ == "__main__":
    main()
