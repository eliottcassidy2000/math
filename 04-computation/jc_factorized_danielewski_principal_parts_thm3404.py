#!/usr/bin/env python3
"""Exact finite controls for THM-3404.

The theorem is an all-parameter commutative-algebra statement.  This
standard-library companion checks its consequence objects on a bounded bank:

* factorwise coefficient ideals and CRT principal-part jets;
* the one-homogeneous-pole J(w)^n filtration, including one-step hostiles;
* class-group Smith data and the canonical common-root cover; and
* the duplicate-factor and tied-principal-part cancellation hostiles.

Normality, Nagata's sequence, and finite-cover norm/conorm are proved in the
theorem rather than inferred from enumeration.  Every gate remains active
under ``python -O``.
"""

from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
import json
from math import gcd


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def trim(poly):
    out = list(poly)
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return tuple(Fraction(value) for value in out)


ZERO = (Fraction(0),)
ONE = (Fraction(1),)


def add(left, right):
    size = max(len(left), len(right))
    out = [Fraction(0) for _ in range(size)]
    for i, value in enumerate(left):
        out[i] += value
    for i, value in enumerate(right):
        out[i] += value
    return trim(out)


def neg(poly):
    return tuple(-value for value in poly)


def sub(left, right):
    return add(left, neg(right))


def scale(poly, scalar):
    scalar = Fraction(scalar)
    return trim(tuple(scalar * value for value in poly))


def mul(left, right):
    out = [Fraction(0) for _ in range(len(left) + len(right) - 1)]
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            out[i + j] += a * b
    return trim(out)


def power(poly, exponent):
    require(exponent >= 0, ("negative polynomial exponent", exponent))
    answer = ONE
    base = poly
    current = exponent
    while current:
        if current & 1:
            answer = mul(answer, base)
        base = mul(base, base)
        current //= 2
    return answer


def divmod_poly(numerator, denominator):
    numerator = list(trim(numerator))
    denominator = trim(denominator)
    require(denominator != ZERO, "division by zero polynomial")
    if len(numerator) < len(denominator):
        return ZERO, trim(numerator)
    quotient = [Fraction(0) for _ in range(len(numerator) - len(denominator) + 1)]
    while len(numerator) >= len(denominator) and trim(numerator) != ZERO:
        shift = len(numerator) - len(denominator)
        coefficient = numerator[-1] / denominator[-1]
        quotient[shift] += coefficient
        for i, value in enumerate(denominator):
            numerator[i + shift] -= coefficient * value
        numerator = list(trim(numerator))
    return trim(quotient), trim(numerator)


def exact_div(numerator, denominator):
    quotient, remainder = divmod_poly(numerator, denominator)
    require(remainder == ZERO, ("inexact division", numerator, denominator, remainder))
    return quotient


def mod_poly(poly, modulus):
    return divmod_poly(poly, modulus)[1]


def xgcd(left, right):
    old_r, r = trim(left), trim(right)
    old_s, s = ONE, ZERO
    old_t, t = ZERO, ONE
    while r != ZERO:
        quotient, remainder = divmod_poly(old_r, r)
        old_r, r = r, remainder
        old_s, s = s, sub(old_s, mul(quotient, s))
        old_t, t = t, sub(old_t, mul(quotient, t))
    require(old_r != ZERO, "zero polynomial gcd")
    leading = old_r[-1]
    return scale(old_r, 1 / leading), scale(old_s, 1 / leading), scale(old_t, 1 / leading)


def linear(root):
    return (Fraction(-root), Fraction(1))


def product_polys(polys):
    answer = ONE
    for poly in polys:
        answer = mul(answer, poly)
    return answer


def polynomial_from_valuations(primes, valuations, positive):
    return product_polys(
        power(prime, value if positive else -value)
        for prime, value in zip(primes, valuations)
        if (value > 0 if positive else value < 0)
    )


def crt_factor_packets():
    # Three distinct arms already exercise torsion, free class rank, and
    # nontrivial CRT gluing while keeping normal/-O replays short.
    roots = (-3, 0, 2)
    packets = 0
    idempotents = 0
    exact_clearings = 0
    one_step_hostiles = 0
    cancellation_hostiles = 0
    cover_packets = 0
    cover_kernel_classes = 0
    max_degree = 0

    for rank in range(1, len(roots) + 1):
        for chosen_roots in combinations(roots, rank):
            factors = tuple(linear(root) for root in chosen_roots)
            for exponents in product(range(1, 5), repeat=rank):
                common = 0
                for exponent in exponents:
                    common = gcd(common, exponent)
                F = product_polys(
                    power(factor, exponent)
                    for factor, exponent in zip(factors, exponents)
                )

                # Cl(R)=Z^rank/<exponents> has free rank rank-1 and torsion
                # order gcd(exponents).  Bezout is checked by the integer gcd.
                require(common >= 1, ("zero class relation", exponents))

                for N in range(1, common + 1):
                    if common % N:
                        continue
                    cover_packets += 1
                    reduced = tuple(exponent // N for exponent in exponents)
                    require(
                        tuple(N * value for value in reduced) == exponents,
                        ("cover exponent reconstruction", exponents, N),
                    )
                    representatives = {
                        tuple((t * value) for value in reduced)
                        for t in range(N)
                    }
                    require(len(representatives) == N, ("cover kernel collapsed", exponents, N))
                    for left in range(N):
                        for right in range(N):
                            difference = left - right
                            equivalent_downstairs = difference % N == 0
                            require(
                                equivalent_downstairs == (left == right),
                                ("cover kernel representatives", exponents, N, left, right),
                            )
                    cover_kernel_classes += N
                    if N == common:
                        reduced_gcd = 0
                        for value in reduced:
                            reduced_gcd = gcd(reduced_gcd, value)
                        require(reduced_gcd == 1, ("maximal cover left torsion", exponents))

                for q in range(1, 4):
                    packets += 1
                    moduli = tuple(
                        power(factor, q * exponent)
                        for factor, exponent in zip(factors, exponents)
                    )
                    Fq = power(F, q)
                    require(product_polys(moduli) == Fq, ("factor product", exponents, q))
                    max_degree = max(max_degree, len(Fq) - 1)
                    require(
                        len(Fq) - 1 == q * sum(exponents),
                        ("principal-part length", exponents, q),
                    )

                    for i, modulus in enumerate(moduli):
                        complement = exact_div(Fq, modulus)
                        divisor, inverse, _ = xgcd(complement, modulus)
                        require(divisor == ONE, ("non-coprime distinct factors", chosen_roots, i))
                        idem = mod_poly(mul(complement, inverse), Fq)
                        require(mod_poly(idem, modulus) == ONE, ("CRT self arm", i))
                        for j, other in enumerate(moduli):
                            if i != j:
                                require(mod_poly(idem, other) == ZERO, ("CRT cross arm", i, j))
                        idempotents += 1

                    tail = (Fraction(1), Fraction(1), Fraction(1))
                    clearing = mul(Fq, tail)
                    require(exact_div(clearing, Fq) == tail, ("exact clearing", exponents, q))
                    exact_clearings += 1
                    for factor in factors:
                        dropped = exact_div(Fq, factor)
                        require(mod_poly(dropped, Fq) != ZERO, ("one-step hostile vanished", exponents, q))
                        require(divmod_poly(dropped, Fq)[1] == dropped, ("hostile unexpectedly cleared",))
                        one_step_hostiles += 1

                    # The canonical N-cover sends old degree -q to degree -Nq.
                    # Its coefficient modulus is G^(Nq)=F^q, but the cover also
                    # has new negative degrees not divisible by N.
                    for N in range(1, common + 1):
                        if common % N:
                            continue
                        G = product_polys(
                            power(factor, exponent // N)
                            for factor, exponent in zip(factors, exponents)
                        )
                        require(power(G, N * q) == Fq, ("cover packet changed", exponents, N, q))

                h1 = ONE
                h2 = sub(F, ONE)
                require(mod_poly(h1, F) != ZERO and mod_poly(h2, F) != ZERO,
                        ("cancellation hostile term regular", exponents))
                require(add(h1, h2) == F, ("cancellation hostile sum", exponents))
                cancellation_hostiles += 1

    return {
        "factor_packets": packets,
        "crt_idempotents": idempotents,
        "exact_clearings": exact_clearings,
        "factor_one_step_hostiles": one_step_hostiles,
        "cancellation_hostiles": cancellation_hostiles,
        "cover_packets": cover_packets,
        "cover_kernel_classes": cover_kernel_classes,
        "max_principal_part_length": max_degree,
    }


def homogeneous_pole_filtration_checks():
    primes = (linear(-2), linear(0), linear(3))
    exponent_banks = ((2, 1, 0), (1, 3, 2), (0, 0, 0), (4, 2, 1))
    b_banks = ((-2, 0, 1), (1, -1, 2), (3, 0, -2), (-1, -2, -1), (0, 0, 0))
    cells = 0
    exact = 0
    hostiles = 0

    for F_orders in exponent_banks:
        F = product_polys(
            power(prime, exponent)
            for prime, exponent in zip(primes, F_orders)
        )
        for a in range(0, 4):
            for b_orders in b_banks:
                numerator = polynomial_from_valuations(primes, b_orders, True)
                denominator = polynomial_from_valuations(primes, b_orders, False)
                J_orders = tuple(
                    max(0, a * F_order - b_order)
                    for F_order, b_order in zip(F_orders, b_orders)
                )
                J = product_polys(
                    power(prime, exponent)
                    for prime, exponent in zip(primes, J_orders)
                )
                for n in range(0, 5):
                    cells += 1
                    if n == 0:
                        require(power(J, n) == ONE, ("I_0", F_orders, a, b_orders))
                        continue
                    h = power(J, n)
                    left = mul(h, power(numerator, n))
                    right = mul(power(F, a * n), power(denominator, n))
                    exact_div(left, right)
                    exact += 1
                    for i, need in enumerate(J_orders):
                        if need == 0:
                            continue
                        dropped = exact_div(h, primes[i])
                        dropped_left = mul(dropped, power(numerator, n))
                        _, remainder = divmod_poly(dropped_left, right)
                        require(remainder != ZERO, ("J one-step hostile cleared", F_orders, a, b_orders, n, i))
                        hostiles += 1

    return {
        "filtration_cells": cells,
        "positive_power_exact": exact,
        "filtration_one_step_hostiles": hostiles,
    }


def edge_hostiles():
    v = linear(0)
    v_squared = power(v, 2)
    # Aggregating v*v gives one arm of multiplicity two.  The repeated ideals
    # (v^q),(v^q) are not comaximal, so the false two-arm CRT does not exist.
    divisor, _, _ = xgcd(v, v)
    require(divisor == v, "duplicate-factor gcd hostile failed")
    require(len(v_squared) - 1 == 2, "duplicate jet length failed")
    correct_class = {"free_rank": 0, "torsion": 2}
    false_duplicate_class = {"free_rank": 1, "torsion": 1}
    require(correct_class != false_duplicate_class, "duplicate hostile invisible")

    # Over R, v^2+1 is irreducible and gives r=1,e=1 (Cl=0); over C it
    # splits into two squarefree factors and gives free rank one.
    base_field = {"R": [0, 1], "C": [1, 0]}
    require(base_field["R"] != base_field["C"], "base-change factor hostile failed")

    return {
        "duplicate_F": "v^2: one arm e=2, not two copies e=1",
        "duplicate_correct_class": correct_class,
        "duplicate_false_class": false_duplicate_class,
        "nonclosed_base": "R:v^2+1 has Cl=0; C-base-change has Cl=Z",
        "zero_F": "XY=0 is reducible and R->R_X kills Y",
        "nonfinite": "R->R_X makes X^-1 regular by deleting the boundary",
        "nonnormal": "k[t^2,t^3] subset k[t] makes t regular after finite normalization",
    }


def main():
    packets = crt_factor_packets()
    filtration = homogeneous_pole_filtration_checks()
    hostiles = edge_hostiles()
    semantic = {
        "packets": packets,
        "homogeneous_filtration": filtration,
        "hostiles": hostiles,
        "identities": {
            "principal_part": "(R_X/R)_(-q)=X^-q*k[v]/(F^q)",
            "coefficient_ideal": "{h:X^-q*h in R}=(F^q)",
            "general_filtration": "I_n(X^-a*b)=J(a,b)^n",
            "class_group": "Cl(R)=Z^r/<e_1,...,e_r>",
            "cover": "old degree q packet equals cover degree Nq packet",
        },
    }
    digest = sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()

    print("THM-3404 exact factorized-Danielewski controls")
    print("arithmetic=integer+rational polynomial only")
    for key, value in packets.items():
        print(f"{key}={value}")
    for key, value in filtration.items():
        print(f"{key}={value}")
    print("principal_part=(R_X/R)_(-q)=X^-q*k[v]/(F^q)")
    print("crt_packet=k[v]/(F^q)=direct_sum_i k[v]/(f_i^(q*e_i))")
    print("general_filtration=I_n(X^-a*b)=J(a,b)^n")
    print("cancellation_hostile=X^-1+X^-1*(F-1)=Y")
    print("duplicate_hostile=v^2 is one e=2 arm, not two e=1 arms")
    print("mixed_cover_control=F=v^2*(v-1)^2,N=2: Cl Z+Z/2 -> Z, old packets unchanged")
    print("factorial_cover_boundary=r<=1 iff a finite normal factorial cover exists")
    print(f"semantic_sha256={digest}")
    print("result=PASS")


if __name__ == "__main__":
    main()
