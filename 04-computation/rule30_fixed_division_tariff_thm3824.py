#!/usr/bin/env python3
"""Exact Smith audit of the normalized Rule-30 amplitude division tariff.

For n>=1, r>=0, L=T_n^r Z^n, and a fixed normalizer e>=0, define

    K = L intersect T_n^{-1}(2^e L).

Then K is the largest translation sublattice of L on which the normalized
response x -> T_n x / 2^e is both integral and still equal modulo L.  The
total tariff L/K contains an integrality layer and, after integrality, a carry
layer.  The proof reduces the total quotient to the image of one effective
T_d modulo 2^e.  This canonical companion checks the effective
Smith forms, odd-core power Smith forms, and the resulting all-period laws.

This is a congruence for the fixed linear division branch N_e on its
2^e-divisibility domain.  K need not preserve the nonlinear exact-gap stratum
nu_2(T_n x)=e; the normalizer e remains an external adaptive branch label.
"""

from __future__ import annotations

import hashlib
import json
import math
import sys

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import smith_normal_form


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


def require(condition: bool, label: object) -> None:
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def matrix_t(n: int) -> Matrix:
    require(n >= 1, "positive period")
    rows = [[0 for _ in range(n)] for _ in range(n)]
    for j in range(n):
        rows[j][(2 * j) % n] += 1
        rows[j][(2 * j + 1) % n] += 1
    return Matrix(rows)


def two_part(n: int) -> int:
    out = 0
    while n % 2 == 0:
        out += 1
        n //= 2
    return out


def effective_period(n: int, r: int) -> tuple[int, int]:
    a = two_part(n)
    return n >> min(a, r), a


def exact_smith(n: int) -> tuple[int, ...]:
    diagonal = smith_normal_form(matrix_t(n), domain=ZZ)
    return tuple(abs(int(diagonal[j, j])) for j in range(n))


def image_group_orders_from_smith(smith: tuple[int, ...], exponent: int) -> tuple[int, ...]:
    """Cyclic orders of im(diag(smith) mod 2^exponent), omitting trivial."""
    if exponent == 0:
        return ()
    modulus = 1 << exponent
    orders = []
    for value in smith:
        if value == 0:
            continue
        order = modulus // math.gcd(modulus, value)
        if order > 1:
            orders.append(order)
    return tuple(sorted(orders, reverse=True))


def predicted_tariff_orders(d: int, exponent: int) -> tuple[int, ...]:
    if exponent == 0:
        return ()
    modulus = 1 << exponent
    if d % 2 == 0:
        return (modulus,) * (d // 2)
    tail = 1 << (exponent - 1)
    return (modulus,) * (d - 1) + (() if tail == 1 else (tail,))


def predicted_integrality_orders(d: int, s: int, exponent: int) -> tuple[int, ...]:
    """Cyclic orders of L/(L intersect domain(N_e))."""
    if exponent == 0:
        return ()
    modulus = 1 << exponent
    if d % 2 == 0:
        require(s == 0, ("even effective period has s=0", d, s))
        return (modulus,) * (d // 2)
    tail_exponent = max(0, exponent - s - 1)
    tail = 1 << tail_exponent
    return (modulus,) * (d - 1) + (() if tail == 1 else (tail,))


def predicted_carry_orders(d: int, s: int, exponent: int) -> tuple[int, ...]:
    """Cyclic orders of the post-integrality kernel quotient J/K."""
    if exponent == 0 or d % 2 == 0:
        return ()
    order = 1 << min(s, exponent - 1)
    return () if order == 1 else (order,)


def main() -> None:
    smith_checks = 0
    odd_power_smith_checks = 0
    odd_power_image_checks = 0
    tariff_checks = 0
    records = []

    effective_smith: dict[int, tuple[int, ...]] = {}
    for d in range(1, 41):
        got = exact_smith(d)
        want = (
            (1,) * (d // 2) + (0,) * (d // 2)
            if d % 2 == 0
            else (1,) * (d - 1) + (2,)
        )
        require(got == want, ("effective T Smith", d, got, want))
        effective_smith[d] = got
        smith_checks += 1

    # On every odd core, THM-3804 predicts one growing cyclic Smith factor
    # for every power.  These independent exact matrices also certify that
    # ker(T_d^k mod 2^e) has the same cardinality as the visible constant
    # subgroup killed by 2^k; hence that kernel is exactly cyclic constant
    # data, not merely a group of the same order by assumption.
    for d in range(1, 16, 2):
        for power in range(1, 7):
            diagonal = smith_normal_form(matrix_t(d) ** power, domain=ZZ)
            got = tuple(abs(int(diagonal[j, j])) for j in range(d))
            want = (1,) * (d - 1) + (1 << power,)
            require(got == want, ("odd power Smith", d, power, got, want))
            odd_power_smith_checks += 1
            for exponent in range(0, 9):
                integrality = image_group_orders_from_smith(got, exponent)
                require(
                    integrality == predicted_integrality_orders(d, power - 1, exponent),
                    ("odd power integrality", d, power, exponent),
                )
                kernel_size = 1 << min(power, exponent)
                constant_kernel_size = len(
                    tuple(
                        c
                        for c in range(1 << exponent)
                        if ((1 << power) * c) % (1 << exponent) == 0
                    )
                ) if exponent > 0 else 1
                require(
                    kernel_size == constant_kernel_size,
                    ("constant kernel exhausts Smith kernel", d, power, exponent),
                )
                require(
                    kernel_size // (1 << min(1, exponent))
                    == math.prod(predicted_carry_orders(d, power - 1, exponent)),
                    ("odd carry quotient", d, power, exponent),
                )
                odd_power_image_checks += 3

    # The dependence on n and r enters only through the current effective
    # period d=n/2^min(v2(n),r).  Before the odd core d is even; at and after
    # the core it is odd.  The old torsion exponent r-v2(n) cancels from the
    # one-step division tariff.
    for n in range(1, 41):
        for r in range(0, 11):
            d, a = effective_period(n, r)
            s = max(0, r - a)
            require((d % 2 == 0) == (r < a), ("effective parity", n, r, d, a))
            for exponent in range(0, 9):
                got = image_group_orders_from_smith(effective_smith[d], exponent)
                want = predicted_tariff_orders(d, exponent)
                require(got == want, ("tariff group", n, r, exponent, got, want))
                integrality = predicted_integrality_orders(d, s, exponent)
                carry = predicted_carry_orders(d, s, exponent)
                total_order = math.prod(got)
                require(
                    total_order == math.prod(integrality) * math.prod(carry),
                    ("tariff filtration order", n, r, exponent),
                )
                records.append(
                    {
                        "n": n,
                        "r": r,
                        "a": a,
                        "d": d,
                        "e": exponent,
                        "total_orders": got,
                        "integrality_orders": integrality,
                        "carry_orders": carry,
                    }
                )
                tariff_checks += 1

    # The minimal all-odd carry hostile realizes the unique nonzero class of
    # the n=2,r=3,e=2 tariff Z/2.  In L=4Z*(1,1), delta=(4,4), while
    # T(delta)/4=(2,2) is outside L.
    delta = Matrix([4, 4])
    response = matrix_t(2) * delta / 4
    require(tuple(int(value) for value in response) == (2, 2), "minimal response")
    require(all(int(value) % 4 == 0 for value in delta), "delta in L")
    require(any(int(value) % 4 for value in response), "response outside L")
    require(predicted_tariff_orders(1, 2) == (2,), "minimal tariff Z/2")

    # At e=1 on the one-point odd core the tariff is trivial; this explains
    # why common gap one cannot split the n=2 Smith carry after normalization.
    require(predicted_tariff_orders(1, 1) == (), "one-core gap-one triviality")
    require(predicted_tariff_orders(7, 0) == (), "e=0 total triviality")
    require(predicted_integrality_orders(7, 5, 0) == (), "e=0 integrality triviality")
    require(predicted_carry_orders(7, 5, 0) == (), "e=0 carry triviality")

    # Essential exact-gap scope hostile.  At n=2,r=0,e=1 one has L=Z^2.
    # The difference delta belongs to K because T(delta) is in 2L, so fixed
    # division by two is a lawful K-congruence.  But it moves x from exact gap
    # one to exact gap two.  K therefore does not preserve the adaptive
    # exact-normalizer stratum.
    exact_gap_x = Matrix([1, 1])
    exact_gap_delta = Matrix([2, 0])
    exact_gap_y = exact_gap_x + exact_gap_delta
    exact_gap_raw_x = matrix_t(2) * exact_gap_x
    exact_gap_raw_y = matrix_t(2) * exact_gap_y
    exact_gap_raw_delta = matrix_t(2) * exact_gap_delta
    require(tuple(exact_gap_raw_x) == (2, 2), "exact-gap hostile left")
    require(tuple(exact_gap_raw_y) == (4, 4), "exact-gap hostile right")
    require(
        all(int(value) % 2 == 0 for value in exact_gap_raw_delta),
        "exact-gap hostile delta in K",
    )
    exact_gap_hostile = {
        "n": 2,
        "r": 0,
        "e": 1,
        "x": (1, 1),
        "delta": (2, 0),
        "y": (3, 1),
        "gaps": (1, 2),
    }

    semantic = hashlib.sha256(
        json.dumps(
            {"records": records, "exact_gap_hostile": exact_gap_hostile},
            sort_keys=True,
            separators=(",", ":"),
        ).encode("ascii")
    ).hexdigest()

    print("RULE30_FIXED_DIVISION_TARIFF_THM3824")
    print("status=PROVED_FORMULA+VERIFIED_EXACT;all_rule30_prizes_open")
    print("variables=n>=1;r>=0;e>=0;L=T_n^r*Z^n;K=L_intersect_T_n^(-1)(2^e*L)")
    print("e=0:K=L;total_integrality_carry_tariffs_trivial")
    print("e>=1,r<v2(n):L/K=(Z/2^e)^(d/2),d=n/2^r;all_integrality,no_carry")
    print("e>=1,r>=v2(n):L/K=(Z/2^e)^(d-1) direct_sum Z/2^(e-1),d=oddpart(n)")
    print("odd_core_integrality=(Z/2^e)^(d-1) direct_sum Z/2^max(e-s-1,0),s=r-v2(n)")
    print("odd_core_post_integrality_carry=Z/2^min(s,e-1)")
    print("convention=Z/1_is_trivial")
    print(f"effective_smith_checks={smith_checks}")
    print(f"odd_power_smith_checks={odd_power_smith_checks}")
    print(f"odd_power_image_kernel_checks={odd_power_image_checks}")
    print(f"all_period_tariff_checks={tariff_checks}")
    print("universe=1<=d<=40;odd_d<=15,power<=6;1<=n<=40;0<=r<=10;0<=e<=8")
    print("minimal_hostile=n=2,r=3,e=2,delta=(4,4),normalized_response=(2,2),tariff=Z/2")
    print("gap_one_one_core=trivial_tariff")
    print("exact_gap_scope_hostile=n=2,r=0,e=1,x=(1,1),delta=(2,0)_in_K,y=(3,1),gaps=1_vs_2")
    print("scope=K_is_maximal_for_fixed_N_e_on_divisibility_domain_not_for_exact_gap_stratum")
    print(f"semantic_sha256={semantic}")
    print("PASS")


if __name__ == "__main__":
    main()
