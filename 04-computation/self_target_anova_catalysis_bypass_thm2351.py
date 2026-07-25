#!/usr/bin/env python3
"""Optimization-safe exact controls for THM-2351.

The script checks the self-target two-token table, its exact uniform ANOVA,
the directional catalysis/bypass inverse, the continuation-row identities,
the boundary cases, equal-prime symmetry, and the Brittenham--Hermiller
parameter family.  The enumerated ledger bank is algebraic and is not
asserted to be realized by knots.
"""

from __future__ import annotations

from fractions import Fraction


Q = Fraction


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def matrix_anova(matrix: tuple[tuple[Q, Q], tuple[Q, Q]]):
    e00, e01 = matrix[0]
    e10, e11 = matrix[1]
    return (
        (e00 + e01 + e10 + e11) / 4,
        (e00 + e01 - e10 - e11) / 4,
        (e00 - e01 + e10 - e11) / 4,
        (e00 - e01 - e10 + e11) / 4,
    )


def reconstruct(mean: Q, row: Q, column: Q, pair: Q):
    return (
        (
            mean + row + column + pair,
            mean + row - column - pair,
        ),
        (
            mean - row + column - pair,
            mean - row - column + pair,
        ),
    )


def self_target_table(
    u_p: int,
    u_q: int,
    u_x: int,
    c_p: int,
    c_q: int,
) -> tuple[tuple[Q, Q], tuple[Q, Q]]:
    return (
        (
            Q(u_x),
            Q(2 * u_p - c_p - u_x),
        ),
        (
            Q(2 * u_q - c_q - u_x),
            Q(-u_x),
        ),
    )


def check_one_ledger(
    u_p: int,
    u_q: int,
    u_x: int,
    c_p: int,
    c_q: int,
) -> tuple[Q, Q, Q, Q]:
    sigma = u_p + u_q - u_x
    b_p = sigma - c_p
    b_q = sigma - c_q
    require(
        min(sigma, c_p, c_q, b_p, b_q) >= 0,
        "ledger is not algebraically admissible",
    )

    table = self_target_table(u_p, u_q, u_x, c_p, c_q)
    mean, row, column, pair = matrix_anova(table)
    require(
        reconstruct(mean, row, column, pair) == table,
        "ANOVA reconstruction failed",
    )
    require(mean == -pair, "mean/pair sign relation failed")
    require(
        mean == Q(b_p + b_q, 4),
        "mean is not one quarter total bypass",
    )
    require(
        pair == -Q(b_p + b_q, 4),
        "pair is not minus one quarter total bypass",
    )
    require(row + column == u_x, "unary sum did not recover u(X)")
    require(
        row - column == Q(u_p - u_q) + Q(c_q - c_p, 2),
        "corrected unary imbalance failed",
    )

    recovered_c_p = (
        sigma
        + 2 * pair
        - row
        + column
        + u_p
        - u_q
    )
    recovered_c_q = (
        sigma
        + 2 * pair
        + row
        - column
        - u_p
        + u_q
    )
    recovered_b_p = (
        -2 * pair + row - column - u_p + u_q
    )
    recovered_b_q = (
        -2 * pair - row + column + u_p - u_q
    )
    require(
        (recovered_c_p, recovered_c_q, recovered_b_p, recovered_b_q)
        == (c_p, c_q, b_p, b_q),
        "directional inverse failed",
    )

    require(-Q(sigma, 2) <= pair <= 0, "pair interval failed")
    require(
        (pair == -Q(sigma, 2)) == (c_p == c_q == 0),
        "pure-bypass boundary iff failed",
    )
    require(
        (pair == 0) == (b_p == b_q == 0),
        "zero-bypass boundary iff failed",
    )
    require(
        (pair > -Q(sigma, 2)) == (c_p + c_q > 0),
        "positive-catalysis test failed",
    )

    corrected_imbalance = row - column - Q(u_p - u_q)
    require(
        corrected_imbalance == Q(c_q - c_p, 2),
        "tournament unary contrast failed",
    )
    require(
        abs(corrected_imbalance)
        <= min(Q(sigma) + 2 * pair, -2 * pair),
        "ANOVA diamond failed",
    )

    mu_u = Q(-sigma)
    mu_x = Q(-sigma + c_p + c_q)
    require(mu_x - mu_u == c_p + c_q, "Mobius catalytic drift failed")
    require(-mu_u - mu_x == b_p + b_q, "Mobius bypass sum failed")
    require(pair == (mu_u + mu_x) / 4, "Mobius/ANOVA bridge failed")
    require(mu_u <= mu_x <= -mu_u, "Mobius interval failed")

    raw = tuple(
        tuple(entry + u_x for entry in row_entries)
        for row_entries in table
    )
    require(
        raw
        == (
            (Q(2 * u_x), Q(2 * u_p - c_p)),
            (Q(2 * u_q - c_q), Q(0)),
        ),
        "raw lift-cost table failed",
    )
    require(
        raw[1][1] < min(raw[0][0], raw[0][1], raw[1][0]),
        "self-target optimum was not unique",
    )
    return mean, row, column, pair


def ledger_bank() -> tuple[int, int, int, int]:
    cases = 0
    pure_bypass = 0
    zero_bypass = 0
    asymmetric = 0
    for u_p in range(1, 7):
        for u_q in range(1, 7):
            lower = max(1, abs(u_p - u_q))
            for u_x in range(lower, u_p + u_q + 1):
                sigma = u_p + u_q - u_x
                for c_p in range(min(u_p, sigma) + 1):
                    for c_q in range(min(u_q, sigma) + 1):
                        check_one_ledger(u_p, u_q, u_x, c_p, c_q)
                        cases += 1
                        pure_bypass += int(c_p == c_q == 0)
                        zero_bypass += int(c_p == c_q == sigma)
                        asymmetric += int(c_p != c_q)
    return cases, pure_bypass, zero_bypass, asymmetric


def boundary_hostiles() -> None:
    pure_bypass_table = self_target_table(2, 2, 3, 0, 0)
    pure_translation_table = self_target_table(2, 2, 3, 1, 1)
    require(
        pure_bypass_table == ((Q(3), Q(1)), (Q(1), Q(-3))),
        "pure-bypass boundary control changed",
    )
    require(
        pure_translation_table == ((Q(3), Q(0)), (Q(0), Q(-3))),
        "pure-translation boundary control changed",
    )
    require(
        matrix_anova(pure_bypass_table)[3] == Q(-1, 2),
        "pure-bypass pair coefficient changed",
    )
    require(
        matrix_anova(pure_translation_table)[3] == 0,
        "pure-translation pair coefficient changed",
    )
    require(
        min(
            (i, j)
            for i in range(2)
            for j in range(2)
            if pure_bypass_table[i][j]
            == min(value for row in pure_bypass_table for value in row)
        )
        == (1, 1),
        "pure-bypass owner control failed",
    )
    require(
        min(
            (i, j)
            for i in range(2)
            for j in range(2)
            if pure_translation_table[i][j]
            == min(value for row in pure_translation_table for value in row)
        )
        == (1, 1),
        "pure-translation owner control failed",
    )


def equal_prime_bank() -> int:
    cases = 0
    for u_p in range(1, 7):
        for u_x in range(1, 2 * u_p + 1):
            sigma = 2 * u_p - u_x
            for common_c in range(min(u_p, sigma) + 1):
                table = self_target_table(
                    u_p, u_p, u_x, common_c, common_c
                )
                mean, row, column, pair = check_one_ledger(
                    u_p, u_p, u_x, common_c, common_c
                )
                require(row == column, "equal-prime unary fields differ")
                require(
                    table[0][1] == table[1][0],
                    "equal-prime split allocations differ",
                )
                common_bypass = sigma - common_c
                require(
                    pair == -Q(common_bypass, 2),
                    "equal-prime quotient pair coefficient failed",
                )
                require(mean == -pair, "equal-prime mean failed")
                cases += 1
    return cases


def brittenham_hermiller_bank() -> int:
    cases = 0
    for u_x in range(1, 6):
        sigma = 6 - u_x
        mean, row, column, pair = check_one_ledger(
            3, 3, u_x, 0, 0
        )
        table = self_target_table(3, 3, u_x, 0, 0)
        require(
            table == (
                (Q(u_x), Q(sigma)),
                (Q(sigma), Q(-u_x)),
            ),
            "Brittenham-Hermiller table failed",
        )
        require(
            (mean, row, column, pair)
            == (
                Q(sigma, 2),
                Q(u_x, 2),
                Q(u_x, 2),
                Q(-sigma, 2),
            ),
            "Brittenham-Hermiller ANOVA failed",
        )
        cases += 1
    return cases


ledger_cases, pure_bypass_cases, zero_bypass_cases, asymmetric_cases = (
    ledger_bank()
)
boundary_hostiles()
equal_prime_cases = equal_prime_bank()
bh_cases = brittenham_hermiller_bank()

print("theorem=THM-2351")
print("status=CLAIMED+VERIFIED-EXACT+UNDER-INDEPENDENT-AUDIT")
print(f"algebraic_ledger_cases={ledger_cases}")
print(f"pure_bypass_boundary_cases={pure_bypass_cases}")
print(f"zero_bypass_boundary_cases={zero_bypass_cases}")
print(f"directionally_asymmetric_cases={asymmetric_cases}")
print(f"equal_prime_quotient_cases={equal_prime_cases}")
print(f"brittenham_hermiller_parameter_cases={bh_cases}")
print("self_target_optimum_is_always_unique=YES")
print("pair_tensor_equals_minus_total_bypass_over_four=YES")
print("root_plus_pair_detects_any_mutual_catalysis=YES")
print("corrected_unaries_recover_both_directions=YES")
print("continuation_rows_recover_catalysis_and_bypass=YES")
print("pair_skew_sector=ZERO")
print("full_anova_stronger_than_directional_pair_ledger=NO")
print("all_checks=PASS")
