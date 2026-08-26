#!/usr/bin/env python3
"""Exact primary certificate for THM-4193.

The proved all-order statement concerns the nontransitive prefixes

    A_n = C3 ▷ P_n

inside transitive contexts ``(P_b,P_c)``.  The certificate also performs the
complete order-at-most-eight singleton census and a deliberately non-load-
bearing arbitrary-context census through factor order seven.

Rows are tournament-class presentations, not child isomorphism classes.
"""

from __future__ import annotations

import hashlib

import tournament_ordinal_cocycle_parity_thm4184 as base
import tournament_source_padding_supermodularity_thm4187 as source


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def transitive_gate(order: int) -> int:
    """THM-4184's exact G_+(P_order), including the order-one endpoint."""
    need(order >= 1, "positive transitive order")
    if order == 1:
        return 0
    numerator = (
        4 * 4**order
        - 9 * (order + 2) * 2**order
        + 24 * order
        + 32
    )
    need(numerator % 18 == 0, "transitive gate integrality")
    return numerator // 18


def transitive_remainder(left: int, right: int) -> int:
    return (
        transitive_gate(left + right)
        - transitive_gate(left)
        - transitive_gate(right)
    )


def alpha(tail: int) -> int:
    """F_{C3▷P_tail}(1,1).  P_0 denotes the empty tail."""
    need(tail >= 0, "nonnegative tail")
    return 12 * (
        2 * 4**tail - (3 * tail + 21) * 2**tail + 1
    )


def transitive_context_formula(tail: int, middle: int, right: int) -> int:
    """F_{C3▷P_tail}(P_middle,P_right) from the normalized cocycle."""
    need(tail >= 0 and middle >= 1 and right >= 1, "context domain")
    cycle_middle = tail + middle
    return (
        9 * transitive_remainder(cycle_middle, right)
        - transitive_remainder(middle, right)
        + source.negative_theta_formula(cycle_middle, right)
    )


def context_increment(tail: int, middle: int, right: int) -> int:
    """F_(tail+1)-F_tail for tail+middle>=2."""
    total = tail + middle
    need(total >= 2 and right >= 1, "increment formula domain")
    bracket = (
        2**total * (4**right - 1)
        - 3
        * (
            2**right * (total + right + 6)
            - (total + 6)
        )
    )
    return 6 * 2**total * bracket


def f5_q_numerator(middle: int, right: int) -> int:
    """18F_5 before the Q_1 endpoint corrections."""
    b = middle
    c = right
    return (
        36860 * 4 ** (b + c)
        - 36860 * 4**b
        - 32 * 4**c
        - 10359 * b * 2 ** (b + c)
        - 10359 * c * 2 ** (b + c)
        - 93294 * 2 ** (b + c)
        + 10359 * b * 2**b
        + 93294 * 2**b
        + 72 * c * 2**c
        + 144 * 2**c
        - 256
    )


def f5_numerator(middle: int, right: int) -> int:
    """The exact integer 18F_5, with T_1=0 rather than Q_1=1."""
    return (
        f5_q_numerator(middle, right)
        + 144 * (right == 1)
        - 18 * (middle == 1)
    )


def direct_prefix(
    cycle: base.TournamentData,
    transitive: tuple[base.TournamentData, ...],
    tail: int,
) -> base.TournamentData:
    if tail == 0:
        return cycle
    return base.ordinal_data(cycle, transitive[tail - 1])


def direct_context_defect(
    prefix: base.TournamentData,
    middle: base.TournamentData,
    right: base.TournamentData,
) -> int:
    return (
        base.remainder(base.ordinal_data(prefix, middle), right)
        - base.remainder(middle, right)
    )


def main() -> None:
    executable = base.find_gentourng()
    base.EXPECTED_CLASSES[7] = 456
    base.EXPECTED_CLASSES[8] = 6880

    one = base.tournament_data(base.transitive(1))
    cycle = base.tournament_data(base.parse("101", 3))
    need(cycle.hamilton == 3 and cycle.gate == 0, "C3 control")

    transitive = []
    current = one
    for order in range(1, 42):
        if order > 1:
            current = base.ordinal_data(current, one)
        need(current.out == base.transitive(order), "transitive recursion")
        need(current.gate == transitive_gate(order), "transitive gate formula")
        transitive.append(current)

    # Exact singleton formula and its sharp crossing.
    singleton_rows = []
    for tail in range(0, 40):
        prefix = direct_prefix(cycle, tuple(transitive), tail)
        direct = direct_context_defect(prefix, one, one)
        need(direct == alpha(tail), "singleton crossing formula")
        singleton_rows.append((tail, direct))
    need([value for _, value in singleton_rows[:5]] == [-216, -468, -900, -1332, -180],
         "negative singleton ledger")
    need(singleton_rows[5] == (5, 10764), "first positive singleton")
    need(all(value < 0 for _, value in singleton_rows[:5]), "negative side")
    need(all(value > 0 for _, value in singleton_rows[5:]), "positive side")

    # Direct ordinal-capacity verification of the all-order closed formula on
    # a substantial exact head.  The theorem proof, not this finite head,
    # proves the unbounded quantifiers.
    direct_formula_checks = 0
    for tail in range(0, 9):
        prefix = direct_prefix(cycle, tuple(transitive), tail)
        for middle in range(1, 9):
            for right in range(1, 9):
                direct = direct_context_defect(
                    prefix,
                    transitive[middle - 1],
                    transitive[right - 1],
                )
                formula = transitive_context_formula(tail, middle, right)
                need(direct == formula, "transitive-context closed formula")
                direct_formula_checks += 1
    need(direct_formula_checks == 9 * 8 * 8, "direct formula check total")

    # Algebraic identity checks for the two proof reductions.  These are
    # coefficient/formula checks; positivity is established symbolically in
    # the theorem by the displayed lower bounds.
    formula_checks = 0
    increment_checks = 0
    f5_checks = 0
    for tail in range(0, 41):
        for middle in range(1, 65):
            for right in range(1, 65):
                value = transitive_context_formula(tail, middle, right)
                if tail >= 5:
                    need(value > 0, "positive transitive context")
                if tail < 40 and tail + middle >= 2:
                    difference = (
                        transitive_context_formula(tail + 1, middle, right)
                        - value
                    )
                    need(
                        difference == context_increment(tail, middle, right),
                        "tail-increment identity",
                    )
                    if tail >= 5:
                        need(difference > 0, "positive tail increment")
                    increment_checks += 1
                if tail == 5:
                    need(18 * value == f5_numerator(middle, right),
                         "F5 expanded numerator")
                    need(f5_numerator(middle, right) > 0, "F5 positivity")
                    f5_checks += 1
                formula_checks += 1

    # Complete class census for the first source-free crossing.  The stream
    # is ordered by order and gentourng label.
    class_banks = {1: (one,)}
    for order in range(2, 9):
        class_banks[order] = tuple(
            base.tournament_data(out) for out in base.classes(executable, order)
        )
    census_digest = hashlib.sha256()
    singleton_census = []
    survivor_rows = []
    for order in range(1, 9):
        source_positive = 0
        source_nonpositive = 0
        source_free_negative = 0
        source_free_nonnegative = 0
        for data in class_banks[order]:
            has_source = any(row.bit_count() == order - 1 for row in data.out)
            value = direct_context_defect(data, one, one)
            if has_source:
                if value > 0:
                    source_positive += 1
                else:
                    source_nonpositive += 1
            else:
                if value < 0:
                    source_free_negative += 1
                else:
                    source_free_nonnegative += 1
                    survivor_rows.append((order, base.label(data.out), value))
            row = (order, base.label(data.out), int(has_source), value)
            census_digest.update(("|".join(map(str, row)) + "\n").encode("ascii"))
        singleton_census.append(
            (
                order,
                len(class_banks[order]),
                source_positive,
                source_nonpositive,
                source_free_negative,
                source_free_nonnegative,
            )
        )
    expected_census = [
        (1, 1, 0, 1, 0, 0),
        (2, 1, 1, 0, 0, 0),
        (3, 2, 1, 0, 1, 0),
        (4, 4, 2, 0, 2, 0),
        (5, 12, 4, 0, 8, 0),
        (6, 56, 12, 0, 44, 0),
        (7, 456, 56, 0, 400, 0),
        (8, 6880, 456, 0, 6423, 1),
    ]
    need(singleton_census == expected_census, "singleton census ledger")
    survivor = direct_prefix(cycle, tuple(transitive), 5)
    expected_survivor = (8, base.label(survivor.out), 10764)
    need(survivor_rows == [expected_survivor], "unique source-free survivor")

    # Non-load-bearing arbitrary-context stress test of the first survivor.
    factors = tuple(
        data
        for order in range(1, 8)
        for data in class_banks[order]
    )
    need(len(factors) == 532, "factor bank through seven")
    remainders = {
        (id(middle), id(right)): base.remainder(middle, right)
        for middle in factors
        for right in factors
    }
    arbitrary_checks = 0
    arbitrary_negative = 0
    arbitrary_minimum = None
    domination_minimum = None
    for middle in factors:
        prefix_middle = base.ordinal_data(survivor, middle)
        for right in factors:
            value = (
                base.remainder(prefix_middle, right)
                - remainders[(id(middle), id(right))]
            )
            record = (
                value,
                len(middle.out),
                base.label(middle.out),
                len(right.out),
                base.label(right.out),
            )
            if arbitrary_minimum is None or record < arbitrary_minimum:
                arbitrary_minimum = record
            arbitrary_negative += value < 0
            domination = value - alpha(5) * middle.hamilton**2 * right.hamilton**2
            domination_record = (
                domination,
                len(middle.out),
                base.label(middle.out),
                len(right.out),
                base.label(right.out),
            )
            if domination_minimum is None or domination_record < domination_minimum:
                domination_minimum = domination_record
            arbitrary_checks += 1
    need(arbitrary_checks == 532**2, "arbitrary-context presentation total")
    need(arbitrary_negative == 0, "arbitrary-context negative row")
    need(arbitrary_minimum == (10764, 1, "", 1, ""),
         "arbitrary-context minimum")
    need(domination_minimum == (0, 1, "", 1, ""),
         "finite singleton-domination minimum")

    print("TOURNAMENT_CYCLE_FIRST_TRANSITIVE_TAIL_CROSSING")
    print("singleton_formula", "12*(2*4^n-(3n+21)*2^n+1)")
    print("singleton_rows_n0_to_n5", singleton_rows[:6])
    print("singleton_sign_crossing", "negative_n_le_4", "positive_n_ge_5")
    print("direct_transitive_context_formula_checks", direct_formula_checks)
    print("algebra_formula_checks", formula_checks)
    print("algebra_increment_checks", increment_checks)
    print("algebra_f5_checks", f5_checks)
    print("first_positive_prefix", "C3▷P5", "value", alpha(5))
    print("singleton_census", singleton_census)
    print("unique_source_free_order_le_8_survivor", expected_survivor)
    print("singleton_census_sha256", census_digest.hexdigest())
    print("arbitrary_context_factor_classes", len(factors))
    print("arbitrary_context_presentations", arbitrary_checks)
    print("arbitrary_context_negative", arbitrary_negative)
    print("arbitrary_context_minimum", arbitrary_minimum)
    print("finite_singleton_domination_minimum", domination_minimum)
    print("PASS")


if __name__ == "__main__":
    main()
