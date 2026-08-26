#!/usr/bin/env python3
"""Exact primary audit for THM-4187.

The all-order proof is symbolic.  This certificate checks the literal padding
defect decomposition on every ordered pair of tournament-class representatives
through factor order six, its transitive-left telescoping consequences, the
quantitative terminal-C3 specialization, and the cycle-first negative family.
Rows are ordered factor presentations, not child isomorphism classes.
"""

from __future__ import annotations

import hashlib

import tournament_ordinal_cocycle_parity_thm4184 as base


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def incoming_mass(data: base.TournamentData, vertex: int) -> int:
    return sum(
        data.capacities[vertex][other]
        for other in range(len(data.out))
        if data.out[other] & (1 << vertex)
    )


def delta(one: base.TournamentData, data: base.TournamentData) -> int:
    """Universal-source padding increment Delta(X)=G(1▷X)-G(X)."""
    return base.ordinal_data(one, data).gate - data.gate


def padding_formula(
    one: base.TournamentData, data: base.TournamentData
) -> tuple[int, tuple[int, ...], tuple[int, ...]]:
    """Return THM-4177's literal nonnegative padding decomposition."""
    child = base.ordinal_data(one, data)
    root_capacities = tuple(child.capacities[0][1 + i] for i in range(len(data.out)))
    degrees = tuple(sum(row) for row in data.capacities)
    defects = []
    value = 0
    for vertex, root_capacity in enumerate(root_capacities):
        incoming_root = sum(
            root_capacities[other]
            for other in range(len(data.out))
            if data.out[other] & (1 << vertex)
        )
        fan_defect = incoming_root - incoming_mass(data, vertex)
        need(fan_defect >= 0, "incoming-fan defect is negative")
        defects.append(fan_defect)
        value += root_capacity * (data.mass - degrees[vertex])
        value += 4 * root_capacity * fan_defect
    need(value == delta(one, data), "padding formula mismatch")
    return value, root_capacities, tuple(defects)


def explicit_supermodularity_decomposition(
    one: base.TournamentData,
    left: base.TournamentData,
    right: base.TournamentData,
) -> tuple[int, int, int]:
    """Compute the exact B-extra and C-vertex terms in the proof."""
    product = base.ordinal_data(left, right)
    product_child = base.ordinal_data(one, product)
    left_child = base.ordinal_data(one, left)
    nleft = len(left.out)
    nright = len(right.out)

    old_root = tuple(left_child.capacities[0][1 + i] for i in range(nleft))
    all_root = tuple(
        product_child.capacities[0][1 + i] for i in range(nleft + nright)
    )
    for i in range(nleft):
        need(
            all_root[i] == right.hamilton * old_root[i],
            "left root-capacity scaling mismatch",
        )

    cross_rows = tuple(
        sum(product.capacities[i][nleft + j] for j in range(nright))
        for i in range(nleft)
    )
    cross_mass = sum(cross_rows)
    left_extra = right.hamilton * sum(
        old_root[i]
        * (left.hamilton * right.mass + cross_mass - cross_rows[i])
        for i in range(nleft)
    )

    degrees = tuple(sum(row) for row in product.capacities)
    right_terms = 0
    for j in range(nright):
        vertex = nleft + j
        root_capacity = all_root[vertex]
        incoming_root = sum(
            all_root[other]
            for other in range(nleft + nright)
            if product.out[other] & (1 << vertex)
        )
        fan_defect = incoming_root - incoming_mass(product, vertex)
        need(fan_defect >= 0, "right-block incoming-fan defect is negative")
        right_terms += root_capacity * (product.mass - degrees[vertex])
        right_terms += 4 * root_capacity * fan_defect

    difference = delta(one, product) - right.hamilton**2 * delta(one, left)
    need(difference == left_extra + right_terms, "supermodularity split mismatch")
    return difference, left_extra, right_terms


def theta(
    left: base.TournamentData,
    middle: base.TournamentData,
    right: base.TournamentData,
) -> int:
    return (
        base.remainder(base.ordinal_data(left, middle), right)
        - left.hamilton**2 * base.remainder(middle, right)
    )


def converse_out(out: tuple[int, ...]) -> tuple[int, ...]:
    order = len(out)
    return tuple(
        sum(
            1 << other
            for other in range(order)
            if other != vertex and not (out[vertex] & (1 << other))
        )
        for vertex in range(order)
    )


def minus_gate(data: base.TournamentData) -> int:
    return base.gate(converse_out(data.out), data.capacities)


def minus_remainder(left: base.TournamentData, right: base.TournamentData) -> int:
    child = base.ordinal_data(left, right)
    return (
        minus_gate(child)
        - right.hamilton**2 * minus_gate(left)
        - left.hamilton**2 * minus_gate(right)
    )


def negative_remainder_formula(order: int) -> int:
    need(order >= 1, "positive transitive-tail order required")
    if order == 1:
        return -72
    return 72 - (27 * order + 126) * 2 ** (order - 1)


def negative_theta_formula(left_tail: int, right_tail: int) -> int:
    need(left_tail >= 1 and right_tail >= 1, "positive tail orders required")
    if left_tail == 1:
        return 144 - (27 * right_tail + 153) * 2**right_tail
    return -2 ** (left_tail - 1) * (
        (27 * (left_tail + right_tail) + 126) * 2**right_tail
        - (27 * left_tail + 126)
    )


def main() -> None:
    executable = base.find_gentourng()
    raw_banks = {1: (base.transitive(1),)}
    raw_banks.update({order: base.classes(executable, order) for order in range(2, 7)})
    banks = {
        order: tuple(base.tournament_data(out) for out in raw_banks[order])
        for order in range(1, 7)
    }
    all_data = tuple(data for order in range(1, 7) for data in banks[order])
    need(len(all_data) == 76, "factor-class total")

    one = base.tournament_data(base.transitive(1))
    cycle = base.tournament_data(base.parse("101", 3))
    need(cycle.hamilton == 3 and cycle.gate == 0, "C3 control")

    # Check THM-4177's source-padding decomposition on all inherited factors.
    padding_checks = 0
    for data in all_data:
        padding_formula(one, data)
        padding_checks += 1

    # Complete class-presentation pair audit of source supermodularity.
    pair_checks = 0
    pair_zero = 0
    pair_minimum_positive = None
    pair_digest = hashlib.sha256()
    for left in all_data:
        for right in all_data:
            difference, left_extra, right_terms = explicit_supermodularity_decomposition(
                one, left, right
            )
            direct_theta = theta(one, left, right)
            need(difference == direct_theta, "Theta/delta identity mismatch")
            need(left_extra >= 0 and right_terms >= 0, "negative proof summand")
            need(difference >= 0, "source supermodularity failure")
            expected_zero = len(left.out) == len(right.out) == 1
            need((difference == 0) == expected_zero, "source equality classification")
            pair_zero += difference == 0
            if difference > 0:
                record = (
                    difference,
                    len(left.out),
                    len(right.out),
                    base.label(left.out),
                    base.label(right.out),
                )
                if pair_minimum_positive is None or record < pair_minimum_positive:
                    pair_minimum_positive = record
            row = (
                len(left.out),
                len(right.out),
                base.label(left.out),
                base.label(right.out),
                difference,
                left_extra,
                right_terms,
            )
            pair_digest.update(("|".join(map(str, row)) + "\n").encode("ascii"))
            pair_checks += 1
    need(pair_checks == 76**2 and pair_zero == 1, "pair census totals")

    # Transitive-left telescoping, including the exact zero boundary.
    transitive_factors = []
    current = one
    for order in range(1, 9):
        if order > 1:
            current = base.ordinal_data(one, current)
        transitive_factors.append(current)
        need(current.out == base.transitive(order), "transitive recursion")

    transitive_remainder_checks = 0
    transitive_remainder_zeros = []
    transitive_theta_checks = 0
    transitive_theta_zeros = []
    for left in transitive_factors:
        left_order = len(left.out)
        for right in all_data:
            value = base.remainder(left, right)
            need(value >= 0, "negative transitive-left remainder")
            if value == 0:
                transitive_remainder_zeros.append(
                    (left_order, len(right.out), base.label(right.out))
                )
            transitive_remainder_checks += 1
        for middle in all_data:
            for right in all_data:
                value = theta(left, middle, right)
                need(value >= 0, "negative transitive-left interaction")
                if value == 0:
                    transitive_theta_zeros.append(
                        (
                            left_order,
                            len(middle.out),
                            len(right.out),
                            base.label(middle.out),
                            base.label(right.out),
                        )
                    )
                transitive_theta_checks += 1
    need(
        transitive_remainder_zeros == [(1, 1, ""), (1, 2, "1"), (2, 1, "")],
        "transitive-left remainder zero boundary",
    )
    need(
        transitive_theta_zeros == [(1, 1, 1, "", "")],
        "transitive-left Theta zero boundary",
    )

    # Converse dual: a transitive right factor is nonnegative for R_-.
    dual_remainder_checks = 0
    dual_remainder_zeros = []
    for right in transitive_factors:
        right_order = len(right.out)
        for left in all_data:
            left_converse = base.tournament_data(converse_out(left.out))
            value = minus_remainder(left, right)
            need(
                value == base.remainder(right, left_converse),
                "converse remainder identity",
            )
            need(value >= 0, "negative transitive-right dual remainder")
            if value == 0:
                dual_remainder_zeros.append(
                    (len(left.out), base.label(left.out), right_order)
                )
            dual_remainder_checks += 1
    need(
        dual_remainder_zeros == [(1, "", 1), (2, "1", 1), (1, "", 2)],
        "transitive-right dual zero boundary",
    )

    # Every no-sink factor in the declared bank is a strict OS+ control.
    no_sink = tuple(data for data in all_data if not base.has_sink(data.out))
    os_checks = 0
    for left in transitive_factors:
        for right in no_sink:
            need(base.remainder(left, right) > 0, "transitive-left OS+ failure")
            os_checks += 1

    # Terminal-C3 exact identity and its sharp 648*a*H(X)^2 consequence.
    c3_identity_checks = 0
    c3_telescoping_checks = 0
    c3_os_checks = 0
    c3_os_equality = []
    for data in all_data:
        source_child = base.ordinal_data(one, data)
        source_capacities = tuple(
            source_child.capacities[0][1 + i] for i in range(len(data.out))
        )
        t_values = tuple(2 * row[0] + 4 * row[1] for row in data.starts)
        aggregate = data.mass + 2 * data.hamilton
        degrees = tuple(sum(row) for row in data.capacities)
        outgoing = tuple(
            sum(
                data.capacities[i][j]
                for j in range(len(data.out))
                if data.out[i] & (1 << j)
            )
            for i in range(len(data.out))
        )
        c3_remainder_formula = (
            15 * sum(t * t for t in t_values)
            + 9
            * sum(
                t_values[i] * (data.mass - degrees[i] + 4 * outgoing[i])
                for i in range(len(data.out))
            )
            - 27 * data.mass**2
            - 108 * data.hamilton * data.mass
            - 120 * data.hamilton**2
        )
        need(
            c3_remainder_formula == base.remainder(data, cycle),
            "terminal-C3 block remainder formula",
        )
        need(sum(source_capacities) == aggregate, "C3 source aggregate")
        need(
            sum(t_values) == 3 * aggregate - 2 * data.hamilton,
            "C3 path-cover aggregate",
        )
        exact = 9 * (
            8 * aggregate * (3 * aggregate - data.hamilton)
            - sum(a * t for a, t in zip(source_capacities, t_values))
        )
        direct = theta(one, data, cycle)
        need(exact == direct, "terminal-C3 exact interaction identity")
        need(direct >= 648 * data.hamilton**2, "terminal-C3 sharp bound")
        need(
            (direct == 648 * data.hamilton**2) == (len(data.out) == 1),
            "terminal-C3 equality classification",
        )
        c3_identity_checks += 1

        for left in transitive_factors:
            left_order = len(left.out)
            local = theta(left, data, cycle)
            bound = 648 * left_order * data.hamilton**2
            need(local >= bound, "transitive terminal-C3 local bound")
            c3_telescoping_checks += 1

            terminal = base.ordinal_data(data, cycle)
            value = base.remainder(left, terminal)
            need(value >= bound, "arbitrary-prefix terminal-C3 OS+ bound")
            if value == bound:
                c3_os_equality.append(
                    (left_order, len(data.out), base.label(data.out), value)
                )
            c3_os_checks += 1
    need(c3_os_equality == [(1, 1, "", 648)], "terminal-C3 OS+ equality")

    # Cycle-first/transitive-tail exact negative family and local interactions.
    transitive_long = []
    current = one
    for order in range(1, 65):
        if order > 1:
            current = base.ordinal_data(current, one)
        transitive_long.append(current)
        need(current.out == base.transitive(order), "long transitive recursion")
        value = base.remainder(cycle, current)
        need(value == negative_remainder_formula(order), "negative remainder formula")
        need(value < 0, "cycle-first remainder is not negative")

    negative_theta_checks = 0
    for left_tail in range(1, 13):
        for right_tail in range(1, 13):
            middle = transitive_long[left_tail - 1]
            right = transitive_long[right_tail - 1]
            value = theta(cycle, middle, right)
            need(
                value == negative_theta_formula(left_tail, right_tail),
                "negative interaction formula",
            )
            need(value < 0, "cycle-first interaction is not negative")
            negative_theta_checks += 1
    need(theta(cycle, one, one) == -216, "minimal inherited hostile")

    print("TOURNAMENT_SOURCE_PADDING_SUPERMODULARITY_PRIMARY")
    print("class_counts", " ".join(f"q{n}={len(banks[n])}" for n in range(1, 7)))
    print("padding_formula_checks", padding_checks)
    print("ordered_pair_supermodularity_checks", pair_checks)
    print("ordered_pair_zero", pair_zero)
    print("ordered_pair_minimum_positive", pair_minimum_positive)
    print("ordered_pair_stream_sha256", pair_digest.hexdigest())
    print("transitive_left_remainder_checks", transitive_remainder_checks)
    print("transitive_left_remainder_zeros", transitive_remainder_zeros)
    print("transitive_left_theta_checks", transitive_theta_checks)
    print("transitive_left_theta_zeros", transitive_theta_zeros)
    print("transitive_right_dual_remainder_checks", dual_remainder_checks)
    print("transitive_right_dual_remainder_zeros", dual_remainder_zeros)
    print("transitive_left_no_sink_os_checks", os_checks)
    print("terminal_c3_identity_checks", c3_identity_checks)
    print("terminal_c3_telescoping_checks", c3_telescoping_checks)
    print("terminal_c3_os_checks", c3_os_checks)
    print("terminal_c3_os_equality", c3_os_equality)
    print("cycle_first_negative_remainder_formula_checks", len(transitive_long))
    print("cycle_first_negative_theta_formula_checks", negative_theta_checks)
    print("minimal_cycle_first_theta_hostile", theta(cycle, one, one))
    print("PASS")


if __name__ == "__main__":
    main()
