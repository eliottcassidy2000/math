#!/usr/bin/env python3
"""Exact primary audit for THM-4208.

The factor universe is every labelled tournament through the declared order.
Factor data are rebuilt by the exact subset/exposed-word engine of THM-4184;
ordinal children use the proved THM-4181/4184 transfer.  The audit checks the
closed form against literal remainders, not merely against its recurrence.
"""

from __future__ import annotations

import hashlib
import itertools
import math

import tournament_ordinal_cocycle_parity_thm4184 as base


def need(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def labelled(order: int) -> tuple[base.TournamentData, ...]:
    width = order * (order - 1) // 2
    answer = []
    for integer in range(1 << width):
        bits = "" if width == 0 else format(integer, f"0{width}b")
        answer.append(base.tournament_data(base.parse(bits, order)))
    return tuple(answer)


def is_transitive(data: base.TournamentData) -> bool:
    n = len(data.out)
    return not any(
        data.out[i] & (1 << j)
        and data.out[j] & (1 << k)
        and data.out[k] & (1 << i)
        for i in range(n)
        for j in range(n)
        for k in range(n)
        if len({i, j, k}) == 3
    )


def is_strong(data: base.TournamentData) -> bool:
    full = (1 << len(data.out)) - 1
    for start in range(len(data.out)):
        reached = 1 << start
        previous = -1
        while reached != previous:
            previous = reached
            for vertex in range(len(data.out)):
                if reached & (1 << vertex):
                    reached |= data.out[vertex]
        if reached != full:
            return False
    return True


def literal_hamilton_endpoint_counts(
    data: base.TournamentData,
) -> tuple[list[int], list[int]]:
    """Count Hamilton starts/ends by literal permutation enumeration."""
    order = len(data.out)
    starts = [0] * order
    ends = [0] * order
    for path in itertools.permutations(range(order)):
        if all(data.out[path[i]] & (1 << path[i + 1]) for i in range(order - 1)):
            starts[path[0]] += 1
            ends[path[-1]] += 1
    return starts, ends


def directed_masses(data: base.TournamentData) -> tuple[list[int], list[int], list[int]]:
    degree = [sum(row) for row in data.capacities]
    outgoing = [
        sum(
            data.capacities[i][j]
            for j in range(len(data.out))
            if data.out[i] & (1 << j)
        )
        for i in range(len(data.out))
    ]
    incoming = [degree[i] - outgoing[i] for i in range(len(data.out))]
    return degree, outgoing, incoming


def signature(data: base.TournamentData) -> tuple[int, ...]:
    h = data.hamilton
    w = data.mass
    degree, _, incoming = directed_masses(data)
    ell = [w - degree[i] - 4 * incoming[i] for i in range(len(data.out))]
    l0 = sum(data.ends[i][0] * ell[i] for i in range(len(data.out)))
    l1 = sum(data.ends[i][1] * ell[i] for i in range(len(data.out)))
    m00 = sum(row[0] * row[0] for row in data.ends)
    m01 = sum(row[0] * row[1] for row in data.ends)
    m11 = sum(row[1] * row[1] for row in data.ends)
    return h, w, l0, l1, m00, m01, m11


def sequence_coefficients(data: base.TournamentData) -> tuple[int, int, int, int]:
    h, w, l0, l1, m00, m01, m11 = signature(data)
    alpha = 117 * (w + h) * (w + 3 * h) - 345 * (m00 + 2 * m01 + m11)
    beta = -18 * h * (w + h)
    gamma = (
        36 * (l0 + l1)
        + 360 * (m00 + m01)
        - 36 * w * w
        - 126 * h * w
        - 108 * h * h
    )
    delta = 12 * (m11 - m01 - 8 * m00) - 18 * (h * w + l0)
    return alpha, beta, gamma, delta


def unary_formula(tail: int, data: base.TournamentData) -> int:
    alpha, beta, gamma, delta = sequence_coefficients(data)
    return 4**tail * alpha + tail * 2**tail * beta + 2**tail * gamma + delta


def endpoint_correction(data: base.TournamentData) -> int:
    h, _, _, _, m00, m01, m11 = signature(data)
    return 3 * (7 * h * h - m00 + 2 * m01 - m11)


def endpoint_data(data: base.TournamentData) -> tuple[int, int, int, int, list[int], list[int]]:
    h = data.hamilton
    values = [row[0] + row[1] for row in data.ends]
    a = sum(values)
    b = data.mass + 2 * h
    moment = sum(value * value for value in values)
    defect = a * (a + 2 * h) - 3 * moment
    slacks = [
        h
        + sum(values[j] for j in range(len(values)) if data.out[j] & (1 << i))
        - values[i]
        for i in range(len(values))
    ]
    return a, b, moment, defect, values, slacks


def startpoint_data(data: base.TournamentData) -> tuple[int, int, int, int]:
    """Converse-dual U-energy data."""
    h = data.hamilton
    values = [row[0] + row[1] for row in data.starts]
    a = sum(values)
    b = data.mass + 2 * h
    moment = sum(value * value for value in values)
    defect = a * (a + 2 * h) - 3 * moment
    return a, b, moment, defect


def left_response_jet(data: base.TournamentData) -> tuple[int, ...]:
    degree, outgoing, _ = directed_masses(data)
    linear = [
        data.mass - degree[i] + 4 * outgoing[i]
        for i in range(len(data.out))
    ]
    s0, s1 = data.mass // 2, data.mass // 2 + data.hamilton
    q00 = sum(row[0] * row[0] for row in data.starts)
    q01 = sum(row[0] * row[1] for row in data.starts)
    q11 = sum(row[1] * row[1] for row in data.starts)
    l0 = sum(data.starts[i][0] * linear[i] for i in range(len(data.out)))
    l1 = sum(data.starts[i][1] * linear[i] for i in range(len(data.out)))
    return (
        data.hamilton * data.mass,
        2 * l0,
        2 * l1,
        2 * data.hamilton * s0,
        2 * data.hamilton * s1,
        2 * s0 * s0 + 6 * q00,
        4 * s0 * s1 + 12 * q01,
        2 * s1 * s1 + 6 * q11,
        2 * q00 - 10 * s0 * s0,
        4 * q01 - 20 * s0 * s1,
        2 * q11 - 10 * s1 * s1,
    )


def right_response_jet(data: base.TournamentData) -> tuple[int, ...]:
    degree, _, incoming = directed_masses(data)
    linear = [
        data.mass - degree[i] - 4 * incoming[i]
        for i in range(len(data.out))
    ]
    s0, s1 = data.mass // 2, data.mass // 2 + data.hamilton
    q00 = sum(row[0] * row[0] for row in data.ends)
    q01 = sum(row[0] * row[1] for row in data.ends)
    q11 = sum(row[1] * row[1] for row in data.ends)
    l0 = sum(data.ends[i][0] * linear[i] for i in range(len(data.out)))
    l1 = sum(data.ends[i][1] * linear[i] for i in range(len(data.out)))
    return (
        data.hamilton * data.mass,
        data.hamilton * s0,
        data.hamilton * s1,
        l0,
        l1,
        s0 * s0,
        s0 * s1,
        s1 * s1,
        q00,
        q01,
        q11,
    )


def dot(left: tuple[int, ...], right: tuple[int, ...]) -> int:
    return sum(x * y for x, y in zip(left, right))


def lollipop_formula(tail: int) -> int:
    return 162 * (36 * 4**tail - (3 * tail + 20) * 2**tail + 4)


def main() -> None:
    factors5 = tuple(data for n in range(1, 6) for data in labelled(n))
    factors4 = tuple(data for n in range(1, 5) for data in labelled(n))
    need(len(factors5) == 1099, "labelled factor count through five")
    need(len(factors4) == 75, "labelled factor count through four")

    one = base.tournament_data(base.transitive(1))
    cycle = base.tournament_data(base.parse("101", 3))
    prefixes = [cycle] + [
        base.ordinal_data(cycle, base.tournament_data(base.transitive(n)))
        for n in range(1, 9)
    ]

    cycle_tail_state_checks = 0
    for tail in range(1, 9):
        data = prefixes[tail]
        power = 2**tail
        degree, outgoing, _ = directed_masses(data)
        q00 = sum(row[0] * row[0] for row in data.starts)
        q01 = sum(row[0] * row[1] for row in data.starts)
        q11 = sum(row[1] * row[1] for row in data.starts)
        linear = [
            sum(
                data.starts[i][parity]
                * (data.mass - degree[i] + 4 * outgoing[i])
                for i in range(len(data.out))
            )
            for parity in range(2)
        ]
        need((data.hamilton, data.mass) == (3, 12 * power - 6), "A_n mass state")
        need(
            (q00, q01, q11)
            == (
                15 * power * power // 2 - 3,
                15 * power * power // 2 - 3,
                15 * power * power // 2 + 6,
            ),
            "A_n rooted Gram state",
        )
        need(
            linear
            == [
                117 * power * power - (9 * tail + 72) * power,
                117 * power * power - (9 * tail + 54) * power - 18,
            ],
            "A_n rooted linear state",
        )
        cycle_tail_state_checks += 1

    digest = hashlib.sha256()
    unary_checks = 0
    moment_checks = 0
    dual_moment_checks = 0
    rooted_chirality_checks = 0
    fan_collapse_checks = 0
    rooted_square_equalities = 0
    zero_defects = 0
    for data in factors5:
        hamilton_starts, hamilton_ends = literal_hamilton_endpoint_counts(data)
        need(sum(hamilton_starts) == data.hamilton, "literal Hamilton start total")
        need(sum(hamilton_ends) == data.hamilton, "literal Hamilton end total")
        need(
            [row[1] - row[0] for row in data.ends] == hamilton_starts,
            "coordinatewise V chirality identity",
        )
        need(
            [row[1] - row[0] for row in data.starts] == hamilton_ends,
            "coordinatewise U chirality identity",
        )
        _, outgoing_fan, incoming_fan = directed_masses(data)
        need(
            [2 * row[0] for row in data.ends] == incoming_fan,
            "V-even coordinate is half the incoming fan",
        )
        need(
            [2 * row[1] for row in data.ends]
            == [incoming_fan[i] + 2 * hamilton_starts[i] for i in range(len(data.out))],
            "V-odd coordinate is incoming fan plus Hamilton starts",
        )
        need(
            [2 * row[0] for row in data.starts] == outgoing_fan,
            "U-even coordinate is half the outgoing fan",
        )
        need(
            [2 * row[1] for row in data.starts]
            == [outgoing_fan[i] + 2 * hamilton_ends[i] for i in range(len(data.out))],
            "U-odd coordinate is outgoing fan plus Hamilton ends",
        )
        need(
            sum((row[0] + row[1]) ** 2 for row in data.ends)
            == sum((incoming_fan[i] + hamilton_starts[i]) ** 2 for i in range(len(data.out))),
            "V-energy fan/endpoint collapse",
        )
        need(
            sum((row[0] + row[1]) ** 2 for row in data.starts)
            == sum((outgoing_fan[i] + hamilton_ends[i]) ** 2 for i in range(len(data.out))),
            "U-energy fan/endpoint collapse",
        )
        start_square = sum(value * value for value in hamilton_starts)
        universal_source = any(
            data.out[i].bit_count() == len(data.out) - 1
            for i in range(len(data.out))
        )
        need(start_square <= data.hamilton**2, "Hamilton-start square bound")
        need(
            (start_square == data.hamilton**2) == universal_source,
            "Hamilton-start square equality boundary",
        )
        need(
            endpoint_correction(data) == 3 * (7 * data.hamilton**2 - start_square),
            "positive n=0 correction from Hamilton starts",
        )
        need(
            endpoint_correction(data) >= 18 * data.hamilton**2,
            "strict n=0 correction bound",
        )
        rooted_chirality_checks += 1
        fan_collapse_checks += 1
        rooted_square_equalities += start_square == data.hamilton**2

        defect_data = endpoint_data(data)
        _, _, _, defect, values, slacks = defect_data
        need(all(value > 0 for value in values), "positive rooted endpoint mass")
        need(all(slack >= 0 for slack in slacks), "coordinate endpoint slack")
        need(
            defect == 2 * sum(v * s for v, s in zip(values, slacks)),
            "defect/slack identity",
        )
        need((defect == 0) == is_transitive(data), "moment equality boundary")
        zero_defects += defect == 0
        moment_checks += 1
        start_a, start_b, _, start_defect = startpoint_data(data)
        need((start_a, start_b) == (defect_data[0], defect_data[1]),
             "dual rooted mass agreement")
        need((start_defect == 0) == is_transitive(data),
             "dual moment equality boundary")
        need(start_defect >= 0, "dual endpoint energy")
        dual_moment_checks += 1

        direct0 = base.remainder(cycle, data)
        need(
            direct0 == unary_formula(0, data) + endpoint_correction(data),
            "unary n=0 correction",
        )
        for tail in range(1, 9):
            direct = base.remainder(prefixes[tail], data)
            need(direct == unary_formula(tail, data), "unary closed form")
            digest.update(f"U|{tail}|{base.label(data.out)}|{direct}\n".encode())
            unary_checks += 1

    context_pairs = 0
    context_rows = 0
    recurrence_checks = 0
    correction_checks = 0
    response_jet_checks = 0
    ordinal_moment_checks = 0
    ordinal_dual_moment_checks = 0
    leading_factor_checks = 0
    for left in factors4:
        left_sequence = sequence_coefficients(left)
        left_endpoint = endpoint_data(left)
        for right in factors4:
            child = base.ordinal_data(left, right)
            child_sequence = sequence_coefficients(child)
            right_sequence = sequence_coefficients(right)
            contextual = [
                child_sequence[i] - right.hamilton**2 * left_sequence[i]
                for i in range(4)
            ]
            contextual[3] += 8 * base.remainder(left, right)
            need(
                endpoint_correction(child)
                == right.hamilton**2 * endpoint_correction(left),
                "n=0 correction is right-homogeneous",
            )
            correction_checks += 1

            values = []
            for tail in range(9):
                direct = (
                    base.remainder(base.ordinal_data(prefixes[tail], left), right)
                    - base.remainder(left, right)
                )
                formula = (
                    4**tail * contextual[0]
                    + tail * 2**tail * contextual[1]
                    + 2**tail * contextual[2]
                    + contextual[3]
                )
                need(direct == formula, "context formula including n=0")
                digest.update(
                    f"F|{tail}|{base.label(left.out)}|{base.label(right.out)}|{direct}\n".encode()
                )
                values.append(direct)
                context_rows += 1
            for index in range(5):
                need(
                    values[index + 4]
                    == 9 * values[index + 3]
                    - 28 * values[index + 2]
                    + 36 * values[index + 1]
                    - 16 * values[index],
                    "C-finite recurrence from n=0",
                )
                recurrence_checks += 1
            context_pairs += 1

            need(
                dot(left_response_jet(left), right_response_jet(right))
                == base.remainder(left, right),
                "exact 11-jet response",
            )
            response_jet_checks += 1

            a_left, b_left, m_left, e_left, _, _ = left_endpoint
            a_right, b_right, m_right, e_right, _, _ = endpoint_data(right)
            a_child, b_child, m_child, e_child, _, _ = endpoint_data(child)
            need(
                a_child == right.hamilton * a_left + b_left * a_right,
                "endpoint mass ordinal law",
            )
            need(b_child == b_left * b_right, "optional mass ordinal law")
            need(
                m_child == right.hamilton**2 * m_left + b_left**2 * m_right,
                "endpoint moment ordinal law",
            )
            need(
                e_child == right.hamilton**2 * e_left + b_left**2 * e_right,
                "endpoint defect ordinal law",
            )
            ordinal_moment_checks += 1
            _, _, mu_left, eu_left = startpoint_data(left)
            _, b_start_right, mu_right, eu_right = startpoint_data(right)
            _, _, mu_child, eu_child = startpoint_data(child)
            need(
                mu_child
                == b_start_right**2 * mu_left + left.hamilton**2 * mu_right,
                "dual endpoint moment ordinal law",
            )
            need(
                eu_child
                == b_start_right**2 * eu_left + left.hamilton**2 * eu_right,
                "dual endpoint defect ordinal law",
            )
            ordinal_dual_moment_checks += 1
            need(
                contextual[0] == b_left**2 * right_sequence[0] > 0,
                "positive contextual 4^n coefficient",
            )
            need(
                contextual[1]
                == -18
                * left.hamilton
                * right.hamilton
                * b_left
                * a_right
                < 0,
                "negative contextual n2^n coefficient",
            )
            leading_factor_checks += 1

    need(
        zero_defects == sum(math.factorial(n) for n in range(1, 6)),
        "transitive labelled equality count",
    )

    lollipop_checks = 0
    prefix = one
    lollipop_first = []
    for tail in range(1, 65):
        data = base.ordinal_data(prefix, cycle)
        expected_mass = 12 * 2**tail - 6
        expected_square = 270 * 4 ** (tail - 1) + 9
        expected_linear = -9 * 4**tail - (27 * tail + 180) * 2**tail + 108
        degree, _, incoming = directed_masses(data)
        actual_square = sum((v0 + 2 * v1) ** 2 for v0, v1 in data.ends)
        actual_linear = sum(
            (data.ends[i][0] + 2 * data.ends[i][1])
            * (data.mass - degree[i] - 4 * incoming[i])
            for i in range(len(data.out))
        )
        need(data.mass == expected_mass, "lollipop mass")
        need(actual_square == expected_square, "lollipop endpoint square")
        need(actual_linear == expected_linear, "lollipop linear moment")
        actual = base.remainder(cycle, data)
        need(actual == lollipop_formula(tail), "C3-lollipop remainder")
        need(36 * 2**tail > 3 * tail + 20, "C3-lollipop positivity")
        if tail <= 6:
            lollipop_first.append(actual)
        prefix = base.ordinal_data(prefix, one)
        lollipop_checks += 1
    need(base.remainder(cycle, cycle) == 3204, "lollipop a=0 boundary")

    hostile_labels = ("110011110111101", "110001101110111")
    hostile = tuple(
        base.tournament_data(base.parse(label, 6)) for label in hostile_labels
    )
    need(
        [(data.hamilton, data.mass, data.gate) for data in hostile]
        == [(43, 338, 22884), (43, 338, 22884)],
        "coarse hostile fibre",
    )
    need(all(not base.has_sink(data.out) and is_strong(data) for data in hostile),
         "hostile is strong and no-sink")
    permutation = (2, 0, 5, 4, 3, 1)
    need(
        all(
            bool(hostile[0].out[i] & (1 << j))
            != bool(hostile[1].out[permutation[i]] & (1 << permutation[j]))
            for i in range(6)
            for j in range(6)
            if i != j
        ),
        "hostile converse isomorphism",
    )
    need(
        all(
            hostile[0].capacities[i][j]
            == hostile[1].capacities[permutation[i]][permutation[j]]
            for i in range(6)
            for j in range(6)
        ),
        "hostile capacity covariance",
    )
    hostile_responses = tuple(base.remainder(cycle, data) for data in hostile)
    need(hostile_responses == (9377256, 9389016), "chirality response split")

    print("TOURNAMENT_CYCLE_PREFIX_ARBITRARY_CONTEXT_PRIMARY")
    print("PASS")
    print("factor_universe_orders=1..5 labelled_count=1099")
    print(f"unary_rows_n=1..8 count={unary_checks}")
    print(f"cycle_tail_state_checks={cycle_tail_state_checks}")
    print(f"moment_rows count={moment_checks} zero_exactly_transitive={zero_defects}")
    print(f"dual_moment_rows count={dual_moment_checks}")
    print(
        "rooted_chirality_rows="
        f"{rooted_chirality_checks} square_equality_universal_source="
        f"{rooted_square_equalities}"
    )
    print(f"fan_endpoint_collapse_rows={fan_collapse_checks}")
    print("context_universe_each_factor_orders=1..4 labelled_count=75")
    print(f"ordered_context_pairs={context_pairs} context_rows_n=0..8={context_rows}")
    print(f"n0_correction_checks={correction_checks} recurrence_checks={recurrence_checks}")
    print(f"response_jet_checks={response_jet_checks}")
    print(f"ordinal_moment_checks={ordinal_moment_checks}")
    print(f"ordinal_dual_moment_checks={ordinal_dual_moment_checks}")
    print(f"leading_factor_checks={leading_factor_checks}")
    print(f"semantic_sha256={digest.hexdigest()}")
    print(f"lollipop_a_range=1..64 checks={lollipop_checks} a0=3204")
    print("lollipop_first_values=" + ",".join(map(str, lollipop_first)))
    print(f"hostile_labels={hostile_labels}")
    print(f"hostile_responses={hostile_responses}")


if __name__ == "__main__":
    main()
