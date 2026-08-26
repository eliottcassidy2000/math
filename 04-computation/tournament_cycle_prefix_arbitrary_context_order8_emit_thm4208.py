#!/usr/bin/env python3
"""Emit exact response jets for the THM-4208 factor-class census."""

from __future__ import annotations

import tournament_ordinal_cocycle_parity_thm4184 as base


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


def left_jet(data: base.TournamentData) -> tuple[int, ...]:
    degree, outgoing, _ = directed_masses(data)
    linear = [data.mass - degree[i] + 4 * outgoing[i] for i in range(len(data.out))]
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


def right_jet(data: base.TournamentData) -> tuple[int, ...]:
    degree, _, incoming = directed_masses(data)
    linear = [data.mass - degree[i] - 4 * incoming[i] for i in range(len(data.out))]
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


def main() -> None:
    base.EXPECTED_CLASSES.update({7: 456, 8: 6880})
    executable = base.find_gentourng()
    banks = {
        1: (base.tournament_data(base.transitive(1)),),
        **{
            order: tuple(
                base.tournament_data(out) for out in base.classes(executable, order)
            )
            for order in range(2, 9)
        },
    }
    expected = (1, 1, 2, 4, 12, 56, 456, 6880)
    if tuple(len(banks[order]) for order in range(1, 9)) != expected:
        raise RuntimeError("class count mismatch")
    cycle = base.tournament_data(base.parse("101", 3))
    a5 = base.ordinal_data(cycle, base.tournament_data(base.transitive(5)))
    bank = tuple((order, data) for order in range(1, 9) for data in banks[order])
    print(len(bank))
    for index, (order, data) in enumerate(bank):
        prefixed = base.ordinal_data(a5, data)
        left = tuple(x - y for x, y in zip(left_jet(prefixed), left_jet(data)))
        right = right_jet(data)
        print(
            index,
            order,
            base.label(data.out) or "-",
            data.hamilton,
            *left,
            *right,
        )


if __name__ == "__main__":
    main()
