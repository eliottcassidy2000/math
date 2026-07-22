#!/usr/bin/env python3
"""Exact-rational referee for THM-2085's Selberg/Hunter constants."""

from fractions import Fraction as F


def certificate(degree: int):
    eps = F(1, degree + 1)
    u_guard = F(2, 7) + eps
    u_danger = F(1, 7) + eps
    u_complement = F(5, 7) + eps

    mixed_lower = (
        u_guard * u_danger
        - 2 * eps * (u_danger + u_guard)
    )
    edge_lower = (
        u_complement * u_danger * u_danger
        - 2 * eps * (
            u_danger * u_danger
            + 2 * u_complement * u_danger
        )
    )
    tree_lower = 6 * edge_lower
    deficit_upper = F(2, 7) - 7 * mixed_lower
    margin = tree_lower - deficit_upper
    return eps, mixed_lower, edge_lower, tree_lower, deficit_upper, margin


def main() -> None:
    first_positive = None
    for degree in range(1, 58):
        if certificate(degree)[-1] > 0:
            first_positive = degree
            break

    row56 = certificate(56)
    row57 = certificate(57)

    assert first_positive == 57
    assert row56[-1] == F(-20416, 21173733)
    assert row57[0] == F(1, 58)
    assert row57[1] == F(5363, 164836)
    assert row57[2] == F(655135, 66923416)
    assert row57[-1] == F(6435, 8365427)
    assert row56[-1] <= 0 < row57[-1]

    print("THM-2085 SELBERG RANK-SEVEN RELATION-GATE REFEREE")
    print(f"degrees 1..56 with positive certificate margin: none")
    print(f"degree 56 margin: {row56[-1]}")
    print(f"degree 57 epsilon: {row57[0]}")
    print(f"degree 57 mixed lower: {row57[1]}")
    print(f"degree 57 restricted-edge lower: {row57[2]}")
    print(f"degree 57 tree lower: {row57[3]}")
    print(f"degree 57 deficit upper: {row57[4]}")
    print(f"degree 57 contradiction margin: {row57[-1]}")
    print("ALL EXACT-RATIONAL ASSERTIONS PASS")


if __name__ == "__main__":
    main()
