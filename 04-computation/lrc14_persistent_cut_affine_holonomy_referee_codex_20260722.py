#!/usr/bin/env python3
"""Exact referee for THM-2089's signed affine-cycle identities."""

from fractions import Fraction as F
from itertools import combinations
from random import Random


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def compose(transports):
    """Return A,B for s_out=A*s_in+B*h."""
    A, B = F(1), F(0)
    for alpha, beta in transports:
        A, B = alpha * A, alpha * B + beta
    return A, B


def cycle_constants(rows):
    D = 1
    N = (-1) ** len(rows)
    for _, b, c in rows:
        D *= c
        N *= b
    R = 0
    for i, (a, _, _) in enumerate(rows):
        left_c = 1
        right_b = 1
        for j in range(i):
            left_c *= rows[j][2]
        for j in range(i + 1, len(rows)):
            right_b *= -rows[j][1]
        R += (-a) * right_b * left_c
    return D, N, R


def directed_transport(edge_data, i, j):
    alpha, beta = edge_data[tuple(sorted((i, j)))]
    if i < j:
        return alpha, beta
    return 1 / alpha, -beta / alpha


def flat_square_count(a_size, b_size, perturb=False):
    left = list(range(a_size))
    right = list(range(a_size, a_size + b_size))
    u = {i: F((i + 2) * (-1 if i % 3 == 0 else 1), i + 1) for i in left + right}
    v = {i: F(i * i - 3 * i + 1, i + 2) for i in left + right}
    data = {}
    for i in left:
        for j in right:
            alpha = u[j] / u[i]
            beta = u[j] * (v[j] - v[i])
            data[(i, j)] = (alpha, beta)
    if perturb and a_size >= 2 and b_size >= 2:
        edge = (left[0], right[0])
        alpha, beta = data[edge]
        data[edge] = (alpha, beta + F(1, 101))

    flat = 0
    nonflat = 0
    for i1, i2 in combinations(left, 2):
        for j1, j2 in combinations(right, 2):
            walk = [(i1, j1), (j1, i2), (i2, j2), (j2, i1)]
            A, B = compose(directed_transport(data, x, y) for x, y in walk)
            if (A, B) == (1, 0):
                flat += 1
            else:
                nonflat += 1
    return flat, nonflat


def main() -> None:
    rng = Random(2089)
    random_checks = 0
    for r in (2, 4, 6):
        for _ in range(500):
            rows = []
            for _ in range(r):
                a = rng.randint(-57, 57)
                b = rng.choice([x for x in range(-57, 58) if x])
                c = rng.choice([x for x in range(-57, 58) if x])
                rows.append((a, b, c))
            transports = [(F(-b, c), F(-a, c)) for a, b, c in rows]
            A, B = compose(transports)
            D, N, R = cycle_constants(rows)
            require(A == F(N, D), "multiplicative clearing mismatch")
            require(B == F(R, D), "offset clearing mismatch")
            random_checks += 1

    flat_25 = flat_square_count(2, 5)
    flat_34 = flat_square_count(3, 4)
    perturbed_34 = flat_square_count(3, 4, perturb=True)
    require(flat_25 == (10, 0), "K2,5 flat square count")
    require(flat_34 == (18, 0), "K3,4 flat square count")
    require(perturbed_34[1] > 0, "perturbed connection stayed flat")
    require(6 * 57**6 == 205_778_683_494, "height constant mismatch")

    print("THM-2089 PERSISTENT-CUT AFFINE-HOLONOMY REFEREE")
    print(f"random cleared-cycle identities: {random_checks}")
    print(f"flat K2,5 squares (flat,nonflat): {flat_25}")
    print(f"flat K3,4 squares (flat,nonflat): {flat_34}")
    print(f"perturbed K3,4 squares (flat,nonflat): {perturbed_34}")
    print(f"universal simple-cycle relation-height invoice: {6 * 57**6}")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
