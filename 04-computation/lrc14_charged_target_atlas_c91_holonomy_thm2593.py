#!/usr/bin/env python3
"""Exact companion for THM-2593.

This is a dependency-free finite check of the additive C91 mapping torus and
the multiplicative unit local system carried by the thirteen THM-2585 slice
factors.  Arithmetic is entirely in F_13[z]/(Phi_7).
"""

from collections import Counter


P = 13
Q = 7
CHECKS = 0

Y_RAW = (
    (0, 9, 9, 0, 0, 4, 4),
    (5, 9, 7, 11, 11, 4, 9),
    (3, 11, 7, 10, 10, 11, 1),
    (9, 11, 4, 2, 10, 5, 10),
    (11, 5, 4, 10, 12, 11, 8),
    (11, 1, 11, 2, 8, 1, 1),
    (6, 9, 12, 8, 4, 4, 7),
    (7, 6, 9, 9, 5, 1, 4),
    (2, 12, 12, 5, 11, 2, 12),
    (2, 5, 2, 1, 3, 9, 8),
    (4, 3, 8, 3, 11, 9, 2),
    (10, 12, 2, 3, 3, 6, 2),
    (8, 4, 9, 2, 2, 6, 4),
)


def check(condition: bool, message: str = "exact check failed") -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(message)


def zeta_power(exponent: int) -> tuple[int, ...]:
    """zeta_7^exponent in the basis 1,zeta,...,zeta^5 over F_13."""
    exponent %= Q
    if exponent < Q - 1:
        return tuple(int(i == exponent) for i in range(Q - 1))
    return tuple(P - 1 for _ in range(Q - 1))


def add(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    return tuple((x + y) % P for x, y in zip(a, b))


def mul(a: tuple[int, ...], b: tuple[int, ...]) -> tuple[int, ...]:
    coefficients = [0] * Q
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            coefficients[(i + j) % Q] += x * y
    top = coefficients[Q - 1]
    return tuple((coefficients[i] - top) % P for i in range(Q - 1))


ONE = zeta_power(0)


def reduce_raw(y: tuple[int, ...]) -> tuple[int, ...]:
    """Reduce seven coefficients modulo Phi_7=1+...+z^6."""
    return tuple((y[i] - y[Q - 1]) % P for i in range(Q - 1))


def multiplication_matrix(a: tuple[int, ...]) -> list[list[int]]:
    columns = [mul(a, zeta_power(j)) for j in range(Q - 1)]
    return [[columns[j][i] for j in range(Q - 1)] for i in range(Q - 1)]


def determinant(matrix: list[list[int]]) -> int:
    a = [[value % P for value in row] for row in matrix]
    out = 1
    for col in range(len(a)):
        pivot = next((row for row in range(col, len(a)) if a[row][col]), None)
        if pivot is None:
            return 0
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
            out = -out
        value = a[col][col] % P
        out = out * value % P
        inverse = pow(value, -1, P)
        for j in range(col, len(a)):
            a[col][j] = a[col][j] * inverse % P
        for row in range(col + 1, len(a)):
            factor = a[row][col]
            for j in range(col, len(a)):
                a[row][j] = (a[row][j] - factor * a[col][j]) % P
    return out % P


def solve(matrix: list[list[int]], target: tuple[int, ...]) -> tuple[int, ...]:
    n = len(matrix)
    a = [[matrix[i][j] % P for j in range(n)] + [target[i] % P] for i in range(n)]
    for col in range(n):
        pivot = next((row for row in range(col, n) if a[row][col]), None)
        check(pivot is not None, "attempted to invert a nonunit")
        if pivot != col:
            a[col], a[pivot] = a[pivot], a[col]
        inverse = pow(a[col][col], -1, P)
        for j in range(col, n + 1):
            a[col][j] = a[col][j] * inverse % P
        for row in range(n):
            if row == col:
                continue
            factor = a[row][col]
            for j in range(col, n + 1):
                a[row][j] = (a[row][j] - factor * a[col][j]) % P
    return tuple(a[i][n] for i in range(n))


def unit_inverse(a: tuple[int, ...]) -> tuple[int, ...]:
    inverse = solve(multiplication_matrix(a), ONE)
    check(mul(a, inverse) == ONE)
    return inverse


def owner_factor(y: tuple[int, ...], kappa: int) -> tuple[int, ...]:
    out = (0,) * (Q - 1)
    for ell, coefficient in enumerate(y):
        term = tuple(coefficient * value % P for value in zeta_power(kappa * ell))
        out = add(out, term)
    return out


def main() -> None:
    y = tuple(reduce_raw(row) for row in Y_RAW)
    check(len(set(y)) == P)

    determinants = tuple(determinant(multiplication_matrix(row)) for row in y)
    check(determinants == (7, 6, 5, 2, 7, 10, 7, 7, 10, 7, 2, 5, 6))
    inverses = tuple(unit_inverse(row) for row in y)

    support_histogram = Counter()
    for q in range(P):
        for kappa in range(1, Q):
            factor = owner_factor(Y_RAW[q], kappa)
            check(any(factor))
            support = 12 * sum(value != 0 for value in factor)
            check(support in (48, 60, 72))
            support_histogram[support] += 1
    check(support_histogram == Counter({48: 8, 60: 32, 72: 38}))

    single_cycles = 0
    additive_edges = 0
    multiplicative_edges = 0
    base_monodromies = 0
    partial_residuals = 0

    for delta in range(1, P):
        # The skew successor is one C_91 cycle.
        state = (0, 0)
        seen = []
        for _ in range(P * Q):
            seen.append(state)
            clock, sheet = state
            next_state = ((clock + 1) % Q, (sheet + delta) % P)

            # h=-sheet kills the pulled-back additive transition delta.
            h_here = -sheet % P
            h_next = -next_state[1] % P
            check((delta + h_next - h_here) % P == 0)
            additive_edges += 1

            # u_sheet=Y_(sheet+delta)Y_sheet^-1 transports the charged section.
            transport = mul(y[next_state[1]], inverses[sheet])
            check(mul(transport, y[sheet]) == y[next_state[1]])
            # The gauge Y_sheet^-1 trivializes that transition on the cover.
            normalized = mul(inverses[next_state[1]], mul(transport, y[sheet]))
            check(normalized == ONE)
            multiplicative_edges += 1
            state = next_state

        check(state == (0, 0))
        check(len(set(seen)) == P * Q)
        single_cycles += 1

        # One base loop has exact order 13 in both local systems.
        for start_sheet in range(P):
            for turns in range(1, P + 1):
                end_sheet = (start_sheet + turns * Q * delta) % P
                additive_holonomy = turns * Q * delta % P
                multiplicative_holonomy = mul(y[end_sheet], inverses[start_sheet])
                expected_trivial = turns == P
                check((additive_holonomy == 0) == expected_trivial)
                check((multiplicative_holonomy == ONE) == expected_trivial)
                base_monodromies += 1

        # If only t of seven edge corrections -delta are supplied, a proper
        # subset cannot kill the old-head cyclic sum because 7<13.
        for t in range(Q + 1):
            residual = (Q - t) * delta % P
            check((residual == 0) == (t == Q))
            partial_residuals += 1
        check(4 * delta % P != 0)  # the three-edge conditional boundary

    # The 78 charged sections repeat over seven clock charts on C_91.
    lifted_histogram = {support: Q * count for support, count in support_histogram.items()}
    check(lifted_histogram == {48: 56, 60: 224, 72: 266})
    check(sum(lifted_histogram.values()) == 546)

    print("THM-2593 charged target-section atlas / C91 holonomy exact companion")
    print("slice_units=13 slice_factors_distinct=13")
    print(f"mapping_tori={single_cycles} single_cycle_size=91")
    print(f"additive_gauge_edges={additive_edges} multiplicative_gauge_edges={multiplicative_edges}")
    print(f"base_loop_monodromy_checks={base_monodromies} exact_order=13")
    print(f"partial_edge_residual_checks={partial_residuals} three_edge_residual=4*delta_nonzero")
    print("charged_cover_vertices=546 support_histogram=48:56,60:224,72:266")
    print("physical_common_carrier=NOT_CONSTRUCTED")
    print(f"ALL_CHECKS=PASS checks={CHECKS}")


if __name__ == "__main__":
    main()
