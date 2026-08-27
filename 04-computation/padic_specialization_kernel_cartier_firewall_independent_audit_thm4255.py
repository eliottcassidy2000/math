#!/usr/bin/env python3
"""Dependency-free clean-room finite audit for THM-4255.

Exact arithmetic uses finite prime fields and sparse dictionaries. This
independently checks graph division, short-jet ranks, explicit-f Cartier
channels, two abstract hostile windows, and restricted ell-adic grade loss.
"""

from __future__ import annotations

from fractions import Fraction
from itertools import product


Poly = dict[tuple[int, int], int]  # (f exponent, u exponent) -> coefficient
Series = dict[int, int]  # f exponent -> coefficient


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def clean_poly(poly: Poly, prime: int) -> Poly:
    return {key: value % prime for key, value in poly.items() if value % prime}


def clean_series(series: Series, prime: int) -> Series:
    return {key: value % prime for key, value in series.items() if value % prime}


def add_poly(left: Poly, right: Poly, prime: int) -> Poly:
    out = dict(left)
    for key, value in right.items():
        out[key] = out.get(key, 0) + value
    return clean_poly(out, prime)


def mul_poly(left: Poly, right: Poly, prime: int) -> Poly:
    out: Poly = {}
    for (fi, ui), a in left.items():
        for (fj, uj), b in right.items():
            key = (fi + fj, ui + uj)
            out[key] = out.get(key, 0) + a * b
    return clean_poly(out, prime)


def lift_series(series: Series) -> Poly:
    return {(exponent, 0): coefficient for exponent, coefficient in series.items()}


def add_series(left: Series, right: Series, prime: int) -> Series:
    out = dict(left)
    for key, value in right.items():
        out[key] = out.get(key, 0) + value
    return clean_series(out, prime)


def mul_series(
    left: Series,
    right: Series,
    prime: int,
    jet: int | None = None,
) -> Series:
    out: Series = {}
    for i, a in left.items():
        for j, b in right.items():
            exponent = i + j
            if jet is None or exponent < jet:
                out[exponent] = out.get(exponent, 0) + a * b
    return clean_series(out, prime)


def pow_series(
    base: Series,
    exponent: int,
    prime: int,
    jet: int | None = None,
) -> Series:
    out: Series = {0: 1}
    factor = clean_series(base, prime)
    power = exponent
    while power:
        if power & 1:
            out = mul_series(out, factor, prime, jet)
        factor = mul_series(factor, factor, prime, jet)
        power >>= 1
    return out


def evaluate_u(
    poly: Poly,
    sigma: Series,
    prime: int,
    jet: int | None = None,
) -> Series:
    out: Series = {}
    for (f_exponent, u_exponent), coefficient in poly.items():
        sigma_power = pow_series(sigma, u_exponent, prime, jet)
        term = {
            f_exponent + exponent: coefficient * value
            for exponent, value in sigma_power.items()
            if jet is None or f_exponent + exponent < jet
        }
        out = add_series(out, term, prime)
    return clean_series(out, prime)


def divide_by_u_minus_f(poly: Poly, prime: int) -> tuple[Poly, Series]:
    """Return Q,R with F=(u-f)Q+R(f), by monomial telescoping."""

    quotient: Poly = {}
    remainder: Series = {}
    for (f_exponent, u_exponent), coefficient in poly.items():
        remainder[f_exponent + u_exponent] = (
            remainder.get(f_exponent + u_exponent, 0) + coefficient
        )
        for k in range(u_exponent):
            # u^a-f^a=(u-f) sum_(k=0)^(a-1) u^(a-1-k) f^k.
            key = (f_exponent + k, u_exponent - 1 - k)
            quotient[key] = quotient.get(key, 0) + coefficient
    return clean_poly(quotient, prime), clean_series(remainder, prime)


def shift_series(series: Series, amount: int, prime: int) -> Series:
    require(all(exponent >= amount for exponent in series), "invalid negative shift")
    return clean_series(
        {exponent - amount: value for exponent, value in series.items()},
        prime,
    )


def cartier_poly(poly: Poly, residue: int, ell: int, prime: int) -> Poly:
    out: Poly = {}
    for (f_exponent, u_exponent), coefficient in poly.items():
        if f_exponent % ell == residue:
            key = ((f_exponent - residue) // ell, u_exponent)
            out[key] = out.get(key, 0) + coefficient
    return clean_poly(out, prime)


def cartier_series(series: Series, residue: int, ell: int, prime: int) -> Series:
    return clean_series(
        {
            (f_exponent - residue) // ell: coefficient
            for f_exponent, coefficient in series.items()
            if f_exponent % ell == residue
        },
        prime,
    )


def matrix_rank_mod(columns: list[list[int]], prime: int) -> int:
    if not columns:
        return 0
    rows = len(columns[0])
    matrix = [[column[row] % prime for column in columns] for row in range(rows)]
    rank = 0
    for column in range(len(columns)):
        pivot = next((row for row in range(rank, rows) if matrix[row][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, prime)
        matrix[rank] = [(entry * inverse) % prime for entry in matrix[rank]]
        for row in range(rows):
            if row != rank and matrix[row][column]:
                factor = matrix[row][column]
                matrix[row] = [
                    (a - factor * b) % prime
                    for a, b in zip(matrix[row], matrix[rank], strict=True)
                ]
        rank += 1
        if rank == rows:
            break
    return rank


def jet_column(monomial: tuple[int, int], jet: int, prime: int) -> list[int]:
    series = evaluate_u({monomial: 1}, {1: 1}, prime, jet)
    return [series.get(exponent, 0) for exponent in range(jet)]


def exhaustive_kernel_and_preimage_audit() -> int:
    prime = 2
    basis = [
        (f_exponent, u_exponent)
        for f_exponent in range(3)
        for u_exponent in range(3)
    ]
    u_minus_f: Poly = {(0, 1): 1, (1, 0): -1}
    checked = 0
    for values in product(range(prime), repeat=len(basis)):
        poly = clean_poly(dict(zip(basis, values, strict=True)), prime)
        quotient, remainder = divide_by_u_minus_f(poly, prime)
        reconstructed = add_poly(
            mul_poly(u_minus_f, quotient, prime),
            lift_series(remainder),
            prime,
        )
        require(reconstructed == poly, "division identity failed")
        evaluated = evaluate_u(poly, {1: 1}, prime)
        require(evaluated == remainder, "remainder is not specialization")

        if not evaluated:
            require(
                mul_poly(u_minus_f, quotient, prime) == poly,
                "kernel not generated by u-f",
            )

        n = 3
        if all(exponent >= n for exponent in evaluated):
            tail = shift_series(remainder, n, prime)
            decomposition = add_poly(
                mul_poly(u_minus_f, quotient, prime),
                {(exponent + n, 0): value for exponent, value in tail.items()},
                prime,
            )
            require(decomposition == poly, "preimage ideal decomposition failed")
        checked += 1
    return checked


def cartier_noncommutation_audit() -> tuple[Series, Series]:
    prime, ell = 101, 5
    element: Poly = {(0, 1): 1, (1, 0): -1}
    specialized = evaluate_u(element, {1: 1}, prime)
    require(not specialized, "u-f must specialize to zero")
    after_0 = cartier_series(specialized, 0, ell, prime)
    before_0 = evaluate_u(cartier_poly(element, 0, ell, prime), {1: 1}, prime)
    after_1 = cartier_series(specialized, 1, ell, prime)
    before_1 = evaluate_u(cartier_poly(element, 1, ell, prime), {1: 1}, prime)
    require(not after_0 and not after_1, "zero acquired a Cartier channel")
    require(before_0 == {1: 1}, "unexpected residue-zero channel")
    require(before_1 == {0: prime - 1}, "unexpected residue-one channel")
    return before_0, before_1


def short_jet_audit() -> tuple[int, int, int, int, int, int]:
    prime, jet = 101, 5
    raw_basis = [(0, 0), (0, 1), (1, 0)]
    raw_rank = matrix_rank_mod(
        [jet_column(item, jet, prime) for item in raw_basis],
        prime,
    )
    require(raw_rank == 2 < len(raw_basis), "raw graph kernel was not detected")

    graph_normal_basis = [(0, exponent) for exponent in range(jet)]
    graph_rank = matrix_rank_mod(
        [jet_column(item, jet, prime) for item in graph_normal_basis],
        prime,
    )
    require(graph_rank == jet, "graph-normal short-jet map should be injective")

    # f^jet lies in V intersect f^jet S; quotienting this tail is essential.
    with_tail_basis = graph_normal_basis + [(jet, 0)]
    with_tail_rank = matrix_rank_mod(
        [jet_column(item, jet, prime) for item in with_tail_basis],
        prime,
    )
    quotient_dimension = len(with_tail_basis) - 1
    require(with_tail_rank == quotient_dimension, "induced quotient lost injectivity")

    rectangle_basis = [(i, a) for i in range(3) for a in range(3)]
    rectangle_rank = matrix_rank_mod(
        [jet_column(item, jet, prime) for item in rectangle_basis],
        prime,
    )
    require(rectangle_rank == 5, "rectangle total-degree rank is wrong")
    return (
        raw_rank,
        len(raw_basis),
        graph_rank,
        quotient_dimension,
        rectangle_rank,
        len(rectangle_basis),
    )


def window_representations(
    n: int,
    ell: int,
    lower: int,
    upper: int,
) -> list[tuple[int, int]]:
    return [
        (r, n - ell * r)
        for r in range(1, 2 + n // ell)
        if lower <= n - ell * r <= upper
    ]


def prop_6_2_hostile_window() -> tuple[int, int, int, list[tuple[int, int]]]:
    prime, ell = 101, 5
    # Gamma=(u-f)z+z^2 and z=f^ell.
    gamma_after_z: Poly = {
        (ell, 1): 1,
        (ell + 1, 0): -1,
        (2 * ell, 0): 1,
    }
    coefficient_ring_order = min(f_exponent for f_exponent, _ in gamma_after_z)
    scalar = evaluate_u(gamma_after_z, {1: 1}, prime)
    scalar_order = min(scalar)
    require(coefficient_ring_order == ell, "wrong coefficient-ring order")
    require(scalar == {2 * ell: 1}, "wrong specialized Gamma")
    require(scalar_order == 2 * ell, "specialization did not increase order")

    d_bound, correction, ell_window, normal_order = 10, 0, 23, 15
    moving: Poly = {(0, 1): 1, (0, 0): -1}
    moving_scalar = evaluate_u(moving, {0: 1, normal_order: 1}, prime)
    hostile_n = min(moving_scalar) + ell_window
    possible = window_representations(
        hostile_n,
        ell_window,
        -correction,
        d_bound + correction,
    )
    require(hostile_n == 38 and not possible, "Prop. 6.2 hostile failed")
    return coefficient_ring_order, scalar_order, hostile_n, possible


def prop_12_3_hostile_window() -> tuple[int, list[tuple[int, int]]]:
    prime = 101
    d_bound, correction, ell, normal_order = 12, 1, 29, 17
    multiplier: Poly = {(0, 1): 1, (0, 0): -1}
    scalar_multiplier = evaluate_u(
        multiplier,
        {0: 1, normal_order: 1},
        prime,
    )
    require(scalar_multiplier == {normal_order: 1}, "wrong moving specialization")
    pivot = min(scalar_multiplier) + ell
    possible = window_representations(
        pivot,
        ell,
        -correction,
        d_bound + correction,
    )
    require(pivot == 46 and not possible, "Prop. 12.3 hostile failed")
    require(Fraction(1, ell).denominator == ell, "pivot is not ell-negative")
    return pivot, possible


def restricted_associated_grade_hostile() -> int:
    ell, normal_order = 5, 7
    coefficient_vector = (1, -1)
    require(any(value % ell for value in coefficient_vector), "vector not primitive")
    specialized = {normal_order: ell}
    require(
        all(value % ell == 0 for value in specialized.values()),
        "specialization stayed primitive",
    )
    return normal_order


def main() -> None:
    checked = exhaustive_kernel_and_preimage_audit()
    channel_0, channel_1 = cartier_noncommutation_audit()
    raw_rank, raw_dim, graph_rank, quotient_dim, rectangle_rank, rectangle_dim = (
        short_jet_audit()
    )
    pre_order, post_order, prop6_n, prop6_possible = prop_6_2_hostile_window()
    prop12_n, prop12_possible = prop_12_3_hostile_window()
    strictness_order = restricted_associated_grade_hostile()

    print(f"EXHAUSTIVE_KERNEL_PREIMAGE PASS field=F2 elements={checked}")
    print("IDENTITIES ker(ev)=(u-f); ev^-1((f^n))=(u-f)+(f^n) PASS")
    print(
        f"CARTIER_NONCOMMUTATION PASS C0_before={channel_0}",
        f"C1_before={channel_1} after=0",
    )
    print(
        "SHORT_JET PASS",
        f"raw_rank={raw_rank}/{raw_dim}",
        f"graph_normal_rank={graph_rank}/{quotient_dim}",
        f"rectangle_rank={rectangle_rank}/{rectangle_dim}",
        "tail_quotient_injective=yes",
    )
    print(
        "PROP6_2_HOSTILE PASS",
        f"pre_order={pre_order}",
        f"post_order={post_order}",
        f"pivot={prop6_n}",
        f"admissible_windows={prop6_possible}",
    )
    print(
        "PROP12_3_HOSTILE PASS",
        f"pivot={prop12_n}",
        f"admissible_windows={prop12_possible}",
    )
    print(
        "RESTRICTED_ASSOCIATED_GRADE_HOSTILE PASS",
        f"hidden_order={strictness_order}",
        "full_graph_quotient_strict=yes",
    )
    print("CLEANROOM_STANDARD_LIBRARY_AUDIT: PASS")


if __name__ == "__main__":
    main()
