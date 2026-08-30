#!/usr/bin/env python3
"""Exact finite controls for THM-4289.

The formal proof lives in the theorem file.  This dependency-free audit checks
the function and ambient-Kahler exact sequences on every A23 resolution triple
D_s, the finite cyclic conductor/Jacobian-ideal quotients, the two-branch
endpoint, differentiation compatibility, and the inherited four-channel
observer kernels.  Every gate remains active under ``python -O``.
"""

from __future__ import annotations

from collections.abc import Iterable


P = 1_000_000_007


def check(condition: bool, message: str) -> None:
    """Optimization-safe assertion."""
    if not condition:
        raise RuntimeError(message)


def rank_mod(rows: Iterable[list[int]], prime: int = P) -> int:
    matrix = [[entry % prime for entry in row] for row in rows]
    if not matrix:
        return 0
    width = len(matrix[0])
    pivot_row = 0
    for column in range(width):
        pivot = next(
            (row for row in range(pivot_row, len(matrix)) if matrix[row][column]),
            None,
        )
        if pivot is None:
            continue
        matrix[pivot_row], matrix[pivot] = matrix[pivot], matrix[pivot_row]
        inverse = pow(matrix[pivot_row][column], prime - 2, prime)
        matrix[pivot_row] = [value * inverse % prime for value in matrix[pivot_row]]
        for row in range(len(matrix)):
            if row == pivot_row or matrix[row][column] == 0:
                continue
            scale = matrix[row][column]
            matrix[row] = [
                (left - scale * right) % prime
                for left, right in zip(matrix[row], matrix[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == len(matrix):
            break
    return pivot_row


def in_span(vector: list[int], rows: list[list[int]]) -> bool:
    return rank_mod(rows + [vector]) == rank_mod(rows)


def add_entry(vector: list[int], branch: int, jets: int, degree: int, value: int) -> None:
    if 0 <= degree < jets:
        vector[branch * jets + degree] = (
            vector[branch * jets + degree] + value
        ) % P


def component_product(left: list[int], right: list[int], branches: int, jets: int) -> list[int]:
    out = [0] * (branches * jets)
    for branch in range(branches):
        offset = branch * jets
        for i in range(jets):
            if left[offset + i] == 0:
                continue
            for j in range(jets - i):
                if right[offset + j]:
                    out[offset + i + j] = (
                        out[offset + i + j]
                        + left[offset + i] * right[offset + j]
                    ) % P
    return out


def triple_function(s: int, jets: int, b_degree: int, z_degree: int) -> list[int]:
    """Restriction of b^i z^j to E:b=0, R:z=0, C:z=b^s."""
    out = [0] * (3 * jets)
    if b_degree == 0:
        add_entry(out, 0, jets, z_degree, 1)
    if z_degree == 0:
        add_entry(out, 1, jets, b_degree, 1)
    add_entry(out, 2, jets, b_degree + s * z_degree, 1)
    return out


def triple_db(s: int, jets: int, b_degree: int, z_degree: int) -> list[int]:
    out = [0] * (3 * jets)
    if z_degree == 0:
        add_entry(out, 1, jets, b_degree, 1)
    add_entry(out, 2, jets, b_degree + s * z_degree, 1)
    return out


def triple_dz(s: int, jets: int, b_degree: int, z_degree: int) -> list[int]:
    out = [0] * (3 * jets)
    if b_degree == 0:
        add_entry(out, 0, jets, z_degree, 1)
    add_entry(out, 2, jets, b_degree + s * z_degree + s - 1, s)
    return out


def standard_vector(branches: int, jets: int, branch: int, degree: int, value: int = 1) -> list[int]:
    out = [0] * (branches * jets)
    add_entry(out, branch, jets, degree, value)
    return out


def equal_value_basis(jets: int) -> list[list[int]]:
    """Basis of triples whose three constant values agree."""
    constant = [0] * (3 * jets)
    for branch in range(3):
        add_entry(constant, branch, jets, 0, 1)
    return [constant] + [
        standard_vector(3, jets, branch, degree)
        for branch in range(3)
        for degree in range(1, jets)
    ]


def phi_obstruction(vector: list[int], s: int, jets: int) -> list[int]:
    """Coefficients b^1,...,b^s of c-r-e(b^s)+e(0)."""
    out: list[int] = []
    for degree in range(1, s + 1):
        value = vector[2 * jets + degree] - vector[jets + degree]
        if degree % s == 0:
            value -= vector[degree // s]
        out.append(value % P)
    return out


def psi_obstruction(vector: list[int], s: int, jets: int) -> list[int]:
    """Coefficients b^0,...,b^(s-1) of c-r-s*b^(s-1)e(b^s)."""
    out: list[int] = []
    for degree in range(s):
        value = vector[2 * jets + degree] - vector[jets + degree]
        shifted = degree - (s - 1)
        if shifted >= 0 and shifted % s == 0:
            value -= s * vector[shifted // s]
        out.append(value % P)
    return out


def branchwise_derivative(vector: list[int], branches: int, jets: int) -> list[int]:
    out = [0] * (branches * jets)
    for branch in range(branches):
        offset = branch * jets
        for degree in range(1, jets):
            out[offset + degree - 1] = degree * vector[offset + degree] % P
    return out


def conductor_functional(vector: list[int], s: int, jets: int) -> list[int]:
    """The explicit c_s/J_f functional in theorem equation (21)."""
    out: list[int] = []
    for degree in range(s + 1, 2 * s + 1):
        value = vector[2 * jets + degree] + vector[jets + degree]
        if degree % s == 0:
            value += s * vector[degree // s]
        out.append(value % P)
    return out


def audit_triple(s: int) -> tuple[int, int, int, int]:
    jets = 2 * s + 5
    monomials = [
        triple_function(s, jets, i, j)
        for i in range(jets)
        for j in range(jets)
    ]
    forms = [
        vector
        for i in range(jets)
        for j in range(jets)
        for vector in (triple_db(s, jets, i, j), triple_dz(s, jets, i, j))
    ]

    function_rank = rank_mod(monomials)
    form_rank = rank_mod(forms)
    check(function_rank == 3 * jets - 2 - s, f"D_{s} function rank")
    check(form_rank == 3 * jets - s, f"D_{s} form rank")

    # Construct the two obstruction maps themselves. Every ambient
    # restriction is killed, both maps are onto, and the kernel dimensions
    # equal the corresponding restriction ranks.
    values = equal_value_basis(jets)
    phi_matrix = [phi_obstruction(vector, s, jets) for vector in values]
    check(rank_mod(phi_matrix) == s, f"D_{s} Phi not onto")
    check(
        all(not any(phi_obstruction(vector, s, jets)) for vector in monomials),
        f"D_{s} ambient function escaped ker Phi",
    )
    check(
        function_rank == len(values) - rank_mod(phi_matrix),
        f"D_{s} function exactness",
    )

    normalized_forms = [
        standard_vector(3, jets, branch, degree)
        for branch in range(3)
        for degree in range(jets)
    ]
    psi_matrix = [psi_obstruction(vector, s, jets) for vector in normalized_forms]
    check(rank_mod(psi_matrix) == s, f"D_{s} Psi not onto")
    check(
        all(not any(psi_obstruction(vector, s, jets)) for vector in forms),
        f"D_{s} ambient form escaped ker Psi",
    )
    check(
        form_rank == len(normalized_forms) - rank_mod(psi_matrix),
        f"D_{s} form exactness",
    )

    # On every basis vector of the equal-value space, Psi(dv) is the
    # coefficientwise derivative of Phi(v).
    for vector in values:
        phi = phi_obstruction(vector, s, jets)
        psi = psi_obstruction(branchwise_derivative(vector, 3, jets), s, jets)
        expected = [(degree + 1) * phi[degree] % P for degree in range(s)]
        check(psi == expected, f"D_{s} differentiation compatibility")

    # The conductor consists of independent branch numerators divisible by
    # z^2 on E and b^(s+1) on R,C.
    conductor = [
        standard_vector(3, jets, 0, degree)
        for degree in range(2, jets)
    ] + [
        standard_vector(3, jets, branch, degree)
        for branch in (1, 2)
        for degree in range(s + 1, jets)
    ]
    conductor_rank = rank_mod(conductor)
    check(conductor_rank == 3 * jets - 2 * s - 4, f"D_{s} conductor rank")

    # f_b restricts as (z^2,0,-s b^(2s)); f_z as
    # (0,-b^(s+1),b^(s+1)).  Multiplication by every ambient monomial spans
    # the normalized image of the Jacobian ideal J_f=(f_b,f_z).
    f_b = [0] * (3 * jets)
    add_entry(f_b, 0, jets, 2, 1)
    add_entry(f_b, 2, jets, 2 * s, -s)
    f_z = [0] * (3 * jets)
    add_entry(f_z, 1, jets, s + 1, -1)
    add_entry(f_z, 2, jets, s + 1, 1)
    jacobian = [
        component_product(generator, monomial, 3, jets)
        for generator in (f_b, f_z)
        for monomial in monomials
    ]
    jacobian_rank = rank_mod(jacobian)
    check(
        conductor_rank - jacobian_rank == s,
        f"D_{s} conductor/Jacobian length",
    )
    check(all(in_span(row, conductor) for row in jacobian), f"D_{s} J not in conductor")

    # bz has normalized numerator (0,0,b^(s+1)). Its first s multiples form
    # a basis modulo J, while b^s*bz is already in J. This certifies the
    # finite truncation of c/J = (R/b^s)[bz], not only its dimension.
    bz = standard_vector(3, jets, 2, s + 1)
    cyclic = [standard_vector(3, jets, 2, s + 1 + degree) for degree in range(s)]
    check(in_span(bz, conductor), f"D_{s} hostile outside conductor")
    check(not in_span(bz, jacobian), f"D_{s} hostile entered Jacobian ideal")
    check(
        rank_mod(jacobian + cyclic) == conductor_rank,
        f"D_{s} cyclic quotient does not span conductor/J",
    )
    check(
        rank_mod(jacobian + cyclic) - jacobian_rank == s,
        f"D_{s} cyclic quotient basis dependent",
    )
    check(
        in_span(standard_vector(3, jets, 2, 2 * s + 1), jacobian),
        f"D_{s} b^s*bz not killed modulo J",
    )
    check(
        all(not any(conductor_functional(vector, s, jets)) for vector in jacobian),
        f"D_{s} quotient functional did not kill J",
    )
    check(
        [conductor_functional(vector, s, jets) for vector in cyclic]
        == [[1 if row == column else 0 for row in range(s)] for column in range(s)],
        f"D_{s} quotient functional not dual to cyclic basis",
    )
    hostile_form = standard_vector(3, jets, 2, 0)
    check(not in_span(hostile_form, forms), f"D_{s} hostile form became Kahler")

    # Check every nonzero truncated numerator arising from a regular
    # normalized branch-form basis vector. Terms shifted beyond the
    # truncation are omitted instead of being counted as vacuous zero checks.
    dualizing_checks = 0
    for branch in range(3):
        shift = 2 if branch == 0 else s + 1
        sign = -1 if branch in (0, 1) else 1
        for degree in range(jets - shift):
            numerator = [0] * (3 * jets)
            add_entry(numerator, branch, jets, degree + shift, sign)
            check(in_span(numerator, monomials), f"D_{s} dualizing lift failed")
            dualizing_checks += 1
    check(dualizing_checks == 3 * jets - 2 * s - 4, f"D_{s} dualizing check count")

    return jets, function_rank, form_rank, conductor_rank - jacobian_rank


def pair_function(m: int, jets: int, b_degree: int, q_degree: int) -> list[int]:
    out = [0] * (2 * jets)
    if q_degree == 0:
        add_entry(out, 0, jets, b_degree, 1)
    add_entry(out, 1, jets, b_degree + m * q_degree, 1)
    return out


def pair_db(m: int, jets: int, b_degree: int, q_degree: int) -> list[int]:
    out = [0] * (2 * jets)
    if q_degree == 0:
        add_entry(out, 0, jets, b_degree, 1)
    add_entry(out, 1, jets, b_degree + m * q_degree, 1)
    return out


def pair_dq(m: int, jets: int, b_degree: int, q_degree: int) -> list[int]:
    out = [0] * (2 * jets)
    add_entry(out, 1, jets, b_degree + m * q_degree + m - 1, m)
    return out


def pair_conductor_functional(vector: list[int], m: int, jets: int) -> list[int]:
    """Endpoint analogue [h_C+h_R]/b^m modulo b^(m-1)."""
    return [
        (vector[degree] + vector[jets + degree]) % P
        for degree in range(m, 2 * m - 1)
    ]


def audit_pair_endpoint(m: int = 12) -> tuple[int, int, int]:
    jets = 2 * m + 5
    monomials = [
        pair_function(m, jets, i, j)
        for i in range(jets)
        for j in range(jets)
    ]
    forms = [
        vector
        for i in range(jets)
        for j in range(jets)
        for vector in (pair_db(m, jets, i, j), pair_dq(m, jets, i, j))
    ]
    function_rank = rank_mod(monomials)
    form_rank = rank_mod(forms)
    check(function_rank == 2 * jets - m, "A_m function rank")
    check(form_rank == 2 * jets - (m - 1), "A_m form rank")

    conductor = [
        standard_vector(2, jets, branch, degree)
        for branch in (0, 1)
        for degree in range(m, jets)
    ]
    f_q = [0] * (2 * jets)
    add_entry(f_q, 0, jets, m, -1)
    add_entry(f_q, 1, jets, m, 1)
    f_b = [0] * (2 * jets)
    add_entry(f_b, 1, jets, 2 * m - 1, -m)
    jacobian = [
        component_product(generator, monomial, 2, jets)
        for generator in (f_b, f_q)
        for monomial in monomials
    ]
    quotient_length = rank_mod(conductor) - rank_mod(jacobian)
    check(quotient_length == m - 1, "A_m conductor/Jacobian length")
    hostile = standard_vector(2, jets, 1, m)
    cyclic = [standard_vector(2, jets, 1, m + degree) for degree in range(m - 1)]
    check(in_span(hostile, conductor), "A_m hostile outside conductor")
    check(not in_span(hostile, jacobian), "A_m hostile entered Jacobian ideal")
    check(
        rank_mod(jacobian + cyclic) == rank_mod(conductor),
        "A_m cyclic quotient does not span conductor/J",
    )
    check(
        rank_mod(jacobian + cyclic) - rank_mod(jacobian) == m - 1,
        "A_m cyclic quotient basis dependent",
    )
    check(
        in_span(standard_vector(2, jets, 1, 2 * m - 1), jacobian),
        "A_m b^(m-1)*q not killed modulo J",
    )
    check(
        all(not any(pair_conductor_functional(vector, m, jets)) for vector in jacobian),
        "A_m quotient functional did not kill J",
    )
    check(
        [pair_conductor_functional(vector, m, jets) for vector in cyclic]
        == [
            [1 if row == column else 0 for row in range(m - 1)]
            for column in range(m - 1)
        ],
        "A_m quotient functional not dual to cyclic basis",
    )
    hostile_form = standard_vector(2, jets, 1, 0)
    check(not in_span(hostile_form, forms), "A_m hostile form became Kahler")
    return function_rank, form_rank, quotient_length


def dot(row: list[int], vector: list[int], prime: int) -> int:
    return sum(left * right for left, right in zip(row, vector)) % prime


def first_observed_order(
    matrix: list[list[int]], vector: list[int], labels: tuple[int, ...], prime: int
) -> int:
    hits = [label for label, row in zip(labels, matrix) if dot(row, vector, prime)]
    check(bool(hits), "nonzero observer kernel vector vanished in full observer")
    return min(hits)


def audit_observers() -> tuple[int, int]:
    # A split finite-field control for the exact matrix in THM-4280.  The
    # identities are algebraic; this checks ranks and the two stated kernels.
    prime = 97
    kappa = next(
        value
        for value in range(2, prime)
        if pow(value, 12, prime) == 1
        and all(pow(value, divisor, prime) != 1 for divisor in (1, 2, 3, 4, 6))
    )
    omega = (-kappa * kappa) % prime
    inv_two = pow(2, prime - 2, prime)
    matrix = [
        [0, 1, -kappa, (omega * omega - kappa) * inv_two],
        [1, 0, 0, 0],
        [0, 0, 0, inv_two],
        [0, 1, kappa, (omega * omega + kappa) * inv_two],
    ]
    matrix = [[entry % prime for entry in row] for row in matrix]
    check(rank_mod(matrix, prime) == 4, "four-channel rank")

    kernel_124 = [0, kappa, 1, 0]
    kernel_247 = [0, kappa, -1, 0]
    check(all(dot(matrix[index], kernel_124, prime) == 0 for index in (0, 1, 2)), "q124 kernel")
    check(dot(matrix[3], kernel_124, prime) != 0, "q124 missing channel")
    check(all(dot(matrix[index], kernel_247, prime) == 0 for index in (1, 2, 3)), "q247 kernel")
    check(dot(matrix[0], kernel_247, prime) != 0, "q247 missing channel")
    labels = (1, 2, 4, 7)
    return (
        first_observed_order(matrix, kernel_124, labels, prime),
        first_observed_order(matrix, kernel_247, labels, prime),
    )


def main() -> None:
    print("THM4289_A23_BLOWDOWN_KAHLER_DUALIZING_OBSERVER_AUDIT_V1")
    print("FIELD F_1000000007 TRIPLE_TRUNCATION N=2S+5")
    for s in range(1, 12):
        jets, function_rank, form_rank, quotient_length = audit_triple(s)
        print(
            f"TRIPLE S {s} N {jets} FUNCTION_RANK {function_rank} "
            f"EXPECTED {3 * jets - 2 - s} FORM_RANK {form_rank} "
            f"EXPECTED {3 * jets - s} C_OVER_J_LENGTH {quotient_length}"
        )

    function_rank, form_rank, quotient_length = audit_pair_endpoint()
    print(
        f"PAIR M 12 N 29 FUNCTION_RANK {function_rank} FORM_RANK {form_rank} "
        f"C_OVER_J_LENGTH {quotient_length}"
    )
    order_124, order_247 = audit_observers()
    check((order_124, order_247) == (7, 1), "observer order regression")
    print(
        f"OBSERVER q124 KERNEL_ORDER {order_124} PRESERVES_P_LE_{order_124} "
        f"FIRST_FAIL_P_{order_124 + 1} CURVE_FIRST_FAIL_D_{order_124}"
    )
    print(
        f"OBSERVER q247 KERNEL_ORDER {order_247} PRESERVES_P_LE_{order_247} "
        f"FIRST_FAIL_P_{order_247 + 1} CURVE_FIRST_FAIL_D_{order_247}"
    )
    print("DEPENDENCY_STATE_SET ARBITRARY 1,2,3,4,5,6,7,8,9,10,11,12")
    print("DEPENDENCY_STATE_SET COMPLEX_HOM 1,2,4,7,12")
    print("DEPENDENCY_STATE_SET ACTUAL_O_HOM 1,2,4,12")
    print("DEPENDENCY_STATE_SET DEGREE_34_42 1")
    for ell in (1, 2, 4):
        channels = 12 - ell
        differential_length = min(channels, 2 * ell)
        print(
            f"DEPENDENCY_FORMULA PARTIAL ELL {ell} REQUIRED_CHANNELS {channels} "
            f"OMEGA_LENGTH {differential_length}"
        )
    print("HOSTILES bz*eta=(0,0,db) AND q*eta=(0,db) PASS")
    print("VERDICT PASS FINITE_EXACT_CONTROLS; JC2 OPEN")


if __name__ == "__main__":
    main()
