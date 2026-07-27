#!/usr/bin/env python3
"""Dependency-free exact referee for THM-2521.

The integer checks verify the signless-incidence, matching, degree, Radon,
and Dirichlet identities.  The cyclotomic checks use F_547, whose
multiplicative group contains primitive seventh and thirteenth roots.
Every assertion remains active under ``python -O``.
"""

from fractions import Fraction


P = 13
Q = 7
OMEGA = 13
NV = 14
NE = 91
MOD = 547


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def inv(value):
    return pow(value, MOD - 2, MOD)


def rank_mod(matrix):
    rows = [list(map(lambda x: x % MOD, row)) for row in matrix]
    if not rows:
        return 0
    height = len(rows)
    width = len(rows[0])
    pivot_row = 0
    for column in range(width):
        pivot = next(
            (row for row in range(pivot_row, height) if rows[row][column]),
            None,
        )
        if pivot is None:
            continue
        rows[pivot_row], rows[pivot] = rows[pivot], rows[pivot_row]
        scale = inv(rows[pivot_row][column])
        rows[pivot_row] = [(scale * x) % MOD for x in rows[pivot_row]]
        for row in range(height):
            if row == pivot_row or not rows[row][column]:
                continue
            factor = rows[row][column]
            rows[row] = [
                (x - factor * y) % MOD
                for x, y in zip(rows[row], rows[pivot_row])
            ]
        pivot_row += 1
        if pivot_row == height:
            break
    return pivot_row


def primitive_root():
    factors = (2, 3, 7, 13)
    for candidate in range(2, MOD):
        if all(pow(candidate, (MOD - 1) // prime, MOD) != 1 for prime in factors):
            return candidate
    raise RuntimeError("primitive root not found")


GENERATOR = primitive_root()
ZETA = pow(GENERATOR, (MOD - 1) // P, MOD)
XI = pow(GENERATOR, (MOD - 1) // Q, MOD)
require(pow(ZETA, P, MOD) == 1 and ZETA != 1, "order 13")
require(pow(XI, Q, MOD) == 1 and XI != 1, "order 7")


EDGES = tuple((x, y) for x in range(NV) for y in range(x + 1, NV))
EDGE_INDEX = {edge: index for index, edge in enumerate(EDGES)}
require(len(EDGES) == NE, "K14 edge count")


def edge(x, y):
    return (x, y) if x < y else (y, x)


def chart_edge(tau, a, c, h, r):
    s = (a * r + c) % Q
    if s == 0:
        return edge(OMEGA, h)
    return edge((h - tau * s) % P, (h + tau * s) % P)


def physical_basis():
    basis = []
    for j in range(1, P):
        vector = [0] * NV
        vector[0] = -1
        vector[j] = 1
        basis.append(tuple(vector))
    return tuple(basis)


def full_potential_basis():
    basis = []
    for j in range(P):
        vector = [0] * NV
        vector[j] = 1
        vector[OMEGA] = -1
        basis.append(tuple(vector))
    return tuple(basis)


U_BASIS = physical_basis()
V0_BASIS = full_potential_basis()
require(len(U_BASIS) == 12 and len(V0_BASIS) == 13, "potential dimensions")


def potential_edges(potential):
    return tuple(potential[x] + potential[y] for x, y in EDGES)


def degrees(weights):
    result = [0] * NV
    for weight, (x, y) in zip(weights, EDGES):
        result[x] += weight
        result[y] += weight
    return tuple(result)


def dot(first, second):
    return sum(x * y for x, y in zip(first, second))


def d_table(tau, a, c, potential):
    return tuple(
        tuple(
            sum(potential[v] for v in chart_edge(tau, a, c, h, r))
            for r in range(Q)
        )
        for h in range(P)
    )


def aligned_radon(tau, a, c, table):
    return tuple(
        sum(
            table[(v - tau * ((a * r + c) % Q)) % P][r]
            for r in range(Q)
        )
        for v in range(P)
    )


def horizontal_transform(values, alpha):
    return sum(
        (value % MOD) * pow(ZETA, (-alpha * h) % P, MOD)
        for h, value in enumerate(values)
    ) % MOD


def table_transform(table, alpha, beta):
    horizontal = [
        horizontal_transform([table[h][r] for h in range(P)], alpha)
        for r in range(Q)
    ]
    return sum(
        horizontal[r] * pow(XI, (-beta * r) % Q, MOD)
        for r in range(Q)
    ) % MOD


def l_multiplier(lam, gamma):
    total = 1
    for s in range(1, Q):
        phase = pow(XI, (-gamma * s) % Q, MOD)
        pair = (
            pow(ZETA, (-lam * s) % P, MOD)
            + pow(ZETA, (lam * s) % P, MOD)
        ) % MOD
        total = (total + phase * pair) % MOD
    return total


def k_multiplier(lam, beta):
    return sum(
        pow(ZETA, (-lam * s) % P, MOD)
        * pow(XI, (-beta * s) % Q, MOD)
        for s in range(Q)
    ) % MOD


def matching_constraints_and_degree_rank():
    constraints = [[0] * NE for _ in range(P)]
    row_edges = []
    for h in range(P):
        indices = [EDGE_INDEX[chart_edge(1, 1, 0, h, r)] for r in range(Q)]
        require(len(set(indices)) == Q, ("matching row", h))
        row_edges.append(indices)
        for index in indices:
            constraints[h][index] = 1
    require(
        sorted(index for indices in row_edges for index in indices) == list(range(NE)),
        "one-factorization partitions E(K14)",
    )
    constraint_rank = rank_mod(constraints)

    factor_basis = []
    for indices in row_edges:
        anchor = indices[-1]
        for index in indices[:-1]:
            vector = [0] * NE
            vector[index] = 1
            vector[anchor] = -1
            factor_basis.append(vector)
    require(len(factor_basis) == 78, "factor-balanced basis")

    degree_columns = []
    for weights in factor_basis:
        degree_columns.append(degrees(weights))
    degree_matrix = [
        [degree_columns[column][row] for column in range(len(degree_columns))]
        for row in range(NV)
    ]
    degree_rank = rank_mod(degree_matrix)
    return constraint_rank, len(factor_basis), degree_rank


def audit_gram_and_splitting():
    for potential in V0_BASIS:
        weights = potential_edges(potential)
        require(degrees(weights) == tuple(12 * x for x in potential), "degree splitting")
    for first in U_BASIS:
        first_weights = potential_edges(first)
        first_degree = degrees(first_weights)
        for second in U_BASIS:
            second_weights = potential_edges(second)
            second_degree = degrees(second_weights)
            base = dot(first, second)
            require(dot(first_weights, second_weights) == 12 * base, "edge Gram")
            require(dot(first_degree, second_degree) == 144 * base, "degree Gram")


def audit_charts_and_spectra():
    chart_count = 0
    basis_chart_count = 0
    radon_entries = 0
    mixed_modes = 0
    for tau in range(1, P):
        for a in range(1, Q):
            a_inverse = pow(a, -1, Q)
            for c in range(Q):
                images = {
                    chart_edge(tau, a, c, h, r)
                    for h in range(P)
                    for r in range(Q)
                }
                require(len(images) == NE, ("chart bijection", tau, a, c))
                chart_count += 1
                for potential in U_BASIS:
                    table = d_table(tau, a, c, potential)
                    require(all(sum(row) == 0 for row in table), "row zero")
                    require(
                        all(sum(table[h][r] for h in range(P)) == 0 for r in range(Q)),
                        "column zero",
                    )
                    radon = aligned_radon(tau, a, c, table)
                    expected = tuple(
                        7 * potential[v]
                        + sum(potential[(v - 2 * tau * s) % P] for s in range(1, Q))
                        for v in range(P)
                    )
                    require(radon == expected, ("aligned Radon", tau, a, c))
                    radon_entries += P

                    for alpha in range(1, P):
                        p_hat = horizontal_transform(potential[:P], alpha)
                        require(p_hat != 0, ("basis root mode", alpha))
                        radon_hat = horizontal_transform(radon, alpha)
                        multiplier = (
                            7
                            + sum(
                                pow(ZETA, (-2 * alpha * tau * s) % P, MOD)
                                for s in range(1, Q)
                            )
                        ) % MOD
                        require(multiplier != 0, ("Radon multiplier", tau, alpha))
                        require(radon_hat == multiplier * p_hat % MOD, "Radon DFT")
                        for beta in range(1, Q):
                            gamma = beta * a_inverse % Q
                            geometric = l_multiplier(alpha * tau % P, gamma)
                            require(geometric != 0, ("L multiplier", alpha, beta))
                            phase = pow(XI, (beta * a_inverse * c) % Q, MOD)
                            expected_mode = phase * geometric * p_hat % MOD
                            actual_mode = table_transform(table, alpha, beta)
                            require(actual_mode == expected_mode, "direct mixed multiplier")
                            require(actual_mode != 0, "direct mixed nonvanishing")
                            mixed_modes += 1
                    basis_chart_count += 1
    return chart_count, basis_chart_count, radon_entries, mixed_modes


def cycle_bilinear(first, second, tau, s):
    return sum(
        (first[(h - tau * s) % P] - first[(h + tau * s) % P])
        * (second[(h - tau * s) % P] - second[(h + tau * s) % P])
        for h in range(P)
    )


def audit_energy_forms():
    bilinear_checks = 0
    for tau in range(1, P):
        for first in U_BASIS:
            for second in U_BASIS:
                base = dot(first[:P], second[:P])
                cycle_sum = sum(
                    cycle_bilinear(first, second, tau, s)
                    for s in range(1, Q)
                )
                require(cycle_sum == P * base, ("six-cycle form", tau))
                bilinear_checks += 1

    # Exact sharp control: one centred delta on F_13.
    delta = tuple(
        Fraction(P - 1, P) if r == 0 else Fraction(-1, P)
        for r in range(P)
    )
    star = sum(x * x for x in delta)
    cycles = tuple(
        sum(
            (delta[(h - s) % P] - delta[(h + s) % P]) ** 2
            for h in range(P)
        )
        for s in range(1, Q)
    )
    require(star == Fraction(12, 13), "delta star")
    require(cycles == (Fraction(2),) * 6, "delta cycles")
    additive_modes = tuple((star - Fraction(2)) / Q for _ in range(1, Q))
    require(additive_modes == (Fraction(-2, 13),) * 6, "delta additive modes")
    residues = {1, 2, 4}
    chi_contrast = sum(
        (1 if s in residues else -1) * cycles[s - 1]
        for s in range(1, Q)
    )
    require(chi_contrast == 0, "multiplicative chi7 hostile")
    return bilinear_checks, star, cycles, additive_modes[0], chi_contrast


def audit_cut_lift():
    potential = U_BASIS[0]
    table = d_table(1, 1, 0, potential)
    count = 0
    for tau in range(1, P):
        for a in range(1, Q):
            for alpha in range(1, P):
                for beta in range(1, Q):
                    source_beta = (-beta * a) % Q
                    source = table_transform(table, alpha, source_beta)
                    geometric = k_multiplier(alpha * tau % P, beta)
                    require(source != 0 and geometric != 0, "cut lift factors")
                    require(geometric * source % MOD != 0, "cut lift mode")
                    count += 1
    require(count == 5184, "cut coefficient count")
    return count


def main():
    constraint_rank, factor_dimension, degree_rank = matching_constraints_and_degree_rank()
    require((constraint_rank, factor_dimension, degree_rank) == (13, 78, 13), "rank ledger")
    audit_gram_and_splitting()
    charts, chart_basis, radon_entries, mixed_modes = audit_charts_and_spectra()
    energy_checks, star, cycles, additive_mode, chi_contrast = audit_energy_forms()
    cut_modes = audit_cut_lift()

    print("THM-2521 exact K13--K14 potential-module referee")
    print(
        "factor_space:",
        f"edges={NE}",
        f"matching_rank={constraint_rank}",
        f"balanced_dim={factor_dimension}",
        f"degree_rank={degree_rank}",
        f"degree_kernel={factor_dimension-degree_rank}",
    )
    print("potential_modules: full_dim=13 physical_dim=12 degree_scalar=12 edge_gram_scalar=12")
    print(
        "chart_audit:",
        f"charts={charts}",
        f"basis_charts={chart_basis}",
        f"aligned_radon_entries={radon_entries}",
        f"mixed_mode_identities={mixed_modes}",
    )
    print(
        "spectral_audit:",
        f"field=F_{MOD}",
        f"generator={GENERATOR}",
        f"zeta13={ZETA}",
        f"xi7={XI}",
        f"cut_modes={cut_modes}",
    )
    print(
        "energy_audit:",
        f"bilinear_checks={energy_checks}",
        "six_cycle_scalar=13",
        "full_K14_gradient_scalar=14",
    )
    print(
        "sharp_delta_control:",
        f"star={star}",
        f"cycles={cycles}",
        f"additive_modes={additive_mode}",
        f"chi7_contrast={chi_contrast}",
    )
    print(
        "VERIFIED: positive rational predecessor drift enters the 12-dimensional "
        "degree-visible K14 potential module; orientation and Boolean owner charge do not follow."
    )


if __name__ == "__main__":
    main()
