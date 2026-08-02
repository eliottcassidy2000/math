#!/usr/bin/env python3
"""Exact order-five transport-circuit boundary behind THM-3110.

The 28-mass forward row-majorization transport in THM-3110 leaves signed
cores of mass 9 and 11.  Each core admits a complete *reverse* transport.
Combining the forward and reverse edge currents gives a unique minimal
Schur circuit through degree four.  This script verifies the four chamber
universes, the exact coefficient-flag ranks, the degree-five exit, and the
smallest raw TP2/ECT hostile.

All arithmetic is integral or rational.  SymPy is used only for exact rank,
small polynomial identities, and no numerical approximation.
"""

import importlib.util
from collections import Counter, defaultdict
from functools import lru_cache
from math import factorial
from pathlib import Path

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


HERE = Path(__file__).resolve().parent
UPSTREAM = HERE / "gmc_product_gamma_arbitrary_anchored_schur_thm3110.py"
SPEC = importlib.util.spec_from_file_location("thm3110_schur", UPSTREAM)
THM3110 = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(THM3110)
BANKS = THM3110.BANKS


# The chamber-independent signed cores left by a compatible 28-mass
# forward transport.  Linear forms are encoded as alpha*a+beta*b.
RESIDUAL_BANKS = (
    (
        (2, ((2, 0), (2, 0), (1, 0), (0, 1), (0, 1), (0, 1))),
        (2, ((2, 1), (1, 0), (1, 0), (1, 0), (0, 1), (0, 1))),
        (5, ((3, 0), (1, 1), (1, 0), (0, 1), (0, 1))),
        (-1, ((2, 0), (1, 0), (1, 0), (1, 0), (0, 3))),
        (-3, ((2, 0), (2, 0), (1, 2), (0, 1))),
        (-3, ((2, 1), (2, 0), (1, 0), (0, 2))),
        (-2, ((3, 0), (1, 1), (1, 1), (0, 1))),
    ),
    (
        (3, ((2, 0), (1, 1), (1, 0), (1, 0), (0, 1), (0, 1), (0, 1))),
        (3, ((3, 0), (1, 0), (1, 0), (0, 2), (0, 1), (0, 1))),
        (4, ((3, 0), (1, 1), (1, 1), (0, 1), (0, 1))),
        (1, ((3, 0), (2, 0), (0, 1), (0, 1), (0, 1), (0, 1))),
        (-3, ((2, 0), (1, 2), (1, 0), (1, 0), (0, 2))),
        (-2, ((2, 0), (2, 0), (1, 0), (0, 3), (0, 1))),
        (-1, ((2, 1), (2, 0), (1, 1), (0, 1), (0, 1))),
        (-4, ((3, 0), (1, 1), (1, 0), (0, 2), (0, 1))),
        (-1, ((3, 0), (2, 0), (0, 2), (0, 1), (0, 1))),
    ),
)


def signed_parts(bank):
    positive = tuple((coefficient, row) for coefficient, row in bank
                     if coefficient > 0)
    negative = tuple((-coefficient, row) for coefficient, row in bank
                     if coefficient < 0)
    return positive, negative


def transport_edges(left, right, supply, demand, chamber):
    """Return one deterministic full row-majorization transport."""

    mass = sum(supply)
    require(mass == sum(demand), "transport mass imbalance")
    size = 2 + len(left) + len(right)
    source = size - 2
    target = size - 1
    flow = THM3110.Dinic(size)
    for index, capacity in enumerate(supply):
        flow.add_edge(source, index, capacity)
    for index, capacity in enumerate(demand):
        flow.add_edge(len(left) + index, target, capacity)

    references = []
    for left_index, (_, high) in enumerate(left):
        for right_index, (_, low) in enumerate(right):
            if THM3110.rows_majorize_in_chamber(high, low, chamber):
                edge_index = len(flow.graph[left_index])
                flow.add_edge(left_index, len(left) + right_index, mass)
                references.append((left_index, right_index, edge_index))

    require(flow.maximum_flow(source, target) == mass,
            "claimed transport is not complete")
    answer = []
    for left_index, right_index, edge_index in references:
        used = mass - flow.graph[left_index][edge_index][1]
        if used:
            answer.append((left_index, right_index, used))
    return tuple(answer)


def row_counter(bank):
    answer = defaultdict(int)
    for coefficient, row in bank:
        answer[row] += coefficient
    return Counter({row: coefficient for row, coefficient in answer.items()
                    if coefficient})


def circuit_universe(invariant, chamber):
    positive, negative = signed_parts(BANKS[invariant])
    residual_positive, residual_negative = signed_parts(
        RESIDUAL_BANKS[invariant])

    positive_residue = row_counter(residual_positive)
    negative_residue = row_counter(residual_negative)
    forward_supply = tuple(
        coefficient - positive_residue[row] for coefficient, row in positive
    )
    forward_demand = tuple(
        coefficient - negative_residue[row] for coefficient, row in negative
    )
    require(sum(forward_supply) == sum(forward_demand) == 28,
            "forward mass drift")
    forward = transport_edges(
        positive, negative, forward_supply, forward_demand, chamber)

    reverse = transport_edges(
        residual_negative,
        residual_positive,
        tuple(coefficient for coefficient, _ in residual_negative),
        tuple(coefficient for coefficient, _ in residual_positive),
        chamber,
    )
    residual_mass = (9, 11)[invariant]
    require(sum(edge[2] for edge in reverse) == residual_mass,
            "reverse mass drift")

    edges = []
    for positive_index, negative_index, capacity in forward:
        edges.append((
            "F", capacity,
            positive[positive_index][1], negative[negative_index][1],
        ))
    for negative_index, positive_index, capacity in reverse:
        edges.append((
            "R", capacity,
            residual_negative[negative_index][1],
            residual_positive[positive_index][1],
        ))

    # Expanding the signed edge current recovers the original signed bank:
    # forward differences enter positively, reverse differences negatively.
    expanded = defaultdict(int)
    for kind, capacity, high, low in edges:
        sign = 1 if kind == "F" else -1
        expanded[high] += sign * capacity
        expanded[low] -= sign * capacity
    require(
        Counter({row: coefficient for row, coefficient in expanded.items()
                 if coefficient}) == row_counter(BANKS[invariant]),
        "edge circuit does not recover the original bank",
    )
    positive_lookup = {row: index for index, (_, row) in enumerate(positive)}
    negative_lookup = {row: index for index, (_, row) in enumerate(negative)}
    reverse_in_original_indices = tuple(
        (
            negative_lookup[residual_negative[negative_index][1]],
            positive_lookup[residual_positive[positive_index][1]],
            capacity,
        )
        for negative_index, positive_index, capacity in reverse
    )
    return tuple(edges), forward, reverse_in_original_indices


def roots(row, first, second):
    return tuple(
        root
        for alpha, beta in row
        for root in range(alpha * first + beta * second)
    )


def forward_heads(values):
    answer = []
    values = list(values)
    while values:
        answer.append(values[0])
        values = [values[index + 1] - values[index]
                  for index in range(len(values) - 1)]
    return tuple(answer)


def rectangular_newton(values, degree):
    """Exact coefficients in binom(a,i) binom(b,j)."""

    first_heads = {
        (first_order, second): coefficient
        for second in range(degree + 1)
        for first_order, coefficient in enumerate(forward_heads(
            values[first, second] for first in range(degree + 1)
        ))
    }
    return {
        (first_order, second_order): coefficient
        for first_order in range(degree + 1)
        for second_order, coefficient in enumerate(forward_heads(
            first_heads[first_order, second]
            for second in range(degree + 1)
        ))
        if coefficient
    }


@lru_cache(maxsize=None)
def symmetric_data(row, first, second):
    return THM3110.complete_and_elementary(roots(row, first, second), 5)


def coefficient_flag_matrix(edges):
    """Flag the edge currents by every Schur coefficient through degree 4."""

    rows = []
    for degree in range(1, 5):
        polynomial_degree = 2 * degree
        for partition in THM3110.partitions(degree):
            edge_polynomials = []
            for _, _, high, low in edges:
                values = {}
                for first in range(polynomial_degree + 1):
                    for second in range(polynomial_degree + 1):
                        values[first, second] = (
                            THM3110.schur_value(
                                symmetric_data(high, first, second), partition)
                            - THM3110.schur_value(
                                symmetric_data(low, first, second), partition)
                        )
                edge_polynomials.append(
                    rectangular_newton(values, polynomial_degree))

            flags = sorted(
                set().union(*(polynomial.keys()
                              for polynomial in edge_polynomials)),
                key=lambda exponent: (sum(exponent), exponent[0]),
                reverse=True,
            )
            for exponent in flags:
                vector = tuple(
                    polynomial.get(exponent, 0)
                    for polynomial in edge_polynomials
                )
                if any(vector):
                    rows.append((degree, partition, exponent, vector))
    return tuple(rows)


def hook_dimension(partition):
    denominator = 1
    for row, length in enumerate(partition):
        for column in range(length):
            denominator *= (
                length - column
                + sum(other > column for other in partition[row + 1:])
            )
    return factorial(sum(partition)) // denominator


def kappa(partition):
    return sum(
        part * (part - 2 * row + 1)
        for row, part in enumerate(partition, 1)
    )


def degree_five_exit_audit():
    points = {
        "tight": ((2, 3), (3, 4), (3, 5)),
        "wide": ((1, 2), (1, 3), (2, 4)),
    }
    for invariant in (0, 1):
        for chamber, supports in points.items():
            for first, second in supports:
                data = THM3110.atom_symmetric_data(
                    invariant, first, second, 5)
                base = THM3110.forced_divisor(invariant, first, second)
                linear = 3 * first + (5 if invariant == 0 else 4) * second
                for partition in THM3110.partitions(5):
                    observed = THM3110.phi_from_data(data, partition)
                    expected = (
                        base * hook_dimension(partition)
                        * (4 * linear + kappa(partition)) // 8
                    )
                    require(observed == expected,
                            "degree-five circuit exit drift")


def p1_polynomial(row, a, b):
    return sp.expand(sum(
        length * (length - 1) / 2
        for length in (alpha * a + beta * b for alpha, beta in row)
    ))


def tp2_hostile():
    """Return the canonical raw 2x2 coefficient-flag sign conflict."""

    positive, negative = signed_parts(BANKS[0])
    residual_positive, residual_negative = signed_parts(RESIDUAL_BANKS[0])
    # At (a,b)=(2,3), these are canonical columns 0, 2, and 4 when
    # ordered by high Ferrers partition and then low Ferrers partition.
    columns = (
        (positive[5][1], negative[2][1]),
        (positive[0][1], negative[0][1]),
        (residual_negative[1][1], residual_positive[0][1]),
    )
    a, b = sp.symbols("a b")
    polynomials = tuple(
        sp.factor(p1_polynomial(high, a, b) - p1_polynomial(low, a, b))
        for high, low in columns
    )
    expected = (a**2, -a * b + 2 * b**2, 2 * a * b + b**2)
    require(all(sp.expand(observed - target) == 0
                for observed, target in zip(polynomials, expected)),
            "TP2 column polynomial drift")

    # Rows are the raw monomial coefficient flags [a^2] and [ab].
    flagged = tuple(
        (sp.Poly(polynomial, a, b).coeff_monomial(a**2),
         sp.Poly(polynomial, a, b).coeff_monomial(a * b))
        for polynomial in polynomials
    )
    negative_minor = flagged[0][0] * flagged[1][1] - flagged[1][0] * flagged[0][1]
    positive_minor = flagged[0][0] * flagged[2][1] - flagged[2][0] * flagged[0][1]
    require((negative_minor, positive_minor) == (-1, 2),
            "TP2 hostile drift")
    return negative_minor, positive_minor


def main():
    expected = {
        (0, "tight"): (330, 21, 20, 16, 5),
        (0, "wide"): (330, 20, 19, 15, 5),
        (1, "tight"): (328, 17, 16, 12, 5),
        (1, "wide"): (328, 20, 19, 14, 6),
    }
    census = []
    orders = []
    for invariant in (0, 1):
        for chamber in ("tight", "wide"):
            edges, forward_order, reverse_order = circuit_universe(
                invariant, chamber)
            flags = coefficient_flag_matrix(edges)
            matrix = sp.Matrix([flag[3] for flag in flags])
            circuit = sp.Matrix([
                capacity if kind == "F" else -capacity
                for kind, capacity, _, _ in edges
            ])
            require(matrix * circuit == sp.zeros(matrix.rows, 1),
                    "low-degree circuit annihilation drift")
            rank = matrix.rank()
            forward_edges = sum(kind == "F" for kind, _, _, _ in edges)
            reverse_edges = sum(kind == "R" for kind, _, _, _ in edges)
            observed = (
                len(flags), len(edges), rank, forward_edges, reverse_edges)
            require(observed == expected[invariant, chamber],
                    "edge/rank census drift")
            require(rank == len(edges) - 1,
                    "circuit is not the unique truncated kernel")
            require(all(entry != 0 for entry in circuit),
                    "circuit lost full support")
            census.append((invariant, chamber, *observed))
            orders.append((invariant, chamber, forward_order, reverse_order))

    degree_five_exit_audit()
    negative_minor, positive_minor = tp2_hostile()

    print("transport=forward_mass:28;reverse_core:I1:9/9;I2:11/11")
    print("circuit_census=" + ";".join(
        f"I{invariant + 1}-{chamber}:matrix{rows}x{edges}:rank{rank}:nullity1:"
        f"F{forward}:R{reverse}"
        for invariant, chamber, rows, edges, rank, forward, reverse in census
    ))
    for invariant, chamber, forward, reverse in orders:
        forward_text = ",".join(
            f"P{left}>N{right}x{capacity}"
            for left, right, capacity in forward)
        reverse_text = ",".join(
            f"N{left}>P{right}x{capacity}"
            for left, right, capacity in reverse)
        print(f"edge_order_I{invariant + 1}_{chamber}=F:{forward_text}|R:{reverse_text}")
    print("flag_basis=binom(a,i)binom(b,j);Schur_degrees=1..4")
    print("low_degree=annihilation_exact;kernel_dimension=1;full_support=PASS")
    print("degree5_exit=Phi_Ij(s_lambda)=base_j*f_lambda/2*(L_j+kappa_lambda/4)")
    print("tp2_columns=a^2,2b^2-ab,b^2+2ab")
    print(f"tp2_raw_minors={negative_minor},{positive_minor}:SIGN_CONFLICT")
    print("boundary=no_ordinary_TP_positroid_ECT_or_coefficientwise_Peano_kernel")
    print("needed_sidecar=chamber_adapted_Newton_cone_or_rooted_rank4_flag")
    print("all_exact_checks=PASS")


if __name__ == "__main__":
    main()
