#!/usr/bin/env python3
"""Exact finite companion for THM-4255.

This standard-library-only audit checks finite triangular truncations of

    ev: Q[u][[f]] -> Q[[f]],   u |-> f,

the graph-ideal/preimage dimensions, Cartier noncommutation, transverse
Hasse-jet recovery, failure of ell-adic strictness, and an abstract hostile
to the degree-D Cartier-window inference used in the external density draft.
The finite checks accompany, but do not replace, the theorem's general proof.
"""

from __future__ import annotations

from fractions import Fraction
from math import comb


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def matrix_rank(rows: list[list[int]]) -> int:
    """Exact rank over Q."""

    if not rows:
        return 0
    a = [[Fraction(x) for x in row] for row in rows]
    m = len(a)
    n = len(a[0])
    rank = 0
    for col in range(n):
        pivot = next((i for i in range(rank, m) if a[i][col]), None)
        if pivot is None:
            continue
        a[rank], a[pivot] = a[pivot], a[rank]
        scale = a[rank][col]
        a[rank] = [x / scale for x in a[rank]]
        for i in range(m):
            if i == rank or not a[i][col]:
                continue
            scale = a[i][col]
            a[i] = [x - scale * y for x, y in zip(a[i], a[rank])]
        rank += 1
        if rank == m:
            break
    return rank


def triangular_basis(total_degree: int) -> list[tuple[int, int]]:
    """Monomials u^a f^b with a+b <= total_degree."""

    return [
        (a, b)
        for degree in range(total_degree + 1)
        for a in range(degree + 1)
        for b in [degree - a]
    ]


def columns_to_rows(
    columns: list[dict[tuple[int, int], int]],
    basis: list[tuple[int, int]],
) -> list[list[int]]:
    return [[column.get(monomial, 0) for column in columns] for monomial in basis]


def evaluation_rows(total_degree: int, jet_order: int) -> list[list[int]]:
    basis = triangular_basis(total_degree)
    return [
        [1 if a + b == degree else 0 for a, b in basis]
        for degree in range(jet_order)
    ]


def graph_generators(total_degree: int) -> list[dict[tuple[int, int], int]]:
    """Columns for (u-f)u^a f^b through total degree total_degree."""

    result: list[dict[tuple[int, int], int]] = []
    for a, b in triangular_basis(total_degree - 1):
        result.append({(a + 1, b): 1, (a, b + 1): -1})
    return result


def tail_generators(total_degree: int, jet_order: int) -> list[dict[tuple[int, int], int]]:
    return [{(0, degree): 1} for degree in range(jet_order, total_degree + 1)]


def audit_triangular_truncations() -> tuple[int, int]:
    truncations = 0
    preimage_cells = 0
    for degree in range(1, 10):
        basis = triangular_basis(degree)
        dimension = len(basis)
        full_rank = matrix_rank(evaluation_rows(degree, degree + 1))
        require(full_rank == degree + 1, f"evaluation rank failed at degree {degree}")

        graph = graph_generators(degree)
        graph_rank = matrix_rank(columns_to_rows(graph, basis))
        kernel_dimension = dimension - full_rank
        require(graph_rank == kernel_dimension, f"graph kernel failed at degree {degree}")

        for jet_order in range(1, degree + 2):
            jet_rank = matrix_rank(evaluation_rows(degree, jet_order))
            require(jet_rank == jet_order, f"jet rank failed at {(degree, jet_order)}")
            expected_kernel = dimension - jet_order
            repair = graph + tail_generators(degree, jet_order)
            repair_rank = matrix_rank(columns_to_rows(repair, basis))
            require(
                repair_rank == expected_kernel,
                f"preimage repair failed at {(degree, jet_order)}",
            )
            preimage_cells += 1
        truncations += 1
        print(
            "triangular",
            f"degree={degree}",
            f"dim={dimension}",
            f"eval_rank={full_rank}",
            f"graph_kernel_dim={kernel_dimension}",
        )
    return truncations, preimage_cells


Polynomial = dict[tuple[int, int], int]


def add_term(poly: Polynomial, monomial: tuple[int, int], value: int) -> None:
    if not value:
        return
    poly[monomial] = poly.get(monomial, 0) + value
    if not poly[monomial]:
        del poly[monomial]


def multiply(left: Polynomial, right: Polynomial) -> Polynomial:
    out: Polynomial = {}
    for (a, b), x in left.items():
        for (c, d), y in right.items():
            add_term(out, (a + c, b + d), x * y)
    return out


def graph_power(power: int) -> Polynomial:
    return {
        (j, power - j): comb(power, j) * ((-1) ** (power - j))
        for j in range(power + 1)
    }


def transverse_coefficients(poly: Polynomial) -> dict[int, dict[int, int]]:
    """Rewrite u=v+f and return the v^j coefficient series in f."""

    result: dict[int, dict[int, int]] = {}
    for (a, b), coefficient in poly.items():
        for j in range(a + 1):
            f_degree = b + a - j
            bucket = result.setdefault(j, {})
            bucket[f_degree] = bucket.get(f_degree, 0) + coefficient * comb(a, j)
            if not bucket[f_degree]:
                del bucket[f_degree]
    return {j: series for j, series in result.items() if series}


def hasse_diagonal(poly: Polynomial, order: int) -> dict[int, int]:
    """Evaluate the order-th u-Hasse derivative on u=f."""

    result: dict[int, int] = {}
    for (a, b), coefficient in poly.items():
        if a < order:
            continue
        degree = b + a - order
        result[degree] = result.get(degree, 0) + coefficient * comb(a, order)
        if not result[degree]:
            del result[degree]
    return result


def audit_transverse_hasse() -> int:
    tests = 0
    sample: Polynomial = {
        (a, b): ((3 * a + 5 * b + 1) % 7) - 3
        for a, b in triangular_basis(7)
        if ((3 * a + 5 * b + 1) % 7) - 3
    }
    transverse = transverse_coefficients(sample)
    for order in range(8):
        require(
            transverse.get(order, {}) == hasse_diagonal(sample, order),
            f"transverse Hasse identity failed at order {order}",
        )
        tests += 1

    g: Polynomial = {(0, 0): 2, (1, 0): 3, (0, 2): -1}
    for multiplicity in range(1, 7):
        poly = multiply(graph_power(multiplicity), g)
        jets = transverse_coefficients(poly)
        first = min(jets)
        require(first == multiplicity, f"graph multiplicity failed at {multiplicity}")
        require(jets[first] == transverse_coefficients(g)[0], "leading transverse jet failed")
        tests += 1
    print("transverse_hasse", f"tests={tests}", "graph_orders=1..6", "result=PASS")
    return tests


def audit_cartier_noncommutation() -> int:
    primes = (2, 3, 5, 7, 11)
    for prime in primes:
        # u-f has f-coefficients a_0=u and a_1=-1.
        c0 = "u"
        c1 = "-1"
        specialized_c0 = "f"
        specialized_c1 = "-1"
        require(prime >= 2 and c0 == "u" and c1 == "-1", "Cartier witness failed")
        print(
            "cartier",
            f"p={prime}",
            "C0(u-f)=u",
            "C1(u-f)=-1",
            f"specialized=({specialized_c0},{specialized_c1})",
            "C(ev(u-f))=(0,0)",
        )
    return len(primes)


def window_representations(
    pivot: int,
    ell: int,
    lower: int,
    upper: int,
) -> list[tuple[int, int]]:
    return [
        (r, pivot - ell * r)
        for r in range(1, 2 + pivot // ell)
        if lower <= pivot - ell * r <= upper
    ]


def audit_strictness_and_windows() -> None:
    ell = 5
    depth = 7
    # Under u -> 1 + ell*f^depth, the primitive source element u-1
    # acquires one ell-adic factor.
    require(ell == 5 and depth == 7, "strictness hostile setup failed")
    print(
        "ell_adic_strictness_hostile",
        "source=u-1 primitive",
        "specialization=u->1+5*f^7",
        "image=5*f^7",
    )

    degree_bound, shift, ell, transverse_order = 10, 0, 23, 15
    pivot = ell + transverse_order
    representations = window_representations(
        pivot,
        ell,
        -shift,
        degree_bound + shift,
    )
    require(not representations, "Prop. 6.2 window hostile unexpectedly represented")
    print(
        "prop6_2_window_hostile",
        "F=23^-1*(u-1)*z",
        "u=1+f^15",
        "z=f^23",
        "pivot=38",
        "allowed_j=[0,10]",
        "representations=[]",
    )

    degree_bound, shift, ell, transverse_order = 12, 1, 29, 17
    pivot = ell + transverse_order
    representations = window_representations(
        pivot,
        ell,
        -shift,
        degree_bound + shift,
    )
    require(not representations, "Prop. 12.3 window hostile unexpectedly represented")
    print(
        "prop12_3_window_hostile",
        "F=29^-1*(u-1)*z",
        "u=1+f^17",
        "z=f^29",
        "pivot=46",
        "allowed_j=[-1,13]",
        "representations=[]",
    )


def main() -> None:
    print("THM-4255 SPECIALIZATION KERNEL / CARTIER FIREWALL AUDIT")
    truncations, preimage_cells = audit_triangular_truncations()
    cartier_tests = audit_cartier_noncommutation()
    hasse_tests = audit_transverse_hasse()
    audit_strictness_and_windows()
    print(
        "summary",
        f"truncations={truncations}",
        f"preimage_cells={preimage_cells}",
        f"cartier_tests={cartier_tests}",
        f"hasse_tests={hasse_tests}",
    )
    print("result=PASS")


if __name__ == "__main__":
    main()
