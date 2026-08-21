#!/usr/bin/env python3
"""Finite-exact quadratic-projection and collision audit for THM-3556.

The explicit cusp-square packet is Z=(L,T,U_*,S).  This script checks an
exact collision over Q(sqrt(5)) and then studies every target polynomial of
total degree at most two.  Constants do not affect a Jacobian, so the
nonconstant monomial pullbacks form a 14-dimensional source-polynomial space
W.  The coefficient
matrix of all 91 brackets on wedge^2(W) gives a linear relaxation of the
Pluecker/decomposability conditions.  Even that relaxation cannot produce a
nonzero constant.
"""

from __future__ import annotations

from itertools import combinations

import sympy as sp
from sympy.polys.domains import GF
from sympy.polys.matrices import DomainMatrix


def require(condition: bool, label: str) -> None:
    """Keep every truth-bearing check active under ``python -O``."""
    if not bool(condition):
        raise ArithmeticError(f"FAILED: {label}")


def exact_rank(matrix: sp.Matrix) -> int:
    """Exact rank over Q using Sympy's fraction-field row reduction."""
    return len(DomainMatrix.from_Matrix(matrix).to_field().rref()[1])


def modular_rank(matrix: sp.Matrix, prime: int) -> int:
    """Independent rank control after coefficientwise reduction mod prime."""
    def reduce_entry(entry: sp.Expr) -> int:
        rational = sp.Rational(entry)
        numerator = int(rational.p) % prime
        denominator = int(rational.q) % prime
        require(denominator != 0, f"denominator survives modulo {prime}")
        return numerator * pow(denominator, prime - 2, prime) % prime

    reduced = matrix.applyfunc(reduce_entry)
    domain_matrix = DomainMatrix.from_Matrix(reduced).convert_to(GF(prime))
    return len(domain_matrix.rref()[1])


def main() -> None:
    v, y = sp.symbols("v y")
    U = sp.expand(
        1 + y - sp.Rational(1, 2) * y**2
        - sp.Rational(3, 2) * v * y * (y - 3)
    )
    T = sp.expand(y**2 - 6 * v * U)
    S = sp.expand(y**3 - 9 * v * U * y)
    L = sp.expand(v**2 * (8 * v * U - y**2))
    packet = [L, T, U, S]

    # The two conjugate source points are distinct over Q(sqrt(5)) and have
    # exactly the same four packet coordinates.
    root_five = sp.sqrt(5)
    point_plus = (
        sp.Rational(8, 27) + sp.Rational(4, 9) * root_five,
        sp.Rational(3, 2) + sp.Rational(3, 5) * root_five,
    )
    point_minus = (
        sp.Rational(8, 27) - sp.Rational(4, 9) * root_five,
        sp.Rational(3, 2) - sp.Rational(3, 5) * root_five,
    )
    expected_collision = [
        -sp.Rational(6724, 3645),
        sp.Rational(57, 20),
        sp.Rational(27, 40),
        sp.Rational(27, 40),
    ]
    collision_ideal = sp.groebner(
        [27 * v - 20 * y + 22,
         y**2 - 3 * y + sp.Rational(9, 20)],
        v, y,
        order="lex",
    )
    require(sp.discriminant(y**2 - 3 * y + sp.Rational(9, 20), y)
            == sp.Rational(36, 5),
            "collision scheme has two distinct geometric points")
    for entry, expected in zip(packet, expected_collision):
        require(collision_ideal.reduce(entry - expected)[1] == 0,
                "packet coordinate is constant on collision scheme")
    require(
        sp.simplify(point_plus[0] - point_minus[0]) != 0
        and sp.simplify(point_plus[1] - point_minus[1]) != 0,
        "collision points are distinct",
    )
    for sign, point in [("plus", point_plus), ("minus", point_minus)]:
        values = [sp.simplify(entry.subs({v: point[0], y: point[1]}))
                  for entry in packet]
        require(values == expected_collision,
                f"direct four-coordinate collision substitution ({sign})")

    # Positive first-jet control: the collision does not itself make the two
    # tangent-plane unit equations inconsistent.  This affine projection has
    # Jacobian one at both sheets, although it is not globally Keller.
    local_B = (
        sp.Rational(9627, 304384) * T
        - sp.Rational(205687, 1369728) * U
    )
    local_jacobian = sp.expand(
        sp.diff(L, v) * sp.diff(local_B, y)
        - sp.diff(L, y) * sp.diff(local_B, v)
    )
    require(collision_ideal.reduce(local_jacobian - 1)[1] == 0,
            "simultaneous unit first-jet control at both collision sheets")

    # These are exactly the nonconstant target monomials of coordinate
    # total degree <=2, pulled back along Z, in deterministic order.
    names = ["L", "T", "U", "S"]
    basis = list(packet)
    basis_names = list(names)
    for first in range(4):
        for second in range(first, 4):
            basis.append(sp.expand(packet[first] * packet[second]))
            basis_names.append(names[first] + names[second])
    require(len(basis) == 14, "quadratic pullback basis length")

    composition_monomials = sorted(
        set().union(*[set(sp.Poly(entry, v, y).monoms()) for entry in basis]),
        key=lambda exponent: (sum(exponent), exponent),
    )
    composition_matrix = sp.Matrix([
        [sp.Poly(entry, v, y).coeff_monomial(monomial) for entry in basis]
        for monomial in composition_monomials
    ])
    composition_rank = exact_rank(composition_matrix)
    require(composition_rank == 14,
            "14 target monomial pullbacks are linearly independent")

    # If P=sum a_i f_i and Q=sum b_i f_i, then
    # Jac(P,Q)=sum_(i<j)(a_i b_j-a_j b_i) Jac(f_i,f_j).  We allow the 91
    # wedge coefficients to vary independently.  This forgets every
    # Pluecker and integrability constraint and is therefore a relaxation.
    pairs = list(combinations(range(14), 2))
    derivatives_v = [sp.Poly(sp.diff(entry, v), v, y) for entry in basis]
    derivatives_y = [sp.Poly(sp.diff(entry, y), v, y) for entry in basis]
    brackets = [
        derivatives_v[first] * derivatives_y[second]
        - derivatives_y[first] * derivatives_v[second]
        for first, second in pairs
    ]
    require(len(brackets) == 91, "wedge-square column count")
    zero_bracket_names = [
        basis_names[first] + "^" + basis_names[second]
        for (first, second), bracket in zip(pairs, brackets)
        if bracket.is_zero
    ]
    require(zero_bracket_names == ["L^LL", "T^TT", "U^UU", "S^SS"],
            "four tautological square-bracket zero columns")

    bracket_monomials = sorted(
        set().union(*[set(entry.monoms()) for entry in brackets]) | {(0, 0)},
        key=lambda exponent: (sum(exponent), exponent),
    )
    coefficient_matrix = sp.Matrix([
        [entry.coeff_monomial(monomial) for entry in brackets]
        for monomial in bracket_monomials
    ])
    require(coefficient_matrix.shape == (139, 91),
            "Jacobian coefficient matrix shape")
    constant_index = bracket_monomials.index((0, 0))
    nonconstant_matrix = coefficient_matrix.copy()
    nonconstant_matrix.row_del(constant_index)

    full_rank = exact_rank(coefficient_matrix)
    nonconstant_rank = exact_rank(nonconstant_matrix)
    require(full_rank == 67, "full exact Jacobian coefficient rank")
    require(nonconstant_rank == 67,
            "nonconstant-row exact Jacobian coefficient rank")

    modular_controls = []
    for prime in (1000003, 1000033):
        modular_full = modular_rank(coefficient_matrix, prime)
        modular_nonconstant = modular_rank(nonconstant_matrix, prime)
        require(modular_full == 67 and modular_nonconstant == 67,
                f"modular rank control at {prime}")
        modular_controls.append((prime, modular_full, modular_nonconstant))

    # Equal ranks mean that the constant row lies in the span of the
    # nonconstant rows.  Thus any wedge coefficient vector killing every
    # nonconstant source monomial also kills the constant coefficient.
    require(full_rank == nonconstant_rank,
            "constant row is redundant modulo nonconstant rows")

    # A still broader relaxation starts directly from the six packet minors.
    # For quadratic A and B, every ambient coefficient A_i B_j-A_j B_i has
    # total target degree <=2.  Forget Pluecker decomposability and polynomial-
    # potential realizability completely: allow each of the six minors an
    # arbitrary coefficient in span{1,W}.  Even these 15*6 descending columns
    # cannot produce a nonzero constant.
    packet_v = [sp.Poly(sp.diff(entry, v), v, y) for entry in packet]
    packet_y = [sp.Poly(sp.diff(entry, y), v, y) for entry in packet]
    packet_minors = [
        packet_v[first] * packet_y[second]
        - packet_y[first] * packet_v[second]
        for first, second in combinations(range(4), 2)
    ]
    coefficient_pullbacks = [sp.Poly(1, v, y)] + [sp.Poly(entry, v, y)
                                                   for entry in basis]
    descended_columns = [
        coefficient * minor
        for coefficient in coefficient_pullbacks
        for minor in packet_minors
    ]
    require(len(descended_columns) == 90,
            "quadratic descending-coefficient relaxation column count")
    descended_monomials = sorted(
        set().union(*[set(entry.monoms()) for entry in descended_columns])
        | {(0, 0)},
        key=lambda exponent: (sum(exponent), exponent),
    )
    descended_matrix = sp.Matrix([
        [entry.coeff_monomial(monomial) for entry in descended_columns]
        for monomial in descended_monomials
    ])
    require(descended_matrix.shape == (139, 90),
            "descending-coefficient matrix shape")
    descended_constant_index = descended_monomials.index((0, 0))
    descended_nonconstant = descended_matrix.copy()
    descended_nonconstant.row_del(descended_constant_index)
    descended_full_rank = exact_rank(descended_matrix)
    descended_nonconstant_rank = exact_rank(descended_nonconstant)
    require(descended_full_rank == 80 and descended_nonconstant_rank == 80,
            "descending-coefficient exact ranks")
    for prime in (1000003, 1000033):
        require(
            modular_rank(descended_matrix, prime) == 80
            and modular_rank(descended_nonconstant, prime) == 80,
            f"descending-coefficient modular rank control at {prime}",
        )

    # Weight ledger used by the cited GGHV degree-floor transfer.  At target
    # cap 18, only L^18 can reach source degree 108.  The common-leading-base
    # monomial v^48*y^24 can likewise arise only from L^12.
    source_degrees = tuple(sp.Poly(entry, v, y).total_degree()
                           for entry in packet)
    require(source_degrees == (6, 4, 3, 5),
            "packet source total-degree weights")
    leading_exponents = ((4, 2), (2, 2), (1, 2), (2, 3))
    leading_terms = []
    for entry, degree in zip(packet, source_degrees):
        polynomial = sp.Poly(entry, v, y)
        leading_terms.append(sp.expand(sum(
            coefficient * v**monomial[0] * y**monomial[1]
            for monomial, coefficient in polynomial.terms()
            if sum(monomial) == degree
        )))
    require(
        leading_terms == [
            -12 * v**4 * y**2,
            9 * v**2 * y**2,
            -sp.Rational(3, 2) * v * y**2,
            sp.Rational(27, 2) * v**2 * y**3,
        ],
        "packet leading monomials",
    )
    cap_eighteen = [
        (a, b, c, d)
        for a in range(19)
        for b in range(19 - a)
        for c in range(19 - a - b)
        for d in range(19 - a - b - c)
    ]
    height_108 = [
        exponent for exponent in cap_eighteen
        if sum(weight * power for weight, power
               in zip(source_degrees, exponent)) == 108
    ]
    require(height_108 == [(18, 0, 0, 0)],
            "unique cap-18 source-height-108 target monomial")
    common_top_72 = [
        exponent for exponent in cap_eighteen
        if (
            sum(pair[0] * power for pair, power
                in zip(leading_exponents, exponent)) == 48
            and sum(pair[1] * power for pair, power
                    in zip(leading_exponents, exponent)) == 24
        )
    ]
    require(common_top_72 == [(12, 0, 0, 0)],
            "unique cap-18 pullback with leading monomial v^48*y^24")

    print("THM-3556 QUADRATIC TARGET-PROJECTION AUDIT")
    print("base_field=Q; extension_scope=every characteristic-zero field")
    print(f"collision_points={point_plus};{point_minus}")
    print(f"collision_packet_L_T_U_S={tuple(expected_collision)}")
    print("collision_first_jet_control=Jac(L,(9627/304384)T"
          "-(205687/1369728)U)=1_on_both_sheets")
    print(f"target_degree_nonconstant_basis={basis_names}")
    print(f"composition_dimension={composition_rank}")
    print("wedge_columns=binomial(14,2)=91")
    print(f"zero_wedge_columns={zero_bracket_names}")
    print(f"coefficient_matrix_shape={coefficient_matrix.shape}")
    print(f"exact_ranks_full_nonconstant={full_rank},{nonconstant_rank}")
    print(f"modular_rank_controls={modular_controls}")
    print("constant_row_in_nonconstant_row_span=True")
    print(f"descent_relaxation_shape={descended_matrix.shape}")
    print("descent_relaxation_exact_ranks_full_nonconstant="
          f"{descended_full_rank},{descended_nonconstant_rank}")
    print(f"packet_source_degrees={source_degrees}")
    print(f"cap18_height108_exponents={height_108}")
    print(f"cap18_common_top_v48_y24_exponents={common_top_72}")
    print("VERDICT: no target-total-degree-<=2 pair has nonzero constant Jacobian.")
    print("SCOPE: the total-degree-3-and-higher nonlinear projection problem remains open here.")


if __name__ == "__main__":
    main()
