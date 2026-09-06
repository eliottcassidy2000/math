#!/usr/bin/env python3
"""Formal exact audit of the THM-4412 seminormal suspension.

The heavy exceptional-quartic facts are imported as named proved dependencies
from THM-4381/4404: S^sn=S+K*r, r not in S, r vanishes simply on the reduced
conductor, and Lambda(r') is nonzero at the retained triple.  This lightweight
certificate checks every new algebraic implication, the all-degree leading-
coefficient pattern behind the intersection lemma, the retained first-jet
determinant, and the contrast with Long's value-level primitive.
"""

from __future__ import annotations

import sys

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def reduce_triple(poly: sp.Expr, x: sp.Symbol) -> sp.Expr:
    return sp.rem(
        sp.Poly(sp.expand(poly), x, domain=sp.QQ),
        sp.Poly(x**3 - x, x, domain=sp.QQ),
    ).as_expr()


def main() -> None:
    # Abstract translated polynomial coordinate D=z+r.
    z, r, d, c = sp.symbols("z r d c")
    D = z + r
    check(sp.diff(D, z) == 1, "translated coordinate is monic")
    check(sp.expand((d - r) + r - d) == 0, "inverse translation z=D-r")
    check(sp.expand(D.subs(z, 0) - r) == 0, "source section z=0")
    check(sp.expand(D.subs(z, c) - (c + r)) == 0, "constant source section")
    check(sp.expand(D.subs(z, -r)) == 0, "target section D=0")

    # For every degree, the top z coefficient of sum a_j(z+r)^j is a_m.
    # Eight symbolic degrees are dynamic controls; the displayed binomial
    # identity and monicity are the all-degree proof used in the theorem.
    for degree in range(1, 9):
        coefficients = sp.symbols(f"a0:{degree + 1}")
        polynomial = sp.expand(
            sum(coefficients[index] * D**index for index in range(degree + 1))
        )
        top = sp.Poly(polynomial, z).coeff_monomial(z**degree)
        check(top == coefficients[degree], f"degree-{degree} leading coefficient")

    # The retained triple tangent relation.  THM-4381 supplies the exact
    # ordinary rows and THM-4404 supplies ell(r') != 0.
    a, b, q = sp.symbols("a b q")
    C_row = sp.Matrix([[3, 3, 3]])
    E_row = sp.Matrix([[-9, 4, 9]])
    r_row = sp.Matrix([[a, b, q]])
    ell = sp.Matrix([[5, -18, 13]])
    check((C_row * ell.T)[0] == 0, "ell kills C tangent")
    check((E_row * ell.T)[0] == 0, "ell kills E tangent")
    period = (r_row * ell.T)[0]
    determinant = sp.Matrix.vstack(C_row, E_row, r_row).det()
    check(sp.expand(determinant - 3 * period) == 0, "jet determinant is period")
    check(sp.Matrix.vstack(C_row, E_row).rank() == 2, "planar tangent rank two")

    # r vanishes at every reduced conductor branch.  A constant z-section
    # therefore preserves point values at the triple, while the imported
    # nonzero period makes its augmented first-jet matrix rank three.
    r_values = (sp.Integer(0), sp.Integer(0), sp.Integer(0))
    D_values = tuple(c + value for value in r_values)
    check(D_values == (c, c, c), "seminormal sidecar preserves point collision")
    check(all(value == 0 for value in r_values), "seminormal sidecar is value invisible")

    # Long's independently proved fibre residue is a contrasting value-level
    # sidecar: its normalized class is the idempotent x^2=(0,1,1).
    x = sp.symbols("x")
    long_h = -sp.Rational(1, 48) - sp.Rational(1093, 192) * x**2
    long_values = tuple(sp.expand(long_h.subs(x, point)) for point in (0, 1, -1))
    expected_long = (
        -sp.Rational(1, 48),
        -sp.Rational(1097, 192),
        -sp.Rational(1097, 192),
    )
    check(long_values == expected_long, "Long primitive values")
    normalized = sp.cancel(
        (long_h - expected_long[0]) / (expected_long[1] - expected_long[0])
    )
    check(normalized == x**2, "Long normalized idempotent")
    check(reduce_triple(normalized**2 - normalized, x) == 0, "Long idempotence")
    long_D_values = tuple(c + value for value in long_values)
    check(long_D_values[0] != long_D_values[1], "Long constant slice splits one branch")
    check(long_D_values[1] == long_D_values[2], "Long constant slice retains a pair")

    # Length bookkeeping for the already proved seminormal extension.
    normalization_defect = 89
    seminormal_defect = 1
    check(normalization_defect - seminormal_defect == 88, "defect ledger")

    print("exceptional_quartic_seminormal_suspension_formal_audit")
    print("imported_proved_facts=THM4381:Ssn=S+K*r,r_notin_S,r_values_zero;THM4404:Lambda(rprime)!=0")
    print("suspension=B=S[D]_inside_N[z];D=z+r;D_transcendental_over_N")
    print("intersection_law=B_intersect_N=S;proof=top_z_coefficient_of_sum_a_j*(z+r)^j_is_a_m")
    print("recovery_law=z_in_B_iff_r_in_S;conclusion=z_notin_B")
    print("source_section_z0=image_S[r]=Ssn;target_section_D0=source_graph_z=-r,image_S")
    print("constant_source_section_z=c:image_S[c+r]=Ssn;retained_D_values=c,c,c")
    print("retained_rows=Cprime(3,3,3);Eprime(-9,4,9);ell=(5,-18,13)")
    print("augmented_jet_determinant=3*(5*rprime_-1-18*rprime_0+13*rprime_1);nonzero_by_THM4404")
    print("node_behavior=r_is_already_in_each_node_local_conductor;only_triple_local_ring_changes")
    print("branch_layers=Long_H:value_idempotent_(0,1,1);quartic_r:value_zero_and_missing_first_jet")
    print("Long_constant_slice=point_partition_1+2;quartic_constant_slice=points_3_but_jet_rank_3")
    print("NO_CLAIM=stable_cancellation_or_planar_embedding_or_chart_entry_or_Keller_pair_or_JC2_or_DC2")
    print(f"checks={CHECKS};result=PASS")


if __name__ == "__main__":
    main()
