#!/usr/bin/env python3
"""Primary exact certificate for THM-4395.

Continue THM-4390's complete fixed residual-weight-at-most-fourteen family
through row ten.  The row-ten bracket and projected P_2/P_3 depth equations
form a global triangular affine graph, including at Phi=0.
"""

from __future__ import annotations

from pathlib import Path
import sys

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import jc2_source_normal_weight14_row9_face_absorption_thm4390 as R9  # noqa: E402


R8 = R9.R8
x = R9.x
CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def exact(value: sp.Expr) -> sp.Expr:
    return R9.exact(value)


def proportional(left: sp.Expr, right: sp.Expr) -> bool:
    return R9.proportional(left, right)


def build_row_nine_state() -> tuple[
    list[sp.Expr],
    list[sp.Expr],
    tuple[sp.Symbol, ...],
    dict[sp.Symbol, sp.Expr],
    dict[int, sp.Expr],
    dict[sp.Symbol, sp.Expr],
    sp.Expr,
]:
    arows, crows, bracket_subs, grows = R9.build_through_row_eight()
    theta8 = tuple(sp.symbols("w14_theta8_0:9"))
    depth8, _, _ = R9.projected_depth_residuals(arows, crows, 8)
    terminal8, depth8_conditions, _, _ = R9.eliminate_linear_fibre(depth8, theta8)
    gate = {
        R8.Delta: sp.Rational(896, 15),
        R8.Theta: sp.Rational(512, 75),
        R8.zeta3: -sp.Rational(3, 2) * R8.Phi,
    }
    check(len(depth8_conditions) == 3, "inherited row-eight gate")
    check(all(exact(value.subs(gate)) == 0 for value in depth8_conditions),
          "row-eight gate vanishes")

    a8 = R9.impose(arows, gate, terminal8, gate)
    c8 = R9.impose(crows, gate, terminal8, gate)
    g9 = exact(grows[9].subs(bracket_subs).subs(gate))
    select9, conditions9, matrix9, pivots9 = R9.solve_bracket_fibre(
        a8, c8, 9, g9, theta8[:7]
    )
    check((matrix9.shape, matrix9.rank(), pivots9) == ((10, 7), 7, tuple(range(7))),
          "row-nine selector")
    check(len(conditions9) == 1, "row-nine condition count")
    face14 = 3339765000 * R9.c70 + 898128000 * R9.c42 + 216513000 * R9.c14
    row9_equation = sp.expand(R9.R9.E9_GATE - face14)
    check(proportional(conditions9[0], row9_equation), "THM-4390 row-nine graph")
    c14_graph = exact(sp.solve(row9_equation, R9.c14)[0])

    a8_selected = R9.impose(a8, select9, {R9.c14: c14_graph})
    c8_selected = R9.impose(c8, select9, {R9.c14: c14_graph})
    a9, c9, theta9 = R9.append_tangent(a8_selected, c8_selected, 9, "w14r10_theta9")
    depth9, _, _ = R9.projected_depth_residuals(a9, c9, 9)
    terminal9, depth9_conditions, matrix_depth9, pivots_depth9 = R9.eliminate_linear_fibre(
        depth9, theta9
    )
    check(depth9_conditions == [], "row-nine depth automatic")
    check((matrix_depth9.shape, matrix_depth9.rank(), pivots_depth9)
          == ((28, 10), 3, (7, 8, 9)), "row-nine terminal depth")
    return (
        R9.impose(a9, terminal9),
        R9.impose(c9, terminal9),
        theta9,
        bracket_subs,
        grows,
        gate,
        c14_graph,
    )


def main() -> None:
    a9, c9, theta9, bracket_subs, grows, gate, c14_graph = build_row_nine_state()
    g10 = exact(
        grows[10]
        .subs(bracket_subs)
        .subs(gate)
        .subs(R9.c14, c14_graph)
    )
    select10, bracket10, matrix10, pivots10 = R9.solve_bracket_fibre(
        a9, c9, 10, g10, theta9[:7]
    )
    check((matrix10.shape, matrix10.rank(), pivots10)
          == ((11, 7), 7, tuple(range(7))), "row-ten bracket selector")
    check(len(bracket10) == 2, "two row-ten bracket conditions")

    xi_equation = (
        13365000 * R8.Phi**2
        + 15035625 * R8.Phi * R8.eta
        + 6014250 * R9.c42
        + 50787000 * R9.c70
        + 57672000 * R8.xi10
        - 964604821504
    )
    xi_candidates = [value for value in bracket10 if proportional(value, xi_equation)]
    check(len(xi_candidates) == 1, "global xi equation")
    remaining_candidates = [value for value in bracket10 if value is not xi_candidates[0]]
    check(len(remaining_candidates) == 1, "remaining bracket equation")
    xi_graph = exact(sp.solve(xi_equation, R8.xi10)[0])
    check(sp.fraction(sp.together(xi_graph))[1] == 57672000,
          "constant xi denominator")

    bracket_remaining = (
        -104916222000 * R8.Phi**2
        + 122625090000 * R8.Phi * R8.alpha11
        + 19707603750 * R8.Phi * R8.beta11
        + 20802470625 * R8.Phi * R9.c51
        - 246422138625 * R8.Phi * R8.eta
        + 20802470625 * R8.alpha11 * R8.eta
        - 89131914000 * R9.c42
        - 194981256000 * R9.c70
        + 61312545000 * R8.eta**2
        + 2707389207937024
    )
    check(
        proportional(remaining_candidates[0].subs(R8.xi10, xi_graph), bracket_remaining),
        "remaining bracket equation after xi",
    )

    selected_a9 = R9.impose(a9, select10, {R8.xi10: xi_graph})
    selected_c9 = R9.impose(c9, select10, {R8.xi10: xi_graph})
    a10, c10, theta10 = R9.append_tangent(selected_a9, selected_c9, 10, "w14r10_theta10")
    depth10, amatrix10, cmatrix10 = R9.projected_depth_residuals(a10, c10, 10)
    terminal10, depth10_conditions, matrix_depth10, pivots_depth10 = R9.eliminate_linear_fibre(
        depth10, theta10
    )
    check((amatrix10.shape, amatrix10.rank()) == ((88, 193), 68),
          "row-ten P2 universe")
    check((cmatrix10.shape, cmatrix10.rank()) == ((99, 304), 83),
          "row-ten P3 universe")
    check((matrix_depth10.shape, matrix_depth10.rank(), pivots_depth10)
          == ((36, 11), 3, (8, 9, 10)), "row-ten terminal depth")
    check(len(depth10_conditions) == 1, "one row-ten depth condition")

    beta_equation = -91 * R8.Phi + 15 * R8.beta11 + 18 * R8.eta
    check(proportional(depth10_conditions[0], beta_equation), "global beta depth equation")
    beta_graph = exact(sp.solve(beta_equation, R8.beta11)[0])
    check(beta_graph == (91 * R8.Phi - 18 * R8.eta) / 15,
          "beta graph formula")

    final_bracket = exact(bracket_remaining.subs(R8.beta11, beta_graph))
    check(sp.diff(final_bracket, R9.c42) == -89131914000,
          "global c42 pivot")
    c42_graph = exact(sp.solve(final_bracket, R9.c42)[0])
    check(sp.fraction(sp.together(c42_graph))[1] == 89131914000,
          "constant c42 denominator")

    xi_final = exact(xi_graph.subs(R9.c42, c42_graph))
    c14_final = exact(
        c14_graph.subs(R8.xi10, xi_final).subs(R9.c42, c42_graph)
    )
    final_graph = {
        R8.beta11: beta_graph,
        R9.c42: c42_graph,
        R8.xi10: xi_final,
        R9.c14: c14_final,
    }
    check(all(variable not in expression.free_symbols for variable, expression in final_graph.items()),
          "triangular graph has no self-dependence")
    check(all(not sp.fraction(sp.together(expression))[1].free_symbols
              for expression in final_graph.values()),
          "all graph denominators are constant units")
    check(all(R9.c23 not in expression.free_symbols for expression in final_graph.values()),
          "c23 remains free")

    difference10 = exact(g10 - R8.predicted_G(10, a9, c9))
    check(
        all(
            exact(
                R9.xcoeff(difference10, degree)
                .subs(select10)
                .subs(R8.xi10, xi_final)
                .subs(R8.beta11, beta_graph)
                .subs(R9.c42, c42_graph)
            )
            == 0
            for degree in range(11)
        ),
        "all row-ten bracket coefficients vanish",
    )
    check(
        all(
            exact(
                value.subs(terminal10)
                .subs(R8.beta11, beta_graph)
                .subs(R9.c42, c42_graph)
            )
            == 0
            for value in depth10
        ),
        "all row-ten depth equations vanish",
    )

    free_control = {
        R8.Phi: 1,
        R8.eta: 2,
        R8.alpha11: 3,
        R9.c51: 5,
        R9.c23: 7,
        R9.c70: 11,
    }
    control_beta = exact(beta_graph.subs(free_control))
    control_c42 = exact(c42_graph.subs(free_control))
    control_xi = exact(xi_final.subs(free_control))
    control_c14 = exact(c14_final.subs(free_control))
    check(all(value != 0 for value in (control_beta, control_c42, control_xi, control_c14)),
          "dense rational graph control")
    old_elliptic = (
        7231154026500 * R8.Phi**3
        + 50541940696500 * R8.Phi**2 * R8.eta
        + 6793915500000 * R8.Phi * R8.eta**2
        + 353642000625 * R8.eta**3
        - 631918028977864704 * R8.Phi
        - 91584545734393856 * R8.eta
    )
    old_weight13 = old_elliptic - 707284001250 * R8.Phi**2 * R9.c51
    check(exact(old_elliptic.subs(free_control)) != 0, "hostile to old elliptic carrier")
    check(exact(old_weight13.subs(free_control)) != 0,
          "hostile to old weight-thirteen row-ten depth carrier")

    print("THM-4395 weight-fourteen row-ten global affine absorption")
    print("imports=THM4390_primary_helpers;inherited_operators=THM4308_R8_and_THM4315_R9")
    print("universe=complete_fixed_residual_weight_at_most_14")
    print("row9=c14_global_graph;row9_depth=automatic")
    print("row10_bracket=shape(11,7),rank7,pivots(0..6),conditions2")
    print("row10_xi_equation=" + sp.sstr(xi_equation))
    print("row10_depth=shape(36,11),rank3,pivots(8,9,10),condition=-91*Phi+15*beta11+18*eta")
    print("row10_beta_graph=" + sp.sstr(beta_graph))
    print("row10_remaining_bracket=" + sp.sstr(final_bracket))
    print("row10_c42_graph=" + sp.sstr(c42_graph))
    print("global_source_projection=A6_coordinates=(Phi,eta,alpha11,c51,c23,c70)")
    print("terminal_fibre=A8;boundary_localization=none;Phi_zero_included=yes")
    print(
        "rational_control=(Phi,eta,alpha11,c51,c23,c70,beta11,c42,xi10,c14)="
        + sp.sstr(
            (1, 2, 3, 5, 7, 11, control_beta, control_c42, control_xi, control_c14)
        )
    )
    print("old_elliptic_and_weight13_carriers=both_nonzero_at_control")
    print("scope=fixed_source_normal_chart_through_row10_only;row11_and_entry_open")
    print(f"checks={CHECKS};result=PASS")


if __name__ == "__main__":
    main()
