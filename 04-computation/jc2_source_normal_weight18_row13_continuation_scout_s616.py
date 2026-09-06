#!/usr/bin/env python3
"""Exact row-thirteen scout beyond THM-4403.

Continue the restricted source

    H14 + lambda18*p**3*y**4 + nu18*y**6

along the global row-eleven and row-twelve source graphs.  Compute the full
row-thirteen bracket quotient, then audit every exact-valuation-thirteen
source monomial p**(13-2*b)*y**b, b=0,...,6.  If the bracket equations admit
a constant-pivot source continuation, append the complete row-thirteen
tangent and compute its projected P_2/P_3 depth quotient.

This imports the THM-4399/4403 primary implementations only as frozen row
operators and inherited state builders.  It is a restricted-family scout,
not a complete residual-weight-at-most-twenty computation.
"""

from __future__ import annotations

from pathlib import Path
import sys

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import jc2_source_normal_row11_late_weight18_response_thm4399 as Q11  # noqa: E402


R9 = Q11.R9
R8 = Q11.R8
x = R9.x
CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def exact(value: sp.Expr) -> sp.Expr:
    return R9.exact(value)


def primitive(value: sp.Expr) -> sp.Expr:
    return R9.primitive(value)


def proportional(left: sp.Expr, right: sp.Expr) -> bool:
    return R9.proportional(left, right)


def resolve(value: sp.Expr, mapping: dict[sp.Symbol, sp.Expr]) -> sp.Expr:
    answer = exact(value)
    for _ in range(12):
        updated = exact(answer.subs(mapping, simultaneous=False))
        if updated == answer:
            return answer
        answer = updated
    raise AssertionError("source graph did not resolve")


def unique_primitives(values: list[sp.Expr]) -> list[sp.Expr]:
    answer: list[sp.Expr] = []
    for value in values:
        value = exact(value)
        if value == 0:
            continue
        candidate = primitive(value)
        if all(not proportional(candidate, old) for old in answer):
            answer.append(candidate)
    return answer


def build_row_eleven_state():
    """Rebuild THM-4399's globally corrected row-eleven target state."""
    (
        a10,
        c10,
        free10,
        g11,
        grows,
        bracket_subs,
        gate,
        graph10,
    ) = Q11.row_ten_state()
    _, base_conditions, _, _ = R9.solve_bracket_fibre(
        a10, c10, 11, g11, free10
    )
    check(len(base_conditions) == 1, "inherited row-eleven principal obstruction")
    F = primitive(base_conditions[0])
    lambda18 = sp.symbols("lambda18")
    g11_late = exact(g11 + lambda18 * x**4)
    select11, conditions11, matrix11, pivots11 = R9.solve_bracket_fibre(
        a10, c10, 11, g11_late, free10
    )
    check(
        (matrix11.shape, matrix11.rank(), pivots11, len(conditions11))
        == ((12, 8), 8, tuple(range(8)), 1),
        "inherited row-eleven late bracket geometry",
    )
    corrected_F = primitive(conditions11[0])
    lambda_response = exact(sp.diff(corrected_F, lambda18))
    check(lambda_response != 0 and not lambda_response.free_symbols,
          "inherited lambda18 constant pivot")
    lambda_graph = exact(sp.solve(corrected_F, lambda18)[0])
    a10_selected = R9.impose(a10, select11, {lambda18: lambda_graph})
    c10_selected = R9.impose(c10, select11, {lambda18: lambda_graph})
    a11_trial, c11_trial, theta11 = R9.append_tangent(
        a10_selected, c10_selected, 11, "row13_scout_theta11"
    )
    depth11, _, _ = R9.projected_depth_residuals(a11_trial, c11_trial, 11)
    terminal11, conditions_depth11, matrix_depth11, pivots_depth11 = (
        R9.eliminate_linear_fibre(depth11, theta11)
    )
    check(not conditions_depth11, "inherited row-eleven depth automatic")
    check(
        (matrix_depth11.shape, matrix_depth11.rank(), pivots_depth11)
        == ((45, 12), 4, (8, 9, 10, 11)),
        "inherited row-eleven depth geometry",
    )
    a11 = R9.impose(a11_trial, terminal11, {lambda18: lambda_graph})
    c11 = R9.impose(c11_trial, terminal11, {lambda18: lambda_graph})
    free11 = tuple(
        symbol for index, symbol in enumerate(theta11)
        if index not in pivots_depth11
    )
    return (
        a11,
        c11,
        free11,
        F,
        lambda18,
        lambda_graph,
        grows,
        bracket_subs,
        gate,
        graph10,
    )


def build_row_twelve_state():
    (
        a11,
        c11,
        free11,
        F,
        lambda18,
        lambda_graph,
        grows,
        bracket_subs,
        gate,
        graph10,
    ) = build_row_eleven_state()
    check(len(free11) == 8, "row-eleven terminal affine eight")

    nu18 = sp.symbols("nu18")
    p3y4_row12 = exact(R9.tcoeff(R8.p**3 * R8.y**4, 12))
    y6_row12 = exact(R9.tcoeff(R8.y**6, 12))
    g12_base = exact(grows[12].subs(bracket_subs).subs(gate).subs(graph10))
    g12 = exact(g12_base + lambda_graph * p3y4_row12 + nu18 * y6_row12)
    select12, conditions12, matrix12, pivots12 = R9.solve_bracket_fibre(
        a11, c11, 12, g12, free11
    )
    check(
        (matrix12.shape, matrix12.rank(), pivots12, len(conditions12))
        == ((13, 8), 8, tuple(range(8)), 2),
        "inherited row-twelve bracket geometry",
    )
    primitive_conditions = [primitive(value) for value in conditions12]
    pivot_candidates = [
        value for value in primitive_conditions if sp.diff(value, nu18) != 0
    ]
    side_candidates = [
        value for value in primitive_conditions if sp.diff(value, nu18) == 0
    ]
    check(
        len(pivot_candidates) == len(side_candidates) == 1,
        "row-twelve triangular conditions",
    )
    side = side_candidates[0]
    check(
        sp.diff(side, R9.c23) == -675,
        "row-twelve inherited c23 constant pivot",
    )
    c23_graph = exact(sp.solve(side, R9.c23)[0])
    pivot_on_side = exact(pivot_candidates[0].subs(R9.c23, c23_graph))
    nu_graph = exact(sp.solve(pivot_on_side, nu18)[0])

    a11_selected = R9.impose(
        a11, select12, {R9.c23: c23_graph}, {nu18: nu_graph}
    )
    c11_selected = R9.impose(
        c11, select12, {R9.c23: c23_graph}, {nu18: nu_graph}
    )
    a12_trial, c12_trial, theta12 = R9.append_tangent(
        a11_selected, c11_selected, 12, "row13_scout_theta12"
    )
    depth12, amatrix12, cmatrix12 = R9.projected_depth_residuals(
        a12_trial, c12_trial, 12
    )
    terminal12, depth_conditions12, depth_matrix12, depth_pivots12 = (
        R9.eliminate_linear_fibre(depth12, theta12)
    )
    check(
        (amatrix12.shape, amatrix12.rank(), cmatrix12.shape, cmatrix12.rank())
        == ((117, 267), 87, (130, 424), 105),
        "row-twelve projected universes",
    )
    check(
        (
            depth_matrix12.shape,
            depth_matrix12.rank(),
            depth_pivots12,
            len(depth_conditions12),
        )
        == ((55, 13), 4, (9, 10, 11, 12), 1),
        "row-twelve terminal geometry",
    )
    depth_condition = primitive(depth_conditions12[0])
    check(
        sp.diff(depth_condition, R9.c70) == -368130186720000,
        "row-twelve inherited c70 constant pivot",
    )
    c70_graph = exact(sp.solve(depth_condition, R9.c70)[0])

    source_graph = {
        R9.c23: c23_graph,
        R9.c70: c70_graph,
        lambda18: lambda_graph,
        nu18: nu_graph,
    }
    source_graph = {
        symbol: resolve(value, source_graph) for symbol, value in source_graph.items()
    }
    base = {R8.Phi, R8.eta, R8.alpha11, R9.c51}
    check(
        set().union(*(value.free_symbols for value in source_graph.values())) <= base,
        "row-twelve graph resolves over affine four",
    )
    check(
        all(not sp.fraction(sp.together(value))[1].free_symbols
            for value in source_graph.values()),
        "row-twelve graph has constant denominators",
    )

    a12 = R9.impose(a12_trial, terminal12, {R9.c70: c70_graph})
    c12 = R9.impose(c12_trial, terminal12, {R9.c70: c70_graph})
    a12 = [resolve(value, source_graph) for value in a12]
    c12 = [resolve(value, source_graph) for value in c12]
    free12 = tuple(
        symbol for index, symbol in enumerate(theta12)
        if index not in depth_pivots12
    )
    check(len(free12) == 9, "row-twelve terminal affine nine")
    target_symbols = set().union(*(value.free_symbols for value in a12 + c12))
    unexpected = target_symbols - base - set(free12) - {x}
    check(
        not unexpected,
        "row-twelve target state has only base and terminal coordinates: "
        + ",".join(map(str, sorted(unexpected, key=str))),
    )

    g13_base = exact(grows[13].subs(bracket_subs).subs(gate).subs(graph10))
    g13 = exact(
        g13_base
        + lambda18 * R9.tcoeff(R8.p**3 * R8.y**4, 13)
        + nu18 * R9.tcoeff(R8.y**6, 13)
    )
    g13 = resolve(g13, source_graph)
    return (
        a12,
        c12,
        free12,
        g13,
        source_graph,
        F,
        lambda18,
        nu18,
    )


def main() -> None:
    (
        a12,
        c12,
        free12,
        g13,
        source_graph,
        F,
        lambda18,
        nu18,
    ) = build_row_twelve_state()

    difference13 = exact(g13 - R8.predicted_G(13, a12, c12))
    equations13 = [R9.xcoeff(difference13, degree) for degree in range(14)]
    check(
        exact(difference13 - sum(equations13[j] * x**j for j in range(14))) == 0,
        "row-thirteen coefficient universe exhaustive",
    )
    select13, conditions13, matrix13, pivots13 = R9.solve_bracket_fibre(
        a12, c12, 13, g13, free12
    )
    reduced13 = [exact(value.subs(select13)) for value in equations13]
    primitives13 = unique_primitives(reduced13)
    raw_cokernel13 = matrix13.T.nullspace()
    check(len(raw_cokernel13) == 5, "row-thirteen raw cokernel dimension five")
    check(len(primitives13) == 1, "row-thirteen selected ideal is principal")

    # Exact-valuation-thirteen diagonal.  Every channel is firewalled from
    # rows <=12 and has leading row x**b.
    rho = sp.symbols("rho13_0:7")
    monomials13 = tuple(R8.p ** (13 - 2*b) * R8.y**b for b in range(7))
    weights13 = tuple(26 - b for b in range(7))
    leading13 = tuple(exact(R9.tcoeff(value, 13)) for value in monomials13)
    check(leading13 == tuple(x**b for b in range(7)),
          "valuation-thirteen leading diagonal")
    check(
        all(exact(R9.tcoeff(value, row)) == 0
            for value in monomials13 for row in range(13)),
        "valuation-thirteen firewall",
    )
    perturbed_g13 = exact(g13 + sum(rho[b] * leading13[b] for b in range(7)))
    select_all, conditions_all, matrix_all, pivots_all = R9.solve_bracket_fibre(
        a12, c12, 13, perturbed_g13, free12
    )
    check(matrix_all == matrix13 and pivots_all == pivots13,
          "same row-thirteen tangent quotient")
    difference_all = exact(perturbed_g13 - R8.predicted_G(13, a12, c12))
    equations_all = [R9.xcoeff(difference_all, degree) for degree in range(14)]
    reduced_all = [exact(value.subs(select_all)) for value in equations_all]
    response_full = sp.Matrix([
        [exact(sp.diff(value, rho[b])) for b in range(7)]
        for value in reduced_all
    ])

    # The row-specific Student observer gives the scalar source response.  It
    # is compared to the full selected quotient to expose extra sidecars.
    student13 = R9.R9.primitive_student_row(13)
    student_base = exact(sum(student13[j] * equations13[j] for j in range(14)))
    student_responses = tuple(exact(sum(
        student13[j] * sp.diff(equations_all[j], rho[b]) for j in range(14)
    )) for b in range(7))
    check(response_full.rank() == 1,
          "valuation-thirteen full response has rank one")
    check(tuple(response_full[:, b].is_zero_matrix for b in range(7))
          == (False, True, False, True, False, True, False),
          "valuation-thirteen parity boundary")
    check(proportional(primitives13[0], primitive(student_base)),
          "active bracket condition is the Student coordinate")
    response_ratios = tuple(exact(value / student_responses[6])
                            for value in student_responses)
    for b in range(7):
        check(response_full[:, b] == response_ratios[b] * response_full[:, 6],
              f"full response b={b} equals Student ratio")

    # Use the minimum-weight Student-visible channel p*y^6 (weight 20) to
    # pay the unique bracket condition, then inspect the complete projected
    # depth quotient.  The other six channel variables are set to zero in
    # this first triangular continuation.
    kappa20 = rho[6]
    kappa_source = exact(g13 + kappa20 * x**6)
    select_kappa, conditions_kappa, matrix_kappa, pivots_kappa = (
        R9.solve_bracket_fibre(a12, c12, 13, kappa_source, free12)
    )
    check(matrix_kappa == matrix13 and pivots_kappa == pivots13,
          "p*y6 uses row-thirteen tangent quotient")
    check(len(conditions_kappa) == 1,
          "p*y6 leaves one bracket condition before response")
    kappa_equation = primitive(conditions_kappa[0])
    kappa_response = exact(sp.diff(kappa_equation, kappa20))
    check(kappa_response != 0 and not kappa_response.free_symbols,
          "p*y6 is a constant bracket pivot")
    kappa_graph = exact(sp.solve(kappa_equation, kappa20)[0])
    check(not sp.fraction(sp.together(kappa_graph))[1].free_symbols,
          "p*y6 bracket graph has constant denominator")
    a12_selected = R9.impose(a12, select_kappa, {kappa20: kappa_graph})
    c12_selected = R9.impose(c12, select_kappa, {kappa20: kappa_graph})
    a13_trial, c13_trial, theta13 = R9.append_tangent(
        a12_selected, c12_selected, 13, "row13_scout_theta13"
    )
    depth13, amatrix13, cmatrix13 = R9.projected_depth_residuals(
        a13_trial, c13_trial, 13
    )
    terminal13, depth_conditions13, depth_matrix13, depth_pivots13 = (
        R9.eliminate_linear_fibre(depth13, theta13)
    )
    depth_primitives13 = unique_primitives(depth_conditions13)
    check(
        (amatrix13.shape, amatrix13.rank(), cmatrix13.shape, cmatrix13.rank())
        == ((133, 308), 97, (147, 491), 117),
        "row-thirteen projected universes",
    )
    check(
        (depth_matrix13.shape, depth_matrix13.rank(), depth_pivots13)
        == ((66, 14), 5, (9, 10, 11, 12, 13)),
        "row-thirteen projected depth selector",
    )
    check(not depth_primitives13,
          "row-thirteen projected depth has no inherited sidecar")
    bracket_after_graph = [
        exact(value.subs(select_kappa).subs(kappa20, kappa_graph))
        for value in [
            R9.xcoeff(kappa_source - R8.predicted_G(13, a12, c12), degree)
            for degree in range(14)
        ]
    ]
    depth_after_terminal = [exact(value.subs(terminal13)) for value in depth13]
    check(all(value == 0 for value in bracket_after_graph),
          "row-thirteen bracket graph is literally sufficient")
    check(all(value == 0 for value in depth_after_terminal),
          "row-thirteen depth terminal is literally sufficient")

    off_graph_residual = exact(kappa_equation.subs(kappa20, kappa_graph + 1))
    check(off_graph_residual == kappa_response,
          "off-kappa hostile preserves the bracket obstruction")

    controls = {
        "dense": {R8.Phi: 1, R8.eta: 2, R8.alpha11: 3, R9.c51: 5},
        "Phi_eta_zero": {
            R8.Phi: 0, R8.eta: 0, R8.alpha11: 1, R9.c51: 2
        },
    }
    control_values: dict[str, sp.Expr] = {}
    for label, control in controls.items():
        value = exact(kappa_graph.subs(control))
        check(all(exact(residual.subs(control)) == 0
                  for residual in bracket_after_graph + depth_after_terminal),
              f"{label} literal row-thirteen control")
        control_values[label] = value

    print("row13_restricted_weight18_continuation_scout")
    print("inherited_source_projection=A4(Phi,eta,alpha11,c51)")
    print("inherited_solved=(c23,c70,lambda18,nu18);denominators=constant")
    print(
        f"row13_bracket_matrix={matrix13.shape};rank={matrix13.rank()};"
        f"pivots={pivots13}"
    )
    print(f"row13_raw_cokernel_dimension={len(raw_cokernel13)}")
    print(f"row13_selected_nonzero={sum(value != 0 for value in reduced13)}")
    print(f"row13_primitive_condition_count={len(primitives13)}")
    for index, value in enumerate(primitives13):
        variables = tuple(sorted(value.free_symbols, key=str))
        poly = sp.Poly(value, *variables, domain=sp.QQ) if variables else None
        print(
            f"row13_condition_{index}_terms="
            f"{len(poly.terms()) if poly is not None else 1};"
            f"degree={poly.total_degree() if poly is not None else 0};"
            f"symbols={tuple(map(str, variables))}"
        )
        print(f"row13_condition_{index}={sp.sstr(value)}")
    print(f"valuation13_weights={weights13}")
    print(f"valuation13_full_response_rank={response_full.rank()}")
    print(f"valuation13_student_responses={sp.sstr(student_responses)}")
    print(f"valuation13_ratios_to_py6={sp.sstr(response_ratios)}")
    for b in range(7):
        column = response_full[:, b]
        print(
            f"valuation13_channel=b{b},monomial=p^{13-2*b}*y^{b},"
            f"weight={weights13[b]},student={sp.sstr(student_responses[b])},"
            f"full_zero={column.is_zero_matrix}"
        )
    print(f"row13_student_base={sp.sstr(primitive(student_base))}")
    print(f"row13_conditions_with_channels={len(conditions_all)}")
    for index, value in enumerate(unique_primitives(reduced_all)):
        print(f"row13_channel_condition_{index}={sp.sstr(value)}")
    print(f"kappa20_response={sp.sstr(kappa_response)}")
    print(f"kappa20_graph={sp.sstr(kappa_graph)}")
    print(
        "row13_depth_universe="
        f"P2{amatrix13.shape}/rank{amatrix13.rank()},"
        f"P3{cmatrix13.shape}/rank{cmatrix13.rank()}"
    )
    print(
        "row13_depth_fibre="
        f"shape{depth_matrix13.shape},rank{depth_matrix13.rank()},"
        f"pivots{depth_pivots13},raw_conditions{len(depth_conditions13)},"
        f"primitive_conditions{len(depth_primitives13)},"
        f"dimension{len(theta13)-depth_matrix13.rank()}"
    )
    for index, value in enumerate(depth_primitives13):
        variables = tuple(sorted(value.free_symbols, key=str))
        poly = sp.Poly(value, *variables, domain=sp.QQ) if variables else None
        print(
            f"row13_depth_condition_{index}_terms="
            f"{len(poly.terms()) if poly is not None else 1};"
            f"degree={poly.total_degree() if poly is not None else 0};"
            f"symbols={tuple(map(str, variables))}"
        )
        print(f"row13_depth_condition_{index}={sp.sstr(value)}")
    print("row13_depth_sidecars=none")
    print("source_projection_after_kappa20=global_A4(Phi,eta,alpha11,c51)")
    print("solved_source_coordinates=(c23,c70,lambda18,nu18,kappa20)")
    print("row13_terminal_fibre=A9")
    for label, value in control_values.items():
        print(f"{label}_control=kappa20:{sp.sstr(value)};PASS")
    print(f"off_kappa_graph_control=residual:{sp.sstr(off_graph_residual)};PASS")
    print("scope=restricted_H14_plus_two_weight18_channels_plus_kappa20_py6_through_row13;not_complete_weight20;not_JC2")
    print(f"checks={CHECKS};exploratory_result=PASS")


if __name__ == "__main__":
    main()
