#!/usr/bin/env python3
"""Audit the THM-4426 companion for a bounded row-15 result.

The dependency is hash-pinned and extended only by constructing G_15. This
artifact is FINITE-EXACT evidence; it does not promote a theorem by itself.
"""

from __future__ import annotations

import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CANON = ROOT / "04-computation" / (
    "jc2_source_normal_row14_weight18_two_memory_"
    "rational_conic_companion_thm4426.py"
)
EXPECTED = "c647bdf2c894d60413b041197cff446884191bf1469372e6b4a6a55ae96fd52a"


def main() -> None:
    raw = CANON.read_bytes()
    actual = hashlib.sha256(raw).hexdigest()
    if actual != EXPECTED:
        raise AssertionError(f"canonical dependency hash changed: {actual}")
    text = raw.decode("utf-8")
    old_range = "grows = {row: tcoeff(source, row) for row in range(4, 15)}"
    new_range = "grows = {row: tcoeff(source, row) for row in range(4, 16)}"
    if text.count(old_range) != 1:
        raise AssertionError("unique row range anchor")
    text = text.replace(old_range, new_range)

    anchor = "    graph_order = (\n"
    if text.count(anchor) != 1:
        raise AssertionError("unique row-fifteen injection anchor")
    injected = r'''
    # Bounded row-fifteen relative-response scout.  Work first over the
    # unspecialized rational boundary conic Q(z,h18)=0, retaining the ten
    # terminal row-fourteen tangent coordinates until the row-fifteen bracket.
    free14 = tuple(theta14[j] for j in range(15) if j not in pivots_depth14)
    a14_boundary = apply(a14_trial, terminal14)
    c14_boundary = apply(c14_trial, terminal14)
    inherited15 = {
        **bracket_subs,
        **gate,
        **resolved13,
        **boundary,
        c51: fixed_c51,
        rho22: rho_parameter.subs(parameter, parameter),
        h18: h18,
        z: z,
    }
    # rho_parameter is a parametrized expression; use the unspecialized graph
    # here instead, then specialize to the G_m below.
    inherited15[rho22] = rho_graph
    g15_boundary = resolve(grows[15], inherited15)
    eq15_boundary = bracket_equations(
        a14_boundary, c14_boundary, 15, g15_boundary
    )
    select15, conditions15, raw15, matrix15, pivots15, row_pivots15 = (
        eliminate_linear(eq15_boundary, free14)
    )
    print("row15_relative_response_scout")
    print(
        f"unpaid_rank={matrix15.shape}_rank{matrix15.rank()}_"
        f"pivots{pivots15}_conditions{len(conditions15)}"
    )
    print("unpaid_condition_hashes=" + sp.sstr(tuple(expression_hash(v) for v in conditions15)))
    print("unpaid_condition_factors=" + sp.sstr(tuple(sp.factor(v) for v in conditions15)))

    response15 = tuple(sp.symbols("row15_r0:8"))
    exact_row15_packet = sum(
        response15[b] * R8.p**(15 - 2*b) * R8.y**b for b in range(8)
    )
    response_leading = tcoeff(exact_row15_packet, 15)
    check(
        response_leading == sum(response15[b] * x**b for b in range(8)),
        "row-fifteen exact-valuation packet",
    )
    eq15_response = bracket_equations(
        a14_boundary, c14_boundary, 15, g15_boundary + response_leading
    )
    select15_response, conditions15_response, raw15_response, matrix15_response, pivots15_response, row_pivots15_response = eliminate_linear(
        eq15_response, free14
    )
    response_matrix, response_rhs = sp.linear_eq_to_matrix(
        conditions15_response, response15
    )
    response_pivots = tuple(response_matrix.rref()[1])
    response_row_pivots = tuple(response_matrix.T.rref()[1])
    if response_pivots:
        response_minor = sp.factor(
            response_matrix.extract(response_row_pivots, response_pivots).det()
        )
    else:
        response_minor = sp.Integer(1)
    all_matrix, all_rhs = sp.linear_eq_to_matrix(
        eq15_response, free14 + response15
    )
    all_pivots = tuple(all_matrix.rref()[1])
    print(
        f"deduplicated_response_after_tangent={response_matrix.shape}_rank{response_matrix.rank()}_"
        f"pivots{response_pivots}_rowpivots{response_row_pivots}"
    )
    print("deduplicated_response_pivot_minor=" + sp.sstr(response_minor))
    print(
        f"joint_bracket={all_matrix.shape}_rank{all_matrix.rank()}_"
        f"pivots{all_pivots}_free{all_matrix.cols-all_matrix.rank()}"
    )
    response_graph, response_debts, response_raw, _, _, _ = eliminate_linear(
        conditions15_response, response15
    )
    print("deduplicated_response_debts_not_used=" + sp.sstr(tuple(response_debts)))
    print(
        "deduplicated_response_graph_hashes_not_used="
        + sp.sstr(tuple((str(key), expression_hash(value)) for key, value in response_graph.items()))
    )
    print(
        "deduplicated_response_graph_factors_not_used="
        + sp.sstr(tuple((str(key), sp.factor(value)) for key, value in response_graph.items()))
    )

    # Do not use the proportionality-deduplicated conditions to construct the
    # response graph: specialization can split formerly proportional raw
    # cokernel rows.  Instead use all six exact left-nullspace compatibilities.
    tangent_matrix, tangent_rhs = sp.linear_eq_to_matrix(
        eq15_response, free14
    )
    left15 = tuple(tangent_matrix.T.nullspace())
    compat15 = tuple(exact((ell.T * tangent_rhs)[0]) for ell in left15)
    raw_response_matrix, raw_response_rhs = sp.linear_eq_to_matrix(
        compat15, response15
    )
    raw_response_pivots = tuple(raw_response_matrix.rref()[1])
    raw_response_row_pivots = tuple(raw_response_matrix.T.rref()[1])
    check(
        len(left15) == 6 and raw_response_matrix.rank() == 1
        and raw_response_pivots == (0,),
        "raw row-fifteen relative response rank one",
    )
    response_row = raw_response_row_pivots[0]
    r0_graph = exact(sp.solve(compat15[response_row], response15[0])[0])
    relative_free = response15[1:]
    raw_compat_after = tuple(
        exact(value.subs(response15[0], r0_graph)) for value in compat15
    )
    raw_compat_on_gm = tuple(
        exact(value.subs({h18: h_parameter, z: z_parameter}))
        for value in raw_compat_after
    )
    check(
        all(value == 0 for value in raw_compat_on_gm),
        "all raw compatibilities vanish on Gm after response",
    )
    print(
        f"raw_cokernel={len(left15)};raw_response_rank={raw_response_matrix.rank()};"
        f"pivots={raw_response_pivots};relative_free={len(relative_free)}"
    )
    print("correct_r0_graph=" + sp.sstr(sp.factor(r0_graph)))
    print(
        "correct_r0_relative_coefficients="
        + sp.sstr(tuple(exact(sp.diff(r0_graph, value)) for value in relative_free))
    )
    response_constant = sp.Integer(
        10852621164972710686787843667734315747451565056000000000000000
    )
    inherited_quartic = conditions15[0]
    expected_r0_graph = exact(
        -inherited_quartic/response_constant
        - sp.Rational(6, 29)*response15[2]
        - sp.Rational(4, 29)*response15[4]
        - sp.Rational(24, 145)*response15[6]
    )
    check(r0_graph == expected_r0_graph, "exact row-fifteen quotient observable")
    r6_graph = exact(-sp.Rational(145, 24)*inherited_quartic/response_constant)
    check(constant_denominator(r6_graph), "weight-twenty-four response graph")
    r6_only_relative = {value: sp.Integer(0) for value in relative_free}
    r6_only_relative[response15[6]] = r6_graph
    check(
        exact(r0_graph.subs(r6_only_relative)) == 0,
        "weight-twenty-four response pays the bracket class",
    )
    print("quotient_observable=145*r0+30*r2+20*r4+24*r6")
    print("least_weight_visible_response=p^3*y^6;weight=24;coefficient=" + sp.sstr(r6_graph))
    print("weight23_hostile=p*y^7_is_bracket_invisible")
    print(
        "raw_compat_after_r0_hashes="
        + sp.sstr(tuple(expression_hash(value) for value in raw_compat_after))
    )

    response_section = {response15[0]: r0_graph}
    bracket_after_response = tuple(
        resolve(value, {**select15_response, **response_section})
        for value in eq15_response
    )
    bracket_after_response_gm = tuple(
        exact(value.subs({h18: h_parameter, z: z_parameter}))
        for value in bracket_after_response
    )
    check(
        all(value == 0 for value in bracket_after_response_gm),
        "sixteen row-fifteen bracket coefficients vanish on Gm",
    )
    tangent_minor15 = sp.factor(
        tangent_matrix.extract(row_pivots15_response, pivots15_response).det()
    )
    print("row15_tangent_pivot_minor=" + sp.sstr(tangent_minor15))
    print("row15_coeffwise_on_Gm=16_of_16_zero")

    inherited_pullback = exact(
        inherited_quartic.subs({h18: h_parameter, z: z_parameter})
    )
    pull_numerator, pull_denominator = sp.fraction(inherited_pullback)
    primitive_pullback = sp.Poly(
        pull_numerator, parameter, domain=sp.QQ
    ).primitive()[1].as_expr()
    primitive_pullback = sp.expand(primitive_pullback)
    if sp.LC(sp.Poly(primitive_pullback, parameter)) < 0:
        primitive_pullback = -primitive_pullback
    pull_factor_list = sp.factor_list(primitive_pullback)
    check(
        sp.degree(primitive_pullback, parameter) == 4
        and len(pull_factor_list[1]) == 1
        and pull_factor_list[1][0][0] == primitive_pullback
        and pull_factor_list[1][0][1] == 1,
        "inherited Gm hostile quartic is QQ-irreducible",
    )
    check(
        sp.Poly(pull_denominator, parameter).monoms() == [(2,)],
        "inherited Gm hostile has only the excluded s-zero pole",
    )
    print("inherited_Gm_quartic_primitive=" + sp.sstr(primitive_pullback))
    print("inherited_Gm_quartic_factor_list=" + sp.sstr(pull_factor_list))
    print(
        f"inherited_Gm_quartic_hash={expression_hash(primitive_pullback)};"
        f"denominator={sp.sstr(pull_denominator)}"
    )
    print("rational_hostile=without_weight24_response_every_s_in_Qstar_fails_row15_bracket")

    # Evaluate the exact generic data on three rational G_m controls.  The
    # symbolic minor is also pulled back to s to detect poles or rank drops.
    response_minor_s = exact(response_minor.subs({h18: h_parameter, z: z_parameter}))
    print("deduplicated_response_minor_on_Gm=" + sp.sstr(sp.factor(response_minor_s)))
    for control_s in (sp.Integer(1), sp.Integer(2), -sp.Integer(1)):
        control = {parameter: control_s}
        specialized_conditions = tuple(
            exact(value.subs({h18: h_parameter, z: z_parameter}).subs(control))
            for value in conditions15
        )
        print(
            f"Gm_control_s={control_s};unpaid_nonzero="
            f"{sum(value != 0 for value in specialized_conditions)};"
            f"hashes={sp.sstr(tuple(expression_hash(value) for value in specialized_conditions))}"
        )

    # Project one step deeper, keeping the seven bracket-neutral source
    # responses.  First take exact ranks at three rational controls; this is a
    # bounded scout, not yet a global rank assertion over Q(s).
    a14_response = apply(
        a14_boundary, select15_response, {response15[0]: r0_graph}
    )
    c14_response = apply(
        c14_boundary, select15_response, {response15[0]: r0_graph}
    )
    a15_trial, c15_trial, theta15 = append_tangent(
        a14_response, c14_response, 15, "gm_theta15"
    )
    depth15, ap15, cp15 = depth_equations(a15_trial, c15_trial, 15)
    terminal15, depth_conditions15, raw_depth15, matrix_depth15, pivots_depth15, row_pivots_depth15 = eliminate_linear(
        depth15, theta15
    )
    print(
        f"row15_depth_universes=A{ap15.shape}_rank{ap15.rank()};"
        f"C{cp15.shape}_rank{cp15.rank()}"
    )
    print(
        f"row15_depth_selector={matrix_depth15.shape}_rank{matrix_depth15.rank()}_"
        f"pivots{pivots_depth15}_conditions{len(depth_conditions15)}"
    )
    print(
        "row15_depth_condition_hashes="
        + sp.sstr(tuple(expression_hash(value) for value in depth_conditions15))
    )
    check(len(depth_conditions15) == 1, "one row-fifteen depth representative")
    depth15_condition = depth_conditions15[0]
    print("row15_depth_condition_factor=" + sp.sstr(sp.factor(depth15_condition)))
    print(
        "row15_depth_condition_symbols="
        + sp.sstr(tuple(sorted(map(str, depth15_condition.free_symbols))))
    )
    depth15_pivot_minor = sp.factor(
        matrix_depth15.extract(row_pivots_depth15, pivots_depth15).det()
    )
    print("row15_depth_pivot_minor=" + sp.sstr(depth15_pivot_minor))
    print(f"row15_terminal_fibre=A{len(theta15)-len(pivots_depth15)}")

    # This is the decisive exact generic test: pull every raw depth residual
    # to Q(s)[r1,...,r7], not merely to finitely many rational controls.
    raw_depth15_gm = tuple(
        exact(value.subs({h18: h_parameter, z: z_parameter}))
        for value in raw_depth15
    )
    check(
        all(value == 0 for value in raw_depth15_gm),
        "all row-fifteen depth residuals vanish over rational Gm",
    )
    print("row15_depth_coeffwise_on_Gm=91_of_91_zero_for_all_seven_relative_parameters")
    r6_only_bracket = tuple(
        exact(value.subs(r6_only_relative)) for value in bracket_after_response_gm
    )
    r6_only_depth = tuple(
        exact(value.subs(r6_only_relative)) for value in raw_depth15_gm
    )
    check(
        all(value == 0 for value in r6_only_bracket)
        and all(value == 0 for value in r6_only_depth),
        "least-weight response kills all row-fifteen bracket and depth coefficients",
    )
    print("weight24_section_verification=16_of_16_bracket_and_91_of_91_depth_zero")
    for control_s in (sp.Integer(1), sp.Integer(2), -sp.Integer(1)):
        control_rows = tuple(
            exact(value.subs({h18: h_parameter, z: z_parameter}).subs(parameter, control_s))
            for value in raw_depth15
        )
        control_matrix, control_rhs = sp.linear_eq_to_matrix(
            control_rows, relative_free
        )
        control_augmented = control_matrix.row_join(control_rhs)
        print(
            f"row15_depth_control_s={control_s};shape={control_matrix.shape};"
            f"response_rank={control_matrix.rank()};augmented_rank={control_augmented.rank()};"
            f"nonzero_rows={sum(value != 0 for value in control_rows)}"
        )
    print("row15_field=every_characteristic_zero_field_for_the_parametric_lift;QQ_irreducibility_only_for_the_rational_hostile")
    print("row15_scope=finite_boundary_Gm_row15_bracket_and_projected_depth_with_exact_valuation15_packet;not_complete_weights15_to22_or_full_source_or_termination_or_JC2")

'''
    text = text.replace(anchor, injected + anchor)
    namespace = {
        "__name__": "__main__",
        "__file__": str(CANON),
        "__package__": None,
    }
    exec(compile(text, str(CANON), "exec"), namespace)


if __name__ == "__main__":
    main()
