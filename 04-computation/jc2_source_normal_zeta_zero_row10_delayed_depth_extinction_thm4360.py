#!/usr/bin/env python3
"""Primary exact certificate for THM-4360.

Restrict THM-4308's finite source-normal response to the complete
zeta_3=Z=0 source plane from THM-4357, retaining beta_11 symbolically.
The certificate reconstructs the general THM-4315 row-nine bracket
equation, independently reruns the row-nine projected P_2/P_3 test on
this non-corner slice, and then applies the row-ten Student compatibility.
The scalar row-ten gate has two apparent one-parameter components, but the
full G_10 equation selects a row-nine tangent which violates an already-
required P_2 functional by a fixed nonzero rational number.

This is a finite projected-depth and necessary row-ten extinction result,
uniformly for arbitrary beta_11 in the residual-weight-at-most-twelve
universe.  It is not an all-row B_2 lift, termination, seam-entry, JC(2),
or DC(2) theorem.
"""

from __future__ import annotations

from math import gcd
from pathlib import Path
import sys

import sympy as sp


# Python's Windows text layer otherwise turns each printed LF into CRLF.
# Pin stdout so normal, optimized, and frozen raw byte streams really match.
if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import jc2_source_normal_student_stein_row9_thm4315 as R9  # noqa: E402


R8 = R9.R8
x = R8.x
t = R8.t
eta = R8.eta
alpha = R8.alpha11
beta = R8.beta11

XI_Z = sp.Rational(1563264759296, 115860375)
U_S = -sp.Rational(5200877686784, 1042743375)
W_S = sp.Rational(2539122786304, 115860375)
K_S = -sp.Rational(32, 5)
LAMBDA_S = sp.Rational(17651227389952, 1042743375)

F = sp.Poly(
    2393096045625 * eta**2 - 415832184456871936,
    eta,
    domain=sp.QQ,
)
C_ETA = sp.Rational(415832184456871936, 2393096045625)
ALPHA_ETA = -sp.Rational(619241095293435904, 1010418330375)
ALPHA_OVER_ETA = sp.cancel(ALPHA_ETA / C_ETA)

CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def exact_zero(expr: sp.Expr) -> bool:
    return sp.cancel(sp.together(sp.expand(expr))) == 0


def check_zero(expr: sp.Expr, label: str) -> None:
    check(exact_zero(expr), label)


def reduce_f(expr: sp.Expr) -> sp.Expr:
    """Reduce coefficientwise in QQ(other symbols)[eta]/(F)."""

    value = sp.cancel(sp.together(sp.expand(expr)))
    numerator, denominator = sp.fraction(value)
    domain = sp.EX
    num_poly = sp.Poly(numerator, eta, domain=domain)
    den_poly = sp.Poly(denominator, eta, domain=domain)
    f_poly = sp.Poly(F.as_expr(), eta, domain=domain)
    den_rem = den_poly.rem(f_poly)
    check(
        den_rem.degree() == 0 and den_rem.as_expr() != 0,
        "F-reduction denominator unit",
    )
    return sp.cancel(num_poly.rem(f_poly).as_expr() / den_rem.as_expr())


def primitive_integer_functional(vector: sp.Matrix) -> list[int]:
    denominator = 1
    for value in vector:
        denominator = sp.ilcm(denominator, int(sp.denom(value)))
    entries = [int(value * denominator) for value in vector]
    common = 0
    for value in entries:
        common = gcd(common, abs(value))
    entries = [value // common for value in entries]
    first = next(value for value in entries if value)
    return entries if first > 0 else [-value for value in entries]


def main() -> None:
    # Reconstruct THM-4308's exact row-eight bracket/depth family.
    arows, crows, bracket_subs = R8.build_bracket_rows()
    amatrix8, cmatrix8, theta8_symbols, terminal_subs = R8.hasse_audit(
        arows, crows, bracket_subs
    )
    check(
        (amatrix8.rows, amatrix8.cols, amatrix8.rank()) == (63, 131, 51),
        "row8 P2",
    )
    check(
        (cmatrix8.rows, cmatrix8.cols, cmatrix8.rank()) == (72, 204, 63),
        "row8 P3",
    )

    s_subs = {R8.Phi: 0, R8.xi10: XI_Z}
    response = {
        symbol: sp.cancel(value.subs(terminal_subs).subs(s_subs))
        for symbol, value in bracket_subs.items()
    }
    check_zero(response[R8.U] - U_S, "zeta-zero U")
    check_zero(response[R8.W] - W_S, "zeta-zero W")
    check_zero(response[R8.Z], "zeta-zero Z")
    check_zero(R8.K.subs(terminal_subs) - K_S, "zeta-zero K")
    check_zero(response[R8.U] + response[R8.W] - LAMBDA_S, "zeta-zero Lambda")

    # THM-4315's named E9 gate becomes a squarefree quadratic on the
    # whole zeta_3=Z=0 plane.  Distinguish that normalization from the
    # literal primitive Student evaluation.
    e9_s = sp.factor(R9.E9_GATE.subs(s_subs))
    expected_e9_s = -sp.Rational(11, 102987) * F.as_expr()
    check_zero(e9_s - expected_e9_s, "restricted E9")
    difference9_full = sp.expand(
        R8.tcoeff(R8.H, 9).subs(bracket_subs)
        - R8.predicted_G(9, arows, crows)
    )
    row9_student = R9.primitive_student_row(9)
    e9_literal = sp.factor(
        sum(
            row9_student[j] * R8.xcoeff(difference9_full, j)
            for j in range(10)
        ).subs(terminal_subs).subs(s_subs)
    )
    expected_e9_literal = sp.Rational(143, 56308142250) * F.as_expr()
    check_zero(e9_literal - expected_e9_literal, "literal Student E9")
    check_zero(
        56308142250 * e9_literal - 143 * F.as_expr(),
        "literal E9 primitive normalization",
    )
    check(
        sp.gcd(F, F.diff()) == sp.Poly(1, eta, domain=sp.QQ),
        "F squarefree",
    )
    check(C_ETA > 0, "two real eta roots")
    check(
        not sp.ntheory.primetest.is_square(int(sp.numer(C_ETA))),
        "eta numerator nonsquare",
    )
    check(
        not sp.ntheory.primetest.is_square(int(sp.denom(C_ETA))),
        "eta denominator nonsquare",
    )
    check(F.degree() == 2 and F.LC() != 0 and F.TC() != 0, "F reduced two-root gate")
    check_zero(sp.diff(e9_s, beta), "named E9 beta-independent")

    # The literal G9 equation selects all seven old terminal directions.
    g9 = sp.expand(R8.tcoeff(R8.H, 9))
    expected_g9 = (
        (20 * R8.U + 10 * R8.W + 4 * R8.Z) * x**6
        + (10 * R8.alpha11 + 6 * R8.beta11) * x**7
        + (5 * R8.upsilon5 + 4 * R8.xi10) * x**8
        + (R8.eta + R8.zeta3) * x**9
    )
    check_zero(g9 - expected_g9, "literal G9")
    check_zero(sp.diff(g9, beta) - 6 * x**7, "literal G9 beta channel")
    fixed_a, fixed_c, old_tangent, _ = R9.solve_row_eight_for_g9(
        arows,
        crows,
        bracket_subs,
        terminal_subs,
        theta8_symbols,
        g9,
    )
    check(old_tangent.rank() == 7, "old terminal rank seven")
    check(old_tangent.cols == 7, "old terminal unique selection")

    fixed_a = [
        sp.expand(row.subs(bracket_subs).subs(terminal_subs).subs(s_subs))
        for row in fixed_a
    ]
    fixed_c = [
        sp.expand(row.subs(bracket_subs).subs(terminal_subs).subs(s_subs))
        for row in fixed_c
    ]
    g9_s = sp.expand(
        g9.subs(bracket_subs).subs(terminal_subs).subs(s_subs)
    )
    bracket_defect = sp.expand(R8.predicted_G(9, fixed_a, fixed_c) - g9_s)
    check(reduce_f(bracket_defect) == 0, "G9 bracket in quotient")

    # Add row nine with all ten fresh determinant-tangent coordinates.
    b9 = R8.B_row(9, fixed_a, fixed_c)
    a9base, c9base = R8.particular_row(9, b9)
    theta9_symbols = list(sp.symbols("theta9_0:10"))
    theta9 = sum(theta9_symbols[j] * x**j for j in range(10))
    a9 = sp.expand(a9base + theta9 * sp.diff(R8.A0, x))
    c9 = sp.expand(c9base + theta9 * sp.diff(R8.C0, x))
    row9_a = fixed_a + [a9]
    row9_c = fixed_c + [c9]

    # Rebuild both projected-depth universes on the whole plane; THM-4315's
    # row-nine depth conclusion itself was scoped only to the cubic corner.
    acoords, amatrix = R9.depth_matrix(2, 9)
    ccoords, cmatrix = R9.depth_matrix(3, 9)
    check(
        (len(acoords), amatrix.cols, amatrix.rank()) == (75, 160, 59),
        "row9 P2 universe",
    )
    check(
        (len(ccoords), cmatrix.cols, cmatrix.rank()) == (85, 251, 73),
        "row9 P3 universe",
    )
    anull = amatrix.T.nullspace()
    cnull = cmatrix.T.nullspace()
    check((len(anull), len(cnull)) == (16, 12), "row9 left nullities")
    avec = sp.Matrix(
        [R8.xcoeff(row9_a[n], r) for n, r in acoords]
    )
    cvec = sp.Matrix(
        [R8.xcoeff(row9_c[n], r) for n, r in ccoords]
    )
    residuals = [sp.expand((row.T * avec)[0]) for row in anull]
    residuals += [sp.expand((row.T * cvec)[0]) for row in cnull]
    depth_tangent, _ = sp.linear_eq_to_matrix(residuals, theta9_symbols)
    check(depth_tangent.rank() == 3, "row9 terminal tangent rank")
    check(depth_tangent.rref()[1] == (7, 8, 9), "row9 terminal pivots")

    zero_free = {symbol: 0 for symbol in theta9_symbols[:7]}
    reduced_residuals = [
        sp.expand(value.subs(zero_free)) for value in residuals
    ]
    pivot_solutions = sp.solve(
        reduced_residuals,
        theta9_symbols[7:],
        dict=True,
    )
    check(len(pivot_solutions) == 1, "row9 pivot solution")
    terminal9 = {**zero_free, **pivot_solutions[0]}
    expected_terminal9 = {
        theta9_symbols[7]: -sp.Rational(317075357581312, 347581125),
        theta9_symbols[8]: -(125 * alpha + 100 * beta + 198 * eta) / 15,
        theta9_symbols[9]: -sp.Rational(553237643264, 23172075),
    }
    for symbol, expected in expected_terminal9.items():
        check(
            reduce_f(terminal9[symbol] - expected) == 0,
            f"row9 particular {symbol}",
        )
    check(
        reduce_f(
            sp.diff(terminal9[theta9_symbols[8]], beta)
            + sp.Rational(20, 3)
        )
        == 0,
        "row9 beta address moves",
    )
    depth_obstructions = [
        reduce_f(value.subs(terminal9)) for value in residuals
    ]
    nonzero_depth = [
        (i, value)
        for i, value in enumerate(depth_obstructions)
        if value != 0
    ]
    print(f"DEPTH_NONZERO_COUNT={len(nonzero_depth)}")
    for index, value in nonzero_depth:
        print(f"DEPTH_NONZERO[{index}]={sp.sstr(value)}")
    check(
        not nonzero_depth,
        "no extra row9 projected-depth equation modulo F",
    )
    check(
        10 - depth_tangent.rank() == 7,
        "new terminal affine dimension seven",
    )

    # The row-ten Student scalar first leaves two apparent source survivors.
    g10 = sp.expand(R8.tcoeff(R8.H, 10))
    expected_g10 = (
        (15 * R8.U + 10 * R8.W + 6 * R8.Z) * x**8
        + (5 * R8.alpha11 + 4 * R8.beta11) * x**9
        + (R8.upsilon5 + R8.xi10) * x**10
    )
    check_zero(g10 - expected_g10, "literal G10")
    check_zero(sp.diff(g10, beta) - 4 * x**9, "literal G10 beta channel")
    row10 = R9.primitive_student_row(10)
    check(
        row10
        == [46189, 0, 14586, 0, 15444, 0, 30888, 0, 99792, 0, 489888],
        "row10 Student row",
    )
    tangent_a = fixed_a + [a9]
    tangent_c = fixed_c + [c9]
    base_a = fixed_a + [a9base]
    base_c = fixed_c + [c9base]
    tangent_change10 = sp.expand(
        R8.predicted_G(10, tangent_a, tangent_c)
        - R8.predicted_G(10, base_a, base_c)
    )
    tangent_scalar10 = sum(
        row10[j] * R8.xcoeff(tangent_change10, j) for j in range(11)
    )
    check_zero(
        tangent_scalar10,
        "row10 obstruction independent of theta9",
    )
    g10_s = sp.expand(
        g10.subs(bracket_subs).subs(terminal_subs).subs(s_subs)
    )
    difference10 = sp.expand(
        g10_s - R8.predicted_G(10, base_a, base_c)
    )
    obstruction10 = sp.factor(
        sum(
            row10[j] * R8.xcoeff(difference10, j)
            for j in range(11)
        )
    )
    obstruction10_mod_f = sp.factor(reduce_f(obstruction10))
    print(f"ROW10_OBSTRUCTION_RAW={sp.sstr(obstruction10)}")
    print(f"ROW10_OBSTRUCTION_MOD_F={sp.sstr(obstruction10_mod_f)}")
    expected_obstruction10_mod_f = sp.Rational(143, 3128230125) * (
        1010418330375 * alpha * eta + 619241095293435904
    )
    expected_obstruction10 = sp.Rational(143, 84462213375) * (
        27281294920125 * alpha * eta
        + 241702700608125 * eta**2
        - 25279541057221296128
    )
    check_zero(obstruction10 - expected_obstruction10, "row10 raw obstruction")
    check_zero(
        obstruction10_mod_f - expected_obstruction10_mod_f,
        "row10 reduced obstruction",
    )
    check_zero(sp.diff(obstruction10, beta), "row10 obstruction beta cancellation")
    alpha_to_eta = {alpha: ALPHA_OVER_ETA * eta}
    check(
        reduce_f((alpha * eta - ALPHA_ETA).subs(alpha_to_eta)) == 0,
        "row10 source parametrization",
    )

    # The full G10 equation selects all ten theta9 coordinates.  Test the
    # selected predecessor against every already-required row-nine depth row.
    difference10_tangent = sp.expand(
        g10_s - R8.predicted_G(10, tangent_a, tangent_c)
    )
    equations10_general = [
        R8.xcoeff(difference10_tangent, j) for j in range(11)
    ]
    matrix10_general, rhs10_general = sp.linear_eq_to_matrix(
        equations10_general,
        theta9_symbols,
    )
    check(matrix10_general.rank() == 10, "general G10 tangent rank ten")
    check(
        matrix10_general.row_join(rhs10_general).rank() == 11,
        "generic G10 augmented rank eleven",
    )
    row_pivots10_general = matrix10_general.T.rref()[1]
    selected10_general = (
        matrix10_general[list(row_pivots10_general), :].inv()
        * rhs10_general[list(row_pivots10_general), :]
    )
    selected10_general_subs = {
        theta9_symbols[j]: sp.cancel(selected10_general[j])
        for j in range(10)
    }
    beta_translation = tuple(
        sp.factor(sp.diff(selected10_general_subs[symbol], beta))
        for symbol in theta9_symbols
    )
    expected_beta_translation = (
        sp.Rational(2384, 945),
        0,
        -sp.Rational(504928, 8505),
        0,
        sp.Rational(9484, 105),
        0,
        -sp.Rational(568, 7),
        0,
        -sp.Rational(20, 3),
        0,
    )
    check(beta_translation == expected_beta_translation, "selected beta translation")
    selected_depth_general = [
        sp.factor(value.subs(selected10_general_subs)) for value in residuals
    ]
    check(
        all(exact_zero(sp.diff(value, beta)) for value in selected_depth_general),
        "all selected prior-depth rows beta-invariant",
    )
    depth12_after_g10 = sp.factor(
        residuals[12].subs(selected10_general_subs)
    )
    depth12_after_g10_mod_f = sp.factor(reduce_f(depth12_after_g10))
    print(
        f"ROW10_DEPTH12_AFTER_G10_RAW={sp.sstr(depth12_after_g10)}"
    )
    print(
        "ROW10_DEPTH12_AFTER_G10_MOD_F="
        f"{sp.sstr(depth12_after_g10_mod_f)}"
    )
    selected_depth_mod_f_general = [
        sp.factor(reduce_f(value)) for value in selected_depth_general
    ]
    nonzero_general = [
        (index, value)
        for index, value in enumerate(selected_depth_mod_f_general)
        if value != 0
    ]
    check([index for index, _ in nonzero_general] == [12, 14, 26], "general depth indices")
    e10_numerator = 1010418330375 * alpha * eta + 619241095293435904
    depth_numerator = 5052091651875 * alpha * eta + 3193963923683016704
    expected_general = {
        12: depth_numerator / sp.Integer(14189651847000),
        14: 13 * e10_numerator / sp.Integer(153248239947600),
        26: -13 * e10_numerator / sp.Integer(204330986596800),
    }
    for index, value in nonzero_general:
        check_zero(value - expected_general[index], f"general depth residual {index}")

    equations10 = [
        value.subs(alpha_to_eta) for value in equations10_general
    ]
    matrix10, rhs10 = sp.linear_eq_to_matrix(
        equations10,
        theta9_symbols,
    )
    check(matrix10.rank() == 10, "G10 tangent rank ten")
    row_pivots10 = matrix10.T.rref()[1]
    check(len(row_pivots10) == 10, "G10 pivot rows")
    selected10 = (
        matrix10[list(row_pivots10), :].inv()
        * rhs10[list(row_pivots10), :]
    )
    selected10_subs = {
        theta9_symbols[j]: sp.cancel(selected10[j]) for j in range(10)
    }
    g10_residuals = [
        reduce_f(value.subs(selected10_subs)) for value in equations10
    ]
    check(
        all(value == 0 for value in g10_residuals),
        "G10 full selection modulo F",
    )
    selected_depth = [
        sp.factor(
            reduce_f(value.subs(alpha_to_eta).subs(selected10_subs))
        )
        for value in residuals
    ]
    nonzero_selected_depth = [
        (i, value)
        for i, value in enumerate(selected_depth)
        if value != 0
    ]
    print(
        "ROW10_SELECTED_DEPTH_NONZERO_COUNT="
        f"{len(nonzero_selected_depth)}"
    )
    for index, value in nonzero_selected_depth:
        print(f"ROW10_SELECTED_DEPTH_NONZERO[{index}]={sp.sstr(value)}")
    check(
        len(nonzero_selected_depth) == 1,
        "unique selected depth obstruction",
    )
    obstructing_index, obstructing_value = nonzero_selected_depth[0]
    check(
        obstructing_index == 12,
        "unique obstruction is A functional 12",
    )
    primitive_obstruction_row = primitive_integer_functional(
        anull[obstructing_index]
    )
    primitive_terms = [
        (coordinate, coefficient)
        for coordinate, coefficient in zip(
            acoords,
            primitive_obstruction_row,
        )
        if coefficient
    ]
    check(
        sp.Matrix([primitive_obstruction_row]) * amatrix
        == sp.zeros(1, amatrix.cols),
        "displayed primitive functional annihilates pi9(P2)",
    )
    scale_to_primitive = next(
        sp.Rational(
            primitive_obstruction_row[j],
            anull[obstructing_index][j],
        )
        for j in range(len(primitive_obstruction_row))
        if primitive_obstruction_row[j]
    )
    primitive_obstruction_value = sp.factor(
        obstructing_value * scale_to_primitive
    )
    print(f"ROW10_DEPTH_OBSTRUCTION_FUNCTIONAL={primitive_terms}")
    print(
        "ROW10_DEPTH_OBSTRUCTION_PRIMITIVE_VALUE="
        f"{primitive_obstruction_value}"
    )
    check(
        primitive_obstruction_value != 0,
        "primitive row9 depth obstruction",
    )

    # Compact exact unit-ideal witness for the delayed-depth extinction.
    incompatibility = sp.expand(depth_numerator - 5 * e10_numerator)
    check(
        incompatibility == 97758447215837184,
        "E10/depth unit-ideal witness",
    )
    check_zero(
        sp.Rational(incompatibility, 14189651847000)
        - primitive_obstruction_value,
        "normalized E10/depth contradiction",
    )

    # Exhibit one actual projected row-nine point coefficientwise.
    row9_a[-1] = sp.expand(row9_a[-1].subs(terminal9))
    row9_c[-1] = sp.expand(row9_c[-1].subs(terminal9))
    aseries = sum(row9_a[n] * t**n for n in range(10))
    cseries = sum(row9_c[n] * t**n for n in range(10))
    source_s = sp.expand(
        R8.G.subs(bracket_subs).subs(terminal_subs).subs(s_subs)
    )
    jac = R8.jacobian(aseries, cseries)
    check(
        reduce_f(R8.tcoeff(jac, 0) - 1) == 0,
        "Jacobian row zero",
    )
    for n in range(1, 9):
        check(reduce_f(R8.tcoeff(jac, n)) == 0, f"Jacobian row {n}")
    for n in range(10):
        check(
            reduce_f(
                R8.tcoeff(R8.P(aseries, cseries) - source_s, n)
            )
            == 0,
            f"P/G row {n}",
        )

    print("THM4360_ZETA_ZERO_ROW10_DELAYED_DEPTH_EXTINCTION=PASS")
    print(
        "ZETA_ZERO_PLANE: Phi=0 xi=1563264759296/115860375; "
        "eta,alpha,beta free"
    )
    print(f"E9_NAMED={sp.sstr(e9_s)}")
    print(f"E9_LITERAL_STUDENT={sp.sstr(e9_literal)}")
    print(f"F_ETA={sp.sstr(F.as_expr())}")
    print(f"ETA_SQUARED={C_ETA}")
    print(
        "ROW9_GEOMETRY closure=two_A2(alpha,beta) "
        "D(beta)=two_(A1xGm); no Q-rational eta"
    )
    print(
        "ROW8_TERMINAL old_affine_dimension=7 G9_rank=7 "
        "selected_image_dimension=0"
    )
    print(
        "ROW9_DEPTH pi9(P2)=75x160/rank59 "
        "pi9(P3)=85x251/rank73"
    )
    print(
        "ROW9_TERMINAL tangent_rank=3 affine_dimension=7 "
        "extra_source_equations=0_mod_F"
    )
    print(
        "ROW9_DEPTH_PARTICULAR q7="
        f"{sp.sstr(reduce_f(terminal9[theta9_symbols[7]]))}; q8="
        f"{sp.sstr(reduce_f(terminal9[theta9_symbols[8]]))}; q9="
        f"{sp.sstr(reduce_f(terminal9[theta9_symbols[9]]))}"
    )
    print(
        f"ROW10_NECESSARY_MOD_F={sp.sstr(obstruction10_mod_f)}=0"
    )
    print(
        f"ROW10_SOURCE alpha*eta={ALPHA_ETA}; "
        f"alpha=({ALPHA_OVER_ETA})*eta modulo F"
    )
    print(
        "ROW10_SCALAR_GEOMETRY closure=two_A1(beta) "
        "D(beta)=two_Gm"
    )
    print(f"BETA_TRANSLATION={beta_translation}")
    print(f"ROW10_GENERAL_DEPTH_INDICES={[index for index, _ in nonzero_general]}")
    print(
        "ROW10_SELECTED_ROW9_DEPTH_NONZERO="
        f"{len(nonzero_selected_depth)}"
    )
    print(
        "ROW10_SELECTED_ROW9_DEPTH_PRIMITIVE_VALUE="
        f"{primitive_obstruction_value}"
    )
    print(
        "ROW10_UNIT_WITNESS=depth_numerator-5*E10_numerator="
        f"{incompatibility}"
    )
    print(f"CHECKS={CHECKS}")
    print(
        "SCOPE finite row9 bracket/projected P2/P3 and row10 necessary "
        "extinction on the zeta3=Z=0 source plane only; no infinite-depth "
        "membership, all-row B2 lift, termination, seam entry, JC2, or DC2"
    )


if __name__ == "__main__":
    main()
