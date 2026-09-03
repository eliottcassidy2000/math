#!/usr/bin/env python3
"""Independent exact audit for THM-4388.

This file does not import the THM-4388 primary certificate.  It rebuilds the
literal weight-13 source face from the audited THM-4308/4315 row operators,
uses an independently organized linear-fibre elimination, treats Phi=0 and
Phi!=0 separately, and checks both row-12 extinction mechanisms.
"""

from __future__ import annotations

from pathlib import Path
import sys

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import jc2_source_normal_bracket_hasse_rows8_thm4308 as R8  # noqa: E402
import jc2_source_normal_student_stein_row9_thm4315 as R9  # noqa: E402


x, t = R8.x, R8.t
Phi, eta = R8.Phi, R8.eta
c51, c23 = sp.symbols("audit_c51 audit_c23")
Y, ratio = sp.symbols("audit_Y audit_ratio")
CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def exact(value: sp.Expr) -> sp.Expr:
    return sp.cancel(sp.together(sp.expand(value)))


def xcoeff(value: sp.Expr, degree: int) -> sp.Expr:
    return sp.expand(value).coeff(x, degree)


def tcoeff(value: sp.Expr, degree: int) -> sp.Expr:
    return sp.expand(value).coeff(t, degree)


def primitive(value: sp.Expr) -> sp.Expr:
    numerator, _ = sp.fraction(exact(value))
    variables = tuple(sorted(numerator.free_symbols, key=str))
    if not variables:
        return sp.signsimp(numerator)
    return sp.Poly(numerator, *variables, domain=sp.QQ).primitive()[1].as_expr()


def proportional(left: sp.Expr, right: sp.Expr, *variables: sp.Symbol) -> bool:
    lp = sp.Poly(sp.fraction(exact(left))[0], *variables, domain=sp.QQ)
    rp = sp.Poly(sp.fraction(exact(right))[0], *variables, domain=sp.QQ)
    return lp.monic() == rp.monic()


def apply(rows: list[sp.Expr], *maps: dict[sp.Symbol, sp.Expr]) -> list[sp.Expr]:
    answer = rows
    for mapping in maps:
        answer = [exact(row.subs(mapping)) for row in answer]
    return answer


def append_tangent(
    arows: list[sp.Expr],
    crows: list[sp.Expr],
    row: int,
    prefix: str,
) -> tuple[list[sp.Expr], list[sp.Expr], tuple[sp.Symbol, ...]]:
    bvalue = R8.B_row(row, arows, crows)
    abase, cbase = R8.particular_row(row, bvalue)
    symbols = tuple(sp.symbols(f"{prefix}_0:{row + 1}"))
    theta = sum(symbols[j] * x**j for j in range(row + 1))
    return (
        arows + [sp.expand(abase + theta * sp.diff(R8.A0, x))],
        crows + [sp.expand(cbase + theta * sp.diff(R8.C0, x))],
        symbols,
    )


def bracket_equations(
    arows: list[sp.Expr],
    crows: list[sp.Expr],
    row: int,
    source: sp.Expr,
) -> list[sp.Expr]:
    difference = exact(source - R8.predicted_G(row, arows, crows))
    check(
        difference == 0 or sp.Poly(difference, x, domain=sp.EX).degree() <= row,
        f"row {row} bracket degree cap",
    )
    return [xcoeff(difference, degree) for degree in range(row + 1)]


def depth_equations(
    arows: list[sp.Expr],
    crows: list[sp.Expr],
    row: int,
) -> tuple[list[sp.Expr], sp.Matrix, sp.Matrix]:
    acoords, amatrix = R9.depth_matrix(2, row)
    ccoords, cmatrix = R9.depth_matrix(3, row)
    avec = sp.Matrix([xcoeff(arows[n], degree) for n, degree in acoords])
    cvec = sp.Matrix([xcoeff(crows[n], degree) for n, degree in ccoords])
    equations = [sp.expand((left.T * avec)[0]) for left in amatrix.T.nullspace()]
    equations += [sp.expand((left.T * cvec)[0]) for left in cmatrix.T.nullspace()]
    return equations, amatrix, cmatrix


def eliminate_linear(
    equations: list[sp.Expr],
    variables: tuple[sp.Symbol, ...],
) -> tuple[
    dict[sp.Symbol, sp.Expr],
    list[sp.Expr],
    list[sp.Expr],
    sp.Matrix,
    tuple[int, ...],
]:
    """Eliminate a lexicographic maximal-rank linear fibre.

    The returned raw residuals retain their literal scaling; ``conditions``
    are primitive and deduplicated up to sign.
    """

    matrix, rhs = sp.linear_eq_to_matrix(equations, variables)
    pivots = tuple(matrix.rref()[1])
    row_pivots = tuple(matrix.T.rref()[1])
    if pivots:
        square = matrix.extract(row_pivots, pivots)
        check(square.rows == square.cols and square.det() != 0, "linear pivot minor")
        values = square.inv() * rhs.extract(row_pivots, (0,))
        substitutions = {
            variables[column]: exact(values[index])
            for index, column in enumerate(pivots)
        }
    else:
        substitutions = {}
    raw = [exact(value.subs(substitutions)) for value in equations]
    conditions: list[sp.Expr] = []
    for value in raw:
        if value == 0:
            continue
        candidate = primitive(value)
        if all(
            exact(candidate - old) != 0 and exact(candidate + old) != 0
            for old in conditions
        ):
            conditions.append(candidate)
    return substitutions, conditions, raw, matrix, pivots


def reduce_mod(value: sp.Expr, modulus: sp.Poly) -> sp.Expr:
    """Reduce a rational expression modulo a squarefree univariate modulus."""

    numerator, denominator = sp.fraction(exact(value))
    variable = modulus.gens[0]
    coefficient_symbols = tuple(
        sorted((numerator.free_symbols | denominator.free_symbols) - {variable}, key=str)
    )
    domain = sp.QQ.frac_field(*coefficient_symbols) if coefficient_symbols else sp.QQ
    mod = sp.Poly(modulus.as_expr(), variable, domain=domain)
    num = sp.Poly(numerator, variable, domain=domain)
    den = sp.Poly(denominator, variable, domain=domain)
    check(sp.gcd(den, mod).degree() == 0, "localized denominator is a unit")
    return sp.factor((num * sp.invert(den, mod)).rem(mod).as_expr())


def quotient_polynomial(value: sp.Expr, parity: int) -> sp.Poly:
    polynomial = sp.Poly(value, Phi, eta, domain=sp.QQ)
    answer = sp.Integer(0)
    for (pdegree, edegree), coefficient in polynomial.terms():
        total = pdegree + edegree
        check(total % 2 == parity, "row-12 sign parity")
        answer += coefficient * ratio**edegree * Y ** ((total - parity) // 2)
    return sp.Poly(answer, Y, ratio, domain=sp.QQ).primitive()[1]


def mod_value(value: sp.Expr, substitutions: dict[sp.Symbol, int], prime: int) -> int:
    scalar = exact(value.subs(substitutions))
    check(not scalar.free_symbols, "finite-control scalar")
    numerator, denominator = sp.fraction(scalar)
    den = int(denominator) % prime
    check(den != 0, "finite-control denominator")
    return int(numerator) % prime * pow(den, -1, prime) % prime


def build_through_row_eight() -> tuple[
    list[sp.Expr],
    list[sp.Expr],
    dict[sp.Symbol, sp.Expr],
    dict[int, sp.Expr],
    tuple[sp.Symbol, ...],
]:
    h13 = sp.expand(R8.H + c51 * R8.p**5 * R8.y + c23 * R8.p**2 * R8.y**3)
    grows = {row: tcoeff(h13, row) for row in range(4, 14)}

    # Independently verify the complete weight-13 face and its lossless
    # THM-4298 leading channels before using any bracket equation.
    face_indices = [
        (i, j)
        for j in range(14)
        for i in range(14)
        if 2 * i + 3 * j == 13
    ]
    check(face_indices == [(5, 1), (2, 3)], "complete weight-13 face")
    deltas = [exact(grows[row] - tcoeff(R8.H, row)) for row in range(7, 14)]
    expected = [
        c51 * x,
        (6 * c51 + c23) * x**3,
        (15 * c51 + 5 * c23) * x**5,
        (20 * c51 + 10 * c23) * x**7,
        (15 * c51 + 10 * c23) * x**9,
        (6 * c51 + 5 * c23) * x**11,
        (c51 + c23) * x**13,
    ]
    check(all(exact(a - b) == 0 for a, b in zip(deltas, expected)), "weight-13 face rows")

    arows = [R8.A0, R8.A1, R8.A2, R8.A3]
    crows = [R8.C0, R8.C1, R8.C2, R8.C3]
    substitutions: dict[sp.Symbol, sp.Expr] = {}
    solve_symbols = {5: R8.upsilon5, 6: R8.U, 7: R8.W, 8: R8.Z}

    for n in range(4, 8):
        abase, cbase = R8.particular_row(n, R8.B_row(n, arows, crows))
        trial_a = arows + [abase]
        trial_c = crows + [cbase]
        m = n + 1
        difference = sp.expand(
            grows[m].subs(substitutions) - R8.predicted_G(m, trial_a, trial_c)
        )
        moment = R9.primitive_student_row(m)
        obstruction = exact(sum(moment[j] * xcoeff(difference, j) for j in range(m + 1)))
        answers = sp.solve(obstruction, solve_symbols[m])
        check(len(answers) == 1, f"row {m} scalar response")
        substitutions[solve_symbols[m]] = exact(answers[0])
        target = sp.expand(difference.subs(solve_symbols[m], answers[0]))
        theta = R8.tangent_solve(m, target)
        arows.append(sp.expand(abase + theta * sp.diff(R8.A0, x)))
        crows.append(sp.expand(cbase + theta * sp.diff(R8.C0, x)))
        check(
            exact(R8.predicted_G(m, arows, crows) - grows[m].subs(substitutions)) == 0,
            f"row {m} reconstructed bracket",
        )

    theta8 = tuple(sp.symbols("audit_theta8_0:9"))
    tangent8 = sum(theta8[j] * x**j for j in range(9))
    abase8, cbase8 = R8.particular_row(8, R8.B_row(8, arows, crows))
    arows.append(sp.expand(abase8 + tangent8 * sp.diff(R8.A0, x)))
    crows.append(sp.expand(cbase8 + tangent8 * sp.diff(R8.C0, x)))
    return arows, crows, substitutions, grows, theta8


def main() -> None:
    arows, crows, bracket_subs, grows, theta8 = build_through_row_eight()

    depth8, ap8, cp8 = depth_equations(arows, crows, 8)
    terminal8, conditions8, _, _, pivots8 = eliminate_linear(depth8, theta8)
    check((ap8.shape, ap8.rank()) == ((63, 131), 51), "row-8 P2 universe")
    check((cp8.shape, cp8.rank()) == ((72, 204), 63), "row-8 P3 universe")
    check(pivots8 == (7, 8) and len(conditions8) == 3, "row-8 terminal fibre")
    gate = {
        R8.Delta: sp.Rational(896, 15),
        R8.Theta: sp.Rational(512, 75),
        R8.zeta3: -sp.Rational(3, 2) * Phi,
    }
    check(all(exact(value.subs(gate)) == 0 for value in conditions8), "inherited Hasse gate")

    a8 = apply(arows, gate, terminal8, gate)
    c8 = apply(crows, gate, terminal8, gate)
    g9 = exact(grows[9].subs(bracket_subs).subs(gate))
    row9_select, row9_conditions, _, row9_matrix, _ = eliminate_linear(
        bracket_equations(a8, c8, 9, g9),
        theta8[:7],
    )
    check(row9_matrix.rank() == 7 and len(row9_conditions) == 1, "row-9 bracket gate")
    e9 = row9_conditions[0]

    a8_selected = apply(a8, row9_select)
    c8_selected = apply(c8, row9_select)
    a9_trial, c9_trial, theta9 = append_tangent(a8_selected, c8_selected, 9, "audit_theta9")
    depth9, ap9, cp9 = depth_equations(a9_trial, c9_trial, 9)
    terminal9, depth9_conditions, _, _, pivots9 = eliminate_linear(depth9, theta9)
    check((ap9.shape, ap9.rank()) == ((75, 160), 59), "row-9 P2 universe")
    check((cp9.shape, cp9.rank()) == ((85, 251), 73), "row-9 P3 universe")
    check(pivots9 == (7, 8, 9) and not depth9_conditions, "row-9 depth automatic")

    # ---- Phi=0 branch -------------------------------------------------
    xi_zero_answers = sp.solve(e9.subs(Phi, 0), R8.xi10)
    check(len(xi_zero_answers) == 1, "Phi=0 xi graph")
    zero9 = {Phi: 0, R8.xi10: exact(xi_zero_answers[0])}
    a9z = apply(a9_trial, zero9, terminal9, zero9)
    c9z = apply(c9_trial, zero9, terminal9, zero9)
    g10z = exact(grows[10].subs(bracket_subs).subs(gate).subs(zero9))
    select10z, conditions10z, _, matrix10z, _ = eliminate_linear(
        bracket_equations(a9z, c9z, 10, g10z),
        theta9[:7],
    )
    check(matrix10z.rank() == 7 and len(conditions10z) == 2, "Phi=0 row-10 bracket")
    fzero = sp.Poly(
        18612736875 * eta**2 - 4820239249178624,
        eta,
        domain=sp.QQ,
    )
    check(fzero.degree() == 2 and sp.gcd(fzero, fzero.diff()).degree() == 0, "Phi=0 reduced boundary")
    alpha_condition = next(value for value in conditions10z if R8.alpha11 in value.free_symbols)
    alpha_zero_answers = sp.solve(alpha_condition, R8.alpha11)
    check(len(alpha_zero_answers) == 1, "Phi=0 alpha graph")
    alpha_zero = exact(alpha_zero_answers[0])
    other10z = next(value for value in conditions10z if value != alpha_condition)
    check(reduce_mod(other10z, fzero) == 0, "Phi=0 boundary quadratic")
    zero10 = {R8.alpha11: alpha_zero}

    a9z_selected = apply(a9z, select10z, zero10)
    c9z_selected = apply(c9z, select10z, zero10)
    a10z_trial, c10z_trial, theta10z = append_tangent(
        a9z_selected, c9z_selected, 10, "audit_zero_theta10"
    )
    depth10z, ap10z, cp10z = depth_equations(a10z_trial, c10z_trial, 10)
    terminal10z, conditions_depth10z, _, matrix_depth10z, pivots10z = eliminate_linear(
        depth10z, theta10z
    )
    reduced10z = [reduce_mod(value, fzero) for value in conditions_depth10z]
    reduced10z = [value for value in reduced10z if value != 0]
    check(pivots10z == (8, 9, 10) and len(reduced10z) == 1, "Phi=0 row-10 depth")
    beta_zero_answers = sp.solve(reduced10z[0], R8.beta11)
    check(len(beta_zero_answers) == 1, "Phi=0 beta graph")
    beta_zero = exact(beta_zero_answers[0])
    check(reduce_mod(beta_zero + sp.Rational(6, 5) * eta, fzero) == 0, "Phi=0 beta value")
    beta_zero_subs = {R8.beta11: beta_zero}

    a10z = apply(a10z_trial, beta_zero_subs, terminal10z, beta_zero_subs)
    c10z = apply(c10z_trial, beta_zero_subs, terminal10z, beta_zero_subs)
    g11z = grows[11]
    for mapping in (bracket_subs, gate, zero9, zero10, beta_zero_subs):
        g11z = exact(g11z.subs(mapping))
    free10z = tuple(theta10z[j] for j in range(len(theta10z)) if j not in pivots10z)
    select11z, conditions11z, raw11z, matrix11z, _ = eliminate_linear(
        bracket_equations(a10z, c10z, 11, g11z),
        free10z,
    )
    reduced11z = [reduce_mod(value, fzero) for value in conditions11z]
    reduced11z = [value for value in reduced11z if value != 0]
    literal11z = [reduce_mod(value, fzero) for value in raw11z]
    literal11z = [value for value in literal11z if value != 0]
    check(matrix11z.rank() == 8 and len(reduced11z) == len(literal11z) == 1, "Phi=0 row-11 bracket")
    check(exact(sp.diff(literal11z[0], c51) - sp.Rational(323, 432) * eta) == 0, "Phi=0 literal c51 coefficient")
    c51_zero_answers = sp.solve(reduced11z[0], c51)
    check(len(c51_zero_answers) == 1, "Phi=0 c51 graph")
    c51_zero = reduce_mod(c51_zero_answers[0], fzero)
    c51_zero_subs = {c51: c51_zero}

    a10z_selected = apply(a10z, select11z, c51_zero_subs)
    c10z_selected = apply(c10z, select11z, c51_zero_subs)
    a11z_trial, c11z_trial, theta11z = append_tangent(
        a10z_selected, c10z_selected, 11, "audit_zero_theta11"
    )
    depth11z, ap11z, cp11z = depth_equations(a11z_trial, c11z_trial, 11)
    terminal11z, conditions_depth11z, _, matrix_depth11z, pivots11z = eliminate_linear(
        depth11z, theta11z
    )
    check(pivots11z == (8, 9, 10, 11), "Phi=0 row-11 depth pivots")
    check(all(reduce_mod(value, fzero) == 0 for value in conditions_depth11z), "Phi=0 row-11 depth automatic")

    a11z = apply(a11z_trial, terminal11z)
    c11z = apply(c11z_trial, terminal11z)
    g12z = grows[12]
    for mapping in (bracket_subs, gate, zero9, zero10, beta_zero_subs, c51_zero_subs):
        g12z = exact(g12z.subs(mapping))
    free11z = tuple(theta11z[j] for j in range(len(theta11z)) if j not in pivots11z)
    _, conditions12z, raw12z, matrix12z, _ = eliminate_linear(
        bracket_equations(a11z, c11z, 12, g12z),
        free11z,
    )
    reduced12z = [reduce_mod(value, fzero) for value in conditions12z]
    reduced12z = [value for value in reduced12z if value != 0]
    literal12z = [reduce_mod(value, fzero) for value in raw12z]
    literal12z = [value for value in literal12z if value != 0]
    check(matrix12z.rank() == 8 and len(reduced12z) == len(literal12z) == 2, "Phi=0 row-12 bracket")
    coefficients23 = [exact(sp.diff(value, c23)) for value in literal12z]
    check(
        coefficients23 == [sp.Rational(323, 360) * eta, -sp.Rational(5, 7)],
        "Phi=0 literal c23 coefficients",
    )
    a0, a1 = coefficients23
    b0 = exact(literal12z[0].subs(c23, 0))
    b1 = exact(literal12z[1].subs(c23, 0))
    incompatibility = reduce_mod(a0 * b1 - a1 * b0, fzero)
    incompatibility_num = sp.fraction(exact(incompatibility))[0]
    incompatibility_poly = sp.Poly(incompatibility_num, eta, domain=sp.QQ)
    check(incompatibility != 0 and sp.gcd(incompatibility_poly, fzero).degree() == 0, "Phi=0 row-12 extinction")

    # ---- Phi!=0 branch ------------------------------------------------
    alpha_answers = sp.solve(e9, R8.alpha11)
    check(len(alpha_answers) == 1, "Phi!=0 alpha graph")
    alpha_graph = exact(alpha_answers[0])
    alpha_den = sp.Poly(sp.fraction(alpha_graph)[1], Phi, eta, domain=sp.QQ)
    check(len(alpha_den.terms()) == 1 and alpha_den.terms()[0][0] == (1, 0), "Phi!=0 alpha localization")
    alpha_subs = {R8.alpha11: alpha_graph}
    a9 = apply(a9_trial, alpha_subs, terminal9, alpha_subs)
    c9 = apply(c9_trial, alpha_subs, terminal9, alpha_subs)
    g10 = exact(grows[10].subs(bracket_subs).subs(gate).subs(alpha_subs))
    select10, conditions10, _, matrix10, _ = eliminate_linear(
        bracket_equations(a9, c9, 10, g10),
        theta9[:7],
    )
    check(matrix10.rank() == 7 and len(conditions10) == 2, "Phi!=0 row-10 bracket")
    responses10 = sp.solve(conditions10, (R8.xi10, R8.beta11), dict=True, simplify=False)
    check(len(responses10) == 1, "Phi!=0 row-10 response graph")
    response10 = {symbol: exact(value) for symbol, value in responses10[0].items()}
    beta_den = sp.Poly(sp.fraction(response10[R8.beta11])[1], Phi, eta, domain=sp.QQ)
    check(len(beta_den.terms()) == 1 and beta_den.terms()[0][0] == (2, 0), "Phi!=0 beta localization")

    a9_selected = apply(a9, select10, response10)
    c9_selected = apply(c9, select10, response10)
    a10_trial, c10_trial, theta10 = append_tangent(
        a9_selected, c9_selected, 10, "audit_theta10"
    )
    depth10, ap10, cp10 = depth_equations(a10_trial, c10_trial, 10)
    terminal10, conditions_depth10, _, matrix_depth10, pivots10 = eliminate_linear(
        depth10, theta10
    )
    check(pivots10 == (8, 9, 10) and len(conditions_depth10) == 1, "Phi!=0 row-10 depth")
    c51_answers = sp.solve(conditions_depth10[0], c51)
    check(len(c51_answers) == 1, "Phi!=0 transverse c51 graph")
    c51_graph = exact(c51_answers[0])
    old_elliptic = (
        7231154026500 * Phi**3
        + 50541940696500 * Phi**2 * eta
        + 6793915500000 * Phi * eta**2
        - 631918028977864704 * Phi
        + 353642000625 * eta**3
        - 91584545734393856 * eta
    )
    check(
        exact(c51_graph - old_elliptic / (707284001250 * Phi**2)) == 0,
        "c51 absorbs the old elliptic cubic",
    )
    c51_subs = {c51: c51_graph}

    a10 = apply(a10_trial, c51_subs, terminal10, c51_subs)
    c10 = apply(c10_trial, c51_subs, terminal10, c51_subs)
    g11 = grows[11]
    for mapping in (bracket_subs, gate, alpha_subs, response10, c51_subs):
        g11 = exact(g11.subs(mapping))
    free10 = tuple(theta10[j] for j in range(len(theta10)) if j not in pivots10)
    select11, conditions11, raw11, matrix11, _ = eliminate_linear(
        bracket_equations(a10, c10, 11, g11),
        free10,
    )
    raw11_nonzero = [(index, value) for index, value in enumerate(raw11) if value != 0]
    check(matrix11.rank() == 8 and len(conditions11) == len(raw11_nonzero) == 1, "Phi!=0 row-11 bracket")
    check(raw11_nonzero[0][0] == 8, "Phi!=0 row-11 residual degree")
    check(
        exact(sp.diff(raw11_nonzero[0][1], c23) - sp.Rational(323, 504) * Phi) == 0,
        "Phi!=0 literal c23 coefficient",
    )
    c23_answers = sp.solve(conditions11[0], c23)
    check(len(c23_answers) == 1, "Phi!=0 c23 graph")
    c23_graph = exact(c23_answers[0])
    c23_den = sp.Poly(sp.fraction(c23_graph)[1], Phi, eta, domain=sp.QQ)
    check(len(c23_den.terms()) == 1 and c23_den.terms()[0][0] == (3, 0), "Phi!=0 c23 localization")
    c23_subs = {c23: c23_graph}

    a10_selected = apply(a10, select11, c23_subs)
    c10_selected = apply(c10, select11, c23_subs)
    a11_trial, c11_trial, theta11 = append_tangent(
        a10_selected, c10_selected, 11, "audit_theta11"
    )
    depth11, ap11, cp11 = depth_equations(a11_trial, c11_trial, 11)
    terminal11, conditions_depth11, _, matrix_depth11, pivots11 = eliminate_linear(
        depth11, theta11
    )
    check(pivots11 == (8, 9, 10, 11) and not conditions_depth11, "Phi!=0 row-11 depth")

    a11 = apply(a11_trial, terminal11)
    c11 = apply(c11_trial, terminal11)
    g12 = grows[12]
    for mapping in (bracket_subs, gate, alpha_subs, response10, c51_subs, c23_subs):
        g12 = exact(g12.subs(mapping))
    free11 = tuple(theta11[j] for j in range(len(theta11)) if j not in pivots11)
    select12, conditions12, _, matrix12, _ = eliminate_linear(
        bracket_equations(a11, c11, 12, g12),
        free11,
    )
    check(matrix12.rank() == 8 and len(conditions12) == 2, "Phi!=0 row-12 bracket")
    parities = []
    for condition in conditions12:
        totals = {(i + j) % 2 for (i, j), _ in sp.Poly(condition, Phi, eta).terms()}
        check(len(totals) == 1, "row-12 condition has a sign character")
        parities.append(next(iter(totals)))
    odd_condition = conditions12[parities.index(1)]
    even_condition = conditions12[parities.index(0)]
    qodd = quotient_polynomial(odd_condition, 1)
    qeven = quotient_polynomial(even_condition, 0)
    check((qodd.degree(Y), qodd.degree(ratio)) == (3, 5), "row-12 odd bidegree")
    check((qeven.degree(Y), qeven.degree(ratio)) == (3, 4), "row-12 even bidegree")

    resultant = sp.Poly(
        sp.resultant(qodd.as_expr(), qeven.as_expr(), Y),
        ratio,
        domain=sp.QQ,
    ).primitive()[1]
    factor_data = [
        (sp.Poly(factor, ratio, domain=sp.QQ), power)
        for factor, power in sp.factor_list(resultant.as_expr())[1]
    ]
    check(sorted((factor.degree(), power) for factor, power in factor_data) == [(1, 1), (13, 1)], "row-12 resultant factors")
    linear_resultant = next(factor for factor, _ in factor_data if factor.degree() == 1)
    check(linear_resultant.monic() == sp.Poly(9 * ratio + 8, ratio).monic(), "extraneous resultant ratio")
    special_gcd = sp.gcd(
        sp.Poly(qodd.as_expr().subs(ratio, -sp.Rational(8, 9)), Y, domain=sp.QQ),
        sp.Poly(qeven.as_expr().subs(ratio, -sp.Rational(8, 9)), Y, domain=sp.QQ),
    )
    check(special_gcd.degree() == 0, "extraneous ratio has no affine Y")

    basis = sp.groebner([qodd.as_expr(), qeven.as_expr()], Y, ratio, order="lex", domain=sp.QQ)
    check(basis.is_zero_dimensional and len(basis.polys) == 2, "row-12 quotient Groebner basis")
    linear_graph = sp.Poly(basis.polys[0].as_expr(), Y, ratio, domain=sp.QQ)
    k13 = sp.Poly(basis.polys[1].as_expr(), ratio, domain=sp.QQ).primitive()[1]
    check(linear_graph.degree(Y) == 1 and k13.degree() == 13, "row-12 quotient graph")
    graph_a = sp.Poly(sp.expand(linear_graph.as_expr()).coeff(Y, 1), ratio, domain=sp.QQ)
    graph_b = sp.Poly(linear_graph.as_expr().subs(Y, 0), ratio, domain=sp.QQ)
    check(sp.gcd(k13, k13.diff()).degree() == 0, "K13 squarefree")
    check(sp.gcd(k13, graph_a).degree() == 0, "Y graph unique")
    check(sp.gcd(k13, graph_b).degree() == 0, "Y graph nonzero")
    check(
        [(sp.Poly(factor, ratio).degree(), power) for factor, power in sp.factor_list(k13.as_expr())[1]]
        == [(13, 1)],
        "K13 irreducible over Q",
    )

    # Positive bracket-carrier control only; its nonzero depth value is a
    # hostile control.  It is not used to infer a negative characteristic-zero
    # fibre census.
    control = {ratio: 6, Y: 9}
    check(mod_value(qodd.as_expr(), control, 23) == 0, "F23 odd carrier control")
    check(mod_value(qeven.as_expr(), control, 23) == 0, "F23 even carrier control")

    a11_selected = apply(a11, select12)
    c11_selected = apply(c11, select12)
    a12_trial, c12_trial, theta12 = append_tangent(
        a11_selected, c11_selected, 12, "audit_theta12"
    )
    depth12, ap12, cp12 = depth_equations(a12_trial, c12_trial, 12)
    _, conditions_depth12, raw_depth12, matrix_depth12, pivots12 = eliminate_linear(
        depth12, theta12
    )
    check((ap12.shape, ap12.rank()) == ((117, 267), 87), "row-12 P2 universe")
    check((cp12.shape, cp12.rank()) == ((130, 424), 105), "row-12 P3 universe")
    check(pivots12 == (9, 10, 11, 12) and len(conditions_depth12) == 1, "row-12 terminal depth")
    n12 = 272008125 * Phi**2 - 43740000 * Phi * eta + 10213932924928
    check(proportional(conditions_depth12[0], n12, Phi, eta), "row-12 depth equation N")
    raw_nonzero = [(index, value) for index, value in enumerate(raw_depth12) if value != 0]
    check(len(raw_nonzero) == 1 and raw_nonzero[0][0] == 50, "row-12 literal depth support")
    check(exact(raw_nonzero[0][1] + n12 / 288360000) == 0, "row-12 literal depth scale")

    nq_slope = 272008125 - 43740000 * ratio
    nq_constant = sp.Integer(10213932924928)
    check(nq_constant != 0, "N slope-zero branch is empty")
    substituted_odd = sp.Poly(
        sp.together(qodd.as_expr().subs(Y, -nq_constant / nq_slope) * nq_slope**3),
        ratio,
        domain=sp.QQ,
    ).primitive()[1]
    substituted_even = sp.Poly(
        sp.together(qeven.as_expr().subs(Y, -nq_constant / nq_slope) * nq_slope**3),
        ratio,
        domain=sp.QQ,
    ).primitive()[1]
    check((substituted_odd.degree(), substituted_even.degree()) == (6, 5), "depth eliminant degrees")
    check(sp.gcd(substituted_odd, substituted_even).degree() == 0, "depth disjoint over Q")
    raw_f0 = sp.Poly(3 * ratio**6 - 4 * ratio**5 - 2 * ratio**2 - 2 * ratio - 2, ratio, modulus=11)
    raw_f1 = sp.Poly(5 * ratio**5 - 5 * ratio**3 + 5 * ratio**2 + 2 * ratio + 5, ratio, modulus=11)
    check(sp.Poly(substituted_odd.as_expr(), ratio, modulus=11).monic() == raw_f0.monic(), "mod-11 odd reduction")
    check(sp.Poly(substituted_even.as_expr(), ratio, modulus=11).monic() == raw_f1.monic(), "mod-11 even reduction")
    bezout_s = 2 * ratio**4 + 4 * ratio**3 + 3 * ratio**2 - 5 * ratio + 2
    bezout_t = ratio**5 - 3 * ratio**4 - 2 * ratio**3 - 3 * ratio**2 + 5 * ratio + 1
    check(
        sp.Poly(bezout_s * raw_f0.as_expr() + bezout_t * raw_f1.as_expr() - 1, ratio, modulus=11).is_zero,
        "mod-11 Bezout identity",
    )
    check(mod_value((272008125 - 43740000 * ratio) * Y + nq_constant, control, 23) == 4, "F23 hostile depth control")

    # The exact THM-4130 overlap: alpha^2=A5 is the sign sheet.  This is a
    # coefficient-atlas map, not a map from the physical source P1_X.
    A5, alpha_root, f21, f31 = sp.symbols("A5 alpha_root f21 f31", nonzero=True)
    phi_hat = -alpha_root * A5**3 * f21 / 2
    eta_hat = -alpha_root * A5**4 * f31 / 2
    check(exact(eta_hat / phi_hat - A5 * f31 / f21) == 0, "normalization quotient r")
    check(
        exact(phi_hat**2).subs(alpha_root**2, A5) == A5**7 * f21**2 / 4,
        "normalization sign square",
    )
    check(phi_hat.subs(alpha_root, -alpha_root) == -phi_hat, "normalization Phi sign")
    check(eta_hat.subs(alpha_root, -alpha_root) == -eta_hat, "normalization eta sign")

    print("THM-4388 independent weight-13 row-12 audit")
    print("imports=audited_THM4308_R8_and_THM4315_R9; primary_THM4388_import=no")
    print("field=algebraically_closed_characteristic_zero")
    print("face=c51*p^5*y+c23*p^2*y^3 channels=(c51,6*c51+c23) lossless=yes")
    print("Phi_zero=F_degree2_squarefree; row11_c51_unit=323*eta/432; row11_depth=automatic")
    print("Phi_zero_row12=c23_coefficients=(323*eta/360,-5/7); compatibility_mod_F=unit; survivors=0")
    print("Phi_nonzero_row10=D-707284001250*Phi^2*c51; old_elliptic_carrier_absorbed=yes")
    print("Phi_nonzero_row11=c23_coefficient=323*Phi/504; localized_graph=lawful")
    print("row12_bracket=13_reduced_quotient_points,26_signed_geometric_lifts")
    print("row12_depth=N=(272008125-43740000*r)*Y+10213932924928")
    print("row12_extinction=substituted_degrees=(6,5); mod11_Bezout=1; survivors=0")
    print("finite_field_roles=F23_positive_carrier_and_hostile_depth_control; F11_good_reduction_Bezout_for_Q_coprimality")
    print("normalization_channels=f21=[p^2*y](R/gamma),f31=[p^3*y](R/gamma)")
    print("normalization_bridge=Phi_hat=-alpha*A5^3*f21/2 eta_hat=-alpha*A5^4*f31/2 alpha^2=A5")
    print("type_note=coefficient_atlas_mu2_only; no_physical_P1_X_to_Ebar_map")
    print(f"checks={CHECKS} result=PASS")


if __name__ == "__main__":
    main()
