#!/usr/bin/env python3
"""Primary exact certificate for THM-4388.

The declared extension of THM-4308 is

    c51*p**5*y + c23*p**2*y**3.

The certificate reconstructs the bracket rows from the literal source,
applies the exact projected P_2/P_3 depth matrices, treats Phi=0 and Phi!=0
separately, and proves that the complete weight-at-most-thirteen extension
is extinct by row twelve.  It imports only proved row operators and does not
modify their globals.
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
c51, c23 = sp.symbols("c51 c23")
Y, ratio = sp.symbols("Y ratio")
H13 = sp.expand(R8.H + c51 * R8.p**5 * R8.y + c23 * R8.p**2 * R8.y**3)


def exact(value: sp.Expr) -> sp.Expr:
    return sp.cancel(sp.together(sp.expand(value)))


def xcoeff(value: sp.Expr, degree: int) -> sp.Expr:
    return sp.expand(value).coeff(x, degree)


def tcoeff(value: sp.Expr, degree: int) -> sp.Expr:
    return sp.expand(value).coeff(t, degree)


def primitive(expression: sp.Expr, *variables: sp.Symbol) -> sp.Expr:
    numerator, _ = sp.fraction(exact(expression))
    return sp.Poly(numerator, *variables, domain=sp.QQ).primitive()[1].as_expr()


def selected_solution(
    matrix: sp.Matrix,
    rhs: sp.Matrix,
    columns: tuple[int, ...],
) -> tuple[sp.Expr, ...]:
    row_indices = tuple(matrix.T.rref()[1])
    square = matrix.extract(row_indices, columns)
    if square.rows != square.cols or square.det() == 0:
        raise AssertionError("selected block is not invertible")
    answer = square.inv() * rhs.extract(row_indices, (0,))
    return tuple(exact(value) for value in answer)


def source_rows(first: int, last: int) -> dict[int, sp.Expr]:
    return {row: sp.expand(tcoeff(H13, row)) for row in range(first, last + 1)}


def build_bracket_rows() -> tuple[
    list[sp.Expr],
    list[sp.Expr],
    dict[sp.Symbol, sp.Expr],
    dict[int, sp.Expr],
]:
    grows = source_rows(4, 13)
    arows = [R8.A0, R8.A1, R8.A2, R8.A3]
    crows = [R8.C0, R8.C1, R8.C2, R8.C3]
    substitutions: dict[sp.Symbol, sp.Expr] = {}
    solve_symbols = {5: R8.upsilon5, 6: R8.U, 7: R8.W, 8: R8.Z}

    if exact(R8.predicted_G(4, arows, crows) - grows[4]) != 0:
        raise AssertionError("inherited G4 mismatch")

    for n in range(4, 8):
        bvalue = R8.B_row(n, arows, crows)
        abase, cbase = R8.particular_row(n, bvalue)
        trial_a = arows + [abase]
        trial_c = crows + [cbase]
        m = n + 1
        difference = sp.expand(
            grows[m].subs(substitutions)
            - R8.predicted_G(m, trial_a, trial_c)
        )
        moment = R9.primitive_student_row(m)
        obstruction = exact(
            sum(moment[j] * xcoeff(difference, j) for j in range(m + 1))
        )
        symbol = solve_symbols[m]
        answer = sp.solve(obstruction, symbol)
        if len(answer) != 1:
            raise AssertionError(f"row {m} does not solve uniquely for {symbol}")
        substitutions[symbol] = exact(answer[0])
        target = sp.expand(difference.subs(symbol, substitutions[symbol]))
        theta = R8.tangent_solve(m, target)
        arows.append(sp.expand(abase + theta * sp.diff(R8.A0, x)))
        crows.append(sp.expand(cbase + theta * sp.diff(R8.C0, x)))
        if exact(R8.predicted_G(m, arows, crows) - grows[m].subs(substitutions)) != 0:
            raise AssertionError(f"row {m} bracket mismatch")
        print(f"BRACKET_E{m}={sp.sstr(primitive(obstruction, *sorted(obstruction.free_symbols, key=str)))}")
        print(f"BRACKET_{symbol}={sp.sstr(substitutions[symbol])}")

    b8 = R8.B_row(8, arows, crows)
    a8base, c8base = R8.particular_row(8, b8)
    theta8_symbols = tuple(sp.symbols("w13_theta8_0:9"))
    theta8 = sum(theta8_symbols[j] * x**j for j in range(9))
    arows.append(sp.expand(a8base + theta8 * sp.diff(R8.A0, x)))
    crows.append(sp.expand(c8base + theta8 * sp.diff(R8.C0, x)))
    return arows, crows, substitutions, grows


def projected_depth_residuals(
    arows: list[sp.Expr],
    crows: list[sp.Expr],
    row: int,
) -> tuple[list[sp.Expr], sp.Matrix, sp.Matrix]:
    acoords, amatrix = R9.depth_matrix(2, row)
    ccoords, cmatrix = R9.depth_matrix(3, row)
    avec = sp.Matrix([xcoeff(arows[n], degree) for n, degree in acoords])
    cvec = sp.Matrix([xcoeff(crows[n], degree) for n, degree in ccoords])
    residuals = [sp.expand((left.T * avec)[0]) for left in amatrix.T.nullspace()]
    residuals += [sp.expand((left.T * cvec)[0]) for left in cmatrix.T.nullspace()]
    return residuals, amatrix, cmatrix


def eliminate_linear_fibre(
    residuals: list[sp.Expr],
    fibre_symbols: tuple[sp.Symbol, ...],
    tag: str,
    print_conditions: bool = True,
) -> tuple[dict[sp.Symbol, sp.Expr], list[sp.Expr], sp.Matrix]:
    matrix, rhs = sp.linear_eq_to_matrix(residuals, fibre_symbols)
    rank = matrix.rank()
    pivots = tuple(matrix.rref()[1])
    row_pivots = tuple(matrix.T.rref()[1])
    pivot_matrix = matrix.extract(row_pivots, pivots)
    solution_vector = pivot_matrix.inv() * rhs.extract(row_pivots, (0,))
    substitutions = {
        fibre_symbols[column]: exact(solution_vector[index])
        for index, column in enumerate(pivots)
    }
    reduced = [exact(value.subs(substitutions)) for value in residuals]
    nonzero = []
    for value in reduced:
        if value == 0:
            continue
        candidate = primitive(value, *sorted(value.free_symbols, key=str))
        if all(exact(candidate - old) != 0 and exact(candidate + old) != 0 for old in nonzero):
            nonzero.append(candidate)
    print(f"{tag}_FIBRE_MATRIX shape={matrix.shape} rank={rank} pivots={pivots}")
    if print_conditions:
        for index, value in enumerate(nonzero):
            print(f"{tag}_CONSISTENCY_{index}={sp.sstr(sp.factor(value))}")
    return substitutions, nonzero, matrix


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
    anew = sp.expand(abase + theta * sp.diff(R8.A0, x))
    cnew = sp.expand(cbase + theta * sp.diff(R8.C0, x))
    return arows + [anew], crows + [cnew], symbols


def solve_bracket_fibre(
    arows: list[sp.Expr],
    crows: list[sp.Expr],
    row: int,
    source: sp.Expr,
    fibre_symbols: tuple[sp.Symbol, ...],
    tag: str,
    print_conditions: bool = True,
) -> tuple[dict[sp.Symbol, sp.Expr], list[sp.Expr], sp.Matrix]:
    difference = exact(source - R8.predicted_G(row, arows, crows))
    equations = [xcoeff(difference, degree) for degree in range(row + 1)]
    return eliminate_linear_fibre(
        equations,
        fibre_symbols,
        tag,
        print_conditions=print_conditions,
    )


def impose_substitutions(
    rows: list[sp.Expr],
    *substitutions: dict[sp.Symbol, sp.Expr],
) -> list[sp.Expr]:
    result = rows
    for mapping in substitutions:
        result = [exact(row.subs(mapping)) for row in result]
    return result


def sign_quotient_polynomial(expression: sp.Expr, parity: int) -> sp.Poly:
    """Substitute eta=ratio*Phi and divide by Phi**parity, then Phi^2=Y."""

    polynomial = sp.Poly(expression, R8.Phi, R8.eta, domain=sp.QQ)
    answer = sp.Integer(0)
    for (pdegree, edegree), coefficient in polynomial.terms():
        total = pdegree + edegree
        if total % 2 != parity:
            raise AssertionError("mixed sign parity in row-twelve condition")
        answer += coefficient * ratio**edegree * Y ** ((total - parity) // 2)
    return sp.Poly(answer, Y, ratio, domain=sp.QQ).primitive()[1]


def quotient_components(expression: sp.Expr) -> tuple[sp.Expr, sp.Expr]:
    """Return even and odd-P components after eta=ratio*Phi, ignoring units."""

    numerator, denominator = sp.fraction(exact(expression))
    den_poly = sp.Poly(denominator, R8.Phi, R8.eta, domain=sp.QQ)
    if len(den_poly.terms()) != 1 or den_poly.terms()[0][0][1] != 0:
        raise AssertionError(f"unexpected localized denominator: {denominator}")
    polynomial = sp.Poly(numerator, R8.Phi, R8.eta, domain=sp.QQ)
    components = [sp.Integer(0), sp.Integer(0)]
    for (pdegree, edegree), coefficient in polynomial.terms():
        total = pdegree + edegree
        parity = total % 2
        components[parity] += (
            coefficient * ratio**edegree * Y ** ((total - parity) // 2)
        )
    return tuple(sp.expand(component) for component in components)  # type: ignore[return-value]


def quotient_vanishes(expression: sp.Expr, basis: sp.GroebnerBasis) -> bool:
    return all(basis.reduce(component)[1] == 0 for component in quotient_components(expression))


def mod_value(expression: sp.Expr, substitutions: dict[sp.Symbol, int], prime: int) -> int:
    value = exact(expression.subs(substitutions))
    if value.free_symbols:
        raise AssertionError(f"non-scalar modular value: {value}")
    numerator, denominator = sp.fraction(value)
    den = int(denominator) % prime
    if den == 0:
        raise ZeroDivisionError
    return int(numerator) % prime * pow(den, -1, prime) % prime


def reduce_mod_univariate(expression: sp.Expr, modulus: sp.Poly) -> sp.Expr:
    value = exact(expression)
    numerator, denominator = sp.fraction(value)
    variable = modulus.gens[0]
    coefficient_symbols = sorted(
        (numerator.free_symbols | denominator.free_symbols) - {variable},
        key=str,
    )
    domain = sp.QQ.frac_field(*coefficient_symbols) if coefficient_symbols else sp.QQ
    num_poly = sp.Poly(numerator, variable, domain=domain)
    den_poly = sp.Poly(denominator, variable, domain=domain)
    mod_poly = sp.Poly(modulus.as_expr(), variable, domain=domain)
    inverse = sp.invert(den_poly, mod_poly)
    return sp.factor((num_poly * inverse).rem(mod_poly).as_expr())


def main() -> None:
    weight13_pairs = [
        (a, b)
        for a in range(7)
        for b in range(5)
        if 2 * a + 3 * b == 13
    ]
    if weight13_pairs != [(2, 3), (5, 1)]:
        raise AssertionError("residual-weight-13 face enumeration changed")
    print("WEIGHT13_FACE monomials=(p^5*y,p^2*y^3) complete=yes")
    print("WEIGHT13_CHANNEL_SIDECAR h7=c51 h8=6*c51+c23 inverse_c23=h8-6*h7")
    grows = source_rows(4, 13)
    for row in range(4, 14):
        delta = exact(grows[row] - tcoeff(R8.H, row))
        if delta != 0:
            print(f"SOURCE_DELTA_G{row}={sp.sstr(delta)}")

    arows, crows, bracket_subs, grows = build_bracket_rows()
    theta8 = tuple(sp.symbols("w13_theta8_0:9"))
    depth8, amatrix8, cmatrix8 = projected_depth_residuals(arows, crows, 8)
    print(
        "ROW8_DEPTH_UNIVERSE "
        f"P2={amatrix8.shape}/{amatrix8.rank()} "
        f"P3={cmatrix8.shape}/{cmatrix8.rank()}"
    )
    terminal8, depth8_conditions, _ = eliminate_linear_fibre(
        depth8,
        theta8,
        "ROW8_DEPTH",
    )
    if len(depth8_conditions) != 3:
        raise AssertionError("row-eight depth did not retain three conditions")
    gate = {
        R8.Delta: sp.Rational(896, 15),
        R8.Theta: sp.Rational(512, 75),
        R8.zeta3: -sp.Rational(3, 2) * R8.Phi,
    }
    if any(exact(condition.subs(gate)) != 0 for condition in depth8_conditions):
        raise AssertionError("inherited Hasse gate does not solve depth conditions")

    a8 = impose_substitutions(arows, gate, terminal8, gate)
    c8rows = impose_substitutions(crows, gate, terminal8, gate)
    free8 = theta8[:7]
    g9 = exact(grows[9].subs(bracket_subs).subs(gate))
    row9_select, row9_conditions, row9_matrix = solve_bracket_fibre(
        a8,
        c8rows,
        9,
        g9,
        free8,
        "ROW9_BRACKET",
    )
    print(f"ROW9_BRACKET_RANK={row9_matrix.rank()}")
    print(f"ROW9_CONDITION_COUNT={len(row9_conditions)}")

    a8_selected = impose_substitutions(a8, row9_select)
    c8_selected = impose_substitutions(c8rows, row9_select)
    a9_trial, c9_trial, theta9 = append_tangent(
        a8_selected,
        c8_selected,
        9,
        "w13_theta9",
    )
    depth9, amatrix9, cmatrix9 = projected_depth_residuals(
        a9_trial,
        c9_trial,
        9,
    )
    terminal9, depth9_conditions, _ = eliminate_linear_fibre(
        depth9,
        theta9,
        "ROW9_DEPTH",
    )
    print(
        "ROW9_DEPTH_UNIVERSE "
        f"P2={amatrix9.shape}/{amatrix9.rank()} "
        f"P3={cmatrix9.shape}/{cmatrix9.rank()}"
    )
    print(f"ROW9_DEPTH_CONDITION_COUNT={len(depth9_conditions)}")

    # Boundary branch Phi=0.  This must be compiled separately because the
    # generic alpha graph divides by Phi.
    e9 = row9_conditions[0]
    xi0_answers = sp.solve(e9.subs(R8.Phi, 0), R8.xi10)
    if len(xi0_answers) != 1:
        raise AssertionError("Phi-zero E9 xi graph missing")
    xi0 = exact(xi0_answers[0])
    zero9_subs = {R8.Phi: 0, R8.xi10: xi0}
    a9_zero = impose_substitutions(a9_trial, zero9_subs, terminal9, zero9_subs)
    c9_zero = impose_substitutions(c9_trial, zero9_subs, terminal9, zero9_subs)
    g10_zero = grows[10]
    for mapping in (bracket_subs, gate, zero9_subs):
        g10_zero = exact(g10_zero.subs(mapping))
    row10_zero_select, row10_zero_conditions, row10_zero_matrix = solve_bracket_fibre(
        a9_zero,
        c9_zero,
        10,
        g10_zero,
        theta9[:7],
        "PHI_ZERO_ROW10_BRACKET",
    )
    print(
        "PHI_ZERO_ROW10_BRACKET "
        f"matrix={row10_zero_matrix.shape}/{row10_zero_matrix.rank()} "
        f"conditions={len(row10_zero_conditions)}"
    )
    fzero = sp.Poly(
        18612736875 * R8.eta**2 - 4820239249178624,
        R8.eta,
        domain=sp.QQ,
    )
    if sp.gcd(fzero, fzero.diff()).degree() != 0 or fzero.eval(0) == 0:
        raise AssertionError("Phi-zero boundary quadratic should be reduced and nonzero")
    alpha_zero_answers = sp.solve(row10_zero_conditions[0], R8.alpha11)
    if len(alpha_zero_answers) != 1:
        raise AssertionError("Phi-zero row-ten alpha graph missing")
    alpha_zero = exact(alpha_zero_answers[0])
    if reduce_mod_univariate(row10_zero_conditions[1], fzero) != 0:
        raise AssertionError("Phi-zero row-ten boundary equation changed")
    zero10_subs = {R8.alpha11: alpha_zero}
    a9_zero_selected = impose_substitutions(
        a9_zero,
        row10_zero_select,
        zero10_subs,
    )
    c9_zero_selected = impose_substitutions(
        c9_zero,
        row10_zero_select,
        zero10_subs,
    )
    a10_zero_trial, c10_zero_trial, theta10_zero = append_tangent(
        a9_zero_selected,
        c9_zero_selected,
        10,
        "w13_zero_theta10",
    )
    depth10_zero, amatrix10_zero, cmatrix10_zero = projected_depth_residuals(
        a10_zero_trial,
        c10_zero_trial,
        10,
    )
    terminal10_zero, depth10_zero_conditions, depth10_zero_matrix = eliminate_linear_fibre(
        depth10_zero,
        theta10_zero,
        "PHI_ZERO_ROW10_DEPTH",
        print_conditions=False,
    )
    reduced_depth10_zero = [
        reduce_mod_univariate(value, fzero)
        for value in depth10_zero_conditions
    ]
    reduced_depth10_zero = [value for value in reduced_depth10_zero if value != 0]
    if len(reduced_depth10_zero) != 1:
        raise AssertionError("Phi-zero row-ten depth should leave beta graph")
    beta_zero_answers = sp.solve(reduced_depth10_zero[0], R8.beta11)
    if len(beta_zero_answers) != 1:
        raise AssertionError("Phi-zero row-ten beta graph missing")
    beta_zero = exact(beta_zero_answers[0])
    if reduce_mod_univariate(beta_zero + sp.Rational(6, 5) * R8.eta, fzero) != 0:
        raise AssertionError("Phi-zero beta graph changed")
    print(
        "PHI_ZERO_ROW10_DEPTH "
        f"P2={amatrix10_zero.shape}/{amatrix10_zero.rank()} "
        f"P3={cmatrix10_zero.shape}/{cmatrix10_zero.rank()} "
        f"terminal={depth10_zero_matrix.shape}/{depth10_zero_matrix.rank()} "
        "beta=-6*eta/5"
    )

    beta_zero_subs = {R8.beta11: beta_zero}
    a10_zero = impose_substitutions(
        a10_zero_trial,
        beta_zero_subs,
        terminal10_zero,
        beta_zero_subs,
    )
    c10_zero = impose_substitutions(
        c10_zero_trial,
        beta_zero_subs,
        terminal10_zero,
        beta_zero_subs,
    )
    g11_zero = grows[11]
    for mapping in (
        bracket_subs,
        gate,
        zero9_subs,
        zero10_subs,
        beta_zero_subs,
    ):
        g11_zero = exact(g11_zero.subs(mapping))
    free10_zero = tuple(
        symbol
        for index, symbol in enumerate(theta10_zero)
        if index not in depth10_zero_matrix.rref()[1]
    )
    row11_zero_select, row11_zero_conditions, row11_zero_matrix = solve_bracket_fibre(
        a10_zero,
        c10_zero,
        11,
        g11_zero,
        free10_zero,
        "PHI_ZERO_ROW11_BRACKET",
        print_conditions=False,
    )
    reduced_row11_zero = [
        reduce_mod_univariate(value, fzero)
        for value in row11_zero_conditions
    ]
    reduced_row11_zero = [value for value in reduced_row11_zero if value != 0]
    if len(reduced_row11_zero) != 1:
        raise AssertionError("Phi-zero row-eleven should leave one c51 graph")
    c51_zero_answers = sp.solve(reduced_row11_zero[0], c51)
    if len(c51_zero_answers) != 1:
        raise AssertionError("Phi-zero row-eleven c51 graph missing")
    c51_zero = reduce_mod_univariate(c51_zero_answers[0], fzero)

    difference11_zero = exact(g11_zero - R8.predicted_G(11, a10_zero, c10_zero))
    equations11_zero = [xcoeff(difference11_zero, degree) for degree in range(12)]
    matrix11_zero_raw, rhs11_zero_raw = sp.linear_eq_to_matrix(equations11_zero, free10_zero)
    row_pivots11_zero = tuple(matrix11_zero_raw.T.rref()[1])
    solution11_zero_raw = (
        matrix11_zero_raw.extract(row_pivots11_zero, tuple(range(len(free10_zero)))).inv()
        * rhs11_zero_raw.extract(row_pivots11_zero, (0,))
    )
    select11_zero_raw = {
        free10_zero[j]: exact(solution11_zero_raw[j])
        for j in range(len(free10_zero))
    }
    literal11_zero = [
        reduce_mod_univariate(value.subs(select11_zero_raw), fzero)
        for value in equations11_zero
    ]
    literal11_zero = [value for value in literal11_zero if value != 0]
    if len(literal11_zero) != 1:
        raise AssertionError("Phi-zero row-eleven literal residual support")
    if exact(sp.diff(literal11_zero[0], c51) - sp.Rational(323, 432) * R8.eta) != 0:
        raise AssertionError("Phi-zero row-eleven c51 coefficient changed")
    print(
        "PHI_ZERO_ROW11_BRACKET "
        f"matrix={row11_zero_matrix.shape}/{row11_zero_matrix.rank()} "
        f"literal_c51_coefficient=323*eta/432 c51_graph={sp.sstr(c51_zero)}"
    )

    c51_zero_subs = {c51: c51_zero}
    a10_zero_selected = impose_substitutions(
        a10_zero,
        row11_zero_select,
        c51_zero_subs,
    )
    c10_zero_selected = impose_substitutions(
        c10_zero,
        row11_zero_select,
        c51_zero_subs,
    )
    a11_zero_trial, c11_zero_trial, theta11_zero = append_tangent(
        a10_zero_selected,
        c10_zero_selected,
        11,
        "w13_zero_theta11",
    )
    depth11_zero, amatrix11_zero, cmatrix11_zero = projected_depth_residuals(
        a11_zero_trial,
        c11_zero_trial,
        11,
    )
    terminal11_zero, depth11_zero_conditions, depth11_zero_matrix = eliminate_linear_fibre(
        depth11_zero,
        theta11_zero,
        "PHI_ZERO_ROW11_DEPTH",
        print_conditions=False,
    )
    reduced_depth11_zero = [
        reduce_mod_univariate(value, fzero)
        for value in depth11_zero_conditions
    ]
    if any(value != 0 for value in reduced_depth11_zero):
        raise AssertionError("Phi-zero row-eleven depth should be automatic")
    print(
        "PHI_ZERO_ROW11_DEPTH "
        f"P2={amatrix11_zero.shape}/{amatrix11_zero.rank()} "
        f"P3={cmatrix11_zero.shape}/{cmatrix11_zero.rank()} "
        f"terminal={depth11_zero_matrix.shape}/{depth11_zero_matrix.rank()} automatic_mod_F=yes"
    )

    a11_zero = impose_substitutions(a11_zero_trial, terminal11_zero)
    c11_zero = impose_substitutions(c11_zero_trial, terminal11_zero)
    g12_zero = grows[12]
    for mapping in (
        bracket_subs,
        gate,
        zero9_subs,
        zero10_subs,
        beta_zero_subs,
        c51_zero_subs,
    ):
        g12_zero = exact(g12_zero.subs(mapping))
    free11_zero = tuple(
        symbol
        for index, symbol in enumerate(theta11_zero)
        if index not in depth11_zero_matrix.rref()[1]
    )
    row12_zero_select, row12_zero_conditions, row12_zero_matrix = solve_bracket_fibre(
        a11_zero,
        c11_zero,
        12,
        g12_zero,
        free11_zero,
        "PHI_ZERO_ROW12_BRACKET",
        print_conditions=False,
    )
    reduced_row12_zero = [
        reduce_mod_univariate(value, fzero)
        for value in row12_zero_conditions
    ]
    reduced_row12_zero = [value for value in reduced_row12_zero if value != 0]
    print(
        "PHI_ZERO_ROW12_BRACKET "
        f"matrix={row12_zero_matrix.shape}/{row12_zero_matrix.rank()} "
        f"reduced_conditions={len(reduced_row12_zero)}"
    )
    for index, value in enumerate(reduced_row12_zero):
        print(f"PHI_ZERO_ROW12_REDUCED_{index}={sp.sstr(value)}")
    difference12_zero = exact(g12_zero - R8.predicted_G(12, a11_zero, c11_zero))
    equations12_zero = [xcoeff(difference12_zero, degree) for degree in range(13)]
    matrix12_zero_raw, rhs12_zero_raw = sp.linear_eq_to_matrix(equations12_zero, free11_zero)
    row_pivots12_zero = tuple(matrix12_zero_raw.T.rref()[1])
    solution12_zero_raw = (
        matrix12_zero_raw.extract(row_pivots12_zero, tuple(range(len(free11_zero)))).inv()
        * rhs12_zero_raw.extract(row_pivots12_zero, (0,))
    )
    select12_zero_raw = {
        free11_zero[j]: exact(solution12_zero_raw[j])
        for j in range(len(free11_zero))
    }
    literal12_zero = [
        reduce_mod_univariate(value.subs(select12_zero_raw), fzero)
        for value in equations12_zero
    ]
    literal12_zero = [value for value in literal12_zero if value != 0]
    if len(literal12_zero) != 2:
        raise AssertionError("Phi-zero row-twelve literal residual support")
    c23_coefficients = [exact(sp.diff(value, c23)) for value in literal12_zero]
    if c23_coefficients != [
        sp.Rational(323, 360) * R8.eta,
        -sp.Rational(5, 7),
    ]:
        raise AssertionError(f"Phi-zero row-twelve c23 coefficients changed: {c23_coefficients}")
    a0, a1 = c23_coefficients
    b0 = exact(literal12_zero[0].subs(c23, 0))
    b1 = exact(literal12_zero[1].subs(c23, 0))
    incompatibility = reduce_mod_univariate(a0 * b1 - a1 * b0, fzero)
    incompatibility_numerator, _ = sp.fraction(exact(incompatibility))
    incompatibility_poly = sp.Poly(incompatibility_numerator, R8.eta, domain=sp.QQ).primitive()[1]
    if sp.gcd(incompatibility_poly, fzero).degree() != 0:
        raise AssertionError("Phi-zero row-twelve equations unexpectedly compatible")
    print(
        "PHI_ZERO_ROW12_EXTINCTION "
        "literal_c23_coefficients=(323*eta/360,-5/7) "
        f"compatibility={sp.sstr(incompatibility)} gcd(F,compatibility)=1 survivors=0"
    )

    alpha_answers = sp.solve(e9, R8.alpha11)
    if len(alpha_answers) != 1:
        raise AssertionError("generic Phi-nonzero row-nine alpha graph missing")
    alpha9 = exact(alpha_answers[0])
    print(f"ROW9_ALPHA_GENERIC={sp.sstr(alpha9)}")

    alpha_subs = {R8.alpha11: alpha9}
    a9 = impose_substitutions(a9_trial, alpha_subs, terminal9, alpha_subs)
    c9rows = impose_substitutions(c9_trial, alpha_subs, terminal9, alpha_subs)
    g10 = exact(grows[10].subs(bracket_subs).subs(gate).subs(alpha_subs))
    row10_select, row10_conditions, row10_matrix = solve_bracket_fibre(
        a9,
        c9rows,
        10,
        g10,
        theta9[:7],
        "ROW10_BRACKET",
    )
    print(f"ROW10_BRACKET_RANK={row10_matrix.rank()}")
    print(f"ROW10_CONDITION_COUNT={len(row10_conditions)}")
    solve10 = sp.solve(
        row10_conditions,
        (R8.xi10, R8.beta11),
        dict=True,
        simplify=False,
    )
    if len(solve10) != 1:
        raise AssertionError(f"row-ten response expected one graph, got {len(solve10)}")
    response10 = {symbol: exact(value) for symbol, value in solve10[0].items()}
    for symbol in (R8.xi10, R8.beta11):
        print(f"ROW10_{symbol}={sp.sstr(response10[symbol])}")

    a9_selected = impose_substitutions(a9, row10_select, response10)
    c9_selected = impose_substitutions(c9rows, row10_select, response10)
    a10_trial, c10_trial, theta10 = append_tangent(
        a9_selected,
        c9_selected,
        10,
        "w13_theta10",
    )
    depth10, amatrix10,chyd = projected_depth_residuals(
        a10_trial,
        c10_trial,
        10,
    )
    terminal10, depth10_conditions, _ = eliminate_linear_fibre(
        depth10,
        theta10,
        "ROW10_DEPTH",
    )
    print(
        "ROW10_DEPTH_UNIVERSE "
        f"P2={amatrix10.shape}/{amatrix10.rank()} "
        f"P3={chyd.shape}/{chyd.rank()}"
    )
    print(f"ROW10_DEPTH_CONDITION_COUNT={len(depth10_conditions)}")
    if len(depth10_conditions) != 1:
        raise AssertionError("row-ten joint depth should have one condition")
    c51_answers = sp.solve(depth10_conditions[0], c51)
    if len(c51_answers) != 1:
        raise AssertionError("row-ten transverse c51 graph missing")
    c51_graph = exact(c51_answers[0])
    print(f"ROW10_c51_TRANSVERSE_GRAPH={sp.sstr(c51_graph)}")
    old_d = (
        7231154026500 * R8.Phi**3
        + 50541940696500 * R8.Phi**2 * R8.eta
        + 6793915500000 * R8.Phi * R8.eta**2
        - 631918028977864704 * R8.Phi
        + 353642000625 * R8.eta**3
        - 91584545734393856 * R8.eta
    )
    if exact(c51_graph - old_d / (707284001250 * R8.Phi**2)) != 0:
        raise AssertionError("unexpected c51 transverse coefficient")

    c51_subs = {c51: c51_graph}
    a10 = impose_substitutions(
        a10_trial,
        c51_subs,
        terminal10,
        c51_subs,
    )
    c10rows = impose_substitutions(
        c10_trial,
        c51_subs,
        terminal10,
        c51_subs,
    )
    g11 = grows[11]
    for mapping in (bracket_subs, gate, alpha_subs, response10, c51_subs):
        g11 = exact(g11.subs(mapping))
    row11_select, row11_conditions, row11_matrix = solve_bracket_fibre(
        a10,
        c10rows,
        11,
        g11,
        theta10[:8],
        "ROW11_BRACKET",
        print_conditions=False,
    )
    print(f"ROW11_BRACKET_RANK={row11_matrix.rank()}")
    print(f"ROW11_CONDITION_COUNT={len(row11_conditions)}")
    if len(row11_conditions) != 1:
        raise AssertionError("row-eleven bracket should have one new condition")
    c23_answers = sp.solve(row11_conditions[0], c23)
    if len(c23_answers) != 1:
        raise AssertionError("row-eleven c23 graph missing")
    c23_graph = exact(c23_answers[0])
    print(f"ROW11_c23_GRAPH={sp.sstr(c23_graph)}")

    # Normalize the literal residual itself, rather than only the primitive
    # polynomial above, to expose the coefficient requested by the hostile.
    row11_difference = exact(g11 - R8.predicted_G(11, a10, c10rows))
    row11_equations = [xcoeff(row11_difference, degree) for degree in range(12)]
    matrix11_raw, rhs11_raw = sp.linear_eq_to_matrix(row11_equations, theta10[:8])
    pivots11_raw = tuple(matrix11_raw.T.rref()[1])
    solution11_raw = matrix11_raw.extract(pivots11_raw, tuple(range(8))).inv() * rhs11_raw.extract(pivots11_raw, (0,))
    select11_raw = {theta10[j]: exact(solution11_raw[j]) for j in range(8)}
    raw11 = [exact(value.subs(select11_raw)) for value in row11_equations]
    raw11_nonzero = [(index, value) for index, value in enumerate(raw11) if value != 0]
    if len(raw11_nonzero) != 1:
        raise AssertionError("generic row-eleven literal residual support")
    if exact(sp.diff(raw11_nonzero[0][1], c23) - sp.Rational(323, 504) * R8.Phi) != 0:
        raise AssertionError("generic row-eleven c23 coefficient changed")
    print(
        "ROW11_LITERAL_RESIDUAL "
        f"x_degree={raw11_nonzero[0][0]} c23_coefficient=323*Phi/504"
    )

    c23_subs = {c23: c23_graph}
    a10_selected = impose_substitutions(a10, row11_select, c23_subs)
    c10_selected = impose_substitutions(c10rows, row11_select, c23_subs)
    a11_trial, c11_trial, theta11 = append_tangent(
        a10_selected,
        c10_selected,
        11,
        "w13_theta11",
    )
    depth11, amatrix11, cmatrix11 = projected_depth_residuals(
        a11_trial,
        c11_trial,
        11,
    )
    terminal11, depth11_conditions, _ = eliminate_linear_fibre(
        depth11,
        theta11,
        "ROW11_DEPTH",
    )
    print(
        "ROW11_DEPTH_UNIVERSE "
        f"P2={amatrix11.shape}/{amatrix11.rank()} "
        f"P3={cmatrix11.shape}/{cmatrix11.rank()}"
    )
    print(f"ROW11_DEPTH_CONDITION_COUNT={len(depth11_conditions)}")

    a11 = impose_substitutions(a11_trial, terminal11)
    c11rows = impose_substitutions(c11_trial, terminal11)
    g12 = grows[12]
    for mapping in (
        bracket_subs,
        gate,
        alpha_subs,
        response10,
        c51_subs,
        c23_subs,
    ):
        g12 = exact(g12.subs(mapping))
    row12_select, row12_conditions, row12_matrix = solve_bracket_fibre(
        a11,
        c11rows,
        12,
        g12,
        theta11[:8],
        "ROW12_BRACKET",
        print_conditions=False,
    )
    print(f"ROW12_BRACKET_RANK={row12_matrix.rank()}")
    print(f"ROW12_CONDITION_COUNT={len(row12_conditions)}")
    if len(row12_conditions) != 2:
        raise AssertionError("row-twelve should leave two conditions")
    qodd = sign_quotient_polynomial(row12_conditions[0], 1)
    qeven = sign_quotient_polynomial(row12_conditions[1], 0)
    n12 = (
        272008125 * R8.Phi**2
        - 43740000 * R8.Phi * R8.eta
        + 10213932924928
    )
    nq = sp.Poly(
        (272008125 - 43740000 * ratio) * Y + 10213932924928,
        Y,
        ratio,
        domain=sp.QQ,
    )
    nq_slope = sp.Poly(272008125 - 43740000 * ratio, ratio, domain=sp.QQ)
    nq_constant = sp.Integer(10213932924928)
    print(f"ROW12_QUOTIENT_ODD={sp.sstr(qodd.as_expr())}")
    print(f"ROW12_QUOTIENT_EVEN={sp.sstr(qeven.as_expr())}")
    print(
        "ROW12_QUOTIENT_BIDEGREES "
        f"odd=({qodd.degree(Y)},{qodd.degree(ratio)}) "
        f"even=({qeven.degree(Y)},{qeven.degree(ratio)})"
    )

    odd_lc = sp.Poly(
        sp.expand(qodd.as_expr()).coeff(Y, qodd.degree(Y)),
        ratio,
        domain=sp.QQ,
    )
    even_lc = sp.Poly(
        sp.expand(qeven.as_expr()).coeff(Y, qeven.degree(Y)),
        ratio,
        domain=sp.QQ,
    )
    odd_at_zero = sp.Poly(qodd.eval(Y, 0), ratio, domain=sp.QQ)
    even_at_zero = sp.Poly(qeven.eval(Y, 0), ratio, domain=sp.QQ)
    print(
        "ROW12_BOUNDARY_GCDS "
        f"Y_infinity_common_leading={sp.sstr(sp.gcd(odd_lc, even_lc).primitive()[1].as_expr())} "
        f"Y_zero={sp.gcd(odd_at_zero, even_at_zero).degree()}"
    )
    if sp.gcd(odd_lc, even_lc).monic() != sp.Poly(9 * ratio + 8, ratio).monic():
        raise AssertionError("row-twelve common infinity factor changed")
    resultant = sp.Poly(
        sp.resultant(qodd.as_expr(), qeven.as_expr(), Y),
        ratio,
        domain=sp.QQ,
    ).primitive()[1]
    resultant_factor = sp.factor_list(resultant.as_expr())
    print(
        "ROW12_RESULTANT "
        f"degree={resultant.degree()} squarefree={sp.gcd(resultant, resultant.diff()).degree() == 0} "
        f"factor_degrees={[sp.Poly(factor, ratio).degree() for factor, _ in resultant_factor[1]]} "
        f"multiplicities={[power for _, power in resultant_factor[1]]}"
    )
    groebner12 = sp.groebner(
        [qodd.as_expr(), qeven.as_expr()],
        Y,
        ratio,
        order="lex",
        domain=sp.QQ,
    )
    print(
        "ROW12_GROEBNER "
        f"zero_dimensional={groebner12.is_zero_dimensional} "
        f"basis_degrees={[(sp.Poly(poly.as_expr(), Y, ratio).degree(Y), sp.Poly(poly.as_expr(), Y, ratio).degree(ratio)) for poly in groebner12.polys]}"
    )
    if not groebner12.is_zero_dimensional or len(groebner12.polys) != 2:
        raise AssertionError("unexpected row-twelve quotient basis")
    linear12 = sp.Poly(groebner12.polys[0].as_expr(), Y, ratio, domain=sp.QQ)
    k13 = sp.Poly(groebner12.polys[1].as_expr(), ratio, domain=sp.QQ)
    if linear12.degree(Y) != 1 or k13.degree() != 13:
        raise AssertionError("row-twelve graph/septimic replacement shape")
    linear_expr = linear12.as_expr()
    graph_a = sp.Poly(sp.expand(linear_expr).coeff(Y, 1), ratio, domain=sp.QQ)
    graph_b = sp.Poly(sp.expand(linear_expr).subs(Y, 0), ratio, domain=sp.QQ)
    k13_primitive = k13.primitive()[1]
    checks = {
        "K_squarefree": sp.gcd(k13_primitive, k13_primitive.diff()).degree() == 0,
        "Y_graph_denominator": sp.gcd(k13_primitive, graph_a).degree() == 0,
        "Y_graph_nonzero": sp.gcd(k13_primitive, graph_b).degree() == 0,
    }
    if not all(checks.values()):
        raise AssertionError(f"row-twelve reduced/nonzero audit failed: {checks}")
    factors = sp.factor_list(k13_primitive.as_expr())[1]
    if [(sp.Poly(factor, ratio).degree(), power) for factor, power in factors] != [(13, 1)]:
        raise AssertionError("K13 should be irreducible over Q")
    linear_factor = sp.Poly(resultant_factor[1][0][0], ratio, domain=sp.QQ).primitive()[1]
    if linear_factor.monic() != sp.Poly(9 * ratio + 8, ratio, domain=sp.QQ).monic():
        raise AssertionError("unexpected extraneous resultant factor")
    special_odd = sp.Poly(qodd.as_expr().subs(ratio, -sp.Rational(8, 9)), Y, domain=sp.QQ)
    special_even = sp.Poly(qeven.as_expr().subs(ratio, -sp.Rational(8, 9)), Y, domain=sp.QQ)
    special_gcd = sp.gcd(special_odd, special_even)
    if special_gcd.degree() != 0:
        raise AssertionError("extraneous ratio unexpectedly lies on affine carrier")
    print(
        "ROW12_REDUCED_CARRIER "
        "K13_degree=13 irreducible=yes squarefree=yes Y_graph_unique=yes "
        "Y_nonzero=yes quotient_points=13 signed_source_points=26 "
        "extraneous_resultant_ratio=-8/9 specialized_gcd_Y=1"
    )

    finite_control: tuple[int, int, int, int] | None = None
    for prime in (17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71):
        try:
            for rvalue in range(prime):
                if mod_value(k13_primitive.as_expr(), {ratio: rvalue}, prime) != 0:
                    continue
                avalue = mod_value(graph_a.as_expr(), {ratio: rvalue}, prime)
                bvalue = mod_value(graph_b.as_expr(), {ratio: rvalue}, prime)
                if avalue == 0:
                    continue
                yvalue = -bvalue * pow(avalue, -1, prime) % prime
                if yvalue == 0:
                    continue
                pvalue = next((candidate for candidate in range(1, prime) if candidate * candidate % prime == yvalue), None)
                if pvalue is None:
                    continue
                finite_control = (prime, rvalue, yvalue, pvalue)
                break
        except ZeroDivisionError:
            continue
        if finite_control is not None:
            break
    if finite_control is None:
        raise AssertionError("no finite-field signed positive control found")
    prime, rvalue, yvalue, pvalue = finite_control
    if mod_value(qodd.as_expr(), {ratio: rvalue, Y: yvalue}, prime) != 0:
        raise AssertionError("finite control misses odd quotient equation")
    if mod_value(qeven.as_expr(), {ratio: rvalue, Y: yvalue}, prime) != 0:
        raise AssertionError("finite control misses even quotient equation")
    finite_phi_eta = {R8.Phi: pvalue, R8.eta: rvalue * pvalue % prime}
    finite_c51 = mod_value(c51_graph, finite_phi_eta, prime)
    finite_c23 = mod_value(
        c23_graph.subs(c51_subs),
        finite_phi_eta,
        prime,
    )
    finite_nq = mod_value(
        nq.as_expr(),
        {ratio: rvalue, Y: yvalue},
        prime,
    )
    if finite_nq == 0:
        raise AssertionError("finite bracket control should be hostile to row-twelve depth")
    print(
        "ROW12_FINITE_CONTROL "
        f"F{prime} ratio={rvalue} Y={yvalue} Phi={pvalue} eta={rvalue * pvalue % prime} "
        f"c51={finite_c51} c23={finite_c23} both_quotient_equations=yes "
        f"sign_sheet_nonzero=yes depth_Nq={finite_nq}_nonzero"
    )

    # Continue the exact terminal fibre over the reduced 13-point quotient.
    a11_selected = impose_substitutions(a11, row12_select)
    c11_selected = impose_substitutions(c11rows, row12_select)
    a12_trial, c12_trial, theta12 = append_tangent(
        a11_selected,
        c11_selected,
        12,
        "w13_theta12",
    )
    depth12, amatrix12, cmatrix12 = projected_depth_residuals(
        a12_trial,
        c12_trial,
        12,
    )
    terminal12, depth12_conditions, depth12_matrix = eliminate_linear_fibre(
        depth12,
        theta12,
        "ROW12_DEPTH",
        print_conditions=False,
    )
    depth12_remainders = [value for value in depth12_conditions if not quotient_vanishes(value, groebner12)]
    print(
        "ROW12_DEPTH_UNIVERSE "
        f"P2={amatrix12.shape}/{amatrix12.rank()} "
        f"P3={cmatrix12.shape}/{cmatrix12.rank()} "
        f"terminal_rank={depth12_matrix.rank()} pivots={depth12_matrix.rref()[1]} "
        f"raw_conditions={len(depth12_conditions)} nonzero_mod_carrier={len(depth12_remainders)}"
    )

    if depth12_remainders:
        print("ROW12_DEPTH_CARRIER_STATUS=extra_conditions")
        for index, value in enumerate(depth12_remainders):
            even_component, odd_component = quotient_components(value)
            even_remainder = sp.factor(groebner12.reduce(even_component)[1])
            odd_remainder = sp.factor(groebner12.reduce(odd_component)[1])
            print(
                f"ROW12_DEPTH_REMAINDER_{index} "
                f"even_degree={sp.Poly(even_remainder, ratio).degree()} "
                f"even_nonzero={even_remainder != 0} odd_nonzero={odd_remainder != 0}"
            )
    else:
        print("ROW12_DEPTH_CARRIER_STATUS=automatic")

    raw_depth12 = [exact(value.subs(terminal12)) for value in depth12]
    raw_depth12_nonzero = [(index, value) for index, value in enumerate(raw_depth12) if value != 0]
    if len(raw_depth12_nonzero) != 1:
        raise AssertionError("row-twelve depth should expose one literal mismatch")
    raw_depth_index, raw_depth_value = raw_depth12_nonzero[0]
    if exact(raw_depth_value + n12 / 288360000) != 0:
        raise AssertionError("row-twelve joint-system literal depth mismatch constant")
    print(
        "ROW12_DEPTH_LITERAL "
        f"joint_residual_index={raw_depth_index} value=-({sp.sstr(n12)})/288360000"
    )

    aleft12 = len(amatrix12.T.nullspace())
    ares12 = depth12[:aleft12]
    cres12 = depth12[aleft12:]
    amat12, arhs12 = sp.linear_eq_to_matrix(ares12, theta12)
    cmat12, crhs12 = sp.linear_eq_to_matrix(cres12, theta12)
    apivots12 = tuple(amat12.rref()[1])
    cpivots12 = tuple(cmat12.rref()[1])
    asolution12 = selected_solution(amat12, arhs12, apivots12)
    csolution12 = selected_solution(cmat12, crhs12, cpivots12)
    amap12 = {column: asolution12[index] for index, column in enumerate(apivots12)}
    cmap12 = {column: csolution12[index] for index, column in enumerate(cpivots12)}
    shared12 = sorted(set(apivots12) & set(cpivots12))
    separate_mismatches = {
        column: exact(amap12[column] - cmap12[column])
        for column in shared12
        if exact(amap12[column] - cmap12[column]) != 0
    }
    if list(separate_mismatches.values()) != [n12 / 108135000]:
        raise AssertionError(f"unexpected separate-depth mismatch: {separate_mismatches}")
    print(
        "ROW12_DEPTH_SEPARATE "
        f"P2_rank={amat12.rank()} pivots={apivots12} "
        f"P3_rank={cmat12.rank()} pivots={cpivots12} "
        f"shared_mismatch=({sp.sstr(n12)})/108135000"
    )

    substituted_odd = sp.Poly(
        sp.together(
            qodd.as_expr().subs(Y, -nq_constant / nq_slope.as_expr())
            * nq_slope.as_expr() ** 3
        ),
        ratio,
        domain=sp.QQ,
    ).primitive()[1]
    substituted_even = sp.Poly(
        sp.together(
            qeven.as_expr().subs(Y, -nq_constant / nq_slope.as_expr())
            * nq_slope.as_expr() ** 3
        ),
        ratio,
        domain=sp.QQ,
    ).primitive()[1]
    if (substituted_odd.degree(), substituted_even.degree()) != (6, 5):
        raise AssertionError("unexpected row-twelve depth elimination degrees")
    if sp.gcd(substituted_odd, substituted_even).degree() != 0:
        raise AssertionError("row-twelve depth should be disjoint from bracket carrier")
    f0_mod11 = sp.Poly(substituted_odd.as_expr(), ratio, modulus=11).monic()
    f1_mod11 = sp.Poly(substituted_even.as_expr(), ratio, modulus=11).monic()
    expected_f0_mod11 = sp.Poly(
        3 * ratio**6 - 4 * ratio**5 - 2 * ratio**2 - 2 * ratio - 2,
        ratio,
        modulus=11,
    ).monic()
    expected_f1_mod11 = sp.Poly(
        5 * ratio**5 - 5 * ratio**3 + 5 * ratio**2 + 2 * ratio + 5,
        ratio,
        modulus=11,
    ).monic()
    if f0_mod11 != expected_f0_mod11 or f1_mod11 != expected_f1_mod11:
        raise AssertionError("mod-11 depth eliminants changed")
    bezout_s = 2 * ratio**4 + 4 * ratio**3 + 3 * ratio**2 - 5 * ratio + 2
    bezout_t = ratio**5 - 3 * ratio**4 - 2 * ratio**3 - 3 * ratio**2 + 5 * ratio + 1
    # Use the unnormalized reductions displayed by the compact certificate.
    raw_f0_mod11 = sp.Poly(
        3 * ratio**6 - 4 * ratio**5 - 2 * ratio**2 - 2 * ratio - 2,
        ratio,
        modulus=11,
    )
    raw_f1_mod11 = sp.Poly(
        5 * ratio**5 - 5 * ratio**3 + 5 * ratio**2 + 2 * ratio + 5,
        ratio,
        modulus=11,
    )
    bezout_check = sp.Poly(
        bezout_s * raw_f0_mod11.as_expr()
        + bezout_t * raw_f1_mod11.as_expr()
        - 1,
        ratio,
        modulus=11,
    )
    if not bezout_check.is_zero:
        raise AssertionError("mod-11 depth Bezout changed")
    print(
        "ROW12_DEPTH_EXTINCTION "
        "Nq=(272008125-43740000*ratio)*Y+10213932924928; "
        "Nq_slope_zero_implies_constant_nonzero; substituted_degrees=(6,5); "
        "mod11_Bezout=1; Phi_nonzero_survivors=0"
    )

    print("WEIGHT13_FACE_RESULT=EXTINCT_BY_ROW12_PROJECTED_DEPTH")
    print(
        "MECHANISM Phi_nonzero: c51 absorbs old elliptic normal at row10; "
        "c23 absorbs row11 bracket; 26 row12 bracket points are depth-disjoint"
    )
    print(
        "MECHANISM Phi_zero: c51 absorbs row11 bracket; two row12 c23 equations "
        "are incompatible on the reduced quadratic boundary"
    )
    print(
        "SCOPE complete fixed THM-4308 source-normal residual-weight<=13 face; "
        "finite projected P2/P3 rows only; no chart/seam entry, all-weight lift, "
        "Keller pair, JC2, or DC2"
    )
    print("RESULT=PASS")

    # Only a depth-compatible carrier can lawfully be advanced to row 13.
    if not depth12_remainders:
        a12 = impose_substitutions(a12_trial, terminal12)
        c12rows = impose_substitutions(c12_trial, terminal12)
        g13 = grows[13]
        for mapping in (
            bracket_subs,
            gate,
            alpha_subs,
            response10,
            c51_subs,
            c23_subs,
        ):
            g13 = exact(g13.subs(mapping))
        free12 = tuple(
            symbol
            for index, symbol in enumerate(theta12)
            if index not in depth12_matrix.rref()[1]
        )
        _, row13_conditions, row13_matrix = solve_bracket_fibre(
            a12,
            c12rows,
            13,
            g13,
            free12,
            "ROW13_BRACKET",
            # type: ignore[call-arg]
        )
        row13_remainders = [value for value in row13_conditions if not quotient_vanishes(value, groebner12)]
        print(
            "ROW13_BRACKET_CARRIER "
            f"matrix={row13_matrix.shape}/{row13_matrix.rank()} free={len(free12)} "
            f"raw_conditions={len(row13_conditions)} nonzero_mod_carrier={len(row13_remainders)}"
        )


if __name__ == "__main__":
    main()
