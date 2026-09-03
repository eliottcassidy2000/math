#!/usr/bin/env python3
"""Primary exact certificate for THM-4380.

This continues the complete THM-4308 source-normal residual-weight-at-most-
twelve family.  It proves that the full row-nine bracket hypersurface has
automatic projected P_2/P_3 depth, treats Phi=0 and Phi!=0 separately, and
shows that every surviving finite jet dies by the row-twelve bracket.  The
calculation is confined to the declared finite source-normal chart; it says
nothing about chart/seam entry or JC(2).
"""

from __future__ import annotations

from math import gcd
from pathlib import Path
import sys

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")

sys.path.insert(0, str(Path(__file__).resolve().parent))
import jc2_source_normal_bracket_hasse_rows8_thm4308 as R8  # noqa: E402
import jc2_source_normal_student_stein_row9_thm4315 as R9  # noqa: E402


CHECKS = 0
x = R8.x
P = R8.Phi
eta = R8.eta
xi = R8.xi10
alpha = R8.alpha11
beta = R8.beta11
Y, r = sp.symbols("Y r")


def check(condition: bool, label: str) -> None:
    global CHECKS
    if not condition:
        raise AssertionError(label)
    CHECKS += 1


def exact(expression: sp.Expr) -> sp.Expr:
    return sp.cancel(sp.together(expression))


def is_zero(expression: sp.Expr) -> bool:
    return exact(sp.expand(expression)) == 0


def check_zero(expression: sp.Expr, label: str) -> None:
    check(is_zero(expression), label)


def selected_solution(
    matrix: sp.Matrix,
    rhs: sp.Matrix,
    columns: tuple[int, ...],
) -> tuple[sp.Expr, ...]:
    """Solve the inherited lexicographic full-rank square subsystem."""

    row_indices = tuple(matrix.T.rref()[1])
    square = matrix.extract(row_indices, columns)
    check(square.rows == square.cols, "selected block square")
    check(square.det() != 0, "selected block invertible")
    answer = square.inv() * rhs.extract(row_indices, (0,))
    return tuple(exact(value) for value in answer)


def primitive_numerator(expression: sp.Expr) -> sp.Poly:
    numerator, _ = sp.fraction(exact(expression))
    return sp.Poly(numerator, P, eta, domain=sp.QQ).primitive()[1]


def reduce_univariate(expression: sp.Expr, variable: sp.Symbol, modulus: sp.Poly) -> sp.Expr:
    value = exact(expression)
    numerator, denominator = sp.fraction(value)
    mod = sp.Poly(modulus.as_expr(), variable, domain=sp.EX)
    num = sp.Poly(numerator, variable, domain=sp.EX)
    den = sp.Poly(denominator, variable, domain=sp.EX)
    inverse = sp.invert(den, mod)
    return sp.factor((num * inverse).rem(mod).as_expr())


def reduce_ideal_numerator(expression: sp.Expr, basis: sp.GroebnerBasis) -> sp.Expr:
    numerator, denominator = sp.fraction(exact(expression))
    den_poly = sp.Poly(denominator, P, eta, domain=sp.QQ)
    # Every rational chart formula used below is localized only at Phi.
    check(all(term[0][1] == 0 for term in den_poly.terms()), "denominator eta-free")
    check(len(den_poly.terms()) == 1, "denominator is a Phi monomial")
    return sp.factor(basis.reduce(numerator)[1])


def modular_value(expression: sp.Expr, substitutions: dict[sp.Symbol, int], prime: int) -> int:
    value = exact(expression.subs(substitutions))
    check(not value.free_symbols, "modular value scalar")
    numerator, denominator = sp.fraction(value)
    den = int(denominator) % prime
    check(den != 0, "modular denominator nonzero")
    return int(numerator) % prime * pow(den, -1, prime) % prime


def append_tangent(
    arows: list[sp.Expr],
    crows: list[sp.Expr],
    row: int,
    prefix: str,
) -> tuple[list[sp.Expr], list[sp.Expr], tuple[sp.Symbol, ...]]:
    bvalue = R8.B_row(row, arows, crows)
    abase, cbase = R8.particular_row(row, bvalue)
    symbols = tuple(sp.symbols(f"{prefix}0:{row + 1}"))
    theta = sum(symbols[j] * x**j for j in range(row + 1))
    anew = sp.expand(abase + theta * sp.diff(R8.A0, x))
    cnew = sp.expand(cbase + theta * sp.diff(R8.C0, x))
    return arows + [anew], crows + [cnew], symbols


def depth_data(
    arows: list[sp.Expr],
    crows: list[sp.Expr],
    row: int,
    symbols: tuple[sp.Symbol, ...],
) -> tuple[
    list[sp.Expr],
    list[sp.Expr],
    sp.Matrix,
    sp.Matrix,
    sp.Matrix,
    sp.Matrix,
    tuple[list[tuple[int, int]], sp.Matrix],
    tuple[list[tuple[int, int]], sp.Matrix],
]:
    acoords, amatrix = R9.depth_matrix(2, row)
    ccoords, cmatrix = R9.depth_matrix(3, row)
    avec = sp.Matrix([R8.xcoeff(arows[n], degree) for n, degree in acoords])
    cvec = sp.Matrix([R8.xcoeff(crows[n], degree) for n, degree in ccoords])
    aresiduals = [sp.expand((left.T * avec)[0]) for left in amatrix.T.nullspace()]
    cresiduals = [sp.expand((left.T * cvec)[0]) for left in cmatrix.T.nullspace()]
    amat, arhs = sp.linear_eq_to_matrix(aresiduals, symbols)
    cmat, crhs = sp.linear_eq_to_matrix(cresiduals, symbols)
    return (
        aresiduals,
        cresiduals,
        amat,
        arhs,
        cmat,
        crhs,
        (acoords, amatrix),
        (ccoords, cmatrix),
    )


def bracket_select(
    arows: list[sp.Expr],
    crows: list[sp.Expr],
    row: int,
    source_row: sp.Expr,
    free_symbols: tuple[sp.Symbol, ...],
) -> tuple[list[sp.Expr], sp.Matrix, tuple[sp.Expr, ...], dict[sp.Symbol, sp.Expr]]:
    difference = exact(source_row - R8.predicted_G(row, arows, crows))
    # In the inherited degree-capped source chart both sides have x-degree at
    # most ``row``.  Thus coefficient comparison in degrees 0,...,row is the
    # entire bracket equation, not a selected projection of it.
    check(
        difference == 0 or sp.Poly(difference, x, domain=sp.EX).degree() <= row,
        f"row-{row} bracket residual exhausts degree cap",
    )
    equations = [
        exact(R8.xcoeff(difference, degree))
        for degree in range(row + 1)
    ]
    matrix, rhs = sp.linear_eq_to_matrix(equations, free_symbols)
    solution = selected_solution(matrix, rhs, tuple(range(len(free_symbols))))
    substitutions = {symbol: solution[j] for j, symbol in enumerate(free_symbols)}
    residuals = [sp.factor(exact(value.subs(substitutions))) for value in equations]
    return residuals, matrix, solution, substitutions


def invariant_even(poly: sp.Poly) -> sp.Poly:
    """Replace eta=r*Phi in an even-total-degree polynomial and Phi^2 by Y."""

    answer = sp.Integer(0)
    for (pdegree, edegree), coefficient in poly.terms():
        total = pdegree + edegree
        check(total % 2 == 0, "even invariant total degree")
        answer += coefficient * Y ** (total // 2) * r**edegree
    return sp.Poly(answer, Y, r, domain=sp.QQ)


def invariant_odd_div_phi(poly: sp.Poly) -> sp.Poly:
    """Replace eta=r*Phi after dividing an odd-total-degree polynomial by Phi."""

    answer = sp.Integer(0)
    for (pdegree, edegree), coefficient in poly.terms():
        total = pdegree + edegree
        check(total % 2 == 1, "odd invariant total degree")
        answer += coefficient * Y ** ((total - 1) // 2) * r**edegree
    return sp.Poly(answer, Y, r, domain=sp.QQ)


def build_row_eight() -> tuple[
    list[sp.Expr],
    list[sp.Expr],
    dict[sp.Symbol, sp.Expr],
    dict[sp.Symbol, sp.Expr],
    list[sp.Symbol],
    list[sp.Expr],
    list[sp.Expr],
    sp.Expr,
]:
    arows, crows, bracket_subs = R8.build_bracket_rows()
    amatrix8, cmatrix8, theta8, terminal8 = R8.hasse_audit(arows, crows, bracket_subs)
    check((amatrix8.shape, amatrix8.rank()) == ((63, 131), 51), "row-eight P2 inherited")
    check((cmatrix8.shape, cmatrix8.rank()) == ((72, 204), 63), "row-eight P3 inherited")
    g9 = sp.expand(R8.tcoeff(R8.H, 9))
    fixed_a, fixed_c, selector, _ = R9.solve_row_eight_for_g9(
        arows,
        crows,
        bracket_subs,
        terminal8,
        theta8,
        g9,
    )
    check(selector.rank() == 7, "row-nine selects old affine-seven fibre")
    return arows, crows, bracket_subs, terminal8, theta8, fixed_a, fixed_c, g9


def row_nine_depth(
    fixed_a: list[sp.Expr],
    fixed_c: list[sp.Expr],
    bracket_subs: dict[sp.Symbol, sp.Expr],
    terminal8: dict[sp.Symbol, sp.Expr],
    g9: sp.Expr,
    source_subs: dict[sp.Symbol, sp.Expr],
    tag: str,
) -> tuple[list[sp.Expr], list[sp.Expr], tuple[sp.Symbol, ...], dict[sp.Symbol, sp.Expr]]:
    arows = [exact(row.subs(bracket_subs).subs(terminal8).subs(source_subs)) for row in fixed_a]
    crows = [exact(row.subs(bracket_subs).subs(terminal8).subs(source_subs)) for row in fixed_c]
    source9 = exact(g9.subs(bracket_subs).subs(terminal8).subs(source_subs))
    check_zero(R8.predicted_G(9, arows, crows) - source9, f"{tag} row-nine bracket")
    arows, crows, symbols = append_tangent(arows, crows, 9, f"{tag}q")
    ares, cres, amat, arhs, cmat, crhs, adata, cdata = depth_data(
        arows,
        crows,
        9,
        symbols,
    )
    check((len(adata[0]), adata[1].cols, adata[1].rank(), len(adata[1].T.nullspace())) == (75, 160, 59, 16), f"{tag} row-nine P2 universe")
    check((len(cdata[0]), cdata[1].cols, cdata[1].rank(), len(cdata[1].T.nullspace())) == (85, 251, 73, 12), f"{tag} row-nine P3 universe")
    joint = amat.col_join(cmat)
    joint_rhs = arhs.col_join(crhs)
    check(amat.rank() == 3 and cmat.rank() == 2 and joint.rank() == 3, f"{tag} row-nine terminal ranks")
    check(amat.rref()[1] == (7, 8, 9) and cmat.rref()[1] == (8, 9), f"{tag} row-nine separate pivots")
    check(joint.rref()[1] == (7, 8, 9), f"{tag} row-nine pivots")
    check(joint[:, :7] == sp.zeros(joint.rows, 7), f"{tag} row-nine free columns")
    solution = selected_solution(joint, joint_rhs, (7, 8, 9))
    substitutions = {symbols[7 + j]: solution[j] for j in range(3)}
    check(all(is_zero(value.subs(substitutions)) for value in ares + cres), f"{tag} row-nine depth automatic")
    arows[-1] = exact(arows[-1].subs(substitutions))
    crows[-1] = exact(crows[-1].subs(substitutions))
    return arows, crows, symbols, substitutions


def phi_zero_branch(
    fixed_a: list[sp.Expr],
    fixed_c: list[sp.Expr],
    bracket_subs: dict[sp.Symbol, sp.Expr],
    terminal8: dict[sp.Symbol, sp.Expr],
    g9: sp.Expr,
) -> None:
    xi0 = sp.factor(sp.solve(R9.E9_GATE.subs(P, 0), xi)[0])
    expected_xi0 = (255605625 * eta**2 + 46483785515008) / 6736896000
    check_zero(xi0 - expected_xi0, "Phi-zero E9 xi graph")
    source9_subs = {P: 0, xi: xi0}
    a9, c9, theta9, _ = row_nine_depth(
        fixed_a,
        fixed_c,
        bracket_subs,
        terminal8,
        g9,
        source9_subs,
        "z",
    )

    g10 = exact(R8.tcoeff(R8.H, 10).subs(bracket_subs).subs(terminal8).subs(source9_subs))
    residuals10, matrix10, _, select10 = bracket_select(a9, c9, 10, g10, theta9[:7])
    check(matrix10.rank() == 7 and matrix10.T.rref()[1] == (0, 1, 2, 3, 4, 5, 7), "Phi-zero row-ten selector")
    expected6 = (
        278844730740000 * alpha * eta
        + 854861882349375 * eta**2
        + 27743789318253707264
    ) / 186472018080000
    expected8 = 13 * (18612736875 * eta**2 - 4820239249178624) / 2302123680000
    check_zero(residuals10[6] - expected6, "Phi-zero row-ten residual six")
    check_zero(residuals10[8] - expected8, "Phi-zero row-ten residual eight")
    check(all(value == 0 or index in (6, 8) for index, value in enumerate(residuals10)), "Phi-zero row-ten residual support")

    f0 = sp.Poly(18612736875 * eta**2 - 4820239249178624, eta, domain=sp.QQ)
    check(f0.degree() == 2 and f0.TC() != 0 and sp.gcd(f0, f0.diff()) == sp.Poly(1, eta, domain=sp.QQ), "Phi-zero two reduced eta roots")
    check(sp.gcd(f0, sp.Poly(eta, eta, domain=sp.QQ)) == sp.Poly(1, eta, domain=sp.QQ), "Phi-zero eta roots nonzero")
    alpha0 = sp.factor(sp.solve(expected6, alpha)[0])
    expected_alpha0 = -(
        854861882349375 * eta**2 + 27743789318253707264
    ) / (278844730740000 * eta)
    check_zero(alpha0 - expected_alpha0, "Phi-zero alpha graph")
    source10_subs = {alpha: alpha0}
    check(reduce_univariate(expected6.subs(source10_subs), eta, f0) == 0, "Phi-zero row-ten residual six on F")
    check(reduce_univariate(expected8, eta, f0) == 0, "Phi-zero row-ten residual eight on F")
    a9 = [exact(row.subs(select10).subs(source10_subs)) for row in a9]
    c9 = [exact(row.subs(select10).subs(source10_subs)) for row in c9]

    a10, c10, theta10 = append_tangent(a9, c9, 10, "zr")
    ares, cres, amat, arhs, cmat, crhs, adata, cdata = depth_data(a10, c10, 10, theta10)
    check((len(adata[0]), adata[1].cols, adata[1].rank(), len(adata[1].T.nullspace())) == (88, 193, 68, 20), "Phi-zero row-ten P2 universe")
    check((len(cdata[0]), cdata[1].cols, cdata[1].rank(), len(cdata[1].T.nullspace())) == (99, 304, 83, 16), "Phi-zero row-ten P3 universe")
    check(amat.rank() == 3 and cmat.rank() == 3 and amat.col_join(cmat).rank() == 3, "Phi-zero row-ten depth ranks")
    check(amat.rref()[1] == (8, 9, 10) and cmat.rref()[1] == (8, 9, 10), "Phi-zero row-ten depth pivots")
    asolution = selected_solution(amat, arhs, (8, 9, 10))
    csolution = selected_solution(cmat, crhs, (8, 9, 10))
    check_zero(asolution[0] - csolution[0] + sp.Rational(4, 15) * (5 * beta + 6 * eta), "Phi-zero joint-depth beta equation")
    check_zero(asolution[1] - csolution[1], "Phi-zero shared coordinate nine")
    check_zero(asolution[2] - csolution[2], "Phi-zero shared coordinate ten")
    beta0 = -sp.Rational(6, 5) * eta
    asubs = {theta10[8 + j]: asolution[j] for j in range(3)}
    check(all(reduce_univariate(value.subs(asubs).subs(beta, beta0), eta, f0) == 0 for value in ares + cres), "Phi-zero row-ten joint depth on F")
    a10 = [exact(row.subs(asubs).subs(beta, beta0)) for row in a10]
    c10 = [exact(row.subs(asubs).subs(beta, beta0)) for row in c10]

    response_subs = {**source9_subs, **source10_subs, beta: beta0}
    response = {
        symbol: reduce_univariate(value.subs(terminal8).subs(response_subs), eta, f0)
        for symbol, value in bracket_subs.items()
    }
    check(response[R8.U] == -sp.Rational(54752794624, 8110125), "Phi-zero U nonzero")
    check(response[R8.W] == sp.Rational(6445195264, 180225), "Phi-zero W nonzero")
    check(response[R8.Z] == -sp.Rational(307951616, 11125), "Phi-zero Z nonzero")

    g11 = exact(R8.tcoeff(R8.H, 11).subs(bracket_subs).subs(terminal8).subs(response_subs))
    residuals11, matrix11, _, _ = bracket_select(a10, c10, 11, g11, theta10[:8])
    check(matrix11.rank() == 8 and matrix11.T.rref()[1] == tuple(range(8)), "Phi-zero row-eleven selector")
    reduced11 = [reduce_univariate(value, eta, f0) for value in residuals11]
    obstruction = sp.Rational(
        12353759608682044232849562614853632,
        137875982427855174564984375,
    )
    check(reduced11[8] == obstruction and obstruction != 0, "Phi-zero row-eleven obstruction")
    check(all(value == 0 or index == 8 for index, value in enumerate(reduced11)), "Phi-zero row-eleven residual support")

    print(f"PHI_ZERO_E9_XI={sp.sstr(xi0)}")
    print(f"PHI_ZERO_ROW10_F={sp.sstr(f0.as_expr())}")
    print(f"PHI_ZERO_ROW10_ALPHA={sp.sstr(alpha0)}")
    print(f"PHI_ZERO_ROW10_BETA={sp.sstr(beta0)}")
    print("PHI_ZERO_ROW10 source_points=2 terminal_fibre=A8 U,W,Z_nonzero=yes")
    print(f"PHI_ZERO_ROW11_OBSTRUCTION={obstruction}")


def phi_nonzero_branch(
    fixed_a: list[sp.Expr],
    fixed_c: list[sp.Expr],
    bracket_subs: dict[sp.Symbol, sp.Expr],
    terminal8: dict[sp.Symbol, sp.Expr],
    g9: sp.Expr,
) -> None:
    alpha9 = sp.factor(sp.solve(R9.E9_GATE, alpha)[0])
    expected_alpha9 = (
        613527750 * P**2
        - 3154140000 * P * eta
        - 255605625 * eta**2
        + 6736896000 * xi
        - 46483785515008
    ) / (511211250 * P)
    check_zero(alpha9 - expected_alpha9, "Phi-nonzero E9 alpha graph")
    source9_subs = {alpha: alpha9}
    a9, c9, theta9, depth9_subs = row_nine_depth(
        fixed_a,
        fixed_c,
        bracket_subs,
        terminal8,
        g9,
        source9_subs,
        "n",
    )
    expected_q7 = 2 * (
        439374375 * P**2
        + 240570000 * P * eta
        - 301968000 * xi
        - 8270057119744
    ) / 27064125
    expected_q8 = (
        1440859860 * P**2
        - 408969000 * P * beta
        + 2344381380 * P * eta
        + 255605625 * eta**2
        - 6736896000 * xi
        + 46483785515008
    ) / (61345350 * P)
    expected_q9 = -4 * (2025 * xi - 731648) / 4455
    check_zero(depth9_subs[theta9[7]] - expected_q7, "Phi-nonzero row-nine q7")
    check_zero(depth9_subs[theta9[8]] - expected_q8, "Phi-nonzero row-nine q8")
    check_zero(depth9_subs[theta9[9]] - expected_q9, "Phi-nonzero row-nine q9")

    g10 = exact(R8.tcoeff(R8.H, 10).subs(bracket_subs).subs(terminal8).subs(source9_subs))
    residuals10, matrix10, _, select10 = bracket_select(a9, c9, 10, g10, theta9[:7])
    check(matrix10.rank() == 7 and matrix10.T.rref()[1] == (0, 1, 2, 3, 4, 5, 7), "Phi-nonzero row-ten selector")
    expected6 = (
        722244600000 * P**3
        + 248448667500 * P**2 * beta
        - 12116672914500 * P**2 * eta
        - 1618073820000 * P * eta**2
        + 21190476096000 * P * xi
        - 120118953400336384 * P
        - 131125685625 * eta**3
        + 3456027648000 * eta * xi
        - 23846181969199104 * eta
    ) / (175375530000 * P)
    expected8 = 13 * (
        13365000 * P**2
        + 15035625 * P * eta
        + 57672000 * xi
        - 964604821504
    ) / 270641250
    check_zero(residuals10[6] - expected6, "Phi-nonzero row-ten residual six")
    check_zero(residuals10[8] - expected8, "Phi-nonzero row-ten residual eight")
    check(all(value == 0 or index in (6, 8) for index, value in enumerate(residuals10)), "Phi-nonzero row-ten residual support")
    xi10 = -(
        13365000 * P**2 + 15035625 * P * eta - 964604821504
    ) / 57672000
    beta10 = (
        11296175760000 * P**3
        + 49737870463500 * P**2 * eta
        + 6793915500000 * P * eta**2
        - 631918028977864704 * P
        + 353642000625 * eta**3
        - 91584545734393856 * eta
    ) / (670058527500 * P**2)
    source10_subs = {xi: xi10, beta: beta10}
    check_zero(expected8.subs(source10_subs), "Phi-nonzero xi graph")
    check_zero(expected6.subs(source10_subs), "Phi-nonzero beta graph")
    alpha10 = sp.factor(alpha9.subs(source10_subs))
    expected_alpha10 = -(
        69009144750 * P**2
        + 357574500000 * P * eta
        + 18612736875 * eta**2
        - 4820239249178624
    ) / (37225473750 * P)
    check_zero(alpha10 - expected_alpha10, "Phi-nonzero alpha row-ten graph")
    a9 = [exact(row.subs(select10).subs(source10_subs)) for row in a9]
    c9 = [exact(row.subs(select10).subs(source10_subs)) for row in c9]

    a10, c10, theta10 = append_tangent(a9, c9, 10, "nr")
    ares10, cres10, amat10, arhs10, cmat10, crhs10, adata10, cdata10 = depth_data(
        a10,
        c10,
        10,
        theta10,
    )
    check((len(adata10[0]), adata10[1].cols, adata10[1].rank(), len(adata10[1].T.nullspace())) == (88, 193, 68, 20), "Phi-nonzero row-ten P2 universe")
    check((len(cdata10[0]), cdata10[1].cols, cdata10[1].rank(), len(cdata10[1].T.nullspace())) == (99, 304, 83, 16), "Phi-nonzero row-ten P3 universe")
    check(amat10.rank() == 3 and cmat10.rank() == 3 and amat10.col_join(cmat10).rank() == 3, "Phi-nonzero row-ten depth ranks")
    check(amat10.rref()[1] == (8, 9, 10) and cmat10.rref()[1] == (8, 9, 10), "Phi-nonzero row-ten depth pivots")
    asolution10 = selected_solution(amat10, arhs10, (8, 9, 10))
    csolution10 = selected_solution(cmat10, crhs10, (8, 9, 10))
    dpoly = sp.Poly(
        7231154026500 * P**3
        + 50541940696500 * P**2 * eta
        + 6793915500000 * P * eta**2
        - 631918028977864704 * P
        + 353642000625 * eta**3
        - 91584545734393856 * eta,
        P,
        eta,
        domain=sp.QQ,
    )
    check_zero(asolution10[0] - csolution10[0] + dpoly.as_expr() / (502543895625 * P**2), "Phi-nonzero row-ten joint-depth curve")
    check_zero(asolution10[1] - csolution10[1], "Phi-nonzero row-ten coordinate nine")
    check_zero(asolution10[2] - csolution10[2], "Phi-nonzero row-ten coordinate ten")
    asubs10 = {theta10[8 + j]: asolution10[j] for j in range(3)}
    check(all(is_zero(value.subs(asubs10)) for value in ares10), "Phi-nonzero P2 standalone depth")
    csubs10 = {theta10[8 + j]: csolution10[j] for j in range(3)}
    check(all(is_zero(value.subs(csubs10)) for value in cres10), "Phi-nonzero P3 standalone depth")
    a10[-1] = exact(a10[-1].subs(asubs10))
    c10[-1] = exact(c10[-1].subs(asubs10))

    g11 = exact(
        R8.tcoeff(R8.H, 11)
        .subs(bracket_subs)
        .subs(terminal8)
        .subs(source9_subs)
        .subs(source10_subs)
    )
    residuals11, matrix11, _, select11 = bracket_select(a10, c10, 11, g11, theta10[:8])
    check(matrix11.rank() == 8 and matrix11.T.rref()[1] == tuple(range(8)), "Phi-nonzero row-eleven selector")
    epoly = sp.Poly(
        35898624648545625000000 * P**6
        + 80771905459227656250000 * P**5 * eta
        + 45434196820815556640625 * P**4 * eta**2
        - 5071037898252514202802000000 * P**4
        - 5410571375319483376728000000 * P**3 * eta
        + 106595566955131983930000000 * P**2 * eta**2
        + 179564695363704166433073974476800 * P**2
        + 16899675818606082900000000 * P * eta**3
        - 2249959222183587936424427520000 * P * eta
        + 851303752055234564062500 * eta**4
        - 268196544113786432750223360000 * eta**2
        + 12360863815079359770761606339756032,
        P,
        eta,
        domain=sp.QQ,
    )
    check_zero(residuals11[8] - epoly.as_expr() / (1971983693758919400000000 * P**2), "Phi-nonzero row-eleven E residual")
    check_zero(residuals11[9] + 2 * dpoly.as_expr() / (1842660950625 * P**2), "Phi-nonzero row-eleven D residual")
    check(all(value == 0 or index in (8, 9) for index, value in enumerate(residuals11)), "Phi-nonzero row-eleven residual support")

    apoly = sp.Poly(
        7231154026500
        + 50541940696500 * r
        + 6793915500000 * r**2
        + 353642000625 * r**3,
        r,
        domain=sp.QQ,
    )
    bpoly = sp.Poly(631918028977864704 + 91584545734393856 * r, r, domain=sp.QQ)
    check(invariant_odd_div_phi(dpoly) == sp.Poly(apoly.as_expr() * Y - bpoly.as_expr(), Y, r, domain=sp.QQ), "D ratio parametrization")
    ebar = invariant_even(epoly)
    kraw = sp.Poly(
        sp.together(ebar.as_expr().subs(Y, bpoly.as_expr() / apoly.as_expr()) * apoly.as_expr() ** 3),
        r,
        domain=sp.QQ,
    )
    kpoly = kraw.primitive()[1]
    expected_k = sp.Poly(
        21252176198679866250754006276839556755825 * r**7
        + 799311827675117522149997435401077574131600 * r**6
        + 14384863896403857958176347858723924433398460 * r**5
        + 142730433788981669223548142320603830110956220 * r**4
        + 786944231209420107856657052701244375708027892 * r**3
        + 1852564642916723756803328543705267257790149632 * r**2
        - 147733098192443646925107791876239203619548432 * r
        + 3714269896529642422852685702214695613036368,
        r,
        domain=sp.QQ,
    )
    check(kpoly == expected_k, "row-eleven ratio septimic")
    check(kpoly.degree() == 7 and kpoly.TC() != 0, "septimic degree and nonzero constant")
    check(kpoly.LC() != 0, "septimic has no projective root at infinity")
    check(sp.gcd(kpoly, kpoly.diff()) == sp.Poly(1, r, domain=sp.QQ), "septimic squarefree")
    check(sp.gcd(kpoly, apoly) == sp.Poly(1, r, domain=sp.QQ), "septimic avoids A zero")
    check(sp.gcd(kpoly, bpoly) == sp.Poly(1, r, domain=sp.QQ), "septimic avoids B zero")
    check(sp.gcd(apoly, bpoly) == sp.Poly(1, r, domain=sp.QQ), "ratio graph has no exceptional fibre")
    check(2 * kpoly.degree() == 14, "seven ratios times two nonzero sign sheets")

    # The independent scratch audit at 84590342a2 worked in the stricter
    # Phi*U!=0 chart.  Recover U from the inherited response (rather than
    # assuming that extra gate), pull its numerator to the ratio graph, and
    # prove that every one of our fourteen points actually lies there.
    u_on_graph = exact(
        bracket_subs[R8.U]
        .subs(terminal8)
        .subs(source9_subs)
        .subs(source10_subs)
    )
    expected_u = (
        32805000 * P**2 + 36905625 * P * eta - 1752089427968
    ) / 259524000
    check_zero(u_on_graph - expected_u, "inherited U response on row-ten graph")
    ubar = invariant_even(primitive_numerator(u_on_graph))
    uraw = sp.Poly(
        sp.together(ubar.as_expr().subs(Y, bpoly.as_expr() / apoly.as_expr()) * apoly.as_expr()),
        r,
        domain=sp.QQ,
    )
    ucontent, upoly = uraw.primitive()
    expected_upoly = sp.Poly(
        -1165769340615 * r**3
        - 16036650988830 * r**2
        - 117079277403336 * r
        + 15165313804484,
        r,
        domain=sp.QQ,
    )
    check(ucontent == 531505152000 and upoly == expected_upoly, "ratio U carrier")
    check(sp.gcd(kpoly, upoly) == sp.Poly(1, r, domain=sp.QQ), "septimic avoids U zero")

    a10 = [exact(row.subs(select11)) for row in a10]
    c10 = [exact(row.subs(select11)) for row in c10]
    a11, c11, theta11 = append_tangent(a10, c10, 11, "nu")
    ares11, cres11, amat11, arhs11, cmat11, crhs11, adata11, cdata11 = depth_data(
        a11,
        c11,
        11,
        theta11,
    )
    check((len(adata11[0]), adata11[1].cols, adata11[1].rank(), len(adata11[1].T.nullspace())) == (102, 228, 77, 25), "Phi-nonzero row-eleven P2 universe")
    check((len(cdata11[0]), cdata11[1].cols, cdata11[1].rank(), len(cdata11[1].T.nullspace())) == (114, 361, 94, 20), "Phi-nonzero row-eleven P3 universe")
    check(amat11.rank() == 4 and cmat11.rank() == 3 and amat11.col_join(cmat11).rank() == 4, "Phi-nonzero row-eleven depth ranks")
    check(amat11.rref()[1] == (8, 9, 10, 11) and cmat11.rref()[1] == (9, 10, 11), "Phi-nonzero row-eleven depth pivots")
    asolution11 = selected_solution(amat11, arhs11, (8, 9, 10, 11))
    csolution11 = selected_solution(cmat11, crhs11, (9, 10, 11))
    basis = sp.groebner([dpoly.as_expr(), epoly.as_expr()], eta, P, order="grlex", domain=sp.QQ)
    check(basis.is_zero_dimensional, "row-eleven source carrier zero dimensional")
    for j in range(3):
        check(reduce_ideal_numerator(asolution11[j + 1] - csolution11[j], basis) == 0, f"row-eleven shared depth coordinate {9 + j}")
    asubs11 = {theta11[8 + j]: asolution11[j] for j in range(4)}
    csubs11 = {theta11[9 + j]: csolution11[j] for j in range(3)}
    check(all(reduce_ideal_numerator(value.subs(asubs11), basis) == 0 for value in ares11), "row-eleven P2 depth on carrier")
    check(all(reduce_ideal_numerator(value.subs(csubs11), basis) == 0 for value in cres11), "row-eleven P3 depth on carrier")
    a11[-1] = exact(a11[-1].subs(asubs11))
    c11[-1] = exact(c11[-1].subs(asubs11))

    g12 = exact(
        R8.tcoeff(R8.H, 12)
        .subs(bracket_subs)
        .subs(terminal8)
        .subs(source9_subs)
        .subs(source10_subs)
    )
    expected_g12 = -13 * x**12 * (3918375 * P**2 + 1366875 * P * eta - 3318235136) / 32440500
    check_zero(g12 - expected_g12, "literal specialized G12")
    residuals12, matrix12, _, select12 = bracket_select(a11, c11, 12, g12, theta11[:8])
    check(matrix12.rank() == 8 and matrix12.T.rref()[1] == tuple(range(8)), "Phi-nonzero row-twelve selector")
    reduced12 = [reduce_ideal_numerator(value.subs(select12), basis) for value in residuals12]
    check(all(value == 0 or index in (8, 9) for index, value in enumerate(reduced12)), "row-twelve reduced residual support")
    ccore = sp.Poly(
        1245684290400 * P**2
        - 592839494100 * P * eta
        - 18612736875 * eta**2
        + 4820239249178624,
        P,
        eta,
        domain=sp.QQ,
    )
    check_zero(reduced12[9] + 1650 * P * ccore.as_expr(), "row-twelve decisive residual remainder")

    cratio = sp.Poly(1245684290400 - 592839494100 * r - 18612736875 * r**2, r, domain=sp.QQ)
    jraw = sp.Poly(
        bpoly.as_expr() * cratio.as_expr() + 4820239249178624 * apoly.as_expr(),
        r,
        domain=sp.QQ,
    )
    jcontent, jpoly = jraw.primitive()
    expected_j = sp.Poly(
        -4607940726340893 * r**2
        - 2340230825466590 * r
        + 113720641620096958,
        r,
        domain=sp.QQ,
    )
    check(jcontent == 7228470067200 and jpoly == expected_j, "row-twelve ratio quadratic")

    prime = 29
    k29 = sp.Poly(kpoly.as_expr(), r, modulus=prime)
    j29 = sp.Poly(jpoly.as_expr(), r, modulus=prime)
    expected_k29 = sp.Poly(8 * r**7 + 9 * r**6 - 10 * r**5 + 12 * r**4 + 2 * r**3 + 8 * r**2 - 12 * r + 3, r, modulus=prime)
    expected_j29 = sp.Poly(-11 * r**2 + 8 * r + 5, r, modulus=prime)
    check(k29 == expected_k29 and j29 == expected_j29, "mod29 retained degrees")
    check(k29.degree() == kpoly.degree() and j29.degree() == jpoly.degree(), "mod29 leading coefficients survive")
    bezout29 = sp.Poly(
        (11 * r + 11) * k29.as_expr()
        + (8 * r**6 + 7 * r**5 + 13 * r**4 + 12 * r**3 - 3 * r**2 + 8 * r + 11) * j29.as_expr()
        - 1,
        r,
        modulus=prime,
    )
    check(bezout29.is_zero, "mod29 K/J Bezout")
    check(sp.gcd(kpoly, jpoly) == sp.Poly(1, r, domain=sp.QQ), "characteristic-zero K/J coprime")

    # A finite-field carrier point is a positive control for every staged
    # equation through row eleven and a hostile for premature extinction.
    finite = {P: 15, eta: 13}
    check(modular_value(kpoly.as_expr(), {r: 5}, 31) == 0, "F31 K root")
    check(modular_value(apoly.as_expr(), {r: 5}, 31) != 0, "F31 A nonzero")
    check(modular_value(bpoly.as_expr(), {r: 5}, 31) == 8 * modular_value(apoly.as_expr(), {r: 5}, 31) % 31, "F31 Y=B/A")
    check(15 * 15 % 31 == 8 and 13 * pow(15, -1, 31) % 31 == 5, "F31 sign-sheet address")
    check(modular_value(xi10, finite, 31) == 20, "F31 xi")
    check(modular_value(alpha10, finite, 31) == 16, "F31 alpha")
    check(modular_value(beta10, finite, 31) == 1, "F31 beta")
    response31 = {
        symbol: modular_value(value.subs(terminal8).subs(source9_subs).subs(source10_subs), finite, 31)
        for symbol, value in bracket_subs.items()
    }
    check((response31[R8.U], response31[R8.W], response31[R8.Z]) == (26, 22, 5), "F31 strict response")
    check(modular_value(dpoly.as_expr(), finite, 31) == 0, "F31 D carrier")
    check(modular_value(epoly.as_expr(), finite, 31) == 0, "F31 E carrier")
    check(modular_value(ccore.as_expr(), finite, 31) == 14, "F31 row-twelve obstruction")
    check(modular_value(jpoly.as_expr(), {r: 5}, 31) == 8, "F31 J hostile")

    print(f"PHI_NONZERO_E9_ALPHA={sp.sstr(alpha9)}")
    print(f"PHI_NONZERO_ROW10_XI={sp.sstr(xi10)}")
    print(f"PHI_NONZERO_ROW10_ALPHA={sp.sstr(alpha10)}")
    print(f"PHI_NONZERO_ROW10_BETA={sp.sstr(beta10)}")
    print(f"PHI_NONZERO_ROW10_D={sp.sstr(dpoly.as_expr())}")
    print(f"RATIO_A={sp.sstr(apoly.as_expr())}")
    print(f"RATIO_B={sp.sstr(bpoly.as_expr())}")
    print(f"ROW11_K={sp.sstr(kpoly.as_expr())}")
    print("ROW11_CARRIER K_degree=7 squarefree=yes gcd(K,A*B)=1 source_points=14 terminal_fibre=A8")
    print(f"ROW11_U_CARRIER={sp.sstr(upoly.as_expr())} gcd(K,Ucarrier)=1")
    print(f"ROW12_J={sp.sstr(jpoly.as_expr())}")
    print("MOD29_K=8*r^7+9*r^6-10*r^5+12*r^4+2*r^3+8*r^2-12*r+3")
    print("MOD29_J=-11*r^2+8*r+5")
    print("MOD29_BEZOUT=(11*r+11)*K+(8*r^6+7*r^5+13*r^4+12*r^3-3*r^2+8*r+11)*J=1")
    print("F31_CONTROL r=5 Y=8 Phi=15 eta=13 xi=20 alpha=16 beta=1 U=26 W=22 Z=5 D=E=0 J=8")


def main() -> None:
    _, _, bracket_subs, terminal8, _, fixed_a, fixed_c, g9 = build_row_eight()
    phi_zero_branch(fixed_a, fixed_c, bracket_subs, terminal8, g9)
    phi_nonzero_branch(fixed_a, fixed_c, bracket_subs, terminal8, g9)
    print("THM-4380 PRIMARY: SOURCE-NORMAL WEIGHT-TWELVE ROW-TWELVE EXTINCTION PASS")
    print(f"CHECKS={CHECKS}")
    print("ROW9 full E9=0 projected-depth automatic; terminal ranks P2/P3/joint=3/2/3; fibre A7")
    print("PHI_ZERO dies at row11 after two genuine row10 joint-depth source points")
    print("PHI_NONZERO row10 joint depth is ratio-presented affine curve; row11 has 14 genuine points; row11 depth automatic; row12 dies by gcd(K,J)=1")
    print("SCOPE exact finite source-normal residual_weight<=12 projected-depth extinction only; no all-row lift, chart/seam entry, Keller pair, JC2, or DC2")


if __name__ == "__main__":
    main()
