#!/usr/bin/env python3
"""Exact verifier for the first weight-13 source-normal extension (THM-4389).

Rebuilds the complete bracket/projected-depth pipeline through row ten after
adjoining p^5 y and p^2 y^3, then audits the resulting elliptic pencil and a
rational hostile outside the THM-4385 central fibre.
"""
from __future__ import annotations

from pathlib import Path
import sys
import sympy as sp

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))
import jc2_source_normal_bracket_hasse_rows8_thm4308 as R8
import jc2_source_normal_student_stein_row9_thm4315 as R9

x, t = R8.x, R8.t
rho13, sigma13 = sp.symbols("rho13 sigma13")
H13 = sp.expand(R8.H + rho13 * R8.p**5 * R8.y + sigma13 * R8.p**2 * R8.y**3)
CHECKS = 0


def check(condition, label):
    global CHECKS
    if not bool(condition):
        raise AssertionError(label)
    CHECKS += 1


def zero(expr):
    return sp.cancel(sp.together(sp.expand(expr))) == 0


def build():
    grows = {n: sp.expand(R8.tcoeff(H13, n)) for n in range(4, 14)}
    arows = [R8.A0, R8.A1, R8.A2, R8.A3]
    crows = [R8.C0, R8.C1, R8.C2, R8.C3]
    check(zero(R8.predicted_G(4, arows, crows) - grows[4]), "row-four inheritance")
    subs = {}
    solve_symbols = {5: R8.upsilon5, 6: R8.U, 7: R8.W, 8: R8.Z}
    obstructions = {}
    for n in range(4, 8):
        b = R8.B_row(n, arows, crows)
        abase, cbase = R8.particular_row(n, b)
        m = n + 1
        diff = sp.expand(grows[m].subs(subs) - R8.predicted_G(m, arows + [abase], crows + [cbase]))
        obs = sp.factor(R8.compatibility_value(m, diff))
        obstructions[m] = obs
        sym = solve_symbols[m]
        sol = sp.solve(obs, sym, dict=True)[0][sym]
        subs[sym] = sp.cancel(sol)
        diff = sp.expand(diff.subs(sym, subs[sym]))
        theta = R8.tangent_solve(m, diff)
        arows.append(sp.expand(abase + theta * sp.diff(R8.A0, x)))
        crows.append(sp.expand(cbase + theta * sp.diff(R8.C0, x)))
        check(zero(R8.predicted_G(m, arows, crows) - grows[m].subs(subs)), f"row-{m} bracket reconstruction")

    b8 = R8.B_row(8, arows, crows)
    abase, cbase = R8.particular_row(8, b8)
    th = list(sp.symbols("w13theta8_0:9"))
    theta = sum(th[j] * x**j for j in range(9))
    arows.append(sp.expand(abase + theta * sp.diff(R8.A0, x)))
    crows.append(sp.expand(cbase + theta * sp.diff(R8.C0, x)))
    check(zero(R8.predicted_G(8, arows[:8], crows[:8]) - grows[8].subs(subs)), "row-eight compatibility")
    return grows, arows, crows, subs, th, obstructions


def row8_gate(arows, crows, subs, th):
    acoords, _, am = R8.depth_matrix(2)
    ccoords, _, cm = R8.depth_matrix(3)
    avec = sp.Matrix([R8.xcoeff(arows[n], d) for n, d in acoords])
    cvec = sp.Matrix([R8.xcoeff(crows[n], d) for n, d in ccoords])
    residuals = [sp.expand((l.T * avec)[0].subs(subs)) for l in am.T.nullspace()]
    residuals += [sp.expand((l.T * cvec)[0].subs(subs)) for l in cm.T.nullspace()]
    mat, rhs = sp.linear_eq_to_matrix(residuals, th)
    check(mat.rank() == 2, "row-eight terminal depth rank")
    compatibility = [sp.factor((l.T * rhs)[0]) for l in mat.T.nullspace()]
    compatibility = [q for q in compatibility if q != 0]
    # Primitive numerator / monic normalization up to rational units.
    uniq = []
    for q in compatibility:
        num = sp.Poly(sp.fraction(sp.cancel(q))[0], R8.Delta, R8.Theta, R8.zeta3, rho13, sigma13,
                      R8.Phi, R8.eta, R8.xi10, domain=sp.QQ)
        prim = num.primitive()[1]
        if prim.LC() < 0:
            prim = -prim
        expr = sp.factor(prim.as_expr())
        if expr not in uniq:
            uniq.append(expr)
    return residuals, mat, rhs, compatibility, uniq


def selected_solution(matrix, rhs, columns):
    row_indices = tuple(matrix.T.rref()[1])
    square = matrix.extract(row_indices, columns)
    check(square.rows == square.cols and square.det() != 0, "selected subsystem invertible")
    return tuple(sp.cancel(v) for v in square.inv() * rhs.extract(row_indices, (0,)))


def bracket_select(arows, crows, row, source_row, free_symbols):
    difference = sp.cancel(sp.together(sp.expand(source_row - R8.predicted_G(row, arows, crows))))
    check(difference == 0 or sp.Poly(difference, x, domain=sp.EX).degree() <= row,
          f"row-{row} full degree cap")
    equations = [sp.cancel(sp.together(R8.xcoeff(difference, d))) for d in range(row + 1)]
    matrix, rhs = sp.linear_eq_to_matrix(equations, free_symbols)
    solution = selected_solution(matrix, rhs, tuple(range(len(free_symbols))))
    substitutions = {symbol: solution[j] for j, symbol in enumerate(free_symbols)}
    residuals = [sp.factor(sp.cancel(sp.together(value.subs(substitutions)))) for value in equations]
    return residuals, matrix, substitutions


def row8_terminal(arows, crows, subs, th):
    residuals, mat, rhs, _, _ = row8_gate(arows, crows, subs, th)
    gates = {R8.Delta: sp.Rational(896, 15), R8.Theta: sp.Rational(512, 75), R8.zeta3: -3 * R8.Phi / 2}
    gated = [sp.expand(q.subs(gates)) for q in residuals]
    sol = sp.solve(gated, th[-2:], dict=True)
    check(len(sol) == 1, "row-eight terminal solution unique in pivot coordinates")
    terminal = {**gates, **sol[0]}
    check(all(zero(q.subs(terminal)) for q in residuals), "all row-eight depth residuals")
    return terminal


def row9_bracket(grows, arows, crows, subs, th, terminal):
    fixed_a = [sp.expand(row.subs(subs).subs(terminal)) for row in arows]
    fixed_c = [sp.expand(row.subs(subs).subs(terminal)) for row in crows]
    residuals, matrix, tangent = bracket_select(
        fixed_a, fixed_c, 9, sp.expand(grows[9].subs(subs).subs(terminal)), tuple(th[:7])
    )
    nonzero = [(i, q) for i, q in enumerate(residuals) if q != 0]
    check(matrix.rank() == 7 and len(nonzero) == 1, "row-nine selector and single gate")
    num = sp.Poly(sp.fraction(sp.cancel(nonzero[0][1]))[0],
                  R8.Phi, R8.eta, R8.xi10, R8.alpha11, R8.beta11, rho13, sigma13,
                  domain=sp.QQ).primitive()[1]
    return fixed_a, fixed_c, tangent, nonzero[0][0], sp.factor(num.as_expr())


def append_tangent(arows, crows, row, prefix):
    b = R8.B_row(row, arows, crows)
    abase, cbase = R8.particular_row(row, b)
    symbols = tuple(sp.symbols(f"{prefix}0:{row+1}"))
    theta = sum(symbols[j] * x**j for j in range(row + 1))
    return (
        arows + [sp.expand(abase + theta * sp.diff(R8.A0, x))],
        crows + [sp.expand(cbase + theta * sp.diff(R8.C0, x))],
        symbols,
    )


def depth_data(arows, crows, row, symbols):
    acoords, am = R9.depth_matrix(2, row)
    ccoords, cm = R9.depth_matrix(3, row)
    avec = sp.Matrix([R8.xcoeff(arows[n], d) for n, d in acoords])
    cvec = sp.Matrix([R8.xcoeff(crows[n], d) for n, d in ccoords])
    ares = [sp.expand((l.T * avec)[0]) for l in am.T.nullspace()]
    cres = [sp.expand((l.T * cvec)[0]) for l in cm.T.nullspace()]
    amat, arhs = sp.linear_eq_to_matrix(ares, symbols)
    cmat, crhs = sp.linear_eq_to_matrix(cres, symbols)
    return ares, cres, amat, arhs, cmat, crhs


def row9_depth_and_row10(grows, fixed_a, fixed_c, tangent9, e9, bracket_subs, terminal8):
    # Use the globally valid xi graph for the unchanged row-nine gate.
    xi_graph = sp.solve(e9, R8.xi10, dict=True)[0][R8.xi10]
    source9 = {R8.xi10: sp.cancel(xi_graph)}
    a9base = [sp.expand(row.subs(tangent9).subs(source9)) for row in fixed_a]
    c9base = [sp.expand(row.subs(tangent9).subs(source9)) for row in fixed_c]
    a9, c9, q9 = append_tangent(a9base, c9base, 9, "w13q9_")
    ares, cres, amat, arhs, cmat, crhs = depth_data(a9, c9, 9, q9)
    joint = amat.col_join(cmat)
    rhs = arhs.col_join(crhs)
    check((amat.rank(), cmat.rank(), joint.rank()) == (3, 2, 3), "row-nine depth ranks")
    pivots = joint.rref()[1]
    check(pivots == (7, 8, 9), "row-nine joint pivots")
    sol = selected_solution(joint, rhs, pivots)
    terminal9 = {q9[pivots[j]]: sol[j] for j in range(3)}
    compat = [sp.factor(sp.cancel(sp.together(q.subs(terminal9)))) for q in ares + cres]
    compat_nonzero = [q for q in compat if q != 0]
    a9fixed = [sp.expand(row.subs(terminal9)) for row in a9]
    c9fixed = [sp.expand(row.subs(terminal9)) for row in c9]
    g10 = sp.expand(grows[10].subs(bracket_subs).subs(terminal8).subs(source9))
    res10, mat10, tangent10 = bracket_select(a9fixed, c9fixed, 10, g10, tuple(q9[:7]))
    nonzero10 = [(i, q) for i, q in enumerate(res10) if q != 0]
    prim10 = []
    vars10 = (R8.Phi, R8.eta, R8.alpha11, R8.beta11, rho13, sigma13)
    for i, q in nonzero10:
        poly = sp.Poly(sp.fraction(sp.cancel(q))[0], *vars10, domain=sp.QQ).primitive()[1]
        if poly.LC() < 0:
            poly = -poly
        prim10.append((i, sp.factor(poly.as_expr())))
    return source9, (amat.rank(), cmat.rank(), joint.rank()), compat_nonzero, terminal9, mat10.rank(), prim10, a9fixed, c9fixed, q9, tangent10


def row10_depth_p_nonzero(a9fixed, c9fixed, q9, tangent10, gates10):
    equations = {i: q for i, q in gates10}
    # The position-eight gate solves alpha on P!=0; position six then solves beta.
    alpha_graph = sp.factor(sp.solve(equations[8], R8.alpha11, dict=True)[0][R8.alpha11])
    beta_graph = sp.factor(
        sp.solve(equations[6].subs(R8.alpha11, alpha_graph), R8.beta11, dict=True)[0][R8.beta11]
    )
    source10 = {R8.alpha11: alpha_graph, R8.beta11: beta_graph}
    check(zero(equations[8].subs(source10)) and zero(equations[6].subs(source10)),
          "row-ten bracket source graphs")
    abase = [sp.expand(row.subs(tangent10).subs(source10)) for row in a9fixed]
    cbase = [sp.expand(row.subs(tangent10).subs(source10)) for row in c9fixed]
    a10, c10, q10 = append_tangent(abase, cbase, 10, "w13q10_")
    ares, cres, amat, arhs, cmat, crhs = depth_data(a10, c10, 10, q10)
    joint = amat.col_join(cmat)
    rhs = arhs.col_join(crhs)
    check((amat.rank(), cmat.rank(), joint.rank()) == (3, 3, 3), "row-ten nonzero-P depth ranks")
    compat = [sp.factor(sp.cancel(sp.together((l.T * rhs)[0]))) for l in joint.T.nullspace()]
    compat = [q for q in compat if q != 0]
    nums = []
    vars4 = (R8.Phi, R8.eta, rho13, sigma13)
    for q in compat:
        num = sp.Poly(sp.fraction(q)[0], *vars4, domain=sp.QQ).primitive()[1]
        if num.LC() < 0:
            num = -num
        nums.append(num)
    common = nums[0]
    for poly in nums[1:]:
        common = sp.gcd(common, poly)
    return (
        source10,
        (amat.rank(), cmat.rank(), joint.rank()),
        len(compat),
        sp.factor(common.as_expr()),
        nums,
        joint,
        rhs,
        a10,
        c10,
        q10,
        ares,
        cres,
    )


def reduce_mod_univariate(expr, var, modulus):
    val = sp.cancel(sp.together(expr))
    num, den = sp.fraction(val)
    mod = sp.Poly(modulus, var, domain=sp.EX)
    np = sp.Poly(num, var, domain=sp.EX)
    dp = sp.Poly(den, var, domain=sp.EX)
    inv = sp.invert(dp, mod)
    return sp.factor((np * inv).rem(mod).as_expr())


def row10_depth_p_zero(a9fixed, c9fixed, q9, tangent10, gates10):
    equations = {i: q for i, q in gates10}
    f0 = sp.Poly(equations[8].subs(R8.Phi, 0), R8.eta, domain=sp.QQ).primitive()[1]
    alpha0 = sp.factor(
        sp.solve(equations[6].subs(R8.Phi, 0), R8.alpha11, dict=True)[0][R8.alpha11]
    )
    source0 = {R8.Phi: 0, R8.alpha11: alpha0}
    abase = [sp.expand(row.subs(tangent10).subs(source0)) for row in a9fixed]
    cbase = [sp.expand(row.subs(tangent10).subs(source0)) for row in c9fixed]
    a10, c10, q10 = append_tangent(abase, cbase, 10, "w13zq10_")
    ares, cres, amat, arhs, cmat, crhs = depth_data(a10, c10, 10, q10)
    joint = amat.col_join(cmat)
    check((amat.rank(), cmat.rank(), joint.rank()) == (3, 3, 3), "row-ten P-zero depth ranks")
    check(amat.rref()[1] == (8, 9, 10) and cmat.rref()[1] == (8, 9, 10),
          "row-ten P-zero pivots")
    asol = selected_solution(amat, arhs, (8, 9, 10))
    csol = selected_solution(cmat, crhs, (8, 9, 10))
    differences = [reduce_mod_univariate(asol[i] - csol[i], R8.eta, f0.as_expr()) for i in range(3)]
    return f0, alpha0, differences


def binary_quartic_invariants(poly, variable):
    coeffs = sp.Poly(poly, variable, domain=sp.EX).all_coeffs()
    check(len(coeffs) == 5, "quartic coefficient count")
    a4, a3, a2, a1, a0 = coeffs
    inv_i = sp.expand(12 * a4 * a0 - 3 * a3 * a1 + a2**2)
    inv_j = sp.expand(
        72 * a4 * a2 * a0
        + 9 * a3 * a2 * a1
        - 27 * a4 * a1**2
        - 27 * a3**2 * a0
        - 2 * a2**3
    )
    delta_flat = sp.expand(4 * inv_i**3 - inv_j**2)
    return inv_i, inv_j, delta_flat


def elliptic_pencil_audit(d13):
    P, e = R8.Phi, R8.eta
    rr, aa = sp.symbols("ratio a13")
    normal = sp.Integer(707284001250)
    d0 = (
        7231154026500 * P**3
        + 50541940696500 * P**2 * e
        + 6793915500000 * P * e**2
        - 631918028977864704 * P
        + 353642000625 * e**3
        - 91584545734393856 * e
    )
    expected = sp.expand(d0 - normal * P**2 * rho13)
    check(zero(d13 - expected), "weight-thirteen normal equation")
    check(sp.diff(expected, rho13) == -normal * P**2, "p5y transverse derivative")
    check(sp.diff(expected, sigma13) == 0, "p2y3 tangent derivative")

    ar = (
        353642000625 * rr**3
        + 6793915500000 * rr**2
        + 50541940696500 * rr
        + 7231154026500
    )
    br = 91584545734393856 * rr + 631918028977864704
    ratio_equation = sp.expand(ar * P**2 - normal * aa * P - br)
    check(zero(expected.subs({e: rr * P, rho13: aa}) - P * ratio_equation),
          "ratio quadratic pencil")
    sheet = 2 * ar * P - normal * aa
    quartic = sp.expand(4 * ar * br + normal**2 * aa**2)
    check(zero(sheet**2 - quartic - 4 * ar * ratio_equation), "completed-square double cover")

    inv_i, inv_j, delta = binary_quartic_invariants(quartic, rr)
    q0_mod7 = sp.Poly(quartic.subs(aa, 0), rr, modulus=7)
    q1_mod7 = sp.Poly(quartic.subs(aa, 1), rr, modulus=7)
    check(tuple(int(v) % 7 for v in q0_mod7.all_coeffs()) == (1, 1, 3, 2, 2),
          "central quartic mod seven")
    check(tuple(int(v) % 7 for v in q1_mod7.all_coeffs()) == (1, 1, 3, 2, 3),
          "unit transverse quartic mod seven")
    inv_data = []
    for aval in (0, 1):
        iv = int(inv_i.subs(aa, aval)) % 7
        jv = int(inv_j.subs(aa, aval)) % 7
        dv = int(delta.subs(aa, aval)) % 7
        jflat = iv**3 * pow(dv, -1, 7) % 7
        inv_data.append((iv, jv, dv, jflat))
    check(inv_data == [(6, 4, 1, 6), (4, 4, 2, 4)], "distinct smooth j-invariant fibres mod seven")

    delta_mod7 = sp.Poly(delta, aa, modulus=7)
    expected_delta_mod7 = sp.Poly(3 * aa**6 + aa**4 - 3 * aa**2 + 1, aa, modulus=7)
    check(delta_mod7 == expected_delta_mod7, "pencil discriminant mod seven")
    check(all(int(delta_mod7.eval(v)) % 7 != 0 for v in range(7)),
          "no rational residue singular fibre mod seven")
    return normal, d0, ar, br, quartic, inv_data, delta_mod7


def hostile_audit(
    source9,
    source10,
    d13,
    gates10,
    joint10,
    rhs10,
    bracket_subs,
    a10,
    c10,
    q10,
    ares10,
    cres10,
):
    P, e = R8.Phi, R8.eta
    hostile = {
        P: sp.Integer(1),
        e: sp.Integer(0),
        rho13: -sp.Rational(35106155434657678, 39293555625),
        sigma13: sp.Integer(1),
    }
    alpha_value = sp.factor(source10[R8.alpha11].subs(hostile))
    beta_value = sp.factor(source10[R8.beta11].subs(hostile))
    hostile[R8.alpha11] = alpha_value
    hostile[R8.beta11] = beta_value
    xi_value = sp.factor(source9[R8.xi10].subs(hostile))
    hostile[R8.xi10] = xi_value
    check(alpha_value == sp.Rational(2410085120016937, 18612736875), "hostile alpha")
    check(beta_value == sp.Rational(91, 15), "hostile beta")
    check(xi_value == sp.Rational(120573932063, 7209000), "hostile xi")
    check(zero(d13.subs(hostile)), "hostile lies on weight-thirteen carrier")
    check(zero(R9.E9_GATE.subs(hostile)), "hostile row-nine source gate")

    d0 = sp.expand(d13.subs(rho13, 0))
    check(d0.subs(hostile) == -631910797823838204, "hostile misses central elliptic fibre")
    check(all(zero(q.subs(hostile)) for _, q in gates10), "hostile row-ten bracket gates")
    check(joint10.row_join(rhs10.subs(hostile)).rank() == joint10.rank() == 3,
          "hostile row-ten joint depth solvable")

    hmatrix = joint10.subs(hostile)
    hrhs = rhs10.subs(hostile)
    terminal_values = selected_solution(hmatrix, hrhs, (8, 9, 10))
    terminal = {q10[j]: 0 for j in range(8)}
    terminal.update({q10[8 + j]: terminal_values[j] for j in range(3)})
    full_point = {**hostile, **terminal}
    check(all(zero(q.subs(full_point)) for q in ares10 + cres10),
          "hostile all row-ten depth residuals")

    direct_source = {
        R8.Delta: sp.Rational(896, 15),
        R8.Theta: sp.Rational(512, 75),
        R8.zeta3: -sp.Rational(3, 2) * hostile[P],
        **hostile,
    }
    for symbol, value in bracket_subs.items():
        direct_source[symbol] = sp.cancel(value.subs(direct_source))
    full_point.update(direct_source)
    aa = sum(a10[n].subs(full_point) * t**n for n in range(11))
    cc = sum(c10[n].subs(full_point) * t**n for n in range(11))
    source_g = sp.expand(-R8.u / 2 + H13).subs(full_point)
    pdefect = sp.expand(R8.P(aa, cc) - source_g)
    jdefect = sp.expand(R8.jacobian(aa, cc))
    pvals = [sp.factor(R8.tcoeff(pdefect, n)) for n in range(11)]
    jvals = [sp.factor(R8.tcoeff(jdefect, n) - (1 if n == 0 else 0)) for n in range(10)]
    if any(not zero(v) for v in pvals):
        print("DEBUG_PDEFECT", [(i, sp.sstr(v)) for i, v in enumerate(pvals) if not zero(v)])
    if any(not zero(v) for v in jvals):
        print("DEBUG_JDEFECT", [(i, sp.sstr(v)) for i, v in enumerate(jvals) if not zero(v)])
    check(all(zero(v) for v in pvals),
          "hostile direct P/G rows zero through ten")
    check(zero(jvals[0]), "hostile direct Jacobian row zero")
    check(all(zero(jvals[n]) for n in range(1, 10)),
          "hostile direct Jacobian rows one through nine")
    check(
        all(sp.degree(a10[n].subs(full_point), x) <= n + 1 for n in range(4, 11))
        and all(sp.degree(c10[n].subs(full_point), x) <= n + 2 for n in range(4, 11)),
        "hostile row degree caps",
    )

    response_subs = {
        R8.Delta: sp.Rational(896, 15),
        R8.Theta: sp.Rational(512, 75),
        R8.zeta3: -3 * P / 2,
        **hostile,
    }
    uval = sp.factor(bracket_subs[R8.U].subs(response_subs))
    wval = sp.factor(bracket_subs[R8.W].subs(response_subs))
    zval = sp.factor(bracket_subs[R8.Z].subs(response_subs))
    check(uval == -sp.Rational(219007077871, 32440500), "hostile U")
    check(wval == sp.Rational(25779284581, 720900), "hostile W")
    check(zval == -sp.Rational(1231789589, 44500), "hostile Z")
    lam = sp.factor(uval + wval + zval)
    discr = sp.factor(wval**2 - 4 * uval * zval)
    check(lam == sp.Rational(43086117893, 32440500), "hostile Lambda")
    check(discr == sp.Rational(6902544853907445578341, 12992420250000), "hostile response discriminant")
    check(hostile[rho13] * hostile[sigma13] * uval * zval * (hostile[rho13] + hostile[sigma13]) != 0,
          "hostile dense weight-thirteen face units")
    return hostile, uval, wval, zval, lam, discr


def main():
    grows, arows, crows, subs, th, obs = build()
    print("WEIGHT13_SOURCE_ROWS")
    expected_dg = {
        4: 0,
        5: 0,
        6: 0,
        7: rho13 * x,
        8: (6 * rho13 + sigma13) * x**3,
        9: 5 * (3 * rho13 + sigma13) * x**5,
        10: 10 * (2 * rho13 + sigma13) * x**7,
        11: 5 * (3 * rho13 + 2 * sigma13) * x**9,
        12: (6 * rho13 + 5 * sigma13) * x**11,
        13: (rho13 + sigma13) * x**13,
    }
    for n in range(4, 14):
        delta = sp.factor(grows[n] - R8.tcoeff(R8.H, n))
        check(zero(delta - expected_dg[n]), f"weight-thirteen source row {n}")
        print(f"dG{n}={sp.sstr(delta)}")
    expected_solutions = {
        R8.upsilon5: R8.upsilon_solution,
        R8.U: R8.U_solution,
        R8.W: R8.W_solution,
        R8.Z: R8.Z_solution,
    }
    for symbol, expected in expected_solutions.items():
        check(zero(subs[symbol] - expected), f"unchanged bracket solution {symbol}")
    print("BRACKET_OBSTRUCTIONS")
    for m in sorted(obs):
        print(f"obs{m}={sp.sstr(obs[m])}")
        print(f"solution_{m}={sp.sstr(subs[{5:R8.upsilon5,6:R8.U,7:R8.W,8:R8.Z}[m]])}")
    residuals, mat, rhs, compat, uniq = row8_gate(arows, crows, subs, th)
    print(f"ROW8_DEPTH tangent_shape={mat.shape} rank={mat.rank()} left_nullity={len(mat.T.nullspace())}")
    print(f"ROW8_COMPAT_NONZERO={len(compat)} UNIQUE={len(uniq)}")
    check(
        {sp.sstr(q) for q in uniq}
        == {
            sp.sstr(15 * R8.Delta - 896),
            sp.sstr(3 * R8.Phi + 2 * R8.zeta3),
            sp.sstr(3030 * R8.Delta + 225 * R8.Theta - 182528),
        },
        "unchanged row-eight source gate ideal",
    )
    for i, q in enumerate(uniq):
        print(f"gate{i}={sp.sstr(q)}")
    terminal = row8_terminal(arows, crows, subs, th)
    print("ROW8_TERMINAL_LAST_TWO")
    for q in th[-2:]:
        print(f"{q}={sp.sstr(sp.factor(terminal[q]))}")
    fixed_a, fixed_c, tangent9, pos9, e9 = row9_bracket(grows, arows, crows, subs, th, terminal)
    check(zero(e9 + R9.E9_GATE), "unchanged row-nine source gate")
    print(f"ROW9 selector_rank=7 residual_position={pos9}")
    print(f"ROW9_GATE={sp.sstr(e9)}")
    print("ROW9_TANGENT")
    for q in th[:7]:
        print(f"{q}={sp.sstr(sp.factor(tangent9[q]))}")
    source9, ranks9, compat9, terminal9, rank10, gates10, a9fixed, c9fixed, q9, tangent10 = row9_depth_and_row10(
        grows, fixed_a, fixed_c, tangent9, e9, subs, terminal
    )
    print(f"ROW9_XI_GRAPH={sp.sstr(source9[R8.xi10])}")
    print(f"ROW9_DEPTH_RANKS={ranks9} residuals_after_solution={len(compat9)}")
    print(f"ROW10_SELECTOR_RANK={rank10}")
    expected_gate6 = (
        1283973857930250 * R8.Phi**2
        - 1709723764698750 * R8.Phi * R8.alpha11
        - 264168692280000 * R8.Phi * R8.beta11
        + 2669098150944000 * R8.Phi * R8.eta
        - 278844730740000 * R8.Phi * rho13
        - 278844730740000 * R8.alpha11 * R8.eta
        - 854861882349375 * R8.eta**2
        - 27743789318253707264
    )
    expected_gate8 = (
        69009144750 * R8.Phi**2
        + 37225473750 * R8.Phi * R8.alpha11
        + 357574500000 * R8.Phi * R8.eta
        + 18612736875 * R8.eta**2
        - 4820239249178624
    )
    check(dict(gates10).keys() == {6, 8}, "row-ten residual support")
    check(zero(dict(gates10)[6] - expected_gate6), "row-ten gate six")
    check(zero(dict(gates10)[8] - expected_gate8), "row-ten gate eight")
    for i, q in gates10:
        print(f"ROW10_GATE_POS{i}={sp.sstr(q)}")
    (
        source10,
        ranks10,
        ncompat10,
        d13,
        nums10,
        joint10,
        rhs10,
        a10,
        c10,
        q10,
        ares10,
        cres10,
    ) = row10_depth_p_nonzero(
        a9fixed, c9fixed, q9, tangent10, gates10
    )
    print(f"ROW10_ALPHA_GRAPH={sp.sstr(source10[R8.alpha11])}")
    print(f"ROW10_BETA_GRAPH={sp.sstr(source10[R8.beta11])}")
    print(f"ROW10_DEPTH_RANKS={ranks10} nonzero_compat={ncompat10}")
    print(f"ROW10_DEPTH_GCD={sp.sstr(d13)}")
    quotients = []
    common = sp.Poly(d13, R8.Phi, R8.eta, rho13, sigma13, domain=sp.QQ)
    for poly in nums10:
        q, rem = sp.div(poly, common)
        check(rem == 0, "row-ten compatibility generated by D13")
        expr = sp.factor(q.as_expr())
        if expr not in quotients:
            quotients.append(expr)
    print(f"ROW10_DEPTH_QUOTIENTS={len(quotients)}")
    for i, q in enumerate(quotients):
        print(f"q{i}={sp.sstr(q)}")
    normal, d0, ar, br, quartic, inv_data, delta_mod7 = elliptic_pencil_audit(d13)
    print(f"NORMAL_COEFFICIENT={normal}")
    print(f"RATIO_PENCIL=A(r)P^2-{normal}*a13*P-B(r)")
    print("QUARTIC_PENCIL=y^2=4*A(r)*B(r)+NORMAL_COEFFICIENT^2*a13^2")
    print(f"MOD7_INVARIANTS_a0_a1={inv_data}")
    print(f"MOD7_DISCRIMINANT={sp.sstr(delta_mod7.as_expr())}")
    f0, alpha0, differences0 = row10_depth_p_zero(a9fixed, c9fixed, q9, tangent10, gates10)
    expected_f0 = 18612736875 * R8.eta**2 - 4820239249178624
    check(zero(f0.as_expr() - expected_f0), "unchanged P-zero boundary polynomial")
    expected_differences0 = [-sp.Rational(4, 15) * (5 * R8.beta11 + 6 * R8.eta), 0, 0]
    check(all(zero(actual - expected) for actual, expected in zip(differences0, expected_differences0)),
          "unchanged P-zero joint depth gate")
    check(zero(d13.subs(R8.Phi, 0) - 19 * R8.eta * expected_f0), "boundary gluing factorization")
    print(f"ROW10_P0_F={sp.sstr(f0.as_expr())}")
    print(f"ROW10_P0_ALPHA={sp.sstr(alpha0)}")
    for i, q in enumerate(differences0):
        print(f"ROW10_P0_DEPTH_DIFF{i}={sp.sstr(q)}")
    hostile, uval, wval, zval, lam, discr = hostile_audit(
        source9,
        source10,
        d13,
        gates10,
        joint10,
        rhs10,
        subs,
        a10,
        c10,
        q10,
        ares10,
        cres10,
    )
    print("RATIONAL_HOSTILE")
    print(f"Phi={hostile[R8.Phi]} eta={hostile[R8.eta]}")
    print(f"a13={hostile[rho13]} b13={hostile[sigma13]}")
    print(f"xi10={hostile[R8.xi10]} alpha11={hostile[R8.alpha11]} beta11={hostile[R8.beta11]}")
    print(f"U={uval} W={wval} Z={zval} Lambda={lam} response_D={discr}")
    print("HOSTILE survives row10, has dense weight13 units, and misses D0=0")
    print(f"CHECKS={CHECKS}")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
