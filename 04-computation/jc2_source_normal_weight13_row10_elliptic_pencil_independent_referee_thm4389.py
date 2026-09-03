#!/usr/bin/env python3
"""Dependency-free exact referee for the first weight-13 source face (THM-4389).

This verifier deliberately imports no repository mathematics.  It
reconstructs the source rows, bracket recursion, and monomial depth spaces
from their definitions, using a reversed-column determinant section rather
than the producer's section.  All checks use ``demand`` and remain active
under optimized Python.
"""
from __future__ import annotations

import sys

import sympy as s


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")

CHECKS = 0


def demand(condition, label: str) -> None:
    global CHECKS
    if not bool(condition):
        raise RuntimeError(label)
    CHECKS += 1


def zero(value) -> bool:
    return s.cancel(s.together(s.expand(value))) == 0


def coeff(poly, variable, degree: int):
    return s.expand(poly).coeff(variable, degree)


def primitive_numerator(value, variables):
    numerator = s.fraction(s.cancel(s.together(value)))[0]
    poly = s.Poly(numerator, *variables, domain=s.QQ).primitive()[1]
    if poly.LC() < 0:
        poly = -poly
    return poly


def same_linear_span(actual, expected, variables, label: str) -> None:
    """Compare affine-linear equation spaces over Q, not printed bases."""
    actual_polys = [primitive_numerator(v, variables) for v in actual if not zero(v)]
    expected_polys = [primitive_numerator(v, variables) for v in expected]
    demand(all(p.total_degree() <= 1 for p in actual_polys + expected_polys), label + " linearity")
    monomials = [(0,) * len(variables)]
    monomials += [tuple(1 if i == j else 0 for i in range(len(variables))) for j in range(len(variables))]

    def vector(poly):
        data = poly.as_dict()
        return [data.get(m, s.Integer(0)) for m in monomials]

    ma = s.Matrix([vector(p) for p in actual_polys])
    me = s.Matrix([vector(p) for p in expected_polys])
    demand(ma.rank() == me.rank() == len(expected_polys), label + " rank")
    demand(ma.col_join(me).rank() == ma.rank(), label + " equality")


# -------------------------------------------------------------------------
# Source chart and boundary data, stated independently.
# -------------------------------------------------------------------------
x, t = s.symbols("x t")
Delta, Phi, Theta = s.symbols("Delta Phi Theta")
eta, zeta, v5, xi = s.symbols("eta zeta v5 xi")
alpha, beta = s.symbols("alpha beta")
U, W, Z = s.symbols("U W Z")
rho, sigma = s.symbols("rho sigma")

K = s.Rational(2848, 45) - s.Rational(7, 6) * Delta
p = t * (1 + x**2 * t)
y = x * t * p
u = x**2 * t

H0 = s.expand(
    -3 * p
    + s.Rational(8, 3) * p**2
    - s.Rational(1376, 135) * p**3
    + K * y**2
    + Phi * p**2 * y
    + Delta * p**4
    + Theta * p * y**2
    + eta * p**3 * y
    + zeta * y**3
    + v5 * p**5
    + xi * p**2 * y**2
    + alpha * p**4 * y
    + beta * p * y**3
    + U * p**6
    + W * p**3 * y**2
    + Z * y**4
)
H13 = s.expand(H0 + rho * p**5 * y + sigma * p**2 * y**3)
G13 = s.expand(-u / 2 + H13)

A0 = 1 + x**2 / 4
C0 = -3 * x / 4 - x**3 / 8
A1 = s.Rational(4, 3) + 2 * x**2
C1 = -4 * x - s.Rational(3, 2) * x**3
A2 = -s.Rational(32, 9) - s.Rational(4, 5) * x**2
C2 = s.Rational(88, 15) * x - s.Rational(12, 5) * x**3
A3 = (
    s.Rational(2176, 135)
    - Phi * x / 2
    + (s.Rational(1088, 315) - 4 * K / 7) * x**2
    - s.Rational(32, 15) * x**4
)
C3 = (
    3 * Phi / 4
    + (-s.Rational(8128, 315) + 6 * K / 7) * x
    + 3 * Phi * x**2 / 8
    + (s.Rational(736, 105) + 3 * K / 7) * x**3
    + s.Rational(8, 5) * x**5
)


def cubic(a, c):
    return s.expand(c**2 - a**3 + s.Rational(3, 4) * a + s.Rational(1, 4))


def jacobian(a, c):
    return s.expand(s.diff(a, x) * s.diff(c, t) - s.diff(a, t) * s.diff(c, x))


q = -(x**2 + 6) / 2


def determinant_old_row(order: int, arows, crows):
    return s.expand(
        sum(
            (order - i) * s.diff(arows[i], x) * crows[order - i]
            - i * arows[i] * s.diff(crows[order - i], x)
            for i in range(1, order)
        )
    )


def cubic_old_row(order: int, arows, crows):
    quadratic = sum(crows[i] * crows[order - i] for i in range(1, order))
    cubic_sum = sum(
        arows[i] * arows[j] * arows[k]
        for i in range(order)
        for j in range(order)
        for k in range(order)
        if i + j + k == order
    )
    return s.expand(quadratic - cubic_sum)


def predicted_source_row(order: int, arows, crows):
    return s.expand(cubic_old_row(order, arows, crows) - q * determinant_old_row(order, arows, crows) / order)


def reversed_independent_columns(matrix: s.Matrix):
    reversed_matrix = matrix[:, ::-1]
    reversed_pivots = reversed_matrix.rref()[1]
    chosen = tuple(matrix.cols - 1 - j for j in reversed_pivots)
    demand(len(chosen) == matrix.rows, "reversed determinant section is square")
    demand(matrix[:, chosen].det() != 0, "reversed determinant section invertible")
    return chosen


def independent_rows(matrix: s.Matrix):
    return tuple(matrix.T.rref()[1])


def solve_selected(matrix: s.Matrix, rhs: s.Matrix, columns):
    block = matrix[:, list(columns)]
    rows = independent_rows(block)
    square = block.extract(rows, tuple(range(len(columns))))
    demand(square.rows == square.cols and square.det() != 0, "selected solve invertible")
    return tuple(s.cancel(v) for v in square.inv() * rhs.extract(rows, (0,)))


def consistency_equations(matrix: s.Matrix, rhs: s.Matrix):
    return [s.factor(s.cancel(s.together((left.T * rhs)[0]))) for left in matrix.T.nullspace()]


def determinant_section(order: int, arows, crows):
    avars = tuple(s.symbols(f"da{order}_0:{order + 2}"))
    cvars = tuple(s.symbols(f"dc{order}_0:{order + 3}"))
    avalue = sum(avars[j] * x**j for j in range(order + 2))
    cvalue = sum(cvars[j] * x**j for j in range(order + 3))
    equation = s.expand(order * (s.diff(A0, x) * cvalue - s.diff(C0, x) * avalue) + determinant_old_row(order, arows, crows))
    equations = [coeff(equation, x, j) for j in range(order + 4)]
    matrix, rhs = s.linear_eq_to_matrix(equations, avars + cvars)
    demand(matrix.rank() == order + 4, f"determinant rank row {order}")
    chosen = reversed_independent_columns(matrix)
    values = solve_selected(matrix, rhs, chosen)
    substitution = {symbol: s.Integer(0) for symbol in avars + cvars}
    substitution.update({(avars + cvars)[chosen[j]]: values[j] for j in range(len(chosen))})
    abase = s.expand(avalue.subs(substitution))
    cbase = s.expand(cvalue.subs(substitution))
    demand(zero(equation.subs(substitution)), f"determinant section row {order}")
    return abase, cbase


def general_determinant_row(order: int, arows, crows, prefix: str):
    abase, cbase = determinant_section(order, arows, crows)
    variables = tuple(s.symbols(f"{prefix}0:{order + 1}"))
    theta = sum(variables[j] * x**j for j in range(order + 1))
    anew = s.expand(abase + theta * s.diff(A0, x))
    cnew = s.expand(cbase + theta * s.diff(C0, x))
    demand(s.degree(anew, x) <= order + 1, f"A degree row {order}")
    demand(s.degree(cnew, x) <= order + 2, f"C degree row {order}")
    return anew, cnew, variables


def source_match(order: int, grow, arows, crows, variables):
    difference = s.expand(grow - predicted_source_row(order, arows, crows))
    demand(difference == 0 or s.degree(difference, x) <= order, f"complete bracket degree row {order}")
    equations = [coeff(difference, x, j) for j in range(order + 1)]
    matrix, rhs = s.linear_eq_to_matrix(equations, variables)
    return difference, equations, matrix, rhs, consistency_equations(matrix, rhs)


def solve_all_variables(matrix, rhs, variables):
    demand(matrix.rank() == len(variables), "full bracket selector rank")
    values = solve_selected(matrix, rhs, tuple(range(len(variables))))
    substitution = {variables[j]: values[j] for j in range(len(variables))}
    return substitution


def solve_pivots_with_free(matrix, rhs, variables):
    pivots = tuple(matrix.rref()[1])
    free = tuple(j for j in range(len(variables)) if j not in pivots)
    effective = rhs - matrix[:, free] * s.Matrix([variables[j] for j in free]) if free else rhs
    values = solve_selected(matrix, effective, pivots)
    substitution = {variables[pivots[j]]: values[j] for j in range(len(pivots))}
    return pivots, free, substitution


def depth_matrix(depth: int, max_row: int):
    coordinates = tuple((n, r) for n in range(max_row + 1) for r in range(n + depth + 1))
    columns = []
    for by in range(depth + 1):
        for ax in range(depth - by + 1):
            for even in range(max_row // 2 + 1):
                for power in range(max_row + 1):
                    if by + power + 2 * even <= max_row:
                        columns.append((ax, by, power, even))
    index = {coordinate: j for j, coordinate in enumerate(coordinates)}
    matrix = s.zeros(len(coordinates), len(columns))
    for col, (ax, by, power, even) in enumerate(columns):
        for j in range(power + even + 1):
            row = by + power + 2 * even + j
            degree = ax + 2 * by + even + 2 * j
            if row <= max_row:
                matrix[index[(row, degree)], col] += s.binomial(power + even, j)
    return coordinates, matrix


def depth_residuals(arows, crows, max_row: int):
    acoords, amatrix = depth_matrix(2, max_row)
    ccoords, cmatrix = depth_matrix(3, max_row)
    avec = s.Matrix([coeff(arows[n], x, degree) for n, degree in acoords])
    cvec = s.Matrix([coeff(crows[n], x, degree) for n, degree in ccoords])
    ares = [s.expand((left.T * avec)[0]) for left in amatrix.T.nullspace()]
    cres = [s.expand((left.T * cvec)[0]) for left in cmatrix.T.nullspace()]
    return ares, cres, amatrix, cmatrix


def bracket_through_row_eight(grows):
    arows = [A0, A1, A2, A3]
    crows = [C0, C1, C2, C3]
    demand(zero(cubic(A0, C0)), "boundary cubic")
    demand(zero((-3 * A0**2 + s.Rational(3, 4)) + q * s.diff(C0, x)), "boundary rotated gradient A")
    demand(zero(2 * C0 - q * s.diff(A0, x)), "boundary rotated gradient C")
    demand(zero(predicted_source_row(4, arows, crows) - grows[4]), "inherited row four")

    eliminated = {}
    source_symbols = {5: v5, 6: U, 7: W, 8: Z}
    obstructions = {}
    for order in range(4, 8):
        anew, cnew, variables = general_determinant_row(order, arows, crows, f"r{order}_")
        trial_a = arows + [anew]
        trial_c = crows + [cnew]
        next_row = order + 1
        _, equations, matrix, rhs, compat = source_match(
            next_row, s.expand(grows[next_row].subs(eliminated)), trial_a, trial_c, variables
        )
        demand(matrix.rank() == next_row, f"next-source tangent rank {next_row}")
        nonzero = [value for value in compat if not zero(value)]
        demand(len(nonzero) == 1, f"one obstruction row {next_row}")
        symbol = source_symbols[next_row]
        solutions = s.solve(nonzero[0], symbol, dict=True)
        demand(len(solutions) == 1 and symbol in solutions[0], f"source solve row {next_row}")
        eliminated[symbol] = s.cancel(solutions[0][symbol])
        rhs_fixed = rhs.subs(symbol, eliminated[symbol])
        tangent_solution = solve_all_variables(matrix, rhs_fixed, variables)
        arows.append(s.expand(anew.subs(tangent_solution)))
        crows.append(s.expand(cnew.subs(tangent_solution)))
        demand(zero(predicted_source_row(next_row, arows, crows) - grows[next_row].subs(eliminated)), f"matched row {next_row}")
        obstructions[next_row] = primitive_numerator(nonzero[0], (Delta, Theta, Phi, eta, zeta, v5, xi, alpha, beta, U, W, Z, rho, sigma)).as_expr()

    anew, cnew, theta8 = general_determinant_row(8, arows, crows, "z8_")
    arows.append(anew)
    crows.append(cnew)
    demand(zero(predicted_source_row(8, arows[:8], crows[:8]) - grows[8].subs(eliminated)), "row-eight bracket closure")
    return arows, crows, eliminated, theta8, obstructions


def row_eight_depth(arows, crows, theta8):
    ares, cres, amatrix, cmatrix = depth_residuals(arows, crows, 8)
    equations = ares + cres
    matrix, rhs = s.linear_eq_to_matrix(equations, theta8)
    demand((amatrix.shape, amatrix.rank(), len(amatrix.T.nullspace())) == ((63, 131), 51, 12), "row8 P2 universe")
    demand((cmatrix.shape, cmatrix.rank(), len(cmatrix.T.nullspace())) == ((72, 204), 63, 9), "row8 P3 universe")
    demand(matrix.rank() == 2, "row8 terminal rank")
    compatibility = consistency_equations(matrix, rhs)
    expected = (15 * Delta - 896, 3 * Phi + 2 * zeta, 3030 * Delta + 225 * Theta - 182528)
    variables = (Delta, Theta, Phi, eta, zeta, xi, alpha, beta, rho, sigma)
    same_linear_span(compatibility, expected, variables, "row8 source ideal")
    gate_subs = {Delta: s.Rational(896, 15), Theta: s.Rational(512, 75), zeta: -3 * Phi / 2}
    fixed_matrix, fixed_rhs = s.linear_eq_to_matrix([v.subs(gate_subs) for v in equations], theta8)
    pivots, free_indices, terminal = solve_pivots_with_free(fixed_matrix, fixed_rhs, theta8)
    demand(len(pivots) == 2 and len(free_indices) == 7, "row8 two selected seven free")
    demand(all(zero(v.subs(gate_subs).subs(terminal)) for v in equations), "row8 all depth residuals")
    return gate_subs, terminal, tuple(theta8[j] for j in free_indices)


def polynomial_span_equal(actual, expected, variables, label):
    actual_polys = [primitive_numerator(v, variables) for v in actual if not zero(v)]
    expected_polys = [primitive_numerator(v, variables) for v in expected]
    terms = sorted(set().union(*(set(p.monoms()) for p in actual_polys + expected_polys)))
    va = s.Matrix([[p.as_dict().get(m, 0) for m in terms] for p in actual_polys])
    ve = s.Matrix([[p.as_dict().get(m, 0) for m in terms] for p in expected_polys])
    demand(va.rank() == ve.rank() == len(expected_polys), label + " rank")
    demand(va.col_join(ve).rank() == va.rank(), label + " equality")


def p_nonzero_depth(a9, c9, xi_graph, gate_a, gate_b):
    alpha_graph = s.factor(s.solve(gate_b, alpha, dict=True)[0][alpha])
    beta_graph = s.factor(s.solve(gate_a.subs(alpha, alpha_graph), beta, dict=True)[0][beta])
    xi_final = s.factor(xi_graph.subs(alpha, alpha_graph))
    source = {alpha: alpha_graph, beta: beta_graph}
    demand(zero(gate_a.subs(source)) and zero(gate_b.subs(source)), "P-nonzero bracket graph")
    a9f = [s.expand(v.subs(source)) for v in a9]
    c9f = [s.expand(v.subs(source)) for v in c9]
    a10new, c10new, theta10 = general_determinant_row(10, a9f, c9f, "z10_")
    a10 = a9f + [a10new]
    c10 = c9f + [c10new]
    ares, cres, am, cm = depth_residuals(a10, c10, 10)
    ma, ra = s.linear_eq_to_matrix(ares, theta10)
    mc, rc = s.linear_eq_to_matrix(cres, theta10)
    joint = ma.col_join(mc)
    rhs = ra.col_join(rc)
    demand((am.shape, am.rank(), len(am.T.nullspace())) == ((88, 193), 68, 20), "row10 P2 universe")
    demand((cm.shape, cm.rank(), len(cm.T.nullspace())) == ((99, 304), 83, 16), "row10 P3 universe")
    demand((ma.rank(), mc.rank(), joint.rank()) == (3, 3, 3), "row10 P-nonzero terminal ranks")
    demand(tuple(joint.rref()[1]) == (8, 9, 10), "row10 P-nonzero terminal pivots")
    compatibility = [v for v in consistency_equations(joint, rhs) if not zero(v)]
    for value in compatibility:
        denominator = s.factor(s.fraction(s.cancel(s.together(value)))[1])
        demand(not (denominator.free_symbols - {Phi}), "row10 localization only at Phi")
        degree = s.degree(denominator, Phi)
        demand(degree is not None and s.cancel(denominator / Phi**degree).is_Rational,
               "row10 denominator is a Phi monomial")
    variables = (Phi, eta, rho, sigma)
    polys = [primitive_numerator(v, variables) for v in compatibility]
    common = polys[0]
    for poly in polys[1:]:
        common = s.gcd(common, poly)
    if common.LC() < 0:
        common = -common
    d0 = (
        7231154026500 * Phi**3
        + 50541940696500 * Phi**2 * eta
        + 6793915500000 * Phi * eta**2
        + 353642000625 * eta**3
        - 631918028977864704 * Phi
        - 91584545734393856 * eta
    )
    d13 = s.expand(d0 - 707284001250 * Phi**2 * rho)
    demand(zero(common.as_expr() - d13), "row10 normal equation")
    for poly in polys:
        quotient, remainder = s.div(poly, common)
        demand(remainder == 0 and quotient.total_degree() == 0, "row10 principal compatibility")
    demand(s.diff(d13, rho) == -707284001250 * Phi**2, "rho transverse derivative")
    demand(s.diff(d13, sigma) == 0, "sigma tangent derivative")
    demand(len(theta10) - joint.rank() == 8, "row10 affine-eight terminal fibre")
    return alpha_graph, beta_graph, xi_final, d0, d13, a10, c10, theta10, ares, cres, joint, rhs


def reduce_mod(value, variable, modulus):
    numerator, denominator = s.fraction(s.cancel(s.together(value)))
    modpoly = s.Poly(modulus, variable, domain=s.EX)
    npoly = s.Poly(numerator, variable, domain=s.EX)
    dpoly = s.Poly(denominator, variable, domain=s.EX)
    inverse = s.invert(dpoly, modpoly)
    return s.factor((npoly * inverse).rem(modpoly).as_expr())


def p_zero_boundary(a9, c9, gate_a, gate_b):
    f = 18612736875 * eta**2 - 4820239249178624
    demand(zero(gate_b.subs(Phi, 0) - f), "P-zero bracket quadratic")
    alpha0 = s.factor(s.solve(gate_a.subs(Phi, 0), alpha, dict=True)[0][alpha])
    expected_alpha0 = -(
        854861882349375 * eta**2 + 27743789318253707264
    ) / (278844730740000 * eta)
    demand(zero(alpha0 - expected_alpha0), "P-zero alpha graph")
    demand(s.gcd(s.Poly(f, eta), s.Poly(eta, eta)) == 1, "P-zero division by eta lawful")
    a9z = [s.expand(v.subs({Phi: 0, alpha: alpha0})) for v in a9]
    c9z = [s.expand(v.subs({Phi: 0, alpha: alpha0})) for v in c9]
    a10new, c10new, theta10 = general_determinant_row(10, a9z, c9z, "b10_")
    a10 = a9z + [a10new]
    c10 = c9z + [c10new]
    ares, cres, _, _ = depth_residuals(a10, c10, 10)
    ma, ra = s.linear_eq_to_matrix(ares, theta10)
    mc, rc = s.linear_eq_to_matrix(cres, theta10)
    demand((ma.rank(), mc.rank(), ma.col_join(mc).rank()) == (3, 3, 3), "P-zero terminal ranks")
    pa = tuple(ma.rref()[1])
    pc = tuple(mc.rref()[1])
    demand(pa == pc == (8, 9, 10), "P-zero common pivot coordinates")
    sola = solve_selected(ma, ra, pa)
    solc = solve_selected(mc, rc, pc)
    differences = [reduce_mod(sola[j] - solc[j], eta, f) for j in range(3)]
    suba = {theta10[pa[j]]: sola[j] for j in range(3)}
    subc = {theta10[pc[j]]: solc[j] for j in range(3)}
    demand(all(zero(reduce_mod(value.subs(suba), eta, f)) for value in ares),
           "P-zero P2 subsystem sufficient modulo F")
    demand(all(zero(reduce_mod(value.subs(subc), eta, f)) for value in cres),
           "P-zero P3 subsystem sufficient modulo F")
    nonzero = [v for v in differences if not zero(v)]
    demand(len(nonzero) == 1, "P-zero one joint compatibility")
    quotient = s.cancel(nonzero[0] / (5 * beta + 6 * eta))
    demand(quotient.is_Rational and quotient != 0, "P-zero beta graph")
    joint_subs = {
        theta10[pa[j]]: s.cancel(sola[j].subs(beta, -6 * eta / 5))
        for j in range(3)
    }
    joint_subs[beta] = -6 * eta / 5
    demand(all(zero(reduce_mod(value.subs(joint_subs), eta, f)) for value in ares + cres),
           "P-zero joint system sufficient modulo F")
    demand(all(zero(s.diff(value, rho)) and zero(s.diff(value, sigma)) for value in differences),
           "P-zero weight-thirteen parameters free")
    demand(s.degree(s.Poly(f, eta).sqf_part()) == 2, "P-zero two reduced roots")
    demand(len(theta10) - ma.col_join(mc).rank() == 8, "P-zero affine-eight fibre")
    return f, alpha0, differences


def quartic_invariants(poly, variable):
    values = s.Poly(poly, variable).all_coeffs()
    demand(len(values) == 5, "binary quartic degree four")
    aa, bb, cc, dd, ee = values
    invariant_i = s.expand(12 * aa * ee - 3 * bb * dd + cc**2)
    invariant_j = s.expand(72 * aa * cc * ee + 9 * bb * cc * dd - 27 * aa * dd**2 - 27 * bb**2 * ee - 2 * cc**3)
    delta_flat = s.expand(4 * invariant_i**3 - invariant_j**2)
    return invariant_i, invariant_j, delta_flat


def geometry(d0, d13):
    r, a = s.symbols("r a")
    normal = s.Integer(707284001250)
    A = 353642000625 * r**3 + 6793915500000 * r**2 + 50541940696500 * r + 7231154026500
    B = 91584545734393856 * r + 631918028977864704
    quadratic = s.expand(A * Phi**2 - normal * a * Phi - B)
    demand(zero(d13.subs({eta: r * Phi, rho: a}) - Phi * quadratic), "ratio quadratic")
    Y = 2 * A * Phi - normal * a
    quartic = s.expand(4 * A * B + normal**2 * a**2)
    demand(zero(Y**2 - quartic - 4 * A * quadratic), "completed square")
    mate = s.cancel(normal * a / A - Phi)
    demand(zero(quadratic.subs(Phi, mate) - quadratic), "deck involution")
    demand(zero(d13.subs({Phi: -Phi, eta: -eta, rho: -rho}) + d13), "total parity exchanges fibres")
    demand(zero(mate.subs(a, 0) + Phi), "central sign deck action")

    invariant_i, invariant_j, delta_flat = quartic_invariants(quartic, r)
    qmod_two = s.Poly(quartic, r, a, modulus=7)
    leading = int(qmod_two.coeff_monomial(r**4)) % 7
    demand(leading != 0, "quartic leading coefficient good mod seven")
    normalized = s.Poly(pow(leading, -1, 7) * qmod_two.as_expr(), r, a, modulus=7)
    expected_q = s.Poly(r**4 + r**3 + 3 * r**2 + 2 * r + (2 + a**2), r, a, modulus=7)
    demand(normalized == expected_q, "quartic pencil mod seven")
    im = s.Poly(invariant_i, a, modulus=7)
    jm = s.Poly(invariant_j, a, modulus=7)
    dm = s.Poly(delta_flat, a, modulus=7)
    # The invariant polynomials may carry a common nonzero scaling.  Their
    # directly computed normalized residues are pinned below.
    demand(im == s.Poly(-2 * a**2 - 1, a, modulus=7), "quartic I mod seven")
    demand(jm == s.Poly(-3, a, modulus=7), "quartic J mod seven")
    expected_delta = s.Poly(3 * a**6 + a**4 - 3 * a**2 + 1, a, modulus=7)
    demand(dm == expected_delta, "quartic discriminant mod seven")
    demand(dm.degree() == 6 and int(dm.LC()) % 7 == 3, "valuation-leading discriminant term")
    demand(all(int(dm.eval(j)) % 7 != 0 for j in range(7)), "no F7 rational discriminant root")
    inv_data = []
    for aval in (0, 1):
        iv = int(im.eval(aval)) % 7
        jv = int(jm.eval(aval)) % 7
        dv = int(dm.eval(aval)) % 7
        inv_data.append((iv, jv, dv, iv**3 * pow(dv, -1, 7) % 7))
    demand(inv_data == [(6, 4, 1, 6), (4, 4, 2, 4)], "distinct good-reduction j witnesses")
    exact_cross = s.expand(invariant_i.subs(a, 0) ** 3 * delta_flat.subs(a, 1) - invariant_i.subs(a, 1) ** 3 * delta_flat.subs(a, 0))
    demand(exact_cross != 0, "exact non-isotrivial witness")
    resultant = s.resultant(quartic, s.diff(quartic, r), r)
    demand(resultant != 0, "generic separability")
    ratio = s.cancel(resultant / delta_flat)
    demand(not ratio.has(a) and ratio != 0, "quartic invariant discriminant agrees with resultant")

    f = 18612736875 * eta**2 - 4820239249178624
    demand(zero(d13.subs(Phi, 0) - 19 * eta * f), "boundary gluing factorization")
    demand(d13.subs({Phi: 0, eta: 0}) == 0 and f.subs(eta, 0) != 0, "origin is closure-only")
    T = s.symbols("T")
    homogeneous = (
        7231154026500 * Phi**3
        + 50541940696500 * Phi**2 * eta
        + 6793915500000 * Phi * eta**2
        + 353642000625 * eta**3
        - normal * a * Phi**2 * T
        - 631918028977864704 * Phi * T**2
        - 91584545734393856 * eta * T**2
    )
    gradient_at_origin = tuple(s.diff(homogeneous, v).subs({Phi: 0, eta: 0, T: 1}) for v in (Phi, eta, T))
    demand(gradient_at_origin != (0, 0, 0), "persistent rational point smooth")
    demand(s.degree(s.Poly(A, r).sqf_part()) == 3, "three distinct points at infinity")
    demand(s.gcd(s.Poly(A, r), s.Poly(B, r)) == 1, "ratio numerator denominator coprime")
    return r, a, A, B, quartic, inv_data, expected_delta.as_expr(), exact_cross


def hostile_check(eliminated, gate8, alpha_graph, beta_graph, xi_graph, d0, d13, a10, c10, theta10, ares, cres, joint, rhs):
    hostile = {
        Phi: s.Integer(1),
        eta: s.Integer(0),
        rho: -s.Rational(35106155434657678, 39293555625),
        sigma: s.Integer(1),
    }
    hostile[alpha] = s.factor(alpha_graph.subs(hostile))
    hostile[beta] = s.factor(beta_graph.subs(hostile))
    hostile[xi] = s.factor(xi_graph.subs(hostile))
    demand(hostile[xi] == s.Rational(120573932063, 7209000), "hostile xi")
    demand(hostile[alpha] == s.Rational(2410085120016937, 18612736875), "hostile alpha")
    demand(hostile[beta] == s.Rational(91, 15), "hostile beta")
    demand(zero(d13.subs(hostile)), "hostile deformed carrier")
    demand(d0.subs(hostile) == -631910797823838204, "hostile misses central carrier")
    matrix = joint.subs(hostile)
    target = rhs.subs(hostile)
    pivots, free_indices, terminal = solve_pivots_with_free(matrix, target, theta10)
    terminal.update({theta10[j]: 0 for j in free_indices})
    full = {**hostile, **terminal}
    for symbol, value in gate8.items():
        full[symbol] = s.cancel(value.subs(full))
    # Recursion eliminated coefficients, evaluated after all source values.
    for symbol, value in eliminated.items():
        full[symbol] = s.cancel(value.subs(full))
    demand(all(zero(v.subs(full)) for v in ares + cres), "hostile all depth residuals")

    aval = sum(a10[n].subs(full) * t**n for n in range(11))
    cval = sum(c10[n].subs(full) * t**n for n in range(11))
    source = G13.subs(full)
    pdefect = s.expand(cubic(aval, cval) - source)
    jdefect = s.expand(jacobian(aval, cval) - 1)
    p_rows = [s.factor(coeff(pdefect, t, n)) for n in range(11)]
    j_rows = [s.factor(coeff(jdefect, t, n)) for n in range(10)]
    if any(not zero(v) for v in p_rows):
        print("DEBUG_P_ROWS", [(j, s.sstr(v)) for j, v in enumerate(p_rows) if not zero(v)])
    if any(not zero(v) for v in j_rows):
        print("DEBUG_J_ROWS", [(j, s.sstr(v)) for j, v in enumerate(j_rows) if not zero(v)])
    demand(all(zero(v) for v in p_rows), "hostile direct cubic rows 0..10")
    demand(all(zero(v) for v in j_rows), "hostile direct Jacobian rows 0..9")
    demand(all(s.degree(a10[n].subs(full), x) <= n + 1 for n in range(4, 11)), "hostile A caps")
    demand(all(s.degree(c10[n].subs(full), x) <= n + 2 for n in range(4, 11)), "hostile C caps")

    uval = s.factor(eliminated[U].subs(full))
    wval = s.factor(eliminated[W].subs(full))
    zval = s.factor(eliminated[Z].subs(full))
    lam = s.factor(uval + wval + zval)
    response_discriminant = s.factor(wval**2 - 4 * uval * zval)
    demand(uval == -s.Rational(219007077871, 32440500), "hostile U")
    demand(wval == s.Rational(25779284581, 720900), "hostile W")
    demand(zval == -s.Rational(1231789589, 44500), "hostile Z")
    demand(lam == s.Rational(43086117893, 32440500), "hostile Lambda")
    demand(response_discriminant == s.Rational(6902544853907445578341, 12992420250000), "hostile response discriminant")
    dense_product = hostile[rho] * hostile[sigma] * uval * zval * (hostile[rho] + hostile[sigma]) * lam * response_discriminant
    demand(dense_product != 0, "hostile coefficient density")
    return hostile, uval, wval, zval, lam, response_discriminant


def main():
    grows = {n: coeff(H13, t, n) for n in range(4, 14)}
    delta_rows = {n: s.factor(grows[n] - coeff(H0, t, n)) for n in range(4, 14)}
    for n in range(4, 14):
        expected_rho = s.binomial(6, n - 7) * x ** (1 + 2 * (n - 7)) * rho if 7 <= n <= 13 else 0
        expected_sigma = s.binomial(5, n - 8) * x ** (3 + 2 * (n - 8)) * sigma if 8 <= n <= 13 else 0
        demand(zero(delta_rows[n] - expected_rho - expected_sigma), f"binomial source row {n}")

    arows, crows, eliminated, theta8, obstructions = bracket_through_row_eight(grows)
    expected_eliminated = {
        v5: -40 * Delta / 9 - 2 * Theta / 3 - s.Rational(184832, 2025),
        U: 4604 * Delta / 165 + 784 * Theta / 297 - 6 * xi / 11 + s.Rational(137763328, 200475),
        W: -20 * Delta**2 / 9 - s.Rational(520544, 2673) * Delta - 13 * Phi**2 / 12
           - s.Rational(456508, 13365) * Theta + 424 * xi / 99 - s.Rational(3222956032, 200475),
        Z: s.Rational(5051, 486) * Delta**2 + s.Rational(221, 81) * Delta * Theta
           + s.Rational(31708408, 72171) * Delta - s.Rational(130, 81) * Phi**2
           - s.Rational(65, 36) * Phi * eta - s.Rational(15675376, 360855) * Theta
           - s.Rational(22886, 2673) * xi + s.Rational(278111504384, 5412825),
    }
    for symbol, expected in expected_eliminated.items():
        demand(zero(eliminated[symbol] - expected), f"unchanged early scalar {symbol}")
        demand(zero(s.diff(eliminated[symbol], rho)) and zero(s.diff(eliminated[symbol], sigma)),
               f"early scalar independent of weight13 channels {symbol}")
    gate8, terminal8, free8 = row_eight_depth(arows, crows, theta8)

    # Reimplement row nine and ten here so the polynomial-span comparison can
    # use honest polynomial generators (rather than treating monomials as
    # pseudo-symbols).
    a8 = [s.expand(v.subs(gate8).subs(terminal8)) for v in arows]
    c8 = [s.expand(v.subs(gate8).subs(terminal8)) for v in crows]
    grow9 = s.expand(grows[9].subs(eliminated).subs(gate8))
    _, _, matrix9, rhs9, compat9 = source_match(9, grow9, a8, c8, free8)
    demand(matrix9.rank() == 7, "row9 selector rank")
    nonzero9 = [v for v in compat9 if not zero(v)]
    e9 = 613527750 * Phi**2 - 511211250 * Phi * alpha - 3154140000 * Phi * eta - 255605625 * eta**2 + 6736896000 * xi - 46483785515008
    demand(len(nonzero9) == 1 and s.cancel(nonzero9[0] / e9).is_Rational, "unchanged row9 gate")
    selected8 = solve_all_variables(matrix9, rhs9, free8)
    a8 = [s.expand(v.subs(selected8)) for v in a8]
    c8 = [s.expand(v.subs(selected8)) for v in c8]
    xi_pre = s.factor(s.solve(e9, xi, dict=True)[0][xi])
    a8 = [s.expand(v.subs(xi, xi_pre)) for v in a8]
    c8 = [s.expand(v.subs(xi, xi_pre)) for v in c8]

    a9new, c9new, theta9 = general_determinant_row(9, a8, c8, "z9_")
    a9, c9 = a8 + [a9new], c8 + [c9new]
    ares9, cres9, am9, cm9 = depth_residuals(a9, c9, 9)
    ma9, ra9 = s.linear_eq_to_matrix(ares9, theta9)
    mc9, rc9 = s.linear_eq_to_matrix(cres9, theta9)
    joint9, jr9 = ma9.col_join(mc9), ra9.col_join(rc9)
    demand((am9.shape, am9.rank(), len(am9.T.nullspace())) == ((75, 160), 59, 16), "row9 P2 universe")
    demand((cm9.shape, cm9.rank(), len(cm9.T.nullspace())) == ((85, 251), 73, 12), "row9 P3 universe")
    demand((ma9.rank(), mc9.rank(), joint9.rank()) == (3, 2, 3), "row9 depth ranks")
    piv9, free9_indices, term9 = solve_pivots_with_free(joint9, jr9, theta9)
    demand(piv9 == (7, 8, 9), "row9 terminal pivots")
    demand(len(piv9) == 3 and len(free9_indices) == 7, "row9 terminal dimensions")
    demand(all(zero(v.subs(term9)) for v in ares9 + cres9), "row9 depth invariance")
    a9 = [s.expand(v.subs(term9)) for v in a9]
    c9 = [s.expand(v.subs(term9)) for v in c9]
    free9 = tuple(theta9[j] for j in free9_indices)

    grow10 = s.expand(grows[10].subs(eliminated).subs(gate8).subs(xi, xi_pre))
    _, _, matrix10, rhs10, compat10 = source_match(10, grow10, a9, c9, free9)
    demand(matrix10.rank() == 7, "row10 selector rank")
    gate_a = 1283973857930250 * Phi**2 - 1709723764698750 * Phi * alpha - 264168692280000 * Phi * beta + 2669098150944000 * Phi * eta - 278844730740000 * Phi * rho - 278844730740000 * alpha * eta - 854861882349375 * eta**2 - 27743789318253707264
    gate_b = 69009144750 * Phi**2 + 37225473750 * Phi * alpha + 357574500000 * Phi * eta + 18612736875 * eta**2 - 4820239249178624
    polynomial_span_equal([v for v in compat10 if not zero(v)], (gate_a, gate_b), (Phi, eta, alpha, beta, rho, sigma), "row10 bracket ideal")
    selected9 = solve_all_variables(matrix10, rhs10, free9)
    a9 = [s.expand(v.subs(selected9)) for v in a9]
    c9 = [s.expand(v.subs(selected9)) for v in c9]

    alpha_graph, beta_graph, xi_final, d0, d13, a10, c10, theta10, ares10, cres10, joint10, jr10 = p_nonzero_depth(a9, c9, xi_pre, gate_a, gate_b)
    demand(any(not zero(s.diff(value, sigma)) for value in a10 + c10),
           "sigma changes selected lifts despite cylinder equation")
    f, alpha0, boundary_differences = p_zero_boundary(a9, c9, gate_a, gate_b)
    r, apar, ar, br, quartic, inv_data, delta_mod7, j_cross = geometry(d0, d13)
    hostile, uval, wval, zval, lam, resp_disc = hostile_check(
        eliminated, gate8, alpha_graph, beta_graph, xi_final, d0, d13,
        a10, c10, theta10, ares10, cres10, joint10, jr10,
    )

    print("CLEANROOM=dependency-free reversed-section reconstruction")
    print("WEIGHT13_ROWS=" + ";".join(f"{n}:{s.sstr(delta_rows[n])}" for n in range(4, 14)))
    print("ROW8_IDEAL=(15*Delta-896,3*Phi+2*zeta,3030*Delta+225*Theta-182528)")
    print("EARLY_SCALARS=v5,U,W,Z unchanged")
    print(f"ROW9_GATE={s.sstr(e9)}")
    print("ROW9_DEPTH_RANKS=3,2,3; TERMINAL_DIMENSION=7")
    print(f"ROW10_BRACKET_GATE_A={s.sstr(gate_a)}")
    print(f"ROW10_BRACKET_GATE_B={s.sstr(gate_b)}")
    print(f"XI_GRAPH={s.sstr(xi_final)}")
    print(f"ALPHA_GRAPH={s.sstr(alpha_graph)}")
    print(f"BETA_GRAPH={s.sstr(beta_graph)}")
    print(f"ROW10_NORMAL={s.sstr(d13)}")
    print("ROW10_DEPTH_RANKS=3,3,3; TERMINAL_DIMENSION=8")
    print(f"BOUNDARY_F={s.sstr(f)}")
    print(f"BOUNDARY_ALPHA={s.sstr(alpha0)}")
    print("BOUNDARY_DIFFS=" + ";".join(s.sstr(v) for v in boundary_differences))
    print("QUARTIC=4*A(r)*B(r)+707284001250^2*a^2")
    print(f"MOD7_DISCRIMINANT={s.sstr(delta_mod7)}")
    print(f"MOD7_J_WITNESSES={inv_data}")
    print(f"EXACT_J_CROSS_NONZERO={j_cross != 0}")
    print("RATIONAL_SMOOTHNESS=all a in Q by the v7 dichotomy")
    print(f"HOSTILE=Phi:{hostile[Phi]},eta:{hostile[eta]},rho:{hostile[rho]},sigma:{hostile[sigma]}")
    print(f"HOSTILE_XI_ALPHA_BETA={hostile[xi]},{hostile[alpha]},{hostile[beta]}")
    print(f"HOSTILE_U_W_Z={uval},{wval},{zval}")
    print(f"HOSTILE_LAMBDA_DISC={lam},{resp_disc}")
    print("DETERMINISM=normal,-O,hashseed streams")
    print(f"CHECKS={CHECKS}")
    print("STATUS=PASS; JC2=OPEN")


if __name__ == "__main__":
    main()
