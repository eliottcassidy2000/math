#!/usr/bin/env python3
"""Private exact discovery scout for split degrees at least 26."""

from __future__ import annotations

import time

import sympy as sp


def coefficients(degree: int, p: sp.Expr, q: sp.Expr, r: sp.Expr, extra: int = 3):
    exponent = sp.Rational(degree, 4)
    cs = [sp.Integer(1)]
    quartic = {2: p, 3: q, 4: r}
    for index in range(1, degree + extra + 1):
        cs.append(sp.factor(sum(
            quartic[step]
            * ((exponent+1)*step-index)
            * cs[index-step]
            for step in range(2, min(4, index)+1)
        )/index))
    return cs


def observables(degree: int, d: sp.Symbol, q: sp.Symbol, s: sp.Symbol):
    cs = coefficients(degree, 2*d, q, d**2-s)
    return (
        sp.factor(4*cs[degree+1]),
        sp.factor(4*cs[degree+2]),
        sp.factor(4*cs[degree+3]+2*d*cs[degree+1]),
        cs,
    )


def invariant_form(
    expression: sp.Expr,
    d: sp.Symbol,
    q: sp.Symbol,
    s: sp.Symbol,
    rho: sp.Symbol,
    z: sp.Symbol,
):
    transformed = sp.Poly(sp.expand(expression.subs({q: 1, d: rho*s**2})), rho, s)
    minimum = min(monomial[1] for monomial, _ in transformed.terms())
    residues = {(monomial[1]-minimum) % 3 for monomial, _ in transformed.terms()}
    if residues != {0}:
        raise RuntimeError(f"noninvariant exponents: min={minimum}, residues={residues}")
    result = sum(
        coefficient*rho**monomial[0]*z**((monomial[1]-minimum)//3)
        for monomial, coefficient in transformed.terms()
    )
    return minimum, sp.factor(sp.Poly(result, rho, z).primitive()[1].as_expr())


def unit_ideal(expressions: list[sp.Expr], *variables: sp.Symbol) -> bool:
    return sp.groebner(expressions, *variables, order="grevlex").contains(sp.Integer(1))


d, q, s = sp.symbols("d q s")
rho, z, omega = sp.symbols("rho z omega")

for degree in (50, 54, 58, 62, 66, 70, 74):
    started = time.time()
    phi, psi, response, cs = observables(degree, d, q, s)
    jacobian = sp.factor(sp.det(sp.Matrix([
        [sp.diff(phi.subs(q, 1), d), sp.diff(phi.subs(q, 1), s)],
        [sp.diff(psi.subs(q, 1), d), sp.diff(psi.subs(q, 1), s)],
    ])))
    f_power, f = invariant_form(phi, d, q, s, rho, z)
    g_power, g = invariant_form(psi, d, q, s, rho, z)
    t_power, t = invariant_form(response, d, q, s, rho, z)
    j_power, jac = invariant_form(jacobian, d, q, s, rho, z)
    print(f"\nM={degree};k={(degree+2)//4}")
    print(f"s_powers=F{f_power},G{g_power},R{t_power},J{j_power}")
    print(f"degrees=F{sp.Poly(f,rho,z).total_degree()},G{sp.Poly(g,rho,z).total_degree()},"
          f"R{sp.Poly(t,rho,z).total_degree()},J{sp.Poly(jac,rho,z).total_degree()}")
    basis = sp.groebner([f, g], rho, z, order="lex")
    print(f"zero_dimensional={basis.is_zero_dimensional};basis_length={len(basis.polys)}")
    if len(basis.polys) <= 3:
        for item in basis.polys:
            expression = item.as_expr()
            print(f"basis_poly_degrees_rho_z={sp.degree(expression,rho)},{sp.degree(expression,z)};terms={len(item.terms())}")
    print(f"z_unit={unit_ideal([f,g,z],rho,z)}")
    print(f"response_unit={unit_ideal([f,g,t],rho,z)}")
    print(f"jacobian_unit={unit_ideal([f,g,jac],rho,z)}")
    phi_prefix = sp.factor(phi.subs({d: -omega**2, q: 1, s: omega}))
    psi_prefix = sp.factor(psi.subs({d: -omega**2, q: 1, s: omega}))
    prefix_gcd = sp.factor(sp.gcd(sp.Poly(phi_prefix,omega),sp.Poly(psi_prefix,omega)).as_expr())
    print(f"prefix_gcd_q1={prefix_gcd}")
    print(f"rootzero_response={sp.factor(response.subs({d:0,q:1,s:0}))}")
    print(f"seconds={time.time()-started:.3f}")
