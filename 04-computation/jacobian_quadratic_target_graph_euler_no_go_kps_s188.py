#!/usr/bin/env python3
"""Exact discriminant/root-support audit for THM-3571."""

from __future__ import annotations

import sympy as sp


def require(condition: bool, label: str) -> None:
    if not condition:
        raise RuntimeError(f"failed truth gate: {label}")


def coefficient(expr: sp.Expr, variable: sp.Symbol, degree: int) -> sp.Expr:
    return sp.Poly(sp.expand(expr), variable).coeff_monomial(variable**degree)


def numerator(expr: sp.Expr) -> sp.Expr:
    return sp.primitive(sp.Poly(sp.cancel(expr).as_numer_denom()[0]))[1].as_expr()


def distinct_root_count(expr: sp.Expr, variable: sp.Symbol) -> int:
    poly = sp.Poly(sp.expand(expr), variable, domain=sp.QQ)
    if poly.degree() <= 0:
        return 0
    return int(poly.sqf_part().degree())


a, b, c, s, p, z = sp.symbols("a b c s p z")
A, B, C, D, E = sp.symbols("A B C D E")

# Every total-degree-at-most-two phi whose graph c+phi=0 contains the
# collision value (-1/4,0,0).
phi = A * a**2 + B * a * b + C * b**2 + D * a + E * b + D / 4 - A / 16
require(sp.expand(phi.subs({a: -sp.Rational(1, 4), b: 0})) == 0, "collision compatibility")

# Normalize the Jelonek surface with a=s(2b-s)/12 and take its first sheet
# c=4(-b+2s)/(3s^2).  The graph section has polynomial strict transform G.
a_sheet = s * (2 * b - s) / 12
c_sheet = 4 * (-b + 2 * s) / (3 * s**2)
G_from_sheet = sp.expand(sp.cancel(sp.Rational(3, 4) * s**2 * (c_sheet + phi.subs(a, a_sheet))))

q2 = s**2 * (A * s**2 + 6 * B * s + 36 * C) / 48
q1 = -(A * s**5 + 3 * B * s**4 - 6 * D * s**3 - 36 * E * s**2 + 48) / 48
q0 = s * (A * s**5 - 9 * A * s - 12 * D * s**3 + 36 * D * s + 384) / 192
G = sp.expand(q2 * b**2 + q1 * b + q0)
require(sp.expand(G_from_sheet - G) == 0, "strict-transform coefficients")
require(q2.subs(s, 0) == 0 and q1.subs(s, 0) == -1, "universal vertical-fibre contribution")

# The target Jelonek equation itself splits into the two conjugate sheets.
L = 27 * a**2 * c**2 - 18 * a * b * c + b**3 * c + 16 * a - b**2
sheet_factorization = (
    (4 * b + 3 * c * s**2 - 8 * s)
    * (12 * b**2 * c - 12 * b * c * s - 12 * b + 3 * c * s**2 + 8 * s)
    / 48
)
require(sp.expand(L.subs(a, a_sheet) - sheet_factorization) == 0, "Jelonek sheet factorization")

delta = sp.expand(q1**2 - 4 * q2 * q0)
delta_tilde = sp.expand(
    256
    - 384 * E * s**2
    - 64 * (24 * C + D) * s**3
    + (36 * A * C - 224 * B - 144 * C * D + 144 * E**2) * s**4
    + (6 * A * B - 32 * A - 24 * B * D + 48 * D * E) * s**5
    + (A**2 - 4 * A * D - 24 * B * E + 48 * C * D + 4 * D**2) * s**6
    + 4 * (B * D - 2 * A * E) * s**7
    + (B**2 - 4 * A * C) * s**8
)
require(sp.expand(256 * delta - delta_tilde) == 0, "strict-transform discriminant")
require(delta_tilde.subs(s, 0) == 256, "nonzero discriminant constant")
require(sp.diff(delta_tilde, s).subs(s, 0) == 0, "zero discriminant linear coefficient")

# The omitted target curve is a=b^2/12, c=4/(3b).  Its graph intersection
# is the reduced support of f_tilde; f_tilde(0)=192 keeps b=0 out.
f_tilde = sp.expand(
    A * b**5
    + 12 * B * b**4
    + (144 * C + 12 * D) * b**3
    + 144 * E * b**2
    + (-9 * A + 36 * D) * b
    + 192
)
omitted_graph = sp.cancel(sp.Rational(4, 3) / b + phi.subs(a, b**2 / 12))
require(sp.cancel(144 * b * omitted_graph - f_tilde) == 0, "omitted-curve polynomial")
require(f_tilde.subs(b, 0) == 192 and sp.degree(f_tilde, b) <= 5, "omitted support bound")

# A constant delta_tilde first forces E=0, D=-24C, and
# B=9C(A+96C)/56.  The residual coefficient ideal contains A and C^2.
constant_equations = [coefficient(delta_tilde - 256, s, j) for j in range(2, 9)]
E_constant = sp.Integer(0)
D_constant = -24 * C
B_constant = 9 * C * (A + 96 * C) / 56
constant_residual = [
    numerator(eq.subs({E: E_constant, D: D_constant, B: B_constant}))
    for eq in constant_equations[3:]
]
constant_gb = sp.groebner(constant_residual, A, C, order="lex", method="f5b")
require(constant_gb.reduce(A)[1] == 0, "constant discriminant forces A")
require(constant_gb.reduce(C**2)[1] == 0, "constant discriminant forces C^2")
require(
    all(sp.expand(eq.subs({A: 0, B: 0, C: 0, D: 0, E: 0})) == 0 for eq in constant_equations),
    "constant discriminant affine boundary",
)

# If a nonconstant polynomial with constant 256 and zero linear coefficient
# has exactly two distinct roots, then, after ordering their multiplicities,
# it has the displayed form.  Exhaust all 16 pairs m<=n, m+n<=8, saturating
# by p!=0.  Sequential elimination keeps the exact Groebner computations tiny.
pairs = [(m, n) for m in range(1, 5) for n in range(m, 9 - m)]
require(len(pairs) == 16, "complete multiplicity-pair atlas")
pair_rows: list[tuple[int, int, str, tuple[str, ...]]] = []

for m, n in pairs:
    model = sp.expand(256 * (1 + p * s) ** m * (1 - sp.Rational(m, n) * p * s) ** n)
    equations = [coefficient(delta_tilde - model, s, j) for j in range(2, 9)]

    E_solutions = sp.solve(equations[0], E)
    require(len(E_solutions) == 1, f"unique E elimination ({m},{n})")
    E_sol = E_solutions[0]
    eq3 = sp.expand(equations[1].subs(E, E_sol))
    D_solutions = sp.solve(eq3, D)
    require(len(D_solutions) == 1, f"unique D elimination ({m},{n})")
    D_sol = D_solutions[0]
    eq4 = sp.expand(equations[2].subs({E: E_sol, D: D_sol}))
    B_solutions = sp.solve(eq4, B)
    require(len(B_solutions) == 1, f"unique B elimination ({m},{n})")
    B_sol = B_solutions[0]

    substitutions = {E: E_sol, D: D_sol, B: B_sol}
    residual = [numerator(eq.subs(substitutions)) for eq in equations[3:]]
    gb = sp.groebner(residual + [z * p - 1], z, A, C, p, order="lex", method="f5b")
    basis = tuple(str(poly.as_expr()) for poly in gb.polys)

    if (m, n) == (2, 2):
        require(gb.reduce(A)[1] == 0 and gb.reduce(C)[1] == 0, "(2,2) affine exceptional cell")
        require(sp.expand(B_sol.subs({A: 0, C: 0})) == 0, "(2,2) forces B=0")
        verdict = "AFFINE"
    elif (m, n) == (2, 4):
        require(gb.reduce(C**2)[1] == 0, "(2,4) nilpotent C gate")
        require(gb.reduce(7 * A + 144 * C)[1] == 0, "(2,4) A/C gate")
        require(sp.expand(B_sol.subs({A: 0, C: 0})) == 0, "(2,4) forces B=0 on points")
        verdict = "AFFINE"
    else:
        require(len(gb.polys) == 1 and gb.polys[0].as_expr() == 1, f"unit ideal ({m},{n})")
        verdict = "EMPTY"

    pair_rows.append((m, n, verdict, basis))

# Independent reduced-root controls.  For G=q2*b^2+q1*b+q0,
# nu=rho(q2)-rho(gcd(q2,q1,q0)).
controls = [
    ("pure-C", {A: 0, B: 0, C: 1, D: 0, E: 0}, (3, 1, 3, -2, 4)),
    ("generic-11100", {A: 1, B: 1, C: 1, D: 0, E: 0}, (8, 3, 5, -9, 16)),
    ("generic-12345", {A: 1, B: 2, C: 3, D: 4, E: 5}, (8, 3, 5, -9, 16)),
]
control_rows = []
for name, values, expected in controls:
    q2v = sp.Poly(sp.expand(q2.subs(values)), s, domain=sp.QQ)
    q1v = sp.Poly(sp.expand(q1.subs(values)), s, domain=sp.QQ)
    q0v = sp.Poly(sp.expand(q0.subs(values)), s, domain=sp.QQ)
    common = sp.gcd(sp.gcd(q2v, q1v), q0v)
    rho_delta = distinct_root_count(delta_tilde.subs(values), s)
    nu = distinct_root_count(q2v.as_expr(), s) - distinct_root_count(common.as_expr(), s)
    rho_f = distinct_root_count(f_tilde.subs(values), b)
    chi_D = 2 - rho_delta - nu
    chi_X = -1 + 2 * rho_delta + 2 * nu - rho_f
    row = (rho_delta, nu, rho_f, chi_D, chi_X)
    require(row == expected, f"Euler control {name}")
    control_rows.append((name,) + row)

print("THM-3571 quadratic target-graph Euler audit")
print("G=q2*b^2+q1*b+q0; 256*disc_b(G)=delta_tilde")
print("chi(D_phi)=2-rho(delta_tilde)-nu")
print("chi(X_phi)=-1+2*rho(delta_tilde)+2*nu-rho(f_tilde)")
print("constant-delta Groebner basis:", [str(poly.as_expr()) for poly in constant_gb.polys])
print("two-root multiplicity atlas (m,n,verdict,basis):")
for row in pair_rows:
    print(row)
print("controls (name,rho_delta,nu,rho_f,chi_D,chi_X):")
for row in control_rows:
    print(row)
print("genuine quadratic consequence: rho(delta_tilde)>=3, nu>=1, rho(f_tilde)<=5")
print("therefore chi(X_phi)>=2; no complete quadratic pullback is A2")
print("all active truth gates passed")
