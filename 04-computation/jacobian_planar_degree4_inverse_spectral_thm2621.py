#!/usr/bin/env python3
"""Exact symbolic companion for THM-2621.

The script checks the universal quartic reduction formulas, the rational
degree-four symplectic hostile, the distinction between Jelonek-leading and
primitive-projection divisors, and the trace-Liouville sidecar.  The geometric
generic-coordinate and resultant valuation arguments are proved in the
theorem; finite fixtures here are controls, not extrapolations.
"""

import sympy as sp


checks = 0


def require(condition, message):
    global checks
    if not bool(condition):
        raise RuntimeError(message)
    checks += 1


def is_zero(expression):
    return sp.factor(sp.cancel(expression)) == 0


T, u, v = sp.symbols("T u v")
a0, a1, a2, a3 = sp.symbols("a0 a1 a2 a3")
D = sp.symbols("D0:7")

f_generic = T**4 + a3*T**3 + a2*T**2 + a1*T + a0
raw = sum(D[k] * T**k for k in range(7))
remainder = sp.Poly(raw, T, domain="EX").rem(
    sp.Poly(f_generic, T, domain="EX")
).as_expr()

rho = (
    D[0] - a0*D[4] + a3*a0*D[5] + a0*(a2-a3**2)*D[6],
    D[1] - a1*D[4] + (a3*a1-a0)*D[5]
    + (a0*a3+a1*a2-a1*a3**2)*D[6],
    D[2] - a2*D[4] + (a3*a2-a1)*D[5]
    + (-a0+a1*a3+a2**2-a2*a3**2)*D[6],
    D[3] - a3*D[4] + (a3**2-a2)*D[5]
    + (-a1+2*a2*a3-a3**3)*D[6],
)
for j in range(4):
    require(sp.expand(remainder).coeff(T, j) == sp.expand(rho[j]),
            f"quartic remainder coefficient rho_{j} changed")
require(sp.Poly(remainder, T).degree() <= 3,
        "quartic reduction retained degree at least four")

# Rational C4 hostile on G_m x A^1.
X, Y = sp.symbols("X Y", nonzero=True)
u_xy = Y**4
v_xy = (X-Y)/(4*Y**3)
jacobian = sp.det(sp.Matrix([
    [sp.diff(u_xy, X), sp.diff(u_xy, Y)],
    [sp.diff(v_xy, X), sp.diff(v_xy, Y)],
]))
require(sp.simplify(jacobian) == -1, "rational hostile Jacobian changed")

f = T**4 - 16*u*v*T**2 - 256*u**3*v**4 + 32*u**2*v**2 - u
den_b = 1 - 256*u**2*v**4
b = ((1+48*u*v**2)*T - 4*v*T**3)/den_b

resultant = sp.factor(sp.resultant(Y**4-u, T-Y-4*v*Y**3, Y))
require(sp.expand(resultant-f) == 0, "hostile elimination did not give f")
require(is_zero(f.subs({T: X, u: u_xy, v: v_xy})),
        "hostile primitive coordinate did not satisfy f")
require(is_zero(b.subs({T: X, u: u_xy, v: v_xy})-Y),
        "hostile companion b(X) did not recover Y")

# All partial derivatives below are coefficientwise at fixed T.
pde = sp.together(sp.diff(f, v)*sp.diff(b, u)
                  - sp.diff(f, u)*sp.diff(b, v)
                  + sp.diff(f, T))
pde_numerator = sp.factor(pde.as_numer_denom()[0])
pde_remainder = sp.Poly(pde_numerator, T, domain="EX").rem(
    sp.Poly(f, T, domain="EX")
).as_expr()
require(sp.expand(pde_remainder) == 0,
        "hostile failed f_v b_u-f_u b_v congruent -f_T")

disc = sp.factor(sp.discriminant(f, T))
expected_disc = -256*u**3*(16*u*v**2-1)**2*(16*u*v**2+1)**4
require(sp.expand(disc-expected_disc) == 0,
        "hostile primitive-coordinate discriminant changed")
require(sp.factor(den_b-(1-16*u*v**2)*(1+16*u*v**2)) == 0,
        "b denominator factorization changed")
require(sp.Poly(f, T).LC() == 1,
        "hostile acquired a Jelonek-leading divisor")

# A genuinely D4 abstract inverse-spectral hostile with arbitrary
# trace--Liouville residue.  It satisfies the PDE and the sheet-defect pole
# law, but it comes from rational symplectic coordinates rather than a
# polynomial A2 realization.
s, t, lam = sp.symbols("s t lambda")
phi = s**4 + s**2
phi_prime = sp.diff(phi, s)
u_st = phi
v_st = t/phi_prime
X_st = 1/(s-1)
Y_st = -t*(s-1)**2 + lam*(s-1)

jacobian_XY = sp.det(sp.Matrix([
    [sp.diff(X_st, s), sp.diff(X_st, t)],
    [sp.diff(Y_st, s), sp.diff(Y_st, t)],
]))
jacobian_uv = sp.det(sp.Matrix([
    [sp.diff(u_st, s), sp.diff(u_st, t)],
    [sp.diff(v_st, s), sp.diff(v_st, t)],
]))
require(sp.simplify(jacobian_XY) == 1,
        "D4 hostile source coordinates lost symplectic Jacobian")
require(sp.simplify(jacobian_uv) == 1,
        "D4 hostile target coordinates lost symplectic Jacobian")

R_d4 = (2-u)*T**4 + 6*T**3 + 7*T**2 + 4*T + 1
f_d4 = sp.cancel(R_d4/(2-u))
b_d4_coefficients = (
    2*(-2*lam+10*u*v-3*v),
    -7*lam+26*u*v-26*v,
    2*(-3*lam+11*u*v-16*v),
    -(u-2)*(-lam+4*u*v-6*v),
)
b_d4 = sum(b_d4_coefficients[j]*T**j for j in range(4))

require(is_zero(R_d4.subs({T: X_st, u: u_st})),
        "D4 hostile primitive coordinate did not satisfy R")
require(is_zero(b_d4.subs({T: X_st, u: u_st, v: v_st})-Y_st),
        "D4 hostile companion did not recover Y")

pde_d4 = sp.together(sp.diff(f_d4, v)*sp.diff(b_d4, u)
                     - sp.diff(f_d4, u)*sp.diff(b_d4, v)
                     - sp.diff(f_d4, T))
pde_d4_numerator = sp.factor(pde_d4.as_numer_denom()[0])
pde_d4_remainder = sp.Poly(pde_d4_numerator, T, domain="EX").rem(
    sp.Poly(R_d4, T, domain="EX")
).as_expr()
require(sp.expand(pde_d4_remainder) == 0,
        "D4 hostile failed f_v b_u-f_u b_v congruent f_T")

disc_R_d4 = sp.factor(sp.discriminant(R_d4, T))
disc_f_d4 = sp.factor(sp.discriminant(f_d4, T))
require(sp.expand(disc_R_d4+16*u*(4*u+1)**2) == 0,
        "D4 hostile raw discriminant changed")
require(sp.expand(disc_f_d4
                  +16*u*(4*u+1)**2/(u-2)**6) == 0,
        "D4 hostile normalized discriminant changed")
require(sp.Poly(R_d4.subs(u, 2), T).degree() == 3,
        "D4 hostile lost its one-sheet defect")
require(sp.cancel(sp.Poly(f_d4, T).coeff_monomial(T**3)-6/(2-u)) == 0,
        "D4 hostile first surviving coefficient lost its full pole")
require(all(sp.denom(coefficient) == 1
            for coefficient in b_d4_coefficients),
        "D4 hostile companion acquired a denominator")

# The squared-pair resolvent U(U^2+2U+1+4u) has one rational root and an
# irreducible generic quadratic.  Together with the non-square opposite-root
# radicand this is the D4, rather than V4/C4, lane.
U = sp.symbols("U")
resolvent_d4 = U*(U**2+2*U+1+4*u)
require(sp.factor(sp.discriminant(U**2+2*U+1+4*u, U)) == -16*u,
        "D4 hostile resolvent quadratic changed")
require(sp.gcd(sp.Poly(s**2+1, s), sp.Poly(2*s, s)).degree() == 0,
        "D4 hostile opposite-root radicand lost its simple zeros")

# Reconstruct alpha from the quartic power sums and companion coefficients.
poly_d4 = sp.Poly(f_d4, T)
aa3 = sp.cancel(poly_d4.coeff_monomial(T**3))
aa2 = sp.cancel(poly_d4.coeff_monomial(T**2))
aa1 = sp.cancel(poly_d4.coeff_monomial(T))
aa0 = sp.cancel(poly_d4.coeff_monomial(1))
power_sums_d4 = (
    -aa3,
    aa3**2-2*aa2,
    -aa3**3+3*aa3*aa2-3*aa1,
    aa3**4-4*aa3**2*aa2+2*aa2**2+4*aa3*aa1-4*aa0,
)

def trace_liouville_coefficient(variable):
    return sp.factor(sp.cancel(sum(
        power_sums_d4[j]*sp.diff(b_d4_coefficients[j], variable)
        + sp.Rational(j, j+1)*b_d4_coefficients[j]
          * sp.diff(power_sums_d4[j], variable)
        for j in range(4)
    )))


alpha_d4_u = trace_liouville_coefficient(u)
alpha_d4_v = trace_liouville_coefficient(v)
require(sp.cancel(alpha_d4_u-(-20*v-lam/(2-u))) == 0,
        "D4 hostile alpha du coefficient changed")
require(sp.cancel(alpha_d4_v+4*(4*u+1)) == 0,
        "D4 hostile alpha dv coefficient changed")
omega_d4_u = alpha_d4_u
omega_d4_v = sp.factor(alpha_d4_v-4*u)
require(sp.cancel(omega_d4_v+20*u+4) == 0,
        "D4 hostile exact dv component changed")
require(sp.cancel(sp.diff(omega_d4_v, u)-sp.diff(omega_d4_u, v)) == 0,
        "D4 hostile trace--Liouville form is not closed")
residue_d4 = sp.simplify(
    sp.limit((2-u)*(-omega_d4_u), u, 2)
)
require(residue_d4 == lam,
        "D4 hostile trace--Liouville residue is not lambda")

# The D4 deck involution in the (X,Y) chart.  Simultaneous substitution is
# load-bearing because both images contain X.
tau_X = -X/(2*X+1)
tau_Y = (2*X+1)*(2*lam-(2*X+1)*Y)
tau_squared_X = sp.cancel(tau_X.subs(
    {X: tau_X, Y: tau_Y}, simultaneous=True
))
tau_squared_Y = sp.cancel(tau_Y.subs(
    {X: tau_X, Y: tau_Y}, simultaneous=True
))
require(sp.cancel(tau_squared_X-X) == 0,
        "D4 hostile deck map is not involutive on X")
require(sp.cancel(tau_squared_Y-Y) == 0,
        "D4 hostile deck map is not involutive on Y")
tau_jacobian = sp.factor(sp.det(sp.Matrix([
    [sp.diff(tau_X, X), sp.diff(tau_X, Y)],
    [sp.diff(tau_Y, X), sp.diff(tau_Y, Y)],
])))
require(tau_jacobian == 1,
        "D4 hostile deck map lost symplectic Jacobian")

s_X = 1+1/X
t_XY = lam*X-Y*X**2
u_XY = sp.factor(s_X**4+s_X**2)
v_XY = sp.cancel(t_XY/(4*s_X**3+2*s_X))
tau_u = sp.cancel(u_XY.subs(
    {X: tau_X, Y: tau_Y}, simultaneous=True
))
tau_v = sp.cancel(v_XY.subs(
    {X: tau_X, Y: tau_Y}, simultaneous=True
))
require(sp.cancel(tau_u-u_XY) == 0,
        "D4 hostile deck map does not fix u")
require(sp.cancel(tau_v-v_XY) == 0,
        "D4 hostile deck map does not fix v")
require(sp.cancel(u_XY.subs(X, sp.Rational(-1, 2))-2) == 0,
        "D4 hostile deck pole does not lie over u=2")

# Exact valuation fixtures for every persisting-sheet degree k and pole order
# e: the first coefficient surviving modulo z has normalized pole order e.
z = sp.symbols("z")
valuation_controls = 0
for e in (1, 2, 3):
    for k in range(4):
        fixture = z**e*T**4
        for j in range(4):
            coefficient = (j+2) if j <= k else z*(j+2)
            fixture += coefficient*T**j
        reduced = sp.Poly(fixture.subs(z, 0), T)
        require(reduced.degree() == k,
                "sheet-defect fixture specialized to the wrong degree")
        normalized = sp.cancel(sp.Poly(fixture, T).coeff_monomial(T**k)/z**e)
        numerator, denominator = sp.fraction(normalized)
        require(sp.degree(denominator, z) == e
                and sp.simplify(numerator.subs(z, 0)) != 0,
                "first surviving coefficient lost its full pole")
        valuation_controls += 1

# Trace alpha=Tr(x dy) for Y^4=u.  Multiplication by Y in the basis
# 1,Y,Y^2,Y^3 has trace(Y^-2)=0, hence alpha=4v du.
mult_y = sp.Matrix([
    [0, 0, 0, u],
    [1, 0, 0, 0],
    [0, 1, 0, 0],
    [0, 0, 1, 0],
])
trace_y_minus_two = sp.simplify(sp.trace(mult_y**2)/u)
require(trace_y_minus_two == 0, "trace(Y^-2) changed")
alpha_du = sp.simplify(trace_y_minus_two/4 + 4*v)
require(alpha_du == 4*v, "hostile trace form is not alpha=4v du")
require(sp.diff(alpha_du, v) == 4,
        "hostile trace curvature magnitude changed")

# Affine SL2 source changes alter x dy only by an exact quadratic differential.
A, B, C, E, x, y = sp.symbols("A B C E x y")
xp = A*x+B*y
yp = C*x+E*y
potential = A*C*x**2/2 + B*C*x*y + B*E*y**2/2
dx_coefficient = sp.expand(xp*sp.diff(yp, x)-sp.diff(potential, x))
dy_coefficient = sp.expand(xp*sp.diff(yp, y)-x-sp.diff(potential, y))
require(dx_coefficient == 0, "source Liouville dx coefficient is not exact")
require(sp.expand(dy_coefficient-x*(A*E-B*C-1)) == 0,
        "source Liouville dy coefficient missed determinant one")

print("THM-2621 planar degree-four inverse-spectral exact controls")
print("quartic_remainder_coefficients=4 PASS")
print("rational_C4_hostile: Jacobian=-1 resultant=f companion=b PDE=PASS")
print(f"projection_discriminant={disc}")
print(f"companion_denominator={sp.factor(den_b)} Jelonek_lead=1")
print("rational_D4_hostile: Jacobians=1 resultant_lead=2-u "
      "surviving_sheets=3 companion_denominator=1 PDE=PASS")
print(f"D4_hostile_raw_discriminant={disc_R_d4}")
print("D4_hostile_trace_Liouville=exact+lambda*dlog(2-u) "
      f"residue={residue_d4}")
print("D4_hostile_deck: tau^2=1 Jacobian=1 fixes=(u,v) pole=(2X+1=0)->u=2")
print(f"sheet_defect_valuation_controls={valuation_controls}")
print("trace_form=4*v*du curvature=-4*du^dv for kappa=-1")
print("trace_Liouville_SL2_change=exact_quadratic")
print(f"exact assertions passed: {checks}")
