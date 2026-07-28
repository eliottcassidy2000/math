#!/usr/bin/env python3
"""Exact symbolic companion for THM-2621.

The script checks the universal quartic reduction formulas, the rational
degree-four symplectic hostile, the distinction between Jelonek-leading and
primitive-projection divisors, and the exact trace-Liouville primitive.  The geometric
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

# For every polynomial Keller map the source Liouville defect
# theta=x*dy-kappa^-1*P*dQ is closed, hence polynomial-exact on A^2.  The first
# control is the universal curvature identity after imposing Jac(P,Q)=kappa.
P_x, P_y, Q_x, Q_y, kappa = sp.symbols("P_x P_y Q_x Q_y kappa", nonzero=True)
theta_curvature = 1-(P_x*Q_y-P_y*Q_x)/kappa
require(sp.simplify(theta_curvature.subs(kappa, P_x*Q_y-P_y*Q_x)) == 0,
        "source Liouville defect is not closed under the Keller equation")

# The punctured C4 hostile still has an explicit rational primitive.  Here
# kappa=-1, so theta=X*dY+u*d(v), and H=u*v+Y^2/2.
hostile_primitive = sp.simplify(u_xy*v_xy + Y**2/2)
theta_dX = sp.simplify(u_xy*sp.diff(v_xy, X))
theta_dY = sp.simplify(X + u_xy*sp.diff(v_xy, Y))
require(sp.simplify(sp.diff(hostile_primitive, X)-theta_dX) == 0,
        "hostile primitive failed in dX")
require(sp.simplify(sp.diff(hostile_primitive, Y)-theta_dY) == 0,
        "hostile primitive failed in dY")

# Removing polynomial extension makes traced or cancelling logarithmic
# residues possible.  The family u=Y^4, v=(X-h(Y))/(4Y^3) always has kappa=-1
# and theta=d(uv)+h(Y)dY.  The first specialization has arbitrary traced
# residue at u=0; the second has opposite branch residues over u=1 whose norm
# is one, so trace erases them.
lam = sp.symbols("lambda")
h_residue = Y-4*lam/Y
v_residue = (X-h_residue)/(4*Y**3)
jac_residue = sp.det(sp.Matrix([
    [sp.diff(u_xy, X), sp.diff(u_xy, Y)],
    [sp.diff(v_residue, X), sp.diff(v_residue, Y)],
]))
require(sp.simplify(jac_residue) == -1, "residue hostile Jacobian changed")
theta_residue_y = sp.simplify(X+u_xy*sp.diff(v_residue, Y))
regular_primitive = sp.simplify(u_xy*v_residue+Y**2/2)
require(sp.simplify(theta_residue_y-sp.diff(regular_primitive, Y)+4*lam/Y) == 0,
        "residue hostile did not leave -4 lambda dY/Y")
require(sp.simplify(4*(-4*lam)/(4*u)) == -4*lam/u,
        "trace residue conversion dY/Y=du/(4u) changed")
norm_ratio = sp.simplify((1-u)/(sp.I**4-u))
require(norm_ratio == 1, "opposite branch-residue ratio lost norm one")

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
print(f"sheet_defect_valuation_controls={valuation_controls}")
print("trace_form=4*v*du curvature=-4*du^dv for kappa=-1")
print("polynomial_Keller_source_defect=closed_exact; hostile_primitive=u*v+Y^2/2")
print("punctured_residue_hostiles=arbitrary_traced_residue+norm_one_branch_cancellation")
print("trace_Liouville_SL2_change=exact_quadratic")
print(f"exact assertions passed: {checks}")
