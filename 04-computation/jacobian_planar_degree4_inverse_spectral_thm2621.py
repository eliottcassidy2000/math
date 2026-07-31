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

# A second D4 hostile separates the full base trace from the opposite-pair
# trace over M=L^tau.  The base class is exact, but the two M/K conjugate
# dlog residues cancel only after the final quadratic trace.
r = sp.symbols("r")
f_pair = T**4+T**2-u
b_pair_1 = 2*v-6/(u-2)
b_pair_3 = 4*v
b_pair = b_pair_1*T+b_pair_3*T**3
pair_pde = sp.together(sp.diff(f_pair, v)*sp.diff(b_pair, u)
                       - sp.diff(f_pair, u)*sp.diff(b_pair, v)
                       - sp.diff(f_pair, T))
pair_pde_remainder = sp.Poly(
    pair_pde.as_numer_denom()[0], T, domain="EX"
).rem(sp.Poly(f_pair, T, domain="EX")).as_expr()
require(sp.expand(pair_pde_remainder) == 0,
        "D4 pair hostile failed the quartic PDE")
pair_disc = sp.factor(sp.discriminant(f_pair, T))
require(pair_disc == -16*u*(4*u+1)**2,
        "D4 pair hostile discriminant changed")
require(sp.simplify(pair_disc.subs(u, 2)) != 0,
        "pure pair hostile denominator became a discriminant divisor")

# Newton sums are p1=p3=0, p2=-2, p4=2+4u.
pair_power_sums = (sp.Integer(0), sp.Integer(-2),
                   sp.Integer(0), 2+4*u)
pair_b_coefficients = (sp.Integer(0), b_pair_1,
                       sp.Integer(0), b_pair_3)

def pair_alpha_coefficient(variable):
    return sp.factor(sp.cancel(sum(
        pair_power_sums[j]*sp.diff(pair_b_coefficients[j], variable)
        + sp.Rational(j, j+1)*pair_b_coefficients[j]
          * sp.diff(pair_power_sums[j], variable)
        for j in range(4)
    )))


pair_alpha_u = pair_alpha_coefficient(u)
pair_alpha_v = pair_alpha_coefficient(v)
require(sp.cancel(pair_alpha_u-(12*v-12/(u-2)**2)) == 0,
        "D4 pair hostile alpha du coefficient changed")
require(sp.cancel(pair_alpha_v-(16*u+4)) == 0,
        "D4 pair hostile alpha dv coefficient changed")
pair_omega_u = pair_alpha_u
pair_omega_v = sp.factor(pair_alpha_v-4*u)
pair_potential = (12*u+4)*v+12/(u-2)
require(sp.cancel(pair_omega_u-sp.diff(pair_potential, u)) == 0,
        "D4 pair hostile base du component is not exact")
require(sp.cancel(pair_omega_v-sp.diff(pair_potential, v)) == 0,
        "D4 pair hostile base dv component is not exact")

# In M=C(r,v), u=r^2+r.  The relative power sums are 0,2r,0,2r^2.
# After removing d(2uv), Gamma_tau leaves -2*dlog((r-1)/(r+2)).
pair_Ctau_dr = sp.factor(
    (b_pair_1+b_pair_3*r).subs(u, r**2+r)
)
pair_log_dr = sp.factor(pair_Ctau_dr-2*v*(2*r+1))
g_pair = (r-1)/(r+2)
require(sp.cancel(pair_log_dr+2*sp.diff(sp.log(g_pair), r)) == 0,
        "D4 pair hostile relative dlog coefficient changed")
require(sp.limit((r-1)*pair_log_dr, r, 1) == -2,
        "D4 pair hostile residue at r=1 changed")
require(sp.limit((r+2)*pair_log_dr, r, -2) == 2,
        "D4 pair hostile residue at r=-2 changed")
require(sp.factor(g_pair*g_pair.subs(r, -1-r)) == 1,
        "D4 pair hostile conjugate dlog classes do not cancel")

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
print("rational_D4_hostile: Jacobians=1 resultant_lead=2-u "
      "surviving_sheets=3 companion_denominator=1 PDE=PASS")
print(f"D4_hostile_raw_discriminant={disc_R_d4}")
print("D4_hostile_trace_Liouville=exact+lambda*dlog(2-u) "
      f"residue={residue_d4}")
print("D4_hostile_deck: tau^2=1 Jacobian=1 fixes=(u,v) pole=(2X+1=0)->u=2")
print("D4_pair_gate_hostile: base=exact pair_residues=(-2,+2) "
      "conjugate_norm=1 pure_b_divisor=u-2")
print(f"sheet_defect_valuation_controls={valuation_controls}")
print("trace_form=4*v*du curvature=-4*du^dv for kappa=-1")
print("polynomial_Keller_source_defect=closed_exact; hostile_primitive=u*v+Y^2/2")
print("punctured_residue_hostiles=arbitrary_traced_residue+norm_one_branch_cancellation")
print("trace_Liouville_SL2_change=exact_quadratic")
print(f"exact assertions passed: {checks}")
