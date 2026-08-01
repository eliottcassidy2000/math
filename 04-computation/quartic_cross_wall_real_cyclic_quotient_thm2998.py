#!/usr/bin/env python3
"""Exact companion for THM-2998.

The companion verifies the derivative-square reciprocal cross wall of a
depressed quartic, its affine (but not projective) covariance, exact real
sign chambers, the smooth ordered star/triangle plane cubic, and the
degree-three cyclic quotient to the standard ``X_0(14)`` elliptic model.
Every truth-bearing gate is an explicit ``require`` call, so ``python -O``
remains a genuine verification run.
"""

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def quartic_invariants(p_value, q_value, r_value):
    """Return Delta, K, H for T^4+pT^2+qT+r."""
    substitution = {
        a: 2 * p_value,
        b: p_value**2 - 4 * r_value,
        c: -(q_value**2),
    }
    quartic = T**4 + p_value * T**2 + q_value * T + r_value
    return (
        sp.factor(sp.discriminant(quartic, T)),
        sp.factor(K.subs(substitution)),
        sp.factor(H.subs(substitution)),
    )


def weierstrass_invariants(a1, a2, a3, a4, a6):
    """Return c4,c6,Delta for a generalized Weierstrass equation."""
    b2 = a1**2 + 4 * a2
    b4 = 2 * a4 + a1 * a3
    b6 = a3**2 + 4 * a6
    b8 = a1**2 * a6 + 4 * a2 * a6 - a1 * a3 * a4 + a2 * a3**2 - a4**2
    c4 = b2**2 - 24 * b4
    c6 = -(b2**3) + 36 * b2 * b4 - 216 * b6
    discriminant = -(b2**2) * b8 - 8 * b4**3 - 27 * b6**2 + 9 * b2 * b4 * b6
    return sp.factor(c4), sp.factor(c6), sp.factor(discriminant)


# ---------------------------------------------------------------------------
# The coefficient wall and derivative-square norm.
# ---------------------------------------------------------------------------

T, U, Z = sp.symbols("T U Z")
p, q, r = sp.symbols("p q r")
a, b, c = sp.symbols("a b c")

S = U**3 + a * U**2 + b * U + c
D = sp.expand(sp.discriminant(S, U))
K = sp.expand(b**3 - a**3 * c)
H = sp.expand(
    a**6 * c**2
    - 2 * a**4 * b * c**2
    - 2 * a**3 * b**3 * c
    - 26 * a**3 * c**3
    + 29 * a**2 * b**2 * c**2
    - 2 * a * b**4 * c
    - 18 * a * b * c**3
    + b**6
    - 26 * b**3 * c**2
    + 189 * c**4
)

f = T**4 + p * T**2 + q * T + r
D_quartic = sp.expand(sp.discriminant(f, T))
abc_from_pqr = {a: 2 * p, b: p**2 - 4 * r, c: -(q**2)}
require(sp.expand(D.subs(abc_from_pqr) - D_quartic) == 0, "matching/quartic discriminants differ")

norm_derivative = sp.factor(sp.resultant(f, sp.diff(f, T) ** 4 - D_quartic, T))
require(
    sp.expand(norm_derivative - 2**8 * D_quartic**2 * H.subs(abc_from_pqr)) == 0,
    "derivative-fourth-power norm wall changed",
)
require(sp.factor(H) == H, "H is no longer irreducible over Q")
require(sp.Poly(H, a, b, c).primitive()[0] == 1, "H is no longer primitive")

# Weighted affine covariance: matching-cubic weights are 2,4,6 and H has
# weight 24.  Translation disappears on depressing the quartic.
lam = sp.symbols("lam", nonzero=True)
require(
    sp.expand(H.subs({a: lam**2 * a, b: lam**4 * b, c: lam**6 * c}) - lam**24 * H) == 0,
    "affine weight-24 covariance changed",
)

# A simple smooth point prevents an unnoticed proper power on the wall.
wall_abc = {a: -14, b: 49, c: -49}
require(H.subs(wall_abc) == 0, "wall control left H=0")
wall_gradient = tuple(sp.diff(H, variable).subs(wall_abc) for variable in (a, b, c))
require(any(value != 0 for value in wall_gradient), "wall control is no longer smooth")

# The wall is not PGL2-invariant.  Start with f=T^4-7T^2+7T on H=0 and
# apply U=(T+14)/(1-2T), whose determinant is 29 and whose pole is not a
# root.  The transformed monic quartic is already depressed and integral.
mobius_U = sp.symbols("mobius_U")
f_wall = T**4 - 7 * T**2 + 7 * T
require(f_wall.subs(T, sp.Rational(1, 2)) != 0, "PGL2 hostile pole became a root")
mobius_inverse = (mobius_U - 14) / (2 * mobius_U + 1)
mobius_raw = sp.cancel((2 * mobius_U + 1) ** 4 * f_wall.subs(T, mobius_inverse))
mobius_monic = sp.expand(mobius_raw / sp.Poly(mobius_raw, mobius_U).LC())
require(
    mobius_monic == mobius_U**4 - 161 * mobius_U**2 - 581 * mobius_U + 1274,
    "integral PGL2 hostile quartic changed",
)
p_inv, q_inv, r_inv = -161, -581, 1274
D_inv, K_inv, H_inv = quartic_invariants(p_inv, q_inv, r_inv)
require(D_inv == 1428170793721, "PGL2 hostile discriminant changed")
require(K_inv == -2238496245503, "PGL2 hostile K changed")
require(H_inv == -641584271212871414216880, "PGL2 hostile did not leave H=0")


# ---------------------------------------------------------------------------
# Exact real-root chambers.
# ---------------------------------------------------------------------------

# Zero-real proof in normalized variables.  Roots are u+/-ib and +/-id
# after a real translation.  If f'(u+ib)^4 is positive real then f' is real
# or purely imaginary.  The two displayed factors show equality with Delta
# forces b^2=d^2, which is inseparable in either case.
u = sp.symbols("u", real=True)
beta, delta = sp.symbols("beta delta", real=True, positive=True)
A_real = u**2 - beta**2 + delta**2
A_imag = 2 * beta * u
derivative_complex = sp.expand(2 * sp.I * beta * (A_real + sp.I * A_imag))
D_zero_real = sp.expand(16 * beta**2 * delta**2 * (A_real**2 + A_imag**2) ** 2)
require(
    sp.factor((derivative_complex**4 - D_zero_real).subs(u, 0))
    == 16 * beta**2 * (beta - delta) ** 5 * (beta + delta) ** 5,
    "same-centre zero-real obstruction changed",
)
real_axis_difference = sp.factor(
    sp.rem(
        sp.Poly(sp.expand(derivative_complex**4 - D_zero_real), u),
        sp.Poly(u**2 - beta**2 + delta**2, u),
    ).as_expr()
)
require(
    real_axis_difference == 256 * beta**6 * (beta - delta) ** 3 * (beta + delta) ** 3,
    "real-derivative zero-real obstruction changed",
)

# Two-real equality classification.  Normalize the conjugate pair to +/-i.
# Let the real roots be y-h and y+h.  Magnitudes force h=1; phase then gives
# y^2-2=+/-2y and exactly four oriented affine representatives.
y, h = sp.symbols("y h", real=True)
A_two = sp.expand((sp.I - (y - h)) * (sp.I - (y + h)))
A_two_bar = sp.conjugate(A_two)
fprime_fourth_two = sp.expand(16 * A_two**4)
D_two_general = sp.expand(-16 * h**2 * A_two**2 * A_two_bar**2)
require(
    sp.factor(
        fprime_fourth_two * sp.conjugate(fprime_fourth_two)
        - D_two_general**2
        - 2**8 * (1 - h**4) * (A_two * A_two_bar) ** 4
    )
    == 0,
    "two-real magnitude equation no longer forces h=1",
)
require(sp.expand(A_two.subs(h, 1) - (y**2 - 2 - 2 * sp.I * y)) == 0, "two-real normalized product changed")
four_y = (
    -1 - sp.sqrt(3),
    -1 + sp.sqrt(3),
    1 - sp.sqrt(3),
    1 + sp.sqrt(3),
)
phase_polynomial = sp.expand((y**2 - 2) ** 2 - 4 * y**2)
require(
    sp.factor(phase_polynomial) == (y**2 - 2 * y - 2) * (y**2 + 2 * y - 2),
    "two-real phase equation changed",
)
require(all(sp.simplify(phase_polynomial.subs(y, value)) == 0 for value in four_y), "two-real class list changed")

# Exact algebraic representative and its reflected partners.
sqrt3 = sp.sqrt(3)
p_two = -2 - sqrt3
q_two = -2 * (1 + sqrt3)
r_two = (3 + 4 * sqrt3) / 4
D_two, K_two, H_two = quartic_invariants(p_two, q_two, r_two)
require(sp.simplify(D_two + 4096 * (7 + 4 * sqrt3)) == 0, "two-real Delta changed")
require(sp.simplify(K_two + 512 * (12 + 7 * sqrt3)) == 0, "two-real K changed")
require(H_two == 0, "two-real equality representative left H=0")

# Four-real gaps.  rho_i=f'(r_i)^2/sqrt(Delta), product rho_i=1.
alpha, beta_gap, gamma = sp.symbols("alpha beta_gap gamma", positive=True)
total = alpha + beta_gap + gamma
rho = (
    alpha * (alpha + beta_gap) * total / (beta_gap * gamma * (beta_gap + gamma)),
    alpha * beta_gap * (beta_gap + gamma) / (gamma * (alpha + beta_gap) * total),
    beta_gap * gamma * (alpha + beta_gap) / (alpha * (beta_gap + gamma) * total),
    gamma * (beta_gap + gamma) * total / (alpha * beta_gap * (alpha + beta_gap)),
)
require(sp.cancel(sp.prod(rho) - 1) == 0, "four-real rho product changed")
gap_equations = tuple(sp.factor(sp.together(value - 1).as_numer_denom()[0]) for value in rho)
require(
    gap_equations
    == tuple(
        sp.expand(expression)
        for expression in (
            alpha * (alpha + beta_gap) * total - beta_gap * gamma * (beta_gap + gamma),
            alpha * beta_gap * (beta_gap + gamma) - gamma * (alpha + beta_gap) * total,
            beta_gap * gamma * (alpha + beta_gap) - alpha * (beta_gap + gamma) * total,
            gamma * (beta_gap + gamma) * total - alpha * beta_gap * (alpha + beta_gap),
        )
    ),
    "four-real collision equations changed",
)

# Exact four-real wall, negative chamber, and positive chamber controls.
require(quartic_invariants(-7, 7, 0) == (2401, -16807, 0), "four-real wall hostile changed")
require(quartic_invariants(-4, 2, 1) == (592, -320, -154368), "four-real negative control changed")
require(quartic_invariants(-10, 1, 0) == (3973, 992000, 980121756189), "four-real positive control changed")
for quartic_control in (T**4 - 7 * T**2 + 7 * T, T**4 - 4 * T**2 + 2 * T + 1, T**4 - 10 * T**2 + T):
    require(sp.Poly(quartic_control, T).count_roots(-sp.oo, sp.oo) == 4, "four-real control changed chamber")


# ---------------------------------------------------------------------------
# Ordered plane cubic, free C3 quotient, and the X_0(14) model.
# ---------------------------------------------------------------------------

x, y_root, z = sp.symbols("x y_root z")
vandermonde = (x - y_root) * (x - z) * (y_root - z)
C_plus = sp.expand(x * y_root * z - vandermonde)
C_minus = sp.expand(x * y_root * z + vandermonde)

# An odd transposition exchanges the two sign sheets.
require(
    sp.expand(C_plus.subs({y_root: z, z: y_root}, simultaneous=True) - C_minus) == 0,
    "odd permutation no longer exchanges C+ and C-",
)

# Smoothness of C+ in all three projective charts; C- follows by the swap.
partials = [sp.diff(C_plus, variable) for variable in (x, y_root, z)]
for chart, variables in (({x: 1}, (y_root, z)), ({y_root: 1}, (x, z)), ({z: 1}, (x, y_root))):
    chart_polynomials = [sp.expand(polynomial.subs(chart)) for polynomial in (C_plus, *partials)]
    groebner = sp.groebner(chart_polynomials, *variables, order="lex")
    require(groebner.polys == [sp.Poly(1, *variables)], "ordered cubic is singular in a projective chart")

cyclic = {x: y_root, y_root: z, z: x}
require(sp.expand(C_plus.subs(cyclic, simultaneous=True) - C_plus) == 0, "C3 no longer preserves C+")

# A projective C3 fixed point is [1:l:l^2] with l^3=1.  None lies on C+.
ell = sp.symbols("ell")
fixed_value = sp.expand(C_plus.subs({x: 1, y_root: ell, z: ell**2}))
require(sp.gcd(fixed_value, ell**3 - 1) == 1, "C3 action acquired a fixed point")

e1 = x + y_root + z
e2 = x * y_root + y_root * z + z * x
e3 = x * y_root * z
for symmetric in (e1, e2, e3):
    require(sp.expand(symmetric.subs(cyclic, simultaneous=True) - symmetric) == 0, "quotient coordinate lost C3 invariance")

# The orbit of [1:0:0] has exactly three points with one symmetric image.
orbit = ((1, 0, 0), (0, 0, 1), (0, 1, 0))
require(len(set(orbit)) == 3, "C3 orbit control collapsed")
require(all(C_plus.subs({x: point[0], y_root: point[1], z: point[2]}) == 0 for point in orbit), "C3 orbit left C+")
require(
    len({tuple(symmetric.subs({x: point[0], y_root: point[1], z: point[2]}) for symmetric in (e1, e2, e3)) for point in orbit}) == 1,
    "symmetric quotient does not collapse one C3 orbit",
)

# Projection of the ordered cubic from [1:0:0] gives a genus-one quartic.
tau = sp.symbols("tau")
projected = sp.factor(C_plus.subs({x: 1, z: tau * y_root}) / y_root)
qa = tau**2 - tau
qb = -(tau**2) + tau + 1
qc = tau - 1
require(sp.expand(projected - (qa * y_root**2 + qb * y_root + qc)) == 0, "ordered projection quadratic changed")
ordered_quartic = sp.factor(qb**2 - 4 * qa * qc)
require(ordered_quartic == tau**4 - 6 * tau**3 + 7 * tau**2 - 2 * tau + 1, "ordered genus-one quartic changed")
I_ordered = 12 * 1 * 1 - 3 * (-6) * (-2) + 7**2
J_ordered = 72 * 1 * 7 * 1 + 9 * (-6) * 7 * (-2) - 27 * 1 * (-2) ** 2 - 27 * (-6) ** 2 * 1 - 2 * 7**3
disc_ordered = sp.discriminant(ordered_quartic, tau)
j_ordered = sp.Rational(256 * I_ordered**3, disc_ordered)
require((I_ordered, J_ordered, disc_ordered, j_ordered) == (25, -506, -7168, -sp.Rational(15625, 28)), "ordered invariants changed")

# On e1!=0, t=e2/e1^2 and s=e3/e1^3.  Squaring C+ gives the exact
# symmetric quotient equation.  Its projective completion is the C3 quotient.
t, s, v = sp.symbols("t s v")
quotient_wall = sp.expand(t**2 - 4 * t**3 - 4 * s + 18 * t * s - 28 * s**2)
require(
    sp.expand(
        sp.discriminant(U**3 - U**2 + t * U - s, U) - s**2 - quotient_wall
    )
    == 0,
    "symmetric discriminant quotient changed",
)
v_definition = 28 * s + 2 - 9 * t
require(
    sp.expand(v_definition**2 - (2 - 7 * t) * (16 * t**2 - 11 * t + 2) + 28 * quotient_wall) == 0,
    "quotient double-cover equation changed",
)

x_e, y_e = sp.symbols("x_e y_e")
x_definition = 32 - 112 * t
y_definition = -112 * v_definition
require(
    sp.expand(y_definition**2 - x_definition * (x_definition**2 + 13 * x_definition + 128) + 351232 * quotient_wall) == 0,
    "quotient elliptic model changed",
)
t_inverse = (32 - x_e) / 112
s_inverse = (64 - 9 * x_e - y_e) / 3136
require(sp.expand(t_inverse.subs(x_e, x_definition) - t) == 0, "quotient inverse t changed")
require(
    sp.expand(s_inverse.subs({x_e: x_definition, y_e: y_definition}) - s) == 0,
    "quotient inverse s changed",
)

# Exact isomorphism to the standard X_0(14) model.
X, Y0 = sp.symbols("X Y0")
x_change = 4 * X - 4
y_change = 8 * Y0 + 4 * X + 4
X_inverse = (x_e + 4) / 4
Y_inverse = (y_e - x_e - 8) / 8
require(sp.expand(X_inverse.subs(x_e, x_change) - X) == 0, "X0 inverse X changed")
require(
    sp.expand(Y_inverse.subs({x_e: x_change, y_e: y_change}) - Y0) == 0,
    "X0 inverse Y changed",
)
E0_equation = sp.expand(y_change**2 - x_change * (x_change**2 + 13 * x_change + 128))
X014_equation = Y0**2 + X * Y0 + Y0 - X**3 - 4 * X + 6
require(sp.expand(E0_equation - 64 * X014_equation) == 0, "X0(14) change of variables changed")

E0_invariants = weierstrass_invariants(0, 13, 0, 128, 0)
X014_invariants = weierstrass_invariants(1, 0, 1, 4, -6)
require(E0_invariants == (-3440, 338624, -(2**18) * 7**3), "quotient invariants changed")
require(X014_invariants == (-215, 5291, -(2**6) * 7**3), "X0(14) invariants changed")
require(
    E0_invariants == (2**4 * X014_invariants[0], 2**6 * X014_invariants[1], 2**12 * X014_invariants[2]),
    "Weierstrass scaling under the X0(14) isomorphism changed",
)
j_quotient = sp.factor(X014_invariants[0] ** 3 / X014_invariants[2])
require(j_quotient == sp.Rational(9938375, 21952), "quotient j-invariant changed")
require(j_ordered != j_quotient, "ordered cubic became birational to its quotient")


print("THM2998 exact quartic cross-wall real/cyclic quotient companion")
print("norm=2^8*Delta^2*H irreducible=1 affine_weight=24")
print("PGL2_hostile_H=-641584271212871414216880")
print("real_chambers=zero:H>0 two:H>=0 four:gap-sign")
print("two_real_oriented_classes=-1-sqrt3,-1+sqrt3,1-sqrt3,1+sqrt3")
print("four_real_gap_walls=4 sign_controls=negative/positive")
print("ordered_plane_cubics=Cplus/Cminus smooth=1 genus=1")
print("C3_action=free quotient_degree=3")
print("ordered_j=-15625/28")
print("quotient_E0=y^2-x*(x^2+13*x+128)")
print("X0_14=Y^2+X*Y+Y-X^3-4*X+6")
print("X0_change=x:4X-4,y:8Y+4X+4")
print("quotient_j=9938375/21952 ordered_not_birational=1")
print("wall_hostile_D_K_H=2401,-16807,0")
print("all_exact_checks=PASS")
