#!/usr/bin/env python3
"""Independent hostile exact audit for THM-3850.

This checker does not import the primary companion.  It freezes a second
derivation of the quadratic branch, the conjugate-numerator valuation law,
finite-root and repeated-discriminant controls, and a rational parameter on
the minimal reducible residual conic.
"""

from __future__ import annotations

import hashlib
import sys

import sympy as sp


if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(newline="\n")


GATES = 0


def gate(condition: bool, label: str) -> None:
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(f"gate failed: {label}")


def eq(left: sp.Expr, right: sp.Expr, label: str) -> None:
    gate(sp.factor(left - right) == 0, label)


def valuation_at_zero(expr: sp.Expr, variable: sp.Symbol) -> int:
    expr = sp.expand(expr)
    poly = sp.Poly(expr, variable)
    terms = poly.terms()
    gate(bool(terms), "valuation input is nonzero")
    return min(monomial[0] for monomial, coefficient in terms if coefficient != 0)


def squarefree_part(poly: sp.Expr, variable: sp.Symbol) -> sp.Expr:
    return sp.factor(sp.Poly(poly, variable).sqf_part().as_expr())


# A second derivation of the branch polynomial from the original packet.
a, x, q, z = sp.symbols("a x q z")
p = sp.Rational(3, 2) + a * x
u = 1 + a * x + a**2 * q
branch = sp.factor((8 * p**3 - 27 * u**2) / a**2)
expected = -27 * a**2 * q**2 + (8 * x**3 - 54 * x * q) * a + 9 * x**2 - 54 * q
eq(branch, expected, "branch recovered from cubic packet")

poly_a = sp.Poly(branch, a)
qa = poly_a.coeff_monomial(a**2)
qb = poly_a.coeff_monomial(a)
qc = poly_a.coeff_monomial(1)
ell_raw = 9 * q - 2 * x**2
eq(qa, -27 * q**2, "quadratic leading coefficient")
eq(qb, 8 * x**3 - 54 * x * q, "quadratic middle coefficient")
eq(qc, 9 * x**2 - 54 * q, "quadratic constant coefficient")
eq(sp.discriminant(poly_a.as_expr(), a), -8 * ell_raw**3, "quadratic discriminant")

# Put z^2=-8L.  The two quadratic numerators are B-Lz and B+Lz.
# Their product is exactly 4*a_2*a_0.  At a root of q this identity
# forces the cancelling numerator to have valuation exactly 2 ord(q),
# because its conjugate and x^2-6q are units.  This closes the possible
# high-multiplicity cancellation loophole directly, without interpolation.
nminus = qb - ell_raw * z
nplus = qb + ell_raw * z
conjugate_product = sp.expand(nminus * nplus).subs(z**2, -8 * ell_raw)
eq(conjugate_product, 4 * qa * qc, "conjugate numerator product")
eq(conjugate_product, -972 * q**2 * (x**2 - 6 * q), "valuation factorization")

c = sp.symbols("c", nonzero=True)
eq(ell_raw.subs({x: c, q: 0}), -2 * c**2, "profile-root L value")
eq(nminus.subs({x: c, q: 0, z: -4 * c}), 0, "cancelling sheet numerator")
eq(nplus.subs({x: c, q: 0, z: -4 * c}), 16 * c**3, "unit conjugate numerator")
eq((-18 * (x**2 - 6 * q) / nplus).subs({x: c, q: 0, z: -4 * c}), -sp.Rational(9, 8) / c,
   "rationalized finite value")
eq(branch.subs({x: c, q: 0}), c**2 * (8 * a * c + 9), "finite specialization")
eq(sp.diff(branch, a).subs({x: c, q: 0, a: -sp.Rational(9, 8) / c}), 8 * c**3,
   "finite point smoothness")

# Explicit local series for multiplicities one through eight independently
# verifies exact 2m cancellation and the same finite normalized value.
t = sp.symbols("t")
local_rows: list[tuple[int, int, sp.Expr]] = []
for multiplicity in range(1, 9):
    local_q = t**multiplicity * (2 + t)
    local_x = 1 + t
    local_l = sp.expand(9 * local_q - 2 * local_x**2)
    local_b = sp.expand(8 * local_x**3 - 54 * local_x * local_q)
    order = 2 * multiplicity + 3
    # z=-4 sqrt(L/(-2)) is the sheet with z(0)=-4 and z^2=-8L.
    root_series = sp.series(sp.sqrt(local_l / -2), t, 0, order).removeO()
    cancelling = sp.series(local_b + 4 * local_l * root_series, t, 0, order).removeO()
    cancel_order = valuation_at_zero(cancelling, t)
    gate(cancel_order == 2 * multiplicity, f"multiplicity {multiplicity}: exact cancellation order")
    regular_value = sp.series(cancelling / (54 * local_q**2), t, 0, 1).removeO()
    eq(regular_value, -sp.Rational(9, 8), f"multiplicity {multiplicity}: finite A value")
    local_rows.append((multiplicity, cancel_order, regular_value))

# Finite exact controls with several distinct roots and multiplicities.
# For each one, factorization over QQ checks irreducibility of the branch,
# while the root packet checks one present and one deleted normalization point.
profiles = [
    ("simple", x + 1, [-1]),
    ("double", (x - 1) ** 2, [1]),
    ("mixed", (x - 1) ** 3 * (x + 2) ** 2, [1, -2]),
    ("three_roots", (x - 1) * (x - 2) * (x + 3), [1, 2, -3]),
    ("uneven", (x - 1) ** 4 * (x + 2) * (x - 3) ** 2, [1, -2, 3]),
]
profile_rows: list[tuple[str, int, int, int, int]] = []
for name, profile, roots in profiles:
    gate(sp.expand(profile).subs(x, 0) != 0, f"{name}: origin is not a profile root")
    specialized = sp.Poly(branch.subs(q, profile), a, x)
    factors = sp.factor_list(specialized.as_expr())[1]
    gate(len(factors) == 1 and factors[0][1] == 1, f"{name}: branch irreducible over QQ")
    L_profile = sp.factor(ell_raw.subs(q, profile))
    ell_profile = squarefree_part(L_profile, x)
    r = int(sp.degree(ell_profile, x))
    gate(r >= 1, f"{name}: nonconstant squarefree part")
    radical = sp.factor(sp.Poly(profile, x).sqf_part().as_expr())
    s = int(sp.degree(radical, x))
    gate(s == len(roots), f"{name}: radical degree")
    for root in roots:
        root = sp.Integer(root)
        eq(L_profile.subs(x, root), -2 * root**2, f"{name}: L at root {root}")
        eq(qb.subs({q: 0, x: root}), 8 * root**3, f"{name}: B at root {root}")
        eq(nminus.subs({q: 0, x: root, z: -4 * root}), 0,
           f"{name}: cancelling sheet at root {root}")
        eq(nplus.subs({q: 0, x: root, z: -4 * root}), 16 * root**3,
           f"{name}: pole sheet at root {root}")
    nu_infinity = 1 if r % 2 else 2
    genus = (r - 1) // 2
    total = s + nu_infinity
    gate(total >= 2, f"{name}: at least two deleted points")
    profile_rows.append((name, s, r, genus, total))

# A deliberately repeated discriminant divisor: L=(x-1)^2(x+1), but q(0)
# is nonzero and L is nonsquare.  The two points over x=1 and the ramified
# point over x=-1 are retained because q is a unit there.  This guards
# against mistaking singular resolution points for punctures.
h = x - 1
ell = x + 1
repeated_L = sp.expand(h**2 * ell)
repeated_q = sp.factor((repeated_L + 2 * x**2) / 9)
eq(9 * repeated_q - 2 * x**2, repeated_L, "repeated-L profile construction")
gate(repeated_q.subs(x, 0) != 0, "repeated-L origin unit")
gate(repeated_q.subs(x, 1) != 0, "even discriminant root retained")
gate(repeated_q.subs(x, -1) != 0, "odd discriminant root retained")
repeated_factors = sp.factor_list(branch.subs(q, repeated_q))[1]
gate(len(repeated_factors) == 1 and repeated_factors[0][1] == 1,
     "repeated-L branch remains irreducible")

# Square-L hostile: q(0) is a unit but L is a square, so the branch really
# splits into two distinct linear factors.  The nonsquare hypothesis is not
# cosmetic, and irreducibility excludes this packet.
square_h = x + 1
square_q = sp.factor((square_h**2 + 2 * x**2) / 9)
square_branch = sp.factor(branch.subs(q, square_q))
# The base field is algebraically closed; adjoin sqrt(-2) so the exact
# factorization is performed over a field where -8 is a square.
square_factors = sp.factor_list(square_branch, extension=sp.sqrt(-2))[1]
eq(9 * square_q - 2 * x**2, square_h**2, "square-L hostile construction")
gate(square_q.subs(x, 0) != 0, "square-L hostile has origin unit")
gate(len(square_factors) == 2 and all(exponent == 1 for factor, exponent in square_factors),
     "square-L hostile splits into two reduced factors")

# General primitive-reducible packet from the incoming all-profile extension.
# If L=h^2, write u=h-sigma*x and v=h+sigma*x with sigma^2=-2.
# Reducing independently modulo sigma^2+2 gives the two denominator graphs.
sigma, hh = sp.symbols("sigma hh")
u_factor = hh - sigma * x
v_factor = hh + sigma * x
square_profile = (hh**2 + 2 * x**2) / 9
E_v = v_factor**2 * a - 3 * (x - sigma * hh)
E_u = u_factor**2 * a - 3 * (x + sigma * hh)
factor_identity = sp.together(branch.subs(q, square_profile) + sp.Rational(1, 3) * E_u * E_v)
factor_numerator = sp.expand(factor_identity.as_numer_denom()[0])
factor_remainder = sp.rem(factor_numerator, sigma**2 + 2, sigma)
eq(factor_remainder, 0, "general primitive reducible factorization")
product_remainder = sp.rem(sp.expand(u_factor * v_factor - 9 * square_profile), sigma**2 + 2, sigma)
eq(product_remainder, 0, "primitive factor product is 9q")

# At a zero of either denominator the other term is a unit.  Since h(0) is
# nonzero in the primitive case, such a root is nonzero.  This gives the
# exact Laurent rings k[x,u^-1] and k[x,v^-1].
root_v_value = sp.rem((x - sigma * hh).subs({x: c, hh: -sigma * c}), sigma**2 + 2, sigma)
root_u_value = sp.rem((x + sigma * hh).subs({x: c, hh: sigma * c}), sigma**2 + 2, sigma)
eq(root_v_value, -c, "v-root numerator unit")
eq(root_u_value, -c, "u-root numerator unit")

# Positive A1 control: u is a nonzero constant and v is linear.  Its mate is
# exactly a two-place Laurent line.  This is the general source of the
# degree-one A1+G_m specialization checked below.
alpha = sp.symbols("alpha", nonzero=True)
h_linear_sigma = sigma * x + alpha
u_linear_sigma = sp.rem(u_factor.subs(hh, h_linear_sigma), sigma**2 + 2, sigma)
v_linear_sigma = sp.rem(v_factor.subs(hh, h_linear_sigma), sigma**2 + 2, sigma)
eq(u_linear_sigma, alpha, "primitive positive control constant factor")
eq(v_linear_sigma, 2 * sigma * x + alpha, "primitive positive control linear factor")
E_u_linear_sigma = sp.rem(sp.expand(E_u.subs(hh, h_linear_sigma)), sigma**2 + 2, sigma)
eq(E_u_linear_sigma, alpha**2 * a + 3 * x - 3 * sigma * alpha,
   "primitive positive A1 graph")
mate_root = alpha * sigma / 4
eq(sp.rem(sp.expand(v_linear_sigma.subs(x, mate_root)), sigma**2 + 2, sigma), 0,
   "primitive mate finite denominator root")
mate_numerator = sp.rem((x - sigma * h_linear_sigma).subs(x, mate_root), sigma**2 + 2, sigma)
gate(mate_numerator != 0, "primitive mate really omits its denominator root")

# The phrase "the first reducible boundary" needs care: there is a second
# degree-one reducible locus.  In fact all degree-one profiles can be closed.
# For q=k*x+beta with k != 0 and beta != 0, the branch is primitive and
# reducible exactly when beta=-9*k^2/8.  At that value it is an A1 graph plus
# a G_m graph, not the disjoint A1 plus three-puncture conic of q=k*x.
k, beta = sp.symbols("k beta", nonzero=True)
linear_q = k * x + beta
linear_L = sp.factor(9 * linear_q - 2 * x**2)
eq(sp.discriminant(linear_L, x), 9 * (9 * k**2 + 8 * beta),
   "degree-one square criterion")
special_linear_q = k * x - sp.Rational(9, 8) * k**2
special_linear_branch = sp.factor(branch.subs(q, special_linear_q))
linear_graph = a * k**2 - sp.Rational(8, 27) * x + sp.Rational(4, 3) * k
linear_d = sp.Rational(8, 9) * x - k
linear_residual = a * linear_d**2 + sp.Rational(8, 9) * x - sp.Rational(4, 3) * k
eq(special_linear_branch, -sp.Rational(2187, 64) * linear_graph * linear_residual,
   "second degree-one reducible factorization")
eq(9 * special_linear_q - 2 * x**2, -sp.Rational(1, 8) * (4 * x - 9 * k) ** 2,
   "second degree-one square L")
eq(linear_residual, a * linear_d**2 + linear_d - k / 3,
   "linear residual Laurent presentation")
# In the residual coordinate ring, d*(3/k)*(1+a*d)=1; hence its ring is
# exactly k[d,d^{-1}], proving that component is G_m rather than merely
# birational to it.
eq(linear_d * (sp.Rational(3, 1) / k) * (1 + a * linear_d) - 1,
   (sp.Rational(3, 1) / k) * linear_residual,
   "linear residual inverse certificate")
linear_hole = sp.Rational(9, 8) * k
eq((sp.Rational(8, 9) * x - k).subs(x, linear_hole), 0,
   "linear residual denominator root")
eq((sp.Rational(8, 9) * x - sp.Rational(4, 3) * k).subs(x, linear_hole),
   -sp.Rational(1, 3) * k, "linear residual genuinely omits its finite root")
linear_meeting = sp.Rational(9, 4) * k
linear_meeting_a = -sp.Rational(2, 3) / k
eq(linear_graph.subs({x: linear_meeting, a: linear_meeting_a}), 0,
   "degree-one components meet: graph side")
eq(linear_residual.subs({x: linear_meeting, a: linear_meeting_a}), 0,
   "degree-one components meet: residual side")

# Minimal reducible boundary q=k*x.  Work independently with k=1 after the
# universal factorization.  Its conic admits the parameter
#   C=9/(T^2+2), V=9T/(T^2+2).
# The omitted parameter values are infinity (the absent C=0 ramification
# point) and the two roots of T^2+2 (the two C=infinity points).
kappa = sp.symbols("kappa", nonzero=True)
v, T, eta = sp.symbols("v T eta")
minimal = sp.factor(branch.subs(q, kappa * x))
residual = sp.factor(minimal / x)
eq(minimal, x * residual, "minimal branch factorization")
eq(residual.subs(x, 0), -54 * kappa, "vertical and residual components are comaximal")
eq(sp.discriminant(residual, a), -8 * x * (9 * kappa - 2 * x) ** 3,
   "minimal residual discriminant")

residual_one = residual.subs(kappa, 1)
conic = v**2 - x * (9 - 2 * x)
x_param = sp.factor(9 / (T**2 + 2))
v_param = sp.factor(9 * T / (T**2 + 2))
eq(conic.subs({x: x_param, v: v_param}), 0, "residual conic parameterization")

ra = sp.Poly(residual_one, a).coeff_monomial(a**2)
rb = sp.Poly(residual_one, a).coeff_monomial(a)
M = 9 - 2 * x
a_map = sp.factor((-rb + eta * M * v) / (2 * ra))
mapped_numerator = sp.together(residual_one.subs(a, a_map)).as_numer_denom()[0]
mapped_remainder = sp.rem(sp.Poly(mapped_numerator, eta), sp.Poly(eta**2 + 8, eta)).as_expr()
mapped_remainder = sp.factor(mapped_remainder.subs(v**2, x * M))
eq(mapped_remainder, 0, "quadratic normalization map")
a_param = sp.factor(a_map.subs({x: x_param, v: v_param}))
eq(a_param, -(T**3 * eta + 6 * T**2 + 4) / (6 * (T**2 + 2)),
   "parameterized residual A coordinate")
gate(sp.degree(sp.together(a_param).as_numer_denom()[0], T) == 3,
     "A has a pole at parameter infinity")
gate(sp.degree(sp.together(a_param).as_numer_denom()[1], T) == 2,
     "finite denominator is T^2+2")
gate(sp.discriminant(T**2 + 2, T) != 0, "two distinct parameter roots at C infinity")
gate(sp.limit(x_param, T, sp.oo) == 0 and sp.limit(v_param, T, sp.oo) == 0,
     "parameter infinity is the missing C=V=0 point")
gate(residual_one.subs(x, 0) != 0, "residual affine curve omits C=0")

# All origin-vanishing profiles.  The common C-order is read coefficientwise.
# The independent table includes both possible cancellations at m=2 and
# checks the A-infinity corner after the exact vertical power is removed.
def c_order(poly: sp.Expr) -> int:
    return valuation_at_zero(poly, x)


def vertical_packet(profile: sp.Expr, expected_order: int, label: str) -> tuple[int, int, int, int]:
    specialized = sp.expand(branch.subs(q, profile))
    coeffs = sp.Poly(specialized, a)
    coeff_orders = tuple(c_order(coeffs.coeff_monomial(a**degree)) for degree in (2, 1, 0))
    actual_order = min(coeff_orders)
    gate(actual_order == expected_order, f"{label}: exact common vertical order")
    residual_local = sp.expand(specialized / x**actual_order)
    residual_coeffs = sp.Poly(residual_local, a)
    d2 = sp.expand(residual_coeffs.coeff_monomial(a**2))
    d1 = sp.expand(residual_coeffs.coeff_monomial(a))
    d0 = sp.expand(residual_coeffs.coeff_monomial(1))
    gate(d2.subs(x, 0) == 0, f"{label}: residual contains A-infinity corner")
    gate(d2 != 0, f"{label}: projective A-infinity divisor is not a component")
    gate(any(coefficient.subs(x, 0) != 0 for coefficient in (d2, d1, d0)),
         f"{label}: residual is not vertically divisible")
    projective_chart = sp.expand(d2 + d1 * t + d0 * t**2)
    eq(projective_chart.subs({x: 0, t: 0}), 0, f"{label}: corner lies on residual closure")
    return (actual_order, *coeff_orders)


vertical_rows = [
    ("m1", x * (2 + 3 * x), 1),
    ("m2_generic", x**2 * (2 + x), 2),
    ("m2_middle_cancel", x**2 * (sp.Rational(4, 27) + 5 * x), 2),
    ("m2_exceptional", x**2 * (sp.Rational(1, 6) + 5 * x + 7 * x**2), 3),
    ("m3", x**3 * (2 + x), 2),
    ("m4", x**4 * (3 + 2 * x), 2),
    ("m7", x**7 * (5 + x + x**2), 2),
]
vertical_audit_rows = [vertical_packet(profile, expected_order, label)
                       for label, profile, expected_order in vertical_rows]

# Symbolic leading coefficients freeze why the exceptional m=2 order is
# exactly three: the A coefficient is -C^3 at a(0)=1/6, regardless of higher
# profile jets.
a0, a1 = sp.symbols("a0 a1", nonzero=True)
m1_symbolic = sp.Poly(sp.expand(branch.subs(q, x * (a0 + a1 * x))), a)
eq(sp.expand(m1_symbolic.coeff_monomial(1)).coeff(x, 1), -54 * a0,
   "m=1 leading vertical coefficient")
m2_symbolic = sp.Poly(sp.expand(branch.subs(q, x**2 * (a0 + a1 * x))), a)
eq(sp.expand(m2_symbolic.coeff_monomial(1)).coeff(x, 2), 9 - 54 * a0,
   "m=2 generic constant leading coefficient")
eq(sp.expand(m2_symbolic.coeff_monomial(a)).coeff(x, 3), 8 - 54 * a0,
   "m=2 A leading coefficient")
eq(sp.expand(m2_symbolic.coeff_monomial(a)).coeff(x, 3).subs(a0, sp.Rational(1, 6)), -1,
   "m=2 exceptional order-three certificate")

semantic_packet = (
    "THM-3850 independent audit",
    "irreducible nonconstant profile branch",
    "projective normalization W^2=squarefree(9b-2C^2)",
    "deleted points are one per distinct profile root plus all infinity points",
    "conjugate numerator product gives exact 2*ord(b) cancellation",
    "minimal b=kappa*C residual is P1 minus three points",
    "all degree-one profiles classified; the other reducible locus is A1 plus G_m",
    "primitive reducible profiles split into Laurent denominator graphs",
    "origin-vanishing profiles force a residual component through finite-C/A-infinity and C-infinity",
    "every nonconstant profile has a nonpolynomial branch component",
    "notation should distinguish the affine hyperelliptic chart from its projective completion",
)

print("THM3850_INDEPENDENT_AUDIT PASS")
print("CONJUGATE_PRODUCT", sp.factor(conjugate_product))
print("LOCAL_MULTIPLICITY_ROWS", local_rows)
print("PROFILE_ROWS", profile_rows)
print("REPEATED_L_CONTROL", "L=(C-1)^2(C+1); all finite discriminant points retained")
print("SQUARE_L_HOSTILE_FACTORS", len(square_factors))
print("DEGREE_ONE_BONUS", "beta=-9*kappa^2/8 gives A1 plus G_m; all other beta!=0 are irreducible")
print("DEGREE_ONE_COMPONENT_MEETING", "C=9*kappa/4, A=-2/(3*kappa)")
print("MINIMAL_CONIC_PARAMETER", "C=9/(T^2+2), V=9T/(T^2+2)")
print("MINIMAL_A_PARAMETER", a_param)
print("MINIMAL_DELETED_PARAMETERS", "T=infinity and the two roots of T^2+2")
print("PRIMITIVE_REDUCIBLE_FACTOR", "Delta=-(E_u*E_v)/3; components k[C,u^-1], k[C,v^-1]")
print("VERTICAL_AUDIT_ROWS", vertical_audit_rows)
print("ALL_PROFILE_CONCLUSION", "at least one component has at least two punctures")
print("AUDIT_VERDICT", "all-profile mathematics passes; repair H_b/projective-completion notation and 'the first' wording")
print("SEMANTIC_SHA256", hashlib.sha256(repr(semantic_packet).encode()).hexdigest())
print("GATES", GATES)
