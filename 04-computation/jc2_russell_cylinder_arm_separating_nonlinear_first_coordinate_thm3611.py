#!/usr/bin/env python3
"""Exact controls for provisional THM-3611 central-transverse first coordinates.

The universal statement is proof-driven.  This companion checks the formal
arm division and jets, the arbitrary-P transport determinant, direct source
specializations, the boundary/first-coefficient degree mechanisms, a sharp
nonlinear formal endpoint, the localized quotient guardrail, algebraic
descent, and both sharp controls (vanishing derivative and collapsed arms).
"""

import sympy as sp


CHECKS = 0


def require(label, condition):
    """Record one deterministic exact gate."""
    global CHECKS
    if condition != True:
        raise RuntimeError(f"FAIL {label}: {condition}")
    CHECKS += 1


def zero(expression):
    """Exact rational-function zero test."""
    return sp.cancel(expression) == 0


def same(left, right):
    """Exact rational-function equality test."""
    return zero(left - right)


def jacobian(first, second, first_var, second_var):
    """Ordinary two-variable Jacobian in the displayed variable order."""
    return sp.expand(
        sp.diff(first, first_var) * sp.diff(second, second_var)
        - sp.diff(first, second_var) * sp.diff(second, first_var)
    )


print("THM-3611 exact companion -- provisional central-transverse nonlinear coordinate rigidity")
print("status=finite exact controls; universal formal/localized proof is proof-driven")


print("SECTION Russell compiler and polynomial graph coordinates")
x, q, w = sp.symbols("x q w")
D = 1 + x**2 * q
a = q / D**2
c = x * D * (D + 2)
b = (D - 1) * (D + 2) ** 2
e = q * (D + 3)
A_graph = sp.cancel(a * c + w)
B_graph = sp.expand(b + c * w)
S_graph = sp.expand(((b + 2) * (e + 3 * w**2) + c * w * (3 * e + w**2)) / 8)

require("compiler b=ac^2", same(b, a * c**2))
require("compiler e=a(b+4)", same(e, a * (b + 4)))
require("compiler determinant", same(jacobian(a, c, x, q), -3))
require("B root-zero arm", same(B_graph, c * A_graph))
require(
    "S root-zero arm",
    same(S_graph, a + sp.Rational(3, 4) * A_graph**2 + c * A_graph**3 / 8),
)
require("B graph polynomial", sp.denom(sp.cancel(B_graph)) == 1)
require("S graph polynomial", sp.denom(sp.cancel(S_graph)) == 1)
print("PASS compiler_gates=7 determinant=-3")


print("SECTION nonlinear P division and central jets with mixed arm labels")
T, C, Avec = sp.symbols("T C Avec")
P_templates = (
    T + T**2,
    3 + 2 * T - C + T**3 + C * T**2 + 2 * C**2,
    -2 - 3 * T + 4 * C + C * T + C**2 * T**2 + T**2,
    5 + T - 2 * C + 2 * T**2 - C * T**2 + T**4,
    2 + T * (T + 4) * (1 + T) + C * (1 + T**2 + C * T),
)
Xi_templates = []
for index, P_template in enumerate(P_templates):
    p0_value = P_template.subs({T: 0, C: 0})
    alpha_value = sp.diff(P_template, T).subs({T: 0, C: 0})
    r_value = sp.diff(P_template, C).subs({T: 0, C: 0})
    numerator = sp.expand(P_template.subs({T: C * Avec}) - p0_value)
    Xi_value = sp.cancel(numerator / C)
    require(f"Xi exact division row={index}", same(numerator, C * Xi_value))
    require(f"Xi polynomial row={index}", sp.denom(Xi_value) == 1)
    require(f"central alpha nonzero row={index}", alpha_value != 0)
    require(
        f"central affine jet row={index}",
        sp.expand(Xi_value.subs(C, 0)) == alpha_value * Avec + r_value,
    )
    require(f"nonlinear first coordinate row={index}", sp.degree(P_template, T) >= 2)
    Xi_templates.append(sp.expand(Xi_value))

require(
    "separated-arm control",
    P_templates[0].subs({T: -4, C: 0}) != P_templates[0].subs({T: 0, C: 0}),
)
require(
    "identified-arm control",
    P_templates[-1].subs({T: -4, C: 0}) == P_templates[-1].subs({T: 0, C: 0}),
)
print(f"PASS nonlinear_P_rows={len(P_templates)} Xi_division=yes arm_labels=mixed")


print("SECTION abstract psi transport determinant and K_U cancellation")
cc, ZZ, ZZ_a, ZZ_c = sp.symbols("cc ZZ ZZ_a ZZ_c")
AA, psi_Z, psi_c = sp.symbols("AA psi_Z psi_c")
K_U, K_c = sp.symbols("K_U K_c")
QQ = sp.Rational(3, 2) * AA + 3 * cc * AA**2 / 8
LL_a = cc * ZZ_a
LL_c = ZZ + cc * ZZ_c
AA_a = psi_Z * ZZ_a
AA_c = psi_Z * ZZ_c + psi_c
SS_a = 1 + QQ * AA_a
SS_c = AA**3 / 8 + QQ * AA_c
MM_a = SS_a + K_U * LL_a
MM_c = SS_c + K_U * LL_c + K_c
det_ac = sp.expand(LL_a * MM_c - LL_c * MM_a)
VV = sp.expand(ZZ * QQ * psi_Z - cc * (AA**3 / 8 + QQ * psi_c + K_c))
require("abstract arbitrary-P determinant", same(det_ac, -(ZZ + cc * ZZ_c + VV * ZZ_a)))
require("abstract K_U cancellation", not sp.cancel(det_ac).has(K_U))

alpha, rr = sp.symbols("alpha rr", nonzero=True)
boundary_A = (ZZ - rr) / alpha
boundary_Q = sp.Rational(3, 2) * boundary_A
boundary_V = (ZZ * boundary_Q / alpha)
require(
    "abstract boundary velocity",
    same(boundary_V, 3 * ZZ * (ZZ - rr) / (2 * alpha**2)),
)
print("PASS abstract_transport_gates=3 sign=plus3 K_U_cancelled=yes")


print("SECTION direct polynomial-source arbitrary-P rows")
U = sp.symbols("U")
K_templates = (
    U**2 + C * U + C**2,
    U**3 - 2 * C * U + C,
    2 * U**2 * C + U - C**3,
    U**2 + C * U - 2 * C,
)
source_rows = (
    (P_templates[0], Xi_templates[0], K_templates[0], sp.Integer(0)),
    (P_templates[1], Xi_templates[1], K_templates[1], x + q),
    (P_templates[2], Xi_templates[2], K_templates[2], x * q + x**2),
    (P_templates[-1], Xi_templates[-1], K_templates[3], q),
)

for index, (P_template, Xi_template, K_template, h_source) in enumerate(source_rows):
    p0_value = P_template.subs({T: 0, C: 0})
    A_source = sp.cancel(A_graph.subs(w, h_source))
    B_source = sp.expand(B_graph.subs(w, h_source))
    S_source = sp.expand(S_graph.subs(w, h_source))
    L_source = sp.expand(P_template.subs({T: B_source, C: c}))
    M_source = sp.expand(S_source + K_template.subs({U: L_source, C: c}))
    Z_source = sp.cancel(Xi_template.subs({Avec: A_source, C: c}))
    require(
        f"localized quotient identity row={index}",
        same(c * Z_source, L_source - p0_value),
    )

    Xi_A_source = sp.diff(Xi_template, Avec).subs({Avec: A_source, C: c})
    Xi_c_source = sp.diff(Xi_template, C).subs({Avec: A_source, C: c})
    psi_Z_source = sp.cancel(1 / Xi_A_source)
    psi_c_source = sp.cancel(-Xi_c_source / Xi_A_source)
    Z_a_source = sp.cancel(-jacobian(Z_source, c, x, q) / 3)
    Z_c_source = sp.cancel(-jacobian(a, Z_source, x, q) / 3)
    Q_source = sp.cancel(sp.Rational(3, 2) * A_source + 3 * c * A_source**2 / 8)
    K_c_source = sp.diff(K_template, C).subs({U: L_source, C: c})
    V_source = sp.cancel(
        Z_source * Q_source * psi_Z_source
        - c * (A_source**3 / 8 + Q_source * psi_c_source + K_c_source)
    )
    transported = sp.cancel(3 * (Z_source + c * Z_c_source + V_source * Z_a_source))
    require(
        f"direct transport row={index}",
        same(jacobian(L_source, M_source, x, q), transported),
    )
    require(f"direct L polynomial row={index}", sp.denom(sp.cancel(L_source)) == 1)
    require(f"direct M polynomial row={index}", sp.denom(sp.cancel(M_source)) == 1)

print(f"PASS direct_source_rows={len(source_rows)} localized_Z_handled=yes")


print("SECTION boundary and first-formal-coefficient degree gates")
av = sp.symbols("av")
boundary_rows = 0
for degree in range(1, 10):
    f_poly = 2 * av**degree + (av ** (degree - 1) if degree > 1 else 0) - 3
    alpha_value = sp.Rational(degree + 1, degree + 2)
    r_value = sp.Rational(2 - degree, degree + 3)
    boundary = sp.Poly(
        f_poly
        + 3 * f_poly * (f_poly - r_value) * sp.diff(f_poly, av) / (2 * alpha_value**2),
        av,
    )
    require(f"boundary degree d={degree}", boundary.degree() == 3 * degree - 1)
    require(
        f"boundary leading coefficient d={degree}",
        boundary.LC() == 3 * degree * sp.Integer(2) ** 3 / (2 * alpha_value**2),
    )
    boundary_rows += 1

formal_rows = 0
for order in range(1, 7):
    for degree in range(0, 7):
        z_poly = 5 * av**degree + (3 * av ** (degree - 1) if degree > 0 else 0)
        velocity = sp.Rational(2 * order - 5, order + 4)
        row = sp.Poly((order + 1) * z_poly + velocity * sp.diff(z_poly, av), av)
        require(f"formal degree N={order} d={degree}", row.degree() == degree)
        require(
            f"formal leading coefficient N={order} d={degree}",
            row.LC() == 5 * (order + 1),
        )
        formal_rows += 1

print(f"PASS boundary_rows={boundary_rows} formal_first_rows={formal_rows}")


print("SECTION nonlinear sharp formal endpoint P=T+T^2")
aa, c0, tt = sp.symbols("aa c0 tt")
A_endpoint = sp.cancel((sp.sqrt(1 + 4 * c0 * tt) - 1) / (2 * c0))
require("endpoint algebraic equation", same(A_endpoint + c0 * A_endpoint**2, tt))
endpoint_series = sp.series(A_endpoint, c0, 0, 6).removeO().expand()
expected_series = (
    tt - c0 * tt**2 + 2 * c0**2 * tt**3 - 5 * c0**3 * tt**4
    + 14 * c0**4 * tt**5 - 42 * c0**5 * tt**6
)
require("endpoint Catalan series", endpoint_series == expected_series)
L_endpoint = c0 * tt
M_endpoint = (
    aa
    + sp.Rational(3, 4) * A_endpoint**2
    + c0 * A_endpoint**3 / 8
    + L_endpoint**2
    + c0 * L_endpoint
)
require("endpoint ac determinant", same(jacobian(L_endpoint, M_endpoint, aa, c0), -tt))
require("endpoint source determinant", same(-3 * jacobian(L_endpoint, M_endpoint, aa, c0), 3 * tt))
require("endpoint nonlinear P value", (T + T**2).subs(T, -4) == 12)
print("PASS nonlinear_formal_endpoint=exact Catalan_terms=6 determinant=3t")


print("SECTION localized quotient, algebraic descent, and two-arm contradiction")
P_control = T + T**2
p0_control = 0
Xi_control = Avec + C * Avec**2
A_control = sp.cancel(a * c)
Z_control = sp.cancel(Xi_control.subs({Avec: A_control, C: c}))
L_control = sp.expand(P_control.subs({T: b, C: c}))
require("control quotient clearing", same(c * Z_control, L_control - p0_control))
require("control Z not global polynomial", sp.denom(Z_control) != 1)
require("control L global polynomial", sp.denom(sp.cancel(L_control)) == 1)

# Birational field coordinates: v=x^2 q and c=x(v+1)(v+3).
vv, cs = sp.symbols("vv cs")
p_v = sp.expand((vv + 1) * (vv + 3))
x_inverse = sp.cancel(cs / p_v)
q_inverse = sp.cancel(vv * p_v**2 / cs**2)
v_recovered = sp.cancel(x_inverse**2 * q_inverse)
c_recovered = sp.cancel(x_inverse * (v_recovered + 1) * (v_recovered + 3))
require("birational inverse recovers v", same(v_recovered, vv))
require("birational inverse recovers c", same(c_recovered, cs))
v_source = x**2 * q
require("source c=xp(v)", same(c, x * ((v_source + 1) * (v_source + 3))))

# Finite algebraic/elimination control for P(T)=T+T^2.
Balg = sp.symbols("Balg")
algebraic_relation = sp.Poly(Balg**2 + Balg - tt * cs, Balg)
require("finite algebraic relation degree", algebraic_relation.degree() == 2)
require(
    "finite algebraic discriminant",
    sp.discriminant(algebraic_relation.as_expr(), Balg) == 1 + 4 * tt * cs,
)
B_endpoint = sp.cancel(c0 * A_endpoint)
require("formal endpoint solves B algebraic relation", same(B_endpoint**2 + B_endpoint, tt * c0))

# Exact Bezout controls model C[x,q] intersect C(c)=C[c].
zv = sp.symbols("zv")
bezout_rows = (
    (zv**2 + 1, (zv - 2) * (zv + 3)),
    (zv**3 + 2 * zv + 2, zv**2 - zv + 1),
    (2 * zv**2 + 3 * zv + 5, zv**3 - 2),
)
for index, (f_poly, g_poly) in enumerate(bezout_rows):
    bezout_f, bezout_g, gcd_value = sp.gcdex(f_poly, g_poly, zv)
    require(f"Bezout gcd row={index}", gcd_value == 1)
    require(
        f"Bezout identity row={index}",
        sp.expand(bezout_f * f_poly + bezout_g * g_poly) == 1,
    )
    require(
        f"Bezout substitution row={index}",
        same(
            (bezout_f * f_poly + bezout_g * g_poly).subs(zv, c),
            1,
        ),
    )

# The descended polynomial s(c) would have incompatible values over c=0.
central_point = {x: 0, q: 0}
hostile_point = {x: 1, q: -1}
require("central point c=0", c.subs(central_point) == 0)
require("central point B=0", B_graph.subs(central_point) == 0)
require("hostile point D=0", D.subs(hostile_point) == 0)
require("hostile point c=0", c.subs(hostile_point) == 0)
require("hostile point B=-4", B_graph.subs(hostile_point) == -4)
require("two-arm value mismatch", B_graph.subs(hostile_point) - B_graph.subs(central_point) == -4)
print(
    "PASS quotient_guardrail=localized_only birational_inverse=yes "
    f"Bezout_rows={len(bezout_rows)} arm_values=0,-4"
)


print("SECTION vanishing-central-derivative formal hostile")
phi, cv = sp.symbols("phi cv")
v_of_phi = -phi / 2 - sp.Rational(27, 32) * phi**2 - phi**3 / 8
a_hostile = v_of_phi / cv**2
A_hostile = phi / cv
Z_hostile = 1 + phi
L_hostile = cv + cv**2 * A_hostile
S_hostile = a_hostile + sp.Rational(3, 4) * A_hostile**2 + cv * A_hostile**3 / 8
j_ls = jacobian(L_hostile, S_hostile, phi, cv)
j_ac = jacobian(a_hostile, cv, phi, cv)
hostile_factor = (6 * phi**2 + 27 * phi + 8) / (16 * cv**2)
require("alpha-zero hostile P_T", sp.diff(C * T + C, T).subs({T: 0, C: 0}) == 0)
require("alpha-zero hostile Z nonconstant", sp.diff(Z_hostile, phi) == 1)
require("alpha-zero hostile numerator determinant", same(j_ls, hostile_factor))
require("alpha-zero hostile coordinate determinant", same(j_ac, -hostile_factor))
require("alpha-zero hostile Keller determinant", same(-3 * j_ls / j_ac, 3))
print("PASS alpha_zero_method_hostile=nonconstant_formal_Z determinant=3")


print("SECTION collapsed-arm coordinate and complete square no-go")
Bv, Cv, p0v, nu, mu, kval = sp.symbols("Bv Cv p0v nu mu kval")
P_star = Cv + Bv * (Bv + 4)
require("collapsed coordinate inverse", sp.expand(P_star - Bv * (Bv + 4)) == Cv)
require("collapsed central derivative", sp.diff(P_star, Bv).subs({Bv: 0, Cv: 0}) == 4)
require("collapsed central arm value", P_star.subs({Bv: 0, Cv: 0}) == 0)
require("collapsed hostile arm value", P_star.subs({Bv: -4, Cv: 0}) == 0)

Y_graph = sp.expand(c * e + (2 * b + 4) * w + c * w**2)
require("B(B+4)=cY", same(B_graph * (B_graph + 4), c * Y_graph))
require("collapsed L=c(1+Y)", same(c + B_graph * (B_graph + 4), c * (1 + Y_graph)))
require(
    "square identity",
    same((c * w + b + 2) ** 2, 4 + c * Y_graph),
)
require("ordinary degree c=7", sp.Poly(c, x, q).total_degree() == 7)
leading_c = sp.Poly(c, x, q).terms()[0]
require("ordinary leading monomial c", leading_c[0] == (5, 2) and leading_c[1] == 1)
require("nonzero-k square degree odd", sp.Poly(4 + kval * c, x, q).total_degree() == 7)
require("k=0 first branch fails modulo D", same(b.subs(q, -1 / x**2), -4))
require("k=0 second branch fails modulo x", (b + 4).subs(x, 0) == 4)

P_family_graph = p0v + nu * c + mu * B_graph * (B_graph + 4)
require(
    "collapsed affine family factorization",
    same(P_family_graph, p0v + c * (nu + mu * Y_graph)),
)
require("collapsed family central derivative", sp.diff(p0v + nu * Cv + mu * Bv * (Bv + 4), Bv).subs({Bv: 0, Cv: 0}) == 4 * mu)
require("collapsed family arm equality", same((p0v + nu * Cv + mu * Bv * (Bv + 4)).subs({Bv: -4, Cv: 0}), p0v))
print("PASS collapsed_arm_coordinate=independent_control square_no_go=complete affine_family=closed")


print(f"PASS total_exact_gates={CHECKS}")
print("RESULT provisional_controls_passed; independent hostile audit still required")
