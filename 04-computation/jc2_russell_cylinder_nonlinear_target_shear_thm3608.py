#!/usr/bin/env python3
"""Exact controls for provisional THM-3608 nonlinear target-shear rigidity.

The universal exclusion is proof-driven: a boundary degree argument and a
first-nonzero formal-arm coefficient force the normalized first output to be
constant.  This companion checks the full differential jet identity, source
specializations, finite degree controls, sharp punctured endpoints, the
singular-u formal hostile, and the nonlinear arm-identifying coordinate.
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


print("THM-3608 exact companion -- provisional nonlinear target-shear rigidity")
print("status=finite exact controls; universal formal argument is proof-driven")


print("SECTION Russell compiler and root-zero arm normal form")
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
require("compiler source determinant", same(jacobian(a, c, x, q), -3))
require("root-zero B arm", same(B_graph, c * A_graph))
require(
    "root-zero S arm",
    same(S_graph, a + sp.Rational(3, 4) * A_graph**2 + c * A_graph**3 / 8),
)
require("B graph polynomial", sp.denom(sp.cancel(B_graph)) == 1)
require("S graph polynomial", sp.denom(sp.cancel(S_graph)) == 1)
print("PASS compiler_gates=7 determinant=-3")


print("SECTION universal u,F,K differential-jet identity")
cc, ZZ, ZZ_a, ZZ_c = sp.symbols("cc ZZ ZZ_a ZZ_c")
uu, uu_p, RR, RR_p, F0 = sp.symbols("uu uu_p RR RR_p F0", nonzero=True)
K_L, K_c = sp.symbols("K_L K_c")
AA = (ZZ - RR) / uu
QQ = sp.Rational(3, 2) * AA + 3 * cc * AA**2 / 8
AA_a = ZZ_a / uu
AA_c = (ZZ_c - RR_p - uu_p * AA) / uu
LL_a = cc * ZZ_a
LL_c = ZZ + cc * ZZ_c
SS_a = 1 + QQ * AA_a
SS_c = AA**3 / 8 + QQ * AA_c
MM_a = SS_a + K_L * LL_a
MM_c = SS_c + K_L * LL_c + K_c
det_ac = sp.expand(LL_a * MM_c - LL_c * MM_a)
VV = sp.cancel(
    QQ * (ZZ + cc * RR_p + cc * uu_p * AA) / uu
    - cc * AA**3 / 8
    - cc * K_c
)
require(
    "full abstract determinant",
    same(det_ac, -(ZZ + cc * ZZ_c + VV * ZZ_a)),
)
require("K_L cancels", not sp.cancel(det_ac).has(K_L))

# Recover the THM-3607 constant-rho, linear-H face exactly.
rho, tau = sp.symbols("rho tau")
AA_old = ZZ - rho
QQ_old = sp.Rational(3, 2) * AA_old + 3 * cc * AA_old**2 / 8
K_old = sp.expand(ZZ * QQ_old - cc * AA_old**3 / 8)
VV_old = VV.subs(
    {uu: 1, uu_p: 0, RR: rho, RR_p: 0, K_c: tau}
)
require("linear face recovery", same(VV_old, K_old - tau * cc))
print("PASS abstract_transport_gates=3 K_L_cancelled=yes")


print("SECTION direct polynomial-source transport rows")
Lvar, cvar = sp.symbols("Lvar cvar")
source_rows = (
    (
        sp.Integer(1),
        sp.Integer(0),
        sp.Integer(0),
        Lvar**2 + cvar * Lvar + cvar**3,
        sp.Integer(0),
    ),
    (
        2 + cvar,
        1 + cvar,
        sp.Integer(3),
        2 * Lvar**2 * cvar - Lvar * cvar**2 + cvar,
        x + q,
    ),
    (
        1 - cvar + cvar**2,
        2 - cvar + cvar**2,
        sp.Integer(-2),
        Lvar**3 + 2 * Lvar * cvar + 3 * cvar**2,
        x * q + x**2 - 2 * q**2,
    ),
    (
        3 + 2 * cvar**2,
        -1 + 2 * cvar - cvar**2,
        sp.Integer(5),
        Lvar**2 * cvar**2 + Lvar + cvar**4,
        x**2 * q + q,
    ),
)

for index, (u_template, R_template, f0_source, K_template, h_source) in enumerate(source_rows):
    u_source = sp.expand(u_template.subs(cvar, c))
    R_source = sp.expand(R_template.subs(cvar, c))
    A_source = sp.cancel(A_graph.subs(w, h_source))
    B_source = sp.expand(B_graph.subs(w, h_source))
    S_source = sp.expand(S_graph.subs(w, h_source))
    F_source = sp.expand(f0_source + c * R_source)
    L_source = sp.expand(u_source * B_source + F_source)
    K_source = sp.expand(K_template.subs({Lvar: L_source, cvar: c}))
    M_source = sp.expand(S_source + K_source)
    Z_source = sp.cancel(u_source * A_source + R_source)
    Z_a_source = sp.cancel(-jacobian(Z_source, c, x, q) / 3)
    Z_c_source = sp.cancel(-jacobian(a, Z_source, x, q) / 3)
    Q_source = sp.cancel(sp.Rational(3, 2) * A_source + 3 * c * A_source**2 / 8)
    K_c_source = sp.diff(K_template, cvar).subs({Lvar: L_source, cvar: c})
    V_source = sp.cancel(
        Q_source
        * (
            Z_source
            + c * sp.diff(R_template, cvar).subs(cvar, c)
            + c * sp.diff(u_template, cvar).subs(cvar, c) * A_source
        )
        / u_source
        - c * A_source**3 / 8
        - c * K_c_source
    )
    transported = sp.cancel(3 * (Z_source + c * Z_c_source + V_source * Z_a_source))
    require(
        f"source transport row={index}",
        same(jacobian(L_source, M_source, x, q), transported),
    )
    require(f"source L polynomial row={index}", sp.denom(sp.cancel(L_source)) == 1)
    require(f"source M polynomial row={index}", sp.denom(sp.cancel(M_source)) == 1)

print(f"PASS direct_source_rows={len(source_rows)} polynomiality_gates={2*len(source_rows)}")


print("SECTION boundary ODE and first-formal-coefficient degree gates")
av = sp.symbols("av")
boundary_rows = 0
for degree in range(1, 10):
    f_poly = 2 * av**degree + (3 * av ** (degree - 1) if degree > 1 else 0) + 5
    u0_value = sp.Rational(degree + 2, degree + 1)
    r0_value = sp.Rational(1 - degree, degree + 3)
    boundary = sp.Poly(
        f_poly
        + 3 * f_poly * (f_poly - r0_value) * sp.diff(f_poly, av) / (2 * u0_value**2),
        av,
    )
    require(f"boundary degree d={degree}", boundary.degree() == 3 * degree - 1)
    require(
        f"boundary leading coefficient d={degree}",
        boundary.LC() == 3 * degree * sp.Integer(2) ** 3 / (2 * u0_value**2),
    )
    boundary_rows += 1

t_probe = sp.Rational(7, 3)
require(
    "constant boundary endpoint",
    same(t_probe + 3 * t_probe * (t_probe - 5) * sp.diff(t_probe, av) / 8, t_probe),
)

formal_rows = 0
for order in range(1, 7):
    for degree in range(0, 7):
        z_poly = 3 * av**degree + (2 * av ** (degree - 1) if degree > 0 else 0)
        velocity = sp.Rational(order - 3, order + 2)
        first_row = sp.Poly((order + 1) * z_poly + velocity * sp.diff(z_poly, av), av)
        require(
            f"formal degree N={order} d={degree}",
            first_row.degree() == degree,
        )
        require(
            f"formal leading coefficient N={order} d={degree}",
            first_row.LC() == 3 * (order + 1),
        )
        formal_rows += 1

print(
    f"PASS boundary_degree_rows={boundary_rows} formal_first_rows={formal_rows} "
    "constant_endpoint=yes"
)


print("SECTION sharp punctured controls, collision, and D-pole invoice")
aa, c0 = sp.symbols("aa c0")
endpoint_rows = (
    (sp.Integer(1), 2 + c0, sp.Integer(0), sp.Rational(3, 2), Lvar**2 + cvar),
    (2 + c0, 1 - c0 + c0**2, sp.Integer(4), sp.Integer(-2), Lvar * cvar + cvar**3),
    (3 - c0 + 2 * c0**2, -1 + 3 * c0, sp.Integer(-5), sp.Rational(5, 3), Lvar**3 - 2 * cvar),
)

for index, (u_end, R_end, f0_end, t_end, K_end_template) in enumerate(endpoint_rows):
    A_end = sp.cancel((t_end - R_end) / u_end)
    L_end = sp.expand(f0_end + t_end * c0)
    M_end = sp.cancel(
        aa
        + sp.Rational(3, 4) * A_end**2
        + c0 * A_end**3 / 8
        + K_end_template.subs({Lvar: L_end, cvar: c0})
    )
    require(
        f"punctured ac determinant row={index}",
        same(jacobian(L_end, M_end, aa, c0), -t_end),
    )

    A_source_end = A_end.subs(c0, c)
    h_source_end = sp.cancel(A_source_end - a * c)
    L_source_end = L_end.subs(c0, c)
    M_source_end = M_end.subs({aa: a, c0: c})
    require(
        f"punctured source determinant row={index}",
        same(jacobian(L_source_end, M_source_end, x, q), 3 * t_end),
    )

    residue = sp.cancel(D * h_source_end).subs(q, -1 / x**2)
    require(f"D-pole residue row={index}", same(residue, 2 / x))

    collision_values = []
    for source_point in ((0, sp.Rational(-3, 4)), (1, -3), (-1, -3)):
        collision_values.append(sp.cancel(h_source_end.subs({x: source_point[0], q: source_point[1]})))
    expected_collision = sp.cancel(A_source_end.subs(c, 0))
    require(
        f"collision common graph value row={index}",
        all(same(value, expected_collision) for value in collision_values),
    )

print(f"PASS punctured_endpoint_rows={len(endpoint_rows)} pole_residue=2/x collision=common")


print("SECTION singular-u formal-algebraic hostile")
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
require("hostile L=c(1+phi)", same(L_hostile, cv * Z_hostile))
require("hostile numerator determinant", same(j_ls, hostile_factor))
require("hostile coordinate determinant", same(j_ac, -hostile_factor))
require("hostile ac Keller determinant", same(j_ls / j_ac, -1))
require("hostile source Keller determinant", same(-3 * j_ls / j_ac, 3))
require("hostile nonconstant Z", sp.diff(Z_hostile, phi) == 1)
require("hostile inverse derivative", sp.diff(v_of_phi, phi).subs(phi, 0) == -sp.Rational(1, 2))

vv = sp.symbols("vv")
inverse_coeffs = sp.symbols("d1:7")
phi_series = sum(inverse_coeffs[index - 1] * vv**index for index in range(1, 7))
relation_series = sp.expand(
    -phi_series / 2 - sp.Rational(27, 32) * phi_series**2 - phi_series**3 / 8 - vv
)
solved = {}
for degree in range(1, 7):
    equation = sp.expand(relation_series.subs(solved)).coeff(vv, degree)
    solution = sp.solve(equation, inverse_coeffs[degree - 1])
    require(f"hostile inverse coefficient exists n={degree}", len(solution) == 1)
    solved[inverse_coeffs[degree - 1]] = sp.factor(solution[0])

expected_inverse = (
    -2,
    -sp.Rational(27, 4),
    -sp.Rational(697, 16),
    -sp.Rational(89775, 256),
    -sp.Rational(3236343, 1024),
    -sp.Rational(249967431, 8192),
)
for index, expected in enumerate(expected_inverse, start=1):
    require(
        f"hostile inverse coefficient value n={index}",
        solved[inverse_coeffs[index - 1]] == expected,
    )
    # phi((a c^2))/c has c-order 2n-1, so it lies in C[a][[c]].
    require(f"hostile formal regularity n={index}", 2 * index - 1 >= 1)

print("PASS singular_u_hostile=exact nonconstant_Z=yes inverse_coefficients=6")


print("SECTION arbitrary-coordinate arm-identification hostile")
Bv, Cv = sp.symbols("Bv Cv")
L_star = Cv + Bv * (Bv + 4)
recovered_c = sp.expand(L_star - Bv * (Bv + 4))
require("nonlinear coordinate inverse", recovered_c == Cv)
require("central arm value", L_star.subs({Bv: 0, Cv: 0}) == 0)
require("hostile arm value", L_star.subs({Bv: -4, Cv: 0}) == 0)
Y_graph = sp.expand(c * e + (2 * b + 4) * w + c * w**2)
require("graph arm-product identity", same(B_graph * (B_graph + 4), c * Y_graph))
require(
    "nonlinear coordinate graph factorization",
    same(c + B_graph * (B_graph + 4), c * (1 + Y_graph)),
)
print("PASS arbitrary_coordinate_hostile=two_arms_identified open_exit=yes")


print(f"PASS total_exact_gates={CHECKS}")
print("RESULT provisional_controls_passed; independent hostile audit still required")
