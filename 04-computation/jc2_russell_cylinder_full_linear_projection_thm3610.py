#!/usr/bin/env python3
"""Exact controls for proved and hostile-audited THM-3610.

The proof is all-degree and uses a formal-arm argument plus the cited
injectivity-on-one-line theorem.  This companion verifies every displayed
algebraic identity, the complete finite row-space trichotomy over a hostile
coefficient box, the boundary and first-formal-coefficient invoices, the
three collision points, and the collision-losing coordinate-plane controls.
"""

from itertools import product

import sympy as sp


CHECKS = 0


def require(label, condition):
    """Record one exact gate and fail with a stable label."""
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
    """Ordinary two-variable Jacobian in the displayed order."""
    return sp.expand(
        sp.diff(first, first_var) * sp.diff(second, second_var)
        - sp.diff(first, second_var) * sp.diff(second, first_var)
    )


print("THM-3610 exact companion -- proved full linear-projection gate")
print("status=finite exact controls; all-degree formal and one-line steps are proof-driven")


print("SECTION compiler and four Russell arm-coordinate identities")
x, q, w = sp.symbols("x q w")
D = 1 + x**2 * q
a = q / D**2
c = x * D * (D + 2)
b = (D - 1) * (D + 2) ** 2
e = q * (D + 3)

B = sp.expand(b + c * w)
C = c
Y = sp.expand(c * e + (2 * b + 4) * w + c * w**2)
S = sp.expand(((b + 2) * (e + 3 * w**2) + c * w * (3 * e + w**2)) / 8)
A = sp.cancel(a * c + w)

require("compiler b=ac^2", same(b, a * c**2))
require("compiler e=a(b+4)", same(e, a * (b + 4)))
require("compiler arm determinant", same(jacobian(a, c, x, q), -3))
require("Russell B=cA", same(B, c * A))
require("Russell C=c", same(C, c))
require("Russell Y=4A+cA^2", same(Y, 4 * A + c * A**2))
require("Russell S arm form", same(S, a + 3 * A**2 / 4 + c * A**3 / 8))
for label, expression in (("B", B), ("C", C), ("Y", Y), ("S", S)):
    require(f"{label} graph polynomial", sp.denom(sp.cancel(expression)) == 1)
print("PASS compiler_and_four_coordinate_gates=11 determinant=-3")


print("SECTION abstract r=1 determinant identities")
aa, cc = sp.symbols("aa cc")
ZZ, ZZ_a, ZZ_c = sp.symbols("ZZ ZZ_a ZZ_c")
rho, sigma, tau = sp.symbols("rho sigma tau")
AA = ZZ - rho

g_y = 4 * AA + cc * AA**2 + tau * cc
g_y_Z = sp.diff(g_y, ZZ)
g_y_c = sp.diff(g_y, cc)
L_a = cc * ZZ_a
L_c = ZZ + cc * ZZ_c
MY_a = g_y_Z * ZZ_a
MY_c = g_y_Z * ZZ_c + g_y_c
det_y = sp.expand(L_a * MY_c - L_c * MY_a)
require("Y-only ac determinant", same(det_y, (cc * g_y_c - ZZ * g_y_Z) * ZZ_a))
require("Y-only xq sign", same(-3 * det_y, 3 * (ZZ * g_y_Z - cc * g_y_c) * ZZ_a))

g_s = (
    3 * AA**2 / 4
    + cc * AA**3 / 8
    + sigma * (4 * AA + cc * AA**2)
    + tau * cc
)
g_s_Z = sp.diff(g_s, ZZ)
g_s_c = sp.diff(g_s, cc)
MS_a = 1 + g_s_Z * ZZ_a
MS_c = g_s_Z * ZZ_c + g_s_c
det_s = sp.expand(L_a * MS_c - L_c * MS_a)
E_s = sp.expand(ZZ + cc * ZZ_c + (ZZ * g_s_Z - cc * g_s_c) * ZZ_a)
require("S+sigmaY ac determinant", same(det_s, -E_s))
require("S+sigmaY xq sign", same(-3 * det_s, 3 * E_s))
print("PASS abstract_determinants=4")


print("SECTION direct source transport identities")
source_graphs = (
    sp.Integer(0),
    sp.Integer(2),
    x,
    q,
    x * q + x**2 - 3 * q**2,
    x**5 * q**2,
)
transport_rows = 0
for index, h in enumerate(source_graphs):
    A_h = sp.cancel(a * c + h)
    Z_h = sp.cancel(A_h + rho)
    B_h = sp.expand(B.subs(w, h))
    Y_h = sp.expand(Y.subs(w, h))
    S_h = sp.expand(S.subs(w, h))
    L_h = sp.expand(B_h + rho * c)

    M_y_h = sp.expand(Y_h + tau * c)
    g_y_h = sp.cancel((4 * (Z_h - rho) + c * (Z_h - rho) ** 2 + tau * c))
    K_y_h = sp.cancel(
        Z_h * sp.diff(g_y, ZZ).subs({ZZ: Z_h, cc: c})
        - c * sp.diff(g_y, cc).subs({ZZ: Z_h, cc: c})
    )
    transported_y = sp.cancel(-K_y_h * jacobian(Z_h, c, x, q))
    require(
        f"direct Y transport row={index}",
        same(jacobian(L_h, M_y_h, x, q), transported_y),
    )
    require(f"direct Y g row={index}", same(M_y_h, g_y_h))

    M_s_h = sp.expand(S_h + sigma * Y_h + tau * c)
    g_s_h = sp.cancel(g_s.subs({ZZ: Z_h, cc: c}))
    K_s_h = sp.cancel(
        Z_h * sp.diff(g_s, ZZ).subs({ZZ: Z_h, cc: c})
        - c * sp.diff(g_s, cc).subs({ZZ: Z_h, cc: c})
    )
    transported_s = sp.cancel(
        3 * Z_h
        + c * jacobian(Z_h, a, x, q)
        - K_s_h * jacobian(Z_h, c, x, q)
    )
    require(
        f"direct S+sigmaY transport row={index}",
        same(jacobian(L_h, M_s_h, x, q), transported_s),
    )
    require(f"direct S+sigmaY g row={index}", same(M_s_h, a + g_s_h))
    transport_rows += 4
print(f"PASS direct_transport_gates={transport_rows}")


print("SECTION x=0 boundary equations and formal coefficient invoice")
arm_a, arm_c, t, lead = sp.symbols("arm_a arm_c t lead", nonzero=True)
boundary_rows = 0
for degree in range(1, 13):
    f = lead * arm_a**degree
    y_ode = sp.expand(4 * f * sp.diff(f, arm_a) - t)
    require(
        f"Y boundary leading degree={degree}",
        sp.Poly(y_ode, arm_a).coeff_monomial(arm_a ** (2 * degree - 1))
        == 4 * degree * lead**2,
    )

    s_ode = sp.expand(
        f
        + f * (sp.Rational(3, 2) * (f - rho) + 4 * sigma) * sp.diff(f, arm_a)
        - t
    )
    require(
        f"S boundary leading degree={degree}",
        sp.Poly(s_ode, arm_a).coeff_monomial(arm_a ** (3 * degree - 1))
        == sp.Rational(3, 2) * degree * lead**3,
    )
    boundary_rows += 2

formal_rows = 0
for order in range(1, 9):
    for degree in range(0, 10):
        z = lead * arm_a**degree
        Z_formal = t + z * arm_c**order
        A_formal = Z_formal - rho
        g_formal = (
            3 * A_formal**2 / 4
            + arm_c * A_formal**3 / 8
            + sigma * (4 * A_formal + arm_c * A_formal**2)
            + tau * arm_c
        )
        # Differentiate the abstract polynomial before substituting Z_formal.
        K_formal = sp.expand(
            (ZZ * g_s_Z - cc * g_s_c).subs({ZZ: Z_formal, cc: arm_c})
        )
        E_formal = sp.expand(
            Z_formal
            + arm_c * sp.diff(Z_formal, arm_c)
            + K_formal * sp.diff(Z_formal, arm_a)
            - t
        )
        coefficient = sp.expand(E_formal).coeff(arm_c, order)
        k0 = t * (sp.Rational(3, 2) * (t - rho) + 4 * sigma)
        expected = sp.expand((order + 1) * z + k0 * sp.diff(z, arm_a))
        require(
            f"formal first coefficient order={order} degree={degree}",
            same(coefficient, expected),
        )
        if degree >= 1:
            require(
                f"formal leading coefficient order={order} degree={degree}",
                sp.Poly(expected, arm_a).coeff_monomial(arm_a**degree)
                == (order + 1) * lead,
            )
        formal_rows += 1 + (degree >= 1)
print(f"PASS boundary_degree_rows={boundary_rows} formal_rows={formal_rows}")


print("SECTION complete {-1,0,1} rank-two row-space census")
rank_two = 0
r_counts = {0: 0, 1: 0, 2: 0}
r1_fixed_c = 0
r1_y_only = 0
r1_s_mixed = 0
for entries in product((-1, 0, 1), repeat=8):
    row_one = entries[:4]
    row_two = entries[4:]
    minors = []
    for i in range(4):
        for j in range(i + 1, 4):
            minors.append(row_one[i] * row_two[j] - row_one[j] * row_two[i])
    if all(value == 0 for value in minors):
        continue

    u1 = row_one[2:]
    u2 = row_two[2:]
    u_det = u1[0] * u2[1] - u1[1] * u2[0]
    if u1 == (0, 0) and u2 == (0, 0):
        projection_rank = 0
    elif u_det == 0:
        projection_rank = 1
    else:
        projection_rank = 2
    r_value = 2 - projection_rank
    require(f"row-space r range matrix={entries}", r_value in (0, 1, 2))
    r_counts[r_value] += 1

    if r_value == 2:
        require(
            f"r2 equals BC plane matrix={entries}",
            u1 == (0, 0) and u2 == (0, 0),
        )
    elif r_value == 0:
        require(f"r0 transverse U determinant matrix={entries}", u_det != 0)
    else:
        if u1[0] != 0 or u2[0] != 0:
            kernel = (u2[0], -u1[0])
        else:
            kernel = (u2[1], -u1[1])
        intersection = (
            kernel[0] * row_one[0] + kernel[1] * row_two[0],
            kernel[0] * row_one[1] + kernel[1] * row_two[1],
        )
        require(f"r1 nonzero BC intersection matrix={entries}", intersection != (0, 0))
        require(
            f"r1 kernel kills U matrix={entries}",
            kernel[0] * u1[0] + kernel[1] * u2[0] == 0
            and kernel[0] * u1[1] + kernel[1] * u2[1] == 0,
        )
        if intersection[0] == 0:
            r1_fixed_c += 1
        else:
            u_direction = u1 if u1 != (0, 0) else u2
            if u_direction[1] == 0:
                r1_y_only += 1
            else:
                r1_s_mixed += 1
    rank_two += 1
print(
    f"PASS rank_two_matrices={rank_two} r0={r_counts[0]} r1={r_counts[1]} "
    f"r2={r_counts[2]} r1_fixedC={r1_fixed_c} r1_Y={r1_y_only} "
    f"r1_Smixed={r1_s_mixed}"
)


print("SECTION collision, line injectivity, and collision-losing controls")
eta = sp.symbols("eta")
collision_points = (
    (sp.Integer(0), -sp.Rational(3, 4)),
    (sp.Integer(1), sp.Integer(-3)),
    (sp.Integer(-1), sp.Integer(-3)),
)
collision_rows = 0
for index, (x_value, q_value) in enumerate(collision_points):
    substitutions = {x: x_value, q: q_value, w: eta}
    require(f"collision a row={index}", same(a.subs(substitutions), -sp.Rational(3, 4)))
    require(f"collision b row={index}", same(b.subs(substitutions), 0))
    require(f"collision c row={index}", same(c.subs(substitutions), 0))
    require(f"collision e row={index}", same(e.subs(substitutions), -3))
    require(f"collision B row={index}", same(B.subs(substitutions), 0))
    require(f"collision C row={index}", same(C.subs(substitutions), 0))
    require(f"collision Y row={index}", same(Y.subs(substitutions), 4 * eta))
    require(
        f"collision S row={index}",
        same(S.subs(substitutions), -sp.Rational(3, 4) + 3 * eta**2 / 4),
    )
    collision_rows += 8

phi, line_first, line_second = sp.symbols("phi line_first line_second")
line_first = 4 * phi
line_second = q + 3 * phi**2 / 4
require("line reconstruct phi", same(line_first / 4, phi))
require("line reconstruct q", same(line_second - 3 * line_first**2 / 64, q))

plane_q, plane_w = sp.symbols("plane_q plane_w")
Y_x0 = 4 * plane_w
S_x0 = plane_q + 3 * plane_w**2 / 4
require("x=0 positive plane Jacobian", jacobian(Y_x0, S_x0, plane_q, plane_w) == -4)

plane_x = sp.symbols("plane_x")
q0_sub = {x: plane_x, q: 0, w: plane_w}
B_q0 = sp.expand(B.subs(q0_sub))
C_q0 = sp.expand(C.subs(q0_sub))
Y_q0 = sp.expand(Y.subs(q0_sub))
S_q0 = sp.expand(S.subs(q0_sub))
W_inverse = sp.expand(Y_q0 * (B_q0 + 2) / 8 - C_q0 * S_q0)
require("q=0 inverse cylinder coordinate", same(W_inverse, plane_w))
require("q=0 positive plane Jacobian", jacobian(C_q0, W_inverse, plane_x, plane_w) == 3)
print(
    f"PASS collision_rows={collision_rows} line_reconstruction=2 "
    "positive_plane_controls=3"
)


print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
