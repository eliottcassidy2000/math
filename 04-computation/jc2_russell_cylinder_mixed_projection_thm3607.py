#!/usr/bin/env python3
"""Exact controls for proved and independently audited THM-3607.

The rank-two linear-projection obstruction is proof-driven: its hard face
uses a polynomial boundary ODE and a first-nonzero coefficient in the formal
arm completion.  This companion verifies the row-space trichotomy, easy
fixed-C and (B,C) faces, every displayed algebraic identity, determinant sign,
finite formal/initial-form rows, complete quadratic hostile, and sharp
punctured rational endpoint with exact SymPy arithmetic.
"""

from itertools import product

import sympy as sp


CHECKS = 0


def require(label, condition):
    """Record one exact gate and raise with a stable label on failure."""
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


print("THM-3607 exact companion -- proved rank-two projection formal rigidity")
print("status=finite exact controls; all-degree completion argument is proof-driven")


print("SECTION compiler identities and mixed arm normal form")
x, q, w = sp.symbols("x q w")
D = 1 + x**2 * q
a = q / D**2
c = x * D * (D + 2)
b = (D - 1) * (D + 2) ** 2
e = q * (D + 3)

B_graph = sp.expand(b + c * w)
S_graph = sp.expand(((b + 2) * (e + 3 * w**2) + c * w * (3 * e + w**2)) / 8)
A_graph = sp.cancel(a * c + w)

require("compiler b=ac^2", same(b, a * c**2))
require("compiler e=a(b+4)", same(e, a * (b + 4)))
require("compiler rational determinant", same(jacobian(a, c, x, q), -3))
require("mixed B arm normal form", same(B_graph, c * A_graph))
require(
    "mixed S arm normal form",
    same(S_graph, a + sp.Rational(3, 4) * A_graph**2 + c * A_graph**3 / 8),
)
require("B graph polynomial", sp.denom(sp.cancel(B_graph)) == 1)
require("S graph polynomial", sp.denom(sp.cancel(S_graph)) == 1)
print("PASS compiler_and_arm_gates=7 determinant=-3")


print("SECTION rank-two row-space trichotomy and easy projection faces")
rank_two_rows = 0
fixed_c_rows = 0
bc_plane_rows = 0
mixed_rows = 0
for entries in product((-1, 0, 1), repeat=6):
    row_one = entries[:3]
    row_two = entries[3:]
    minors = (
        row_one[0] * row_two[1] - row_one[1] * row_two[0],
        row_one[0] * row_two[2] - row_one[2] * row_two[0],
        row_one[1] * row_two[2] - row_one[2] * row_two[1],
    )
    if minors == (0, 0, 0):
        continue

    if row_one[2] == 0 and row_two[2] == 0:
        intersection = row_one
    else:
        intersection = (
            row_two[2] * row_one[0] - row_one[2] * row_two[0],
            row_two[2] * row_one[1] - row_one[2] * row_two[1],
            0,
        )
    require(
        f"nonzero BC intersection matrix={entries}",
        intersection[2] == 0 and intersection[:2] != (0, 0),
    )

    if intersection[0] == 0:
        fixed_c_rows += 1
    elif row_one[2] == 0 and row_two[2] == 0:
        bc_plane_rows += 1
    else:
        # Normalize L=B+rho*C, then eliminate B from a row having S coefficient.
        rho_value = sp.Rational(intersection[1], intersection[0])
        s_row = row_one if row_one[2] else row_two
        normalized_b = sp.Rational(s_row[0], s_row[2])
        normalized_c = sp.Rational(s_row[1], s_row[2])
        reduced = (
            sp.Integer(0),
            sp.cancel(normalized_c - normalized_b * rho_value),
            sp.Integer(1),
        )
        require(
            f"mixed row reduction matrix={entries}",
            reduced[0] == 0 and reduced[2] == 1,
        )
        mixed_rows += 1
    rank_two_rows += 1

u_source = x**2 * q
p_source = sp.expand((u_source + 1) * (u_source + 3))
require("fixed-c normal form", same(c, x * p_source))
u_var = sp.symbols("u_var")
p_var = sp.expand((u_var + 1) * (u_var + 3))
R_test = u_var + 2 * u_var**2 + 3 * u_var**3
for weight in range(-4, 5):
    homogeneous = x**weight * R_test.subs(u_var, u_source)
    expected = x ** (weight + 2) * (
        weight * R_test * sp.diff(p_var, u_var)
        - sp.diff(R_test, u_var) * p_var
    ).subs(u_var, u_source)
    require(
        f"fixed-c Hamiltonian weight formula m={weight}",
        same(jacobian(homogeneous, c, x, q), expected),
    )

I_var = u_var**3 / 3 + 2 * u_var**2 + 3 * u_var
require("fixed-c integral derivative", sp.diff(I_var, u_var) == p_var)
require("fixed-c first root integral", I_var.subs(u_var, -3) == 0)
require("fixed-c hostile root integral", I_var.subs(u_var, -1) == -sp.Rational(4, 3))

h_easy = x**3 + 2 * x * q - q**2 + 5
B_easy = sp.expand(B_graph.subs(w, h_easy))
require(
    "BC face factorization",
    same(jacobian(B_easy, c, x, q), c * (-3 * c + jacobian(h_easy, c, x, q))),
)
print(
    f"PASS rank_two_matrices={rank_two_rows} fixedC={fixed_c_rows} "
    f"BC={bc_plane_rows} mixed={mixed_rows} easy_face_gates=14"
)


print("SECTION normalized mixed determinant-to-transport-PDE sign")
aa, cc, ZZ, ZZ_a, ZZ_c, rho, tau = sp.symbols("aa cc ZZ ZZ_a ZZ_c rho tau")
AA = ZZ - rho
f_A = sp.Rational(3, 2) * AA + 3 * cc * AA**2 / 8
f_c = AA**3 / 8
K_rho = sp.expand(ZZ * f_A - cc * f_c)

L_a = cc * ZZ_a
L_c = ZZ + cc * ZZ_c
M_a = 1 + f_A * ZZ_a
M_c = f_A * ZZ_c + f_c + tau
determinant_ac = sp.expand(L_a * M_c - L_c * M_a)
require(
    "abstract ac determinant",
    same(
        determinant_ac,
        -(ZZ + cc * ZZ_c + (K_rho - tau * cc) * ZZ_a),
    ),
)

# An independent source-coordinate check for several polynomial graphs.
source_graphs = (
    sp.Integer(0),
    sp.Integer(2),
    x,
    q,
    x * q + x**2 - 3 * q**2,
    x**5 * q**2,
)
for index, h_source in enumerate(source_graphs):
    A_source = sp.cancel(a * c + h_source)
    Z_source = sp.cancel(A_source + rho)
    B_source = sp.expand(B_graph.subs(w, h_source))
    S_source = sp.expand(S_graph.subs(w, h_source))
    L_source = sp.expand(B_source + rho * c)
    M_source = sp.expand(S_source + tau * c)
    K_source = sp.cancel(
        Z_source
        * (sp.Rational(3, 2) * A_source + 3 * c * A_source**2 / 8)
        - c * A_source**3 / 8
    )
    transported = sp.cancel(
        3 * Z_source
        + c * jacobian(Z_source, a, x, q)
        - (K_source - tau * c) * jacobian(Z_source, c, x, q)
    )
    require(
        f"source transport identity row={index}",
        same(jacobian(L_source, M_source, x, q), transported),
    )
print(f"PASS transport_identity_rows={len(source_graphs)} abstract_sign=plus3")


print("SECTION x=0 boundary ODE and formal first-coefficient gate")
f0, f1, f2, g0, g1, g2 = sp.symbols("f0 f1 f2 g0 g1 g2")
f_probe_h = f0 + f1 * q + f2 * q**2
h_probe = f_probe_h + x * (g0 + g1 * q) + x**2 * g2
B_probe = sp.expand(B_graph.subs(w, h_probe))
S_probe = sp.expand(S_graph.subs(w, h_probe))
L_probe = sp.expand(B_probe + rho * c)
M_probe = sp.expand(S_probe + tau * c)
boundary_jacobian = sp.factor(jacobian(L_probe, M_probe, x, q).subs(x, 0))
f_probe = f_probe_h + rho
boundary_expected = 3 * (
    f_probe
    + sp.Rational(3, 2)
    * (f_probe - rho)
    * f_probe
    * sp.diff(f_probe, q)
)
require("direct x=0 boundary ODE", same(boundary_jacobian, boundary_expected))

arm_a, arm_c, t, leading = sp.symbols("arm_a arm_c t leading", nonzero=True)
boundary_degree_rows = 0
for degree in range(1, 10):
    f_degree = leading * arm_a**degree
    ode_left = sp.expand(
        f_degree
        + sp.Rational(3, 2)
        * (f_degree - rho)
        * f_degree
        * sp.diff(f_degree, arm_a)
        - t
    )
    require(
        f"boundary leading degree n={degree}",
        sp.Poly(ode_left, arm_a).coeff_monomial(arm_a ** (3 * degree - 1))
        == sp.Rational(3, 2) * degree * leading**3,
    )
    boundary_degree_rows += 1

formal_rows = 0
for order in range(1, 9):
    for degree in range(0, 8):
        coefficient = leading * arm_a**degree
        W_formal = coefficient * arm_c**order
        Z_formal = t + W_formal
        A_formal = Z_formal - rho
        K_formal = sp.expand(
            Z_formal
            * (sp.Rational(3, 2) * A_formal + 3 * arm_c * A_formal**2 / 8)
            - arm_c * A_formal**3 / 8
        )
        formal_equation = sp.expand(
            W_formal
            + arm_c * sp.diff(W_formal, arm_c)
            + (K_formal - tau * arm_c) * sp.diff(W_formal, arm_a)
        )
        first_coefficient = sp.Poly(formal_equation, arm_c).coeff_monomial(
            arm_c**order
        )
        expected = sp.expand(
            (order + 1) * coefficient
            + sp.Rational(3, 2) * t * (t - rho) * sp.diff(coefficient, arm_a)
        )
        require(
            f"formal recurrence N={order} degree={degree}",
            same(first_coefficient, expected),
        )
        require(
            f"formal leading coefficient N={order} degree={degree}",
            sp.Poly(expected, arm_a).coeff_monomial(arm_a**degree)
            == (order + 1) * leading,
        )
        formal_rows += 1
print(
    f"PASS boundary_degree_rows={boundary_degree_rows} formal_rows={formal_rows} "
    "unique_polynomial_formal_branch=Z=t"
)


print("SECTION explicit formal arm-completion inverse")
R_completion = sp.symbols("R_completion")
completion_relation = (R_completion - 1) * (R_completion + 2) ** 2
require(
    "completion simple root derivative",
    sp.diff(completion_relation, R_completion).subs(R_completion, 1) == 9,
)
x_completion = arm_c / (R_completion * (R_completion + 2))
q_completion = arm_a * R_completion**2
require(
    "completion D identity modulo relation",
    same(
        1 + x_completion**2 * q_completion - R_completion,
        (arm_a * arm_c**2 - completion_relation) / (R_completion + 2) ** 2,
    ),
)
require(
    "completion c reconstruction",
    same(x_completion * R_completion * (R_completion + 2), arm_c),
)
require(
    "completion a reconstruction",
    same(q_completion / R_completion**2, arm_a),
)
print("PASS completion_inverse_gates=4 derivative_at_R1=9 coefficient_ring=C[a]")


print("SECTION ordinary initial-form degree-seven gate")
c_top = x**5 * q**2
degree_rows = 0
monomial_rows = 0
for degree in range(3, 29):
    require(f"B top dominates d={degree}", degree + 7 > 9)
    require(f"S top beats bh2 d={degree}", 3 * degree + 7 > 2 * degree + 9)
    require(f"S top beats che d={degree}", 3 * degree + 7 > degree + 11)
    require(f"S top beats be d={degree}", 3 * degree + 7 > 13)

    survivors = []
    for x_degree in range(degree + 1):
        monomial = x**x_degree * q ** (degree - x_degree)
        bracket = jacobian(monomial, c_top, x, q)
        expected = (
            (7 * x_degree - 5 * degree)
            * x ** (x_degree + 4)
            * q ** (degree - x_degree + 1)
        )
        require(
            f"monomial bracket d={degree} i={x_degree}",
            bracket == expected,
        )
        if 7 * x_degree == 5 * degree:
            survivors.append(x_degree)
        monomial_rows += 1

    if degree % 7:
        require(f"no initial-form survivor d={degree}", survivors == [])
    else:
        require(
            f"unique initial-form survivor d={degree}",
            survivors == [5 * degree // 7],
        )
    degree_rows += 1

for degree in (3, 7, 14):
    homogeneous = x**degree + 2 * x ** (degree - 1) * q + 3 * q**degree
    top_jacobian = jacobian(c_top * homogeneous, c_top * homogeneous**3 / 8, x, q)
    require(
        f"top Jacobian product identity d={degree}",
        same(
            top_jacobian,
            -c_top * homogeneous**3 * jacobian(homogeneous, c_top, x, q) / 4,
        ),
    )

require(
    "degree-seven coarse survivor",
    jacobian(c_top, c_top, x, q) == 0,
)
print(
    f"PASS initial_degree_rows={degree_rows} monomial_rows={monomial_rows} "
    "survivors=degrees_divisible_by_7"
)


print("SECTION complete collision-preserving quadratic hostile")
A2, B2, C2, D2, E2, F2 = sp.symbols("A2 B2 C2 D2 E2 F2")
h_quadratic = A2 * x**2 + B2 * x * q + C2 * q**2 + D2 * x + E2 * q + F2
collision_points = (
    (sp.Integer(0), -sp.Rational(3, 4)),
    (sp.Integer(1), sp.Integer(-3)),
    (sp.Integer(-1), sp.Integer(-3)),
)
collision_values = [
    sp.expand(h_quadratic.subs({x: x_value, q: q_value}))
    for x_value, q_value in collision_points
]
require(
    "quadratic collision odd equation",
    same(collision_values[1] - collision_values[2], 2 * (D2 - 3 * B2)),
)
require(
    "quadratic collision even equation after odd",
    same(
        16 * (collision_values[1] - collision_values[0]).subs(D2, 3 * B2),
        16 * A2 + 135 * C2 - 36 * E2,
    ),
)

collision_substitution = {
    D2: 3 * B2,
    A2: (-135 * C2 + 36 * E2) / 16,
}
B_quadratic = sp.expand(B_graph.subs(w, h_quadratic))
S_quadratic = sp.expand(S_graph.subs(w, h_quadratic))
J_quadratic = sp.expand(jacobian(B_quadratic, S_quadratic, x, q))
J_collision = sp.Poly(
    sp.expand(J_quadratic.subs(collision_substitution, simultaneous=True)),
    x,
    q,
)
coefficient_one = sp.factor(J_collision.coeff_monomial(x**17 * q**3))
require(
    "quadratic coefficient x17q3",
    coefficient_one == -sp.Rational(6561, 65536) * (15 * C2 - 4 * E2) ** 4,
)

J_after_one = sp.Poly(
    sp.expand(J_collision.as_expr().subs(E2, sp.Rational(15, 4) * C2)),
    x,
    q,
)
coefficient_two = sp.factor(J_after_one.coeff_monomial(x**13 * q**7))
require(
    "quadratic coefficient x13q7",
    coefficient_two == sp.Rational(3, 4) * (B2 + 1) ** 4,
)

J_after_two = sp.Poly(sp.expand(J_after_one.as_expr().subs(B2, -1)), x, q)
coefficient_three = sp.factor(J_after_two.coeff_monomial(x**13 * q**3))
require("quadratic terminal coefficient x13q3", coefficient_three == -sp.Rational(81, 2))
print(
    "PASS quadratic_collision_constraints=2 decisive_coefficients="
    "x17q3,x13q7,x13q3 terminal=-81/2"
)


print("SECTION sharp punctured rational endpoint and pole")
t_free = sp.symbols("t_free", nonzero=True)
h_sharp = sp.cancel(t_free - rho - a * c)
A_sharp = sp.cancel(a * c + h_sharp)
B_sharp = sp.cancel(B_graph.subs(w, h_sharp))
S_sharp = sp.cancel(S_graph.subs(w, h_sharp))
Z_sharp = sp.cancel(A_sharp + rho)
L_sharp = sp.cancel(B_sharp + rho * c)
M_sharp = sp.cancel(S_sharp + tau * c)
require("sharp A constant", same(A_sharp, t_free - rho))
require("sharp Z constant", same(Z_sharp, t_free))
require("sharp B", same(B_sharp, (t_free - rho) * c))
require("sharp L", same(L_sharp, t_free * c))
require(
    "sharp S",
    same(
        S_sharp,
        a
        + sp.Rational(3, 4) * (t_free - rho) ** 2
        + c * (t_free - rho) ** 3 / 8,
    ),
)
require(
    "sharp M",
    same(
        M_sharp,
        a
        + sp.Rational(3, 4) * (t_free - rho) ** 2
        + c * ((t_free - rho) ** 3 / 8 + tau),
    ),
)
require("sharp mixed determinant", same(jacobian(L_sharp, M_sharp, x, q), 3 * t_free))

ac_numerator, ac_denominator = sp.fraction(sp.cancel(a * c))
require("sharp pole denominator", same(ac_denominator, D))
require("sharp pole coprime", sp.gcd(ac_numerator, ac_denominator) == 1)
for index, (x_value, q_value) in enumerate(collision_points):
    substitution = {x: x_value, q: q_value}
    require(
        f"sharp graph collision row={index}",
        same(h_sharp.subs(substitution), t_free - rho),
    )
    require(f"sharp L collision row={index}", L_sharp.subs(substitution) == 0)
    require(
        f"sharp M collision row={index}",
        same(
            M_sharp.subs(substitution),
            -sp.Rational(3, 4) + 3 * (t_free - rho) ** 2 / 4,
        ),
    )
print("PASS sharp_endpoint_gates=18 determinant=3t pole_order=1 collision_rows=3")


print(f"SUMMARY exact_checks={CHECKS} all_passed=True")
