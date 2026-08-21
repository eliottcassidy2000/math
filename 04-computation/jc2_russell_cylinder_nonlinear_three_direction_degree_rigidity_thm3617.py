#!/usr/bin/env python3
"""Exact controls for provisional THM-3617 nonlinear target rigidity.

The universal no-go is proof-driven through injectivity on one line and
Jung--van der Kulk degree divisibility.  This companion verifies the Russell
compiler, the line inverse, the quotient-safe three-direction certificate,
the postcomposition chain rule, all target-monomial degree invoices, the three
graph degree profiles, the exceptional Jacobian trace, and sharp method
hostiles.  It uses exact SymPy arithmetic and no assertion gates.
"""

from itertools import product

import sympy as sp


CHECKS = 0


def require(label, condition):
    """Record one exact deterministic gate."""
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


def homogeneous_component(expression, degree, x_var, q_var, scale_var):
    """Extract one ordinary homogeneous component by exact scaling."""
    scaled = sp.expand(
        expression.subs(
            {x_var: scale_var * x_var, q_var: scale_var * q_var},
            simultaneous=True,
        )
    )
    return sp.expand(scaled.coeff(scale_var, degree))


def target_exponents(maximum_degree):
    """All exponent pairs of target total degree at most the cap."""
    return tuple(
        (first, second)
        for first in range(maximum_degree + 1)
        for second in range(maximum_degree + 1 - first)
    )


print("THM-3617 exact companion -- provisional nonlinear three-direction rigidity")
print("status=verified exact controls; independent hostile audit pending; cited all-degree steps are proof-driven")


print("SECTION Russell compiler and graph polynomiality")
x, q, w = sp.symbols("x q w")
D = 1 + x**2 * q
a = q / D**2
c = x * D * (D + 2)
b = (D - 1) * (D + 2) ** 2
e = q * (D + 3)

B = sp.expand(b + c * w)
Ccoord = c
Y = sp.expand(c * e + (2 * b + 4) * w + c * w**2)
S = sp.expand(((b + 2) * (e + 3 * w**2) + c * w * (3 * e + w**2)) / 8)
Aarm = sp.cancel(a * c + w)

require("compiler b=ac^2", same(b, a * c**2))
require("compiler e=a(b+4)", same(e, a * (b + 4)))
require("compiler arm determinant", same(jacobian(a, c, x, q), -3))
require("Russell B=cA", same(B, c * Aarm))
require("Russell C=c", same(Ccoord, c))
require("Russell Y arm", same(Y, 4 * Aarm + c * Aarm**2))
require("Russell S arm", same(S, a + 3 * Aarm**2 / 4 + c * Aarm**3 / 8))
for label, expression in (("B", B), ("C", Ccoord), ("Y", Y), ("S", S)):
    require(f"{label} graph polynomial", sp.denom(sp.cancel(expression)) == 1)
print("PASS compiler_gates=11 determinant=-3")


print("SECTION arbitrary target coefficient universe and x=0 reconstruction")
T, U = sp.symbols("T U")
theta, lam, eta, nu = sp.symbols("theta lambda eta nu", nonzero=True)
r00, r10, r01, r20, r11, r02, r21, r12, r03 = sp.symbols(
    "r00 r10 r01 r20 r11 r02 r21 r12 r03"
)
v10, v01, v11, v02 = sp.symbols("v10 v01 v11 v02")
q00, q10, q01, q20, q11, q02, q30, q21, q12, q03 = sp.symbols(
    "q00 q10 q01 q20 q11 q02 q30 q21 q12 q03"
)

Rpoly = (
    r00
    + r10 * T
    + r01 * U
    + r20 * T**2
    + r11 * T * U
    + r02 * U**2
    + lam * T**3
    + r21 * T**2 * U
    + r12 * T * U**2
    + r03 * U**3
)
Vpoly = nu + v10 * T + v01 * U + eta * T**2 + v11 * T * U + v02 * U**2
Qpoly = (
    q00
    + q10 * T
    + q01 * U
    + q20 * T**2
    + q11 * T * U
    + q02 * U**2
    + q30 * T**3
    + q21 * T**2 * U
    + q12 * T * U**2
    + q03 * U**3
)

require("R target degree", sp.Poly(Rpoly, T, U).total_degree() == 3)
require("R B3 coefficient", sp.Poly(Rpoly, T, U).coeff_monomial(T**3) == lam)
require("V target degree", sp.Poly(Vpoly, T, U).total_degree() == 2)
require("V B2 coefficient", sp.Poly(Vpoly, T, U).coeff_monomial(T**2) == eta)
require("V origin value", Vpoly.subs({T: 0, U: 0}) == nu)
require("Q target degree cap", sp.Poly(Qpoly, T, U).total_degree() == 3)

phi = sp.symbols("phi")
line_B = B.subs({x: 0, w: phi})
line_C = Ccoord.subs(x, 0)
line_Y = Y.subs({x: 0, w: phi})
line_S = S.subs({x: 0, w: phi})
line_F = sp.expand(theta * line_Y + Rpoly.subs({T: line_B, U: line_C}))
line_G = sp.expand(
    Vpoly.subs({T: line_B, U: line_C}) * line_S
    + Qpoly.subs({T: line_B, U: line_C})
)
recovered_phi = sp.cancel((line_F - r00) / (4 * theta))
recovered_q = sp.cancel((line_G - q00) / nu - 3 * recovered_phi**2 / 4)

require("line B=0", line_B == 0)
require("line C=0", line_C == 0)
require("line Y=4phi", line_Y == 4 * phi)
require("line S=q+3phi2/4", same(line_S, q + 3 * phi**2 / 4))
require("line F", same(line_F, 4 * theta * phi + r00))
require("line G", same(line_G, nu * (q + 3 * phi**2 / 4) + q00))
require("line reconstruct phi", same(recovered_phi, phi))
require("line reconstruct q", same(recovered_q, q))
print("PASS target_coefficient_gates=6 line_restriction_and_reconstruction=8")


print("SECTION quotient-safe three-direction certificate and postcomposition chain rule")
Ytarget, Starget = sp.symbols("Ytarget Starget")
f_restricted = sp.expand(theta * Ytarget + Rpoly.subs(U, 0))
g_restricted = sp.expand(Vpoly.subs(U, 0) * Starget + Qpoly.subs(U, 0))
grad_f = tuple(sp.diff(f_restricted, variable) for variable in (T, Ytarget, Starget))
grad_g = tuple(sp.diff(g_restricted, variable) for variable in (T, Ytarget, Starget))

first_vector = tuple(sp.Poly(entry, T, Starget).coeff_monomial(T**2) for entry in grad_f)
second_vector = tuple(sp.Poly(entry, T, Starget).coeff_monomial(1) for entry in grad_f)
third_vector = tuple(sp.Poly(entry, T, Starget).coeff_monomial(1) for entry in grad_g)
direction_matrix = sp.Matrix((first_vector, second_vector, third_vector))
require("first gradient coefficient vector", first_vector == (3 * lam, 0, 0))
require("second gradient coefficient vector", second_vector == (r10, theta, 0))
require("third gradient coefficient vector", third_vector == (q10, 0, nu))
require("three-direction determinant", direction_matrix.det() == 3 * lam * theta * nu)

# The determinant above is an ancillary ambient control.  The load-bearing
# certificate works in the Russell quotient.  A polynomial H(U,V) killed by
# ell*d_U+m*d_V is a polynomial in mU-ell V; exact powers supply controls for
# that characteristic-zero kernel.  If the resulting S-free linear form has
# a Y coefficient, its highest pure Y power survives the quotient normal
# form with leading monomial C*Y.  If it has no Y coefficient, no Y appears.
ell, mm, UU, VV = sp.symbols("ell mm UU VV")
kernel_coordinate = mm * UU - ell * VV
directional_kernel_rows = 0
for degree in range(9):
    packet = sp.expand(kernel_coordinate**degree)
    require(
        f"directional kernel power {degree}",
        zero(ell * sp.diff(packet, UU) + mm * sp.diff(packet, VV)),
    )
    directional_kernel_rows += 1

Bt, Ct, Yt = sp.symbols("Bt Ct Yt")
n0, nB, nC, nY = sp.symbols("n0 nB nC nY")
quotient_basis = sp.groebner(
    [Ct * Yt - Bt * (Bt + 4)],
    Ct,
    Yt,
    Bt,
    order="lex",
    domain=sp.QQ[n0, nB, nC, nY],
)
linear_form = n0 + nB * Bt + nC * Ct + nY * Yt
pure_y_rows = 0
for degree in range(2, 9):
    normal = sp.expand(quotient_basis.reduce(sp.expand(linear_form**degree))[1])
    pure_coefficient = sp.Poly(normal, Ct, Yt, Bt).coeff_monomial(Yt**degree)
    require(f"quotient pure-Y survivor degree={degree}", pure_coefficient == nY**degree)
    pure_y_rows += 1

no_y_rows = 0
for degree in range(1, 9):
    no_y_packet = sp.expand((linear_form.subs(nY, 0)) ** degree)
    require(f"quotient no-Y branch degree={degree}", sp.diff(no_y_packet, Yt) == 0)
    no_y_rows += 1
require("linear one-form has no B3", sp.Poly(linear_form, Bt, Ct, Yt).coeff_monomial(Bt**3) == 0)

H1u, H1v, H2u, H2v = sp.symbols("H1u H1v H2u H2v")
Lx, Lq, Mx, Mq = sp.symbols("Lx Lq Mx Mq")
composed_jacobian = sp.expand(
    (H1u * Lx + H1v * Mx) * (H2u * Lq + H2v * Mq)
    - (H1u * Lq + H1v * Mq) * (H2u * Lx + H2v * Mx)
)
factored_jacobian = sp.expand((H1u * H2v - H1v * H2u) * (Lx * Mq - Lq * Mx))
require("postcomposition chain rule", composed_jacobian == factored_jacobian)
print(
    "PASS ambient_direction_gates=4 determinant=3*lambda*theta*nu "
    f"quotient_directional_kernel={directional_kernel_rows} pure_Y={pure_y_rows} "
    f"no_Y={no_y_rows} linear_B3=1 chain_rule_gates=1"
)


print("SECTION complete target-monomial degree invoices")
R_exponents = target_exponents(3)
V_exponents = target_exponents(2)
Q_exponents = target_exponents(3)
require("R exponent count", len(R_exponents) == 10)
require("V exponent count", len(V_exponents) == 6)
require("Q exponent count", len(Q_exponents) == 10)


def degree_invoice(label, b_degree, c_degree, y_degree, s_degree, expected_f, expected_g):
    """Check unique target leading monomials for one source degree profile."""
    r_weights = {(i, j): i * b_degree + j * c_degree for i, j in R_exponents}
    v_weights = {(i, j): i * b_degree + j * c_degree + s_degree for i, j in V_exponents}
    q_weights = {(i, j): i * b_degree + j * c_degree for i, j in Q_exponents}
    f_weights = {("Y", 0): y_degree}
    f_weights.update({("R", exponent): weight for exponent, weight in r_weights.items()})
    g_weights = {("VS", exponent): weight for exponent, weight in v_weights.items()}
    g_weights.update({("Q", exponent): weight for exponent, weight in q_weights.items()})
    f_max = max(f_weights.values())
    g_max = max(g_weights.values())
    f_winners = tuple(key for key, value in f_weights.items() if value == f_max)
    g_winners = tuple(key for key, value in g_weights.items() if value == g_max)
    require(f"{label} F maximum", f_max == expected_f)
    require(f"{label} G maximum", g_max == expected_g)
    require(f"{label} F unique B3", f_winners == (("R", (3, 0)),))
    require(f"{label} G unique B2S", g_winners == (("VS", (2, 0)),))
    for exponent, weight in r_weights.items():
        if exponent != (3, 0):
            require(f"{label} R lower exponent={exponent}", weight < expected_f)
    require(f"{label} Y lower", y_degree < expected_f)
    for exponent, weight in v_weights.items():
        if exponent != (2, 0):
            require(f"{label} VS lower exponent={exponent}", weight < expected_g)
    for exponent, weight in q_weights.items():
        require(f"{label} Q lower exponent={exponent}", weight < expected_g)


invoice_rows = 0
for degree in range(3, 101):
    degree_invoice(
        f"high d={degree}",
        degree + 7,
        7,
        2 * degree + 7,
        3 * degree + 7,
        3 * degree + 21,
        5 * degree + 21,
    )
    invoice_rows += 1
degree_invoice("low", 9, 7, 11, 13, 27, 31)
degree_invoice("generic quadratic", 9, 7, 11, 13, 27, 31)
degree_invoice("cancelled quadratic", 8, 7, 9, 10, 24, 26)
invoice_rows += 3
print(
    f"PASS target_exponents=R{len(R_exponents)},V{len(V_exponents)},Q{len(Q_exponents)} "
    f"degree_invoice_profiles={invoice_rows}"
)


print("SECTION exact source leading forms for high graph degrees")
scale = sp.symbols("scale")
c7 = x**5 * q**2
k = x * q
Fcore = sp.expand(theta * Y + lam * B**3)
Gcore = sp.expand((nu + eta * B**2) * S)
high_rows = 0
for degree in range(3, 8):
    h_top = x**degree + 2 * x ** (degree - 1) * q + 3 * q**degree
    h_source = h_top + 5 * x + 7 * q + 11
    B_source = sp.expand(B.subs(w, h_source))
    Y_source = sp.expand(Y.subs(w, h_source))
    S_source = sp.expand(S.subs(w, h_source))
    F_source = sp.expand(Fcore.subs(w, h_source))
    G_source = sp.expand(Gcore.subs(w, h_source))
    expected_B = sp.expand(c7 * h_top)
    expected_Y = sp.expand(c7 * h_top**2)
    expected_S = sp.expand(c7 * h_top**3 / 8)
    expected_F = sp.expand(lam * c7**3 * h_top**3)
    expected_G = sp.expand(eta * c7**3 * h_top**5 / 8)
    require(
        f"high B top d={degree}",
        homogeneous_component(B_source, degree + 7, x, q, scale) == expected_B,
    )
    require(
        f"high Y top d={degree}",
        homogeneous_component(Y_source, 2 * degree + 7, x, q, scale) == expected_Y,
    )
    require(
        f"high S top d={degree}",
        homogeneous_component(S_source, 3 * degree + 7, x, q, scale) == expected_S,
    )
    require(
        f"high F top d={degree}",
        homogeneous_component(F_source, 3 * degree + 21, x, q, scale) == expected_F,
    )
    require(
        f"high G top d={degree}",
        homogeneous_component(G_source, 5 * degree + 21, x, q, scale) == expected_G,
    )
    require(f"high F degree d={degree}", sp.Poly(F_source, x, q).total_degree() == 3 * degree + 21)
    require(f"high G degree d={degree}", sp.Poly(G_source, x, q).total_degree() == 5 * degree + 21)
    high_rows += 7
print(f"PASS high_source_rows={high_rows} tested_d=3..7")


print("SECTION low and quadratic source leading forms")
low_graphs = (sp.Integer(0), sp.Integer(5), 2 * x - 3 * q + 7)
low_rows = 0
for index, h_source in enumerate(low_graphs):
    F_source = sp.expand(Fcore.subs(w, h_source))
    G_source = sp.expand(Gcore.subs(w, h_source))
    require(
        f"low F top row={index}",
        homogeneous_component(F_source, 27, x, q, scale) == sp.expand(lam * c7**3 * k**3),
    )
    require(
        f"low G top row={index}",
        homogeneous_component(G_source, 31, x, q, scale) == sp.expand(eta * c7**3 * k**5 / 8),
    )
    require(f"low F degree row={index}", sp.Poly(F_source, x, q).total_degree() == 27)
    require(f"low G degree row={index}", sp.Poly(G_source, x, q).total_degree() == 31)
    low_rows += 4

A2, B2, C2, D1, E1, n = sp.symbols("A2 B2 C2 D1 E1 n")
h2 = A2 * x**2 + B2 * x * q + C2 * q**2
h1 = D1 * x + E1 * q
h_quadratic = h2 + h1 + n
B_quadratic = sp.expand(B.subs(w, h_quadratic))
S_quadratic = sp.expand(S.subs(w, h_quadratic))
F_quadratic = sp.expand(Fcore.subs(w, h_quadratic))
G_quadratic = sp.expand(Gcore.subs(w, h_quadratic))
require(
    "generic quadratic B9",
    homogeneous_component(B_quadratic, 9, x, q, scale) == sp.expand(c7 * (k + h2)),
)
require(
    "generic quadratic S13",
    homogeneous_component(S_quadratic, 13, x, q, scale) == sp.expand(c7 * (k + h2) ** 3 / 8),
)
require(
    "generic quadratic F27",
    homogeneous_component(F_quadratic, 27, x, q, scale)
    == sp.expand(lam * c7**3 * (k + h2) ** 3),
)
require(
    "generic quadratic G31",
    homogeneous_component(G_quadratic, 31, x, q, scale)
    == sp.expand(eta * c7**3 * (k + h2) ** 5 / 8),
)

h_cancel_linear = -k + h1 + n
B_cancel = sp.expand(B.subs(w, h_cancel_linear))
S_cancel = sp.expand(S.subs(w, h_cancel_linear))
F_cancel = sp.expand(Fcore.subs(w, h_cancel_linear))
G_cancel = sp.expand(Gcore.subs(w, h_cancel_linear))
require(
    "cancelled quadratic B8",
    homogeneous_component(B_cancel, 8, x, q, scale) == sp.expand(c7 * h1),
)
require(
    "cancelled quadratic S10",
    homogeneous_component(S_cancel, 10, x, q, scale) == sp.expand(c7 * h1**3 / 8),
)
require(
    "cancelled quadratic F24",
    homogeneous_component(F_cancel, 24, x, q, scale) == sp.expand(lam * c7**3 * h1**3),
)
require(
    "cancelled quadratic G26",
    homogeneous_component(G_cancel, 26, x, q, scale) == sp.expand(eta * c7**3 * h1**5 / 8),
)

quadratic_census = 0
quadratic_cancel = 0
for coefficients in product((-1, 0, 1), repeat=3):
    top_value = h2.subs({A2: coefficients[0], B2: coefficients[1], C2: coefficients[2]})
    is_cancelled = sp.expand(k + top_value) == 0
    require(
        f"quadratic cancellation tuple={coefficients}",
        is_cancelled == (coefficients == (0, -1, 0)),
    )
    quadratic_census += 1
    quadratic_cancel += is_cancelled

linear_census = 0
linear_zero = 0
for coefficients in product((-1, 0, 1), repeat=2):
    linear_value = h1.subs({D1: coefficients[0], E1: coefficients[1]})
    is_zero = sp.expand(linear_value) == 0
    require(f"linear zero tuple={coefficients}", is_zero == (coefficients == (0, 0)))
    linear_census += 1
    linear_zero += is_zero

require("quadratic census size", quadratic_census == 27)
require("quadratic unique cancellation", quadratic_cancel == 1)
require("linear census size", linear_census == 9)
require("linear unique zero", linear_zero == 1)
print(
    f"PASS low_rows={low_rows} quadratic_form_gates=8 "
    f"quadratic_census={quadratic_census} linear_census={linear_census}"
)


print("SECTION Jung--van der Kulk divisibility grid")
divisibility_rows = 0
for degree in range(3, 301):
    degree_F = 3 * degree + 21
    degree_G = 5 * degree + 21
    require(f"high ordered d={degree}", degree_F < degree_G)
    require(f"high difference d={degree}", degree_G - degree_F == 2 * degree)
    require(f"high less than double d={degree}", degree_G < 2 * degree_F)
    require(f"high nondivisibility d={degree}", degree_G % degree_F != 0)
    divisibility_rows += 4
for label, degree_F, degree_G in (
    ("low", 27, 31),
    ("generic quadratic", 27, 31),
    ("cancelled quadratic", 24, 26),
):
    require(f"{label} ordered", degree_F < degree_G)
    require(f"{label} less than double", degree_G < 2 * degree_F)
    require(f"{label} nondivisibility", degree_G % degree_F != 0)
    divisibility_rows += 3
print(f"PASS divisibility_rows={divisibility_rows} high_grid=3..300 boundary_profiles=3")


print("SECTION exceptional h=-xq+n arbitrary-coefficient trace")
h_exceptional = -x * q + n
B_exceptional = sp.expand(B.subs(w, h_exceptional))
C_exceptional = Ccoord
Y_exceptional = sp.expand(Y.subs(w, h_exceptional))
S_exceptional = sp.expand(S.subs(w, h_exceptional))

B_line = sp.expand(B_exceptional.subs(x, 0))
C_line = sp.expand(C_exceptional.subs(x, 0))
Y_line = sp.expand(Y_exceptional.subs(x, 0))
S_line = sp.expand(S_exceptional.subs(x, 0))
Bx_line = sp.expand(sp.diff(B_exceptional, x).subs(x, 0))
Cx_line = sp.expand(sp.diff(C_exceptional, x).subs(x, 0))
Yx_line = sp.expand(sp.diff(Y_exceptional, x).subs(x, 0))
Bq_line = sp.expand(sp.diff(B_exceptional, q).subs(x, 0))
Cq_line = sp.expand(sp.diff(C_exceptional, q).subs(x, 0))
Yq_line = sp.expand(sp.diff(Y_exceptional, q).subs(x, 0))
Sq_line = sp.expand(sp.diff(S_exceptional, q).subs(x, 0))

require("exceptional B line", B_line == 0)
require("exceptional C line", C_line == 0)
require("exceptional Y line", Y_line == 4 * n)
require("exceptional S line", S_line == q + 3 * n**2 / 4)
require("exceptional Bx line", Bx_line == 3 * n)
require("exceptional Cx line", Cx_line == 3)
require("exceptional Yx line", Yx_line == 8 * q + 3 * n**2)
require("exceptional Bq line", Bq_line == 0)
require("exceptional Cq line", Cq_line == 0)
require("exceptional Yq line", Yq_line == 0)
require("exceptional Sq line", Sq_line == 1)

R_T_origin = sp.diff(Rpoly, T).subs({T: 0, U: 0})
R_U_origin = sp.diff(Rpoly, U).subs({T: 0, U: 0})
V_origin = Vpoly.subs({T: 0, U: 0})
F_line_exceptional = theta * Y_line + Rpoly.subs({T: B_line, U: C_line})
Fx_line_exceptional = sp.expand(
    theta * Yx_line + R_T_origin * Bx_line + R_U_origin * Cx_line
)
Fq_line_exceptional = sp.expand(
    theta * Yq_line + R_T_origin * Bq_line + R_U_origin * Cq_line
)
Gq_line_exceptional = sp.expand(
    (sp.diff(Vpoly, T).subs({T: 0, U: 0}) * Bq_line
     + sp.diff(Vpoly, U).subs({T: 0, U: 0}) * Cq_line) * S_line
    + V_origin * Sq_line
    + sp.diff(Qpoly, T).subs({T: 0, U: 0}) * Bq_line
    + sp.diff(Qpoly, U).subs({T: 0, U: 0}) * Cq_line
)
J_line_exceptional = sp.expand(Fx_line_exceptional * Gq_line_exceptional)

expected_F_line = 4 * theta * n + r00
expected_Fx_line = theta * (8 * q + 3 * n**2) + 3 * n * r10 + 3 * r01
expected_J_line = sp.expand(nu * expected_Fx_line)
require("exceptional F line arbitrary target", same(F_line_exceptional, expected_F_line))
require("exceptional Fx line arbitrary target", same(Fx_line_exceptional, expected_Fx_line))
require("exceptional Fq line arbitrary target", Fq_line_exceptional == 0)
require("exceptional Gq line arbitrary target", Gq_line_exceptional == nu)
require("exceptional Jacobian line arbitrary target", J_line_exceptional == expected_J_line)
require("exceptional q coefficient", sp.Poly(J_line_exceptional, q).coeff_monomial(q) == 8 * theta * nu)
for symbol in (eta, v10, v01, v11, v02, q10, q01, q20, q11, q02, q30, q21, q12, q03):
    require(f"exceptional trace independent of {symbol}", symbol not in J_line_exceptional.free_symbols)
print("PASS exceptional_source_jets=11 target_trace_gates=20 q_coefficient=8*theta*nu")


print("SECTION sharp coefficient and degree-cap method hostiles")
h_degree_seven = x**7 + 2 * q**7 + x + 1
h_degree_three = x**3 + 2 * q**3 + x + 1
B7 = sp.expand(B.subs(w, h_degree_seven))
Y7 = sp.expand(Y.subs(w, h_degree_seven))
S7 = sp.expand(S.subs(w, h_degree_seven))
B3 = sp.expand(B.subs(w, h_degree_three))
Y3 = sp.expand(Y.subs(w, h_degree_three))
S3 = sp.expand(S.subs(w, h_degree_three))

lambda_zero_F = sp.expand(Y7 + B7**2)
lambda_zero_G = sp.expand((1 + B7**2) * S7)
require("lambda-zero hostile F degree", sp.Poly(lambda_zero_F, x, q).total_degree() == 28)
require("lambda-zero hostile G degree", sp.Poly(lambda_zero_G, x, q).total_degree() == 56)
require("lambda-zero hostile divisible profile", 56 == 2 * 28)

eta_zero_F = sp.expand(Y3 + B3**3)
eta_zero_G = sp.expand(S3 + B3**3)
require("eta-zero hostile F degree", sp.Poly(eta_zero_F, x, q).total_degree() == 30)
require("eta-zero hostile G degree", sp.Poly(eta_zero_G, x, q).total_degree() == 30)
require("eta-zero hostile equal profile", sp.Poly(eta_zero_F, x, q).total_degree() == sp.Poly(eta_zero_G, x, q).total_degree())

degree_four_R_F = sp.expand(Y7 + B7**4 + B7**3)
degree_four_R_G = sp.expand((1 + B7**2) * S7)
require("degree-four R hostile F degree", sp.Poly(degree_four_R_F, x, q).total_degree() == 56)
require("degree-four R hostile G degree", sp.Poly(degree_four_R_G, x, q).total_degree() == 56)
require("degree-four R hostile equal profile", sp.Poly(degree_four_R_F, x, q).total_degree() == sp.Poly(degree_four_R_G, x, q).total_degree())

degree_six_Q_F = sp.expand(Y3 + B3**3)
degree_six_Q_G = sp.expand((1 + B3**2) * S3 + B3**6)
require("degree-six Q hostile F degree", sp.Poly(degree_six_Q_F, x, q).total_degree() == 30)
require("degree-six Q hostile G degree", sp.Poly(degree_six_Q_G, x, q).total_degree() == 60)
require("degree-six Q hostile divisible profile", 60 == 2 * 30)

theta_zero_line_F = sp.expand(line_F.subs(theta, 0))
nu_zero_line_G = sp.expand(line_G.subs(nu, 0))
require("theta-zero first line coordinate constant", theta_zero_line_F == r00)
require("theta-zero trace coefficient", (8 * theta * nu).subs(theta, 0) == 0)
require("nu-zero second line coordinate constant", nu_zero_line_G == q00)
require("nu-zero trace coefficient", (8 * theta * nu).subs(nu, 0) == 0)
print("PASS hostile_gates=16 lambda_profile=28,56 eta_profile=30,30 R4_profile=56,56 Q6_profile=30,60")


print(f"PASS total_exact_gates={CHECKS}")
print("RESULT PASS -- provisional theorem candidate verified-exact; independent hostile audit pending")
