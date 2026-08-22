#!/usr/bin/env python3
"""Exact controls for proved and independently hostile-audited THM-3641.

The companion verifies the full ordinary-triple projective curvature-debt
atlas for decomposable Russell-cylinder quadratic-fold pairs.  It includes
the generic finite rho chart, both exceptional finite charts, the infinite
chart, the principal THM-3624 boundary, and exact polynomial controls.

The result is deliberately conditional: zero debt removes one retained J2
quotient and does not construct an actual J0 witness or a Keller pair.
"""

import ast
from pathlib import Path

import sympy as sp


CHECKS = 0


def require(label, condition):
    """Record one optimization-safe exact gate."""
    global CHECKS
    if condition != True:
        raise RuntimeError(f"FAIL {label}: {condition}")
    CHECKS += 1


def hadamard(left, right):
    """Coordinatewise product of equal-shaped column matrices."""
    return left.multiply_elementwise(right)


def verify_chart(label, h, p, q, m, s_curve, lambda_row, mu_row, expected_debt):
    """Verify moments and the exhaustive rank-one/rank-two jet identity."""
    ones = sp.ones(3, 1)
    h2 = hadamard(h, h)
    hp = hadamard(h, p)
    hq = hadamard(h, q)
    require(f"{label} lambda kills constants", sp.factor((lambda_row * ones)[0]) == 0)
    require(f"{label} lambda kills tangents", sp.factor((lambda_row * h)[0]) == 0)
    symplectic = q - hp
    require(
        f"{label} transformed symplectic constant",
        sp.factor(symplectic[0] - symplectic[1]) == 0
        and sp.factor(symplectic[1] - symplectic[2]) == 0,
    )
    require(
        f"{label} moment constant",
        sp.factor((mu_row * ones)[0] - (lambda_row * p)[0]) == 0,
    )
    require(
        f"{label} moment linear",
        sp.factor(
            (mu_row * h)[0]
            - (lambda_row * (3 * q - 2 * hp))[0]
        )
        == 0,
    )
    require(
        f"{label} moment quadratic",
        sp.factor((mu_row * h2)[0] - (lambda_row * hq)[0]) == 0,
    )
    require(
        f"{label} mixed moment",
        sp.factor(2 * (mu_row * h)[0] - (lambda_row * (q + hp))[0]) == 0,
    )

    rank_parameter = sp.symbols(f"rank_parameter_{label}")
    u_yy, u_yh, u_hh = sp.symbols(f"u_yy_{label} u_yh_{label} u_hh_{label}")
    a_y, a_h = sp.symbols(f"a_y_{label} a_h_{label}")
    b_y, b_h = sp.symbols(f"b_y_{label} b_h_{label}")

    retained_residue = (
        m
        + u_yy * p
        + u_yh * (q + hp)
        + u_hh * hq
        + b_y * p
        + b_h * (3 * q - 2 * hp)
        + rank_parameter * a_y * (2 * q - 3 * hp)
        - rank_parameter * a_h * hq
    )
    determinant_derivative = (
        s_curve
        + u_yy * ones
        + 2 * u_yh * h
        + u_hh * h2
        + b_y * ones
        + b_h * h
        - rank_parameter * a_y * h
        - rank_parameter * a_h * h2
    )
    direct_debt = sp.factor((lambda_row * m)[0] - (mu_row * s_curve)[0])
    require(
        f"{label} displayed debt",
        sp.factor(direct_debt - expected_debt) == 0,
    )
    require(
        f"{label} full all-cell identity",
        sp.factor(
            (lambda_row * retained_residue)[0]
            - (mu_row * determinant_derivative)[0]
            - expected_debt
        )
        == 0,
    )


print("THM-3641 exact companion -- proved curvature-debt atlas")
print("status=PROVED_VERIFIED_EXACT_INDEPENDENTLY_HOSTILE_AUDITED")


print("SECTION stable decomposable coefficients")
t = sp.symbols("t")
U1, A0, A1, C0, C1, L0 = sp.symbols("U1 A0 A1 C0 C1 L0")
V1, B0, B1, E0, E1, H0 = sp.symbols("V1 B0 B1 E0 E1 H0")
f_x = U1 + t * A1 + t**2 * C1
f_t = A0 + 2 * t * C0 + 3 * t**2 * L0
g_x = V1 + t * B1 + t**2 * E1
g_t = B0 + 2 * t * E0 + 3 * t**2 * H0
jacobian = sp.Poly(sp.expand(f_x * g_t - f_t * g_x), t)
require("stable J0", jacobian.nth(0) == U1 * B0 - A0 * V1)
require(
    "stable J1",
    jacobian.nth(1) == 2 * U1 * E0 + A1 * B0 - A0 * B1 - 2 * C0 * V1,
)
require(
    "stable J2",
    jacobian.nth(2)
    == 3 * U1 * H0 + 2 * A1 * E0 + C1 * B0 - A0 * E1 - 2 * C0 * B1 - 3 * L0 * V1,
)
print("PASS stable_coefficients=J0,J1,J2")


print("SECTION universal compiler jets")
x, q_symbol, xi = sp.symbols("x q_symbol xi")
D_general = 1 + x**2 * q_symbol
c_general = sp.expand(x * D_general * (D_general + 2))
z_general = sp.expand(q_symbol * (D_general + 3) + 3)
delta_c_general = sp.diff(c_general, q_symbol)
delta_z_general = sp.diff(z_general, q_symbol)

v_minus, u, v_plus = sp.symbols("v_minus u v_plus")
r_minus, r_zero, r_plus = sp.symbols("r_minus r_zero r_plus")
tau_c = sp.Matrix([2 * (v_minus + 6), 3, 2 * (6 - v_plus)])
tau_z = sp.Matrix([-2 * (v_minus + 9), 4 * u, 2 * (9 - v_plus)])
n_c = sp.Matrix([2, 0, -2])
n_z = sp.Matrix([-2, 4, -2])
m_c = sp.Matrix([-18 - 2 * v_minus, 0, -18 + 2 * v_plus])
m_z = sp.Matrix([12 + 2 * v_minus, 0, -12 + 2 * v_plus])
c_second = sp.Matrix(
    [
        2 * (r_minus - v_minus**2 - 18 * v_minus - 54),
        0,
        2 * (v_plus**2 - 18 * v_plus + 54 - r_plus),
    ]
)
z_second = sp.Matrix(
    [
        2 * (v_minus**2 + 12 * v_minus + 9 - r_minus),
        4 * r_zero + sp.Rational(9, 8),
        2 * (v_plus**2 - 12 * v_plus + 9 - r_plus),
    ]
)


def local_compiler_jet(point, q_value, slope, curvature):
    """Return tangent, vertical value/derivative, and curve curvature."""
    q_series = q_value + slope * xi + curvature * xi**2 / 2
    substitution = {x: point + xi, q_symbol: q_series}
    curve_c = sp.expand(c_general.subs(substitution))
    curve_z = sp.expand(z_general.subs(substitution))
    vertical_c = sp.expand(delta_c_general.subs(substitution))
    vertical_z = sp.expand(delta_z_general.subs(substitution))
    return (
        sp.diff(curve_c, xi).subs(xi, 0),
        sp.diff(curve_z, xi).subs(xi, 0),
        vertical_c.subs(xi, 0),
        vertical_z.subs(xi, 0),
        sp.diff(vertical_c, xi).subs(xi, 0),
        sp.diff(vertical_z, xi).subs(xi, 0),
        sp.diff(curve_c, xi, 2).subs(xi, 0),
        sp.diff(curve_z, xi, 2).subs(xi, 0),
    )


branch_packets = (
    (-1, -3, v_minus, r_minus),
    (0, -sp.Rational(3, 4), u, r_zero),
    (1, -3, v_plus, r_plus),
)
for index, packet in enumerate(branch_packets):
    actual = local_compiler_jet(*packet)
    expected = (
        tau_c[index],
        tau_z[index],
        n_c[index],
        n_z[index],
        m_c[index],
        m_z[index],
        c_second[index],
        z_second[index],
    )
    for jet_index, (left, right) in enumerate(zip(actual, expected)):
        require(
            f"compiler branch={index} jet={jet_index}",
            sp.factor(left - right) == 0,
        )

for index in range(3):
    require(
        f"compiler symplectic determinant branch={index}",
        sp.factor(tau_c[index] * n_z[index] - tau_z[index] * n_c[index]) == 12,
    )
print("PASS compiler_jets=tangent,vertical,vertical_prime,curvature")
print("PASS branchwise_symplectic_determinant=12")


print("SECTION generic finite rho chart")
rho = sp.symbols("rho")
M = 3 - u * rho
K = 9 - u * rho
generic_v_minus = -(18 + rho * (2 * u + 9)) / (4 + rho)
generic_v_plus = (18 + rho * (2 * u - 9)) / (4 - rho)
generic_substitution = {v_minus: generic_v_minus, v_plus: generic_v_plus}
slope_surface = 4 * u * (v_minus + v_plus) + 4 * v_minus * v_plus - 27 * v_minus + 27 * v_plus - 162
require(
    "generic finite slope surface",
    sp.factor(slope_surface.subs(generic_substitution)) == 0,
)
require(
    "generic finite common tangent",
    all(
        sp.factor(value) == 0
        for value in (
            4 * tau_c.subs(generic_substitution)
            - rho * tau_z.subs(generic_substitution)
            - 4 * M * sp.ones(3, 1)
        )
    ),
)

generic_h = sp.simplify(tau_z.subs(generic_substitution))
generic_p = sp.simplify((4 * n_c - rho * n_z) / (4 * M))
generic_q = n_z
generic_m = sp.simplify((4 * m_c.subs(generic_substitution) - rho * m_z.subs(generic_substitution)) / (4 * M))
generic_s = sp.simplify((4 * c_second.subs(generic_substitution) - rho * z_second.subs(generic_substitution)) / (4 * M))
expected_h = sp.Matrix([-4 * K / (4 + rho), 4 * u, 4 * K / (4 - rho)])
require(
    "generic finite normalized tangents",
    all(sp.factor(generic_h[i] - expected_h[i]) == 0 for i in range(3)),
)
require(
    "generic ordinary left-middle factor",
    sp.factor(generic_h[1] - generic_h[0] - 4 * (4 * u + 9) / (rho + 4)) == 0,
)
require(
    "generic ordinary middle-right factor",
    sp.factor(generic_h[2] - generic_h[1] - 4 * (9 - 4 * u) / (4 - rho)) == 0,
)
require(
    "generic ordinary endpoint factor",
    sp.factor(generic_h[2] - generic_h[0] - 32 * K / (16 - rho**2)) == 0,
)

generic_lambda = sp.Matrix(
    [[
        (rho + 4) * (9 - 4 * u) / (8 * K),
        -1,
        (4 - rho) * (9 + 4 * u) / (8 * K),
    ]]
)
generic_mu = sp.Matrix(
    [[
        (rho + 4) ** 2 * (9 - 4 * u) / (16 * M * K),
        rho / M,
        -(rho - 4) ** 2 * (9 + 4 * u) / (16 * M * K),
    ]]
)
generic_numerator = (
    (rho + 4) ** 3 * (9 - 4 * u) * r_minus
    - 32 * rho**2 * K * r_zero
    + (4 - rho) ** 3 * (9 + 4 * u) * r_plus
    + 27 * (576 - 320 * rho * u + 225 * rho**2 - 9 * rho**3 * u)
)
generic_debt = -generic_numerator / (32 * K * M**2)
verify_chart(
    "generic_finite",
    generic_h,
    generic_p,
    generic_q,
    generic_m,
    generic_s,
    generic_lambda,
    generic_mu,
    generic_debt,
)
require(
    "generic zero-debt curvature coefficient nonzero",
    sp.factor(
        sp.diff(generic_numerator, r_minus)
        - (rho + 4) ** 3 * (9 - 4 * u)
    )
    == 0,
)
print("PASS generic_finite_identity_and_zero_debt_hyperplane")
print("PASS generic_ordinary=K*(4u-9)*(4u+9)!=0")


print("SECTION exceptional finite charts")
free_v = sp.symbols("free_v")

rho4_substitution = {
    v_minus: -9,
    u: sp.Rational(9, 4),
    v_plus: free_v,
}
require("rho=4 slope surface", sp.factor(slope_surface.subs(rho4_substitution)) == 0)
require(
    "rho=4 common tangent",
    all(
        sp.factor(value) == 0
        for value in (
            4 * tau_c.subs(rho4_substitution)
            - 4 * tau_z.subs(rho4_substitution)
            + 24 * sp.ones(3, 1)
        )
    ),
)
rho4_h = tau_z.subs(rho4_substitution)
rho4_p = (n_z - n_c) / 6
rho4_q = n_z
rho4_m = (m_z.subs(rho4_substitution) - m_c.subs(rho4_substitution)) / 6
rho4_s = (z_second.subs(rho4_substitution) - c_second.subs(rho4_substitution)) / 6
require(
    "rho=4 ordinary factors",
    sp.factor(rho4_h[2] - rho4_h[0] - 2 * (9 - free_v)) == 0
    and sp.factor(rho4_h[2] - rho4_h[1] - (9 - 2 * free_v)) == 0,
)
rho4_lambda = sp.Matrix(
    [[
        (2 * free_v - 9) / (2 * (free_v - 9)),
        -1,
        -9 / (2 * (free_v - 9)),
    ]]
)
rho4_mu = sp.Matrix(
    [[-(2 * free_v - 9) / (3 * (free_v - 9)), -sp.Rational(2, 3), 0]]
)
rho4_debt = (
    32 * (free_v - 9) * r_zero
    - 16 * (2 * free_v - 9) * r_minus
    + 27 * (117 - 29 * free_v)
) / (72 * (free_v - 9))
verify_chart(
    "rho_plus_4",
    rho4_h,
    rho4_p,
    rho4_q,
    rho4_m,
    rho4_s,
    rho4_lambda,
    rho4_mu,
    rho4_debt,
)

rhom4_substitution = {
    v_minus: free_v,
    u: -sp.Rational(9, 4),
    v_plus: 9,
}
require("rho=-4 slope surface", sp.factor(slope_surface.subs(rhom4_substitution)) == 0)
require(
    "rho=-4 common tangent",
    all(
        sp.factor(value) == 0
        for value in (
            4 * tau_c.subs(rhom4_substitution)
            + 4 * tau_z.subs(rhom4_substitution)
            + 24 * sp.ones(3, 1)
        )
    ),
)
rhom4_h = tau_z.subs(rhom4_substitution)
rhom4_p = -(n_c + n_z) / 6
rhom4_q = n_z
rhom4_m = -(m_c.subs(rhom4_substitution) + m_z.subs(rhom4_substitution)) / 6
rhom4_s = -(c_second.subs(rhom4_substitution) + z_second.subs(rhom4_substitution)) / 6
require(
    "rho=-4 ordinary factors",
    sp.factor(rhom4_h[2] - rhom4_h[0] - 2 * (free_v + 9)) == 0
    and sp.factor(rhom4_h[1] - rhom4_h[0] - (2 * free_v + 9)) == 0,
)
rhom4_lambda = sp.Matrix(
    [[
        9 / (2 * (free_v + 9)),
        -1,
        (2 * free_v + 9) / (2 * (free_v + 9)),
    ]]
)
rhom4_mu = sp.Matrix(
    [[0, sp.Rational(2, 3), (2 * free_v + 9) / (3 * (free_v + 9))]]
)
rhom4_debt = (
    32 * (free_v + 9) * r_zero
    - 16 * (2 * free_v + 9) * r_plus
    - 27 * (29 * free_v + 117)
) / (72 * (free_v + 9))
verify_chart(
    "rho_minus_4",
    rhom4_h,
    rhom4_p,
    rhom4_q,
    rhom4_m,
    rhom4_s,
    rhom4_lambda,
    rhom4_mu,
    rhom4_debt,
)
print("PASS exceptional_rho=+4,-4 identities")
print("PASS exceptional_ordinary=vplus_notin(9,9/2),vminus_notin(-9,-9/2)")


print("SECTION infinite projective chart")
infinite_substitution = {v_minus: -9 - 2 * u, v_plus: 9 - 2 * u}
require(
    "infinite slope surface",
    sp.factor(slope_surface.subs(infinite_substitution)) == 0,
)
require(
    "infinite common tangent",
    all(
        sp.factor(value) == 0
        for value in (tau_z.subs(infinite_substitution) - 4 * u * sp.ones(3, 1))
    ),
)
infinite_h = tau_c.subs(infinite_substitution)
infinite_p = n_z / (4 * u)
infinite_q = n_c
infinite_m = m_z.subs(infinite_substitution) / (4 * u)
infinite_s = z_second.subs(infinite_substitution) / (4 * u)
require(
    "infinite ordinary factors",
    sp.factor(infinite_h[1] - infinite_h[0] - (4 * u + 9)) == 0
    and sp.factor(infinite_h[2] - infinite_h[1] - (4 * u - 9)) == 0
    and sp.factor(infinite_h[2] - infinite_h[0] - 8 * u) == 0,
)
infinite_lambda = sp.Matrix(
    [[(4 * u - 9) / (8 * u), -1, (4 * u + 9) / (8 * u)]]
)
infinite_mu = sp.Matrix(
    [[
        (9 - 4 * u) / (16 * u**2),
        -1 / u,
        -(4 * u + 9) / (16 * u**2),
    ]]
)
infinite_debt = (
    (9 - 4 * u) * r_minus
    + 32 * u * r_zero
    - (9 + 4 * u) * r_plus
    - 243 * u
) / (32 * u**3)
verify_chart(
    "rho_infinity",
    infinite_h,
    infinite_p,
    infinite_q,
    infinite_m,
    infinite_s,
    infinite_lambda,
    infinite_mu,
    infinite_debt,
)
print("PASS infinite_chart_identity")
print("PASS infinite_ordinary=u_notin(0,+9/4,-9/4)")


print("SECTION principal recovery and polynomial controls")
principal_lambda = sp.Matrix(
    [[(9 - 4 * u) / 18, -1, (9 + 4 * u) / 18]]
)
principal_mu = sp.Matrix(
    [[(9 - 4 * u) / 27, 0, -(9 + 4 * u) / 27]]
)
principal_debt = -sp.Rational(2, 81) * (
    (9 - 4 * u) * r_minus + (9 + 4 * u) * r_plus + 243
)
require(
    "principal lambda specialization",
    all(
        sp.factor(generic_lambda.subs(rho, 0)[i] - principal_lambda[i]) == 0
        for i in range(3)
    ),
)
require(
    "principal mu specialization",
    all(
        sp.factor(generic_mu.subs(rho, 0)[i] - principal_mu[i]) == 0
        for i in range(3)
    ),
)
require(
    "principal debt specialization",
    sp.factor(generic_debt.subs(rho, 0) - principal_debt) == 0,
)
print("PASS principal_boundary=THM3624_invoice")

Q1 = sp.Poly(
    x**5
    + sp.Rational(9, 2) * x**4
    - 2 * x**3
    - sp.Rational(27, 4) * x**2
    + x
    - sp.Rational(3, 4),
    x,
    domain=sp.QQ,
)
Qh = sp.Poly(
    (
        44069 * x**11
        + 112059 * x**10
        - 154749 * x**9
        - 406377 * x**8
        + 188107 * x**7
        + 513081 * x**6
        - 82835 * x**5
        - 230931 * x**4
        + 5408 * x
        - 4056
    )
    / 5408,
    x,
    domain=sp.QQ,
)
Q6 = sp.Poly(
    Q1.as_expr() - sp.Rational(259, 36) * x**2 * (x**2 - 1) ** 2,
    x,
    domain=sp.QQ,
)
Q_even = sp.Poly(
    -sp.Rational(27, 4) * x**6
    + 18 * x**4
    - sp.Rational(27, 2) * x**2
    - sp.Rational(3, 4),
    x,
    domain=sp.QQ,
)
Qstar = sp.Poly(Q_even.as_expr() + x * (1 - x**2) ** 3, x, domain=sp.QQ)
points = (-1, 0, 1)


def polynomial_packet(poly):
    """Return retained values, slopes, and curvatures."""
    return tuple(
        tuple(sp.diff(poly.as_expr(), x, order).subs(x, point) for point in points)
        for order in range(3)
    )


expected_packets = {
    "Q1": (
        (-3, -sp.Rational(3, 4), -3),
        (-sp.Rational(9, 2), 1, sp.Rational(9, 2)),
        (sp.Rational(65, 2), -sp.Rational(27, 2), sp.Rational(97, 2)),
    ),
    "Qh": (
        (-3, -sp.Rational(3, 4), -3),
        (-sp.Rational(9, 2), 1, sp.Rational(9, 2)),
        (0, 0, -sp.Rational(243, 13)),
    ),
    "Q6": (
        (-3, -sp.Rational(3, 4), -3),
        (-sp.Rational(9, 2), 1, sp.Rational(9, 2)),
        (-sp.Rational(451, 18), -sp.Rational(251, 9), -sp.Rational(163, 18)),
    ),
    "Qstar": (
        (-3, -sp.Rational(3, 4), -3),
        (-sp.Rational(9, 2), 1, sp.Rational(9, 2)),
        (-sp.Rational(27, 2), -27, -sp.Rational(27, 2)),
    ),
}
for name, poly in (("Q1", Q1), ("Qh", Qh), ("Q6", Q6), ("Qstar", Qstar)):
    require(f"{name} retained packet", polynomial_packet(poly) == expected_packets[name])
    packet = expected_packets[name]
    debt_value = principal_debt.subs(
        {u: packet[1][1], r_minus: packet[2][0], r_plus: packet[2][2]}
    )
    expected_value = -sp.Rational(2072, 81) if name == "Q1" else 0
    require(f"{name} principal debt", sp.factor(debt_value - expected_value) == 0)

require("Q6 exact perturbation", Q6 == Q1 - sp.Poly(sp.Rational(259, 36) * x**2 * (x**2 - 1) ** 2, x))
require("Q6 degree", Q6.degree() == 6)
require("Qstar exact odd perturbation", Qstar == Q_even + sp.Poly(x * (1 - x**2) ** 3, x))
require("Qstar degree", Qstar.degree() == 7)

control_tau_c = tau_c.subs({v_minus: -sp.Rational(9, 2), v_plus: sp.Rational(9, 2)})
control_tau_z = tau_z.subs({v_minus: -sp.Rational(9, 2), u: 1, v_plus: sp.Rational(9, 2)})
control_determinants = tuple(
    sp.det(
        sp.Matrix(
            [
                [control_tau_c[i], control_tau_z[i]],
                [control_tau_c[j], control_tau_z[j]],
            ]
        )
    )
    for i, j in ((0, 1), (1, 2), (0, 2))
)
require("control ordinary triple", control_determinants == (39, 15, 54))
print("PASS controls=Q1(-2072/81),Qh(0),Q6(0),Qstar(0)")
print("PASS control_tangent_determinants=(39,15,54)")


print("SECTION Hermite realizability and minimality")
hermite_rows = []
for derivative_order in range(3):
    for point in points:
        hermite_rows.append(
            [
                sp.diff(x**degree, x, derivative_order).subs(x, point)
                for degree in range(9)
            ]
        )
hermite_matrix = sp.Matrix(hermite_rows)
require("degree-eight Hermite rank", hermite_matrix.rank() == 9)
require("degree-eight Hermite determinant", hermite_matrix.det() != 0)

value_slope_rows = []
for derivative_order in range(2):
    for point in points:
        value_slope_rows.append(
            [
                sp.diff(x**degree, x, derivative_order).subs(x, point)
                for degree in range(6)
            ]
        )
value_slope_matrix = sp.Matrix(value_slope_rows)
require("degree-five value-slope rank", value_slope_matrix.rank() == 6)
require("degree-five value-slope determinant", value_slope_matrix.det() != 0)
Q1_coefficient_column = sp.Matrix([Q1.nth(degree) for degree in range(6)])
Q1_value_slope_packet = sp.Matrix(
    list(expected_packets["Q1"][0]) + list(expected_packets["Q1"][1])
)
require(
    "degree-five unique packet is Q1",
    value_slope_matrix * Q1_coefficient_column == Q1_value_slope_packet,
)
print("PASS every_value_slope_curvature_packet_has_degree<=8_control")
print("PASS Q6_minimal_degree_in_fixed_Q1_slope_packet=6")


print("SECTION active-gate audit")
source_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(source_tree))
require("assert-free AST", assert_nodes == 0)
print("PASS ast_assert_nodes=0")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
