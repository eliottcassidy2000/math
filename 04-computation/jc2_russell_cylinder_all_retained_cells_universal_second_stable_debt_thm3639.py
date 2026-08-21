#!/usr/bin/env python3
"""Exact controls for proved and independently hostile-audited THM-3639.

The companion derives the quadratic-fold stable coefficients, reconstructs
the minimal ordinary-triple compiler jets, checks representative and
higher-lift invariance, proves the exhaustive rank-one/rank-two SL_2 normal
form, and verifies the division-free endpoint-curvature identity

    lambda(J2) = -2072/81 + (5/27) J0'(-1) - (13/27) J0'(1).

All arithmetic is exact.  Every gate stays active under ``python -O``.
"""

import ast
from pathlib import Path

import sympy as sp


CHECKS = 0


def require(label, condition):
    """Record one active exact gate."""
    global CHECKS
    if condition != True:
        raise RuntimeError(f"FAIL {label}: {condition}")
    CHECKS += 1


def hadamard(left, right):
    """Coordinatewise product of equal-shaped column matrices."""
    return left.multiply_elementwise(right)


def det2(first, second):
    """Determinant of two two-component columns."""
    return first[0] * second[1] - first[1] * second[0]


print("THM-3639 exact companion -- proved all-retained-cell debt")
print("status=PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED")


print("SECTION stable Jacobian coefficients")
t = sp.symbols("t")
U1, A0s, A1, C0, C1, D0 = sp.symbols("U1 A0s A1 C0 C1 D0")
V1, B0s, B1, E0, E1, H0 = sp.symbols("V1 B0s B1 E0 E1 H0")
f_x = U1 + t * A1 + t**2 * C1
f_t = A0s + 2 * t * C0 + 3 * t**2 * D0
g_x = V1 + t * B1 + t**2 * E1
g_t = B0s + 2 * t * E0 + 3 * t**2 * H0
jacobian = sp.Poly(sp.expand(f_x * g_t - f_t * g_x), t)
expected_j0 = U1 * B0s - A0s * V1
expected_j1 = 2 * U1 * E0 + A1 * B0s - A0s * B1 - 2 * C0 * V1
expected_j2 = (
    3 * U1 * H0
    + 2 * A1 * E0
    + C1 * B0s
    - A0s * E1
    - 2 * C0 * B1
    - 3 * D0 * V1
)
require("stable coefficient J0", jacobian.nth(0) == expected_j0)
require("stable coefficient J1", jacobian.nth(1) == expected_j1)
require("stable coefficient J2", jacobian.nth(2) == expected_j2)
print("PASS stable_coefficients=J0,J1,J2")


print("SECTION minimal ordinary-triple compiler")
x, q = sp.symbols("x q")
Q = sp.Poly(
    x**5
    + sp.Rational(9, 2) * x**4
    - 2 * x**3
    - sp.Rational(27, 4) * x**2
    + x
    - sp.Rational(3, 4),
    x,
    domain=sp.QQ,
)
D_general = 1 + x**2 * q
b_general = sp.expand((D_general - 1) * (D_general + 2) ** 2)
c_general = sp.expand(x * D_general * (D_general + 2))
e_general = sp.expand(q * (D_general + 3))
b = sp.Poly(b_general.subs(q, Q.as_expr()), x, domain=sp.QQ)
c = sp.Poly(c_general.subs(q, Q.as_expr()), x, domain=sp.QQ)
e = sp.Poly(e_general.subs(q, Q.as_expr()), x, domain=sp.QQ)
delta_c = sp.Poly(sp.diff(c_general, q).subs(q, Q.as_expr()), x, domain=sp.QQ)
delta_z = sp.Poly(sp.diff(e_general, q).subs(q, Q.as_expr()), x, domain=sp.QQ)
points = (-1, 0, 1)

require("compiler surface relation", c**2 * e == b * (b + 4))
require(
    "retained target triple",
    [(b.eval(point), c.eval(point), e.eval(point)) for point in points]
    == [(0, 0, -3)] * 3,
)
tau_c = sp.Matrix([c.diff().eval(point) for point in points])
tau_z = sp.Matrix([e.diff().eval(point) for point in points])
c_second = sp.Matrix([c.diff().diff().eval(point) for point in points])
normal_c = sp.Matrix([delta_c.eval(point) for point in points])
normal_z = sp.Matrix([delta_z.eval(point) for point in points])
normal_c_prime = sp.Matrix([delta_c.diff().eval(point) for point in points])
require("tangent c", tau_c == sp.Matrix([3, 3, 3]))
require("tangent z", tau_z == sp.Matrix([-9, 4, 9]))
require(
    "tangent pair determinants",
    tuple(
        sp.det(sp.Matrix([[tau_c[i], tau_z[i]], [tau_c[j], tau_z[j]]]))
        for i, j in ((0, 1), (1, 2), (0, 2))
    )
    == (39, 15, 54),
)
require(
    "retained c second derivative",
    c_second == sp.Matrix([sp.Rational(157, 2), 0, -sp.Rational(221, 2)]),
)
require("vertical c values", normal_c == sp.Matrix([2, 0, -2]))
require("vertical z values", normal_z == sp.Matrix([-2, 4, -2]))
require("vertical c derivatives", normal_c_prime == sp.Matrix([-9, 0, -9]))
lambda_row = sp.Matrix([[sp.Rational(5, 18), -1, sp.Rational(13, 18)]])
require("lambda kills tangent c", (lambda_row * tau_c)[0] == 0)
require("lambda kills tangent z", (lambda_row * tau_z)[0] == 0)
print("PASS triple_tangents=(39,15,54) lambda=(5,-18,13)/18")
print("PASS vertical_jets=nc(2,0,-2)_nz(-2,4,-2)_ncprime(-9,0,-9)")


print("SECTION representative and unrestricted-quotient invariance")
linear_evaluation = sp.Matrix([[tau_c[i], tau_z[i]] for i in range(3)])
quadratic_evaluation = sp.Matrix(
    [[tau_c[i] ** 2, tau_c[i] * tau_z[i], tau_z[i] ** 2] for i in range(3)]
)
require("three-line linear rank", linear_evaluation.rank() == 2)
require("three-line quadratic determinant", quadratic_evaluation.det() == 31590)

# A cubic target germ is the first possible graph-ideal initial form.  One
# vertical derivative leaves branch order two, so neither its retained value
# nor first branch derivative can alter the quotient.
xi, epsilon = sp.symbols("xi epsilon")
h30, h21, h12, h03 = sp.symbols("h30 h21 h12 h03")
yc, yz = sp.symbols("yc yz")
general_cubic = h30 * yc**3 + h21 * yc**2 * yz + h12 * yc * yz**2 + h03 * yz**3
for index in range(3):
    branch_cubic = general_cubic.subs(
        {
            yc: tau_c[index] * xi + normal_c[index] * epsilon,
            yz: tau_z[index] * xi + normal_z[index] * epsilon,
        }
    )
    vertical_derivative = sp.Poly(
        sp.diff(branch_cubic, epsilon).subs(epsilon, 0), xi
    )
    require(
        f"cubic vertical branch order index={index}",
        vertical_derivative.nth(0) == 0 and vertical_derivative.nth(1) == 0,
    )

# Verify the determinant form of the unrestricted residue branch by branch.
du, dv, da, db = sp.symbols("du dv da db")
nu_u, nu_v, nu_a, nu_b = sp.symbols("nu_u nu_v nu_a nu_b")
dnu_u, dnu_v = sp.symbols("dnu_u dnu_v")
a0, b0 = sp.symbols("a0 b0")
unrestricted_scalar = (
    3 * du * nu_b
    + 2 * da * nu_v
    + b0 * dnu_u
    - a0 * dnu_v
    - 2 * nu_u * db
    - 3 * nu_a * dv
)
vertical_det_derivative = (
    dnu_u * b0 - a0 * dnu_v + nu_u * db - da * nu_v
)
determinant_scalar = vertical_det_derivative + 3 * (
    det2((du, dv), (nu_a, nu_b)) - det2((nu_u, nu_v), (da, db))
)
require("unrestricted determinant formula", sp.expand(unrestricted_scalar - determinant_scalar) == 0)
print("PASS ker_gamma_order>=3 and delta_ker_gamma_branch_order>=2")
print("PASS unrestricted_R2_depends_only_on_zero_2jets_first_1jets")


print("SECTION exhaustive SL2 normal forms")
g11, g12, g21, g22 = sp.symbols("g11 g12 g21 g22")
gradient_matrix = sp.Matrix([[g11, g12], [g21, g22]])
r_det = sp.expand(gradient_matrix.det())
rank_two_change = sp.Matrix(
    [[g22 / r_det, -g12 / r_det], [-g21, g11]]
)
require("rank-two change determinant", sp.simplify(rank_two_change.det()) == 1)
require(
    "rank-two diagonalization",
    sp.simplify(rank_two_change * gradient_matrix)
    == sp.diag(1, r_det),
)
normal_a0, normal_b0, r_symbol = sp.symbols("normal_a0 normal_b0 r_symbol")
normal_j0_values = [
    3 * normal_b0 - r_symbol * tau_z[index] * normal_a0
    for index in range(3)
]
rank_two_stable_solution = sp.solve(
    [value - 1 for value in normal_j0_values],
    (normal_a0, normal_b0),
    dict=True,
)
require(
    "rank-two stable value normal form",
    rank_two_stable_solution
    == [{normal_a0: 0, normal_b0: sp.Rational(1, 3)}],
)

# For a rank-one factor G=w*alpha^T, the common J0 values are
# Delta*(3 alpha_c + tau_z alpha_z).  Their differences force alpha_z=0,
# because any one value equal to one forces Delta nonzero.
Delta, alpha_c, alpha_z = sp.symbols("Delta alpha_c alpha_z")
rank_one_values = [
    sp.expand(Delta * (3 * alpha_c + tau_z[index] * alpha_z))
    for index in range(3)
]
require(
    "rank-one endpoint difference",
    sp.expand(rank_one_values[0] - rank_one_values[2]) == -18 * Delta * alpha_z,
)
require(
    "rank-one left-middle difference",
    sp.expand(rank_one_values[0] - rank_one_values[1]) == -13 * Delta * alpha_z,
)
rank_one_gradient = sp.Matrix([[1, 0], [0, 0]])
rank_one_free_a = sp.symbols("rank_one_free_a")
rank_one_stable = sp.Matrix([rank_one_free_a, sp.Rational(1, 3)])
rank_one_shear = sp.Matrix([[1, -3 * rank_one_free_a], [0, 1]])
require("rank-one stabilizer determinant", rank_one_shear.det() == 1)
require(
    "rank-one stabilizer preserves gradient",
    rank_one_shear * rank_one_gradient == rank_one_gradient,
)
require(
    "rank-one stabilizer kills A0",
    rank_one_shear * rank_one_stable == sp.Matrix([0, sp.Rational(1, 3)]),
)
print("PASS rank2_M=diag(1,detG)*Ginverse lies_in_SL2")
print("PASS rank1_factor_forces_alpha_z=0 and_SL2_shear_kills_A0")
print("PASS exhaustive_normal_form=gradU(1,0)_gradV(0,r)_A0=0_B0=1/3")


print("SECTION division-free universal curvature identity")
r = sp.symbols("r")
u_cc, u_cz, u_zz = sp.symbols("u_cc u_cz u_zz")
a_c, a_z, b_c, b_z = sp.symbols("a_c a_z b_c b_z")
H_u = sp.Matrix([[u_cc, u_cz], [u_cz, u_zz]])
dU = tau_c
dV = r * tau_z
dA = a_c * tau_c + a_z * tau_z
dB = b_c * tau_c + b_z * tau_z
nu_U = normal_c
nu_V = r * normal_z
nu_A = a_c * normal_c + a_z * normal_z
nu_B = b_c * normal_c + b_z * normal_z
U_second = []
dot_nu_U = []
for index in range(3):
    tangent = sp.Matrix([tau_c[index], tau_z[index]])
    normal = sp.Matrix([normal_c[index], normal_z[index]])
    U_second.append(c_second[index] + (tangent.T * H_u * tangent)[0])
    dot_nu_U.append(normal_c_prime[index] + (tangent.T * H_u * normal)[0])
U_second = sp.Matrix(U_second)
dot_nu_U = sp.Matrix(dot_nu_U)

E = U_second / 3 + hadamard(dU, dB) - hadamard(dA, dV)
R2 = (
    3 * hadamard(dU, nu_B)
    + 2 * hadamard(dA, nu_V)
    + dot_nu_U / 3
    - 2 * hadamard(nu_U, dB)
    - 3 * hadamard(nu_A, dV)
)
debt = sp.factor((lambda_row * R2)[0])
expected_raw_debt = (
    162 * r * a_c
    + 216 * r * a_z
    - 24 * b_c
    - 162 * b_z
    - 8 * u_cc
    - 108 * u_cz
    - 72 * u_zz
    - 27
) / 9
require("raw retained debt", sp.expand(debt - expected_raw_debt) == 0)
curvature_identity = sp.factor(
    debt
    + sp.Rational(2072, 81)
    - sp.Rational(5, 27) * E[0]
    + sp.Rational(13, 27) * E[2]
)
require("division-free curvature identity", curvature_identity == 0)

# Higher stable restrictions have a common value and tangent derivative.
# Insert all such free data into J2 and verify that lambda removes it.
x0, y0, z0, w0 = sp.symbols("x0 y0 z0 w0")
x_c, x_z, y_c, y_z = sp.symbols("x_c x_z y_c y_z")
dX = x_c * tau_c + x_z * tau_z
dY = y_c * tau_c + y_z * tau_z
ones = sp.ones(3, 1)
C_values = nu_U + x0 * ones
E_values = nu_V + y0 * ones
D_values = nu_A + z0 * ones
H_values = nu_B + w0 * ones
C_prime = dot_nu_U + dX
full_J2_values = (
    3 * hadamard(dU, H_values)
    + 2 * hadamard(dA, E_values)
    + C_prime / 3
    - 2 * hadamard(C_values, dB)
    - 3 * hadamard(D_values, dV)
)
require(
    "all higher lifts annihilated",
    sp.expand((lambda_row * (full_J2_values - R2))[0]) == 0,
)

# Rank-two and rank-one eliminations are independent controls on the unified
# identity rather than prerequisites for it.
rank_two_solution = sp.solve(list(E), (a_c, a_z, b_c), dict=True)
expected_rank_two_solution = {
    a_c: (6 * b_z + 4 * u_cz - 7) / (6 * r),
    a_z: (65 * u_zz - 58) / (195 * r),
    b_c: -u_cc / 3 - sp.Rational(3658, 1755),
}
require("rank-two affine solution count", len(rank_two_solution) == 1)
for variable in (a_c, a_z, b_c):
    require(
        f"rank-two solution {variable}",
        sp.factor(rank_two_solution[0][variable] - expected_rank_two_solution[variable])
        == 0,
    )
require(
    "rank-two debt",
    sp.factor(debt.subs(rank_two_solution[0])) == -sp.Rational(2072, 81),
)

rank_one_solution = sp.solve(
    list(E.subs(r, 0)), (u_cc, u_cz, u_zz), dict=True
)
expected_rank_one_solution = {
    u_cc: -3 * b_c - sp.Rational(3658, 585),
    u_cz: sp.Rational(7, 4) - sp.Rational(3, 2) * b_z,
    u_zz: sp.Rational(58, 65),
}
require("rank-one affine solution count", len(rank_one_solution) == 1)
for variable in (u_cc, u_cz, u_zz):
    require(
        f"rank-one solution {variable}",
        sp.factor(rank_one_solution[0][variable] - expected_rank_one_solution[variable])
        == 0,
    )
require(
    "rank-one debt",
    sp.factor(debt.subs(r, 0).subs(rank_one_solution[0]))
    == -sp.Rational(2072, 81),
)

# The n_c' contribution is load-bearing: deleting it changes the forced
# constant by exactly three.
R2_without_normal_derivative = R2 - normal_c_prime / 3
mutated_debt = sp.factor(
    (lambda_row * R2_without_normal_derivative)[0]
    .subs(rank_two_solution[0])
)
require("normal derivative mutation active", mutated_debt == -sp.Rational(1829, 81))
print("PASS raw_lambda_J2=(162rac+216raz-24bc-162bz-8ucc-108ucz-72uzz-27)/9")
print("PASS lambda_J2=-2072/81+(5/27)J0prime_minus-(13/27)J0prime_plus")
print("PASS r_nonzero_rank2_and_r_zero_rank1_both_debt=-2072/81")
print("PASS all_F2_G2_F3_G3_lifts_annihilated_by_lambda")
print("PASS active_mutation_without_ncprime=-1829/81")


print("SECTION J1 implication and local boundary controls")
J1_values = (
    2 * hadamard(dU, E_values)
    + dA / 3
    - 2 * hadamard(C_values, dV)
)
rank_two_j1 = sp.simplify(J1_values.subs(rank_two_solution[0]))
rank_two_x0 = (65 * u_zz - 58) / (1170 * r**2)
rank_two_y0 = (-6 * b_z - 144 * r**2 - 4 * u_cz + 7) / (36 * r)
require(
    "rank-two retained J1 local values",
    sp.simplify(rank_two_j1.subs({x0: rank_two_x0, y0: rank_two_y0}))
    == sp.zeros(3, 1),
)
rank_one_j1 = sp.simplify(
    J1_values.subs(r, 0).subs(rank_one_solution[0]).subs({a_z: 0, y0: -a_c / 6})
)
require("rank-one retained J1 boundary", rank_one_j1 == sp.zeros(3, 1))
require("forced debt nonzero", -sp.Rational(2072, 81) != 0)
print("PASS J1_is_not_needed_for_debt; if_J1=0_then_J2_cannot_vanish")
print("PASS retained_local_J1_controls_rank2_and_rank1")
print("PASS scope=fixed_Q1_compiler_quadratic_fold_JC2_OPEN")


print("SECTION source AST gate")
source_text = Path(__file__).read_text(encoding="utf-8")
tree = ast.parse(source_text)
assertion_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
require("no assertion statements", assertion_nodes == 0)
print(f"PASS ast_assertion_nodes={assertion_nodes}")

print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
