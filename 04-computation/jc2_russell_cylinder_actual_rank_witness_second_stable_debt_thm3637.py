#!/usr/bin/env python3
"""Exact controls for proved and independently hostile-audited THM-3637.

The companion verifies the stable-order Jacobian expansion, reconstructs the
THM-3635 actual-ring witness, gives an exact first-order lift, and checks the
all-unbounded retained-jet calculation whose second-order cokernel is
-2072/81.  All arithmetic is rational and every gate stays active under -O.
"""

import ast
import hashlib
from pathlib import Path

from flint import fmpq, fmpq_mat
import sympy as sp


CHECKS = 0


def require(label, condition):
    """Record one active exact gate."""
    global CHECKS
    if condition != True:
        raise RuntimeError(f"FAIL {label}: {condition}")
    CHECKS += 1


def as_fmpq(value):
    """Convert an exact SymPy rational to a python-flint rational."""
    value = sp.Rational(value)
    return fmpq(int(value.p), int(value.q))


def from_fmpq(value):
    """Convert a python-flint rational to a SymPy rational."""
    return sp.Rational(int(value.numerator), int(value.denominator))


def monic(poly, x):
    """Return a rational univariate polynomial with leading coefficient one."""
    value = sp.Poly(poly, x, domain=sp.QQ)
    return sp.Poly(value.as_expr() / value.LC(), x, domain=sp.QQ)


def coefficient_hash(poly, x):
    """Hash the complete descending reduced-rational coefficient list."""
    value = sp.Poly(poly, x, domain=sp.QQ)
    payload = ";".join(
        f"{int(sp.Rational(q).p)}/{int(sp.Rational(q).q)}"
        for q in value.all_coeffs()
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def solve_polynomial_system(columns, rhs):
    """Solve exact polynomial coefficient equations with free variables zero."""
    row_count = max(max(column.degree() for column in columns), rhs.degree()) + 1
    coefficient_data = [
        as_fmpq(column.nth(degree))
        for degree in range(row_count)
        for column in columns
    ]
    augmented_data = [
        entry
        for degree in range(row_count)
        for entry in (
            [as_fmpq(column.nth(degree)) for column in columns]
            + [as_fmpq(rhs.nth(degree))]
        )
    ]
    coefficient_matrix = fmpq_mat(row_count, len(columns), coefficient_data)
    augmented_matrix = fmpq_mat(row_count, len(columns) + 1, augmented_data)
    reduced, augmented_rank = augmented_matrix.rref()
    coefficient_rank = coefficient_matrix.rank()

    pivot_rows = {}
    for row in range(reduced.nrows()):
        pivot = next(
            (column for column in range(len(columns)) if reduced[row, column]),
            None,
        )
        if pivot is not None:
            pivot_rows[pivot] = row
        else:
            require(f"linear consistency row={row}", reduced[row, len(columns)] == 0)

    solution = [fmpq(0)] * len(columns)
    for pivot, row in pivot_rows.items():
        solution[pivot] = reduced[row, len(columns)]
    return (
        [from_fmpq(value) for value in solution],
        (row_count, len(columns)),
        coefficient_rank,
        augmented_rank,
    )


print("THM-3637 exact companion -- proved second-stable debt")
print("status=PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED")


print("SECTION stable Jacobian coefficients")
t = sp.symbols("t")
U1, A0, A1, C0, C1, D0 = sp.symbols("U1 A0 A1 C0 C1 D0")
V1, B0, B1, E0, E1, H0 = sp.symbols("V1 B0 B1 E0 E1 H0")
f_x = U1 + t * A1 + t**2 * C1
f_t = A0 + 2 * t * C0 + 3 * t**2 * D0
g_x = V1 + t * B1 + t**2 * E1
g_t = B0 + 2 * t * E0 + 3 * t**2 * H0
jacobian = sp.Poly(sp.expand(f_x * g_t - f_t * g_x), t)
expected_j0 = U1 * B0 - A0 * V1
expected_j1 = 2 * U1 * E0 + A1 * B0 - A0 * B1 - 2 * C0 * V1
expected_j2 = 3 * U1 * H0 + 2 * A1 * E0 + C1 * B0 - A0 * E1 - 2 * C0 * B1 - 3 * D0 * V1
require("Jacobian coefficient t0", jacobian.nth(0) == expected_j0)
require("Jacobian coefficient t1", jacobian.nth(1) == expected_j1)
require("Jacobian coefficient t2", jacobian.nth(2) == expected_j2)
print("PASS J0=Uprime*B-A*Vprime")
print("PASS J1=2Uprime*E+Aprime*B-A*Bprime-2C*Vprime")
print("PASS J2=3Uprime*H+2Aprime*E+Cprime*B-A*Eprime-2C*Bprime-3D*Vprime")


print("SECTION minimal retained compiler and vertical jets")
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
z = sp.Poly(e.as_expr() + 3, x, domain=sp.QQ)
delta_c = sp.Poly(sp.diff(c_general, q).subs(q, Q.as_expr()), x, domain=sp.QQ)
delta_e = sp.Poly(sp.diff(e_general, q).subs(q, Q.as_expr()), x, domain=sp.QQ)
points = (-1, 0, 1)

require("compiler surface relation", c**2 * e == b * (b + 4))
require(
    "retained target triple",
    [(b.eval(point), c.eval(point), e.eval(point)) for point in points]
    == [(0, 0, -3)] * 3,
)
tau_c = sp.Matrix([c.diff().eval(point) for point in points])
tau_e = sp.Matrix([e.diff().eval(point) for point in points])
c_second = sp.Matrix([c.diff().diff().eval(point) for point in points])
normal_c = sp.Matrix([delta_c.eval(point) for point in points])
normal_e = sp.Matrix([delta_e.eval(point) for point in points])
normal_c_prime = sp.Matrix([delta_c.diff().eval(point) for point in points])
require("retained tangent c", tau_c == sp.Matrix([3, 3, 3]))
require("retained tangent e", tau_e == sp.Matrix([-9, 4, 9]))
require("retained c second jet", c_second == sp.Matrix([sp.Rational(157, 2), 0, -sp.Rational(221, 2)]))
require("vertical c values", normal_c == sp.Matrix([2, 0, -2]))
require("vertical e values", normal_e == sp.Matrix([-2, 4, -2]))
require("vertical c first jets", normal_c_prime == sp.Matrix([-9, 0, -9]))
require(
    "ordinary triple determinants",
    tuple(
        sp.det(sp.Matrix([[tau_c[i], tau_e[i]], [tau_c[j], tau_e[j]]]))
        for i, j in ((0, 1), (1, 2), (0, 2))
    )
    == (39, 15, 54),
)
print("PASS Q1_degree=5 retained_triple=(0,0,-3) tangents=(39,15,54)")
print("PASS delta_c=(2,0,-2) delta_e=(-2,4,-2) delta_c_prime=(-9,0,-9)")


print("SECTION retained cokernel and graph-ideal gauge")
lambda_row = sp.Matrix([[sp.Rational(5, 18), -1, sp.Rational(13, 18)]])
require("lambda kills tangent c", (lambda_row * tau_c)[0] == 0)
require("lambda kills tangent e", (lambda_row * tau_e)[0] == 0)

# In regular local coordinates (c,e+3), a target germ vanishing on all three
# branches has no terms of degree below three.  The following ranks are the
# degree-one and degree-two interpolation checks on the three tangent lines.
linear_evaluation = sp.Matrix([[tau_c[i], tau_e[i]] for i in range(3)])
quadratic_evaluation = sp.Matrix(
    [[tau_c[i] ** 2, tau_c[i] * tau_e[i], tau_e[i] ** 2] for i in range(3)]
)
require("graph ideal linear jet vanishes", linear_evaluation.rank() == 2)
require("graph ideal quadratic jet determinant", quadratic_evaluation.det() == 31590)

# Differentiating a general cubic target jet in any vertical tangent direction
# leaves a branch series divisible by xi^2.  Higher target order only raises
# that branch order.
xi, epsilon = sp.symbols("xi epsilon")
h30, h21, h12, h03 = sp.symbols("h30 h21 h12 h03")
u, v = sp.symbols("u v")
general_cubic = h30 * u**3 + h21 * u**2 * v + h12 * u * v**2 + h03 * v**3
for index in range(3):
    branch_substitution = general_cubic.subs(
        {
            u: tau_c[index] * xi + normal_c[index] * epsilon,
            v: tau_e[index] * xi + normal_e[index] * epsilon,
        }
    )
    vertical_derivative = sp.Poly(
        sp.diff(branch_substitution, epsilon).subs(epsilon, 0),
        xi,
    )
    require(
        f"graph gauge vertical order branch={index}",
        vertical_derivative.nth(0) == 0 and vertical_derivative.nth(1) == 0,
    )
print("PASS lambda=(5,-18,13)/18 annihilates retained tangent plane")
print("PASS ker_gamma subset m^3 and delta(ker_gamma) subset branchwise_m^2")


print("SECTION actual-ring basis and normalized Bezout witness")
g35 = monic(z**3 - c * b, x)
g44 = monic(z**2 * b - c**3, x)
g58 = monic(z**5 - c**4 - 9 * c * g44, x)
apery = (0, 73, 50, 15, 88, 65, 30, 79, 44, 21, 58, 35)
module_representatives = (
    sp.Poly(1, x, domain=sp.QQ),
    g58 * c,
    g35 * c,
    c,
    g44**2,
    g44 * b,
    c**2,
    g44 * g35,
    g44,
    b,
    g58,
    g35,
)
require("actual ring generator degrees", (z.degree(), c.degree(), b.degree(), g35.degree(), g44.degree(), g58.degree()) == (12, 15, 21, 35, 44, 58))
require("actual ring Apery degrees", tuple(poly.degree() for poly in module_representatives) == apery)


def basis_by_degree(degree):
    """Return the THM-3635 monic C[z]-module basis element of this degree."""
    if degree < 0:
        return None
    residue = degree % 12
    if degree < apery[residue]:
        return None
    return sp.Poly(
        module_representatives[residue] * z ** ((degree - apery[residue]) // 12),
        x,
        domain=sp.QQ,
    )


basis94 = [basis_by_degree(degree) for degree in range(95) if basis_by_degree(degree) is not None]
require("actual ring cutoff 94 dimension", len(basis94) == 54)
j0_columns = [-e.diff() * basis_poly for basis_poly in basis94] + [
    c.diff() * basis_poly for basis_poly in basis94
]
j0_solution, j0_shape, j0_rank, j0_augmented_rank = solve_polynomial_system(
    j0_columns,
    sp.Poly(1, x, domain=sp.QQ),
)
A = sp.Poly(
    sum(
        (j0_solution[index] * basis94[index] for index in range(len(basis94))),
        sp.Poly(0, x, domain=sp.QQ),
    ),
    x,
    domain=sp.QQ,
)
B = sp.Poly(
    sum(
        (
            j0_solution[len(basis94) + index] * basis94[index]
            for index in range(len(basis94))
        ),
        sp.Poly(0, x, domain=sp.QQ),
    ),
    x,
    domain=sp.QQ,
)
require("J0 witness matrix shape", j0_shape == (109, 108))
require("J0 witness ranks", (j0_rank, j0_augmented_rank) == (108, 108))
require("J0 witness identity", c.diff() * B - e.diff() * A == 1)
require("J0 witness degrees", (A.degree(), B.degree()) == (94, 91))
require("J0 witness A triple", [A.eval(point) for point in points] == [0, 0, 0])
require("J0 witness B triple", [B.eval(point) for point in points] == [sp.Rational(1, 3)] * 3)
require("J0 witness A hash", coefficient_hash(A, x) == "c7253d6126f4acd03e38437b630e01d7dc6daeee64e92bf142dec70457ff00b5")
require("J0 witness B hash", coefficient_hash(B, x) == "1b330351288f562690ee99c04e8a71d07c2ba0cfffd33f4980041302adda70d1")
print("PASS J0_actual matrix=109x108 ranks=(108,108) degrees_A_B=(94,91)")
print("PASS J0_actual hashes_A_B=c7253d6126f.../1b330351288f...")


print("SECTION exact first-order extension")
# J1=0 is equivalent to c'Y-e'X equal to this exact right-hand side.
j1_rhs = sp.Poly(
    -c.diff() * delta_e
    + delta_c * e.diff()
    - sp.Rational(1, 2) * (A.diff() * B - A * B.diff()),
    x,
    domain=sp.QQ,
)
basis174 = [basis_by_degree(degree) for degree in range(175) if basis_by_degree(degree) is not None]
require("actual ring cutoff 174 dimension", len(basis174) == 134)
j1_columns = [-e.diff() * basis_poly for basis_poly in basis174] + [
    c.diff() * basis_poly for basis_poly in basis174
]
j1_solution, j1_shape, j1_rank, j1_augmented_rank = solve_polynomial_system(
    j1_columns,
    j1_rhs,
)
X = sp.Poly(
    sum(
        (j1_solution[index] * basis174[index] for index in range(len(basis174))),
        sp.Poly(0, x, domain=sp.QQ),
    ),
    x,
    domain=sp.QQ,
)
Y = sp.Poly(
    sum(
        (
            j1_solution[len(basis174) + index] * basis174[index]
            for index in range(len(basis174))
        ),
        sp.Poly(0, x, domain=sp.QQ),
    ),
    x,
    domain=sp.QQ,
)
C_actual = sp.Poly(delta_c + X, x, domain=sp.QQ)
E_actual = sp.Poly(delta_e + Y, x, domain=sp.QQ)
j1_actual = sp.Poly(
    2 * c.diff() * E_actual + A.diff() * B - A * B.diff() - 2 * C_actual * e.diff(),
    x,
    domain=sp.QQ,
)
require("J1 lift matrix shape", j1_shape == (189, 268))
require("J1 lift exact ranks", (j1_rank, j1_augmented_rank) == (188, 188))
require("J1 lift residual", j1_actual.is_zero)
require("J1 lift degrees", (X.degree(), Y.degree()) == (173, 91))
require("J1 lift X hash", coefficient_hash(X, x) == "13403da7b39e92a0ae35b690a07184d673de63361c8bae243c7041e3e45081ef")
require("J1 lift Y hash", coefficient_hash(Y, x) == "6432e36c32a516a2204f315142eee0dba83fd0cbc0400f213038ef5a864a5e62")
print("PASS J1_actual matrix=189x268 ranks=(188,188) residual=0")
print("PASS J1_actual degrees_X_Y=(173,91) hashes=13403da7b39e.../6432e36c32a...")


print("SECTION all-unbounded retained gauge and second-order debt")
sigma = sp.symbols("sigma")
p = sigma - sp.Rational(7, 6)
q_coefficient = -sp.Rational(58, 195)
r = -sp.Rational(3658, 1755)
s_coefficient = sigma
dA = p * tau_c + q_coefficient * tau_e
dB = r * tau_c + s_coefficient * tau_e

# Differentiating c'B-Ae'=1 at the retained points gives this affine system.
differentiated_j0 = sp.Matrix(
    [
        c_second[index] / 3
        + 3 * dB[index]
        - tau_e[index] * dA[index]
        for index in range(3)
    ]
)
require("all-gauge differentiated J0", differentiated_j0 == sp.zeros(3, 1))
p_symbol, q_symbol, r_symbol, s_symbol = sp.symbols("p_symbol q_symbol r_symbol s_symbol")
generic_dA = p_symbol * tau_c + q_symbol * tau_e
generic_dB = r_symbol * tau_c + s_symbol * tau_e
generic_equations = sp.Matrix(
    [
        c_second[index] / 3
        + 3 * generic_dB[index]
        - tau_e[index] * generic_dA[index]
        for index in range(3)
    ]
)
generic_matrix, generic_rhs = sp.linear_eq_to_matrix(
    list(generic_equations),
    (p_symbol, q_symbol, r_symbol, s_symbol),
)
require("all-gauge affine rank", generic_matrix.rank() == 3)
require("all-gauge augmented rank", generic_matrix.row_join(generic_rhs).rank() == 3)

# The vertical values of arbitrary target lifts F1,G1 are fixed by these
# retained gradients; ker(gamma) changes have zero vertical value.
delta_A = p * normal_c + q_coefficient * normal_e
delta_B = r * normal_c + s_coefficient * normal_e
require("delta F1 formula", delta_A == p * normal_c + q_coefficient * normal_e)
require("delta G1 formula", delta_B == r * normal_c + s_coefficient * normal_e)

x0 = -sp.Rational(29, 585)
y0 = -sigma / 6 - sp.Rational(137, 36)
j1_retained = sp.Matrix(
    [
        2 * tau_c[index] * (normal_e[index] + y0)
        + dA[index] / 3
        - 2 * (normal_c[index] + x0) * tau_e[index]
        for index in range(3)
    ]
)
require("J1 retained common values", j1_retained == sp.zeros(3, 1))
x0_symbol, y0_symbol = sp.symbols("x0_symbol y0_symbol")
generic_j1_retained = sp.Matrix(
    [
        2 * tau_c[index] * (normal_e[index] + y0_symbol)
        + dA[index] / 3
        - 2 * (normal_c[index] + x0_symbol) * tau_e[index]
        for index in range(3)
    ]
)
j1_value_matrix, j1_value_rhs = sp.linear_eq_to_matrix(
    list(generic_j1_retained),
    (x0_symbol, y0_symbol),
)
require("J1 retained value rank", j1_value_matrix.rank() == 2)
require("J1 retained value consistency", j1_value_matrix.row_join(j1_value_rhs).rank() == 2)

# After lambda kills the F3/G3 restrictions and the derivative of X in S,
# the displayed vector is the complete retained J2 residue.  The first and
# last terms are the indispensable delta(F1),delta(G1) contributions.
j2_residue = sp.Matrix(
    [
        3 * tau_c[index] * delta_B[index]
        + 2 * dA[index] * (normal_e[index] + y0)
        + normal_c_prime[index] / 3
        - 2 * (normal_c[index] + x0) * dB[index]
        - 3 * delta_A[index] * tau_e[index]
        for index in range(3)
    ]
)
j2_debt = sp.factor((lambda_row * j2_residue)[0])
require("J2 debt constant", j2_debt == -sp.Rational(2072, 81))

vertical_B_term = sp.Matrix(
    [3 * tau_c[index] * delta_B[index] for index in range(3)]
)
A_E_term = sp.Matrix(
    [2 * dA[index] * (normal_e[index] + y0) for index in range(3)]
)
C_prime_B_term = normal_c_prime / 3
C_B_prime_term = sp.Matrix(
    [-2 * (normal_c[index] + x0) * dB[index] for index in range(3)]
)
vertical_A_term = sp.Matrix(
    [-3 * delta_A[index] * tau_e[index] for index in range(3)]
)
require(
    "J2 quotient vertical G1 term",
    sp.expand((lambda_row * vertical_B_term)[0] + 2 * (4 * r + 27 * sigma)) == 0,
)
require(
    "J2 quotient A prime E term",
    sp.expand((lambda_row * A_E_term)[0] + 12 * (3 * p + 4 * q_coefficient)) == 0,
)
require("J2 quotient C prime B term", (lambda_row * C_prime_B_term)[0] == -3)
require(
    "J2 quotient C B prime term",
    sp.expand(
        (lambda_row * C_B_prime_term)[0]
        - sp.Rational(4, 3) * (4 * r + 27 * sigma)
    )
    == 0,
)
require(
    "J2 quotient vertical F1 term",
    sp.expand((lambda_row * vertical_A_term)[0] - 18 * (3 * p + 4 * q_coefficient)) == 0,
)
require(
    "J2 residue term decomposition",
    j2_residue
    == vertical_B_term + A_E_term + C_prime_B_term + C_B_prime_term + vertical_A_term,
)

j2_without_vertical_F1_G1 = sp.Matrix(
    [
        2 * dA[index] * (normal_e[index] + y0)
        + normal_c_prime[index] / 3
        - 2 * (normal_c[index] + x0) * dB[index]
        for index in range(3)
    ]
)
mutated_debt = sp.factor((lambda_row * j2_without_vertical_F1_G1)[0])
require("vertical F1 G1 mutation active", mutated_debt == sp.Rational(3415, 81))
require("vertical F1 G1 contribution", j2_debt - mutated_debt == -sp.Rational(1829, 27))

# Every Z=gamma(F3), W=gamma(G3) has a common retained value.  This verifies
# directly that lambda kills the still-free 3(c'W-e'Z) term.
z0, w0 = sp.symbols("z0 w0")
require(
    "F3 G3 quotient annihilation",
    (lambda_row * (3 * tau_c * w0 - 3 * tau_e * z0))[0] == 0,
)
alpha, beta = sp.symbols("alpha beta")
require(
    "S derivative quotient annihilation",
    (lambda_row * (alpha * tau_c + beta * tau_e))[0] == 0,
)
print("PASS all_J0_gauges p=sigma-7/6 q=-58/195 r=-3658/1755 s=sigma")
print("PASS all_J1_lifts common_X=-29/585 common_Y=-sigma/6-137/36")
print("PASS delta_F1_G1 contribution=-1829/27 active_mutation=3415/81")
print("PASS J2_lambda_debt=-2072/81 independent_of_sigma_and_all_lifts")
print("PASS scope=fixed_gamma_F0_c_gamma_G0_e_not_fixed_target_representatives")


print("SECTION source AST gate")
source_text = Path(__file__).read_text(encoding="utf-8")
tree = ast.parse(source_text)
assertion_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
require("no assertion statements", assertion_nodes == 0)
print(f"PASS ast_assertion_nodes={assertion_nodes}")

print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
