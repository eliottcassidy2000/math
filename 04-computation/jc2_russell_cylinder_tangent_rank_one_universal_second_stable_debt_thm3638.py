#!/usr/bin/env python3
"""Exact controls for proved and independently hostile-audited THM-3638.

The main calculation is symbolic and unbounded: after the retained-tangent
rank-one cell is normalized by constant output changes, the differentiated
zero-stable determinant identity forces a Hessian invoice whose substitution
in the second-stable retained quotient is always -2072/81.  A finite exact
sidecar reconstructs the fixed actual-ring control U=c, V=(e+3)^2 and its
sharp symmetric degree-101/102 boundary.
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


def monic(poly, variable):
    """Normalize a rational univariate polynomial to leading coefficient one."""
    value = sp.Poly(poly, variable, domain=sp.QQ)
    return sp.Poly(value.as_expr() / value.LC(), variable, domain=sp.QQ)


def coefficient_hash(poly, variable):
    """Hash the complete descending reduced-rational coefficient list."""
    value = sp.Poly(poly, variable, domain=sp.QQ)
    payload = ";".join(
        f"{int(sp.Rational(entry).p)}/{int(sp.Rational(entry).q)}"
        for entry in value.all_coeffs()
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def rational_hash(value):
    """Hash one reduced rational in p/q form."""
    value = sp.Rational(value)
    payload = f"{int(value.p)}/{int(value.q)}"
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def matrix_rank(columns, rhs=None):
    """Exact rank of polynomial coefficient columns, optionally augmented."""
    maximum_degree = max(column.degree() for column in columns)
    rows = []
    for degree in range(maximum_degree + 1):
        row = [as_fmpq(column.nth(degree)) for column in columns]
        if rhs is not None:
            row.append(as_fmpq(rhs.nth(degree)))
        rows.extend(row)
    width = len(columns) + (1 if rhs is not None else 0)
    matrix = fmpq_mat(maximum_degree + 1, width, rows)
    return (maximum_degree + 1, len(columns)), matrix.rank()


print("THM-3638 exact companion -- tangent-rank-one second-stable debt")
print("status=PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED")


print("SECTION stable expansion")
t = sp.symbols("t")
U_p, V_p, A, A_p, B, B_p = sp.symbols("U_p V_p A A_p B B_p")
C, C_p, E, E_p, D, D_p, H, H_p = sp.symbols(
    "C C_p E E_p D D_p H H_p"
)
f_x = U_p + t * A_p + t**2 * C_p + t**3 * D_p
f_t = A + 2 * t * C + 3 * t**2 * D
g_x = V_p + t * B_p + t**2 * E_p + t**3 * H_p
g_t = B + 2 * t * E + 3 * t**2 * H
jacobian_series = sp.Poly(sp.expand(f_x * g_t - f_t * g_x), t)
J0 = U_p * B - A * V_p
J1 = 2 * U_p * E + A_p * B - A * B_p - 2 * C * V_p
J2 = 3 * U_p * H + 2 * A_p * E + C_p * B - A * E_p - 2 * C * B_p - 3 * D * V_p
require("stable coefficient J0", jacobian_series.nth(0) == J0)
require("stable coefficient J1", jacobian_series.nth(1) == J1)
require("stable coefficient J2", jacobian_series.nth(2) == J2)
print("PASS stable_coefficients=J0,J1,J2")


print("SECTION compiler and retained jets")
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
D_q = 1 + x**2 * q
b_q = sp.expand((D_q - 1) * (D_q + 2) ** 2)
c_q = sp.expand(x * D_q * (D_q + 2))
z_q = sp.expand(q * (D_q + 3) + 3)


def restrict(expression):
    """Restrict a compiler expression to q=Q(x)."""
    return sp.Poly(sp.expand(expression.subs(q, Q.as_expr())), x, domain=sp.QQ)


b = restrict(b_q)
c = restrict(c_q)
z = restrict(z_q)
delta_b = restrict(sp.diff(b_q, q))
delta_c = restrict(sp.diff(c_q, q))
delta_z = restrict(sp.diff(z_q, q))
points = (-1, 0, 1)


def values(poly):
    """Evaluate a polynomial at the retained triple."""
    poly = sp.Poly(poly, x, domain=sp.QQ)
    return sp.Matrix([poly.eval(point) for point in points])


def lam(vector):
    """The retained tangent-cokernel row."""
    vector = sp.Matrix(vector)
    return sp.factor((5 * vector[0] - 18 * vector[1] + 13 * vector[2]) / 18)


ones = sp.ones(3, 1)
tau_c = values(c.diff())
tau_z = values(z.diff())
n_c = values(delta_c)
n_z = values(delta_z)
n_c_prime = values(delta_c.diff())
n_z_prime = values(delta_z.diff())
c_second = values(c.diff().diff())
require("retained c tangent", tau_c == sp.Matrix([3, 3, 3]))
require("retained z tangent", tau_z == sp.Matrix([-9, 4, 9]))
require("vertical c values", n_c == sp.Matrix([2, 0, -2]))
require("vertical z values", n_z == sp.Matrix([-2, 4, -2]))
require("vertical c derivatives", n_c_prime == sp.Matrix([-9, 0, -9]))
require("vertical z derivatives", n_z_prime == sp.Matrix([3, 0, -3]))
require(
    "retained c second derivatives",
    c_second == sp.Matrix([sp.Rational(157, 2), 0, -sp.Rational(221, 2)]),
)
require("lambda annihilates ones", lam(ones) == 0)
require("lambda annihilates tau c", lam(tau_c) == 0)
require("lambda annihilates tau z", lam(tau_z) == 0)
tangent_rows = [(3, -9), (3, 4), (3, 9)]
require(
    "ordinary triple tangent determinants",
    tuple(
        sp.det(sp.Matrix([tangent_rows[i], tangent_rows[j]]))
        for i, j in ((0, 1), (1, 2), (0, 2))
    )
    == (39, 15, 54),
)
print("PASS tangents=tau_c,tau_z verticals=n_c,n_z lambda_annihilator=true")
print("PASS ordinary_triple_determinants=(39,15,54)")


print("SECTION graph-ideal representative invariance")
slopes = (-9, 4, 9)
linear_line_matrix = sp.Matrix([[3, slope] for slope in slopes])
quadratic_line_matrix = sp.Matrix(
    [[9, 3 * slope, slope**2] for slope in slopes]
)
require("three-line linear initial form rank", linear_line_matrix.rank() == 2)
require("three-line quadratic initial form rank", quadratic_line_matrix.det() != 0)
print("PASS ker_gamma_at_ordinary_triple subset m^3")


print("SECTION normalized rank-one Hessian invoice")
u_cc, u_cz, u_zz, b_c, b_z = sp.symbols("u_cc u_cz u_zz b_c b_z")
H_U = sp.Matrix([[u_cc, u_cz], [u_cz, u_zz]])
branch_tangents = [sp.Matrix([3, slope]) for slope in slopes]
vertical_vectors = [
    sp.Matrix([n_c[index], n_z[index]]) for index in range(3)
]
vertical_derivatives = [
    sp.Matrix([n_c_prime[index], n_z_prime[index]]) for index in range(3)
]
U_second = sp.Matrix(
    [
        c_second[index]
        + (branch_tangents[index].T * H_U * branch_tangents[index])[0]
        for index in range(3)
    ]
)
B_prime = sp.Matrix([3 * b_c + slope * b_z for slope in slopes])
J0_prime_invoice = sp.simplify(U_second / 3 + 3 * B_prime)
normalized_hessian_solution = {
    u_cc: -3 * b_c - sp.Rational(3658, 585),
    u_cz: sp.Rational(7, 4) - sp.Rational(3, 2) * b_z,
    u_zz: sp.Rational(58, 65),
}
require(
    "normalized J0 prime Hessian solution",
    J0_prime_invoice.subs(normalized_hessian_solution) == sp.zeros(3, 1),
)
require(
    "normalized Hessian system full rank",
    sp.linear_eq_to_matrix(list(J0_prime_invoice), (u_cc, u_cz, u_zz))[0].rank()
    == 3,
)

delta_B = sp.Matrix(
    [
        b_c * vertical_vectors[index][0] + b_z * vertical_vectors[index][1]
        for index in range(3)
    ]
)
delta_U_prime = sp.Matrix(
    [
        (branch_tangents[index].T * H_U * vertical_vectors[index])[0]
        + vertical_derivatives[index][0]
        for index in range(3)
    ]
)
second_stable_debt = sp.factor(
    9 * lam(delta_B)
    + sp.Rational(1, 3) * lam(delta_U_prime)
    - 2 * lam(sp.matrix_multiply_elementwise(n_c, B_prime))
)
normalized_debt = sp.factor(second_stable_debt.subs(normalized_hessian_solution))
require("universal normalized debt", normalized_debt == -sp.Rational(2072, 81))
require("normalized debt independent of B gradient", not normalized_debt.free_symbols)
mutated_delta_U_prime = delta_U_prime + sp.Matrix([1, 0, 0])
mutated_debt = sp.factor(
    (
        9 * lam(delta_B)
        + sp.Rational(1, 3) * lam(mutated_delta_U_prime)
        - 2 * lam(sp.matrix_multiply_elementwise(n_c, B_prime))
    ).subs(normalized_hessian_solution)
)
require("active vertical-derivative mutation", mutated_debt == -sp.Rational(4129, 162))
print("PASS J0prime forces ucc=-3bc-3658/585 ucz=7/4-3bz/2 uzz=58/65")
print("PASS lambda_J2=-2072/81 independent_of_bc_bz_and_higher_lifts=true")
print("PASS hostile_vertical_derivative_mutation=-4129/162")


print("SECTION unnormalized shear-covariance cross-check")
v_cc, v_cz, v_zz, a_value = sp.symbols("v_cc v_cz v_zz a_value")
H_V = sp.Matrix([[v_cc, v_cz], [v_cz, v_zz]])
V_second = sp.Matrix(
    [
        (branch_tangents[index].T * H_V * branch_tangents[index])[0]
        for index in range(3)
    ]
)
unnormalized_B_prime = sp.simplify((a_value * V_second - U_second / 3) / 3)
unnormalized_compatibility = sp.factor(lam(unnormalized_B_prime))
require(
    "unnormalized J0 compatibility factor",
    sp.factor(
        unnormalized_compatibility
        - (195 * a_value * v_zz - 65 * u_zz + 58) / 9
    )
    == 0,
)
unnormalized_b_z = sp.factor(
    (unnormalized_B_prime[2] - unnormalized_B_prime[0]) / 18
)
unnormalized_b_c = sp.factor(
    (unnormalized_B_prime[1] - 4 * unnormalized_b_z) / 3
)
unnormalized_delta_B = sp.Matrix(
    [
        unnormalized_b_c * vertical_vectors[index][0]
        + unnormalized_b_z * vertical_vectors[index][1]
        for index in range(3)
    ]
)
delta_V_prime = sp.Matrix(
    [
        (branch_tangents[index].T * H_V * vertical_vectors[index])[0]
        for index in range(3)
    ]
)
unnormalized_debt = sp.factor(
    9 * lam(unnormalized_delta_B)
    + sp.Rational(1, 3) * lam(delta_U_prime)
    - a_value * lam(delta_V_prime)
    - 2
    * lam(sp.matrix_multiply_elementwise(n_c, unnormalized_B_prime))
)
require(
    "unnormalized debt formula",
    sp.factor(
        unnormalized_debt
        - 40 * (117 * a_value * v_zz - 39 * u_zz - 17) / 81
    )
    == 0,
)
require(
    "compatibility forces universal debt",
    sp.factor(
        unnormalized_debt
        + sp.Rational(2072, 81)
        + sp.Rational(8, 27) * (65 * u_zz - 58 - 195 * a_value * v_zz)
    )
    == 0,
)
print("PASS unnormalized_factor=J0compatibility*(-8/27) around -2072/81")


print("SECTION actual-ring rank-one t0 control")
g35 = monic(z**3 - c * b, x)
g44 = monic(z**2 * b - c**3, x)
g58 = monic(z**5 - c**4 - 9 * c * g44, x)
require("actual control g35 scale", g35 == sp.Poly(sp.Rational(2, 9) * (z**3 - c * b), x))
require("actual control g44 scale", g44 == sp.Poly(sp.Rational(2, 9) * (z**2 * b - c**3), x))
require(
    "actual control g58 scale",
    g58 == sp.Poly(sp.Rational(4, 81) * (z**5 - c**4 - 9 * c * g44), x),
)
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


def basis_by_degree(degree):
    """The THM-3635 monic C[z]-module basis element of a semigroup degree."""
    if degree < 0:
        return None
    residue = degree % 12
    if degree < apery[residue]:
        return None
    power = (degree - apery[residue]) // 12
    return sp.Poly(module_representatives[residue] * z**power, x, domain=sp.QQ)


actual_U = c
actual_V = z**2
constant_one = sp.Poly(1, x, domain=sp.QQ)


def determinant_columns(cutoff):
    """Columns for U'B-A V' on the complete symmetric filtered piece."""
    degrees = [
        degree for degree in range(cutoff + 1) if basis_by_degree(degree) is not None
    ]
    columns = []
    for degree in degrees:
        basis_poly = basis_by_degree(degree)
        columns.extend((actual_U.diff() * basis_poly, -actual_V.diff() * basis_poly))
    return degrees, columns


degrees101, columns101 = determinant_columns(101)
shape101, rank101 = matrix_rank(columns101)
_, augmented_rank101 = matrix_rank(columns101, constant_one)
require("actual cutoff 101 dimensions", len(degrees101) == 61 and shape101 == (125, 122))
require("actual cutoff 101 no solution", (rank101, augmented_rank101) == (122, 123))

degrees102, columns102 = determinant_columns(102)
shape102, rank102 = matrix_rank(columns102)
_, augmented_rank102 = matrix_rank(columns102, constant_one)
require("actual cutoff 102 dimensions", len(degrees102) == 62 and shape102 == (126, 124))
require("actual cutoff 102 unique solution", (rank102, augmented_rank102) == (124, 124))

augmented_rows = []
for degree in range(shape102[0]):
    augmented_rows.extend(
        [as_fmpq(column.nth(degree)) for column in columns102]
        + [as_fmpq(1 if degree == 0 else 0)]
    )
augmented_matrix = fmpq_mat(shape102[0], len(columns102) + 1, augmented_rows)
reduced_augmented, reduced_rank = augmented_matrix.rref()
require("actual cutoff 102 RREF rank", reduced_rank == len(columns102))
solution = [sp.Rational(0)] * len(columns102)
pivot_columns = set()
for row in range(reduced_augmented.nrows()):
    pivot = next(
        (
            column
            for column in range(len(columns102))
            if reduced_augmented[row, column]
        ),
        None,
    )
    if pivot is not None:
        pivot_columns.add(pivot)
        solution[pivot] = from_fmpq(
            reduced_augmented[row, len(columns102)]
            / reduced_augmented[row, pivot]
        )
require("actual cutoff 102 no free columns", len(pivot_columns) == len(columns102))

actual_A = sp.Poly(0, x, domain=sp.QQ)
actual_B = sp.Poly(0, x, domain=sp.QQ)
for index, degree in enumerate(degrees102):
    actual_B += solution[2 * index] * basis_by_degree(degree)
    actual_A += solution[2 * index + 1] * basis_by_degree(degree)
require(
    "actual rank-one determinant",
    actual_U.diff() * actual_B - actual_A * actual_V.diff() == constant_one,
)
require("actual witness degrees", (actual_A.degree(), actual_B.degree()) == (93, 102))
require(
    "actual A triple value",
    [actual_A.eval(point) for point in points] == [-sp.Rational(29, 195)] * 3,
)
require(
    "actual B triple value",
    [actual_B.eval(point) for point in points] == [sp.Rational(1, 3)] * 3,
)
require(
    "actual A hash",
    coefficient_hash(actual_A, x)
    == "a17fb9620b7d2ee308f8ca1f311205ac619cc4e70886e641d7ed4b0863d5c73a",
)
require(
    "actual B hash",
    coefficient_hash(actual_B, x)
    == "d43b93e89fcf1aa61bcfce9b3f497ff2d5c723cadec9c7313e52a1238a2ddf26",
)

actual_A_prime = values(actual_A.diff())
actual_B_prime = values(actual_B.diff())
actual_A_z_gradient = sp.factor((actual_A_prime[2] - actual_A_prime[0]) / 18)
actual_B_z_gradient = sp.factor((actual_B_prime[2] - actual_B_prime[0]) / 18)
required_A_z_gradient = sp.factor(
    3 * actual_A.eval(0) * actual_B_z_gradient
)
first_stable_local_mismatch = sp.factor(
    actual_A_z_gradient - required_A_z_gradient
)
require("actual B z gradient", actual_B_z_gradient == sp.Rational(7, 6))
require("actual J1 required A z gradient", required_A_z_gradient == -sp.Rational(203, 390))
require("actual minimal witness fails J1 locally", first_stable_local_mismatch != 0)
require(
    "actual J1 mismatch hash",
    rational_hash(first_stable_local_mismatch)
    == "a48a97b7b9d3d930c9472c610004657e5728cf0dfa0297429a9738482860decc",
)
print("PASS actual_U=c actual_V=z^2 cutoff101_ranks=(122,123)")
print("PASS cutoff102_ranks=(124,124) unique_AB_degrees=(93,102)")
print("PASS actual_AB_hashes=a17fb9620b7d... d43b93e89fcf...")
print("PASS minimal_actual_witness_J1=LOCAL_FAIL mismatch_hash=a48a97b7b9d...")


print("SECTION source AST gate")
source_text = Path(__file__).read_text(encoding="utf-8")
tree = ast.parse(source_text)
assertion_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
require("no assertion statements", assertion_nodes == 0)
print(f"PASS ast_assertion_nodes={assertion_nodes}")

print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
