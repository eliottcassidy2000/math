#!/usr/bin/env python3
"""Exact controls for proved THM-3635.

This companion verifies the retained jet plane, the minimal Hermite curve,
an exact finite module/SAGBI presentation of its actual restriction ring,
normalization and conductor certificates, and the exhaustive fixed-pair
degree-93/94 determinant boundary.  It constructs the witness by a rational
RREF in the full 41-dimensional normalization quotient; no floating point or
bounded target-monomial heuristic is used.
"""

import ast
import hashlib
import math
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


def monic(poly, x):
    """Return a rational univariate polynomial with leading coefficient one."""
    value = sp.Poly(poly, x, domain=sp.QQ)
    return sp.Poly(value.as_expr() / value.LC(), x, domain=sp.QQ)


def as_fmpq(value):
    """Convert an exact SymPy rational to a python-flint rational."""
    value = sp.Rational(value)
    return fmpq(int(value.p), int(value.q))


def from_fmpq(value):
    """Convert a python-flint rational to a SymPy rational."""
    return sp.Rational(int(value.numerator), int(value.denominator))


def coefficient_hash(poly, x):
    """Hash the complete descending reduced-rational coefficient list."""
    value = sp.Poly(poly, x, domain=sp.QQ)
    payload = ";".join(
        f"{int(sp.Rational(q).p)}/{int(sp.Rational(q).q)}"
        for q in value.all_coeffs()
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def matrix_rank(columns, rhs=None):
    """Exact rank of polynomial coefficient columns, optionally augmented."""
    max_degree = max(column.degree() for column in columns)
    rows = []
    for degree in range(max_degree + 1):
        row = [as_fmpq(column.nth(degree)) for column in columns]
        if rhs is not None:
            row.append(as_fmpq(rhs.nth(degree)))
        rows.extend(row)
    width = len(columns) + (1 if rhs is not None else 0)
    matrix = fmpq_mat(max_degree + 1, width, rows)
    return (max_degree + 1, len(columns)), matrix.rank()


print("THM-3635 exact companion -- proved retained-curve actual-rank witness")
print("status=PROVED + VERIFIED-EXACT + HOSTILE-AUDITED")


print("SECTION retained one-jet plane")
v_minus, u, v_plus = sp.symbols("v_minus u v_plus")
tau_c = sp.Matrix([2 * (v_minus + 6), 3, 2 * (6 - v_plus)])
tau_e = sp.Matrix([-2 * (v_minus + 9), 4 * u, 2 * (9 - v_plus)])
ones = sp.ones(3, 1)
slope_debt = (
    4 * u * (v_minus + v_plus)
    + 4 * v_minus * v_plus
    - 27 * v_minus
    + 27 * v_plus
    - 162
)
require(
    "slope determinant",
    sp.expand(sp.Matrix.hstack(ones, tau_c, tau_e).det() - 2 * slope_debt) == 0,
)

x = sp.symbols("x")
L = sp.Poly(x * (x**2 - 1), x, domain=sp.QQ)
hostile_U = sp.expand(L.as_expr() * (3 * x**2 - 2) / 2)
hostile_V = L.as_expr()
hostile_A = sp.expand(
    L.as_expr() * (sp.Rational(225, 8) * x**3 - sp.Rational(75, 4) * x)
)
hostile_B = sp.expand(1 + sp.Rational(45, 4) * x * L.as_expr())


def derivative_column(poly):
    return sp.Matrix([sp.diff(poly, x).subs(x, point) for point in (-1, 0, 1)])


hostile_derivatives = sp.Matrix.hstack(
    derivative_column(hostile_U),
    derivative_column(hostile_V),
    derivative_column(hostile_A),
    derivative_column(hostile_B),
)
require("enlarged hostile derivative rank", hostile_derivatives.rank() == 3)
require(
    "enlarged hostile rank-three minor",
    hostile_derivatives[:, :3].det() == -sp.Rational(225, 2),
)
print("PASS jet_plane=span(tau_c,tau_e) slope_debt=det/2")
print("PASS enlarged_hostile derivative_rank=3 minor_UVA=-225/2")


print("SECTION minimal Hermite polynomial and compiler")
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
hermite_rows = []
for point in (-1, 0, 1):
    hermite_rows.append([sp.Rational(point) ** degree for degree in range(6)])
for point in (-1, 0, 1):
    hermite_rows.append(
        [0 if degree == 0 else degree * sp.Rational(point) ** (degree - 1) for degree in range(6)]
    )
hermite_matrix = sp.Matrix(hermite_rows)
require("Hermite interpolation invertible", hermite_matrix.det() != 0)
require(
    "Hermite values",
    [Q.eval(point) for point in (-1, 0, 1)]
    == [-3, -sp.Rational(3, 4), -3],
)
require(
    "Hermite derivatives",
    [Q.diff().eval(point) for point in (-1, 0, 1)]
    == [-sp.Rational(9, 2), 1, sp.Rational(9, 2)],
)
require("minimal Hermite degree", Q.degree() == 5 and Q.LC() != 0)

D = sp.Poly(1 + x**2 * Q.as_expr(), x, domain=sp.QQ)
b = sp.Poly((D.as_expr() - 1) * (D.as_expr() + 2) ** 2, x, domain=sp.QQ)
c = sp.Poly(x * D.as_expr() * (D.as_expr() + 2), x, domain=sp.QQ)
e = sp.Poly(Q.as_expr() * (D.as_expr() + 3), x, domain=sp.QQ)
z = sp.Poly(e.as_expr() + 3, x, domain=sp.QQ)

require("restricted compiler relation", c**2 * e == b * (b + 4))
for index, point in enumerate((-1, 0, 1)):
    require(
        f"retained compiler point {index}",
        (b.eval(point), c.eval(point), e.eval(point)) == (0, 0, -3),
    )
require(
    "minimal tangent c",
    [c.diff().eval(point) for point in (-1, 0, 1)] == [3, 3, 3],
)
require(
    "minimal tangent e",
    [e.diff().eval(point) for point in (-1, 0, 1)] == [-9, 4, 9],
)
tangent_rows = [(3, -9), (3, 4), (3, 9)]
tangent_determinants = (
    sp.det(sp.Matrix([tangent_rows[0], tangent_rows[1]])),
    sp.det(sp.Matrix([tangent_rows[1], tangent_rows[2]])),
    sp.det(sp.Matrix([tangent_rows[0], tangent_rows[2]])),
)
require("ordinary triple tangent determinants", tangent_determinants == (39, 15, 54))
print("PASS Q1=unique_degree5 rho0_u1 tangents=(39,15,54)")


print("SECTION finite module and SAGBI certificate")
g35_raw = z**3 - c * b
g44_raw = z**2 * b - c**3
require("g35 cancellation", g35_raw.degree() == 35 and g35_raw.LC() == sp.Rational(9, 2))
require("g44 cancellation", g44_raw.degree() == 44 and g44_raw.LC() == sp.Rational(9, 2))
g35 = monic(g35_raw, x)
g44 = monic(g44_raw, x)
g58_raw = z**5 - c**4 - 9 * c * g44
require("g58 cancellation", g58_raw.degree() == 58 and g58_raw.LC() == sp.Rational(81, 4))
g58 = monic(g58_raw, x)
sagbi_generators = (z, c, b, g35, g44, g58)
require(
    "SAGBI generator degrees",
    tuple(poly.degree() for poly in sagbi_generators) == (12, 15, 21, 35, 44, 58),
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
require(
    "Apery representative degrees",
    tuple(poly.degree() for poly in module_representatives) == apery,
)
require(
    "Apery representatives monic",
    all(poly.LC() == 1 for poly in module_representatives),
)


def basis_by_degree(degree):
    """The monic C[z]-module basis element of a semigroup degree."""
    if degree < 0:
        return None
    residue = degree % 12
    if degree < apery[residue]:
        return None
    power = (degree - apery[residue]) // 12
    return sp.Poly(
        module_representatives[residue] * z**power,
        x,
        domain=sp.QQ,
    )


gaps = tuple(degree for degree in range(1, 77) if basis_by_degree(degree) is None)
expected_gaps = (
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 14, 16, 17, 18, 19, 20,
    22, 23, 25, 26, 28, 29, 31, 32, 34, 37, 38, 40, 41, 43, 46, 49,
    52, 53, 55, 61, 64, 67, 76,
)
require("degree gaps", gaps == expected_gaps)
require("semigroup conductor", all(basis_by_degree(degree) is not None for degree in range(77, 101)))


def normal_form(poly):
    """Reduce modulo the vector subspace S, returning gap coefficients."""
    remainder = sp.Poly(poly, x, domain=sp.QQ)
    gap_part = {}
    while not remainder.is_zero:
        degree = remainder.degree()
        coefficient = remainder.LC()
        basis_poly = basis_by_degree(degree)
        if basis_poly is None:
            gap_part[degree] = gap_part.get(degree, 0) + coefficient
            remainder -= sp.Poly(coefficient * x**degree, x, domain=sp.QQ)
        else:
            remainder -= coefficient * basis_poly
    return {degree: sp.factor(value) for degree, value in gap_part.items() if value != 0}


for residue, representative in enumerate(module_representatives):
    for name, generator in (("z", z), ("c", c), ("b", b)):
        require(
            f"module closure residue={residue} generator={name}",
            normal_form(representative * generator) == {},
        )

require("coprime normalization pole orders", math.gcd(z.degree(), g35.degree()) == 1)
print("PASS SAGBI_degrees=(12,15,21,35,44,58) module_closure=36")
print("PASS apery=(0,73,50,15,88,65,30,79,44,21,58,35) gaps=41 conductor_semigroup=77")
print("PASS normalization=C[x] fraction_field_pole_gcd=1")


print("SECTION global conductor")


def conductor_matrix(cutoff):
    degrees = [degree for degree in range(cutoff + 1) if basis_by_degree(degree) is not None]
    rows = []
    for power in range(1, 12):
        normal_forms = [normal_form(basis_by_degree(degree) * x**power) for degree in degrees]
        for gap in gaps:
            rows.extend(as_fmpq(value.get(gap, 0)) for value in normal_forms)
    matrix = fmpq_mat(11 * len(gaps), len(degrees), rows)
    return degrees, matrix


degrees81, conductor81 = conductor_matrix(81)
degrees82, conductor82 = conductor_matrix(82)
require("conductor cutoff dimensions", (len(degrees81), len(degrees82)) == (41, 42))
require("no conductor through degree 81", conductor81.rank() == 41)
reduced82, conductor_rank82 = conductor82.rref()
require("degree 82 conductor rank", conductor_rank82 == 41)

pivot_rows = {}
for row in range(reduced82.nrows()):
    pivot = next(
        (column for column in range(reduced82.ncols()) if reduced82[row, column]),
        None,
    )
    if pivot is not None:
        pivot_rows[pivot] = row
free_columns = [column for column in range(reduced82.ncols()) if column not in pivot_rows]
require("degree 82 conductor nullity", len(free_columns) == 1)
free_column = free_columns[0]
conductor_vector = [fmpq(0)] * reduced82.ncols()
conductor_vector[free_column] = fmpq(1)
for pivot, row in pivot_rows.items():
    conductor_vector[pivot] = -reduced82[row, free_column]
h = sp.Poly(0, x, domain=sp.QQ)
for value, degree in zip(conductor_vector, degrees82):
    h += from_fmpq(value) * basis_by_degree(degree)
h = monic(h, x)
require("conductor degree", h.degree() == 82)
for power in range(12):
    require(f"conductor module gate power={power}", normal_form(h * x**power) == {})

H76_quotient, H76_remainder = sp.div(h, L**2, domain=sp.QQ)
require("conductor retained factor", H76_remainder.is_zero)
H76 = monic(H76_quotient, x)
require("conductor remote degree", H76.degree() == 76)
require("conductor remote squarefree", sp.gcd(H76, H76.diff()).degree() == 0)
require("conductor factors disjoint", sp.gcd(H76, L).degree() == 0)
require(
    "conductor hash",
    coefficient_hash(h, x) == "23956683350b2a13a5b07b99ec75ac26840a9933b0163e1696c687cdef1124f2",
)
require(
    "H76 hash",
    coefficient_hash(H76, x) == "08b3d82e2ae493e8ffa079c2913a6c2095ac24c1157e8c78c9b2fb2476b4e3b8",
)
print("PASS global_conductor=h*C[x] degree=82 factor=L^2*H76 H76_squarefree_coprimeL=true")
print("PASS hashes h=23956683350b... H76=08b3d82e2ae4...")


print("SECTION actual-ring Bezout correction")
bezout_B, bezout_minus_A, bezout_gcd = sp.gcdex(c.diff(), e.diff())
B0 = sp.Poly(bezout_B, x, domain=sp.QQ)
A0 = sp.Poly(-bezout_minus_A, x, domain=sp.QQ)
require("derivative gcd", sp.Poly(bezout_gcd, x, domain=sp.QQ) == sp.Poly(1, x, domain=sp.QQ))
require("reduced Bezout identity", c.diff() * B0 - e.diff() * A0 == 1)
require("reduced Bezout degrees", A0.degree() < c.diff().degree() and B0.degree() < e.diff().degree())

quotient_columns = []
for degree in range(82):
    c_part = normal_form(c.diff() * x**degree)
    e_part = normal_form(e.diff() * x**degree)
    quotient_columns.append(
        [c_part.get(gap, 0) for gap in gaps]
        + [e_part.get(gap, 0) for gap in gaps]
    )
quotient_rhs = (
    [-normal_form(A0).get(gap, 0) for gap in gaps]
    + [-normal_form(B0).get(gap, 0) for gap in gaps]
)

coefficient_data = [
    as_fmpq(quotient_columns[column][row])
    for row in range(82)
    for column in range(82)
]
augmented_data = [
    entry
    for row in range(82)
    for entry in (
        [as_fmpq(quotient_columns[column][row]) for column in range(82)]
        + [as_fmpq(quotient_rhs[row])]
    )
]
quotient_matrix = fmpq_mat(82, 82, coefficient_data)
quotient_augmented = fmpq_mat(82, 83, augmented_data)
require("quotient correction rank", quotient_matrix.rank() == 81)
reduced_augmented, augmented_rank = quotient_augmented.rref()
require("quotient correction consistency", augmented_rank == 81)

correction_pivots = {}
for row in range(reduced_augmented.nrows()):
    pivot = next(
        (column for column in range(82) if reduced_augmented[row, column]),
        None,
    )
    if pivot is not None:
        correction_pivots[pivot] = row
correction_free = [column for column in range(82) if column not in correction_pivots]
require("correction free column", correction_free == [81])

T_coefficients = [fmpq(0)] * 82
for pivot, row in correction_pivots.items():
    T_coefficients[pivot] = reduced_augmented[row, 82]
T = sp.Poly(
    sum(from_fmpq(value) * x**degree for degree, value in enumerate(T_coefficients)),
    x,
    domain=sp.QQ,
)

reduced_homogeneous, homogeneous_rank = quotient_matrix.rref()
require("homogeneous quotient rank", homogeneous_rank == 81)
homogeneous_pivots = {}
for row in range(reduced_homogeneous.nrows()):
    pivot = next(
        (column for column in range(82) if reduced_homogeneous[row, column]),
        None,
    )
    if pivot is not None:
        homogeneous_pivots[pivot] = row
kernel_vector = [fmpq(0)] * 82
kernel_vector[81] = fmpq(1)
for pivot, row in homogeneous_pivots.items():
    kernel_vector[pivot] = -reduced_homogeneous[row, 81]
K = monic(
    sum(from_fmpq(value) * x**degree for degree, value in enumerate(kernel_vector)),
    x,
)
expected_K = sp.Poly(
    L.as_expr() * (27 * x**2 - 4 * x - 36) * H76.as_expr() / 27,
    x,
    domain=sp.QQ,
)
require("kernel factorization", K == expected_K)
require("kernel c membership", normal_form(c.diff() * K) == {})
require("kernel e membership", normal_form(e.diff() * K) == {})
require(
    "kernel hash",
    coefficient_hash(K, x) == "7f4fe2af8b57163cb7ea49760b46c77feb56df38e96e67b08a0c5b989b6e6e59",
)

A = sp.Poly(A0 + c.diff() * T, x, domain=sp.QQ)
B = sp.Poly(B0 + e.diff() * T, x, domain=sp.QQ)
require("actual A membership", normal_form(A) == {})
require("actual B membership", normal_form(B) == {})
require("actual witness determinant", c.diff() * B - e.diff() * A == 1)
require("actual witness degrees", (T.degree(), A.degree(), B.degree()) == (80, 94, 91))
require("actual A triple value", [A.eval(point) for point in (-1, 0, 1)] == [0, 0, 0])
require(
    "actual B triple value",
    [B.eval(point) for point in (-1, 0, 1)]
    == [sp.Rational(1, 3)] * 3,
)
require(
    "T hash",
    coefficient_hash(T, x) == "b93aa43aec99bfbc94fc7fefd2dd5c57a815d5eb7e1bca1fc7ee2f48fd3fdd66",
)
require(
    "A hash",
    coefficient_hash(A, x) == "c7253d6126f4acd03e38437b630e01d7dc6daeee64e92bf142dec70457ff00b5",
)
require(
    "B hash",
    coefficient_hash(B, x) == "1b330351288f562690ee99c04e8a71d07c2ba0cfffd33f4980041302adda70d1",
)
actual_derivatives = sp.Matrix.hstack(
    derivative_column(c.as_expr()),
    derivative_column(e.as_expr()),
    derivative_column(A.as_expr()),
    derivative_column(B.as_expr()),
)
require("actual derivative plane rank", actual_derivatives.rank() == 2)
print("PASS quotient_correction matrix_rank=81 augmented_rank=81 free_column=x^81")
print("PASS kernel=K degree=81 factor=L*(27x^2-4x-36)*H76/27")
print("PASS T_degree=80 A_degree=94 B_degree=91 triple_values_A=0 B=1/3")
print("PASS determinant=cprime*B-A*eprime=1 actual_derivative_rank=2")
print("PASS hashes T=b93aa43aec99... A=c7253d6126f... B=1b330351288f...")


print("SECTION exhaustive fixed-pair degree boundary")


def filtered_columns(cutoff):
    degrees = [degree for degree in range(cutoff + 1) if basis_by_degree(degree) is not None]
    columns = []
    for degree in degrees:
        basis_poly = basis_by_degree(degree)
        columns.extend((c.diff() * basis_poly, -e.diff() * basis_poly))
    return degrees, columns


constant_one = sp.Poly(1, x, domain=sp.QQ)
degrees93, columns93 = filtered_columns(93)
shape93, rank93 = matrix_rank(columns93)
shape93_augmented, rank93_augmented = matrix_rank(columns93, constant_one)
require("degree 93 dimensions", len(degrees93) == 53 and shape93 == (108, 106))
require("degree 93 augmented shape", shape93_augmented == shape93)
require("degree 93 no solution", (rank93, rank93_augmented) == (106, 107))

degrees94, columns94 = filtered_columns(94)
shape94, rank94 = matrix_rank(columns94)
shape94_augmented, rank94_augmented = matrix_rank(columns94, constant_one)
require("degree 94 dimensions", len(degrees94) == 54 and shape94 == (109, 108))
require("degree 94 augmented shape", shape94_augmented == shape94)
require("degree 94 unique solution", (rank94, rank94_augmented) == (108, 108))
require("constructed witness lies in cutoff 94", A.degree() <= 94 and B.degree() <= 94)
print("PASS cutoff93 dimS=53 matrix=108x106 ranks=(106,107) no_solution")
print("PASS cutoff94 dimS=54 matrix=109x108 ranks=(108,108) unique_solution")
print("PASS minimality_scope=fixed_Q1_fixed_Uc_Ve_symmetric_source_x_degree_only")


print("SECTION source AST gate")
source_text = Path(__file__).read_text(encoding="utf-8")
tree = ast.parse(source_text)
assertion_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
require("no assertion statements", assertion_nodes == 0)
print(f"PASS ast_assertion_nodes={assertion_nodes}")

print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
