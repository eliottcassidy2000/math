"""Complete retained order-six two-form cokernel for THM-3681."""

import ast
import hashlib
from pathlib import Path

import sympy as sp
from flint import fmpq, fmpq_mat


CUTOFF = 6
EXPECTED_RELATION_SHA256 = (
    "204b836a68a19c64404a409d7a6298aa370a1ae9a8607b3100d288e2e384c982"
)
x, q, t, X = sp.symbols("x q t X")
POINTS = (-1, 0, 1)
gates = 0


def require(condition, label):
    global gates
    if not condition:
        raise RuntimeError(label)
    gates += 1


Q = sp.expand(
    (
        22868 * x**8
        - 89583 * x**6
        + 2916 * x**5
        + 123684 * x**4
        - 5832 * x**3
        - 63530 * x**2
        + 2916 * x
        - 2187
    )
    / 2916
)

D = 1 + x**2 * q
b = sp.expand((D - 1) * (D + 2) ** 2)
c = sp.expand(x * D * (D + 2))
e = sp.expand(q * (D + 3))
y = c / 3
z = e + 3

require(sp.Poly(Q, x).degree() == 8, "Q_dagger degree")
require(sp.expand(c**2 * e - b * (b + 4)) == 0, "Russell surface relation")
require(tuple(Q.subs(x, point) for point in POINTS) == (-3, sp.Rational(-3, 4), -3), "retained values")
require(
    tuple(sp.diff(Q, x).subs(x, point) for point in POINTS)
    == (sp.Rational(-9, 2), 1, sp.Rational(9, 2)),
    "retained slopes",
)
require(
    5 * sp.diff(Q, x, 2).subs(x, -1)
    + 13 * sp.diff(Q, x, 2).subs(x, 1)
    + 243
    == 0,
    "zero second debt",
)


def multiply(left, right):
    """Multiply in Q[X,t]/(monomials of total degree above six)."""
    answer = {}
    for (i, j), a in left.items():
        for (u, v), b_value in right.items():
            if i + j + u + v > CUTOFF:
                continue
            key = (i + u, j + v)
            answer[key] = answer.get(key, 0) + a * b_value
    return {key: sp.cancel(value) for key, value in answer.items() if value}


def power(value, exponent):
    answer = {(0, 0): sp.S.One}
    for _ in range(exponent):
        answer = multiply(answer, value)
    return answer


def jet(expression, point):
    shifted = sp.Poly(sp.expand(expression.subs(x, X + point)), X, t, domain=sp.QQ)
    return {
        (source_degree, stable_degree): coefficient
        for (source_degree, stable_degree), coefficient in shifted.terms()
        if source_degree + stable_degree <= CUTOFF
    }


pulled_y = sp.expand(y.subs(q, Q + t**2))
pulled_z = sp.expand(z.subs(q, Q + t**2))
y_x = sp.diff(pulled_y, x)
y_t = sp.diff(pulled_y, t)
z_x = sp.diff(pulled_z, x)
z_t = sp.diff(pulled_z, t)
area = sp.expand(y_x * z_t - y_t * z_x)

packets = {}
for point in POINTS:
    y_jet = jet(pulled_y, point)
    z_jet = jet(pulled_z, point)
    require(y_jet.get((0, 0), 0) == 0, f"y vanishes at branch {point}")
    require(z_jet.get((0, 0), 0) == 0, f"z vanishes at branch {point}")
    packets[point] = (
        [power(y_jet, degree) for degree in range(CUTOFF + 1)],
        [power(z_jet, degree) for degree in range(CUTOFF + 1)],
        (jet(area, point), jet(y_x, point), jet(z_x, point)),
    )


def monomial_pullback(kind, y_degree, z_degree, w_degree, point):
    y_powers, z_powers, bases = packets[point]
    value = multiply(y_powers[y_degree], z_powers[z_degree])
    value = multiply(value, {(0, w_degree): sp.S.One})
    return multiply(value, bases[kind])


# Taylor coefficient rows: J_(2j) through x-degree 3-j at each branch.
ROWS = tuple(
    (stable_degree, point, source_degree)
    for stable_degree, max_source_degree in ((0, 3), (2, 2), (4, 1), (6, 0))
    for source_degree in range(max_source_degree + 1)
    for point in POINTS
)

COLUMNS = []
LABELS = []
for kind in range(3):
    for y_degree in range(CUTOFF + 1):
        for z_degree in range(CUTOFF + 1 - y_degree):
            for w_degree in range(CUTOFF + 1 - y_degree - z_degree):
                branch = {
                    point: monomial_pullback(kind, y_degree, z_degree, w_degree, point)
                    for point in POINTS
                }
                COLUMNS.append(
                    tuple(
                        branch[point].get((source_degree, stable_degree), 0)
                        for stable_degree, point, source_degree in ROWS
                    )
                )
                LABELS.append((kind, y_degree, z_degree, w_degree))

require(len(ROWS) == 30, "order-six row count")
require(len(COLUMNS) == 252, "complete order-six two-form monomial count")
require(len(set(LABELS)) == 252, "two-form labels distinct")


def fq(value):
    value = sp.Rational(value)
    return fmpq(int(value.p), int(value.q))


def left_nullspace(source):
    reduced, source_rank = source.transpose().rref()
    pivots = []
    for row in range(source_rank):
        pivot = next(column for column in range(source.nrows()) if reduced[row, column])
        pivots.append(pivot)
    pivot_set = set(pivots)
    free = [column for column in range(source.nrows()) if column not in pivot_set]
    answer = []
    for free_column in free:
        vector = [fmpq(0)] * source.nrows()
        vector[free_column] = fmpq(1)
        for row, pivot in enumerate(pivots):
            vector[pivot] = -reduced[row, free_column]
        column = fmpq_mat(source.nrows(), 1, vector)
        require(
            source.transpose() * column == fmpq_mat(source.ncols(), 1),
            f"left-nullspace vector {free_column}",
        )
        answer.append(vector)
    return source_rank, answer


matrix = fmpq_mat(
    len(ROWS),
    len(COLUMNS),
    [fq(COLUMNS[column][row]) for row in range(len(ROWS)) for column in range(len(COLUMNS))],
)
rank, left_kernel = left_nullspace(matrix)
require(rank == 26, "order-six rank")
require(len(left_kernel) == 4, "order-six left nullity")

constant = [
    fmpq(1) if stable_degree == 0 and source_degree == 0 else fmpq(0)
    for stable_degree, _point, source_degree in ROWS
]
responses = [
    sum((vector[row] * constant[row] for row in range(len(ROWS))), fmpq(0))
    for vector in left_kernel
]
active = [index for index, response in enumerate(responses) if response]
require(len(active) == 1, "unique constant-active left-cokernel line")

vector = left_kernel[active[0]]
j6_indices = [index for index, row in enumerate(ROWS) if row[0] == 6]
lambda_row = (fmpq(5, 18), fmpq(-1), fmpq(13, 18))
scale = lambda_row[0] / vector[j6_indices[0]]
relation = [value * scale for value in vector]
require(tuple(relation[index] for index in j6_indices) == lambda_row, "J6 block is Lambda")

for column, label in zip(COLUMNS, LABELS):
    require(
        sum((relation[row] * fq(column[row]) for row in range(len(ROWS))), fmpq(0))
        == 0,
        f"universal relation on monomial {label}",
    )

payload = ";".join(
    f"{ROWS[index]}:{relation[index]}"
    for index in range(len(ROWS))
    if relation[index]
)
relation_sha256 = hashlib.sha256(payload.encode("ascii")).hexdigest()
require(relation_sha256 == EXPECTED_RELATION_SHA256, "relation hash")

constant_response = responses[active[0]] * scale
require(constant_response == fmpq(326763520, 1594323), "sixth debt value")

row_index = {row: index for index, row in enumerate(ROWS)}
require(
    sum((relation[row_index[(6, point, 0)]] for point in POINTS), fmpq(0)) == 0,
    "Lambda kills stable-only constants",
)
require(
    sum((relation[row_index[(4, point, 0)]] for point in POINTS), fmpq(0)) == 0,
    "R4 value block kills stable-only constants",
)
require(
    sum((relation[row_index[(2, point, 0)]] for point in POINTS), fmpq(0)) == 0,
    "R2 value block kills stable-only constants",
)
require(
    tuple(3 - stable_degree // 2 for stable_degree in (0, 2, 4, 6))
    == (3, 2, 1, 0),
    "quadratic alpha weights",
)

# Active mutation: altering any nonzero row coefficient destroys the identity
# on at least one complete-universe column.
mutated = list(relation)
mutated[0] += 1
require(
    any(
        sum((mutated[row] * fq(column[row]) for row in range(len(ROWS))), fmpq(0))
        != 0
        for column in COLUMNS
    ),
    "active relation mutation",
)

# Sharp predecessor control: the complete order-four window reproduces
# Q_dagger's zero fourth debt, so no constant-active row exists there.
ROWS4 = tuple(
    (stable_degree, point, source_degree)
    for stable_degree, max_source_degree in ((0, 2), (2, 1), (4, 0))
    for source_degree in range(max_source_degree + 1)
    for point in POINTS
)
row4_indices = [ROWS.index(row) for row in ROWS4]
column4_indices = [
    index
    for index, (_kind, y_degree, z_degree, w_degree) in enumerate(LABELS)
    if y_degree + z_degree + w_degree <= 4
]
matrix4 = fmpq_mat(
    len(ROWS4),
    len(column4_indices),
    [
        fq(COLUMNS[column][row])
        for row in row4_indices
        for column in column4_indices
    ],
)
rank4, kernel4 = left_nullspace(matrix4)
constant4 = [
    fmpq(1) if stable_degree == 0 and source_degree == 0 else fmpq(0)
    for stable_degree, _point, source_degree in ROWS4
]
responses4 = [
    sum((vector[row] * constant4[row] for row in range(len(ROWS4))), fmpq(0))
    for vector in kernel4
]
require((len(ROWS4), len(column4_indices), rank4, len(kernel4)) == (18, 105, 15, 3), "order-four control shape")
require(not any(responses4), "order-four constant response vanishes")

source = Path(__file__).read_text(encoding="utf-8")
require(
    sum(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))) == 0,
    "assertion-free AST",
)

print("Qdagger_retained_values_slopes_and_D2=PASS")
print("order4_control=rows18_columns105_rank15_leftnullity3_constantactive0")
print("order6_complete=rows30_columns252_rank26_leftnullity4_constantactive1")
print(f"lambda_relation_sha256={relation_sha256}")
print("forced_lambda_J6=-326763520/1594323=-326763520/3^13")
print(f"RESULT=PASS;gates={gates}")
