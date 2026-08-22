"""Symbolic QQ(r) sixth-debt certificate for THM-3683."""

import ast
import hashlib
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix


CUTOFF = 6
x, q, t, X, r = sp.symbols("x q t X r")
POINTS = (-1, 0, 1)
GATES = 0


def require(condition, label):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(label)


Q1 = (
    x**5
    + sp.Rational(9, 2) * x**4
    - 2 * x**3
    - sp.Rational(27, 4) * x**2
    + x
    - sp.Rational(3, 4)
)
P = x**2 * (x**2 - 1) ** 2
Q6 = sp.expand(Q1 - sp.Rational(259, 36) * P)
R1 = sp.expand(P * (1 - x**2))
R2 = sp.expand(P * (4 - 9 * x))
p_parabola = (
    sp.Rational(520, 9) * r**2
    - sp.Rational(1688, 81) * r
    - sp.Rational(5717, 729)
)
Q = sp.expand(Q6 + p_parabola * R1 + r * R2)

D = 1 + x**2 * q
b = sp.expand((D - 1) * (D + 2) ** 2)
c = sp.expand(x * D * (D + 2))
e = sp.expand(q * (D + 3))
y = c / 3
z = e + 3

require(sp.expand(c**2 * e - b * (b + 4)) == 0, "Russell surface relation")
require(sp.Poly(Q, x).degree() == 8, "family degree")
require(
    tuple(sp.factor(Q.subs(x, point)) for point in POINTS)
    == (-3, sp.Rational(-3, 4), -3),
    "retained values",
)
require(
    tuple(sp.factor(sp.diff(Q, x).subs(x, point)) for point in POINTS)
    == (sp.Rational(-9, 2), 1, sp.Rational(9, 2)),
    "retained slopes",
)
require(
    sp.factor(
        5 * sp.diff(Q, x, 2).subs(x, -1)
        + 13 * sp.diff(Q, x, 2).subs(x, 1)
        + 243
    )
    == 0,
    "zero second debt",
)
fourth_debt = sp.factor(
    sp.Rational(64, 6561)
    * (729 * p_parabola - 42120 * r**2 + 15192 * r + 5717)
)
require(fourth_debt == 0, "zero fourth debt parabola")


def multiply(left, right):
    """Multiply in QQ(r)[X,t] through total degree six."""
    answer = {}
    for (i, j), left_value in left.items():
        for (u, v), right_value in right.items():
            if i + j + u + v > CUTOFF:
                continue
            key = (i + u, j + v)
            answer[key] = answer.get(key, 0) + left_value * right_value
    return {key: sp.factor(value) for key, value in answer.items() if value != 0}


def power(value, exponent):
    answer = {(0, 0): sp.S.One}
    for _ in range(exponent):
        answer = multiply(answer, value)
    return answer


def jet(expression, point):
    shifted = sp.Poly(
        sp.expand(expression.subs(x, X + point)), X, t, domain=sp.QQ[r]
    )
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
    require(y_jet.get((0, 0), 0) == 0, f"y retained zero {point}")
    require(z_jet.get((0, 0), 0) == 0, f"z retained zero {point}")
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
                    point: monomial_pullback(
                        kind, y_degree, z_degree, w_degree, point
                    )
                    for point in POINTS
                }
                COLUMNS.append(
                    tuple(
                        branch[point].get((source_degree, stable_degree), 0)
                        for stable_degree, point, source_degree in ROWS
                    )
                )
                LABELS.append((kind, y_degree, z_degree, w_degree))

require(len(ROWS) == 30, "row count")
require(len(COLUMNS) == 252, "column count")
require(len(set(LABELS)) == 252, "distinct monomial labels")

# COLUMNS are the rows of the transpose of the retained pullback matrix.
matrix = DomainMatrix.from_list_sympy(
    len(COLUMNS), len(ROWS), [list(column) for column in COLUMNS]
).to_field()
require(matrix.domain == sp.QQ.frac_field(r), "QQ(r) domain")
require(matrix.rank() == 26, "symbolic order-six rank")
kernel = matrix.nullspace(divide_last=True).to_Matrix()
require(kernel.shape == (4, 30), "symbolic order-six nullity")

j6_indices = [index for index, row in enumerate(ROWS) if row[0] == 6]
constant_indices = [
    index for index, row in enumerate(ROWS) if row[0] == 0 and row[2] == 0
]
active_rows = [
    row
    for row in range(kernel.rows)
    if any(kernel[row, column] != 0 for column in j6_indices)
]
require(len(active_rows) == 1, "one active kernel basis row")
active = active_rows[0]
scale = sp.Rational(5, 18) / kernel[active, j6_indices[0]]
relation = [sp.factor(kernel[active, column] * scale) for column in range(kernel.cols)]
require(
    tuple(relation[column] for column in j6_indices)
    == (sp.Rational(5, 18), -1, sp.Rational(13, 18)),
    "Lambda J6 block",
)

for column, label in zip(COLUMNS, LABELS):
    require(
        sp.factor(sum(relation[row] * column[row] for row in range(len(ROWS))))
        == 0,
        f"symbolic identity {label}",
    )

nonzero = [(ROWS[index], relation[index]) for index in range(len(ROWS)) if relation[index] != 0]
require(len(nonzero) == 23, "23-entry symbolic relation")
for stable_degree in (2, 4, 6):
    require(
        sp.factor(
            sum(
                relation[index]
                for index, row in enumerate(ROWS)
                if row[0] == stable_degree and row[2] == 0
            )
        )
        == 0,
        f"J{stable_degree} value block kills constants",
    )

debt = sp.factor(sum(relation[index] for index in constant_indices))
debt_quartic = (
    72783360 * r**4
    - 77822208 * r**3
    - 28419741 * r**2
    + 7849770 * r
    - 1276420
)
expected_debt = sp.factor(-sp.Rational(256, 1594323) * debt_quartic)
require(sp.factor(debt - expected_debt) == 0, "sixth debt quartic")
require(
    debt.subs(r, 0) == sp.Rational(326763520, 1594323),
    "Qdagger specialization",
)

# The inherited order-four cokernel is exactly the J6-zero subspace.
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
matrix4 = DomainMatrix.from_list_sympy(
    len(column4_indices),
    len(ROWS4),
    [
        [COLUMNS[column][row] for row in row4_indices]
        for column in column4_indices
    ],
).to_field()
require(matrix4.rank() == 15, "symbolic order-four rank")
kernel4 = matrix4.nullspace(divide_last=True).to_Matrix()
require(kernel4.shape == (3, 18), "symbolic order-four nullity")
embedded = []
for local_kernel_row in range(kernel4.rows):
    extension = [sp.S.Zero for _ in ROWS]
    for local_row, global_row in enumerate(row4_indices):
        extension[global_row] = kernel4[local_kernel_row, local_row]
    require(
        all(
            sp.factor(
                sum(extension[row] * column[row] for row in range(len(ROWS)))
            )
            == 0
            for column in COLUMNS
        ),
        f"order-four embedding {local_kernel_row}",
    )
    require(
        sp.factor(sum(extension[index] for index in constant_indices)) == 0,
        f"order-four constant-inactive {local_kernel_row}",
    )
    embedded.append(extension)
require(sp.Matrix([*embedded, relation]).rank() == 4, "active quotient dimension one")

quartic = sp.Poly(debt_quartic, r, domain=sp.QQ)
require(quartic.is_irreducible, "debt quartic irreducible")
require(quartic.count_roots(-sp.oo, sp.oo) == 2, "two real exceptional roots")
require(quartic.count_roots(sp.Rational(-1, 2), sp.Rational(-2, 5)) == 1, "negative root isolation")
require(quartic.count_roots(sp.Rational(13, 10), sp.Rational(4, 3)) == 1, "positive root isolation")
discriminant = sp.factor(sp.discriminant(debt_quartic, r))
require(discriminant < 0, "two nonreal exceptional roots")

# Work once in the irreducible quartic field.  Every complex root is an
# embedding of this field, so rank and constant-response statements transfer
# to all four exceptional specializations.
alpha = sp.CRootOf(debt_quartic, 0)
exceptional_field = sp.QQ.algebraic_field(alpha)
exceptional_rows = [
    [
        exceptional_field.from_sympy(sp.sympify(value).subs(r, alpha))
        for value in column
    ]
    for column in COLUMNS
]
exceptional_matrix = DomainMatrix(
    exceptional_rows,
    (len(exceptional_rows), len(exceptional_rows[0])),
    exceptional_field,
    fmt="dense",
)
require(exceptional_matrix.rank() == 26, "exceptional-field order-six rank")
exceptional_kernel = exceptional_matrix.nullspace(divide_last=True)
require(exceptional_kernel.shape == (4, 30), "exceptional-field nullity")
exceptional_kernel_rows = exceptional_kernel.to_list()
exceptional_responses = [
    sum((kernel_row[index] for index in constant_indices), exceptional_field.zero)
    for kernel_row in exceptional_kernel_rows
]
require(not any(exceptional_responses), "exceptional-field constant response vanishes")

payload = ";".join(f"{row}:{sp.sstr(value)}" for row, value in nonzero)
relation_sha256 = hashlib.sha256(payload.encode("ascii")).hexdigest()
quartic_sha256 = hashlib.sha256(sp.sstr(debt_quartic).encode("ascii")).hexdigest()

source = Path(__file__).read_text(encoding="utf-8")
require(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "assertion-free AST",
)

print("family=THM3677_zero_fourth_parabola;domain=QQ(r)")
print("order4=rows18_columns105_rank15_nullity3_constantinactive3")
print("order6=rows30_columns252_rank26_nullity4_activequotient1_nonzero23")
print(f"sixth_debt={expected_debt}")
print("exceptional_roots=irreducible_quartic;real=2;nonreal=2;isolating=(-1/2,-2/5),(13/10,4/3)")
print("exceptional_field=rank26_nullity4_constantactive0")
print(f"discriminant={discriminant}")
print(f"relation_sha256={relation_sha256}")
print(f"quartic_sha256={quartic_sha256}")
print(f"RESULT=PASS;gates={GATES}")
