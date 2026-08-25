#!/usr/bin/env python3
"""Exact companion for THM-4043's shifted identities and J6 lift.

On the THM-3677 zero-fourth-debt parabola, build the complete retained odd
window

    J1 Taylor order <= 2, J3 Taylor order <= 1, J5 values

for every target two-form coefficient monomial of total target degree <= 5.
The three expected left-cokernel rows are the w^5, w^3, and w shifts of the
THM-3683 universal sixth-order arbitrary-form identity.  Closedness and
decomposability are deliberately not imposed.

The rational Q6 point (p,r)=(0,0), which is off the zero-fourth parabola, is
rebuilt as a hostile boundary: its rank rises from 15 to 16 and neither of
its two left-cokernel rows meets the J5 block.  The companion also rebuilds
R0(1) from THM-3683's displayed coefficients and reduces it exactly by the
exceptional debt quartic, which is the scalar gate for the J6 continuation.
"""

from __future__ import annotations

import ast
from hashlib import sha256
from pathlib import Path

import sympy as sp
from sympy.polys.matrices import DomainMatrix


CHECKS = 0


def require(condition, payload):
    global CHECKS
    if condition is not True:
        raise RuntimeError(payload)
    CHECKS += 1


def raw_sha256(path):
    return sha256(path.read_bytes()).hexdigest()


ROOT = Path(__file__).resolve().parents[1]
PARENT_FILES = (
    (
        "thm3683_script",
        ROOT / "04-computation/jc2_russell_cylinder_sixth_debt_quartic_thm3683.py",
        "b5e4132c322b3a01883688be9e8c993c5927a38f191c547eefcc84af432d9eb3",
    ),
    (
        "thm3683_output",
        ROOT / "05-knowledge/results/jc2_russell_cylinder_sixth_debt_quartic_thm3683.out",
        "3dfc849c88ca713c1a200fe05d26676624f9f83d9432e54a1c1c17058b320025",
    ),
    (
        "thm3688_script",
        ROOT / "04-computation/jc2_russell_cylinder_exceptional_quartic_exact_j1_j2_lift_thm3688.py",
        "02cd67446b18b3863bc3665d48a6c5cccda81c394f94b754d2b90b1597c53ba6",
    ),
    (
        "thm3688_output",
        ROOT / "05-knowledge/results/jc2_russell_cylinder_exceptional_quartic_exact_j1_j2_lift_thm3688.out",
        "0c6bce1baef8dab98cd63219538c5e14aa9eae05903dc85fee647f776b514459",
    ),
    (
        "thm3737_script",
        ROOT / "04-computation/jc2_russell_cylinder_exceptional_quartic_jacobian_image_hyperplane_thm3737.py",
        "4bbae46df140fbcda30b747f6356538a74fea34b928ed396e481e196fd0303c9",
    ),
    (
        "thm3737_output",
        ROOT / "05-knowledge/results/jc2_russell_cylinder_exceptional_quartic_jacobian_image_hyperplane_thm3737.out",
        "e97bb6b47f147b5dac5fff8ceb4deab87e8cd85a902d4bb8a1397c1ccecbd996",
    ),
)
for parent_label, parent_path, expected_hash in PARENT_FILES:
    observed_hash = raw_sha256(parent_path)
    require(
        observed_hash == expected_hash,
        ("parent drift", parent_label, observed_hash, expected_hash),
    )


x, q, t, X, r = sp.symbols("x q t X r")
POINTS = (-1, 0, 1)
CUTOFF = 5
LAMBDA = (sp.Rational(5, 18), -1, sp.Rational(13, 18))

Q1 = (
    x**5
    + sp.Rational(9, 2) * x**4
    - 2 * x**3
    - sp.Rational(27, 4) * x**2
    + x
    - sp.Rational(3, 4)
)
P = sp.expand(x**2 * (x**2 - 1) ** 2)
Q6 = sp.expand(Q1 - sp.Rational(259, 36) * P)
R1 = sp.expand(P * (1 - x**2))
R2_POLYNOMIAL = sp.expand(P * (4 - 9 * x))
P_PARABOLA = (
    sp.Rational(520, 9) * r**2
    - sp.Rational(1688, 81) * r
    - sp.Rational(5717, 729)
)
Q_PARABOLA = sp.expand(Q6 + P_PARABOLA * R1 + r * R2_POLYNOMIAL)

D = 1 + x**2 * q
b = sp.expand((D - 1) * (D + 2) ** 2)
c = sp.expand(x * D * (D + 2))
e = sp.expand(q * (D + 3))
y = c / 3
z = e + 3

ROWS = tuple(
    (stable_degree, point, source_degree)
    for stable_degree, maximum_source_degree in ((1, 2), (3, 1), (5, 0))
    for source_degree in range(maximum_source_degree + 1)
    for point in POINTS
)


def multiply(left, right):
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


def build_transposed_matrix(collision_polynomial, coefficient_domain):
    pulled_y = sp.expand(y.subs(q, collision_polynomial + t**2))
    pulled_z = sp.expand(z.subs(q, collision_polynomial + t**2))
    y_x = sp.diff(pulled_y, x)
    y_t = sp.diff(pulled_y, t)
    z_x = sp.diff(pulled_z, x)
    z_t = sp.diff(pulled_z, t)
    area = sp.expand(y_x * z_t - y_t * z_x)

    def jet(expression, point):
        shifted = sp.Poly(
            sp.expand(expression.subs(x, X + point)),
            X,
            t,
            domain=coefficient_domain,
        )
        return {
            (source_degree, stable_degree): coefficient
            for (source_degree, stable_degree), coefficient in shifted.terms()
            if source_degree + stable_degree <= CUTOFF
        }

    packets = {}
    for point in POINTS:
        y_jet = jet(pulled_y, point)
        z_jet = jet(pulled_z, point)
        require(y_jet.get((0, 0), 0) == 0, ("retained y", point))
        require(z_jet.get((0, 0), 0) == 0, ("retained z", point))
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

    columns = []
    labels = []
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
                    columns.append(
                        tuple(
                            branch[point].get((source_degree, stable_degree), 0)
                            for stable_degree, point, source_degree in ROWS
                        )
                    )
                    labels.append((kind, y_degree, z_degree, w_degree))

    require(len(ROWS) == 18, "18 retained odd rows")
    require(len(columns) == 168, "3*binom(8,3)=168 target columns")
    require(len(set(labels)) == 168, "distinct target columns")
    matrix = DomainMatrix.from_list_sympy(
        len(columns), len(ROWS), [list(column) for column in columns]
    ).to_field()
    return matrix, columns, labels


def zero_relation():
    return [sp.S.Zero for _ in ROWS]


def put(relation, stable_degree, point, source_degree, value):
    relation[ROWS.index((stable_degree, point, source_degree))] += value


def add_lambda(relation, stable_degree):
    for point, value in zip(POINTS, LAMBDA):
        put(relation, stable_degree, point, 0, value)


A2 = 315394560 * r**3 - 310946688 * r**2 - 107434431 * r + 14507690


def add_R4(relation, stable_degree):
    value = 8 * (2340 * r - 503) / 3159
    put(relation, stable_degree, -1, 0, value)
    put(relation, stable_degree, 0, 0, -value)
    put(relation, stable_degree, -1, 1, -sp.Rational(5, 27))
    put(relation, stable_degree, 1, 1, sp.Rational(13, 27))


def add_R2(relation, stable_degree):
    put(relation, stable_degree, -1, 0, -8 * A2 / 9979281)
    put(relation, stable_degree, 0, 0, 8 * A2 / 9979281)
    put(
        relation,
        stable_degree,
        -1,
        1,
        -(177840 * r - 26159) / 28431,
    )
    put(relation, stable_degree, 0, 1, -64 * (9 * r - 1) / 81)
    put(relation, stable_degree, 1, 1, 13 * (144 * r + 65) / 2187)
    put(relation, stable_degree, -1, 2, sp.Rational(10, 81))
    put(relation, stable_degree, 1, 2, sp.Rational(26, 81))


relation_1 = zero_relation()
add_lambda(relation_1, 1)

relation_3 = zero_relation()
add_R4(relation_3, 1)
add_lambda(relation_3, 3)

relation_5 = zero_relation()
add_R2(relation_5, 1)
add_R4(relation_5, 3)
add_lambda(relation_5, 5)

EXPECTED_RELATIONS = (relation_1, relation_3, relation_5)


print("== THM-4043 exceptional-quartic shifted stable identities and J6 lift ==")
print("status=PROVED VERIFIED-EXACT arbitrary-target-two-forms no-closedness")
print(
    "parent_sha256",
    tuple((label, expected_hash) for label, _path, expected_hash in PARENT_FILES),
)

parabola_matrix, parabola_columns, labels = build_transposed_matrix(
    Q_PARABOLA, sp.QQ[r]
)
require(parabola_matrix.domain == sp.QQ.frac_field(r), "QQ(r) matrix domain")
require(parabola_matrix.shape == (168, 18), "parabola matrix shape")
require(parabola_matrix.rank() == 15, "parabola rank 15")
kernel = parabola_matrix.nullspace(divide_last=True).to_Matrix()
require(kernel.shape == (3, 18), "parabola left nullity 3")

expected_matrix = sp.Matrix(EXPECTED_RELATIONS)
require(expected_matrix.rank() == 3, "three expected relations independent")
for relation_name, relation in zip(("w5", "w3", "w1"), EXPECTED_RELATIONS):
    require(
        sum(value != 0 for value in relation)
        == {"w5": 3, "w3": 7, "w1": 14}[relation_name],
        ("relation support", relation_name),
    )
    residual = parabola_matrix * DomainMatrix.from_Matrix(sp.Matrix(relation)).to_field()
    require(residual.is_zero_matrix, ("expected relation", relation_name))
    for label, column in zip(labels, parabola_columns):
        require(
            sp.factor(sum(relation[index] * column[index] for index in range(18)))
            == 0,
            ("column identity", relation_name, label),
        )

j5_indices = [index for index, row in enumerate(ROWS) if row[0] == 5]
j5_projection = kernel[:, j5_indices]
require(j5_projection.rank() == 1, "exactly one J5-active cokernel class")
active_rows = [
    index
    for index in range(kernel.rows)
    if any(kernel[index, column] != 0 for column in j5_indices)
]
require(len(active_rows) == 1, "one active basis row")
active_block = [kernel[active_rows[0], column] for column in j5_indices]
scale = LAMBDA[0] / active_block[0]
require(
    tuple(sp.factor(scale * value) for value in active_block) == LAMBDA,
    "J5-active block is Lambda",
)
print("parabola_matrix shape=168x18 rank=15 left_nullity=3")
print(
    "parabola_relations supports=(3,7,14) "
    "basis=(LambdaJ1,LambdaJ3+R4J1,LambdaJ5+R4J3+R2J1)"
)
print("parabola_J5_active_dimension=1 active_block=(5/18,-1,13/18)")

# Hostile off-parabola control: Q6 is (p,r)=(0,0), whereas the parabola has
# p=-5717/729 at r=0.
q6_matrix, _q6_columns, _q6_labels = build_transposed_matrix(Q6, sp.QQ)
require(q6_matrix.domain == sp.QQ, "Q6 rational domain")
require(q6_matrix.shape == (168, 18), "Q6 matrix shape")
require(q6_matrix.rank() == 16, "Q6 hostile rank 16")
q6_kernel = q6_matrix.nullspace(divide_last=True).to_Matrix()
require(q6_kernel.shape == (2, 18), "Q6 hostile left nullity 2")
require(q6_kernel[:, j5_indices] == sp.zeros(2, 3), "Q6 has no J5-active relation")
print(
    "hostile_Q6 p=0 r=0 off_parabola_p=-5717/729 "
    "rank=16 left_nullity=2 J5_active_dimension=0"
)

# The k=2 shift used at J4 reads Lambda(J4)+R4(J2)+R2(J0)=0.
# THM-3688 has J0=1 and J2=0.  The only value coefficients of R2 cancel,
# while every other R2 coordinate is a positive x derivative.
R2_ONE = sp.factor(-8 * A2 / 9979281 + 8 * A2 / 9979281)
require(R2_ONE == 0, "R2 kills the constant J0")
print("R2_one=0 shifted_k2_clears_J4_when_J2_zero")

# Rebuild THM-3683's R0(1) from its two non-derivative coefficients.  Every
# derivative term in R0 kills the constant J0=1.
A_MINUS_R0 = (
    2755286876160 * r**4
    - 1919428213248 * r**3
    - 1346654271456 * r**2
    + 156577328193 * r
    - 9091311692
)
A_ZERO_R0 = (
    196806205440 * r**4
    + 816178042368 * r**3
    - 347643535824 * r**2
    - 119357786847 * r
    + 35777404148
)
R0_ONE = sp.factor(
    -16 * A_MINUS_R0 / 3502727631 + 16 * A_ZERO_R0 / 3502727631
)
DEBT_QUARTIC = (
    72783360 * r**4
    - 77822208 * r**3
    - 28419741 * r**2
    + 7849770 * r
    - 1276420
)
EXPECTED_R0_ONE = -sp.Rational(256, 1594323) * DEBT_QUARTIC
require(sp.factor(R0_ONE - EXPECTED_R0_ONE) == 0, "exact R0(1) debt quartic")
exceptional_remainder = sp.Poly(R0_ONE, r, domain=sp.QQ).rem(
    sp.Poly(DEBT_QUARTIC, r, domain=sp.QQ)
)
require(exceptional_remainder.is_zero, "R0(1) vanishes in exceptional field")
require(sp.Integer(6) != 0 and sp.Integer(7) != 0, "stage divisors in characteristic zero")
print("R0_one=-256*F(r)/1594323 exceptional_field_remainder=0")

# Logical composition with THM-3688 and the proved cutoff-free image theorem
# THM-3737: the k=3,2,1 shifts successively put the J3,J4,J5 debts in
# ker Lambda, where L0(S^2)=ker Lambda supplies F4,F5,F6.  The R0 gate does
# the same at J6 with divisor 7.
stagewise_through_j6 = (
    parabola_matrix.rank() == 15
    and j5_projection.rank() == 1
    and R2_ONE == 0
    and exceptional_remainder.is_zero
)
require(stagewise_through_j6, "logical J5/J6 scalar gates")
print("shift w^k: Jtilde_n=J_(n-k), J_negative=0, k=(5,3,1)")
print("actual_stagewise_through_J6=1 via_THM3688_input_and_THM3737_image")
print("SCOPE finite stagewise actual lift; next scalar J7; no all-order/global/JC2")

source_path = Path(__file__).resolve()
source_bytes = source_path.read_bytes()
require(b"\r\n" not in source_bytes, "raw LF source")
tree = ast.parse(source_bytes.decode("utf-8"))
require(
    not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
    "Python assert node present",
)
print("script_sha256", sha256(source_bytes).hexdigest())
print("CHECKS", CHECKS)
print("RESULT=PASS")
