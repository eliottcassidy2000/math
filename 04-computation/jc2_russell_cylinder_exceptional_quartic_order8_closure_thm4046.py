#!/usr/bin/env python3
"""Exact exceptional-quartic retained J8 identity and J7 shift probe.

Universe: all three coefficient slots of an arbitrary target two-form, with
all y^a z^b w^c coefficient monomials a+b+c <= 8.  The pullback density rows
are J0 through x-degree 4, J2 through degree 3, J4 through degree 2, J6
through degree 1, and the three J8 values.  Thus the complete retained matrix
has 45 rows and 3*binom(11,3)=495 columns.  No closedness or decomposability
condition is imposed.

The exceptional calculation is over K=QQ(alpha), where alpha obeys the
THM-3683 irreducible debt quartic.  Local coordinate jets are built through
total degree 9 before differentiating, so the order-eight top face is not
lost.  Exact K arithmetic uses the power basis (1,alpha,alpha^2,alpha^3).
"""

from __future__ import annotations

import ast
from hashlib import sha256
from math import comb
from pathlib import Path

from flint import fmpq, fmpq_mat


CHECKS = 0
PRIME = 137
CUTOFF = 8
COORDINATE_CUTOFF = 9
POINTS = (-1, 0, 1)
PARENT_SHA256 = (
    ("THM3673_script", "81d8237c9fefae07176d82bffb7d9b84763be3c13bb72755da2e0530d90db573"),
    ("THM3673_output", "e07c9dd7986e3c4e45b6bff4bbf6947043b3fd110d804bf93931d777408170ae"),
    ("THM3675_script", "543dd998a41924e16c1f6d1cf1439432fe050cef84e7c9301ecd25fb977bfd00"),
    ("THM3675_output", "699d1494ad63078014daec6e9efd03244a9f5afbd28435a0d6738728af5115f8"),
    ("THM3683_script", "b5e4132c322b3a01883688be9e8c993c5927a38f191c547eefcc84af432d9eb3"),
    ("THM3683_output", "3dfc849c88ca713c1a200fe05d26676624f9f83d9432e54a1c1c17058b320025"),
    ("THM3737_script", "4bbae46df140fbcda30b747f6356538a74fea34b928ed396e481e196fd0303c9"),
    ("THM3737_output", "e97bb6b47f147b5dac5fff8ceb4deab87e8cd85a902d4bb8a1397c1ccecbd996"),
    ("THM4043_script", "ef76a4aefe213c63ff4ce40d97fb57f6a9cf1b6ea90a7a4032b57c6e9c462de3"),
    ("THM4043_output", "7566341e1eff4664ddd19ae2332455b10375d7ae730d19da530097f2a50055a9"),
)


def require(condition, label):
    global CHECKS
    if not condition:
        raise RuntimeError(label)
    CHECKS += 1


QZERO = fmpq(0)
QONE = fmpq(1)
RELATION = (
    fmpq(1276420, 72783360),
    fmpq(-7849770, 72783360),
    fmpq(28419741, 72783360),
    fmpq(77822208, 72783360),
)
KZERO = (QZERO, QZERO, QZERO, QZERO)
KONE = (QONE, QZERO, QZERO, QZERO)
ALPHA = (QZERO, QONE, QZERO, QZERO)


def q(value, denominator=1):
    return fmpq(value, denominator)


def kconstant(value):
    if isinstance(value, fmpq):
        return (value, QZERO, QZERO, QZERO)
    return (q(value), QZERO, QZERO, QZERO)


def kadd(left, right):
    return tuple(a + b for a, b in zip(left, right))


def kneg(value):
    return tuple(-item for item in value)


def ksub(left, right):
    return tuple(a - b for a, b in zip(left, right))


def kmul(left, right):
    temporary = [QZERO] * 7
    for i, a in enumerate(left):
        if not a:
            continue
        for j, b in enumerate(right):
            if b:
                temporary[i + j] += a * b
    for degree in range(6, 3, -1):
        coefficient = temporary[degree]
        if coefficient:
            for offset, relation in enumerate(RELATION):
                temporary[degree - 4 + offset] += coefficient * relation
    return tuple(temporary[:4])


def kpow(value, exponent):
    answer = KONE
    base = value
    while exponent:
        if exponent & 1:
            answer = kmul(answer, base)
        base = kmul(base, base)
        exponent //= 2
    return answer


def kmul_q(value, scalar):
    scalar = scalar if isinstance(scalar, fmpq) else q(scalar)
    return tuple(item * scalar for item in value)


def multiplication_matrix(value):
    basis = (
        KONE,
        (QZERO, QONE, QZERO, QZERO),
        (QZERO, QZERO, QONE, QZERO),
        (QZERO, QZERO, QZERO, QONE),
    )
    columns = [kmul(value, item) for item in basis]
    return fmpq_mat(4, 4, [columns[column][row] for row in range(4) for column in range(4)])


def kinv(value):
    require(value != KZERO, "inverse of zero")
    matrix = multiplication_matrix(value)
    solution = matrix.solve(fmpq_mat(4, 1, [QONE, QZERO, QZERO, QZERO]))
    answer = tuple(solution[row, 0] for row in range(4))
    require(kmul(value, answer) == KONE, "field inverse residual")
    return answer


def ksum(values):
    answer = KZERO
    for value in values:
        answer = kadd(answer, value)
    return answer


def kdot(left, right):
    answer = KZERO
    for a, b in zip(left, right):
        if a != KZERO and b != KZERO:
            answer = kadd(answer, kmul(a, b))
    return answer


def poly_add(left, right):
    answer = dict(left)
    for key, value in right.items():
        answer[key] = kadd(answer.get(key, KZERO), value)
        if answer[key] == KZERO:
            del answer[key]
    return answer


def poly_neg(value):
    return {key: kneg(coefficient) for key, coefficient in value.items()}


def poly_sub(left, right):
    return poly_add(left, poly_neg(right))


def poly_scale(value, scalar):
    if not isinstance(scalar, tuple):
        scalar = kconstant(scalar)
    answer = {}
    for key, coefficient in value.items():
        product = kmul(coefficient, scalar)
        if product != KZERO:
            answer[key] = product
    return answer


def poly_mul(left, right, cutoff):
    answer = {}
    for (i, j), a in left.items():
        for (u, v), b in right.items():
            if i + j + u + v > cutoff:
                continue
            key = (i + u, j + v)
            answer[key] = kadd(answer.get(key, KZERO), kmul(a, b))
    return {key: value for key, value in answer.items() if value != KZERO}


def poly_power(value, exponent, cutoff):
    answer = {(0, 0): KONE}
    for _ in range(exponent):
        answer = poly_mul(answer, value, cutoff)
    return answer


def poly_diff(value, axis):
    answer = {}
    for (h_degree, t_degree), coefficient in value.items():
        exponent = (h_degree, t_degree)[axis]
        if not exponent:
            continue
        key = (h_degree - (axis == 0), t_degree - (axis == 1))
        answer[key] = kmul_q(coefficient, exponent)
    return answer


def poly_truncate(value, cutoff):
    return {key: coefficient for key, coefficient in value.items() if sum(key) <= cutoff}


def univariate_add(left, right):
    answer = dict(left)
    for degree, value in right.items():
        answer[degree] = kadd(answer.get(degree, KZERO), value)
        if answer[degree] == KZERO:
            del answer[degree]
    return answer


def univariate_scale(value, scalar):
    if not isinstance(scalar, tuple):
        scalar = kconstant(scalar)
    return {
        degree: product
        for degree, coefficient in value.items()
        if (product := kmul(coefficient, scalar)) != KZERO
    }


def univariate_shift_degree(value, amount):
    return {degree + amount: coefficient for degree, coefficient in value.items()}


def shifted_univariate(value, point, cutoff):
    answer = {}
    for degree, coefficient in value.items():
        for h_degree in range(min(degree, cutoff) + 1):
            multiplier = comb(degree, h_degree) * point ** (degree - h_degree)
            term = kmul_q(coefficient, multiplier)
            key = (h_degree, 0)
            answer[key] = kadd(answer.get(key, KZERO), term)
    return {key: coefficient for key, coefficient in answer.items() if coefficient != KZERO}


Q1 = {
    0: kconstant(q(-3, 4)),
    1: kconstant(1),
    2: kconstant(q(-27, 4)),
    3: kconstant(-2),
    4: kconstant(q(9, 2)),
    5: kconstant(1),
}
P = {2: KONE, 4: kconstant(-2), 6: KONE}
Q6 = univariate_add(Q1, univariate_scale(P, q(-259, 36)))
R1 = univariate_add(P, univariate_scale(univariate_shift_degree(P, 2), -1))
R2 = univariate_add(univariate_scale(P, 4), univariate_scale(univariate_shift_degree(P, 1), -9))


def parabola_parameter(r_value):
    return ksum(
        (
            kmul_q(kmul(r_value, r_value), q(520, 9)),
            kmul_q(r_value, q(-1688, 81)),
            kconstant(q(-5717, 729)),
        )
    )


def collision_coefficients(r_value, p_value):
    return univariate_add(
        Q6,
        univariate_add(univariate_scale(R1, p_value), univariate_scale(R2, r_value)),
    )


ROWS = tuple(
    (stable_degree, point, source_degree)
    for stable_degree, maximum_source_degree in (
        (0, 4),
        (2, 3),
        (4, 2),
        (6, 1),
        (8, 0),
    )
    for source_degree in range(maximum_source_degree + 1)
    for point in POINTS
)
require(len(ROWS) == 45, "45 retained rows")
J8_INDICES = tuple(index for index, row in enumerate(ROWS) if row[0] == 8)
LOWER_INDICES = tuple(index for index, row in enumerate(ROWS) if row[0] != 8)
CONSTANT_INDICES = tuple(
    index for index, row in enumerate(ROWS) if row[0] == 0 and row[2] == 0
)


def build_matrix(r_value, p_value):
    collision = collision_coefficients(r_value, p_value)
    packets = {}
    for point in POINTS:
        q_jet = shifted_univariate(collision, point, COORDINATE_CUTOFF)
        q_jet = poly_add(q_jet, {(0, 2): KONE})
        x_jet = {(0, 0): kconstant(point), (1, 0): KONE}
        x_squared_q = poly_mul(poly_mul(x_jet, x_jet, COORDINATE_CUTOFF), q_jet, COORDINATE_CUTOFF)
        D_jet = poly_add({(0, 0): KONE}, x_squared_q)
        D_plus_2 = poly_add(D_jet, {(0, 0): kconstant(2)})
        D_plus_3 = poly_add(D_jet, {(0, 0): kconstant(3)})
        y_jet = poly_scale(
            poly_mul(poly_mul(x_jet, D_jet, COORDINATE_CUTOFF), D_plus_2, COORDINATE_CUTOFF),
            q(1, 3),
        )
        z_jet = poly_add(
            poly_mul(q_jet, D_plus_3, COORDINATE_CUTOFF),
            {(0, 0): kconstant(3)},
        )
        require(y_jet.get((0, 0), KZERO) == KZERO, ("retained y zero", point))
        require(z_jet.get((0, 0), KZERO) == KZERO, ("retained z zero", point))
        y_h = poly_diff(y_jet, 0)
        y_t = poly_diff(y_jet, 1)
        z_h = poly_diff(z_jet, 0)
        z_t = poly_diff(z_jet, 1)
        area = poly_sub(
            poly_mul(y_h, z_t, CUTOFF),
            poly_mul(y_t, z_h, CUTOFF),
        )
        y_coefficients = poly_truncate(y_jet, CUTOFF)
        z_coefficients = poly_truncate(z_jet, CUTOFF)
        packets[point] = (
            [poly_power(y_coefficients, degree, CUTOFF) for degree in range(CUTOFF + 1)],
            [poly_power(z_coefficients, degree, CUTOFF) for degree in range(CUTOFF + 1)],
            (area, poly_truncate(y_h, CUTOFF), poly_truncate(z_h, CUTOFF)),
        )

    def monomial_pullback(kind, y_degree, z_degree, w_degree, point):
        y_powers, z_powers, bases = packets[point]
        value = poly_mul(y_powers[y_degree], z_powers[z_degree], CUTOFF)
        value = poly_mul(value, {(0, w_degree): KONE}, CUTOFF)
        return poly_mul(value, bases[kind], CUTOFF)

    matrix = []
    labels = []
    for kind in range(3):
        for y_degree in range(CUTOFF + 1):
            for z_degree in range(CUTOFF + 1 - y_degree):
                for w_degree in range(CUTOFF + 1 - y_degree - z_degree):
                    branches = {
                        point: monomial_pullback(kind, y_degree, z_degree, w_degree, point)
                        for point in POINTS
                    }
                    matrix.append(
                        [
                            branches[point].get((source_degree, stable_degree), KZERO)
                            for stable_degree, point, source_degree in ROWS
                        ]
                    )
                    labels.append((kind, y_degree, z_degree, w_degree))
    require(len(matrix) == 495, "495 target monomial columns")
    require(len(set(labels)) == 495, "495 distinct target monomial labels")
    return matrix, labels


def q_mod(value, prime=PRIME):
    return int(value.numerator) % prime * pow(int(value.denominator) % prime, -1, prime) % prime


def k_mod(value, root, prime=PRIME):
    answer = 0
    power = 1
    for coefficient in value:
        answer = (answer + q_mod(coefficient, prime) * power) % prime
        power = power * root % prime
    return answer


def matrix_mod(matrix, root, prime=PRIME):
    return [[k_mod(value, root, prime) for value in row] for row in matrix]


def mod_rref(matrix, prime=PRIME):
    rows = [[value % prime for value in row] for row in matrix]
    pivot_columns = []
    pivot_row = 0
    for column in range(len(rows[0])):
        selected = next(
            (index for index in range(pivot_row, len(rows)) if rows[index][column]),
            None,
        )
        if selected is None:
            continue
        rows[pivot_row], rows[selected] = rows[selected], rows[pivot_row]
        inverse = pow(rows[pivot_row][column], -1, prime)
        rows[pivot_row] = [value * inverse % prime for value in rows[pivot_row]]
        for index in range(len(rows)):
            if index == pivot_row:
                continue
            coefficient = rows[index][column]
            if coefficient:
                rows[index] = [
                    (value - coefficient * pivot) % prime
                    for value, pivot in zip(rows[index], rows[pivot_row])
                ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == len(rows):
            break
    return rows, pivot_columns


def mod_rank(matrix, prime=PRIME):
    return len(mod_rref(matrix, prime)[1])


def independent_row_indices_mod(matrix, prime=PRIME):
    transpose = [list(column) for column in zip(*matrix)]
    return mod_rref(transpose, prime)[1]


def k_rref_nullspace(matrix):
    rows = [list(row) for row in matrix]
    pivot_columns = []
    pivot_product = KONE
    pivot_row = 0
    for column in range(len(rows[0])):
        selected = next(
            (index for index in range(pivot_row, len(rows)) if rows[index][column] != KZERO),
            None,
        )
        if selected is None:
            continue
        rows[pivot_row], rows[selected] = rows[selected], rows[pivot_row]
        pivot_value = rows[pivot_row][column]
        pivot_product = kmul(pivot_product, pivot_value)
        inverse = kinv(pivot_value)
        rows[pivot_row] = [kmul(value, inverse) for value in rows[pivot_row]]
        for index in range(len(rows)):
            if index == pivot_row:
                continue
            coefficient = rows[index][column]
            if coefficient != KZERO:
                rows[index] = [
                    ksub(value, kmul(coefficient, pivot))
                    for value, pivot in zip(rows[index], rows[pivot_row])
                ]
        pivot_columns.append(column)
        pivot_row += 1
        if pivot_row == len(rows):
            break
    free_columns = [
        column for column in range(len(rows[0])) if column not in pivot_columns
    ]
    basis = []
    for free_column in free_columns:
        vector = [KZERO] * len(rows[0])
        vector[free_column] = KONE
        for index, pivot_column in enumerate(pivot_columns):
            vector[pivot_column] = kneg(rows[index][free_column])
        basis.append(vector)
    return rows, pivot_columns, free_columns, basis, pivot_product


def canonical_k(value):
    return "(" + ",".join(str(coordinate) for coordinate in value) + ")"


def canonical_relation(relation):
    return ";".join(
        f"J{stable}@{point}:d{source}={canonical_k(relation[index])}"
        for index, (stable, point, source) in enumerate(ROWS)
        if relation[index] != KZERO
    )


def alpha_polynomial(*coefficients):
    return ksum(
        kmul_q(kpow(ALPHA, degree), coefficient)
        for degree, coefficient in enumerate(coefficients)
    )


def expected_predecessor_blocks():
    """THM-3683's displayed R4, R2, R0, placed at J6, J4, J2."""
    expected = {}

    def put(stable, point, source, value):
        expected[(stable, point, source)] = value

    r4_value = kmul_q(alpha_polynomial(-503, 2340), q(8, 3159))
    put(6, -1, 0, r4_value)
    put(6, 0, 0, kneg(r4_value))
    put(6, -1, 1, kconstant(q(-5, 27)))
    put(6, 1, 1, kconstant(q(13, 27)))

    a2 = alpha_polynomial(14507690, -107434431, -310946688, 315394560)
    put(4, -1, 0, kmul_q(a2, q(-8, 9979281)))
    put(4, 0, 0, kmul_q(a2, q(8, 9979281)))
    put(4, -1, 1, kmul_q(alpha_polynomial(26159, -177840), q(1, 28431)))
    put(4, 0, 1, kmul_q(alpha_polynomial(1, -9), q(64, 81)))
    put(4, 1, 1, kmul_q(alpha_polynomial(65, 144), q(13, 2187)))
    put(4, -1, 2, kconstant(q(10, 81)))
    put(4, 1, 2, kconstant(q(26, 81)))

    a_minus = alpha_polynomial(
        -9091311692,
        156577328193,
        -1346654271456,
        -1919428213248,
        2755286876160,
    )
    a_zero = alpha_polynomial(
        35777404148,
        -119357786847,
        -347643535824,
        816178042368,
        196806205440,
    )
    b_minus = alpha_polynomial(374589757, -2408051880, -11632937472, 7569469440)
    put(2, -1, 0, kmul_q(a_minus, q(-16, 3502727631)))
    put(2, 0, 0, kmul_q(a_zero, q(16, 3502727631)))
    put(2, -1, 1, kmul_q(b_minus, q(2, 89813529)))
    put(2, 0, 1, kmul_q(alpha_polynomial(-514, 1233, 58968), q(-256, 85293)))
    put(2, 1, 1, kmul_q(alpha_polynomial(-10889, 182880, 331776), q(26, 531441)))
    put(2, -1, 2, kmul_q(alpha_polynomial(-14087, 121680), q(4, 85293)))
    put(2, 1, 2, kmul_q(alpha_polynomial(65, 144), q(52, 6561)))
    put(2, -1, 3, kconstant(q(-20, 243)))
    put(2, 1, 3, kconstant(q(52, 243)))
    return expected


def rational_matrix_rank(matrix, columns=None):
    if columns is None:
        columns = tuple(range(len(matrix[0])))
    flat = []
    for row in matrix:
        for index in columns:
            require(row[index][1:] == (QZERO, QZERO, QZERO), "rational specialization drift")
            flat.append(row[index][0])
    return fmpq_mat(len(matrix), len(columns), flat).rank()


print("== exceptional-quartic retained order-eight identity probe ==")
print("status=FINITE-EXACT arbitrary-target-two-forms no-closedness")
print("universe=45 retained rows x 495 complete target coefficient 8-jet columns")
print("coordinate_jet_cutoff=9 before differentiation; density_cutoff=8")
print(f"parent_sha256={PARENT_SHA256}")
print("logical_dependencies=THM3673,THM3675,THM3683,THM3737,THM4043")

# Exact exceptional matrix over K=QQ(alpha).
p_alpha = parabola_parameter(ALPHA)
exceptional_matrix, labels = build_matrix(ALPHA, p_alpha)
require((len(exceptional_matrix), len(exceptional_matrix[0])) == (495, 45), "matrix shape")

quartic_roots = tuple(
    value
    for value in range(PRIME)
    if (
        72783360 * value**4
        - 77822208 * value**3
        - 28419741 * value**2
        + 7849770 * value
        - 1276420
    )
    % PRIME
    == 0
)
require(quartic_roots == (44, 82, 92, 134), "split quartic roots mod 137")

root_ranks = []
root_lower_ranks = []
root_matrices = {}
for root in quartic_roots:
    reduction = matrix_mod(exceptional_matrix, root)
    root_matrices[root] = reduction
    root_ranks.append(mod_rank(reduction))
    root_lower_ranks.append(mod_rank([[row[index] for index in LOWER_INDICES] for row in reduction]))
require(tuple(root_ranks) == (40, 40, 40, 40), "four exceptional modular ranks")
require(tuple(root_lower_ranks) == (38, 38, 38, 38), "four exceptional lower-block ranks")
print(f"split_prime=137 roots={quartic_roots} ranks={tuple(root_ranks)} lower42_ranks={tuple(root_lower_ranks)}")

selected_indices = independent_row_indices_mod(root_matrices[44])
require(len(selected_indices) == 40, "40 modularly independent target rows")
selected = [exceptional_matrix[index] for index in selected_indices]
_rref, pivot_columns, free_columns, nullspace, pivot_product = k_rref_nullspace(selected)
require(len(pivot_columns) == 40, "exact rank-40 minor")
require(len(nullspace) == 5, "exact selected nullity five")
pivot_product_norm = multiplication_matrix(pivot_product).det()
require(pivot_product != KZERO and pivot_product_norm != 0, "nonzero exact pivot product")
require(
    tuple((index, ROWS[index]) for index in free_columns)
    == (
        (2, (0, 1, 0)),
        (17, (2, 1, 0)),
        (29, (4, 1, 0)),
        (38, (6, 1, 0)),
        (44, (8, 1, 0)),
    ),
    "canonical free-coordinate gauge",
)
for basis_index, relation in enumerate(nullspace):
    require(
        all(kdot(row, relation) == KZERO for row in exceptional_matrix),
        ("all-495 exact annihilation", basis_index),
    )
require(
    all(all(relation[index] == KZERO for index in J8_INDICES) for relation in nullspace[:4]),
    "four J8-inactive relations",
)
require(nullspace[4][J8_INDICES[0]] != KZERO, "one J8-active relation")

lambda_block = (kconstant(q(5, 18)), kconstant(-1), kconstant(q(13, 18)))
scale = kmul(lambda_block[0], kinv(nullspace[4][J8_INDICES[0]]))
active_relation = [kmul(value, scale) for value in nullspace[4]]
require(tuple(active_relation[index] for index in J8_INDICES) == lambda_block, "J8 block Lambda")
require(all(kdot(row, active_relation) == KZERO for row in exceptional_matrix), "active all-495 residual")
predecessor_blocks = expected_predecessor_blocks()
for index, row in enumerate(ROWS):
    if row[0] in (2, 4, 6):
        require(
            active_relation[index] == predecessor_blocks.get(row, KZERO),
            ("literal THM3683 predecessor block", row),
        )
support = tuple(index for index, value in enumerate(active_relation) if value != KZERO)
require(len(support) == 35, "35-entry active relation")
value_block_sums = {
    stable: ksum(
        active_relation[index]
        for index, row in enumerate(ROWS)
        if row[0] == stable and row[2] == 0
    )
    for stable in (0, 2, 4, 6, 8)
}
require(all(value_block_sums[stable] == KZERO for stable in (2, 4, 6, 8)), "positive value blocks kill constants")

kappa = value_block_sums[0]
require(kappa != KZERO, "nonzero exact constant response")
kappa_norm = multiplication_matrix(kappa).det()
require(kappa_norm != 0, "nonzero kappa norm")
kappa_inverse = kinv(kappa)
require(kmul(kappa, kappa_inverse) == KONE, "kappa inverse certificate")

canonical = canonical_relation(active_relation)
relation_hash = sha256(canonical.encode()).hexdigest()
minor_label_payload = ";".join(str(labels[index]) for index in selected_indices)
minor_labels_hash = sha256(minor_label_payload.encode()).hexdigest()
print("exceptional_exact_rank=40 nullity=5 J8_active_dimension=1")
print("canonical_gauge=free value rows at x=1 for J0,J2,J4,J6,J8; first four set zero")
print("J8_block=((5/18,0,0,0),(-1,0,0,0),(13/18,0,0,0))")
print("lower_blocks_literal_THM3683_shift=J6:R4,J4:R2,J2:R0")
print(f"active_support_count={len(support)}")
print(f"relation_canonical={canonical}")
print(f"relation_sha256={relation_hash}")
print(f"rank40_minor_target_labels_sha256={minor_labels_hash}")
print(f"rank40_pivot_product_coordinates={canonical_k(pivot_product)}")
print(f"rank40_pivot_product_norm={pivot_product_norm}")
print(f"kappa_coordinates={canonical_k(kappa)}")
print(f"kappa_norm={kappa_norm}")
print(f"kappa_inverse_coordinates={canonical_k(kappa_inverse)}")
print("positive_stable_value_block_sums=J2:0,J4:0,J6:0,J8:0")

kappa_reductions = tuple(k_mod(kappa, root) for root in quartic_roots)
require(kappa_reductions == (105, 71, 89, 8), "four nonzero modular constant responses")
relation_reduction_hashes = []
for root in quartic_roots:
    reduced_support = [
        (ROWS[index], k_mod(value, root))
        for index, value in enumerate(active_relation)
        if k_mod(value, root)
    ]
    relation_reduction_hashes.append(
        sha256(repr(reduced_support).encode()).hexdigest()
    )
require(all(value for value in kappa_reductions), "four nonzero split responses")
print(f"kappa_mod137_at_roots={kappa_reductions}")
print(f"relation_mod137_hashes={tuple(relation_reduction_hashes)}")

# Exact hostile controls.  The zero-fourth parabola away from the exceptional
# quartic loses the J8-active class.  Leaving the parabola raises rank again.
r_zero = KZERO
parabola_r0_matrix, _ = build_matrix(r_zero, parabola_parameter(r_zero))
parabola_r0_rank = rational_matrix_rank(parabola_r0_matrix)
parabola_r0_lower_rank = rational_matrix_rank(parabola_r0_matrix, LOWER_INDICES)
require((parabola_r0_rank, parabola_r0_lower_rank) == (41, 38), "exact off-quartic parabola hostile")
require(3 - parabola_r0_rank + parabola_r0_lower_rank == 0, "off-quartic J8 projection zero")

q6_matrix, _ = build_matrix(r_zero, KZERO)
q6_rank = rational_matrix_rank(q6_matrix)
q6_lower_rank = rational_matrix_rank(q6_matrix, LOWER_INDICES)
require((q6_rank, q6_lower_rank) == (42, 39), "exact off-parabola Q6 hostile")
require(3 - q6_rank + q6_lower_rank == 0, "Q6 J8 projection zero")

off_parabola_alpha_matrix, _ = build_matrix(ALPHA, kadd(p_alpha, KONE))
off_parabola_mod44 = matrix_mod(off_parabola_alpha_matrix, 44)
off_parabola_mod44_rank = mod_rank(off_parabola_mod44)
off_parabola_mod44_lower_rank = mod_rank(
    [[row[index] for index in LOWER_INDICES] for row in off_parabola_mod44]
)
require((off_parabola_mod44_rank, off_parabola_mod44_lower_rank) == (42, 39), "off-parabola alpha hostile")
print("hostile_parabola_r0_exact rank=41 nullity=4 lower42_rank=38 J8_active_dimension=0")
print("hostile_Q6_exact rank=42 nullity=3 lower42_rank=39 J8_active_dimension=0")
print("hostile_alpha_pplus1_mod137_root44 rank=42 nullity=3 lower42_rank=39 J8_active_dimension=0")

# Structural consequences.  The exact active relation prolongs THM-3683 as
#   Lambda(J8)+R4(J6)+R2(J4)+R0(J2)+S0(J0)=0.
# Applying it to w*omega uses w=t and gives Jtilde_n=J_(n-1), hence
#   Lambda(J7)+R4(J5)+R2(J3)+R0(J1)=0.
# If an inherited stagewise lift has J1=...=J6=0, its provisional seventh
# debt D7 obeys Lambda(D7)=0.  THM-3737 then solves
#   D7+8 L0(F8,G8)=0.
# With J7 killed, the unshifted identity gives Lambda(D8)=-kappa.  The next
# correction is 9 L0(F9,G9), whose Lambda is zero, so J8 cannot be killed.
print("identity=Lambda(J8)+THM3683_R4(J6)+THM3683_R2(J4)+THM3683_R0(J2)+S0(J0)=0")
print("w_shift=Lambda(J7)+THM3683_R4(J5)+THM3683_R2(J3)+THM3683_R0(J1)=0")
print("stage7_if_J1_through_J6_zero: Lambda(D7)=0; solve D7+8*L0(F8,G8)=0 via THM3737")
print("stage8_if_J2_J4_J6_zero_and_J0_one: Lambda(D8)=-kappa!=0")
print("stage8_correction=D8+9*L0(F9,G9); Lambda(L0)=0; no J8 extension")
print("character_dilation_H=eta*t^k,k>=2: blocks J0,Jk,J2k,J3k,J4k weights eta^4,eta^3,eta^2,eta,1")
print("formal_conjugacy_invoice_for_constant_c: 0=eta^4*c*kappa; positive stable constant-series blocks vanish")
print("CONSEQUENCE via THM3673/3675 mechanism: obstruction transfers to every 0!=H in t^2*C[t]")
print("SCOPE retained finite-jet obstruction for this Russell compiler; no global JC2 claim")

source_path = Path(__file__).resolve()
source_bytes = source_path.read_bytes()
require(b"\r\n" not in source_bytes, "raw LF source")
tree = ast.parse(source_bytes.decode("utf-8"))
require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "Python assert node")
print(f"script_sha256={sha256(source_bytes).hexdigest()}")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
