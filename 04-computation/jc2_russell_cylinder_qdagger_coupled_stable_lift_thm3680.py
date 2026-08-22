"""Exact Q_dagger target-ring lift through J1=J2=0 (THM-3680)."""

import contextlib
import gc
import hashlib
import io
import runpy
import sys
from pathlib import Path

import sympy as sp
from flint import fmpq, fmpq_mat, nmod_mat, nmod_poly


sys.set_int_max_str_digits(0)
MODULUS = 1000003
BASE_PATH = Path(
    "04-computation/jc2_russell_cylinder_qdagger_actual_j0_lift_thm3678.py"
)
BASE_SHA256 = "63f59e111e298652073bfeb7187eebaeb99c8d5e8def6efbdb5420b099d92d46"
EXPECTED = {
    "pivot_columns": "0aae92931398ae5428dd9d229bc804a46ea838ef1312d7fc5d139152b2635c37",
    "pivot_rows": "6f65529dfc366c2fce2de7dbbd54f86d37b07de6021a9372e749f718e349315f",
    "F2_target": "ab540e283912a0965a521f9e0bfe5a4b1c6b3074c15d88e22c4518b5160d144e",
    "F2_restriction": "f0cddc4c6fd0b94fbb3f02413fb9fa5f1ae74efc902406b873c671d6598689ff",
    "G2_target": "340a5c25b2a3eac6d44d915ea5db7c0effa359e7cc8bf0079928de8e1d7d7c66",
    "G2_restriction": "5c46e799104a50ef818ed125bdc4fbb830f8dddbb43f4baab500a34d34df1536",
    "F3_target": "e51511c60ce30902a5a16937ffab87ca9d90f4d8f3acc5488b214360c92eb750",
    "F3_restriction": "b4eb6e645727231e05bc031bf2e697574915a36fa31644e6fb0122b73b29b731",
    "G3_target": "b42158e37532fcbea6f074be62b8cfe97acab0ec6397f8f41ac18229437479ae",
    "G3_restriction": "bf1ebe51fbd16660623d3977e731e4cc021f50d79f628a09cad77e99847164f3",
}

gates = 0


def require(condition, label):
    global gates
    if not condition:
        raise RuntimeError(label)
    gates += 1


require(BASE_PATH.is_file(), "THM-3678 companion missing")
require(hashlib.sha256(BASE_PATH.read_bytes()).hexdigest() == BASE_SHA256, "THM-3678 companion hash")
base_transcript = io.StringIO()
with contextlib.redirect_stdout(base_transcript):
    base = runpy.run_path(str(BASE_PATH))
require(base_transcript.getvalue().endswith("RESULT=PASS;gates=238\n"), "THM-3678 exact replay")

x = base["x"]
Bq, Cq, Eq = base["B"], base["C"], base["E"]
F1q, G1q = base["F"], base["G"]
delta_F1q, delta_G1q = base["delta_F"], base["delta_G"]
delta_cq, delta_eq = base["dC"], base["dE"]
WEIGHTS = base["WEIGHTS"]


def modular_rational(value):
    value = sp.Rational(value)
    return int(value.p) * pow(int(value.q), -1, MODULUS) % MODULUS


def to_nmod(polynomial):
    polynomial = sp.Poly(polynomial, x, domain=sp.QQ)
    if polynomial.is_zero:
        return nmod_poly([], MODULUS)
    coefficients = [0] * (polynomial.degree() + 1)
    for (degree,), coefficient in polynomial.terms():
        coefficients[degree] = modular_rational(coefficient)
    return nmod_poly(coefficients, MODULUS)


B, C, E = (to_nmod(polynomial) for polynomial in (Bq, Cq, Eq))
F1, G1 = (to_nmod(polynomial) for polynomial in (F1q, G1q))
delta_F1, delta_G1 = (to_nmod(polynomial) for polynomial in (delta_F1q, delta_G1q))
delta_c, delta_e = (to_nmod(polynomial) for polynomial in (delta_cq, delta_eq))


def modular_packet(cutoff):
    generators = (B, C, E)
    powers = []
    for generator, weight in zip(generators, WEIGHTS):
        row = [nmod_poly([1], MODULUS)]
        for _ in range(cutoff // weight):
            row.append(row[-1] * generator)
        powers.append(row)
    answer = []
    for i in range(cutoff // WEIGHTS[0] + 1):
        for j in range(cutoff // WEIGHTS[1] + 1):
            for k in range(cutoff // WEIGHTS[2] + 1):
                if i * WEIGHTS[0] + j * WEIGHTS[1] + k * WEIGHTS[2] <= cutoff:
                    answer.append(powers[0][i] * powers[1][j] * powers[2][k])
    return answer


known_j1 = (
    2 * C.derivative() * delta_e
    + F1.derivative() * G1
    - F1 * G1.derivative()
    - 2 * delta_c * E.derivative()
)
known_j2 = (
    3 * C.derivative() * delta_G1
    + 2 * F1.derivative() * delta_e
    + delta_c.derivative() * G1
    - F1 * delta_e.derivative()
    - 2 * delta_c * G1.derivative()
    - 3 * delta_F1 * E.derivative()
)
require((known_j1.degree(), known_j2.degree()) == (386, 204), "known-term degrees")


def mod_coefficient(polynomial, degree):
    return int(polynomial[degree]) if degree <= polynomial.degree() else 0


def modular_system(cutoff):
    monomials = modular_packet(cutoff)
    zero = nmod_poly([], MODULUS)
    columns = []
    for polynomial in monomials:
        columns.append(
            (
                -2 * polynomial * E.derivative(),
                polynomial.derivative() * G1 - 2 * polynomial * G1.derivative(),
            )
        )
    for polynomial in monomials:
        columns.append(
            (
                2 * C.derivative() * polynomial,
                2 * F1.derivative() * polynomial - F1 * polynomial.derivative(),
            )
        )
    for polynomial in monomials:
        columns.append((zero, -3 * polynomial * E.derivative()))
    for polynomial in monomials:
        columns.append((zero, 3 * C.derivative() * polynomial))

    rows_j1 = max([known_j1.degree(), *(pair[0].degree() for pair in columns)]) + 1
    rows_j2 = max([known_j2.degree(), *(pair[1].degree() for pair in columns)]) + 1
    operator_data = []
    augmented_data = []
    for row in range(rows_j1):
        values = [mod_coefficient(pair[0], row) for pair in columns]
        operator_data.extend(values)
        augmented_data.extend([*values, -mod_coefficient(known_j1, row) % MODULUS])
    for row in range(rows_j2):
        values = [mod_coefficient(pair[1], row) for pair in columns]
        operator_data.extend(values)
        augmented_data.extend([*values, -mod_coefficient(known_j2, row) % MODULUS])
    operator = nmod_mat(rows_j1 + rows_j2, len(columns), operator_data, MODULUS)
    augmented = nmod_mat(
        rows_j1 + rows_j2, len(columns) + 1, augmented_data, MODULUS
    )
    return monomials, columns, rows_j1, rows_j2, operator, augmented


# A nearby hostile control.  This is only a modular cutoff signal, not a
# characteristic-zero impossibility theorem.
lower_packet, _lower_columns, lower_rows_j1, lower_rows_j2, lower_operator, lower_augmented = modular_system(366)
lower_operator_rank = lower_operator.rank()
lower_augmented_rank = lower_augmented.rank()
require(len(lower_packet) == 953, "cutoff-366 packet")
require((lower_rows_j1, lower_rows_j2, lower_operator.ncols()) == (387, 561, 3812), "cutoff-366 shape")
require((lower_operator_rank, lower_augmented_rank) == (941, 942), "cutoff-366 modular ranks")
del lower_packet, _lower_columns, lower_operator, lower_augmented
gc.collect()

# The feasible modular system selects a deterministic characteristic-zero
# square.  Modular rank itself is not asserted as rational full-operator rank.
packet_mod, columns_mod, rows_j1, rows_j2, operator_mod, augmented_mod = modular_system(369)
reduced_augmented, augmented_rank = augmented_mod.rref()
operator_rank = operator_mod.rank()
require(len(packet_mod) == 973, "cutoff-369 packet")
require((rows_j1, rows_j2, operator_mod.ncols()) == (390, 564, 3892), "cutoff-369 shape")
require((operator_rank, augmented_rank) == (947, 947), "cutoff-369 modular ranks")

pivot_columns = []
for row in range(augmented_rank):
    pivot = next(
        (column for column in range(len(columns_mod) + 1) if reduced_augmented[row, column]),
        None,
    )
    require(pivot is not None and pivot < len(columns_mod), f"modular pivot row {row}")
    pivot_columns.append(pivot)
require(len(pivot_columns) == 947, "pivot column count")

selected_data = []
for row in range(rows_j1 + rows_j2):
    selected_data.extend(int(operator_mod[row, column]) for column in pivot_columns)
selected_mod = nmod_mat(rows_j1 + rows_j2, len(pivot_columns), selected_data, MODULUS)
transpose_reduced, transpose_rank = selected_mod.transpose().rref()
require(transpose_rank == 947, "selected modular column independence")
pivot_rows = []
for row in range(transpose_rank):
    pivot_rows.append(
        next(
            column
            for column in range(rows_j1 + rows_j2)
            if transpose_reduced[row, column]
        )
    )
require(len(pivot_rows) == 947, "pivot row count")


def integer_hash(values):
    return hashlib.sha256(",".join(str(value) for value in values).encode("ascii")).hexdigest()


pivot_column_hash = integer_hash(pivot_columns)
pivot_row_hash = integer_hash(pivot_rows)
require(pivot_column_hash == EXPECTED["pivot_columns"], "pivot column hash")
require(pivot_row_hash == EXPECTED["pivot_rows"], "pivot row hash")
del packet_mod, columns_mod, operator_mod, augmented_mod, reduced_augmented, selected_mod, transpose_reduced
gc.collect()

# Rebuild the selected 947-square over Q.
monomials = base["packet"](369)
restrictions = [item[1] for item in monomials]
require(len(monomials) == 973, "exact cutoff-369 packet")
zero_q = sp.Poly(0, x, domain=sp.QQ)


def exact_column(column):
    block, index = divmod(column, len(restrictions))
    polynomial = restrictions[index]
    if block == 0:
        return (
            -2 * polynomial * Eq.diff(),
            polynomial.diff() * G1q - 2 * polynomial * G1q.diff(),
        )
    if block == 1:
        return (
            2 * Cq.diff() * polynomial,
            2 * F1q.diff() * polynomial - F1q * polynomial.diff(),
        )
    if block == 2:
        return zero_q, -3 * polynomial * Eq.diff()
    if block == 3:
        return zero_q, 3 * Cq.diff() * polynomial
    raise RuntimeError(f"bad block {block}")


selected_exact_columns = [exact_column(column) for column in pivot_columns]
known_j1_q = (
    2 * Cq.diff() * delta_eq
    + F1q.diff() * G1q
    - F1q * G1q.diff()
    - 2 * delta_cq * Eq.diff()
)
known_j2_q = (
    3 * Cq.diff() * delta_G1q
    + 2 * F1q.diff() * delta_eq
    + delta_cq.diff() * G1q
    - F1q * delta_eq.diff()
    - 2 * delta_cq * G1q.diff()
    - 3 * delta_F1q * Eq.diff()
)


def as_fmpq(value):
    value = sp.Rational(value)
    return fmpq(int(value.p), int(value.q))


def as_rational(value):
    return sp.Rational(int(value.numerator), int(value.denominator))


square_data = []
rhs_data = []
for row in pivot_rows:
    if row < rows_j1:
        square_data.extend(as_fmpq(pair[0].nth(row)) for pair in selected_exact_columns)
        rhs_data.append(as_fmpq(-known_j1_q.nth(row)))
    else:
        degree = row - rows_j1
        square_data.extend(as_fmpq(pair[1].nth(degree)) for pair in selected_exact_columns)
        rhs_data.append(as_fmpq(-known_j2_q.nth(degree)))

square = fmpq_mat(947, 947, square_data)
rhs = fmpq_mat(947, 1, rhs_data)
solution = square.solve(rhs)
require(square * solution == rhs, "selected exact square residual")

target_vectors = [[fmpq(0)] * len(restrictions) for _ in range(4)]
lift_restrictions = [sp.Poly(0, x, domain=sp.QQ) for _ in range(4)]
for local_index, column in enumerate(pivot_columns):
    coefficient = solution[local_index, 0]
    if not coefficient:
        continue
    block, index = divmod(column, len(restrictions))
    target_vectors[block][index] = coefficient
    lift_restrictions[block] += as_rational(coefficient) * restrictions[index]

F2q, G2q, F3q, G3q = lift_restrictions
j1_residual = known_j1_q - 2 * F2q * Eq.diff() + 2 * Cq.diff() * G2q
j2_residual = (
    known_j2_q
    + F2q.diff() * G1q
    - 2 * F2q * G1q.diff()
    + 2 * F1q.diff() * G2q
    - F1q * G2q.diff()
    - 3 * F3q * Eq.diff()
    + 3 * Cq.diff() * G3q
)
require(j1_residual.is_zero, "full rational J1 residual")
require(j2_residual.is_zero, "full rational J2 residual")


def target_hash(values):
    payload = ";".join(
        f"{i},{j},{k}:{value}"
        for value, ((i, j, k), _restriction, _delta) in zip(values, monomials)
        if value
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def polynomial_hash(polynomial):
    polynomial = sp.Poly(polynomial, x, domain=sp.QQ)
    payload = ";".join(
        f"{degree}:{polynomial.nth(degree)}" for degree in range(polynomial.degree() + 1)
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


expected_lifts = {
    "F2": (281, 369, EXPECTED["F2_target"], EXPECTED["F2_restriction"]),
    "G2": (278, 366, EXPECTED["G2_target"], EXPECTED["G2_restriction"]),
    "F3": (281, 369, EXPECTED["F3_target"], EXPECTED["F3_restriction"]),
    "G3": (104, 369, EXPECTED["G3_target"], EXPECTED["G3_restriction"]),
}
lift_results = {}
for label, vector, polynomial in zip(
    ("F2", "G2", "F3", "G3"), target_vectors, lift_restrictions
):
    result = (
        sum(bool(value) for value in vector),
        polynomial.degree(),
        target_hash(vector),
        polynomial_hash(polynomial),
    )
    require(result == expected_lifts[label], f"{label} frozen certificate")
    lift_results[label] = result

print("base_THM3678_replay=PASS;J0=1")
print(
    f"modulus={MODULUS};cutoff=366;monomials=953;rows=387+561;columns=3812;"
    f"rank={lower_operator_rank};augmented={lower_augmented_rank};feasible=False"
)
print(
    f"modulus={MODULUS};cutoff=369;monomials=973;rows=390+564;columns=3892;"
    f"rank={operator_rank};augmented={augmented_rank};feasible=True"
)
print(
    f"exact_skeleton=947x947;column_hash={pivot_column_hash};row_hash={pivot_row_hash}"
)
for label in ("F2", "G2", "F3", "G3"):
    support, degree, target_digest, restriction_digest = lift_results[label]
    print(
        f"{label}_support={support};target_hash={target_digest};degree={degree};"
        f"restriction_hash={restriction_digest}"
    )
print(f"RESULT=PASS;full_QQ_J0=1_J1=J2=0;gates={gates}")
