#!/usr/bin/env python3
"""Exact reconstruction companion for THM-3651.

The script has two logically separate halves.

1.  Over K=Q(s), s^2=94, it reconstructs one actual target-ring pair for
    the plus degree-seven double-zero fold with J0=1 and J1=J2=0.  Modular
    arithmetic selects K-column/K-row skeletons only.  Faithful rational
    2-by-2 block matrices solve the selected squares exactly, and the proof
    ends with full K[x] residual checks.  Coefficient hashes bind the chosen
    target representatives and their restrictions.

2.  Independently, it builds the complete arbitrary-target-two-form local
    pullback matrices through total source orders 4, 5, and 6.  Target
    coordinates are retained through degree 7 before differentiation.  The
    order-six inconsistency is certified by a pinned 25-support left row and
    an active hostile deletion.

Conjugating s sends every plus-fold statement to the minus fold.  Nothing in
this file asserts an actual lift beyond J2, an integrable/decomposable
realization of the arbitrary two-form order-five survivor, a Keller pair, or
a consequence for JC(2).
"""

import ast
import hashlib
from math import comb
from pathlib import Path

from flint import fmpq, fmpq_mat, fmpq_poly, fmpz, nmod_mat, nmod_poly


RADICAND = 94
PIVOT_PRIME = 1000033
PIVOT_ROOT = 20857
J0_CUTOFF = 160
LIFT_CUTOFF = 301
CHECKS = 0


def require(label, condition):
    """Optimization-safe exact gate."""
    global CHECKS
    if condition != True:
        raise RuntimeError(f"FAIL {label}: {condition}")
    CHECKS += 1


def q(value=0):
    return value if isinstance(value, fmpq) else fmpq(value)


def qpair(value):
    value = q(value)
    return int(value.numerator), int(value.denominator)


def qtext(value):
    numerator, denominator = qpair(value)
    return f"{numerator}/{denominator}"


def sha256_text(payload):
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def integer_vector_hash(values):
    return sha256_text(",".join(str(int(value)) for value in values))


# K elements are pairs (a,b) representing a+s*b, with s^2=94.
def ka(value=0, radical=0):
    return q(value), q(radical)


KZERO = ka()
KONE = ka(1)


def kadd(left, right):
    return left[0] + right[0], left[1] + right[1]


def kneg(value):
    return -value[0], -value[1]


def ksub(left, right):
    return left[0] - right[0], left[1] - right[1]


def kmul(left, right):
    return (
        left[0] * right[0] + RADICAND * left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def kscale(value, scalar):
    scalar = q(scalar)
    return value[0] * scalar, value[1] * scalar


def kconjugate(value):
    return value[0], -value[1]


def kiszero(value):
    return not value[0] and not value[1]


def kequal(left, right):
    return left[0] == right[0] and left[1] == right[1]


def ktext(value):
    return f"{qtext(value[0])}+s*{qtext(value[1])}"


def knorm(value):
    return value[0] * value[0] - RADICAND * value[1] * value[1]


def kpower(value, exponent):
    answer = KONE
    for _ in range(exponent):
        answer = kmul(answer, value)
    return answer


# K[x] elements are pairs of fmpq_poly objects.
def kp(real=0, radical=0):
    real_poly = real if isinstance(real, fmpq_poly) else fmpq_poly(real)
    radical_poly = radical if isinstance(radical, fmpq_poly) else fmpq_poly(radical)
    return real_poly, radical_poly


KPZERO = kp()
KPONE = kp([1])


def kpadd(left, right):
    return left[0] + right[0], left[1] + right[1]


def kpneg(value):
    return -value[0], -value[1]


def kpsub(left, right):
    return left[0] - right[0], left[1] - right[1]


def kpmul(left, right):
    return (
        left[0] * right[0] + RADICAND * left[1] * right[1],
        left[0] * right[1] + left[1] * right[0],
    )


def kpscale(value, scalar):
    scalar = q(scalar)
    return value[0] * scalar, value[1] * scalar


def kpmul_element(value, scalar):
    return (
        value[0] * scalar[0] + RADICAND * value[1] * scalar[1],
        value[0] * scalar[1] + value[1] * scalar[0],
    )


def kpderivative(value):
    return value[0].derivative(), value[1].derivative()


def kpdegree(value):
    return max(value[0].degree(), value[1].degree())


def kpcoeff(value, degree):
    return value[0][degree], value[1][degree]


def kpiszero(value):
    return value[0].is_zero() and value[1].is_zero()


def kpequal(left, right):
    return left[0] == right[0] and left[1] == right[1]


def kppow(value, exponent):
    answer = KPONE
    for _ in range(exponent):
        answer = kpmul(answer, value)
    return answer


def kptranslate(value, point):
    substitution = fmpq_poly([point, 1])
    return value[0](substitution), value[1](substitution)


def kpolynomial_hash(value):
    degree = kpdegree(value)
    payload = ";".join(ktext(kpcoeff(value, index)) for index in range(degree, -1, -1))
    return sha256_text(payload)


def ktarget_hash(values, packet):
    entries = []
    for value, (metadata, _restriction, _delta) in zip(values, packet):
        if kiszero(value):
            continue
        entries.append(f"{metadata[0]},{metadata[1]},{metadata[2]}:{ktext(value)}")
    return sha256_text(";".join(entries))


def modular_q(value, prime=PIVOT_PRIME):
    numerator, denominator = qpair(value)
    return numerator * pow(denominator, -1, prime) % prime


def modular_k(value, prime=PIVOT_PRIME, root=PIVOT_ROOT):
    return (modular_q(value[0], prime) + root * modular_q(value[1], prime)) % prime


def kp_to_nmod(value, prime=PIVOT_PRIME, root=PIVOT_ROOT):
    degree = kpdegree(value)
    if degree < 0:
        return nmod_poly([], prime)
    return nmod_poly([modular_k(kpcoeff(value, index), prime, root) for index in range(degree + 1)], prime)


def selected_k_solve(matrix_entries, rhs):
    """Solve one square K-system via its faithful rational block matrix."""
    size = len(rhs)
    require("selected K square row count", len(matrix_entries) == size)
    require("selected K square column count", all(len(row) == size for row in matrix_entries))
    block_data = []
    for row in matrix_entries:
        block_data.extend(value[0] for value in row)
        block_data.extend(RADICAND * value[1] for value in row)
    for row in matrix_entries:
        block_data.extend(value[1] for value in row)
        block_data.extend(value[0] for value in row)
    block_rhs = [value[0] for value in rhs] + [value[1] for value in rhs]
    rational_matrix = fmpq_mat(2 * size, 2 * size, block_data)
    rational_rhs = fmpq_mat(2 * size, 1, block_rhs)
    rational_solution = rational_matrix.solve(rational_rhs)
    require("selected rational block residual", rational_matrix * rational_solution == rational_rhs)
    return [
        (rational_solution[index, 0], rational_solution[size + index, 0])
        for index in range(size)
    ]


print("THM-3651 reconstruction -- degree-seven double-zero sixth-order closure")
print("status=VERIFIED_EXACT_RECONSTRUCTION")
print("field=Q(s);s^2=94;fold=plus;minus=coefficientwise_conjugate")
require("94 is nonsquare over Q", not fmpz(RADICAND).is_square())
require("pivot modulus is prime", fmpz(PIVOT_PRIME).is_prime() == 1)
require(
    "pivot root squares to 94",
    PIVOT_ROOT * PIVOT_ROOT % PIVOT_PRIME == RADICAND,
)


# ---------------------------------------------------------------------------
# Compiler and the two conjugate folds.
# ---------------------------------------------------------------------------
xpoly = kp([0, 1])
x2 = kpmul(xpoly, xpoly)
Q1 = kp([
    q(-3) / 4,
    1,
    q(-27) / 4,
    -2,
    q(9) / 2,
    1,
])
P = kp([0, 0, 1, 0, -2, 0, 1])
beta = ka(q(-211) / 130, q(99) / 260)
alpha = ka(q(-1683) / 260, q(-11) / 65)
require("alpha-beta fourth-zero line", kequal(alpha, kscale(kadd(ka(259), kscale(beta, 16)), q(-1) / 36)))
require("beta quadratic", kiszero(kadd(kadd(kscale(kmul(beta, beta), 520), kscale(beta, 1688)), ka(-5717))))
Q = kpadd(Q1, kpmul(P, kp([alpha[0], beta[0]], [alpha[1], beta[1]])))
require("plus fold degree seven", kpdegree(Q) == 7)
Q6 = kpadd(Q1, kpscale(P, q(-259) / 36))
r_parameter = kscale(beta, q(-1) / 9)
Q_beta_line = kpadd(
    Q6,
    kpmul_element(kpmul(P, kpsub(xpoly, kp([q(4) / 9]))), beta),
)
require("Qbeta line identity", kpequal(Q, Q_beta_line))
p_of_r = kadd(
    kadd(
        kscale(kmul(r_parameter, r_parameter), q(520) / 9),
        kscale(r_parameter, q(-1688) / 81),
    ),
    ka(q(-5717) / 729),
)
require("THM3677 parabola parameter vanishes", kiszero(p_of_r))
R1 = kpmul(P, kpsub(KPONE, x2))
R2 = kpmul(P, kp([4, -9]))
Q_parabola_specialization = kpadd(
    Q6,
    kpadd(kpmul_element(R1, p_of_r), kpmul_element(R2, r_parameter)),
)
require(
    "Qbeta equals THM3677 p=0 r=-beta/9 specialization",
    kpequal(Q, Q_parabola_specialization),
)

D = kpadd(KPONE, kpmul(x2, Q))
Dminus1 = kpsub(D, KPONE)
Dplus2 = kpadd(D, kp([2]))
Dplus3 = kpadd(D, kp([3]))
b_restricted = kpmul(Dminus1, kpmul(Dplus2, Dplus2))
c_restricted = kpmul(xpoly, kpmul(D, Dplus2))
e_restricted = kpmul(Q, Dplus3)
restricted = (b_restricted, c_restricted, e_restricted)
restriction_degrees = tuple(kpdegree(value) for value in restricted)
require("restriction degrees", restriction_degrees == (27, 19, 16))
require(
    "restricted Russell surface relation",
    kpiszero(
        kpsub(
            kpmul(c_restricted, kpmul(c_restricted, e_restricted)),
            kpmul(b_restricted, kpadd(b_restricted, kp([4]))),
        )
    ),
)

# Partial-q derivatives of the compiler, restricted at q=Q.
delta_b = kpscale(kpmul(x2, kpmul(D, Dplus2)), 3)
delta_c = kpscale(kpmul(kpmul(xpoly, x2), kpadd(D, KPONE)), 2)
delta_e = kpadd(kpscale(D, 2), kp([2]))
deltas = (delta_b, delta_c, delta_e)
require(
    "differentiated Russell surface relation",
    kpiszero(
        kpsub(
            kpadd(
                kpscale(
                    kpmul(c_restricted, kpmul(delta_c, e_restricted)),
                    2,
                ),
                kpmul(kpmul(c_restricted, c_restricted), delta_e),
            ),
            kpmul(kpadd(kpscale(b_restricted, 2), kp([4])), delta_b),
        )
    ),
)


def eval_kpoly(value, point):
    return value[0](point), value[1](point)


points = (-1, 0, 1)
fold_values = tuple(eval_kpoly(Q, point) for point in points)
fold_slopes = tuple(eval_kpoly(kpderivative(Q), point) for point in points)
require("retained fold values", fold_values == (ka(-3), ka(q(-3) / 4), ka(-3)))
require("retained fold slopes", fold_slopes == (ka(q(-9) / 2), ka(1), ka(q(9) / 2)))
curvatures = tuple(eval_kpoly(kpderivative(kpderivative(Q)), point) for point in points)
second_debt = kadd(kadd(kscale(curvatures[0], 5), kscale(curvatures[2], 13)), ka(243))
require("zero second debt", kiszero(second_debt))
print("PASS folds=Q1+P*(alpha_plus+beta_plus*x),conjugate;degree=7;D2=D4=0")


def exact_monomials(cutoff):
    powers = []
    for generator, degree in zip(restricted, restriction_degrees):
        row = [KPONE]
        for _ in range(cutoff // degree):
            row.append(kpmul(row[-1], generator))
        powers.append(row)

    answer = []
    for b_exponent in range(cutoff // restriction_degrees[0] + 1):
        for c_exponent in range(cutoff // restriction_degrees[1] + 1):
            for e_exponent in range(cutoff // restriction_degrees[2] + 1):
                metadata = (b_exponent, c_exponent, e_exponent)
                weighted_degree = sum(
                    exponent * degree
                    for exponent, degree in zip(metadata, restriction_degrees)
                )
                if weighted_degree > cutoff:
                    continue
                monomial = kpmul(
                    powers[0][b_exponent],
                    kpmul(powers[1][c_exponent], powers[2][e_exponent]),
                )
                delta_monomial = KPZERO
                if b_exponent:
                    delta_monomial = kpadd(
                        delta_monomial,
                        kpscale(
                            kpmul(
                                kpmul(
                                    powers[0][b_exponent - 1],
                                    powers[1][c_exponent],
                                ),
                                kpmul(powers[2][e_exponent], deltas[0]),
                            ),
                            b_exponent,
                        ),
                    )
                if c_exponent:
                    delta_monomial = kpadd(
                        delta_monomial,
                        kpscale(
                            kpmul(
                                kpmul(
                                    powers[0][b_exponent],
                                    powers[1][c_exponent - 1],
                                ),
                                kpmul(powers[2][e_exponent], deltas[1]),
                            ),
                            c_exponent,
                        ),
                    )
                if e_exponent:
                    delta_monomial = kpadd(
                        delta_monomial,
                        kpscale(
                            kpmul(
                                kpmul(
                                    powers[0][b_exponent],
                                    powers[1][c_exponent],
                                ),
                                kpmul(powers[2][e_exponent - 1], deltas[2]),
                            ),
                            e_exponent,
                        ),
                    )
                answer.append((metadata, monomial, delta_monomial))
    return answer


def modular_pivot_skeleton(operator, augmented):
    reduced_augmented, augmented_rank = augmented.rref()
    operator_rank = operator.rank()
    require("modular operator/augmented rank equality", operator_rank == augmented_rank)
    pivot_columns = []
    for row in range(augmented_rank):
        pivot = next(column for column in range(operator.ncols() + 1) if reduced_augmented[row, column])
        require(f"operator pivot row {row}", pivot < operator.ncols())
        pivot_columns.append(pivot)
    selected = nmod_mat(
        operator.nrows(),
        len(pivot_columns),
        [operator[row, column] for row in range(operator.nrows()) for column in pivot_columns],
        PIVOT_PRIME,
    )
    transpose_reduced, transpose_rank = selected.transpose().rref()
    require("selected modular column rank", transpose_rank == augmented_rank)
    pivot_rows = [
        next(column for column in range(operator.nrows()) if transpose_reduced[row, column])
        for row in range(transpose_rank)
    ]
    return operator_rank, augmented_rank, pivot_columns, pivot_rows


print("SECTION actual_target_ring_lift_through_J2")
print("monomial_order=nested_increasing_b,c,e weighted_by=27,19,16")
print("J0_blocks=G1,F1 coupled_blocks=F2,G2,F3,G3 rows=J1_then_J2_ascending_x")

base_packet = exact_monomials(J0_CUTOFF)
base_polynomials = [item[1] for item in base_packet]
require("base monomial count", len(base_packet) == 141)
j0_columns = [kpmul(kpderivative(restricted[1]), polynomial) for polynomial in base_polynomials]
j0_columns += [kpneg(kpmul(polynomial, kpderivative(restricted[2]))) for polynomial in base_polynomials]
j0_row_count = max(kpdegree(column) for column in j0_columns) + 1
require("J0 row envelope", j0_row_count == 179)
j0_operator = nmod_mat(
    j0_row_count,
    len(j0_columns),
    [modular_k(kpcoeff(j0_columns[column], row)) for row in range(j0_row_count) for column in range(len(j0_columns))],
    PIVOT_PRIME,
)
j0_augmented = nmod_mat(
    j0_row_count,
    len(j0_columns) + 1,
    [
        value
        for row in range(j0_row_count)
        for value in (
            *[modular_k(kpcoeff(column, row)) for column in j0_columns],
            1 if row == 0 else 0,
        )
    ],
    PIVOT_PRIME,
)
j0_rank, j0_aug_rank, j0_pivot_columns, j0_pivot_rows = modular_pivot_skeleton(j0_operator, j0_augmented)
require("J0 modular ranks", (j0_rank, j0_aug_rank) == (178, 178))
require("J0 pivot column hash", integer_vector_hash(j0_pivot_columns) == "7784507c6b1b77d78bfbd33fa2385eaa0c5a0c4e1a457c50faa6c27737a07768")
require("J0 pivot row hash", integer_vector_hash(j0_pivot_rows) == "2113e21fc204fdc9fc3be762bc1d022a2388bb14a8390ce4b8603ac661de7454")
print(
    f"J0_cutoff={J0_CUTOFF} monomials={len(base_packet)} shape={j0_row_count}x{len(j0_columns)} "
    f"rank={j0_rank} selected_K_square=178 "
    f"pivot_column_hash={integer_vector_hash(j0_pivot_columns)} "
    f"pivot_row_hash={integer_vector_hash(j0_pivot_rows)}"
)

j0_square = [
    [kpcoeff(j0_columns[column], row) for column in j0_pivot_columns]
    for row in j0_pivot_rows
]
j0_rhs = [KONE if row == 0 else KZERO for row in j0_pivot_rows]
j0_selected_solution = selected_k_solve(j0_square, j0_rhs)
j0_solution = [KZERO for _ in j0_columns]
for column, value in zip(j0_pivot_columns, j0_selected_solution):
    j0_solution[column] = value
G1_vector = j0_solution[: len(base_packet)]
F1_vector = j0_solution[len(base_packet) :]


def linear_combination(values, polynomials):
    answer = KPZERO
    for value, polynomial in zip(values, polynomials):
        if not kiszero(value):
            answer = kpadd(answer, kpmul_element(polynomial, value))
    return answer


F1 = linear_combination(F1_vector, base_polynomials)
G1 = linear_combination(G1_vector, base_polynomials)
delta_F1 = linear_combination(F1_vector, [item[2] for item in base_packet])
delta_G1 = linear_combination(G1_vector, [item[2] for item in base_packet])
j0_residual = kpsub(kpmul(kpderivative(restricted[1]), G1), kpmul(F1, kpderivative(restricted[2])))
require("full K[x] J0 identity", kpequal(j0_residual, KPONE))

expected_base_data = {
    "F1": (89, "0ee493528fa3b1bf9ca6a824af0a6b7cdfd1588c2a99f358b171a7cea0f219d8", 160, "89178122fa5660177bc2ec8f2c097232fee112efa5ccba496f5ec1a2574e13bc", 153, "ce48af3830f332279a6f9c6e5ad26883baf59d9b3c33ce488901f6d1e247f019"),
    "G1": (86, "6b02b7b4af6d9a50ec526e0367802a811ebfbc03ab263aa660cd67e772346a38", 157, "1e21a0e832129ef01b121061019923d04df1d7faae2f7440601dc516e5049d3f", 150, "a53ae91d4966f864882344580697213481687c81da1ee13c96f840a1c4e75d26"),
}
for label, vector, polynomial, delta in (
    ("F1", F1_vector, F1, delta_F1),
    ("G1", G1_vector, G1, delta_G1),
):
    support = sum(not kiszero(value) for value in vector)
    target_digest = ktarget_hash(vector, base_packet)
    restriction_digest = kpolynomial_hash(polynomial)
    delta_digest = kpolynomial_hash(delta)
    expected = expected_base_data[label]
    require(f"{label} support", support == expected[0])
    require(f"{label} target hash", target_digest == expected[1])
    require(f"{label} restriction degree", kpdegree(polynomial) == expected[2])
    require(f"{label} restriction hash", restriction_digest == expected[3])
    require(f"{label} delta degree", kpdegree(delta) == expected[4])
    require(f"{label} delta hash", delta_digest == expected[5])
    print(
        f"{label}_support={support} target_hash={target_digest} "
        f"restriction_degree={kpdegree(polynomial)} restriction_hash={restriction_digest} "
        f"delta_degree={kpdegree(delta)} delta_hash={delta_digest}"
    )

lift_packet = exact_monomials(LIFT_CUTOFF)
lift_polynomials = [item[1] for item in lift_packet]
require("lift monomial count", len(lift_packet) == 744)


def modular_lift_system():
    monomials = [kp_to_nmod(polynomial) for polynomial in lift_polynomials]
    c = kp_to_nmod(restricted[1])
    e = kp_to_nmod(restricted[2])
    delta_c_mod = kp_to_nmod(deltas[1])
    delta_e_mod = kp_to_nmod(deltas[2])
    f1 = kp_to_nmod(F1)
    g1 = kp_to_nmod(G1)
    delta_f1 = kp_to_nmod(delta_F1)
    delta_g1 = kp_to_nmod(delta_G1)
    zero = nmod_poly([], PIVOT_PRIME)
    columns = []
    for polynomial in monomials:
        columns.append((-2 * polynomial * e.derivative(), polynomial.derivative() * g1 - 2 * polynomial * g1.derivative()))
    for polynomial in monomials:
        columns.append((2 * c.derivative() * polynomial, 2 * f1.derivative() * polynomial - f1 * polynomial.derivative()))
    for polynomial in monomials:
        columns.append((zero, -3 * polynomial * e.derivative()))
    for polynomial in monomials:
        columns.append((zero, 3 * c.derivative() * polynomial))
    known_j1 = 2 * c.derivative() * delta_e_mod + f1.derivative() * g1 - f1 * g1.derivative() - 2 * delta_c_mod * e.derivative()
    known_j2 = 3 * c.derivative() * delta_g1 + 2 * f1.derivative() * delta_e_mod + delta_c_mod.derivative() * g1 - f1 * delta_e_mod.derivative() - 2 * delta_c_mod * g1.derivative() - 3 * delta_f1 * e.derivative()
    row_count_j1 = max([known_j1.degree()] + [pair[0].degree() for pair in columns]) + 1
    row_count_j2 = max([known_j2.degree()] + [pair[1].degree() for pair in columns]) + 1
    operator = nmod_mat(
        row_count_j1 + row_count_j2,
        len(columns),
        [
            value
            for block in (0, 1)
            for row in range((row_count_j1, row_count_j2)[block])
            for value in [int(pair[block][row]) for pair in columns]
        ],
        PIVOT_PRIME,
    )
    augmented = nmod_mat(
        row_count_j1 + row_count_j2,
        len(columns) + 1,
        [
            value
            for block, known in ((0, known_j1), (1, known_j2))
            for row in range((row_count_j1, row_count_j2)[block])
            for value in (*[int(pair[block][row]) for pair in columns], int(-known[row]))
        ],
        PIVOT_PRIME,
    )
    return operator, augmented, row_count_j1, row_count_j2


lift_operator, lift_augmented, row_count_j1, row_count_j2 = modular_lift_system()
require("lift row envelopes", (row_count_j1, row_count_j2) == (320, 461))
lift_rank, lift_aug_rank, lift_pivot_columns, lift_pivot_rows = modular_pivot_skeleton(lift_operator, lift_augmented)
require("lift modular ranks", (lift_rank, lift_aug_rank) == (776, 776))
require("lift pivot column hash", integer_vector_hash(lift_pivot_columns) == "a9697a7a900ca6ddb98a14981d01847456e04cfb48932e8c300d07065c7084d6")
require("lift pivot row hash", integer_vector_hash(lift_pivot_rows) == "449e1c433c85fe550e6dda35664a1584b70fe5c4da4846003173ade228fbbfef")
print(
    f"pivot_prime={PIVOT_PRIME} sqrt94={PIVOT_ROOT} operator_shape={lift_operator.nrows()}x{lift_operator.ncols()} "
    f"row_envelopes={row_count_j1}+{row_count_j2} operator_rank={lift_rank} augmented_rank={lift_aug_rank} "
    f"pivot_column_hash={integer_vector_hash(lift_pivot_columns)} "
    f"pivot_row_hash={integer_vector_hash(lift_pivot_rows)}"
)

known_j1_exact = kpadd(
    kpsub(
        kpadd(kpscale(kpmul(kpderivative(restricted[1]), deltas[2]), 2), kpmul(kpderivative(F1), G1)),
        kpmul(F1, kpderivative(G1)),
    ),
    kpneg(kpscale(kpmul(deltas[1], kpderivative(restricted[2])), 2)),
)
known_j2_exact = kpadd(
    kpadd(
        kpadd(
            kpscale(kpmul(kpderivative(restricted[1]), delta_G1), 3),
            kpscale(kpmul(kpderivative(F1), deltas[2]), 2),
        ),
        kpmul(kpderivative(deltas[1]), G1),
    ),
    kpneg(
        kpadd(
            kpadd(kpmul(F1, kpderivative(deltas[2])), kpscale(kpmul(deltas[1], kpderivative(G1)), 2)),
            kpscale(kpmul(delta_F1, kpderivative(restricted[2])), 3),
        )
    ),
)


def exact_lift_column(column):
    block, index = divmod(column, len(lift_polynomials))
    polynomial = lift_polynomials[index]
    if block == 0:
        return (
            kpneg(kpscale(kpmul(polynomial, kpderivative(restricted[2])), 2)),
            kpsub(kpmul(kpderivative(polynomial), G1), kpscale(kpmul(polynomial, kpderivative(G1)), 2)),
        )
    if block == 1:
        return (
            kpscale(kpmul(kpderivative(restricted[1]), polynomial), 2),
            kpsub(kpscale(kpmul(kpderivative(F1), polynomial), 2), kpmul(F1, kpderivative(polynomial))),
        )
    if block == 2:
        return KPZERO, kpneg(kpscale(kpmul(polynomial, kpderivative(restricted[2])), 3))
    if block == 3:
        return KPZERO, kpscale(kpmul(kpderivative(restricted[1]), polynomial), 3)
    raise RuntimeError(f"bad lift block {block}")


selected_exact_columns = [exact_lift_column(column) for column in lift_pivot_columns]
lift_square = []
lift_rhs = []
for row in lift_pivot_rows:
    if row < row_count_j1:
        lift_square.append([kpcoeff(pair[0], row) for pair in selected_exact_columns])
        lift_rhs.append(kneg(kpcoeff(known_j1_exact, row)))
    else:
        degree = row - row_count_j1
        lift_square.append([kpcoeff(pair[1], degree) for pair in selected_exact_columns])
        lift_rhs.append(kneg(kpcoeff(known_j2_exact, degree)))

lift_selected_solution = selected_k_solve(lift_square, lift_rhs)
lift_solution = [KZERO for _ in range(4 * len(lift_polynomials))]
for column, value in zip(lift_pivot_columns, lift_selected_solution):
    lift_solution[column] = value
target_vectors = [
    lift_solution[block * len(lift_polynomials) : (block + 1) * len(lift_polynomials)]
    for block in range(4)
]
restrictions = [linear_combination(vector, lift_polynomials) for vector in target_vectors]
F2, G2, F3, G3 = restrictions
j1_residual = kpadd(
    known_j1_exact,
    kpadd(
        kpneg(kpscale(kpmul(F2, kpderivative(restricted[2])), 2)),
        kpscale(kpmul(kpderivative(restricted[1]), G2), 2),
    ),
)
j2_residual = kpadd(
    known_j2_exact,
    kpadd(
        kpadd(
            kpsub(kpmul(kpderivative(F2), G1), kpscale(kpmul(F2, kpderivative(G1)), 2)),
            kpsub(kpscale(kpmul(kpderivative(F1), G2), 2), kpmul(F1, kpderivative(G2))),
        ),
        kpadd(
            kpneg(kpscale(kpmul(F3, kpderivative(restricted[2])), 3)),
            kpscale(kpmul(kpderivative(restricted[1]), G3), 3),
        ),
    ),
)
require("full K[x] J1 identity", kpiszero(j1_residual))
require("full K[x] J2 identity", kpiszero(j2_residual))
print("exact_selected_K_square=776 rational_block=1552 full_K[x]_J0=1_J1=J2=0")

expected_lift_data = {
    "F2": (230, "7e38ac68a6cbd4702e1c50ba566849fa9d0f9198b058523111f1a6c0070d863c", 301, "7e7c7072644f5870d9798a1b157759f92bfb4406a4d8deaf560c76e3f53b1483"),
    "G2": (228, "0acaf6b4655dc7cea05fae327c4cb89de446c3c618b6768c1f290be22ce4f2b9", 298, "884f6b51f22cd3f22a1834f9663a9e7146823b24fd04cbeb8b28e08fe684d830"),
    "F3": (230, "c4007f9663bbd876cb663afbb796fd7e7904f346e8c6746e60cec8175fe5987e", 301, "d066e6a37d2d878dc648f937c2467863d6561985d6d1038323263045d79ad62c"),
    "G3": (86, "d7b163e2ec6949316a0dcccdbf105302dc69e36bbcc05ec67d61a2979b044ccb", 300, "757b39f18070c1bdb46a9eb199edc42dd0efaf7b4f1d0aba1866f2c1b1dd3a92"),
}
for label, vector, polynomial in zip(("F2", "G2", "F3", "G3"), target_vectors, restrictions):
    support = sum(not kiszero(value) for value in vector)
    target_digest = ktarget_hash(vector, lift_packet)
    restriction_digest = kpolynomial_hash(polynomial)
    expected = expected_lift_data[label]
    require(f"{label} support", support == expected[0])
    require(f"{label} target hash", target_digest == expected[1])
    require(f"{label} restriction degree", kpdegree(polynomial) == expected[2])
    require(f"{label} restriction hash", restriction_digest == expected[3])
    print(
        f"{label}_support={support} target_hash={target_digest} "
        f"restriction_degree={kpdegree(polynomial)} restriction_hash={restriction_digest}"
    )
print("PASS actual_lift=plus_J0=1_J1=J2=0;minus=conjugate;cutoffs_chosen_not_minimal")


# ---------------------------------------------------------------------------
# Complete local arbitrary-two-form matrices through total order six.
# ---------------------------------------------------------------------------
print("SECTION arbitrary_two_form_total_order_matrices")


def jet_add(left, right):
    answer = dict(left)
    for key, value in right.items():
        answer[key] = kadd(answer.get(key, KZERO), value)
        if kiszero(answer[key]):
            del answer[key]
    return answer


def jet_scale(value, scalar):
    answer = {}
    for key, coefficient in value.items():
        product = kmul(coefficient, scalar) if isinstance(scalar, tuple) else kscale(coefficient, scalar)
        if not kiszero(product):
            answer[key] = product
    return answer


def jet_multiply(left, right, cutoff=7):
    answer = {}
    for (left_h, left_t), left_value in left.items():
        for (right_h, right_t), right_value in right.items():
            key = (left_h + right_h, left_t + right_t)
            if sum(key) > cutoff:
                continue
            answer[key] = kadd(answer.get(key, KZERO), kmul(left_value, right_value))
    return {key: value for key, value in answer.items() if not kiszero(value)}


def jet_power(value, exponent, cutoff=7):
    answer = {(0, 0): KONE}
    for _ in range(exponent):
        answer = jet_multiply(answer, value, cutoff)
    return answer


def jet_derivative(value, coordinate):
    answer = {}
    for (h_degree, t_degree), coefficient in value.items():
        degree = (h_degree, t_degree)[coordinate]
        if not degree:
            continue
        key = (h_degree - (coordinate == 0), t_degree - (coordinate == 1))
        answer[key] = kscale(coefficient, degree)
    return answer


def polynomial_as_h_jet(value, point, cutoff=7):
    translated = kptranslate(value, point)
    return {
        (degree, 0): kpcoeff(translated, degree)
        for degree in range(min(cutoff, kpdegree(translated)) + 1)
        if not kiszero(kpcoeff(translated, degree))
    }


local_packets = {}
for point in points:
    local_x = {(0, 0): ka(point), (1, 0): KONE}
    local_q = jet_add(polynomial_as_h_jet(Q, point), {(0, 2): KONE})
    local_D = jet_add({(0, 0): KONE}, jet_multiply(jet_power(local_x, 2), local_q))
    local_y = jet_scale(jet_multiply(local_x, jet_multiply(local_D, jet_add(local_D, {(0, 0): ka(2)}))), q(1) / 3)
    local_z = jet_add(jet_multiply(local_q, jet_add(local_D, {(0, 0): ka(3)})), {(0, 0): ka(3)})
    require(f"local y vanishes {point}", kiszero(local_y.get((0, 0), KZERO)))
    require(f"local z vanishes {point}", kiszero(local_z.get((0, 0), KZERO)))
    y_h = jet_derivative(local_y, 0)
    y_t = jet_derivative(local_y, 1)
    z_h = jet_derivative(local_z, 0)
    z_t = jet_derivative(local_z, 1)
    area = jet_add(jet_multiply(y_h, z_t, 6), jet_scale(jet_multiply(y_t, z_h, 6), -1))
    local_packets[point] = (
        [jet_power(local_y, exponent, 7) for exponent in range(7)],
        [jet_power(local_z, exponent, 7) for exponent in range(7)],
        (area, y_h, z_h),
    )

expected_regression = ka(q(1253248661) / 4225, q(477739493) / 25350)
require("degree-seven-before-differentiation regression", kequal(local_packets[-1][2][1].get((6, 0), KZERO), expected_regression))


def local_column(kind, y_degree, z_degree, w_degree, point, cutoff):
    y_powers, z_powers, bases = local_packets[point]
    value = jet_multiply(y_powers[y_degree], z_powers[z_degree], cutoff)
    value = jet_multiply(value, {(0, w_degree): KONE}, cutoff)
    return jet_multiply(value, bases[kind], cutoff)


def row_labels(order):
    return tuple(
        (point, h_degree, total_degree - h_degree)
        for point in points
        for total_degree in range(order + 1)
        for h_degree in range(total_degree + 1)
    )


def column_labels(order):
    return tuple(
        (kind, y_degree, z_degree, w_degree)
        for kind in range(3)
        for y_degree in range(order + 1)
        for z_degree in range(order + 1 - y_degree)
        for w_degree in range(order + 1 - y_degree - z_degree)
    )


def local_matrix(order):
    rows = row_labels(order)
    columns = column_labels(order)
    data = []
    for kind, y_degree, z_degree, w_degree in columns:
        branch_values = {
            point: local_column(kind, y_degree, z_degree, w_degree, point, order)
            for point in points
        }
        data.append(tuple(branch_values[point].get((h_degree, t_degree), KZERO) for point, h_degree, t_degree in rows))
    # Use the historical denominator-clearing normalization 12.  Division by
    # 12 gives the identical solvability/obstruction statement for unit
    # density, but the pinned certificate pairing below is against this RHS.
    target = tuple(ka(12) if h_degree == 0 and t_degree == 0 else KZERO for _point, h_degree, t_degree in rows)
    return rows, columns, data, target


def rational_block_rank(columns, row_count, extra=()):
    all_columns = [*columns, *extra]
    column_count = len(all_columns)
    real_rows = []
    radical_rows = []
    for row in range(row_count):
        real_rows.extend(value[row][0] for value in all_columns)
        real_rows.extend(RADICAND * value[row][1] for value in all_columns)
        radical_rows.extend(value[row][1] for value in all_columns)
        radical_rows.extend(value[row][0] for value in all_columns)
    matrix = fmpq_mat(2 * row_count, 2 * column_count, real_rows + radical_rows)
    rank = matrix.rank()
    require("rational block rank even", rank % 2 == 0)
    return rank // 2


matrix_packets = {}
for order, expected in ((4, (45, 105, 40, 40)), (5, (63, 168, 57, 57)), (6, (84, 252, 77, 78))):
    rows, columns, data, target = local_matrix(order)
    operator_rank = rational_block_rank(data, len(rows))
    augmented_rank = rational_block_rank(data, len(rows), (target,))
    require(f"N{order} shape/ranks", (len(rows), len(columns), operator_rank, augmented_rank) == expected)
    matrix_packets[order] = (rows, columns, data, target)
    print(f"N{order}=rows{len(rows)}_columns{len(columns)}_rank{operator_rank}_augmented{augmented_rank}")


certificate_table = {
    (-1, 0, 0): ka(q(-416265780017) / 23562825, q(129267792) / 54925),
    (-1, 1, 0): ka(q(2464348103) / 966680, q(75921) / 10985),
    (-1, 0, 2): ka(q(-19437583987) / 5800080, q(163017) / 21970),
    (-1, 2, 0): ka(q(6975) / 143, q(20) / 3),
    (-1, 1, 2): ka(q(-47295) / 572, q(-95) / 13),
    (-1, 3, 0): ka(q(25) / 11),
    (-1, 0, 4): ka(q(39555) / 286, q(90) / 13),
    (-1, 2, 2): ka(q(-75) / 22),
    (-1, 1, 4): ka(q(225) / 44),
    (-1, 0, 6): ka(q(-675) / 88),
    (0, 0, 0): ka(q(-160054491789) / 7854275),
    (0, 1, 0): ka(q(1342208) / 1859, q(10944) / 169),
    (0, 0, 2): ka(q(46861679809) / 4833400, q(-163017) / 21970),
    (0, 1, 2): ka(q(1944) / 143, q(-108) / 13),
    (0, 0, 4): ka(q(-69093) / 143, q(-90) / 13),
    (0, 0, 6): ka(q(1215) / 44),
    (1, 1, 0): ka(q(-1521051607) / 371800),
    (1, 0, 2): ka(q(-14152473763) / 2230800),
    (1, 2, 0): ka(q(1467) / 11, q(4) / 3),
    (1, 1, 2): ka(q(9459) / 44, q(1)),
    (1, 3, 0): ka(q(-65) / 11),
    (1, 0, 4): ka(q(7587) / 22),
    (1, 2, 2): ka(q(-195) / 22),
    (1, 1, 4): ka(q(-585) / 44),
    (1, 0, 6): ka(q(-1755) / 88),
}

rows6, columns6, data6, target6 = matrix_packets[6]
certificate = tuple(certificate_table.get(row, KZERO) for row in rows6)
require("sixth certificate support", sum(not kiszero(value) for value in certificate) == 25)
for column, label in zip(data6, columns6):
    pairing = KZERO
    for weight, entry in zip(certificate, column):
        pairing = kadd(pairing, kmul(weight, entry))
    require(f"N6 certificate annihilates {label}", kiszero(pairing))

target_pairing = KZERO
for weight, entry in zip(certificate, target6):
    target_pairing = kadd(target_pairing, kmul(weight, entry))
expected_pairing = ka(q(-275824386272) / 604175, q(1551213504) / 54925)
require("sixth certificate target pairing", kequal(target_pairing, expected_pairing))
require("sixth certificate norm", knorm(target_pairing) == q(177370060592178176) / 1329185)
unit_target_pairing = kscale(target_pairing, q(1) / 12)
require(
    "unit target pairing",
    kequal(
        unit_target_pairing,
        ka(q(-68956096568) / 1812525, q(129267792) / 54925),
    ),
)


def certificate_hash(values):
    payload = ";".join(
        f"{point},{h_degree},{t_degree}:{ktext(value)}"
        for (point, h_degree, t_degree), value in zip(rows6, values)
        if not kiszero(value)
    )
    return sha256_text(payload)


plus_certificate_hash = certificate_hash(certificate)
minus_certificate_hash = certificate_hash(tuple(kconjugate(value) for value in certificate))
require("plus certificate hash", plus_certificate_hash == "acd75c71e93c95a163fdc91b5d241f3041a7ea98ef6118fee237e974a224a753")
require("minus certificate hash", minus_certificate_hash == "599a55ad0e15634dfde922c17c8e98a397a7f2b141a51b8b5ee36707730cb4dc")

# THM-3683 specialization cross-control.  The quartic formula is inherited
# from that independently hostile-audited package, but its specialization is
# recomputed here.
quartic_value = ka()
for degree, coefficient in enumerate(
    (-1276420, 7849770, -28419741, -77822208, 72783360)
):
    quartic_value = kadd(
        quartic_value,
        kscale(kpower(r_parameter, degree), coefficient),
    )
quartic_reduction = kscale(
    kadd(kscale(beta, 21544632), ka(-97639283)),
    q(2187) / 33800,
)
require("quartic beta reduction", kequal(quartic_value, quartic_reduction))
sixth_debt = kscale(quartic_value, q(-256) / 1594323)
require(
    "THM3683 specialization cross-control value",
    kequal(
        sixth_debt,
        ka(q(275824386272) / 200201625, q(-210658624) / 2471625),
    ),
)
require("certificate/debt scale", kequal(target_pairing, kscale(sixth_debt, q(-3645) / 11)))

mutated = list(certificate)
deleted_index = rows6.index((1, 0, 6))
mutated[deleted_index] = KZERO
nonzero_mutations = []
for label, column in zip(columns6, data6):
    pairing = KZERO
    for weight, entry in zip(mutated, column):
        pairing = kadd(pairing, kmul(weight, entry))
    if not kiszero(pairing):
        nonzero_mutations.append((label, pairing))
require("hostile deletion nonzero count", len(nonzero_mutations) == 47)
hostile_label = (0, 0, 0, 1)  # (dy^dz)*w
hostile_value = dict(nonzero_mutations)[hostile_label]
require("hostile deletion named residual", kequal(hostile_value, ka(q(-5265) / 22)))

print(
    f"N6_certificate_support=25 hash={plus_certificate_hash} "
    f"pairing_against_RHS12={ktext(target_pairing)} norm={qtext(knorm(target_pairing))}"
)
print(f"pairing_against_unit_RHS={ktext(unit_target_pairing)}")
print(f"N6_conjugate_certificate_hash={minus_certificate_hash}")
print("hostile_delete=(branch+1,h0,t6);nonzero_columns=47/252;(dy^dz)*w_residual=-5265/22")
print(
    "retained_degree7_before_differentiation_regression="
    f"[h^6]y_h@-1={ktext(expected_regression)}"
)
print(f"THM3683_specialization_cross_control={ktext(sixth_debt)};certificate_scale=-3645/11")
print("PASS arbitrary_two_form_survives_total_order5_but_fails_total_order6")
print("PASS separation=actual_pair_only_through_J2;order5_survivor_arbitrary_two_form_only")

source_text = Path(__file__).read_text(encoding="utf-8")
assertion_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source_text)))
require("assertion-free AST", assertion_nodes == 0)
print(f"PASS ast_assertion_nodes={assertion_nodes}")
print("PASS scope=fixed_conjugate_degree7_folds_quadratic_fold_JC2_OPEN")
print(f"CHECKS={CHECKS}")
print("RESULT=PASS")
