"""Exact quartic-field J1=J2 lift via the 89-gap quotient compiler.

K[X] arithmetic is represented by four FLINT rational polynomials in the
power basis (1,alpha,alpha^2,alpha^3).  This avoids slow symbolic convolution
in the number field while retaining an exact rational verification path.
"""

import contextlib
import gc
import hashlib
import importlib.util
import io
import math
import os
from pathlib import Path
import subprocess
import sys
import time
import types

import sympy as sp
from flint import fmpq, fmpq_mat, fmpq_poly, fmpz_mat


J0_PATH = Path("04-computation/jc2_russell_cylinder_exceptional_quartic_exact_j0_lift_thm3687.py")
START_TIME = time.monotonic()


def progress(label):
    if os.environ.get("EXACT_COUPLED_PROGRESS"):
        print(f"progress={label};elapsed={time.monotonic() - START_TIME:.1f}s", file=sys.stderr, flush=True)


def load_module(path, name, quiet=False):
    if path.is_file():
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        runner = lambda: spec.loader.exec_module(module)
    else:
        source = subprocess.check_output(
            ["git", "show", f"HEAD:{path.as_posix()}"], text=True
        )
        module = types.ModuleType(name)
        runner = lambda: exec(compile(source, path.as_posix(), "exec"), module.__dict__)
    if quiet:
        with contextlib.redirect_stdout(io.StringIO()):
            runner()
    else:
        runner()
    return module


exact = load_module(J0_PATH, "quartic_exact_j0", quiet=True)
progress("J0-replay")


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def as_fmpq(value):
    value = sp.Rational(value)
    return fmpq(int(value.p), int(value.q))


ZERO_Q = fmpq(0)
ONE_Q = fmpq(1)
ZERO_POLY = fmpq_poly([])
RELATION = (
    fmpq(1276420, 72783360),
    fmpq(-7849770, 72783360),
    fmpq(28419741, 72783360),
    fmpq(77822208, 72783360),
)


def elem_zero():
    return (ZERO_Q, ZERO_Q, ZERO_Q, ZERO_Q)


def elem_one():
    return (ONE_Q, ZERO_Q, ZERO_Q, ZERO_Q)


def elem_add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def elem_neg(value):
    return tuple(-item for item in value)


def elem_sub(left, right):
    return tuple(a - b for a, b in zip(left, right))


def elem_mul(left, right):
    temporary = [ZERO_Q] * 7
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


def elem_from_rational_rows(matrix, row, column):
    return tuple(matrix[4 * row + coordinate, column] for coordinate in range(4))


def multiplication_block(value):
    basis = (
        (ONE_Q, ZERO_Q, ZERO_Q, ZERO_Q),
        (ZERO_Q, ONE_Q, ZERO_Q, ZERO_Q),
        (ZERO_Q, ZERO_Q, ONE_Q, ZERO_Q),
        (ZERO_Q, ZERO_Q, ZERO_Q, ONE_Q),
    )
    columns = [elem_mul(value, item) for item in basis]
    return [[columns[column][row] for column in range(4)] for row in range(4)]


def elem_inv(value):
    require(any(value), "field inverse of zero")
    matrix = fmpq_mat(4, 4, [entry for row in multiplication_block(value) for entry in row])
    solution = matrix.solve(fmpq_mat(4, 1, [ONE_Q, ZERO_Q, ZERO_Q, ZERO_Q]))
    answer = tuple(solution[row, 0] for row in range(4))
    require(elem_mul(value, answer) == elem_one(), "field inverse residual")
    return answer


def kpoly_zero():
    return (ZERO_POLY, ZERO_POLY, ZERO_POLY, ZERO_POLY)


def kpoly_add(left, right):
    return tuple(a + b for a, b in zip(left, right))


def kpoly_neg(value):
    return tuple(-item for item in value)


def kpoly_sub(left, right):
    return tuple(a - b for a, b in zip(left, right))


def kpoly_mul(left, right):
    temporary = [fmpq_poly([]) for _ in range(7)]
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


def kpoly_mul_elem(polynomial, element):
    constant = tuple(fmpq_poly([item]) for item in element)
    return kpoly_mul(polynomial, constant)


def kpoly_scale(polynomial, scalar):
    scalar = fmpq(scalar)
    return tuple(item * scalar for item in polynomial)


def kpoly_diff(polynomial):
    return tuple(item.derivative() for item in polynomial)


def kpoly_shift(polynomial, degree):
    return tuple(item.left_shift(degree) for item in polynomial)


def kpoly_coeff(polynomial, degree):
    return tuple(item[degree] if degree <= item.degree() else ZERO_Q for item in polynomial)


def kpoly_degree(polynomial):
    degrees = [item.degree() for item in polynomial if item]
    return max(degrees) if degrees else -1


def kpoly_is_zero(polynomial):
    return not any(polynomial)


def kpoly_divmod(numerator, denominator):
    """Euclidean division in K[X], with an exact reconstruction gate."""
    require(not kpoly_is_zero(denominator), "K[X] division by zero")
    denominator_degree = kpoly_degree(denominator)
    denominator_lc_inverse = elem_inv(kpoly_coeff(denominator, denominator_degree))
    quotient = kpoly_zero()
    remainder = numerator
    while not kpoly_is_zero(remainder) and kpoly_degree(remainder) >= denominator_degree:
        shift = kpoly_degree(remainder) - denominator_degree
        coefficient = elem_mul(
            kpoly_coeff(remainder, kpoly_degree(remainder)), denominator_lc_inverse
        )
        term = kpoly_shift(kpoly_constant(coefficient), shift)
        quotient = kpoly_add(quotient, term)
        remainder = kpoly_sub(
            remainder,
            kpoly_shift(kpoly_mul_elem(denominator, coefficient), shift),
        )
    # Every caller immediately substitutes the quotient into either a full
    # L0 identity or the full high-J2 residual.  Those end-to-end K[X] gates
    # certify the division more cheaply than remultiplying enormous quotient
    # coefficients here.
    return quotient, remainder


def from_sympy_poly(polynomial):
    if not polynomial:
        return kpoly_zero()
    degree = polynomial.degree()
    rows = [[ZERO_Q] * (degree + 1) for _ in range(4)]
    for power in range(degree + 1):
        coordinates = exact.field_coordinates(
            polynomial.get((power,), exact.field.zero)
        )
        for coordinate, value in enumerate(coordinates):
            rows[coordinate][power] = as_fmpq(value)
    return tuple(fmpq_poly(row) for row in rows)


B = from_sympy_poly(exact.B)
C = from_sympy_poly(exact.C)
E = from_sympy_poly(exact.E)
dC = from_sympy_poly(exact.dC)
dE = from_sympy_poly(exact.dE)
F1 = from_sympy_poly(exact.F1)
G1 = from_sympy_poly(exact.G1)
delta_F1 = from_sympy_poly(exact.delta_F1)
delta_G1 = from_sympy_poly(exact.delta_G1)
C_prime = kpoly_diff(C)
E_prime = kpoly_diff(E)
F1_prime = kpoly_diff(F1)
G1_prime = kpoly_diff(G1)


four_poly = tuple(fmpq_poly([item]) for item in (fmpq(4), ZERO_Q, ZERO_Q, ZERO_Q))
require(
    kpoly_is_zero(
        kpoly_sub(
            kpoly_mul(kpoly_mul(C, C), E),
            kpoly_mul(B, kpoly_add(B, four_poly)),
        )
    ),
    "native quartic-field compiler relation",
)
require((kpoly_degree(B), kpoly_degree(C), kpoly_degree(E)) == (30, 21, 18), "native generator degrees")


def kpoly_constant(value):
    return tuple(fmpq_poly([item]) for item in value)


def kpoly_monic(polynomial, label):
    require(not kpoly_is_zero(polynomial), f"canonical nonzero {label}")
    answer = kpoly_mul_elem(
        polynomial,
        elem_inv(kpoly_coeff(polynomial, kpoly_degree(polynomial))),
    )
    require(kpoly_coeff(answer, kpoly_degree(answer)) == elem_one(), f"canonical monic {label}")
    return answer


Z = kpoly_monic(kpoly_add(E, kpoly_constant((fmpq(3), ZERO_Q, ZERO_Q, ZERO_Q))), "Z")
C0 = kpoly_monic(C, "C0")
B0 = kpoly_monic(B, "B0")
G71 = kpoly_monic(kpoly_sub(kpoly_mul(kpoly_mul(kpoly_mul(Z, Z), Z), Z), kpoly_mul(kpoly_mul(C0, C0), B0)), "G71")
G83 = kpoly_monic(kpoly_sub(kpoly_mul(kpoly_mul(kpoly_mul(C0, C0), C0), C0), kpoly_mul(kpoly_mul(kpoly_mul(Z, Z), Z), B0)), "G83")
G124 = kpoly_monic(kpoly_sub(kpoly_mul(kpoly_mul(kpoly_mul(G71, Z), Z), Z), kpoly_mul(kpoly_mul(G83, C0), C0)), "G124")
module_basis = (
    kpoly_constant(elem_one()),
    kpoly_mul(C0, G124),
    kpoly_mul(C0, G71),
    C0,
    kpoly_mul(G83, G83),
    kpoly_mul(B0, G83),
    kpoly_mul(C0, C0),
    kpoly_mul(C0, kpoly_mul(G83, G83)),
    kpoly_mul(kpoly_mul(C0, B0), G83),
    kpoly_mul(kpoly_mul(C0, C0), C0),
    kpoly_mul(G71, G83),
    G83,
    B0,
    kpoly_mul(kpoly_mul(C0, G71), G83),
    kpoly_mul(C0, G83),
    kpoly_mul(C0, B0),
    G124,
    G71,
)
canonical_items = []
for residue, basis in enumerate(module_basis):
    power = 0
    polynomial = basis
    while kpoly_degree(polynomial) <= 375:
        canonical_items.append(((residue, power), polynomial))
        power += 1
        polynomial = kpoly_mul(polynomial, Z)
canonical_items.sort(key=lambda item: kpoly_degree(item[1]))
monomials = [item[0] for item in canonical_items]
restrictions = {index: item[1] for index, item in enumerate(canonical_items)}
size = len(monomials)
require(size == 287, "canonical packet size")
require(
    [kpoly_degree(item[1]) for item in canonical_items]
    == sorted(kpoly_degree(item[1]) for item in canonical_items),
    "canonical degree order",
)
progress("canonical-packet")


known_j1 = kpoly_sub(
    kpoly_add(
        kpoly_scale(kpoly_mul(C_prime, dE), 2),
        kpoly_mul(F1_prime, G1),
    ),
    kpoly_add(
        kpoly_mul(F1, G1_prime),
        kpoly_scale(kpoly_mul(dC, E_prime), 2),
    ),
)
known_j2 = kpoly_add(
    kpoly_add(
        kpoly_scale(kpoly_mul(C_prime, delta_G1), 3),
        kpoly_scale(kpoly_mul(F1_prime, dE), 2),
    ),
    kpoly_add(
        kpoly_mul(kpoly_diff(dC), G1),
        kpoly_neg(kpoly_mul(F1, kpoly_diff(dE))),
    ),
)
known_j2 = kpoly_add(
    known_j2,
    kpoly_add(
        kpoly_scale(kpoly_mul(dC, G1_prime), -2),
        kpoly_scale(kpoly_mul(delta_F1, E_prime), -3),
    ),
)
require((kpoly_degree(known_j1), kpoly_degree(known_j2)) == (392, 207), "known debt degrees")


def l0(F, G):
    return kpoly_sub(kpoly_mul(C_prime, G), kpoly_mul(E_prime, F))


def m_operator(F, G):
    return kpoly_add(
        kpoly_sub(
            kpoly_mul(kpoly_diff(F), G1),
            kpoly_scale(kpoly_mul(F, G1_prime), 2),
        ),
        kpoly_sub(
            kpoly_scale(kpoly_mul(F1_prime, G), 2),
            kpoly_mul(F1, kpoly_diff(G)),
        ),
    )


# THM-3703 gives one monic canonical restriction in each non-gap degree.
# Descending reduction therefore exposes the exact 89-dimensional quotient
# K[X]/S without a Groebner or raw target-monomial solve.
canonical_by_degree = {
    kpoly_degree(polynomial): (index, polynomial)
    for index, polynomial in restrictions.items()
}
require(len(canonical_by_degree) == size, "one canonical basis element per degree")
gap_degrees = [degree for degree in range(376) if degree not in canonical_by_degree]
require(len(gap_degrees) == 89 and gap_degrees[-1] == 169, "exact 89-gap quotient")


def canonical_reduce(polynomial, want_coordinates=False):
    require(kpoly_degree(polynomial) <= 375, "canonical reduction cutoff")
    remainder = polynomial
    coordinates = [elem_zero() for _ in range(size)]
    for degree in range(kpoly_degree(remainder), -1, -1):
        if degree not in canonical_by_degree:
            continue
        coefficient = kpoly_coeff(remainder, degree)
        if any(coefficient):
            index, basis = canonical_by_degree[degree]
            coordinates[index] = elem_add(coordinates[index], coefficient)
            remainder = kpoly_sub(remainder, kpoly_mul_elem(basis, coefficient))
    if want_coordinates:
        return coordinates, remainder
    return remainder


def membership_vector(F, G):
    F_remainder = canonical_reduce(F)
    G_remainder = canonical_reduce(G)
    return [kpoly_coeff(F_remainder, degree) for degree in gap_degrees] + [
        kpoly_coeff(G_remainder, degree) for degree in gap_degrees
    ]


def vector_polynomial(vector):
    return tuple(
        fmpq_poly([value[coordinate] for value in vector])
        for coordinate in range(4)
    )


# Kernel pairs are (C'U,E'U).  For deg(U)<=355, membership in S^2 is a
# 178-by-356 gap problem of exact rank 177.  The first 177 monomials and first
# 177 quotient coordinates form the canonical square; U-degrees 177..355 are
# its 179 free directions.  This is the quotient compiler behind the former
# 395-row raw restriction solve.
membership_columns = {}
for degree in range(177):
    membership_columns[degree] = vector_polynomial(
        membership_vector(kpoly_shift(C_prime, degree), kpoly_shift(E_prime, degree))
    )
    if degree % 60 == 0:
        progress(f"gap-column-{degree}")
membership_pivots = list(range(177))
membership_rows = list(range(177))
omitted_membership_row = 177


def integer_square_solve(row_indices, column_indices, columns, vectors):
    """Solve a K-square by expanding K/Q and clearing left rows exactly."""
    require(len(row_indices) == len(column_indices), "integer square K shape")
    require(vectors and all(len(vector) == len(row_indices) for vector in vectors), "integer RHS shape")
    dimension = 4 * len(row_indices)
    rhs_width = len(vectors)
    square = fmpz_mat(dimension, dimension)
    scaled_rhs = fmpq_mat(dimension, rhs_width)
    for local_row, quotient_row in enumerate(row_indices):
        if local_row % 60 == 0:
            progress(f"gap-square-row-{local_row}")
        blocks = [
            multiplication_block(kpoly_coeff(columns[column], quotient_row))
            for column in column_indices
        ]
        for coordinate_row in range(4):
            expanded_row = 4 * local_row + coordinate_row
            left = [entry for block in blocks for entry in block[coordinate_row]]
            right = [vector[local_row][coordinate_row] for vector in vectors]
            denominator = 1
            for value in left:
                denominator = denominator // math.gcd(denominator, int(value.denominator)) * int(value.denominator)
            for column, value in enumerate(left):
                square[expanded_row, column] = int(value.numerator) * (
                    denominator // int(value.denominator)
                )
            for column, value in enumerate(right):
                scaled_rhs[expanded_row, column] = value * denominator
    rhs_column_denominators = []
    for column in range(rhs_width):
        denominator = 1
        for row in range(dimension):
            denominator = denominator // math.gcd(
                denominator, int(scaled_rhs[row, column].denominator)
            ) * int(scaled_rhs[row, column].denominator)
        rhs_column_denominators.append(denominator)
    rhs = fmpz_mat(dimension, rhs_width)
    for row in range(dimension):
        for column, denominator in enumerate(rhs_column_denominators):
            value = scaled_rhs[row, column]
            rhs[row, column] = int(value.numerator) * (
                denominator // int(value.denominator)
            )
    solution = square.solve(rhs)
    for row in range(dimension):
        for column, denominator in enumerate(rhs_column_denominators):
            solution[row, column] /= denominator
    # The constructed pairs are substituted into all 178 gap coordinates and
    # both full Jacobian polynomials below.  Those K[X] identities are a much
    # cheaper and stronger end-to-end certificate than multiplying this dense
    # denominator-cleared Q-square back by every intermediate RHS here.
    del square, rhs, scaled_rhs
    gc.collect()
    return solution


def bezout_particular(right_hand_side):
    quotient, F = kpoly_divmod(kpoly_mul(F1, right_hand_side), C_prime)
    U = kpoly_neg(quotient)
    G = kpoly_add(kpoly_mul(G1, right_hand_side), kpoly_mul(E_prime, U))
    require(l0(F, G) == right_hand_side, "Bezout particular residual")
    return F, G


def negative_membership_rhs(F, G):
    return [elem_neg(value) for value in membership_vector(F, G)[:177]]


def upoly_from_solution(solution, scenario, free_degree=None):
    coordinate_rows = [[ZERO_Q] * 356 for _ in range(4)]
    for row, degree in enumerate(membership_pivots):
        value = elem_from_rational_rows(solution, row, scenario)
        for coordinate in range(4):
            coordinate_rows[coordinate][degree] = value[coordinate]
    if free_degree is not None:
        coordinate_rows[0][free_degree] += ONE_Q
    return tuple(fmpq_poly(row) for row in coordinate_rows)


half = fmpq(1, 2)
h1 = kpoly_scale(known_j1, -half)
F2_particular, G2_particular = bezout_particular(h1)
selected_u_degrees = list(range(182, 356))
require(len(selected_u_degrees) == 174, "selected U directions")
solve_u_base = integer_square_solve(
    membership_rows,
    membership_pivots,
    membership_columns,
    [negative_membership_rhs(F2_particular, G2_particular)],
)
U_base = upoly_from_solution(solve_u_base, 0)
del solve_u_base
gc.collect()
F2_base = kpoly_add(F2_particular, kpoly_mul(C_prime, U_base))
G2_base = kpoly_add(G2_particular, kpoly_mul(E_prime, U_base))
require(kpoly_is_zero(canonical_reduce(F2_base)), "J1 base F2 membership")
require(kpoly_is_zero(canonical_reduce(G2_base)), "J1 base G2 membership")
require(l0(F2_base, G2_base) == h1, "J1 base L0 residual")
progress("J1-gap-solutions")


# On an L0-kernel pair `(C'U,E'U)`, the entire transverse operator collapses:
#
#     M(C'U,E'U)=U' + A U,
#     A=C''G1-2C'G1'+2F1'E'-F1E''.
#
# The U' coefficient is exactly `C'G1-F1E'=1`.  Moreover deg A=214 and
# lc(A)=162*lc(E)*lc(F1).  Thus monic U-degrees 182..355 give a triangular
# response on degrees 396..569 with one constant nonzero diagonal.
base_pre_j2 = kpoly_add(known_j2, m_operator(F2_base, G2_base))
kernel_potential = kpoly_add(
    kpoly_sub(
        kpoly_mul(kpoly_diff(C_prime), G1),
        kpoly_scale(kpoly_mul(C_prime, G1_prime), 2),
    ),
    kpoly_sub(
        kpoly_scale(kpoly_mul(F1_prime, E_prime), 2),
        kpoly_mul(F1, kpoly_diff(E_prime)),
    ),
)


require(kpoly_degree(kernel_potential) == 214, "kernel-potential degree")
expected_diagonal = tuple(
    fmpq(162) * item
    for item in elem_mul(
        kpoly_coeff(E, kpoly_degree(E)),
        kpoly_coeff(F1, kpoly_degree(F1)),
    )
)
require(
    kpoly_coeff(kernel_potential, 214) == expected_diagonal,
    "kernel-potential leading coefficient",
)

# On degrees >=396 the derivative U' and the membership corrections of
# degrees <=176 are invisible.  Therefore the entire 174-step triangular
# solve is just the high quotient of `-base_pre_j2` by A.  Long division also
# gates the three nominal rows 570..572: successful support has no quotient
# term above degree 355.
transverse_quotient, _transverse_remainder = kpoly_divmod(
    kpoly_neg(base_pre_j2), kernel_potential
)
require(kpoly_degree(transverse_quotient) <= 355, "transverse top-three compatibility")
parameters = [
    kpoly_coeff(transverse_quotient, degree) for degree in selected_u_degrees
]


free_coordinate_rows = [[ZERO_Q] * 356 for _ in range(4)]
for degree, parameter in zip(selected_u_degrees, parameters):
    for coordinate in range(4):
        free_coordinate_rows[coordinate][degree] = parameter[coordinate]
U_free = tuple(fmpq_poly(row) for row in free_coordinate_rows)
F2_seed = kpoly_add(F2_particular, kpoly_mul(C_prime, U_free))
G2_seed = kpoly_add(G2_particular, kpoly_mul(E_prime, U_free))
solve_u_final = integer_square_solve(
    membership_rows,
    membership_pivots,
    membership_columns,
    [negative_membership_rhs(F2_seed, G2_seed)],
)
U_final = kpoly_add(U_free, upoly_from_solution(solve_u_final, 0))
del solve_u_final
F2 = kpoly_add(F2_particular, kpoly_mul(C_prime, U_final))
G2 = kpoly_add(G2_particular, kpoly_mul(E_prime, U_final))
require(all(value == elem_zero() for value in membership_vector(F2, G2)), "final J1 target membership")
pre_j2 = kpoly_add(known_j2, m_operator(F2, G2))
require(
    all(kpoly_coeff(pre_j2, degree) == elem_zero() for degree in range(396, kpoly_degree(pre_j2) + 1)),
    "all high J2 rows vanish",
)
require(kpoly_degree(pre_j2) <= 395, "J2 reduced to L0 range")
progress("transverse-solve")


third = fmpq(1, 3)
h2 = kpoly_scale(pre_j2, -third)
F3_particular, G3_particular = bezout_particular(h2)
solve_v_base = integer_square_solve(
    membership_rows,
    membership_pivots,
    membership_columns,
    [negative_membership_rhs(F3_particular, G3_particular)],
)
V_base = upoly_from_solution(solve_v_base, 0)
del solve_v_base
F3 = kpoly_add(F3_particular, kpoly_mul(C_prime, V_base))
G3 = kpoly_add(G3_particular, kpoly_mul(E_prime, V_base))
final_membership = membership_vector(F3, G3)
require(final_membership[omitted_membership_row] == elem_zero(), "omitted gap compatibility")
require(all(value == elem_zero() for value in final_membership), "J2 target membership")
require(l0(F3, G3) == h2, "J2 L0 residual")

j1_residual = kpoly_add(known_j1, kpoly_scale(l0(F2, G2), 2))
j2_residual = kpoly_add(pre_j2, kpoly_scale(l0(F3, G3), 3))
require(kpoly_is_zero(j1_residual), "full J1 residual")
require(kpoly_is_zero(j2_residual), "full J2 residual")
progress("full-residuals")


F2_values, F2_remainder = canonical_reduce(F2, want_coordinates=True)
G2_values, G2_remainder = canonical_reduce(G2, want_coordinates=True)
F3_values, F3_remainder = canonical_reduce(F3, want_coordinates=True)
G3_values, G3_remainder = canonical_reduce(G3, want_coordinates=True)
require(
    all(kpoly_is_zero(value) for value in (F2_remainder, G2_remainder, F3_remainder, G3_remainder)),
    "canonical target reconstruction",
)


def elem_text(value):
    return ",".join(str(item) for item in value)


def elem_hash(value):
    return hashlib.sha256(elem_text(value).encode("ascii")).hexdigest()


def target_hash(values):
    payload = ";".join(
        f"{residue},{power}:{elem_text(value)}"
        for value, (residue, power) in zip(values, monomials)
        if any(value)
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def index_hash(values):
    return hashlib.sha256(",".join(map(str, values)).encode("ascii")).hexdigest()


print("field=THM3683_exceptional_quartic")
print("compiler=THM3703_89_gap_quotient;membership_K_square=177;membership_Q_square=708")
print(
    f"gap_hash={index_hash(gap_degrees)};membership_pivots_hash={index_hash(membership_pivots)};"
    f"selected_U_degrees_hash={index_hash(selected_u_degrees)}"
)
print(
    f"transverse_K_triangular=174;constant_diagonal_nonzero={int(any(expected_diagonal))};"
    f"constant_diagonal_hash={elem_hash(expected_diagonal)}"
)
for label, values, polynomial in (
    ("F2", F2_values, F2),
    ("G2", G2_values, G2),
    ("F3", F3_values, F3),
    ("G3", G3_values, G3),
):
    print(
        f"{label}_support={sum(any(value) for value in values)};"
        f"{label}_degree={kpoly_degree(polynomial)};{label}_target_hash={target_hash(values)}"
    )
print("J0=1;J1=0;J2=0;all_four_embeddings=uniform;full_quartic_field_residuals=0")
print("NO_CLAIM=J3_or_J4_or_global_pair_or_Keller_counterexample")
print("RESULT=PASS")
