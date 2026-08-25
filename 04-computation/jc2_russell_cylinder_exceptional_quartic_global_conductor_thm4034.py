"""Exact global conductor certificate for the THM-3703 restriction algebra.

The positive modular minor selects an exact square only.  The conductor
candidate is solved over the quartic field and every normalization-basis
product is reduced exactly against the free K[Z] Apéry basis.
"""

import contextlib
import gc
import hashlib
import importlib.util
import io
import math
import subprocess
import types
from pathlib import Path

from flint import fmpq, fmpz_mat, nmod_mat, nmod_poly
import sympy as sp

SAGBI_PATH = Path(
    "04-computation/jc2_russell_cylinder_exceptional_quartic_sagbi_module_thm3703.py"
)
PRIME = 137
ALPHA_MOD = 44
MAX_DEGREE = 195
MONIC_DEGREE = 178


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def load_module(path, name):
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
    with contextlib.redirect_stdout(io.StringIO()):
        runner()
    return module


sagbi = load_module(SAGBI_PATH, "quartic_sagbi_thm3703")

field = sagbi.field
alpha = sagbi.alpha
ring = sagbi.polynomial_ring
X = sagbi.X
Z = sagbi.Z
zero = field.zero
one = field.one
gaps = tuple(sagbi.EXPECTED_GAPS)
require(len(gaps) == 89 and gaps[-1] == 169, "THM-3703 gap boundary")
require(Z.degree() == 18 and Z.LC == one, "normalization basis degree")
require(72783360 % PRIME != 0, "good reduction quartic leading coefficient")
require(
    (
        72783360 * ALPHA_MOD**4
        - 77822208 * ALPHA_MOD**3
        - 28419741 * ALPHA_MOD**2
        + 7849770 * ALPHA_MOD
        - 1276420
    )
    % PRIME
    == 0,
    "good reduction quartic root",
)


def field_coordinates(value):
    return sagbi.field_coordinates(value)


def rational_mod(value):
    value = sp.Rational(value)
    require(int(value.q) % PRIME != 0, "good reduction denominator")
    return int(value.p) % PRIME * pow(int(value.q) % PRIME, -1, PRIME) % PRIME


def field_mod(value):
    return sum(
        rational_mod(coordinate) * pow(ALPHA_MOD, power, PRIME)
        for power, coordinate in enumerate(field_coordinates(value))
    ) % PRIME


canonical = {}
for basis in sagbi.module_basis:
    polynomial = basis
    while polynomial.degree() <= MAX_DEGREE:
        require(polynomial.degree() not in canonical, "unique canonical degree")
        canonical[polynomial.degree()] = polynomial
        polynomial *= Z
require(all(degree in canonical for degree in range(170, MAX_DEGREE + 1)), "semigroup tail")


def normal_form(polynomial):
    remainder = polynomial
    for degree in range(remainder.degree(), -1, -1):
        basis = canonical.get(degree)
        if basis is None:
            continue
        coefficient = remainder.get((degree,), zero)
        if coefficient:
            remainder -= coefficient * basis
    require(remainder.degree() <= 169 if remainder else True, "gap remainder degree")
    return remainder


powers = [ring.one]
for _ in range(MAX_DEGREE):
    powers.append(powers[-1] * X)
remainders = [normal_form(polynomial) for polynomial in powers]
require(all(normal_form(basis) == ring.zero for basis in canonical.values()), "basis membership")


def coefficient(polynomial, degree):
    return polynomial.get((degree,), zero)


row_keys = tuple((shift, gap) for shift in range(18) for gap in gaps)
matrix_data = []
for shift, gap in row_keys:
    matrix_data.extend(
        field_mod(coefficient(remainders[degree + shift], gap))
        for degree in range(MONIC_DEGREE)
    )
matrix_mod = nmod_mat(len(row_keys), MONIC_DEGREE, matrix_data, PRIME)
transpose_reduced, rank = matrix_mod.transpose().rref()
require(rank == MONIC_DEGREE, "positive full-column rank")
selected_rows = []
for row in range(rank):
    selected_rows.append(next(column for column in range(len(row_keys)) if transpose_reduced[row, column]))
require(len(set(selected_rows)) == MONIC_DEGREE, "distinct selected rows")
selected_mod = nmod_mat(
    MONIC_DEGREE,
    MONIC_DEGREE,
    [int(matrix_mod[row, column]) for row in selected_rows for column in range(MONIC_DEGREE)],
    PRIME,
)
det_mod = int(selected_mod.det())
require(det_mod != 0, "positive selected determinant")


alpha_powers = (one, alpha, alpha**2, alpha**3)


def multiplication_block(value):
    columns = [field_coordinates(value * power) for power in alpha_powers]
    return [[columns[column][row] for column in range(4)] for row in range(4)]


dimension = 4 * MONIC_DEGREE
square = fmpz_mat(dimension, dimension)
rhs_scaled = [None] * dimension
for local_row, selected_row in enumerate(selected_rows):
    shift, gap = row_keys[selected_row]
    blocks = [
        multiplication_block(coefficient(remainders[degree + shift], gap))
        for degree in range(MONIC_DEGREE)
    ]
    right_value = -coefficient(remainders[MONIC_DEGREE + shift], gap)
    right_coordinates = field_coordinates(right_value)
    for coordinate_row in range(4):
        expanded_row = 4 * local_row + coordinate_row
        left = [entry for block in blocks for entry in block[coordinate_row]]
        right = sp.Rational(right_coordinates[coordinate_row])
        denominator = int(right.q)
        for value in left:
            rational = sp.Rational(value)
            denominator = denominator // math.gcd(denominator, int(rational.q)) * int(rational.q)
        for column, value in enumerate(left):
            rational = sp.Rational(value)
            square[expanded_row, column] = int(rational.p) * (denominator // int(rational.q))
        rhs_scaled[expanded_row] = fmpq(int(right.p) * (denominator // int(right.q)), 1)
rhs_denominator = 1
for value in rhs_scaled:
    rhs_denominator = rhs_denominator // math.gcd(rhs_denominator, int(value.denominator)) * int(value.denominator)
rhs = fmpz_mat(dimension, 1)
for row, value in enumerate(rhs_scaled):
    rhs[row, 0] = int(value.numerator) * (rhs_denominator // int(value.denominator))
solution = square.solve(rhs)
for row in range(dimension):
    solution[row, 0] /= rhs_denominator


def field_from_solution(index):
    answer = zero
    for coordinate, power in enumerate(alpha_powers):
        value = solution[4 * index + coordinate, 0]
        answer += field.convert(sp.Rational(int(value.numerator), int(value.denominator))) * power
    return answer


conductor = X**MONIC_DEGREE
for degree in range(MONIC_DEGREE):
    value = field_from_solution(degree)
    if value:
        conductor += value * X**degree
del square, rhs, solution
gc.collect()
require(conductor.degree() == MONIC_DEGREE and conductor.LC == one, "monic conductor candidate")
for shift in range(18):
    require(normal_form(conductor * X**shift) == ring.zero, f"conductor membership {shift}")

retained_factor = X**2 * (X**2 - 1) ** 2
quotient, remainder = divmod(conductor, retained_factor)
require(not remainder and quotient.degree() == 172, "retained triple factor")
for point in (-1, 0, 1):
    require(quotient.evaluate(X, field.convert(point)) != zero, f"retained multiplicity {point}")


def polynomial_mod(polynomial):
    return nmod_poly(
        [field_mod(coefficient(polynomial, degree)) for degree in range(polynomial.degree() + 1)],
        PRIME,
    )


quotient_mod = polynomial_mod(quotient)
retained_factor_mod = polynomial_mod(retained_factor)
require(quotient_mod.gcd(quotient_mod.derivative()).degree() == 0, "squarefree quotient reduction")
require(quotient_mod.gcd(retained_factor_mod).degree() == 0, "retained-factor coprimality reduction")

# The three divided-difference resultants give an intrinsic conductor formula.
# The exact proof uses trace duals to place every resultant in the conductor;
# this good reduction supplies the sharp common-divisor degree bound.  Keeping
# all three pairs is essential: the first two alone retain a degree-eight
# spurious factor at the distinguished triple fibre.
x_symbol, y_symbol = sp.symbols("x y")


def modular_expression(polynomial, variable):
    return sum(
        int(polynomial[degree]) * variable**degree
        for degree in range(polynomial.degree() + 1)
    )


def divided_difference(polynomial):
    numerator = sp.Poly(
        modular_expression(polynomial, x_symbol)
        - modular_expression(polynomial, y_symbol),
        x_symbol,
        y_symbol,
        modulus=PRIME,
    )
    quotient_polynomial, remainder_polynomial = sp.div(
        numerator,
        sp.Poly(x_symbol - y_symbol, x_symbol, y_symbol, modulus=PRIME),
    )
    require(remainder_polynomial.is_zero, "divided-difference remainder")
    return quotient_polynomial.as_expr()


B_mod, C_mod, E_mod = map(polynomial_mod, (sagbi.B, sagbi.C, sagbi.E))
delta_B, delta_C, delta_E = map(divided_difference, (B_mod, C_mod, E_mod))
resultants = tuple(
    sp.Poly(sp.resultant(left, right, y_symbol), x_symbol, modulus=PRIME).monic()
    for left, right in (
        (delta_B, delta_C),
        (delta_B, delta_E),
        (delta_C, delta_E),
    )
)
require(tuple(item.degree() for item in resultants) == (578, 488, 338), "resultant degrees")
triple_gcd = resultants[0].gcd(resultants[1]).gcd(resultants[2]).monic()
conductor_sympy = sp.Poly(
    sum(
        field_mod(coefficient(conductor, degree)) * x_symbol**degree
        for degree in range(conductor.degree() + 1)
    ),
    x_symbol,
    modulus=PRIME,
).monic()
require(triple_gcd == conductor_sympy, "triple-resultant conductor reduction")
pair_gcd = resultants[0].gcd(resultants[1]).monic()
pair_extra, pair_remainder = sp.div(pair_gcd, triple_gcd)
require(pair_remainder.is_zero, "pair-resultant quotient")
require(
    pair_extra.monic()
    == sp.Poly(x_symbol**2 * (x_symbol**2 - 1) ** 3, x_symbol, modulus=PRIME).monic(),
    "two-resultant hostile factor",
)


hostile_failed_shifts = tuple(
    shift for shift in range(18) if remainders[170 + shift]
)
require(hostile_failed_shifts, "x^170 is not a conductor generator")
require(bool(remainders[170]) and hostile_failed_shifts[0] == 0, "x^170 is not in S")


def index_hash(values):
    return hashlib.sha256(",".join(map(str, values)).encode("ascii")).hexdigest()


def coefficient_hash(polynomial):
    payload = ";".join(
        f"{degree}:" + ",".join(str(item) for item in field_coordinates(coefficient(polynomial, degree)))
        for degree in range(polynomial.degree() + 1)
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


print("field=THM3683_exceptional_quartic;restriction_algebra=THM3703")
print(
    f"condition_shape={len(row_keys)}x{MONIC_DEGREE};rank_mod_{PRIME}={rank};"
    f"selected_det={det_mod};selected_rows_hash={index_hash(selected_rows)}"
)
print(
    f"conductor_degree={conductor.degree()};conductor_support="
    f"{sum(bool(coefficient(conductor, degree)) for degree in range(conductor.degree() + 1))};"
    f"conductor_hash={coefficient_hash(conductor)}"
)
print(
    f"conductor_products=18;all_exact_membership=True;"
    f"normalization_basis=1,x,...,x^17_over_K[Z]"
)
print(
    f"retained_factor=x^2*(x^2-1)^2;quotient_degree={quotient.degree()};"
    f"quotient_hash={coefficient_hash(quotient)};retained_multiplicity_exact=2,2,2;"
    f"off_retained_squarefree=True"
)
print(
    f"degree_delta=89;global_conductor_colength=178;S_mod_conductor_length=89;"
    f"numerical_semigroup_conductor=170;x170_first_failed_shift={hostile_failed_shifts[0]}"
)
print("resultant_degrees=578,488,338;triple_gcd_degree=178;intrinsic_conductor=True")
print("two_resultant_hostile_mod137_degree=186;extra=x^2*(x^2-1)^3")
print("NO_CLAIM=singularity_classification_or_higher_stable_lift_or_Keller_or_JC2")
print("RESULT=PASS")
