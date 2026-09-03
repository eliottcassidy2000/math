"""Exact scratch scout for the THM-4381 seminormal class.

This is deliberately not a theorem companion.  It computes the image in
N/S of the fixed-x first normal sidecar of the Russell target compiler.  The
map is representative-level and is not identified here with THM-4067's
fixed-c first-output derivative.
"""

import contextlib
import hashlib
import io
import subprocess
import types
from pathlib import Path


CONDUCTOR_PATH = Path(
    "04-computation/jc2_russell_cylinder_exceptional_quartic_global_conductor_thm4034.py"
)
CONDUCTOR_CORE_MARKER = (
    "# The three divided-difference resultants give an intrinsic conductor formula."
)


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def load_conductor_core(path, name):
    if path.is_file():
        source = path.read_text()
    else:
        source = subprocess.check_output(
            ["git", "show", f"HEAD:{path.as_posix()}"], text=True
        )
    require(source.count(CONDUCTOR_CORE_MARKER) == 1, "unique core marker")
    source = source.split(CONDUCTOR_CORE_MARKER, 1)[0]
    module = types.ModuleType(name)
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(source, path.as_posix(), "exec"), module.__dict__)
    return module


certificate = load_conductor_core(CONDUCTOR_PATH, "quartic_conductor_s615")
sagbi = certificate.sagbi
field = certificate.field
ring = certificate.ring
X = certificate.X
one = certificate.one
zero = certificate.zero
Z = sagbi.Z
Q = sagbi.Q
B = sagbi.B
C = sagbi.C
E = sagbi.E
gaps = tuple(sagbi.EXPECTED_GAPS)
require(len(gaps) == 89, "normalization quotient dimension")


def coefficient(polynomial, degree):
    return polynomial.get((degree,), zero)


def field_coordinates(value):
    return certificate.field_coordinates(value)


def serialize_field(value):
    return ",".join(str(item) for item in field_coordinates(value))


def vector_hash(vectors):
    payload = ";".join(
        ":".join(serialize_field(value) for value in vector) for vector in vectors
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def integer_vector_hash(vectors):
    payload = ";".join(":".join(map(str, vector)) for vector in vectors)
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


# Build the monic filtered S-basis far enough for all products below.  This
# is the same Apery reduction as THM-3703/4034, extended from degree 195 to
# the actual product ceiling rather than silently truncating it.
conductor = certificate.conductor
h172 = certificate.quotient
L = X * (X**2 - 1)
r = L * h172
coefficient_basis = tuple(
    certificate.canonical[degree]
    for degree in sorted(certificate.canonical)
    if degree < conductor.degree()
)
require(len(coefficient_basis) == 89, "S/conductor coefficient basis")

D = 1 + X**2 * Q
normal_generators = (
    3 * X**2 * D * (D + 2),
    2 * X**3 * (D + 1),
    2 * (D + 1),
)
generator_labels = ("gamma1(B)", "gamma1(C)", "gamma1(E)")
require(
    2 * C * E * normal_generators[1] + C**2 * normal_generators[2]
    == (2 * B + 4) * normal_generators[0],
    "normal derivative respects the Russell relation",
)

product_ceiling = max(
    basis.degree() + generator.degree()
    for basis in coefficient_basis
    for generator in normal_generators
)
canonical = {}
for basis in sagbi.module_basis:
    polynomial = basis
    while polynomial.degree() <= product_ceiling:
        require(polynomial.degree() not in canonical, "unique filtered degree")
        canonical[polynomial.degree()] = polynomial
        polynomial *= Z
require(
    all(degree in canonical for degree in range(170, product_ceiling + 1)),
    "filtered tail coverage",
)


def normal_form(polynomial):
    remainder = polynomial
    for degree in range(remainder.degree(), -1, -1):
        basis = canonical.get(degree)
        if basis is None:
            continue
        value = coefficient(remainder, degree)
        if value:
            remainder -= value * basis
    require(not remainder or remainder.degree() <= 169, "gap normal-form boundary")
    return remainder


def checked_gap_vector(polynomial):
    remainder = normal_form(polynomial)
    gap_set = set(gaps)
    require(
        all(
            coefficient(remainder, degree) == zero
            for degree in range(170)
            if degree not in gap_set
        ),
        "normal form supported on gaps",
    )
    return tuple(coefficient(remainder, gap) for gap in gaps)


columns = []
column_labels = []
for generator_label, generator in zip(generator_labels, normal_generators):
    for index, basis in enumerate(coefficient_basis):
        columns.append(checked_gap_vector(basis * generator))
        column_labels.append(f"{generator_label}*s{index}")


PRIME = certificate.PRIME


def vector_mod(vector):
    return tuple(certificate.field_mod(value) for value in vector)


def independent_rows_mod(rows):
    echelon = []
    pivots = []
    selected = []
    for source_index, source_row in enumerate(rows):
        row = list(vector_mod(source_row))
        for pivot, old_row in zip(pivots, echelon):
            value = row[pivot]
            if value:
                row = [
                    (left - value * right) % PRIME
                    for left, right in zip(row, old_row)
                ]
        pivot = next((index for index, value in enumerate(row) if value), None)
        if pivot is None:
            continue
        inverse = pow(row[pivot], -1, PRIME)
        row = [(value * inverse) % PRIME for value in row]
        pivots.append(pivot)
        echelon.append(row)
        selected.append(source_index)
    return tuple(echelon), tuple(pivots), tuple(selected)


r_vector = checked_gap_vector(r)
require(any(r_vector), "THM-4381 class is nonzero in N/S")
echelon, pivots, selected = independent_rows_mod(columns)
rank = len(echelon)
augmented_rank = len(independent_rows_mod(columns + [r_vector])[0])

# Individual generator channels show which first normal derivative already
# sees how much of N/S.  The retained functional is also printed as the
# local triple-point coordinate, but it is not used to infer global equality.
individual_ranks = []
individual_augmented_ranks = []
for offset in range(3):
    block = columns[offset * 89 : (offset + 1) * 89]
    individual_ranks.append(len(independent_rows_mod(block)[0]))
    individual_augmented_ranks.append(len(independent_rows_mod(block + [r_vector])[0]))

points = (-1, 0, 1)
weights = (field.convert(5) / field.convert(18), -one, field.convert(13) / field.convert(18))


def evaluate(polynomial, point):
    return polynomial.evaluate(X, field.convert(point))


def retained_lambda(polynomial):
    return sum(
        (weight * evaluate(polynomial, point) for weight, point in zip(weights, points)),
        zero,
    )


retained_values = tuple(
    tuple(evaluate(generator, point) for point in points)
    for generator in normal_generators
)
retained_responses = tuple(retained_lambda(generator) for generator in normal_generators)
r_derivative_response = retained_lambda(r.diff(X))
require(bool(r_derivative_response), "r spans retained first-jet defect")


# The correctly differentiated special-fibre quotient is
#
#     E_x = N / dS,       d=d/dx,
#
# not N/S.  Its monic complement has degrees g-1 for the 89 gaps g of S.
# Search actual target monomials H(b,c,e) for the positive statement
# r' in dS + gamma1(H).  A full-rank good reduction is a sound
# characteristic-zero certificate; failure to reach full modular rank would
# not be used as a characteristic-zero obstruction.
derivative_gaps = tuple(gap - 1 for gap in gaps)
require(len(set(derivative_gaps)) == 89, "derivative quotient rows")

TARGET_CUTOFF = 240
DERIVATIVE_PRIME = 421
DERIVATIVE_ALPHA = 126
require(DERIVATIVE_PRIME > TARGET_CUTOFF + 1, "derivative good prime exceeds degrees")
require(
    (
        72783360 * DERIVATIVE_ALPHA**4
        - 77822208 * DERIVATIVE_ALPHA**3
        - 28419741 * DERIVATIVE_ALPHA**2
        + 7849770 * DERIVATIVE_ALPHA
        - 1276420
    )
    % DERIVATIVE_PRIME
    == 0,
    "derivative reduction is a quartic root",
)


def rational_mod_at(value, prime):
    numerator, denominator = map(int, value.as_numer_denom())
    require(denominator % prime != 0, "good derivative denominator")
    return numerator % prime * pow(denominator % prime, -1, prime) % prime


def field_mod_at(value, prime, alpha_value):
    return sum(
        rational_mod_at(coordinate, prime) * pow(alpha_value, power, prime)
        for power, coordinate in enumerate(field_coordinates(value))
    ) % prime


def polynomial_mod_at(polynomial, prime, alpha_value):
    return certificate.nmod_poly(
        [
            field_mod_at(coefficient(polynomial, degree), prime, alpha_value)
            for degree in range(polynomial.degree() + 1)
        ],
        prime,
    )


def coefficient_mod(polynomial, degree):
    return int(polynomial[degree]) if degree <= polynomial.degree() else 0


B_mod = polynomial_mod_at(B, DERIVATIVE_PRIME, DERIVATIVE_ALPHA)
C_mod = polynomial_mod_at(C, DERIVATIVE_PRIME, DERIVATIVE_ALPHA)
E_mod = polynomial_mod_at(E, DERIVATIVE_PRIME, DERIVATIVE_ALPHA)
Z_mod = polynomial_mod_at(Z, DERIVATIVE_PRIME, DERIVATIVE_ALPHA)
r_mod = polynomial_mod_at(r, DERIVATIVE_PRIME, DERIVATIVE_ALPHA)
normal_generators_mod = tuple(
    polynomial_mod_at(generator, DERIVATIVE_PRIME, DERIVATIVE_ALPHA)
    for generator in normal_generators
)
require(
    tuple(generator.degree() for generator in (B_mod, C_mod, E_mod)) == (30, 21, 18),
    "target degrees survive derivative reduction",
)

target_powers_mod = []
for generator, degree in zip((B_mod, C_mod, E_mod), (30, 21, 18)):
    powers = [certificate.nmod_poly([1], DERIVATIVE_PRIME)]
    for _ in range(TARGET_CUTOFF // degree):
        powers.append(powers[-1] * generator)
    target_powers_mod.append(tuple(powers))

target_records_mod = []
for i, B_power in enumerate(target_powers_mod[0]):
    for j, C_power in enumerate(target_powers_mod[1]):
        for k, E_power in enumerate(target_powers_mod[2]):
            degree = 30 * i + 21 * j + 18 * k
            if degree > TARGET_CUTOFF:
                continue
            normal = certificate.nmod_poly([], DERIVATIVE_PRIME)
            if i:
                normal += i * target_powers_mod[0][i - 1] * C_power * E_power * normal_generators_mod[0]
            if j:
                normal += j * B_power * target_powers_mod[1][j - 1] * E_power * normal_generators_mod[1]
            if k:
                normal += k * B_power * C_power * target_powers_mod[2][k - 1] * normal_generators_mod[2]
            target_records_mod.append((degree, f"B^{i}C^{j}E^{k}", normal))
target_records_mod.sort(key=lambda item: (item[0], item[1]))
max_density_degree = max(item[2].degree() for item in target_records_mod if item[2])
require(max_density_degree < DERIVATIVE_PRIME - 1, "derivative leading coefficients are units")

derivative_canonical_mod = {}
for basis in sagbi.module_basis:
    polynomial = polynomial_mod_at(basis, DERIVATIVE_PRIME, DERIVATIVE_ALPHA)
    while polynomial.degree() <= max_density_degree + 1:
        if polynomial.degree() > 0:
            degree_before = polynomial.degree()
            derivative = polynomial.derivative() * pow(
                degree_before, -1, DERIVATIVE_PRIME
            )
            require(derivative.degree() == degree_before - 1, "derivative degree survives")
            require(coefficient_mod(derivative, derivative.degree()) == 1, "monic derivative basis")
            degree = derivative.degree()
            require(degree not in derivative_canonical_mod, "unique derivative degree")
            derivative_canonical_mod[degree] = derivative
        polynomial *= Z_mod
require(
    all(
        degree in derivative_canonical_mod
        for degree in range(169, max_density_degree + 1)
    ),
    "derivative filtered tail coverage",
)


def derivative_normal_form_mod(polynomial):
    remainder = polynomial
    for degree in range(remainder.degree(), -1, -1):
        basis = derivative_canonical_mod.get(degree)
        if basis is None:
            continue
        value = coefficient_mod(remainder, degree)
        if value:
            remainder -= value * basis
    require(
        not remainder or remainder.degree() <= derivative_gaps[-1],
        "derivative normal-form boundary",
    )
    derivative_gap_set = set(derivative_gaps)
    require(
        all(
            coefficient_mod(remainder, degree) == 0
            for degree in range(derivative_gaps[-1] + 1)
            if degree not in derivative_gap_set
        ),
        "derivative normal form supported on shifted gaps",
    )
    return remainder


def derivative_gap_vector_mod(polynomial):
    remainder = derivative_normal_form_mod(polynomial)
    return tuple(coefficient_mod(remainder, degree) for degree in derivative_gaps)


def independent_integer_rows_mod(rows, prime):
    echelon = []
    pivots = []
    selected = []
    for source_index, source_row in enumerate(rows):
        row = [value % prime for value in source_row]
        for pivot, old_row in zip(pivots, echelon):
            value = row[pivot]
            if value:
                row = [
                    (left - value * right) % prime
                    for left, right in zip(row, old_row)
                ]
        pivot = next((index for index, value in enumerate(row) if value), None)
        if pivot is None:
            continue
        inverse = pow(row[pivot], -1, prime)
        row = [(value * inverse) % prime for value in row]
        pivots.append(pivot)
        echelon.append(row)
        selected.append(source_index)
    return tuple(echelon), tuple(pivots), tuple(selected)


target_vectors = []
target_labels = []
target_degrees = []
for degree, label, normal in target_records_mod:
    if not normal:
        continue
    target_degrees.append(degree)
    target_labels.append(label)
    target_vectors.append(derivative_gap_vector_mod(normal))

r_prime_vector = derivative_gap_vector_mod(r_mod.derivative())
require(any(r_prime_vector), "seminormal derivative class is nonzero in E_x")
target_echelon, target_pivots, target_selected = independent_integer_rows_mod(
    target_vectors, DERIVATIVE_PRIME
)
target_rank = len(target_echelon)
target_augmented_rank = len(
    independent_integer_rows_mod(target_vectors + [r_prime_vector], DERIVATIVE_PRIME)[0]
)

rank_ladder = []
for cutoff in (30, 60, 90, 120, 150, 180, 210, 240):
    prefix = [
        vector for degree, vector in zip(target_degrees, target_vectors) if degree <= cutoff
    ]
    prefix_rank = len(independent_integer_rows_mod(prefix, DERIVATIVE_PRIME)[0])
    prefix_augmented = len(
        independent_integer_rows_mod(prefix + [r_prime_vector], DERIVATIVE_PRIME)[0]
    )
    rank_ladder.append((cutoff, len(prefix), prefix_rank, prefix_augmented))

print("scope=representative_level_fixed_x_first_normal_sidecar;not_full_THM4067")
print("rings=A=K[b,c,e]/(c^2e-b(b+4));S=gamma0(A);N=K[x];N/S_dimension=89")
print("map=gamma1:A->N;gamma1(H)=d/ds H(B(x,Q+s),C(x,Q+s),E(x,Q+s))|s=0")
print("normal_generator_degrees=" + ",".join(str(item.degree()) for item in normal_generators))
print("product_ceiling=" + str(product_ceiling))
print("retained_values=" + ";".join(
    ",".join(serialize_field(value) for value in row) for row in retained_values
))
print("retained_Lambda_responses=" + ",".join(serialize_field(value) for value in retained_responses))
print("r_derivative_Lambda_nonzero=True")
print("r_derivative_Lambda_hash=" + hashlib.sha256(
    serialize_field(r_derivative_response).encode("ascii")
).hexdigest())
print(
    "r_derivative_Lambda_mod421="
    + str(field_mod_at(r_derivative_response, DERIVATIVE_PRIME, DERIVATIVE_ALPHA))
)
print("individual_S_module_ranks=" + ",".join(map(str, individual_ranks)))
print("individual_augmented_ranks_with_r=" + ",".join(map(str, individual_augmented_ranks)))
print(f"NS_combined_modular_rank={rank};augmented_rank_with_r={augmented_rank};r_hit={augmented_rank == rank}")
print("selected_column_count=" + str(len(selected)))
print("selected_column_hash=" + hashlib.sha256(
    ";".join(column_labels[index] for index in selected).encode("ascii")
).hexdigest())
print("column_vector_hash=" + vector_hash(columns))
print("r_vector_hash=" + vector_hash((r_vector,)))
print("special_fibre_density_quotient=E_x=N/dS;dimension=89;basis_degrees=" + ",".join(map(str, derivative_gaps)))
print(f"derivative_good_reduction=p{DERIVATIVE_PRIME},alpha{DERIVATIVE_ALPHA}")
print("target_monomial_cutoff=" + str(TARGET_CUTOFF) + ";target_column_count=" + str(len(target_vectors)))
print("target_rank_ladder=" + ";".join(
    f"{cutoff}:{count}:{local_rank}:{local_augmented}"
    for cutoff, count, local_rank, local_augmented in rank_ladder
))
print(
    f"E_x_target_modular_rank={target_rank};augmented_rank_with_rprime={target_augmented_rank};"
    f"rprime_hit={target_augmented_rank == target_rank}"
)
print("E_x_rank_logic=positive_full_minor_implies_characteristic_zero_surjectivity")
print("E_x_selected_target_count=" + str(len(target_selected)))
print("E_x_selected_target_hash=" + hashlib.sha256(
    ";".join(target_labels[index] for index in target_selected).encode("ascii")
).hexdigest())
print("E_x_target_vector_hash=" + integer_vector_hash(target_vectors))
print("rprime_vector_hash=" + integer_vector_hash((r_prime_vector,)))
print("RESULT=PASS")
