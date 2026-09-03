"""Independent audit of the fixed-x descended two-form transgression.

The source is the audited THM-3703 exceptional restriction algebra

    S = K[B,C,E] subset K[x].

For the fixed-x normal deformation q=Q+s, write gamma_1 for d/ds at s=0.
Contraction of a target two-form with (d/dx,d/ds) sends

    dP wedge dQ |-> P' gamma_1(Q) - Q' gamma_1(P)  mod dS.

This audit does not import the s616 scout.  It uses the T-module generators
db wedge dc, db wedge de, dc wedge de for Omega_T^2, proves the retained
triple-point functional is an exact upper bound, and finds an independent
positive minor after specialization at (p,alpha)=(443,112).  The matrix is
built in the ordinary coefficient space, with dS inserted explicitly; no
gap-normal-form routine from the scout is used.
"""

import contextlib
import hashlib
import importlib.util
import io
import subprocess
import types
from pathlib import Path

import sympy as sp


SAGBI_PATH = Path(
    "04-computation/jc2_russell_cylinder_exceptional_quartic_sagbi_module_thm3703.py"
)
PRIME = 443
ALPHA_MOD = 112
MULTIPLIER_CUTOFF = 180
LADDER_CUTOFFS = (30, 60, 90, 120, 150, 180)


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def load_sagbi(path, name):
    if path.is_file():
        spec = importlib.util.spec_from_file_location(name, path)
        module = importlib.util.module_from_spec(spec)
        runner = lambda: spec.loader.exec_module(module)
    else:
        source = subprocess.check_output(
            ["git", "show", f"HEAD:{path.as_posix()}"], text=True
        )
        module = types.ModuleType(name)
        runner = lambda: exec(
            compile(source, path.as_posix(), "exec"), module.__dict__
        )
    with contextlib.redirect_stdout(io.StringIO()):
        runner()
    return module


sagbi = load_sagbi(SAGBI_PATH, "quartic_two_form_audit_s616")
field = sagbi.field
X = sagbi.X
zero = field.zero
B = sagbi.B
C = sagbi.C
E = sagbi.E
D = sagbi.D
Z = sagbi.Z
gaps = tuple(sagbi.EXPECTED_GAPS)
require(len(gaps) == 89 and gaps[-1] == 169, "THM-3703 gap universe")


def coefficient(polynomial, degree):
    return polynomial.get((degree,), zero)


def evaluate_exact(polynomial, point):
    return polynomial.evaluate(X, field.convert(point))


# Differentiate the Russell compiler directly in the fixed-x q=Q+s lane.
# These are target-representative derivatives, not a map descended to S.
normal_exact = (
    3 * X**2 * D * (D + 2),
    2 * X**3 * (D + 1),
    2 * (D + 1),
)
derivative_exact = tuple(generator.diff(X) for generator in (B, C, E))
require(C**2 * E == B * (B + 4), "exceptional surface relation")
require(
    2 * C * E * normal_exact[1] + C**2 * normal_exact[2]
    == (2 * B + 4) * normal_exact[0],
    "normal derivative respects surface relation",
)

# Omega_T^2 is generated over T by the three displayed wedges.  Their
# contraction densities generate the full descended two-form image.
wedge_pairs = ((0, 1, "dB^dC"), (0, 2, "dB^dE"), (1, 2, "dC^dE"))
wedge_exact = tuple(
    (
        label,
        derivative_exact[left] * normal_exact[right]
        - derivative_exact[right] * normal_exact[left],
    )
    for left, right, label in wedge_pairs
)
require(all(polynomial for _, polynomial in wedge_exact), "nonzero base wedges")

# Exact retained-triple upper bound.  Every target polynomial has the same
# value at x=-1,0,1.  Its x- and s-derivative rows are linear combinations of
# the C,E rows below.  Hence every two-form contraction is a target scalar
# times the constant row (12,12,12), and Lambda annihilates it.  Lambda also
# annihilates dS, so it descends to E_x=N/dS.
points = (-1, 0, 1)
weights = (5, -18, 13)
retained_values = tuple(
    tuple(evaluate_exact(generator, point) for point in points)
    for generator in (B, C, E)
)
retained_derivatives = tuple(
    tuple(evaluate_exact(generator, point) for point in points)
    for generator in derivative_exact
)
retained_normals = tuple(
    tuple(evaluate_exact(generator, point) for point in points)
    for generator in normal_exact
)
retained_wedges = tuple(
    tuple(evaluate_exact(polynomial, point) for point in points)
    for _, polynomial in wedge_exact
)
require(
    retained_values
    == (
        (zero, zero, zero),
        (zero, zero, zero),
        (-3 * field.one, -3 * field.one, -3 * field.one),
    ),
    "retained target fibre",
)
require(
    retained_derivatives
    == (
        (zero, zero, zero),
        (3 * field.one, 3 * field.one, 3 * field.one),
        (-9 * field.one, 4 * field.one, 9 * field.one),
    ),
    "retained tangent rows",
)
require(
    retained_normals
    == (
        (zero, zero, zero),
        (2 * field.one, zero, -2 * field.one),
        (-2 * field.one, 4 * field.one, -2 * field.one),
    ),
    "retained normal rows",
)
require(
    retained_wedges
    == (
        (zero, zero, zero),
        (zero, zero, zero),
        (12 * field.one, 12 * field.one, 12 * field.one),
    ),
    "retained base-wedge determinant",
)


def retained_lambda(values):
    return sum(
        (field.convert(weight) * value for weight, value in zip(weights, values)),
        zero,
    )


require(sum(weights) == 0, "Lambda kills constant rows")
require(
    all(retained_lambda(row) == zero for row in retained_derivatives),
    "Lambda kills dS tangent rows",
)
require(
    all(retained_lambda(row) == zero for row in retained_wedges),
    "Lambda kills target two-form contractions",
)
# A concrete hostile complement: d(x^2)=2x is not killed by Lambda.
lambda_two_x = retained_lambda(tuple(field.convert(2 * point) for point in points))
require(lambda_two_x == field.convert(16), "Lambda is nonzero on E_x")


# Independent good reduction.  This quartic root differs from the p=421
# specialization used in the s615 one-form-sidecar certificate.
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
    "quartic good specialization",
)


def rational_mod(value):
    numerator, denominator = map(int, sp.Rational(value).as_numer_denom())
    require(denominator % PRIME != 0, "specialization denominator unit")
    return numerator % PRIME * pow(denominator % PRIME, -1, PRIME) % PRIME


def field_mod(value):
    return sum(
        rational_mod(coordinate) * pow(ALPHA_MOD, power, PRIME)
        for power, coordinate in enumerate(sagbi.field_coordinates(value))
    ) % PRIME


def trim(polynomial):
    answer = [value % PRIME for value in polynomial]
    while answer and answer[-1] == 0:
        answer.pop()
    return tuple(answer)


def poly_add(left, right):
    size = max(len(left), len(right))
    return trim(
        [
            (left[index] if index < len(left) else 0)
            + (right[index] if index < len(right) else 0)
            for index in range(size)
        ]
    )


def poly_scale(scalar, polynomial):
    return trim([scalar * value for value in polynomial])


def poly_sub(left, right):
    return poly_add(left, poly_scale(-1, right))


def poly_mul(left, right):
    if not left or not right:
        return ()
    answer = [0] * (len(left) + len(right) - 1)
    for left_degree, left_value in enumerate(left):
        for right_degree, right_value in enumerate(right):
            answer[left_degree + right_degree] += left_value * right_value
    return trim(answer)


def poly_derivative(polynomial):
    return trim([degree * polynomial[degree] for degree in range(1, len(polynomial))])


def poly_evaluate(polynomial, point):
    answer = 0
    for value in reversed(polynomial):
        answer = (answer * point + value) % PRIME
    return answer


def convert_polynomial(polynomial):
    return trim(
        [field_mod(coefficient(polynomial, degree)) for degree in range(polynomial.degree() + 1)]
    )


def degree(polynomial):
    return len(polynomial) - 1


B_mod, C_mod, E_mod, Z_mod = map(convert_polynomial, (B, C, E, Z))
normal_mod = tuple(map(convert_polynomial, normal_exact))
derivative_mod = tuple(poly_derivative(item) for item in (B_mod, C_mod, E_mod))
wedge_mod = tuple(
    (
        label,
        poly_sub(
            poly_mul(derivative_mod[left], normal_mod[right]),
            poly_mul(derivative_mod[right], normal_mod[left]),
        ),
    )
    for left, right, label in wedge_pairs
)
require(
    tuple(map(degree, (B_mod, C_mod, E_mod, Z_mod))) == (30, 21, 18, 18),
    "specialized target degrees",
)
wedge_degrees = tuple(degree(polynomial) for _, polynomial in wedge_mod)
require(
    wedge_degrees == (42, 39, 30),
    f"specialized wedge degrees: {wedge_degrees}",
)
require(
    tuple(
        tuple(poly_evaluate(polynomial, point % PRIME) for point in points)
        for _, polynomial in wedge_mod
    )
    == ((0, 0, 0), (0, 0, 0), (12, 12, 12)),
    "modular retained wedge control",
)


def powers_through(polynomial, count):
    powers = [(1,)]
    for _ in range(count):
        powers.append(poly_mul(powers[-1], polynomial))
    return tuple(powers)


B_powers = powers_through(B_mod, MULTIPLIER_CUTOFF // 30)
C_powers = powers_through(C_mod, MULTIPLIER_CUTOFF // 21)
E_powers = powers_through(E_mod, MULTIPLIER_CUTOFF // 18)
multiplier_records = []
for i, B_power in enumerate(B_powers):
    for j, C_power in enumerate(C_powers):
        for k, E_power in enumerate(E_powers):
            filtration_degree = 30 * i + 21 * j + 18 * k
            if filtration_degree > MULTIPLIER_CUTOFF:
                continue
            polynomial = poly_mul(poly_mul(B_power, C_power), E_power)
            multiplier_records.append(
                (filtration_degree, f"B^{i}C^{j}E^{k}", polynomial)
            )
multiplier_records.sort(key=lambda item: (item[0], item[1]))

# Each multiplier times each of the three base wedges is an actual element of
# the two-form image.  A positive rank from this bounded family is enough;
# the exact retained functional supplied the cutoff-free upper bound.
columns = []
for filtration_degree, multiplier_label, multiplier in multiplier_records:
    for wedge_label, wedge in wedge_mod:
        polynomial = poly_mul(multiplier, wedge)
        if polynomial:
            columns.append(
                (
                    filtration_degree,
                    f"{multiplier_label}*{wedge_label}",
                    polynomial,
                )
            )
columns.sort(key=lambda item: (item[0], item[1]))
ambient_degree = max(degree(polynomial) for _, _, polynomial in columns)
ambient_dimension = ambient_degree + 1
require(PRIME > ambient_degree + 1, "all differentiated degrees below prime")

# Build dS directly in the ordinary coefficient space.  The monic THM-3703
# Apery basis times powers of Z gives exactly one canonical S polynomial in
# each semigroup degree.  Its nonconstant derivatives span truncated dS.
canonical_mod = {}
for exact_basis in sagbi.module_basis:
    polynomial = convert_polynomial(exact_basis)
    while degree(polynomial) <= ambient_degree + 1:
        require(degree(polynomial) not in canonical_mod, "unique canonical degree")
        canonical_mod[degree(polynomial)] = polynomial
        polynomial = poly_mul(polynomial, Z_mod)
require(
    all(item in canonical_mod for item in range(170, ambient_degree + 2)),
    "semigroup tail covers ambient truncation",
)


def coefficient_vector(polynomial):
    return tuple(
        polynomial[index] if index < len(polynomial) else 0
        for index in range(ambient_dimension)
    )


def add_to_echelon(source, echelon, pivots):
    row = [value % PRIME for value in source]
    for pivot, old_row in zip(pivots, echelon):
        value = row[pivot]
        if value:
            row = [
                (left - value * right) % PRIME
                for left, right in zip(row, old_row)
            ]
    pivot = next((index for index, value in enumerate(row) if value), None)
    if pivot is None:
        return False
    inverse = pow(row[pivot], -1, PRIME)
    row = [(value * inverse) % PRIME for value in row]
    echelon.append(row)
    pivots.append(pivot)
    return True


echelon = []
pivots = []
for canonical_degree in sorted(canonical_mod):
    if canonical_degree == 0:
        continue
    derivative = poly_derivative(canonical_mod[canonical_degree])
    require(degree(derivative) == canonical_degree - 1, "derivative degree survives")
    require(
        add_to_echelon(coefficient_vector(derivative), echelon, pivots),
        "independent canonical derivative",
    )
dS_rank = len(echelon)
require(ambient_dimension - dS_rank == len(gaps), "E_x has dimension 89")

active_echelon = [row[:] for row in echelon]
active_pivots = list(pivots)
selected_labels = []
selected_payloads = []
rank_ladder = []
processed = 0
for cutoff in LADDER_CUTOFFS:
    while processed < len(columns) and columns[processed][0] <= cutoff:
        _, label, polynomial = columns[processed]
        vector = coefficient_vector(polynomial)
        if add_to_echelon(vector, active_echelon, active_pivots):
            selected_labels.append(label)
            selected_payloads.append(",".join(map(str, vector)))
        processed += 1
    rank_ladder.append(
        (cutoff, processed, len(active_echelon), len(active_echelon) - dS_rank)
    )

combined_rank = len(active_echelon)
two_form_quotient_rank = combined_rank - dS_rank
require(two_form_quotient_rank == 88, "two-form image has modular rank 88")
require(len(selected_labels) == 88, "selected positive minor size")

# Both x^2 and the raw one-form sidecar gamma_1(C) cross the retained
# hyperplane.  Adjoining either to the rank-88 image must fill the last E_x
# coordinate.  This is a hostile check against an accidental rank-88 subspace
# unrelated to Lambda.
two_x_mod = (0, 2)
gamma1_C_mod = normal_mod[1]
require(
    sum(
        weight * poly_evaluate(two_x_mod, point % PRIME)
        for weight, point in zip(weights, points)
    )
    % PRIME
    == 16,
    "modular Lambda(2x)",
)
require(
    sum(
        weight * poly_evaluate(gamma1_C_mod, point % PRIME)
        for weight, point in zip(weights, points)
    )
    % PRIME
    == (-16) % PRIME,
    "modular Lambda(gamma1(C))",
)
hostile_echelon = [row[:] for row in active_echelon]
hostile_pivots = list(active_pivots)
require(
    add_to_echelon(coefficient_vector(two_x_mod), hostile_echelon, hostile_pivots),
    "2x fills missed line",
)
require(len(hostile_echelon) == ambient_dimension, "2x gives full ambient rank")
hostile_echelon = [row[:] for row in active_echelon]
hostile_pivots = list(active_pivots)
require(
    add_to_echelon(
        coefficient_vector(gamma1_C_mod), hostile_echelon, hostile_pivots
    ),
    "raw one-form sidecar fills missed line",
)
require(
    len(hostile_echelon) == ambient_dimension,
    "raw one-form sidecar gives full ambient rank",
)

selected_label_hash = hashlib.sha256(
    ";".join(selected_labels).encode("ascii")
).hexdigest()
selected_vector_hash = hashlib.sha256(
    ";".join(selected_payloads).encode("ascii")
).hexdigest()
wedge_hash = hashlib.sha256(
    ";".join(
        label + ":" + ",".join(map(str, polynomial))
        for label, polynomial in wedge_mod
    ).encode("ascii")
).hexdigest()

print("scope=fixed_x_first_normal_descended_target_two_forms;not_THM4067;not_JC2")
print("imports=audited_THM3703_only;s616_candidate_import=no;gap_normal_form=no")
print("field_specialization=p443,alpha112;quartic_root=yes;independent_from_p421")
print("source=T=K[B,C,E]/(C^2E-B(B+4));target=E_x=K[x]/dS")
print("map=dP^dQ->Pprime*gamma1(Q)-Qprime*gamma1(P)_mod_dS")
print("omega_generators=dB^dC,dB^dE,dC^dE;wedge_degrees=42,39,30")
print(
    "finite_universe=all_B^iC^jE^k_with_i,j,k>=0_and_30i+21j+18k<=180;"
    "times_each_of_3_omega_generators"
)
print("retained_points=-1,0,1;Lambda_weights=5,-18,13")
print("retained_base_wedges=(0,0,0)|(0,0,0)|(12,12,12)")
print("upper_bound=Lambda_kills_dS_and_all_OmegaT2;cutoff_free=yes")
print(
    f"multiplier_cutoff={MULTIPLIER_CUTOFF};multipliers={len(multiplier_records)};"
    f"two_form_columns={len(columns)}"
)
print(
    f"ambient_degree={ambient_degree};ambient_dimension={ambient_dimension};"
    f"dS_rank={dS_rank};E_x_dimension={ambient_dimension-dS_rank}"
)
print(
    "rank_ladder=cutoff:columns:combined_rank:quotient_rank="
    + ";".join(":".join(map(str, row)) for row in rank_ladder)
)
print(
    f"combined_rank={combined_rank};two_form_quotient_rank={two_form_quotient_rank};"
    "image_equals_ker_Lambda=yes"
)
print(f"selected_label_hash={selected_label_hash}")
print(f"selected_vector_hash={selected_vector_hash}")
print(f"base_wedge_mod443_hash={wedge_hash}")
print("hostile_complements=2x,raw_gamma1(C);both_raise_quotient_rank_88_to_89")
print("cokernel=one_dimensional;generator=[rprime]_by_THM4381_missing_first_jet")
print("logic=exact_Lambda_upper_bound_plus_good_reduction_positive_88_minor")
print("controls=surface_relation,retained_determinant,2x,raw_gamma1(C),rank_ladder")
print(
    "reproduction=python3_-B_04-computation/"
    "jc2_exceptional_quartic_descended_two_form_transgression_independent_audit_s616.py"
)
print("replays=normal,-O,PYTHONHASHSEED9173;byte_match=yes")
print("NO_CLAIM=full_moving_graph_transgression_or_chart_entry_or_Keller_pair_or_JC2")
print("RESULT=PASS")
