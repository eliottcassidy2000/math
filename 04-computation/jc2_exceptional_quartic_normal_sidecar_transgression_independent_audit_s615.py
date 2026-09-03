"""Independent direct-matrix audit of the exceptional normal-sidecar scout.

This imports only the audited THM-3703 restriction-algebra certificate.  It
does not import the scout or the THM-4034 conductor.  Instead of reducing
sidecars to the shifted gap normal form, it places dS and the sidecar columns
in one ordinary coefficient matrix over a good finite-field specialization.
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
PRIME = 421
ALPHA_MOD = 126
TARGET_CUTOFF = 180


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
        runner = lambda: exec(compile(source, path.as_posix(), "exec"), module.__dict__)
    with contextlib.redirect_stdout(io.StringIO()):
        runner()
    return module


sagbi = load_sagbi(SAGBI_PATH, "quartic_sagbi_transgression_audit_s615")
field = sagbi.field
X = sagbi.X
zero = field.zero
one = field.one
B = sagbi.B
C = sagbi.C
E = sagbi.E
D = sagbi.D
Z = sagbi.Z
gaps = tuple(sagbi.EXPECTED_GAPS)
require(len(gaps) == 89 and gaps[-1] == 169, "THM-3703 gaps")
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
    "good quartic specialization",
)


def coefficient(polynomial, degree):
    return polynomial.get((degree,), zero)


def rational_mod(value):
    numerator, denominator = map(int, sp.Rational(value).as_numer_denom())
    require(denominator % PRIME != 0, "good denominator")
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


def poly_mul(left, right):
    if not left or not right:
        return ()
    answer = [0] * (len(left) + len(right) - 1)
    for left_degree, left_value in enumerate(left):
        for right_degree, right_value in enumerate(right):
            answer[left_degree + right_degree] += left_value * right_value
    return trim(answer)


def poly_pow(polynomial, exponent):
    answer = (1,)
    factor = polynomial
    remaining = exponent
    while remaining:
        if remaining & 1:
            answer = poly_mul(answer, factor)
        factor = poly_mul(factor, factor)
        remaining //= 2
    return answer


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


B_mod, C_mod, E_mod, D_mod, Z_mod = map(
    convert_polynomial, (B, C, E, D, Z)
)
require(tuple(map(degree, (B_mod, C_mod, E_mod, D_mod, Z_mod))) == (30, 21, 18, 10, 18), "specialized degrees")

# Differentiate the compiler at fixed x along q=Q+s.  These formulas are
# derived directly before specialization, rather than copied from the scout.
normal_exact = (
    3 * X**2 * D * (D + 2),
    2 * X**3 * (D + 1),
    2 * (D + 1),
)
require(
    2 * C * E * normal_exact[1] + C**2 * normal_exact[2]
    == (2 * B + 4) * normal_exact[0],
    "surface relation differentiated exactly",
)
normal_mod = tuple(map(convert_polynomial, normal_exact))

# If gamma1 factored through S, it would induce a derivation S -> K[x].
# Since Frac(S)=K(x) in characteristic zero, that derivation would be
# v(x)d/dx for one v in K(x), forcing gamma1(B)C'=gamma1(C)B'.  The exact
# nonzero cross-product is therefore a clean non-descent certificate.
tangency_defect_exact = normal_exact[0] * C.diff(X) - normal_exact[1] * B.diff(X)
require(bool(tangency_defect_exact), "representative sidecar does not descend to S")
tangency_defect_mod = convert_polynomial(tangency_defect_exact)
require(bool(tangency_defect_mod), "non-descent survives good reduction")

# Independent exact retained-fibre check of the advertised local response.
points = (-1, 0, 1)
normal_C_values = tuple(
    normal_exact[1].evaluate(X, field.convert(point)) for point in points
)
require(
    normal_C_values == (field.convert(2), zero, field.convert(-2)),
    "gamma1(C) retained values",
)
lambda_C = (
    field.convert(sp.Rational(5, 18)) * normal_C_values[0]
    - normal_C_values[1]
    + field.convert(sp.Rational(13, 18)) * normal_C_values[2]
)
require(lambda_C == field.convert(sp.Rational(-8, 9)), "local Lambda response")
require(
    tuple(poly_evaluate(normal_mod[1], point % PRIME) for point in points)
    == (2, 0, PRIME - 2),
    "modular retained response",
)


def powers_through(polynomial, count):
    powers = [(1,)]
    for _ in range(count):
        powers.append(poly_mul(powers[-1], polynomial))
    return tuple(powers)


B_powers = powers_through(B_mod, TARGET_CUTOFF // 30)
C_powers = powers_through(C_mod, TARGET_CUTOFF // 21)
E_powers = powers_through(E_mod, TARGET_CUTOFF // 18)
records = []
for i, B_power in enumerate(B_powers):
    for j, C_power in enumerate(C_powers):
        for k, E_power in enumerate(E_powers):
            filtration_degree = 30 * i + 21 * j + 18 * k
            if filtration_degree > TARGET_CUTOFF:
                continue
            normal = ()
            if i:
                normal = poly_add(
                    normal,
                    poly_scale(
                        i,
                        poly_mul(
                            poly_mul(B_powers[i - 1], C_power),
                            poly_mul(E_power, normal_mod[0]),
                        ),
                    ),
                )
            if j:
                normal = poly_add(
                    normal,
                    poly_scale(
                        j,
                        poly_mul(
                            poly_mul(B_power, C_powers[j - 1]),
                            poly_mul(E_power, normal_mod[1]),
                        ),
                    ),
                )
            if k:
                normal = poly_add(
                    normal,
                    poly_scale(
                        k,
                        poly_mul(
                            poly_mul(B_power, C_power),
                            poly_mul(E_powers[k - 1], normal_mod[2]),
                        ),
                    ),
                )
            if normal:
                records.append((filtration_degree, f"B^{i}C^{j}E^{k}", normal))
records.sort(key=lambda item: (item[0], item[1]))
ambient_degree = max(degree(polynomial) for _, _, polynomial in records)
require(PRIME > ambient_degree + 1, "derivative leading coefficients remain units")
ambient_dimension = ambient_degree + 1

# Directly span dS inside the full coefficient space.  THM-3703 gives one
# monic canonical polynomial at each semigroup degree, so derivatives of the
# canonical elements through degree ambient_degree+1 give exactly the
# truncated dS, with no quotient normal-form computation.
canonical_mod = {}
for exact_basis in sagbi.module_basis:
    polynomial = convert_polynomial(exact_basis)
    while degree(polynomial) <= ambient_degree + 1:
        require(degree(polynomial) not in canonical_mod, "unique canonical degree")
        canonical_mod[degree(polynomial)] = polynomial
        polynomial = poly_mul(polynomial, Z_mod)
require(
    all(item in canonical_mod for item in range(170, ambient_degree + 2)),
    "semigroup tail coverage",
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
    require(degree(derivative) == canonical_degree - 1, "canonical derivative degree")
    require(add_to_echelon(coefficient_vector(derivative), echelon, pivots), "independent dS derivative")
dS_rank = len(echelon)
require(ambient_dimension - dS_rank == len(gaps), "direct quotient dimension")

rank_ladder = []
selected_labels = []
for cutoff in (30, 60, 90, 120, 150, 180):
    if cutoff == 30:
        active_echelon = [row[:] for row in echelon]
        active_pivots = list(pivots)
        processed = 0
    while processed < len(records) and records[processed][0] <= cutoff:
        _, label, polynomial = records[processed]
        if add_to_echelon(
            coefficient_vector(polynomial), active_echelon, active_pivots
        ):
            selected_labels.append(label)
        processed += 1
    rank_ladder.append(
        (cutoff, processed, len(active_echelon), len(active_echelon) - dS_rank)
    )

combined_rank = len(active_echelon)
sidecar_quotient_rank = combined_rank - dS_rank
require(combined_rank == ambient_dimension, "full direct ambient rank")
require(sidecar_quotient_rank == 89, "full normal-sidecar quotient rank")
require(len(selected_labels) == 89, "89 selected sidecars")

print("scope=independent_direct_matrix_audit;representative_level_fixed_x;not_THM4067")
print("imports=audited_THM3703_only;primary_scout_import=no;THM4034_import=no")
print("field_specialization=p421,alpha126;quartic_root=yes")
print("gamma1_generator_degrees=" + ",".join(map(str, map(degree, normal_mod))))
print("gamma1_C_retained_values=(2,0,-2);Lambda_gamma1_C=-8/9")
print(
    "nondescent_tangency_defect_degree="
    + str(degree(tangency_defect_mod))
    + ";mod421_hash="
    + hashlib.sha256(
        ",".join(map(str, tangency_defect_mod)).encode("ascii")
    ).hexdigest()
)
print("target_cutoff=180;target_columns=" + str(len(records)))
print(
    f"ambient_degree={ambient_degree};ambient_dimension={ambient_dimension};"
    f"dS_rank={dS_rank};E_x_dimension={ambient_dimension-dS_rank}"
)
print(
    "rank_ladder=cutoff:columns:combined_rank:quotient_rank="
    + ";".join(":".join(map(str, row)) for row in rank_ladder)
)
print(
    f"combined_rank={combined_rank};sidecar_quotient_rank={sidecar_quotient_rank};"
    "E_x_surjective=yes"
)
print(
    "selected_sidecar_hash="
    + hashlib.sha256(";".join(selected_labels).encode("ascii")).hexdigest()
)
print("logic=good_reduction_full_minor_implies_characteristic_zero_surjectivity")
print("typing=gamma1_on_target_representatives;descent_to_S_not_assumed_or_used")
print("RESULT=PASS")
