"""Finite-exact scout for descended two-form transgression at the quartic triple.

This is a bridge computation, not a theorem companion.  It starts from the
independently audited fixed-x sidecar data of the s615 session.  For an actual
target two-form dP wedge dQ it forms the contracted density

    tau(P,Q) = P' gamma1(Q) - Q' gamma1(P)

and computes its image in E_x=K[x]/dS.  A good-reduction rank is used only in
the sound positive-minor direction.  The characteristic-zero upper bound is
the exact retained-triple identity proved below, not modular rank failure.
"""

import contextlib
import hashlib
import importlib.util
import io
import subprocess
import types
from pathlib import Path


AUDIT_PATH = Path(
    "04-computation/"
    "jc2_exceptional_quartic_normal_sidecar_transgression_independent_audit_s615.py"
)

CHECK_COUNT = 0


def require(condition, label):
    global CHECK_COUNT
    CHECK_COUNT += 1
    if not condition:
        raise RuntimeError(label)


def load_audit(path, name):
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


audit = load_audit(AUDIT_PATH, "quartic_two_form_transgression_s616")
PRIME = audit.PRIME
TARGET_CUTOFF = audit.TARGET_CUTOFF
require((PRIME, audit.ALPHA_MOD, TARGET_CUTOFF) == (421, 126, 180), "frozen good fibre")


# At the retained fibre B=0,C=0,E=-3, the Russell surface is smooth because
# d(C^2 E-B(B+4))/dB=-4.  Thus C,E are regular local coordinates.  The
# ordinary and fixed-x normal tangent rows are exact field-valued rows.
points = (-1, 0, 1)
field = audit.field
zero = audit.zero
one = audit.one


def evaluate_exact(polynomial, point):
    return polynomial.evaluate(audit.X, field.convert(point))


C_tangent = tuple(evaluate_exact(audit.C.diff(audit.X), point) for point in points)
E_tangent = tuple(evaluate_exact(audit.E.diff(audit.X), point) for point in points)
C_normal = tuple(evaluate_exact(audit.normal_exact[1], point) for point in points)
E_normal = tuple(evaluate_exact(audit.normal_exact[2], point) for point in points)
B_tangent = tuple(evaluate_exact(audit.B.diff(audit.X), point) for point in points)
B_normal = tuple(evaluate_exact(audit.normal_exact[0], point) for point in points)

require(B_tangent == (zero, zero, zero), "B ordinary tangent vanishes")
require(B_normal == (zero, zero, zero), "B normal tangent vanishes")
require(C_tangent == tuple(field.convert(item) for item in (3, 3, 3)), "C tangent row")
require(E_tangent == tuple(field.convert(item) for item in (-9, 4, 9)), "E tangent row")
require(C_normal == tuple(field.convert(item) for item in (2, 0, -2)), "C normal row")
require(E_normal == tuple(field.convert(item) for item in (-2, 4, -2)), "E normal row")

base_wedge = tuple(
    C_tangent[index] * E_normal[index]
    - E_tangent[index] * C_normal[index]
    for index in range(3)
)
require(base_wedge == tuple(field.convert(12) for _ in range(3)), "constant base wedge")

weights = (
    field.convert(5) / field.convert(18),
    -one,
    field.convert(13) / field.convert(18),
)
require(sum(weights, zero) == zero, "Lambda kills constants")
require(sum((w * v for w, v in zip(weights, C_tangent)), zero) == zero, "Lambda dC")
require(sum((w * v for w, v in zip(weights, E_tangent)), zero) == zero, "Lambda dE")


# Rebuild every target monomial together with its restriction and its normal
# sidecar.  This is deliberately independent of the s615 scout's gap normal
# form; the imported audit uses a direct ordinary coefficient matrix.
records = []
for i, B_power in enumerate(audit.B_powers):
    for j, C_power in enumerate(audit.C_powers):
        for k, E_power in enumerate(audit.E_powers):
            filtration_degree = 30 * i + 21 * j + 18 * k
            if filtration_degree > TARGET_CUTOFF:
                continue
            restriction = audit.poly_mul(audit.poly_mul(B_power, C_power), E_power)
            normal = ()
            if i:
                normal = audit.poly_add(
                    normal,
                    audit.poly_scale(
                        i,
                        audit.poly_mul(
                            audit.poly_mul(audit.B_powers[i - 1], C_power),
                            audit.poly_mul(E_power, audit.normal_mod[0]),
                        ),
                    ),
                )
            if j:
                normal = audit.poly_add(
                    normal,
                    audit.poly_scale(
                        j,
                        audit.poly_mul(
                            audit.poly_mul(B_power, audit.C_powers[j - 1]),
                            audit.poly_mul(E_power, audit.normal_mod[1]),
                        ),
                    ),
                )
            if k:
                normal = audit.poly_add(
                    normal,
                    audit.poly_scale(
                        k,
                        audit.poly_mul(
                            audit.poly_mul(B_power, C_power),
                            audit.poly_mul(audit.E_powers[k - 1], audit.normal_mod[2]),
                        ),
                    ),
                )
            records.append(
                (filtration_degree, f"B^{i}C^{j}E^{k}", restriction, normal)
            )
records.sort(key=lambda item: (item[0], item[1]))

carriers = (
    ("B", audit.B_mod, audit.normal_mod[0]),
    ("C", audit.C_mod, audit.normal_mod[1]),
    ("E", audit.E_mod, audit.normal_mod[2]),
)
columns = []
for carrier_label, carrier, carrier_normal in carriers:
    carrier_derivative = audit.poly_derivative(carrier)
    for filtration_degree, target_label, restriction, normal in records:
        tau = audit.poly_add(
            audit.poly_mul(carrier_derivative, normal),
            audit.poly_scale(
                -1,
                audit.poly_mul(audit.poly_derivative(restriction), carrier_normal),
            ),
        )
        if tau:
            columns.append(
                (filtration_degree, carrier_label, target_label, tau)
            )

ambient_degree = max(audit.degree(item[3]) for item in columns)
require(PRIME > ambient_degree + 1, "good prime exceeds every density degree")
ambient_dimension = ambient_degree + 1


def coefficient_vector(polynomial):
    return tuple(
        polynomial[index] if index < len(polynomial) else 0
        for index in range(ambient_dimension)
    )


# Span dS in the common ambient coefficient space.  THM-3703's monic Apery
# basis has one canonical element in every semigroup degree, so its
# derivatives are independent and give the exact truncated dS.
canonical_mod = {}
for exact_basis in audit.sagbi.module_basis:
    polynomial = audit.convert_polynomial(exact_basis)
    while audit.degree(polynomial) <= ambient_degree + 1:
        require(audit.degree(polynomial) not in canonical_mod, "unique canonical degree")
        canonical_mod[audit.degree(polynomial)] = polynomial
        polynomial = audit.poly_mul(polynomial, audit.Z_mod)
require(
    all(degree in canonical_mod for degree in range(170, ambient_degree + 2)),
    "semigroup tail coverage",
)

base_echelon = []
base_pivots = []
for canonical_degree in sorted(canonical_mod):
    if canonical_degree == 0:
        continue
    derivative = audit.poly_derivative(canonical_mod[canonical_degree])
    require(audit.degree(derivative) == canonical_degree - 1, "derivative degree survives")
    require(
        audit.add_to_echelon(
            coefficient_vector(derivative), base_echelon, base_pivots
        ),
        "independent dS derivative",
    )
dS_rank = len(base_echelon)
require(ambient_dimension - dS_rank == 89, "E_x has THM-3703 dimension 89")


def lambda_mod(polynomial):
    inverse_18 = pow(18, -1, PRIME)
    modular_weights = (5 * inverse_18 % PRIME, PRIME - 1, 13 * inverse_18 % PRIME)
    return sum(
        weight * audit.poly_evaluate(polynomial, point % PRIME)
        for weight, point in zip(modular_weights, points)
    ) % PRIME


require(all(lambda_mod(item[3]) == 0 for item in columns), "all sampled wedge responses lie in ker Lambda")

individual_ranks = []
for carrier_label, _, _ in carriers:
    echelon = [row[:] for row in base_echelon]
    pivots = list(base_pivots)
    for _, label, _, polynomial in columns:
        if label == carrier_label:
            audit.add_to_echelon(coefficient_vector(polynomial), echelon, pivots)
    individual_ranks.append(len(echelon) - dS_rank)

echelon = [row[:] for row in base_echelon]
pivots = list(base_pivots)
selected_labels = []
for filtration_degree, carrier_label, target_label, polynomial in sorted(
    columns, key=lambda item: (item[0], item[1], item[2])
):
    if audit.add_to_echelon(coefficient_vector(polynomial), echelon, pivots):
        selected_labels.append(f"{carrier_label}|{target_label}")

wedge_quotient_rank = len(echelon) - dS_rank
require(wedge_quotient_rank == 88, "positive rank-88 minor")
require(len(selected_labels) == 88, "selected rank-88 columns")

print("scope=descended_target_two_forms;fixed_x_first_normal_slice;not_full_THM4067")
print("source=Omega^2_T;T=K[b,c,e]/(c^2e-b(b+4));target=E_x=K[x]/dS")
print("map=tau(dP_wedge_dQ)=Pprime*gamma1(Q)-Qprime*gamma1(P) mod dS")
print("retained_base_wedge=Cprime*gamma1(E)-Eprime*gamma1(C)=(12,12,12)")
print("exact_upper_bound=Lambda(tau(Omega))=0_for_every_descended_target_two_form")
print(f"good_reduction=p{PRIME},alpha{audit.ALPHA_MOD};target_cutoff={TARGET_CUTOFF}")
print(
    f"ambient_degree={ambient_degree};ambient_dimension={ambient_dimension};"
    f"dS_rank={dS_rank};E_x_dimension={ambient_dimension-dS_rank}"
)
print("individual_BCE_quotient_ranks=" + ",".join(map(str, individual_ranks)))
print(
    f"combined_wedge_quotient_rank={wedge_quotient_rank};"
    "positive_minor_plus_exact_upper_bound=ker_Lambda"
)
print("selected_column_count=" + str(len(selected_labels)))
print(
    "selected_column_hash="
    + hashlib.sha256(";".join(selected_labels).encode("ascii")).hexdigest()
)
print("seminormal_class_rprime=missed_because_Lambda(rprime)_nonzero")
print("destroyed_information=representative_level_gamma1_without_wedge_pairing")
print("needed_sidecar=non_descended_primitive_or_independent_target_coordinate")
print("scout_check_count=" + str(CHECK_COUNT))
print("RESULT=PASS")
