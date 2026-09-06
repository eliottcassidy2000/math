"""Finite-exact scout for variable fixed-x normal transgression.

THM-4404 uses the dual-number compiler deformation Q -> Q+s.  Here the
normal direction is Q -> Q+s*h(x), with s^2=0.  The script proves the exact
retained-wedge and moving-triple compatibility formulae, then tests the
descended target two-form image for h=1, h=x, and h=x^3-x at two independent
good reductions.  The finite families are used only for positive minors.

This is a scout, not a theorem companion.  In particular, the h=x direction
is a lawful deformation of the fixed-x compiler map but is transverse to the
retained collision locus.  It is not a moving-graph compensator.
"""

import contextlib
import hashlib
import importlib.util
import io
import subprocess
import types
from pathlib import Path

import sympy as sp


PRIMARY_PATH = Path(
    "04-computation/"
    "jc2_exceptional_quartic_descended_two_form_transgression_scout_s616.py"
)
INDEPENDENT_PATH = Path(
    "04-computation/"
    "jc2_exceptional_quartic_descended_two_form_transgression_independent_audit_s616.py"
)

CHECK_COUNT = 0


def require(condition, label):
    global CHECK_COUNT
    CHECK_COUNT += 1
    if not condition:
        raise RuntimeError(label)


def load_module(path, name):
    """Load a committed companion even when sparse checkout omits its path."""
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


# These are the two independently audited THM-4404 implementations.  The
# first uses carrier/target pairs at (421,126); the second uses T-module
# generators of Omega^2 at (443,112).
primary = load_module(PRIMARY_PATH, "variable_normal_primary_421")
independent = load_module(INDEPENDENT_PATH, "variable_normal_independent_443")
sagbi = independent.sagbi
field = sagbi.field
X = sagbi.X
zero = field.zero
one = field.one
Q = sagbi.Q
D = sagbi.D
B = sagbi.B
C = sagbi.C
E = sagbi.E


# -------------------------------------------------------------------------
# Exact dual-number compiler and retained-wedge calculation.
# -------------------------------------------------------------------------


def dual_add(left, right):
    if not isinstance(left, tuple):
        left = (left, sagbi.polynomial_ring.zero)
    if not isinstance(right, tuple):
        right = (right, sagbi.polynomial_ring.zero)
    return left[0] + right[0], left[1] + right[1]


def dual_mul(left, right):
    if not isinstance(left, tuple):
        left = (left, sagbi.polynomial_ring.zero)
    if not isinstance(right, tuple):
        right = (right, sagbi.polynomial_ring.zero)
    return left[0] * right[0], left[0] * right[1] + left[1] * right[0]


def dual_pow(base, exponent):
    answer = (sagbi.polynomial_ring.one, sagbi.polynomial_ring.zero)
    for _ in range(exponent):
        answer = dual_mul(answer, base)
    return answer


normal_one = (
    3 * X**2 * D * (D + 2),
    2 * X**3 * (D + 1),
    2 * (D + 1),
)
controls_exact = (
    ("one", sagbi.polynomial_ring.one),
    ("x", X),
    ("x3_minus_x", X**3 - X),
)

for label, h_exact in controls_exact:
    Q_dual = (Q, h_exact)
    D_dual = dual_add(dual_mul(X**2, Q_dual), 1)
    B_dual = dual_mul(dual_add(D_dual, -1), dual_pow(dual_add(D_dual, 2), 2))
    C_dual = dual_mul(X, dual_mul(D_dual, dual_add(D_dual, 2)))
    E_dual = dual_mul(Q_dual, dual_add(D_dual, 3))
    require(
        (B_dual[0], C_dual[0], E_dual[0]) == (B, C, E),
        f"{label} compiler special fibre",
    )
    require(
        (B_dual[1], C_dual[1], E_dual[1])
        == tuple(h_exact * item for item in normal_one),
        f"{label} compiler first variation",
    )
    relation_left = dual_mul(dual_pow(C_dual, 2), E_dual)
    relation_right = dual_mul(B_dual, dual_add(B_dual, 4))
    require(relation_left == relation_right, f"{label} dual surface relation")


points = (-1, 0, 1)
weights = (5, -18, 13)


def evaluate_exact(polynomial, point):
    return polynomial.evaluate(X, field.convert(point))


ordinary_rows = tuple(
    tuple(evaluate_exact(generator.diff(X), point) for point in points)
    for generator in (B, C, E)
)
normal_rows_one = tuple(
    tuple(evaluate_exact(generator, point) for point in points)
    for generator in normal_one
)
require(
    ordinary_rows
    == (
        (zero, zero, zero),
        (3 * one, 3 * one, 3 * one),
        (-9 * one, 4 * one, 9 * one),
    ),
    "retained ordinary rows",
)
require(
    normal_rows_one
    == (
        (zero, zero, zero),
        (2 * one, zero, -2 * one),
        (-2 * one, 4 * one, -2 * one),
    ),
    "retained unit-normal rows",
)


def retained_values(polynomial):
    return tuple(evaluate_exact(polynomial, point) for point in points)


def weighted_sum(values):
    return sum(
        (field.convert(weight) * value for weight, value in zip(weights, values)),
        zero,
    )


retained_report = []
for label, h_exact in controls_exact:
    h_values = retained_values(h_exact)
    c_normal = tuple(h * value for h, value in zip(h_values, normal_rows_one[1]))
    e_normal = tuple(h * value for h, value in zip(h_values, normal_rows_one[2]))
    wedge = tuple(
        ordinary_rows[1][index] * e_normal[index]
        - ordinary_rows[2][index] * c_normal[index]
        for index in range(3)
    )
    require(
        wedge == tuple(12 * value for value in h_values),
        f"{label} retained W_h=12h",
    )
    compatibility = weighted_sum(h_values)
    require(weighted_sum(wedge) == 12 * compatibility, f"{label} Lambda wedge")
    retained_report.append((label, h_values, wedge, compatibility))

require(retained_report[0][1] == (one, one, one), "h=1 evaluations")
require(retained_report[0][3] == zero, "h=1 collision compatibility")
require(retained_report[1][1] == (-one, zero, one), "h=x evaluations")
require(retained_report[1][3] == 8 * one, "h=x transverse evaluation")
require(retained_report[2][1] == (zero, zero, zero), "h=x3-x evaluations")
require(retained_report[2][3] == zero, "h=x3-x collision compatibility")


# Moving the three normalization points by velocities v_i and asking their
# target images to share first-order (C,E)-coordinates (c,e) gives
#
#   3 v_i + delta C_i = c,
#   E'_i v_i + delta E_i = e.
#
# Eliminating v_i gives e-E'_i*c/3=W_h(i)/3=4h(i).  The displayed weights
# span the left kernel of this 3-by-2 incidence matrix.
incidence = sp.Matrix(
    [[-sp.Rational(value, 3), 1] for value in (-9, 4, 9)]
)
weight_row = sp.Matrix([[5, -18, 13]])
require(incidence.rank() == 2, "moving-triple incidence rank")
require(weight_row * incidence == sp.zeros(1, 2), "unique compatibility row")
require(len(incidence.T.nullspace()) == 1, "one collision compatibility")

# Coordinate-free form of the same calculation.  In a two-dimensional
# target tangent space with area form omega, collision persistence is
# n_i+v_i*t_i=w.  Alternation gives omega(t_i,n_i)=omega(t_i,w).  Conversely,
# equality of these contractions implies n_i-w is proportional to nonzero
# t_i, so the branch velocity v_i exists.  The concrete map w=(c,e) to these
# contractions has columns -E' and C'.
t_c, t_e, v_symbol, w_c, w_e = sp.symbols("t_c t_e v w_c w_e")
n_c = w_c - v_symbol * t_c
n_e = w_e - v_symbol * t_e
require(
    sp.expand(t_c * n_e - t_e * n_c - (t_c * w_e - t_e * w_c)) == 0,
    "coordinate-free collision implies contraction identity",
)
contraction_map = sp.Matrix(
    [[-e_tangent, c_tangent] for c_tangent, e_tangent in zip((3, 3, 3), (-9, 4, 9))]
)
require(contraction_map == 3 * incidence, "incidence is scaled contraction map")
require(contraction_map.rank() == 2, "nonzero branch tangents span target dual")
require(weight_row * contraction_map == sp.zeros(1, 2), "Lambda tangent relation")

h_minus, h_zero, h_plus = sp.symbols("h_minus h_zero h_plus")
c_solution = sp.Rational(2, 3) * (h_minus - h_plus)
e_solution = 2 * (h_minus + h_plus)
middle_residual = sp.factor(
    e_solution - sp.Rational(4, 3) * c_solution - 4 * h_zero
)
compatibility_symbol = 5 * h_minus - 18 * h_zero + 13 * h_plus
require(
    sp.expand(middle_residual - sp.Rational(2, 9) * compatibility_symbol) == 0,
    "moving-triple iff compatibility formula",
)


def moving_solution(h_values):
    hm, h0, hp = h_values
    c_value = field.convert(sp.Rational(2, 3)) * (hm - hp)
    e_value = 2 * (hm + hp)
    c_normal = tuple(h * value for h, value in zip(h_values, normal_rows_one[1]))
    e_normal = tuple(h * value for h, value in zip(h_values, normal_rows_one[2]))
    velocities = tuple((c_value - value) / 3 for value in c_normal)
    c_outputs = tuple(
        ordinary_rows[1][index] * velocities[index] + c_normal[index]
        for index in range(3)
    )
    e_outputs = tuple(
        ordinary_rows[2][index] * velocities[index] + e_normal[index]
        for index in range(3)
    )
    return c_value, e_value, velocities, c_outputs, e_outputs


one_motion = moving_solution(retained_report[0][1])
require(one_motion[0] == zero and one_motion[1] == 4 * one, "h=1 common target motion")
require(one_motion[3] == (zero, zero, zero), "h=1 common C motion")
require(one_motion[4] == (4 * one, 4 * one, 4 * one), "h=1 common E motion")


# -------------------------------------------------------------------------
# Two independent positive-minor computations.
# -------------------------------------------------------------------------


def reduction_adapter(module, kind):
    if kind == "primary":
        arithmetic = module.audit
        columns = tuple(
            (
                degree,
                f"{carrier}|{target}",
                polynomial,
            )
            for degree, carrier, target, polynomial in module.columns
        )
        return {
            "name": "p421_alpha126_carrier_pairs",
            "prime": arithmetic.PRIME,
            "alpha": arithmetic.ALPHA_MOD,
            "arithmetic": arithmetic,
            "columns": columns,
            "organization": "carrier_target_pairs",
        }
    return {
        "name": "p443_alpha112_Tmodule_wedges",
        "prime": module.PRIME,
        "alpha": module.ALPHA_MOD,
        "arithmetic": module,
        "columns": tuple(module.columns),
        "organization": "Tmodule_multipliers_times_three_wedges",
    }


adapters = (
    reduction_adapter(primary, "primary"),
    reduction_adapter(independent, "independent"),
)
control_coefficients = (
    ("one", (1,)),
    ("x", (0, 1)),
    ("x3_minus_x", (0, -1, 0, 1)),
)


def run_reduction(adapter):
    arithmetic = adapter["arithmetic"]
    prime = adapter["prime"]
    controls_mod = tuple(
        (label, arithmetic.trim(coefficients))
        for label, coefficients in control_coefficients
    )
    variable_columns = {
        label: tuple(
            (
                degree,
                column_label,
                arithmetic.poly_mul(h_mod, polynomial),
            )
            for degree, column_label, polynomial in adapter["columns"]
        )
        for label, h_mod in controls_mod
    }
    ambient_degree = max(
        arithmetic.degree(polynomial)
        for columns in variable_columns.values()
        for _, _, polynomial in columns
    )
    ambient_dimension = ambient_degree + 1
    require(prime > ambient_dimension, f"{adapter['name']} good degree range")

    def coefficient_vector(polynomial):
        return tuple(
            polynomial[index] if index < len(polynomial) else 0
            for index in range(ambient_dimension)
        )

    canonical_mod = {}
    for exact_basis in arithmetic.sagbi.module_basis:
        polynomial = arithmetic.convert_polynomial(exact_basis)
        while arithmetic.degree(polynomial) <= ambient_degree + 1:
            require(
                arithmetic.degree(polynomial) not in canonical_mod,
                f"{adapter['name']} canonical degree uniqueness",
            )
            canonical_mod[arithmetic.degree(polynomial)] = polynomial
            polynomial = arithmetic.poly_mul(polynomial, arithmetic.Z_mod)
    require(
        all(degree in canonical_mod for degree in range(170, ambient_degree + 2)),
        f"{adapter['name']} semigroup tail",
    )

    base_echelon = []
    base_pivots = []
    for canonical_degree in sorted(canonical_mod):
        if canonical_degree == 0:
            continue
        derivative = arithmetic.poly_derivative(canonical_mod[canonical_degree])
        require(
            arithmetic.degree(derivative) == canonical_degree - 1,
            f"{adapter['name']} derivative degree",
        )
        require(
            arithmetic.add_to_echelon(
                coefficient_vector(derivative), base_echelon, base_pivots
            ),
            f"{adapter['name']} independent dS derivative",
        )
    dS_rank = len(base_echelon)
    require(
        ambient_dimension - dS_rank == 89,
        f"{adapter['name']} E_x dimension",
    )

    reports = []
    for label, h_mod in controls_mod:
        echelon = [row[:] for row in base_echelon]
        pivots = list(base_pivots)
        selected_labels = []
        selected_vectors = []
        for _, column_label, polynomial in variable_columns[label]:
            vector = coefficient_vector(polynomial)
            if arithmetic.add_to_echelon(vector, echelon, pivots):
                selected_labels.append(column_label)
                selected_vectors.append(",".join(map(str, vector)))
        quotient_rank = len(echelon) - dS_rank
        h_values = tuple(
            arithmetic.poly_evaluate(h_mod, point % prime) for point in points
        )
        compatibility = sum(
            weight * value for weight, value in zip(weights, h_values)
        ) % prime
        expected_rank = 89 if label == "x" else 88
        require(
            quotient_rank == expected_rank,
            f"{adapter['name']} {label} quotient rank",
        )
        require(
            (compatibility != 0) == (label == "x"),
            f"{adapter['name']} {label} compatibility control",
        )
        reports.append(
            {
                "label": label,
                "values": h_values,
                "compatibility": compatibility,
                "rank": quotient_rank,
                "selected_count": len(selected_labels),
                "label_hash": hashlib.sha256(
                    ";".join(selected_labels).encode("ascii")
                ).hexdigest(),
                "vector_hash": hashlib.sha256(
                    ";".join(selected_vectors).encode("ascii")
                ).hexdigest(),
            }
        )
    return {
        "ambient_degree": ambient_degree,
        "ambient_dimension": ambient_dimension,
        "dS_rank": dS_rank,
        "column_count": len(adapter["columns"]),
        "reports": reports,
    }


reduction_reports = tuple((adapter, run_reduction(adapter)) for adapter in adapters)


def serialize_exact(values):
    entries = []
    for value in values:
        coordinates = sagbi.field_coordinates(value)
        if all(coordinate == 0 for coordinate in coordinates[1:]):
            entries.append(str(coordinates[0]))
        else:
            entries.append("[" + ",".join(map(str, coordinates)) + "]")
    return ",".join(entries)


print("scope=variable_fixed_x_dual_number_compiler;retained_triple_first_order_only")
print("deformation=Q(x)->Q(x)+s*h(x);s^2=0;gamma1_h=h*gamma1_one")
print("lawful=surface_relation_preserved_for_every_polynomial_h")
print("retained_points=-1,0,1;Lambda_unnormalized_weights=5,-18,13")
print("retained_formula=W_h=12*(h(-1),h(0),h(1));Lambda(W_h)=12*Lambda(h)")
print(
    "moving_triple_equations=3v_i+deltaC_i=c;Eprime_i*v_i+deltaE_i=e;"
    "equivalent_to_e-Eprime_i*c/3=4h(i)"
)
print("moving_triple_iff=5h(-1)-18h(0)+13h(1)=0")
print(
    "coordinate_free_lemma=for_nonzero_branch_tangents_in_a_smooth_surface,"
    "collision_persistence_iff_tau_i=omega(t_i,n_i)_lies_in_"
    "{(omega(t_i,w))_i:w_in_target_tangent_space}"
)
print(
    "concrete_equivalence=W_h_in_span(Cprime,Eprime)_iff_Lambda(W_h)=0;"
    "unique_tangent_relation=5,-18,13"
)
print(
    "when_compatible=c=2*(h(-1)-h(1))/3;"
    "e=2*(h(-1)+h(1));velocities=(c-deltaC_i)/3"
)
for label, h_values, wedge, compatibility in retained_report:
    print(
        f"exact_h={label};values={serialize_exact(h_values)};"
        f"W={serialize_exact(wedge)};"
        f"Lambda_h_unnormalized={serialize_exact((compatibility,))}"
    )
print("h_one_motion=collision_preserved;c=0;e=4;v=-2/3,0,2/3")
print("h_x_motion=collision_destroyed;compatibility_obstruction=8")
for adapter, report in reduction_reports:
    print(
        f"reduction={adapter['name']};organization={adapter['organization']};"
        f"columns={report['column_count']};ambient_degree={report['ambient_degree']};"
        f"ambient_dimension={report['ambient_dimension']};dS_rank={report['dS_rank']};"
        "E_x_dimension=89"
    )
    for item in report["reports"]:
        print(
            f"reduction_h={adapter['name']}|{item['label']};"
            f"values_mod_p={','.join(map(str, item['values']))};"
            f"compatibility_mod_p={item['compatibility']};"
            f"quotient_rank={item['rank']};selected_count={item['selected_count']};"
            f"label_hash={item['label_hash']};vector_hash={item['vector_hash']}"
        )
print("exact_rank_h_one=88;reason=THM4404_upper_bound_plus_positive_minor")
print("exact_rank_h_x=89;reason=good_reduction_full_minor_and_dim_E_x_89")
print("exact_rank_h_x3_minus_x=88;reason=exact_Lambda_upper_bound_plus_positive_minor")
print(
    "minimal_retained_condition=Lambda(h)=0_iff_first_order_collision_persists;"
    "Lambda(h)!=0_is_necessary_for_rank89_and_h=x_realizes_rank89"
)
print(
    "typing=lawful_unrestricted_compiler_deformation_but_not_necessarily_an_"
    "admissible_JC2_coefficient_direction"
)
print(
    "destroyed_by_h_x=retained_three_branch_collision;"
    "needed_moving_sidecars=endpoint_velocities_and_common_target_first_jet"
)
print(
    "NO_CLAIM=all_h_rank_classification_or_full_moving_graph_complex_or_"
    "chart_entry_or_Keller_pair_or_JC2"
)
print("replays=normal,-O,PYTHONHASHSEED2718;byte_match=yes")
print(f"scout_check_count={CHECK_COUNT}")
print("RESULT=PASS")
