"""Exact positive-minor saturation for the exceptional-quartic Jacobian image.

The restriction algebra data are loaded from THM-3703 and converted to four
FLINT rational polynomials, representing K[X] in the power basis
(1, alpha, alpha^2, alpha^3).  A good-reduction determinant is used only in
the sound direction: a nonzero reduction certifies a nonzero exact K-minor.
"""

import contextlib
import hashlib
import importlib.util
import io
from pathlib import Path
import subprocess
import types

import sympy as sp
from flint import fmpq, fmpq_mat, fmpq_poly, nmod_mat


SAGBI_PATH = Path(
    "04-computation/jc2_russell_cylinder_exceptional_quartic_sagbi_module_thm3703.py"
)


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


def kpoly_constant(value):
    return tuple(fmpq_poly([item]) for item in value)


def kpoly_add(left, right):
    return tuple(a + b for a, b in zip(left, right))


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
    temporary = [fmpq_poly([]) for _ in range(7)]
    for i, coefficient_polynomial in enumerate(polynomial):
        if not coefficient_polynomial:
            continue
        for j, scalar in enumerate(element):
            if scalar:
                temporary[i + j] += coefficient_polynomial * scalar
    for degree in range(6, 3, -1):
        coefficient_polynomial = temporary[degree]
        if coefficient_polynomial:
            for offset, relation in enumerate(RELATION):
                temporary[degree - 4 + offset] += coefficient_polynomial * relation
    return tuple(temporary[:4])


def kpoly_diff(polynomial):
    return tuple(item.derivative() for item in polynomial)


def kpoly_coeff(polynomial, degree):
    return tuple(item[degree] if degree <= item.degree() else ZERO_Q for item in polynomial)


def kpoly_degree(polynomial):
    degrees = [item.degree() for item in polynomial if item]
    return max(degrees) if degrees else -1


def kpoly_is_zero(polynomial):
    return not any(polynomial)


def kpoly_monic(polynomial):
    require(not kpoly_is_zero(polynomial), "monic zero")
    answer = kpoly_mul_elem(
        polynomial, elem_inv(kpoly_coeff(polynomial, kpoly_degree(polynomial)))
    )
    require(kpoly_coeff(answer, kpoly_degree(answer)) == elem_one(), "monic residual")
    return answer


def kpoly_evaluate(polynomial, point):
    answer = elem_zero()
    for degree in range(kpoly_degree(polynomial), -1, -1):
        answer = elem_add(tuple(point * item for item in answer), kpoly_coeff(polynomial, degree))
    return answer


def from_sympy_poly(polynomial):
    if not polynomial:
        return kpoly_zero()
    rows = [[ZERO_Q] * (polynomial.degree() + 1) for _ in range(4)]
    for degree in range(polynomial.degree() + 1):
        coordinates = sagbi.field_coordinates(
            polynomial.get((degree,), sagbi.field.zero)
        )
        for coordinate, value in enumerate(coordinates):
            rows[coordinate][degree] = as_fmpq(value)
    return tuple(fmpq_poly(row) for row in rows)


Z = from_sympy_poly(sagbi.Z)
C = from_sympy_poly(sagbi.C)
E = from_sympy_poly(sagbi.E)
C_prime = kpoly_diff(C)
E_prime = kpoly_diff(E)
module_basis = tuple(from_sympy_poly(item) for item in sagbi.module_basis)
require(tuple(kpoly_degree(item) for item in module_basis) == sagbi.EXPECTED_APERY, "S Apéry basis")


# THM-3687 gives 1=C'G1-E'F1, so each p_r belongs to the image module.
# Together with C'p_r and E'p_r these 54 elements generate
# I=C'S+E'S over K[Z].  Interleaving by residue is a deterministic compiler
# choice; it saves ten reductions relative to the three-block ordering.
initial_generators = []
for residue, basis in enumerate(module_basis):
    initial_generators.extend(
        (
            (f"S:{residue}", basis),
            (f"C:{residue}", kpoly_mul(C_prime, basis)),
            (f"E:{residue}", kpoly_mul(E_prime, basis)),
        )
    )
require(len(initial_generators) == 54, "initial generator count")


z_powers = {0: kpoly_constant(elem_one()), 1: Z}


def z_power(power):
    if power not in z_powers:
        z_powers[power] = kpoly_mul(z_power(power - 1), Z)
    return z_powers[power]


points = (fmpq(-1), fmpq(0), fmpq(1))


def retained_lambda(polynomial):
    values = tuple(kpoly_evaluate(polynomial, point) for point in points)
    weighted = elem_zero()
    for weight, value in zip((fmpq(5, 18), fmpq(-1), fmpq(13, 18)), values):
        weighted = elem_add(weighted, tuple(weight * item for item in value))
    return weighted


require(
    tuple(kpoly_evaluate(Z, point) for point in points)
    == (elem_zero(), elem_zero(), elem_zero()),
    "Z retained roots",
)
three = (fmpq(3), ZERO_Q, ZERO_Q, ZERO_Q)
require(
    tuple(kpoly_evaluate(C_prime, point) for point in points)
    == (three, three, three),
    "C-prime retained values",
)
require(
    tuple(kpoly_evaluate(E_prime, point) for point in points)
    == tuple((fmpq(value), ZERO_Q, ZERO_Q, ZERO_Q) for value in (-9, 4, 9)),
    "E-prime retained values",
)
require(
    all(retained_lambda(generator) == elem_zero() for _, generator in initial_generators),
    "image generators lie in retained kernel",
)


X_POLY = (fmpq_poly([0, 1]), ZERO_POLY, ZERO_POLY, ZERO_POLY)
require(retained_lambda(kpoly_constant(elem_one())) == elem_zero(), "Lambda one")
require(
    retained_lambda(X_POLY) == (fmpq(4, 9), ZERO_Q, ZERO_Q, ZERO_Q),
    "Lambda X",
)

x_powers = {1: X_POLY}


def x_power(power):
    if power == 0:
        return kpoly_constant(elem_one())
    if power not in x_powers:
        x_powers[power] = kpoly_mul(x_power(power - 1), X_POLY)
    return x_powers[power]


# Since Lambda(X^r)=1 for positive even r and 4/9 for odd r, the
# K[Z]-kernel has this transparent basis in residue order modulo 18.
explicit_basis = [kpoly_constant(elem_one()), kpoly_mul(X_POLY, Z)]
for residue in range(2, 18):
    coefficient = fmpq(9, 4) if residue % 2 == 0 else ONE_Q
    explicit_basis.append(
        kpoly_sub(
            x_power(residue),
            kpoly_mul_elem(X_POLY, (coefficient, ZERO_Q, ZERO_Q, ZERO_Q)),
        )
    )

image_apery = (0, 19, *range(2, 18))
require(tuple(kpoly_degree(item) for item in explicit_basis) == image_apery, "explicit basis degrees")
require(all(retained_lambda(item) == elem_zero() for item in explicit_basis), "explicit kernel basis")
require(
    tuple(retained_lambda(x_power(power)) for power in range(1, 18))
    == tuple(
        (fmpq(1) if power % 2 == 0 else fmpq(4, 9), ZERO_Q, ZERO_Q, ZERO_Q)
        for power in range(1, 18)
    ),
    "Lambda power pattern",
)


# Filter the 54 K[Z]-generators through degree 150.  A positive determinant
# in one good split fibre is a certificate that the corresponding exact
# determinant in K is nonzero.  This is not a characteristic-zero inference
# from modular rank failure.
CUTOFF = 150
PRIME = 137
ROOT = 44
F_QUARTIC_AT_ROOT = (
    72783360 * ROOT**4
    - 77822208 * ROOT**3
    - 28419741 * ROOT**2
    + 7849770 * ROOT
    - 1276420
)
require(72783360 % PRIME != 0, "good reduction leading coefficient")
require(F_QUARTIC_AT_ROOT % PRIME == 0, "good reduction split root")

filtered_columns = []
filtered_metadata = []
for label, generator in initial_generators:
    power = 0
    polynomial = generator
    while kpoly_degree(polynomial) <= CUTOFF:
        filtered_columns.append(polynomial)
        filtered_metadata.append((label, power))
        power += 1
        polynomial = kpoly_mul(polynomial, Z)
require(len(filtered_columns) == 165, "filtered column count")


def rational_mod(value, prime):
    numerator = int(value.numerator)
    denominator = int(value.denominator)
    require(denominator % prime != 0, "good reduction denominator")
    return numerator % prime * pow(denominator % prime, -1, prime) % prime


def element_mod(value, prime=PRIME, root=ROOT):
    return sum(
        rational_mod(coordinate, prime) * pow(root, power, prime)
        for power, coordinate in enumerate(value)
    ) % prime


def matrix_mod(columns, cutoff):
    data = []
    for degree in range(cutoff + 1):
        data.extend(element_mod(kpoly_coeff(column, degree)) for column in columns)
    return nmod_mat(cutoff + 1, len(columns), data, PRIME)


full_matrix_mod = matrix_mod(filtered_columns, CUTOFF)
_full_reduced, full_rank = full_matrix_mod.rref()
require(full_rank == 150, "good reduction full rank")

selected_columns = tuple(range(150))
selected_rows = tuple(range(150))
selected_square_mod = nmod_mat(
    150,
    150,
    [
        int(full_matrix_mod[row, column])
        for row in selected_rows
        for column in selected_columns
    ],
    PRIME,
)
minor_mod = int(selected_square_mod.det())
require(minor_mod == 73, "good reduction selected determinant")
require(
    all(retained_lambda(filtered_columns[column]) == elem_zero() for column in selected_columns),
    "selected exact columns lie in retained kernel",
)

# Hostile insufficient-cutoff control: the same generator grammar has not yet
# saturated the retained hyperplane at degree 100.
HOSTILE_CUTOFF = 100
hostile_columns = [
    column for column in filtered_columns if kpoly_degree(column) <= HOSTILE_CUTOFF
]
hostile_matrix_mod = matrix_mod(hostile_columns, HOSTILE_CUTOFF)
_hostile_reduced, hostile_rank = hostile_matrix_mod.rref()
require((len(hostile_columns), hostile_rank) == (69, 69), "insufficient cutoff control")
require(HOSTILE_CUTOFF + 1 - hostile_rank == 32, "insufficient cutoff codimension")

filtered_kernel_dimension = CUTOFF
require(full_rank == filtered_kernel_dimension, "matching filtered dimensions")
codimension = sum((degree - residue) // 18 for residue, degree in enumerate(image_apery))
require(codimension == 1, "global kernel codimension")


def elem_text(value):
    return ",".join(str(item) for item in value)


def polynomial_hash(polynomial):
    payload = ";".join(
        f"{degree}:{elem_text(kpoly_coeff(polynomial, degree))}"
        for degree in range(kpoly_degree(polynomial) + 1)
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def packet_hash(polynomials):
    return hashlib.sha256(
        ";".join(polynomial_hash(item) for item in polynomials).encode("ascii")
    ).hexdigest()


def index_hash(indices):
    return hashlib.sha256(",".join(map(str, indices)).encode("ascii")).hexdigest()


def exact_square_hash():
    payload = ";".join(
        elem_text(kpoly_coeff(filtered_columns[column], row))
        for row in selected_rows
        for column in selected_columns
    )
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


print("field=THM3683_exceptional_quartic")
print("restriction_module=THM3703_direct_sum_KZ_rank18;degZ=18")
print(
    f"image_generators_KZ=54;filtered_cutoff={CUTOFF};filtered_columns={len(filtered_columns)};"
    f"selected_K_square=150"
)
print(
    f"good_reduction=p{PRIME}:alpha{ROOT};full_rank={full_rank};"
    f"selected_det_mod_p={minor_mod};omitted_row=150"
)
print(
    f"selected_columns_hash={index_hash(selected_columns)};"
    f"selected_rows_hash={index_hash(selected_rows)};exact_square_hash={exact_square_hash()}"
)
print(
    f"hostile_cutoff={HOSTILE_CUTOFF};hostile_columns={len(hostile_columns)};"
    f"hostile_rank={hostile_rank};hostile_codimension={HOSTILE_CUTOFF + 1 - hostile_rank}"
)
print("image_apery_mod18=" + ",".join(map(str, image_apery)))
print(f"image_codimension={codimension};filtered_kernel_dimension={filtered_kernel_dimension}")
print("retained_points=-1,0,1;Cprime_values=3,3,3;Eprime_values=-9,4,9")
print("Lambda=5*ev(-1)/18-ev(0)+13*ev(1)/18;Lambda(1)=0;Lambda(X)=4/9")
print("explicit_kernel_basis=1,XZ,X^r-X(r_odd_3_to_17),X^r-9X/4(r_even_2_to_16)")
print(f"initial_generator_hash={packet_hash(item for _, item in initial_generators)}")
print(f"explicit_kernel_basis_hash={packet_hash(explicit_basis)}")
print("image_L0=kernel_Lambda;all_four_embeddings=uniform;exact_characteristic_zero=1")
print("stagewise_gate=J_n_debt_in_image_iff_Lambda(debt)=0;degree_cutoff=none")
print("NO_CLAIM=all_order_solution_or_algebraization_or_global_pair_or_Keller_counterexample")
print("RESULT=PASS")
