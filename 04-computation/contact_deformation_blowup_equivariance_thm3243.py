"""Exact companion for THM-3243's contact-deformation blowup theorem."""

from collections import deque
from itertools import product
import ast
import hashlib
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCIES = {
    ROOT
    / (
        "01-canon/theorems/"
        "THM-3234-singer-owner-compactification-and-pointed-heisenberg-carrier-gate.md"
    ): "ef77a1f8fce16eb851eb38d5110a61ab73aa693f2d0ee9e11a912aa4fc302c87",
    ROOT
    / (
        "01-canon/theorems/"
        "THM-3240-exact-address-heisenberg-clutch-on-carrier-imbalance.md"
    ): "7d23f2920adfecb17d8a149aada08a8e34215111546eac63b40570e898e14f14",
    ROOT
    / (
        "01-canon/theorems/"
        "THM-3241-finite-field-contact-Singer-realization-and-order-gate.md"
    ): "05ba274b802f36272416f8d8ede18e2ebc5b5073cf459a27daf53409b7c75c81",
}


def lf_sha256(path):
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


for dependency, expected_hash in DEPENDENCIES.items():
    require(lf_sha256(dependency) == expected_hash, "dependency hash: %s" % dependency.name)

syntax_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax_tree))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax_tree)
)
require(assert_nodes == 0, "assert statements are optimization-sensitive")
require(float_literals == 0, "floating literals are forbidden")


def matrix_multiply(left, right, prime):
    a, b, c, d = left
    e, f, g, h = right
    return (
        (a * e + b * g) % prime,
        (a * f + b * h) % prime,
        (c * e + d * g) % prime,
        (c * f + d * h) % prime,
    )


def matrix_vector(matrix, vector, prime):
    a, b, c, d = matrix
    x, y = vector
    return ((a * x + b * y) % prime, (c * x + d * y) % prime)


def determinant(matrix, prime):
    a, b, c, d = matrix
    return (a * d - b * c) % prime


def matrix_inverse(matrix, prime):
    a, b, c, d = matrix
    det = determinant(matrix, prime)
    require(det != 0, "singular matrix")
    inverse_det = pow(det, -1, prime)
    return (
        d * inverse_det % prime,
        -b * inverse_det % prime,
        -c * inverse_det % prime,
        a * inverse_det % prime,
    )


def matrix_power(matrix, exponent, prime):
    result = (1, 0, 0, 1)
    base = matrix
    while exponent:
        if exponent % 2:
            result = matrix_multiply(result, base, prime)
        base = matrix_multiply(base, base, prime)
        exponent //= 2
    return result


def matrix_order(matrix, prime):
    result = (1, 0, 0, 1)
    for exponent in range(1, prime**2 + 1):
        result = matrix_multiply(result, matrix, prime)
        if result == (1, 0, 0, 1):
            return exponent
    raise RuntimeError("matrix order exceeds finite-field bound")


def scalar_matrix_value(matrix, prime):
    a, b, c, d = matrix
    if b % prime == 0 and c % prime == 0 and a % prime == d % prime:
        return a % prime
    return None


def projective_direction(vector, prime):
    x, y = vector
    require((x, y) != (0, 0), "zero has no projective direction")
    if x % prime:
        inverse_x = pow(x, -1, prime)
        return (1, y * inverse_x % prime)
    return (0, 1)


def projective_order(matrix, prime):
    for exponent in range(1, prime**2 + 1):
        if scalar_matrix_value(matrix_power(matrix, exponent, prime), prime) is not None:
            return exponent
    raise RuntimeError("projective order exceeds finite-field bound")


def scalar_order(value, prime):
    require(value % prime != 0, "zero scalar")
    current = 1
    for exponent in range(1, prime):
        current = current * value % prime
        if current == 1:
            return exponent
    raise RuntimeError("scalar order exceeds base-field bound")


def generated_matrix_group(generators, prime):
    identity = (1, 0, 0, 1)
    moves = list(generators) + [matrix_inverse(g, prime) for g in generators]
    seen = {identity}
    queue = deque([identity])
    while queue:
        current = queue.popleft()
        for move in moves:
            following = matrix_multiply(current, move, prime)
            if following not in seen:
                seen.add(following)
                queue.append(following)
    return seen


def heisenberg_action(element, vector, prime):
    x, y, z = element
    r, w = vector
    return ((r + x) % prime, (w + z - y * r) % prime)


def direction_action(element, direction, prime):
    _, y, _ = element
    derivative = (1, 0, -y % prime, 1)
    representative = (0, 1) if direction == (0, 1) else direction
    return projective_direction(matrix_vector(derivative, representative, prime), prime)


def flag_action(element, flag, prime):
    vector, direction = flag
    return (
        heisenberg_action(element, vector, prime),
        direction_action(element, direction, prime),
    )


def poly_clean(poly, prime):
    return {monomial: coefficient % prime for monomial, coefficient in poly.items() if coefficient % prime}


def poly_add(left, right, prime):
    out = dict(left)
    for monomial, coefficient in right.items():
        out[monomial] = (out.get(monomial, 0) + coefficient) % prime
    return poly_clean(out, prime)


def poly_scale(poly, scalar, prime):
    return poly_clean({monomial: scalar * coefficient for monomial, coefficient in poly.items()}, prime)


def poly_multiply(left, right, prime):
    out = {}
    for (left_x, left_y), left_coefficient in left.items():
        for (right_x, right_y), right_coefficient in right.items():
            monomial = (left_x + right_x, left_y + right_y)
            out[monomial] = (
                out.get(monomial, 0) + left_coefficient * right_coefficient
            ) % prime
    return poly_clean(out, prime)


def poly_power(poly, exponent, prime):
    out = {(0, 0): 1}
    base = poly
    while exponent:
        if exponent % 2:
            out = poly_multiply(out, base, prime)
        base = poly_multiply(base, base, prime)
        exponent //= 2
    return out


def affine_component(matrix, translation, row, prime):
    if row == 0:
        left, right = matrix[0], matrix[1]
    else:
        left, right = matrix[2], matrix[3]
    return poly_clean(
        {(1, 0): left, (0, 1): right, (0, 0): translation[row]},
        prime,
    )


def centre_ideal_identity(matrix, translation, row, prime):
    component = affine_component(matrix, translation, row, prime)
    left = poly_add(poly_power(component, prime, prime), poly_scale(component, -1, prime), prime)
    if row == 0:
        coefficients = (matrix[0], matrix[1])
    else:
        coefficients = (matrix[2], matrix[3])
    x_frobenius = {(prime, 0): 1, (1, 0): -1}
    y_frobenius = {(0, prime): 1, (0, 1): -1}
    right = poly_add(
        poly_scale(x_frobenius, coefficients[0], prime),
        poly_scale(y_frobenius, coefficients[1], prime),
        prime,
    )
    return left == right


prime_rows = []
for prime in (3, 5, 7, 13):
    elements = list(product(range(prime), repeat=3))
    points = set(product(range(prime), repeat=2))
    origin_orbit = {heisenberg_action(element, (0, 0), prime) for element in elements}
    origin_stabilizer = [
        element for element in elements if heisenberg_action(element, (0, 0), prime) == (0, 0)
    ]
    require(origin_orbit == points, "Heisenberg origin orbit")
    require(
        origin_stabilizer == [(0, y, 0) for y in range(prime)],
        "Heisenberg origin stabilizer",
    )

    vertical = ((0, 0), (0, 1))
    nonvertical = ((0, 0), (1, 0))
    vertical_orbit = {flag_action(element, vertical, prime) for element in elements}
    nonvertical_orbit = {flag_action(element, nonvertical, prime) for element in elements}
    all_flags = {
        (point, direction)
        for point in points
        for direction in [(0, 1)] + [(1, slope) for slope in range(prime)]
    }
    require(len(vertical_orbit) == prime**2, "vertical flag orbit")
    require(len(nonvertical_orbit) == prime**3, "nonvertical flag orbit")
    require(vertical_orbit.isdisjoint(nonvertical_orbit), "flag orbit disjointness")
    require(vertical_orbit | nonvertical_orbit == all_flags, "full flag decomposition")
    require(
        sum(flag_action(element, vertical, prime) == vertical for element in elements) == prime,
        "vertical stabilizer",
    )
    require(
        sum(flag_action(element, nonvertical, prime) == nonvertical for element in elements) == 1,
        "nonvertical stabilizer",
    )
    prime_rows.append(
        (prime, len(origin_orbit), len(origin_stabilizer), len(vertical_orbit), len(nonvertical_orbit))
    )


prime = 13
identity = (1, 0, 0, 1)
linear_generators = [
    (1, 1, 0, 1),
    (0, 1, 1, 0),
    (2, 0, 0, 1),
]
generated_gl = generated_matrix_group(linear_generators, prime)
expected_gl_order = (prime**2 - 1) * (prime**2 - prime)
require(len(generated_gl) == expected_gl_order, "GL2 generator closure")
require(all(determinant(matrix, prime) for matrix in generated_gl), "GL2 nonsingularity")

directions = {(0, 1)} | {(1, slope) for slope in range(prime)}
direction_orbit = {
    projective_direction(matrix_vector(matrix, (1, 0), prime), prime)
    for matrix in generated_gl
}
flag_stabilizer = sum(
    projective_direction(matrix_vector(matrix, (1, 0), prime), prime) == (1, 0)
    for matrix in generated_gl
)
require(direction_orbit == directions, "GL2 direction transitivity")
require(flag_stabilizer == prime * (prime - 1) ** 2, "AGL flag stabilizer")

affine_generators = [
    (identity, (1, 0)),
    (identity, (0, 1)),
    (linear_generators[0], (0, 0)),
    (linear_generators[1], (0, 0)),
    (linear_generators[2], (0, 0)),
]
ideal_checks = 0
for matrix, translation in affine_generators:
    for row in (0, 1):
        require(
            centre_ideal_identity(matrix, translation, row, prime),
            "full rational-centre ideal identity",
        )
        ideal_checks += 1

singer = (1, 4, 2, 1)
require(determinant(singer, prime) != 0, "Singer matrix singular")
singer_order = matrix_order(singer, prime)
singer_projective_order = projective_order(singer, prime)
radial_matrix = matrix_power(singer, prime + 1, prime)
radial_scalar = scalar_matrix_value(radial_matrix, prime)
require(singer_order == prime**2 - 1, "Singer order")
require(singer_projective_order == prime + 1, "Singer projective order")
require(radial_scalar is not None, "Singer radial power is not scalar")
radial_order = scalar_order(radial_scalar, prime)
require(radial_order == prime - 1, "Singer radial order")

singer_direction_orbit = []
direction = (1, 0)
for _ in range(prime + 1):
    require(direction not in singer_direction_orbit, "Singer direction repeats early")
    singer_direction_orbit.append(direction)
    direction = projective_direction(matrix_vector(singer, direction, prime), prime)
require(direction == (1, 0), "Singer direction orbit does not close")
require(set(singer_direction_orbit) == directions, "Singer misses a rational direction")

# THM-3241's p=13 order-two chart: c=(2x)H modulo x^2-2.
contact_classes = set()
for helper_constant, helper_linear in product(range(prime), repeat=2):
    contact_constant = 4 * helper_linear % prime
    contact_linear = 2 * helper_constant % prime
    contact_classes.add((contact_constant, contact_linear))
require(len(contact_classes) == prime**2, "contact chart is not bijective")
require((4 * 10 % prime, 2) == (1, 2), "explicit Singer helper")

full_centre_points = prime**2
rational_flags = prime**2 * (prime + 1)
require(rational_flags == 2366, "p=13 rational flag count")
require(prime**2 + prime**3 == rational_flags, "Heisenberg flag split")
require(expected_gl_order * prime**2 // rational_flags == 1872, "p=13 flag stabilizer")

print("THM-3243 exact companion")
print("dependency_hashes=%d assert_nodes=%d float_literals=%d" % (len(DEPENDENCIES), assert_nodes, float_literals))
for row in prime_rows:
    print(
        "p=%d point_orbit=%d point_stabilizer=%d vertical_flags=%d nonvertical_flags=%d"
        % row
    )
print(
    "p13_gl2=%d generators=%d directions=%d flag_stabilizer=%d"
    % (len(generated_gl), len(linear_generators), len(direction_orbit), flag_stabilizer)
)
print(
    "p13_centre_ideal_checks=%d centre_points=%d rational_flags=%d"
    % (ideal_checks, full_centre_points, rational_flags)
)
print(
    "p13_singer_order=%d projective_order=%d radial_scalar=%d radial_order=%d"
    % (singer_order, singer_projective_order, radial_scalar, radial_order)
)
print(
    "p13_contact_classes=%d delayed=1 exact=%d helper=(1,10) contact=(1,2)"
    % (len(contact_classes), len(contact_classes) - 1)
)
print("THM-3243 PASS")
