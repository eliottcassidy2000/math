"""Exact companion for THM-3241's finite-field contact Singer gate."""

from math import comb
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
        "THM-3232-root-free-contact-stratum-norm-and-discriminant-power.md"
    ): "dcf42b70cad1bc940aaeea13d7d0a4b981678b8d6ca87e828b06624c242d77a3",
    ROOT
    / (
        "01-canon/theorems/"
        "THM-3234-singer-owner-compactification-and-pointed-heisenberg-carrier-gate.md"
    ): "ef77a1f8fce16eb851eb38d5110a61ab73aa693f2d0ee9e11a912aa4fc302c87",
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


def trim(poly, prime):
    out = [int(value) % prime for value in poly]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def add(left, right, prime):
    degree = max(len(left), len(right))
    out = [0] * degree
    for index in range(degree):
        out[index] = (
            (left[index] if index < len(left) else 0)
            + (right[index] if index < len(right) else 0)
        ) % prime
    return trim(out, prime)


def sub(left, right, prime):
    return add(left, [-value for value in right], prime)


def scale(poly, scalar, prime):
    return trim([scalar * value for value in poly], prime)


def mul(left, right, prime):
    out = [0] * (len(left) + len(right) - 1)
    for i, left_i in enumerate(left):
        for j, right_j in enumerate(right):
            out[i + j] = (out[i + j] + left_i * right_j) % prime
    return trim(out, prime)


def power(poly, exponent, prime):
    out = [1]
    for _ in range(exponent):
        out = mul(out, poly, prime)
    return out


def derivative(poly, prime):
    if len(poly) <= 1:
        return [0]
    return trim([index * poly[index] for index in range(1, len(poly))], prime)


def hasse(poly, order, prime):
    if order >= len(poly):
        return [0]
    return trim(
        [poly[index] * comb(index, order) for index in range(order, len(poly))],
        prime,
    )


def divmod_poly(dividend, divisor, prime):
    dividend = trim(dividend, prime)
    divisor = trim(divisor, prime)
    require(divisor != [0], "division by zero polynomial")
    if len(dividend) < len(divisor):
        return [0], dividend
    quotient = [0] * (len(dividend) - len(divisor) + 1)
    remainder = list(dividend)
    inverse_lead = pow(divisor[-1], -1, prime)
    while remainder != [0] and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1] * inverse_lead % prime
        quotient[shift] = (quotient[shift] + coefficient) % prime
        remainder = sub(
            remainder,
            [0] * shift + [coefficient * value for value in divisor],
            prime,
        )
    return trim(quotient, prime), trim(remainder, prime)


def remainder(poly, modulus, prime):
    return divmod_poly(poly, modulus, prime)[1]


def inverse_mod(poly, modulus, prime):
    old_r, new_r = trim(poly, prime), trim(modulus, prime)
    old_s, new_s = [1], [0]
    while new_r != [0]:
        quotient, next_r = divmod_poly(old_r, new_r, prime)
        old_r, new_r = new_r, next_r
        old_s, new_s = new_s, sub(old_s, mul(quotient, new_s, prime), prime)
    require(len(old_r) == 1 and old_r[0] != 0, "class is not invertible")
    return remainder(scale(old_s, pow(old_r[0], -1, prime), prime), modulus, prime)


def class_mul(left, right, modulus, prime):
    return remainder(mul(left, right, prime), modulus, prime)


def class_power(poly, exponent, modulus, prime):
    out = [1]
    base = remainder(poly, modulus, prime)
    while exponent:
        if exponent % 2:
            out = class_mul(out, base, modulus, prime)
        base = class_mul(base, base, modulus, prime)
        exponent //= 2
    return out


def multiplication_matrix(poly, modulus, prime):
    degree = len(modulus) - 1
    columns = []
    for shift in range(degree):
        basis = [0] * shift + [1]
        value = class_mul(poly, basis, modulus, prime)
        columns.append(value + [0] * (degree - len(value)))
    return [[columns[column][row] for column in range(degree)] for row in range(degree)]


def omega(poly_f, poly_g, order, prime):
    return sub(
        mul(hasse(poly_f, order, prime), derivative(poly_g, prime), prime),
        mul(hasse(poly_g, order, prime), derivative(poly_f, prime), prime),
        prime,
    )


def contact_numerator_denominator(poly_f, poly_g, stratum, order, prime):
    numerator = remainder(scale(omega(poly_f, poly_g, order, prime), -1, prime), stratum, prime)
    denominator = remainder(
        mul(derivative(poly_f, prime), derivative(poly_g, prime), prime),
        stratum,
        prime,
    )
    contact = class_mul(numerator, inverse_mod(denominator, stratum, prime), stratum, prime)
    return numerator, denominator, contact


def determinant_mod(matrix, prime):
    size = len(matrix)
    work = [[value % prime for value in row] for row in matrix]
    result = 1
    sign = 1
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column] % prime),
            None,
        )
        if pivot is None:
            return 0
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign = -sign
        pivot_value = work[column][column] % prime
        result = result * pivot_value % prime
        inverse_pivot = pow(pivot_value, -1, prime)
        for row in range(column + 1, size):
            factor = work[row][column] * inverse_pivot % prime
            for entry in range(column, size):
                work[row][entry] = (
                    work[row][entry] - factor * work[column][entry]
                ) % prime
    return sign * result % prime


def resultant(poly_f, poly_g, prime):
    poly_f = trim(poly_f, prime)
    poly_g = trim(poly_g, prime)
    degree_f = len(poly_f) - 1
    degree_g = len(poly_g) - 1
    require(degree_f >= 1 and degree_g >= 0, "resultant degree boundary")
    if degree_g == 0:
        return pow(poly_g[0], degree_f, prime)
    width = degree_f + degree_g
    rows = []
    descending_f = list(reversed(poly_f))
    descending_g = list(reversed(poly_g))
    for shift in range(degree_g):
        row = [0] * width
        row[shift : shift + degree_f + 1] = descending_f
        rows.append(row)
    for shift in range(degree_f):
        row = [0] * width
        row[shift : shift + degree_g + 1] = descending_g
        rows.append(row)
    return determinant_mod(rows, prime)


def is_irreducible_degree_two_or_three(poly, prime):
    degree = len(trim(poly, prime)) - 1
    require(degree in (2, 3), "root test only handles degree two or three")
    return all(
        sum(coefficient * pow(value, index, prime) for index, coefficient in enumerate(poly))
        % prime
        != 0
        for value in range(prime)
    )


def realization(prime, stratum, prescribed, order):
    require(is_irreducible_degree_two_or_three(stratum, prime), "irreducible stratum")
    stratum_derivative = derivative(stratum, prime)
    derivative_power = class_power(stratum_derivative, order - 1, stratum, prime)
    helper = class_mul(
        prescribed,
        inverse_mod(derivative_power, stratum, prime),
        stratum,
        prime,
    )
    poly_f = stratum
    poly_g = add(stratum, mul(power(stratum, order, prime), helper, prime), prime)
    require(
        remainder(derivative(poly_g, prime), stratum, prime) == stratum_derivative,
        "derivative agrees on stratum",
    )
    for lower in range(2, order):
        _, _, lower_contact = contact_numerator_denominator(
            poly_f, poly_g, stratum, lower, prime
        )
        require(lower_contact == [0], "lower contact vanishes")
    _, _, contact = contact_numerator_denominator(poly_f, poly_g, stratum, order, prime)
    require(contact == remainder(prescribed, stratum, prime), "prescribed contact realized")
    return helper, poly_f, poly_g


def multiplicative_order(poly, modulus, prime):
    group_order = prime ** (len(modulus) - 1) - 1
    value = [1]
    for order in range(1, group_order + 1):
        value = class_mul(value, poly, modulus, prime)
        if value == [1]:
            return order
    raise RuntimeError("unit order not found")


def prime_divisors(number):
    divisors = []
    divisor = 2
    while divisor * divisor <= number:
        if number % divisor == 0:
            divisors.append(divisor)
            while number % divisor == 0:
                number //= divisor
        divisor += 1
    if number > 1:
        divisors.append(number)
    return divisors


# Uniform prescribed-unit realization bank in degrees two and three.
realization_bank = (
    (5, [2, 0, 1], [1, 1], 2),
    (5, [2, 0, 1], [1, 1], 3),
    (5, [2, 0, 1], [1, 1], 4),
    (7, [1, 0, 1], [2, 3], 2),
    (7, [1, 0, 1], [2, 3], 4),
    (13, [11, 0, 1], [1, 2], 2),
    (13, [11, 0, 1], [1, 2], 3),
    (5, [1, 1, 0, 1], [1, 1], 3),
)
realization_checks = 0
for prime, stratum, prescribed, order in realization_bank:
    realization(prime, stratum, prescribed, order)
    realization_checks += 1


# The p=13 degree-two Singer contact is exactly THM-3234's multiplier.
prime = 13
singer_s = [11, 0, 1]
alpha = [1, 2]
alpha_helper, alpha_f, alpha_g = realization(prime, singer_s, alpha, 2)
require(alpha_helper == [1, 10], "explicit Singer helper H=1+10x")
alpha_matrix = multiplication_matrix(alpha, singer_s, prime)
require(alpha_matrix == [[1, 4], [2, 1]], "exact THM-3234 Singer matrix")
alpha_order = multiplicative_order(alpha, singer_s, prime)
require(alpha_order == 168, "alpha has Singer order")
group_order = prime ** 2 - 1
order_primes = prime_divisors(group_order)
require(order_primes == [2, 3, 7], "factorization of 168")

alpha_numerator, alpha_denominator, alpha_contact = contact_numerator_denominator(
    alpha_f, alpha_g, singer_s, 2, prime
)
require(alpha_contact == alpha, "physical pair contact equals alpha")
gate_residues = []
gate_resultants = []
for divisor in order_primes:
    exponent = group_order // divisor
    numerator_power = class_power(alpha_numerator, exponent, singer_s, prime)
    denominator_power = class_power(alpha_denominator, exponent, singer_s, prime)
    gate = sub(numerator_power, denominator_power, prime)
    gate_residues.append(gate)
    gate_resultant = resultant(singer_s, gate, prime)
    gate_resultants.append(gate_resultant)
    require(gate != [0] and gate_resultant != 0, "root-free Singer order gate")

alpha_norm = class_power(alpha, 14, singer_s, prime)
require(alpha_norm == [6], "Singer norm is primitive base unit 6")
require(multiplicative_order(alpha_norm, singer_s, prime) == 12, "norm order 12")


# The full truncated deformation slice H -> c_2=S'H is a 169-state affine
# carrier; alpha is compatible with this linear identification.
deformation_contacts = set()
helper_by_contact = {}
deformation_intertwining_checks = 0
singer_derivative = derivative(singer_s, prime)
for constant in range(prime):
    for linear in range(prime):
        helper = trim([constant, linear], prime)
        contact = class_mul(singer_derivative, helper, singer_s, prime)
        padded_contact = tuple(contact + [0] * (2 - len(contact)))
        deformation_contacts.add(padded_contact)
        helper_by_contact[padded_contact] = tuple(helper + [0] * (2 - len(helper)))
        shifted_helper = class_mul(alpha, helper, singer_s, prime)
        shifted_contact = class_mul(singer_derivative, shifted_helper, singer_s, prime)
        require(
            shifted_contact == class_mul(alpha, contact, singer_s, prime),
            "Singer action intertwines helper and contact slices",
        )
        deformation_intertwining_checks += 1
require(len(deformation_contacts) == 169, "contact deformation atlas is bijective")
require((0, 0) in deformation_contacts, "zero is delayed-contact class")


# Pull back THM-3234's full affine Heisenberg action through the deformation
# isomorphism.  This is a moduli reparametrization, not a physical LRC map.
heisenberg_pullback_checks = 0
for shift_r in range(prime):
    for shear in range(prime):
        for shift_w in range(prime):
            for coordinate_r, coordinate_w in deformation_contacts:
                moved_contact = (
                    (coordinate_r + shift_r) % prime,
                    (coordinate_w + shift_w - shear * coordinate_r) % prime,
                )
                require(
                    moved_contact in helper_by_contact,
                    "Heisenberg action pulls back through contact atlas",
                )
                heisenberg_pullback_checks += 1
require(heisenberg_pullback_checks == 13 ** 5, "full Heisenberg deformation census")


# Same norm and open algebra-generation gate do not imply Singer order.
beta = class_power(alpha, 49, singer_s, prime)
beta_order = multiplicative_order(beta, singer_s, prime)
beta_norm = class_power(beta, 14, singer_s, prime)
require(beta_order == 24, "beta order 24")
require(beta_norm == alpha_norm, "beta has same norm as alpha")
require(class_power(alpha, prime, singer_s, prime) != alpha, "alpha generates F169")
require(class_power(beta, prime, singer_s, prime) != beta, "beta generates F169")
beta_helper, _, _ = realization(prime, singer_s, beta, 2)


# Orbit counts under multiplication pin the carrier interpretation.
require(group_order // alpha_order == 1, "alpha is one punctured orbit")
require(group_order // beta_order == 7, "beta has seven punctured orbits")


def format_poly(poly):
    return "[" + ",".join(str(value) for value in poly) + "]"


print("dependency_thm3232_sha256=%s" % DEPENDENCIES[next(iter(DEPENDENCIES))])
print("dependency_thm3234_sha256=%s" % list(DEPENDENCIES.values())[1])
print("assert_nodes=%d,float_literals=%d" % (assert_nodes, float_literals))
print("prescribed_contact_realizations=%d" % realization_checks)
print("p13_stratum=x2-2,alpha=%s,H=%s" % (format_poly(alpha), format_poly(alpha_helper)))
print("alpha_matrix=1,4;2,1")
print("alpha_order=%d,norm=%s,norm_order=12" % (alpha_order, format_poly(alpha_norm)))
print("deformation_atlas=169,nonzero_exact_contact=168,intertwining_checks=%d" % deformation_intertwining_checks)
print("heisenberg_moduli_pullback_checks=%d" % heisenberg_pullback_checks)
print("singer_prime_gates=%s" % ",".join(str(value) for value in order_primes))
print("singer_gate_residues=%s" % ";".join(format_poly(value) for value in gate_residues))
print("singer_gate_resultants=%s" % ",".join(str(value) for value in gate_resultants))
print("beta=alpha^49=%s,H_beta=%s" % (format_poly(beta), format_poly(beta_helper)))
print("beta_order=%d,beta_norm=%s,punctured_orbits=7" % (beta_order, format_poly(beta_norm)))
print("hierarchy=unit_then_algebra-generator_then-multiplicative-Singer")
print("scope=algebraic-contact-carrier-not-lawful-LRC-owner-packet-map")
print("all_exact_checks=PASS")
