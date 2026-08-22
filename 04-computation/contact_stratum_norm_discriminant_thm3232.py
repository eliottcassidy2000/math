"""Exact companion for THM-3232's contact-stratum norm."""

from fractions import Fraction
from math import comb
import ast
import hashlib
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCY = ROOT / (
    "01-canon/theorems/"
    "THM-3229-hasse-pluecker-simple-root-contact-gcd-flag-and-degree-termination.md"
)
DEPENDENCY_SHA256 = "316f13d5a896d185124555ce5e32bdd041cd391497517c000538a03b5ddf0aab"


def lf_sha256(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    return hashlib.sha256(payload).hexdigest()


require(lf_sha256(DEPENDENCY) == DEPENDENCY_SHA256, "THM-3229 dependency hash")
syntax_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax_tree))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax_tree)
)
require(assert_nodes == 0, "assert statements are optimization-sensitive")
require(float_literals == 0, "floating literals are forbidden")


def trim(poly):
    out = [Fraction(value) for value in poly]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def add(left, right):
    degree = max(len(left), len(right))
    out = [Fraction(0)] * degree
    for index in range(degree):
        out[index] = (left[index] if index < len(left) else 0) + (
            right[index] if index < len(right) else 0
        )
    return trim(out)


def sub(left, right):
    return add(left, [-Fraction(value) for value in right])


def scale(poly, scalar):
    return trim([Fraction(scalar) * value for value in poly])


def mul(left, right):
    out = [Fraction(0)] * (len(left) + len(right) - 1)
    for i, left_i in enumerate(left):
        for j, right_j in enumerate(right):
            out[i + j] += left_i * right_j
    return trim(out)


def power(poly, exponent):
    out = [Fraction(1)]
    for _ in range(exponent):
        out = mul(out, poly)
    return out


def monic(poly):
    poly = trim(poly)
    return scale(poly, Fraction(1, 1) / poly[-1])


def derivative(poly):
    if len(poly) <= 1:
        return [Fraction(0)]
    return trim([index * poly[index] for index in range(1, len(poly))])


def divmod_poly(dividend, divisor):
    dividend = trim(dividend)
    divisor = trim(divisor)
    require(divisor != [0], "division by zero polynomial")
    if len(dividend) < len(divisor):
        return [Fraction(0)], dividend
    quotient = [Fraction(0)] * (len(dividend) - len(divisor) + 1)
    remainder = list(dividend)
    while remainder != [0] and len(remainder) >= len(divisor):
        shift = len(remainder) - len(divisor)
        coefficient = remainder[-1] / divisor[-1]
        quotient[shift] += coefficient
        remainder = sub(
            remainder,
            [Fraction(0)] * shift + [coefficient * value for value in divisor],
        )
    return trim(quotient), trim(remainder)


def remainder(poly, modulus):
    return divmod_poly(poly, modulus)[1]


def quotient_exact(dividend, divisor):
    quotient, rem = divmod_poly(dividend, divisor)
    require(rem == [0], "expected exact polynomial division")
    return quotient


def gcd_poly(left, right):
    left = trim(left)
    right = trim(right)
    while right != [0]:
        left, right = right, remainder(left, right)
    return monic(left)


def inverse_mod(poly, modulus):
    old_r, new_r = trim(poly), trim(modulus)
    old_s, new_s = [Fraction(1)], [Fraction(0)]
    while new_r != [0]:
        quotient, next_r = divmod_poly(old_r, new_r)
        old_r, new_r = new_r, next_r
        old_s, new_s = new_s, sub(old_s, mul(quotient, new_s))
    require(len(old_r) == 1 and old_r[0] != 0, "class is not invertible")
    return remainder(scale(old_s, Fraction(1, 1) / old_r[0]), modulus)


def eval_poly(poly, point):
    value = Fraction(0)
    point = Fraction(point)
    for coefficient in reversed(poly):
        value = value * point + coefficient
    return value


def hasse_poly(poly, order):
    if order >= len(poly):
        return [Fraction(0)]
    return trim(
        [
            Fraction(poly[index]) * comb(index, order)
            for index in range(order, len(poly))
        ]
    )


def omega(poly_f, poly_g, order):
    return sub(
        mul(hasse_poly(poly_f, order), derivative(poly_g)),
        mul(hasse_poly(poly_g, order), derivative(poly_f)),
    )


def determinant(matrix):
    size = len(matrix)
    require(all(len(row) == size for row in matrix), "square matrix required")
    work = [[Fraction(value) for value in row] for row in matrix]
    value = Fraction(1)
    sign = 1
    for column in range(size):
        pivot = next(
            (row for row in range(column, size) if work[row][column] != 0),
            None,
        )
        if pivot is None:
            return Fraction(0)
        if pivot != column:
            work[column], work[pivot] = work[pivot], work[column]
            sign = -sign
        pivot_value = work[column][column]
        value *= pivot_value
        for row in range(column + 1, size):
            factor = work[row][column] / pivot_value
            for entry in range(column, size):
                work[row][entry] -= factor * work[column][entry]
    return sign * value


def resultant(poly_f, poly_g):
    poly_f = trim(poly_f)
    poly_g = trim(poly_g)
    degree_f = len(poly_f) - 1
    degree_g = len(poly_g) - 1
    require(degree_f >= 1 and degree_g >= 0, "resultant degree boundary")
    if degree_g == 0:
        return poly_g[0] ** degree_f
    width = degree_f + degree_g
    rows = []
    descending_f = list(reversed(poly_f))
    descending_g = list(reversed(poly_g))
    for shift in range(degree_g):
        row = [Fraction(0)] * width
        row[shift : shift + degree_f + 1] = descending_f
        rows.append(row)
    for shift in range(degree_f):
        row = [Fraction(0)] * width
        row[shift : shift + degree_g + 1] = descending_g
        rows.append(row)
    return determinant(rows)


def multiplication_matrix(poly, modulus):
    modulus = monic(modulus)
    degree = len(modulus) - 1
    columns = []
    for shift in range(degree):
        product = remainder(mul(poly, [Fraction(0)] * shift + [Fraction(1)]), modulus)
        columns.append(product + [Fraction(0)] * (degree - len(product)))
    return [[columns[column][row] for column in range(degree)] for row in range(degree)]


def class_norm_trace(numerator, denominator, modulus):
    inverse = inverse_mod(denominator, modulus)
    representative = remainder(mul(numerator, inverse), modulus)
    matrix = multiplication_matrix(representative, modulus)
    norm = determinant(matrix)
    trace = sum(matrix[index][index] for index in range(len(matrix)))
    return representative, norm, trace


def discriminant(monic_poly):
    monic_poly = monic(monic_poly)
    degree = len(monic_poly) - 1
    sign = -1 if degree * (degree - 1) // 2 % 2 else 1
    return sign * resultant(monic_poly, derivative(monic_poly))


def linear_factor(root):
    return [Fraction(-root), Fraction(1)]


def root_product(roots):
    out = [Fraction(1)]
    for root in roots:
        out = mul(out, linear_factor(root))
    return monic(out)


def compose_affine(poly, rho, sigma):
    affine = [Fraction(sigma), Fraction(rho)]
    out = [Fraction(0)]
    for index, coefficient in enumerate(poly):
        out = add(out, scale(power(affine, index), coefficient))
    return trim(out)


def contact_norm(poly_f, poly_g, stratum, order):
    numerator = scale(omega(poly_f, poly_g, order), -1)
    denominator = mul(derivative(poly_f), derivative(poly_g))
    representative, norm, trace = class_norm_trace(
        numerator, denominator, stratum
    )
    degree = len(stratum) - 1
    formula = ((-1) ** degree) * resultant(stratum, omega(poly_f, poly_g, order))
    formula /= resultant(stratum, denominator)
    require(norm == formula, "norm/resultant quotient and sign")
    require(norm != 0, "first-contact norm is nonzero")
    return representative, norm, trace


# Split rational bank with two roots at each of contact depths two and three.
split_roots = (Fraction(-3), Fraction(-1), Fraction(2), Fraction(4))
split_contacts = (2, 2, 3, 3)
split_p = root_product(split_roots)
split_q = [Fraction(1)]
for root, order in zip(split_roots, split_contacts):
    split_q = mul(split_q, power(linear_factor(root), order - 1))
split_f = split_p
split_g = mul(split_p, add([Fraction(1)], split_q))
current = split_p
split_stratum_checks = 0
split_branch_checks = 0
for order in (2, 3):
    next_flag = gcd_poly(current, omega(split_f, split_g, order))
    stratum = monic(quotient_exact(current, next_flag))
    expected_roots = [
        root
        for root, contact_order in zip(split_roots, split_contacts)
        if contact_order == order
    ]
    require(stratum == root_product(expected_roots), "exact split contact stratum")
    _, norm, _ = contact_norm(split_f, split_g, stratum, order)
    branch_product = Fraction(1)
    for root in expected_roots:
        first_f = eval_poly(derivative(split_f), root)
        first_g = eval_poly(derivative(split_g), root)
        branch_value = -eval_poly(omega(split_f, split_g, order), root) / (
            first_f * first_g
        )
        require(branch_value != 0, "split branch coefficient")
        branch_product *= branch_value
        split_branch_checks += 1
    require(norm == branch_product, "norm equals unordered branch product")
    split_stratum_checks += 1
    current = next_flag
require(current == [1], "split contact flag terminates")


# Irreducible diagonal-power bank: g=P+P^m and c_m=(P')^(m-1).
irreducible_p = [Fraction(-2), Fraction(0), Fraction(1)]
irreducible_resultant = resultant(irreducible_p, derivative(irreducible_p))
require(irreducible_resultant == -8, "quadratic derivative resultant sign")
power_norms = []
power_family_checks = 0
for order in range(2, 7):
    poly_f = irreducible_p
    poly_g = add(irreducible_p, power(irreducible_p, order))
    for lower_order in range(2, order):
        require(
            remainder(omega(poly_f, poly_g, lower_order), irreducible_p) == [0],
            "lower Hasse--Pluecker vanishing in power family",
        )
    representative, norm, trace = contact_norm(
        poly_f, poly_g, irreducible_p, order
    )
    expected_class = remainder(
        power(derivative(irreducible_p), order - 1), irreducible_p
    )
    require(representative == expected_class, "power-family contact class")
    require(norm == irreducible_resultant ** (order - 1), "discriminant power")
    if order == 2:
        require(trace == 0 and norm == -8, "trace-cancellation hostile")
    power_norms.append(norm)
    power_family_checks += 1


# Affine covariance and its discriminant-balanced scalar invariant.
cubic_p = [Fraction(0), Fraction(-1), Fraction(0), Fraction(1)]
cubic_order = 4
cubic_f = cubic_p
cubic_g = add(cubic_p, power(cubic_p, cubic_order))
_, cubic_norm, _ = contact_norm(cubic_f, cubic_g, cubic_p, cubic_order)
cubic_disc = discriminant(cubic_p)
require(cubic_disc == 4 and cubic_norm == -64, "cubic sign controls")
rho = Fraction(3)
sigma = Fraction(2)
affine_f = compose_affine(cubic_f, rho, sigma)
affine_g = compose_affine(cubic_g, rho, sigma)
affine_stratum = scale(compose_affine(cubic_p, rho, sigma), rho ** -3)
_, affine_norm, _ = contact_norm(
    affine_f, affine_g, affine_stratum, cubic_order
)
affine_disc = discriminant(affine_stratum)
require(
    affine_norm == rho ** (3 * (cubic_order - 1)) * cubic_norm,
    "affine contact-norm weight",
)
require(
    affine_disc == rho ** (-3 * 2) * cubic_disc,
    "affine discriminant weight",
)
cubic_balanced = cubic_norm ** 2 * cubic_disc ** 3
affine_balanced = affine_norm ** 2 * affine_disc ** 3
require(cubic_balanced == affine_balanced != 0, "balanced affine invariant")


# Constant row rescaling, swapping, and a genuine GL2 target-frame law.
scaled_f = scale(cubic_f, 2)
scaled_g = scale(cubic_g, -3)
_, scaled_norm, _ = contact_norm(scaled_f, scaled_g, cubic_p, cubic_order)
require(scaled_norm == cubic_norm, "constant row rescaling invariance")
_, swapped_norm, _ = contact_norm(cubic_g, cubic_f, cubic_p, cubic_order)
require(swapped_norm == (-1) ** 3 * cubic_norm, "row-swap norm sign")
swapped_balanced = swapped_norm ** 2 * cubic_disc ** 3
require(swapped_balanced == cubic_balanced, "balanced row-swap invariance")

frame_f = add(scale(cubic_f, 2), cubic_g)
frame_g = add(scale(cubic_f, -1), scale(cubic_g, 3))
frame_det = Fraction(7)
_, frame_norm, _ = contact_norm(frame_f, frame_g, cubic_p, cubic_order)
expected_frame_ratio = Fraction(7, 6) ** 3
require(frame_norm == expected_frame_ratio * cubic_norm, "GL2 normalized frame law")
old_wedge_norm = ((-1) ** 3) * resultant(
    cubic_p, omega(cubic_f, cubic_g, cubic_order)
)
new_wedge_norm = ((-1) ** 3) * resultant(
    cubic_p, omega(frame_f, frame_g, cubic_order)
)
require(new_wedge_norm == frame_det ** 3 * old_wedge_norm, "Pluecker norm law")


# Finite-field irreducible Frobenius norm: P=x^2+2 over F5 and c=P'=2x.
def trim_mod(poly, prime):
    out = [int(value) % prime for value in poly]
    while len(out) > 1 and out[-1] == 0:
        out.pop()
    return out


def reduce_mod(poly, modulus, prime):
    out = trim_mod(poly, prime)
    modulus = trim_mod(modulus, prime)
    while len(out) >= len(modulus) and out != [0]:
        shift = len(out) - len(modulus)
        factor = out[-1] * pow(modulus[-1], -1, prime) % prime
        for index, coefficient in enumerate(modulus):
            out[index + shift] = (out[index + shift] - factor * coefficient) % prime
        out = trim_mod(out, prime)
    return out


def mul_mod(left, right, modulus, prime):
    product = [0] * (len(left) + len(right) - 1)
    for i, left_i in enumerate(left):
        for j, right_j in enumerate(right):
            product[i + j] = (product[i + j] + left_i * right_j) % prime
    return reduce_mod(product, modulus, prime)


def power_mod(poly, exponent, modulus, prime):
    out = [1]
    base = reduce_mod(poly, modulus, prime)
    while exponent:
        if exponent % 2:
            out = mul_mod(out, base, modulus, prime)
        base = mul_mod(base, base, modulus, prime)
        exponent //= 2
    return out


finite_prime = 5
finite_p = [2, 0, 1]
finite_c = [0, 2]
finite_norm = int(resultant(finite_p, finite_c)) % finite_prime
finite_frobenius = power_mod(finite_c, 1 + finite_prime, finite_p, finite_prime)
require(finite_norm == 3, "finite-field resultant norm")
require(finite_frobenius == [finite_norm], "finite-field Frobenius orbit norm")


# Exact boundaries: a coarse pre-stratum norm vanishes, and a repeated divisor
# has zero discriminant, so neither is an etale first-contact norm.
boundary_roots = (Fraction(-1), Fraction(1))
boundary_contacts = (2, 3)
boundary_p = root_product(boundary_roots)
boundary_q = mul(linear_factor(-1), power(linear_factor(1), 2))
boundary_f = boundary_p
boundary_g = mul(boundary_p, add([Fraction(1)], boundary_q))
coarse_omega = omega(boundary_f, boundary_g, 2)
require(resultant(boundary_p, coarse_omega) == 0, "coarse mixed-depth norm vanishes")
boundary_g2 = gcd_poly(boundary_p, coarse_omega)
boundary_s2 = monic(quotient_exact(boundary_p, boundary_g2))
_, boundary_norm, _ = contact_norm(boundary_f, boundary_g, boundary_s2, 2)
require(boundary_norm != 0 and len(boundary_s2) == 2, "exact stratum repairs norm")
repeated = power(linear_factor(1), 2)
require(discriminant(repeated) == 0, "non-etale repeated-root boundary")


print("dependency_thm3229_sha256=%s" % DEPENDENCY_SHA256)
print("assert_nodes=%d,float_literals=%d" % (assert_nodes, float_literals))
print("split_contact_strata=2,3")
print("split_stratum_checks=%d,split_branch_checks=%d" % (split_stratum_checks, split_branch_checks))
print("irreducible_power_family_orders=2..6")
print("irreducible_power_norms=%s" % ",".join(str(value) for value in power_norms))
print("trace_hostile=Q[x]/(x2-2):trace(2x)=0,norm(2x)=-8")
print("affine_covariance=r3,m4,rho3:weight9")
print("balanced_invariant=%s" % cubic_balanced)
print("target_frame=diagonal-invariant,swap-balanced,GL2-weight-checked")
print("finite_field_frobenius=F5:x2+2:norm3")
print("boundaries=coarse-mixed-depth-zero,repeated-discriminant-zero")
print("scope=root-free-orbit-product-not-root-selector-or-additive-carry")
print("all_exact_checks=PASS")
