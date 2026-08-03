"""Exact companion for THM-3229's Hasse--Pluecker contact gcd flag."""

from fractions import Fraction
from math import comb, factorial
import ast
import hashlib
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCY = ROOT / (
    "01-canon/theorems/"
    "THM-3221-selected-root-osculating-separation-and-minimal-jet-prime-carry.md"
)
DEPENDENCY_SHA256 = "6ae707482cf73ebbd995722fe056d093df37d3aea0af6f58d6989e18d0b06f84"


def lf_sha256(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    return hashlib.sha256(payload).hexdigest()


require(lf_sha256(DEPENDENCY) == DEPENDENCY_SHA256, "THM-3221 dependency hash")
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
    for i in range(degree):
        out[i] = (left[i] if i < len(left) else 0) + (
            right[i] if i < len(right) else 0
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
        subtractor = [Fraction(0)] * shift + [coefficient * value for value in divisor]
        remainder = sub(remainder, subtractor)
    return trim(quotient), trim(remainder)


def gcd_poly(left, right):
    left = trim(left)
    right = trim(right)
    while right != [0]:
        _, remainder = divmod_poly(left, right)
        left, right = right, remainder
    return monic(left)


def eval_poly(poly, point):
    point = Fraction(point)
    value = Fraction(0)
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


def hasse_eval(poly, order, point):
    return eval_poly(hasse_poly(poly, order), point)


def omega(poly_f, poly_g, order):
    return sub(
        mul(hasse_poly(poly_f, order), hasse_poly(poly_g, 1)),
        mul(hasse_poly(poly_g, order), hasse_poly(poly_f, 1)),
    )


def linear_factor(root):
    return [Fraction(-root), Fraction(1)]


def root_product(roots):
    out = [Fraction(1)]
    for root in roots:
        out = mul(out, linear_factor(root))
    return monic(out)


# A four-root bank with prescribed contacts 2,3,4,5.
roots = (Fraction(-2), Fraction(0), Fraction(1), Fraction(3))
contacts = (2, 3, 4, 5)
poly_p = root_product(roots)
poly_q = [Fraction(1)]
for root, contact_order in zip(roots, contacts):
    poly_q = mul(poly_q, power(linear_factor(root), contact_order - 1))
poly_f = poly_p
poly_g = mul(poly_p, add([Fraction(1)], poly_q))
degree_f = len(poly_f) - 1
degree_g = len(poly_g) - 1
degree_bound = max(degree_f, degree_g)


orientation_checks = 0
prescribed_contact_checks = 0
for root, contact_order in zip(roots, contacts):
    require(eval_poly(poly_f, root) == 0, "f common root")
    require(eval_poly(poly_g, root) == 0, "g common root")
    first_f = hasse_eval(poly_f, 1, root)
    first_g = hasse_eval(poly_g, 1, root)
    require(first_f != 0 and first_g != 0, "simple common root")
    for order in range(2, degree_bound + 1):
        normalized_difference = (
            hasse_eval(poly_g, order, root) / first_g
            - hasse_eval(poly_f, order, root) / first_f
        )
        omega_value = eval_poly(omega(poly_f, poly_g, order), root)
        require(
            normalized_difference == -omega_value / (first_f * first_g),
            "Hasse--Pluecker orientation",
        )
        orientation_checks += 1
        if order < contact_order:
            require(normalized_difference == 0, "prescribed lower contact")
            prescribed_contact_checks += 1
        elif order == contact_order:
            require(normalized_difference != 0, "prescribed first separation")
            prescribed_contact_checks += 1


# The squarefree simple-root gcd flag removes exactly the contact-m stratum.
gcd_flag_checks = 0
current = gcd_poly(poly_f, poly_g)
require(current == poly_p, "initial simple common-root gcd")
for order in range(2, degree_bound + 1):
    current = gcd_poly(current, omega(poly_f, poly_g, order))
    surviving_roots = [
        root for root, contact_order in zip(roots, contacts) if contact_order > order
    ]
    expected = root_product(surviving_roots)
    require(current == expected, "nested Hasse contact gcd flag")
    gcd_flag_checks += 1
require(current == [1], "degree-bound termination")


# Polynomial degree and constant-target scaling laws.
degree_checks = 0
scaling_checks = 0
for order in range(2, degree_bound + 1):
    omega_order = omega(poly_f, poly_g, order)
    if omega_order != [0]:
        require(
            len(omega_order) - 1 <= degree_f + degree_g - order - 1,
            "Omega degree bound",
        )
    scaled = omega(scale(poly_f, 2), scale(poly_g, -3), order)
    require(scaled == scale(omega_order, -6), "constant target scaling")
    degree_checks += 1
    scaling_checks += 1


# Proportional pairs never terminate; this is the exact equality boundary.
proportional_checks = 0
for order in range(2, degree_f + 1):
    require(omega(poly_f, scale(poly_f, 7), order) == [0], "proportional wall")
    proportional_checks += 1


# Hasse derivatives, unlike ordinary derivatives, retain small-characteristic jets.
small_f = [0, 1]
small_g = [0, 1, 0, 1]
p = 3
first_f_mod = int(hasse_eval(small_f, 1, 0)) % p
first_g_mod = int(hasse_eval(small_g, 1, 0)) % p
omega_2_mod = int(eval_poly(omega(small_f, small_g, 2), 0)) % p
omega_3_mod = int(eval_poly(omega(small_f, small_g, 3), 0)) % p
ordinary_third_mod = factorial(3) % p
require(first_f_mod == first_g_mod == 1, "F3 simple root")
require(omega_2_mod == 0 and omega_3_mod == 2, "F3 Hasse contact three")
require((-omega_3_mod) % p == 1, "F3 transition orientation")
require(ordinary_third_mod == 0, "ordinary derivative loses F3 third jet")


# A multiple common root is outside the saturated simple-root flag.
multiple_f = [0, 0, 1]
multiple_g = mul(multiple_f, [1, 1])
require(hasse_eval(multiple_f, 1, 0) == 0, "multiple f derivative")
require(hasse_eval(multiple_g, 1, 0) == 0, "multiple g derivative")


print("dependency_thm3221_sha256=%s" % DEPENDENCY_SHA256)
print("assert_nodes=%d,float_literals=%d" % (assert_nodes, float_literals))
print("prescribed_contact_bank=2,3,4,5")
print("hasse_pluecker_orientation_checks=%d" % orientation_checks)
print("prescribed_contact_checks=%d" % prescribed_contact_checks)
print("nested_gcd_flag_checks=%d" % gcd_flag_checks)
print("degree_bound_termination=D%d" % degree_bound)
print("omega_degree_checks=%d" % degree_checks)
print("constant_scaling_checks=%d" % scaling_checks)
print("proportional_boundary_checks=%d" % proportional_checks)
print("small_characteristic_hostile=F3:x_vs_x+x3_contact3_Hasse_only")
print("multiple_root_boundary=first_derivatives_zero")
print("scope=root_free_simple-divisor-flag_not_root_ordering_or_resonant_moment_PRS")
print("all_exact_checks=PASS")
