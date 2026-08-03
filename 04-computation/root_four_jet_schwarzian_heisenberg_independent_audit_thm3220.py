"""Independent raw-series audit of THM-3220.

This companion does not import the primary verifier and does not take its
log-coordinate implementation as an oracle.  It composes and inverts
truncated series directly, then reconstructs every claimed coordinate law.
"""

from fractions import Fraction
from itertools import product
import ast
import hashlib
from pathlib import Path

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
PRIMARY_SCRIPT = ROOT / (
    "04-computation/root_four_jet_schwarzian_heisenberg_thm3220.py"
)
PRIMARY_OUTPUT = ROOT / (
    "05-knowledge/results/root_four_jet_schwarzian_heisenberg_thm3220.out"
)
PRIMARY_SCRIPT_SHA256 = (
    "f399b6e86b639e9b27859dcbc9d2fcc6b44f4fb24aae51706e205a39e38c9e98"
)
PRIMARY_OUTPUT_SHA256 = (
    "1255b5ceebdf2d6b26ee92d9cbc3901e686728ba7ee79fdf6dbe7be52934e743"
)


def lf_sha256(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    return hashlib.sha256(payload).hexdigest()


require(lf_sha256(PRIMARY_SCRIPT) == PRIMARY_SCRIPT_SHA256, "primary script hash")
require(lf_sha256(PRIMARY_OUTPUT) == PRIMARY_OUTPUT_SHA256, "primary output hash")
syntax_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax_tree))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax_tree)
)
require(assert_nodes == 0, "optimization-sensitive assert statement")
require(float_literals == 0, "floating literal in exact audit")


def series_mul(left, right, order):
    result = [0] * order
    for i, left_coefficient in enumerate(left[:order]):
        for j, right_coefficient in enumerate(right[: order - i]):
            result[i + j] += left_coefficient * right_coefficient
    return tuple(result)


def series_compose(outer, inner, order):
    result = [0] * order
    power = [0] * order
    power[0] = 1
    for outer_coefficient in outer[:order]:
        for degree in range(order):
            result[degree] += outer_coefficient * power[degree]
        power = list(series_mul(tuple(power), inner, order))
    return tuple(result)


def series_inverse(series, order):
    require(series[0] == 0 and series[1] != 0, "invertible pointed series")
    inverse = [0] * order
    inverse[1] = sp.cancel(sp.Integer(1) / series[1])
    for degree in range(2, order):
        inverse[degree] = 0
        residual = series_compose(series, tuple(inverse), order)[degree]
        inverse[degree] = sp.cancel(-residual / series[1])
    candidate = tuple(inverse)
    identity = (0, 1) + (0,) * (order - 2)
    require(
        all(
            sp.simplify(left - right) == 0
            for left, right in zip(series_compose(series, candidate, order), identity)
        ),
        "left formal inverse",
    )
    require(
        all(
            sp.simplify(left - right) == 0
            for left, right in zip(series_compose(candidate, series, order), identity)
        ),
        "right formal inverse",
    )
    return candidate


def log_coordinates(series):
    alpha, beta, gamma = series[2:5]
    return (
        alpha,
        sp.expand(beta - alpha**2),
        sp.expand(gamma - sp.Rational(5, 2) * alpha * beta + sp.Rational(3, 2) * alpha**3),
    )


def raw_from_log(log_vector):
    first, second, central = log_vector
    return (
        0,
        1,
        first,
        sp.expand(second + first**2),
        sp.expand(central + sp.Rational(5, 2) * first * second + first**3),
    )


def transition(target, source):
    return series_compose(target, series_inverse(source, 5), 5)


def same_series(left, right):
    return all(sp.simplify(a - b) == 0 for a, b in zip(left, right))


# Derive the Heisenberg law from raw substitution rather than invoking it.
ai, bi, ci, ao, bo, co = sp.symbols("ai bi ci ao bo co")
inner = (0, 1, ai, bi, ci)
outer = (0, 1, ao, bo, co)
raw_composite = series_compose(outer, inner, 5)
inner_log = log_coordinates(inner)
outer_log = log_coordinates(outer)
actual_log = log_coordinates(raw_composite)
expected_log = (
    inner_log[0] + outer_log[0],
    inner_log[1] + outer_log[1],
    inner_log[2]
    + outer_log[2]
    + sp.Rational(1, 2)
    * (inner_log[0] * outer_log[1] - inner_log[1] * outer_log[0]),
)
require(
    all(sp.simplify(actual - expected) == 0 for actual, expected in zip(actual_log, expected_log)),
    "raw derivation of Heisenberg law",
)
raw_inverse = series_inverse(inner, 5)
require(
    all(
        sp.simplify(actual + expected) == 0
        for actual, expected in zip(log_coordinates(raw_inverse), inner_log)
    ),
    "raw derivation of logarithmic inverse",
)


# Recover the Schwarzian dictionary from independent derivative symbols.
f1, f2, f3, f4 = sp.symbols("f1 f2 f3 f4", nonzero=True)
normalized_root_jet = (
    0,
    1,
    f2 / (2 * f1),
    f3 / (6 * f1),
    f4 / (24 * f1),
)
root_log = log_coordinates(normalized_root_jet)
kappa = f2 / f1
schwarzian = f3 / f1 - sp.Rational(3, 2) * kappa**2
schwarzian_prime = f4 / f1 - 4 * f2 * f3 / f1**2 + 3 * f2**3 / f1**3
require(sp.simplify(root_log[0] - kappa / 2) == 0, "root curvature coordinate")
require(sp.simplify(root_log[1] - schwarzian / 6) == 0, "root Schwarzian coordinate")
require(
    sp.simplify(root_log[2] - (schwarzian_prime - kappa * schwarzian) / 24) == 0,
    "root covariant fourth coordinate",
)


# Direct raw transition triangles reconstruct the half-area invoice.
vertex_logs = (
    (Fraction(0), Fraction(0), Fraction(0)),
    (Fraction(1), Fraction(0), Fraction(0)),
    (Fraction(0), Fraction(2), Fraction(0)),
    (Fraction(2), Fraction(-1), Fraction(3, 2)),
    (Fraction(-1), Fraction(4), Fraction(-3)),
    (Fraction(3, 2), Fraction(5, 2), Fraction(7, 2)),
)
vertices = tuple(raw_from_log(vector) for vector in vertex_logs)
triangle_checks = 0
for source in vertices:
    for middle in vertices:
        for target in vertices:
            first_edge = transition(middle, source)
            second_edge = transition(target, middle)
            direct_edge = transition(target, source)
            require(
                same_series(series_compose(second_edge, first_edge, 5), direct_edge),
                "strict raw transition triangle",
            )
            first_log = log_coordinates(first_edge)
            second_log = log_coordinates(second_edge)
            direct_log = log_coordinates(direct_edge)
            area = (
                first_log[0] * second_log[1] - first_log[1] * second_log[0]
            ) / 2
            require(
                sp.simplify(direct_log[2] - first_log[2] - second_log[2] - area)
                == 0,
                "central half-area transgression",
            )
            triangle_checks += 1

hostile_10 = log_coordinates(transition(vertices[1], vertices[0]))
hostile_21 = log_coordinates(transition(vertices[2], vertices[1]))
hostile_20 = log_coordinates(transition(vertices[2], vertices[0]))
require(hostile_10 == (1, 0, 0), "area hostile first edge")
require(hostile_21 == (-1, 2, -1), "area hostile second edge")
require(hostile_20 == (0, 2, 0), "area hostile direct edge")
require(hostile_10[2] + hostile_21[2] != hostile_20[2], "naive central sum fails")


# A genuinely nonlinear common source change cancels from every transition.
def output_scale(series, scalar):
    return tuple(scalar * coefficient for coefficient in series)


coordinate_change_checks = 0
source_changes = (
    (Fraction(2), Fraction(3, 5), Fraction(-1, 3), Fraction(2, 7)),
    (Fraction(3), Fraction(-2, 5), Fraction(4, 9), Fraction(1, 6)),
)
for rho, quadratic, cubic, quartic in source_changes:
    source_change = (0, 1 / rho, quadratic, cubic, quartic)
    input_rescale = (0, 1 / rho, 0, 0, 0)
    for source in vertices:
        for target in vertices:
            source_y = output_scale(series_compose(source, source_change, 5), rho)
            target_y = output_scale(series_compose(target, source_change, 5), rho)
            transition_y = transition(target_y, source_y)
            transition_x = transition(target, source)
            conjugated = output_scale(
                series_compose(transition_x, input_rescale, 5), rho
            )
            require(same_series(transition_y, conjugated), "nonlinear coordinate cancellation")
            old_log = log_coordinates(transition_x)
            new_log = log_coordinates(transition_y)
            require(
                all(
                    sp.simplify(new_log[index] - old_log[index] / rho ** (index + 1))
                    == 0
                    for index in range(3)
                ),
                "transition tensor weights one two three",
            )
            coordinate_change_checks += 1


# Characteristic two is checked from raw polynomial composition only.
def series_compose_mod(outer, inner, order, prime):
    raw = series_compose(outer, inner, order)
    return tuple(int(coefficient) % prime for coefficient in raw)


f2_elements = tuple((0, 1, aa, bb, cc) for aa, bb, cc in product(range(2), repeat=3))
f2_identity = (0, 1, 0, 0, 0)


def finite_power(element, exponent, prime):
    running = f2_identity
    for _ in range(exponent):
        running = series_compose_mod(element, running, 5, prime)
    return running


f2_orders = []
for element in f2_elements:
    for exponent in range(1, 9):
        if finite_power(element, exponent, 2) == f2_identity:
            f2_orders.append(exponent)
            break
require(sorted(f2_orders) == [1, 2, 2, 2, 2, 2, 4, 4], "D8 order census")
f2_center = tuple(
    element
    for element in f2_elements
    if all(
        series_compose_mod(element, other, 5, 2)
        == series_compose_mod(other, element, 5, 2)
        for other in f2_elements
    )
)
require(
    f2_center == ((0, 1, 0, 0, 0), (0, 1, 0, 0, 1)),
    "D8 center from raw composition",
)
for element in f2_elements:
    aa, bb = element[2], element[3]
    require(
        finite_power(element, 2, 2) == (0, 1, 0, 0, aa * (aa + bb) % 2),
        "D8 square refinement",
    )


# Quadratic root germs directly trace the twisted cubic and its discriminant.
r, s = sp.symbols("r s")
quadratic_r = (0, 1, r, 0, 0)
quadratic_s = (0, 1, s, 0, 0)
require(
    same_series(raw_from_log((r, -r**2, sp.Rational(3, 2) * r**3)), quadratic_r),
    "quadratic twisted-cubic logarithm",
)
commutator = series_inverse(quadratic_s, 5)
commutator = series_compose(series_inverse(quadratic_r, 5), commutator, 5)
commutator = series_compose(quadratic_s, commutator, 5)
commutator = series_compose(quadratic_r, commutator, 5)
require(
    same_series(commutator, (0, 1, 0, 0, r * s * (s - r))),
    "raw quadratic central commutator",
)
T = sp.symbols("T")
require(
    sp.factor(sp.discriminant(T * (T - r) * (T - s), T) - (r * s * (s - r)) ** 2)
    == 0,
    "two-germ oriented discriminant",
)
r0, r1, r2 = sp.symbols("r0 r1 r2")
quadratic_vertices = tuple((0, 1, value, 0, 0) for value in (r0, r1, r2))
edge_10 = log_coordinates(transition(quadratic_vertices[1], quadratic_vertices[0]))
edge_21 = log_coordinates(transition(quadratic_vertices[2], quadratic_vertices[1]))
oriented_area = sp.factor(edge_10[0] * edge_21[1] - edge_10[1] * edge_21[0])
vandermonde = -(r1 - r0) * (r2 - r1) * (r2 - r0)
require(sp.factor(oriented_area - vandermonde) == 0, "three-germ Vandermonde area")
cubic = (T - r0) * (T - r1) * (T - r2)
require(
    sp.factor(sp.discriminant(cubic, T) - oriented_area**2) == 0,
    "three-germ discriminant square",
)


# Raw formal commutators recover the positive-Witt leading coefficient.
witt_checks = 0
for first_degree in range(2, 11):
    for second_degree in range(2, 11):
        if first_degree == second_degree:
            continue
        order = first_degree + second_degree
        first = [Fraction(0)] * order
        second = [Fraction(0)] * order
        first[1] = Fraction(1)
        second[1] = Fraction(1)
        first_coefficient = Fraction(2)
        second_coefficient = Fraction(-3)
        first[first_degree] = first_coefficient
        second[second_degree] = second_coefficient
        bracket = series_inverse(tuple(second), order)
        bracket = series_compose(series_inverse(tuple(first), order), bracket, order)
        bracket = series_compose(tuple(second), bracket, order)
        bracket = series_compose(tuple(first), bracket, order)
        leading_degree = first_degree + second_degree - 1
        identity = (Fraction(0), Fraction(1)) + (Fraction(0),) * (order - 2)
        require(
            bracket[:leading_degree] == identity[:leading_degree],
            "positive-Witt lower-degree cancellation",
        )
        require(
            bracket[leading_degree]
            == (first_degree - second_degree) * first_coefficient * second_coefficient,
            "positive-Witt leading coefficient and sign",
        )
        witt_checks += 1


# Work in F_p[delta]/(delta^2-Delta) to test the deck/Frobenius character.
def quadratic_mul(left, right, delta_square, prime):
    return (
        (left[0] * right[0] + left[1] * right[1] * delta_square) % prime,
        (left[0] * right[1] + left[1] * right[0]) % prime,
    )


def quadratic_power(element, exponent, delta_square, prime):
    running = (1, 0)
    base = element
    remaining = exponent
    while remaining:
        if remaining % 2:
            running = quadratic_mul(running, base, delta_square, prime)
        base = quadratic_mul(base, base, delta_square, prime)
        remaining //= 2
    return running


frobenius_checks = 0
deck_checks = 0
for prime in (3, 5, 7, 11, 13):
    for delta_square in range(1, prime):
        delta = (0, 1)
        inverse_delta = (0, pow(delta_square, -1, prime))
        character = pow(delta_square, (prime - 1) // 2, prime)
        for velocity in range(1, prime):
            inverse_cube = quadratic_power(inverse_delta, 3, delta_square, prime)
            scalar = 2 * velocity**3 % prime
            gamma = (scalar * inverse_cube[0] % prime, scalar * inverse_cube[1] % prime)
            gamma_square = quadratic_mul(gamma, gamma, delta_square, prime)
            expected_norm = 4 * velocity**6 * pow(delta_square, -3, prime) % prime
            require(gamma_square == (expected_norm, 0), "quadratic deck norm")
            deck_gamma = (gamma[0], -gamma[1] % prime)
            require(deck_gamma == (-gamma[0] % prime, -gamma[1] % prime), "deck sign")
            gamma_frobenius = quadratic_power(gamma, prime, delta_square, prime)
            require(
                gamma_frobenius
                == (character * gamma[0] % prime, character * gamma[1] % prime),
                "quadratic Frobenius eigencharacter",
            )
            frobenius_checks += 1
            deck_checks += 1


def fraction_mod(value, prime):
    rational = Fraction(value)
    require(rational.denominator % prime != 0, "p-integral denominator")
    return rational.numerator * pow(rational.denominator, -1, prime) % prime


def raw_from_log_mod(log_vector, prime):
    first, second, central = log_vector
    inverse_two = pow(2, -1, prime)
    return (
        0,
        1,
        first % prime,
        (second + first**2) % prime,
        (central + 5 * inverse_two * first * second + first**3) % prime,
    )


# Exhaustive small-prime order checks use only raw coefficient composition.
order_p_checks = 0
for prime in (3, 5, 7):
    identity = (0, 1, 0, 0, 0)
    for log_vector in product(range(prime), repeat=3):
        if log_vector == (0, 0, 0):
            continue
        generator = raw_from_log_mod(log_vector, prime)
        running = identity
        for exponent in range(1, prime + 1):
            running = series_compose_mod(generator, running, 5, prime)
            if exponent < prime:
                require(running != identity, "premature odd-prime reset")
            else:
                require(running == identity, "odd-prime order-p reset")
        order_p_checks += 1


# Independent raw iteration recovers the full first logarithmic carry.
carry_checks = 0
carry_vectors = (
    (Fraction(1), Fraction(2), Fraction(3)),
    (Fraction(2), Fraction(3), Fraction(5)),
    (Fraction(-3), Fraction(4), Fraction(2)),
    (Fraction(5), Fraction(-2), Fraction(7)),
)
for prime in (3, 5, 7, 11, 13, 17):
    for log_vector in carry_vectors:
        generator = tuple(Fraction(value) for value in raw_from_log(log_vector))
        running = (Fraction(0), Fraction(1), Fraction(0), Fraction(0), Fraction(0))
        for _ in range(prime):
            running = series_compose(generator, running, 5)
        for coefficient, expected in zip(running[2:5], log_vector):
            require(fraction_mod(coefficient, prime) == 0, "p-fold coefficient divisibility")
            require(
                fraction_mod(coefficient / prime, prime) == fraction_mod(expected, prime),
                "primitive logarithmic first carry",
            )
        carry_checks += 1


# Both stated failure boundaries are replayed from raw germs.
u_four = (0, 1, 0, 0, 0)
u_five = (0, 1, 0, 0, 0, 1)
require(u_four == u_five[:5] and u_five[5] == 1, "four-jet recovery hostile")
t = sp.symbols("t")
unit = (1, t, 0, 0, 0)
base_f = (0, 1, 0, 0, 0)
base_g = (0, 1, 1, 0, 0)
scaled_f = series_mul(unit, base_f, 5)
scaled_g = series_mul(unit, base_g, 5)
old_transition = log_coordinates(transition(base_g, base_f))
new_transition = log_coordinates(transition(scaled_g, scaled_f))
require(old_transition[1] == -1, "unscaled common-unit control")
require(sp.simplify(new_transition[1] - (-1 - t)) == 0, "nonconstant-unit hostile")


print("primary_script_sha256=%s" % PRIMARY_SCRIPT_SHA256)
print("primary_output_sha256=%s" % PRIMARY_OUTPUT_SHA256)
print("assert_nodes=%d,float_literals=%d" % (assert_nodes, float_literals))
print("independent_engine=raw_truncated_series_composition_and_inversion")
print("heisenberg_law=derived_symbolically_from_raw_coefficients")
print("selected_root_dictionary=kappa_over_2,Schwarzian_over_6,covariant_fourth_over_24")
print("strict_raw_transition_triangles=%d" % triangle_checks)
print("nonlinear_coordinate_change_checks=%d" % coordinate_change_checks)
print("area_hostile=naive_central_sum_fails_by_one")
print("characteristic_two_boundary=D8_with_quadratic_square_refinement")
print("quadratic_curve=twisted_cubic_with_oriented_discriminant_center")
print("positive_witt_raw_commutators=%d" % witt_checks)
print("deck_sign_checks=%d" % deck_checks)
print("quadratic_extension_frobenius_checks=%d" % frobenius_checks)
print("raw_nonzero_order_p_checks=%d" % order_p_checks)
print("raw_full_first_carry_checks=%d" % carry_checks)
print("hostiles=four_jet_nonrecovery_and_nonconstant_common_unit")
print("scope=no_global_root_selector_no_physical_carrier")
print("independent_exact_audit=PASS")
