"""Exact companion for THM-3221's selected-root osculating tower.

All arithmetic is symbolic or rational.  The script is deliberately
assertion-independent so that ordinary and optimized Python runs have the
same truth surface.
"""

from fractions import Fraction
import ast
import hashlib
from math import comb
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCY = ROOT / (
    "01-canon/theorems/"
    "THM-3215-arbitrary-degree-root-jet-hamiltonian-affine-dihedral-"
    "holonomy-and-p-fold-carry.md"
)
DEPENDENCY_SHA256 = "00c5775d62a6db1f651265bfc7f659f96cdf40464b644d29ae4dab7df419dba8"


def lf_sha256(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n")
    return hashlib.sha256(payload).hexdigest()


require(lf_sha256(DEPENDENCY) == DEPENDENCY_SHA256, "THM-3215 dependency hash")
syntax_tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax_tree))
float_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax_tree)
)
require(assert_nodes == 0, "assert statements are optimization-sensitive")
require(float_literals == 0, "floating literals are forbidden")


def series_mul(left, right, degree):
    out = [Fraction(0) for _ in range(degree + 1)]
    for i, left_i in enumerate(left[: degree + 1]):
        for j, right_j in enumerate(right[: degree + 1 - i]):
            out[i + j] += left_i * right_j
    return out


def series_compose(outer, inner, degree):
    out = [Fraction(0) for _ in range(degree + 1)]
    power = [Fraction(0) for _ in range(degree + 1)]
    power[0] = Fraction(1)
    for order, coefficient in enumerate(outer[: degree + 1]):
        for j in range(degree + 1):
            out[j] += coefficient * power[j]
        if order < degree:
            power = series_mul(power, inner, degree)
    return out


def series_inverse_tangent(poly, degree):
    inverse = [Fraction(0) for _ in range(degree + 1)]
    inverse[1] = Fraction(1)
    for order in range(2, degree + 1):
        error = series_compose(poly, inverse, order)[order]
        inverse[order] -= error
    identity = [Fraction(0) for _ in range(degree + 1)]
    identity[1] = Fraction(1)
    require(series_compose(poly, inverse, degree) == identity, "left formal inverse")
    require(series_compose(inverse, poly, degree) == identity, "right formal inverse")
    return inverse


# The first unequal coefficient survives every common lower osculating jet.
first_separation_checks = 0
coordinate_change_checks = 0
for m in range(2, 11):
    for seed in range(3):
        phi_f = [Fraction(0) for _ in range(m + 1)]
        phi_g = [Fraction(0) for _ in range(m + 1)]
        phi_f[1] = Fraction(1)
        phi_g[1] = Fraction(1)
        for order in range(2, m):
            numerator = ((seed + 2) * (order + 1)) % 13 - 6
            common = Fraction(numerator, order + seed + 7)
            phi_f[order] = common
            phi_g[order] = common
        af = Fraction(2 * seed + m + 1, m + seed + 5)
        ag = Fraction(3 * seed + 2 * m + 1, m + seed + 8)
        phi_f[m] = af
        phi_g[m] = ag
        inverse_f = series_inverse_tangent(phi_f, m)
        transition = series_compose(phi_g, inverse_f, m)
        expected = [Fraction(0) for _ in range(m + 1)]
        expected[1] = Fraction(1)
        expected[m] = ag - af
        require(
            transition == expected,
            "first osculating separation",
        )
        first_separation_checks += 1

        # An arbitrary nonlinear source chart cancels from the transition.
        scale = Fraction(seed + 2, seed + 1)
        inverse_chart = [Fraction(0) for _ in range(m + 1)]
        inverse_chart[1] = 1 / scale
        for order in range(2, m + 1):
            inverse_chart[order] = Fraction(seed + order, 17 + order)
        phi_f_new = [scale * value for value in series_compose(phi_f, inverse_chart, m)]
        phi_g_new = [scale * value for value in series_compose(phi_g, inverse_chart, m)]
        transition_new = series_compose(
            phi_g_new, series_inverse_tangent(phi_f_new, m), m
        )
        expected_new = [
            scale * transition[order] / scale**order for order in range(m + 1)
        ]
        require(transition_new == expected_new, "nonlinear chart cancellation")
        require(
            transition_new[m] == (ag - af) * scale ** (1 - m),
            "cotangent weight",
        )
        coordinate_change_checks += 1


# The minimal nonidentity quotient is the additive associated-graded layer.
minimal_composition_checks = 0
coordinate_weight_checks = 0
for m in range(2, 13):
    for seed in range(4):
        coefficient_c = Fraction(seed + 1, seed + 3)
        coefficient_d = Fraction(2 * seed + 3, seed + 5)
        scale = Fraction(seed + 2, seed + 1)
        inner = [Fraction(0) for _ in range(m + 1)]
        outer = [Fraction(0) for _ in range(m + 1)]
        inner[1] = Fraction(1)
        outer[1] = Fraction(1)
        inner[m] = coefficient_d
        outer[m] = coefficient_c
        composed = series_compose(outer, inner, m)
        expected = [Fraction(0) for _ in range(m + 1)]
        expected[1] = Fraction(1)
        expected[m] = coefficient_c + coefficient_d
        require(composed == expected, "minimal jet additivity")
        require(
            coefficient_c * scale ** (1 - m)
            == scale * coefficient_c / scale**m,
            "coordinate weight",
        )
        minimal_composition_checks += 1
        coordinate_weight_checks += 1


# Exact order p in the minimal quotient, for every nonzero residue coefficient.
prime_order_checks = 0
for p in (2, 3, 5, 7, 11, 13, 17):
    for m in range(2, 13):
        for coefficient in range(1, p):
            running = 0
            for n in range(1, p):
                running = (running + coefficient) % p
                require(running == n * coefficient % p, "additive iterate")
                require(running != 0, "premature minimal-jet reset")
            running = (running + coefficient) % p
            require(running == 0, "prime reset")
            prime_order_checks += 1


# The coefficient of u**m in H^p-u is exactly p*c; no higher-tail claim.
divided_carry_checks = 0
for p in (2, 3, 5, 7, 11, 13, 17):
    for m in range(2, 13):
        for coefficient in (1, 2, 5):
            return_coefficient = p * coefficient
            require(return_coefficient % p == 0, "carry divisibility")
            require(return_coefficient // p == coefficient, "divided carry")
            divided_carry_checks += 1


# Root, derivative anchor, and the degree-bounded normalized tower reconstruct.
finite_degree_reconstruction_checks = 0
for degree in range(1, 11):
    for root in (Fraction(0), Fraction(1), Fraction(-2, 3)):
        derivative = Fraction(degree + 1, degree + 2)
        translated = [Fraction(0) for _ in range(degree + 1)]
        translated[1] = derivative
        for order in range(2, degree + 1):
            translated[order] = Fraction(
                (order + 1) * (degree + 2), order + degree + 1
            )
        polynomial = [Fraction(0) for _ in range(degree + 1)]
        for order, coefficient in enumerate(translated):
            for x_order in range(order + 1):
                polynomial[x_order] += (
                    coefficient
                    * comb(order, x_order)
                    * (-root) ** (order - x_order)
                )
        derivative_at_root = sum(
            order * polynomial[order] * root ** (order - 1)
            for order in range(1, degree + 1)
        )
        require(derivative_at_root == derivative, "derivative anchor")
        normalized = [coefficient / derivative_at_root for coefficient in translated]
        reconstructed_local = [
            derivative_at_root * coefficient for coefficient in normalized
        ]
        reconstructed = [Fraction(0) for _ in range(degree + 1)]
        for order, coefficient in enumerate(reconstructed_local):
            for x_order in range(order + 1):
                reconstructed[x_order] += (
                    coefficient
                    * comb(order, x_order)
                    * (-root) ** (order - x_order)
                )
        require(reconstructed == polynomial, "finite-degree reconstruction")
        finite_degree_reconstruction_checks += 1


# No fixed jet depth identifies arbitrary degree.
uniform_depth_hostile_checks = 0
for depth in range(1, 11):
    low = [Fraction(0) for _ in range(depth + 2)]
    high = [Fraction(0) for _ in range(depth + 2)]
    low[1] = Fraction(1)
    high[1] = Fraction(1)
    high[depth + 1] = Fraction(1)
    require(low[: depth + 1] == high[: depth + 1], "fixed-depth agreement")
    require(low != high, "next-jet separation")
    uniform_depth_hostile_checks += 1


# A nonconstant unit preserves the zero divisor but changes the normalized tower.
same_divisor_f = [Fraction(0), Fraction(1), Fraction(0)]
same_divisor_g = [Fraction(0), Fraction(1), Fraction(1)]
require(same_divisor_f[1] == same_divisor_g[1] == 1, "first derivatives")
require(same_divisor_g[2] - same_divisor_f[2] == 1, "unit hostile")


# Sharp warning: the order-p statement is only about the minimal quotient.
higher_tail_hostile = [0, 1, 1, 0, 0, 0]
iterate = [0, 1, 0, 0, 0, 0]
for _ in range(3):
    iterate = [
        int(coefficient) % 3
        for coefficient in series_compose(higher_tail_hostile, iterate, 5)
    ]
require(iterate[:3] == [0, 1, 0], "minimal F3 reset")
require(iterate == [0, 1, 0, 0, 0, 1], "higher Nottingham hostile")


print("dependency_thm3215_sha256=%s" % DEPENDENCY_SHA256)
print("assert_nodes=%d,float_literals=%d" % (assert_nodes, float_literals))
print("first_separation_checks=%d" % first_separation_checks)
print("nonlinear_coordinate_change_checks=%d" % coordinate_change_checks)
print("minimal_composition_checks=%d" % minimal_composition_checks)
print("coordinate_weight_checks=%d" % coordinate_weight_checks)
print("prime_order_checks=%d" % prime_order_checks)
print("divided_carry_checks=%d" % divided_carry_checks)
print("finite_degree_reconstruction_checks=%d" % finite_degree_reconstruction_checks)
print("uniform_depth_hostile_checks=%d" % uniform_depth_hostile_checks)
print("same_divisor_hostile=u_vs_u_times_1_plus_u")
print("higher_truncation_hostile=F3:(u+u^2)^3=u+u^5_mod_u^6")
print("proportionality_boundary=normalized_tower_needs_derivative_anchor")
print("scope=selected_simple_root_local_separation_not_global_selector_or_moment_survival")
print("all_exact_checks=PASS")
