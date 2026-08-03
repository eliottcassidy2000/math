#!/usr/bin/env python3
"""Exact arbitrary-degree root-jet holonomy controls for THM-3215."""

from __future__ import annotations

import ast
import hashlib
from itertools import product
from pathlib import Path

import sympy as sp


ROOT = Path(__file__).resolve().parents[1]
DEPENDENCY = (
    ROOT
    / "01-canon"
    / "theorems"
    / "THM-3206-heterogeneous-factorial-exterior-reflection-groupoid-and-fixed-plane-holonomy.md"
)
EXPECTED_DEPENDENCY_SHA256 = (
    "6b355bcf8c0098824fb51f3e2019bdfce355924ba5c680266b2a0907b49f672a"
)


def require(condition: bool, message: object) -> None:
    if not condition:
        raise RuntimeError(message)


def lf_sha256(path: Path) -> str:
    data = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(data).hexdigest()


def zero_matrix(matrix: sp.Matrix) -> bool:
    return matrix.applyfunc(sp.factor) == sp.zeros(*matrix.shape)


dependency_hash = lf_sha256(DEPENDENCY)
require(
    dependency_hash == EXPECTED_DEPENDENCY_SHA256,
    ("dependency hash", dependency_hash),
)


# The universal Hamiltonian two-jet and its scalar square.
f_value, f_prime, f_second = sp.symbols("f_value f_prime f_second")
jet = sp.Matrix(
    [[f_prime, f_second], [-2 * f_value, -f_prime]]
)
jet_discriminant = f_prime**2 - 2 * f_value * f_second
require(
    zero_matrix(jet * jet - jet_discriminant * sp.eye(2)),
    "universal jet square",
)


# THM-3206's factorial reflection is exactly this jet in a root-adapted frame.
x, s, v = sp.symbols("x s v")
q = v * x**2 - x + s
q_prime = sp.diff(q, x)
q_second = sp.diff(q, x, 2)
factorial_reflection = sp.Matrix(
    [
        [-(1 + 2 * v), 2 * v],
        [-2 * (s + v + 1), 2 * v + 1],
    ]
)
root_frame = sp.Matrix([[1, 0], [1 + x, 1]])
quadratic_jet = sp.Matrix(
    [[q_prime, q_second], [-2 * q, -q_prime]]
)
require(
    zero_matrix(root_frame.inv() * factorial_reflection * root_frame - quadratic_jet),
    "factorial quadratic jet conjugacy",
)
require(
    sp.factor(q_prime**2 - 2 * q * q_second - (1 - 4 * s * v)) == 0,
    "quadratic discriminant is constant",
)


# Arbitrary-degree controls: only the value and first two derivatives enter.
u = sp.symbols("u")
arbitrary_degree_checks = 0
root_polynomials: list[sp.Expr] = []
for degree in range(1, 9):
    for seed in range(1, 7):
        coefficients = [
            sp.Integer(((-1) ** (seed + index)) * (seed + 2 * index + 1))
            for index in range(1, degree)
        ]
        polynomial = sp.Integer(seed + 1) * u
        polynomial += sum(
            coefficient * u ** (index + 1)
            for index, coefficient in enumerate(coefficients, start=1)
        )
        value = polynomial.subs(u, 0)
        first = sp.diff(polynomial, u).subs(u, 0)
        second = sp.diff(polynomial, u, 2).subs(u, 0)
        local_jet = sp.Matrix([[first, second], [-2 * value, -first]])
        require(value == 0 and first != 0, ("simple root", degree, seed))
        require(
            local_jet * local_jet == first**2 * sp.eye(2),
            ("arbitrary-degree reflection", degree, seed),
        )
        root_polynomials.append(polynomial)
        arbitrary_degree_checks += 1


# Affine-dihedral word law at a common simple root.
def reflection(parameter: sp.Expr) -> sp.Matrix:
    return sp.Matrix([[1, parameter], [0, -1]])


def translation(parameter: sp.Expr) -> sp.Matrix:
    return sp.Matrix([[1, parameter], [0, 1]])


t1, t2, t3, shift = sp.symbols("t1 t2 t3 shift")
require(reflection(t1) ** 2 == sp.eye(2), "reflection square")
require(
    reflection(t2) * reflection(t1) == translation(t1 - t2),
    "two-reflection translation",
)
require(
    translation(t2 - t3) * translation(t1 - t2)
    == translation(t1 - t3),
    "strict transition composition",
)
require(
    reflection(t1) * translation(shift) * reflection(t1)
    == translation(-shift),
    "dihedral conjugation",
)

MAX_SYMBOLIC_WORD = 8
word_parameters = sp.symbols("kappa1:" + str(MAX_SYMBOLIC_WORD + 1))
for length in range(1, MAX_SYMBOLIC_WORD + 1):
    actual = sp.eye(2)
    for parameter in word_parameters[:length]:
        actual = reflection(parameter) * actual
    theta = sum(
        (-1) ** (index + 1) * parameter
        for index, parameter in enumerate(word_parameters[:length], start=1)
    )
    expected = sp.Matrix([[1, theta], [0, (-1) ** length]])
    require(zero_matrix(actual - expected), ("word law", length))


# The pair translation is the next Pluecker jet.
lambda_1, lambda_2, second_1, second_2 = sp.symbols(
    "lambda_1 lambda_2 second_1 second_2", nonzero=True
)
root_jet_1 = sp.Matrix([[lambda_1, second_1], [0, -lambda_1]])
root_jet_2 = sp.Matrix([[lambda_2, second_2], [0, -lambda_2]])
mu = lambda_1 * lambda_2
omega_2 = second_1 * lambda_2 - second_2 * lambda_1
require(
    root_jet_2 * root_jet_1 == mu * sp.eye(2) + omega_2 * sp.Matrix([[0, 1], [0, 0]]),
    "second-jet Pluecker shear",
)


# Coordinate changes give a common affine change of every logarithmic jet.
coordinate_scale, coordinate_second = sp.symbols(
    "coordinate_scale coordinate_second", nonzero=True
)
kappa_1, kappa_2, kappa_3 = sp.symbols("kappa_1 kappa_2 kappa_3")


def transformed_kappa(kappa: sp.Expr) -> sp.Expr:
    return kappa / coordinate_scale - coordinate_second / coordinate_scale**2


require(
    sp.factor(
        (transformed_kappa(kappa_1) - transformed_kappa(kappa_2))
        - (kappa_1 - kappa_2) / coordinate_scale
    )
    == 0,
    "coordinate change of pair one-form",
)
require(
    sp.factor(
        (
            transformed_kappa(kappa_1)
            - transformed_kappa(kappa_2)
            + transformed_kappa(kappa_3)
        )
        - (
            (kappa_1 - kappa_2 + kappa_3) / coordinate_scale
            - coordinate_second / coordinate_scale**2
        )
    )
    == 0,
    "odd word retains connection shift",
)


# The same pair one-form is the quadratic coefficient of the transition germ.
f_first, g_first, f_second_symbol, g_second_symbol = sp.symbols(
    "f_first g_first f_second_symbol g_second_symbol", nonzero=True
)
transition_first = g_first / f_first
transition_second = (
    g_second_symbol * f_first - g_first * f_second_symbol
) / f_first**3
pair_theta = f_second_symbol / f_first - g_second_symbol / g_first
require(
    sp.factor(pair_theta + f_first * transition_second / transition_first) == 0,
    "transition germ quadratic coefficient",
)


# The factorial common-root pair has the explicit rank-one nilpotent of THM-3206.
v_1, v_2 = sp.symbols("v_1 v_2")
s_1 = x - v_1 * x**2
s_2 = x - v_2 * x**2


def factorial_matrix(s_value: sp.Expr, v_value: sp.Expr) -> sp.Matrix:
    return sp.Matrix(
        [
            [-(1 + 2 * v_value), 2 * v_value],
            [-2 * (s_value + v_value + 1), 2 * v_value + 1],
        ]
    )


lambda_q1 = 2 * v_1 * x - 1
lambda_q2 = 2 * v_2 * x - 1
pair_product = factorial_matrix(s_2, v_2) * factorial_matrix(s_1, v_1)
nilpotent = sp.simplify(pair_product - lambda_q1 * lambda_q2 * sp.eye(2))
expected_nilpotent = 2 * (v_1 - v_2) * sp.Matrix(
    [[x + 1, -1], [(x + 1) ** 2, -(x + 1)]]
)
require(zero_matrix(nilpotent - expected_nilpotent), "factorial pair nilpotent")
require(zero_matrix(nilpotent * nilpotent), "factorial nilpotent square")


# Exact binomial carry for a rank-one unipotent return.
carry_exponent = sp.symbols("carry_exponent", integer=True, positive=True)
nilpotent_unit = sp.Matrix([[0, 1], [0, 0]])
for exponent in range(1, 13):
    left = (mu * sp.eye(2) + omega_2 * nilpotent_unit) ** exponent
    right = (
        mu**exponent * sp.eye(2)
        + exponent * mu ** (exponent - 1) * omega_2 * nilpotent_unit
    )
    require(zero_matrix(left - right), ("nilpotent binomial", exponent))


# Pure integer modular helpers for exhaustive finite-field controls.
Matrix2 = tuple[tuple[int, int], tuple[int, int]]
IDENTITY: Matrix2 = ((1, 0), (0, 1))


def matrix_multiply(left: Matrix2, right: Matrix2, prime: int) -> Matrix2:
    return tuple(
        tuple(
            sum(left[row][index] * right[index][column] for index in range(2))
            % prime
            for column in range(2)
        )
        for row in range(2)
    )  # type: ignore[return-value]


def matrix_power(matrix: Matrix2, exponent: int, prime: int) -> Matrix2:
    answer = IDENTITY
    factor = matrix
    remaining = exponent
    while remaining:
        if remaining & 1:
            answer = matrix_multiply(answer, factor, prime)
        factor = matrix_multiply(factor, factor, prime)
        remaining //= 2
    return answer


def scalar_matrix(matrix: Matrix2, prime: int) -> bool:
    return (
        matrix[0][1] % prime == 0
        and matrix[1][0] % prime == 0
        and matrix[0][0] % prime == matrix[1][1] % prime
    )


def factorial_matrix_mod(s_value: int, v_value: int, prime: int) -> Matrix2:
    return (
        (-(1 + 2 * v_value) % prime, 2 * v_value % prime),
        (-2 * (s_value + v_value + 1) % prime, (2 * v_value + 1) % prime),
    )


coordinate_change_checks = 0
for prime in (3, 5, 7):
    for scale in range(1, prime):
        inverse_scale = pow(scale, -1, prime)
        for second in range(prime):
            for left in range(prime):
                for right in range(prime):
                    left_new = (
                        left * inverse_scale
                        - second * inverse_scale * inverse_scale
                    ) % prime
                    right_new = (
                        right * inverse_scale
                        - second * inverse_scale * inverse_scale
                    ) % prime
                    require(
                        (left_new - right_new) % prime
                        == (left - right) * inverse_scale % prime,
                        ("finite coordinate law", prime, scale, second, left, right),
                    )
                    coordinate_change_checks += 1


factorial_pair_checks = 0
nontrivial_pair_order_checks = 0
factorial_word_checks = 0
for prime in (3, 5, 7, 11, 13):
    for root in range(prime):
        bank: list[tuple[int, int, int, int]] = []
        for v_value in range(prime):
            s_value = (root - v_value * root * root) % prime
            derivative = (2 * v_value * root - 1) % prime
            if s_value != 0 and derivative != 0:
                kappa = (2 * v_value * pow(derivative, -1, prime)) % prime
                bank.append((v_value, s_value, derivative, kappa))
        for left in bank:
            for right in bank:
                v_left, s_left, lambda_left, kappa_left = left
                v_right, s_right, lambda_right, kappa_right = right
                product_matrix = matrix_multiply(
                    factorial_matrix_mod(s_right, v_right, prime),
                    factorial_matrix_mod(s_left, v_left, prime),
                    prime,
                )
                require(
                    scalar_matrix(product_matrix, prime)
                    == (kappa_left == kappa_right)
                    == (v_left == v_right),
                    ("factorial pair classification", prime, root, left, right),
                )
                if v_left != v_right:
                    require(
                        scalar_matrix(matrix_power(product_matrix, prime, prime), prime),
                        ("p-fold scalar return", prime, root, left, right),
                    )
                    for exponent in range(1, prime):
                        require(
                            not scalar_matrix(
                                matrix_power(product_matrix, exponent, prime), prime
                            ),
                            ("exact projective order", prime, exponent, left, right),
                        )
                    nontrivial_pair_order_checks += 1
                factorial_pair_checks += 1

        if prime <= 7:
            for length in (3, 4, 5):
                for word in product(bank, repeat=length):
                    product_matrix = IDENTITY
                    theta = 0
                    for index, (v_value, s_value, _, kappa) in enumerate(word, start=1):
                        product_matrix = matrix_multiply(
                            factorial_matrix_mod(s_value, v_value, prime),
                            product_matrix,
                            prime,
                        )
                        theta += (1 if index % 2 else -1) * kappa
                    theta %= prime
                    require(
                        scalar_matrix(product_matrix, prime)
                        == (length % 2 == 0 and theta == 0),
                        ("factorial word classification", prime, root, length, word),
                    )
                    factorial_word_checks += 1


# Heterogeneous four-reflection scalar return: pairwise sharpness does not extend.
hostile_prime = 7
hostile_root = 1
hostile_v = (0, 2, 3, 5)
hostile_rows: list[tuple[int, int, int, int]] = []
hostile_product = IDENTITY
hostile_theta = 0
for index, v_value in enumerate(hostile_v, start=1):
    s_value = (hostile_root - v_value * hostile_root**2) % hostile_prime
    derivative = (2 * v_value * hostile_root - 1) % hostile_prime
    kappa = 2 * v_value * pow(derivative, -1, hostile_prime) % hostile_prime
    hostile_rows.append((s_value, v_value, derivative, kappa))
    hostile_product = matrix_multiply(
        factorial_matrix_mod(s_value, v_value, hostile_prime),
        hostile_product,
        hostile_prime,
    )
    hostile_theta += (1 if index % 2 else -1) * kappa
require(hostile_theta % hostile_prime == 0, "heterogeneous hostile theta")
require(scalar_matrix(hostile_product, hostile_prime), "heterogeneous scalar hostile")
require(len(set(hostile_v)) == 4, "heterogeneous hostile distinctness")


# Exact characteristic-zero first-carry controls for p>=5.
padic_carry_checks = 0
integer_q1 = factorial_matrix(1, 0)
integer_q2 = factorial_matrix(-1, 2)
integer_pair = integer_q2 * integer_q1
integer_mu = sp.Integer(-3)
integer_omega = sp.Integer(4)
integer_pair_root_frame = sp.Matrix([[1, 0], [2, 1]])
require(
    integer_pair_root_frame.inv() * integer_pair * integer_pair_root_frame
    == integer_mu * sp.eye(2) + integer_omega * nilpotent_unit,
    "integer carry root frame",
)
for prime in (5, 7, 11, 13, 17):
    difference = integer_pair**prime - integer_mu**prime * sp.eye(2)
    require(
        all(int(entry) % prime == 0 for entry in difference),
        ("carry divisibility", prime),
    )
    first_carry = difference.applyfunc(lambda entry: int(entry) // prime)
    expected_carry = (
        integer_mu ** (prime - 1)
        * integer_omega
        * integer_pair_root_frame
        * nilpotent_unit
        * integer_pair_root_frame.inv()
    )
    require(
        all(
            int(first_carry[row, column] - expected_carry[row, column]) % prime == 0
            for row in range(2)
            for column in range(2)
        ),
        ("first carry residue", prime),
    )
    require(
        any(int(entry) % prime != 0 for entry in first_carry),
        ("primitive first carry", prime),
    )
    padic_carry_checks += 1


# Boundaries: higher jets, multiple roots, and the separate infinity chart.
higher_jet_f = u
higher_jet_g = u + u**3
higher_f_matrix = sp.Matrix(
    [
        [sp.diff(higher_jet_f, u).subs(u, 0), sp.diff(higher_jet_f, u, 2).subs(u, 0)],
        [0, -sp.diff(higher_jet_f, u).subs(u, 0)],
    ]
)
higher_g_matrix = sp.Matrix(
    [
        [sp.diff(higher_jet_g, u).subs(u, 0), sp.diff(higher_jet_g, u, 2).subs(u, 0)],
        [0, -sp.diff(higher_jet_g, u).subs(u, 0)],
    ]
)
require(higher_jet_f != higher_jet_g, "higher-jet hostile distinctness")
require(higher_f_matrix == higher_g_matrix, "two-jet invisibility hostile")
require(higher_g_matrix * higher_f_matrix == sp.eye(2), "two-jet scalar hostile")

multiple_root = u**2
multiple_jet = sp.Matrix(
    [
        [sp.diff(multiple_root, u).subs(u, 0), sp.diff(multiple_root, u, 2).subs(u, 0)],
        [-2 * multiple_root.subs(u, 0), -sp.diff(multiple_root, u).subs(u, 0)],
    ]
)
require(multiple_jet.rank() == 1, "multiple-root rank")
require(multiple_jet**2 == sp.zeros(2), "multiple-root nilpotence")

s_infinity_1, s_infinity_2 = sp.symbols("s_infinity_1 s_infinity_2")
infinity_1 = factorial_matrix(s_infinity_1, 0)
infinity_2 = factorial_matrix(s_infinity_2, 0)
require(
    infinity_2 * infinity_1
    == sp.eye(2) + sp.Matrix([[0, 0], [2 * (s_infinity_2 - s_infinity_1), 0]]),
    "quadratic infinity chart",
)


source = Path(__file__).read_text(encoding="utf-8")
require(
    not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
    "assert-sensitive test",
)
require(
    not any(expression.has(sp.Float) for expression in (
        tuple(jet) + tuple(quadratic_jet) + tuple(nilpotent)
        + tuple(integer_pair) + tuple(multiple_jet)
    )),
    "floating atom entered exact evidence",
)


print("dependency_thm3206_sha256=" + dependency_hash)
print("universal_two_jet_square=J(f)^2=(fprime^2-2*f*fsecond)I")
print("factorial_quadratic_conjugacy=P_x^-1*F(s,v)*P_x=J_x(q)")
print("arbitrary_degree_simple_root_checks=" + str(arbitrary_degree_checks))
print("symbolic_affine_dihedral_words=1.." + str(MAX_SYMBOLIC_WORD))
print("flat_transition_cocycle=T_jk*T_ij=T_ik")
print("coordinate_change_checks=" + str(coordinate_change_checks))
print("pair_shear=Omega2=fsecond*gprime-gsecond*fprime")
print("factorial_common_root_pair_checks=" + str(factorial_pair_checks))
print("factorial_nontrivial_order_p_checks=" + str(nontrivial_pair_order_checks))
print("factorial_common_root_word_checks=" + str(factorial_word_checks))
print("heterogeneous_scalar_four_word_rows=" + repr(hostile_rows))
print("heterogeneous_scalar_four_word_matrix=" + repr(hostile_product))
print("padic_first_carry_checks=" + str(padic_carry_checks))
print("higher_jet_hostile=(u,u+u^3):same_two_jet_scalar_return")
print("multiple_root_hostile=u^2:rank1_nilpotent_not_reflection")
print("infinity_boundary=separate_homogeneous_chart_required")
print("sympy_float_atoms=0")
print("scope=common_selected_simple_root_local_jet_not_global_root_selector_or_GMC2")
print("all_exact_checks=PASS")
