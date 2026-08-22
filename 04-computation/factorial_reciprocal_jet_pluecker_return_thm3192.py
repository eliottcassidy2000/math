#!/usr/bin/env python3
"""Exact reciprocal-jet and PRS chart controls for THM-3192."""

from hashlib import sha256
from math import comb, factorial
from pathlib import Path

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def lf_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEPENDENCIES = (
    (HERE / "factorial_quadratic_recurrence_census_thm3124.py",
     "575fd6a8ae2bbedb4a6991a66382d2772033ea666424a1b9b57244f01deb1b7c"),
    (ROOT / "05-knowledge/results/factorial_quadratic_recurrence_census_thm3124.out",
     "2383aa94f5ae03d5a4009fe596b3e4cc675b557020c77828e15ece45b885e05a"),
    (HERE / "factorial_six_step_prime_resonance_third_euclidean_newton_thm3176.py",
     "a988112c15153b1925fe5c4fb18b7132fdd1eed998247850fc4b691fbc0c51f7"),
    (ROOT / "05-knowledge/results/factorial_six_step_prime_resonance_third_euclidean_newton_thm3176.out",
     "71b73fdcf2247d23ef602d1879b55206846687a03c1891efa560f3beb975cd66"),
    (HERE / "factorial_exterior_path_convolution_thm3186.py",
     "5bacf9a0d7da1e467f9f22fc0d21d0c0b1968dab955736440e9010fb00e21eff"),
    (ROOT / "05-knowledge/results/factorial_exterior_path_convolution_thm3186.out",
     "66d600435bdb58b0fd0b3c82c49c7f9b2154b51e624a12e8b6bb3424ccac5776"),
)
for dependency, expected_hash in DEPENDENCIES:
    require(lf_hash(dependency) == expected_hash,
            ("dependency hash drift", str(dependency)))


def compound(matrix):
    pairs = ((0, 1), (0, 2), (1, 2))
    return sp.Matrix([
        [sp.det(matrix.extract(rows, columns)) for columns in pairs]
        for rows in pairs
    ])


# Exact reciprocal gauge.
n = sp.symbols("n", integer=True, positive=True)
d, z = sp.symbols("d z", nonzero=True)
A_n = 2 * (n + 1) * (2 * n + 1)
B_n = n * (n + 1) * z * (z - 4 * d)
C_n = (d - n - 1) * z
J_n = sp.Matrix([
    [A_n, B_n, C_n],
    [1, 0, 0],
    [0, 0, d * z],
])

# Entrywise degree shifts from
# D_(n+1) S_n(1/z) D_n^(-1), where
# D_n=diag(z^n,z^(n-1),z^n).
S_reciprocal = sp.Matrix([
    [2 * (n + 1) * (2 * n + 1) / z,
     n * (n + 1) * (1 - 4 * d / z),
     d - n - 1],
    [1, 0, 0],
    [0, 0, d],
])
degree_gauge = sp.Matrix([
    [z, z**2, z],
    [1, z, 1],
    [z, z**2, z],
])
gauged = S_reciprocal.multiply_elementwise(degree_gauge)
require((gauged - J_n).applyfunc(sp.cancel) == sp.zeros(3),
        "reciprocal gauge identity")


# z-adic Smith and exterior layers over Q(d,n)[[z]].
determinant = sp.factor(sp.det(J_n))
require(determinant == -d * n * (n + 1) * z**2 * (z - 4 * d),
        "reciprocal determinant")
unit_minor = sp.factor(J_n.extract((0, 2), (0, 2)).det())
require(sp.factor(unit_minor - A_n * d * z) == 0,
        "valuation-one minor")
exterior = compound(J_n)
expected_exterior = sp.Matrix([
    [-B_n, -C_n, 0],
    [0, A_n * d * z, B_n * d * z],
    [0, d * z, 0],
])
require((exterior - expected_exterior).applyfunc(sp.expand) == sp.zeros(3),
        "reciprocal exterior reconstruction")
normalized_exterior = exterior.applyfunc(lambda entry: sp.cancel(entry / z))
residue = normalized_exterior.subs(z, 0)
expected_residue = sp.Matrix([
    [4 * d * n * (n + 1), n + 1 - d, 0],
    [0, A_n * d, 0],
    [0, d, 0],
])
require((residue - expected_residue).applyfunc(sp.factor) == sp.zeros(3)
        and residue.rank() == 2,
        "normalized exterior residue")
right_kernel = sp.Matrix([0, 0, 1])
left_conormal = sp.Matrix([[0, 1, -A_n]])
require(residue * right_kernel == sp.zeros(3, 1),
        "residue right kernel")
require((left_conormal * residue).applyfunc(sp.factor) == sp.zeros(1, 3),
        "residue left conormal")


def first_layer(index):
    scalar_a = 2 * (index + 1) * (2 * index + 1)
    return sp.Matrix([
        [4 * d * index * (index + 1), index + 1 - d, 0],
        [0, scalar_a * d, 0],
        [0, d, 0],
    ])


def inclusion(index):
    scalar_a = 2 * (index + 1) * (2 * index + 1)
    return sp.Matrix([[1, 0], [0, 1], [0, 1 / scalar_a]])


projection = sp.Matrix([[1, 0, 0], [0, 1, 0]])


def internal_flag(index):
    scalar_a = 2 * (index + 1) * (2 * index + 1)
    return sp.Matrix([
        [4 * d * index * (index + 1), index + 1 - d],
        [0, scalar_a * d],
    ])


require((projection * inclusion(n) - sp.eye(2)).applyfunc(sp.factor)
        == sp.zeros(2), "flag splitting")
require((first_layer(n) - inclusion(n) * internal_flag(n) * projection)
        .applyfunc(sp.factor) == sp.zeros(3), "first-layer flag factorization")
FIRST_LAYER_WORD_CHECKS = 0
FIRST_LAYER_CONVOLUTION_CHECKS = 0
for length in range(1, 6):
    direct_word = sp.eye(3)
    internal_word = sp.eye(2)
    for shift in range(length):
        direct_word = first_layer(n + shift) * direct_word
        internal_word = internal_flag(n + shift) * internal_word
    expected_word = inclusion(n + length - 1) * internal_word * projection
    require((direct_word - expected_word).applyfunc(sp.factor) == sp.zeros(3),
            ("first-layer word factorization", length))
    FIRST_LAYER_WORD_CHECKS += 1
    gate_u = sum(
        sp.prod(4 * (n + shift) * (n + shift + 1)
                for shift in range(j_shift + 1, length))
        * (n + j_shift + 1)
        * sp.prod(2 * (n + shift + 1) * (2 * (n + shift) + 1)
                  for shift in range(j_shift))
        for j_shift in range(length)
    )
    gate_v = sum(
        sp.prod(4 * (n + shift) * (n + shift + 1)
                for shift in range(j_shift + 1, length))
        * sp.prod(2 * (n + shift + 1) * (2 * (n + shift) + 1)
                  for shift in range(j_shift))
        for j_shift in range(length)
    )
    require(sp.factor(internal_word[0, 1]
                      - d**(length - 1) * (gate_u - d * gate_v)) == 0,
            ("first-layer single gate", length))
    FIRST_LAYER_CONVOLUTION_CHECKS += 1
require(sp.factor(internal_flag(n)[0, 0] / internal_flag(n)[1, 1]
                  - sp.Rational(2) * n / (2 * n + 1)) == 0,
        "first-layer diagonal ratio")

second_layer_return = (exterior * right_kernel).applyfunc(
    lambda entry: sp.cancel(entry / z**2)
).subs(z, 0)
require((second_layer_return
         - sp.Matrix([0, -4 * d**2 * n * (n + 1), 0]))
        .applyfunc(sp.factor) == sp.zeros(3, 1),
        "hidden column second-layer return")


# Direct moment polynomials, reciprocal coefficients, and truncation.
v = sp.symbols("v")


def moment_polynomial(order, scalar_d):
    return sp.Poly(sum(
        sp.Integer(comb(order, index)) * v**index * sum(
            sp.Integer(comb(order - index, ell))
            * scalar_d**(order - index - ell)
            * (-1)**ell
            * sp.Integer(factorial(2 * index + ell))
            for ell in range(order - index + 1)
        )
        for index in range(order + 1)
    ), v)


moments = [moment_polynomial(order, d) for order in range(9)]
reciprocals = [sp.Poly(sp.expand(z**order
                                * moments[order].as_expr().subs(v, 1 / z)), z)
               for order in range(9)]
JET_COEFFICIENT_CHECKS = 0
TRUNCATION_CHECKS = 0
for order in range(9):
    for jet_index in range(order + 1):
        require(reciprocals[order].nth(jet_index)
                == moments[order].nth(order - jet_index),
                ("reciprocal coefficient identity", order, jet_index))
        JET_COEFFICIENT_CHECKS += 1
for order in range(1, 8):
    recurrence_residual = sp.expand(
        reciprocals[order + 1].as_expr()
        - 2 * (order + 1) * (2 * order + 1)
          * reciprocals[order].as_expr()
        - order * (order + 1) * z * (z - 4 * d)
          * reciprocals[order - 1].as_expr()
        - (d - order - 1) * d**order * z**(order + 1)
    )
    require(recurrence_residual == 0, ("reciprocal recurrence", order))
    for height in range(order + 1):
        modulus = z**(height + 1)
        boundary_remainder = sp.rem(
            sp.Poly((d - order - 1) * d**order * z**(order + 1), z),
            sp.Poly(modulus, z),
        )
        homogeneous_remainder = sp.rem(
            sp.Poly(recurrence_residual
                    + (d - order - 1) * d**order * z**(order + 1), z),
            sp.Poly(modulus, z),
        )
        require(boundary_remainder.is_zero and homogeneous_remainder.is_zero,
                ("jet truncation", order, height))
        TRUNCATION_CHECKS += 1


# Exact offset-six reciprocal-jet projection to the maintained PRS charts.
p = sp.symbols("p", integer=True, positive=True)
offset = p + 6


def factorial_ratio_from_2p(argument):
    """Return argument!/(2p)! for a fixed integral argument-2p."""

    shift = sp.expand(argument - 2 * p)
    require(shift.is_Integer, "nonintegral factorial shift")
    shift = int(shift)
    if shift >= 0:
        return sp.prod(2 * p + index for index in range(1, shift + 1))
    return sp.cancel(1 / sp.prod(2 * p - index
                                for index in range(0, -shift)))


def normalized_top_coefficient(degree_shift, coefficient_shift):
    """[v^(p+coefficient_shift)]A_(p+degree_shift)/(2p)! exactly."""

    codimension = degree_shift - coefficient_shift
    require(codimension >= 0, "negative top codimension")
    choose = (sp.prod(p + degree_shift - index
                      for index in range(codimension))
              * sp.Rational(1, factorial(codimension)))
    total = sum(
        sp.binomial(codimension, ell)
        * offset**(codimension - ell)
        * (-1)**ell
        * factorial_ratio_from_2p(2 * (p + coefficient_shift) + ell)
        for ell in range(codimension + 1)
    )
    return sp.factor(choose * total)


A = {shift: normalized_top_coefficient(4, shift)
     for shift in range(-1, 5)}
B = {shift: normalized_top_coefficient(5, shift)
     for shift in range(0, 6)}
quotient_lead = (2 * p + 10) * (2 * p + 9)
R = {
    shift: sp.factor(
        (2 * p + 7) * B[shift]
        - (2 * p + 7) * quotient_lead * A[shift - 1]
        + 2 * (p + 5) * (p + 6) * A[shift]
    )
    for shift in range(4)
}
H = 24 * p + 109
D_2 = (p + 5) * (2 * p + 3) * H**2
N_1 = -2 * (2 * p + 5) * (2 * p + 7) * (2 * p + 3) * H
N_0 = 4 * (p + 6) * (2 * p + 5) * (28 * p + 123)
S = {
    shift: sp.factor(D_2 * A[shift]
                     - N_1 * R[shift - 1]
                     - N_0 * R[shift])
    for shift in (1, 2)
}
J_wall = (256 * p**4 - 27648 * p**3 - 365600 * p**2
          - 1528800 * p - 2096649)
K_wall = (5120 * p**5 - 810240 * p**4 - 14447872 * p**3
          - 92004672 * p**2 - 256323456 * p - 265142241)


def toeplitz_chart(first, second):
    f0, f1, f2 = first
    g0, g1, g2 = second
    matrix = sp.Matrix([[f0, 0, g0],
                        [f1, f0, g1],
                        [f2, f1, g2]])
    return sp.factor(matrix.det())


def wedge_chart(first, second):
    f0, f1 = first
    g0, g1 = second
    return sp.factor(f0 * g1 - f1 * g0)


P_H = toeplitz_chart((A[4], A[3], A[2]),
                     (B[5], B[4], B[3]))
P_J = toeplitz_chart((R[3], R[2], R[1]),
                     (A[4], A[3], A[2]))
P_K = wedge_chart((S[2], S[1]), (R[3], R[2]))
U_H = (256 * sp.prod((p + index)**2 for index in range(1, 5))
       * (2 * p + 1)**2 * (2 * p + 3)**2 * (2 * p + 5)**2
       * (2 * p + 7))
U_J = (64 * sp.prod((p + index)**2 for index in range(1, 5))
       * (p + 5) * (2 * p + 1)**2 * (2 * p + 3))
U_K = (64 * sp.prod((p + index)**2 for index in range(1, 6))
       * (p + 6) * (2 * p + 1) * (2 * p + 7) / (2 * p - 1))
require(sp.factor(P_H - U_H * R[3]) == 0, "H Toeplitz chart")
require(sp.factor(P_J - U_J * S[2]) == 0, "J Toeplitz chart")
require(sp.factor(P_K - U_K * K_wall) == 0, "K wedge chart")

expected_R_lead = (-8 * sp.prod(p + index for index in range(1, 6))
                   * (2 * p + 1) * (2 * p + 3) * H)
expected_S_lead = (4 * sp.prod(p + index for index in range(1, 6))
                   * (2 * p + 7) * J_wall)
require(sp.factor(R[3] - expected_R_lead) == 0, "H leading row")
require(sp.factor(S[2] - expected_S_lead) == 0, "J leading row")


def unit_for_every_prime_at_least_197(expression):
    numerator, denominator = sp.fraction(sp.cancel(expression))
    constants = (int(numerator.subs(p, 0)),
                 int(denominator.subs(p, 0)))
    require(all(value != 0 for value in constants),
            "candidate unit has zero residue polynomial")
    return all(all(prime_factor < 197
                   for prime_factor in sp.factorint(abs(value)))
               for value in constants)


UNIT_EXPRESSIONS = (U_H, U_J, U_K,
                    sp.factor(R[3] / H), sp.factor(S[2] / J_wall))
require(all(unit_for_every_prime_at_least_197(expression)
            for expression in UNIT_EXPRESSIONS),
        "offset-six p-unit range")
D_3 = (2 * p - 1) * (2 * p + 7) * J_wall**2
P_1 = (-2 * (2 * p + 1) * (2 * p + 3) * H
       * (2 * p - 1) * J_wall)
P_0 = 4 * (p + 6) * (2 * p + 1) * K_wall
require(sp.factor(D_3 * R[3] - P_1 * S[2]) == 0,
        "third quotient leading cancellation")
require(sp.factor(D_3 * R[2] - P_1 * S[1] - P_0 * S[2]) == 0,
        "third quotient connection cancellation")


# Direct integer controls for the rational top-jet formulas.
DIRECT_PRS_COEFFICIENT_CHECKS = 0
for prime in (5, 7, 11):
    scalar_d = prime + 6
    denominator = factorial(2 * prime)
    direct_A = moment_polynomial(prime + 4, scalar_d)
    direct_B = moment_polynomial(prime + 5, scalar_d)
    direct_a = {shift: sp.Rational(direct_A.nth(prime + shift), denominator)
                for shift in A}
    direct_b = {shift: sp.Rational(direct_B.nth(prime + shift), denominator)
                for shift in B}
    direct_r = {
        shift: ((2 * prime + 7) * direct_b[shift]
                - (2 * prime + 7) * ((2 * prime + 10) * (2 * prime + 9))
                  * direct_a[shift - 1]
                + 2 * (prime + 5) * (prime + 6) * direct_a[shift])
        for shift in R
    }
    h_value = 24 * prime + 109
    d2_value = (prime + 5) * (2 * prime + 3) * h_value**2
    n1_value = (-2 * (2 * prime + 5) * (2 * prime + 7)
                * (2 * prime + 3) * h_value)
    n0_value = (4 * (prime + 6) * (2 * prime + 5)
                * (28 * prime + 123))
    direct_s = {
        shift: (d2_value * direct_a[shift]
                - n1_value * direct_r[shift - 1]
                - n0_value * direct_r[shift])
        for shift in S
    }
    for symbolic, direct, label in ((A, direct_a, "A"),
                                    (B, direct_b, "B"),
                                    (R, direct_r, "R"),
                                    (S, direct_s, "S")):
        for shift, expression in symbolic.items():
            require(direct[shift] == expression.subs(p, prime),
                    ("direct PRS top jet", label, prime, shift))
            DIRECT_PRS_COEFFICIENT_CHECKS += 1


# Two complementary chart hostiles.
def reciprocal_transfer(index, scalar_d):
    return sp.Matrix([
        [2 * (index + 1) * (2 * index + 1),
         index * (index + 1) * z * (z - 4 * scalar_d),
         (scalar_d - index - 1) * z],
        [1, 0, 0],
        [0, 0, scalar_d * z],
    ])


hostile_product = sp.eye(3)
for index in (1, 2, 3):
    hostile_product = compound(reciprocal_transfer(index, 5)) * hostile_product
hostile_vector = (hostile_product * right_kernel).applyfunc(sp.factor)
expected_hostile = sp.Matrix([
    60 * z**4 * (z - 20) * (4 * z - 105),
    3000 * z**4 * (z - 20) * (z**2 - 20 * z + 140),
    7500 * z**4 * (z - 20),
])
require(hostile_vector == expected_hostile, "reciprocal cancellation hostile")
hostile_point = sp.Rational(105, 4)
hostile_value = hostile_vector.subs(z, hostile_point)
require(hostile_value[0] == 0
        and hostile_value[1] != 0 and hostile_value[2] != 0,
        "selected chart wall killed full reciprocal wedge")
require(all(sp.det(reciprocal_transfer(index, 5)).subs(z, hostile_point) != 0
            for index in (1, 2, 3)),
        "hostile used a singular local transfer")

P_AB = wedge_chart((A[4], A[3]), (B[5], B[4]))
H_ROOT = sp.Rational(-109, 24)
require(sp.factor(P_H.subs(p, H_ROOT)) == 0
        and sp.factor(P_AB.subs(p, H_ROOT)) != 0,
        "H chart wall lost every pair chart")

EXACT_OBJECTS = (tuple(J_n), tuple(exterior), tuple(residue),
                 tuple(A.values()), tuple(B.values()), tuple(R.values()),
                 tuple(S.values()), (P_H, P_J, P_K, U_H, U_J, U_K),
                 tuple(hostile_vector))
require(not any(expression.has(sp.Float)
                for group in EXACT_OBJECTS for expression in group),
        "floating atom entered exact evidence")


print("THM-3192 RECIPROCAL COEFFICIENT-JET / PRS PLUECKER RETURN")
print("dependency_hash_checks=" + repr(len(DEPENDENCIES)))
print("reciprocal_gauge=J_n=D_(n+1)S_n(1/z)D_n^-1")
print("z_smith_J=(1,z,z)")
print("z_smith_Lambda2J=(z,z,z^2)")
print("normalized_exterior_residue_rank=2")
print("first_layer_flag_word_checks=" + repr(FIRST_LAYER_WORD_CHECKS))
print("first_layer_single_gate_checks="
      + repr(FIRST_LAYER_CONVOLUTION_CHECKS))
print("first_layer_diagonal_ratio=2n/(2n+1)")
print("hidden_e12_return_layer=z^2_to_-4d^2n(n+1)e02")
print("reciprocal_coefficient_checks=" + repr(JET_COEFFICIENT_CHECKS))
print("homogeneous_truncation_checks=" + repr(TRUNCATION_CHECKS))
print("direct_prs_top_coefficient_checks="
      + repr(DIRECT_PRS_COEFFICIENT_CHECKS))
print("offset6_chart_map=(P2(A,B)=U_H*R3,P2(R,A)=U_J*S2,P1(S,R)=U_K*K)")
print("offset6_unit_range=p>=197")
print("offset6_symbolic_p_unit_checks=" + repr(len(UNIT_EXPRESSIONS)))
print("reciprocal_path_hostile=(d=5,z=105/4,e01=0,e02e12_nonzero)")
print("H_algebraic_chart_hostile=(p=-109/24,P2=0,P1_nonzero)")
print("sympy_float_atoms=0")
print("scope=exact_jet_carrier_and_HJK_chart_map_not_arbitrary_offset_or_GMC2")
print("all_exact_checks=PASS")
