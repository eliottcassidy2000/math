#!/usr/bin/env python3
"""Exact companion for the THM-3365 Vieta/Vandermonde FC(3) benchmark."""

from __future__ import annotations

import ast
import hashlib
import itertools
from fractions import Fraction
from math import comb, factorial
from pathlib import Path

import sympy as sp


EXPECTED_SEMANTIC_SHA256 = "5e4b4741681568dc3359d9441c3cbb8022241fba938f1f6356f7d4eade33e18f"


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def multiply_sparse(left, right):
    product = {}
    for alpha, coefficient_left in left.items():
        for beta, coefficient_right in right.items():
            exponent = tuple(alpha[index] + beta[index] for index in range(3))
            product[exponent] = (
                product.get(exponent, 0) + coefficient_left * coefficient_right
            )
    return {exponent: coefficient for exponent, coefficient in product.items()
            if coefficient}


def factorial_functional(polynomial) -> int:
    return sum(
        coefficient
        * factorial(exponent[0])
        * factorial(exponent[1])
        * factorial(exponent[2])
        for exponent, coefficient in polynomial.items()
    )


def closed_even_moment(index: int) -> int:
    return (
        factorial(2 * index) ** 2
        * factorial(3 * index)
        // factorial(index)
    )


def chamber_sum(index: int) -> Fraction:
    return sum(
        Fraction(
            comb(2 * index, split)
            * factorial(2 * index + split)
            * factorial(4 * index - split),
            2 ** (2 * index + split),
        )
        for split in range(2 * index + 1)
    )


def odd_product(index: int) -> int:
    product = 1
    for value in range(1, index + 1):
        product *= 2 * value - 1
    return product


def beta_integral_by_expansion(index: int) -> Fraction:
    return sum(
        Fraction((-1) ** split * comb(2 * index, split),
                 2 * index + 2 * split + 1)
        for split in range(2 * index + 1)
    )


def beta_integral_by_half_gamma(index: int) -> Fraction:
    return Fraction(
        2 ** (4 * index + 1)
        * factorial(2 * index) ** 2
        * factorial(3 * index + 1),
        factorial(index) * factorial(6 * index + 2),
    )


def permutation_sign(permutation) -> int:
    inversions = sum(
        permutation[left] > permutation[right]
        for left in range(3)
        for right in range(left + 1, 3)
    )
    return -1 if inversions % 2 else 1


source = Path(__file__).read_text(encoding="utf-8")
syntax = ast.parse(source)
assertion_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(syntax))
floating_literals = sum(
    isinstance(node, ast.Constant) and isinstance(node.value, float)
    for node in ast.walk(syntax)
)
require(assertion_nodes == 0, "companion contains a Python assertion node")
require(floating_literals == 0, "companion contains a floating literal")


# The nonlinear transfer is the Vieta map E=(e1,e2,e3).
x, y, z = sp.symbols("x y z")
variables = (x, y, z)
e1 = x + y + z
e2 = x * y + y * z + z * x
e3 = x * y * z
vieta = sp.Matrix((e1, e2, e3))
jacobian = vieta.jacobian(variables)
determinant = sp.expand(jacobian.det())
vandermonde = sp.expand((x - y) * (y - z) * (z - x))
require(sp.expand(determinant + vandermonde) == 0,
        "Vieta Jacobian determinant has the wrong Vandermonde sign")
require(sp.factor(determinant) == (x - y) * (x - z) * (y - z),
        "Vieta determinant did not factor into the three difference forms")


# A sparse implementation independent of SymPy evaluates L(V^n) directly.
vandermonde_sparse = {
    (2, 1, 0): -1,
    (2, 0, 1): 1,
    (1, 2, 0): 1,
    (1, 0, 2): -1,
    (0, 2, 1): -1,
    (0, 1, 2): 1,
}
require(
    sp.Poly(vandermonde, variables).terms()
    == sorted(vandermonde_sparse.items(), reverse=True),
    "sparse and symbolic Vandermonde expansions disagree",
)

# Exact toric audit of the ordered six-column exponent configuration.
toric_columns = (
    (2, 1, 0),
    (2, 0, 1),
    (1, 2, 0),
    (1, 0, 2),
    (0, 2, 1),
    (0, 1, 2),
)
coefficient_vector = (-1, 1, 1, -1, -1, 1)
require(tuple(vandermonde_sparse) == toric_columns,
        "toric column order disagrees with the sparse Vandermonde")
require(tuple(vandermonde_sparse[column] for column in toric_columns)
        == coefficient_vector,
        "toric coefficient vector disagrees with the Vandermonde signs")

toric_matrix = sp.Matrix([
    [column[coordinate] for column in toric_columns]
    for coordinate in range(3)
])
coefficient_column = sp.Matrix(coefficient_vector)
require(toric_matrix.rank() == 3, "A_V does not have rank three")
require(toric_matrix * coefficient_column == sp.zeros(3, 1),
        "c_V is not in the integer kernel of A_V")

# Every conformal subvector is a subset of the three positive and three
# negative unit coordinates of c_V. Exhausting that cube proves Graver
# indecomposability without importing a toric-algebra implementation.
positive_positions = tuple(
    position for position, value in enumerate(coefficient_vector) if value > 0
)
negative_positions = tuple(
    position for position, value in enumerate(coefficient_vector) if value < 0
)
proper_conformal_checks = 0
proper_conformal_relations = []
for positive_mask in range(1 << len(positive_positions)):
    for negative_mask in range(1 << len(negative_positions)):
        relation = [0] * 6
        for bit, position in enumerate(positive_positions):
            if positive_mask & (1 << bit):
                relation[position] = 1
        for bit, position in enumerate(negative_positions):
            if negative_mask & (1 << bit):
                relation[position] = -1
        relation = tuple(relation)
        if relation == (0, 0, 0, 0, 0, 0):
            continue
        if relation == coefficient_vector:
            continue
        proper_conformal_checks += 1
        if toric_matrix * sp.Matrix(relation) == sp.zeros(3, 1):
            proper_conformal_relations.append(relation)

require(proper_conformal_checks == 62,
        "conformal-subvector universe has the wrong size")
require(not proper_conformal_relations,
        "c_V has a proper nonzero conformal subrelation")

# A smaller support relation proves that this Graver vector is not a circuit.
four_column_relation = (1, -1, 0, 0, -1, 1)
require(toric_matrix * sp.Matrix(four_column_relation) == sp.zeros(3, 1),
        "displayed four-column vector is not a relation")
require(sum(value != 0 for value in four_column_relation) == 4,
        "displayed non-circuit witness does not have support four")
require(
    tuple(
        toric_columns[0][coordinate] + toric_columns[5][coordinate]
        for coordinate in range(3)
    )
    == tuple(
        toric_columns[1][coordinate] + toric_columns[4][coordinate]
        for coordinate in range(3)
    )
    == (2, 2, 2),
    "210+012=201+021 relation failed",
)

gale_rows = (
    (2, 2, 3),
    (-2, -1, -2),
    (-1, -2, -2),
    (1, 0, 0),
    (0, 1, 0),
    (0, 0, 1),
)
gale_matrix = sp.Matrix(gale_rows)
require(gale_matrix.rank() == 3, "displayed Gale matrix lacks full rank")
require(toric_matrix * gale_matrix == sp.zeros(3, 3),
        "displayed Gale matrix is not killed by A_V")
require(all(any(entry != 0 for entry in row) for row in gale_rows),
        "displayed Gale matrix has a zero row")
proportional_gale_pairs = []
for left in range(6):
    for right in range(left + 1, 6):
        if sp.Matrix((gale_rows[left], gale_rows[right])).rank() < 2:
            proportional_gale_pairs.append((left + 1, right + 1))
require(not proportional_gale_pairs,
        "displayed Gale matrix has a proportional row pair")
bouquet_partition = tuple((index,) for index in range(1, 7))

all_plus_sparse = {column: 1 for column in toric_columns}
signed_first_moment = factorial_functional(vandermonde_sparse)
all_plus_first_moment = factorial_functional(all_plus_sparse)
require(signed_first_moment == 0 and all_plus_first_moment == 12,
        "same-support coefficient hostile has the wrong factorial moments")

power = {(0, 0, 0): 1}
direct_moments = []
direct_rows = []
for exponent in range(13):
    value = factorial_functional(power)
    expected = 0 if exponent % 2 else closed_even_moment(exponent // 2)
    require(value == expected,
            f"direct factorial moment mismatch at exponent {exponent}")
    direct_moments.append(value)
    direct_rows.append((exponent, len(power), value))
    power = multiply_sparse(power, vandermonde_sparse)


# Independent ordered-chamber, Beta-integral, and duplication controls.
chamber_rows = []
beta_rows = []
simplex_rows = []
for index in range(65):
    closed = closed_even_moment(index)
    chamber = chamber_sum(index)
    require(chamber.denominator == 1 and chamber.numerator == closed,
            f"ordered-chamber sum mismatch at m={index}")

    beta_expanded = beta_integral_by_expansion(index)
    beta_half_gamma = beta_integral_by_half_gamma(index)
    require(beta_expanded == beta_half_gamma,
            f"Beta integral mismatch at m={index}")

    # Rational forms of the two half-integer Gamma duplication identities.
    require(
        factorial(2 * index) * 2 ** index
        == 4 ** index * factorial(index) * odd_product(index),
        f"first duplication cross-product failed at m={index}",
    )
    triple_index = 3 * index + 1
    require(
        factorial(2 * triple_index)
        == 2 ** triple_index
        * factorial(triple_index)
        * odd_product(triple_index),
        f"second duplication cross-product failed at m={index}",
    )

    chamber_from_beta = (
        Fraction(factorial(6 * index + 1), 2 ** (4 * index))
        * beta_expanded
    )
    require(chamber_from_beta.denominator == 1,
            f"Beta chamber value is nonintegral at m={index}")
    require(chamber_from_beta.numerator == closed,
            f"Beta chamber formula mismatch at m={index}")

    simplex = Fraction(closed, factorial(6 * index + 2))
    normalized_simplex = 2 * simplex
    require(factorial(6 * index + 2) * simplex == closed,
            f"simplex polar normalization failed at m={index}")

    chamber_rows.append((index, closed))
    beta_rows.append((index, beta_expanded.numerator,
                      beta_expanded.denominator))
    simplex_rows.append((index, simplex.numerator, simplex.denominator,
                         normalized_simplex.numerator,
                         normalized_simplex.denominator))

require(simplex_rows[0] == (0, 1, 2, 1, 1),
        "simplex area normalization at m=0 is wrong")
require(simplex_rows[1] == (1, 1, 1680, 1, 840),
        "simplex Vandermonde-square normalization at m=1 is wrong")


# Alternation under all S3 variable permutations is the odd-moment mechanism.
permutation_rows = []
for permutation in itertools.permutations(range(3)):
    image = sp.expand(vandermonde.xreplace(
        {variables[index]: variables[permutation[index]] for index in range(3)}
    ))
    sign = permutation_sign(permutation)
    require(sp.expand(image - sign * vandermonde) == 0,
            f"alternation failed for permutation {permutation}")
    permutation_rows.append((permutation, sign))


# Exact scalar conventions: positivity is meaningful only for real scalars.
scalar_rows = []
for scalar in (sp.Integer(-3), sp.Integer(2), sp.I):
    moments = []
    for exponent in range(7):
        base = 0 if exponent % 2 else closed_even_moment(exponent // 2)
        moments.append(sp.expand(scalar ** exponent * base))
    if scalar.is_real:
        require(all(moments[2 * index] > 0 for index in range(1, 4)),
                f"real nonzero scalar lost even-moment positivity: {scalar}")
        scalar_type = "real-positive-even-moments"
    else:
        require(moments[2] == -closed_even_moment(1),
                "complex scalar control failed to reverse the second moment")
        require(all(moments[2 * index] != 0 for index in range(1, 4)),
                "complex nonzero scalar acquired a vanishing even moment")
        scalar_type = "complex-nonzero-not-ordered"
    scalar_rows.append((str(scalar), scalar_type,
                        tuple(str(value) for value in moments)))


# Hostile: three linear factors do not imply alternation or odd cancellation.
xyz_moments = tuple(factorial(exponent) ** 3 for exponent in range(1, 13))
require(all(value > 0 for value in xyz_moments),
        "xyz hostile unexpectedly acquired a vanishing moment")
require(xyz_moments[0] == 1 and xyz_moments[1] == 8,
        "xyz hostile normalization failed")


semantic_payload = (
    str(determinant),
    str(vandermonde),
    toric_columns,
    coefficient_vector,
    proper_conformal_checks,
    tuple(proper_conformal_relations),
    four_column_relation,
    gale_rows,
    tuple(proportional_gale_pairs),
    bouquet_partition,
    (signed_first_moment, all_plus_first_moment),
    tuple(direct_rows),
    tuple(chamber_rows),
    tuple(beta_rows),
    tuple(simplex_rows),
    tuple(permutation_rows),
    tuple(scalar_rows),
    xyz_moments,
)
semantic_sha256 = hashlib.sha256(
    repr(semantic_payload).encode("utf-8")
).hexdigest()
require(
    semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
    f"semantic digest mismatch: {semantic_sha256}",
)


print("THM-3365 FACTORIAL VANDERMONDE ALL-POWER BENCHMARK")
print("assertion_nodes", assertion_nodes, "floating_literals", floating_literals)
print("Vieta_map", (str(e1), str(e2), str(e3)))
print("Jacobian_determinant", str(determinant))
print("factorization", "-(x-y)(y-z)(z-x)")
print("TORIC BOUQUET AUDIT")
print("A_V_columns", toric_columns, "rank", toric_matrix.rank())
print("c_V", coefficient_vector, "A_V_c_V", (0, 0, 0))
print("proper_conformal_checks", proper_conformal_checks,
      "proper_conformal_relations", tuple(proper_conformal_relations),
      "Graver", True)
print("four_column_relation", four_column_relation,
      "support", 4, "c_V_is_circuit", False)
print("Gale_rows", gale_rows, "rank", gale_matrix.rank(),
      "A_V_times_Gale_zero", True)
print("proportional_Gale_pairs", tuple(proportional_gale_pairs),
      "bouquet_partition", bouquet_partition)
print("same_support_L_signed_all_plus",
      (signed_first_moment, all_plus_first_moment))
print("direct_exponents", (0, 12), "direct_moments", tuple(direct_moments))
print("chamber_beta_closed_checks", len(chamber_rows), "m_range", (0, 64))
print("even_formula", "(2m)!^2(3m)!/m!", "odd_formula", 0)
print("simplex_m0_unnormalized_normalized", ((1, 2), (1, 1)))
print("simplex_m1_unnormalized_normalized", ((1, 1680), (1, 840)))
print("permutation_signs", tuple(permutation_rows))
print("scalar_controls", tuple(scalar_rows))
print("xyz_hostile_n1_to_n12", xyz_moments)
print("semantic_sha256", semantic_sha256)
print("status=ALL EXACT CHECKS PASSED; alternating determinant benchmark only, not FC(3)")
