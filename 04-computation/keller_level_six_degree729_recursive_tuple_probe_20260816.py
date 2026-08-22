#!/usr/bin/env python3
"""FINITE-EXACT candidate for the fixed Keller sixth eliminant.

At the target (1,1,1) over F_733, this script builds the complete five-stage
inverse tower of dimensions 3,9,27,81,243.  All inversions are performed by
recursive 3 by 3 adjugates in nested coefficient tuples; neither level-five
THM-3525 script is imported.  The sixth cubic core is then normed *as a
polynomial* by five successive polynomial-valued 3 by 3 determinants.

The resulting degree-729 polynomial is checked for squarefreeness with FLINT,
against independently evaluated scalar recursive norms, and once against the
literal flattened 243 by 243 regular determinant.  This is a good-reduction
generic degree/separability candidate only.  It is not an image equation,
irreducibility result, all-level law, arbitrary-map theorem, or general
Jacobian-conjecture claim.
"""

from __future__ import annotations

import hashlib
import sys
import types
from pathlib import Path

from flint import nmod_poly


ROOT = Path(__file__).resolve().parents[1]
SUPPORT_PATH = ROOT / "04-computation/keller_level_four_degree81_finite_field_probe_20260816.py"
SUPPORT_SHA256 = "4039b4081c9f0d95b197d2e3a7581c66433382e53dac3b95fa2526c3a4ba4f2e"
SUPPORT_SENTINEL = "\nbase = PrimeField()\n"
MODULUS = 733
TARGET = (1, 1, 1)
EXPECTED_COEFFICIENT_SHA256 = "7aba23e306b00b14b8c60c34f9762ba8b35aecac111065058dfe9d4b3f1ecd51"
EXPECTED_SEMANTIC_SHA256 = "8009c86f1c8f290829df2ba8332dc2b09929b08cbd55376f48a13acd8c2c427c"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


support_raw = SUPPORT_PATH.read_bytes()
support_hash = hashlib.sha256(support_raw.replace(b"\r\n", b"\n")).hexdigest()
require(support_hash == SUPPORT_SHA256, ("THM-3498 support hash drift", support_hash))
support_text = support_raw.decode("utf-8").replace("\r\n", "\n")
require(support_text.count(SUPPORT_SENTINEL) == 1, "degree-81 definition boundary changed")
support_module = types.ModuleType("level_six_recursive_tuple_support")
support_module.__file__ = str(SUPPORT_PATH)
sys.modules[support_module.__name__] = support_module
exec(
    compile(support_text.split(SUPPORT_SENTINEL, 1)[0], str(SUPPORT_PATH), "exec"),
    support_module.__dict__,
)

# The audited generic engine reads these globals dynamically.
support_module.MODULUS = MODULUS
support_module.TARGET = TARGET

PrimeField = support_module.PrimeField
Cubic = support_module.Cubic
sub = support_module.sub
fmap = support_module.fmap
l_value = support_module.l_value
poly_evaluate = support_module.poly_evaluate


def determinant_three(base, matrix):
    """Determinant of a 3 by 3 matrix over any commutative base algebra."""

    positive = base.add(
        base.mul(matrix[0][0], base.mul(matrix[1][1], matrix[2][2])),
        base.add(
            base.mul(matrix[0][1], base.mul(matrix[1][2], matrix[2][0])),
            base.mul(matrix[0][2], base.mul(matrix[1][0], matrix[2][1])),
        ),
    )
    negative = base.add(
        base.mul(matrix[0][2], base.mul(matrix[1][1], matrix[2][0])),
        base.add(
            base.mul(matrix[0][1], base.mul(matrix[1][0], matrix[2][2])),
            base.mul(matrix[0][0], base.mul(matrix[1][2], matrix[2][1])),
        ),
    )
    return base.add(positive, base.neg(negative))


def relative_norm(ring: Cubic, value):
    base = ring.base
    zero, one = base.scalar(0), base.scalar(1)
    basis = ((one, zero, zero), (zero, one, zero), (zero, zero, one))
    columns = [ring.mul(value, vector) for vector in basis]
    matrix = [[columns[column][row] for column in range(3)] for row in range(3)]
    return determinant_three(base, matrix)


def absolute_norm(ring, value) -> int:
    if isinstance(ring, PrimeField):
        return value % MODULUS
    return absolute_norm(ring.base, relative_norm(ring, value))


def recursive_inverse(ring, value):
    """Invert by transitive cubic adjugates, never by a flattened solve."""

    if isinstance(ring, PrimeField):
        return ring.inv(value)
    base = ring.base
    zero, one = base.scalar(0), base.scalar(1)
    basis = ((one, zero, zero), (zero, one, zero), (zero, zero, one))
    columns = [ring.mul(value, vector) for vector in basis]
    matrix = [[columns[column][row] for column in range(3)] for row in range(3)]
    norm = determinant_three(base, matrix)
    inverse_norm = recursive_inverse(base, norm)

    # Column zero of M(value)^(-1) is C_(0,i)/det(M), i=0,1,2.
    cofactors = (
        base.add(
            base.mul(matrix[1][1], matrix[2][2]),
            base.neg(base.mul(matrix[1][2], matrix[2][1])),
        ),
        base.neg(
            base.add(
                base.mul(matrix[1][0], matrix[2][2]),
                base.neg(base.mul(matrix[1][2], matrix[2][0])),
            )
        ),
        base.add(
            base.mul(matrix[1][0], matrix[2][1]),
            base.neg(base.mul(matrix[1][1], matrix[2][0])),
        ),
    )
    result = tuple(base.mul(cofactor, inverse_norm) for cofactor in cofactors)
    require(ring.mul(value, result) == ring.one, ("recursive inverse failed", ring.label))
    return result


def unit_norm(ring, value, label: str) -> int:
    norm = absolute_norm(ring, value)
    require(norm != 0, ("nonunit gate", label))
    inverse = recursive_inverse(ring, value)
    one = ring.scalar(1) if isinstance(ring, PrimeField) else ring.one
    require(ring.mul(value, inverse) == one, ("unit product", label))
    return norm


def make_extension_recursive(base, a, b, c, label: str):
    """Create one inverse cubic and return its leading/derivative unit norms."""

    leading = l_value(base, a, b, c)
    leading_norm = unit_norm(base, leading, f"{label}:L")
    trisection = sub(base, base.scalar(4), base.mul(base.scalar(3), base.mul(b, c)))
    inverse_leading = recursive_inverse(base, leading)
    p = base.neg(base.mul(trisection, inverse_leading))
    q = base.mul(base.mul(base.scalar(2), c), inverse_leading)
    extension = Cubic(base, p, q, label)
    derivative = sub(
        extension,
        extension.mul(extension.scalar(3), extension.power(extension.theta, 2)),
        extension.embed(p),
    )
    derivative_norm = unit_norm(extension, derivative, f"{label}:derivative")
    return extension, leading_norm, derivative_norm


def inverse_coordinates_recursive(algebra, a, b, c, x, label: str):
    """Evaluate the inverse graph while exposing both chart denominators."""

    two, three, nine, twelve = map(algebra.scalar, (2, 3, 9, 12))
    denominator = algebra.add(
        algebra.add(
            algebra.mul(
                sub(algebra, algebra.mul(twelve, a), algebra.power(b, 2)),
                algebra.power(x, 2),
            ),
            algebra.mul(b, x),
        ),
        two,
    )
    denominator_norm = unit_norm(algebra, denominator, f"{label}:y-denominator")
    y_numerator = algebra.mul(
        algebra.mul(three, algebra.mul(a, x)),
        algebra.add(
            algebra.mul(sub(algebra, algebra.mul(nine, algebra.mul(a, c)), b), x),
            two,
        ),
    )
    y = sub(
        algebra,
        b,
        algebra.mul(y_numerator, recursive_inverse(algebra, denominator)),
    )
    x_cube = algebra.power(x, 3)
    x_cube_norm = unit_norm(algebra, x_cube, f"{label}:x^3")
    z_numerator = sub(
        algebra,
        sub(algebra, algebra.mul(two, x), c),
        algebra.mul(three, algebra.mul(algebra.power(x, 2), y)),
    )
    z = algebra.mul(z_numerator, recursive_inverse(algebra, x_cube))
    return (x, y, z), denominator_norm, x_cube_norm


def poly_trim(ring, polynomial):
    zero = ring.scalar(0)
    result = list(polynomial)
    while len(result) > 1 and result[-1] == zero:
        result.pop()
    return result


def poly_add(ring, left, right):
    zero = ring.scalar(0)
    result = [zero for _ in range(max(len(left), len(right)))]
    for index in range(len(result)):
        a = left[index] if index < len(left) else zero
        b = right[index] if index < len(right) else zero
        result[index] = ring.add(a, b)
    return poly_trim(ring, result)


def poly_neg(ring, polynomial):
    return poly_trim(ring, [ring.neg(value) for value in polynomial])


def poly_mul(ring, left, right):
    zero = ring.scalar(0)
    result = [zero for _ in range(len(left) + len(right) - 1)]
    for i, a in enumerate(left):
        if a == zero:
            continue
        for j, b in enumerate(right):
            if b != zero:
                result[i + j] = ring.add(result[i + j], ring.mul(a, b))
    return poly_trim(ring, result)


def polynomial_determinant_three(ring, matrix):
    positive = poly_add(
        ring,
        poly_mul(ring, matrix[0][0], poly_mul(ring, matrix[1][1], matrix[2][2])),
        poly_add(
            ring,
            poly_mul(ring, matrix[0][1], poly_mul(ring, matrix[1][2], matrix[2][0])),
            poly_mul(ring, matrix[0][2], poly_mul(ring, matrix[1][0], matrix[2][1])),
        ),
    )
    negative = poly_add(
        ring,
        poly_mul(ring, matrix[0][2], poly_mul(ring, matrix[1][1], matrix[2][0])),
        poly_add(
            ring,
            poly_mul(ring, matrix[0][1], poly_mul(ring, matrix[1][0], matrix[2][2])),
            poly_mul(ring, matrix[0][0], poly_mul(ring, matrix[1][2], matrix[2][1])),
        ),
    )
    return poly_add(ring, positive, poly_neg(ring, negative))


def relative_polynomial_norm(ring: Cubic, polynomial):
    """Norm a polynomial-valued element across one cubic tuple layer."""

    base = ring.base
    zero, one = base.scalar(0), base.scalar(1)
    basis = ((one, zero, zero), (zero, one, zero), (zero, zero, one))
    matrix = [[None for _ in range(3)] for _ in range(3)]
    for column, vector in enumerate(basis):
        products = [ring.mul(coefficient, vector) for coefficient in polynomial]
        for row in range(3):
            matrix[row][column] = poly_trim(base, [value[row] for value in products])
    return polynomial_determinant_three(base, matrix)


# Build the lawful inverse tower and freeze every gate before taking a norm.
base = PrimeField()
rings = [base]
targets = [tuple(base.scalar(value) for value in TARGET)]
inverse_points = []
unit_gate_ledger = []
for level in range(1, 6):
    current = rings[-1]
    target = targets[-1]
    extension, leading_norm, derivative_norm = make_extension_recursive(
        current, *target, f"K{level}"
    )
    embedded_target = tuple(extension.embed(value) for value in target)
    source, denominator_norm, x_cube_norm = inverse_coordinates_recursive(
        extension, *embedded_target, extension.theta, f"K{level}"
    )
    require(fmap(extension, *source) == embedded_target, ("inverse graph", level))
    rings.append(extension)
    targets.append(source)
    inverse_points.append(source)
    unit_gate_ledger.append(
        (level, extension.dim, leading_norm, derivative_norm, denominator_norm, x_cube_norm)
    )

tower_dimensions = tuple(ring.dim for ring in rings[1:])
require(tower_dimensions == (3, 9, 27, 81, 243), "inverse tower dimensions changed")

K5 = rings[-1]
q5 = inverse_points[-1]
L5 = l_value(K5, *q5)
terminal_leading_norm = unit_norm(K5, L5, "sixth-core:L5")
T5 = sub(K5, K5.scalar(4), K5.mul(K5.scalar(3), K5.mul(q5[1], q5[2])))
C5 = K5.neg(K5.mul(K5.scalar(2), q5[2]))


def sixth_core_value(value: int):
    return K5.add(
        K5.mul(L5, K5.scalar(pow(value, 3, MODULUS))),
        K5.add(K5.mul(T5, K5.scalar(value)), C5),
    )


# The primary route never interpolates values: it norms the cubic polynomial
# itself through polynomial-valued 3 by 3 determinants.
current_ring = K5
current_polynomial = [C5, T5, K5.scalar(0), L5]
degree_ledger = [len(current_polynomial) - 1]
while not isinstance(current_ring, PrimeField):
    current_polynomial = relative_polynomial_norm(current_ring, current_polynomial)
    current_ring = current_ring.base
    degree_ledger.append(len(current_polynomial) - 1)

P6 = [coefficient % MODULUS for coefficient in current_polynomial]
require(tuple(degree_ledger) == (3, 9, 27, 81, 243, 729), "relative norm degree ledger")
require(P6[-1] == terminal_leading_norm, "top coefficient is not Norm(L5)")
require(P6[0] == absolute_norm(K5, C5), "constant coefficient is not Norm(C5)")

flint_polynomial = nmod_poly(P6, MODULUS)
require(flint_polynomial.degree() == 729, "sixth eliminant lost degree")
derivative_gcd = flint_polynomial.gcd(flint_polynomial.derivative())
require(derivative_gcd == nmod_poly([1], MODULUS), "sixth eliminant is not squarefree")
_factor_unit, factors = flint_polynomial.factor()
factor_degree_exponents = tuple(
    sorted((factor.degree(), exponent) for factor, exponent in factors)
)
require(
    sum(degree * exponent for degree, exponent in factor_degree_exponents) == 729,
    "factor degrees do not sum to 729",
)
require(
    all(exponent == 1 for _degree, exponent in factor_degree_exponents),
    "factorization contains a repeated factor",
)

# Held-out scalar norms commute with the polynomial-valued norm.  The last
# three points lie beyond the 0..729 additive interpolation window.
check_points = (0, 1, 2, 729, 730, 731, 732)
recursive_value_checks = tuple(
    (value, absolute_norm(K5, sixth_core_value(value))) for value in check_points
)
for value, expected in recursive_value_checks:
    require(poly_evaluate(P6, value) == expected, ("polynomial/scalar norm", value))

# One deliberately expensive hostile forgets every cubic layer and computes
# the literal 243 by 243 regular determinant at a held-out point.
flat_check_point = 730
flat_determinant = K5.norm(sixth_core_value(flat_check_point))
require(
    flat_determinant == dict(recursive_value_checks)[flat_check_point],
    "flattened 243 by 243 determinant disagrees",
)

# Odd rank 243 forces the constant-sign flip to negate the X=0 norm.
wrong_sign_zero = absolute_norm(K5, K5.neg(C5))
require(wrong_sign_zero == (-P6[0]) % MODULUS, "constant-sign parity hostile")
require(wrong_sign_zero != P6[0], "constant-sign hostile did not fire")

# Deleting the certified degree-729 term must fail at a nonzero held-out point.
truncated_at_730 = poly_evaluate(P6[:-1], 730)
require(truncated_at_730 != poly_evaluate(P6, 730), "degree-729 truncation hostile")

coefficient_ledger = "\n".join(
    f"{index}:{coefficient}" for index, coefficient in enumerate(P6)
)
coefficient_sha256 = hashlib.sha256(coefficient_ledger.encode("ascii")).hexdigest()
if EXPECTED_COEFFICIENT_SHA256:
    require(coefficient_sha256 == EXPECTED_COEFFICIENT_SHA256, "coefficient ledger changed")

semantic_lines = [
    f"support={support_hash}",
    f"prime={MODULUS};target={TARGET};dimensions={tower_dimensions}",
    f"units={tuple(unit_gate_ledger)};terminal_L={terminal_leading_norm}",
    f"degrees={tuple(degree_ledger)};gcd_degree={derivative_gcd.degree()}",
    f"factors={factor_degree_exponents}",
    f"values={recursive_value_checks};flat=({flat_check_point},{flat_determinant})",
    f"wrong_sign_zero={wrong_sign_zero};truncated_at_730={truncated_at_730}",
    f"coefficients={coefficient_sha256}",
]
semantic_sha256 = hashlib.sha256("\n".join(semantic_lines).encode("ascii")).hexdigest()
if EXPECTED_SEMANTIC_SHA256:
    require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256, "semantic ledger changed")

print("== fixed Keller level-six degree-729 recursive-tuple probe ==")
print(f"support_sha256={support_hash}")
print(f"prime={MODULUS}; target={TARGET}; tower_dimensions={tower_dimensions}")
print("unit_gate_columns=(level,dimension,Norm(L),Norm(cubic_derivative),Norm(y_denominator),Norm(x^3))")
print(f"unit_gate_ledger={tuple(unit_gate_ledger)}")
print("five inverse graph identities: PASS; every listed chart/derivative gate is a unit")
print(f"sixth_core_terminal_L_norm={terminal_leading_norm}")
print(f"relative_polynomial_degree_ledger={tuple(degree_ledger)}")
print(f"sixth_core_norm_degree={flint_polynomial.degree()}; derivative_gcd_degree={derivative_gcd.degree()}")
print(f"recursive_scalar_checks={recursive_value_checks}")
print(f"flat_243x243_determinant_check=({flat_check_point},{flat_determinant})")
print(f"hostiles=(constant_sign_zero={wrong_sign_zero},truncated_at_730={truncated_at_730})")
print(f"factor_degree_exponent_ledger={factor_degree_exponents}")
print(f"ascending_coefficient_sha256={coefficient_sha256}")
print(f"semantic_sha256={semantic_sha256}")
print("FINITE-EXACT candidate verdict: full degree 729 and squarefree on one lawful good fibre")
print("scope: fixed-map generic sixth x-eliminant degree/separability candidate only; no R6/R7 image, irreducibility, all-level, arbitrary-map, or general JC claim")
print("all exact checks passed")
