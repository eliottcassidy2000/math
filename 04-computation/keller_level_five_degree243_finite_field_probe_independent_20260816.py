#!/usr/bin/env python3
"""Independent degree-243 separability probe for the fixed Keller fifth iterate.

This extends the proved THM-3498 degree-81 construction by one inverse cubic
stage.  At one good target over F_251 it builds the complete nested inverse
algebra of dimensions 3, 9, 27, and 81.  It then forms the norm of the fifth
inverse cubic core as a polynomial in a fresh x-coordinate.

The primary norm route contracts one cubic layer at a time.  Selected values
are checked against the determinant of the literal 81 by 81 multiplication
matrix.  A full degree-243 squarefree result is a good-reduction witness for
generic degree and separability only; it is not an image equation,
irreducibility theorem, or general Jacobian-conjecture claim.
"""

from __future__ import annotations

import hashlib
import sys
import types
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SUPPORT_PATH = ROOT / "04-computation/keller_level_four_degree81_finite_field_probe_20260816.py"
SUPPORT_SHA256 = "4039b4081c9f0d95b197d2e3a7581c66433382e53dac3b95fa2526c3a4ba4f2e"
SUPPORT_SENTINEL = "\nbase = PrimeField()\n"
EXPECTED_LEDGER_SHA256 = "912f32ec0b9b375d9db2ba71d7fdf224456c86862871ae0f8f92bac5038f00ab"
EXPECTED_SEMANTIC_SHA256 = "1397fb3a3173e8dfbe867f4fe4c2d527ef33da3c94fb60fcb227f8c60b9d15b7"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


support_raw = SUPPORT_PATH.read_bytes()
support_hash = hashlib.sha256(support_raw.replace(b"\r\n", b"\n")).hexdigest()
if SUPPORT_SHA256:
    require(support_hash == SUPPORT_SHA256, ("THM-3498 support hash drift", support_hash))
support_text = support_raw.decode("utf-8").replace("\r\n", "\n")
require(support_text.count(SUPPORT_SENTINEL) == 1, "degree-81 definition boundary changed")
support_prefix = support_text.split(SUPPORT_SENTINEL, 1)[0]
support_module = types.ModuleType("thm3498_degree81_support")
support_module.__file__ = str(SUPPORT_PATH)
sys.modules[support_module.__name__] = support_module
exec(compile(support_prefix, str(SUPPORT_PATH), "exec"), support_module.__dict__)

# Replace only the finite field and target after loading the audited generic
# engine.  Its methods read MODULUS dynamically from their defining module.
MODULUS = 251
TARGET = (1, 1, 1)
support_module.MODULUS = MODULUS
support_module.TARGET = TARGET

PrimeField = support_module.PrimeField
Cubic = support_module.Cubic
sub = support_module.sub
fmap = support_module.fmap
l_value = support_module.l_value
inverse_coordinates = support_module.inverse_coordinates
make_extension = support_module.make_extension
poly_evaluate = support_module.poly_evaluate
poly_derivative = support_module.poly_derivative
poly_gcd = support_module.poly_gcd
interpolate_consecutive = support_module.interpolate_consecutive


def determinant_three(base, matrix):
    """Determinant of a 3 by 3 matrix over an arbitrary commutative base."""

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


def norm_to_base(ring: Cubic, value):
    base = ring.base
    zero = base.scalar(0)
    one = base.scalar(1)
    basis = ((one, zero, zero), (zero, one, zero), (zero, zero, one))
    columns = [ring.mul(value, vector) for vector in basis]
    matrix = [[columns[column][row] for column in range(3)] for row in range(3)]
    return determinant_three(base, matrix)


def absolute_norm(ring, value) -> int:
    if isinstance(ring, PrimeField):
        return value % MODULUS
    return absolute_norm(ring.base, norm_to_base(ring, value))


def build_inverse_tower():
    base = PrimeField()
    targets = [tuple(base.scalar(value) for value in TARGET)]
    rings = [base]
    inverse_points = []
    for level in range(1, 5):
        current = rings[-1]
        target = targets[-1]
        extension = make_extension(current, *target, f"K{level}")
        embedded_target = tuple(extension.embed(value) for value in target)
        source = inverse_coordinates(extension, *embedded_target, extension.theta)
        require(fmap(extension, *source) == embedded_target, ("inverse graph", level))
        rings.append(extension)
        inverse_points.append(source)
        targets.append(source)
    return rings, inverse_points


rings, inverse_points = build_inverse_tower()
K4 = rings[-1]
q4 = inverse_points[-1]
require(tuple(ring.dim for ring in rings[1:]) == (3, 9, 27, 81), "tower dimensions")

L4 = l_value(K4, *q4)
T4 = sub(K4, K4.scalar(4), K4.mul(K4.scalar(3), K4.mul(q4[1], q4[2])))
C4 = K4.neg(K4.mul(K4.scalar(2), q4[2]))


def fifth_core_value(value: int):
    return K4.add(
        K4.mul(L4, K4.scalar(pow(value, 3, MODULUS))),
        K4.add(K4.mul(T4, K4.scalar(value)), C4),
    )


def norm_value(value: int) -> int:
    return absolute_norm(K4, fifth_core_value(value))


grid_values = [norm_value(value) for value in range(244)]
P5 = interpolate_consecutive(grid_values)
require(len(P5) - 1 == 243, "the fifth eliminant lost degree")
for value, expected in enumerate(grid_values):
    require(poly_evaluate(P5, value) == expected, ("interpolation grid", value))
for value in (244, 245, 249):
    require(poly_evaluate(P5, value) == norm_value(value), ("off-grid recursive norm", value))
require(poly_gcd(P5, poly_derivative(P5)) == [1], "degree-243 eliminant is not squarefree")

# With only 243 nodes, interpolation is forced into degree at most 242 and
# must fail at the omitted node because the true leading coefficient and
# 243! are both units modulo 251.
short_P5 = interpolate_consecutive(grid_values[:243])
require(len(short_P5) - 1 <= 242, "short interpolation hostile degree")
require(
    poly_evaluate(short_P5, 243) != grid_values[243],
    "243-node interpolation hostile did not fire",
)

# A separately represented literal regular determinant checks the transitive
# cubic contraction at one grid point and two held-out points.
flat_checks = []
for value in (0, 244, 249):
    recursive = norm_value(value)
    flat = K4.norm(fifth_core_value(value))
    require(flat == recursive, ("flat/recursive norm", value, flat, recursive))
    flat_checks.append((value, flat))

# The defining inverse core is L X^3+T X-2z.  Flipping the constant sign at
# X=0 negates an odd (81-dimensional) norm, so this is a deterministic sign
# hostile rather than a merely different random polynomial.
wrong_sign_zero = absolute_norm(K4, K4.neg(C4))
require(wrong_sign_zero == (-flat_checks[0][1]) % MODULUS, "constant-sign norm parity")
require(wrong_sign_zero != flat_checks[0][1], "constant-sign hostile did not fire")

ledger = "\n".join(f"{index}:{coefficient}" for index, coefficient in enumerate(P5))
ledger_sha256 = hashlib.sha256(ledger.encode("ascii")).hexdigest()
require(ledger_sha256 == EXPECTED_LEDGER_SHA256, "degree-243 coefficient ledger changed")
semantic_lines = [
    f"support={support_hash}",
    f"prime={MODULUS};target={TARGET}",
    f"degree={len(P5)-1};gcd_degree={len(poly_gcd(P5, poly_derivative(P5)))-1}",
    f"ledger={ledger_sha256}",
    "flat=" + ";".join(f"{value}:{result}" for value, result in flat_checks),
]
semantic_sha256 = hashlib.sha256("\n".join(semantic_lines).encode("ascii")).hexdigest()
require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256, "semantic ledger changed")

print("== fixed Keller level-five degree-243 finite-field probe ==")
print(f"support_sha256={support_hash}")
print(f"prime={MODULUS}; target={TARGET}; tower_dimensions={(3, 9, 27, 81)}")
print("four inverse graphs: PASS; all inverse-chart denominators and cubic derivatives are units")
print(f"fifth_core_norm=degree {len(P5)-1}; derivative_gcd_degree={len(poly_gcd(P5, poly_derivative(P5)))-1}")
print(f"flat_81x81_determinant_checks={tuple(flat_checks)}")
print(f"hostiles=(243_node_interpolation_fails_at_243=True,constant_sign_zero={wrong_sign_zero})")
print(f"ascending_coefficient_sha256={ledger_sha256}")
print(f"semantic_sha256={semantic_sha256}")
print("finite-field verdict: full degree 243 and squarefree")
print("scope: generic fifth x-eliminant degree/separability witness only; no image prime, irreducibility, R5 identity, all-level, or general JC claim")
print("all exact checks passed")
