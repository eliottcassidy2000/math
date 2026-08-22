#!/usr/bin/env python3
"""Independent regular-matrix/Fourier audit of the Keller degree-729 gate.

This script imports only the definition prefix of the audited degree-81
regular-representation engine.  It builds five inverse stages over F_733,
represents every algebra element by its literal FLINT multiplication matrix,
and evaluates the sixth cubic-core norm at all 732 nonzero field elements.
Multiplicative Fourier inversion then reconstructs the complete polynomial.

It neither imports nor reads the recursive-tuple degree-729 candidate.  The
matching coefficient digest therefore compares two disjoint representations:
symbolic transitive cubic norms versus literal 243 by 243 determinants and
Fourier inversion.  The result is only a good-reduction degree/separability
gate, not an image, irreducibility, all-level, arbitrary-map, or general-JC
statement.
"""

from __future__ import annotations

import hashlib
import sys
import types
from pathlib import Path

from flint import nmod_poly


ROOT = Path(__file__).resolve().parents[1]
SUPPORT_PATH = (
    ROOT / "04-computation/jc_level4_degree81_fourier_flint_independent_audit_20260816.py"
)
SUPPORT_SHA256 = "3d0bee9dd97993160fc7275cb4a96e77893013c421029eada4bbd5b46ac5d3e6"
SUPPORT_SENTINEL = "\nK0 = RegularAlgebra()\n"
P = 733
TARGET = (1, 1, 1)
EXPECTED_COEFFICIENT_SHA256 = (
    "7aba23e306b00b14b8c60c34f9762ba8b35aecac111065058dfe9d4b3f1ecd51"
)
EXPECTED_SEMANTIC_SHA256 = (
    "d3eedf7368e7b98681ca41b76529c49403a351f79f2eb858e69df95362dc3518"
)


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


support_raw = SUPPORT_PATH.read_bytes()
support_hash = hashlib.sha256(support_raw.replace(b"\r\n", b"\n")).hexdigest()
require(support_hash == SUPPORT_SHA256, ("degree-81 FLINT support drift", support_hash))
support_text = support_raw.decode("utf-8").replace("\r\n", "\n")
require(support_text.count(SUPPORT_SENTINEL) == 1, "FLINT definition boundary changed")
support_module = types.ModuleType("level_six_degree729_flint_support")
support_module.__file__ = str(SUPPORT_PATH)
sys.modules[support_module.__name__] = support_module
exec(
    compile(
        support_text.split(SUPPORT_SENTINEL, 1)[0],
        str(SUPPORT_PATH),
        "exec",
    ),
    support_module.__dict__,
)

# The inherited engine reads the modulus dynamically from its defining module.
support_module.P = P
support_module.TARGET = TARGET

RegularAlgebra = support_module.RegularAlgebra
sub = support_module.sub
fmap = support_module.fmap
l_value = support_module.l_value
inverse_coordinates = support_module.inverse_coordinates
make_extension = support_module.make_extension


def determinant_nonzero(value) -> int:
    return int(value.det()) % P


# Build a literal regular-representation tower.  Every inverse operation in
# make_extension/inverse_coordinates checks a nonzero determinant first.
K0 = RegularAlgebra()
rings = [K0]
targets = [tuple(K0.scalar(value) for value in TARGET)]
inverse_points = []
unit_ledger = []
for level in range(1, 6):
    current = rings[-1]
    target = targets[-1]
    leading = l_value(current, *target)
    leading_norm = determinant_nonzero(leading)
    require(leading_norm != 0, ("leading L gate", level))
    extension = make_extension(current, *target, f"K{level}")
    embedded_target = tuple(extension.embed(value) for value in target)
    source = inverse_coordinates(
        extension, *embedded_target, extension.theta, f"K{level}"
    )
    require(fmap(extension, *source) == embedded_target, ("inverse graph", level))

    theta_square = extension.theta * extension.theta
    derivative = sub(
        extension,
        extension.mul(extension.scalar(3), theta_square),
        extension.embed(
            current.neg(
                current.mul(
                    sub(
                        current,
                        current.scalar(4),
                        current.mul(current.scalar(3), current.mul(target[1], target[2])),
                    ),
                    current.inverse(leading, f"K{level} audit leading L"),
                )
            )
        ),
    )
    derivative_norm = determinant_nonzero(derivative)
    require(derivative_norm != 0, ("derivative gate", level))

    rings.append(extension)
    targets.append(source)
    inverse_points.append(source)
    unit_ledger.append((level, extension.dim, leading_norm, derivative_norm))

dimensions = tuple(ring.dim for ring in rings[1:])
require(dimensions == (3, 9, 27, 81, 243), "inverse tower dimensions changed")

K5 = rings[-1]
q5 = inverse_points[-1]
L5 = l_value(K5, *q5)
T5 = sub(K5, K5.scalar(4), K5.mul(K5.scalar(3), K5.mul(q5[1], q5[2])))
C5 = K5.neg(K5.mul(K5.scalar(2), q5[2]))
terminal_leading_norm = determinant_nonzero(L5)
require(terminal_leading_norm != 0, "sixth-core leading coefficient is not a unit")


def norm_at(value: int) -> int:
    element = K5.add(
        K5.mul(L5, K5.scalar(pow(value, 3, P))),
        K5.add(K5.mul(T5, K5.scalar(value)), C5),
    )
    return K5.norm(element)


# Since deg <= 3*243=729 < |F_733^*|=732, multiplicative Fourier inversion
# has no aliasing.  Since 1/732=-1 mod 733, the inverse factor is -1.
nonzero_values = tuple((value, norm_at(value)) for value in range(1, P))
coefficients_0_to_731 = []
for exponent in range(P - 1):
    character_sum = sum(
        polynomial_value * pow(value, (-exponent) % (P - 1), P)
        for value, polynomial_value in nonzero_values
    ) % P
    coefficients_0_to_731.append((-character_sum) % P)

fourier_tail = tuple(coefficients_0_to_731[730:])
require(fourier_tail == (0, 0), ("Fourier tail above degree 729", fourier_tail))
coefficients = coefficients_0_to_731[:730]
require(coefficients[-1] == terminal_leading_norm, "degree-729 leading norm mismatch")
constant_control = norm_at(0)
require(coefficients[0] == constant_control, "unused X=0 determinant control failed")

polynomial = nmod_poly(coefficients, P)
require(polynomial.degree() == 729, "degree-729 leading coefficient vanished")
derivative_gcd = polynomial.gcd(polynomial.derivative())
require(derivative_gcd == nmod_poly([1], P), "degree-729 polynomial is not squarefree")
for value, expected in nonzero_values:
    require(int(polynomial(value)) % P == expected, ("Horner mismatch", value))

_factor_unit, factors = polynomial.factor()
factor_degrees = tuple(sorted((factor.degree(), exponent) for factor, exponent in factors))
require(sum(degree * exponent for degree, exponent in factor_degrees) == 729, "factor sum")
require(all(exponent == 1 for _degree, exponent in factor_degrees), "repeated factor")

# Odd dimension forces the sign-flipped constant core to negate its norm.
wrong_sign_zero = K5.norm(K5.neg(C5))
require(wrong_sign_zero == (-constant_control) % P, "constant-sign parity hostile")
require(wrong_sign_zero != constant_control, "constant-sign hostile did not fire")

coefficient_ledger = "\n".join(
    f"{index}:{coefficient}" for index, coefficient in enumerate(coefficients)
)
coefficient_sha256 = hashlib.sha256(coefficient_ledger.encode("ascii")).hexdigest()
require(
    coefficient_sha256 == EXPECTED_COEFFICIENT_SHA256,
    ("recursive-tuple/Fourier coefficient disagreement", coefficient_sha256),
)

sample_values = tuple(nonzero_values[index - 1] for index in (1, 2, 729, 730, 731, 732))
semantic_lines = [
    f"support={support_hash}",
    f"prime={P};target={TARGET};dimensions={dimensions}",
    f"units={tuple(unit_ledger)};terminal_L={terminal_leading_norm}",
    f"degree={polynomial.degree()};tail={fourier_tail};gcd={derivative_gcd.degree()}",
    f"factors={factor_degrees}",
    f"samples={sample_values};zero={constant_control};wrong_zero={wrong_sign_zero}",
    f"coefficients={coefficient_sha256}",
]
semantic_sha256 = hashlib.sha256("\n".join(semantic_lines).encode("ascii")).hexdigest()
if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
    require(
        semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
        ("independent audit semantic drift", semantic_sha256),
    )

print("== independent F_733 level-six Fourier/FLINT audit ==")
print(f"support_sha256={support_hash}")
print(f"target={TARGET}; regular_matrix_dimensions={dimensions}")
print(f"unit_gate_columns=(level,dimension,Norm(L),Norm(cubic_derivative));ledger={tuple(unit_ledger)}")
print("five inverse graphs and all leading/derivative/chart units: PASS")
print(f"732-point multiplicative Fourier inversion: degrees 730..731={fourier_tail};degree={polynomial.degree()}")
print("unused X=0 determinant and all-point Horner controls: PASS")
print(f"sample_determinants={sample_values};constant_sign_zero_hostile={wrong_sign_zero}")
print(f"factor_degree_exponent_ledger={factor_degrees}")
print(f"ascending_coefficient_sha256={coefficient_sha256}")
print(f"semantic_sha256={semantic_sha256}")
print("independent verdict: full degree 729 and squarefree on the lawful F_733 fibre")
print("scope: fixed-map generic sixth-eliminant gate only; no R6/R7 image, irreducibility, all-level, arbitrary-map, or general JC claim")
print("all independent exact checks passed")
