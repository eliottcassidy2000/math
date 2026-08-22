#!/usr/bin/env python3
"""Independent FLINT/Fourier audit of the Keller degree-243 gate.

This script imports only the definition prefix of the already audited
degree-81 FLINT regular-matrix engine.  It changes the field to F_251,
builds a fourth inverse stage, and recovers the fifth eliminant from all 250
nonzero field values by multiplicative Fourier inversion.  It neither imports
nor reads the tuple/interpolation candidate.
"""

from __future__ import annotations

import hashlib
import sys
import types
from pathlib import Path

from flint import nmod_poly


ROOT = Path(__file__).resolve().parents[1]
SUPPORT_PATH = ROOT / "04-computation/jc_level4_degree81_fourier_flint_independent_audit_20260816.py"
SUPPORT_SHA256 = "3d0bee9dd97993160fc7275cb4a96e77893013c421029eada4bbd5b46ac5d3e6"
SUPPORT_SENTINEL = "\nK0 = RegularAlgebra()\n"
P = 251
TARGET = (1, 1, 1)
EXPECTED_LEDGER_SHA256 = "912f32ec0b9b375d9db2ba71d7fdf224456c86862871ae0f8f92bac5038f00ab"
EXPECTED_SEMANTIC_SHA256 = "0f44c226329aedf3a3c232dc44fcd228cfce9587483d92f47090baac2f2be7ea"


def require(condition: bool, payload: object) -> None:
    if not condition:
        raise RuntimeError(payload)


support_raw = SUPPORT_PATH.read_bytes()
support_hash = hashlib.sha256(support_raw.replace(b"\r\n", b"\n")).hexdigest()
require(support_hash == SUPPORT_SHA256, ("degree-81 FLINT support drift", support_hash))
support_text = support_raw.decode("utf-8").replace("\r\n", "\n")
require(support_text.count(SUPPORT_SENTINEL) == 1, "FLINT definition boundary changed")
support_prefix = support_text.split(SUPPORT_SENTINEL, 1)[0]
support_module = types.ModuleType("thm3498_degree81_flint_support")
support_module.__file__ = str(SUPPORT_PATH)
sys.modules[support_module.__name__] = support_module
exec(compile(support_prefix, str(SUPPORT_PATH), "exec"), support_module.__dict__)

# The inherited matrix engine reads P and TARGET from its defining module.
support_module.P = P
support_module.TARGET = TARGET

RegularAlgebra = support_module.RegularAlgebra
sub = support_module.sub
fmap = support_module.fmap
l_value = support_module.l_value
inverse_coordinates = support_module.inverse_coordinates
make_extension = support_module.make_extension


K0 = RegularAlgebra()
t0 = tuple(K0.scalar(value) for value in TARGET)
rings = [K0]
targets = [t0]
inverse_points = []
for level in range(1, 5):
    current = rings[-1]
    target = targets[-1]
    extension = make_extension(current, *target, f"K{level}")
    embedded_target = tuple(extension.embed(value) for value in target)
    source = inverse_coordinates(
        extension, *embedded_target, extension.theta, f"K{level}"
    )
    require(fmap(extension, *source) == embedded_target, ("inverse graph", level))
    rings.append(extension)
    targets.append(source)
    inverse_points.append(source)

require(tuple(ring.dim for ring in rings[1:]) == (3, 9, 27, 81), "tower dimensions")
K4 = rings[-1]
q4 = inverse_points[-1]
L4 = l_value(K4, *q4)
T4 = sub(K4, K4.scalar(4), K4.mul(K4.scalar(3), K4.mul(q4[1], q4[2])))
C4 = K4.neg(K4.mul(K4.scalar(2), q4[2]))


def norm_at(value: int) -> int:
    element = K4.add(
        K4.mul(L4, K4.scalar(pow(value, 3, P))),
        K4.add(K4.mul(T4, K4.scalar(value)), C4),
    )
    return K4.norm(element)


# Since degree <= 3*81=243 < |F_251^*|=250, multiplicative Fourier
# inversion recovers every coefficient without aliasing.  As 1/250=-1 mod
# 251, c_k=-sum_x P(x)x^{-k}.
nonzero_values = {value: norm_at(value) for value in range(1, P)}
coefficients_0_to_249 = []
for exponent in range(P - 1):
    character_sum = sum(
        polynomial_value * pow(value, (-exponent) % (P - 1), P)
        for value, polynomial_value in nonzero_values.items()
    ) % P
    coefficients_0_to_249.append((-character_sum) % P)

require(
    all(value == 0 for value in coefficients_0_to_249[244:]),
    "Fourier tail above degree 243 is nonzero",
)
coefficients = coefficients_0_to_249[:244]
require(coefficients[-1] != 0, "degree-243 leading coefficient vanished")
require(coefficients[0] == norm_at(0), "unused X=0 determinant control failed")

polynomial = nmod_poly(coefficients, P)
require(polynomial.degree() == 243, "FLINT polynomial degree changed")
require(
    polynomial.gcd(polynomial.derivative()) == nmod_poly([1], P),
    "degree-243 polynomial is not squarefree",
)
for value, expected in nonzero_values.items():
    require(int(polynomial(value)) % P == expected, ("Horner mismatch", value))

_factor_unit, factors = polynomial.factor()
factor_degrees = tuple(sorted((factor.degree(), exponent) for factor, exponent in factors))
require(sum(degree * exponent for degree, exponent in factor_degrees) == 243, "factor degree sum")
require(all(exponent == 1 for _degree, exponent in factor_degrees), "factor repeated")

ledger = "\n".join(f"{index}:{coefficient}" for index, coefficient in enumerate(coefficients))
ledger_sha256 = hashlib.sha256(ledger.encode("ascii")).hexdigest()
require(ledger_sha256 == EXPECTED_LEDGER_SHA256, "degree-243 ledger changed")
semantic_lines = [
    f"support={support_hash}",
    f"prime={P};target={TARGET}",
    f"degree={polynomial.degree()};tail={tuple(coefficients_0_to_249[244:])}",
    f"factors={factor_degrees}",
    f"ledger={ledger_sha256}",
]
semantic_sha256 = hashlib.sha256("\n".join(semantic_lines).encode("ascii")).hexdigest()
require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256, "semantic ledger changed")

# At X=0, changing -2z to +2z negates the determinant in odd dimension 81.
wrong_sign_zero = K4.norm(K4.neg(C4))
require(wrong_sign_zero == (-norm_at(0)) % P, "constant-sign determinant parity")
require(wrong_sign_zero != norm_at(0), "constant-sign hostile did not fire")

print("== independent F_251 level-five Fourier/FLINT audit ==")
print(f"support_sha256={support_hash}")
print(f"target={TARGET}; regular_matrix_dimensions={(3, 9, 27, 81)}")
print("inverse graphs and all leading/derivative/chart units: PASS")
print("250-point multiplicative Fourier inversion: degrees 244..249 vanish; degree=243")
print("unused X=0 determinant and all-point Horner controls: PASS")
print(f"factor_degree_exponent_ledger={factor_degrees}")
print(f"constant_sign_zero_hostile={wrong_sign_zero}")
print("FLINT derivative gcd=1; the 243-root norm-product is squarefree")
print(f"ascending_coefficient_sha256={ledger_sha256}")
print(f"semantic_sha256={semantic_sha256}")
print("scope: good-reduction generic fifth-eliminant gate only; no R5 image prime, irreducibility, all-level, or general JC claim")
print("all independent exact checks passed")
